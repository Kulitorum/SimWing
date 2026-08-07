#include "viewer_protocol.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <istream>
#include <limits>
#include <new>
#include <ostream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <unordered_set>
#include <utility>

namespace simwing::viewer {
namespace {

constexpr std::array<std::uint8_t, 4> kFrameMagic{'S', 'W', 'V', 'F'};
constexpr std::array<std::uint8_t, 4> kTraceMagic{'S', 'W', 'V', 'T'};
constexpr std::size_t kFrameEnvelopeBytes = 16;

void clearError(ProtocolError* error) {
    if (error != nullptr) {
        *error = {};
    }
}

bool fail(
    ProtocolError* error,
    ProtocolErrorCode code,
    std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

bool finite(double value) {
    return std::isfinite(value);
}

bool finite(const Vec3d& value) {
    return finite(value.x) && finite(value.y) && finite(value.z);
}

bool validString(
    const std::string& value,
    const char* label,
    bool requireNonempty,
    const ProtocolLimits& limits,
    ProtocolError* error) {
    if (requireNonempty && value.empty()) {
        return fail(error, ProtocolErrorCode::InvalidData,
                    std::string(label) + " must not be empty");
    }
    if (value.size() > limits.maxStringBytes) {
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    std::string(label) + " exceeds the string-size limit");
    }
    if (value.find('\0') != std::string::npos) {
        return fail(error, ProtocolErrorCode::InvalidData,
                    std::string(label) + " contains an embedded NUL");
    }
    return true;
}

bool countWithin(
    std::size_t count,
    std::uint32_t limit,
    const char* label,
    ProtocolError* error) {
    if (count > limit || count > std::numeric_limits<std::uint32_t>::max()) {
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    std::string(label) + " count exceeds the protocol limit");
    }
    return true;
}

std::size_t associatedCount(
    FieldAssociation association,
    const DiagnosticFrame& frame,
    bool* valid) {
    *valid = true;
    switch (association) {
    case FieldAssociation::Global:
        return 1;
    case FieldAssociation::Vertex:
        return frame.vertices.size();
    case FieldAssociation::Triangle:
        return frame.triangles.size();
    case FieldAssociation::Line:
        return frame.lines.size();
    case FieldAssociation::Contact:
        return frame.contacts.size();
    case FieldAssociation::Sealing:
        return frame.sealing.size();
    }
    *valid = false;
    return 0;
}

bool validEntityReference(
    const EntityReference& reference,
    const DiagnosticFrame& frame) {
    switch (reference.kind) {
    case EntityKind::Vertex:
        return reference.index < frame.vertices.size();
    case EntityKind::Triangle:
        return reference.index < frame.triangles.size();
    case EntityKind::Line:
        return reference.index < frame.lines.size();
    }
    return false;
}

template <class Element>
bool uniqueStableIds(
    const std::vector<Element>& elements,
    const char* label,
    ProtocolError* error) {
    std::unordered_set<std::uint64_t> ids;
    ids.reserve(elements.size());
    for (const Element& element : elements) {
        if (element.stableId == 0) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        std::string(label) + " stable ID must be nonzero");
        }
        if (!ids.insert(element.stableId).second) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        std::string(label) + " stable IDs are not unique");
        }
    }
    return true;
}

class BufferWriter {
public:
    BufferWriter(std::vector<std::uint8_t>& bytes, std::size_t maxBytes)
        : bytes_(bytes), maxBytes_(maxBytes) {}

    void u8(std::uint8_t value) {
        if (makeRoom(1)) {
            bytes_.push_back(value);
        }
    }

    void u16(std::uint16_t value) { integer(value); }
    void u32(std::uint32_t value) { integer(value); }
    void u64(std::uint64_t value) { integer(value); }

    void real(double value) {
        static_assert(sizeof(double) == sizeof(std::uint64_t));
        static_assert(std::numeric_limits<double>::is_iec559);
        u64(std::bit_cast<std::uint64_t>(value));
    }

    void vector(const Vec3d& value) {
        real(value.x);
        real(value.y);
        real(value.z);
    }

    void string(const std::string& value) {
        u32(static_cast<std::uint32_t>(value.size()));
        if (makeRoom(value.size())) {
            bytes_.insert(bytes_.end(), value.begin(), value.end());
        }
    }

    [[nodiscard]] bool exceeded() const noexcept { return exceeded_; }

private:
    template <class Integer>
    void integer(Integer value) {
        static_assert(std::is_unsigned_v<Integer>);
        if (makeRoom(sizeof(Integer))) {
            for (std::size_t i = 0; i < sizeof(Integer); ++i) {
                bytes_.push_back(static_cast<std::uint8_t>(value & 0xffU));
                value >>= 8U;
            }
        }
    }

    bool makeRoom(std::size_t count) {
        if (exceeded_ || count > maxBytes_ - bytes_.size()) {
            exceeded_ = true;
            return false;
        }
        return true;
    }

    std::vector<std::uint8_t>& bytes_;
    std::size_t maxBytes_ = 0;
    bool exceeded_ = false;
};

class BufferReader {
public:
    BufferReader(
        std::span<const std::uint8_t> bytes,
        ProtocolError* error,
        const ProtocolLimits& limits)
        : bytes_(bytes), error_(error), limits_(limits) {}

    [[nodiscard]] bool u8(std::uint8_t& value) {
        if (!need(1)) {
            return false;
        }
        value = bytes_[position_++];
        return true;
    }

    [[nodiscard]] bool u16(std::uint16_t& value) { return integer(value); }
    [[nodiscard]] bool u32(std::uint32_t& value) { return integer(value); }
    [[nodiscard]] bool u64(std::uint64_t& value) { return integer(value); }

    [[nodiscard]] bool real(double& value) {
        std::uint64_t bits = 0;
        if (!u64(bits)) {
            return false;
        }
        value = std::bit_cast<double>(bits);
        return true;
    }

    [[nodiscard]] bool vector(Vec3d& value) {
        return real(value.x) && real(value.y) && real(value.z);
    }

    [[nodiscard]] bool string(std::string& value, const char* label) {
        std::uint32_t size = 0;
        if (!u32(size)) {
            return false;
        }
        if (size > limits_.maxStringBytes) {
            return fail(error_, ProtocolErrorCode::LimitExceeded,
                        std::string(label) + " exceeds the string-size limit");
        }
        if (!need(size)) {
            return false;
        }
        value.assign(
            reinterpret_cast<const char*>(bytes_.data() + position_), size);
        position_ += size;
        return true;
    }

    [[nodiscard]] std::size_t remaining() const noexcept {
        return bytes_.size() - position_;
    }

private:
    [[nodiscard]] bool need(std::size_t count) {
        if (count > remaining()) {
            return fail(error_, ProtocolErrorCode::Truncated,
                        "diagnostic frame is truncated");
        }
        return true;
    }

    template <class Integer>
    [[nodiscard]] bool integer(Integer& value) {
        static_assert(std::is_unsigned_v<Integer>);
        if (!need(sizeof(Integer))) {
            return false;
        }
        std::uint64_t decoded = 0;
        for (std::size_t i = 0; i < sizeof(Integer); ++i) {
            decoded |= static_cast<std::uint64_t>(bytes_[position_++])
                       << (8U * i);
        }
        value = static_cast<Integer>(decoded);
        return true;
    }

    std::span<const std::uint8_t> bytes_;
    std::size_t position_ = 0;
    ProtocolError* error_ = nullptr;
    const ProtocolLimits& limits_;
};

void writeFramePayload(const DiagnosticFrame& frame, BufferWriter& writer) {
    writer.string(frame.sceneChecksum);
    writer.string(frame.solverCommit);
    writer.u64(frame.step);
    writer.real(frame.simulationTimeSeconds);
    writer.real(frame.timeStepSeconds);
    writer.u32(frame.couplingIteration);

    writer.real(frame.couplingResiduals.displacementMetres);
    writer.real(frame.couplingResiduals.tractionNewtons);
    writer.real(frame.couplingResiduals.fluid);
    writer.real(frame.couplingResiduals.structure);
    writer.real(frame.couplingResiduals.interfacePowerWatts);

    writer.real(frame.conservation.fluidMassKilograms);
    writer.vector(frame.conservation.totalMomentumNewtonSeconds);
    writer.real(frame.conservation.totalEnergyJoules);
    writer.vector(frame.conservation.interfaceForceResidualNewtons);
    writer.vector(frame.conservation.interfaceMomentResidualNewtonMetres);
    writer.real(frame.conservation.interfacePowerResidualWatts);

    writer.u32(static_cast<std::uint32_t>(frame.vertices.size()));
    for (const DiagnosticVertex& vertex : frame.vertices) {
        writer.u64(vertex.stableId);
        writer.vector(vertex.positionMetres);
    }

    writer.u32(static_cast<std::uint32_t>(frame.triangles.size()));
    for (const DiagnosticTriangle& triangle : frame.triangles) {
        writer.u64(triangle.stableId);
        writer.u32(triangle.vertex0);
        writer.u32(triangle.vertex1);
        writer.u32(triangle.vertex2);
        writer.u64(triangle.negativeRegionId);
        writer.u64(triangle.positiveRegionId);
    }

    writer.u32(static_cast<std::uint32_t>(frame.lines.size()));
    for (const DiagnosticLine& line : frame.lines) {
        writer.u64(line.stableId);
        writer.u32(line.vertex0);
        writer.u32(line.vertex1);
        writer.u32(line.role);
    }

    writer.u32(static_cast<std::uint32_t>(frame.contacts.size()));
    for (const ContactMarker& contact : frame.contacts) {
        writer.u64(contact.stableId);
        writer.u8(static_cast<std::uint8_t>(contact.first.kind));
        writer.u32(contact.first.index);
        writer.u8(static_cast<std::uint8_t>(contact.second.kind));
        writer.u32(contact.second.index);
        writer.vector(contact.positionMetres);
        writer.vector(contact.normal);
        writer.real(contact.signedGapMetres);
    }

    writer.u32(static_cast<std::uint32_t>(frame.sealing.size()));
    for (const SealingMarker& marker : frame.sealing) {
        writer.u64(marker.stableId);
        writer.u64(marker.negativeRegionId);
        writer.u64(marker.positiveRegionId);
        writer.vector(marker.positionMetres);
        writer.vector(marker.normal);
        writer.real(marker.apertureMetres);
        writer.u8(marker.sealed ? 1U : 0U);
    }

    writer.u32(static_cast<std::uint32_t>(frame.scalarFields.size()));
    for (const ScalarField& field : frame.scalarFields) {
        writer.string(field.name);
        writer.string(field.unit);
        writer.u8(static_cast<std::uint8_t>(field.association));
        writer.u32(static_cast<std::uint32_t>(field.values.size()));
        for (double value : field.values) {
            writer.real(value);
        }
    }

    writer.u32(static_cast<std::uint32_t>(frame.vectorFields.size()));
    for (const VectorField& field : frame.vectorFields) {
        writer.string(field.name);
        writer.string(field.unit);
        writer.u8(static_cast<std::uint8_t>(field.association));
        writer.u32(static_cast<std::uint32_t>(field.values.size()));
        for (const Vec3d& value : field.values) {
            writer.vector(value);
        }
    }
}

template <class Element, class ReadElement>
bool readElements(
    BufferReader& reader,
    std::vector<Element>& elements,
    std::uint32_t limit,
    const char* label,
    ProtocolError* error,
    ReadElement readElement) {
    std::uint32_t count = 0;
    if (!reader.u32(count)) {
        return false;
    }
    if (count > limit) {
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    std::string(label) + " count exceeds the protocol limit");
    }
    elements.resize(count);
    for (Element& element : elements) {
        if (!readElement(reader, element)) {
            return false;
        }
    }
    return true;
}

bool readEntityReference(BufferReader& reader, EntityReference& reference) {
    std::uint8_t kind = 0;
    if (!reader.u8(kind) || !reader.u32(reference.index)) {
        return false;
    }
    reference.kind = static_cast<EntityKind>(kind);
    return true;
}

bool readFramePayload(
    std::span<const std::uint8_t> payload,
    DiagnosticFrame& frame,
    ProtocolError* error,
    const ProtocolLimits& limits) {
    BufferReader reader(payload, error, limits);
    if (!reader.string(frame.sceneChecksum, "scene checksum")
        || !reader.string(frame.solverCommit, "solver commit")
        || !reader.u64(frame.step)
        || !reader.real(frame.simulationTimeSeconds)
        || !reader.real(frame.timeStepSeconds)
        || !reader.u32(frame.couplingIteration)
        || !reader.real(frame.couplingResiduals.displacementMetres)
        || !reader.real(frame.couplingResiduals.tractionNewtons)
        || !reader.real(frame.couplingResiduals.fluid)
        || !reader.real(frame.couplingResiduals.structure)
        || !reader.real(frame.couplingResiduals.interfacePowerWatts)
        || !reader.real(frame.conservation.fluidMassKilograms)
        || !reader.vector(frame.conservation.totalMomentumNewtonSeconds)
        || !reader.real(frame.conservation.totalEnergyJoules)
        || !reader.vector(frame.conservation.interfaceForceResidualNewtons)
        || !reader.vector(frame.conservation.interfaceMomentResidualNewtonMetres)
        || !reader.real(frame.conservation.interfacePowerResidualWatts)) {
        return false;
    }

    if (!readElements(reader, frame.vertices, limits.maxVertices, "vertex",
                      error,
                      [](BufferReader& input, DiagnosticVertex& vertex) {
                          return input.u64(vertex.stableId)
                                 && input.vector(vertex.positionMetres);
                      })
        || !readElements(reader, frame.triangles, limits.maxTriangles,
                         "triangle", error,
                         [](BufferReader& input,
                            DiagnosticTriangle& triangle) {
                             return input.u64(triangle.stableId)
                                    && input.u32(triangle.vertex0)
                                    && input.u32(triangle.vertex1)
                                    && input.u32(triangle.vertex2)
                                    && input.u64(triangle.negativeRegionId)
                                    && input.u64(triangle.positiveRegionId);
                         })
        || !readElements(reader, frame.lines, limits.maxLines, "line", error,
                         [](BufferReader& input, DiagnosticLine& line) {
                             return input.u64(line.stableId)
                                    && input.u32(line.vertex0)
                                    && input.u32(line.vertex1)
                                    && input.u32(line.role);
                         })
        || !readElements(reader, frame.contacts, limits.maxContactMarkers,
                         "contact marker", error,
                         [](BufferReader& input, ContactMarker& marker) {
                             return input.u64(marker.stableId)
                                    && readEntityReference(input, marker.first)
                                    && readEntityReference(input, marker.second)
                                    && input.vector(marker.positionMetres)
                                    && input.vector(marker.normal)
                                    && input.real(marker.signedGapMetres);
                         })
        || !readElements(reader, frame.sealing, limits.maxSealingMarkers,
                         "sealing marker", error,
                         [&](BufferReader& input, SealingMarker& marker) {
                             std::uint8_t sealed = 0;
                             if (!input.u64(marker.stableId)
                                 || !input.u64(marker.negativeRegionId)
                                 || !input.u64(marker.positiveRegionId)
                                 || !input.vector(marker.positionMetres)
                                 || !input.vector(marker.normal)
                                 || !input.real(marker.apertureMetres)
                                 || !input.u8(sealed)) {
                                 return false;
                             }
                             if (sealed > 1) {
                                 return fail(error,
                                             ProtocolErrorCode::InvalidData,
                                             "sealing marker has an invalid state");
                             }
                             marker.sealed = sealed != 0;
                             return true;
                         })) {
        return false;
    }

    if (!readElements(reader, frame.scalarFields, limits.maxScalarFields,
                      "scalar field", error,
                      [&](BufferReader& input, ScalarField& field) {
                          std::uint8_t association = 0;
                          std::uint32_t count = 0;
                          if (!input.string(field.name, "field name")
                              || !input.string(field.unit, "field unit")
                              || !input.u8(association)
                              || !input.u32(count)) {
                              return false;
                          }
                          if (count > limits.maxFieldValues) {
                              return fail(error,
                                          ProtocolErrorCode::LimitExceeded,
                                          "scalar field value count exceeds the protocol limit");
                          }
                          field.association =
                              static_cast<FieldAssociation>(association);
                          field.values.resize(count);
                          for (double& value : field.values) {
                              if (!input.real(value)) {
                                  return false;
                              }
                          }
                          return true;
                      })
        || !readElements(reader, frame.vectorFields, limits.maxVectorFields,
                         "vector field", error,
                         [&](BufferReader& input, VectorField& field) {
                             std::uint8_t association = 0;
                             std::uint32_t count = 0;
                             if (!input.string(field.name, "field name")
                                 || !input.string(field.unit, "field unit")
                                 || !input.u8(association)
                                 || !input.u32(count)) {
                                 return false;
                             }
                             if (count > limits.maxFieldValues) {
                                 return fail(error,
                                             ProtocolErrorCode::LimitExceeded,
                                             "vector field value count exceeds the protocol limit");
                             }
                             field.association =
                                 static_cast<FieldAssociation>(association);
                             field.values.resize(count);
                             for (Vec3d& value : field.values) {
                                 if (!input.vector(value)) {
                                     return false;
                                 }
                             }
                             return true;
                         })) {
        return false;
    }

    if (reader.remaining() != 0) {
        return fail(error, ProtocolErrorCode::TrailingData,
                    "diagnostic frame payload has trailing data");
    }
    return true;
}

template <class Integer>
void appendLittle(std::vector<std::uint8_t>& bytes, Integer value) {
    static_assert(std::is_unsigned_v<Integer>);
    for (std::size_t i = 0; i < sizeof(Integer); ++i) {
        bytes.push_back(static_cast<std::uint8_t>(value & 0xffU));
        value >>= 8U;
    }
}

template <class Integer>
Integer littleAt(std::span<const std::uint8_t> bytes, std::size_t offset) {
    static_assert(std::is_unsigned_v<Integer>);
    std::uint64_t decoded = 0;
    for (std::size_t i = 0; i < sizeof(Integer); ++i) {
        decoded |= static_cast<std::uint64_t>(bytes[offset + i]) << (8U * i);
    }
    return static_cast<Integer>(decoded);
}

template <class Integer>
void setLittleAt(
    std::vector<std::uint8_t>& bytes,
    std::size_t offset,
    Integer value) {
    static_assert(std::is_unsigned_v<Integer>);
    for (std::size_t i = 0; i < sizeof(Integer); ++i) {
        bytes[offset + i] = static_cast<std::uint8_t>(value & 0xffU);
        value >>= 8U;
    }
}

bool writeBytes(
    std::ostream& stream,
    const void* data,
    std::size_t size,
    ProtocolError* error) {
    stream.write(static_cast<const char*>(data),
                 static_cast<std::streamsize>(size));
    if (!stream) {
        return fail(error, ProtocolErrorCode::IoError,
                    "failed to write viewer trace");
    }
    return true;
}

template <class Integer>
bool writeLittle(
    std::ostream& stream,
    Integer value,
    ProtocolError* error) {
    std::array<std::uint8_t, sizeof(Integer)> bytes{};
    for (std::size_t i = 0; i < sizeof(Integer); ++i) {
        bytes[i] = static_cast<std::uint8_t>(value & 0xffU);
        value >>= 8U;
    }
    return writeBytes(stream, bytes.data(), bytes.size(), error);
}

bool readExact(
    std::istream& stream,
    void* data,
    std::size_t size,
    ProtocolError* error,
    const char* message) {
    stream.read(static_cast<char*>(data), static_cast<std::streamsize>(size));
    if (stream.gcount() != static_cast<std::streamsize>(size)) {
        return fail(error, ProtocolErrorCode::Truncated, message);
    }
    return true;
}

template <class Integer>
bool readLittle(
    std::istream& stream,
    Integer& value,
    ProtocolError* error,
    const char* message) {
    std::array<std::uint8_t, sizeof(Integer)> bytes{};
    if (!readExact(stream, bytes.data(), bytes.size(), error, message)) {
        return false;
    }
    value = littleAt<Integer>(bytes, 0);
    return true;
}

bool writeStreamString(
    std::ostream& stream,
    const std::string& value,
    ProtocolError* error) {
    return writeLittle(stream, static_cast<std::uint32_t>(value.size()), error)
           && writeBytes(stream, value.data(), value.size(), error);
}

bool readStreamString(
    std::istream& stream,
    std::string& value,
    const char* label,
    const ProtocolLimits& limits,
    ProtocolError* error) {
    std::uint32_t size = 0;
    if (!readLittle(stream, size, error, "viewer trace header is truncated")) {
        return false;
    }
    if (size > limits.maxStringBytes) {
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    std::string(label) + " exceeds the string-size limit");
    }
    value.resize(size);
    return readExact(stream, value.data(), size, error,
                     "viewer trace header is truncated");
}

bool validTraceHeader(
    const TraceHeader& header,
    const ProtocolLimits& limits,
    ProtocolError* error) {
    return validString(header.sceneChecksum, "scene checksum", true, limits,
                       error)
           && validString(header.solverCommit, "solver commit", true, limits,
                          error);
}

} // namespace

bool validateFrame(
    const DiagnosticFrame& frame,
    ProtocolError* error,
    const ProtocolLimits& limits) {
    clearError(error);
    if (!validString(frame.sceneChecksum, "scene checksum", true, limits,
                     error)
        || !validString(frame.solverCommit, "solver commit", true, limits,
                        error)
        || !countWithin(frame.vertices.size(), limits.maxVertices, "vertex",
                        error)
        || !countWithin(frame.triangles.size(), limits.maxTriangles,
                        "triangle", error)
        || !countWithin(frame.lines.size(), limits.maxLines, "line", error)
        || !countWithin(frame.contacts.size(), limits.maxContactMarkers,
                        "contact marker", error)
        || !countWithin(frame.sealing.size(), limits.maxSealingMarkers,
                        "sealing marker", error)
        || !countWithin(frame.scalarFields.size(), limits.maxScalarFields,
                        "scalar field", error)
        || !countWithin(frame.vectorFields.size(), limits.maxVectorFields,
                        "vector field", error)) {
        return false;
    }

    if (!finite(frame.simulationTimeSeconds)
        || frame.simulationTimeSeconds < 0.0
        || !finite(frame.timeStepSeconds) || frame.timeStepSeconds < 0.0) {
        return fail(error, ProtocolErrorCode::InvalidData,
                    "simulation time and time step must be finite and nonnegative");
    }

    const CouplingResiduals& residuals = frame.couplingResiduals;
    const ConservationValues& conservation = frame.conservation;
    if (!finite(residuals.displacementMetres)
        || !finite(residuals.tractionNewtons) || !finite(residuals.fluid)
        || !finite(residuals.structure)
        || !finite(residuals.interfacePowerWatts)
        || !finite(conservation.fluidMassKilograms)
        || conservation.fluidMassKilograms < 0.0
        || !finite(conservation.totalMomentumNewtonSeconds)
        || !finite(conservation.totalEnergyJoules)
        || !finite(conservation.interfaceForceResidualNewtons)
        || !finite(conservation.interfaceMomentResidualNewtonMetres)
        || !finite(conservation.interfacePowerResidualWatts)) {
        return fail(error, ProtocolErrorCode::InvalidData,
                    "residual and conservation values must be finite");
    }

    if (!uniqueStableIds(frame.vertices, "vertex", error)
        || !uniqueStableIds(frame.triangles, "triangle", error)
        || !uniqueStableIds(frame.lines, "line", error)
        || !uniqueStableIds(frame.contacts, "contact marker", error)
        || !uniqueStableIds(frame.sealing, "sealing marker", error)) {
        return false;
    }

    for (const DiagnosticVertex& vertex : frame.vertices) {
        if (!finite(vertex.positionMetres)) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        "vertex position is not finite");
        }
    }
    for (const DiagnosticTriangle& triangle : frame.triangles) {
        if (triangle.vertex0 >= frame.vertices.size()
            || triangle.vertex1 >= frame.vertices.size()
            || triangle.vertex2 >= frame.vertices.size()
            || triangle.vertex0 == triangle.vertex1
            || triangle.vertex1 == triangle.vertex2
            || triangle.vertex2 == triangle.vertex0
            || triangle.negativeRegionId == 0
            || triangle.positiveRegionId == 0
            || triangle.negativeRegionId == triangle.positiveRegionId) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        "triangle has invalid vertices or side-region IDs");
        }
    }
    for (const DiagnosticLine& line : frame.lines) {
        if (line.vertex0 >= frame.vertices.size()
            || line.vertex1 >= frame.vertices.size()
            || line.vertex0 == line.vertex1) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        "line has invalid vertex indices");
        }
    }
    for (const ContactMarker& marker : frame.contacts) {
        if (!validEntityReference(marker.first, frame)
            || !validEntityReference(marker.second, frame)
            || !finite(marker.positionMetres) || !finite(marker.normal)
            || !finite(marker.signedGapMetres)) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        "contact marker contains invalid data or entity indices");
        }
    }
    for (const SealingMarker& marker : frame.sealing) {
        if (!finite(marker.positionMetres) || !finite(marker.normal)
            || !finite(marker.apertureMetres) || marker.apertureMetres < 0.0
            || marker.negativeRegionId == 0
            || marker.positiveRegionId == 0
            || marker.negativeRegionId == marker.positiveRegionId) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        "sealing marker contains invalid or non-finite data");
        }
    }

    std::unordered_set<std::string> fieldNames;
    fieldNames.reserve(frame.scalarFields.size() + frame.vectorFields.size());
    const auto validateField = [&](const auto& field,
                                   const char* kind) -> bool {
        if (!validString(field.name, "field name", true, limits, error)
            || !validString(field.unit, "field unit", false, limits, error)
            || !countWithin(field.values.size(), limits.maxFieldValues,
                            "field value", error)) {
            return false;
        }
        bool validAssociation = false;
        const std::size_t expected =
            associatedCount(field.association, frame, &validAssociation);
        if (!validAssociation || field.values.size() != expected) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        std::string(kind)
                            + " field value count does not match its association");
        }
        if (!fieldNames.insert(field.name).second) {
            return fail(error, ProtocolErrorCode::InvalidData,
                        "diagnostic field names must be unique");
        }
        for (const auto& value : field.values) {
            if (!finite(value)) {
                return fail(error, ProtocolErrorCode::InvalidData,
                            std::string(kind)
                                + " field contains a non-finite value");
            }
        }
        return true;
    };
    for (const ScalarField& field : frame.scalarFields) {
        if (!validateField(field, "scalar")) {
            return false;
        }
    }
    for (const VectorField& field : frame.vectorFields) {
        if (!validateField(field, "vector")) {
            return false;
        }
    }
    return true;
}

bool serializeFrame(
    const DiagnosticFrame& frame,
    std::vector<std::uint8_t>& bytes,
    ProtocolError* error,
    const ProtocolLimits& limits) {
    clearError(error);
    bytes.clear();
    if (!validateFrame(frame, error, limits)) {
        return false;
    }

    try {
        if (limits.maxFrameBytes < kFrameEnvelopeBytes
            || limits.maxFrameBytes
                   > std::numeric_limits<std::size_t>::max()) {
            return fail(error, ProtocolErrorCode::LimitExceeded,
                        "configured diagnostic frame byte limit is invalid");
        }
        bytes.insert(bytes.end(), kFrameMagic.begin(), kFrameMagic.end());
        appendLittle(bytes, kFrameProtocolVersion);
        appendLittle(bytes, static_cast<std::uint16_t>(0));
        appendLittle(bytes, static_cast<std::uint64_t>(0));

        BufferWriter frameWriter(
            bytes, static_cast<std::size_t>(limits.maxFrameBytes));
        writeFramePayload(frame, frameWriter);
        if (frameWriter.exceeded()) {
            bytes.clear();
            return fail(error, ProtocolErrorCode::LimitExceeded,
                        "serialized diagnostic frame exceeds the byte limit");
        }
        const std::uint64_t payloadSize =
            static_cast<std::uint64_t>(bytes.size() - kFrameEnvelopeBytes);
        setLittleAt(bytes, 8, payloadSize);
        return true;
    } catch (const std::bad_alloc&) {
        bytes.clear();
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    "unable to allocate the diagnostic frame buffer");
    } catch (const std::length_error&) {
        bytes.clear();
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    "diagnostic frame size exceeds the platform limit");
    }
}

bool deserializeFrame(
    std::span<const std::uint8_t> bytes,
    DiagnosticFrame& frame,
    ProtocolError* error,
    const ProtocolLimits& limits) {
    clearError(error);
    if (bytes.size() < kFrameEnvelopeBytes) {
        return fail(error, ProtocolErrorCode::Truncated,
                    "diagnostic frame envelope is truncated");
    }
    if (!std::equal(kFrameMagic.begin(), kFrameMagic.end(), bytes.begin())) {
        return fail(error, ProtocolErrorCode::InvalidMagic,
                    "diagnostic frame magic is invalid");
    }
    const std::uint16_t version = littleAt<std::uint16_t>(bytes, 4);
    const std::uint16_t flags = littleAt<std::uint16_t>(bytes, 6);
    if (version != kFrameProtocolVersion) {
        return fail(error, ProtocolErrorCode::UnsupportedVersion,
                    "diagnostic frame version is unsupported");
    }
    if (flags != 0) {
        return fail(error, ProtocolErrorCode::UnsupportedVersion,
                    "diagnostic frame uses unsupported flags");
    }
    const std::uint64_t payloadSize = littleAt<std::uint64_t>(bytes, 8);
    if (limits.maxFrameBytes < kFrameEnvelopeBytes
        || payloadSize > limits.maxFrameBytes
        || payloadSize > std::numeric_limits<std::size_t>::max()
        || payloadSize > limits.maxFrameBytes - kFrameEnvelopeBytes) {
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    "diagnostic frame exceeds the byte limit");
    }
    const std::size_t expectedSize =
        kFrameEnvelopeBytes + static_cast<std::size_t>(payloadSize);
    if (bytes.size() < expectedSize) {
        return fail(error, ProtocolErrorCode::Truncated,
                    "diagnostic frame payload is truncated");
    }
    if (bytes.size() > expectedSize) {
        return fail(error, ProtocolErrorCode::TrailingData,
                    "diagnostic frame has trailing data");
    }

    try {
        DiagnosticFrame candidate;
        if (!readFramePayload(bytes.subspan(kFrameEnvelopeBytes), candidate,
                              error, limits)
            || !validateFrame(candidate, error, limits)) {
            return false;
        }
        frame = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    "unable to allocate diagnostic frame data");
    } catch (const std::length_error&) {
        return fail(error, ProtocolErrorCode::LimitExceeded,
                    "diagnostic frame count exceeds the platform limit");
    }
}

TraceWriter::TraceWriter(std::ostream& output, ProtocolLimits limits)
    : output_(&output), limits_(limits) {}

bool TraceWriter::writeHeader(const TraceHeader& header) {
    error_ = {};
    if (headerWritten_) {
        return fail(&error_, ProtocolErrorCode::InvalidData,
                    "viewer trace header was already written");
    }
    if (!validTraceHeader(header, limits_, &error_)
        || !writeBytes(*output_, kTraceMagic.data(), kTraceMagic.size(),
                       &error_)
        || !writeLittle(*output_, kTraceProtocolVersion, &error_)
        || !writeLittle(*output_, kFrameProtocolVersion, &error_)
        || !writeLittle(*output_, static_cast<std::uint32_t>(0), &error_)
        || !writeStreamString(*output_, header.sceneChecksum, &error_)
        || !writeStreamString(*output_, header.solverCommit, &error_)) {
        return false;
    }
    header_ = header;
    headerWritten_ = true;
    return true;
}

bool TraceWriter::writeFrame(const DiagnosticFrame& frame) {
    error_ = {};
    if (!headerWritten_) {
        return fail(&error_, ProtocolErrorCode::InvalidData,
                    "viewer trace header must be written before frames");
    }
    if (frame.sceneChecksum != header_.sceneChecksum
        || frame.solverCommit != header_.solverCommit) {
        return fail(&error_, ProtocolErrorCode::InvalidData,
                    "frame provenance does not match the viewer trace header");
    }
    std::vector<std::uint8_t> bytes;
    if (!serializeFrame(frame, bytes, &error_, limits_)
        || !writeLittle(*output_, static_cast<std::uint64_t>(bytes.size()),
                        &error_)
        || !writeBytes(*output_, bytes.data(), bytes.size(), &error_)) {
        return false;
    }
    return true;
}

const ProtocolError& TraceWriter::error() const noexcept {
    return error_;
}

TraceReader::TraceReader(std::istream& input, ProtocolLimits limits)
    : input_(&input), limits_(limits) {}

bool TraceReader::readHeader(TraceHeader& header) {
    error_ = {};
    header = {};
    if (headerRead_) {
        return fail(&error_, ProtocolErrorCode::InvalidData,
                    "viewer trace header was already read");
    }
    std::array<std::uint8_t, 4> magic{};
    std::uint16_t traceVersion = 0;
    std::uint16_t frameVersion = 0;
    std::uint32_t flags = 0;
    if (!readExact(*input_, magic.data(), magic.size(), &error_,
                   "viewer trace header is truncated")) {
        return false;
    }
    if (magic != kTraceMagic) {
        return fail(&error_, ProtocolErrorCode::InvalidMagic,
                    "viewer trace magic is invalid");
    }
    if (!readLittle(*input_, traceVersion, &error_,
                    "viewer trace header is truncated")
        || !readLittle(*input_, frameVersion, &error_,
                       "viewer trace header is truncated")
        || !readLittle(*input_, flags, &error_,
                       "viewer trace header is truncated")) {
        return false;
    }
    if (traceVersion != kTraceProtocolVersion
        || frameVersion != kFrameProtocolVersion || flags != 0) {
        return fail(&error_, ProtocolErrorCode::UnsupportedVersion,
                    "viewer trace version or flags are unsupported");
    }
    TraceHeader candidate;
    try {
        if (!readStreamString(*input_, candidate.sceneChecksum,
                              "scene checksum", limits_, &error_)
            || !readStreamString(*input_, candidate.solverCommit,
                                 "solver commit", limits_, &error_)
            || !validTraceHeader(candidate, limits_, &error_)) {
            return false;
        }
    } catch (const std::bad_alloc&) {
        return fail(&error_, ProtocolErrorCode::LimitExceeded,
                    "unable to allocate viewer trace header data");
    }
    header_ = candidate;
    header = std::move(candidate);
    headerRead_ = true;
    return true;
}

TraceReadStatus TraceReader::readNext(DiagnosticFrame& frame) {
    error_ = {};
    if (!headerRead_) {
        fail(&error_, ProtocolErrorCode::InvalidData,
             "viewer trace header must be read before frames");
        return TraceReadStatus::Error;
    }

    std::array<std::uint8_t, sizeof(std::uint64_t)> sizeBytes{};
    input_->read(reinterpret_cast<char*>(sizeBytes.data()),
                 static_cast<std::streamsize>(sizeBytes.size()));
    const std::streamsize received = input_->gcount();
    if (received == 0 && input_->eof()) {
        return TraceReadStatus::End;
    }
    if (received != static_cast<std::streamsize>(sizeBytes.size())) {
        fail(&error_, ProtocolErrorCode::Truncated,
             "viewer trace frame length is truncated");
        return TraceReadStatus::Error;
    }
    const std::uint64_t size = littleAt<std::uint64_t>(sizeBytes, 0);
    if (size > limits_.maxFrameBytes
        || size > std::numeric_limits<std::size_t>::max()) {
        fail(&error_, ProtocolErrorCode::LimitExceeded,
             "viewer trace frame exceeds the byte limit");
        return TraceReadStatus::Error;
    }

    try {
        std::vector<std::uint8_t> bytes(static_cast<std::size_t>(size));
        DiagnosticFrame decoded;
        if (!readExact(*input_, bytes.data(), bytes.size(), &error_,
                       "viewer trace frame is truncated")
            || !deserializeFrame(bytes, decoded, &error_, limits_)) {
            return TraceReadStatus::Error;
        }
        if (decoded.sceneChecksum != header_.sceneChecksum
            || decoded.solverCommit != header_.solverCommit) {
            fail(&error_, ProtocolErrorCode::InvalidData,
                 "frame provenance does not match the viewer trace header");
            return TraceReadStatus::Error;
        }
        frame = std::move(decoded);
    } catch (const std::bad_alloc&) {
        fail(&error_, ProtocolErrorCode::LimitExceeded,
             "unable to allocate viewer trace frame data");
        return TraceReadStatus::Error;
    }

    return TraceReadStatus::Frame;
}

const ProtocolError& TraceReader::error() const noexcept {
    return error_;
}

} // namespace simwing::viewer

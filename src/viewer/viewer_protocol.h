#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <span>
#include <string>
#include <vector>

namespace simwing::viewer {

inline constexpr std::uint16_t kFrameProtocolVersion = 1;
inline constexpr std::uint16_t kTraceProtocolVersion = 1;

// Limits are part of the untrusted-input boundary. Tests may lower them to
// exercise rejection paths without allocating production-sized buffers.
struct ProtocolLimits {
    std::uint64_t maxFrameBytes = 256ULL * 1024ULL * 1024ULL;
    std::uint32_t maxVertices = 5'000'000;
    std::uint32_t maxTriangles = 10'000'000;
    std::uint32_t maxLines = 10'000'000;
    std::uint32_t maxContactMarkers = 2'000'000;
    std::uint32_t maxSealingMarkers = 2'000'000;
    std::uint32_t maxScalarFields = 256;
    std::uint32_t maxVectorFields = 256;
    std::uint32_t maxFieldValues = 10'000'000;
    std::uint32_t maxStringBytes = 1024;
};

enum class ProtocolErrorCode {
    None,
    InvalidData,
    InvalidMagic,
    UnsupportedVersion,
    LimitExceeded,
    Truncated,
    TrailingData,
    IoError
};

struct ProtocolError {
    ProtocolErrorCode code = ProtocolErrorCode::None;
    std::string message;

    [[nodiscard]] explicit operator bool() const noexcept {
        return code != ProtocolErrorCode::None;
    }
};

struct Vec3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct DiagnosticVertex {
    std::uint64_t stableId = 0;
    Vec3d positionMetres;
};

struct DiagnosticTriangle {
    std::uint64_t stableId = 0;
    std::uint32_t vertex0 = 0;
    std::uint32_t vertex1 = 0;
    std::uint32_t vertex2 = 0;
    std::uint64_t negativeRegionId = 0;
    std::uint64_t positiveRegionId = 0;
};

struct DiagnosticLine {
    std::uint64_t stableId = 0;
    std::uint32_t vertex0 = 0;
    std::uint32_t vertex1 = 0;
    std::uint32_t role = 0;
};

enum class EntityKind : std::uint8_t {
    Vertex = 1,
    Triangle = 2,
    Line = 3
};

struct EntityReference {
    EntityKind kind = EntityKind::Vertex;
    std::uint32_t index = 0;
};

struct ContactMarker {
    std::uint64_t stableId = 0;
    EntityReference first;
    EntityReference second;
    Vec3d positionMetres;
    Vec3d normal;
    double signedGapMetres = 0.0;
};

struct SealingMarker {
    std::uint64_t stableId = 0;
    std::uint64_t negativeRegionId = 0;
    std::uint64_t positiveRegionId = 0;
    Vec3d positionMetres;
    Vec3d normal;
    double apertureMetres = 0.0;
    bool sealed = false;
};

enum class FieldAssociation : std::uint8_t {
    Global = 0,
    Vertex = 1,
    Triangle = 2,
    Line = 3,
    Contact = 4,
    Sealing = 5
};

struct ScalarField {
    std::string name;
    std::string unit;
    FieldAssociation association = FieldAssociation::Global;
    std::vector<double> values;
};

struct VectorField {
    std::string name;
    std::string unit;
    FieldAssociation association = FieldAssociation::Global;
    std::vector<Vec3d> values;
};

struct CouplingResiduals {
    double displacementMetres = 0.0;
    double tractionNewtons = 0.0;
    double fluid = 0.0;
    double structure = 0.0;
    double interfacePowerWatts = 0.0;
};

struct ConservationValues {
    double fluidMassKilograms = 0.0;
    Vec3d totalMomentumNewtonSeconds;
    double totalEnergyJoules = 0.0;
    Vec3d interfaceForceResidualNewtons;
    Vec3d interfaceMomentResidualNewtonMetres;
    double interfacePowerResidualWatts = 0.0;
};

// A frame owns all of its data. Publishing a const frame or a shared_ptr to a
// const frame therefore never aliases solver arrays or OpenGL buffers.
struct DiagnosticFrame {
    std::string sceneChecksum;
    std::string solverCommit;
    std::uint64_t step = 0;
    double simulationTimeSeconds = 0.0;
    double timeStepSeconds = 0.0;
    std::uint32_t couplingIteration = 0;
    CouplingResiduals couplingResiduals;
    ConservationValues conservation;
    std::vector<DiagnosticVertex> vertices;
    std::vector<DiagnosticTriangle> triangles;
    std::vector<DiagnosticLine> lines;
    std::vector<ContactMarker> contacts;
    std::vector<SealingMarker> sealing;
    std::vector<ScalarField> scalarFields;
    std::vector<VectorField> vectorFields;
};

struct TraceHeader {
    std::string sceneChecksum;
    std::string solverCommit;
};

[[nodiscard]] bool validateFrame(
    const DiagnosticFrame& frame,
    ProtocolError* error = nullptr,
    const ProtocolLimits& limits = ProtocolLimits{});

[[nodiscard]] bool serializeFrame(
    const DiagnosticFrame& frame,
    std::vector<std::uint8_t>& bytes,
    ProtocolError* error = nullptr,
    const ProtocolLimits& limits = ProtocolLimits{});

[[nodiscard]] bool deserializeFrame(
    std::span<const std::uint8_t> bytes,
    DiagnosticFrame& frame,
    ProtocolError* error = nullptr,
    const ProtocolLimits& limits = ProtocolLimits{});

class TraceWriter {
public:
    explicit TraceWriter(
        std::ostream& output,
        ProtocolLimits limits = ProtocolLimits{});

    [[nodiscard]] bool writeHeader(const TraceHeader& header);
    [[nodiscard]] bool writeFrame(const DiagnosticFrame& frame);
    // Writes an explicit successful end marker. Growing-file readers remain
    // pending at a temporary physical EOF and stop only when this marker is
    // visible. A writer cannot append frames after finishing.
    [[nodiscard]] bool finish();
    [[nodiscard]] const ProtocolError& error() const noexcept;

private:
    std::ostream* output_ = nullptr;
    ProtocolLimits limits_;
    TraceHeader header_;
    ProtocolError error_;
    bool headerWritten_ = false;
    bool finished_ = false;
};

enum class TraceReadStatus {
    Frame,
    // Follow-mode only: the physical file currently ends between records or
    // partway through a record. Calling readNext again resumes transactionally
    // without rereading or discarding the bytes already received.
    Pending,
    End,
    Error
};

enum class TraceReadMode {
    // A natural physical EOF between records is the end of a completed legacy
    // trace. A partial record remains an error.
    Replay,
    // A natural or partial EOF is temporary. Only an explicit end marker is
    // End; callers should poll at a bounded rate.
    Follow
};

class TraceReader {
public:
    explicit TraceReader(
        std::istream& input,
        ProtocolLimits limits = ProtocolLimits{});
    TraceReader(
        std::istream& input,
        TraceReadMode mode,
        ProtocolLimits limits = ProtocolLimits{});

    [[nodiscard]] bool readHeader(TraceHeader& header);
    [[nodiscard]] TraceReadStatus readNext(DiagnosticFrame& frame);
    [[nodiscard]] const ProtocolError& error() const noexcept;

private:
    std::istream* input_ = nullptr;
    ProtocolLimits limits_;
    TraceHeader header_;
    ProtocolError error_;
    bool headerRead_ = false;
    bool ended_ = false;
    TraceReadMode mode_ = TraceReadMode::Replay;
    std::uint8_t lengthBytes_[sizeof(std::uint64_t)]{};
    std::size_t lengthBytesRead_ = 0;
    std::uint64_t pendingFrameSize_ = 0;
    std::vector<std::uint8_t> pendingFrameBytes_;
    std::size_t pendingFrameBytesRead_ = 0;
};

} // namespace simwing::viewer

#include "scene.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <exception>
#include <istream>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::array<char, 8> binaryMagic{
    'S', 'I', 'M', 'W', 'S', 'C', 'N', '2'};
constexpr std::uint32_t maximumEntityCount = 10'000'000;
constexpr std::uint32_t maximumStringBytes = 1'048'576;
constexpr std::size_t maximumDecodedAllocationBytes =
    256ULL * 1024ULL * 1024ULL;

bool finite(double value) {
    return std::isfinite(value);
}

bool finite(const Vec2& value) {
    return finite(value.x) && finite(value.y);
}

bool finite(const Vec3& value) {
    return finite(value.x) && finite(value.y) && finite(value.z);
}

bool exactlyZero(const Vec3& value) {
    return value.x == 0.0 && value.y == 0.0 && value.z == 0.0;
}

bool finite(const Quaternion& value) {
    return finite(value.w) && finite(value.x) && finite(value.y)
           && finite(value.z);
}

Vec3 subtract(const Vec3& left, const Vec3& right) {
    return {left.x - right.x, left.y - right.y, left.z - right.z};
}

Vec3 cross(const Vec3& left, const Vec3& right) {
    return {left.y * right.z - left.z * right.y,
            left.z * right.x - left.x * right.z,
            left.x * right.y - left.y * right.x};
}

double squaredLength(const Vec3& value) {
    return value.x * value.x + value.y * value.y + value.z * value.z;
}

bool valid(RegionKind value) {
    return value == RegionKind::Outside || value == RegionKind::Cell;
}

bool valid(SurfaceRole value) {
    return value == SurfaceRole::Skin || value == SurfaceRole::Rib
           || value == SurfaceRole::Diagonal || value == SurfaceRole::MiniRib;
}

bool valid(OpeningRole value) {
    return value == OpeningRole::Intake || value == OpeningRole::Vent
           || value == OpeningRole::Crossport;
}

bool valid(AttachmentKind value) {
    return value == AttachmentKind::SurfaceVertex
           || value == AttachmentKind::PilotHarness
           || value == AttachmentKind::SuspensionJunction;
}

bool valid(SuspensionLineRole value) {
    return value == SuspensionLineRole::Suspension
           || value == SuspensionLineRole::Riser
           || value == SuspensionLineRole::Brake
           || value == SuspensionLineRole::Harness;
}

void add(ValidationReport& report,
         ValidationCode code,
         EntityKind entityKind,
         StableId entityId,
         std::string message) {
    report.diagnostics.push_back(
        {DiagnosticSeverity::Error, code, entityKind, entityId,
         std::move(message)});
}

template <typename T>
std::vector<const T*> sortedById(const std::vector<T>& values) {
    std::vector<const T*> result;
    result.reserve(values.size());
    for (const T& value : values) {
        result.push_back(&value);
    }
    std::stable_sort(result.begin(), result.end(), [](const T* left,
                                                      const T* right) {
        return left->id < right->id;
    });
    return result;
}

template <typename T>
void validateIds(const std::vector<T>& values,
                 EntityKind kind,
                 std::string_view label,
                 ValidationReport& report) {
    std::map<StableId, std::size_t> counts;
    for (const T& value : values) {
        if (value.id == invalidStableId) {
            add(report, ValidationCode::InvalidId, kind, value.id,
                std::string(label) + " ID must be nonzero");
        } else {
            ++counts[value.id];
        }
    }
    for (const auto& [id, count] : counts) {
        if (count > 1) {
            add(report, ValidationCode::DuplicateId, kind, id,
                std::string("duplicate ") + std::string(label) + " ID");
        }
    }
}

template <typename T>
std::unordered_set<StableId> idSet(const std::vector<T>& values) {
    std::unordered_set<StableId> result;
    result.reserve(values.size());
    for (const T& value : values) {
        if (value.id != invalidStableId) {
            result.insert(value.id);
        }
    }
    return result;
}

const Vertex* findVertex(const std::unordered_map<StableId, const Vertex*>& byId,
                         StableId id) {
    const auto found = byId.find(id);
    return found == byId.end() ? nullptr : found->second;
}

bool degenerate(const Triangle& triangle,
                const std::unordered_map<StableId, const Vertex*>& vertices) {
    const Vertex* a = findVertex(vertices, triangle.vertexIds[0]);
    const Vertex* b = findVertex(vertices, triangle.vertexIds[1]);
    const Vertex* c = findVertex(vertices, triangle.vertexIds[2]);
    if (!a || !b || !c || !finite(a->positionMeters)
        || !finite(b->positionMeters) || !finite(c->positionMeters)) {
        return false;
    }
    const Vec3 ab = subtract(b->positionMeters, a->positionMeters);
    const Vec3 ac = subtract(c->positionMeters, a->positionMeters);
    const Vec3 bc = subtract(c->positionMeters, b->positionMeters);
    const double longestSquared = std::max(
        {squaredLength(ab), squaredLength(ac), squaredLength(bc)});
    if (!(longestSquared > 0.0) || !finite(longestSquared)) {
        return true;
    }
    const double crossSquared = squaredLength(cross(ab, ac));
    // A relative test avoids baking a particular canopy scale into the
    // validator while rejecting coincident and effectively collinear faces.
    return !(crossSquared > 1.0e-24 * longestSquared * longestSquared);
}

std::string diagnosticSummary(const ValidationReport& report) {
    if (report.diagnostics.empty()) {
        return {};
    }
    const ValidationDiagnostic& first = report.diagnostics.front();
    std::ostringstream result;
    result << first.message;
    if (report.diagnostics.size() > 1) {
        result << " (and " << report.diagnostics.size() - 1
               << " more validation error(s))";
    }
    return result.str();
}

Scene canonicalScene(const Scene& scene) {
    Scene result = scene;
    const auto byId = [](const auto& left, const auto& right) {
        return left.id < right.id;
    };
    std::stable_sort(result.regions.begin(), result.regions.end(), byId);
    std::stable_sort(result.vertices.begin(), result.vertices.end(), byId);
    std::stable_sort(result.fabricMaterials.begin(),
                     result.fabricMaterials.end(), byId);
    std::stable_sort(result.seamMaterials.begin(),
                     result.seamMaterials.end(), byId);
    std::stable_sort(result.triangles.begin(), result.triangles.end(), byId);
    std::stable_sort(result.openings.begin(), result.openings.end(), byId);
    std::stable_sort(result.seams.begin(), result.seams.end(), byId);
    std::stable_sort(result.lineMaterials.begin(), result.lineMaterials.end(),
                     byId);
    std::stable_sort(result.pilots.begin(), result.pilots.end(), byId);
    std::stable_sort(result.suspensionJunctions.begin(),
                     result.suspensionJunctions.end(), byId);
    std::stable_sort(result.attachments.begin(), result.attachments.end(), byId);
    std::stable_sort(result.suspensionLines.begin(),
                     result.suspensionLines.end(), byId);
    return result;
}

bool withinBinarySafetyLimits(const Scene& scene, std::string& error) {
    std::size_t remaining = maximumDecodedAllocationBytes;
    const auto claim = [&](std::size_t count,
                           std::size_t elementSize,
                           std::string_view label) {
        if (count > maximumEntityCount
            || (elementSize != 0 && count > remaining / elementSize)) {
            error = std::string(label)
                + " exceeds the scene binary allocation limit";
            return false;
        }
        remaining -= count * elementSize;
        return true;
    };
    const auto claimString = [&](const std::string& value,
                                 std::string_view label) {
        if (value.size() > maximumStringBytes
            || value.size() > remaining) {
            error = std::string(label)
                + " exceeds the scene binary string/allocation limit";
            return false;
        }
        remaining -= value.size();
        return true;
    };

    if (!claim(scene.regions.size(), sizeof(FluidRegion), "region count")
        || !claim(scene.vertices.size(), sizeof(Vertex), "vertex count")
        || !claim(scene.fabricMaterials.size(), sizeof(FabricMaterial),
                  "fabric-material count")
        || !claim(scene.seamMaterials.size(), sizeof(SeamMaterial),
                  "seam-material count")
        || !claim(scene.triangles.size(), sizeof(Triangle), "triangle count")
        || !claim(scene.openings.size(), sizeof(Opening), "opening count")
        || !claim(scene.seams.size(), sizeof(Seam), "seam count")
        || !claim(scene.lineMaterials.size(), sizeof(LineMaterial),
                  "line-material count")
        || !claim(scene.pilots.size(), sizeof(Pilot), "pilot count")
        || !claim(scene.suspensionJunctions.size(),
                  sizeof(SuspensionJunction), "suspension-junction count")
        || !claim(scene.attachments.size(), sizeof(Attachment),
                  "attachment count")
        || !claim(scene.suspensionLines.size(), sizeof(SuspensionLine),
                  "suspension-line count")
        || !claimString(scene.metadata.designChecksum, "design checksum")
        || !claimString(scene.metadata.exporterVersion, "exporter version")) {
        return false;
    }
    for (const FluidRegion& region : scene.regions) {
        if (!claimString(region.name, "region name")) {
            return false;
        }
    }
    for (const FabricMaterial& material : scene.fabricMaterials) {
        if (!claimString(material.name, "fabric-material name")) {
            return false;
        }
    }
    for (const SeamMaterial& material : scene.seamMaterials) {
        if (!claimString(material.name, "seam-material name")) {
            return false;
        }
    }
    for (const Opening& opening : scene.openings) {
        if (!claim(opening.orderedVertexIds.size(), sizeof(StableId),
                   "opening boundary count")) {
            return false;
        }
    }
    for (const Seam& seam : scene.seams) {
        if (!claim(seam.firstOrderedVertexIds.size(), sizeof(StableId),
                   "first seam-chain vertex count")
            || !claim(seam.secondOrderedVertexIds.size(), sizeof(StableId),
                      "second seam-chain vertex count")) {
            return false;
        }
    }
    for (const LineMaterial& material : scene.lineMaterials) {
        if (!claimString(material.name, "line-material name")) {
            return false;
        }
    }
    for (const Pilot& pilot : scene.pilots) {
        if (!claimString(pilot.name, "pilot name")) {
            return false;
        }
    }
    return true;
}

class BinaryWriter {
public:
    explicit BinaryWriter(std::ostream& output) : output_(output) {}

    void bytes(const char* data, std::size_t count) {
        if (!ok_) {
            return;
        }
        output_.write(data, static_cast<std::streamsize>(count));
        if (!output_) {
            fail("failed to write scene stream");
        }
    }

    void u8(std::uint8_t value) {
        const char byte = static_cast<char>(value);
        bytes(&byte, 1);
    }

    void u32(std::uint32_t value) {
        std::array<char, 4> encoded{};
        for (std::size_t i = 0; i < encoded.size(); ++i) {
            encoded[i] = static_cast<char>((value >> (8 * i)) & 0xffu);
        }
        bytes(encoded.data(), encoded.size());
    }

    void u64(std::uint64_t value) {
        std::array<char, 8> encoded{};
        for (std::size_t i = 0; i < encoded.size(); ++i) {
            encoded[i] = static_cast<char>((value >> (8 * i)) & 0xffu);
        }
        bytes(encoded.data(), encoded.size());
    }

    void real(double value) {
        u64(std::bit_cast<std::uint64_t>(value));
    }

    void string(const std::string& value) {
        if (value.size() > maximumStringBytes) {
            fail("scene string exceeds the binary safety limit");
            return;
        }
        u32(static_cast<std::uint32_t>(value.size()));
        bytes(value.data(), value.size());
    }

    void count(std::size_t value) {
        if (value > maximumEntityCount) {
            fail("scene entity collection exceeds the binary safety limit");
            return;
        }
        u32(static_cast<std::uint32_t>(value));
    }

    [[nodiscard]] bool ok() const noexcept { return ok_; }
    [[nodiscard]] const std::string& error() const noexcept { return error_; }

private:
    void fail(std::string message) {
        if (ok_) {
            ok_ = false;
            error_ = std::move(message);
        }
    }

    std::ostream& output_;
    bool ok_ = true;
    std::string error_;
};

class BinaryReader {
public:
    explicit BinaryReader(std::istream& input) : input_(input) {}

    bool bytes(char* data, std::size_t count) {
        if (!ok_) {
            return false;
        }
        input_.read(data, static_cast<std::streamsize>(count));
        if (input_.gcount() != static_cast<std::streamsize>(count)) {
            fail("truncated scene stream");
            return false;
        }
        return true;
    }

    std::uint8_t u8() {
        char value = 0;
        bytes(&value, 1);
        return static_cast<std::uint8_t>(static_cast<unsigned char>(value));
    }

    std::uint32_t u32() {
        std::array<char, 4> encoded{};
        if (!bytes(encoded.data(), encoded.size())) {
            return 0;
        }
        std::uint32_t result = 0;
        for (std::size_t i = 0; i < encoded.size(); ++i) {
            result |= static_cast<std::uint32_t>(
                          static_cast<unsigned char>(encoded[i]))
                      << (8 * i);
        }
        return result;
    }

    std::uint64_t u64() {
        std::array<char, 8> encoded{};
        if (!bytes(encoded.data(), encoded.size())) {
            return 0;
        }
        std::uint64_t result = 0;
        for (std::size_t i = 0; i < encoded.size(); ++i) {
            result |= static_cast<std::uint64_t>(
                          static_cast<unsigned char>(encoded[i]))
                      << (8 * i);
        }
        return result;
    }

    double real() {
        return std::bit_cast<double>(u64());
    }

    std::string string() {
        const std::uint32_t size = u32();
        if (!ok_) {
            return {};
        }
        if (size > maximumStringBytes) {
            fail("scene string exceeds the binary safety limit");
            return {};
        }
        if (!claimAllocation(size, sizeof(char), "scene string")) {
            return {};
        }
        std::string result(size, '\0');
        if (size != 0) {
            bytes(result.data(), size);
        }
        return result;
    }

    std::uint32_t count(std::string_view label) {
        const std::uint32_t value = u32();
        if (ok_ && value > maximumEntityCount) {
            fail(std::string(label) + " count exceeds the binary safety limit");
            return 0;
        }
        return value;
    }

    void fail(std::string message) {
        if (ok_) {
            ok_ = false;
            error_ = std::move(message);
        }
    }

    [[nodiscard]] bool ok() const noexcept { return ok_; }
    [[nodiscard]] const std::string& error() const noexcept { return error_; }

    bool claimAllocation(std::size_t count,
                         std::size_t elementSize,
                         std::string_view label) {
        if (!ok_) {
            return false;
        }
        if (elementSize != 0
            && count > remainingAllocationBytes_ / elementSize) {
            fail(std::string(label)
                 + " exceeds the scene decoded-allocation limit");
            return false;
        }
        remainingAllocationBytes_ -= count * elementSize;
        return true;
    }

private:
    std::istream& input_;
    std::size_t remainingAllocationBytes_ = maximumDecodedAllocationBytes;
    bool ok_ = true;
    std::string error_;
};

void writeVec2(BinaryWriter& writer, const Vec2& value) {
    writer.real(value.x);
    writer.real(value.y);
}

void writeVec3(BinaryWriter& writer, const Vec3& value) {
    writer.real(value.x);
    writer.real(value.y);
    writer.real(value.z);
}

void writeQuaternion(BinaryWriter& writer, const Quaternion& value) {
    writer.real(value.w);
    writer.real(value.x);
    writer.real(value.y);
    writer.real(value.z);
}

Vec2 readVec2(BinaryReader& reader) {
    return {reader.real(), reader.real()};
}

Vec3 readVec3(BinaryReader& reader) {
    return {reader.real(), reader.real(), reader.real()};
}

Quaternion readQuaternion(BinaryReader& reader) {
    return {reader.real(), reader.real(), reader.real(), reader.real()};
}

template <typename T>
void resizeFromCount(std::vector<T>& values,
                     BinaryReader& reader,
                     std::string_view label) {
    const std::uint32_t count = reader.count(label);
    if (reader.claimAllocation(count, sizeof(T), label)) {
        values.resize(count);
    }
}

} // namespace

bool ValidationReport::ok() const noexcept {
    return errorCount() == 0;
}

std::size_t ValidationReport::errorCount() const noexcept {
    return static_cast<std::size_t>(std::count_if(
        diagnostics.begin(), diagnostics.end(),
        [](const ValidationDiagnostic& diagnostic) {
            return diagnostic.severity == DiagnosticSeverity::Error;
        }));
}

ValidationReport validateScene(const Scene& scene) {
    ValidationReport report;

    if (scene.metadata.schemaMajor != sceneSchemaMajor
        || scene.metadata.schemaMinor > sceneSchemaMinor) {
        add(report, ValidationCode::UnsupportedSchema, EntityKind::Scene, 0,
            "unsupported scene schema version");
    }
    if (scene.metadata.lengthUnit != LengthUnit::Metre
        || scene.metadata.handedness != CoordinateHandedness::RightHanded
        || scene.metadata.upAxis != UpAxis::PositiveZ) {
        add(report, ValidationCode::InvalidCoordinateSystem,
            EntityKind::Scene, 0,
            "scene coordinates must be metres, right-handed, and Z-up");
    }
    if (!finite(scene.metadata.sourceLengthToMeters)
        || !(scene.metadata.sourceLengthToMeters > 0.0)) {
        add(report, ValidationCode::InvalidPhysicalValue, EntityKind::Scene, 0,
            "source length conversion must be finite and positive");
    }
    if (scene.metadata.designChecksum.empty()
        || scene.metadata.exporterVersion.empty()) {
        add(report, ValidationCode::MissingMetadata, EntityKind::Scene, 0,
            "design checksum and exporter version are required");
    }

    validateIds(scene.regions, EntityKind::Region, "region", report);
    validateIds(scene.vertices, EntityKind::Vertex, "vertex", report);
    validateIds(scene.fabricMaterials, EntityKind::FabricMaterial,
                "fabric material", report);
    validateIds(scene.seamMaterials, EntityKind::SeamMaterial,
                "seam material", report);
    validateIds(scene.triangles, EntityKind::Triangle, "triangle", report);
    validateIds(scene.openings, EntityKind::Opening, "opening", report);
    validateIds(scene.seams, EntityKind::Seam, "seam", report);
    validateIds(scene.lineMaterials, EntityKind::LineMaterial,
                "line material", report);
    validateIds(scene.pilots, EntityKind::Pilot, "pilot", report);
    validateIds(scene.suspensionJunctions, EntityKind::SuspensionJunction,
                "suspension junction", report);
    validateIds(scene.attachments, EntityKind::Attachment, "attachment", report);
    validateIds(scene.suspensionLines, EntityKind::SuspensionLine,
                "suspension line", report);

    const auto regionIds = idSet(scene.regions);
    const auto vertexIds = idSet(scene.vertices);
    const auto fabricMaterialIds = idSet(scene.fabricMaterials);
    const auto seamMaterialIds = idSet(scene.seamMaterials);
    const auto lineMaterialIds = idSet(scene.lineMaterials);
    const auto pilotIds = idSet(scene.pilots);
    const auto suspensionJunctionIds = idSet(scene.suspensionJunctions);
    const auto attachmentIds = idSet(scene.attachments);

    std::size_t outsideCount = 0;
    for (const FluidRegion* region : sortedById(scene.regions)) {
        if (!valid(region->kind)) {
            add(report, ValidationCode::InvalidEnumValue, EntityKind::Region,
                region->id, "region kind is invalid");
        } else if (region->kind == RegionKind::Outside) {
            ++outsideCount;
        }
        if (region->name.empty()) {
            add(report, ValidationCode::MissingMetadata, EntityKind::Region,
                region->id, "region name is required");
        }
    }
    if (outsideCount != 1) {
        add(report, ValidationCode::InvalidSideRegions, EntityKind::Scene, 0,
            "scene must contain exactly one outside fluid region");
    }

    std::unordered_map<StableId, const Vertex*> verticesById;
    verticesById.reserve(scene.vertices.size());
    std::unordered_map<StableId, std::size_t> vertexIdCounts;
    vertexIdCounts.reserve(scene.vertices.size());
    for (const Vertex& vertex : scene.vertices) {
        ++vertexIdCounts[vertex.id];
    }
    for (const Vertex* vertex : sortedById(scene.vertices)) {
        // A face that names a duplicated vertex ID is already invalid. Do not
        // choose one duplicate arbitrarily for the geometric checks, because
        // that would make diagnostics depend on insertion order.
        if (vertexIdCounts[vertex->id] == 1) {
            verticesById.emplace(vertex->id, vertex);
        }
        if (!finite(vertex->positionMeters)) {
            add(report, ValidationCode::NonFiniteValue, EntityKind::Vertex,
                vertex->id, "vertex position is not finite");
        }
    }

    for (const FabricMaterial* material
         : sortedById(scene.fabricMaterials)) {
        const bool allFinite = finite(material->warpStiffnessNewtonsPerMeter)
            && finite(material->weftStiffnessNewtonsPerMeter)
            && finite(material->shearStiffnessNewtonsPerMeter)
            && finite(material->bendingStiffnessNewtonMeters)
            && finite(material->arealDensityKgPerSquareMeter)
            && finite(material->dampingSeconds)
            && finite(material->porosityFraction)
            && finite(material->permeabilitySquareMeters);
        if (!allFinite) {
            add(report, ValidationCode::NonFiniteValue,
                EntityKind::FabricMaterial, material->id,
                "fabric material contains a non-finite value");
        } else if (!(material->warpStiffnessNewtonsPerMeter > 0.0)
                   || !(material->weftStiffnessNewtonsPerMeter > 0.0)
                   || !(material->shearStiffnessNewtonsPerMeter > 0.0)
                   || material->bendingStiffnessNewtonMeters < 0.0
                   || !(material->arealDensityKgPerSquareMeter > 0.0)
                   || material->dampingSeconds < 0.0
                   || material->porosityFraction < 0.0
                   || material->porosityFraction > 1.0
                   || material->permeabilitySquareMeters < 0.0) {
            add(report, ValidationCode::InvalidPhysicalValue,
                EntityKind::FabricMaterial, material->id,
                "fabric material has an invalid physical value");
        }
        if (material->name.empty()) {
            add(report, ValidationCode::MissingMetadata,
                EntityKind::FabricMaterial, material->id,
                "fabric material name is required");
        }
    }

    for (const SeamMaterial* material : sortedById(scene.seamMaterials)) {
        if (!finite(material->linearDensityKgPerMeter)
            || !finite(material->axialStiffnessNewtons)) {
            add(report, ValidationCode::NonFiniteValue,
                EntityKind::SeamMaterial, material->id,
                "seam material contains a non-finite value");
        } else if (!(material->linearDensityKgPerMeter > 0.0)
                   || !(material->axialStiffnessNewtons > 0.0)) {
            add(report, ValidationCode::InvalidPhysicalValue,
                EntityKind::SeamMaterial, material->id,
                "seam material density and axial stiffness must be positive");
        }
        if (material->name.empty()) {
            add(report, ValidationCode::MissingMetadata,
                EntityKind::SeamMaterial, material->id,
                "seam material name is required");
        }
    }

    std::map<std::array<StableId, 3>, StableId> triangleGeometry;
    for (const Triangle* triangle : sortedById(scene.triangles)) {
        if (!valid(triangle->role)) {
            add(report, ValidationCode::InvalidEnumValue,
                EntityKind::Triangle, triangle->id,
                "triangle surface role is invalid");
        }
        for (const Vec2& coordinate : triangle->materialCoordinates) {
            if (!finite(coordinate)) {
                add(report, ValidationCode::NonFiniteValue,
                    EntityKind::Triangle, triangle->id,
                    "triangle material coordinate is not finite");
                break;
            }
        }
        std::unordered_set<StableId> localVertices;
        for (StableId vertexId : triangle->vertexIds) {
            if (!vertexIds.contains(vertexId)) {
                add(report, ValidationCode::MissingVertexReference,
                    EntityKind::Triangle, triangle->id,
                    "triangle references a missing vertex");
            }
            localVertices.insert(vertexId);
        }
        if (localVertices.size() != 3 || degenerate(*triangle, verticesById)) {
            add(report, ValidationCode::DegenerateTriangle,
                EntityKind::Triangle, triangle->id,
                "triangle has repeated or effectively collinear vertices");
        }
        if (!fabricMaterialIds.contains(triangle->materialId)) {
            add(report, ValidationCode::MissingMaterialReference,
                EntityKind::Triangle, triangle->id,
                "triangle references a missing fabric material");
        }
        if (triangle->sheetId == invalidStableId) {
            add(report, ValidationCode::InvalidId, EntityKind::Triangle,
                triangle->id,
                "triangle fabric sheet ID must be non-zero");
        }
        const bool negativeExists =
            regionIds.contains(triangle->negativeSideRegionId);
        const bool positiveExists =
            regionIds.contains(triangle->positiveSideRegionId);
        if (!negativeExists || !positiveExists) {
            add(report, ValidationCode::MissingRegionReference,
                EntityKind::Triangle, triangle->id,
                "triangle references a missing side region");
        }
        const bool isInternalSheet =
            triangle->role == SurfaceRole::Diagonal
            || triangle->role == SurfaceRole::MiniRib;
        if (triangle->negativeSideRegionId
                == triangle->positiveSideRegionId
            && !isInternalSheet) {
            add(report, ValidationCode::InvalidSideRegions,
                EntityKind::Triangle, triangle->id,
                "skin and rib triangle side regions must be distinct");
        }
        std::array<StableId, 3> key = triangle->vertexIds;
        std::sort(key.begin(), key.end());
        const auto [found, inserted] = triangleGeometry.emplace(key,
                                                                triangle->id);
        if (!inserted) {
            add(report, ValidationCode::DuplicateTriangle,
                EntityKind::Triangle, triangle->id,
                "triangle duplicates an existing geometric face");
        }
    }

    for (const Opening* opening : sortedById(scene.openings)) {
        if (!valid(opening->role)) {
            add(report, ValidationCode::InvalidEnumValue, EntityKind::Opening,
                opening->id, "opening role is invalid");
        }
        bool invalidLoop = opening->orderedVertexIds.size() < 3;
        std::unordered_set<StableId> loopVertices;
        for (StableId vertexId : opening->orderedVertexIds) {
            if (!vertexIds.contains(vertexId)) {
                add(report, ValidationCode::MissingVertexReference,
                    EntityKind::Opening, opening->id,
                    "opening references a missing boundary vertex");
            }
            if (!loopVertices.insert(vertexId).second) {
                invalidLoop = true;
            }
        }
        if (invalidLoop) {
            add(report, ValidationCode::InvalidOpening, EntityKind::Opening,
                opening->id,
                "opening boundary must contain at least three unique vertices");
        }
        if (!regionIds.contains(opening->negativeSideRegionId)
            || !regionIds.contains(opening->positiveSideRegionId)) {
            add(report, ValidationCode::MissingRegionReference,
                EntityKind::Opening, opening->id,
                "opening references a missing side region");
        }
        if (opening->negativeSideRegionId
            == opening->positiveSideRegionId) {
            add(report, ValidationCode::InvalidSideRegions,
                EntityKind::Opening, opening->id,
                "opening side regions must be distinct");
        }
    }

    std::map<std::array<StableId, 2>, StableId> seamVertexPairs;
    for (const Seam* seam : sortedById(scene.seams)) {
        bool invalidChains = seam->firstOrderedVertexIds.size() < 2
            || seam->firstOrderedVertexIds.size()
                   != seam->secondOrderedVertexIds.size();
        std::unordered_set<StableId> allChainVertices;
        const auto validateChain = [&](const std::vector<StableId>& chain,
                                       const char* label) {
            std::unordered_set<StableId> localVertices;
            for (const StableId vertexId : chain) {
                if (!vertexIds.contains(vertexId)) {
                    add(report, ValidationCode::MissingVertexReference,
                        EntityKind::Seam, seam->id,
                        std::string("seam ") + label
                            + " chain references a missing vertex");
                }
                if (!localVertices.insert(vertexId).second
                    || !allChainVertices.insert(vertexId).second) {
                    invalidChains = true;
                }
            }
        };
        validateChain(seam->firstOrderedVertexIds, "first");
        validateChain(seam->secondOrderedVertexIds, "second");
        if (seam->firstOrderedVertexIds.size()
            == seam->secondOrderedVertexIds.size()) {
            for (std::size_t index = 0;
                 index < seam->firstOrderedVertexIds.size(); ++index) {
                std::array<StableId, 2> pair{
                    seam->firstOrderedVertexIds[index],
                    seam->secondOrderedVertexIds[index]};
                if (pair[0] == pair[1]) {
                    invalidChains = true;
                }
                std::ranges::sort(pair);
                if (!seamVertexPairs.emplace(pair, seam->id).second) {
                    invalidChains = true;
                }
            }
        }
        if (invalidChains) {
            add(report, ValidationCode::InvalidSeam, EntityKind::Seam,
                seam->id,
                "seam needs two disjoint equal-length ordered chains with at least two unique vertex pairs");
        }
        if (!seamMaterialIds.contains(seam->materialId)) {
            add(report, ValidationCode::MissingMaterialReference,
                EntityKind::Seam, seam->id,
                "seam references a missing seam material");
        }
    }

    for (const LineMaterial* material : sortedById(scene.lineMaterials)) {
        const bool allFinite = finite(material->diameterMeters)
            && finite(material->linearDensityKgPerMeter)
            && finite(material->axialStiffnessNewtons)
            && finite(material->dragCoefficient);
        if (!allFinite) {
            add(report, ValidationCode::NonFiniteValue,
                EntityKind::LineMaterial, material->id,
                "line material contains a non-finite value");
        } else if (!(material->diameterMeters > 0.0)
                   || !(material->linearDensityKgPerMeter > 0.0)
                   || !(material->axialStiffnessNewtons > 0.0)
                   || material->dragCoefficient < 0.0) {
            add(report, ValidationCode::InvalidPhysicalValue,
                EntityKind::LineMaterial, material->id,
                "line material has an invalid physical value");
        }
        if (material->name.empty()) {
            add(report, ValidationCode::MissingMetadata,
                EntityKind::LineMaterial, material->id,
                "line material name is required");
        }
    }

    for (const Pilot* pilot : sortedById(scene.pilots)) {
        const double quaternionNormSquared = pilot->bodyToWorld.w
                * pilot->bodyToWorld.w
            + pilot->bodyToWorld.x * pilot->bodyToWorld.x
            + pilot->bodyToWorld.y * pilot->bodyToWorld.y
            + pilot->bodyToWorld.z * pilot->bodyToWorld.z;
        if (!finite(pilot->massKg)
            || !finite(pilot->centerOfMassPositionMeters)
            || !finite(pilot->linearVelocityMetersPerSecond)
            || !finite(pilot->bodyToWorld)
            || !finite(pilot->principalInertiaKgSquareMeters)) {
            add(report, ValidationCode::NonFiniteValue, EntityKind::Pilot,
                pilot->id, "pilot contains a non-finite value");
        } else if (!(pilot->massKg > 0.0)
                   || !(pilot->principalInertiaKgSquareMeters.x > 0.0)
                   || !(pilot->principalInertiaKgSquareMeters.y > 0.0)
                   || !(pilot->principalInertiaKgSquareMeters.z > 0.0)
                   || std::abs(quaternionNormSquared - 1.0) > 1.0e-6) {
            add(report, ValidationCode::InvalidPhysicalValue,
                EntityKind::Pilot, pilot->id,
                "pilot mass/inertia must be positive and orientation normalized");
        }
        if (pilot->name.empty()) {
            add(report, ValidationCode::MissingMetadata, EntityKind::Pilot,
                pilot->id, "pilot name is required");
        }
    }

    for (const SuspensionJunction* junction
         : sortedById(scene.suspensionJunctions)) {
        if (!finite(junction->positionMeters) || !finite(junction->massKg)) {
            add(report, ValidationCode::NonFiniteValue,
                EntityKind::SuspensionJunction, junction->id,
                "suspension junction contains a non-finite value");
        } else if (junction->massKg < 0.0
                   || (!junction->fixed && !(junction->massKg > 0.0))) {
            add(report, ValidationCode::InvalidPhysicalValue,
                EntityKind::SuspensionJunction, junction->id,
                "dynamic suspension junction requires positive explicit mass");
        }
        if (vertexIds.contains(junction->id)) {
            add(report, ValidationCode::InvalidId,
                EntityKind::SuspensionJunction, junction->id,
                "suspension junction ID must not collide with a surface vertex ID");
        }
    }

    for (const Attachment* attachment : sortedById(scene.attachments)) {
        if (!valid(attachment->kind)) {
            add(report, ValidationCode::InvalidEnumValue,
                EntityKind::Attachment, attachment->id,
                "attachment kind is invalid");
            continue;
        }
        if (attachment->kind == AttachmentKind::SurfaceVertex) {
            if (!vertexIds.contains(attachment->vertexId)
                || attachment->pilotId != invalidStableId
                || attachment->suspensionJunctionId != invalidStableId
                || !exactlyZero(attachment->pilotLocalPositionMeters)) {
                add(report, ValidationCode::InvalidAttachmentTarget,
                    EntityKind::Attachment, attachment->id,
                    "surface attachment must reference one existing vertex");
            }
        } else if (attachment->kind == AttachmentKind::PilotHarness) {
            if (!pilotIds.contains(attachment->pilotId)) {
                add(report, ValidationCode::MissingPilotReference,
                    EntityKind::Attachment, attachment->id,
                    "pilot attachment references a missing pilot");
            }
            if (attachment->vertexId != invalidStableId
                || attachment->suspensionJunctionId != invalidStableId
                || !finite(attachment->pilotLocalPositionMeters)) {
                add(report, ValidationCode::InvalidAttachmentTarget,
                    EntityKind::Attachment, attachment->id,
                    "pilot attachment target or local position is invalid");
            }
        } else {
            if (!suspensionJunctionIds.contains(
                    attachment->suspensionJunctionId)
                || attachment->vertexId != invalidStableId
                || attachment->pilotId != invalidStableId
                || !exactlyZero(attachment->pilotLocalPositionMeters)) {
                add(report, ValidationCode::InvalidAttachmentTarget,
                    EntityKind::Attachment, attachment->id,
                    "junction attachment must reference one existing suspension junction");
            }
        }
    }

    for (const SuspensionLine* line
         : sortedById(scene.suspensionLines)) {
        if (!valid(line->role)) {
            add(report, ValidationCode::InvalidEnumValue,
                EntityKind::SuspensionLine, line->id,
                "suspension line role is invalid");
        }
        if (!attachmentIds.contains(line->startAttachmentId)
            || !attachmentIds.contains(line->endAttachmentId)) {
            add(report, ValidationCode::DanglingLineAttachment,
                EntityKind::SuspensionLine, line->id,
                "suspension line references a missing attachment");
        }
        if (line->startAttachmentId == line->endAttachmentId
            || !finite(line->restLengthMeters)
            || !(line->restLengthMeters > 0.0)) {
            add(report, ValidationCode::InvalidSuspensionLine,
                EntityKind::SuspensionLine, line->id,
                "suspension line endpoints and rest length are invalid");
        }
        if (!lineMaterialIds.contains(line->materialId)) {
            add(report, ValidationCode::MissingMaterialReference,
                EntityKind::SuspensionLine, line->id,
                "suspension line references a missing line material");
        }
    }

    std::sort(report.diagnostics.begin(), report.diagnostics.end(),
              [](const ValidationDiagnostic& left,
                 const ValidationDiagnostic& right) {
        if (left.entityKind != right.entityKind) {
            return left.entityKind < right.entityKind;
        }
        if (left.entityId != right.entityId) {
            return left.entityId < right.entityId;
        }
        if (left.code != right.code) {
            return left.code < right.code;
        }
        return left.message < right.message;
    });
    return report;
}

bool writeScene(const Scene& scene,
                std::ostream& output,
                std::string* errorMessage) {
    const ValidationReport validation = validateScene(scene);
    if (!validation.ok()) {
        if (errorMessage) {
            *errorMessage = diagnosticSummary(validation);
        }
        return false;
    }

    std::string safetyError;
    if (!withinBinarySafetyLimits(scene, safetyError)) {
        if (errorMessage) {
            *errorMessage = std::move(safetyError);
        }
        return false;
    }

    try {
        const Scene canonical = canonicalScene(scene);
        BinaryWriter writer(output);
        writer.bytes(binaryMagic.data(), binaryMagic.size());
        writer.u32(sceneBinaryVersion);
        writer.u32(canonical.metadata.schemaMajor);
        writer.u32(canonical.metadata.schemaMinor);
        writer.u8(static_cast<std::uint8_t>(canonical.metadata.lengthUnit));
        writer.u8(static_cast<std::uint8_t>(canonical.metadata.handedness));
        writer.u8(static_cast<std::uint8_t>(canonical.metadata.upAxis));
        writer.real(canonical.metadata.sourceLengthToMeters);
        writer.string(canonical.metadata.designChecksum);
        writer.string(canonical.metadata.exporterVersion);

        writer.count(canonical.regions.size());
        for (const FluidRegion& region : canonical.regions) {
            writer.u64(region.id);
            writer.u8(static_cast<std::uint8_t>(region.kind));
            writer.string(region.name);
        }
        writer.count(canonical.vertices.size());
        for (const Vertex& vertex : canonical.vertices) {
            writer.u64(vertex.id);
            writeVec3(writer, vertex.positionMeters);
        }
        writer.count(canonical.fabricMaterials.size());
        for (const FabricMaterial& material : canonical.fabricMaterials) {
            writer.u64(material.id);
            writer.string(material.name);
            writer.real(material.warpStiffnessNewtonsPerMeter);
            writer.real(material.weftStiffnessNewtonsPerMeter);
            writer.real(material.shearStiffnessNewtonsPerMeter);
            writer.real(material.bendingStiffnessNewtonMeters);
            writer.real(material.arealDensityKgPerSquareMeter);
            writer.real(material.dampingSeconds);
            writer.real(material.porosityFraction);
            writer.real(material.permeabilitySquareMeters);
        }
        writer.count(canonical.seamMaterials.size());
        for (const SeamMaterial& material : canonical.seamMaterials) {
            writer.u64(material.id);
            writer.string(material.name);
            writer.real(material.linearDensityKgPerMeter);
            writer.real(material.axialStiffnessNewtons);
        }
        writer.count(canonical.triangles.size());
        for (const Triangle& triangle : canonical.triangles) {
            writer.u64(triangle.id);
            for (StableId vertexId : triangle.vertexIds) {
                writer.u64(vertexId);
            }
            for (const Vec2& coordinate : triangle.materialCoordinates) {
                writeVec2(writer, coordinate);
            }
            writer.u64(triangle.negativeSideRegionId);
            writer.u64(triangle.positiveSideRegionId);
            writer.u64(triangle.materialId);
            writer.u64(triangle.sheetId);
            writer.u8(static_cast<std::uint8_t>(triangle.role));
        }
        writer.count(canonical.openings.size());
        for (const Opening& opening : canonical.openings) {
            writer.u64(opening.id);
            writer.count(opening.orderedVertexIds.size());
            for (StableId vertexId : opening.orderedVertexIds) {
                writer.u64(vertexId);
            }
            writer.u64(opening.negativeSideRegionId);
            writer.u64(opening.positiveSideRegionId);
            writer.u8(static_cast<std::uint8_t>(opening.role));
        }
        writer.count(canonical.seams.size());
        for (const Seam& seam : canonical.seams) {
            writer.u64(seam.id);
            writer.u64(seam.materialId);
            writer.count(seam.firstOrderedVertexIds.size());
            for (const StableId vertexId : seam.firstOrderedVertexIds) {
                writer.u64(vertexId);
            }
            writer.count(seam.secondOrderedVertexIds.size());
            for (const StableId vertexId : seam.secondOrderedVertexIds) {
                writer.u64(vertexId);
            }
        }
        writer.count(canonical.lineMaterials.size());
        for (const LineMaterial& material : canonical.lineMaterials) {
            writer.u64(material.id);
            writer.string(material.name);
            writer.real(material.diameterMeters);
            writer.real(material.linearDensityKgPerMeter);
            writer.real(material.axialStiffnessNewtons);
            writer.real(material.dragCoefficient);
        }
        writer.count(canonical.pilots.size());
        for (const Pilot& pilot : canonical.pilots) {
            writer.u64(pilot.id);
            writer.string(pilot.name);
            writer.real(pilot.massKg);
            writeVec3(writer, pilot.centerOfMassPositionMeters);
            writeVec3(writer, pilot.linearVelocityMetersPerSecond);
            writeQuaternion(writer, pilot.bodyToWorld);
            writeVec3(writer, pilot.principalInertiaKgSquareMeters);
        }
        writer.count(canonical.suspensionJunctions.size());
        for (const SuspensionJunction& junction :
             canonical.suspensionJunctions) {
            writer.u64(junction.id);
            writeVec3(writer, junction.positionMeters);
            writer.real(junction.massKg);
            writer.u8(junction.fixed ? 1 : 0);
        }
        writer.count(canonical.attachments.size());
        for (const Attachment& attachment : canonical.attachments) {
            writer.u64(attachment.id);
            writer.u8(static_cast<std::uint8_t>(attachment.kind));
            writer.u64(attachment.vertexId);
            writer.u64(attachment.pilotId);
            writeVec3(writer, attachment.pilotLocalPositionMeters);
            writer.u64(attachment.suspensionJunctionId);
        }
        writer.count(canonical.suspensionLines.size());
        for (const SuspensionLine& line : canonical.suspensionLines) {
            writer.u64(line.id);
            writer.u64(line.startAttachmentId);
            writer.u64(line.endAttachmentId);
            writer.u64(line.materialId);
            writer.real(line.restLengthMeters);
            writer.u8(static_cast<std::uint8_t>(line.role));
        }
        if (!writer.ok()) {
            if (errorMessage) {
                *errorMessage = writer.error();
            }
            return false;
        }
    } catch (const std::exception& exception) {
        if (errorMessage) {
            *errorMessage = std::string("failed to serialize scene: ")
                + exception.what();
        }
        return false;
    }

    if (errorMessage) {
        errorMessage->clear();
    }
    return true;
}

bool readScene(std::istream& input,
               Scene& scene,
               std::string* errorMessage) {
    try {
        BinaryReader reader(input);
        std::array<char, binaryMagic.size()> magic{};
        reader.bytes(magic.data(), magic.size());
        if (reader.ok() && magic != binaryMagic) {
            reader.fail("scene stream has an invalid magic header");
        }
        const std::uint32_t binaryVersion = reader.u32();
        if (reader.ok() && binaryVersion != sceneBinaryVersion) {
            reader.fail("scene stream uses an unsupported binary version");
        }

        Scene decoded;
        decoded.metadata.schemaMajor = reader.u32();
        decoded.metadata.schemaMinor = reader.u32();
        decoded.metadata.lengthUnit = static_cast<LengthUnit>(reader.u8());
        decoded.metadata.handedness =
            static_cast<CoordinateHandedness>(reader.u8());
        decoded.metadata.upAxis = static_cast<UpAxis>(reader.u8());
        decoded.metadata.sourceLengthToMeters = reader.real();
        decoded.metadata.designChecksum = reader.string();
        decoded.metadata.exporterVersion = reader.string();

        resizeFromCount(decoded.regions, reader, "region");
        for (FluidRegion& region : decoded.regions) {
            region.id = reader.u64();
            region.kind = static_cast<RegionKind>(reader.u8());
            region.name = reader.string();
        }
        resizeFromCount(decoded.vertices, reader, "vertex");
        for (Vertex& vertex : decoded.vertices) {
            vertex.id = reader.u64();
            vertex.positionMeters = readVec3(reader);
        }
        resizeFromCount(decoded.fabricMaterials, reader, "fabric material");
        for (FabricMaterial& material : decoded.fabricMaterials) {
            material.id = reader.u64();
            material.name = reader.string();
            material.warpStiffnessNewtonsPerMeter = reader.real();
            material.weftStiffnessNewtonsPerMeter = reader.real();
            material.shearStiffnessNewtonsPerMeter = reader.real();
            material.bendingStiffnessNewtonMeters = reader.real();
            material.arealDensityKgPerSquareMeter = reader.real();
            material.dampingSeconds = reader.real();
            material.porosityFraction = reader.real();
            material.permeabilitySquareMeters = reader.real();
        }
        resizeFromCount(decoded.seamMaterials, reader, "seam material");
        for (SeamMaterial& material : decoded.seamMaterials) {
            material.id = reader.u64();
            material.name = reader.string();
            material.linearDensityKgPerMeter = reader.real();
            material.axialStiffnessNewtons = reader.real();
        }
        resizeFromCount(decoded.triangles, reader, "triangle");
        for (Triangle& triangle : decoded.triangles) {
            triangle.id = reader.u64();
            for (StableId& vertexId : triangle.vertexIds) {
                vertexId = reader.u64();
            }
            for (Vec2& coordinate : triangle.materialCoordinates) {
                coordinate = readVec2(reader);
            }
            triangle.negativeSideRegionId = reader.u64();
            triangle.positiveSideRegionId = reader.u64();
            triangle.materialId = reader.u64();
            triangle.sheetId = reader.u64();
            triangle.role = static_cast<SurfaceRole>(reader.u8());
        }
        resizeFromCount(decoded.openings, reader, "opening");
        for (Opening& opening : decoded.openings) {
            opening.id = reader.u64();
            resizeFromCount(opening.orderedVertexIds, reader,
                            "opening boundary vertex");
            for (StableId& vertexId : opening.orderedVertexIds) {
                vertexId = reader.u64();
            }
            opening.negativeSideRegionId = reader.u64();
            opening.positiveSideRegionId = reader.u64();
            opening.role = static_cast<OpeningRole>(reader.u8());
        }
        resizeFromCount(decoded.seams, reader, "seam");
        for (Seam& seam : decoded.seams) {
            seam.id = reader.u64();
            seam.materialId = reader.u64();
            resizeFromCount(seam.firstOrderedVertexIds, reader,
                            "first seam-chain vertex");
            for (StableId& vertexId : seam.firstOrderedVertexIds) {
                vertexId = reader.u64();
            }
            resizeFromCount(seam.secondOrderedVertexIds, reader,
                            "second seam-chain vertex");
            for (StableId& vertexId : seam.secondOrderedVertexIds) {
                vertexId = reader.u64();
            }
        }
        resizeFromCount(decoded.lineMaterials, reader, "line material");
        for (LineMaterial& material : decoded.lineMaterials) {
            material.id = reader.u64();
            material.name = reader.string();
            material.diameterMeters = reader.real();
            material.linearDensityKgPerMeter = reader.real();
            material.axialStiffnessNewtons = reader.real();
            material.dragCoefficient = reader.real();
        }
        resizeFromCount(decoded.pilots, reader, "pilot");
        for (Pilot& pilot : decoded.pilots) {
            pilot.id = reader.u64();
            pilot.name = reader.string();
            pilot.massKg = reader.real();
            pilot.centerOfMassPositionMeters = readVec3(reader);
            pilot.linearVelocityMetersPerSecond = readVec3(reader);
            pilot.bodyToWorld = readQuaternion(reader);
            pilot.principalInertiaKgSquareMeters = readVec3(reader);
        }
        resizeFromCount(decoded.suspensionJunctions, reader,
                        "suspension junction");
        for (SuspensionJunction& junction : decoded.suspensionJunctions) {
            junction.id = reader.u64();
            junction.positionMeters = readVec3(reader);
            junction.massKg = reader.real();
            const std::uint8_t fixed = reader.u8();
            if (fixed > 1) {
                reader.fail("suspension junction fixed flag is invalid");
            }
            junction.fixed = fixed != 0;
        }
        resizeFromCount(decoded.attachments, reader, "attachment");
        for (Attachment& attachment : decoded.attachments) {
            attachment.id = reader.u64();
            attachment.kind = static_cast<AttachmentKind>(reader.u8());
            attachment.vertexId = reader.u64();
            attachment.pilotId = reader.u64();
            attachment.pilotLocalPositionMeters = readVec3(reader);
            attachment.suspensionJunctionId = reader.u64();
        }
        resizeFromCount(decoded.suspensionLines, reader, "suspension line");
        for (SuspensionLine& line : decoded.suspensionLines) {
            line.id = reader.u64();
            line.startAttachmentId = reader.u64();
            line.endAttachmentId = reader.u64();
            line.materialId = reader.u64();
            line.restLengthMeters = reader.real();
            line.role = static_cast<SuspensionLineRole>(reader.u8());
        }

        if (!reader.ok()) {
            if (errorMessage) {
                *errorMessage = reader.error();
            }
            return false;
        }
        if (input.peek() != std::char_traits<char>::eof()) {
            if (errorMessage) {
                *errorMessage = "scene stream contains trailing data";
            }
            return false;
        }

        const ValidationReport validation = validateScene(decoded);
        if (!validation.ok()) {
            if (errorMessage) {
                *errorMessage = diagnosticSummary(validation);
            }
            return false;
        }
        scene = std::move(decoded);
    } catch (const std::exception& exception) {
        if (errorMessage) {
            *errorMessage = std::string("failed to deserialize scene: ")
                + exception.what();
        }
        return false;
    }

    if (errorMessage) {
        errorMessage->clear();
    }
    return true;
}

} // namespace simwing::fsi

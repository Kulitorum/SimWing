#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace simwing::fsi {

// Stable IDs are unique within each entity collection. Vertex and
// SuspensionJunction additionally share one ID namespace because both become
// Structure nodes and DiagnosticVertex entities without ID transformation.
using StableId = std::uint64_t;

inline constexpr StableId invalidStableId = 0;
inline constexpr std::uint32_t sceneSchemaMajor = 2;
inline constexpr std::uint32_t sceneSchemaMinor = 1;
inline constexpr std::uint32_t sceneBinaryVersion = 2;

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Quaternion {
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

enum class LengthUnit : std::uint8_t {
    Metre = 1,
};

enum class CoordinateHandedness : std::uint8_t {
    RightHanded = 1,
};

enum class UpAxis : std::uint8_t {
    PositiveZ = 1,
};

struct SceneMetadata {
    std::uint32_t schemaMajor = sceneSchemaMajor;
    std::uint32_t schemaMinor = sceneSchemaMinor;
    LengthUnit lengthUnit = LengthUnit::Metre;
    CoordinateHandedness handedness = CoordinateHandedness::RightHanded;
    UpAxis upAxis = UpAxis::PositiveZ;

    // Geometry in the scene is always SI. These fields retain enough
    // provenance to reproduce the export boundary.
    double sourceLengthToMeters = 1.0;
    std::string designChecksum;
    std::string exporterVersion;
};

enum class RegionKind : std::uint8_t {
    Outside = 1,
    Cell = 2,
};

struct FluidRegion {
    StableId id = invalidStableId;
    RegionKind kind = RegionKind::Cell;
    std::string name;
};

struct Vertex {
    StableId id = invalidStableId;
    Vec3 positionMeters;
};

struct FabricMaterial {
    StableId id = invalidStableId;
    std::string name;
    double warpStiffnessNewtonsPerMeter = 0.0;
    double weftStiffnessNewtonsPerMeter = 0.0;
    double shearStiffnessNewtonsPerMeter = 0.0;
    double bendingStiffnessNewtonMeters = 0.0;
    double arealDensityKgPerSquareMeter = 0.0;
    double dampingSeconds = 0.0;
    double porosityFraction = 0.0;
    double permeabilitySquareMeters = 0.0;
};

// Measured properties of the assembled seam line. They intentionally do not
// prescribe how discrete paired vertices share tributary length or stitch
// load; that belongs to a separately verified structural seam model.
struct SeamMaterial {
    StableId id = invalidStableId;
    std::string name;
    double linearDensityKgPerMeter = 0.0;
    double axialStiffnessNewtons = 0.0;
};

enum class SurfaceRole : std::uint8_t {
    Skin = 1,
    Rib = 2,
    Diagonal = 3,
    MiniRib = 4,
};

struct Triangle {
    StableId id = invalidStableId;
    std::array<StableId, 3> vertexIds{};
    // Undeformed material-chart coordinates in metres. Their signed area is
    // the authoritative fabric area used for mass and membrane assembly.
    std::array<Vec2, 3> materialCoordinates{};
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    StableId materialId = invalidStableId;
    // Stable identity of the authored fabric sheet. Triangles sharing a mesh
    // edge bend together only when this ID also matches. This prevents a
    // welded skin/rib or skin/diagonal attachment edge from acquiring a
    // fictitious fabric hinge.
    StableId sheetId = invalidStableId;
    SurfaceRole role = SurfaceRole::Skin;
};

enum class OpeningRole : std::uint8_t {
    Intake = 1,
    Vent = 2,
    Crossport = 3,
};

struct Opening {
    StableId id = invalidStableId;
    std::vector<StableId> orderedVertexIds;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    OpeningRole role = OpeningRole::Intake;
};

// The chains have equal cardinality and pair element-by-element. Distinct
// vertex IDs preserve the intended sewn coincidence without welding topology
// or assuming a constraint stiffness/discretization rule.
struct Seam {
    StableId id = invalidStableId;
    StableId materialId = invalidStableId;
    std::vector<StableId> firstOrderedVertexIds;
    std::vector<StableId> secondOrderedVertexIds;
};

struct LineMaterial {
    StableId id = invalidStableId;
    std::string name;
    double diameterMeters = 0.0;
    double linearDensityKgPerMeter = 0.0;
    double axialStiffnessNewtons = 0.0;
    double dragCoefficient = 0.0;
};

struct Pilot {
    StableId id = invalidStableId;
    std::string name;
    double massKg = 0.0;
    Vec3 centerOfMassPositionMeters;
    Vec3 linearVelocityMetersPerSecond;
    Quaternion bodyToWorld;
    Vec3 principalInertiaKgSquareMeters;
};

enum class AttachmentKind : std::uint8_t {
    SurfaceVertex = 1,
    PilotHarness = 2,
    SuspensionJunction = 3,
};

struct SuspensionJunction {
    StableId id = invalidStableId;
    Vec3 positionMeters;
    double massKg = 0.0;
    bool fixed = false;
};

struct Attachment {
    StableId id = invalidStableId;
    AttachmentKind kind = AttachmentKind::SurfaceVertex;

    // Exactly one target is active according to kind. A surface attachment
    // identifies an existing scene vertex. A pilot attachment identifies a
    // pilot and locates the harness point in that pilot's body frame.
    StableId vertexId = invalidStableId;
    StableId pilotId = invalidStableId;
    Vec3 pilotLocalPositionMeters;
    StableId suspensionJunctionId = invalidStableId;
};

enum class SuspensionLineRole : std::uint8_t {
    Suspension = 1,
    Riser = 2,
    Brake = 3,
    Harness = 4,
};

struct SuspensionLine {
    StableId id = invalidStableId;
    StableId startAttachmentId = invalidStableId;
    StableId endAttachmentId = invalidStableId;
    StableId materialId = invalidStableId;
    double restLengthMeters = 0.0;
    SuspensionLineRole role = SuspensionLineRole::Suspension;
};

struct Scene {
    SceneMetadata metadata;
    std::vector<FluidRegion> regions;
    std::vector<Vertex> vertices;
    std::vector<FabricMaterial> fabricMaterials;
    std::vector<SeamMaterial> seamMaterials;
    std::vector<Triangle> triangles;
    std::vector<Opening> openings;
    std::vector<Seam> seams;
    std::vector<LineMaterial> lineMaterials;
    std::vector<Pilot> pilots;
    std::vector<SuspensionJunction> suspensionJunctions;
    std::vector<Attachment> attachments;
    std::vector<SuspensionLine> suspensionLines;
};

enum class DiagnosticSeverity : std::uint8_t {
    Error = 1,
};

enum class EntityKind : std::uint8_t {
    Scene = 1,
    Region = 2,
    Vertex = 3,
    FabricMaterial = 4,
    Triangle = 5,
    Opening = 6,
    LineMaterial = 7,
    Pilot = 8,
    Attachment = 9,
    SuspensionLine = 10,
    SeamMaterial = 11,
    Seam = 12,
    SuspensionJunction = 13,
};

enum class ValidationCode : std::uint16_t {
    UnsupportedSchema = 1,
    InvalidCoordinateSystem = 2,
    MissingMetadata = 3,
    DuplicateId = 4,
    InvalidId = 5,
    NonFiniteValue = 6,
    InvalidPhysicalValue = 7,
    DegenerateTriangle = 8,
    DuplicateTriangle = 9,
    MissingVertexReference = 10,
    MissingRegionReference = 11,
    MissingMaterialReference = 12,
    InvalidSideRegions = 13,
    InvalidOpening = 14,
    InvalidEnumValue = 15,
    MissingPilotReference = 16,
    InvalidAttachmentTarget = 17,
    DanglingLineAttachment = 18,
    InvalidSuspensionLine = 19,
    InvalidSeam = 20,
};

struct ValidationDiagnostic {
    DiagnosticSeverity severity = DiagnosticSeverity::Error;
    ValidationCode code = ValidationCode::MissingMetadata;
    EntityKind entityKind = EntityKind::Scene;
    StableId entityId = invalidStableId;
    std::string message;
};

struct ValidationReport {
    std::vector<ValidationDiagnostic> diagnostics;

    [[nodiscard]] bool ok() const noexcept;
    [[nodiscard]] std::size_t errorCount() const noexcept;
};

// Diagnostics are sorted by entity kind, stable ID, code, and message. Thus
// validation output is deterministic even if top-level entity vectors are not
// in ID order.
[[nodiscard]] ValidationReport validateScene(const Scene& scene);

// The compact binary stream is explicitly little-endian and emits every
// top-level entity collection in stable-ID order. writeScene validates before
// writing; readScene validates before replacing its output argument.
[[nodiscard]] bool writeScene(const Scene& scene,
                              std::ostream& output,
                              std::string* errorMessage = nullptr);
[[nodiscard]] bool readScene(std::istream& input,
                             Scene& scene,
                             std::string* errorMessage = nullptr);

} // namespace simwing::fsi

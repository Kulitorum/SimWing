#pragma once

#include "fluid/planar_cut_surface.h"
#include "fluid/moving_interface.h"
#include "transfer.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t planarFaceResolvedBridgeVersion = 4;

enum class PlanarFaceCorrespondenceMode : std::uint8_t {
    FixedMaterial = 0,
    RigidNormalTranslation = 1,
};

struct PlanarFaceResolvedBridgeSettings {
    ConservativeTransferSettings transfer;
    PlanarFaceCorrespondenceMode correspondenceMode =
        PlanarFaceCorrespondenceMode::FixedMaterial;
    double geometryToleranceMeters = 1.0e-12;
    double absoluteVelocityToleranceMetersPerSecond = 1.0e-11;
    double relativeVelocityTolerance = 1.0e-11;
    double minimumNormalAlignment = 1.0 - 1.0e-12;
    double minimumOverlapAreaSquareMeters = 1.0e-14;
    double absoluteAreaToleranceSquareMeters = 1.0e-11;
    double relativeAreaTolerance = 1.0e-11;
    double absoluteForceToleranceNewtons = 1.0e-10;
    double relativeForceTolerance = 1.0e-11;
    double absoluteMomentToleranceNewtonMeters = 1.0e-10;
    double relativeMomentTolerance = 1.0e-11;
    double absolutePowerToleranceWatts = 1.0e-10;
    double relativePowerTolerance = 1.0e-11;
};

struct PlanarFaceResolvedBridgeDiagnostics {
    std::uint32_t version = planarFaceResolvedBridgeVersion;
    std::uint64_t fluidSurfaceStableId = 0;
    std::size_t fluidFaceCount = 0;
    std::size_t structureTriangleCount = 0;
    std::size_t overlapPatchCount = 0;
    double fluidAreaSquareMeters = 0.0;
    double referenceStructureAreaSquareMeters = 0.0;
    double areaResidualSquareMeters = 0.0;
    StructureVector3 fluidPressureForceNewtons;
    StructureVector3 structureSurfaceForceNewtons;
    StructureVector3 forceResidualNewtons;
    double forceResidualNormNewtons = 0.0;
    StructureVector3 fluidPressureMomentNewtonMeters;
    StructureVector3 structureSurfaceMomentNewtonMeters;
    StructureVector3 momentResidualNewtonMeters;
    double momentResidualNormNewtonMeters = 0.0;
    double fluidPressurePowerWatts = 0.0;
    double structureSurfacePowerWatts = 0.0;
    double powerResidualWatts = 0.0;
    double maximumFacePowerResidualWatts = 0.0;
    PlanarFaceCorrespondenceMode correspondenceMode =
        PlanarFaceCorrespondenceMode::FixedMaterial;
    double gridPlaneCoordinateMeters = 0.0;
    double physicalPlaneCoordinateMeters = 0.0;
    double normalTranslationFromReferenceMeters = 0.0;
    double maximumRigidPositionResidualMeters = 0.0;
    double maximumRigidVelocityResidualMetersPerSecond = 0.0;
    bool finite = true;

    bool operator==(
        const PlanarFaceResolvedBridgeDiagnostics&) const = default;
};

class PlanarFaceResolvedTransferResult final {
public:
    [[nodiscard]] const ConservativeTransferResult&
    transferResult() const noexcept;
    [[nodiscard]] const PlanarFaceResolvedBridgeDiagnostics&
    diagnostics() const noexcept;

    bool operator==(
        const PlanarFaceResolvedTransferResult&) const = default;

private:
    friend class PlanarFaceResolvedFluidStructureBridge;
    PlanarFaceResolvedTransferResult(
        ConservativeTransferResult transferResult,
        PlanarFaceResolvedBridgeDiagnostics diagnostics);

    ConservativeTransferResult transferResult_;
    PlanarFaceResolvedBridgeDiagnostics diagnostics_;
};

// Planar bridge for one set of axis-aligned MAC tiles.
// Reference triangles are clipped against every tile once. Each overlap area
// and centroid becomes a material barycentric quadrature patch, so later
// samples may carry nonuniform face traction while exact linear-triangle force,
// moment, and power transfer remains available. FixedMaterial retains the
// original exact grid-face binding. RigidNormalTranslation additionally lets
// the grid plane rebase and the physical surface move along its normal, while
// requiring every structural coupling node to be the matching rigid normal
// translation of its reference point. General transverse deformation,
// curved/Eulerian remap, and nonplanar cut-cell motion remain future operators.
class PlanarFaceResolvedFluidStructureBridge final {
public:
    PlanarFaceResolvedFluidStructureBridge(
        const Structure& target,
        std::uint64_t fluidSurfaceStableId,
        std::vector<CouplingSurfaceNodeDefinition> nodes,
        std::vector<CouplingSurfaceTriangleDefinition> triangles,
        std::vector<fluid::MovingInterfaceFaceDiagnostics> referenceFaces,
        const PlanarFaceResolvedBridgeSettings& settings = {});

    [[nodiscard]] std::uint64_t fluidSurfaceStableId() const noexcept;
    [[nodiscard]] const ConservativeSurfaceTransfer& transfer() const noexcept;
    [[nodiscard]] std::size_t overlapPatchCount() const noexcept;
    [[nodiscard]] double referenceAreaSquareMeters() const noexcept;

    [[nodiscard]] PlanarFaceResolvedTransferResult evaluate(
        const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
        std::span<const CouplingNodeKinematics> nodeKinematics) const;

    [[nodiscard]] PlanarFaceResolvedTransferResult evaluateMovingPlane(
        const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
        std::span<const CouplingNodeKinematics> nodeKinematics,
        double physicalPlaneCoordinateMeters) const;

    [[nodiscard]] PlanarFaceResolvedTransferResult evaluateCutSurface(
        const fluid::PlanarCutSurfacePressureDiagnostics& cutSurface,
        std::span<const CouplingNodeKinematics> nodeKinematics) const;

private:
    struct FaceDefinition {
        std::uint64_t minusRegionStableId = 0;
        std::uint64_t plusRegionStableId = 0;
        fluid::GridFaceAxis axis = fluid::GridFaceAxis::X;
        std::size_t i = 0;
        std::size_t j = 0;
        std::size_t k = 0;
        fluid::Vector3 lowerCornerMeters;
        fluid::Vector3 upperCornerMeters;
        double areaSquareMeters = 0.0;
    };

    struct OverlapPatch {
        std::uint64_t stableId = 0;
        std::uint64_t triangleStableId = 0;
        std::size_t faceIndex = 0;
        std::array<double, 3> barycentricCoordinates{};
        double areaSquareMeters = 0.0;
    };

    [[nodiscard]] PlanarFaceResolvedTransferResult evaluateImpl(
        const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
        std::span<const CouplingNodeKinematics> nodeKinematics,
        std::optional<double> physicalPlaneCoordinateMeters) const;

    std::uint64_t fluidSurfaceStableId_ = 0;
    ConservativeSurfaceTransfer transfer_;
    PlanarFaceResolvedBridgeSettings settings_;
    std::vector<FaceDefinition> faces_;
    std::vector<OverlapPatch> overlaps_;
    std::vector<StructureVector3> referenceNodePositions_;
    fluid::GridFaceAxis axis_ = fluid::GridFaceAxis::X;
    double referencePlaneCoordinateMeters_ = 0.0;
    double referenceAreaSquareMeters_ = 0.0;
};

} // namespace simwing::fsi

#pragma once

#include "structure.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t couplingSurfaceTopologyVersion = 1;

struct CouplingSurfaceNodeDefinition {
    std::uint64_t stableId = 0;
    std::size_t structureNode = 0;

    bool operator==(const CouplingSurfaceNodeDefinition&) const = default;
};

struct CouplingSurfaceTriangleDefinition {
    std::uint64_t stableId = 0;
    std::array<std::uint64_t, 3> nodeStableIds{};

    bool operator==(const CouplingSurfaceTriangleDefinition&) const = default;
};

struct CouplingNodeKinematics {
    std::uint64_t stableId = 0;
    StructureVector3 positionMeters;
    StructureVector3 velocityMetersPerSecond;

    bool operator==(const CouplingNodeKinematics&) const = default;
};

struct CouplingTriangleTraction {
    std::uint64_t stableId = 0;
    StructureVector3 tractionPascals;

    bool operator==(const CouplingTriangleTraction&) const = default;
};

// One constant-traction quadrature patch inside a current linear triangle.
// Barycentric coordinates make the geometric correspondence explicit and
// preserve force, moment, and work exactly for linear triangle kinematics.
// Points must be ordered by triangleStableId, then stableId.
struct CouplingTriangleTractionQuadrature {
    std::uint64_t stableId = 0;
    std::uint64_t triangleStableId = 0;
    std::array<double, 3> barycentricCoordinates{};
    double areaSquareMeters = 0.0;
    StructureVector3 tractionPascals;

    bool operator==(
        const CouplingTriangleTractionQuadrature&) const = default;
};

struct CouplingNodeLoad {
    std::uint64_t stableId = 0;
    std::size_t structureNode = 0;
    StructureVector3 forceNewtons;

    bool operator==(const CouplingNodeLoad&) const = default;
};

struct ConservativeTransferSettings {
    StructureVector3 momentReferenceMeters;
    double minimumTriangleAreaSquareMeters = 1.0e-16;
    double minimumQuadratureAreaSquareMeters = 1.0e-18;
    double barycentricTolerance = 1.0e-12;
};

// Surface and nodal ledgers are accumulated independently. Residuals therefore
// diagnose the actual transfer rather than restating one side of the exchange.
struct ConservativeTransferDiagnostics {
    std::size_t nodeCount = 0;
    std::size_t triangleCount = 0;
    std::size_t quadraturePointCount = 0;
    double surfaceAreaSquareMeters = 0.0;
    StructureVector3 momentReferenceMeters;
    StructureVector3 integratedSurfaceForceNewtons;
    StructureVector3 transferredNodalForceNewtons;
    StructureVector3 forceResidualNewtons;
    double forceResidualNormNewtons = 0.0;
    StructureVector3 integratedSurfaceMomentNewtonMeters;
    StructureVector3 transferredNodalMomentNewtonMeters;
    StructureVector3 momentResidualNewtonMeters;
    double momentResidualNormNewtonMeters = 0.0;
    double integratedSurfacePowerWatts = 0.0;
    double transferredNodalPowerWatts = 0.0;
    double powerResidualWatts = 0.0;
    bool finite = true;

    bool operator==(const ConservativeTransferDiagnostics&) const = default;
};

class ConservativeTransferResult final {
public:
    [[nodiscard]] std::uint64_t surfaceFingerprint() const noexcept;
    [[nodiscard]] std::uint64_t targetDefinitionFingerprint() const noexcept;
    [[nodiscard]] std::span<const CouplingNodeLoad> nodeLoads() const noexcept;
    [[nodiscard]] const ConservativeTransferDiagnostics& diagnostics() const noexcept;

    bool operator==(const ConservativeTransferResult&) const = default;

private:
    friend class ConservativeSurfaceTransfer;
    ConservativeTransferResult() = default;

    std::uint64_t surfaceFingerprint_ = 0;
    std::uint64_t targetDefinitionFingerprint_ = 0;
    std::vector<CouplingNodeLoad> nodeLoads_;
    ConservativeTransferDiagnostics diagnostics_;
};

// A canonical stable-ID view of the structural coupling surface. It accepts a
// uniform world-space traction per current linear triangle and transfers the
// exact integrated wrench/work rate through barycentric nodal loads.
class ConservativeSurfaceTransfer final {
public:
    ConservativeSurfaceTransfer(
        const Structure& target,
        std::vector<CouplingSurfaceNodeDefinition> nodes,
        std::vector<CouplingSurfaceTriangleDefinition> triangles);

    [[nodiscard]] std::uint64_t fingerprint() const noexcept;
    [[nodiscard]] std::uint64_t targetDefinitionFingerprint() const noexcept;
    [[nodiscard]] std::span<const CouplingSurfaceNodeDefinition> nodes() const noexcept;
    [[nodiscard]] std::span<const CouplingSurfaceTriangleDefinition>
    triangles() const noexcept;

    [[nodiscard]] std::vector<CouplingNodeKinematics>
    captureKinematics(const Structure& target) const;

    [[nodiscard]] ConservativeTransferResult evaluate(
        std::span<const CouplingNodeKinematics> nodeKinematics,
        std::span<const CouplingTriangleTraction> triangleTractions,
        const ConservativeTransferSettings& settings = {}) const;

    [[nodiscard]] ConservativeTransferResult evaluateQuadrature(
        std::span<const CouplingNodeKinematics> nodeKinematics,
        std::span<const CouplingTriangleTractionQuadrature> quadrature,
        const ConservativeTransferSettings& settings = {}) const;

    // Adds the immutable result to existing pending structural loads. All
    // topology/result bindings and the bounded activation scale are checked
    // before the first load is mutated.
    void addLoadsTo(Structure& target,
                    const ConservativeTransferResult& result,
                    double activationScale = 1.0) const;

private:
    std::uint64_t targetDefinitionFingerprint_ = 0;
    std::uint64_t fingerprint_ = 0;
    std::vector<CouplingSurfaceNodeDefinition> nodes_;
    std::vector<CouplingSurfaceTriangleDefinition> triangles_;
    std::vector<std::array<std::size_t, 3>> triangleNodeIndices_;
};

} // namespace simwing::fsi

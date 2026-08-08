#include "fluid_frame.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::viewer {
namespace {

Vec3d toViewer(const fsi::fluid::Vector3& value) noexcept {
    return {value.x, value.y, value.z};
}

void addGlobalScalar(DiagnosticFrame& frame,
                     std::string name,
                     std::string unit,
                     const double value) {
    frame.scalarFields.push_back({
        std::move(name), std::move(unit),
        FieldAssociation::Global, {value}});
}

} // namespace

DiagnosticFrame buildPeriodicFluidFrame(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const fsi::fluid::MacVelocityField& velocityMetersPerSecond,
    const fsi::fluid::CellScalarField& pressurePascals,
    const fsi::fluid::PeriodicFlowStrangSubcyclingDiagnostics& diagnostics,
    const PeriodicFluidFrameContext& context) {
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)) {
        throw std::invalid_argument(
            "periodic fluid frame fields do not match their grid");
    }
    if (!fsi::fluid::isFinite(velocityMetersPerSecond)
        || !fsi::fluid::isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "periodic fluid frame fields must be finite");
    }
    if (!diagnostics.accepted
        || diagnostics.completedSubstepCount
            != diagnostics.plannedSubstepCount
        || diagnostics.plannedSubstepCount == 0
        || diagnostics.substeps.size()
            != diagnostics.plannedSubstepCount
        || diagnostics.requestedIntervalSeconds <= 0.0
        || context.step == 0
        || !std::isfinite(context.simulationTimeSeconds)
        || context.simulationTimeSeconds < 0.0
        || !std::isfinite(context.densityKgPerCubicMeter)
        || context.densityKgPerCubicMeter <= 0.0) {
        throw std::invalid_argument(
            "periodic fluid frame requires a committed accepted interval");
    }

    DiagnosticFrame frame;
    frame.sceneChecksum = context.sceneChecksum;
    frame.solverCommit = context.solverCommit;
    frame.step = context.step;
    frame.simulationTimeSeconds = context.simulationTimeSeconds;
    frame.timeStepSeconds = diagnostics.requestedIntervalSeconds;
    frame.couplingResiduals.fluid =
        diagnostics.finalDivergenceL2PerSecond;
    frame.conservation.fluidMassKilograms =
        context.densityKgPerCubicMeter
        * grid.cellVolumeCubicMeters()
        * static_cast<double>(grid.cellCount());
    frame.conservation.totalMomentumNewtonSeconds =
        toViewer(diagnostics.momentumAfterNewtonSeconds);
    frame.conservation.totalEnergyJoules =
        diagnostics.kineticEnergyAfterJoules;

    ScalarField pressure{
        "pressure", "Pa", FieldAssociation::Vertex, {}};
    ScalarField speed{
        "speed", "m/s", FieldAssociation::Vertex, {}};
    VectorField velocity{
        "velocity", "m/s", FieldAssociation::Vertex, {}};
    frame.vertices.reserve(grid.cellCount());
    pressure.values.reserve(grid.cellCount());
    speed.values.reserve(grid.cellCount());
    velocity.values.reserve(grid.cellCount());

    const fsi::fluid::GridCellCounts counts = grid.cellCounts();
    const auto xFaces = velocityMetersPerSecond.xFaces();
    const auto yFaces = velocityMetersPerSecond.yFaces();
    const auto zFaces = velocityMetersPerSecond.zFaces();
    const auto pressureValues = pressurePascals.values();
    for (std::size_t k = 0; k < counts.z; ++k) {
        const std::size_t nextK = k + 1 == counts.z ? 0 : k + 1;
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t nextJ = j + 1 == counts.y ? 0 : j + 1;
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t nextI = i + 1 == counts.x ? 0 : i + 1;
                const std::size_t index = grid.cellIndex(i, j, k);
                const fsi::fluid::Vector3 cellVelocity{
                    0.5 * (xFaces[index]
                           + xFaces[grid.cellIndex(nextI, j, k)]),
                    0.5 * (yFaces[index]
                           + yFaces[grid.cellIndex(i, nextJ, k)]),
                    0.5 * (zFaces[index]
                           + zFaces[grid.cellIndex(i, j, nextK)]),
                };
                frame.vertices.push_back({
                    static_cast<std::uint64_t>(index) + 1,
                    toViewer(grid.cellCenterMeters(i, j, k)),
                });
                pressure.values.push_back(pressureValues[index]);
                speed.values.push_back(std::hypot(
                    cellVelocity.x, cellVelocity.y, cellVelocity.z));
                velocity.values.push_back(toViewer(cellVelocity));
            }
        }
    }
    frame.scalarFields.push_back(std::move(pressure));
    frame.scalarFields.push_back(std::move(speed));
    addGlobalScalar(
        frame, "substeps", "1",
        static_cast<double>(diagnostics.plannedSubstepCount));
    addGlobalScalar(
        frame, "maximum outgoing CFL", "1",
        diagnostics.maximumObservedOutgoingCourantNumber);
    addGlobalScalar(
        frame, "maximum diffusion number", "1",
        diagnostics.maximumObservedDiffusionNumber);
    addGlobalScalar(
        frame, "divergence L2", "1/s",
        diagnostics.finalDivergenceL2PerSecond);
    addGlobalScalar(
        frame, "kinetic energy loss", "J",
        diagnostics.totalEnergyLossJoules);
    frame.vectorFields.push_back(std::move(velocity));
    frame.vectorFields.push_back({
        "momentum residual", "N*s", FieldAssociation::Global,
        {toViewer(diagnostics.momentumResidualNewtonSeconds)}});

    ProtocolError error;
    if (!validateFrame(frame, &error)) {
        throw std::invalid_argument(
            "periodic fluid frame is invalid: " + error.message);
    }
    return frame;
}

} // namespace simwing::viewer

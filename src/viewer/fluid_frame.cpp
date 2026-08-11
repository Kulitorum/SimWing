#include "fluid_frame.h"

#include "cfd_slice.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::viewer {
namespace {

Vec3d toViewer(const fsi::fluid::Vector3& value) noexcept {
    return {value.x, value.y, value.z};
}

bool isFinite(const Vec3d& value) noexcept {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

void addGlobalScalar(DiagnosticFrame& frame,
                     std::string name,
                     std::string unit,
                     const double value) {
    frame.scalarFields.push_back({
        std::move(name), std::move(unit),
        FieldAssociation::Global, {value}});
}

std::size_t previous(const std::size_t index,
                     const std::size_t count) noexcept {
    return index == 0 ? count - 1 : index - 1;
}

std::size_t next(const std::size_t index,
                 const std::size_t count) noexcept {
    return index + 1 == count ? 0 : index + 1;
}

} // namespace

PeriodicFluidCellFields buildPeriodicFluidCellFields(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const fsi::fluid::MacVelocityField& velocityMetersPerSecond) {
    if (!velocityMetersPerSecond.matches(grid)
        || !fsi::fluid::isFinite(velocityMetersPerSecond)) {
        throw std::invalid_argument(
            "periodic fluid cell fields require a matching finite velocity");
    }

    PeriodicFluidCellFields result;
    result.velocityMetersPerSecond.reserve(grid.cellCount());
    result.speedMetersPerSecond.reserve(grid.cellCount());
    result.vorticityPerSecond.reserve(grid.cellCount());
    result.vorticityMagnitudePerSecond.reserve(grid.cellCount());
    fsi::fluid::CellScalarField finiteVolumeDivergence(grid);
    fsi::fluid::computeDivergence(
        grid, velocityMetersPerSecond, finiteVolumeDivergence);
    const auto divergenceValues = finiteVolumeDivergence.values();
    if (!std::ranges::all_of(
            divergenceValues,
            [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "periodic fluid divergence must be finite");
    }
    result.divergencePerSecond.assign(
        divergenceValues.begin(), divergenceValues.end());
    result.maximumAbsoluteDivergencePerSecond =
        fsi::fluid::maximumAbsoluteValue(finiteVolumeDivergence);

    const fsi::fluid::GridCellCounts counts = grid.cellCounts();
    const auto xFaces = velocityMetersPerSecond.xFaces();
    const auto yFaces = velocityMetersPerSecond.yFaces();
    const auto zFaces = velocityMetersPerSecond.zFaces();
    for (std::size_t k = 0; k < counts.z; ++k) {
        const std::size_t nextK = next(k, counts.z);
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t nextJ = next(j, counts.y);
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t nextI = next(i, counts.x);
                const std::size_t index = grid.cellIndex(i, j, k);
                const Vec3d cellVelocity{
                    0.5 * (xFaces[index]
                           + xFaces[grid.cellIndex(nextI, j, k)]),
                    0.5 * (yFaces[index]
                           + yFaces[grid.cellIndex(i, nextJ, k)]),
                    0.5 * (zFaces[index]
                           + zFaces[grid.cellIndex(i, j, nextK)]),
                };
                const double speed = std::hypot(
                    cellVelocity.x, cellVelocity.y, cellVelocity.z);
                if (!isFinite(cellVelocity) || !std::isfinite(speed)) {
                    throw std::invalid_argument(
                        "periodic cell velocity diagnostics must be finite");
                }
                result.velocityMetersPerSecond.push_back(cellVelocity);
                result.speedMetersPerSecond.push_back(speed);
            }
        }
    }

    const fsi::fluid::Vector3 spacing = grid.cellSpacingMeters();
    const auto centeredVelocity = [&](const std::size_t i,
                                      const std::size_t j,
                                      const std::size_t k) -> const Vec3d& {
        return result.velocityMetersPerSecond[grid.cellIndex(i, j, k)];
    };
    for (std::size_t k = 0; k < counts.z; ++k) {
        const std::size_t previousK = previous(k, counts.z);
        const std::size_t nextK = next(k, counts.z);
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t previousJ = previous(j, counts.y);
            const std::size_t nextJ = next(j, counts.y);
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t previousI = previous(i, counts.x);
                const std::size_t nextI = next(i, counts.x);
                const Vec3d curl{
                    (centeredVelocity(i, nextJ, k).z
                     - centeredVelocity(i, previousJ, k).z)
                        / (2.0 * spacing.y)
                        - (centeredVelocity(i, j, nextK).y
                           - centeredVelocity(i, j, previousK).y)
                            / (2.0 * spacing.z),
                    (centeredVelocity(i, j, nextK).x
                     - centeredVelocity(i, j, previousK).x)
                        / (2.0 * spacing.z)
                        - (centeredVelocity(nextI, j, k).z
                           - centeredVelocity(previousI, j, k).z)
                            / (2.0 * spacing.x),
                    (centeredVelocity(nextI, j, k).y
                     - centeredVelocity(previousI, j, k).y)
                        / (2.0 * spacing.x)
                        - (centeredVelocity(i, nextJ, k).x
                           - centeredVelocity(i, previousJ, k).x)
                            / (2.0 * spacing.y),
                };
                const double magnitude = std::hypot(curl.x, curl.y, curl.z);
                if (!isFinite(curl) || !std::isfinite(magnitude)) {
                    throw std::invalid_argument(
                        "periodic vorticity diagnostics must be finite");
                }
                result.vorticityPerSecond.push_back(curl);
                result.vorticityMagnitudePerSecond.push_back(magnitude);
                result.maximumVorticityPerSecond = std::max(
                    result.maximumVorticityPerSecond, magnitude);
            }
        }
    }
    return result;
}

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
    PeriodicFluidCellFields cellFields = buildPeriodicFluidCellFields(
        grid, velocityMetersPerSecond);
    ScalarField speed{
        "speed", "m/s", FieldAssociation::Vertex,
        std::move(cellFields.speedMetersPerSecond)};
    ScalarField divergence{
        "divergence", "1/s", FieldAssociation::Vertex,
        std::move(cellFields.divergencePerSecond)};
    ScalarField vorticityMagnitude{
        "vorticity magnitude", "1/s", FieldAssociation::Vertex,
        std::move(cellFields.vorticityMagnitudePerSecond)};
    VectorField velocity{
        "velocity", "m/s", FieldAssociation::Vertex,
        std::move(cellFields.velocityMetersPerSecond)};
    VectorField vorticity{
        "vorticity", "1/s", FieldAssociation::Vertex,
        std::move(cellFields.vorticityPerSecond)};
    frame.vertices.reserve(grid.cellCount());
    pressure.values.reserve(grid.cellCount());

    const fsi::fluid::GridCellCounts counts = grid.cellCounts();
    const auto pressureValues = pressurePascals.values();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                frame.vertices.push_back({
                    static_cast<std::uint64_t>(index) + 1,
                    toViewer(grid.cellCenterMeters(i, j, k)),
                });
                pressure.values.push_back(pressureValues[index]);
            }
        }
    }
    frame.scalarFields.push_back(std::move(pressure));
    frame.scalarFields.push_back(std::move(speed));
    frame.scalarFields.push_back(std::move(divergence));
    frame.scalarFields.push_back(std::move(vorticityMagnitude));
    addGlobalScalar(
        frame, std::string(cfdGridVertexBeginFieldName), "1", 0.0);
    addGlobalScalar(
        frame, std::string(cfdGridCellCountXFieldName), "1",
        static_cast<double>(counts.x));
    addGlobalScalar(
        frame, std::string(cfdGridCellCountYFieldName), "1",
        static_cast<double>(counts.y));
    addGlobalScalar(
        frame, std::string(cfdGridCellCountZFieldName), "1",
        static_cast<double>(counts.z));
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
        frame, "divergence maximum absolute", "1/s",
        cellFields.maximumAbsoluteDivergencePerSecond);
    addGlobalScalar(
        frame, "vorticity maximum", "1/s",
        cellFields.maximumVorticityPerSecond);
    addGlobalScalar(
        frame, "kinetic energy loss", "J",
        diagnostics.totalEnergyLossJoules);
    frame.vectorFields.push_back(std::move(velocity));
    frame.vectorFields.push_back(std::move(vorticity));
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

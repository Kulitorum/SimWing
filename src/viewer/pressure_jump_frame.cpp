#include "pressure_jump_frame.h"

#include "fluid_frame.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::viewer {
namespace {

using fsi::fluid::GridFaceAxis;
using fsi::fluid::GridFacePressureJump;

Vec3d toViewer(const fsi::fluid::Vector3& value) noexcept {
    return {value.x, value.y, value.z};
}

Vec3d interpolate(const Vec3d& first,
                  const Vec3d& second,
                  const double fraction) noexcept {
    return {
        first.x + fraction * (second.x - first.x),
        first.y + fraction * (second.y - first.y),
        first.z + fraction * (second.z - first.z),
    };
}

bool sameGridFace(const GridFacePressureJump& first,
                  const GridFacePressureJump& second) noexcept {
    return first.axis == second.axis
        && first.i == second.i
        && first.j == second.j
        && first.k == second.k;
}

double wrapCoordinate(const double value,
                      const double lower,
                      const double upper) noexcept {
    if (value < lower) return value + (upper - lower);
    if (value >= upper) return value - (upper - lower);
    return value;
}

std::pair<std::size_t, std::size_t> adjacentCells(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const GridFacePressureJump& face) {
    const auto counts = grid.cellCounts();
    const std::size_t plus = grid.cellIndex(face.i, face.j, face.k);
    switch (face.axis) {
    case GridFaceAxis::X:
        return {
            grid.cellIndex(
                face.i == 0 ? counts.x - 1 : face.i - 1,
                face.j, face.k),
            plus};
    case GridFaceAxis::Y:
        return {
            grid.cellIndex(
                face.i,
                face.j == 0 ? counts.y - 1 : face.j - 1,
                face.k),
            plus};
    case GridFaceAxis::Z:
        return {
            grid.cellIndex(
                face.i, face.j,
                face.k == 0 ? counts.z - 1 : face.k - 1),
            plus};
    }
    throw std::invalid_argument("pressure-jump frame has an invalid axis");
}

double aggregateJump(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const fsi::fluid::SharpPressureJumpField& pressureJumps,
    const GridFacePressureJump& face) {
    const std::size_t index = grid.cellIndex(face.i, face.j, face.k);
    switch (face.axis) {
    case GridFaceAxis::X:
        return pressureJumps.xFaceJumpsPascals()[index];
    case GridFaceAxis::Y:
        return pressureJumps.yFaceJumpsPascals()[index];
    case GridFaceAxis::Z:
        return pressureJumps.zFaceJumpsPascals()[index];
    }
    throw std::invalid_argument("pressure-jump frame has an invalid axis");
}

std::array<Vec3d, 4> crossingQuad(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const GridFacePressureJump& face) {
    const auto spacing = grid.cellSpacingMeters();
    const auto lower = grid.lowerMeters();
    const auto upper = grid.upperMeters();
    switch (face.axis) {
    case GridFaceAxis::X: {
        auto center = grid.xFaceCenterMeters(face.i, face.j, face.k);
        center.x = wrapCoordinate(
            center.x + (face.crossingFraction - 0.5) * spacing.x,
            lower.x, upper.x);
        return {{
            {center.x, center.y - 0.5 * spacing.y,
             center.z - 0.5 * spacing.z},
            {center.x, center.y + 0.5 * spacing.y,
             center.z - 0.5 * spacing.z},
            {center.x, center.y + 0.5 * spacing.y,
             center.z + 0.5 * spacing.z},
            {center.x, center.y - 0.5 * spacing.y,
             center.z + 0.5 * spacing.z},
        }};
    }
    case GridFaceAxis::Y: {
        auto center = grid.yFaceCenterMeters(face.i, face.j, face.k);
        center.y = wrapCoordinate(
            center.y + (face.crossingFraction - 0.5) * spacing.y,
            lower.y, upper.y);
        return {{
            {center.x - 0.5 * spacing.x, center.y,
             center.z - 0.5 * spacing.z},
            {center.x - 0.5 * spacing.x, center.y,
             center.z + 0.5 * spacing.z},
            {center.x + 0.5 * spacing.x, center.y,
             center.z + 0.5 * spacing.z},
            {center.x + 0.5 * spacing.x, center.y,
             center.z - 0.5 * spacing.z},
        }};
    }
    case GridFaceAxis::Z: {
        auto center = grid.zFaceCenterMeters(face.i, face.j, face.k);
        center.z = wrapCoordinate(
            center.z + (face.crossingFraction - 0.5) * spacing.z,
            lower.z, upper.z);
        return {{
            {center.x - 0.5 * spacing.x,
             center.y - 0.5 * spacing.y, center.z},
            {center.x + 0.5 * spacing.x,
             center.y - 0.5 * spacing.y, center.z},
            {center.x + 0.5 * spacing.x,
             center.y + 0.5 * spacing.y, center.z},
            {center.x - 0.5 * spacing.x,
             center.y + 0.5 * spacing.y, center.z},
        }};
    }
    }
    throw std::invalid_argument("pressure-jump frame has an invalid axis");
}

Vec3d axisNormal(const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return {1.0, 0.0, 0.0};
    case GridFaceAxis::Y:
        return {0.0, 1.0, 0.0};
    case GridFaceAxis::Z:
        return {0.0, 0.0, 1.0};
    }
    throw std::invalid_argument("pressure-jump frame has an invalid axis");
}

bool finiteProjectionDiagnostics(
    const fsi::fluid::ProjectionDiagnostics& value) noexcept {
    return std::ranges::all_of(
        std::array{
            value.compatibilityDivergencePerSecond,
            value.initialResidualPascalsPerSquareMeter,
            value.finalResidualPascalsPerSquareMeter,
            value.divergenceL2BeforePerSecond,
            value.divergenceL2AfterPerSecond,
            value.divergenceMaximumBeforePerSecond,
            value.divergenceMaximumAfterPerSecond,
            value.kineticEnergyBeforeJoules,
            value.kineticEnergyAfterJoules,
            value.pressureMeanPascals,
            value.pressureJumpSourceCompatibilityPascalsPerSquareMeter,
        },
        [](const double sample) { return std::isfinite(sample); });
}

void addGlobalScalar(DiagnosticFrame& frame,
                     const char* name,
                     const char* unit,
                     const double value) {
    frame.scalarFields.push_back({
        name, unit, FieldAssociation::Global, {value}});
}

} // namespace

DiagnosticFrame buildPressureJumpFrame(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const fsi::fluid::MacVelocityField& velocityMetersPerSecond,
    const fsi::fluid::CellScalarField& pressurePascals,
    const fsi::fluid::SharpPressureJumpField& pressureJumps,
    const fsi::fluid::ProjectionDiagnostics& diagnostics,
    const PressureJumpFrameContext& context) {
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)
        || !pressureJumps.matches(grid)
        || !fsi::fluid::isFinite(velocityMetersPerSecond)
        || !fsi::fluid::isFinite(pressurePascals)
        || !diagnostics.converged
        || diagnostics.pressureJumpFaceCount != pressureJumps.faceCount()
        || !finiteProjectionDiagnostics(diagnostics)
        || context.step == 0
        || !std::isfinite(context.simulationTimeSeconds)
        || context.simulationTimeSeconds < 0.0
        || !std::isfinite(context.timeStepSeconds)
        || !(context.timeStepSeconds > 0.0)
        || !std::isfinite(context.densityKgPerCubicMeter)
        || !(context.densityKgPerCubicMeter > 0.0)) {
        throw std::invalid_argument(
            "pressure-jump frame requires matching finite accepted state");
    }
    const ProtocolLimits protocolLimits;
    if (grid.cellCount() > protocolLimits.maxVertices
        || pressureJumps.faceCount()
            > (protocolLimits.maxVertices - grid.cellCount()) / 4
        || pressureJumps.faceCount()
            > protocolLimits.maxTriangles / 2) {
        throw std::length_error(
            "pressure-jump frame geometry exceeds protocol limits");
    }

    DiagnosticFrame frame;
    frame.sceneChecksum = context.sceneChecksum;
    frame.solverCommit = context.solverCommit;
    frame.step = context.step;
    frame.simulationTimeSeconds = context.simulationTimeSeconds;
    frame.timeStepSeconds = context.timeStepSeconds;
    frame.couplingResiduals.fluid =
        diagnostics.divergenceL2AfterPerSecond;
    frame.conservation.fluidMassKilograms =
        context.densityKgPerCubicMeter
        * grid.cellVolumeCubicMeters()
        * static_cast<double>(grid.cellCount());
    frame.conservation.totalEnergyJoules =
        diagnostics.kineticEnergyAfterJoules;

    const PeriodicFluidCellFields cellFields =
        buildPeriodicFluidCellFields(grid, velocityMetersPerSecond);
    ScalarField pressure{
        "pressure sample", "Pa", FieldAssociation::Vertex, {}};
    ScalarField speed{
        "speed", "m/s", FieldAssociation::Vertex, {}};
    ScalarField divergence{
        "divergence", "1/s", FieldAssociation::Vertex, {}};
    VectorField velocity{
        "velocity", "m/s", FieldAssociation::Vertex, {}};
    ScalarField jump{
        "pressure jump", "Pa", FieldAssociation::Triangle, {}};
    ScalarField fraction{
        "crossing fraction", "1", FieldAssociation::Triangle, {}};
    VectorField normal{
        "interface normal", "1", FieldAssociation::Triangle, {}};
    const std::size_t expectedVertices = grid.cellCount()
        + 4 * pressureJumps.faceCount();
    const std::size_t expectedTriangles = 2 * pressureJumps.faceCount();
    frame.vertices.reserve(expectedVertices);
    frame.triangles.reserve(expectedTriangles);
    pressure.values.reserve(expectedVertices);
    speed.values.reserve(expectedVertices);
    divergence.values.reserve(expectedVertices);
    velocity.values.reserve(expectedVertices);
    jump.values.reserve(expectedTriangles);
    fraction.values.reserve(expectedTriangles);
    normal.values.reserve(expectedTriangles);

    Vec3d momentum;
    const double cellMass = context.densityKgPerCubicMeter
        * grid.cellVolumeCubicMeters();
    const auto pressureValues = pressurePascals.values();
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                frame.vertices.push_back({
                    static_cast<std::uint64_t>(index) + 1,
                    toViewer(grid.cellCenterMeters(i, j, k)),
                });
                pressure.values.push_back(pressureValues[index]);
                speed.values.push_back(
                    cellFields.speedMetersPerSecond[index]);
                divergence.values.push_back(
                    cellFields.divergencePerSecond[index]);
                velocity.values.push_back(
                    cellFields.velocityMetersPerSecond[index]);
                momentum.x += cellMass
                    * cellFields.velocityMetersPerSecond[index].x;
                momentum.y += cellMass
                    * cellFields.velocityMetersPerSecond[index].y;
                momentum.z += cellMass
                    * cellFields.velocityMetersPerSecond[index].z;
            }
        }
    }
    frame.conservation.totalMomentumNewtonSeconds = {
        momentum.x, momentum.y, momentum.z};

    const auto faces = pressureJumps.faces();
    std::size_t uniqueFaceCount = 0;
    double cumulativeJump = 0.0;
    double smoothMinusPressure = 0.0;
    double smoothDelta = 0.0;
    std::size_t minusCell = 0;
    std::size_t plusCell = 0;
    for (std::size_t crossing = 0; crossing < faces.size(); ++crossing) {
        const auto& face = faces[crossing];
        if (crossing == 0 || !sameGridFace(faces[crossing - 1], face)) {
            ++uniqueFaceCount;
            cumulativeJump = 0.0;
            const auto adjacent = adjacentCells(grid, face);
            minusCell = adjacent.first;
            plusCell = adjacent.second;
            smoothMinusPressure = pressureValues[minusCell];
            smoothDelta = pressureValues[plusCell]
                - aggregateJump(grid, pressureJumps, face)
                - smoothMinusPressure;
        }
        const double surfacePressure = smoothMinusPressure
            + face.crossingFraction * smoothDelta
            + cumulativeJump + 0.5 * face.pressureJumpPascals;
        const Vec3d surfaceVelocity = interpolate(
            cellFields.velocityMetersPerSecond[minusCell],
            cellFields.velocityMetersPerSecond[plusCell],
            face.crossingFraction);
        const double surfaceSpeed = std::hypot(
            surfaceVelocity.x, surfaceVelocity.y, surfaceVelocity.z);
        const double surfaceDivergence =
            cellFields.divergencePerSecond[minusCell]
            + face.crossingFraction
                * (cellFields.divergencePerSecond[plusCell]
                   - cellFields.divergencePerSecond[minusCell]);
        const auto quad = crossingQuad(grid, face);
        const std::uint32_t firstVertex = static_cast<std::uint32_t>(
            frame.vertices.size());
        for (const Vec3d& position : quad) {
            frame.vertices.push_back({
                static_cast<std::uint64_t>(frame.vertices.size()) + 1,
                position});
            pressure.values.push_back(surfacePressure);
            speed.values.push_back(surfaceSpeed);
            divergence.values.push_back(surfaceDivergence);
            velocity.values.push_back(surfaceVelocity);
        }
        const std::uint64_t firstTriangle =
            static_cast<std::uint64_t>(frame.triangles.size()) + 1;
        frame.triangles.push_back({
            firstTriangle,
            firstVertex, firstVertex + 1, firstVertex + 2,
            face.minusRegionStableId, face.plusRegionStableId});
        frame.triangles.push_back({
            firstTriangle + 1,
            firstVertex, firstVertex + 2, firstVertex + 3,
            face.minusRegionStableId, face.plusRegionStableId});
        for (std::size_t half = 0; half < 2; ++half) {
            jump.values.push_back(face.pressureJumpPascals);
            fraction.values.push_back(face.crossingFraction);
            normal.values.push_back(axisNormal(face.axis));
        }
        cumulativeJump += face.pressureJumpPascals;
    }

    frame.scalarFields.push_back(std::move(pressure));
    frame.scalarFields.push_back(std::move(speed));
    frame.scalarFields.push_back(std::move(divergence));
    frame.scalarFields.push_back(std::move(jump));
    frame.scalarFields.push_back(std::move(fraction));
    addGlobalScalar(frame, "authored crossings", "1",
                    static_cast<double>(pressureJumps.faceCount()));
    addGlobalScalar(frame, "crossed grid faces", "1",
                    static_cast<double>(uniqueFaceCount));
    addGlobalScalar(
        frame, "jump source compatibility", "Pa/m^2",
        diagnostics
            .pressureJumpSourceCompatibilityPascalsPerSquareMeter);
    addGlobalScalar(
        frame, "divergence L2", "1/s",
        diagnostics.divergenceL2AfterPerSecond);
    frame.vectorFields.push_back(std::move(velocity));
    frame.vectorFields.push_back(std::move(normal));

    ProtocolError error;
    if (!validateFrame(frame, &error)) {
        throw std::invalid_argument(
            "pressure-jump frame is invalid: " + error.message);
    }
    return frame;
}

} // namespace simwing::viewer

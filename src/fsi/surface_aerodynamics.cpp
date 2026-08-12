#include "surface_aerodynamics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace simwing::fsi {
namespace {

constexpr std::size_t noCell = std::numeric_limits<std::size_t>::max();
constexpr double minimumArea = 1.0e-14;
constexpr double minimumVolume = 1.0e-12;

[[nodiscard]] Vec3 add(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

[[nodiscard]] Vec3 subtract(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

[[nodiscard]] Vec3 scale(const Vec3& value, const double factor) {
    return {value.x * factor, value.y * factor, value.z * factor};
}

[[nodiscard]] double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]] Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

[[nodiscard]] double norm(const Vec3& value) {
    return std::sqrt(dot(value, value));
}

[[nodiscard]] Vec3 normalized(const Vec3& value, const Vec3& fallback) {
    const double magnitude = norm(value);
    return magnitude > 1.0e-12 ? scale(value, 1.0 / magnitude) : fallback;
}

[[nodiscard]] StructureVector3 toStructure(const Vec3& value) {
    return {value.x, value.y, value.z};
}

[[nodiscard]] Vec3 fromStructure(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

[[nodiscard]] bool finite(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

struct TriangleGeometry {
    Vec3 areaVector;
    Vec3 normal;
    double area = 0.0;
};

[[nodiscard]] TriangleGeometry triangleGeometry(
    const SceneFluidSurfaceTriangle& triangle,
    const SceneFluidSurfaceState& state) {
    const Vec3& a = state.vertices[triangle.vertexIndices[0]].positionMeters;
    const Vec3& b = state.vertices[triangle.vertexIndices[1]].positionMeters;
    const Vec3& c = state.vertices[triangle.vertexIndices[2]].positionMeters;
    const Vec3 areaVector = scale(cross(subtract(b, a), subtract(c, a)), 0.5);
    const double area = norm(areaVector);
    if (!(area > minimumArea) || !std::isfinite(area)) {
        throw std::runtime_error(
            "surface aerodynamics encountered a degenerate moving triangle");
    }
    return {areaVector, scale(areaVector, 1.0 / area), area};
}

[[nodiscard]] double tetrahedronContribution(
    const std::array<std::size_t, 3>& indices,
    const SceneFluidSurfaceState& state) {
    const Vec3& a = state.vertices[indices[0]].positionMeters;
    const Vec3& b = state.vertices[indices[1]].positionMeters;
    const Vec3& c = state.vertices[indices[2]].positionMeters;
    return dot(a, cross(b, c)) / 6.0;
}

[[nodiscard]] std::vector<double> signedRegionVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state) {
    std::vector<double> volumes(surface.regions.size(), 0.0);
    for (const SceneFluidSurfaceTriangle& triangle : surface.triangles) {
        const double contribution = tetrahedronContribution(
            triangle.vertexIndices, state);
        volumes[triangle.negativeSideRegionIndex] += contribution;
        volumes[triangle.positiveSideRegionIndex] -= contribution;
    }
    // Openings are not material, but their authored caps close the control
    // volumes used by the gas law. No cap force is applied to the structure.
    for (const SceneFluidSurfaceOpening& opening : surface.openings) {
        for (const auto& cap : opening.capTriangleVertexIndices) {
            const double contribution = tetrahedronContribution(cap, state);
            volumes[opening.negativeSideRegionIndex] += contribution;
            volumes[opening.positiveSideRegionIndex] -= contribution;
        }
    }
    return volumes;
}

[[nodiscard]] SceneFluidSurfaceState referenceSurfaceState(
    const SceneFluidSurfaceDefinition& surface) {
    SceneFluidSurfaceState state;
    state.definitionFingerprint = surface.fingerprint;
    state.vertices.reserve(surface.vertices.size());
    for (const SceneFluidSurfaceVertex& vertex : surface.vertices) {
        state.vertices.push_back({vertex.id, vertex.referencePositionMeters, {}});
    }
    return state;
}

void validateSettings(const SurfaceAerodynamicsSettings& settings) {
    const auto positive = [](const double value) {
        return std::isfinite(value) && value > 0.0;
    };
    if (!positive(settings.timeStepSeconds)
        || !positive(settings.airDensityKgPerCubicMeter)
        || !positive(settings.ambientPressurePascals)
        || !positive(settings.airTemperatureKelvin)
        || !positive(settings.specificGasConstantJoulesPerKilogramKelvin)
        || !finite(fromStructure(settings.targetWindMetersPerSecond))
        || !std::isfinite(settings.windRampSeconds)
        || settings.windRampSeconds < 0.0
        || !std::isfinite(settings.trimIncidenceRadians)
        || !std::isfinite(settings.liftCurveSlopePerRadian)
        || !std::isfinite(settings.minimumLiftCoefficient)
        || !std::isfinite(settings.maximumLiftCoefficient)
        || settings.minimumLiftCoefficient > settings.maximumLiftCoefficient
        || !std::isfinite(settings.zeroLiftDragCoefficient)
        || settings.zeroLiftDragCoefficient < 0.0
        || !std::isfinite(settings.inducedDragFactor)
        || settings.inducedDragFactor < 0.0
        || !std::isfinite(settings.initialCellPressureDynamicFraction)
        || settings.initialCellPressureDynamicFraction < 0.0
        || !positive(settings.maximumCellPressureDynamicFraction)
        || !std::isfinite(settings.intakeDischargeCoefficient)
        || settings.intakeDischargeCoefficient < 0.0
        || !std::isfinite(settings.crossportDischargeCoefficient)
        || settings.crossportDischargeCoefficient < 0.0
        || !std::isfinite(settings.maximumMassChangeFractionPerStep)
        || settings.maximumMassChangeFractionPerStep <= 0.0
        || settings.maximumMassChangeFractionPerStep > 1.0
        || !std::isfinite(settings.minimumVolumeFraction)
        || settings.minimumVolumeFraction <= 0.0
        || settings.minimumVolumeFraction > 1.0) {
        throw std::invalid_argument("surface aerodynamic settings are invalid");
    }
}

[[nodiscard]] Vec3 averageEdgePosition(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const bool xAxis,
    const bool high) {
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (const SceneFluidSurfaceVertex& vertex : surface.vertices) {
        const double coordinate = xAxis ? vertex.referencePositionMeters.x
                                        : vertex.referencePositionMeters.y;
        minimum = std::min(minimum, coordinate);
        maximum = std::max(maximum, coordinate);
    }
    const double extent = maximum - minimum;
    const double threshold = high ? maximum - 0.08 * extent
                                  : minimum + 0.08 * extent;
    Vec3 sum{};
    std::size_t count = 0;
    for (std::size_t index = 0; index < surface.vertices.size(); ++index) {
        const double coordinate = xAxis
            ? surface.vertices[index].referencePositionMeters.x
            : surface.vertices[index].referencePositionMeters.y;
        if ((high && coordinate >= threshold)
            || (!high && coordinate <= threshold)) {
            sum = add(sum, state.vertices[index].positionMeters);
            ++count;
        }
    }
    return count != 0 ? scale(sum, 1.0 / static_cast<double>(count)) : Vec3{};
}

[[nodiscard]] std::pair<double, Vec3> openingAreaVector(
    const SceneFluidSurfaceOpening& opening,
    const SceneFluidSurfaceState& state) {
    Vec3 areaVector{};
    for (const auto& cap : opening.capTriangleVertexIndices) {
        const Vec3& a = state.vertices[cap[0]].positionMeters;
        const Vec3& b = state.vertices[cap[1]].positionMeters;
        const Vec3& c = state.vertices[cap[2]].positionMeters;
        areaVector = add(
            areaVector,
            scale(cross(subtract(b, a), subtract(c, a)), 0.5));
    }
    return {norm(areaVector), areaVector};
}

} // namespace

SurfaceAerodynamicsModel::SurfaceAerodynamicsModel(
    const SceneFluidSurfaceDefinition& surface,
    SurfaceAerodynamicsSettings settings)
    : surface_(&surface), settings_(settings) {
    validateSceneFluidSurfaceDefinition(surface);
    validateSettings(settings_);

    regionToCell_.assign(surface.regions.size(), noCell);
    const SceneFluidSurfaceState reference = referenceSurfaceState(surface);
    const std::vector<double> volumes = signedRegionVolumes(surface, reference);
    for (std::size_t region = 0; region < surface.regions.size(); ++region) {
        if (surface.regions[region].kind != RegionKind::Cell) {
            continue;
        }
        const double signedVolume = volumes[region];
        if (!std::isfinite(signedVolume)
            || std::abs(signedVolume) <= minimumVolume) {
            throw std::invalid_argument(
                "surface aerodynamic cell has no closed reference volume");
        }
        regionToCell_[region] = cells_.size();
        cells_.push_back(
            {region, surface.regions[region].id, signedVolume,
             signedVolume >= 0.0 ? 1.0 : -1.0});
    }
    if (cells_.empty()) {
        throw std::invalid_argument(
            "surface aerodynamic model needs at least one authored cell");
    }

    externalTriangles_.reserve(surface.triangles.size());
    triangleOrientationSigns_.reserve(surface.triangles.size());
    for (const SceneFluidSurfaceTriangle& triangle : surface.triangles) {
        const bool negativeOutside = surface.regions[
            triangle.negativeSideRegionIndex].kind == RegionKind::Outside;
        const bool positiveOutside = surface.regions[
            triangle.positiveSideRegionIndex].kind == RegionKind::Outside;
        externalTriangles_.push_back(
            triangle.role == SurfaceRole::Skin
            && negativeOutside != positiveOutside);
        const std::size_t negativeCell =
            regionToCell_[triangle.negativeSideRegionIndex];
        const std::size_t positiveCell =
            regionToCell_[triangle.positiveSideRegionIndex];
        double orientationSign = 1.0;
        if (negativeCell != noCell) {
            orientationSign = cells_[negativeCell].orientationSign;
        }
        if (positiveCell != noCell) {
            if (negativeCell != noCell
                && cells_[positiveCell].orientationSign != orientationSign) {
                throw std::invalid_argument(
                    "surface aerodynamic adjacent cells disagree on winding orientation");
            }
            orientationSign = cells_[positiveCell].orientationSign;
        }
        triangleOrientationSigns_.push_back(orientationSign);
    }
}

const SurfaceAerodynamicsSettings& SurfaceAerodynamicsModel::settings() const
    noexcept {
    return settings_;
}

std::uint64_t SurfaceAerodynamicsModel::surfaceDefinitionFingerprint() const
    noexcept {
    return surface_->fingerprint;
}

SurfaceAerodynamicsState SurfaceAerodynamicsModel::initialState(
    const SceneFluidSurfaceState& surfaceState) const {
    validateSceneFluidSurfaceState(*surface_, surfaceState);
    const std::vector<double> signedVolumes = signedRegionVolumes(
        *surface_, surfaceState);
    const double gasFactor = settings_.specificGasConstantJoulesPerKilogramKelvin
        * settings_.airTemperatureKelvin;

    SurfaceAerodynamicsState state;
    state.surfaceDefinitionFingerprint = surface_->fingerprint;
    state.acceptedStepCount = surfaceState.acceptedStepCount;
    state.simulationTimeSeconds = surfaceState.simulationTimeSeconds;
    state.cells.reserve(cells_.size());
    for (const CellTopology& topology : cells_) {
        const double referenceVolume = std::abs(
            topology.referenceSignedVolumeCubicMeters);
        const double volume = topology.orientationSign
            * signedVolumes[topology.regionIndex];
        if (!(volume >= referenceVolume * settings_.minimumVolumeFraction)
            || !std::isfinite(volume)) {
            throw std::invalid_argument(
                "surface aerodynamic initial cell volume is collapsed or inverted");
        }
        const double absolutePressure = settings_.ambientPressurePascals;
        state.cells.push_back(
            {topology.regionId, referenceVolume, volume,
             absolutePressure * volume / gasFactor, 0.0});
    }
    return state;
}

SurfaceAerodynamicsCandidate SurfaceAerodynamicsModel::advance(
    const SurfaceAerodynamicsState& current,
    const SceneFluidSurfaceState& surfaceState) const {
    validateSceneFluidSurfaceState(*surface_, surfaceState);
    if (current.surfaceDefinitionFingerprint != surface_->fingerprint
        || current.acceptedStepCount != surfaceState.acceptedStepCount
        || current.simulationTimeSeconds != surfaceState.simulationTimeSeconds
        || current.cells.size() != cells_.size()) {
        throw std::invalid_argument(
            "surface aerodynamic state does not match the moving surface epoch");
    }
    for (std::size_t index = 0; index < cells_.size(); ++index) {
        if (current.cells[index].regionId != cells_[index].regionId
            || !(current.cells[index].airMassKilograms > 0.0)
            || !std::isfinite(current.cells[index].airMassKilograms)) {
            throw std::invalid_argument(
                "surface aerodynamic cell state identity or mass is invalid");
        }
    }

    const double ramp = settings_.windRampSeconds == 0.0
        ? 1.0
        : std::clamp((current.simulationTimeSeconds + settings_.timeStepSeconds)
                         / settings_.windRampSeconds,
                     0.0, 1.0);
    const Vec3 wind = scale(
        fromStructure(settings_.targetWindMetersPerSecond), ramp);
    Vec3 meanVelocity{};
    for (const SceneFluidSurfaceVertexState& vertex : surfaceState.vertices) {
        meanVelocity = add(meanVelocity, vertex.velocityMetersPerSecond);
    }
    meanVelocity = scale(
        meanVelocity, 1.0 / static_cast<double>(surfaceState.vertices.size()));
    const Vec3 relativeWind = subtract(wind, meanVelocity);
    const double relativeSpeed = norm(relativeWind);
    const Vec3 dragDirection = normalized(relativeWind, {0.0, 1.0, 0.0});

    Vec3 spanDirection = normalized(
        subtract(averageEdgePosition(*surface_, surfaceState, true, true),
                 averageEdgePosition(*surface_, surfaceState, true, false)),
        {1.0, 0.0, 0.0});
    Vec3 chordDirection = subtract(
        averageEdgePosition(*surface_, surfaceState, false, true),
        averageEdgePosition(*surface_, surfaceState, false, false));
    chordDirection = subtract(
        chordDirection, scale(spanDirection, dot(chordDirection, spanDirection)));
    chordDirection = normalized(chordDirection, {0.0, 1.0, 0.0});
    Vec3 referenceLiftDirection = normalized(
        cross(spanDirection, chordDirection), {0.0, 0.0, 1.0});
    if (referenceLiftDirection.z < 0.0) {
        referenceLiftDirection = scale(referenceLiftDirection, -1.0);
    }
    Vec3 liftDirection = normalized(
        cross(spanDirection, dragDirection), referenceLiftDirection);
    if (dot(liftDirection, referenceLiftDirection) < 0.0) {
        liftDirection = scale(liftDirection, -1.0);
    }

    const double angleOfAttack = std::atan2(
        dot(relativeWind, referenceLiftDirection),
        dot(relativeWind, chordDirection)) + settings_.trimIncidenceRadians;
    const double liftCoefficient = std::clamp(
        settings_.liftCurveSlopePerRadian * angleOfAttack,
        settings_.minimumLiftCoefficient,
        settings_.maximumLiftCoefficient);
    const double dragCoefficient = settings_.zeroLiftDragCoefficient
        + settings_.inducedDragFactor * liftCoefficient * liftCoefficient;
    const double dynamicPressure = 0.5
        * settings_.airDensityKgPerCubicMeter * relativeSpeed * relativeSpeed;

    std::vector<TriangleGeometry> geometry;
    geometry.reserve(surface_->triangles.size());
    double projectedAreaSum = 0.0;
    double externalAreaSum = 0.0;
    for (std::size_t index = 0; index < surface_->triangles.size(); ++index) {
        geometry.push_back(triangleGeometry(surface_->triangles[index], surfaceState));
        if (externalTriangles_[index]) {
            projectedAreaSum += geometry.back().area
                * std::abs(dot(geometry.back().normal, liftDirection));
            externalAreaSum += geometry.back().area;
        }
    }
    if (!(externalAreaSum > minimumArea)) {
        throw std::runtime_error(
            "surface aerodynamics found no external skin area");
    }
    const double planformArea = std::max(
        0.5 * projectedAreaSum, 0.25 * externalAreaSum);
    const double lift = dynamicPressure * planformArea * liftCoefficient;
    const double drag = dynamicPressure * planformArea * dragCoefficient;
    const Vec3 aerodynamicForce = add(
        scale(liftDirection, lift), scale(dragDirection, drag));

    SurfaceAerodynamicsCandidate candidate;
    candidate.nextState = current;
    candidate.nextState.acceptedStepCount = current.acceptedStepCount + 1;
    candidate.nextState.simulationTimeSeconds =
        current.simulationTimeSeconds + settings_.timeStepSeconds;

    const std::vector<double> signedVolumes = signedRegionVolumes(
        *surface_, surfaceState);
    const double gasFactor = settings_.specificGasConstantJoulesPerKilogramKelvin
        * settings_.airTemperatureKelvin;
    std::vector<double> absolutePressures(cells_.size(), 0.0);
    std::vector<double> geometricVolumeFractions(cells_.size(), 1.0);
    for (std::size_t index = 0; index < cells_.size(); ++index) {
        auto& next = candidate.nextState.cells[index];
        const double geometricVolume = cells_[index].orientationSign
            * signedVolumes[cells_[index].regionIndex];
        const double minimumCellVolume = next.referenceVolumeCubicMeters
            * settings_.minimumVolumeFraction;
        if (!std::isfinite(geometricVolume)) {
            throw std::runtime_error(
                "surface aerodynamic cell "
                + std::to_string(cells_[index].regionId)
                + " has a non-finite volume at step "
                + std::to_string(current.acceptedStepCount));
        }
        geometricVolumeFractions[index] = geometricVolume
            / next.referenceVolumeCubicMeters;
        // A collapsed ram-air cell remains a valid reduced-order state. The
        // gas law uses a small positive regularized volume so pressure can
        // drive recovery; the signed geometric fraction remains diagnostic.
        next.volumeCubicMeters = std::max(geometricVolume, minimumCellVolume);
        absolutePressures[index] = next.airMassKilograms * gasFactor
            / next.volumeCubicMeters;
    }

    std::vector<double> massFlowRates(cells_.size(), 0.0);
    for (const SceneFluidSurfaceOpening& opening : surface_->openings) {
        const std::size_t negativeCell = regionToCell_[
            opening.negativeSideRegionIndex];
        const std::size_t positiveCell = regionToCell_[
            opening.positiveSideRegionIndex];
        const bool negativeOutside = surface_->regions[
            opening.negativeSideRegionIndex].kind == RegionKind::Outside;
        const bool positiveOutside = surface_->regions[
            opening.positiveSideRegionIndex].kind == RegionKind::Outside;
        const auto [area, areaVector] = openingAreaVector(opening, surfaceState);
        if (!(area > minimumArea)) {
            continue;
        }
        if (negativeOutside != positiveOutside) {
            const std::size_t cell = negativeOutside ? positiveCell : negativeCell;
            if (cell == noCell) {
                continue;
            }
            Vec3 inward = scale(areaVector, 1.0 / area);
            if (!negativeOutside) {
                inward = scale(inward, -1.0);
            }
            inward = scale(inward, cells_[cell].orientationSign);
            const double captureSpeed = std::max(0.0, dot(relativeWind, inward));
            const double reservoirPressure = settings_.ambientPressurePascals
                + 0.5 * settings_.airDensityKgPerCubicMeter
                    * captureSpeed * captureSpeed;
            const double pressureDifference = reservoirPressure
                - absolutePressures[cell];
            const double rateMagnitude = settings_.intakeDischargeCoefficient
                * area * std::sqrt(2.0 * settings_.airDensityKgPerCubicMeter
                                   * std::abs(pressureDifference));
            massFlowRates[cell] += std::copysign(rateMagnitude, pressureDifference);
            candidate.diagnostics.totalOpeningMassFlowKilogramsPerSecond
                += rateMagnitude;
        } else if (negativeCell != noCell && positiveCell != noCell
                   && negativeCell != positiveCell) {
            const double pressureDifference = absolutePressures[negativeCell]
                - absolutePressures[positiveCell];
            const double rateMagnitude = settings_.crossportDischargeCoefficient
                * area * std::sqrt(2.0 * settings_.airDensityKgPerCubicMeter
                                   * std::abs(pressureDifference));
            const double rate = std::copysign(rateMagnitude, pressureDifference);
            massFlowRates[negativeCell] -= rate;
            massFlowRates[positiveCell] += rate;
            candidate.diagnostics.totalOpeningMassFlowKilogramsPerSecond
                += rateMagnitude;
        }
    }

    const Vec3 targetWind = fromStructure(settings_.targetWindMetersPerSecond);
    const double targetDynamicPressure = 0.5
        * settings_.airDensityKgPerCubicMeter * dot(targetWind, targetWind);
    const double bootstrapGauge =
        settings_.initialCellPressureDynamicFraction
        * targetDynamicPressure * ramp;
    const double maximumGauge = std::max(
        bootstrapGauge,
        settings_.maximumCellPressureDynamicFraction
            * std::max(dynamicPressure, 1.0e-3));
    candidate.diagnostics.minimumCellGaugePressurePascals =
        std::numeric_limits<double>::infinity();
    candidate.diagnostics.maximumCellGaugePressurePascals =
        -std::numeric_limits<double>::infinity();
    candidate.diagnostics.minimumCellVolumeFraction =
        std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < cells_.size(); ++index) {
        auto& next = candidate.nextState.cells[index];
        const double proposedChange = massFlowRates[index]
            * settings_.timeStepSeconds;
        const double maximumChange = settings_.maximumMassChangeFractionPerStep
            * next.airMassKilograms;
        next.airMassKilograms += std::clamp(
            proposedChange, -maximumChange, maximumChange);
        double absolutePressure = next.airMassKilograms * gasFactor
            / next.volumeCubicMeters;
        absolutePressure = std::clamp(
            absolutePressure,
            settings_.ambientPressurePascals + bootstrapGauge,
            settings_.ambientPressurePascals + maximumGauge);
        next.airMassKilograms = absolutePressure * next.volumeCubicMeters
            / gasFactor;
        next.gaugePressurePascals = absolutePressure
            - settings_.ambientPressurePascals;
        absolutePressures[index] = absolutePressure;
        candidate.diagnostics.minimumCellGaugePressurePascals = std::min(
            candidate.diagnostics.minimumCellGaugePressurePascals,
            next.gaugePressurePascals);
        candidate.diagnostics.maximumCellGaugePressurePascals = std::max(
            candidate.diagnostics.maximumCellGaugePressurePascals,
            next.gaugePressurePascals);
        candidate.diagnostics.minimumCellVolumeFraction = std::min(
            candidate.diagnostics.minimumCellVolumeFraction,
            geometricVolumeFractions[index]);
    }

    candidate.triangleTractions.reserve(surface_->triangles.size());
    candidate.trianglePressureJumpPascals.reserve(surface_->triangles.size());
    candidate.triangleExternalTractionPascals.reserve(surface_->triangles.size());
    Vec3 pressureForce{};
    const double distributionDenominator = projectedAreaSum > minimumArea
        ? projectedAreaSum : externalAreaSum;
    for (std::size_t index = 0; index < surface_->triangles.size(); ++index) {
        const SceneFluidSurfaceTriangle& triangle = surface_->triangles[index];
        const std::size_t negativeCell = regionToCell_[
            triangle.negativeSideRegionIndex];
        const std::size_t positiveCell = regionToCell_[
            triangle.positiveSideRegionIndex];
        const double negativeGauge = negativeCell == noCell ? 0.0
            : absolutePressures[negativeCell] - settings_.ambientPressurePascals;
        const double positiveGauge = positiveCell == noCell ? 0.0
            : absolutePressures[positiveCell] - settings_.ambientPressurePascals;
        const double pressureJump = negativeGauge - positiveGauge;
        const Vec3 pressureTraction = scale(
            geometry[index].normal,
            pressureJump * triangleOrientationSigns_[index]);
        Vec3 externalTraction{};
        if (externalTriangles_[index]) {
            const double weight = projectedAreaSum > minimumArea
                ? geometry[index].area
                    * std::abs(dot(geometry[index].normal, liftDirection))
                : geometry[index].area;
            externalTraction = scale(
                aerodynamicForce,
                weight / (distributionDenominator * geometry[index].area));
        }
        const Vec3 traction = add(pressureTraction, externalTraction);
        candidate.triangleTractions.push_back(
            {triangle.id, toStructure(traction)});
        candidate.trianglePressureJumpPascals.push_back(pressureJump);
        candidate.triangleExternalTractionPascals.push_back(norm(externalTraction));
        pressureForce = add(
            pressureForce, scale(pressureTraction, geometry[index].area));
        if (!finite(traction)) {
            candidate.diagnostics.finite = false;
        }
    }

    candidate.diagnostics.windMetersPerSecond = toStructure(wind);
    candidate.diagnostics.relativeWindMetersPerSecond = toStructure(relativeWind);
    candidate.diagnostics.aerodynamicForceNewtons = toStructure(aerodynamicForce);
    candidate.diagnostics.pressureForceNewtons = toStructure(pressureForce);
    candidate.diagnostics.dynamicPressurePascals = dynamicPressure;
    candidate.diagnostics.planformAreaSquareMeters = planformArea;
    candidate.diagnostics.angleOfAttackRadians = angleOfAttack;
    candidate.diagnostics.liftCoefficient = liftCoefficient;
    candidate.diagnostics.dragCoefficient = dragCoefficient;
    candidate.diagnostics.liftNewtons = lift;
    candidate.diagnostics.dragNewtons = drag;
    candidate.diagnostics.finite = candidate.diagnostics.finite
        && finite(aerodynamicForce) && finite(pressureForce)
        && std::isfinite(candidate.diagnostics.minimumCellGaugePressurePascals)
        && std::isfinite(candidate.diagnostics.maximumCellGaugePressurePascals)
        && std::isfinite(candidate.diagnostics.minimumCellVolumeFraction);
    if (!candidate.diagnostics.finite) {
        throw std::runtime_error(
            "surface aerodynamic candidate produced non-finite diagnostics");
    }
    return candidate;
}

} // namespace simwing::fsi

#include "amr_external_flow_case.h"

#include "cfd_slice.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace simwing::fsi::amr {
namespace {

std::size_t cellCount(const fluid::GridCellCounts counts) {
    return counts.x * counts.y * counts.z;
}

std::size_t cellIndex(const fluid::GridCellCounts counts,
                      const std::size_t i,
                      const std::size_t j,
                      const std::size_t k) {
    return i + counts.x * (j + counts.y * k);
}

fluid::Vector3 cellSpacing(
    const WindTunnelProjectedCoarseGrid& grid) {
    return {
        (grid.upperMeters.x - grid.lowerMeters.x)
            / static_cast<double>(grid.cellCounts.x),
        (grid.upperMeters.y - grid.lowerMeters.y)
            / static_cast<double>(grid.cellCounts.y),
        (grid.upperMeters.z - grid.lowerMeters.z)
            / static_cast<double>(grid.cellCounts.z),
    };
}

void validateSettings(const ExternalFlowTransportSettings& settings) {
    if (!std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || !std::isfinite(settings.markerPulsePeriodSeconds)
        || !(settings.markerPulsePeriodSeconds > 0.0)
        || !std::isfinite(settings.markerPulseDurationSeconds)
        || !(settings.markerPulseDurationSeconds > 0.0)
        || settings.markerPulseDurationSeconds
            > settings.markerPulsePeriodSeconds) {
        throw std::invalid_argument(
            "AMR external-flow transport settings are invalid");
    }
}

void addGlobalScalar(viewer::DiagnosticFrame& frame,
                     std::string name,
                     std::string unit,
                     const double value) {
    frame.scalarFields.push_back({
        std::move(name), std::move(unit),
        viewer::FieldAssociation::Global, {value}});
}

viewer::Vec3d toViewer(const fluid::Vector3 value) {
    return {value.x, value.y, value.z};
}

} // namespace

struct ExternalFlowTransportCase::Implementation {
    ExternalFlowTransportSettings settings;
    WindTunnelProjectedCoarseGrid projected;
    std::vector<double> marker;
    ExternalFlowTransportDiagnostics diagnostics;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;

    explicit Implementation(ExternalFlowTransportSettings requestedSettings)
        : settings(std::move(requestedSettings)) {
        validateSettings(settings);
        projected = evaluateWindTunnelProjectedCoarseGrid(
            settings.projection);
        if (!projected.diagnostics.accepted
            || projected.velocityMetersPerSecond.size()
                != cellCount(projected.cellCounts)) {
            throw std::runtime_error(
                "AMR external-flow projection did not produce an accepted grid");
        }
        marker.assign(cellCount(projected.cellCounts), 0.0);
        diagnostics.projection = projected.diagnostics;
    }

    double markerValue(const std::vector<double>& source,
                       const long long i,
                       const long long j,
                       const long long k,
                       const double inletMarker) const {
        const auto counts = projected.cellCounts;
        if (j < 0) {
            return inletMarker;
        }
        if (j >= static_cast<long long>(counts.y)) {
            return source[cellIndex(
                counts,
                static_cast<std::size_t>(std::clamp<long long>(
                    i, 0, static_cast<long long>(counts.x) - 1)),
                counts.y - 1,
                static_cast<std::size_t>(std::clamp<long long>(
                    k, 0, static_cast<long long>(counts.z) - 1)))];
        }
        if (i < 0 || i >= static_cast<long long>(counts.x)
            || k < 0 || k >= static_cast<long long>(counts.z)) {
            return 0.0;
        }
        return source[cellIndex(
            counts, static_cast<std::size_t>(i),
            static_cast<std::size_t>(j), static_cast<std::size_t>(k))];
    }

    void advanceMarker() {
        const auto counts = projected.cellCounts;
        const fluid::Vector3 spacing = cellSpacing(projected);
        const double phase = std::fmod(
            simulationTimeSeconds, settings.markerPulsePeriodSeconds);
        const double inletMarker =
            phase < settings.markerPulseDurationSeconds ? 1.0 : 0.0;
        std::vector<double> next(marker.size(), 0.0);
        double maximumCourant = 0.0;
        for (std::size_t k = 0; k < counts.z; ++k) {
            for (std::size_t j = 0; j < counts.y; ++j) {
                for (std::size_t i = 0; i < counts.x; ++i) {
                    const std::size_t index = cellIndex(counts, i, j, k);
                    const fluid::Vector3 velocity =
                        projected.velocityMetersPerSecond[index];
                    const double c = marker[index];
                    const auto sample = [&](const long long di,
                                            const long long dj,
                                            const long long dk) {
                        return markerValue(
                            marker,
                            static_cast<long long>(i) + di,
                            static_cast<long long>(j) + dj,
                            static_cast<long long>(k) + dk,
                            inletMarker);
                    };
                    const double derivativeX = velocity.x >= 0.0
                        ? (c - sample(-1, 0, 0)) / spacing.x
                        : (sample(1, 0, 0) - c) / spacing.x;
                    const double derivativeY = velocity.y >= 0.0
                        ? (c - sample(0, -1, 0)) / spacing.y
                        : (sample(0, 1, 0) - c) / spacing.y;
                    const double derivativeZ = velocity.z >= 0.0
                        ? (c - sample(0, 0, -1)) / spacing.z
                        : (sample(0, 0, 1) - c) / spacing.z;
                    next[index] = c - settings.timeStepSeconds
                        * (velocity.x * derivativeX
                           + velocity.y * derivativeY
                           + velocity.z * derivativeZ);
                    maximumCourant = std::max(
                        maximumCourant,
                        settings.timeStepSeconds
                            * (std::abs(velocity.x) / spacing.x
                               + std::abs(velocity.y) / spacing.y
                               + std::abs(velocity.z) / spacing.z));
                }
            }
        }
        if (!(maximumCourant <= 1.0)) {
            throw std::runtime_error(
                "AMR external-flow marker CFL exceeds the donor-cell bound");
        }

        ExternalFlowTransportDiagnostics nextDiagnostics;
        nextDiagnostics.projection = projected.diagnostics;
        nextDiagnostics.maximumOutgoingCourantNumber = maximumCourant;
        nextDiagnostics.minimumMarker =
            *std::ranges::min_element(next);
        nextDiagnostics.maximumMarker =
            *std::ranges::max_element(next);
        const double cellVolume = spacing.x * spacing.y * spacing.z;
        for (const double value : next) {
            nextDiagnostics.markerIntegralCubicMeters += value * cellVolume;
        }
        nextDiagnostics.finite =
            std::isfinite(nextDiagnostics.maximumOutgoingCourantNumber)
            && std::isfinite(nextDiagnostics.minimumMarker)
            && std::isfinite(nextDiagnostics.maximumMarker)
            && std::isfinite(nextDiagnostics.markerIntegralCubicMeters);
        nextDiagnostics.accepted = nextDiagnostics.projection.accepted
            && nextDiagnostics.finite
            && nextDiagnostics.maximumOutgoingCourantNumber <= 1.0
            && nextDiagnostics.minimumMarker >= -1.0e-12
            && nextDiagnostics.maximumMarker <= 1.0 + 1.0e-12;
        if (!nextDiagnostics.accepted) {
            throw std::runtime_error(
                "AMR external-flow marker advance was rejected");
        }
        marker = std::move(next);
        diagnostics = std::move(nextDiagnostics);
    }

    viewer::DiagnosticFrame frame() const {
        if (!diagnostics.accepted || acceptedStepCount == 0) {
            throw std::logic_error(
                "AMR external-flow frame requires an accepted step");
        }
        const auto counts = projected.cellCounts;
        const fluid::Vector3 spacing = cellSpacing(projected);
        const std::size_t count = cellCount(counts);

        viewer::DiagnosticFrame result;
        result.sceneChecksum = externalFlowTransportCaseChecksum;
        result.solverCommit = externalFlowTransportCaseSolverId;
        result.step = acceptedStepCount;
        result.simulationTimeSeconds = simulationTimeSeconds;
        result.timeStepSeconds = settings.timeStepSeconds;
        result.couplingResiduals.fluid =
            diagnostics.projection.projectedMaximumDivergencePerSecond;
        const double density =
            settings.projection.airDensityKilogramsPerCubicMeter;
        const double domainVolume =
            (projected.upperMeters.x - projected.lowerMeters.x)
            * (projected.upperMeters.y - projected.lowerMeters.y)
            * (projected.upperMeters.z - projected.lowerMeters.z);
        result.conservation.fluidMassKilograms = density * domainVolume;

        std::vector<double> speed(count, 0.0);
        std::vector<double> refinementMask(count, 0.0);
        result.vertices.reserve(count);
        double kineticEnergy = 0.0;
        fluid::Vector3 momentum;
        const double cellVolume = spacing.x * spacing.y * spacing.z;
        for (std::size_t k = 0; k < counts.z; ++k) {
            for (std::size_t j = 0; j < counts.y; ++j) {
                for (std::size_t i = 0; i < counts.x; ++i) {
                    const std::size_t index = cellIndex(counts, i, j, k);
                    const fluid::Vector3 velocity =
                        projected.velocityMetersPerSecond[index];
                    speed[index] = std::hypot(
                        velocity.x, velocity.y, velocity.z);
                    const bool refined = i >= counts.x / 4
                        && i < 3 * counts.x / 4
                        && j >= counts.y / 4
                        && j < 3 * counts.y / 4
                        && k >= counts.z / 4
                        && k < 3 * counts.z / 4;
                    refinementMask[index] = refined ? 1.0 : 0.0;
                    const fluid::Vector3 center{
                        projected.lowerMeters.x
                            + (static_cast<double>(i) + 0.5) * spacing.x,
                        projected.lowerMeters.y
                            + (static_cast<double>(j) + 0.5) * spacing.y,
                        projected.lowerMeters.z
                            + (static_cast<double>(k) + 0.5) * spacing.z,
                    };
                    result.vertices.push_back({
                        static_cast<std::uint64_t>(index) + 1,
                        toViewer(center),
                    });
                    const double cellMass = density * cellVolume;
                    momentum.x += cellMass * velocity.x;
                    momentum.y += cellMass * velocity.y;
                    momentum.z += cellMass * velocity.z;
                    kineticEnergy += 0.5 * cellMass
                        * (velocity.x * velocity.x
                           + velocity.y * velocity.y
                           + velocity.z * velocity.z);
                }
            }
        }
        result.conservation.totalMomentumNewtonSeconds = toViewer(momentum);
        result.conservation.totalEnergyJoules = kineticEnergy;
        result.scalarFields.push_back({
            "external_flow.marker", "1",
            viewer::FieldAssociation::Vertex, marker});
        result.scalarFields.push_back({
            "pressure", "Pa", viewer::FieldAssociation::Vertex,
            projected.pressurePascals});
        result.scalarFields.push_back({
            "speed", "m/s", viewer::FieldAssociation::Vertex,
            std::move(speed)});
        result.scalarFields.push_back({
            "divergence", "1/s", viewer::FieldAssociation::Vertex,
            projected.divergencePerSecond});
        result.scalarFields.push_back({
            "AMR refined coarse cell", "1",
            viewer::FieldAssociation::Vertex, std::move(refinementMask)});
        addGlobalScalar(
            result, std::string(viewer::cfdGridVertexBeginFieldName),
            "1", 0.0);
        addGlobalScalar(
            result, std::string(viewer::cfdGridCellCountXFieldName),
            "1", static_cast<double>(counts.x));
        addGlobalScalar(
            result, std::string(viewer::cfdGridCellCountYFieldName),
            "1", static_cast<double>(counts.y));
        addGlobalScalar(
            result, std::string(viewer::cfdGridCellCountZFieldName),
            "1", static_cast<double>(counts.z));
        addGlobalScalar(
            result, "maximum outgoing CFL", "1",
            diagnostics.maximumOutgoingCourantNumber);
        addGlobalScalar(
            result, "projected divergence maximum", "1/s",
            diagnostics.projection.projectedMaximumDivergencePerSecond);
        addGlobalScalar(
            result, "projection iterations", "1",
            static_cast<double>(diagnostics.projection.solverIterations));
        addGlobalScalar(
            result, "marker integral", "m^3",
            diagnostics.markerIntegralCubicMeters);
        result.vectorFields.push_back({
            "velocity", "m/s", viewer::FieldAssociation::Vertex, {}});
        result.vectorFields.back().values.reserve(count);
        for (const fluid::Vector3 value :
             projected.velocityMetersPerSecond) {
            result.vectorFields.back().values.push_back(toViewer(value));
        }

        viewer::ProtocolError error;
        if (!viewer::validateFrame(result, &error)) {
            throw std::runtime_error(
                "AMR external-flow frame is invalid: " + error.message);
        }
        return result;
    }
};

ExternalFlowTransportCase::ExternalFlowTransportCase(
    ExternalFlowTransportSettings settings)
    : implementation_(
          std::make_unique<Implementation>(std::move(settings))) {}

ExternalFlowTransportCase::~ExternalFlowTransportCase() = default;

viewer::TraceHeader ExternalFlowTransportCase::traceHeader() const {
    return {externalFlowTransportCaseChecksum,
            externalFlowTransportCaseSolverId};
}

viewer::DiagnosticFrame ExternalFlowTransportCase::advance() {
    implementation_->advanceMarker();
    ++implementation_->acceptedStepCount;
    implementation_->simulationTimeSeconds +=
        implementation_->settings.timeStepSeconds;
    return implementation_->frame();
}

const ExternalFlowTransportSettings&
ExternalFlowTransportCase::stepSettings() const noexcept {
    return implementation_->settings;
}

const ExternalFlowTransportDiagnostics&
ExternalFlowTransportCase::diagnostics() const noexcept {
    return implementation_->diagnostics;
}

std::uint64_t ExternalFlowTransportCase::acceptedStepCount() const noexcept {
    return implementation_->acceptedStepCount;
}

double ExternalFlowTransportCase::simulationTimeSeconds() const noexcept {
    return implementation_->simulationTimeSeconds;
}

} // namespace simwing::fsi::amr

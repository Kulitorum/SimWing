#pragma once

#include "fluid/planar_region_fragment_opening_pressure_projection.h"
#include "fluid/planar_region_fragment_opening_resistance.h"

#include <cstddef>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

struct PlanarPressureRegionFragmentOpeningPressureStepSettings {
    PlanarPressureRegionFragmentOpeningPressureProjectionSettings projection;
    bool useAuthoredPressureDrive = false;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureStepSettings&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPressureStepLimits {
    PlanarPressureRegionFragmentOpeningResistanceLimits resistanceLimits;
    PlanarPressureRegionFragmentOpeningPressureProjectionLimits
        projectionLimits;
    std::size_t maximumWorkingBytes = 4096ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionFragmentOpeningPressureStepDiagnostics {
    bool accepted = false;
    bool finite = false;
    bool energyAccepted = false;
    std::uint64_t sourceOpeningFluxFingerprint = 0;
    std::uint64_t resultOpeningFluxFingerprint = 0;
    std::size_t workingStorageBytes = 0;
    PlanarPressureRegionFragmentOpeningResistanceDiagnostics resistance;
    PlanarPressureRegionFragmentOpeningPressureProjectionDiagnostics
        projection;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    double geometryPressureWorkJoules = 0.0;
    double authoredPressureWorkJoules = 0.0;
    double correctionKineticEnergyJoules = 0.0;
    double dissipatedEnergyJoules = 0.0;
    double energyResidualJoules = 0.0;
    double energyToleranceJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureStepDiagnostics&)
        const = default;
};

// One opt-in constrained aperture update. Calibrated resistance first advances
// the predicted material-relative opening samples. When explicitly enabled,
// the same midpoint stage also reroutes each wall's fixed authored pressure
// jump into the corresponding opening-fluid degree. The augmented pressure
// projection then corrects Cartesian and opening velocities together so the
// published endpoint closes moving continuity. The aggregate identity is
//
//   deltaK = authoredPressureWork + geometryPressureWork
//            - correctionKineticEnergy
//            - resistanceDissipation.
//
// All four caller-owned fields are copied before either nested stage and commit
// together only after nested and aggregate acceptance. The step owns no
// authored static-pressure evolution/relaxation, traction application,
// topology rebase, scene coefficient mapping, or production selection.
[[nodiscard]] PlanarPressureRegionFragmentOpeningPressureStepDiagnostics
advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    std::vector<double>& orientedTopologyLinkVelocityMetersPerSecond,
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>&
        openingVelocitySamples,
    PlanarPressureRegionFragmentOpeningFluxState& openingFluxState,
    std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings = {},
    const PlanarPressureRegionFragmentOpeningPressureStepLimits& limits = {});

} // namespace simwing::fsi::fluid

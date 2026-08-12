#pragma once

#include "scene_fluid_surface.h"
#include "transfer.h"

#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr char surfaceAerodynamicsModelId[] =
    "simwing-surface-aerodynamics-v1";

// Fast interactive aerodynamic closure for the moving scene surface. It is
// deliberately a reduced-order model: a bounded wing polar supplies the net
// external force, while explicit cell gas states and authored openings supply
// two-sided pressure traction. It does not claim to be a Navier-Stokes solve.
struct SurfaceAerodynamicsSettings {
    double timeStepSeconds = 1.0 / 60.0;
    double airDensityKgPerCubicMeter = 1.225;
    double ambientPressurePascals = 101325.0;
    double airTemperatureKelvin = 288.15;
    double specificGasConstantJoulesPerKilogramKelvin = 287.05;
    StructureVector3 targetWindMetersPerSecond{0.0, 10.0, 0.0};
    double windRampSeconds = 1.0;

    double trimIncidenceRadians = 0.14;
    double liftCurveSlopePerRadian = 4.8;
    double minimumLiftCoefficient = -0.45;
    double maximumLiftCoefficient = 1.25;
    double zeroLiftDragCoefficient = 0.075;
    double inducedDragFactor = 0.12;

    // This pressure target rises with the wind ramp. It is the explicit
    // bootstrap that gives the membrane a flying shape while intake flow
    // starts, without applying a full-pressure impulse on the first frame.
    double initialCellPressureDynamicFraction = 0.65;
    double maximumCellPressureDynamicFraction = 2.0;
    double intakeDischargeCoefficient = 0.65;
    double crossportDischargeCoefficient = 0.72;
    double maximumMassChangeFractionPerStep = 0.10;
    double minimumVolumeFraction = 0.03;
};

struct SurfacePneumaticCellState {
    StableId regionId = invalidStableId;
    double referenceVolumeCubicMeters = 0.0;
    double volumeCubicMeters = 0.0;
    double airMassKilograms = 0.0;
    double gaugePressurePascals = 0.0;

    bool operator==(const SurfacePneumaticCellState&) const = default;
};

struct SurfaceAerodynamicsState {
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::vector<SurfacePneumaticCellState> cells;

    bool operator==(const SurfaceAerodynamicsState&) const = default;
};

struct SurfaceAerodynamicsDiagnostics {
    StructureVector3 windMetersPerSecond;
    StructureVector3 relativeWindMetersPerSecond;
    StructureVector3 aerodynamicForceNewtons;
    StructureVector3 pressureForceNewtons;
    double dynamicPressurePascals = 0.0;
    double planformAreaSquareMeters = 0.0;
    double angleOfAttackRadians = 0.0;
    double liftCoefficient = 0.0;
    double dragCoefficient = 0.0;
    double liftNewtons = 0.0;
    double dragNewtons = 0.0;
    double minimumCellGaugePressurePascals = 0.0;
    double maximumCellGaugePressurePascals = 0.0;
    double minimumCellVolumeFraction = 1.0;
    double totalOpeningMassFlowKilogramsPerSecond = 0.0;
    bool finite = true;

    bool operator==(const SurfaceAerodynamicsDiagnostics&) const = default;
};

struct SurfaceAerodynamicsCandidate {
    SurfaceAerodynamicsState nextState;
    std::vector<CouplingTriangleTraction> triangleTractions;
    std::vector<double> trianglePressureJumpPascals;
    std::vector<double> triangleExternalTractionPascals;
    SurfaceAerodynamicsDiagnostics diagnostics;

    bool operator==(const SurfaceAerodynamicsCandidate&) const = default;
};

class SurfaceAerodynamicsModel final {
public:
    explicit SurfaceAerodynamicsModel(
        const SceneFluidSurfaceDefinition& surface,
        SurfaceAerodynamicsSettings settings = {});

    [[nodiscard]] const SurfaceAerodynamicsSettings& settings() const noexcept;
    [[nodiscard]] std::uint64_t surfaceDefinitionFingerprint() const noexcept;
    [[nodiscard]] SurfaceAerodynamicsState initialState(
        const SceneFluidSurfaceState& surfaceState) const;
    [[nodiscard]] SurfaceAerodynamicsCandidate advance(
        const SurfaceAerodynamicsState& current,
        const SceneFluidSurfaceState& surfaceState) const;

private:
    struct CellTopology {
        std::size_t regionIndex = 0;
        StableId regionId = invalidStableId;
        double referenceSignedVolumeCubicMeters = 0.0;
        double orientationSign = 1.0;
    };

    const SceneFluidSurfaceDefinition* surface_ = nullptr;
    SurfaceAerodynamicsSettings settings_;
    std::vector<CellTopology> cells_;
    std::vector<std::size_t> regionToCell_;
    std::vector<bool> externalTriangles_;
    std::vector<double> triangleOrientationSigns_;
};

} // namespace simwing::fsi

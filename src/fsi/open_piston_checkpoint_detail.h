#pragma once

#include "open_piston_case.h"

namespace simwing::fsi {

struct OpenPistonCaseCheckpoint::Detail {
    StructureCheckpoint structure;
    fluid::MovingInterfaceFluidCheckpoint fluid;
    fluid::PlanarControlVolumeDiagnostics controlVolumeDiagnostics;
    fluid::PlanarControlVolumeRebaseDiagnostics lastRebaseDiagnostics;
    fluid::PlanarCutSurfacePressureDiagnostics cutSurfaceDiagnostics;
    PlanarFaceResolvedBridgeDiagnostics bridgeDiagnostics;
    OpenPistonConservationDiagnostics conservationDiagnostics;
    double lastRebaseVelocityResidualMetersPerSecond = 0.0;
    std::uint64_t topologyRebaseCount = 0;
    std::size_t movingPlaneCoordinate = 0;
    double surfaceOffsetMeters = 0.0;
};

} // namespace simwing::fsi

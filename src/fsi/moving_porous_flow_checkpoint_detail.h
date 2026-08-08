#pragma once

#include "moving_porous_flow_case.h"

#include <stdexcept>

namespace simwing::fsi {

struct MovingPorousFlowCaseCheckpoint::Detail {
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t scalarSampleCount = 0;
    std::size_t pressureJumpCount = 0;
    fluid::MacVelocityField velocityMetersPerSecond;
    fluid::CellScalarField pressurePascals;
    fluid::SharpPressureJumpField pressureJumps;
    fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics diagnostics;
    fluid::MovingPorousFaceTopology porousTopology;
    double sheetPositionMeters = 0.0;
    double sheetVelocityMetersPerSecond = 0.0;
    std::uint64_t topologyRebaseCount = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
};

struct MovingPorousFlowCaseCheckpointCodecAccess {
    [[nodiscard]] static const MovingPorousFlowCaseCheckpoint::Detail& detail(
        const MovingPorousFlowCaseCheckpoint& checkpoint) {
        if (!checkpoint.detail) {
            throw std::invalid_argument(
                "moving porous-flow checkpoint payload is absent");
        }
        return *checkpoint.detail;
    }
};

} // namespace simwing::fsi

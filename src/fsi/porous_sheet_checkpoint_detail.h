#pragma once

#include "porous_sheet_case.h"

#include <memory>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {

struct CoupledPorousSheetCheckpoint::Detail {
    StructureCheckpoint structure;
    fluid::MacVelocityField velocity;
    fluid::CellScalarField pressure;
    CoupledPorousSheetStepDiagnostics diagnostics;

    Detail(StructureCheckpoint structureValue,
           fluid::MacVelocityField velocityValue,
           fluid::CellScalarField pressureValue,
           CoupledPorousSheetStepDiagnostics diagnosticsValue)
        : structure(std::move(structureValue)),
          velocity(std::move(velocityValue)),
          pressure(std::move(pressureValue)),
          diagnostics(std::move(diagnosticsValue)) {}
};

struct CoupledPorousSheetCheckpointCodecAccess {
    [[nodiscard]] static const CoupledPorousSheetCheckpoint::Detail& detail(
        const CoupledPorousSheetCheckpoint& checkpoint) {
        if (!checkpoint.detail) {
            throw std::invalid_argument(
                "coupled porous sheet checkpoint payload is absent");
        }
        return *checkpoint.detail;
    }
};

} // namespace simwing::fsi

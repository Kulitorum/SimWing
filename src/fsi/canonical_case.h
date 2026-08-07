#pragma once

#include "structure.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr char canonicalCaseChecksum[] =
    "sha256:simwing-canonical-structural-case-v1";
inline constexpr char canonicalCaseSolverId[] =
    "simwing-fsi-canonical-worker-v1";

// A deliberately small, deterministic structural case used by the first FSI
// worker/viewer vertical slice. It has one orthotropic membrane plus distance,
// cable, and suspension-tie constraints. Its loads are analytic constants;
// there is no inherited Playground or aerodynamic input.
class CanonicalStructuralCase final {
public:
    CanonicalStructuralCase();

    CanonicalStructuralCase(const CanonicalStructuralCase&) = delete;
    CanonicalStructuralCase& operator=(const CanonicalStructuralCase&) = delete;
    CanonicalStructuralCase(CanonicalStructuralCase&&) = delete;
    CanonicalStructuralCase& operator=(CanonicalStructuralCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;

private:
    Structure structure_;
    viewer::StructureFrameMapping frameMapping_;
    StructureStepSettings stepSettings_;
};

} // namespace simwing::fsi

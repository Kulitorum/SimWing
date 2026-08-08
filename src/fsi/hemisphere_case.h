#pragma once

#include "structure.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstddef>

namespace simwing::fsi {

inline constexpr char anchoredHemisphereCaseChecksum[] =
    "sha256:simwing-three-point-anchored-hemisphere-case-v2";
inline constexpr char anchoredHemisphereCaseSolverId[] =
    "simwing-fsi-three-point-anchored-hemisphere-worker-v2";
inline constexpr std::size_t hemisphereRadialSegments = 36;
inline constexpr std::size_t hemisphereLatitudeSegments = 10;
inline constexpr std::size_t hemisphereAnchorCount = 3;
static_assert(hemisphereRadialSegments % hemisphereAnchorCount == 0);

// A pressure-loaded fabric hemisphere held by three equally spaced points on
// its equator. The free rim exposes the dynamic membrane modes hidden by the
// original fully clamped boundary. Pressure traction is integrated on the
// accepted surface geometry and conservatively shared by each triangle's three
// nodes; it remains an analytic structural case rather than a CFD claim.
class AnchoredHemisphereCase final {
public:
    AnchoredHemisphereCase();

    AnchoredHemisphereCase(const AnchoredHemisphereCase&) = delete;
    AnchoredHemisphereCase& operator=(const AnchoredHemisphereCase&) = delete;
    AnchoredHemisphereCase(AnchoredHemisphereCase&&) = delete;
    AnchoredHemisphereCase& operator=(AnchoredHemisphereCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] double pressurePascals() const noexcept;
    [[nodiscard]] double apexRadialDisplacementMeters() const;
    [[nodiscard]] double maximumFreeRimDisplacementMeters() const;

private:
    Structure structure_;
    viewer::StructureFrameMapping frameMapping_;
    StructureStepSettings stepSettings_;
    double pressurePascals_ = 0.0;
};

} // namespace simwing::fsi

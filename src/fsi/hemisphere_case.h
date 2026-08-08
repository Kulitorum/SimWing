#pragma once

#include "structure.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <array>
#include <cstddef>

namespace simwing::fsi {

inline constexpr char anchoredHemisphereCaseChecksum[] =
    "sha256:simwing-three-point-flapping-hemisphere-case-v4";
inline constexpr char anchoredHemisphereCaseSolverId[] =
    "simwing-fsi-three-point-flapping-hemisphere-worker-v4";
inline constexpr std::size_t hemisphereRadialSegments = 36;
inline constexpr std::size_t hemisphereLatitudeSegments = 10;
inline constexpr std::array<std::size_t, 3> hemisphereAnchorLongitudes{
    0, 12, 24};
inline constexpr std::size_t hemisphereAnchorCount =
    hemisphereAnchorLongitudes.size();
// A finite 2 kN/m segment stiffness keeps the diagnostic outline visible
// without imposing the old exact-length (zero-compliance) rim.
inline constexpr double hemisphereRimComplianceMetersPerNewton = 5.0e-4;

[[nodiscard]] constexpr bool isHemisphereAnchorLongitude(
    const std::size_t longitude) noexcept {
    for (const std::size_t anchor : hemisphereAnchorLongitudes) {
        if (longitude % hemisphereRadialSegments == anchor) {
            return true;
        }
    }
    return false;
}

// A fabric hemisphere held by three points on its compliant rim. A spatially
// alternating follower-pressure mode drives actual membrane deformation while
// the non-collinear supports remove the rigid-body rotation admitted by only
// two positional pins. Pressure traction is integrated on the accepted surface
// geometry and conservatively shared by each triangle's three nodes; it
// remains an analytic structural case rather than a CFD claim.
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

#pragma once

#include "structure.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <array>
#include <cstddef>

namespace simwing::fsi {

inline constexpr char anchoredHemisphereCaseChecksum[] =
    "sha256:simwing-two-point-anchored-hemisphere-case-v3";
inline constexpr char anchoredHemisphereCaseSolverId[] =
    "simwing-fsi-two-point-anchored-hemisphere-worker-v3";
inline constexpr std::size_t hemisphereRadialSegments = 36;
inline constexpr std::size_t hemisphereLatitudeSegments = 10;
// Retain two points from the previous three-point set and release the third.
// The unequal arcs deliberately leave an asymmetric, visibly flapping case.
inline constexpr std::array<std::size_t, 2> hemisphereAnchorLongitudes{0, 12};
inline constexpr std::size_t hemisphereAnchorCount =
    hemisphereAnchorLongitudes.size();
// A finite 10 kN/m segment stiffness keeps the diagnostic outline visible
// without imposing the old exact-length (zero-compliance) rim.
inline constexpr double hemisphereRimComplianceMetersPerNewton = 1.0e-4;

[[nodiscard]] constexpr bool isHemisphereAnchorLongitude(
    const std::size_t longitude) noexcept {
    for (const std::size_t anchor : hemisphereAnchorLongitudes) {
        if (longitude % hemisphereRadialSegments == anchor) {
            return true;
        }
    }
    return false;
}

// A pressure-loaded fabric hemisphere held by two points from the original
// three-point set. Its compliant free rim exposes the dynamic membrane modes
// hidden by the original fully clamped boundary. Pressure traction is
// integrated on the accepted surface geometry and conservatively shared by
// each triangle's three nodes; it remains an analytic structural case rather
// than a CFD claim.
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

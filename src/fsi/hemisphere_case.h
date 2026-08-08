#pragma once

#include "structure.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstddef>

namespace simwing::fsi {

inline constexpr char clampedHemisphereCaseChecksum[] =
    "sha256:simwing-clamped-hemisphere-case-v1";
inline constexpr char clampedHemisphereCaseSolverId[] =
    "simwing-fsi-clamped-hemisphere-worker-v1";
inline constexpr std::size_t clampedHemisphereRadialSegments = 32;
inline constexpr std::size_t clampedHemisphereLatitudeSegments = 10;

// A pressure-loaded fabric hemisphere with its complete equatorial ring fixed.
// The case exercises a curved membrane, signed bending hinges, live follower
// pressure, and a many-element immutable viewer frame without introducing a
// synthetic fluid claim. Pressure traction is integrated on the accepted
// surface geometry and conservatively shared by each triangle's three nodes.
class ClampedHemisphereCase final {
public:
    ClampedHemisphereCase();

    ClampedHemisphereCase(const ClampedHemisphereCase&) = delete;
    ClampedHemisphereCase& operator=(const ClampedHemisphereCase&) = delete;
    ClampedHemisphereCase(ClampedHemisphereCase&&) = delete;
    ClampedHemisphereCase& operator=(ClampedHemisphereCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] double pressurePascals() const noexcept;
    [[nodiscard]] double apexRadialDisplacementMeters() const;

private:
    Structure structure_;
    viewer::StructureFrameMapping frameMapping_;
    StructureStepSettings stepSettings_;
    double pressurePascals_ = 0.0;
};

} // namespace simwing::fsi

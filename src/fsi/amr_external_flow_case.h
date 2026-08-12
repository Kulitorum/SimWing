#pragma once

#include "amr_external_flow_projection.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi::amr {

inline constexpr char externalFlowTransportCaseChecksum[] =
    "sha256:simwing-amr-external-flow-transport-v1";
inline constexpr char externalFlowTransportCaseSolverId[] =
    "simwing-amr-external-flow-transport-worker-v1";
inline constexpr char staticWingExternalFlowSolverId[] =
    "simwing-amr-static-wing-direct-forcing-worker-v1";
inline constexpr char staticWingSlabExternalFlowSolverId[] =
    "simwing-amr-static-wing-slab-direct-forcing-worker-v1";
inline constexpr char staticWingFastSlabExternalFlowSolverId[] =
    "simwing-amr-static-wing-fast-slab-preview-worker-v1";

struct ExternalFlowTransportSettings {
    WindTunnelProjectionSettings projection;
    double timeStepSeconds = 0.005;
    double markerPulsePeriodSeconds = 1.2;
    double markerPulseDurationSeconds = 0.3;
    bool clipStaticWingToWindTunnel = false;
    bool approximateFastPreview = false;

    bool operator==(const ExternalFlowTransportSettings&) const = default;
};

inline constexpr std::size_t maximumSpanwiseSlabResolutionScale = 4;

// At scale one, three coarse X cells at the production 0.375 m spacing are
// centred on one authored span station and only the middle X cell is refined.
// Higher scales subdivide the same physical X/Y/Z domain and refinement box,
// and conservatively reduce dt with scale squared for the sharper interface
// forcing velocities. The X faces remain explicit far-field boundaries, so
// this is a rapid local plumbing/interface probe rather than an aerodynamic-
// load result.
[[nodiscard]] ExternalFlowTransportSettings
makeSpanwiseSlabSettings(
    double centreXMeters = 0.0,
    std::size_t resolutionScale = 1);

// Explicitly approximate interactive preset. It retains the requested slab
// grid but performs six forcing/projection passes instead of twelve and uses
// 1.25 times the conservative slab timestep. Its pressure relative tolerance
// is 1e-6 instead of 1e-10. It is a qualitative transport and viewer probe,
// not a pressure/load/convergence result.
[[nodiscard]] ExternalFlowTransportSettings
makeFastSpanwiseSlabPreviewSettings(
    double centreXMeters = 0.0,
    std::size_t resolutionScale = 1);

struct ExternalFlowTransportDiagnostics {
    WindTunnelProjectionDiagnostics projection;
    WindTunnelMomentumStepDiagnostics momentum;
    double maximumOutgoingCourantNumber = 0.0;
    double minimumMarker = 0.0;
    double maximumMarker = 0.0;
    double markerIntegralCubicMeters = 0.0;
    bool finite = false;
    bool accepted = false;

    bool operator==(const ExternalFlowTransportDiagnostics&) const = default;
};

// Visible production-grid bridge. Persistent two-level momentum and a passive
// marker are advanced and published on the averaged-down coarse diagnostic
// grid. The scene overload additionally retains the fixed authoritative
// surface and activates the explicitly labelled cut-cell direct-forcing
// approximation; it does not claim sharp pressure jumps, viscosity,
// turbulence, aerodynamic loads, or structural coupling.
class ExternalFlowTransportCase final {
public:
    explicit ExternalFlowTransportCase(
        ExternalFlowTransportSettings settings = {});
    ExternalFlowTransportCase(
        const Scene& staticWingScene,
        ExternalFlowTransportSettings settings = {});
    ~ExternalFlowTransportCase();

    ExternalFlowTransportCase(const ExternalFlowTransportCase&) = delete;
    ExternalFlowTransportCase& operator=(const ExternalFlowTransportCase&) =
        delete;
    ExternalFlowTransportCase(ExternalFlowTransportCase&&) = delete;
    ExternalFlowTransportCase& operator=(ExternalFlowTransportCase&&) =
        delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();
    [[nodiscard]] const ExternalFlowTransportSettings&
    stepSettings() const noexcept;
    [[nodiscard]] const ExternalFlowTransportDiagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

} // namespace simwing::fsi::amr

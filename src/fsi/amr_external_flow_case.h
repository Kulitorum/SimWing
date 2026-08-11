#pragma once

#include "amr_external_flow_projection.h"
#include "viewer_protocol.h"

#include <cstdint>
#include <memory>

namespace simwing::fsi::amr {

inline constexpr char externalFlowTransportCaseChecksum[] =
    "sha256:simwing-amr-external-flow-transport-v1";
inline constexpr char externalFlowTransportCaseSolverId[] =
    "simwing-amr-external-flow-transport-worker-v1";

struct ExternalFlowTransportSettings {
    WindTunnelProjectionSettings projection;
    double timeStepSeconds = 0.01;
    double markerPulsePeriodSeconds = 1.2;
    double markerPulseDurationSeconds = 0.3;

    bool operator==(const ExternalFlowTransportSettings&) const = default;
};

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

// A first visible production-grid bridge. The two-level AMR velocity is
// projected once, then an explicitly labelled passive marker is advanced on
// its averaged-down coarse diagnostic grid and published as immutable CFD
// slices. Momentum advection/diffusion and wing coupling are deliberately not
// claimed by this case.
class ExternalFlowTransportCase final {
public:
    explicit ExternalFlowTransportCase(
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

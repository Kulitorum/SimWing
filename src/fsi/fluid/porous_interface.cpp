#include "fluid/porous_interface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

void validateResistance(const DarcyForchheimerResistance& resistance) {
    if (!std::isfinite(resistance.linearPascalSecondsPerMeter)
        || !std::isfinite(
            resistance.quadraticPascalSecondsSquaredPerSquareMeter)
        || resistance.linearPascalSecondsPerMeter < 0.0
        || resistance.quadraticPascalSecondsSquaredPerSquareMeter < 0.0
        || (resistance.linearPascalSecondsPerMeter == 0.0
            && resistance.quadraticPascalSecondsSquaredPerSquareMeter
                == 0.0)) {
        throw std::invalid_argument(
            "porous resistance requires a finite nonnegative active coefficient");
    }
}

std::uint8_t axisOrdinal(const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return 0;
    case GridFaceAxis::Y:
        return 1;
    case GridFaceAxis::Z:
        return 2;
    }
    throw std::invalid_argument("porous crossing has an invalid axis");
}

auto canonicalKey(const GridFacePressureJump& crossing) {
    return std::tuple{
        axisOrdinal(crossing.axis),
        crossing.k, crossing.j, crossing.i,
        crossing.crossingFraction,
        crossing.surfaceStableId,
        crossing.minusRegionStableId,
        crossing.plusRegionStableId};
}

double faceNormalVelocity(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& velocity,
    const PorousGridFaceCrossing& crossing) {
    const auto counts = grid.cellCounts();
    static_cast<void>(axisOrdinal(crossing.axis));
    if (crossing.i >= counts.x
        || crossing.j >= counts.y
        || crossing.k >= counts.z) {
        throw std::invalid_argument("porous crossing index is out of range");
    }
    const std::size_t index = grid.cellIndex(
        crossing.i, crossing.j, crossing.k);
    switch (crossing.axis) {
    case GridFaceAxis::X:
        return velocity.xFaces()[index];
    case GridFaceAxis::Y:
        return velocity.yFaces()[index];
    case GridFaceAxis::Z:
        return velocity.zFaces()[index];
    }
    throw std::invalid_argument("porous crossing has an invalid axis");
}

double indexedFaceNormalVelocity(
    const MacVelocityField& velocity,
    const GridFaceAxis axis,
    const std::size_t index) {
    switch (axis) {
    case GridFaceAxis::X:
        return velocity.xFaces()[index];
    case GridFaceAxis::Y:
        return velocity.yFaces()[index];
    case GridFaceAxis::Z:
        return velocity.zFaces()[index];
    }
    throw std::invalid_argument("pressure jump has an invalid axis");
}

double faceArea(const PeriodicCartesianGrid& grid,
                const GridFaceAxis axis) {
    const auto spacing = grid.cellSpacingMeters();
    switch (axis) {
    case GridFaceAxis::X:
        return spacing.y * spacing.z;
    case GridFaceAxis::Y:
        return spacing.x * spacing.z;
    case GridFaceAxis::Z:
        return spacing.x * spacing.y;
    }
    throw std::invalid_argument("porous crossing has an invalid axis");
}

void validateIterationSettings(const PorousIterationSettings& iteration) {
    switch (iteration.constitutiveEvaluation) {
    case PorousConstitutiveEvaluation::Endpoint:
    case PorousConstitutiveEvaluation::Midpoint:
        break;
    default:
        throw std::invalid_argument(
            "porous projection constitutive time is invalid");
    }
    if (!std::isfinite(
            iteration.absoluteNormalVelocityToleranceMetersPerSecond)
        || iteration.absoluteNormalVelocityToleranceMetersPerSecond < 0.0
        || !std::isfinite(iteration.relativeNormalVelocityTolerance)
        || iteration.relativeNormalVelocityTolerance < 0.0
        || (iteration.absoluteNormalVelocityToleranceMetersPerSecond == 0.0
            && iteration.relativeNormalVelocityTolerance == 0.0)) {
        throw std::invalid_argument(
            "porous projection velocity tolerances are invalid");
    }
    if (!std::isfinite(iteration.absolutePressureJumpTolerancePascals)
        || iteration.absolutePressureJumpTolerancePascals < 0.0
        || !std::isfinite(iteration.relativePressureJumpTolerance)
        || iteration.relativePressureJumpTolerance < 0.0
        || (iteration.absolutePressureJumpTolerancePascals == 0.0
            && iteration.relativePressureJumpTolerance == 0.0)) {
        throw std::invalid_argument(
            "porous projection pressure-jump tolerances are invalid");
    }
    if (!std::isfinite(iteration.relaxation)
        || !(iteration.relaxation > 0.0)
        || iteration.relaxation > 1.0
        || iteration.maximumNonlinearIterations == 0) {
        throw std::invalid_argument(
            "porous projection relaxation and iteration bound are invalid");
    }
}

SharpPressureJumpField combinedPressureJumps(
    const PeriodicCartesianGrid& grid,
    const PorousPressureJumpField& porous,
    const SharpPressureJumpField* prescribed) {
    if (prescribed == nullptr || prescribed->empty()) {
        return porous.pressureJumps();
    }
    std::vector<GridFacePressureJump> faces;
    faces.reserve(
        porous.pressureJumps().faceCount() + prescribed->faceCount());
    faces.insert(
        faces.end(),
        porous.pressureJumps().faces().begin(),
        porous.pressureJumps().faces().end());
    faces.insert(
        faces.end(), prescribed->faces().begin(), prescribed->faces().end());
    return SharpPressureJumpField(grid, std::move(faces));
}

void relaxVelocity(const MacVelocityField& candidate,
                   const double relaxation,
                   MacVelocityField& iterate) {
    const double retained = 1.0 - relaxation;
    for (std::size_t face = 0; face < iterate.xFaces().size(); ++face) {
        iterate.xFaces()[face] =
            retained * iterate.xFaces()[face]
            + relaxation * candidate.xFaces()[face];
        iterate.yFaces()[face] =
            retained * iterate.yFaces()[face]
            + relaxation * candidate.yFaces()[face];
        iterate.zFaces()[face] =
            retained * iterate.zFaces()[face]
            + relaxation * candidate.zFaces()[face];
    }
}

void moveToConstitutiveTime(
    const MacVelocityField& original,
    const PorousConstitutiveEvaluation evaluation,
    MacVelocityField& endpointOrEvaluation) {
    if (evaluation == PorousConstitutiveEvaluation::Endpoint) {
        return;
    }
    for (std::size_t face = 0;
         face < endpointOrEvaluation.xFaces().size(); ++face) {
        endpointOrEvaluation.xFaces()[face] = 0.5
            * (original.xFaces()[face]
               + endpointOrEvaluation.xFaces()[face]);
        endpointOrEvaluation.yFaces()[face] = 0.5
            * (original.yFaces()[face]
               + endpointOrEvaluation.yFaces()[face]);
        endpointOrEvaluation.zFaces()[face] = 0.5
            * (original.zFaces()[face]
               + endpointOrEvaluation.zFaces()[face]);
    }
}

void accumulateJumpLedgers(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& constitutiveVelocity,
    const std::span<const GridFacePressureJump> pressureJumps,
    PorousProjectionDiagnostics& diagnostics) {
    for (const auto& jump : pressureJumps) {
        const double force = jump.pressureJumpPascals
            * faceArea(grid, jump.axis);
        const auto index = grid.cellIndex(jump.i, jump.j, jump.k);
        diagnostics.totalPressureJumpPowerToFluidWatts += force
            * indexedFaceNormalVelocity(
                constitutiveVelocity, jump.axis, index);
        switch (jump.axis) {
        case GridFaceAxis::X:
            diagnostics.totalPressureJumpForceOnFluidNewtons.x += force;
            break;
        case GridFaceAxis::Y:
            diagnostics.totalPressureJumpForceOnFluidNewtons.y += force;
            break;
        case GridFaceAxis::Z:
            diagnostics.totalPressureJumpForceOnFluidNewtons.z += force;
            break;
        }
    }
}

void setJumpLedgers(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& constitutiveVelocity,
    const PorousPressureJumpField* porous,
    const SharpPressureJumpField* prescribed,
    const double timeStepSeconds,
    PorousProjectionDiagnostics& diagnostics) {
    diagnostics.totalDissipationWatts = porous == nullptr
        ? 0.0 : porous->totalDissipationWatts();
    diagnostics.totalPorousDissipationJoules =
        diagnostics.totalDissipationWatts * timeStepSeconds;
    diagnostics.totalPressureJumpForceOnFluidNewtons = {};
    diagnostics.totalPressureJumpPowerToFluidWatts = 0.0;
    if (porous != nullptr) {
        accumulateJumpLedgers(
            grid, constitutiveVelocity, porous->pressureJumps().faces(),
            diagnostics);
    }
    if (prescribed != nullptr) {
        accumulateJumpLedgers(
            grid, constitutiveVelocity, prescribed->faces(), diagnostics);
    }
    diagnostics.totalPressureJumpImpulseOnFluidNewtonSeconds = {
        diagnostics.totalPressureJumpForceOnFluidNewtons.x * timeStepSeconds,
        diagnostics.totalPressureJumpForceOnFluidNewtons.y * timeStepSeconds,
        diagnostics.totalPressureJumpForceOnFluidNewtons.z * timeStepSeconds,
    };
    diagnostics.totalPressureJumpWorkToFluidJoules =
        diagnostics.totalPressureJumpPowerToFluidWatts * timeStepSeconds;
}

} // namespace

double porousPressureJumpPascals(
    const DarcyForchheimerResistance& resistance,
    const double relativeNormalVelocityMetersPerSecond) {
    validateResistance(resistance);
    if (!std::isfinite(relativeNormalVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "porous relative normal velocity must be finite");
    }
    const double magnitude = std::abs(
        relativeNormalVelocityMetersPerSecond);
    const double jump = -relativeNormalVelocityMetersPerSecond
        * (resistance.linearPascalSecondsPerMeter
           + resistance.quadraticPascalSecondsSquaredPerSquareMeter
               * magnitude);
    if (!std::isfinite(jump)) {
        throw std::overflow_error("porous pressure jump is not finite");
    }
    return jump;
}

double porousRelativeNormalVelocityMetersPerSecond(
    const DarcyForchheimerResistance& resistance,
    const double pressureJumpPascals) {
    validateResistance(resistance);
    if (!std::isfinite(pressureJumpPascals)) {
        throw std::invalid_argument("porous pressure jump must be finite");
    }
    const double pressureDrop = std::abs(pressureJumpPascals);
    if (pressureDrop == 0.0) return 0.0;

    double speed = 0.0;
    if (resistance.quadraticPascalSecondsSquaredPerSquareMeter == 0.0) {
        speed = pressureDrop / resistance.linearPascalSecondsPerMeter;
    } else if (resistance.linearPascalSecondsPerMeter == 0.0) {
        speed = std::sqrt(pressureDrop)
            / std::sqrt(
                resistance.quadraticPascalSecondsSquaredPerSquareMeter);
    } else {
        const double maximum = std::numeric_limits<double>::max();
        if (resistance.quadraticPascalSecondsSquaredPerSquareMeter
            <= maximum / pressureDrop / 4.0) {
            const double quadraticRoot = 2.0
                * std::sqrt(
                    resistance
                        .quadraticPascalSecondsSquaredPerSquareMeter)
                * std::sqrt(pressureDrop);
            const double root = std::hypot(
                resistance.linearPascalSecondsPerMeter,
                quadraticRoot);
            speed = (pressureDrop / root)
                / (0.5
                   * (1.0
                      + resistance.linearPascalSecondsPerMeter / root));
        } else {
            const double linearOverPressure =
                resistance.linearPascalSecondsPerMeter / pressureDrop;
            const double quadraticOverPressure =
                resistance.quadraticPascalSecondsSquaredPerSquareMeter
                / pressureDrop;
            const double rootOverPressure = std::hypot(
                linearOverPressure,
                2.0 * std::sqrt(quadraticOverPressure));
            speed = (1.0 / rootOverPressure)
                / (0.5 * (1.0
                           + linearOverPressure / rootOverPressure));
        }
    }
    if (!std::isfinite(speed) || !(speed > 0.0)) {
        throw std::overflow_error(
            "porous relative normal velocity is not finite");
    }
    return std::copysign(speed, -pressureJumpPascals);
}

PorousPressureJumpField::PorousPressureJumpField(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& fluidVelocityMetersPerSecond,
    std::vector<PorousGridFaceCrossing> crossings)
    : pressureJumps_(grid) {
    if (!fluidVelocityMetersPerSecond.matches(grid)
        || !isFinite(fluidVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "porous crossings require a matching finite MAC velocity");
    }

    samples_.reserve(crossings.size());
    std::vector<GridFacePressureJump> pressureJumps;
    pressureJumps.reserve(crossings.size());
    for (const auto& crossing : crossings) {
        validateResistance(crossing.resistance);
        if (!std::isfinite(crossing.surfaceNormalVelocityMetersPerSecond)) {
            throw std::invalid_argument(
                "porous surface normal velocity must be finite");
        }
        const double fluidVelocity = faceNormalVelocity(
            grid, fluidVelocityMetersPerSecond, crossing);
        const double relativeVelocity = fluidVelocity
            - crossing.surfaceNormalVelocityMetersPerSecond;
        const double jump = porousPressureJumpPascals(
            crossing.resistance, relativeVelocity);
        const double area = faceArea(grid, crossing.axis);
        const double flow = relativeVelocity * area;
        const double dissipation = -jump * flow;
        if (!std::isfinite(flow)
            || !std::isfinite(dissipation)
            || dissipation < 0.0) {
            throw std::overflow_error(
                "porous flux or dissipation is not finite");
        }
        GridFacePressureJump pressureJump{
            crossing.surfaceStableId,
            crossing.minusRegionStableId,
            crossing.plusRegionStableId,
            crossing.axis,
            crossing.i,
            crossing.j,
            crossing.k,
            jump,
            crossing.crossingFraction,
        };
        pressureJumps.push_back(pressureJump);
        samples_.push_back({
            pressureJump,
            fluidVelocity,
            crossing.surfaceNormalVelocityMetersPerSecond,
            relativeVelocity,
            area,
            flow,
            dissipation,
        });
        totalDissipationWatts_ += dissipation;
        if (!std::isfinite(totalDissipationWatts_)) {
            throw std::overflow_error(
                "total porous dissipation is not finite");
        }
    }

    SharpPressureJumpField canonical(
        grid, std::move(pressureJumps));
    std::sort(samples_.begin(), samples_.end(),
              [](const auto& first, const auto& second) {
                  return canonicalKey(first.pressureJump)
                      < canonicalKey(second.pressureJump);
              });
    const auto canonicalFaces = canonical.faces();
    for (std::size_t index = 0; index < samples_.size(); ++index) {
        if (samples_[index].pressureJump != canonicalFaces[index]) {
            throw std::logic_error(
                "porous crossing canonicalization is inconsistent");
        }
    }
    pressureJumps_ = std::move(canonical);
}

const SharpPressureJumpField&
PorousPressureJumpField::pressureJumps() const noexcept {
    return pressureJumps_;
}

std::span<const PorousGridFaceSample>
PorousPressureJumpField::samples() const noexcept {
    return samples_;
}

double PorousPressureJumpField::totalDissipationWatts() const noexcept {
    return totalDissipationWatts_;
}

namespace {

struct PorousProjectionCandidate {
    ProjectionDiagnostics projection;
    bool finite = false;
};

template<typename Projector>
PorousProjectionDiagnostics projectVelocityWithPorousInterfacesCore(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField* prescribedPressureJumps,
    const PorousIterationSettings& iterationSettings,
    const double timeStepSeconds,
    Projector&& projectCandidate) {
    validateIterationSettings(iterationSettings);
    if (!predictedVelocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)
        || !isFinite(predictedVelocityMetersPerSecond)
        || !isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "porous projection requires matching finite fields");
    }
    if (prescribedPressureJumps != nullptr
        && !prescribedPressureJumps->matches(grid)) {
        throw std::invalid_argument(
            "prescribed porous-projection jumps do not match their grid");
    }

    PorousProjectionDiagnostics diagnostics;
    diagnostics.constitutiveEvaluation =
        iterationSettings.constitutiveEvaluation;
    diagnostics.porousCrossingCount = porousCrossings.size();
    const MacVelocityField originalVelocity =
        predictedVelocityMetersPerSecond;
    if (porousCrossings.empty()) {
        const PorousProjectionCandidate candidate = projectCandidate(
            predictedVelocityMetersPerSecond, pressurePascals,
            prescribedPressureJumps);
        diagnostics.projection = candidate.projection;
        diagnostics.accepted = diagnostics.projection.converged;
        if (diagnostics.accepted) {
            MacVelocityField constitutiveVelocity =
                predictedVelocityMetersPerSecond;
            moveToConstitutiveTime(
                originalVelocity,
                iterationSettings.constitutiveEvaluation,
                constitutiveVelocity);
            setJumpLedgers(
                grid, constitutiveVelocity, nullptr,
                prescribedPressureJumps, timeStepSeconds, diagnostics);
        }
        diagnostics.finite = isFinite(predictedVelocityMetersPerSecond)
            && isFinite(pressurePascals)
            && std::isfinite(diagnostics.totalPressureJumpForceOnFluidNewtons.x)
            && std::isfinite(diagnostics.totalPressureJumpForceOnFluidNewtons.y)
            && std::isfinite(diagnostics.totalPressureJumpForceOnFluidNewtons.z)
            && std::isfinite(diagnostics.totalPressureJumpPowerToFluidWatts)
            && std::isfinite(diagnostics.totalPressureJumpWorkToFluidJoules);
        return diagnostics;
    }

    const CellScalarField originalPressure = pressurePascals;
    MacVelocityField iterateVelocity = originalVelocity;
    CellScalarField pressureWarmStart = originalPressure;

    for (std::size_t iteration = 0;
         iteration < iterationSettings.maximumNonlinearIterations;
         ++iteration) {
        MacVelocityField sampledVelocity = iterateVelocity;
        moveToConstitutiveTime(
            originalVelocity, iterationSettings.constitutiveEvaluation,
            sampledVelocity);
        const PorousPressureJumpField sampled(
            grid, sampledVelocity, porousCrossings);
        const SharpPressureJumpField pressureJumps = combinedPressureJumps(
            grid, sampled, prescribedPressureJumps);
        MacVelocityField candidateVelocity = originalVelocity;
        CellScalarField candidatePressure = pressureWarmStart;
        const PorousProjectionCandidate candidate = projectCandidate(
            candidateVelocity, candidatePressure, &pressureJumps);
        diagnostics.projection = candidate.projection;
        diagnostics.nonlinearIterationCount = iteration + 1;
        if (!diagnostics.projection.converged) {
            diagnostics.samples.assign(
                sampled.samples().begin(), sampled.samples().end());
            setJumpLedgers(
                grid, sampledVelocity, &sampled,
                prescribedPressureJumps, timeStepSeconds, diagnostics);
            diagnostics.finite = candidate.finite
                && isFinite(candidateVelocity)
                && isFinite(candidatePressure);
            return diagnostics;
        }

        MacVelocityField candidateConstitutiveVelocity = candidateVelocity;
        moveToConstitutiveTime(
            originalVelocity, iterationSettings.constitutiveEvaluation,
            candidateConstitutiveVelocity);
        const PorousPressureJumpField endpoint(
            grid, candidateConstitutiveVelocity, porousCrossings);
        const auto sampledFaces = sampled.samples();
        const auto endpointFaces = endpoint.samples();
        double maximumVelocityResidual = 0.0;
        double maximumJumpResidual = 0.0;
        double velocityScale = 0.0;
        double jumpScale = 0.0;
        for (std::size_t face = 0; face < endpointFaces.size(); ++face) {
            maximumVelocityResidual = std::max(
                maximumVelocityResidual,
                std::abs(
                    endpointFaces[face]
                        .relativeNormalVelocityMetersPerSecond
                    - sampledFaces[face]
                        .relativeNormalVelocityMetersPerSecond));
            maximumJumpResidual = std::max(
                maximumJumpResidual,
                std::abs(
                    endpointFaces[face].pressureJump.pressureJumpPascals
                    - sampledFaces[face].pressureJump.pressureJumpPascals));
            velocityScale = std::max({
                velocityScale,
                std::abs(endpointFaces[face]
                             .relativeNormalVelocityMetersPerSecond),
                std::abs(sampledFaces[face]
                             .relativeNormalVelocityMetersPerSecond),
            });
            jumpScale = std::max({
                jumpScale,
                std::abs(endpointFaces[face]
                             .pressureJump.pressureJumpPascals),
                std::abs(sampledFaces[face]
                             .pressureJump.pressureJumpPascals),
            });
        }
        if (iteration == 0) {
            diagnostics
                .initialMaximumNormalVelocityResidualMetersPerSecond =
                    maximumVelocityResidual;
        }
        diagnostics.finalMaximumNormalVelocityResidualMetersPerSecond =
            maximumVelocityResidual;
        diagnostics.finalMaximumPressureJumpResidualPascals =
            maximumJumpResidual;
        diagnostics.samples.assign(
            endpointFaces.begin(), endpointFaces.end());
        setJumpLedgers(
            grid, candidateConstitutiveVelocity, &endpoint,
            prescribedPressureJumps, timeStepSeconds, diagnostics);
        diagnostics.finite = candidate.finite
            && std::isfinite(maximumVelocityResidual)
            && std::isfinite(maximumJumpResidual)
            && std::isfinite(diagnostics.totalDissipationWatts)
            && std::isfinite(diagnostics.totalPorousDissipationJoules)
            && std::isfinite(diagnostics.totalPressureJumpForceOnFluidNewtons.x)
            && std::isfinite(diagnostics.totalPressureJumpForceOnFluidNewtons.y)
            && std::isfinite(diagnostics.totalPressureJumpForceOnFluidNewtons.z)
            && std::isfinite(
                diagnostics.totalPressureJumpImpulseOnFluidNewtonSeconds.x)
            && std::isfinite(
                diagnostics.totalPressureJumpImpulseOnFluidNewtonSeconds.y)
            && std::isfinite(
                diagnostics.totalPressureJumpImpulseOnFluidNewtonSeconds.z)
            && std::isfinite(diagnostics.totalPressureJumpPowerToFluidWatts)
            && std::isfinite(diagnostics.totalPressureJumpWorkToFluidJoules)
            && isFinite(candidateVelocity) && isFinite(candidatePressure);
        const double velocityTolerance =
            iterationSettings
                .absoluteNormalVelocityToleranceMetersPerSecond
            + iterationSettings.relativeNormalVelocityTolerance
                * velocityScale;
        const double jumpTolerance =
            iterationSettings.absolutePressureJumpTolerancePascals
            + iterationSettings.relativePressureJumpTolerance
                * jumpScale;
        if (diagnostics.finite
            && maximumVelocityResidual <= velocityTolerance
            && maximumJumpResidual <= jumpTolerance) {
            diagnostics.accepted = true;
            predictedVelocityMetersPerSecond = std::move(candidateVelocity);
            pressurePascals = std::move(candidatePressure);
            return diagnostics;
        }
        if (!diagnostics.finite) {
            return diagnostics;
        }
        relaxVelocity(
            candidateVelocity, iterationSettings.relaxation,
            iterateVelocity);
        pressureWarmStart = std::move(candidatePressure);
    }
    return diagnostics;
}

PorousProjectionDiagnostics projectVelocityWithPorousInterfacesImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField* prescribedPressureJumps,
    const PorousProjectionSettings& settings) {
    const auto projectCandidate =
        [&](MacVelocityField& candidateVelocity,
            CellScalarField& candidatePressure,
            const SharpPressureJumpField* pressureJumps) {
            const ProjectionDiagnostics projection =
                pressureJumps == nullptr || pressureJumps->empty()
                ? projectVelocity(
                    grid, candidateVelocity, candidatePressure,
                    settings.projection)
                : projectVelocityWithPressureJumps(
                    grid, candidateVelocity, candidatePressure,
                    *pressureJumps, settings.projection);
            return PorousProjectionCandidate{
                .projection = projection,
                .finite = projection.converged,
            };
        };
    return projectVelocityWithPorousInterfacesCore(
        grid, predictedVelocityMetersPerSecond, pressurePascals,
        porousCrossings, prescribedPressureJumps, settings.iteration,
        settings.projection.timeStepSeconds, projectCandidate);
}

MovingPorousProjectionDiagnostics
projectVelocityWithMovingAndPorousInterfacesImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& movingInterfaces,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField* prescribedPressureJumps,
    const MovingPorousProjectionSettings& settings) {
    validateIterationSettings(settings.iteration);
    if (!movingInterfaces.matches(grid)
        || !predictedVelocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)
        || !isFinite(predictedVelocityMetersPerSecond)
        || !isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "moving porous projection requires matching finite inputs");
    }
    if (prescribedPressureJumps != nullptr
        && !prescribedPressureJumps->matches(grid)) {
        throw std::invalid_argument(
            "prescribed moving-porous jumps do not match their grid");
    }

    MovingPorousProjectionDiagnostics diagnostics;
    const auto projectCandidate =
        [&](MacVelocityField& candidateVelocity,
            CellScalarField& candidatePressure,
            const SharpPressureJumpField* pressureJumps) {
            diagnostics.movingInterface =
                pressureJumps == nullptr || pressureJumps->empty()
                ? projectVelocityWithMovingInterfaces(
                    grid, candidateVelocity, candidatePressure,
                    movingInterfaces, settings.movingProjection)
                : projectVelocityWithMovingInterfacesAndPressureJumps(
                    grid, candidateVelocity, candidatePressure,
                    movingInterfaces, *pressureJumps,
                    settings.movingProjection);
            return PorousProjectionCandidate{
                .projection = diagnostics.movingInterface.projection,
                .finite = diagnostics.movingInterface.finite,
            };
        };
    diagnostics.porous = projectVelocityWithPorousInterfacesCore(
        grid, predictedVelocityMetersPerSecond, pressurePascals,
        porousCrossings, prescribedPressureJumps, settings.iteration,
        settings.movingProjection.projection.timeStepSeconds,
        projectCandidate);
    diagnostics.finite = diagnostics.porous.finite
        && diagnostics.movingInterface.finite;
    diagnostics.accepted = diagnostics.porous.accepted
        && diagnostics.finite;
    diagnostics.porous.finite = diagnostics.finite;
    diagnostics.porous.accepted = diagnostics.accepted;
    return diagnostics;
}

} // namespace

PorousProjectionDiagnostics projectVelocityWithPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const PorousProjectionSettings& settings) {
    return projectVelocityWithPorousInterfacesImpl(
        grid, predictedVelocityMetersPerSecond, pressurePascals,
        porousCrossings, nullptr, settings);
}

PorousProjectionDiagnostics projectVelocityWithPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField& prescribedPressureJumps,
    const PorousProjectionSettings& settings) {
    return projectVelocityWithPorousInterfacesImpl(
        grid, predictedVelocityMetersPerSecond, pressurePascals,
        porousCrossings, &prescribedPressureJumps, settings);
}

MovingPorousProjectionDiagnostics
projectVelocityWithMovingAndPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& movingInterfaces,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const MovingPorousProjectionSettings& settings) {
    return projectVelocityWithMovingAndPorousInterfacesImpl(
        grid, predictedVelocityMetersPerSecond, pressurePascals,
        movingInterfaces, porousCrossings, nullptr, settings);
}

MovingPorousProjectionDiagnostics
projectVelocityWithMovingAndPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& movingInterfaces,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField& prescribedPressureJumps,
    const MovingPorousProjectionSettings& settings) {
    return projectVelocityWithMovingAndPorousInterfacesImpl(
        grid, predictedVelocityMetersPerSecond, pressurePascals,
        movingInterfaces, porousCrossings, &prescribedPressureJumps,
        settings);
}

namespace {

bool finiteVector(const Vector3& value) noexcept {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Vector3 addVector(
    const Vector3& first,
    const Vector3& second) noexcept {
    return {
        first.x + second.x,
        first.y + second.y,
        first.z + second.z,
    };
}

Vector3 scaleVector(
    const Vector3& value,
    const double scale) noexcept {
    return {value.x * scale, value.y * scale, value.z * scale};
}

Vector3 axialVector(
    const GridFaceAxis axis,
    const double value) {
    switch (axis) {
    case GridFaceAxis::X:
        return {value, 0.0, 0.0};
    case GridFaceAxis::Y:
        return {0.0, value, 0.0};
    case GridFaceAxis::Z:
        return {0.0, 0.0, value};
    }
    throw std::invalid_argument("porous traction has an invalid face axis");
}

bool closeValue(
    const double first,
    const double second) noexcept {
    return std::abs(first - second)
        <= 1.0e-10
            + 1.0e-11 * std::max(std::abs(first), std::abs(second));
}

bool closeVector(
    const Vector3& first,
    const Vector3& second) noexcept {
    return closeValue(first.x, second.x)
        && closeValue(first.y, second.y)
        && closeValue(first.z, second.z);
}

double wrapCoordinate(
    const double value,
    const double lower,
    const double upper) noexcept {
    if (value < lower) return value + (upper - lower);
    if (value >= upper) return value - (upper - lower);
    return value;
}

std::pair<Vector3, Vector3> crossingBounds(
    const PeriodicCartesianGrid& grid,
    const GridFacePressureJump& crossing) {
    const Vector3 spacing = grid.cellSpacingMeters();
    const Vector3 lower = grid.lowerMeters();
    const Vector3 upper = grid.upperMeters();
    switch (crossing.axis) {
    case GridFaceAxis::X: {
        Vector3 center = grid.xFaceCenterMeters(
            crossing.i, crossing.j, crossing.k);
        center.x = wrapCoordinate(
            center.x + (crossing.crossingFraction - 0.5) * spacing.x,
            lower.x, upper.x);
        return {
            {center.x, center.y - 0.5 * spacing.y,
             center.z - 0.5 * spacing.z},
            {center.x, center.y + 0.5 * spacing.y,
             center.z + 0.5 * spacing.z},
        };
    }
    case GridFaceAxis::Y: {
        Vector3 center = grid.yFaceCenterMeters(
            crossing.i, crossing.j, crossing.k);
        center.y = wrapCoordinate(
            center.y + (crossing.crossingFraction - 0.5) * spacing.y,
            lower.y, upper.y);
        return {
            {center.x - 0.5 * spacing.x, center.y,
             center.z - 0.5 * spacing.z},
            {center.x + 0.5 * spacing.x, center.y,
             center.z + 0.5 * spacing.z},
        };
    }
    case GridFaceAxis::Z: {
        Vector3 center = grid.zFaceCenterMeters(
            crossing.i, crossing.j, crossing.k);
        center.z = wrapCoordinate(
            center.z + (crossing.crossingFraction - 0.5) * spacing.z,
            lower.z, upper.z);
        return {
            {center.x - 0.5 * spacing.x,
             center.y - 0.5 * spacing.y, center.z},
            {center.x + 0.5 * spacing.x,
             center.y + 0.5 * spacing.y, center.z},
        };
    }
    }
    throw std::invalid_argument("porous traction has an invalid face axis");
}

PorousSurfaceTractionDiagnostics evaluatePorousSurfaceTractionImpl(
    const PeriodicCartesianGrid& grid,
    const PorousProjectionDiagnostics& source,
    const double timeStepSeconds) {
    if (!std::isfinite(timeStepSeconds) || !(timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "porous traction requires a finite positive time step");
    }
    if (!source.accepted || !source.finite
        || !source.projection.converged
        || source.porousCrossingCount != source.samples.size()
        || !std::isfinite(source.totalDissipationWatts)
        || !std::isfinite(source.totalPorousDissipationJoules)
        || source.totalDissipationWatts < 0.0
        || source.totalPorousDissipationJoules < 0.0
        || !closeValue(
            source.totalPorousDissipationJoules,
            source.totalDissipationWatts * timeStepSeconds)) {
        throw std::invalid_argument(
            "porous traction requires accepted consistent diagnostics");
    }

    PorousSurfaceTractionDiagnostics result;
    result.timeStepSeconds = timeStepSeconds;
    result.faces.reserve(source.samples.size());
    std::map<std::uint64_t, PorousSurfaceTractionAggregate> surfaces;
    const GridCellCounts counts = grid.cellCounts();
    for (const auto& sample : source.samples) {
        const auto& crossing = sample.pressureJump;
        if (crossing.surfaceStableId == 0
            || crossing.minusRegionStableId == 0
            || crossing.plusRegionStableId == 0
            || crossing.minusRegionStableId == crossing.plusRegionStableId
            || crossing.i >= counts.x
            || crossing.j >= counts.y
            || crossing.k >= counts.z
            || !std::isfinite(crossing.crossingFraction)
            || !(crossing.crossingFraction > 0.0)
            || !(crossing.crossingFraction < 1.0)
            || !std::isfinite(crossing.pressureJumpPascals)
            || !std::isfinite(
                sample.fluidNormalVelocityMetersPerSecond)
            || !std::isfinite(
                sample.surfaceNormalVelocityMetersPerSecond)
            || !std::isfinite(
                sample.relativeNormalVelocityMetersPerSecond)
            || !std::isfinite(sample.faceAreaSquareMeters)
            || !(sample.faceAreaSquareMeters > 0.0)
            || !std::isfinite(
                sample.volumeFlowRateCubicMetersPerSecond)
            || !std::isfinite(sample.dissipationWatts)
            || sample.dissipationWatts < 0.0
            || !closeValue(
                sample.faceAreaSquareMeters,
                faceArea(grid, crossing.axis))
            || !closeValue(
                sample.relativeNormalVelocityMetersPerSecond,
                sample.fluidNormalVelocityMetersPerSecond
                    - sample.surfaceNormalVelocityMetersPerSecond)
            || !closeValue(
                sample.volumeFlowRateCubicMetersPerSecond,
                sample.relativeNormalVelocityMetersPerSecond
                    * sample.faceAreaSquareMeters)
            || !closeValue(
                sample.dissipationWatts,
                -crossing.pressureJumpPascals
                    * sample.volumeFlowRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "porous traction sample is inconsistent");
        }

        const double axialFluidForce =
            crossing.pressureJumpPascals
            * sample.faceAreaSquareMeters;
        const Vector3 fluidForce = axialVector(
            crossing.axis, axialFluidForce);
        const Vector3 fluidImpulse = scaleVector(
            fluidForce, timeStepSeconds);
        const double fluidPower = axialFluidForce
            * sample.fluidNormalVelocityMetersPerSecond;
        const double fluidWork = fluidPower * timeStepSeconds;
        const Vector3 surfaceForce = scaleVector(fluidForce, -1.0);
        const Vector3 surfaceImpulse = scaleVector(
            surfaceForce, timeStepSeconds);
        const double surfacePower = -axialFluidForce
            * sample.surfaceNormalVelocityMetersPerSecond;
        const double surfaceWork = surfacePower * timeStepSeconds;
        const double dissipatedEnergy =
            sample.dissipationWatts * timeStepSeconds;
        const double energyResidual =
            fluidWork + surfaceWork + dissipatedEnergy;
        const auto bounds = crossingBounds(grid, crossing);

        PorousFaceTractionDiagnostics face{
            crossing.surfaceStableId,
            crossing.minusRegionStableId,
            crossing.plusRegionStableId,
            crossing.axis,
            crossing.i,
            crossing.j,
            crossing.k,
            crossing.crossingFraction,
            bounds.first,
            bounds.second,
            sample.faceAreaSquareMeters,
            crossing.pressureJumpPascals,
            sample.fluidNormalVelocityMetersPerSecond,
            sample.surfaceNormalVelocityMetersPerSecond,
            sample.relativeNormalVelocityMetersPerSecond,
            fluidForce,
            fluidImpulse,
            fluidPower,
            fluidWork,
            surfaceForce,
            surfaceImpulse,
            surfacePower,
            surfaceWork,
            sample.dissipationWatts,
            dissipatedEnergy,
            energyResidual,
        };
        if (!finiteVector(face.pressureForceOnFluidNewtons)
            || !finiteVector(face.pressureImpulseOnFluidNewtonSeconds)
            || !finiteVector(face.pressureForceOnSurfaceNewtons)
            || !finiteVector(face.pressureImpulseOnSurfaceNewtonSeconds)
            || !std::isfinite(face.pressurePowerToFluidWatts)
            || !std::isfinite(face.pressureWorkToFluidJoules)
            || !std::isfinite(face.pressurePowerToSurfaceWatts)
            || !std::isfinite(face.pressureWorkToSurfaceJoules)
            || !std::isfinite(face.dissipatedEnergyJoules)
            || !std::isfinite(face.energyResidualJoules)) {
            throw std::overflow_error(
                "porous traction ledger is not finite");
        }
        result.faces.push_back(face);

        auto [surfaceIterator, inserted] =
            surfaces.try_emplace(crossing.surfaceStableId);
        auto& surface = surfaceIterator->second;
        if (inserted) {
            surface.stableId = crossing.surfaceStableId;
        }
        ++surface.faceCount;
        surface.areaSquareMeters += face.areaSquareMeters;
        surface.pressureForceOnFluidNewtons = addVector(
            surface.pressureForceOnFluidNewtons,
            face.pressureForceOnFluidNewtons);
        surface.pressureImpulseOnFluidNewtonSeconds = addVector(
            surface.pressureImpulseOnFluidNewtonSeconds,
            face.pressureImpulseOnFluidNewtonSeconds);
        surface.pressurePowerToFluidWatts +=
            face.pressurePowerToFluidWatts;
        surface.pressureWorkToFluidJoules +=
            face.pressureWorkToFluidJoules;
        surface.pressureForceOnSurfaceNewtons = addVector(
            surface.pressureForceOnSurfaceNewtons,
            face.pressureForceOnSurfaceNewtons);
        surface.pressureImpulseOnSurfaceNewtonSeconds = addVector(
            surface.pressureImpulseOnSurfaceNewtonSeconds,
            face.pressureImpulseOnSurfaceNewtonSeconds);
        surface.pressurePowerToSurfaceWatts +=
            face.pressurePowerToSurfaceWatts;
        surface.pressureWorkToSurfaceJoules +=
            face.pressureWorkToSurfaceJoules;
        surface.dissipationWatts += face.dissipationWatts;
        surface.dissipatedEnergyJoules += face.dissipatedEnergyJoules;
        surface.energyResidualJoules += face.energyResidualJoules;

        result.totalPressureForceOnFluidNewtons = addVector(
            result.totalPressureForceOnFluidNewtons,
            face.pressureForceOnFluidNewtons);
        result.totalPressureImpulseOnFluidNewtonSeconds = addVector(
            result.totalPressureImpulseOnFluidNewtonSeconds,
            face.pressureImpulseOnFluidNewtonSeconds);
        result.totalPressurePowerToFluidWatts +=
            face.pressurePowerToFluidWatts;
        result.totalPressureWorkToFluidJoules +=
            face.pressureWorkToFluidJoules;
        result.totalPressureForceOnSurfaceNewtons = addVector(
            result.totalPressureForceOnSurfaceNewtons,
            face.pressureForceOnSurfaceNewtons);
        result.totalPressureImpulseOnSurfaceNewtonSeconds = addVector(
            result.totalPressureImpulseOnSurfaceNewtonSeconds,
            face.pressureImpulseOnSurfaceNewtonSeconds);
        result.totalPressurePowerToSurfaceWatts +=
            face.pressurePowerToSurfaceWatts;
        result.totalPressureWorkToSurfaceJoules +=
            face.pressureWorkToSurfaceJoules;
        result.totalDissipationWatts += face.dissipationWatts;
        result.totalDissipatedEnergyJoules +=
            face.dissipatedEnergyJoules;
        result.maximumAbsoluteFaceEnergyResidualJoules = std::max(
            result.maximumAbsoluteFaceEnergyResidualJoules,
            std::abs(face.energyResidualJoules));
    }

    result.surfaces.reserve(surfaces.size());
    bool surfacesAccepted = true;
    for (auto& [stableId, surface] : surfaces) {
        static_cast<void>(stableId);
        const double surfaceEnergyScale =
            std::abs(surface.pressureWorkToFluidJoules)
            + std::abs(surface.pressureWorkToSurfaceJoules)
            + std::abs(surface.dissipatedEnergyJoules);
        surfacesAccepted = surfacesAccepted
            && surface.faceCount > 0
            && std::isfinite(surface.areaSquareMeters)
            && surface.areaSquareMeters > 0.0
            && finiteVector(surface.pressureForceOnFluidNewtons)
            && finiteVector(
                surface.pressureImpulseOnFluidNewtonSeconds)
            && finiteVector(surface.pressureForceOnSurfaceNewtons)
            && finiteVector(
                surface.pressureImpulseOnSurfaceNewtonSeconds)
            && std::isfinite(surface.pressurePowerToFluidWatts)
            && std::isfinite(surface.pressureWorkToFluidJoules)
            && std::isfinite(surface.pressurePowerToSurfaceWatts)
            && std::isfinite(surface.pressureWorkToSurfaceJoules)
            && std::isfinite(surface.dissipationWatts)
            && std::isfinite(surface.dissipatedEnergyJoules)
            && std::isfinite(surface.energyResidualJoules)
            && std::abs(surface.energyResidualJoules)
                <= 1.0e-10 + 1.0e-11 * surfaceEnergyScale
            && closeVector(
                surface.pressureForceOnSurfaceNewtons,
                scaleVector(
                    surface.pressureForceOnFluidNewtons, -1.0))
            && closeVector(
                surface.pressureImpulseOnSurfaceNewtonSeconds,
                scaleVector(
                    surface.pressureImpulseOnFluidNewtonSeconds,
                    -1.0));
        result.surfaces.push_back(surface);
    }
    result.energyResidualJoules =
        result.totalPressureWorkToFluidJoules
        + result.totalPressureWorkToSurfaceJoules
        + result.totalDissipatedEnergyJoules;
    result.finite =
        finiteVector(result.totalPressureForceOnFluidNewtons)
        && finiteVector(
            result.totalPressureImpulseOnFluidNewtonSeconds)
        && finiteVector(result.totalPressureForceOnSurfaceNewtons)
        && finiteVector(
            result.totalPressureImpulseOnSurfaceNewtonSeconds)
        && std::isfinite(result.totalPressurePowerToFluidWatts)
        && std::isfinite(result.totalPressureWorkToFluidJoules)
        && std::isfinite(result.totalPressurePowerToSurfaceWatts)
        && std::isfinite(result.totalPressureWorkToSurfaceJoules)
        && std::isfinite(result.totalDissipationWatts)
        && std::isfinite(result.totalDissipatedEnergyJoules)
        && std::isfinite(result.energyResidualJoules)
        && std::isfinite(
            result.maximumAbsoluteFaceEnergyResidualJoules);
    const double energyScale =
        std::abs(result.totalPressureWorkToFluidJoules)
        + std::abs(result.totalPressureWorkToSurfaceJoules)
        + std::abs(result.totalDissipatedEnergyJoules);
    result.accepted = result.finite && surfacesAccepted
        && closeValue(
            result.totalDissipationWatts,
            source.totalDissipationWatts)
        && closeValue(
            result.totalDissipatedEnergyJoules,
            source.totalPorousDissipationJoules)
        && std::abs(result.energyResidualJoules)
            <= 1.0e-10 + 1.0e-11 * energyScale
        && result.maximumAbsoluteFaceEnergyResidualJoules
            <= 1.0e-10 + 1.0e-11 * energyScale
        && closeVector(
            result.totalPressureForceOnSurfaceNewtons,
            scaleVector(result.totalPressureForceOnFluidNewtons, -1.0))
        && closeVector(
            result.totalPressureImpulseOnSurfaceNewtonSeconds,
            scaleVector(
                result.totalPressureImpulseOnFluidNewtonSeconds, -1.0));
    return result;
}

} // namespace

PorousSurfaceTractionDiagnostics evaluatePorousSurfaceTraction(
    const PeriodicCartesianGrid& grid,
    const PorousProjectionDiagnostics& diagnostics,
    const double timeStepSeconds) {
    return evaluatePorousSurfaceTractionImpl(
        grid, diagnostics, timeStepSeconds);
}

PorousSurfaceTractionDiagnostics evaluatePorousSurfaceTraction(
    const PeriodicCartesianGrid& grid,
    const MovingPorousProjectionDiagnostics& diagnostics,
    const double timeStepSeconds) {
    if (!diagnostics.accepted || !diagnostics.finite
        || !diagnostics.porous.accepted || !diagnostics.porous.finite
        || !diagnostics.movingInterface.finite
        || diagnostics.porous.projection
            != diagnostics.movingInterface.projection) {
        throw std::invalid_argument(
            "porous traction requires accepted moving-porous diagnostics");
    }
    return evaluatePorousSurfaceTractionImpl(
        grid, diagnostics.porous, timeStepSeconds);
}

} // namespace simwing::fsi::fluid

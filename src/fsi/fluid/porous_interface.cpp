#include "fluid/porous_interface.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

void validateProjectionSettings(const PorousProjectionSettings& settings) {
    switch (settings.constitutiveEvaluation) {
    case PorousConstitutiveEvaluation::Endpoint:
    case PorousConstitutiveEvaluation::Midpoint:
        break;
    default:
        throw std::invalid_argument(
            "porous projection constitutive time is invalid");
    }
    if (!std::isfinite(
            settings.absoluteNormalVelocityToleranceMetersPerSecond)
        || settings.absoluteNormalVelocityToleranceMetersPerSecond < 0.0
        || !std::isfinite(settings.relativeNormalVelocityTolerance)
        || settings.relativeNormalVelocityTolerance < 0.0
        || (settings.absoluteNormalVelocityToleranceMetersPerSecond == 0.0
            && settings.relativeNormalVelocityTolerance == 0.0)) {
        throw std::invalid_argument(
            "porous projection velocity tolerances are invalid");
    }
    if (!std::isfinite(settings.absolutePressureJumpTolerancePascals)
        || settings.absolutePressureJumpTolerancePascals < 0.0
        || !std::isfinite(settings.relativePressureJumpTolerance)
        || settings.relativePressureJumpTolerance < 0.0
        || (settings.absolutePressureJumpTolerancePascals == 0.0
            && settings.relativePressureJumpTolerance == 0.0)) {
        throw std::invalid_argument(
            "porous projection pressure-jump tolerances are invalid");
    }
    if (!std::isfinite(settings.relaxation)
        || !(settings.relaxation > 0.0)
        || settings.relaxation > 1.0
        || settings.maximumNonlinearIterations == 0) {
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

PorousProjectionDiagnostics projectVelocityWithPorousInterfacesImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField* prescribedPressureJumps,
    const PorousProjectionSettings& settings) {
    validateProjectionSettings(settings);
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
    diagnostics.constitutiveEvaluation = settings.constitutiveEvaluation;
    diagnostics.porousCrossingCount = porousCrossings.size();
    if (porousCrossings.empty()) {
        if (prescribedPressureJumps == nullptr
            || prescribedPressureJumps->empty()) {
            diagnostics.projection = projectVelocity(
                grid, predictedVelocityMetersPerSecond, pressurePascals,
                settings.projection);
        } else {
            diagnostics.projection = projectVelocityWithPressureJumps(
                grid, predictedVelocityMetersPerSecond, pressurePascals,
                *prescribedPressureJumps, settings.projection);
        }
        diagnostics.accepted = diagnostics.projection.converged;
        diagnostics.finite = isFinite(predictedVelocityMetersPerSecond)
            && isFinite(pressurePascals);
        return diagnostics;
    }

    const MacVelocityField originalVelocity =
        predictedVelocityMetersPerSecond;
    const CellScalarField originalPressure = pressurePascals;
    MacVelocityField iterateVelocity = originalVelocity;
    CellScalarField pressureWarmStart = originalPressure;

    for (std::size_t iteration = 0;
         iteration < settings.maximumNonlinearIterations; ++iteration) {
        MacVelocityField sampledVelocity = iterateVelocity;
        moveToConstitutiveTime(
            originalVelocity, settings.constitutiveEvaluation,
            sampledVelocity);
        const PorousPressureJumpField sampled(
            grid, sampledVelocity, porousCrossings);
        const SharpPressureJumpField pressureJumps = combinedPressureJumps(
            grid, sampled, prescribedPressureJumps);
        MacVelocityField candidateVelocity = originalVelocity;
        CellScalarField candidatePressure = pressureWarmStart;
        diagnostics.projection = projectVelocityWithPressureJumps(
            grid, candidateVelocity, candidatePressure,
            pressureJumps, settings.projection);
        diagnostics.nonlinearIterationCount = iteration + 1;
        if (!diagnostics.projection.converged) {
            diagnostics.samples.assign(
                sampled.samples().begin(), sampled.samples().end());
            diagnostics.totalDissipationWatts =
                sampled.totalDissipationWatts();
            return diagnostics;
        }

        MacVelocityField candidateConstitutiveVelocity = candidateVelocity;
        moveToConstitutiveTime(
            originalVelocity, settings.constitutiveEvaluation,
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
        diagnostics.totalDissipationWatts =
            endpoint.totalDissipationWatts();
        diagnostics.samples.assign(
            endpointFaces.begin(), endpointFaces.end());
        diagnostics.finite = std::isfinite(maximumVelocityResidual)
            && std::isfinite(maximumJumpResidual)
            && std::isfinite(diagnostics.totalDissipationWatts)
            && isFinite(candidateVelocity) && isFinite(candidatePressure);
        const double velocityTolerance =
            settings.absoluteNormalVelocityToleranceMetersPerSecond
            + settings.relativeNormalVelocityTolerance * velocityScale;
        const double jumpTolerance =
            settings.absolutePressureJumpTolerancePascals
            + settings.relativePressureJumpTolerance * jumpScale;
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
            candidateVelocity, settings.relaxation, iterateVelocity);
        pressureWarmStart = std::move(candidatePressure);
    }
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

} // namespace simwing::fsi::fluid

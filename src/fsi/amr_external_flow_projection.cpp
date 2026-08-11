#include "amr_external_flow_projection.h"

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace simwing::fsi::amr {
namespace {

using FaceFields = std::array<std::unique_ptr<amrex::MultiFab>, 3>;

struct Hierarchy {
    amrex::Vector<amrex::Geometry> geometry;
    amrex::Vector<amrex::BoxArray> boxes;
    amrex::Vector<amrex::DistributionMapping> distribution;
    amrex::Vector<FaceFields> velocity;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> divergence;
};

void validateProjectionSettings(
    const WindTunnelProjectionSettings& settings) {
    if (!std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || !std::isfinite(settings.airDensityKilogramsPerCubicMeter)
        || !(settings.airDensityKilogramsPerCubicMeter > 0.0)
        || !std::isfinite(settings.relativeTolerance)
        || !(settings.relativeTolerance > 0.0)
        || !(settings.relativeTolerance < 1.0)
        || settings.maximumIterations == 0
        || settings.maximumIterations
            > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(
            "AMR wind-tunnel projection settings are invalid");
    }
}

amrex::IntVect cellCounts(const fluid::GridCellCounts counts) {
    if (counts.x > static_cast<std::size_t>(
                       std::numeric_limits<int>::max())
        || counts.y > static_cast<std::size_t>(
                          std::numeric_limits<int>::max())
        || counts.z > static_cast<std::size_t>(
                          std::numeric_limits<int>::max())) {
        throw std::overflow_error("AMR grid count exceeds integer indexing");
    }
    return {
        static_cast<int>(counts.x),
        static_cast<int>(counts.y),
        static_cast<int>(counts.z),
    };
}

FaceFields makeFaceFields(
    const amrex::BoxArray& cells,
    const amrex::DistributionMapping& distribution,
    const fluid::Vector3 freestream) {
    const std::array<double, 3> values{
        freestream.x, freestream.y, freestream.z};
    FaceFields result;
    for (int axis = 0; axis < 3; ++axis) {
        const amrex::BoxArray faces = amrex::convert(
            cells, amrex::IntVect::TheDimensionVector(axis));
        result[static_cast<std::size_t>(axis)] =
            std::make_unique<amrex::MultiFab>(
                faces, distribution, 1, 1);
        result[static_cast<std::size_t>(axis)]->setVal(values[axis]);
    }
    return result;
}

void initializeDivergentWake(
    FaceFields& fields,
    const amrex::Box& levelDomain,
    const fluid::Vector3 freestream) {
    amrex::MultiFab& streamwise = *fields[1];
    const amrex::Box faceDomain = amrex::convert(
        levelDomain, amrex::IntVect::TheDimensionVector(1));
    const double amplitude = 0.08 * freestream.y;
    const double inverseX = 1.0
        / static_cast<double>(std::max(1, levelDomain.length(0)));
    const double inverseY = 1.0
        / static_cast<double>(std::max(1, faceDomain.length(1) - 1));
    const double inverseZ = 1.0
        / static_cast<double>(std::max(1, levelDomain.length(2)));
    const int xLower = levelDomain.smallEnd(0);
    const int yLower = faceDomain.smallEnd(1);
    const int zLower = levelDomain.smallEnd(2);
    for (amrex::MFIter iterator(streamwise); iterator.isValid(); ++iterator) {
        const amrex::Box box = iterator.validbox();
        auto velocity = streamwise.array(iterator);
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     const int i, const int j, const int k) noexcept {
                const double x =
                    (static_cast<double>(i - xLower) + 0.5) * inverseX;
                const double y =
                    static_cast<double>(j - yLower) * inverseY;
                const double z =
                    (static_cast<double>(k - zLower) + 0.5) * inverseZ;
                const double crossSection =
                    std::sin(3.14159265358979323846 * x)
                    * std::sin(3.14159265358979323846 * z);
                velocity(i, j, k) = freestream.y
                    - amplitude * y * crossSection;
            });
    }
    amrex::Gpu::streamSynchronize();
}

Hierarchy makeHierarchy(const WindTunnelGridSettings& settings) {
    const amrex::IntVect counts = cellCounts(settings.coarseCellCounts);
    const amrex::Box coarseDomain(
        amrex::IntVect(0), counts - amrex::IntVect(1));
    const amrex::RealBox physicalDomain(
        {settings.lowerMeters.x,
         settings.lowerMeters.y,
         settings.lowerMeters.z},
        {settings.upperMeters.x,
         settings.upperMeters.y,
         settings.upperMeters.z});
    const std::array<int, 3> nonPeriodic{0, 0, 0};

    Hierarchy result;
    result.geometry.emplace_back(
        coarseDomain, &physicalDomain, amrex::CoordSys::cartesian,
        nonPeriodic.data());
    result.boxes.emplace_back(coarseDomain);
    result.boxes[0].maxSize(static_cast<int>(settings.maximumGridSize));
    result.distribution.emplace_back(result.boxes[0]);

    const int ratio = static_cast<int>(settings.refinementRatio);
    const amrex::Box refinedDomain = amrex::refine(coarseDomain, ratio);
    result.geometry.emplace_back(
        refinedDomain, &physicalDomain, amrex::CoordSys::cartesian,
        nonPeriodic.data());
    const amrex::IntVect quarter(
        static_cast<int>(settings.coarseCellCounts.x / 4),
        static_cast<int>(settings.coarseCellCounts.y / 4),
        static_cast<int>(settings.coarseCellCounts.z / 4));
    const amrex::IntVect threeQuarter(
        static_cast<int>(3 * settings.coarseCellCounts.x / 4 - 1),
        static_cast<int>(3 * settings.coarseCellCounts.y / 4 - 1),
        static_cast<int>(3 * settings.coarseCellCounts.z / 4 - 1));
    result.boxes.emplace_back(
        amrex::refine(amrex::Box(quarter, threeQuarter), ratio));
    result.boxes[1].maxSize(
        static_cast<int>(settings.maximumGridSize
                         * settings.refinementRatio));
    result.distribution.emplace_back(result.boxes[1]);

    for (std::size_t level = 0; level < result.boxes.size(); ++level) {
        result.velocity.push_back(makeFaceFields(
            result.boxes[level], result.distribution[level],
            settings.freestreamMetersPerSecond));
        initializeDivergentWake(
            result.velocity[level], result.geometry[level].Domain(),
            settings.freestreamMetersPerSecond);
        result.divergence.push_back(std::make_unique<amrex::MultiFab>(
            result.boxes[level], result.distribution[level], 1, 0));
    }
    return result;
}

void computeDivergence(Hierarchy& hierarchy) {
    for (std::size_t level = 0; level < hierarchy.velocity.size(); ++level) {
        const amrex::Array<const amrex::MultiFab*, 3> velocity{
            hierarchy.velocity[level][0].get(),
            hierarchy.velocity[level][1].get(),
            hierarchy.velocity[level][2].get(),
        };
        amrex::computeDivergence(
            *hierarchy.divergence[level], velocity,
            hierarchy.geometry[level]);
    }
}

void synchronizeCoarseFaces(Hierarchy& hierarchy) {
    const amrex::Array<const amrex::MultiFab*, 3> fine{
        hierarchy.velocity[1][0].get(),
        hierarchy.velocity[1][1].get(),
        hierarchy.velocity[1][2].get(),
    };
    const amrex::Array<amrex::MultiFab*, 3> coarse{
        hierarchy.velocity[0][0].get(),
        hierarchy.velocity[0][1].get(),
        hierarchy.velocity[0][2].get(),
    };
    amrex::average_down_faces(
        fine, coarse, amrex::IntVect(2), hierarchy.geometry[0]);
}

double activeMaximumMagnitude(
    const Hierarchy& hierarchy,
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& fields,
    std::size_t& activeCellCount) {
    const amrex::IntVect ratio(2);
    const amrex::iMultiFab coarseMask = amrex::makeFineMask(
        *fields[0], hierarchy.boxes[1], ratio, 1, 0);
    double maximum = 0.0;
    activeCellCount = 0;
    for (std::size_t level = 0; level < fields.size(); ++level) {
        const amrex::MultiFab& field = *fields[level];
        for (amrex::MFIter iterator(field); iterator.isValid(); ++iterator) {
            const amrex::Box box = iterator.validbox();
            const auto values = field.const_array(iterator);
            const auto mask = level == 0
                ? coarseMask.const_array(iterator)
                : amrex::Array4<const int>();
            for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
                for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                    for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                        if (level == 0 && mask(i, j, k) == 0) {
                            continue;
                        }
                        maximum = std::max(
                            maximum, std::abs(values(i, j, k)));
                        ++activeCellCount;
                    }
                }
            }
        }
    }
    return maximum;
}

double maximumActivePressurePascals(
    const Hierarchy& hierarchy,
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& potential,
    const double density) {
    std::size_t ignoredCellCount = 0;
    return density * activeMaximumMagnitude(
        hierarchy, potential, ignoredCellCount);
}

void auditBoundaryNormalVelocity(
    const Hierarchy& hierarchy,
    const fluid::Vector3 freestream,
    WindTunnelProjectionDiagnostics& diagnostics) {
    const amrex::Box& domain = hierarchy.geometry[0].Domain();
    const amrex::MultiFab& streamwise = *hierarchy.velocity[0][1];
    const double amplitude = 0.08 * freestream.y;
    const double inverseX =
        1.0 / static_cast<double>(domain.length(0));
    const double inverseZ =
        1.0 / static_cast<double>(domain.length(2));
    for (amrex::MFIter iterator(streamwise); iterator.isValid(); ++iterator) {
        const amrex::Box box = iterator.validbox();
        const auto values = streamwise.const_array(iterator);
        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                if (box.contains(amrex::IntVect(i, domain.smallEnd(1), k))) {
                    diagnostics.lowerYInflowNormalVelocityErrorMetersPerSecond =
                        std::max(
                            diagnostics
                                .lowerYInflowNormalVelocityErrorMetersPerSecond,
                            std::abs(
                                values(i, domain.smallEnd(1), k)
                                - freestream.y));
                }
                const int outlet = domain.bigEnd(1) + 1;
                if (box.contains(amrex::IntVect(i, outlet, k))) {
                    const double x =
                        (static_cast<double>(i - domain.smallEnd(0)) + 0.5)
                        * inverseX;
                    const double z =
                        (static_cast<double>(k - domain.smallEnd(2)) + 0.5)
                        * inverseZ;
                    const double originalOutletVelocity = freestream.y
                        - amplitude
                            * std::sin(3.14159265358979323846 * x)
                            * std::sin(3.14159265358979323846 * z);
                    diagnostics.upperYOutletNormalVelocityChangeMetersPerSecond =
                        std::max(
                            diagnostics
                                .upperYOutletNormalVelocityChangeMetersPerSecond,
                            std::abs(
                                values(i, outlet, k)
                                - originalOutletVelocity));
                }
            }
        }
    }
}

} // namespace

WindTunnelProjectionDiagnostics evaluateWindTunnelPressureProjection(
    const WindTunnelProjectionSettings& settings) {
    validateProjectionSettings(settings);
    if (!amrex::Initialized()) {
        throw std::logic_error("AMReX runtime is not initialized");
    }

    WindTunnelProjectionDiagnostics diagnostics;
    diagnostics.hierarchy =
        evaluateWindTunnelBoundaryInitialization(settings.grid);
    Hierarchy hierarchy = makeHierarchy(settings.grid);
    synchronizeCoarseFaces(hierarchy);
    computeDivergence(hierarchy);
    std::size_t initialActiveCells = 0;
    diagnostics.initialMaximumDivergencePerSecond =
        activeMaximumMagnitude(
            hierarchy, hierarchy.divergence, initialActiveCells);

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> potential;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rightHandSide;
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3>> gradient;
    for (std::size_t level = 0; level < hierarchy.boxes.size(); ++level) {
        potential.push_back(std::make_unique<amrex::MultiFab>(
            hierarchy.boxes[level], hierarchy.distribution[level], 1, 1));
        rightHandSide.push_back(std::make_unique<amrex::MultiFab>(
            hierarchy.boxes[level], hierarchy.distribution[level], 1, 0));
        gradient.emplace_back();
        potential.back()->setVal(0.0);
        amrex::MultiFab::LinComb(
            *rightHandSide.back(),
            1.0 / settings.timeStepSeconds,
            *hierarchy.divergence[level], 0,
            0.0, *hierarchy.divergence[level], 0,
            0, 1, 0);
        for (int axis = 0; axis < 3; ++axis) {
            const amrex::BoxArray faces = amrex::convert(
                hierarchy.boxes[level],
                amrex::IntVect::TheDimensionVector(axis));
            gradient.back()[static_cast<std::size_t>(axis)] =
                std::make_unique<amrex::MultiFab>(
                    faces, hierarchy.distribution[level], 1, 0);
        }
    }

    amrex::LPInfo information;
    information.setAgglomeration(true);
    information.setConsolidation(true);
    information.setMaxCoarseningLevel(30);
    amrex::MLPoisson poisson(
        hierarchy.geometry, hierarchy.boxes,
        hierarchy.distribution, information);
    poisson.setMaxOrder(2);
    poisson.setDomainBC(
        {amrex::LinOpBCType::Neumann,
         amrex::LinOpBCType::Neumann,
         amrex::LinOpBCType::Neumann},
        {amrex::LinOpBCType::Neumann,
         amrex::LinOpBCType::Dirichlet,
         amrex::LinOpBCType::Neumann});
    for (int level = 0;
         level < static_cast<int>(hierarchy.boxes.size()); ++level) {
        poisson.setLevelBC(level, nullptr);
    }

    amrex::Vector<amrex::MultiFab*> potentialPointers;
    amrex::Vector<const amrex::MultiFab*> rightHandSidePointers;
    amrex::Vector<amrex::Array<amrex::MultiFab*, 3>> gradientPointers;
    for (std::size_t level = 0; level < hierarchy.boxes.size(); ++level) {
        potentialPointers.push_back(potential[level].get());
        rightHandSidePointers.push_back(rightHandSide[level].get());
        gradientPointers.push_back({
            gradient[level][0].get(),
            gradient[level][1].get(),
            gradient[level][2].get(),
        });
    }

    amrex::MLMG multigrid(poisson);
    multigrid.setVerbose(0);
    multigrid.setBottomVerbose(0);
    multigrid.setMaxIter(static_cast<int>(settings.maximumIterations));
    multigrid.setBottomSolver(amrex::MLMG::BottomSolver::bicgstab);
    multigrid.solve(
        potentialPointers, rightHandSidePointers,
        settings.relativeTolerance, 0.0);
    multigrid.getGradSolution(gradientPointers);
    diagnostics.solverIterations =
        static_cast<std::size_t>(multigrid.getNumIters());
    diagnostics.solverFinalResidual = multigrid.getFinalResidual();

    for (std::size_t level = 0; level < hierarchy.velocity.size(); ++level) {
        for (int axis = 0; axis < 3; ++axis) {
            amrex::MultiFab::Saxpy(
                *hierarchy.velocity[level][static_cast<std::size_t>(axis)],
                -settings.timeStepSeconds,
                *gradient[level][static_cast<std::size_t>(axis)],
                0, 0, 1, 0);
        }
    }
    synchronizeCoarseFaces(hierarchy);
    computeDivergence(hierarchy);
    diagnostics.projectedMaximumDivergencePerSecond =
        activeMaximumMagnitude(
            hierarchy, hierarchy.divergence,
            diagnostics.activeCompositeCellCount);
    diagnostics.maximumDivergenceReductionRatio =
        diagnostics.initialMaximumDivergencePerSecond > 0.0
        ? diagnostics.projectedMaximumDivergencePerSecond
            / diagnostics.initialMaximumDivergencePerSecond
        : 0.0;
    diagnostics.maximumPressureCorrectionPascals =
        maximumActivePressurePascals(
            hierarchy, potential,
            settings.airDensityKilogramsPerCubicMeter);
    auditBoundaryNormalVelocity(
        hierarchy, settings.grid.freestreamMetersPerSecond,
        diagnostics);
    diagnostics.pressureOutletReferenceOwned = true;
    diagnostics.finite = diagnostics.hierarchy.finite
        && std::isfinite(diagnostics.solverFinalResidual)
        && std::isfinite(diagnostics.initialMaximumDivergencePerSecond)
        && std::isfinite(diagnostics.projectedMaximumDivergencePerSecond)
        && std::isfinite(diagnostics.maximumDivergenceReductionRatio)
        && std::isfinite(diagnostics.maximumPressureCorrectionPascals)
        && std::isfinite(
            diagnostics.lowerYInflowNormalVelocityErrorMetersPerSecond)
        && std::isfinite(
            diagnostics.upperYOutletNormalVelocityChangeMetersPerSecond);
    diagnostics.accepted = diagnostics.hierarchy.accepted
        && diagnostics.finite
        && diagnostics.activeCompositeCellCount == initialActiveCells
        && diagnostics.solverIterations > 0
        && diagnostics.initialMaximumDivergencePerSecond > 0.0
        && diagnostics.maximumPressureCorrectionPascals > 0.0
        && diagnostics.maximumDivergenceReductionRatio < 1.0e-7
        && diagnostics.lowerYInflowNormalVelocityErrorMetersPerSecond
            < 1.0e-11
        && diagnostics.upperYOutletNormalVelocityChangeMetersPerSecond
            > 0.0
        && diagnostics.pressureOutletReferenceOwned;
    return diagnostics;
}

} // namespace simwing::fsi::amr

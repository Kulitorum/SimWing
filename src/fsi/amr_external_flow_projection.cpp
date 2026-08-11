#include "amr_external_flow_projection.h"

#include <AMReX_Array.H>
#include <AMReX_BCRec.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFIter.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

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
            > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || settings.staticWingForcingProjectionIterations == 0
        || settings.staticWingForcingProjectionIterations > 100
        || !std::isfinite(settings.staticWingDirectForcingRelaxation)
        || !(settings.staticWingDirectForcingRelaxation > 0.0)
        || settings.staticWingDirectForcingRelaxation > 2.0) {
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
    const amrex::IntVect refinementLower(
        static_cast<int>(settings.refinedCoarseCellLower[0]),
        static_cast<int>(settings.refinedCoarseCellLower[1]),
        static_cast<int>(settings.refinedCoarseCellLower[2]));
    const amrex::IntVect refinementUpper(
        static_cast<int>(settings.refinedCoarseCellUpperExclusive[0] - 1),
        static_cast<int>(settings.refinedCoarseCellUpperExclusive[1] - 1),
        static_cast<int>(settings.refinedCoarseCellUpperExclusive[2] - 1));
    result.boxes.emplace_back(
        amrex::refine(
            amrex::Box(refinementLower, refinementUpper), ratio));
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

namespace {

void copyProjectedCoarseGrid(
    const WindTunnelProjectionSettings& settings,
    const Hierarchy& hierarchy,
    const amrex::MultiFab& potential,
    WindTunnelProjectedCoarseGrid& output) {
    const auto counts = settings.grid.coarseCellCounts;
    const std::size_t cellCount = counts.x * counts.y * counts.z;
    output.cellCounts = counts;
    output.lowerMeters = settings.grid.lowerMeters;
    output.upperMeters = settings.grid.upperMeters;
    output.velocityMetersPerSecond.assign(cellCount, {});
    output.pressurePascals.assign(cellCount, 0.0);
    output.divergencePerSecond.assign(cellCount, 0.0);

    const amrex::MultiFab& divergence = *hierarchy.divergence[0];
    for (amrex::MFIter iterator(divergence); iterator.isValid(); ++iterator) {
        const amrex::Box box = iterator.validbox();
        const auto xVelocity = hierarchy.velocity[0][0]->const_array(iterator);
        const auto yVelocity = hierarchy.velocity[0][1]->const_array(iterator);
        const auto zVelocity = hierarchy.velocity[0][2]->const_array(iterator);
        const auto pressurePotential = potential.const_array(iterator);
        const auto divergenceValues = divergence.const_array(iterator);
        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    const std::size_t index = static_cast<std::size_t>(i)
                        + counts.x
                            * (static_cast<std::size_t>(j)
                               + counts.y * static_cast<std::size_t>(k));
                    output.velocityMetersPerSecond[index] = {
                        0.5 * (xVelocity(i, j, k)
                               + xVelocity(i + 1, j, k)),
                        0.5 * (yVelocity(i, j, k)
                               + yVelocity(i, j + 1, k)),
                        0.5 * (zVelocity(i, j, k)
                               + zVelocity(i, j, k + 1)),
                    };
                    output.pressurePascals[index] =
                        settings.airDensityKilogramsPerCubicMeter
                        * pressurePotential(i, j, k);
                    output.divergencePerSecond[index] =
                        divergenceValues(i, j, k);
                }
            }
        }
    }
}

void fillCoarsePhysicalGhosts(
    Hierarchy& hierarchy,
    const fluid::Vector3 freestream) {
    const std::array<double, 3> prescribed{
        freestream.x, freestream.y, freestream.z};
    const amrex::Box& cellDomain = hierarchy.geometry[0].Domain();
    for (int axis = 0; axis < 3; ++axis) {
        amrex::MultiFab& field = *hierarchy.velocity[0][axis];
        field.FillBoundary(hierarchy.geometry[0].periodicity());
        const amrex::Box faceDomain = amrex::convert(
            cellDomain, amrex::IntVect::TheDimensionVector(axis));
        for (amrex::MFIter iterator(field); iterator.isValid(); ++iterator) {
            const amrex::Box box = iterator.fabbox();
            auto values = field.array(iterator);
            for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
                for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                    for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                        const bool outsideX = i < faceDomain.smallEnd(0)
                            || i > faceDomain.bigEnd(0);
                        const bool outsideZ = k < faceDomain.smallEnd(2)
                            || k > faceDomain.bigEnd(2);
                        if (outsideX || outsideZ
                            || j < faceDomain.smallEnd(1)) {
                            values(i, j, k) = prescribed[axis];
                        } else if (j > faceDomain.bigEnd(1)) {
                            values(i, j, k) =
                                values(i, faceDomain.bigEnd(1), k);
                        }
                    }
                }
            }
        }
    }
}

void fillFineVelocityGhosts(Hierarchy& hierarchy) {
    const amrex::IntVect ratio(2);
    const amrex::Array<amrex::MultiFab*, 3> coarse{
        hierarchy.velocity[0][0].get(),
        hierarchy.velocity[0][1].get(),
        hierarchy.velocity[0][2].get(),
    };
    const amrex::Array<amrex::MultiFab*, 3> fine{
        hierarchy.velocity[1][0].get(),
        hierarchy.velocity[1][1].get(),
        hierarchy.velocity[1][2].get(),
    };
    amrex::Vector<amrex::Array<amrex::MultiFab*, 3>> coarseData{coarse};
    amrex::Vector<amrex::Array<amrex::MultiFab*, 3>> fineData{fine};
    const amrex::Vector<amrex::Real> times{0.0};
    amrex::Array<amrex::Vector<amrex::BCRec>, 3> boundaryRecords;
    for (auto& componentRecords : boundaryRecords) {
        componentRecords.resize(1);
        for (int axis = 0; axis < 3; ++axis) {
            componentRecords[0].setLo(axis, amrex::BCType::int_dir);
            componentRecords[0].setHi(axis, amrex::BCType::int_dir);
        }
    }
    amrex::Array<amrex::PhysBCFunctNoOp, 3> physicalBoundaries;
    amrex::Interpolater* interpolater = &amrex::face_cons_linear_interp;
    amrex::FillPatchTwoLevels(
        fine, 0.0, coarseData, times, fineData, times,
        0, 0, 1,
        hierarchy.geometry[0], hierarchy.geometry[1],
        physicalBoundaries, 0, physicalBoundaries, 0,
        ratio, interpolater, boundaryRecords, 0);
}

void fillVelocityGhosts(
    Hierarchy& hierarchy,
    const fluid::Vector3 freestream) {
    fillCoarsePhysicalGhosts(hierarchy, freestream);
    fillFineVelocityGhosts(hierarchy);
}

void enforceCoarseNormalVelocityBoundaries(
    Hierarchy& hierarchy,
    const fluid::Vector3 freestream) {
    const amrex::Box& domain = hierarchy.geometry[0].Domain();
    for (int axis = 0; axis < 3; ++axis) {
        amrex::MultiFab& field = *hierarchy.velocity[0][axis];
        for (amrex::MFIter iterator(field); iterator.isValid(); ++iterator) {
            const amrex::Box box = iterator.validbox();
            auto values = field.array(iterator);
            for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
                for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                    for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                        if (axis == 0
                            && (i == domain.smallEnd(0)
                                || i == domain.bigEnd(0) + 1)) {
                            values(i, j, k) = freestream.x;
                        } else if (axis == 1
                                   && j == domain.smallEnd(1)) {
                            values(i, j, k) = freestream.y;
                        } else if (axis == 1
                                   && j == domain.bigEnd(1) + 1) {
                            values(i, j, k) = values(i, j - 1, k);
                        } else if (axis == 2
                                   && (k == domain.smallEnd(2)
                                       || k == domain.bigEnd(2) + 1)) {
                            values(i, j, k) = freestream.z;
                        }
                    }
                }
            }
        }
    }
}

double advectMomentumPredictor(
    Hierarchy& hierarchy,
    const WindTunnelProjectionSettings& settings) {
    fillVelocityGhosts(
        hierarchy, settings.grid.freestreamMetersPerSecond);
    double maximumCourant = 0.0;
    for (std::size_t level = 0;
         level < hierarchy.velocity.size(); ++level) {
        FaceFields next;
        for (int component = 0; component < 3; ++component) {
            const amrex::BoxArray faces = amrex::convert(
                hierarchy.boxes[level],
                amrex::IntVect::TheDimensionVector(component));
            next[component] = std::make_unique<amrex::MultiFab>(
                faces, hierarchy.distribution[level], 1, 1);
            next[component]->setVal(0.0);
        }
        const auto inverseSpacing =
            hierarchy.geometry[level].InvCellSizeArray();
        for (int component = 0; component < 3; ++component) {
            const amrex::MultiFab& source =
                *hierarchy.velocity[level][component];
            for (amrex::MFIter iterator(source);
                 iterator.isValid(); ++iterator) {
                const amrex::Box box = iterator.validbox();
                const auto u =
                    hierarchy.velocity[level][0]->const_array(iterator);
                const auto v =
                    hierarchy.velocity[level][1]->const_array(iterator);
                const auto w =
                    hierarchy.velocity[level][2]->const_array(iterator);
                const auto q = source.const_array(iterator);
                auto destination = next[component]->array(iterator);
                for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
                    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                        for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                            double advectX = 0.0;
                            double advectY = 0.0;
                            double advectZ = 0.0;
                            if (component == 0) {
                                advectX = u(i, j, k);
                                advectY = 0.25
                                    * (v(i - 1, j, k) + v(i, j, k)
                                       + v(i - 1, j + 1, k)
                                       + v(i, j + 1, k));
                                advectZ = 0.25
                                    * (w(i - 1, j, k) + w(i, j, k)
                                       + w(i - 1, j, k + 1)
                                       + w(i, j, k + 1));
                            } else if (component == 1) {
                                advectX = 0.25
                                    * (u(i, j - 1, k) + u(i + 1, j - 1, k)
                                       + u(i, j, k) + u(i + 1, j, k));
                                advectY = v(i, j, k);
                                advectZ = 0.25
                                    * (w(i, j - 1, k) + w(i, j, k)
                                       + w(i, j - 1, k + 1)
                                       + w(i, j, k + 1));
                            } else {
                                advectX = 0.25
                                    * (u(i, j, k - 1) + u(i + 1, j, k - 1)
                                       + u(i, j, k) + u(i + 1, j, k));
                                advectY = 0.25
                                    * (v(i, j, k - 1) + v(i, j + 1, k - 1)
                                       + v(i, j, k) + v(i, j + 1, k));
                                advectZ = w(i, j, k);
                            }
                            const double center = q(i, j, k);
                            const double derivativeX = inverseSpacing[0]
                                * (advectX >= 0.0
                                       ? center - q(i - 1, j, k)
                                       : q(i + 1, j, k) - center);
                            const double derivativeY = inverseSpacing[1]
                                * (advectY >= 0.0
                                       ? center - q(i, j - 1, k)
                                       : q(i, j + 1, k) - center);
                            const double derivativeZ = inverseSpacing[2]
                                * (advectZ >= 0.0
                                       ? center - q(i, j, k - 1)
                                       : q(i, j, k + 1) - center);
                            const double predicted = center
                                - settings.timeStepSeconds
                                    * (advectX * derivativeX
                                       + advectY * derivativeY
                                       + advectZ * derivativeZ);
                            if (!std::isfinite(center)
                                || !std::isfinite(advectX)
                                || !std::isfinite(advectY)
                                || !std::isfinite(advectZ)
                                || !std::isfinite(derivativeX)
                                || !std::isfinite(derivativeY)
                                || !std::isfinite(derivativeZ)
                                || !std::isfinite(predicted)) {
                                throw std::runtime_error(
                                    "AMR momentum predictor encountered a non-finite stencil at level "
                                    + std::to_string(level)
                                    + " component "
                                    + std::to_string(component)
                                    + " face ("
                                    + std::to_string(i) + ","
                                    + std::to_string(j) + ","
                                    + std::to_string(k) + ") values center="
                                    + std::to_string(center)
                                    + " advect=("
                                    + std::to_string(advectX) + ","
                                    + std::to_string(advectY) + ","
                                    + std::to_string(advectZ)
                                    + ") derivatives=("
                                    + std::to_string(derivativeX) + ","
                                    + std::to_string(derivativeY) + ","
                                    + std::to_string(derivativeZ) + ")");
                            }
                            destination(i, j, k) = predicted;
                            maximumCourant = std::max(
                                maximumCourant,
                                settings.timeStepSeconds
                                    * (std::abs(advectX) * inverseSpacing[0]
                                       + std::abs(advectY) * inverseSpacing[1]
                                       + std::abs(advectZ)
                                           * inverseSpacing[2]));
                        }
                    }
                }
            }
        }
        hierarchy.velocity[level] = std::move(next);
    }
    if (!std::isfinite(maximumCourant) || !(maximumCourant <= 1.0)) {
        throw std::runtime_error(
            "AMR momentum predictor exceeds the donor-cell CFL bound");
    }
    enforceCoarseNormalVelocityBoundaries(
        hierarchy, settings.grid.freestreamMetersPerSecond);
    synchronizeCoarseFaces(hierarchy);
    return maximumCourant;
}

std::size_t interfaceCellIndex(
    const fluid::GridCellCounts counts,
    const int i,
    const int j,
    const int k) {
    return static_cast<std::size_t>(i)
        + counts.x
            * (static_cast<std::size_t>(j)
               + counts.y * static_cast<std::size_t>(k));
}

fluid::Vector3 cellVelocity(
    const amrex::Array4<const double>& xVelocity,
    const amrex::Array4<const double>& yVelocity,
    const amrex::Array4<const double>& zVelocity,
    const int i,
    const int j,
    const int k) {
    return {
        0.5 * (xVelocity(i, j, k) + xVelocity(i + 1, j, k)),
        0.5 * (yVelocity(i, j, k) + yVelocity(i, j + 1, k)),
        0.5 * (zVelocity(i, j, k) + zVelocity(i, j, k + 1)),
    };
}

fluid::Vector3 projectedNormalVelocity(
    const StaticWingInterfaceCell& cell,
    const fluid::Vector3 velocity) {
    const auto& tensor = cell.normalProjectionTensor;
    return {
        tensor[0] * velocity.x + tensor[1] * velocity.y
            + tensor[2] * velocity.z,
        tensor[1] * velocity.x + tensor[3] * velocity.y
            + tensor[4] * velocity.z,
        tensor[2] * velocity.x + tensor[4] * velocity.y
            + tensor[5] * velocity.z,
    };
}

double normalSpeed(const StaticWingInterfaceCell& cell,
                   const fluid::Vector3 velocity) {
    const fluid::Vector3 normal = projectedNormalVelocity(cell, velocity);
    const double squared = velocity.x * normal.x
        + velocity.y * normal.y + velocity.z * normal.z;
    return std::sqrt(std::max(0.0, squared));
}

double maximumStaticWingNormalSpeed(
    const Hierarchy& hierarchy,
    const StaticWingInterface& interface) {
    const amrex::IntVect ratio(2);
    const amrex::iMultiFab coarseMask = amrex::makeFineMask(
        *hierarchy.divergence[0], hierarchy.boxes[1], ratio, 1, 0);
    double maximum = 0.0;
    for (std::size_t level = 0;
         level < hierarchy.velocity.size(); ++level) {
        const auto& interfaceLevel = interface.level(level);
        const auto counts = interfaceLevel.cellCounts;
        const amrex::MultiFab& cells = *hierarchy.divergence[level];
        for (amrex::MFIter iterator(cells); iterator.isValid(); ++iterator) {
            const amrex::Box box = iterator.validbox();
            const auto xVelocity =
                hierarchy.velocity[level][0]->const_array(iterator);
            const auto yVelocity =
                hierarchy.velocity[level][1]->const_array(iterator);
            const auto zVelocity =
                hierarchy.velocity[level][2]->const_array(iterator);
            const auto mask = level == 0
                ? coarseMask.const_array(iterator)
                : amrex::Array4<const int>();
            for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
                for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                    for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                        if (level == 0 && mask(i, j, k) == 0) {
                            continue;
                        }
                        const auto& cell = interfaceLevel.cells.at(
                            interfaceCellIndex(counts, i, j, k));
                        if (!cell.active()) {
                            continue;
                        }
                        maximum = std::max(
                            maximum,
                            normalSpeed(
                                cell,
                                cellVelocity(
                                    xVelocity, yVelocity, zVelocity,
                                    i, j, k)));
                    }
                }
            }
        }
    }
    return maximum;
}

WindTunnelMomentumStepDiagnostics::StaticWingDirectForcingDiagnostics
applyStaticWingDirectForcing(
    Hierarchy& hierarchy,
    const WindTunnelProjectionSettings& settings,
    const StaticWingInterface& interface) {
    WindTunnelMomentumStepDiagnostics::StaticWingDirectForcingDiagnostics
        diagnostics;
    diagnostics.binding = interface.diagnostics();
    diagnostics.active = true;
    diagnostics.maximumSurfaceNormalSpeedBeforeMetersPerSecond =
        maximumStaticWingNormalSpeed(hierarchy, interface);

    for (std::size_t level = 0;
         level < hierarchy.velocity.size(); ++level) {
        const auto& interfaceLevel = interface.level(level);
        const auto counts = interfaceLevel.cellCounts;
        amrex::MultiFab correction(
            hierarchy.boxes[level], hierarchy.distribution[level], 3, 1);
        amrex::MultiFab active(
            hierarchy.boxes[level], hierarchy.distribution[level], 1, 1);
        correction.setVal(0.0);
        active.setVal(0.0);
        for (amrex::MFIter iterator(correction);
             iterator.isValid(); ++iterator) {
            const amrex::Box box = iterator.validbox();
            const auto xVelocity =
                hierarchy.velocity[level][0]->const_array(iterator);
            const auto yVelocity =
                hierarchy.velocity[level][1]->const_array(iterator);
            const auto zVelocity =
                hierarchy.velocity[level][2]->const_array(iterator);
            auto delta = correction.array(iterator);
            auto activeValues = active.array(iterator);
            for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
                for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                    for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                        const auto& cell = interfaceLevel.cells.at(
                            interfaceCellIndex(counts, i, j, k));
                        if (!cell.active()) {
                            continue;
                        }
                        const fluid::Vector3 projected =
                            projectedNormalVelocity(
                                cell,
                                cellVelocity(
                                    xVelocity, yVelocity, zVelocity,
                                    i, j, k));
                        delta(i, j, k, 0) =
                            -settings.staticWingDirectForcingRelaxation
                            * projected.x;
                        delta(i, j, k, 1) =
                            -settings.staticWingDirectForcingRelaxation
                            * projected.y;
                        delta(i, j, k, 2) =
                            -settings.staticWingDirectForcingRelaxation
                            * projected.z;
                        activeValues(i, j, k) = 1.0;
                    }
                }
            }
        }
        correction.FillBoundary(hierarchy.geometry[level].periodicity());
        active.FillBoundary(hierarchy.geometry[level].periodicity());
        const amrex::Box& domain = hierarchy.geometry[level].Domain();
        for (int component = 0; component < 3; ++component) {
            amrex::MultiFab& velocity =
                *hierarchy.velocity[level][component];
            for (amrex::MFIter iterator(velocity);
                 iterator.isValid(); ++iterator) {
                const amrex::Box box = iterator.validbox();
                auto values = velocity.array(iterator);
                const auto delta = correction.const_array(iterator);
                const auto activeValues = active.const_array(iterator);
                for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
                    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                        for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                            amrex::IntVect upper(i, j, k);
                            amrex::IntVect lower = upper;
                            --lower[component];
                            double summedDelta = 0.0;
                            double weight = 0.0;
                            if (domain.contains(lower)) {
                                const double cellWeight = activeValues(
                                    lower[0], lower[1], lower[2]);
                                summedDelta += cellWeight * delta(
                                    lower[0], lower[1], lower[2], component);
                                weight += cellWeight;
                            }
                            if (domain.contains(upper)) {
                                const double cellWeight = activeValues(
                                    upper[0], upper[1], upper[2]);
                                summedDelta += cellWeight * delta(
                                    upper[0], upper[1], upper[2], component);
                                weight += cellWeight;
                            }
                            if (weight > 0.0) {
                                values(i, j, k) += summedDelta / weight;
                            }
                        }
                    }
                }
            }
        }
    }
    enforceCoarseNormalVelocityBoundaries(
        hierarchy, settings.grid.freestreamMetersPerSecond);
    synchronizeCoarseFaces(hierarchy);
    diagnostics.maximumSurfaceNormalSpeedAfterForcingMetersPerSecond =
        maximumStaticWingNormalSpeed(hierarchy, interface);
    diagnostics.finite = diagnostics.binding.finite
        && std::isfinite(
            diagnostics.maximumSurfaceNormalSpeedBeforeMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumSurfaceNormalSpeedAfterForcingMetersPerSecond);
    diagnostics.accepted = diagnostics.binding.accepted
        && diagnostics.finite
        && diagnostics.maximumSurfaceNormalSpeedBeforeMetersPerSecond > 0.0
        && diagnostics.maximumSurfaceNormalSpeedAfterForcingMetersPerSecond
            < diagnostics.maximumSurfaceNormalSpeedBeforeMetersPerSecond;
    return diagnostics;
}

double coarseKineticEnergy(
    const WindTunnelProjectedCoarseGrid& grid,
    const double density) {
    const fluid::Vector3 spacing = {
        (grid.upperMeters.x - grid.lowerMeters.x)
            / static_cast<double>(grid.cellCounts.x),
        (grid.upperMeters.y - grid.lowerMeters.y)
            / static_cast<double>(grid.cellCounts.y),
        (grid.upperMeters.z - grid.lowerMeters.z)
            / static_cast<double>(grid.cellCounts.z),
    };
    const double cellMass = density * spacing.x * spacing.y * spacing.z;
    double result = 0.0;
    for (const fluid::Vector3 velocity : grid.velocityMetersPerSecond) {
        result += 0.5 * cellMass
            * (velocity.x * velocity.x
               + velocity.y * velocity.y
               + velocity.z * velocity.z);
    }
    return result;
}

double maximumCellVelocityChange(
    const WindTunnelProjectedCoarseGrid& before,
    const WindTunnelProjectedCoarseGrid& after) {
    if (before.velocityMetersPerSecond.size()
        != after.velocityMetersPerSecond.size()) {
        throw std::logic_error(
            "AMR momentum snapshots have inconsistent velocity storage");
    }
    double result = 0.0;
    for (std::size_t index = 0;
         index < before.velocityMetersPerSecond.size(); ++index) {
        const fluid::Vector3 difference{
            after.velocityMetersPerSecond[index].x
                - before.velocityMetersPerSecond[index].x,
            after.velocityMetersPerSecond[index].y
                - before.velocityMetersPerSecond[index].y,
            after.velocityMetersPerSecond[index].z
                - before.velocityMetersPerSecond[index].z,
        };
        result = std::max(
            result,
            std::hypot(difference.x, difference.y, difference.z));
    }
    return result;
}

WindTunnelProjectionDiagnostics projectHierarchy(
    const WindTunnelProjectionSettings& settings,
    Hierarchy& hierarchy,
    const WindTunnelBoundaryDiagnostics& hierarchyDiagnostics,
    WindTunnelProjectedCoarseGrid* projectedCoarseGrid,
    const bool requireNontrivialProjection) {
    WindTunnelProjectionDiagnostics diagnostics;
    diagnostics.hierarchy = hierarchyDiagnostics;
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
        && (!requireNontrivialProjection
            || (diagnostics.initialMaximumDivergencePerSecond > 0.0
                && diagnostics.maximumPressureCorrectionPascals > 0.0))
        && (diagnostics.initialMaximumDivergencePerSecond == 0.0
            ? diagnostics.projectedMaximumDivergencePerSecond == 0.0
            : diagnostics.maximumDivergenceReductionRatio < 1.0e-7)
        && diagnostics.lowerYInflowNormalVelocityErrorMetersPerSecond
            < 1.0e-11
        && diagnostics.upperYOutletNormalVelocityChangeMetersPerSecond
            > 0.0
        && diagnostics.pressureOutletReferenceOwned;
    if (projectedCoarseGrid != nullptr) {
        amrex::average_down(
            *potential[1], *potential[0], 0, 1, amrex::IntVect(2));
        copyProjectedCoarseGrid(
            settings, hierarchy, *potential[0], *projectedCoarseGrid);
        projectedCoarseGrid->diagnostics = diagnostics;
    }
    return diagnostics;
}

WindTunnelProjectionDiagnostics evaluateProjection(
    const WindTunnelProjectionSettings& settings,
    WindTunnelProjectedCoarseGrid* projectedCoarseGrid) {
    validateProjectionSettings(settings);
    if (!amrex::Initialized()) {
        throw std::logic_error("AMReX runtime is not initialized");
    }
    const WindTunnelBoundaryDiagnostics hierarchyDiagnostics =
        evaluateWindTunnelBoundaryInitialization(settings.grid);
    Hierarchy hierarchy = makeHierarchy(settings.grid);
    return projectHierarchy(
        settings, hierarchy, hierarchyDiagnostics, projectedCoarseGrid, true);
}

} // namespace

WindTunnelProjectionDiagnostics evaluateWindTunnelPressureProjection(
    const WindTunnelProjectionSettings& settings) {
    return evaluateProjection(settings, nullptr);
}

WindTunnelProjectedCoarseGrid evaluateWindTunnelProjectedCoarseGrid(
    const WindTunnelProjectionSettings& settings) {
    WindTunnelProjectedCoarseGrid result;
    static_cast<void>(evaluateProjection(settings, &result));
    return result;
}

WindTunnelMomentumStepResult evaluateWindTunnelMomentumAdvance(
    const WindTunnelProjectionSettings& settings) {
    validateProjectionSettings(settings);
    if (!amrex::Initialized()) {
        throw std::logic_error("AMReX runtime is not initialized");
    }
    const WindTunnelBoundaryDiagnostics hierarchyDiagnostics =
        evaluateWindTunnelBoundaryInitialization(settings.grid);
    Hierarchy hierarchy = makeHierarchy(settings.grid);
    WindTunnelProjectedCoarseGrid before;
    WindTunnelMomentumStepResult result;
    result.diagnostics.initialProjection = projectHierarchy(
        settings, hierarchy, hierarchyDiagnostics, &before, true);
    result.diagnostics.maximumOutgoingCourantNumber =
        advectMomentumPredictor(hierarchy, settings);
    result.diagnostics.correctedProjection = projectHierarchy(
        settings, hierarchy, hierarchyDiagnostics,
        &result.projectedCoarseGrid, false);
    result.diagnostics.maximumCellVelocityChangeMetersPerSecond =
        maximumCellVelocityChange(before, result.projectedCoarseGrid);
    result.diagnostics.kineticEnergyBeforeJoules = coarseKineticEnergy(
        before, settings.airDensityKilogramsPerCubicMeter);
    result.diagnostics.kineticEnergyAfterJoules = coarseKineticEnergy(
        result.projectedCoarseGrid,
        settings.airDensityKilogramsPerCubicMeter);
    result.diagnostics.finite =
        std::isfinite(
            result.diagnostics.maximumOutgoingCourantNumber)
        && std::isfinite(
            result.diagnostics.maximumCellVelocityChangeMetersPerSecond)
        && std::isfinite(result.diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(result.diagnostics.kineticEnergyAfterJoules)
        && result.diagnostics.staticWing.finite;
    result.diagnostics.accepted =
        result.diagnostics.initialProjection.accepted
        && result.diagnostics.correctedProjection.accepted
        && result.diagnostics.finite
        && result.diagnostics.maximumOutgoingCourantNumber > 0.0
        && result.diagnostics.maximumOutgoingCourantNumber <= 1.0
        && result.diagnostics.maximumCellVelocityChangeMetersPerSecond > 0.0
        && result.diagnostics.kineticEnergyBeforeJoules > 0.0
        && result.diagnostics.kineticEnergyAfterJoules > 0.0;
    return result;
}

namespace {

amrex::Vector<FaceFields> cloneVelocity(const Hierarchy& hierarchy) {
    amrex::Vector<FaceFields> result;
    result.resize(hierarchy.velocity.size());
    for (std::size_t level = 0;
         level < hierarchy.velocity.size(); ++level) {
        for (int axis = 0; axis < 3; ++axis) {
            const amrex::MultiFab& source =
                *hierarchy.velocity[level][axis];
            result[level][axis] = std::make_unique<amrex::MultiFab>(
                source.boxArray(), source.DistributionMap(),
                source.nComp(), source.nGrowVect());
            amrex::MultiFab::Copy(
                *result[level][axis], source,
                0, 0, source.nComp(), source.nGrowVect());
        }
    }
    return result;
}

void finishMomentumDiagnostics(
    const WindTunnelProjectionSettings& settings,
    const WindTunnelProjectedCoarseGrid& before,
    WindTunnelMomentumStepResult& result) {
    result.diagnostics.maximumCellVelocityChangeMetersPerSecond =
        maximumCellVelocityChange(before, result.projectedCoarseGrid);
    result.diagnostics.kineticEnergyBeforeJoules = coarseKineticEnergy(
        before, settings.airDensityKilogramsPerCubicMeter);
    result.diagnostics.kineticEnergyAfterJoules = coarseKineticEnergy(
        result.projectedCoarseGrid,
        settings.airDensityKilogramsPerCubicMeter);
    result.diagnostics.finite =
        std::isfinite(
            result.diagnostics.maximumOutgoingCourantNumber)
        && std::isfinite(
            result.diagnostics.maximumCellVelocityChangeMetersPerSecond)
        && std::isfinite(result.diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(result.diagnostics.kineticEnergyAfterJoules)
        && result.diagnostics.staticWing.finite;
    result.diagnostics.accepted =
        result.diagnostics.initialProjection.accepted
        && result.diagnostics.correctedProjection.accepted
        && result.diagnostics.finite
        && result.diagnostics.maximumOutgoingCourantNumber > 0.0
        && result.diagnostics.maximumOutgoingCourantNumber <= 1.0
        && result.diagnostics.maximumCellVelocityChangeMetersPerSecond > 0.0
        && result.diagnostics.kineticEnergyBeforeJoules > 0.0
        && result.diagnostics.kineticEnergyAfterJoules > 0.0
        && result.diagnostics.staticWing.accepted;
}

} // namespace

struct WindTunnelMomentumState::Implementation {
    WindTunnelProjectionSettings settings;
    WindTunnelBoundaryDiagnostics hierarchyDiagnostics;
    Hierarchy hierarchy;
    WindTunnelProjectedCoarseGrid current;
    std::unique_ptr<StaticWingInterface> staticWing;

    explicit Implementation(WindTunnelProjectionSettings requestedSettings,
                            const Scene* requestedStaticWing,
                            const bool clipStaticWingToWindTunnel = false)
        : settings(std::move(requestedSettings)) {
        validateProjectionSettings(settings);
        if (!amrex::Initialized()) {
            throw std::logic_error("AMReX runtime is not initialized");
        }
        hierarchyDiagnostics =
            evaluateWindTunnelBoundaryInitialization(settings.grid);
        if (requestedStaticWing != nullptr) {
            staticWing = std::make_unique<StaticWingInterface>(
                *requestedStaticWing, settings.grid,
                clipStaticWingToWindTunnel);
        }
        hierarchy = makeHierarchy(settings.grid);
        const auto initial = projectHierarchy(
            settings, hierarchy, hierarchyDiagnostics, &current, true);
        if (!initial.accepted) {
            throw std::runtime_error(
                "AMR momentum state initial projection was rejected");
        }
    }

    WindTunnelMomentumStepResult advance() {
        amrex::Vector<FaceFields> savedVelocity = cloneVelocity(hierarchy);
        const WindTunnelProjectedCoarseGrid before = current;
        try {
            WindTunnelMomentumStepResult result;
            result.diagnostics.initialProjection = before.diagnostics;
            result.diagnostics.maximumOutgoingCourantNumber =
                advectMomentumPredictor(hierarchy, settings);
            if (staticWing) {
                std::vector<double> accumulatedPressure;
                std::size_t accumulatedProjectionIterations = 0;
                for (std::size_t iteration = 0;
                     iteration
                         < settings.staticWingForcingProjectionIterations;
                     ++iteration) {
                    auto forcing = applyStaticWingDirectForcing(
                        hierarchy, settings, *staticWing);
                    if (!forcing.accepted) {
                        throw std::runtime_error(
                            "AMR static-wing direct forcing was rejected");
                    }
                    if (iteration == 0) {
                        result.diagnostics.staticWing = forcing;
                    } else {
                        result.diagnostics.staticWing
                            .maximumSurfaceNormalSpeedAfterForcingMetersPerSecond =
                                forcing
                                    .maximumSurfaceNormalSpeedAfterForcingMetersPerSecond;
                    }
                    result.diagnostics.correctedProjection = projectHierarchy(
                        settings, hierarchy, hierarchyDiagnostics,
                        &result.projectedCoarseGrid, false);
                    if (accumulatedPressure.empty()) {
                        accumulatedPressure.assign(
                            result.projectedCoarseGrid.pressurePascals.size(),
                            0.0);
                    }
                    if (accumulatedPressure.size()
                        != result.projectedCoarseGrid
                               .pressurePascals.size()) {
                        throw std::logic_error(
                            "AMR static-wing pressure snapshots changed shape");
                    }
                    for (std::size_t cell = 0;
                         cell < accumulatedPressure.size(); ++cell) {
                        accumulatedPressure[cell] +=
                            result.projectedCoarseGrid
                                .pressurePascals[cell];
                    }
                    accumulatedProjectionIterations +=
                        result.diagnostics.correctedProjection
                            .solverIterations;
                    result.diagnostics.staticWing
                        .maximumSurfaceNormalSpeedAfterProjectionMetersPerSecond =
                            maximumStaticWingNormalSpeed(
                                hierarchy, *staticWing);
                    ++result.diagnostics.staticWing
                          .forcingProjectionIterations;
                }
                result.projectedCoarseGrid.pressurePascals =
                    std::move(accumulatedPressure);
                result.diagnostics.correctedProjection.solverIterations =
                    accumulatedProjectionIterations;
                result.diagnostics.correctedProjection
                    .maximumPressureCorrectionPascals = 0.0;
                for (const double pressure :
                     result.projectedCoarseGrid.pressurePascals) {
                    result.diagnostics.correctedProjection
                        .maximumPressureCorrectionPascals = std::max(
                            result.diagnostics.correctedProjection
                                .maximumPressureCorrectionPascals,
                            std::abs(pressure));
                }
                result.projectedCoarseGrid.diagnostics =
                    result.diagnostics.correctedProjection;
                result.diagnostics.staticWing.finite =
                    result.diagnostics.staticWing.finite
                    && std::isfinite(
                        result.diagnostics.correctedProjection
                            .maximumPressureCorrectionPascals)
                    && std::isfinite(
                        result.diagnostics.staticWing
                            .maximumSurfaceNormalSpeedAfterProjectionMetersPerSecond);
                result.diagnostics.staticWing.accepted =
                    result.diagnostics.staticWing.accepted
                    && result.diagnostics.staticWing.finite
                    && result.diagnostics.staticWing
                           .forcingProjectionIterations
                        == settings.staticWingForcingProjectionIterations
                    && result.diagnostics.staticWing
                           .maximumSurfaceNormalSpeedAfterProjectionMetersPerSecond
                        < result.diagnostics.staticWing
                              .maximumSurfaceNormalSpeedBeforeMetersPerSecond;
            } else {
                result.diagnostics.correctedProjection = projectHierarchy(
                    settings, hierarchy, hierarchyDiagnostics,
                    &result.projectedCoarseGrid, false);
            }
            finishMomentumDiagnostics(settings, before, result);
            if (!result.diagnostics.accepted) {
                throw std::runtime_error(
                    "AMR momentum state advance was rejected");
            }
            current = result.projectedCoarseGrid;
            return result;
        } catch (...) {
            hierarchy.velocity = std::move(savedVelocity);
            throw;
        }
    }
};

WindTunnelMomentumState::WindTunnelMomentumState(
    WindTunnelProjectionSettings settings)
    : implementation_(
          std::make_unique<Implementation>(std::move(settings), nullptr)) {}

WindTunnelMomentumState::WindTunnelMomentumState(
    WindTunnelProjectionSettings settings,
    const Scene& staticWingScene,
    const bool clipStaticWingToWindTunnel)
    : implementation_(std::make_unique<Implementation>(
          std::move(settings), &staticWingScene,
          clipStaticWingToWindTunnel)) {}

WindTunnelMomentumState::~WindTunnelMomentumState() = default;

WindTunnelMomentumStepResult WindTunnelMomentumState::advance() {
    return implementation_->advance();
}

const WindTunnelProjectedCoarseGrid&
WindTunnelMomentumState::projectedCoarseGrid() const noexcept {
    return implementation_->current;
}

const StaticWingInterface*
WindTunnelMomentumState::staticWingInterface() const noexcept {
    return implementation_->staticWing.get();
}

} // namespace simwing::fsi::amr

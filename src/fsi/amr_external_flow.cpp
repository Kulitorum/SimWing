#include "amr_external_flow.h"

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::amr {
namespace {

using FaceFields = std::array<std::unique_ptr<amrex::MultiFab>, 3>;

bool multiplyBounded(const std::size_t first,
                     const std::size_t second,
                     std::size_t& product) noexcept {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        return false;
    }
    product = first * second;
    return true;
}

void validateSettings(const WindTunnelGridSettings& settings) {
    const auto counts = settings.coarseCellCounts;
    std::size_t xy = 0;
    std::size_t xyz = 0;
    if (counts.x < 8 || counts.y < 8 || counts.z < 8
        || !multiplyBounded(counts.x, counts.y, xy)
        || !multiplyBounded(xy, counts.z, xyz)
        || xyz > 16'777'216
        || settings.maximumGridSize < 4
        || settings.maximumGridSize > 256
        || settings.refinementRatio != 2
        || settings.refinedCoarseCellLower[0] == 0
        || settings.refinedCoarseCellLower[1] == 0
        || settings.refinedCoarseCellLower[2] == 0
        || settings.refinedCoarseCellUpperExclusive[0] >= counts.x
        || settings.refinedCoarseCellUpperExclusive[1] >= counts.y
        || settings.refinedCoarseCellUpperExclusive[2] >= counts.z
        || settings.refinedCoarseCellLower[0]
            >= settings.refinedCoarseCellUpperExclusive[0]
        || settings.refinedCoarseCellLower[1]
            >= settings.refinedCoarseCellUpperExclusive[1]
        || settings.refinedCoarseCellLower[2]
            >= settings.refinedCoarseCellUpperExclusive[2]
        || !std::isfinite(settings.lowerMeters.x)
        || !std::isfinite(settings.lowerMeters.y)
        || !std::isfinite(settings.lowerMeters.z)
        || !std::isfinite(settings.upperMeters.x)
        || !std::isfinite(settings.upperMeters.y)
        || !std::isfinite(settings.upperMeters.z)
        || !(settings.upperMeters.x > settings.lowerMeters.x)
        || !(settings.upperMeters.y > settings.lowerMeters.y)
        || !(settings.upperMeters.z > settings.lowerMeters.z)
        || !std::isfinite(settings.freestreamMetersPerSecond.x)
        || !std::isfinite(settings.freestreamMetersPerSecond.y)
        || !std::isfinite(settings.freestreamMetersPerSecond.z)
        || !(settings.freestreamMetersPerSecond.y > 0.0)) {
        throw std::invalid_argument(
            "AMR wind-tunnel grid settings are invalid");
    }
}

amrex::IntVect intVect(const fluid::GridCellCounts counts) {
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

FaceFields makeFaceFields(const amrex::BoxArray& cells,
                          const amrex::DistributionMapping& distribution,
                          const fluid::Vector3 freestream) {
    FaceFields result;
    const std::array<double, 3> values{
        freestream.x, freestream.y, freestream.z};
    for (int axis = 0; axis < 3; ++axis) {
        amrex::BoxArray faces = amrex::convert(
            cells, amrex::IntVect::TheDimensionVector(axis));
        result[static_cast<std::size_t>(axis)] =
            std::make_unique<amrex::MultiFab>(
                faces, distribution, 1, 1);
        result[static_cast<std::size_t>(axis)]->setVal(values[axis]);
    }
    return result;
}

void initializeInteriorWake(FaceFields& fields,
                            const amrex::Box& cellDomain,
                            const fluid::Vector3 freestream) {
    amrex::MultiFab& streamwise = *fields[1];
    const amrex::Box faceDomain = amrex::convert(
        cellDomain, amrex::IntVect::TheDimensionVector(1));
    const double amplitude = 0.08 * freestream.y;
    const double inverseX = 1.0
        / static_cast<double>(std::max(1, cellDomain.length(0)));
    const double inverseY = 1.0
        / static_cast<double>(std::max(1, faceDomain.length(1) - 1));
    const double inverseZ = 1.0
        / static_cast<double>(std::max(1, cellDomain.length(2)));
    for (amrex::MFIter iterator(streamwise); iterator.isValid(); ++iterator) {
        const amrex::Box box = iterator.validbox();
        auto values = streamwise.array(iterator);
        const int xLower = cellDomain.smallEnd(0);
        const int yLower = faceDomain.smallEnd(1);
        const int zLower = cellDomain.smallEnd(2);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(
                                   const int i,
                                   const int j,
                                   const int k) noexcept {
            const double x = (static_cast<double>(i - xLower) + 0.5)
                * inverseX;
            const double y = static_cast<double>(j - yLower) * inverseY;
            const double z = (static_cast<double>(k - zLower) + 0.5)
                * inverseZ;
            const double crossSection = std::sin(3.14159265358979323846 * x)
                * std::sin(3.14159265358979323846 * z);
            values(i, j, k) = freestream.y
                - amplitude * y * crossSection;
        });
    }
    amrex::Gpu::streamSynchronize();
}

void fillCoarsePhysicalBoundaries(
    FaceFields& fields,
    const amrex::Geometry& geometry,
    const amrex::Box& cellDomain,
    const fluid::Vector3 freestream) {
    const std::array<double, 3> farField{
        freestream.x, freestream.y, freestream.z};
    for (int axis = 0; axis < 3; ++axis) {
        amrex::MultiFab& field = *fields[static_cast<std::size_t>(axis)];
        field.FillBoundary(geometry.periodicity());
        const amrex::Box faceDomain = amrex::convert(
            cellDomain, amrex::IntVect::TheDimensionVector(axis));
        for (amrex::MFIter iterator(field); iterator.isValid(); ++iterator) {
            const amrex::Box allocated = iterator.fabbox();
            auto values = field.array(iterator);
            const int xLower = faceDomain.smallEnd(0);
            const int xUpper = faceDomain.bigEnd(0);
            const int yLower = faceDomain.smallEnd(1);
            const int yUpper = faceDomain.bigEnd(1);
            const int zLower = faceDomain.smallEnd(2);
            const int zUpper = faceDomain.bigEnd(2);
            const double prescribed = farField[axis];
            amrex::ParallelFor(
                allocated, [=] AMREX_GPU_DEVICE(
                               const int i,
                               const int j,
                               const int k) noexcept {
                    const bool outsideX = i < xLower || i > xUpper;
                    const bool outsideZ = k < zLower || k > zUpper;
                    if (outsideX || outsideZ || j < yLower) {
                        values(i, j, k) = prescribed;
                    } else if (j > yUpper) {
                        values(i, j, k) = values(i, yUpper, k);
                    }
                });
        }
    }
    amrex::Gpu::streamSynchronize();
}

std::size_t ownedBytes(const FaceFields& fields) {
    std::size_t result = 0;
    for (const auto& field : fields) {
        for (amrex::MFIter iterator(*field); iterator.isValid(); ++iterator) {
            const std::size_t bytes = (*field)[iterator].nBytesOwned();
            if (bytes > std::numeric_limits<std::size_t>::max() - result) {
                throw std::overflow_error(
                    "AMR velocity storage accounting overflowed");
            }
            result += bytes;
        }
    }
    return result;
}

std::size_t validCellCount(const amrex::BoxArray& boxes) {
    std::size_t result = 0;
    for (int index = 0; index < boxes.size(); ++index) {
        const auto points = boxes[index].numPts();
        if (points < 0
            || static_cast<std::uint64_t>(points)
                > std::numeric_limits<std::size_t>::max() - result) {
            throw std::overflow_error("AMR valid-cell count overflowed");
        }
        result += static_cast<std::size_t>(points);
    }
    return result;
}

void auditCoarseBoundaries(
    const FaceFields& fields,
    const amrex::Box& cellDomain,
    const fluid::Vector3 freestream,
    WindTunnelBoundaryDiagnostics& diagnostics) {
    const std::array<double, 3> expected{
        freestream.x, freestream.y, freestream.z};
    for (int axis = 0; axis < 3; ++axis) {
        const amrex::MultiFab& field =
            *fields[static_cast<std::size_t>(axis)];
        const amrex::Box faceDomain = amrex::convert(
            cellDomain, amrex::IntVect::TheDimensionVector(axis));
        for (amrex::MFIter iterator(field); iterator.isValid(); ++iterator) {
            const auto values = field.const_array(iterator);
            const amrex::Box allocated = iterator.fabbox();
            for (int k = allocated.smallEnd(2);
                 k <= allocated.bigEnd(2); ++k) {
                for (int j = allocated.smallEnd(1);
                     j <= allocated.bigEnd(1); ++j) {
                    for (int i = allocated.smallEnd(0);
                         i <= allocated.bigEnd(0); ++i) {
                        const bool outsideX = i < faceDomain.smallEnd(0)
                            || i > faceDomain.bigEnd(0);
                        const bool outsideZ = k < faceDomain.smallEnd(2)
                            || k > faceDomain.bigEnd(2);
                        if (outsideX || outsideZ) {
                            diagnostics.farFieldMaximumErrorMetersPerSecond =
                                std::max(
                                    diagnostics
                                        .farFieldMaximumErrorMetersPerSecond,
                                    std::abs(values(i, j, k)
                                             - expected[axis]));
                        } else if (j < faceDomain.smallEnd(1)) {
                            diagnostics
                                .lowerYInflowMaximumErrorMetersPerSecond =
                                std::max(
                                    diagnostics
                                        .lowerYInflowMaximumErrorMetersPerSecond,
                                    std::abs(values(i, j, k)
                                             - expected[axis]));
                        } else if (j > faceDomain.bigEnd(1)) {
                            diagnostics
                                .upperYOutflowGradientMaximumErrorMetersPerSecond =
                                std::max(
                                    diagnostics
                                        .upperYOutflowGradientMaximumErrorMetersPerSecond,
                                    std::abs(
                                        values(i, j, k)
                                        - values(
                                            i, faceDomain.bigEnd(1), k)));
                        } else if (axis == 1) {
                            diagnostics
                                .maximumInteriorWakePerturbationMetersPerSecond =
                                std::max(
                                    diagnostics
                                        .maximumInteriorWakePerturbationMetersPerSecond,
                                    std::abs(values(i, j, k)
                                             - freestream.y));
                        }
                    }
                }
            }
        }
    }
}

} // namespace

struct Runtime::Implementation {
    amrex::AMReX* instance = nullptr;
};

Runtime::Runtime(int& argc, char**& argv)
    : implementation_(std::make_unique<Implementation>()) {
    if (amrex::Initialized()) {
        throw std::logic_error(
            "SimWing AMReX runtime is already initialized");
    }
    implementation_->instance = amrex::Initialize(argc, argv, false);
    if (implementation_->instance == nullptr) {
        throw std::runtime_error("AMReX initialization failed");
    }
}

Runtime::~Runtime() {
    if (implementation_ != nullptr
        && implementation_->instance != nullptr) {
        amrex::Finalize(implementation_->instance);
    }
}

WindTunnelBoundaryDiagnostics evaluateWindTunnelBoundaryInitialization(
    const WindTunnelGridSettings& settings) {
    validateSettings(settings);
    if (!amrex::Initialized()) {
        throw std::logic_error("AMReX runtime is not initialized");
    }

    const amrex::IntVect counts = intVect(settings.coarseCellCounts);
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
    const amrex::Geometry coarseGeometry(
        coarseDomain, &physicalDomain, amrex::CoordSys::cartesian,
        nonPeriodic.data());

    amrex::BoxArray coarseBoxes(coarseDomain);
    coarseBoxes.maxSize(static_cast<int>(settings.maximumGridSize));
    const amrex::DistributionMapping coarseDistribution(coarseBoxes);
    FaceFields coarseVelocity = makeFaceFields(
        coarseBoxes, coarseDistribution,
        settings.freestreamMetersPerSecond);
    initializeInteriorWake(
        coarseVelocity, coarseDomain,
        settings.freestreamMetersPerSecond);
    fillCoarsePhysicalBoundaries(
        coarseVelocity, coarseGeometry, coarseDomain,
        settings.freestreamMetersPerSecond);

    const int ratio = static_cast<int>(settings.refinementRatio);
    const amrex::Box refinedDomain = amrex::refine(coarseDomain, ratio);
    const amrex::Geometry refinedGeometry(
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
    amrex::BoxArray refinedBoxes(
        amrex::refine(
            amrex::Box(refinementLower, refinementUpper), ratio));
    refinedBoxes.maxSize(
        static_cast<int>(settings.maximumGridSize
                         * settings.refinementRatio));
    const amrex::DistributionMapping refinedDistribution(refinedBoxes);
    FaceFields refinedVelocity = makeFaceFields(
        refinedBoxes, refinedDistribution,
        settings.freestreamMetersPerSecond);
    initializeInteriorWake(
        refinedVelocity, refinedDomain,
        settings.freestreamMetersPerSecond);
    for (auto& field : refinedVelocity) {
        field->FillBoundary(refinedGeometry.periodicity());
    }

    WindTunnelBoundaryDiagnostics diagnostics;
    diagnostics.coarseCellCounts = settings.coarseCellCounts;
    diagnostics.refinedDomainCellCounts = {
        settings.coarseCellCounts.x * settings.refinementRatio,
        settings.coarseCellCounts.y * settings.refinementRatio,
        settings.coarseCellCounts.z * settings.refinementRatio,
    };
    diagnostics.coarseCellSpacingMeters = {
        (settings.upperMeters.x - settings.lowerMeters.x)
            / static_cast<double>(settings.coarseCellCounts.x),
        (settings.upperMeters.y - settings.lowerMeters.y)
            / static_cast<double>(settings.coarseCellCounts.y),
        (settings.upperMeters.z - settings.lowerMeters.z)
            / static_cast<double>(settings.coarseCellCounts.z),
    };
    diagnostics.physicalBoundaries = positiveYWindTunnelBoundaries;
    diagnostics.coarseBlockCount =
        static_cast<std::size_t>(coarseBoxes.size());
    diagnostics.refinedBlockCount =
        static_cast<std::size_t>(refinedBoxes.size());
    diagnostics.coarseValidCellCount = validCellCount(coarseBoxes);
    diagnostics.refinedValidCellCount = validCellCount(refinedBoxes);
    diagnostics.allocatedVelocityBytes = ownedBytes(coarseVelocity)
        + ownedBytes(refinedVelocity);
    auditCoarseBoundaries(
        coarseVelocity, coarseDomain,
        settings.freestreamMetersPerSecond, diagnostics);
    diagnostics.allPhysicalDirectionsNonPeriodic =
        !coarseGeometry.isPeriodic(0)
        && !coarseGeometry.isPeriodic(1)
        && !coarseGeometry.isPeriodic(2)
        && !refinedGeometry.isPeriodic(0)
        && !refinedGeometry.isPeriodic(1)
        && !refinedGeometry.isPeriodic(2);
    diagnostics.finite =
        std::isfinite(
            diagnostics.maximumInteriorWakePerturbationMetersPerSecond)
        && std::isfinite(
            diagnostics.lowerYInflowMaximumErrorMetersPerSecond)
        && std::isfinite(
            diagnostics.upperYOutflowGradientMaximumErrorMetersPerSecond)
        && std::isfinite(diagnostics.farFieldMaximumErrorMetersPerSecond);
    diagnostics.accepted = diagnostics.finite
        && diagnostics.allPhysicalDirectionsNonPeriodic
        && diagnostics.physicalBoundaries
            == positiveYWindTunnelBoundaries
        && diagnostics.coarseBlockCount > 1
        && diagnostics.refinedBlockCount > 0
        && diagnostics.coarseValidCellCount > 0
        && diagnostics.refinedValidCellCount > 0
        && diagnostics.allocatedVelocityBytes > 0
        && diagnostics.maximumInteriorWakePerturbationMetersPerSecond > 0.0
        && diagnostics.lowerYInflowMaximumErrorMetersPerSecond == 0.0
        && diagnostics.upperYOutflowGradientMaximumErrorMetersPerSecond
            == 0.0
        && diagnostics.farFieldMaximumErrorMetersPerSecond == 0.0;
    return diagnostics;
}

} // namespace simwing::fsi::amr

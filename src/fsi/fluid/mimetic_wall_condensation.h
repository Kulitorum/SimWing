#pragma once

#include "mimetic_local_cell.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t mimeticWallCondensationVersion = 1;

struct MimeticWallCondensationSettings {
    std::size_t maximumHalfFaces = 4096;
    std::size_t maximumOwnedBytes = 256ULL * 1024ULL * 1024ULL;

    bool operator==(
        const MimeticWallCondensationSettings&) const = default;
};

// Exact local Schur data for eliminating prescribed-zero-flux material-wall
// trace equations from one already condensed cell-scalar operator. The local
// trace matrix is diagonal plus seven signed low-rank modes, so the wall block
// inverse normally retains only one equilibrated 7x7 Woodbury inverse plus
// linear face data. If that auxiliary inverse is numerically singular but a
// bounded wall principal block of at most eight faces remains invertible, its
// exact Schur metric is retained instead and later wall solves reconstruct the
// small block from the immutable local operator. Active traces remain in their
// original local indices. A cell with no walls is the exact identity case; an
// all-wall cell rejects because its local constant mode needs a global gauge
// rather than a wall Schur inverse.
struct MimeticWallCondensation {
    std::uint32_t version = mimeticWallCondensationVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t localOperatorFingerprint = 0;
    std::size_t halfFaceCount = 0;
    std::size_t wallHalfFaceCount = 0;
    std::size_t activeHalfFaceCount = 0;
    std::size_t ownedStorageBytes = 0;
    double conservationDenominator = 0.0;
    std::vector<std::uint8_t> wallMask;
    std::vector<double> conservationCoupling;
    std::vector<double> condensedOperatorDiagonal;
    std::array<double, 49> inverseWoodburyCore{};
    std::array<double, 49> wallSchurMetric{};

    bool operator==(const MimeticWallCondensation&) const = default;
};

[[nodiscard]] MimeticWallCondensation buildMimeticWallCondensation(
    const MimeticLocalCellOperator& localOperator,
    std::span<const std::uint8_t> wallMask,
    const MimeticWallCondensationSettings& settings = {});

void validateMimeticWallCondensation(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator);

// Applies the active-active Schur complement. Wall entries of active traces
// must be exactly zero and remain zero in the returned full-local-index vector.
[[nodiscard]] std::vector<double>
applyMimeticWallCondensedTraceOperator(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator,
    std::span<const double> activeTraceValues);

// Condenses a full local trace right-hand side into active rows. Returned wall
// entries are zero.
[[nodiscard]] std::vector<double>
condenseMimeticWallTraceRightHandSide(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator,
    std::span<const double> fullRightHandSide);

// Reconstructs wall trace values after the active Schur system has been solved.
[[nodiscard]] std::vector<double> reconstructMimeticWallTraces(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator,
    std::span<const double> fullRightHandSide,
    std::span<const double> activeTraceValues);

} // namespace simwing::fsi::fluid

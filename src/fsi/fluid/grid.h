#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    bool operator==(const Vector3&) const = default;
};

struct GridCellCounts {
    std::size_t x = 0;
    std::size_t y = 0;
    std::size_t z = 0;

    bool operator==(const GridCellCounts&) const = default;
};

// A uniform right-handed Cartesian cell grid. The first verification backend
// is periodic in all directions; keeping the geometry object separate from
// fields makes the later block/AMR ownership boundary explicit.
class PeriodicCartesianGrid final {
public:
    PeriodicCartesianGrid(GridCellCounts cellCounts,
                          Vector3 lowerMeters,
                          Vector3 upperMeters);

    [[nodiscard]] GridCellCounts cellCounts() const noexcept;
    [[nodiscard]] Vector3 lowerMeters() const noexcept;
    [[nodiscard]] Vector3 upperMeters() const noexcept;
    [[nodiscard]] Vector3 cellSpacingMeters() const noexcept;
    [[nodiscard]] std::size_t cellCount() const noexcept;
    [[nodiscard]] double cellVolumeCubicMeters() const noexcept;

    [[nodiscard]] std::size_t cellIndex(std::size_t i,
                                        std::size_t j,
                                        std::size_t k) const;
    [[nodiscard]] Vector3 cellCenterMeters(std::size_t i,
                                           std::size_t j,
                                           std::size_t k) const;
    [[nodiscard]] Vector3 xFaceCenterMeters(std::size_t i,
                                            std::size_t j,
                                            std::size_t k) const;
    [[nodiscard]] Vector3 yFaceCenterMeters(std::size_t i,
                                            std::size_t j,
                                            std::size_t k) const;
    [[nodiscard]] Vector3 zFaceCenterMeters(std::size_t i,
                                            std::size_t j,
                                            std::size_t k) const;

private:
    GridCellCounts cellCounts_;
    Vector3 lowerMeters_;
    Vector3 upperMeters_;
    Vector3 cellSpacingMeters_;
    std::size_t cellCount_ = 0;
};

// Cell-centred scalar samples. A field is shape-bound but deliberately does
// not own a grid or allocation policy.
class CellScalarField final {
public:
    explicit CellScalarField(const PeriodicCartesianGrid& grid,
                             double initialValue = 0.0);

    [[nodiscard]] GridCellCounts cellCounts() const noexcept;
    [[nodiscard]] bool matches(const PeriodicCartesianGrid& grid) const noexcept;
    [[nodiscard]] std::span<double> values() noexcept;
    [[nodiscard]] std::span<const double> values() const noexcept;

    bool operator==(const CellScalarField&) const = default;

private:
    GridCellCounts cellCounts_;
    std::vector<double> values_;
};

// A periodic staggered (MAC) velocity field. Each component stores one unique
// face per cell; the closing face is the corresponding face at index zero.
class MacVelocityField final {
public:
    explicit MacVelocityField(const PeriodicCartesianGrid& grid,
                              double initialValue = 0.0);

    [[nodiscard]] GridCellCounts cellCounts() const noexcept;
    [[nodiscard]] bool matches(const PeriodicCartesianGrid& grid) const noexcept;
    [[nodiscard]] std::span<double> xFaces() noexcept;
    [[nodiscard]] std::span<const double> xFaces() const noexcept;
    [[nodiscard]] std::span<double> yFaces() noexcept;
    [[nodiscard]] std::span<const double> yFaces() const noexcept;
    [[nodiscard]] std::span<double> zFaces() noexcept;
    [[nodiscard]] std::span<const double> zFaces() const noexcept;

    bool operator==(const MacVelocityField&) const = default;

private:
    GridCellCounts cellCounts_;
    std::array<std::vector<double>, 3> components_;
};

[[nodiscard]] bool isFinite(const CellScalarField& field) noexcept;
[[nodiscard]] bool isFinite(const MacVelocityField& field) noexcept;

// These operators use the same periodic finite-volume stencil as projection.
// Their shared definition is a load-bearing discrete conservation contract.
void computeDivergence(const PeriodicCartesianGrid& grid,
                       const MacVelocityField& velocity,
                       CellScalarField& divergencePerSecond);
void computePressureGradient(const PeriodicCartesianGrid& grid,
                             const CellScalarField& pressurePascals,
                             MacVelocityField& gradientPascalsPerMeter);
void applyNegativeLaplacian(const PeriodicCartesianGrid& grid,
                            const CellScalarField& input,
                            CellScalarField& output);

[[nodiscard]] double mean(const CellScalarField& field);
[[nodiscard]] double l2Norm(const CellScalarField& field);
[[nodiscard]] double maximumAbsoluteValue(const CellScalarField& field);
[[nodiscard]] double kineticEnergyJoules(const PeriodicCartesianGrid& grid,
                                          const MacVelocityField& velocity,
                                          double densityKgPerCubicMeter);

} // namespace simwing::fsi::fluid

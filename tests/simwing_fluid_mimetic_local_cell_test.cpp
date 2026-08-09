#include "fluid/mimetic_local_cell.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using namespace simwing::fsi::fluid;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(const double actual,
               const double expected,
               const double tolerance,
               const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g expected %.17g tolerance %.3g)\n",
                     message, actual, expected, tolerance);
        ++failures;
    }
}

template<typename Callback>
void expectInvalid(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

template<typename Callback>
void expectLimited(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected, message);
}

Vector3 add(const Vector3& first, const Vector3& second) {
    return {
        first.x + second.x,
        first.y + second.y,
        first.z + second.z,
    };
}

Vector3 subtract(const Vector3& first, const Vector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

Vector3 scale(const Vector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double dot(const Vector3& first, const Vector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

Vector3 cross(const Vector3& first, const Vector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double length(const Vector3& value) {
    return std::sqrt(dot(value, value));
}

MimeticLocalCellGeometry cuboid(const Vector3 center,
                                const Vector3 dimensions) {
    const double xArea = dimensions.y * dimensions.z;
    const double yArea = dimensions.x * dimensions.z;
    const double zArea = dimensions.x * dimensions.y;
    MimeticLocalCellGeometry result;
    result.volumeCubicMeters = dimensions.x * dimensions.y * dimensions.z;
    result.centroidMeters = center;
    result.halfFaces = {
        {xArea, {center.x - 0.5 * dimensions.x, center.y, center.z},
         {-1.0, 0.0, 0.0}},
        {xArea, {center.x + 0.5 * dimensions.x, center.y, center.z},
         {1.0, 0.0, 0.0}},
        {yArea, {center.x, center.y - 0.5 * dimensions.y, center.z},
         {0.0, -1.0, 0.0}},
        {yArea, {center.x, center.y + 0.5 * dimensions.y, center.z},
         {0.0, 1.0, 0.0}},
        {zArea, {center.x, center.y, center.z - 0.5 * dimensions.z},
         {0.0, 0.0, -1.0}},
        {zArea, {center.x, center.y, center.z + 0.5 * dimensions.z},
         {0.0, 0.0, 1.0}},
    };
    return result;
}

MimeticLocalCellGeometry tetrahedron() {
    const std::array<Vector3, 4> vertices{{
        {0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
        {0.3, 1.5, 0.0},
        {0.2, 0.4, 1.2},
    }};
    MimeticLocalCellGeometry result;
    for (const auto& vertex : vertices) {
        result.centroidMeters = add(
            result.centroidMeters, scale(vertex, 0.25));
    }
    result.volumeCubicMeters = std::abs(dot(
        subtract(vertices[1], vertices[0]),
        cross(subtract(vertices[2], vertices[0]),
              subtract(vertices[3], vertices[0])))) / 6.0;
    const std::array<std::array<std::size_t, 3>, 4> faces{{
        {1, 2, 3},
        {0, 3, 2},
        {0, 1, 3},
        {0, 2, 1},
    }};
    for (const auto& indices : faces) {
        const Vector3 first = vertices[indices[0]];
        const Vector3 second = vertices[indices[1]];
        const Vector3 third = vertices[indices[2]];
        const Vector3 centroid = scale(
            add(add(first, second), third), 1.0 / 3.0);
        Vector3 areaVector = scale(
            cross(subtract(second, first), subtract(third, first)), 0.5);
        if (dot(areaVector, subtract(result.centroidMeters, centroid)) > 0.0) {
            areaVector = scale(areaVector, -1.0);
        }
        const double area = length(areaVector);
        result.halfFaces.push_back({
            area, centroid, scale(areaVector, 1.0 / area),
        });
    }
    return result;
}

double scalarAt(const double constant,
                const Vector3 gradient,
                const Vector3 position) {
    return constant + dot(gradient, position);
}

bool positiveDefinite(const MimeticLocalCellOperator& localOperator) {
    const std::size_t count = localOperator.halfFaceCount;
    std::vector<double> lower(count * count, 0.0);
    double maximumDiagonal = 0.0;
    for (std::size_t row = 0; row < count; ++row) {
        maximumDiagonal = std::max(
            maximumDiagonal,
            localOperator.inverseFluxInnerProduct[row * count + row]);
    }
    for (std::size_t row = 0; row < count; ++row) {
        for (std::size_t column = 0; column <= row; ++column) {
            double value = localOperator.inverseFluxInnerProduct[
                row * count + column];
            for (std::size_t inner = 0; inner < column; ++inner) {
                value -= lower[row * count + inner]
                    * lower[column * count + inner];
            }
            if (row == column) {
                if (value <= maximumDiagonal * 1.0e-12) {
                    return false;
                }
                lower[row * count + column] = std::sqrt(value);
            } else {
                lower[row * count + column] = value
                    / lower[column * count + column];
            }
        }
    }
    return true;
}

void testCartesianEquivalence() {
    const Vector3 dimensions{2.0, 4.0, 6.0};
    const auto geometry = cuboid({}, dimensions);
    const auto localOperator = buildMimeticLocalCellOperator(geometry);
    const auto translated = buildMimeticLocalCellOperator(
        cuboid({8.0, -4.0, 2.0}, dimensions));
    check(localOperator == translated,
          "Cartesian local operator is exactly translation invariant");
    check(localOperator.version == mimeticLocalCellVersion
              && localOperator.fingerprint != 0
              && localOperator.halfFaceCount == 6
              && positiveDefinite(localOperator),
          "Cartesian local operator is bounded and SPD");
    const double expectedDiagonal = 2.0 / geometry.volumeCubicMeters;
    for (std::size_t row = 0; row < 6; ++row) {
        for (std::size_t column = 0; column < 6; ++column) {
            checkNear(
                localOperator.inverseFluxInnerProduct[row * 6 + column],
                row == column ? expectedDiagonal : 0.0,
                2.0e-16,
                "Cartesian mimetic matrix is the exact diagonal half-cell stencil");
        }
    }

    const double constant = 1.25;
    const Vector3 gradient{2.0, -3.0, 0.5};
    std::vector<double> traces;
    for (const auto& face : geometry.halfFaces) {
        traces.push_back(scalarAt(constant, gradient, face.centroidMeters));
    }
    const auto fluxes = applyMimeticLocalNormalFlux(
        localOperator,
        scalarAt(constant, gradient, geometry.centroidMeters), traces);
    for (std::size_t face = 0; face < fluxes.size(); ++face) {
        checkNear(fluxes[face],
                  -dot(gradient, geometry.halfFaces[face].outwardUnitNormal),
                  2.0e-15,
                  "Cartesian mimetic flux is exact for a linear scalar");
    }

    const auto leftGeometry = cuboid({}, {2.0, 4.0, 6.0});
    const auto rightGeometry = cuboid({3.0, 0.0, 0.0}, {4.0, 4.0, 6.0});
    const auto left = buildMimeticLocalCellOperator(leftGeometry);
    const auto right = buildMimeticLocalCellOperator(rightGeometry);
    const double sharedArea = leftGeometry.halfFaces[1].areaSquareMeters;
    const double leftTraceConductance = sharedArea * sharedArea
        * left.inverseFluxInnerProduct[1 * 6 + 1];
    const double rightTraceConductance = sharedArea * sharedArea
        * right.inverseFluxInnerProduct[0 * 6 + 0];
    const double leftScalar = 5.0;
    const double rightScalar = 1.0;
    const double sharedTrace = (
        leftTraceConductance * leftScalar
        + rightTraceConductance * rightScalar)
        / (leftTraceConductance + rightTraceConductance);
    std::vector<double> leftTraces(6, leftScalar);
    std::vector<double> rightTraces(6, rightScalar);
    leftTraces[1] = sharedTrace;
    rightTraces[0] = sharedTrace;
    const double leftIntegratedFlux = sharedArea
        * applyMimeticLocalNormalFlux(
            left, leftScalar, leftTraces)[1];
    const double rightIntegratedFlux = sharedArea
        * applyMimeticLocalNormalFlux(
            right, rightScalar, rightTraces)[0];
    const double centerDistance = rightGeometry.centroidMeters.x
        - leftGeometry.centroidMeters.x;
    const double expectedFlux = sharedArea
        * (leftScalar - rightScalar) / centerDistance;
    checkNear(leftIntegratedFlux, expectedFlux, 3.0e-15,
              "condensed Cartesian trace equals area-over-centre-distance flux");
    checkNear(leftIntegratedFlux + rightIntegratedFlux, 0.0, 1.0e-14,
              "shared Cartesian half-face fluxes cancel conservatively");
}

void testTetrahedralConsistencyAndBalance() {
    const auto geometry = tetrahedron();
    const auto localOperator = buildMimeticLocalCellOperator(geometry);
    check(localOperator.halfFaceCount == 4
              && localOperator.maximumAreaClosureErrorSquareMeters < 5.0e-16
              && localOperator.maximumDivergenceTheoremErrorCubicMeters
                  < 5.0e-16
              && localOperator.maximumAlgebraicConsistencyError < 2.0e-15
              && positiveDefinite(localOperator),
          "skew tetrahedron builds a closed, consistent SPD local operator");
    const double constant = -0.75;
    const Vector3 gradient{0.8, -1.1, 0.35};
    std::vector<double> traces;
    for (const auto& face : geometry.halfFaces) {
        traces.push_back(scalarAt(constant, gradient, face.centroidMeters));
    }
    const double cellScalar = scalarAt(
        constant, gradient, geometry.centroidMeters);
    const auto fluxes = applyMimeticLocalNormalFlux(
        localOperator, cellScalar, traces);
    for (std::size_t face = 0; face < fluxes.size(); ++face) {
        checkNear(fluxes[face],
                  -dot(gradient, geometry.halfFaces[face].outwardUnitNormal),
                  3.0e-15,
                  "skew tetrahedron reproduces every linear normal flux");
    }
    const auto harmonic = balanceMimeticLocalCell(
        localOperator, traces, 0.0);
    checkNear(harmonic.cellScalar, cellScalar, 3.0e-15,
              "zero-source local balance reproduces the linear cell scalar");
    checkNear(harmonic.integratedOutwardFluxSum, 0.0, 2.0e-15,
              "zero-source tetrahedral flux closes locally");
    checkNear(harmonic.conservationResidual, 0.0, 2.0e-15,
              "zero-source tetrahedral conservation residual is roundoff only");

    const std::array<double, 4> arbitrary{{1.0, -0.25, 0.6, 2.0}};
    const double requestedSource = 0.37;
    const auto forced = balanceMimeticLocalCell(
        localOperator, arbitrary, requestedSource);
    checkNear(forced.integratedOutwardFluxSum, requestedSource, 3.0e-15,
              "forced tetrahedral balance preserves its integrated source");
    checkNear(forced.conservationResidual, 0.0, 3.0e-15,
              "forced tetrahedral conservation residual is roundoff only");

    const std::vector<double> constantTraces(4, 3.25);
    const auto constantBalance = balanceMimeticLocalCell(
        localOperator, constantTraces, 0.0);
    checkNear(constantBalance.cellScalar, 3.25, 2.0e-15,
              "constant face trace remains the exact local null mode");
    check(std::ranges::all_of(
              constantBalance.outwardNormalFluxes,
              [](const double value) { return std::abs(value) < 3.0e-15; }),
          "constant local null mode has no face flux");
}

void testRejectedGeometryAndCorruption() {
    const auto acceptedGeometry = cuboid({}, {2.0, 4.0, 6.0});
    const auto accepted = buildMimeticLocalCellOperator(acceptedGeometry);

    auto open = acceptedGeometry;
    open.halfFaces.pop_back();
    expectInvalid(
        [&] { static_cast<void>(buildMimeticLocalCellOperator(open)); },
        "local mimetic assembly rejects an open half-face shell");
    auto wrongVolume = acceptedGeometry;
    wrongVolume.volumeCubicMeters *= 1.1;
    expectInvalid(
        [&] { static_cast<void>(
            buildMimeticLocalCellOperator(wrongVolume)); },
        "local mimetic assembly rejects inconsistent volume geometry");
    auto wrongNormal = acceptedGeometry;
    wrongNormal.halfFaces.front().outwardUnitNormal.x = -2.0;
    expectInvalid(
        [&] { static_cast<void>(
            buildMimeticLocalCellOperator(wrongNormal)); },
        "local mimetic assembly rejects a non-unit normal");
    auto nonFinite = acceptedGeometry;
    nonFinite.halfFaces.front().centroidMeters.x =
        std::numeric_limits<double>::quiet_NaN();
    expectInvalid(
        [&] { static_cast<void>(
            buildMimeticLocalCellOperator(nonFinite)); },
        "local mimetic assembly rejects non-finite geometry");

    MimeticLocalCellSettings limits;
    limits.maximumHalfFaces = 5;
    expectLimited(
        [&] { static_cast<void>(
            buildMimeticLocalCellOperator(acceptedGeometry, limits)); },
        "local mimetic assembly bounds half-face count");
    limits = {};
    limits.maximumOperatorBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildMimeticLocalCellOperator(acceptedGeometry, limits)); },
        "local mimetic assembly bounds dense matrix storage");

    auto corrupt = accepted;
    corrupt.inverseFluxInnerProduct[1] += 0.01;
    corrupt.inverseFluxInnerProduct[6] += 0.01;
    expectInvalid(
        [&] { validateMimeticLocalCellOperator(corrupt); },
        "local mimetic validation rejects fingerprinted matrix corruption");
    expectInvalid(
        [&] { static_cast<void>(applyMimeticLocalNormalFlux(
            accepted, 0.0, std::vector<double>(5, 0.0))); },
        "local mimetic application rejects a short trace field");
    expectInvalid(
        [&] { static_cast<void>(balanceMimeticLocalCell(
            accepted, std::vector<double>(6, 0.0),
            std::numeric_limits<double>::infinity())); },
        "local mimetic balance rejects a non-finite source");
}

} // namespace

int main() {
    try {
        testCartesianEquivalence();
        testTetrahedralConsistencyAndBalance();
        testRejectedGeometryAndCorruption();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d mimetic local-cell check(s) failed\n", failures);
        return 1;
    }
    std::puts("all mimetic local-cell checks passed");
    return 0;
}

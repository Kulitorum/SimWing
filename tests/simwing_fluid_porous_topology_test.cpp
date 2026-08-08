#include "fluid/porous_topology.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

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
    check(std::abs(actual - expected) <= tolerance, message);
}

template <class Function>
void checkRejected(Function&& function, const char* message) {
    bool rejected = false;
    try {
        function();
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected, message);
}

PeriodicCartesianGrid grid() {
    return {{4, 5, 6}, {-2.0, -5.0, 10.0}, {2.0, 5.0, 16.0}};
}

void testInteriorAndAdjacentSelection() {
    const auto geometry = grid();
    const MovingPorousFaceTopology current{
        movingPorousFaceTopologyVersion, GridFaceAxis::X, 1, 0};
    const auto interior = selectMovingPorousTopology(
        geometry, current, -0.75);
    check(interior.topology == current
              && interior.rebaseDirection
                  == PorousTopologyRebaseDirection::None,
          "porous topology retains an interior X segment");
    checkNear(interior.crossingFraction, 0.75, 1.0e-15,
              "porous topology retains the exact interior fraction");

    const auto positive = selectMovingPorousTopology(
        geometry, current, 0.25);
    check(positive.topology.faceCoordinate == 2
              && positive.topology.periodicImage == 0
              && positive.rebaseDirection
                  == PorousTopologyRebaseDirection::Positive,
          "porous topology advances one positive X segment");
    checkNear(positive.crossingFraction, 0.75, 1.0e-15,
              "porous topology recomputes the positive crossing fraction");

    const auto negative = selectMovingPorousTopology(
        geometry, current, -1.75);
    check(negative.topology.faceCoordinate == 0
              && negative.topology.periodicImage == 0
              && negative.rebaseDirection
                  == PorousTopologyRebaseDirection::Negative,
          "porous topology advances one negative X segment");
    checkNear(negative.crossingFraction, 0.75, 1.0e-15,
              "porous topology recomputes the negative crossing fraction");
}

void testPeriodicWrapAndAxes() {
    const auto geometry = grid();
    const auto positiveWrap = selectMovingPorousTopology(
        geometry,
        {movingPorousFaceTopologyVersion, GridFaceAxis::X, 3, 0},
        2.25);
    check(positiveWrap.topology.faceCoordinate == 0
              && positiveWrap.topology.periodicImage == 1
              && positiveWrap.rebaseDirection
                  == PorousTopologyRebaseDirection::Positive,
          "porous topology preserves unwrapped position across positive wrap");
    checkNear(positiveWrap.crossingFraction, 0.75, 1.0e-15,
              "porous topology preserves fraction across positive wrap");

    const auto negativeWrap = selectMovingPorousTopology(
        geometry,
        {movingPorousFaceTopologyVersion, GridFaceAxis::X, 0, 0},
        -2.75);
    check(negativeWrap.topology.faceCoordinate == 3
              && negativeWrap.topology.periodicImage == -1
              && negativeWrap.rebaseDirection
                  == PorousTopologyRebaseDirection::Negative,
          "porous topology preserves unwrapped position across negative wrap");
    checkNear(negativeWrap.crossingFraction, 0.75, 1.0e-15,
              "porous topology preserves fraction across negative wrap");

    const auto y = selectMovingPorousTopology(
        geometry,
        {movingPorousFaceTopologyVersion, GridFaceAxis::Y, 2, 0},
        0.5);
    check(y.topology.faceCoordinate == 3
              && y.rebaseDirection
                  == PorousTopologyRebaseDirection::Positive,
          "porous topology selects an adjacent Y segment");
    checkNear(y.crossingFraction, 0.25, 1.0e-15,
              "porous topology uses Y spacing");

    const auto z = selectMovingPorousTopology(
        geometry,
        {movingPorousFaceTopologyVersion, GridFaceAxis::Z, 5, 0},
        16.25);
    check(z.topology.faceCoordinate == 0
              && z.topology.periodicImage == 1
              && z.rebaseDirection
                  == PorousTopologyRebaseDirection::Positive,
          "porous topology wraps a Z segment");
    checkNear(z.crossingFraction, 0.75, 1.0e-15,
              "porous topology uses Z spacing and image");
}

void testTransactionalRejection() {
    const auto geometry = grid();
    const MovingPorousFaceTopology current{
        movingPorousFaceTopologyVersion, GridFaceAxis::X, 1, 0};
    checkRejected(
        [&] { static_cast<void>(selectMovingPorousTopology(
            geometry, current, -0.5)); },
        "porous topology rejects exact segment boundaries");
    checkRejected(
        [&] { static_cast<void>(selectMovingPorousTopology(
            geometry, current, 1.0)); },
        "porous topology rejects skipped segments");
    checkRejected(
        [&] { static_cast<void>(selectMovingPorousTopology(
            geometry, current,
            std::numeric_limits<double>::quiet_NaN())); },
        "porous topology rejects non-finite physical position");
    checkRejected(
        [&] { static_cast<void>(movingPorousCrossingFraction(
            geometry,
            {movingPorousFaceTopologyVersion, GridFaceAxis::X, 4, 0},
            0.0)); },
        "porous topology rejects an out-of-grid face coordinate");
    checkRejected(
        [&] { static_cast<void>(movingPorousCrossingFraction(
            geometry,
            {movingPorousFaceTopologyVersion,
             static_cast<GridFaceAxis>(99), 0, 0},
            0.0)); },
        "porous topology rejects an invalid axis");
    checkRejected(
        [&] { static_cast<void>(movingPorousCrossingFraction(
            geometry,
            {movingPorousFaceTopologyVersion + 1,
             GridFaceAxis::X, 1, 0},
            -1.0)); },
        "porous topology rejects an unsupported version");

    const auto first = selectMovingPorousTopology(
        geometry, current, 0.25);
    const auto second = selectMovingPorousTopology(
        geometry, current, 0.25);
    check(first == second,
          "porous topology selection is deterministic and owning");
}

} // namespace

int main() {
    testInteriorAndAdjacentSelection();
    testPeriodicWrapAndAxes();
    testTransactionalRejection();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d moving porous topology check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all moving porous topology checks passed");
    return 0;
}

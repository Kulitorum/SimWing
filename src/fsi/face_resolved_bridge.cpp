#include "face_resolved_bridge.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace simwing::fsi {
namespace {

struct Vector2 {
    double u = 0.0;
    double v = 0.0;
};

[[nodiscard]] bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

[[nodiscard]] bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

[[nodiscard]] StructureVector3 toStructure(const fluid::Vector3& value) {
    return {value.x, value.y, value.z};
}

[[nodiscard]] StructureVector3 add(const StructureVector3& first,
                                   const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

[[nodiscard]] StructureVector3 subtract(const StructureVector3& first,
                                        const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

[[nodiscard]] StructureVector3 scale(const StructureVector3& value,
                                     const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

[[nodiscard]] StructureVector3 cross(const StructureVector3& first,
                                     const StructureVector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

[[nodiscard]] double dot(const StructureVector3& first,
                         const StructureVector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

[[nodiscard]] double length(const StructureVector3& value) {
    return std::hypot(value.x, value.y, value.z);
}

[[nodiscard]] int axisOrdinal(const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return 0;
    case fluid::GridFaceAxis::Y:
        return 1;
    case fluid::GridFaceAxis::Z:
        return 2;
    }
    throw std::invalid_argument("face-resolved bridge has an unknown face axis");
}

[[nodiscard]] StructureVector3 axisUnit(const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return {1.0, 0.0, 0.0};
    case fluid::GridFaceAxis::Y:
        return {0.0, 1.0, 0.0};
    case fluid::GridFaceAxis::Z:
        return {0.0, 0.0, 1.0};
    }
    throw std::invalid_argument("face-resolved bridge has an unknown face axis");
}

[[nodiscard]] double normalCoordinate(
    const StructureVector3& value,
    const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return value.x;
    case fluid::GridFaceAxis::Y:
        return value.y;
    case fluid::GridFaceAxis::Z:
        return value.z;
    }
    throw std::invalid_argument("face-resolved bridge has an unknown face axis");
}

void setNormalCoordinate(StructureVector3& value,
                         const fluid::GridFaceAxis axis,
                         const double coordinate) {
    switch (axis) {
    case fluid::GridFaceAxis::X:
        value.x = coordinate;
        return;
    case fluid::GridFaceAxis::Y:
        value.y = coordinate;
        return;
    case fluid::GridFaceAxis::Z:
        value.z = coordinate;
        return;
    }
    throw std::invalid_argument("face-resolved bridge has an unknown face axis");
}

[[nodiscard]] bool transverseIndicesMatch(
    const fluid::MovingInterfaceFaceDiagnostics& actual,
    const fluid::GridFaceAxis axis,
    const std::size_t expectedI,
    const std::size_t expectedJ,
    const std::size_t expectedK) {
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return actual.j == expectedJ && actual.k == expectedK;
    case fluid::GridFaceAxis::Y:
        return actual.i == expectedI && actual.k == expectedK;
    case fluid::GridFaceAxis::Z:
        return actual.i == expectedI && actual.j == expectedJ;
    }
    throw std::invalid_argument("face-resolved bridge has an unknown face axis");
}

[[nodiscard]] Vector2 project(const StructureVector3& value,
                              const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return {value.y, value.z};
    case fluid::GridFaceAxis::Y:
        return {value.z, value.x};
    case fluid::GridFaceAxis::Z:
        return {value.x, value.y};
    }
    throw std::invalid_argument("face-resolved bridge has an unknown face axis");
}

[[nodiscard]] Vector2 subtract(const Vector2& first, const Vector2& second) {
    return {first.u - second.u, first.v - second.v};
}

[[nodiscard]] double cross2(const Vector2& first, const Vector2& second) {
    return first.u * second.v - first.v * second.u;
}

[[nodiscard]] double signedTwiceArea(
    const std::vector<Vector2>& polygon) {
    double result = 0.0;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const auto& first = polygon[index];
        const auto& second = polygon[(index + 1) % polygon.size()];
        result += cross2(first, second);
    }
    return result;
}

[[nodiscard]] std::pair<double, Vector2> areaAndCentroid(
    const std::vector<Vector2>& polygon) {
    if (polygon.size() < 3) {
        return {};
    }
    const double twiceArea = signedTwiceArea(polygon);
    if (!(twiceArea > 0.0)) {
        return {};
    }
    Vector2 centroid;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const auto& first = polygon[index];
        const auto& second = polygon[(index + 1) % polygon.size()];
        const double weight = cross2(first, second);
        centroid.u += (first.u + second.u) * weight;
        centroid.v += (first.v + second.v) * weight;
    }
    centroid.u /= 3.0 * twiceArea;
    centroid.v /= 3.0 * twiceArea;
    return {0.5 * twiceArea, centroid};
}

template<typename Inside, typename Intersection>
[[nodiscard]] std::vector<Vector2> clipPolygon(
    const std::vector<Vector2>& input,
    const Inside& inside,
    const Intersection& intersection) {
    std::vector<Vector2> result;
    if (input.empty()) {
        return result;
    }
    Vector2 previous = input.back();
    bool previousInside = inside(previous);
    for (const Vector2 current : input) {
        const bool currentInside = inside(current);
        if (currentInside != previousInside) {
            result.push_back(intersection(previous, current));
        }
        if (currentInside) {
            result.push_back(current);
        }
        previous = current;
        previousInside = currentInside;
    }
    return result;
}

[[nodiscard]] std::vector<Vector2> clipToRectangle(
    std::vector<Vector2> polygon,
    const Vector2 lower,
    const Vector2 upper) {
    const auto clipU = [&](const double boundary, const bool keepGreater) {
        polygon = clipPolygon(
            polygon,
            [&](const Vector2 point) {
                return keepGreater ? point.u >= boundary
                                   : point.u <= boundary;
            },
            [&](const Vector2 first, const Vector2 second) {
                const double fraction = (boundary - first.u)
                    / (second.u - first.u);
                return Vector2{
                    boundary,
                    first.v + fraction * (second.v - first.v),
                };
            });
    };
    const auto clipV = [&](const double boundary, const bool keepGreater) {
        polygon = clipPolygon(
            polygon,
            [&](const Vector2 point) {
                return keepGreater ? point.v >= boundary
                                   : point.v <= boundary;
            },
            [&](const Vector2 first, const Vector2 second) {
                const double fraction = (boundary - first.v)
                    / (second.v - first.v);
                return Vector2{
                    first.u + fraction * (second.u - first.u),
                    boundary,
                };
            });
    };
    clipU(lower.u, true);
    clipU(upper.u, false);
    clipV(lower.v, true);
    clipV(upper.v, false);
    return polygon;
}

[[nodiscard]] std::vector<Vector2> clipToConvexPolygon(
    std::vector<Vector2> subject,
    const std::vector<Vector2>& clipper,
    const double tolerance) {
    for (std::size_t edge = 0; edge < clipper.size(); ++edge) {
        const Vector2 first = clipper[edge];
        const Vector2 second = clipper[(edge + 1) % clipper.size()];
        const Vector2 direction = subtract(second, first);
        subject = clipPolygon(
            subject,
            [&](const Vector2 point) {
                return cross2(direction, subtract(point, first)) >= -tolerance;
            },
            [&](const Vector2 start, const Vector2 end) {
                const Vector2 travel = subtract(end, start);
                const double denominator = cross2(direction, travel);
                if (std::abs(denominator) <= tolerance) {
                    return start;
                }
                const double fraction = cross2(
                    direction, subtract(first, start)) / denominator;
                return Vector2{
                    start.u + fraction * travel.u,
                    start.v + fraction * travel.v,
                };
            });
        if (subject.empty()) {
            break;
        }
    }
    return subject;
}

[[nodiscard]] std::array<double, 3> barycentric(
    const Vector2 point,
    const std::array<Vector2, 3>& triangle) {
    const double denominator = cross2(
        subtract(triangle[1], triangle[0]),
        subtract(triangle[2], triangle[0]));
    const double second = cross2(
        subtract(point, triangle[0]),
        subtract(triangle[2], triangle[0])) / denominator;
    const double third = cross2(
        subtract(triangle[1], triangle[0]),
        subtract(point, triangle[0])) / denominator;
    return {1.0 - second - third, second, third};
}

[[nodiscard]] double combinedTolerance(const double absoluteTolerance,
                                       const double relativeTolerance,
                                       const double firstMagnitude,
                                       const double secondMagnitude) {
    return absoluteTolerance
        + relativeTolerance * std::max(firstMagnitude, secondMagnitude);
}

void validateSettings(const PlanarFaceResolvedBridgeSettings& settings) {
    const double nonnegative[] = {
        settings.geometryToleranceMeters,
        settings.absoluteVelocityToleranceMetersPerSecond,
        settings.relativeVelocityTolerance,
        settings.minimumOverlapAreaSquareMeters,
        settings.absoluteAreaToleranceSquareMeters,
        settings.relativeAreaTolerance,
        settings.absoluteForceToleranceNewtons,
        settings.relativeForceTolerance,
        settings.absoluteMomentToleranceNewtonMeters,
        settings.relativeMomentTolerance,
        settings.absolutePowerToleranceWatts,
        settings.relativePowerTolerance,
    };
    if (!std::ranges::all_of(nonnegative, [](const double value) {
            return std::isfinite(value) && value >= 0.0;
        })
        || !(settings.minimumOverlapAreaSquareMeters > 0.0)
        || !finite(settings.transfer.momentReferenceMeters)
        || !std::isfinite(
            settings.transfer.minimumTriangleAreaSquareMeters)
        || !(settings.transfer.minimumTriangleAreaSquareMeters > 0.0)
        || !std::isfinite(
            settings.transfer.minimumQuadratureAreaSquareMeters)
        || !(settings.transfer.minimumQuadratureAreaSquareMeters > 0.0)
        || settings.minimumOverlapAreaSquareMeters
            < settings.transfer.minimumQuadratureAreaSquareMeters
        || !std::isfinite(settings.transfer.barycentricTolerance)
        || settings.transfer.barycentricTolerance < 0.0
        || !std::isfinite(settings.minimumNormalAlignment)
        || settings.minimumNormalAlignment < 0.0
        || settings.minimumNormalAlignment > 1.0
        || (settings.correspondenceMode
                != PlanarFaceCorrespondenceMode::FixedMaterial
            && settings.correspondenceMode
                != PlanarFaceCorrespondenceMode::RigidNormalTranslation)) {
        throw std::invalid_argument(
            "face-resolved bridge settings are invalid");
    }
}

[[nodiscard]] std::tuple<int, std::size_t, std::size_t, std::size_t>
faceKey(const fluid::MovingInterfaceFaceDiagnostics& face) {
    return {axisOrdinal(face.axis), face.k, face.j, face.i};
}

[[nodiscard]] StructureVector3 faceCenter(
    const fluid::MovingInterfaceFaceDiagnostics& face) {
    return scale(add(toStructure(face.lowerCornerMeters),
                     toStructure(face.upperCornerMeters)),
                 0.5);
}

[[nodiscard]] StructureVector3 physicalFaceCenter(
    const fluid::MovingInterfaceFaceDiagnostics& face,
    const fluid::GridFaceAxis axis,
    const std::optional<double> physicalPlaneCoordinateMeters) {
    StructureVector3 center = faceCenter(face);
    if (physicalPlaneCoordinateMeters.has_value()) {
        setNormalCoordinate(
            center, axis, *physicalPlaneCoordinateMeters);
    }
    return center;
}

} // namespace

PlanarFaceResolvedTransferResult::PlanarFaceResolvedTransferResult(
    ConservativeTransferResult transferResult,
    PlanarFaceResolvedBridgeDiagnostics diagnostics)
    : transferResult_(std::move(transferResult)),
      diagnostics_(std::move(diagnostics)) {}

const ConservativeTransferResult&
PlanarFaceResolvedTransferResult::transferResult() const noexcept {
    return transferResult_;
}

const PlanarFaceResolvedBridgeDiagnostics&
PlanarFaceResolvedTransferResult::diagnostics() const noexcept {
    return diagnostics_;
}

PlanarFaceResolvedFluidStructureBridge::
PlanarFaceResolvedFluidStructureBridge(
    const Structure& target,
    const std::uint64_t fluidSurfaceStableId,
    std::vector<CouplingSurfaceNodeDefinition> nodes,
    std::vector<CouplingSurfaceTriangleDefinition> triangles,
    std::vector<fluid::MovingInterfaceFaceDiagnostics> referenceFaces,
    const PlanarFaceResolvedBridgeSettings& settings)
    : fluidSurfaceStableId_(fluidSurfaceStableId),
      transfer_(target, std::move(nodes), std::move(triangles)),
      settings_(settings) {
    validateSettings(settings_);
    if (fluidSurfaceStableId_ == 0) {
        throw std::invalid_argument(
            "face-resolved bridge surface stable ID must be nonzero");
    }

    std::vector<fluid::MovingInterfaceFaceDiagnostics> selectedFaces;
    for (const auto& face : referenceFaces) {
        if (face.surfaceStableId == fluidSurfaceStableId_) {
            selectedFaces.push_back(face);
        }
    }
    if (selectedFaces.empty()) {
        throw std::invalid_argument(
            "reference fluid geometry lacks the bound surface stable ID");
    }
    const auto axis = selectedFaces.front().axis;
    const auto minusRegion = selectedFaces.front().minusRegionStableId;
    const auto plusRegion = selectedFaces.front().plusRegionStableId;
    const double planeCoordinate = normalCoordinate(
        toStructure(selectedFaces.front().lowerCornerMeters), axis);
    axis_ = axis;
    referencePlaneCoordinateMeters_ = planeCoordinate;
    std::tuple<int, std::size_t, std::size_t, std::size_t> previousKey;
    bool havePreviousKey = false;
    for (const auto& face : selectedFaces) {
        const auto key = faceKey(face);
        if ((havePreviousKey && !(previousKey < key))
            || face.axis != axis
            || face.minusRegionStableId == 0
            || face.plusRegionStableId == 0
            || (settings_.correspondenceMode
                    == PlanarFaceCorrespondenceMode::FixedMaterial
                && face.minusRegionStableId == face.plusRegionStableId)
            || face.minusRegionStableId != minusRegion
            || face.plusRegionStableId != plusRegion
            || !finite(face.lowerCornerMeters)
            || !finite(face.upperCornerMeters)
            || !std::isfinite(face.areaSquareMeters)
            || !(face.areaSquareMeters > 0.0)) {
            throw std::invalid_argument(
                "reference fluid faces are invalid, nonplanar, or not canonical");
        }
        havePreviousKey = true;
        previousKey = key;
        const StructureVector3 lower = toStructure(face.lowerCornerMeters);
        const StructureVector3 upper = toStructure(face.upperCornerMeters);
        if (std::abs(normalCoordinate(lower, axis) - planeCoordinate)
                > settings_.geometryToleranceMeters
            || std::abs(normalCoordinate(upper, axis) - planeCoordinate)
                > settings_.geometryToleranceMeters) {
            throw std::invalid_argument(
                "reference fluid faces do not share one plane");
        }
        const Vector2 lower2 = project(lower, axis);
        const Vector2 upper2 = project(upper, axis);
        const double rectangleArea =
            (upper2.u - lower2.u) * (upper2.v - lower2.v);
        const double areaTolerance = combinedTolerance(
            settings_.absoluteAreaToleranceSquareMeters,
            settings_.relativeAreaTolerance,
            rectangleArea, face.areaSquareMeters);
        if (!(upper2.u > lower2.u) || !(upper2.v > lower2.v)
            || std::abs(rectangleArea - face.areaSquareMeters)
                > areaTolerance) {
            throw std::invalid_argument(
                "reference fluid face bounds and area disagree");
        }
        faces_.push_back({
            face.minusRegionStableId,
            face.plusRegionStableId,
            face.axis,
            face.i,
            face.j,
            face.k,
            face.lowerCornerMeters,
            face.upperCornerMeters,
            face.areaSquareMeters,
        });
    }

    for (std::size_t first = 0; first < faces_.size(); ++first) {
        const Vector2 firstLower = project(
            toStructure(faces_[first].lowerCornerMeters), axis);
        const Vector2 firstUpper = project(
            toStructure(faces_[first].upperCornerMeters), axis);
        for (std::size_t second = first + 1;
             second < faces_.size(); ++second) {
            const Vector2 secondLower = project(
                toStructure(faces_[second].lowerCornerMeters), axis);
            const Vector2 secondUpper = project(
                toStructure(faces_[second].upperCornerMeters), axis);
            const double width = std::min(firstUpper.u, secondUpper.u)
                - std::max(firstLower.u, secondLower.u);
            const double height = std::min(firstUpper.v, secondUpper.v)
                - std::max(firstLower.v, secondLower.v);
            if (width > 0.0 && height > 0.0
                && width * height
                    >= settings_.minimumOverlapAreaSquareMeters) {
                throw std::invalid_argument(
                    "reference fluid face tiles overlap");
            }
        }
    }

    std::map<std::uint64_t, StructureVector3> referenceNodePositions;
    referenceNodePositions_.reserve(transfer_.nodes().size());
    for (const auto& node : transfer_.nodes()) {
        const StructureVector3 position =
            target.definition().nodes[node.structureNode].positionMeters;
        referenceNodePositions.emplace(node.stableId, position);
        referenceNodePositions_.push_back(position);
    }
    std::vector<double> faceCoverage(faces_.size(), 0.0);
    std::vector<double> triangleCoverage(
        transfer_.triangles().size(), 0.0);
    std::vector<std::vector<std::vector<Vector2>>> facePolygons(faces_.size());
    std::uint64_t nextOverlapStableId = 1;
    const StructureVector3 normal = axisUnit(axis);
    for (std::size_t triangleIndex = 0;
         triangleIndex < transfer_.triangles().size(); ++triangleIndex) {
        const auto& triangle = transfer_.triangles()[triangleIndex];
        std::array<StructureVector3, 3> positions;
        std::array<Vector2, 3> projected;
        for (std::size_t corner = 0; corner < 3; ++corner) {
            positions[corner] = referenceNodePositions.at(
                triangle.nodeStableIds[corner]);
            projected[corner] = project(positions[corner], axis);
            if (std::abs(normalCoordinate(positions[corner], axis)
                         - planeCoordinate)
                    > settings_.geometryToleranceMeters) {
                throw std::invalid_argument(
                    "reference structural triangle is not on the fluid plane");
            }
        }
        const StructureVector3 areaVector = cross(
            subtract(positions[1], positions[0]),
            subtract(positions[2], positions[0]));
        const double twiceArea = length(areaVector);
        if (!(twiceArea > 0.0)
            || dot(scale(areaVector, 1.0 / twiceArea), normal)
                < settings_.minimumNormalAlignment) {
            throw std::invalid_argument(
                "reference structural triangle orientation does not match the fluid faces");
        }
        const double triangleArea = 0.5 * twiceArea;
        referenceAreaSquareMeters_ += triangleArea;
        const std::vector<Vector2> trianglePolygon{
            projected[0], projected[1], projected[2]};
        for (std::size_t faceIndex = 0;
             faceIndex < faces_.size(); ++faceIndex) {
            const Vector2 lower = project(
                toStructure(faces_[faceIndex].lowerCornerMeters), axis);
            const Vector2 upper = project(
                toStructure(faces_[faceIndex].upperCornerMeters), axis);
            std::vector<Vector2> overlap = clipToRectangle(
                trianglePolygon, lower, upper);
            const auto [overlapArea, centroid] = areaAndCentroid(overlap);
            if (overlapArea
                < settings_.minimumOverlapAreaSquareMeters) {
                continue;
            }
            for (const auto& existing : facePolygons[faceIndex]) {
                const auto duplicate = clipToConvexPolygon(
                    overlap, existing, 0.0);
                if (areaAndCentroid(duplicate).first
                    >= settings_.minimumOverlapAreaSquareMeters) {
                    throw std::invalid_argument(
                        "reference structural triangles overlap on a fluid face");
                }
            }
            facePolygons[faceIndex].push_back(overlap);
            const auto coordinates = barycentric(centroid, projected);
            overlaps_.push_back({
                nextOverlapStableId++,
                triangle.stableId,
                faceIndex,
                coordinates,
                overlapArea,
            });
            faceCoverage[faceIndex] += overlapArea;
            triangleCoverage[triangleIndex] += overlapArea;
        }
        const double triangleTolerance = combinedTolerance(
            settings_.absoluteAreaToleranceSquareMeters,
            settings_.relativeAreaTolerance,
            triangleCoverage[triangleIndex], triangleArea);
        if (std::abs(triangleCoverage[triangleIndex] - triangleArea)
            > triangleTolerance) {
            throw std::invalid_argument(
                "fluid face tiles do not completely cover a structural triangle");
        }
    }
    for (std::size_t faceIndex = 0;
         faceIndex < faces_.size(); ++faceIndex) {
        const double tolerance = combinedTolerance(
            settings_.absoluteAreaToleranceSquareMeters,
            settings_.relativeAreaTolerance,
            faceCoverage[faceIndex], faces_[faceIndex].areaSquareMeters);
        if (std::abs(faceCoverage[faceIndex]
                     - faces_[faceIndex].areaSquareMeters) > tolerance) {
            throw std::invalid_argument(
                "structural triangles do not completely cover a fluid face tile");
        }
    }
    if (overlaps_.empty()) {
        throw std::invalid_argument(
            "face-resolved bridge has no positive-area overlap patches");
    }
}

std::uint64_t
PlanarFaceResolvedFluidStructureBridge::fluidSurfaceStableId() const noexcept {
    return fluidSurfaceStableId_;
}

const ConservativeSurfaceTransfer&
PlanarFaceResolvedFluidStructureBridge::transfer() const noexcept {
    return transfer_;
}

std::size_t
PlanarFaceResolvedFluidStructureBridge::overlapPatchCount() const noexcept {
    return overlaps_.size();
}

double PlanarFaceResolvedFluidStructureBridge::referenceAreaSquareMeters()
    const noexcept {
    return referenceAreaSquareMeters_;
}

PlanarFaceResolvedTransferResult
PlanarFaceResolvedFluidStructureBridge::evaluate(
    const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
    const std::span<const CouplingNodeKinematics> nodeKinematics) const {
    return evaluateImpl(fluidDiagnostics, nodeKinematics, std::nullopt);
}

PlanarFaceResolvedTransferResult
PlanarFaceResolvedFluidStructureBridge::evaluateMovingPlane(
    const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
    const std::span<const CouplingNodeKinematics> nodeKinematics,
    const double physicalPlaneCoordinateMeters) const {
    return evaluateImpl(
        fluidDiagnostics, nodeKinematics, physicalPlaneCoordinateMeters);
}

PlanarFaceResolvedTransferResult
PlanarFaceResolvedFluidStructureBridge::evaluateCutSurface(
    const fluid::PlanarCutSurfacePressureDiagnostics& cutSurface,
    const std::span<const CouplingNodeKinematics> nodeKinematics) const {
    if (settings_.correspondenceMode
            != PlanarFaceCorrespondenceMode::RigidNormalTranslation
        || cutSurface.version
            != fluid::planarCutSurfacePressureVersion
        || cutSurface.sourceInterfaceVersion
            != fluid::faceAlignedMovingInterfaceVersion
        || !cutSurface.finite || !cutSurface.accepted
        || cutSurface.surfaceStableId != fluidSurfaceStableId_
        || cutSurface.axis != axis_
        || cutSurface.faceCount != cutSurface.faces.size()
        || cutSurface.faceCount != faces_.size()
        || !finite(cutSurface.momentReferenceMeters)
        || !std::isfinite(cutSurface.gridPlaneCoordinateMeters)
        || !std::isfinite(cutSurface.physicalPlaneCoordinateMeters)
        || !std::isfinite(cutSurface.periodicPositionResidualMeters)
        || cutSurface.periodicPositionResidualMeters
            > settings_.geometryToleranceMeters
        || !std::isfinite(cutSurface.normalVelocityMetersPerSecond)
        || !std::isfinite(
            cutSurface.maximumNormalVelocitySpreadMetersPerSecond)
        || !std::isfinite(
            cutSurface.reactionSourcePhysicalPlaneCoordinateMeters)
        || !std::isfinite(
            cutSurface.reactionSourceNormalVelocityMetersPerSecond)
        || !std::isfinite(cutSurface.areaSquareMeters)
        || !std::isfinite(cutSurface.sourceAreaSquareMeters)
        || !finite(cutSurface.pressureForceNewtons)
        || !finite(cutSurface.sourcePressureForceNewtons)
        || !finite(cutSurface.pressureMomentNewtonMeters)
        || !std::isfinite(cutSurface.pressurePowerWatts)
        || !std::isfinite(cutSurface.sourcePressurePowerWatts)
        || length(subtract(
               toStructure(cutSurface.momentReferenceMeters),
               settings_.transfer.momentReferenceMeters))
            > settings_.geometryToleranceMeters) {
        throw std::invalid_argument(
            "face-resolved bridge requires a matching accepted cut surface");
    }

    fluid::MovingInterfaceProjectionDiagnostics source;
    source.projection.converged = true;
    source.interfaceVersion = cutSurface.sourceInterfaceVersion;
    source.interfaceFaceCount = cutSurface.faceCount;
    source.fluidRegionCount = 1;
    source.finite = true;
    source.faces.reserve(cutSurface.faces.size());
    for (const auto& face : cutSurface.faces) {
        StructureVector3 expectedLower = toStructure(
            face.gridLowerCornerMeters);
        StructureVector3 expectedUpper = toStructure(
            face.gridUpperCornerMeters);
        setNormalCoordinate(
            expectedLower, axis_, cutSurface.physicalPlaneCoordinateMeters);
        setNormalCoordinate(
            expectedUpper, axis_, cutSurface.physicalPlaneCoordinateMeters);
        if (face.surfaceStableId != fluidSurfaceStableId_
            || face.axis != axis_
            || !finite(face.gridLowerCornerMeters)
            || !finite(face.gridUpperCornerMeters)
            || !finite(face.physicalLowerCornerMeters)
            || !finite(face.physicalUpperCornerMeters)
            || length(subtract(
                   toStructure(face.physicalLowerCornerMeters),
                   expectedLower)) > settings_.geometryToleranceMeters
            || length(subtract(
                   toStructure(face.physicalUpperCornerMeters),
                   expectedUpper)) > settings_.geometryToleranceMeters) {
            throw std::invalid_argument(
                "cut-surface pressure geometry changed before transfer");
        }
        source.faces.push_back({
            face.surfaceStableId,
            face.minusRegionStableId,
            face.plusRegionStableId,
            face.axis,
            face.i,
            face.j,
            face.k,
            face.gridLowerCornerMeters,
            face.gridUpperCornerMeters,
            face.areaSquareMeters,
            face.normalVelocityMetersPerSecond,
            face.pressureTractionPascals,
            face.pressureForceNewtons,
            face.pressurePowerWatts,
        });
    }
    fluid::MovingInterfaceSurfaceDiagnostics surface;
    surface.stableId = cutSurface.surfaceStableId;
    surface.faceCount = cutSurface.faceCount;
    surface.areaSquareMeters = cutSurface.sourceAreaSquareMeters;
    surface.pressureForceNewtons = cutSurface.sourcePressureForceNewtons;
    surface.pressurePowerWatts = cutSurface.sourcePressurePowerWatts;
    source.surfaces.push_back(surface);

    auto result = evaluateImpl(
        source, nodeKinematics,
        cutSurface.physicalPlaneCoordinateMeters);
    const auto& diagnostics = result.diagnostics();
    const double forceTolerance = combinedTolerance(
        settings_.absoluteForceToleranceNewtons,
        settings_.relativeForceTolerance,
        length(diagnostics.fluidPressureForceNewtons),
        length(toStructure(cutSurface.pressureForceNewtons)));
    const double momentTolerance = combinedTolerance(
        settings_.absoluteMomentToleranceNewtonMeters,
        settings_.relativeMomentTolerance,
        length(diagnostics.fluidPressureMomentNewtonMeters),
        length(toStructure(cutSurface.pressureMomentNewtonMeters)));
    const double powerTolerance = combinedTolerance(
        settings_.absolutePowerToleranceWatts,
        settings_.relativePowerTolerance,
        std::abs(diagnostics.fluidPressurePowerWatts),
        std::abs(cutSurface.pressurePowerWatts));
    if (std::abs(diagnostics.gridPlaneCoordinateMeters
                 - cutSurface.gridPlaneCoordinateMeters)
            > settings_.geometryToleranceMeters
        || length(subtract(
               diagnostics.fluidPressureForceNewtons,
               toStructure(cutSurface.pressureForceNewtons)))
            > forceTolerance
        || length(subtract(
               diagnostics.fluidPressureMomentNewtonMeters,
               toStructure(cutSurface.pressureMomentNewtonMeters)))
            > momentTolerance
        || std::abs(diagnostics.fluidPressurePowerWatts
                    - cutSurface.pressurePowerWatts)
            > powerTolerance) {
        throw std::invalid_argument(
            "cut-surface and structural transfer ledgers do not close");
    }
    return result;
}

PlanarFaceResolvedTransferResult
PlanarFaceResolvedFluidStructureBridge::evaluateImpl(
    const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
    const std::span<const CouplingNodeKinematics> nodeKinematics,
    const std::optional<double> physicalPlaneCoordinateMeters) const {
    const bool movingCorrespondence =
        physicalPlaneCoordinateMeters.has_value();
    if (movingCorrespondence
            != (settings_.correspondenceMode
                == PlanarFaceCorrespondenceMode::RigidNormalTranslation)
        || (movingCorrespondence
            && !std::isfinite(*physicalPlaneCoordinateMeters))) {
        throw std::invalid_argument(
            "face-resolved bridge correspondence mode does not match evaluation");
    }
    if (fluidDiagnostics.interfaceVersion
            != fluid::faceAlignedMovingInterfaceVersion
        || !fluidDiagnostics.projection.converged
        || !fluidDiagnostics.finite
        || fluidDiagnostics.interfaceFaceCount
            != fluidDiagnostics.faces.size()) {
        throw std::invalid_argument(
            "face-resolved bridge requires accepted finite fluid diagnostics");
    }
    if (nodeKinematics.size() != transfer_.nodes().size()) {
        throw std::invalid_argument(
            "face-resolved bridge kinematics do not match the coupling surface");
    }
    double maximumRigidPositionResidualMeters = 0.0;
    for (std::size_t index = 0; index < nodeKinematics.size(); ++index) {
        if (nodeKinematics[index].stableId
                != transfer_.nodes()[index].stableId
            || !finite(nodeKinematics[index].positionMeters)
            || !finite(nodeKinematics[index].velocityMetersPerSecond)) {
            throw std::invalid_argument(
                "face-resolved bridge kinematics are non-finite or not canonical");
        }
        if (movingCorrespondence) {
            StructureVector3 expected = referenceNodePositions_[index];
            setNormalCoordinate(
                expected, axis_, *physicalPlaneCoordinateMeters);
            maximumRigidPositionResidualMeters = std::max(
                maximumRigidPositionResidualMeters,
                length(subtract(
                    nodeKinematics[index].positionMeters, expected)));
        }
    }
    if (maximumRigidPositionResidualMeters
        > settings_.geometryToleranceMeters) {
        throw std::invalid_argument(
            "moving face-resolved bridge requires rigid normal translation");
    }
    std::vector<const fluid::MovingInterfaceFaceDiagnostics*> selectedFaces;
    for (const auto& face : fluidDiagnostics.faces) {
        if (face.surfaceStableId == fluidSurfaceStableId_) {
            selectedFaces.push_back(&face);
        }
    }
    if (selectedFaces.size() != faces_.size()) {
        throw std::invalid_argument(
            "current fluid diagnostics do not match the bound face count");
    }
    const double gridPlaneCoordinateMeters = normalCoordinate(
        toStructure(selectedFaces.front()->lowerCornerMeters), axis_);
    const double currentPhysicalPlaneCoordinateMeters =
        movingCorrespondence
        ? *physicalPlaneCoordinateMeters : gridPlaneCoordinateMeters;

    PlanarFaceResolvedBridgeDiagnostics diagnostics;
    diagnostics.fluidSurfaceStableId = fluidSurfaceStableId_;
    diagnostics.fluidFaceCount = faces_.size();
    diagnostics.structureTriangleCount = transfer_.triangles().size();
    diagnostics.overlapPatchCount = overlaps_.size();
    diagnostics.referenceStructureAreaSquareMeters =
        referenceAreaSquareMeters_;
    diagnostics.correspondenceMode = settings_.correspondenceMode;
    diagnostics.gridPlaneCoordinateMeters = gridPlaneCoordinateMeters;
    diagnostics.physicalPlaneCoordinateMeters =
        currentPhysicalPlaneCoordinateMeters;
    diagnostics.normalTranslationFromReferenceMeters =
        currentPhysicalPlaneCoordinateMeters
        - referencePlaneCoordinateMeters_;
    diagnostics.maximumRigidPositionResidualMeters =
        maximumRigidPositionResidualMeters;
    std::vector<double> mappedFacePower(faces_.size(), 0.0);
    std::vector<CouplingTriangleTractionQuadrature> quadrature;
    quadrature.reserve(overlaps_.size());
    for (std::size_t faceIndex = 0;
         faceIndex < faces_.size(); ++faceIndex) {
        const auto& expected = faces_[faceIndex];
        const auto& actual = *selectedFaces[faceIndex];
        const StructureVector3 actualLower =
            toStructure(actual.lowerCornerMeters);
        const StructureVector3 actualUpper =
            toStructure(actual.upperCornerMeters);
        const Vector2 actualLower2 = project(actualLower, axis_);
        const Vector2 actualUpper2 = project(actualUpper, axis_);
        const Vector2 expectedLower2 = project(
            toStructure(expected.lowerCornerMeters), axis_);
        const Vector2 expectedUpper2 = project(
            toStructure(expected.upperCornerMeters), axis_);
        const bool commonGeometryMatches =
            actual.minusRegionStableId == expected.minusRegionStableId
            && actual.plusRegionStableId == expected.plusRegionStableId
            && actual.axis == axis_;
        const bool fixedGeometryMatches =
            actual.i == expected.i && actual.j == expected.j
            && actual.k == expected.k
            && std::abs(actual.lowerCornerMeters.x
                        - expected.lowerCornerMeters.x)
                <= settings_.geometryToleranceMeters
            && std::abs(actual.lowerCornerMeters.y
                        - expected.lowerCornerMeters.y)
                <= settings_.geometryToleranceMeters
            && std::abs(actual.lowerCornerMeters.z
                        - expected.lowerCornerMeters.z)
                <= settings_.geometryToleranceMeters
            && std::abs(actual.upperCornerMeters.x
                        - expected.upperCornerMeters.x)
                <= settings_.geometryToleranceMeters
            && std::abs(actual.upperCornerMeters.y
                        - expected.upperCornerMeters.y)
                <= settings_.geometryToleranceMeters
            && std::abs(actual.upperCornerMeters.z
                        - expected.upperCornerMeters.z)
                <= settings_.geometryToleranceMeters;
        const bool movingGeometryMatches =
            transverseIndicesMatch(
                actual, axis_, expected.i, expected.j, expected.k)
            && std::abs(actualLower2.u - expectedLower2.u)
                <= settings_.geometryToleranceMeters
            && std::abs(actualLower2.v - expectedLower2.v)
                <= settings_.geometryToleranceMeters
            && std::abs(actualUpper2.u - expectedUpper2.u)
                <= settings_.geometryToleranceMeters
            && std::abs(actualUpper2.v - expectedUpper2.v)
                <= settings_.geometryToleranceMeters
            && std::abs(normalCoordinate(actualLower, axis_)
                        - gridPlaneCoordinateMeters)
                <= settings_.geometryToleranceMeters
            && std::abs(normalCoordinate(actualUpper, axis_)
                        - gridPlaneCoordinateMeters)
                <= settings_.geometryToleranceMeters;
        const bool geometryMatches = commonGeometryMatches
            && (movingCorrespondence
                ? movingGeometryMatches : fixedGeometryMatches);
        const double areaTolerance = combinedTolerance(
            settings_.absoluteAreaToleranceSquareMeters,
            settings_.relativeAreaTolerance,
            actual.areaSquareMeters, expected.areaSquareMeters);
        if (!geometryMatches
            || !std::isfinite(actual.areaSquareMeters)
            || std::abs(actual.areaSquareMeters
                        - expected.areaSquareMeters) > areaTolerance
            || !std::isfinite(actual.normalVelocityMetersPerSecond)
            || !finite(actual.pressureTractionPascals)
            || !finite(actual.pressureForceNewtons)
            || !std::isfinite(actual.pressurePowerWatts)) {
            throw std::invalid_argument(
                "current fluid face geometry or pressure ledger changed incompatibly");
        }
        const StructureVector3 traction =
            toStructure(actual.pressureTractionPascals);
        const StructureVector3 reconstructedForce = scale(
            traction, actual.areaSquareMeters);
        const StructureVector3 sourceForce =
            toStructure(actual.pressureForceNewtons);
        const double forceTolerance = combinedTolerance(
            settings_.absoluteForceToleranceNewtons,
            settings_.relativeForceTolerance,
            length(reconstructedForce), length(sourceForce));
        const double reconstructedPower =
            dot(sourceForce, axisUnit(actual.axis))
            * actual.normalVelocityMetersPerSecond;
        const double powerTolerance = combinedTolerance(
            settings_.absolutePowerToleranceWatts,
            settings_.relativePowerTolerance,
            std::abs(reconstructedPower),
            std::abs(actual.pressurePowerWatts));
        if (length(subtract(reconstructedForce, sourceForce))
                > forceTolerance
            || std::abs(reconstructedPower - actual.pressurePowerWatts)
                > powerTolerance) {
            throw std::invalid_argument(
                "current fluid face traction, force, and power are inconsistent");
        }
        diagnostics.fluidAreaSquareMeters += actual.areaSquareMeters;
        diagnostics.fluidPressureForceNewtons = add(
            diagnostics.fluidPressureForceNewtons, sourceForce);
        diagnostics.fluidPressureMomentNewtonMeters = add(
            diagnostics.fluidPressureMomentNewtonMeters,
            cross(subtract(
                      physicalFaceCenter(
                          actual, axis_,
                          physicalPlaneCoordinateMeters),
                           settings_.transfer.momentReferenceMeters),
                  sourceForce));
        diagnostics.fluidPressurePowerWatts += actual.pressurePowerWatts;
    }

    if (movingCorrespondence) {
        const double expectedNormalVelocity =
            selectedFaces.front()->normalVelocityMetersPerSecond;
        for (const auto* face : selectedFaces) {
            diagnostics.maximumRigidVelocityResidualMetersPerSecond =
                std::max(
                    diagnostics.maximumRigidVelocityResidualMetersPerSecond,
                    std::abs(face->normalVelocityMetersPerSecond
                             - expectedNormalVelocity));
        }
        const StructureVector3 expectedVelocity = scale(
            axisUnit(axis_), expectedNormalVelocity);
        for (const auto& node : nodeKinematics) {
            diagnostics.maximumRigidVelocityResidualMetersPerSecond =
                std::max(
                    diagnostics.maximumRigidVelocityResidualMetersPerSecond,
                    length(subtract(
                        node.velocityMetersPerSecond,
                        expectedVelocity)));
        }
        const double velocityTolerance = combinedTolerance(
            settings_.absoluteVelocityToleranceMetersPerSecond,
            settings_.relativeVelocityTolerance,
            std::abs(expectedNormalVelocity),
            std::abs(expectedNormalVelocity)
                + diagnostics.maximumRigidVelocityResidualMetersPerSecond);
        if (diagnostics.maximumRigidVelocityResidualMetersPerSecond
            > velocityTolerance) {
            throw std::invalid_argument(
                "moving face-resolved bridge velocity is not rigid and normal");
        }
    }

    const auto aggregateSurface = std::lower_bound(
        fluidDiagnostics.surfaces.begin(), fluidDiagnostics.surfaces.end(),
        fluidSurfaceStableId_,
        [](const fluid::MovingInterfaceSurfaceDiagnostics& candidate,
           const std::uint64_t stableId) {
            return candidate.stableId < stableId;
        });
    if (aggregateSurface == fluidDiagnostics.surfaces.end()
        || aggregateSurface->stableId != fluidSurfaceStableId_
        || aggregateSurface->faceCount != faces_.size()
        || !std::isfinite(aggregateSurface->areaSquareMeters)
        || !(aggregateSurface->areaSquareMeters > 0.0)
        || !finite(aggregateSurface->pressureForceNewtons)
        || !std::isfinite(aggregateSurface->pressurePowerWatts)) {
        throw std::invalid_argument(
            "fluid surface aggregate is absent or disagrees with its faces");
    }
    const double aggregateAreaTolerance = combinedTolerance(
        settings_.absoluteAreaToleranceSquareMeters,
        settings_.relativeAreaTolerance,
        diagnostics.fluidAreaSquareMeters,
        aggregateSurface->areaSquareMeters);
    const double aggregateForceTolerance = combinedTolerance(
        settings_.absoluteForceToleranceNewtons,
        settings_.relativeForceTolerance,
        length(diagnostics.fluidPressureForceNewtons),
        length(toStructure(aggregateSurface->pressureForceNewtons)));
    const double aggregatePowerTolerance = combinedTolerance(
        settings_.absolutePowerToleranceWatts,
        settings_.relativePowerTolerance,
        std::abs(diagnostics.fluidPressurePowerWatts),
        std::abs(aggregateSurface->pressurePowerWatts));
    if (std::abs(diagnostics.fluidAreaSquareMeters
                 - aggregateSurface->areaSquareMeters)
            > aggregateAreaTolerance
        || length(subtract(
               diagnostics.fluidPressureForceNewtons,
               toStructure(aggregateSurface->pressureForceNewtons)))
            > aggregateForceTolerance
        || std::abs(diagnostics.fluidPressurePowerWatts
                    - aggregateSurface->pressurePowerWatts)
            > aggregatePowerTolerance) {
        throw std::invalid_argument(
            "fluid surface aggregate area, force, or power disagrees with its faces");
    }

    for (const auto& overlap : overlaps_) {
        const auto& face = *selectedFaces[overlap.faceIndex];
        quadrature.push_back({
            overlap.stableId,
            overlap.triangleStableId,
            overlap.barycentricCoordinates,
            overlap.areaSquareMeters,
            toStructure(face.pressureTractionPascals),
        });
        const auto triangle = std::lower_bound(
            transfer_.triangles().begin(), transfer_.triangles().end(),
            overlap.triangleStableId,
            [](const CouplingSurfaceTriangleDefinition& candidate,
               const std::uint64_t stableId) {
                return candidate.stableId < stableId;
            });
        const std::size_t triangleIndex = static_cast<std::size_t>(
            triangle - transfer_.triangles().begin());
        const auto& triangleDefinition = transfer_.triangles()[triangleIndex];
        StructureVector3 velocity;
        for (std::size_t corner = 0; corner < 3; ++corner) {
            const auto node = std::lower_bound(
                transfer_.nodes().begin(), transfer_.nodes().end(),
                triangleDefinition.nodeStableIds[corner],
                [](const CouplingSurfaceNodeDefinition& candidate,
                   const std::uint64_t stableId) {
                    return candidate.stableId < stableId;
                });
            const std::size_t nodeIndex = static_cast<std::size_t>(
                node - transfer_.nodes().begin());
            velocity = add(velocity, scale(
                nodeKinematics[nodeIndex].velocityMetersPerSecond,
                overlap.barycentricCoordinates[corner]));
        }
        mappedFacePower[overlap.faceIndex] += dot(
            scale(toStructure(face.pressureTractionPascals),
                  overlap.areaSquareMeters),
            velocity);
    }
    for (std::size_t faceIndex = 0;
         faceIndex < faces_.size(); ++faceIndex) {
        const double residual = mappedFacePower[faceIndex]
            - selectedFaces[faceIndex]->pressurePowerWatts;
        diagnostics.maximumFacePowerResidualWatts = std::max(
            diagnostics.maximumFacePowerResidualWatts, std::abs(residual));
        const double tolerance = combinedTolerance(
            settings_.absolutePowerToleranceWatts,
            settings_.relativePowerTolerance,
            std::abs(mappedFacePower[faceIndex]),
            std::abs(selectedFaces[faceIndex]->pressurePowerWatts));
        if (std::abs(residual) > tolerance) {
            throw std::invalid_argument(
                "mapped structural velocity does not match one fluid face power ledger");
        }
    }

    ConservativeTransferResult transferred = transfer_.evaluateQuadrature(
        nodeKinematics, quadrature, settings_.transfer);
    const auto& target = transferred.diagnostics();
    diagnostics.structureSurfaceForceNewtons =
        target.integratedSurfaceForceNewtons;
    diagnostics.forceResidualNewtons = subtract(
        target.integratedSurfaceForceNewtons,
        diagnostics.fluidPressureForceNewtons);
    diagnostics.forceResidualNormNewtons = length(
        diagnostics.forceResidualNewtons);
    diagnostics.structureSurfaceMomentNewtonMeters =
        target.integratedSurfaceMomentNewtonMeters;
    diagnostics.momentResidualNewtonMeters = subtract(
        target.integratedSurfaceMomentNewtonMeters,
        diagnostics.fluidPressureMomentNewtonMeters);
    diagnostics.momentResidualNormNewtonMeters = length(
        diagnostics.momentResidualNewtonMeters);
    diagnostics.structureSurfacePowerWatts =
        target.integratedSurfacePowerWatts;
    diagnostics.powerResidualWatts =
        target.integratedSurfacePowerWatts
        - diagnostics.fluidPressurePowerWatts;
    diagnostics.areaResidualSquareMeters =
        target.surfaceAreaSquareMeters - diagnostics.fluidAreaSquareMeters;
    diagnostics.finite =
        std::isfinite(diagnostics.fluidAreaSquareMeters)
        && std::isfinite(diagnostics.areaResidualSquareMeters)
        && finite(diagnostics.fluidPressureForceNewtons)
        && finite(diagnostics.structureSurfaceForceNewtons)
        && finite(diagnostics.forceResidualNewtons)
        && std::isfinite(diagnostics.forceResidualNormNewtons)
        && finite(diagnostics.fluidPressureMomentNewtonMeters)
        && finite(diagnostics.structureSurfaceMomentNewtonMeters)
        && finite(diagnostics.momentResidualNewtonMeters)
        && std::isfinite(diagnostics.momentResidualNormNewtonMeters)
        && std::isfinite(diagnostics.fluidPressurePowerWatts)
        && std::isfinite(diagnostics.structureSurfacePowerWatts)
        && std::isfinite(diagnostics.powerResidualWatts)
        && std::isfinite(diagnostics.maximumFacePowerResidualWatts)
        && std::isfinite(diagnostics.gridPlaneCoordinateMeters)
        && std::isfinite(diagnostics.physicalPlaneCoordinateMeters)
        && std::isfinite(
            diagnostics.normalTranslationFromReferenceMeters)
        && std::isfinite(
            diagnostics.maximumRigidPositionResidualMeters)
        && std::isfinite(
            diagnostics.maximumRigidVelocityResidualMetersPerSecond)
        && target.finite;
    if (!diagnostics.finite) {
        throw std::invalid_argument(
            "face-resolved bridge produced non-finite diagnostics");
    }

    const double areaTolerance = combinedTolerance(
        settings_.absoluteAreaToleranceSquareMeters,
        settings_.relativeAreaTolerance,
        diagnostics.fluidAreaSquareMeters,
        target.surfaceAreaSquareMeters);
    const double forceTolerance = combinedTolerance(
        settings_.absoluteForceToleranceNewtons,
        settings_.relativeForceTolerance,
        length(diagnostics.fluidPressureForceNewtons),
        length(target.integratedSurfaceForceNewtons));
    const double momentTolerance = combinedTolerance(
        settings_.absoluteMomentToleranceNewtonMeters,
        settings_.relativeMomentTolerance,
        length(diagnostics.fluidPressureMomentNewtonMeters),
        length(target.integratedSurfaceMomentNewtonMeters));
    const double powerTolerance = combinedTolerance(
        settings_.absolutePowerToleranceWatts,
        settings_.relativePowerTolerance,
        std::abs(diagnostics.fluidPressurePowerWatts),
        std::abs(target.integratedSurfacePowerWatts));
    if (std::abs(diagnostics.areaResidualSquareMeters) > areaTolerance
        || diagnostics.forceResidualNormNewtons > forceTolerance
        || diagnostics.momentResidualNormNewtonMeters > momentTolerance
        || std::abs(diagnostics.powerResidualWatts) > powerTolerance) {
        throw std::invalid_argument(
            "face-resolved fluid and structural area, force, moment, or power ledgers do not close");
    }
    return {std::move(transferred), std::move(diagnostics)};
}

} // namespace simwing::fsi

#include "hemisphere_case.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

constexpr double radiusMeters = 1.0;
constexpr double fabricArealDensityKgPerSquareMeter = 0.08;
constexpr double basePressurePascals = 55.0;
constexpr double pressureAmplitudePascals = 45.0;
constexpr double pressureFrequencyHertz = 0.4;

[[nodiscard]] StructureVector3 subtract(const StructureVector3& first,
                                        const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

[[nodiscard]] StructureVector3 cross(const StructureVector3& first,
                                     const StructureVector3& second) {
    return {first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
}

[[nodiscard]] double dot(const StructureVector3& first,
                         const StructureVector3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

[[nodiscard]] double length(const StructureVector3& value) {
    return std::sqrt(dot(value, value));
}

[[nodiscard]] StructureVector3 scaled(const StructureVector3& value,
                                      double scale) {
    return {scale * value.x, scale * value.y, scale * value.z};
}

[[nodiscard]] std::size_t ringNode(std::size_t latitude,
                                   std::size_t longitude) {
    return 1 + (latitude - 1) * clampedHemisphereRadialSegments
        + longitude % clampedHemisphereRadialSegments;
}

[[nodiscard]] StructureMembraneMaterial membraneMaterial() {
    StructureMembraneMaterial material;
    material.warpStiffnessNewtonsPerMeter = 9000.0;
    material.weftStiffnessNewtonsPerMeter = 9000.0;
    material.couplingStiffnessNewtonsPerMeter = 2400.0;
    material.shearStiffnessNewtonsPerMeter = 3300.0;
    material.warpPreTensionNewtonsPerMeter = 0.0;
    material.weftPreTensionNewtonsPerMeter = 0.0;
    material.dampingSeconds = 0.015;
    material.compressionStiffnessRatio = 0.2;
    return material;
}

[[nodiscard]] std::array<StructureVector2, 3> intrinsicChart(
    const StructureDefinition& definition,
    const StructureTriangleDefinition& triangle) {
    const StructureVector3& first =
        definition.nodes[triangle.nodes[0]].positionMeters;
    const StructureVector3& second =
        definition.nodes[triangle.nodes[1]].positionMeters;
    const StructureVector3& third =
        definition.nodes[triangle.nodes[2]].positionMeters;
    const StructureVector3 edge = subtract(second, first);
    const StructureVector3 diagonal = subtract(third, first);
    const double edgeLength = length(edge);
    const double thirdX = dot(diagonal, edge) / edgeLength;
    const double thirdYSquared = std::max(
        0.0, dot(diagonal, diagonal) - thirdX * thirdX);
    return {StructureVector2{0.0, 0.0},
            StructureVector2{edgeLength, 0.0},
            StructureVector2{thirdX, std::sqrt(thirdYSquared)}};
}

struct EdgeIncidence {
    std::size_t from = 0;
    std::size_t to = 0;
    std::size_t opposite = 0;
};

void addRestShapeDihedrals(StructureDefinition& definition) {
    using Edge = std::pair<std::size_t, std::size_t>;
    std::map<Edge, std::vector<EdgeIncidence>> incidences;
    for (const StructureTriangleDefinition& triangle : definition.triangles) {
        for (std::size_t edge = 0; edge < 3; ++edge) {
            const std::size_t from = triangle.nodes[edge];
            const std::size_t to = triangle.nodes[(edge + 1) % 3];
            const std::size_t opposite = triangle.nodes[(edge + 2) % 3];
            incidences[{std::min(from, to), std::max(from, to)}].push_back(
                {from, to, opposite});
        }
    }

    for (const auto& [edge, adjacent] : incidences) {
        static_cast<void>(edge);
        if (adjacent.size() == 1) {
            continue;
        }
        if (adjacent.size() != 2
            || adjacent[0].from != adjacent[1].to
            || adjacent[0].to != adjacent[1].from) {
            throw std::logic_error(
                "clamped hemisphere mesh is not an oriented manifold");
        }

        const EdgeIncidence& first = adjacent[0];
        const EdgeIncidence& second = adjacent[1];
        const StructureVector3& a =
            definition.nodes[first.from].positionMeters;
        const StructureVector3& b =
            definition.nodes[first.to].positionMeters;
        const StructureVector3& c =
            definition.nodes[first.opposite].positionMeters;
        const StructureVector3& d =
            definition.nodes[second.opposite].positionMeters;
        const StructureVector3 edgeVector = subtract(b, a);
        const double edgeLength = length(edgeVector);
        const StructureVector3 firstArea = cross(edgeVector, subtract(c, a));
        const StructureVector3 secondArea = cross(subtract(d, a), edgeVector);
        const StructureVector3 firstNormal =
            scaled(firstArea, 1.0 / length(firstArea));
        const StructureVector3 secondNormal =
            scaled(secondArea, 1.0 / length(secondArea));
        const StructureVector3 edgeUnit =
            scaled(edgeVector, 1.0 / edgeLength);
        const double cosine = std::clamp(
            dot(firstNormal, secondNormal), -1.0, 1.0);
        const double sine = dot(
            edgeUnit, cross(firstNormal, secondNormal));
        definition.dihedrals.push_back(
            {{first.from, first.to, first.opposite, second.opposite},
             std::atan2(sine, cosine),
             5.0});
    }
}

[[nodiscard]] StructureDefinition makeDefinition() {
    StructureDefinition definition;
    definition.nodes.reserve(
        1 + clampedHemisphereLatitudeSegments
                * clampedHemisphereRadialSegments);
    definition.nodes.push_back(
        {{0.0, 0.0, radiusMeters}, 0.0, false});
    for (std::size_t latitude = 1;
         latitude <= clampedHemisphereLatitudeSegments;
         ++latitude) {
        const double polarAngle =
            0.5 * std::numbers::pi * static_cast<double>(latitude)
            / static_cast<double>(clampedHemisphereLatitudeSegments);
        const double radialDistance = radiusMeters * std::sin(polarAngle);
        const double z = radiusMeters * std::cos(polarAngle);
        const bool fixed = latitude == clampedHemisphereLatitudeSegments;
        for (std::size_t longitude = 0;
             longitude < clampedHemisphereRadialSegments;
             ++longitude) {
            const double azimuth =
                2.0 * std::numbers::pi * static_cast<double>(longitude)
                / static_cast<double>(clampedHemisphereRadialSegments);
            definition.nodes.push_back(
                {{radialDistance * std::cos(azimuth),
                  radialDistance * std::sin(azimuth),
                  z},
                 0.0,
                 fixed});
        }
    }

    definition.triangles.reserve(
        clampedHemisphereRadialSegments
        * (2 * clampedHemisphereLatitudeSegments - 1));
    for (std::size_t longitude = 0;
         longitude < clampedHemisphereRadialSegments;
         ++longitude) {
        definition.triangles.push_back(
            {{0,
              ringNode(1, longitude),
              ringNode(1, longitude + 1)}});
    }
    for (std::size_t latitude = 1;
         latitude < clampedHemisphereLatitudeSegments;
         ++latitude) {
        for (std::size_t longitude = 0;
             longitude < clampedHemisphereRadialSegments;
             ++longitude) {
            const std::size_t upper = ringNode(latitude, longitude);
            const std::size_t upperNext = ringNode(latitude, longitude + 1);
            const std::size_t lower = ringNode(latitude + 1, longitude);
            const std::size_t lowerNext =
                ringNode(latitude + 1, longitude + 1);
            definition.triangles.push_back(
                {{upper, lower, lowerNext}});
            definition.triangles.push_back(
                {{upper, lowerNext, upperNext}});
        }
    }

    const StructureMembraneMaterial material = membraneMaterial();
    definition.membranes.reserve(definition.triangles.size());
    for (std::size_t triangleIndex = 0;
         triangleIndex < definition.triangles.size();
         ++triangleIndex) {
        const StructureTriangleDefinition& triangle =
            definition.triangles[triangleIndex];
        const StructureVector3& first =
            definition.nodes[triangle.nodes[0]].positionMeters;
        const StructureVector3& second =
            definition.nodes[triangle.nodes[1]].positionMeters;
        const StructureVector3& third =
            definition.nodes[triangle.nodes[2]].positionMeters;
        const double area = 0.5 * length(cross(
            subtract(second, first), subtract(third, first)));
        for (const std::size_t node : triangle.nodes) {
            definition.nodes[node].massKg +=
                fabricArealDensityKgPerSquareMeter * area / 3.0;
        }
        definition.membranes.push_back(
            {triangleIndex,
             intrinsicChart(definition, triangle),
             material,
             StructureMaterialRole::Bulk});
    }

    const std::size_t boundaryLatitude = clampedHemisphereLatitudeSegments;
    definition.constraints.reserve(clampedHemisphereRadialSegments);
    for (std::size_t longitude = 0;
         longitude < clampedHemisphereRadialSegments;
         ++longitude) {
        const std::size_t first = ringNode(boundaryLatitude, longitude);
        const std::size_t second = ringNode(boundaryLatitude, longitude + 1);
        definition.constraints.push_back(
            {StructureConstraintKind::Distance,
             first,
             second,
             length(subtract(definition.nodes[second].positionMeters,
                             definition.nodes[first].positionMeters)),
             0.0});
    }
    addRestShapeDihedrals(definition);
    return definition;
}

[[nodiscard]] viewer::StructureFrameMappingDefinition makeFrameMapping() {
    viewer::StructureFrameMappingDefinition mapping;
    mapping.vertexStableIds.reserve(
        1 + clampedHemisphereLatitudeSegments
                * clampedHemisphereRadialSegments);
    for (std::size_t index = 0;
         index < 1 + clampedHemisphereLatitudeSegments
                       * clampedHemisphereRadialSegments;
         ++index) {
        mapping.vertexStableIds.push_back(100'000 + index);
    }
    const std::size_t triangleCount = clampedHemisphereRadialSegments
        * (2 * clampedHemisphereLatitudeSegments - 1);
    mapping.triangles.reserve(triangleCount);
    for (std::size_t index = 0; index < triangleCount; ++index) {
        mapping.triangles.push_back({200'000 + index, 1, 2});
    }
    mapping.lines.reserve(clampedHemisphereRadialSegments);
    for (std::size_t index = 0;
         index < clampedHemisphereRadialSegments;
         ++index) {
        mapping.lines.push_back(
            {300'000 + index,
             static_cast<std::uint32_t>(
                 StructureConstraintKind::Distance)});
    }
    return mapping;
}

[[nodiscard]] StructureStepSettings makeStepSettings() {
    StructureStepSettings settings;
    settings.timeStepSeconds = 1.0 / 120.0;
    settings.substeps = 4;
    settings.constraintIterations = 20;
    settings.cableConstraintSweepPairs = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.7;
    settings.workerThreads = 0;
    return settings;
}

} // namespace

ClampedHemisphereCase::ClampedHemisphereCase()
    : structure_(makeDefinition()),
      frameMapping_(structure_, makeFrameMapping()),
      stepSettings_(makeStepSettings()) {}

viewer::TraceHeader ClampedHemisphereCase::traceHeader() const {
    return {clampedHemisphereCaseChecksum, clampedHemisphereCaseSolverId};
}

viewer::DiagnosticFrame ClampedHemisphereCase::advance() {
    const double loadTime = structure_.simulationTimeSeconds();
    pressurePascals_ = basePressurePascals
        + pressureAmplitudePascals
            * std::sin(2.0 * std::numbers::pi
                       * pressureFrequencyHertz * loadTime);

    const StructureDefinition& definition = structure_.definition();
    const std::vector<StructureNodeState> states = structure_.nodeStates();
    std::vector<StructureVector3> forces(states.size());
    for (const StructureTriangleDefinition& triangle : definition.triangles) {
        const StructureVector3& first =
            states[triangle.nodes[0]].positionMeters;
        const StructureVector3& second =
            states[triangle.nodes[1]].positionMeters;
        const StructureVector3& third =
            states[triangle.nodes[2]].positionMeters;
        const StructureVector3 nodalForce = scaled(
            cross(subtract(second, first), subtract(third, first)),
            pressurePascals_ / 6.0);
        for (const std::size_t node : triangle.nodes) {
            forces[node].x += nodalForce.x;
            forces[node].y += nodalForce.y;
            forces[node].z += nodalForce.z;
        }
    }
    structure_.setExternalForces(forces);
    const StructureDiagnostics diagnostics = structure_.step(stepSettings_);
    if (!diagnostics.finite) {
        throw std::runtime_error(
            "clamped hemisphere step produced non-finite diagnostics");
    }

    viewer::StructureFrameContext context;
    context.sceneChecksum = clampedHemisphereCaseChecksum;
    context.solverCommit = clampedHemisphereCaseSolverId;
    context.timeStepSeconds = stepSettings_.timeStepSeconds;
    context.couplingIteration = 0;
    context.couplingResiduals.structure =
        diagnostics.maximumMembraneResidual;
    viewer::DiagnosticFrame frame =
        viewer::buildStructureFrame(structure_, frameMapping_, context);

    std::vector<double> radialDisplacement;
    radialDisplacement.reserve(frame.vertices.size());
    for (const viewer::DiagnosticVertex& vertex : frame.vertices) {
        radialDisplacement.push_back(
            std::sqrt(vertex.positionMetres.x * vertex.positionMetres.x
                      + vertex.positionMetres.y * vertex.positionMetres.y
                      + vertex.positionMetres.z * vertex.positionMetres.z)
            - radiusMeters);
    }
    frame.scalarFields.push_back(
        {"dome.radial_displacement",
         "m",
         viewer::FieldAssociation::Vertex,
         std::move(radialDisplacement)});
    frame.scalarFields.push_back(
        {"dome.pressure",
         "Pa",
         viewer::FieldAssociation::Triangle,
         std::vector<double>(frame.triangles.size(), pressurePascals_)});
    std::vector<viewer::Vec3d> pressureForces;
    pressureForces.reserve(forces.size());
    for (const StructureVector3& force : forces) {
        pressureForces.push_back({force.x, force.y, force.z});
    }
    frame.vectorFields.push_back(
        {"dome.pressure_force",
         "N",
         viewer::FieldAssociation::Vertex,
         std::move(pressureForces)});
    viewer::ProtocolError error;
    if (!viewer::validateFrame(frame, &error)) {
        throw std::logic_error(
            "clamped hemisphere frame is invalid: " + error.message);
    }
    return frame;
}

const Structure& ClampedHemisphereCase::structure() const noexcept {
    return structure_;
}

const StructureStepSettings& ClampedHemisphereCase::stepSettings() const
    noexcept {
    return stepSettings_;
}

double ClampedHemisphereCase::pressurePascals() const noexcept {
    return pressurePascals_;
}

double ClampedHemisphereCase::apexRadialDisplacementMeters() const {
    return structure_.nodeStates().front().positionMeters.z - radiusMeters;
}

} // namespace simwing::fsi

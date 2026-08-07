#include "scene_structure.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <string_view>
#include <utility>

namespace simwing::fsi {
namespace {

template <typename T>
std::vector<const T*> sortedById(const std::vector<T>& values) {
    std::vector<const T*> sorted;
    sorted.reserve(values.size());
    for (const T& value : values) {
        sorted.push_back(&value);
    }
    std::ranges::sort(sorted, {}, [](const T* value) { return value->id; });
    return sorted;
}

void addDiagnostic(SceneStructureAssembly& assembly,
                   SceneStructureDiagnosticCode code,
                   EntityKind entityKind,
                   StableId entityId,
                   std::string message,
                   std::optional<ValidationCode> validationCode =
                       std::nullopt) {
    assembly.diagnostics.push_back(
        {SceneStructureDiagnosticSeverity::Error,
         code,
         entityKind,
         entityId,
         validationCode,
         std::move(message)});
}

void sortDiagnostics(SceneStructureAssembly& assembly) {
    std::ranges::sort(
        assembly.diagnostics,
        [](const SceneStructureDiagnostic& left,
           const SceneStructureDiagnostic& right) {
            if (left.entityKind != right.entityKind) {
                return left.entityKind < right.entityKind;
            }
            if (left.entityId != right.entityId) {
                return left.entityId < right.entityId;
            }
            if (left.code != right.code) {
                return left.code < right.code;
            }
            if (left.sceneValidationCode != right.sceneValidationCode) {
                return left.sceneValidationCode < right.sceneValidationCode;
            }
            return left.message < right.message;
        });
}

void reject(SceneStructureAssembly& assembly) {
    assembly.assembled = false;
    assembly.definition = {};
    assembly.mappings = {};
    assembly.totalFabricMassKg = 0.0;
    sortDiagnostics(assembly);
}

bool checkedMappingBytes(std::size_t count,
                         std::size_t elementSize,
                         std::size_t& remaining) {
    if (elementSize != 0 && count > remaining / elementSize) {
        return false;
    }
    remaining -= count * elementSize;
    return true;
}

bool checkBounds(const Scene& scene,
                 const SceneStructureLimits& limits,
                 SceneStructureAssembly& assembly) {
    const bool nodeCountOverflow = scene.suspensionJunctions.size()
        > std::numeric_limits<std::size_t>::max() - scene.vertices.size();
    const std::size_t nodeCount = nodeCountOverflow
        ? std::numeric_limits<std::size_t>::max()
        : scene.vertices.size() + scene.suspensionJunctions.size();
    if (nodeCountOverflow || nodeCount > limits.maximumNodes) {
        addDiagnostic(assembly,
                      SceneStructureDiagnosticCode::MappingOverflow,
                      EntityKind::Scene,
                      invalidStableId,
                      "scene vertex and junction count exceeds the structure node bound");
    }
    if (scene.triangles.size() > limits.maximumTriangles) {
        addDiagnostic(
            assembly,
            SceneStructureDiagnosticCode::MappingOverflow,
            EntityKind::Scene,
            invalidStableId,
            "scene triangle count exceeds the structure triangle bound");
    }
    if (scene.suspensionLines.size() > limits.maximumConstraints) {
        addDiagnostic(
            assembly,
            SceneStructureDiagnosticCode::MappingOverflow,
            EntityKind::Scene,
            invalidStableId,
            "scene suspension-line count exceeds the structure constraint bound");
    }

    std::size_t remaining = limits.maximumMappingBytes;
    const bool fits =
        checkedMappingBytes(scene.vertices.size(), sizeof(StableId), remaining)
        && checkedMappingBytes(scene.suspensionJunctions.size(),
                               sizeof(StableId), remaining)
        && checkedMappingBytes(scene.triangles.size(), sizeof(StableId),
                               remaining)
        && checkedMappingBytes(scene.triangles.size(), sizeof(StableId),
                               remaining)
        && checkedMappingBytes(scene.suspensionLines.size(), sizeof(StableId),
                               remaining)
        // A manifold triangulation has fewer than two hinges per triangle;
        // claiming two pairs is a simple overflow-safe conservative bound.
        && checkedMappingBytes(scene.triangles.size(),
                               2 * sizeof(std::array<StableId, 2>),
                               remaining);
    if (!fits) {
        addDiagnostic(
            assembly,
            SceneStructureDiagnosticCode::MappingOverflow,
            EntityKind::Scene,
            invalidStableId,
            "stable-ID mappings exceed the configured byte bound");
    }
    return assembly.diagnostics.empty();
}

StructureVector3 convert(const Vec3& value) {
    return {value.x, value.y, value.z};
}

StructureVector2 convert(const Vec2& value) {
    return {value.x, value.y};
}

Vec3 subtract(const Vec3& first, const Vec3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

Vec3 cross(const Vec3& first, const Vec3& second) {
    return {first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
}

double dot(const Vec3& first, const Vec3& second) {
    return first.x * second.x + first.y * second.y
        + first.z * second.z;
}

double length(const Vec3& value) {
    return std::sqrt(dot(value, value));
}

Vec3 scaled(const Vec3& value, double scale) {
    return {scale * value.x, scale * value.y, scale * value.z};
}

bool compatibleMaterial(const FabricMaterial& material) {
    // This is the diagonal (zero coupling) specialization of softwing's
    // checked SymmetricMatrix3 inverse. Mirroring its conditioning test here
    // prevents an assembly reported as successful from failing in Structure's
    // constructor.
    const double determinant = material.warpStiffnessNewtonsPerMeter
        * material.weftStiffnessNewtonsPerMeter
        * material.shearStiffnessNewtonsPerMeter;
    const double scale = std::max(
        {material.warpStiffnessNewtonsPerMeter,
         material.weftStiffnessNewtonsPerMeter,
         material.shearStiffnessNewtonsPerMeter});
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon()
        * scale * scale * scale;
    return scale > 0.0 && std::isfinite(determinant)
        && std::isfinite(tolerance) && determinant > tolerance
        && std::isfinite(1.0 / material.warpStiffnessNewtonsPerMeter)
        && std::isfinite(1.0 / material.weftStiffnessNewtonsPerMeter)
        && std::isfinite(1.0 / material.shearStiffnessNewtonsPerMeter);
}

bool compatibleMaterialCoordinates(const Triangle& triangle) {
    const Vec2 first{
        triangle.materialCoordinates[1].x
            - triangle.materialCoordinates[0].x,
        triangle.materialCoordinates[1].y
            - triangle.materialCoordinates[0].y};
    const Vec2 second{
        triangle.materialCoordinates[2].x
            - triangle.materialCoordinates[0].x,
        triangle.materialCoordinates[2].y
            - triangle.materialCoordinates[0].y};
    const double determinant = first.x * second.y - second.x * first.y;
    const Vec2 third{
        triangle.materialCoordinates[2].x
            - triangle.materialCoordinates[1].x,
        triangle.materialCoordinates[2].y
            - triangle.materialCoordinates[1].y};
    const double maximumEdgeLengthSquared = std::max(
        {first.x * first.x + first.y * first.y,
         second.x * second.x + second.y * second.y,
         third.x * third.x + third.y * third.y});
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
        * maximumEdgeLengthSquared;
    // softwing requires positive material-chart area, not merely a
    // nonsingular chart. Matching that contract here makes successful bridge
    // output safe to pass to Structure without a later constructor failure.
    return maximumEdgeLengthSquared > 0.0 && std::isfinite(determinant)
        && std::isfinite(tolerance) && determinant > tolerance;
}

double materialChartArea(const Triangle& triangle) {
    const Vec2 first{
        triangle.materialCoordinates[1].x
            - triangle.materialCoordinates[0].x,
        triangle.materialCoordinates[1].y
            - triangle.materialCoordinates[0].y};
    const Vec2 second{
        triangle.materialCoordinates[2].x
            - triangle.materialCoordinates[0].x,
        triangle.materialCoordinates[2].y
            - triangle.materialCoordinates[0].y};
    return 0.5 * (first.x * second.y - second.x * first.y);
}

struct EdgeIncidence {
    const Triangle* triangle = nullptr;
    StableId from = invalidStableId;
    StableId to = invalidStableId;
    StableId opposite = invalidStableId;
};

void assembleBending(
    const std::vector<const Triangle*>& triangles,
    const std::map<StableId, const Vertex*>& verticesById,
    const std::map<StableId, const FabricMaterial*>& materialsById,
    const std::map<StableId, std::size_t>& nodeIndices,
    SceneStructureAssembly& assembly) {
    // The same spatial edge may be shared by several welded physical sheets
    // (for example skin and rib). Sheet identity is therefore part of the
    // bending key: only adjacent triangles cut from one authored sheet form a
    // dihedral.
    std::map<std::array<StableId, 3>, std::vector<EdgeIncidence>> edges;
    for (const Triangle* triangle : triangles) {
        for (std::size_t corner = 0; corner < 3; ++corner) {
            const StableId from = triangle->vertexIds[corner];
            const StableId to = triangle->vertexIds[(corner + 1) % 3];
            const StableId opposite = triangle->vertexIds[(corner + 2) % 3];
            std::array<StableId, 2> canonicalEdge{from, to};
            std::ranges::sort(canonicalEdge);
            const std::array<StableId, 3> key{
                triangle->sheetId, canonicalEdge[0], canonicalEdge[1]};
            edges[key].push_back({triangle, from, to, opposite});
        }
    }

    assembly.definition.dihedrals.reserve(triangles.size() * 3 / 2);
    assembly.mappings.dihedralTriangleIds.reserve(
        triangles.size() * 3 / 2);
    for (const auto& [sheetEdge, incidences] : edges) {
        const StableId sheetId = sheetEdge[0];
        const std::array<StableId, 2> edge{sheetEdge[1], sheetEdge[2]};
        if (incidences.size() == 1) {
            // A boundary has no adjacent fabric strip and intentionally gets
            // no bending constraint.
            continue;
        }
        if (incidences.size() != 2) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedBendingTopology,
                EntityKind::Triangle,
                incidences.front().triangle->id,
                "non-manifold edge " + std::to_string(edge[0]) + "-"
                    + std::to_string(edge[1])
                    + " in sheet " + std::to_string(sheetId)
                    + " cannot define one fabric bending hinge");
            continue;
        }
        const EdgeIncidence& first = incidences[0];
        const EdgeIncidence& second = incidences[1];
        if (first.from != second.to || first.to != second.from) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedBendingTopology,
                EntityKind::Triangle,
                first.triangle->id,
                "inconsistently oriented edge " + std::to_string(edge[0])
                    + "-" + std::to_string(edge[1])
                    + " in sheet " + std::to_string(sheetId)
                    + " cannot define a signed fabric bending hinge");
            continue;
        }

        const double firstRigidity = materialsById.at(
            first.triangle->materialId)->bendingStiffnessNewtonMeters;
        const double secondRigidity = materialsById.at(
            second.triangle->materialId)->bendingStiffnessNewtonMeters;
        if (!(firstRigidity > 0.0) || !(secondRigidity > 0.0)) {
            // Series strips with a zero-rigidity side form a free hinge.
            continue;
        }

        const Vec3& a = verticesById.at(first.from)->positionMeters;
        const Vec3& b = verticesById.at(first.to)->positionMeters;
        const Vec3& c = verticesById.at(first.opposite)->positionMeters;
        const Vec3& d = verticesById.at(second.opposite)->positionMeters;
        const Vec3 edgeVector = subtract(b, a);
        const double edgeLength = length(edgeVector);
        const Vec3 firstAreaVector = cross(edgeVector, subtract(c, a));
        const Vec3 secondAreaVector = cross(subtract(d, a), edgeVector);
        const double firstAreaLength = length(firstAreaVector);
        const double secondAreaLength = length(secondAreaVector);
        const auto chartCoordinate = [](const Triangle& triangle,
                                        StableId vertexId) {
            for (std::size_t corner = 0; corner < 3; ++corner) {
                if (triangle.vertexIds[corner] == vertexId) {
                    return triangle.materialCoordinates[corner];
                }
            }
            return Vec2{};
        };
        const auto chartStrip = [&](const EdgeIncidence& incidence) {
            const Vec2 chartA = chartCoordinate(
                *incidence.triangle, incidence.from);
            const Vec2 chartB = chartCoordinate(
                *incidence.triangle, incidence.to);
            const Vec2 chartC = chartCoordinate(
                *incidence.triangle, incidence.opposite);
            const Vec2 chartEdge{chartB.x - chartA.x,
                                 chartB.y - chartA.y};
            const Vec2 chartOpposite{chartC.x - chartA.x,
                                     chartC.y - chartA.y};
            const double chartEdgeLength = std::hypot(
                chartEdge.x, chartEdge.y);
            const double chartAltitude = std::abs(
                chartEdge.x * chartOpposite.y
                - chartEdge.y * chartOpposite.x) / chartEdgeLength;
            return std::array<double, 2>{chartEdgeLength, chartAltitude};
        };
        const std::array<double, 2> firstChart = chartStrip(first);
        const std::array<double, 2> secondChart = chartStrip(second);
        const double compliance =
            firstChart[1] / (2.0 * firstRigidity * firstChart[0])
            + secondChart[1] / (2.0 * secondRigidity * secondChart[0]);
        const Vec3 edgeUnit = scaled(edgeVector, 1.0 / edgeLength);
        const Vec3 firstNormal = scaled(firstAreaVector,
                                        1.0 / firstAreaLength);
        const Vec3 secondNormal = scaled(secondAreaVector,
                                         1.0 / secondAreaLength);
        const double cosine = std::clamp(dot(firstNormal, secondNormal),
                                         -1.0, 1.0);
        const double sine = dot(
            edgeUnit, cross(firstNormal, secondNormal));
        const double restAngle = std::atan2(sine, cosine);
        if (!std::isfinite(compliance) || compliance < 0.0
            || !std::isfinite(restAngle)) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::MaterialIncompatible,
                EntityKind::Triangle,
                first.triangle->id,
                "fabric bending hinge is outside the representable XPBD range");
            continue;
        }

        assembly.definition.dihedrals.push_back(
            {{nodeIndices.at(first.from),
              nodeIndices.at(first.to),
              nodeIndices.at(first.opposite),
              nodeIndices.at(second.opposite)},
             restAngle,
             compliance});
        assembly.mappings.dihedralTriangleIds.push_back(
            {first.triangle->id, second.triangle->id});
    }
}

StructureMembraneMaterial convert(const FabricMaterial& material) {
    StructureMembraneMaterial result;
    result.warpStiffnessNewtonsPerMeter =
        material.warpStiffnessNewtonsPerMeter;
    result.weftStiffnessNewtonsPerMeter =
        material.weftStiffnessNewtonsPerMeter;
    result.couplingStiffnessNewtonsPerMeter = 0.0;
    result.shearStiffnessNewtonsPerMeter =
        material.shearStiffnessNewtonsPerMeter;
    result.warpPreTensionNewtonsPerMeter = 0.0;
    result.weftPreTensionNewtonsPerMeter = 0.0;
    result.dampingSeconds = material.dampingSeconds;
    result.compressionStiffnessRatio = 1.0;
    return result;
}

template <typename Vector>
std::optional<std::size_t> indexOf(const Vector& ids,
                                   StableId id) noexcept {
    const auto found = std::lower_bound(ids.begin(), ids.end(), id);
    if (found == ids.end() || *found != id) {
        return std::nullopt;
    }
    return static_cast<std::size_t>(found - ids.begin());
}

} // namespace

std::optional<std::size_t> SceneStructureMappings::nodeIndex(
    StableId vertexId) const noexcept {
    return indexOf(nodeVertexIds, vertexId);
}

std::optional<std::size_t> SceneStructureMappings::triangleIndex(
    StableId triangleId) const noexcept {
    return indexOf(triangleIds, triangleId);
}

std::optional<std::size_t> SceneStructureMappings::junctionNodeIndex(
    StableId junctionId) const noexcept {
    const std::optional<std::size_t> local = indexOf(
        nodeSuspensionJunctionIds, junctionId);
    if (!local) {
        return std::nullopt;
    }
    return nodeVertexIds.size() + *local;
}

std::optional<std::size_t> SceneStructureMappings::membraneIndex(
    StableId triangleId) const noexcept {
    return indexOf(membraneTriangleIds, triangleId);
}

std::optional<std::size_t> SceneStructureMappings::constraintIndex(
    StableId suspensionLineId) const noexcept {
    return indexOf(constraintSuspensionLineIds, suspensionLineId);
}

bool SceneStructureAssembly::ok() const noexcept {
    return assembled && diagnostics.empty();
}

SceneStructureAssembly assembleSceneStructure(
    const Scene& scene,
    const SceneStructureLimits& limits) {
    SceneStructureAssembly assembly;
    try {
        const ValidationReport validation = validateScene(scene);
        if (!validation.ok()) {
            for (const ValidationDiagnostic& diagnostic :
                 validation.diagnostics) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::InvalidScene,
                    diagnostic.entityKind,
                    diagnostic.entityId,
                    diagnostic.message,
                    diagnostic.code);
            }
            reject(assembly);
            return assembly;
        }
        if (!checkBounds(scene, limits, assembly)) {
            reject(assembly);
            return assembly;
        }

        const auto pilots = sortedById(scene.pilots);
        for (const Pilot* pilot : pilots) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedPilot,
                EntityKind::Pilot,
                pilot->id,
                "rigid pilot dynamics are not checkpoint-safe in simwing_structure");
        }

        for (const Seam* seam : sortedById(scene.seams)) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedSeam,
                EntityKind::Seam,
                seam->id,
                "seam topology requires a verified stitch and tributary load-sharing model");
        }

        const auto attachments = sortedById(scene.attachments);
        std::map<StableId, const Attachment*> attachmentsById;
        for (const Attachment* attachment : attachments) {
            attachmentsById.emplace(attachment->id, attachment);
            if (attachment->kind == AttachmentKind::PilotHarness) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::UnsupportedPilotHarness,
                    EntityKind::Attachment,
                    attachment->id,
                    "pilot-harness attachments require rigid-payload checkpoint support");
            }
        }

        const auto lines = sortedById(scene.suspensionLines);
        for (const SuspensionLine* line : lines) {
            const Attachment* start = attachmentsById.at(
                line->startAttachmentId);
            const Attachment* end = attachmentsById.at(line->endAttachmentId);
            if (start->kind == AttachmentKind::PilotHarness
                || end->kind == AttachmentKind::PilotHarness) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                    EntityKind::SuspensionLine,
                    line->id,
                    "pilot harness suspension topology is not checkpoint-safe");
            } else {
                const StableId startNodeId =
                    start->kind == AttachmentKind::SurfaceVertex
                    ? start->vertexId
                    : start->suspensionJunctionId;
                const StableId endNodeId =
                    end->kind == AttachmentKind::SurfaceVertex
                    ? end->vertexId
                    : end->suspensionJunctionId;
                if (start->kind == end->kind
                    && startNodeId == endNodeId) {
                    addDiagnostic(
                        assembly,
                        SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                        EntityKind::SuspensionLine,
                        line->id,
                        "suspension line endpoints must be distinct structural nodes");
                }
            }
        }
        if (!assembly.diagnostics.empty()) {
            reject(assembly);
            return assembly;
        }

        const auto vertices = sortedById(scene.vertices);
        const auto junctions = sortedById(scene.suspensionJunctions);
        const auto triangles = sortedById(scene.triangles);
        const auto materials = sortedById(scene.fabricMaterials);
        const auto lineMaterials = sortedById(scene.lineMaterials);

        std::map<StableId, std::size_t> nodeIndices;
        assembly.definition.nodes.reserve(vertices.size() + junctions.size());
        assembly.mappings.nodeVertexIds.reserve(vertices.size());
        for (const Vertex* vertex : vertices) {
            const std::size_t index = assembly.definition.nodes.size();
            nodeIndices.emplace(vertex->id, index);
            assembly.definition.nodes.push_back(
                {convert(vertex->positionMeters), 0.0, false});
            assembly.mappings.nodeVertexIds.push_back(vertex->id);
        }
        assembly.mappings.nodeSuspensionJunctionIds.reserve(
            junctions.size());
        for (const SuspensionJunction* junction : junctions) {
            const std::size_t index = assembly.definition.nodes.size();
            nodeIndices.emplace(junction->id, index);
            assembly.definition.nodes.push_back(
                {convert(junction->positionMeters),
                 junction->massKg,
                 junction->fixed});
            assembly.mappings.nodeSuspensionJunctionIds.push_back(
                junction->id);
        }

        std::map<StableId, const Vertex*> verticesById;
        for (const Vertex* vertex : vertices) {
            verticesById.emplace(vertex->id, vertex);
        }

        std::map<StableId, const FabricMaterial*> materialsById;
        for (const FabricMaterial* material : materials) {
            materialsById.emplace(material->id, material);
            if (!compatibleMaterial(*material)) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::MaterialIncompatible,
                    EntityKind::FabricMaterial,
                    material->id,
                    "fabric stiffness is too ill-conditioned for the XPBD membrane inverse");
            }
        }

        std::vector<long double> lumpedMasses(vertices.size(), 0.0L);
        long double totalFabricMass = 0.0L;
        assembly.definition.triangles.reserve(triangles.size());
        assembly.definition.membranes.reserve(triangles.size());
        assembly.mappings.triangleIds.reserve(triangles.size());
        assembly.mappings.membraneTriangleIds.reserve(triangles.size());
        for (const Triangle* triangle : triangles) {
            if (!compatibleMaterialCoordinates(*triangle)) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::DegenerateMaterialCoordinates,
                    EntityKind::Triangle,
                    triangle->id,
                    "triangle material coordinates are singular or ill-conditioned");
                continue;
            }

            StructureTriangleDefinition convertedTriangle;
            for (std::size_t corner = 0; corner < 3; ++corner) {
                convertedTriangle.nodes[corner] = nodeIndices.at(
                    triangle->vertexIds[corner]);
            }
            const std::size_t structureTriangle =
                assembly.definition.triangles.size();
            assembly.definition.triangles.push_back(convertedTriangle);
            assembly.mappings.triangleIds.push_back(triangle->id);

            StructureMembraneDefinition membrane;
            membrane.triangle = structureTriangle;
            for (std::size_t corner = 0; corner < 3; ++corner) {
                membrane.materialCoordinates[corner] = convert(
                    triangle->materialCoordinates[corner]);
            }
            membrane.material = convert(*materialsById.at(
                triangle->materialId));
            membrane.role = StructureMaterialRole::Bulk;
            assembly.definition.membranes.push_back(membrane);
            assembly.mappings.membraneTriangleIds.push_back(triangle->id);

            const double area = materialChartArea(*triangle);
            const double mass = area
                * materialsById.at(triangle->materialId)
                      ->arealDensityKgPerSquareMeter;
            if (!(mass > 0.0) || !std::isfinite(mass)) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::MappingOverflow,
                    EntityKind::Triangle,
                    triangle->id,
                    "triangle fabric mass is non-finite or underflows to zero");
                continue;
            }
            const long double share = static_cast<long double>(mass) / 3.0L;
            for (const std::size_t node : convertedTriangle.nodes) {
                lumpedMasses[node] += share;
            }
            totalFabricMass += static_cast<long double>(mass);
        }

        if (assembly.diagnostics.empty()) {
            assembleBending(triangles, verticesById, materialsById,
                            nodeIndices, assembly);
        }

        std::map<StableId, const LineMaterial*> lineMaterialsById;
        for (const LineMaterial* material : lineMaterials) {
            lineMaterialsById.emplace(material->id, material);
        }
        assembly.definition.constraints.reserve(lines.size());
        assembly.mappings.constraintSuspensionLineIds.reserve(lines.size());
        for (const SuspensionLine* line : lines) {
            const Attachment* start = attachmentsById.at(
                line->startAttachmentId);
            const Attachment* end = attachmentsById.at(line->endAttachmentId);
            const LineMaterial* material = lineMaterialsById.at(
                line->materialId);
            const double compliance = line->restLengthMeters
                / material->axialStiffnessNewtons;
            if (!std::isfinite(compliance) || compliance < 0.0) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::MaterialIncompatible,
                    EntityKind::LineMaterial,
                    material->id,
                    "line stiffness cannot be represented as finite XPBD compliance");
                continue;
            }
            assembly.definition.constraints.push_back(
                {StructureConstraintKind::Cable,
                 nodeIndices.at(
                     start->kind == AttachmentKind::SurfaceVertex
                         ? start->vertexId
                         : start->suspensionJunctionId),
                 nodeIndices.at(
                     end->kind == AttachmentKind::SurfaceVertex
                         ? end->vertexId
                         : end->suspensionJunctionId),
                 line->restLengthMeters,
                 compliance});
            assembly.mappings.constraintSuspensionLineIds.push_back(line->id);
        }

        for (std::size_t index = 0; index < lumpedMasses.size(); ++index) {
            const double mass = static_cast<double>(lumpedMasses[index]);
            if (!(mass > 0.0) || !std::isfinite(mass)) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::ZeroMassDynamicNode,
                    EntityKind::Vertex,
                    assembly.mappings.nodeVertexIds[index],
                    "scene vertex has no finite positive lumped fabric mass");
            } else {
                assembly.definition.nodes[index].massKg = mass;
            }
        }
        assembly.totalFabricMassKg = static_cast<double>(totalFabricMass);
        if (!std::isfinite(assembly.totalFabricMassKg)) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::MappingOverflow,
                EntityKind::Scene,
                invalidStableId,
                "total fabric mass exceeds the representable structure range");
        }

        if (!assembly.diagnostics.empty()) {
            reject(assembly);
            return assembly;
        }
        assembly.assembled = true;
        return assembly;
    } catch (const std::bad_alloc&) {
        addDiagnostic(
            assembly,
            SceneStructureDiagnosticCode::MappingOverflow,
            EntityKind::Scene,
            invalidStableId,
            "scene-to-structure assembly exhausted its bounded allocation");
        reject(assembly);
        return assembly;
    }
}

} // namespace simwing::fsi

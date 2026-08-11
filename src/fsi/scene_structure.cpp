#include "scene_structure.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <queue>
#include <set>
#include <stdexcept>
#include <string_view>
#include <tuple>
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
    assembly.settings = {};
    assembly.totalFabricMassKg = 0.0;
    assembly.totalSeamMassKg = 0.0;
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
    std::size_t seamConstraintCount = 0;
    bool constraintCountOverflow = false;
    for (const Seam& seam : scene.seams) {
        const std::size_t count = seam.firstOrderedVertexIds.size();
        if (count > std::numeric_limits<std::size_t>::max() / 3) {
            constraintCountOverflow = true;
            break;
        }
        const std::size_t generated = 3 * count - 2;
        if (generated > std::numeric_limits<std::size_t>::max()
                            - seamConstraintCount) {
            constraintCountOverflow = true;
            break;
        }
        seamConstraintCount += generated;
    }
    if (constraintCountOverflow
        || scene.suspensionLines.size()
            > std::numeric_limits<std::size_t>::max()
                - seamConstraintCount
        || scene.suspensionLines.size() + seamConstraintCount
            > limits.maximumConstraints) {
        addDiagnostic(
            assembly,
            SceneStructureDiagnosticCode::MappingOverflow,
            EntityKind::Scene,
            invalidStableId,
            "scene seam and suspension-line constraints exceed the structure bound");
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
        && checkedMappingBytes(scene.suspensionLines.size(), sizeof(StableId),
                               remaining)
        && checkedMappingBytes(
            scene.seams.size(),
            sizeof(SceneStructureMappings::SeamConstraintRange), remaining)
        && checkedMappingBytes(scene.attachments.size(), sizeof(StableId),
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

void assembleSeams(
    const std::vector<const Seam*>& seams,
    const std::vector<const SeamMaterial*>& materials,
    const std::map<StableId, const Vertex*>& verticesById,
    const std::map<StableId, std::size_t>& nodeIndices,
    const SceneStructureSettings& settings,
    std::vector<long double>& lumpedMasses,
    long double& totalSeamMass,
    SceneStructureAssembly& assembly) {
    std::map<StableId, const SeamMaterial*> materialsById;
    for (const SeamMaterial* material : materials) {
        materialsById.emplace(material->id, material);
    }

    assembly.mappings.constraintSeamRanges.reserve(seams.size());
    for (const Seam* seam : seams) {
        const std::size_t diagnosticCount = assembly.diagnostics.size();
        const SeamMaterial* material = materialsById.at(seam->materialId);
        const std::size_t count = seam->firstOrderedVertexIds.size();
        std::array<std::vector<double>, 2> segmentLengths;
        for (auto& lengths : segmentLengths) {
            lengths.reserve(count - 1);
        }

        for (std::size_t index = 0; index < count; ++index) {
            const Vec3& first = verticesById.at(
                seam->firstOrderedVertexIds[index])->positionMeters;
            const Vec3& second = verticesById.at(
                seam->secondOrderedVertexIds[index])->positionMeters;
            const double separation = length(subtract(second, first));
            if (!std::isfinite(separation)
                || separation > settings.seamCoincidenceToleranceMeters) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::UnsupportedSeam,
                    EntityKind::Seam,
                    seam->id,
                    "paired seam vertices are not coincident within the structural tolerance");
                break;
            }
        }
        for (std::size_t chain = 0; chain < 2; ++chain) {
            const auto& ids = chain == 0
                ? seam->firstOrderedVertexIds
                : seam->secondOrderedVertexIds;
            for (std::size_t index = 0; index + 1 < count; ++index) {
                const Vec3& first = verticesById.at(ids[index])
                                        ->positionMeters;
                const Vec3& second = verticesById.at(ids[index + 1])
                                         ->positionMeters;
                const double segmentLength = length(subtract(second, first));
                const double compliance = 2.0 * segmentLength
                    / material->axialStiffnessNewtons;
                if (!(segmentLength > 0.0)
                    || !std::isfinite(segmentLength)
                    || !std::isfinite(compliance)) {
                    addDiagnostic(
                        assembly,
                        SceneStructureDiagnosticCode::MaterialIncompatible,
                        EntityKind::Seam,
                        seam->id,
                        "seam segment stiffness cannot be represented as finite XPBD compliance");
                    break;
                }
                segmentLengths[chain].push_back(segmentLength);
            }
        }
        if (assembly.diagnostics.size() != diagnosticCount) {
            continue;
        }

        const std::size_t firstConstraint =
            assembly.definition.constraints.size();
        for (std::size_t index = 0; index < count; ++index) {
            if (seam->firstOrderedVertexIds[index]
                == seam->secondOrderedVertexIds[index]) {
                continue;
            }
            assembly.definition.constraints.push_back(
                {StructureConstraintKind::SeamStitch,
                 nodeIndices.at(seam->firstOrderedVertexIds[index]),
                 nodeIndices.at(seam->secondOrderedVertexIds[index]),
                 0.0,
                 0.0});
        }
        for (std::size_t chain = 0; chain < 2; ++chain) {
            const auto& ids = chain == 0
                ? seam->firstOrderedVertexIds
                : seam->secondOrderedVertexIds;
            for (std::size_t index = 0; index + 1 < count; ++index) {
                const double segmentLength = segmentLengths[chain][index];
                assembly.definition.constraints.push_back(
                    {StructureConstraintKind::Distance,
                     nodeIndices.at(ids[index]),
                     nodeIndices.at(ids[index + 1]),
                     segmentLength,
                     2.0 * segmentLength
                         / material->axialStiffnessNewtons});
            }
        }
        assembly.mappings.constraintSeamRanges.push_back(
            {seam->id,
             firstConstraint,
             assembly.definition.constraints.size() - firstConstraint});

        for (std::size_t index = 0; index + 1 < count; ++index) {
            const double centrelineLength = 0.5
                * (segmentLengths[0][index] + segmentLengths[1][index]);
            const double segmentMass = centrelineLength
                * material->linearDensityKgPerMeter;
            const long double nodeShare =
                static_cast<long double>(segmentMass) / 4.0L;
            for (const StableId id : {
                     seam->firstOrderedVertexIds[index],
                     seam->firstOrderedVertexIds[index + 1],
                     seam->secondOrderedVertexIds[index],
                     seam->secondOrderedVertexIds[index + 1]}) {
                lumpedMasses[nodeIndices.at(id)] += nodeShare;
            }
            totalSeamMass += static_cast<long double>(segmentMass);
        }
    }
}

struct SuspensionEndpointKey {
    StructureSuspensionEndpointKind kind =
        StructureSuspensionEndpointKind::SurfaceAttachment;
    StableId stableId = invalidStableId;

    auto operator<=>(const SuspensionEndpointKey&) const = default;
};

SuspensionEndpointKey endpointKey(const Attachment& attachment) {
    switch (attachment.kind) {
    case AttachmentKind::SurfaceVertex:
        return {StructureSuspensionEndpointKind::SurfaceAttachment,
                attachment.id};
    case AttachmentKind::SuspensionJunction:
        return {StructureSuspensionEndpointKind::Junction,
                attachment.suspensionJunctionId};
    case AttachmentKind::PilotHarness:
        return {StructureSuspensionEndpointKind::PilotHarness,
                attachment.id};
    }
    throw std::logic_error("validated attachment has an unknown kind");
}

StructureSuspensionEndpointDefinition convert(
    const SuspensionEndpointKey& endpoint) {
    return {endpoint.kind, endpoint.stableId};
}

bool assemblePilotSuspension(
    const Pilot& pilot,
    const std::vector<const Attachment*>& attachments,
    const std::vector<const SuspensionLine*>& lines,
    const std::vector<const SuspensionJunction*>& junctions,
    const std::map<StableId, std::size_t>& nodeIndices,
    const std::map<StableId, const LineMaterial*>& lineMaterialsById,
    const SceneStructureSettings& settings,
    SceneStructureAssembly& assembly) {
    if (lines.empty()) {
        addDiagnostic(
            assembly, SceneStructureDiagnosticCode::UnsupportedPilot,
            EntityKind::Pilot, pilot.id,
            "rigid pilot requires a connected suspension tree");
        return false;
    }

    std::map<StableId, const Attachment*> attachmentsById;
    for (const Attachment* attachment : attachments) {
        attachmentsById.emplace(attachment->id, attachment);
    }
    std::map<StableId, const SuspensionJunction*> junctionsById;
    for (const SuspensionJunction* junction : junctions) {
        junctionsById.emplace(junction->id, junction);
    }

    using Neighbor = std::pair<SuspensionEndpointKey,
                               const SuspensionLine*>;
    std::map<SuspensionEndpointKey, std::vector<Neighbor>> adjacency;
    std::set<std::pair<SuspensionEndpointKey, SuspensionEndpointKey>> pairs;
    for (const SuspensionLine* line : lines) {
        const SuspensionEndpointKey first = endpointKey(
            *attachmentsById.at(line->startAttachmentId));
        const SuspensionEndpointKey second = endpointKey(
            *attachmentsById.at(line->endAttachmentId));
        auto canonical = std::pair{first, second};
        if (canonical.second < canonical.first) {
            std::swap(canonical.first, canonical.second);
        }
        if (first == second || !pairs.insert(canonical).second) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                EntityKind::SuspensionLine, line->id,
                "suspension tree contains a self-edge or duplicate endpoint pair");
            continue;
        }
        adjacency[first].push_back({second, line});
        adjacency[second].push_back({first, line});
    }
    for (auto& [endpoint, neighbors] : adjacency) {
        std::ranges::sort(
            neighbors, [](const Neighbor& left, const Neighbor& right) {
                return std::tuple{left.first.kind, left.first.stableId,
                                  left.second->id}
                    < std::tuple{right.first.kind, right.first.stableId,
                                 right.second->id};
            });
        if (endpoint.kind
                == StructureSuspensionEndpointKind::SurfaceAttachment
            && neighbors.size() != 1) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                EntityKind::Attachment, endpoint.stableId,
                "surface suspension attachment must be a tree leaf");
        }
        if (endpoint.kind == StructureSuspensionEndpointKind::Junction) {
            const SuspensionJunction* junction =
                junctionsById.at(endpoint.stableId);
            if (junction->fixed || neighbors.size() < 2) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                    EntityKind::SuspensionJunction, junction->id,
                    junction->fixed
                        ? "pilot suspension junction must be dynamic"
                        : "pilot suspension junction cannot be a terminal leaf");
            }
        }
    }
    if (!assembly.diagnostics.empty()) {
        return false;
    }

    std::set<SuspensionEndpointKey> visited;
    std::map<StableId,
             std::pair<SuspensionEndpointKey, SuspensionEndpointKey>>
        orientedLines;
    for (const auto& [componentStart, ignored] : adjacency) {
        static_cast<void>(ignored);
        if (visited.contains(componentStart)) {
            continue;
        }
        std::vector<SuspensionEndpointKey> component;
        std::queue<SuspensionEndpointKey> pending;
        pending.push(componentStart);
        visited.insert(componentStart);
        std::size_t degreeSum = 0;
        while (!pending.empty()) {
            const SuspensionEndpointKey current = pending.front();
            pending.pop();
            component.push_back(current);
            degreeSum += adjacency.at(current).size();
            for (const Neighbor& neighbor : adjacency.at(current)) {
                if (visited.insert(neighbor.first).second) {
                    pending.push(neighbor.first);
                }
            }
        }
        std::vector<SuspensionEndpointKey> roots;
        std::size_t surfaceLeaves = 0;
        for (const SuspensionEndpointKey& endpoint : component) {
            if (endpoint.kind
                == StructureSuspensionEndpointKind::PilotHarness) {
                roots.push_back(endpoint);
            } else if (endpoint.kind
                       == StructureSuspensionEndpointKind::SurfaceAttachment) {
                ++surfaceLeaves;
            }
        }
        if (roots.size() != 1 || surfaceLeaves == 0
            || degreeSum / 2 != component.size() - 1) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                EntityKind::SuspensionLine, lines.front()->id,
                "each suspension component must be an acyclic tree with one harness root and at least one surface leaf");
            continue;
        }

        std::set<SuspensionEndpointKey> rooted;
        std::queue<SuspensionEndpointKey> outward;
        rooted.insert(roots.front());
        outward.push(roots.front());
        while (!outward.empty()) {
            const SuspensionEndpointKey parent = outward.front();
            outward.pop();
            for (const Neighbor& neighbor : adjacency.at(parent)) {
                if (!rooted.insert(neighbor.first).second) {
                    continue;
                }
                // SuspensionSystem segments point from the canopy leaf toward
                // the rigid-payload root.
                orientedLines.emplace(
                    neighbor.second->id,
                    std::pair{neighbor.first, parent});
                outward.push(neighbor.first);
            }
        }
    }
    if (!assembly.diagnostics.empty()
        || orientedLines.size() != lines.size()) {
        if (assembly.diagnostics.empty()) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                EntityKind::SuspensionLine, lines.front()->id,
                "suspension tree could not be oriented toward the pilot");
        }
        return false;
    }

    StructureSuspensionDefinition suspension;
    suspension.pilotStableId = pilot.id;
    suspension.pilotMassKg = pilot.massKg;
    suspension.pilotInitialCenterOfMassWorldMeters = convert(
        pilot.centerOfMassPositionMeters);
    suspension.pilotInitialLinearVelocityMetersPerSecond = convert(
        pilot.linearVelocityMetersPerSecond);
    suspension.pilotInitialBodyToWorld = {
        pilot.bodyToWorld.w, pilot.bodyToWorld.x,
        pilot.bodyToWorld.y, pilot.bodyToWorld.z};
    suspension.pilotPrincipalInertiaKgSquareMeters = convert(
        pilot.principalInertiaKgSquareMeters);
    suspension.solverIterations = settings.suspensionSolverIterations;
    suspension.attachmentTolerance =
        settings.suspensionAttachmentTolerance;
    suspension.minimumLineLengthMeters =
        settings.suspensionMinimumLineLengthMeters;
    suspension.maximumLineResidualMeters =
        settings.suspensionMaximumLineResidualMeters;
    suspension.maximumControlWorkJoules =
        settings.suspensionMaximumControlWorkJoules;

    for (const auto& [endpoint, neighbors] : adjacency) {
        static_cast<void>(neighbors);
        if (endpoint.kind
            == StructureSuspensionEndpointKind::SurfaceAttachment) {
            const Attachment* attachment = attachmentsById.at(
                endpoint.stableId);
            suspension.attachments.push_back(
                {endpoint.stableId, nodeIndices.at(attachment->vertexId)});
        } else if (endpoint.kind
                   == StructureSuspensionEndpointKind::Junction) {
            suspension.junctions.push_back(
                {endpoint.stableId, nodeIndices.at(endpoint.stableId)});
        } else {
            const Attachment* attachment = attachmentsById.at(
                endpoint.stableId);
            if (attachment->pilotId != pilot.id) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::UnsupportedPilotHarness,
                    EntityKind::Attachment, attachment->id,
                    "harness root belongs to a different pilot");
                continue;
            }
            suspension.harnessPoints.push_back(
                {endpoint.stableId,
                 convert(attachment->pilotLocalPositionMeters)});
            assembly.mappings.pilotHarnessAttachmentIds.push_back(
                endpoint.stableId);
        }
    }
    for (const SuspensionLine* line : lines) {
        const auto oriented = orientedLines.at(line->id);
        suspension.segments.push_back(
            {line->id, convert(oriented.first), convert(oriented.second),
             line->restLengthMeters,
             lineMaterialsById.at(line->materialId)
                 ->axialStiffnessNewtons,
             0.0, static_cast<std::uint32_t>(line->role)});
        assembly.mappings.suspensionSegmentLineIds.push_back(line->id);
    }
    if (!assembly.diagnostics.empty()) {
        return false;
    }
    assembly.definition.suspension = std::move(suspension);
    return true;
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

std::optional<SceneStructureMappings::SeamConstraintRange>
SceneStructureMappings::seamConstraintRange(
    StableId seamId) const noexcept {
    const auto found = std::lower_bound(
        constraintSeamRanges.begin(), constraintSeamRanges.end(), seamId,
        [](const SeamConstraintRange& range, StableId id) {
            return range.seamId < id;
        });
    if (found == constraintSeamRanges.end() || found->seamId != seamId) {
        return std::nullopt;
    }
    return *found;
}

std::optional<std::size_t> SceneStructureMappings::suspensionSegmentIndex(
    StableId suspensionLineId) const noexcept {
    return indexOf(suspensionSegmentLineIds, suspensionLineId);
}

std::optional<std::size_t> SceneStructureMappings::pilotHarnessIndex(
    StableId attachmentId) const noexcept {
    return indexOf(pilotHarnessAttachmentIds, attachmentId);
}

bool SceneStructureAssembly::ok() const noexcept {
    return assembled && diagnostics.empty();
}

SceneStructureAssembly assembleSceneStructure(
    const Scene& scene,
    const SceneStructureLimits& limits,
    const SceneStructureSettings& settings) {
    SceneStructureAssembly assembly;
    assembly.settings = settings;
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
        if (!(settings.seamCoincidenceToleranceMeters >= 0.0)
            || !std::isfinite(settings.seamCoincidenceToleranceMeters)) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedSeam,
                EntityKind::Scene, invalidStableId,
                "invalid seam coincidence tolerance");
        }
        const bool invalidSuspensionSettings =
            settings.suspensionSolverIterations <= 0
            || !(settings.suspensionAttachmentTolerance > 0.0)
            || !std::isfinite(settings.suspensionAttachmentTolerance)
            || !(settings.suspensionMinimumLineLengthMeters > 0.0)
            || !std::isfinite(settings.suspensionMinimumLineLengthMeters)
            || settings.suspensionMaximumLineResidualMeters < 0.0
            || !std::isfinite(
                settings.suspensionMaximumLineResidualMeters)
            || settings.suspensionMaximumControlWorkJoules < 0.0
            || !std::isfinite(
                settings.suspensionMaximumControlWorkJoules);
        if (invalidSuspensionSettings) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                EntityKind::Scene, invalidStableId,
                "invalid suspension solver settings");
        }
        if (settings.fabricSelfContact) {
            const StructureFabricContactDefinition& contact =
                *settings.fabricSelfContact;
            if (!(contact.halfThicknessMeters > 0.0)
                || !std::isfinite(contact.halfThicknessMeters)
                || contact.normalComplianceMetersPerNewton < 0.0
                || !std::isfinite(contact.normalComplianceMetersPerNewton)
                || contact.staticFriction < 0.0
                || !std::isfinite(contact.staticFriction)
                || contact.dynamicFriction < 0.0
                || !std::isfinite(contact.dynamicFriction)
                || contact.dynamicFriction > contact.staticFriction) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::MaterialIncompatible,
                    EntityKind::Scene, invalidStableId,
                    "invalid explicit fabric contact settings");
            }
        }

        const auto pilots = sortedById(scene.pilots);
        if (pilots.size() > 1) {
            for (const Pilot* pilot : pilots) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::UnsupportedPilot,
                    EntityKind::Pilot,
                    pilot->id,
                    "one Structure instance supports exactly one rigid pilot");
            }
        }

        std::set<StableId> fabricVertexIds;
        for (const Triangle& triangle : scene.triangles) {
            fabricVertexIds.insert(triangle.vertexIds.begin(),
                                   triangle.vertexIds.end());
        }
        for (const Seam* seam : sortedById(scene.seams)) {
            const auto checkChain = [&](const std::vector<StableId>& chain) {
                for (const StableId vertexId : chain) {
                    if (!fabricVertexIds.contains(vertexId)) {
                        addDiagnostic(
                            assembly,
                            SceneStructureDiagnosticCode::UnsupportedSeam,
                            EntityKind::Seam,
                            seam->id,
                            "seam vertices must belong to fabric triangles");
                        return;
                    }
                }
            };
            checkChain(seam->firstOrderedVertexIds);
            checkChain(seam->secondOrderedVertexIds);
        }

        const auto attachments = sortedById(scene.attachments);
        std::map<StableId, const Attachment*> attachmentsById;
        for (const Attachment* attachment : attachments) {
            attachmentsById.emplace(attachment->id, attachment);
        }

        const auto lines = sortedById(scene.suspensionLines);
        for (const SuspensionLine* line : lines) {
            const Attachment* start = attachmentsById.at(
                line->startAttachmentId);
            const Attachment* end = attachmentsById.at(line->endAttachmentId);
            if (pilots.empty()
                && (start->kind == AttachmentKind::PilotHarness
                    || end->kind == AttachmentKind::PilotHarness)) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::UnsupportedSuspensionTopology,
                    EntityKind::SuspensionLine,
                    line->id,
                    "a pilot harness line requires exactly one rigid pilot");
            } else if (start->kind != AttachmentKind::PilotHarness
                       && end->kind != AttachmentKind::PilotHarness) {
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

        const auto allVertices = sortedById(scene.vertices);
        const auto junctions = sortedById(scene.suspensionJunctions);
        const auto triangles = sortedById(scene.triangles);
        const auto materials = sortedById(scene.fabricMaterials);
        const auto seams = sortedById(scene.seams);
        const auto seamMaterials = sortedById(scene.seamMaterials);
        const auto lineMaterials = sortedById(scene.lineMaterials);
        std::set<StableId> structuralVertexIds;
        for (const Triangle* triangle : triangles) {
            structuralVertexIds.insert(triangle->vertexIds.begin(),
                                       triangle->vertexIds.end());
        }
        for (const Attachment* attachment : attachments) {
            if (attachment->kind == AttachmentKind::SurfaceVertex) {
                structuralVertexIds.insert(attachment->vertexId);
            }
        }
        std::vector<const Vertex*> vertices;
        vertices.reserve(structuralVertexIds.size());
        for (const Vertex* vertex : allVertices) {
            if (structuralVertexIds.contains(vertex->id)) {
                vertices.push_back(vertex);
            }
        }

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
        long double totalSeamMass = 0.0L;
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
        if (pilots.empty()) {
            assembly.definition.constraints.reserve(lines.size());
            assembly.mappings.constraintSuspensionLineIds.reserve(
                lines.size());
            for (const SuspensionLine* line : lines) {
                const Attachment* start = attachmentsById.at(
                    line->startAttachmentId);
                const Attachment* end = attachmentsById.at(
                    line->endAttachmentId);
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
                assembly.mappings.constraintSuspensionLineIds.push_back(
                    line->id);
            }
        } else if (pilots.size() == 1) {
            static_cast<void>(assemblePilotSuspension(
                *pilots.front(), attachments, lines, junctions,
                nodeIndices, lineMaterialsById, settings, assembly));
        }

        if (assembly.diagnostics.empty()) {
            assembleSeams(
                seams, seamMaterials, verticesById, nodeIndices, settings,
                lumpedMasses, totalSeamMass, assembly);
        }

        for (std::size_t index = 0; index < lumpedMasses.size(); ++index) {
            const double mass = static_cast<double>(lumpedMasses[index]);
            if (!(mass > 0.0) || !std::isfinite(mass)) {
                addDiagnostic(
                    assembly,
                    SceneStructureDiagnosticCode::ZeroMassDynamicNode,
                    EntityKind::Vertex,
                    assembly.mappings.nodeVertexIds[index],
                    "scene vertex has no finite positive lumped fabric or seam mass");
            } else {
                assembly.definition.nodes[index].massKg = mass;
            }
        }
        assembly.totalFabricMassKg = static_cast<double>(totalFabricMass);
        assembly.totalSeamMassKg = static_cast<double>(totalSeamMass);
        assembly.definition.fabricSelfContact =
            settings.fabricSelfContact;
        if (!std::isfinite(assembly.totalFabricMassKg)) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::MappingOverflow,
                EntityKind::Scene,
                invalidStableId,
                "total fabric mass exceeds the representable structure range");
        }
        if (!std::isfinite(assembly.totalSeamMassKg)) {
            addDiagnostic(
                assembly,
                SceneStructureDiagnosticCode::MappingOverflow,
                EntityKind::Scene,
                invalidStableId,
                "total seam mass exceeds the representable structure range");
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

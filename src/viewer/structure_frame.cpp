#include "structure_frame.h"

#include "scene_structure.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace simwing::viewer {
namespace {

template <class Entity>
void requireUniqueNonzeroIds(
    const std::vector<Entity>& entities,
    const char* label,
    const auto& idOf) {
    std::unordered_set<std::uint64_t> ids;
    ids.reserve(entities.size());
    for (const Entity& entity : entities) {
        const std::uint64_t id = idOf(entity);
        if (id == 0 || !ids.insert(id).second) {
            throw std::invalid_argument(
                std::string("Structure frame ") + label
                + " IDs must be nonzero and unique");
        }
    }
}

void requireUniqueNonzeroIds(
    const std::vector<std::uint64_t>& ids,
    const char* label) {
    requireUniqueNonzeroIds(
        ids, label, [](std::uint64_t id) { return id; });
}

[[nodiscard]] Vec3d toViewer(const fsi::StructureVector3& value) {
    return {value.x, value.y, value.z};
}

[[nodiscard]] fsi::StructureVector3 subtract(
    const fsi::StructureVector3& first,
    const fsi::StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

[[nodiscard]] double length(const fsi::StructureVector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

[[nodiscard]] double constraintViolation(
    const fsi::StructureConstraintDefinition& constraint,
    double currentLength) {
    const double extension = currentLength - constraint.restLengthMeters;
    switch (constraint.kind) {
    case fsi::StructureConstraintKind::Distance:
    case fsi::StructureConstraintKind::SuspensionTie:
        return std::abs(extension);
    case fsi::StructureConstraintKind::Cable:
        return std::max(0.0, extension);
    }
    throw std::invalid_argument("Unknown structure constraint kind");
}

void addGlobalScalar(
    DiagnosticFrame& frame,
    std::string name,
    std::string unit,
    double value) {
    frame.scalarFields.push_back(
        {std::move(name), std::move(unit), FieldAssociation::Global, {value}});
}

void addGlobalVector(
    DiagnosticFrame& frame,
    std::string name,
    std::string unit,
    const fsi::StructureVector3& value) {
    frame.vectorFields.push_back(
        {std::move(name), std::move(unit), FieldAssociation::Global,
         {toViewer(value)}});
}

[[nodiscard]] bool sameMappings(
    const fsi::SceneStructureMappings& first,
    const fsi::SceneStructureMappings& second) {
    return first.nodeVertexIds == second.nodeVertexIds
        && first.nodeSuspensionJunctionIds
               == second.nodeSuspensionJunctionIds
        && first.triangleIds == second.triangleIds
        && first.membraneTriangleIds == second.membraneTriangleIds
        && first.constraintSuspensionLineIds
               == second.constraintSuspensionLineIds
        && first.suspensionSegmentLineIds
               == second.suspensionSegmentLineIds
        && first.pilotHarnessAttachmentIds
               == second.pilotHarnessAttachmentIds
        && first.dihedralTriangleIds == second.dihedralTriangleIds;
}

[[nodiscard]] bool sameDefinition(
    const fsi::StructureDefinition& first,
    const fsi::StructureDefinition& second) {
    if (first.nodes.size() != second.nodes.size()
        || first.triangles.size() != second.triangles.size()
        || first.constraints.size() != second.constraints.size()
        || first.membranes.size() != second.membranes.size()
        || first.dihedrals.size() != second.dihedrals.size()) {
        return false;
    }
    for (std::size_t index = 0; index < first.nodes.size(); ++index) {
        const fsi::StructureNodeDefinition& left = first.nodes[index];
        const fsi::StructureNodeDefinition& right = second.nodes[index];
        if (left.positionMeters != right.positionMeters
            || left.massKg != right.massKg || left.fixed != right.fixed) {
            return false;
        }
    }
    for (std::size_t index = 0; index < first.triangles.size(); ++index) {
        if (first.triangles[index].nodes != second.triangles[index].nodes) {
            return false;
        }
    }
    for (std::size_t index = 0; index < first.constraints.size(); ++index) {
        const fsi::StructureConstraintDefinition& left =
            first.constraints[index];
        const fsi::StructureConstraintDefinition& right =
            second.constraints[index];
        if (left.kind != right.kind || left.firstNode != right.firstNode
            || left.secondNode != right.secondNode
            || left.restLengthMeters != right.restLengthMeters
            || left.complianceMetersPerNewton
                   != right.complianceMetersPerNewton) {
            return false;
        }
    }
    for (std::size_t index = 0; index < first.membranes.size(); ++index) {
        const fsi::StructureMembraneDefinition& left =
            first.membranes[index];
        const fsi::StructureMembraneDefinition& right =
            second.membranes[index];
        if (left.triangle != right.triangle || left.role != right.role) {
            return false;
        }
        for (std::size_t corner = 0; corner < 3; ++corner) {
            if (left.materialCoordinates[corner].x
                    != right.materialCoordinates[corner].x
                || left.materialCoordinates[corner].y
                       != right.materialCoordinates[corner].y) {
                return false;
            }
        }
        const fsi::StructureMembraneMaterial& leftMaterial = left.material;
        const fsi::StructureMembraneMaterial& rightMaterial = right.material;
        if (leftMaterial.warpStiffnessNewtonsPerMeter
                != rightMaterial.warpStiffnessNewtonsPerMeter
            || leftMaterial.weftStiffnessNewtonsPerMeter
                   != rightMaterial.weftStiffnessNewtonsPerMeter
            || leftMaterial.couplingStiffnessNewtonsPerMeter
                   != rightMaterial.couplingStiffnessNewtonsPerMeter
            || leftMaterial.shearStiffnessNewtonsPerMeter
                   != rightMaterial.shearStiffnessNewtonsPerMeter
            || leftMaterial.warpPreTensionNewtonsPerMeter
                   != rightMaterial.warpPreTensionNewtonsPerMeter
            || leftMaterial.weftPreTensionNewtonsPerMeter
                   != rightMaterial.weftPreTensionNewtonsPerMeter
            || leftMaterial.dampingSeconds
                   != rightMaterial.dampingSeconds
            || leftMaterial.compressionStiffnessRatio
                   != rightMaterial.compressionStiffnessRatio) {
            return false;
        }
    }
    for (std::size_t index = 0; index < first.dihedrals.size(); ++index) {
        const fsi::StructureDihedralDefinition& left =
            first.dihedrals[index];
        const fsi::StructureDihedralDefinition& right =
            second.dihedrals[index];
        if (left.nodes != right.nodes
            || left.restAngleRadians != right.restAngleRadians
            || left.complianceRadiansPerNewtonMeter
                   != right.complianceRadiansPerNewtonMeter) {
            return false;
        }
    }
    return true;
}

} // namespace

StructureFrameMapping::StructureFrameMapping(
    const fsi::Structure& structure,
    StructureFrameMappingDefinition definition)
    : structureFingerprint_(structure.definitionFingerprint()),
      vertexStableIds_(std::move(definition.vertexStableIds)),
      triangles_(std::move(definition.triangles)),
      lines_(std::move(definition.lines)) {
    const fsi::StructureDefinition& topology = structure.definition();
    const std::size_t harnessCount = topology.suspension
        ? topology.suspension->harnessPoints.size() : 0;
    const std::size_t suspensionLineCount = topology.suspension
        ? topology.suspension->segments.size() : 0;
    constexpr std::size_t maximumViewerVertices =
        std::numeric_limits<std::uint32_t>::max();
    if (harnessCount > maximumViewerVertices
        || topology.nodes.size()
               > maximumViewerVertices - harnessCount) {
        throw std::invalid_argument(
            "Structure frame has too many vertices for the viewer protocol");
    }
    if (vertexStableIds_.size() != topology.nodes.size() + harnessCount
        || triangles_.size() != topology.triangles.size()
        || lines_.size()
               != topology.constraints.size() + suspensionLineCount) {
        throw std::invalid_argument(
            "Structure frame mapping sizes do not match structure topology");
    }

    requireUniqueNonzeroIds(vertexStableIds_, "vertex");
    requireUniqueNonzeroIds(
        triangles_, "triangle",
        [](const StructureFrameTriangleMapping& triangle) {
            return triangle.stableId;
        });
    requireUniqueNonzeroIds(
        lines_, "line", [](const StructureFrameLineMapping& line) {
            return line.stableId;
        });

    for (const StructureFrameTriangleMapping& triangle : triangles_) {
        if (triangle.negativeRegionId == 0
            || triangle.positiveRegionId == 0) {
            throw std::invalid_argument(
                "Structure frame triangle needs nonzero side regions");
        }
    }

    for (const fsi::StructureTriangleDefinition& triangle :
         topology.triangles) {
        for (const std::size_t node : triangle.nodes) {
            if (node >= topology.nodes.size()) {
                throw std::invalid_argument(
                    "Structure frame triangle references an unknown node");
            }
        }
    }
    for (const fsi::StructureConstraintDefinition& constraint :
         topology.constraints) {
        if (constraint.firstNode >= topology.nodes.size()
            || constraint.secondNode >= topology.nodes.size()) {
            throw std::invalid_argument(
                "Structure frame line references an unknown node");
        }
    }
    for (std::size_t index = topology.constraints.size();
         index < lines_.size(); ++index) {
        if (lines_[index].vertex0 >= vertexStableIds_.size()
            || lines_[index].vertex1 >= vertexStableIds_.size()
            || lines_[index].vertex0 == lines_[index].vertex1) {
            throw std::invalid_argument(
                "Structure frame suspension line has invalid endpoints");
        }
    }
}

std::uint64_t StructureFrameMapping::structureFingerprint() const noexcept {
    return structureFingerprint_;
}

const std::vector<std::uint64_t>&
StructureFrameMapping::vertexStableIds() const noexcept {
    return vertexStableIds_;
}

const std::vector<StructureFrameTriangleMapping>&
StructureFrameMapping::triangles() const noexcept {
    return triangles_;
}

const std::vector<StructureFrameLineMapping>&
StructureFrameMapping::lines() const noexcept {
    return lines_;
}

StructureFrameMapping makeStructureFrameMapping(
    const fsi::Scene& scene,
    const fsi::SceneStructureAssembly& assembly,
    const fsi::Structure& structure) {
    if (!assembly.ok()) {
        throw std::invalid_argument(
            "Cannot map a failed scene-to-structure assembly");
    }

    // Reassembly is intentional here. It validates the scene and provides a
    // transactional oracle for every conversion rule owned by scene_structure
    // without copying those material, mass, or topology rules into the viewer.
    const fsi::SceneStructureAssembly canonical =
        fsi::assembleSceneStructure(scene, {}, assembly.settings);
    if (!canonical.ok()) {
        throw std::invalid_argument(
            "Cannot map a scene that fails structural assembly");
    }
    if (!sameMappings(assembly.mappings, canonical.mappings)) {
        throw std::invalid_argument(
            "Scene structure stable-ID mappings do not match the scene");
    }

    const fsi::Structure canonicalStructure(canonical.definition);
    const fsi::Structure assemblyStructure(assembly.definition);
    if (!sameDefinition(canonical.definition, assembly.definition)
        || !sameDefinition(assembly.definition, structure.definition())
        || canonicalStructure.definitionFingerprint()
               != assemblyStructure.definitionFingerprint()
        || assemblyStructure.definitionFingerprint()
               != structure.definitionFingerprint()) {
        throw std::invalid_argument(
            "Scene, assembly, and Structure definitions do not match");
    }

    std::map<fsi::StableId, const fsi::Triangle*> trianglesById;
    for (const fsi::Triangle& triangle : scene.triangles) {
        trianglesById.emplace(triangle.id, &triangle);
    }
    std::map<fsi::StableId, const fsi::SuspensionLine*> linesById;
    for (const fsi::SuspensionLine& line : scene.suspensionLines) {
        linesById.emplace(line.id, &line);
    }

    StructureFrameMappingDefinition definition;
    definition.vertexStableIds = canonical.mappings.nodeVertexIds;
    definition.vertexStableIds.insert(
        definition.vertexStableIds.end(),
        canonical.mappings.nodeSuspensionJunctionIds.begin(),
        canonical.mappings.nodeSuspensionJunctionIds.end());
    definition.vertexStableIds.insert(
        definition.vertexStableIds.end(),
        canonical.mappings.pilotHarnessAttachmentIds.begin(),
        canonical.mappings.pilotHarnessAttachmentIds.end());
    definition.triangles.reserve(canonical.mappings.triangleIds.size());
    for (const fsi::StableId id : canonical.mappings.triangleIds) {
        const auto found = trianglesById.find(id);
        if (found == trianglesById.end()) {
            throw std::invalid_argument(
                "Scene structure triangle mapping references an unknown ID");
        }
        definition.triangles.push_back(
            {id,
             found->second->negativeSideRegionId,
             found->second->positiveSideRegionId});
    }
    definition.lines.reserve(
        canonical.mappings.constraintSuspensionLineIds.size()
        + canonical.mappings.suspensionSegmentLineIds.size());
    for (const fsi::StableId id :
         canonical.mappings.constraintSuspensionLineIds) {
        const auto found = linesById.find(id);
        if (found == linesById.end()) {
            throw std::invalid_argument(
                "Scene structure line mapping references an unknown ID");
        }
        definition.lines.push_back(
            {id, static_cast<std::uint32_t>(found->second->role)});
    }
    if (canonical.definition.suspension) {
        const fsi::StructureSuspensionDefinition& suspension =
            *canonical.definition.suspension;
        if (suspension.segments.size()
            != canonical.mappings.suspensionSegmentLineIds.size()) {
            throw std::invalid_argument(
                "Scene suspension segment mapping is incomplete");
        }
        const auto endpointVertex = [&](const auto& endpoint) {
            if (endpoint.kind
                == fsi::StructureSuspensionEndpointKind::SurfaceAttachment) {
                const auto found = std::ranges::find_if(
                    suspension.attachments, [&](const auto& attachment) {
                        return attachment.stableId == endpoint.stableId;
                    });
                if (found != suspension.attachments.end()) {
                    return found->node;
                }
            } else if (endpoint.kind
                       == fsi::StructureSuspensionEndpointKind::Junction) {
                const auto found = std::ranges::find_if(
                    suspension.junctions, [&](const auto& junction) {
                        return junction.stableId == endpoint.stableId;
                    });
                if (found != suspension.junctions.end()) {
                    return found->node;
                }
            } else if (endpoint.kind
                       == fsi::StructureSuspensionEndpointKind::PilotHarness) {
                const std::optional<std::size_t> harness =
                    canonical.mappings.pilotHarnessIndex(endpoint.stableId);
                if (harness) {
                    return canonical.definition.nodes.size() + *harness;
                }
            }
            throw std::invalid_argument(
                "Scene suspension endpoint mapping is incomplete");
        };
        for (std::size_t index = 0; index < suspension.segments.size();
             ++index) {
            const fsi::StableId id =
                canonical.mappings.suspensionSegmentLineIds[index];
            const auto found = linesById.find(id);
            if (found == linesById.end()
                || suspension.segments[index].stableId != id) {
                throw std::invalid_argument(
                    "Scene suspension mapping references an unknown line");
            }
            const std::size_t first = endpointVertex(
                suspension.segments[index].from);
            const std::size_t second = endpointVertex(
                suspension.segments[index].to);
            definition.lines.push_back(
                {id, static_cast<std::uint32_t>(found->second->role),
                 static_cast<std::uint32_t>(first),
                 static_cast<std::uint32_t>(second)});
        }
    }
    return StructureFrameMapping(structure, std::move(definition));
}

DiagnosticFrame buildStructureFrame(
    const fsi::Structure& structure,
    const StructureFrameMapping& mapping,
    const StructureFrameContext& context) {
    if (mapping.structureFingerprint() != structure.definitionFingerprint()) {
        throw std::invalid_argument(
            "Structure frame mapping belongs to a different structure");
    }

    const fsi::StructureDefinition& definition = structure.definition();
    const fsi::StructureCheckpoint checkpoint = structure.checkpoint();
    const fsi::StructureDiagnostics diagnostics = structure.diagnostics();
    const std::optional<fsi::StructureSuspensionState> suspensionState =
        structure.suspensionState();
    const std::size_t harnessCount = suspensionState
        ? suspensionState->harnessPositionsMeters.size() : 0;
    const std::size_t suspensionLineCount = suspensionState
        ? suspensionState->segments.size() : 0;
    if (checkpoint.nodes.size() != definition.nodes.size()
        || checkpoint.pendingExternalForcesNewtons.size()
               != definition.nodes.size()
        || mapping.vertexStableIds().size()
               != definition.nodes.size() + harnessCount
        || mapping.triangles().size() != definition.triangles.size()
        || mapping.lines().size()
               != definition.constraints.size() + suspensionLineCount
        || definition.suspension.has_value()
               != suspensionState.has_value()
        || (definition.suspension
            && (definition.suspension->harnessPoints.size() != harnessCount
                || definition.suspension->segments.size()
                       != suspensionLineCount))) {
        throw std::logic_error(
            "Structure frame state and topology are inconsistent");
    }

    DiagnosticFrame frame;
    frame.sceneChecksum = context.sceneChecksum;
    frame.solverCommit = context.solverCommit;
    frame.step = checkpoint.acceptedStepCount;
    frame.simulationTimeSeconds = checkpoint.simulationTimeSeconds;
    frame.timeStepSeconds = context.timeStepSeconds;
    frame.couplingIteration = context.couplingIteration;
    frame.couplingResiduals = context.couplingResiduals;
    frame.conservation = context.conservation;

    frame.vertices.reserve(mapping.vertexStableIds().size());
    for (std::size_t index = 0; index < checkpoint.nodes.size(); ++index) {
        frame.vertices.push_back(
            {mapping.vertexStableIds()[index],
             toViewer(checkpoint.nodes[index].positionMeters)});
    }
    if (suspensionState) {
        for (std::size_t index = 0; index < harnessCount; ++index) {
            frame.vertices.push_back(
                {mapping.vertexStableIds()[definition.nodes.size() + index],
                 toViewer(
                     suspensionState->harnessPositionsMeters[index])});
        }
    }

    frame.triangles.reserve(definition.triangles.size());
    for (std::size_t index = 0; index < definition.triangles.size(); ++index) {
        const fsi::StructureTriangleDefinition& triangle =
            definition.triangles[index];
        const StructureFrameTriangleMapping& ids = mapping.triangles()[index];
        frame.triangles.push_back(
            {ids.stableId,
             static_cast<std::uint32_t>(triangle.nodes[0]),
             static_cast<std::uint32_t>(triangle.nodes[1]),
             static_cast<std::uint32_t>(triangle.nodes[2]),
             ids.negativeRegionId,
             ids.positiveRegionId});
    }

    std::vector<double> lineLengths;
    std::vector<double> lineViolations;
    frame.lines.reserve(mapping.lines().size());
    lineLengths.reserve(mapping.lines().size());
    lineViolations.reserve(mapping.lines().size());
    for (std::size_t index = 0; index < definition.constraints.size();
         ++index) {
        const fsi::StructureConstraintDefinition& constraint =
            definition.constraints[index];
        const StructureFrameLineMapping& ids = mapping.lines()[index];
        frame.lines.push_back(
            {ids.stableId,
             static_cast<std::uint32_t>(constraint.firstNode),
             static_cast<std::uint32_t>(constraint.secondNode),
             ids.role});
        const double currentLength = length(subtract(
            checkpoint.nodes[constraint.secondNode].positionMeters,
            checkpoint.nodes[constraint.firstNode].positionMeters));
        lineLengths.push_back(currentLength);
        lineViolations.push_back(
            constraintViolation(constraint, currentLength));
    }
    if (suspensionState) {
        for (std::size_t index = 0; index < suspensionLineCount; ++index) {
            const std::size_t mappingIndex =
                definition.constraints.size() + index;
            const StructureFrameLineMapping& ids =
                mapping.lines()[mappingIndex];
            const fsi::StructureSuspensionSegmentState& segment =
                suspensionState->segments[index];
            if (segment.stableId != ids.stableId) {
                throw std::logic_error(
                    "Structure suspension state and frame mapping disagree");
            }
            frame.lines.push_back(
                {ids.stableId, ids.vertex0, ids.vertex1, ids.role});
            const double currentLength = length(subtract(
                {frame.vertices[ids.vertex1].positionMetres.x,
                 frame.vertices[ids.vertex1].positionMetres.y,
                 frame.vertices[ids.vertex1].positionMetres.z},
                {frame.vertices[ids.vertex0].positionMetres.x,
                 frame.vertices[ids.vertex0].positionMetres.y,
                 frame.vertices[ids.vertex0].positionMetres.z}));
            lineLengths.push_back(currentLength);
            lineViolations.push_back(std::max(
                0.0, currentLength - segment.commandedRestLengthMeters));
        }
    }

    std::vector<Vec3d> velocities;
    std::vector<Vec3d> pendingForces;
    velocities.reserve(checkpoint.nodes.size());
    pendingForces.reserve(checkpoint.pendingExternalForcesNewtons.size());
    for (std::size_t index = 0; index < checkpoint.nodes.size(); ++index) {
        velocities.push_back(
            toViewer(checkpoint.nodes[index].velocityMetersPerSecond));
        pendingForces.push_back(
            toViewer(checkpoint.pendingExternalForcesNewtons[index]));
    }
    if (suspensionState) {
        for (std::size_t index = 0; index < harnessCount; ++index) {
            velocities.push_back(toViewer(
                suspensionState->harnessVelocitiesMetersPerSecond[index]));
            pendingForces.push_back({});
        }
    }

    frame.scalarFields.push_back(
        {"structure.constraint_length", "m", FieldAssociation::Line,
         std::move(lineLengths)});
    frame.scalarFields.push_back(
        {"structure.constraint_violation", "m", FieldAssociation::Line,
         std::move(lineViolations)});
    addGlobalScalar(frame, "structure.dynamic_mass", "kg",
                    diagnostics.totalDynamicMassKg);
    addGlobalScalar(frame, "structure.kinetic_energy", "J",
                    diagnostics.kineticEnergyJoules);
    addGlobalScalar(frame, "structure.maximum_distance_error", "m",
                    diagnostics.maximumDistanceErrorMeters);
    addGlobalScalar(frame, "structure.maximum_cable_extension", "m",
                    diagnostics.maximumCableExtensionMeters);
    addGlobalScalar(frame, "structure.maximum_absolute_membrane_strain", "1",
                    diagnostics.maximumAbsoluteMembraneStrain);
    addGlobalScalar(frame, "structure.maximum_membrane_residual", "1",
                    diagnostics.maximumMembraneResidual);
    addGlobalScalar(frame, "structure.maximum_contact_penetration", "m",
                    diagnostics.maximumContactPenetrationMeters);
    addGlobalScalar(frame, "structure.maximum_suspension_residual", "m",
                    diagnostics.maximumSuspensionResidualMeters);

    frame.vectorFields.push_back(
        {"structure.velocity", "m/s", FieldAssociation::Vertex,
         std::move(velocities)});
    frame.vectorFields.push_back(
        {"structure.pending_external_force", "N", FieldAssociation::Vertex,
         std::move(pendingForces)});
    addGlobalVector(frame, "structure.total_pending_external_force", "N",
                    diagnostics.pendingExternalForceNewtons);
    addGlobalVector(frame, "structure.last_applied_external_force", "N",
                    diagnostics.lastAppliedExternalForceNewtons);
    addGlobalVector(frame, "structure.center_of_mass", "m",
                    diagnostics.centerOfMassMeters);
    addGlobalVector(frame, "structure.linear_momentum", "kg*m/s",
                    diagnostics.linearMomentumKgMetersPerSecond);

    ProtocolError error;
    if (!validateFrame(frame, &error)) {
        throw std::invalid_argument(
            "Invalid structure diagnostic frame: " + error.message);
    }
    return frame;
}

} // namespace simwing::viewer

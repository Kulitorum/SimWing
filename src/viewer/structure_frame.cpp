#include "structure_frame.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

} // namespace

StructureFrameMapping::StructureFrameMapping(
    const fsi::Structure& structure,
    StructureFrameMappingDefinition definition)
    : structureFingerprint_(structure.definitionFingerprint()),
      vertexStableIds_(std::move(definition.vertexStableIds)),
      triangles_(std::move(definition.triangles)),
      lines_(std::move(definition.lines)) {
    const fsi::StructureDefinition& topology = structure.definition();
    if (topology.nodes.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument(
            "Structure frame has too many vertices for the viewer protocol");
    }
    if (vertexStableIds_.size() != topology.nodes.size()
        || triangles_.size() != topology.triangles.size()
        || lines_.size() != topology.constraints.size()) {
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
            || triangle.positiveRegionId == 0
            || triangle.negativeRegionId == triangle.positiveRegionId) {
            throw std::invalid_argument(
                "Structure frame triangle needs two distinct nonzero regions");
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
    if (checkpoint.nodes.size() != definition.nodes.size()
        || checkpoint.pendingExternalForcesNewtons.size()
               != definition.nodes.size()
        || mapping.vertexStableIds().size() != definition.nodes.size()
        || mapping.triangles().size() != definition.triangles.size()
        || mapping.lines().size() != definition.constraints.size()) {
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

    frame.vertices.reserve(definition.nodes.size());
    for (std::size_t index = 0; index < checkpoint.nodes.size(); ++index) {
        frame.vertices.push_back(
            {mapping.vertexStableIds()[index],
             toViewer(checkpoint.nodes[index].positionMeters)});
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
    frame.lines.reserve(definition.constraints.size());
    lineLengths.reserve(definition.constraints.size());
    lineViolations.reserve(definition.constraints.size());
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

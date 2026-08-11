#include "structure.h"

#include "structure_checkpoint_detail.h"

#include <softwing/soft_body.h>
#include <softwing/suspension.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace simwing::fsi {
namespace {

[[nodiscard]] bool finite(const StructureVector2& value) {
    return std::isfinite(value.x) && std::isfinite(value.y);
}

[[nodiscard]] bool finite(const StructureVector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

[[nodiscard]] bool finite(const StructureQuaternion& value) {
    return std::isfinite(value.w) && std::isfinite(value.x)
        && std::isfinite(value.y) && std::isfinite(value.z);
}

[[nodiscard]] softwing::Vec3 toSoftwing(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

[[nodiscard]] StructureVector3 fromSoftwing(const softwing::Vec3& value) {
    return {value.x, value.y, value.z};
}

[[nodiscard]] softwing::Quaternion toSoftwing(
    const StructureQuaternion& value) {
    return {value.w, value.x, value.y, value.z};
}

[[nodiscard]] StructureQuaternion fromSoftwing(
    const softwing::Quaternion& value) {
    return {value.w, value.x, value.y, value.z};
}

[[nodiscard]] StructureVector3 add(const StructureVector3& first,
                                   const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

[[nodiscard]] StructureVector3 scaled(const StructureVector3& value,
                                      double factor) {
    return {value.x * factor, value.y * factor, value.z * factor};
}

[[nodiscard]] double squaredLength(const StructureVector3& value) {
    return value.x * value.x + value.y * value.y + value.z * value.z;
}

[[nodiscard]] softwing::MaterialRole toSoftwing(StructureMaterialRole role) {
    switch (role) {
    case StructureMaterialRole::Bulk:
        return softwing::MaterialRole::Bulk;
    case StructureMaterialRole::Seam:
        return softwing::MaterialRole::Seam;
    case StructureMaterialRole::Reinforcement:
        return softwing::MaterialRole::Reinforcement;
    }
    throw std::invalid_argument("Unknown structure material role");
}

[[nodiscard]] softwing::OrthotropicMembraneMaterial toSoftwing(
    const StructureMembraneMaterial& material) {
    softwing::OrthotropicMembraneMaterial result;
    result.warpStiffness = material.warpStiffnessNewtonsPerMeter;
    result.weftStiffness = material.weftStiffnessNewtonsPerMeter;
    result.couplingStiffness = material.couplingStiffnessNewtonsPerMeter;
    result.shearStiffness = material.shearStiffnessNewtonsPerMeter;
    result.warpPreTension = material.warpPreTensionNewtonsPerMeter;
    result.weftPreTension = material.weftPreTensionNewtonsPerMeter;
    result.dampingTime = material.dampingSeconds;
    result.compressionStiffnessRatio = material.compressionStiffnessRatio;
    return result;
}

void requireNodeIndex(std::size_t index,
                      std::size_t nodeCount,
                      const char* entity) {
    if (index >= nodeCount) {
        throw std::invalid_argument(std::string(entity)
                                    + " references an unknown node");
    }
}

void validateDefinition(const StructureDefinition& definition) {
    for (const StructureNodeDefinition& node : definition.nodes) {
        if (!finite(node.positionMeters) || !std::isfinite(node.massKg)
            || node.massKg < 0.0 || (!node.fixed && !(node.massKg > 0.0))) {
            throw std::invalid_argument("Invalid structure node definition");
        }
    }

    for (const StructureTriangleDefinition& triangle : definition.triangles) {
        for (const std::size_t node : triangle.nodes) {
            requireNodeIndex(node, definition.nodes.size(), "Triangle");
        }
        if (triangle.nodes[0] == triangle.nodes[1]
            || triangle.nodes[1] == triangle.nodes[2]
            || triangle.nodes[2] == triangle.nodes[0]) {
            throw std::invalid_argument("Triangle repeats a node");
        }
    }

    for (const StructureConstraintDefinition& constraint :
         definition.constraints) {
        requireNodeIndex(constraint.firstNode,
                         definition.nodes.size(),
                         "Constraint");
        requireNodeIndex(constraint.secondNode,
                         definition.nodes.size(),
                         "Constraint");
        if (constraint.firstNode == constraint.secondNode
            || !std::isfinite(constraint.restLengthMeters)
            || constraint.restLengthMeters < 0.0
            || !std::isfinite(constraint.complianceMetersPerNewton)
            || constraint.complianceMetersPerNewton < 0.0) {
            throw std::invalid_argument("Invalid structure constraint");
        }
        switch (constraint.kind) {
        case StructureConstraintKind::Distance:
        case StructureConstraintKind::Cable:
        case StructureConstraintKind::SuspensionTie:
            break;
        case StructureConstraintKind::SeamStitch:
            if (constraint.restLengthMeters != 0.0
                || constraint.complianceMetersPerNewton != 0.0) {
                throw std::invalid_argument(
                    "Seam stitch must be rigid and zero-rest");
            }
            break;
        default:
            throw std::invalid_argument("Unknown structure constraint kind");
        }
    }

    for (const StructureMembraneDefinition& membrane : definition.membranes) {
        if (membrane.triangle >= definition.triangles.size()) {
            throw std::invalid_argument(
                "Membrane references an unknown triangle");
        }
        if (!std::ranges::all_of(
                membrane.materialCoordinates,
                [](const StructureVector2& value) { return finite(value); })) {
            throw std::invalid_argument(
                "Membrane material coordinates must be finite");
        }
        softwing::validateOrthotropicMembraneMaterial(
            toSoftwing(membrane.material));
        static_cast<void>(toSoftwing(membrane.role));
    }

    for (const StructureDihedralDefinition& dihedral : definition.dihedrals) {
        for (const std::size_t node : dihedral.nodes) {
            requireNodeIndex(node, definition.nodes.size(), "Dihedral");
        }
        std::array<std::size_t, 4> uniqueNodes = dihedral.nodes;
        std::ranges::sort(uniqueNodes);
        if (std::ranges::adjacent_find(uniqueNodes) != uniqueNodes.end()
            || !std::isfinite(dihedral.restAngleRadians)
            || !std::isfinite(dihedral.complianceRadiansPerNewtonMeter)
            || dihedral.complianceRadiansPerNewtonMeter < 0.0) {
            throw std::invalid_argument("Invalid dihedral definition");
        }
    }

    if (definition.fabricSelfContact) {
        const StructureFabricContactDefinition& contact =
            *definition.fabricSelfContact;
        if (definition.triangles.empty()
            || !std::isfinite(contact.halfThicknessMeters)
            || !(contact.halfThicknessMeters > 0.0)
            || !std::isfinite(contact.normalComplianceMetersPerNewton)
            || contact.normalComplianceMetersPerNewton < 0.0
            || !std::isfinite(contact.staticFriction)
            || contact.staticFriction < 0.0
            || !std::isfinite(contact.dynamicFriction)
            || contact.dynamicFriction < 0.0
            || contact.dynamicFriction > contact.staticFriction) {
            throw std::invalid_argument(
                "Invalid structure fabric self-contact definition");
        }
    }

    if (definition.suspension) {
        const StructureSuspensionDefinition& suspension =
            *definition.suspension;
        const auto requireStableIds = [](const auto& values,
                                         const char* label) {
            std::vector<std::uint64_t> ids;
            ids.reserve(values.size());
            for (const auto& value : values) {
                if (value.stableId == 0) {
                    throw std::invalid_argument(
                        std::string(label) + " stable ID must be nonzero");
                }
                ids.push_back(value.stableId);
            }
            std::ranges::sort(ids);
            if (std::ranges::adjacent_find(ids) != ids.end()) {
                throw std::invalid_argument(
                    std::string(label) + " stable IDs must be unique");
            }
        };
        if (suspension.pilotStableId == 0
            || !(suspension.pilotMassKg > 0.0)
            || !std::isfinite(suspension.pilotMassKg)
            || !finite(suspension.pilotCenterOfMassLocalMeters)
            || !finite(suspension.pilotInitialCenterOfMassWorldMeters)
            || !finite(suspension.pilotInitialLinearVelocityMetersPerSecond)
            || !finite(suspension.pilotInitialAngularVelocityRadiansPerSecond)
            || !finite(suspension.pilotInitialBodyToWorld)
            || !(suspension.pilotPrincipalInertiaKgSquareMeters.x > 0.0)
            || !(suspension.pilotPrincipalInertiaKgSquareMeters.y > 0.0)
            || !(suspension.pilotPrincipalInertiaKgSquareMeters.z > 0.0)
            || !finite(suspension.pilotPrincipalInertiaKgSquareMeters)
            || suspension.attachments.empty()
            || suspension.harnessPoints.empty()
            || suspension.segments.empty()
            || suspension.solverIterations <= 0
            || !(suspension.attachmentTolerance > 0.0)
            || !std::isfinite(suspension.attachmentTolerance)
            || !(suspension.minimumLineLengthMeters > 0.0)
            || !std::isfinite(suspension.minimumLineLengthMeters)
            || suspension.maximumLineResidualMeters < 0.0
            || !std::isfinite(suspension.maximumLineResidualMeters)
            || suspension.maximumControlWorkJoules < 0.0
            || !std::isfinite(suspension.maximumControlWorkJoules)) {
            throw std::invalid_argument(
                "Invalid structure suspension definition");
        }
        requireStableIds(suspension.attachments,
                         "Suspension attachment");
        requireStableIds(suspension.junctions, "Suspension junction");
        requireStableIds(suspension.harnessPoints, "Pilot harness");
        requireStableIds(suspension.segments, "Suspension segment");
        for (const StructureSuspensionAttachmentDefinition& attachment :
             suspension.attachments) {
            requireNodeIndex(attachment.node, definition.nodes.size(),
                             "Suspension attachment");
        }
        for (const StructureSuspensionJunctionDefinition& junction :
             suspension.junctions) {
            requireNodeIndex(junction.node, definition.nodes.size(),
                             "Suspension junction");
            if (definition.nodes[junction.node].fixed) {
                throw std::invalid_argument(
                    "Suspension junction must be a dynamic node");
            }
        }
        for (const StructurePilotHarnessDefinition& harness :
             suspension.harnessPoints) {
            if (!finite(harness.localPositionMeters)) {
                throw std::invalid_argument(
                    "Pilot harness position must be finite");
            }
        }
        for (const StructureSuspensionSegmentDefinition& segment :
             suspension.segments) {
            if (segment.from.stableId == 0 || segment.to.stableId == 0
                || (segment.from.kind == segment.to.kind
                    && segment.from.stableId == segment.to.stableId)
                || !(segment.restLengthMeters > 0.0)
                || !std::isfinite(segment.restLengthMeters)
                || !(segment.axialStiffnessNewtons > 0.0)
                || !std::isfinite(segment.axialStiffnessNewtons)
                || segment.axialDampingNewtonSecondsPerMeter < 0.0
                || !std::isfinite(
                    segment.axialDampingNewtonSecondsPerMeter)) {
                throw std::invalid_argument(
                    "Invalid structure suspension segment");
            }
        }
    }
}

class Fingerprint {
public:
    void add(std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= 1099511628211ULL;
            value >>= 8U;
        }
    }

    void add(double value) {
        add(std::bit_cast<std::uint64_t>(value));
    }

    void add(const StructureVector3& value) {
        add(value.x);
        add(value.y);
        add(value.z);
    }

    void add(const StructureQuaternion& value) {
        add(value.w);
        add(value.x);
        add(value.y);
        add(value.z);
    }

    [[nodiscard]] std::uint64_t value() const { return value_; }

private:
    std::uint64_t value_ = 14695981039346656037ULL;
};

[[nodiscard]] std::uint64_t computeDefinitionFingerprint(
    const StructureDefinition& definition) {
    Fingerprint hash;
    hash.add(std::uint64_t{2});
    hash.add(static_cast<std::uint64_t>(definition.nodes.size()));
    for (const auto& node : definition.nodes) {
        hash.add(node.positionMeters.x);
        hash.add(node.positionMeters.y);
        hash.add(node.positionMeters.z);
        hash.add(node.massKg);
        hash.add(node.fixed ? std::uint64_t{1} : std::uint64_t{0});
    }
    hash.add(static_cast<std::uint64_t>(definition.triangles.size()));
    for (const auto& triangle : definition.triangles) {
        for (const auto node : triangle.nodes)
            hash.add(static_cast<std::uint64_t>(node));
    }
    hash.add(static_cast<std::uint64_t>(definition.constraints.size()));
    for (const auto& constraint : definition.constraints) {
        hash.add(static_cast<std::uint64_t>(constraint.kind));
        hash.add(static_cast<std::uint64_t>(constraint.firstNode));
        hash.add(static_cast<std::uint64_t>(constraint.secondNode));
        hash.add(constraint.restLengthMeters);
        hash.add(constraint.complianceMetersPerNewton);
    }
    hash.add(static_cast<std::uint64_t>(definition.membranes.size()));
    for (const auto& membrane : definition.membranes) {
        hash.add(static_cast<std::uint64_t>(membrane.triangle));
        for (const auto& chart : membrane.materialCoordinates) {
            hash.add(chart.x);
            hash.add(chart.y);
        }
        const auto& material = membrane.material;
        hash.add(material.warpStiffnessNewtonsPerMeter);
        hash.add(material.weftStiffnessNewtonsPerMeter);
        hash.add(material.couplingStiffnessNewtonsPerMeter);
        hash.add(material.shearStiffnessNewtonsPerMeter);
        hash.add(material.warpPreTensionNewtonsPerMeter);
        hash.add(material.weftPreTensionNewtonsPerMeter);
        hash.add(material.dampingSeconds);
        hash.add(material.compressionStiffnessRatio);
        hash.add(static_cast<std::uint64_t>(membrane.role));
    }
    hash.add(static_cast<std::uint64_t>(definition.dihedrals.size()));
    for (const auto& dihedral : definition.dihedrals) {
        for (const auto node : dihedral.nodes)
            hash.add(static_cast<std::uint64_t>(node));
        hash.add(dihedral.restAngleRadians);
        hash.add(dihedral.complianceRadiansPerNewtonMeter);
    }
    hash.add(definition.fabricSelfContact ? std::uint64_t{1}
                                          : std::uint64_t{0});
    if (definition.fabricSelfContact) {
        const auto& contact = *definition.fabricSelfContact;
        hash.add(contact.halfThicknessMeters);
        hash.add(contact.normalComplianceMetersPerNewton);
        hash.add(contact.staticFriction);
        hash.add(contact.dynamicFriction);
    }
    hash.add(definition.suspension ? std::uint64_t{1} : std::uint64_t{0});
    if (definition.suspension) {
        const StructureSuspensionDefinition& suspension =
            *definition.suspension;
        hash.add(suspension.pilotStableId);
        hash.add(suspension.pilotMassKg);
        hash.add(suspension.pilotCenterOfMassLocalMeters);
        hash.add(suspension.pilotInitialCenterOfMassWorldMeters);
        hash.add(suspension.pilotInitialLinearVelocityMetersPerSecond);
        hash.add(suspension.pilotInitialAngularVelocityRadiansPerSecond);
        hash.add(suspension.pilotInitialBodyToWorld);
        hash.add(suspension.pilotPrincipalInertiaKgSquareMeters);
        hash.add(static_cast<std::uint64_t>(suspension.attachments.size()));
        for (const auto& attachment : suspension.attachments) {
            hash.add(attachment.stableId);
            hash.add(static_cast<std::uint64_t>(attachment.node));
        }
        hash.add(static_cast<std::uint64_t>(suspension.junctions.size()));
        for (const auto& junction : suspension.junctions) {
            hash.add(junction.stableId);
            hash.add(static_cast<std::uint64_t>(junction.node));
        }
        hash.add(static_cast<std::uint64_t>(suspension.harnessPoints.size()));
        for (const auto& harness : suspension.harnessPoints) {
            hash.add(harness.stableId);
            hash.add(harness.localPositionMeters);
        }
        hash.add(static_cast<std::uint64_t>(suspension.segments.size()));
        for (const auto& segment : suspension.segments) {
            hash.add(segment.stableId);
            hash.add(static_cast<std::uint64_t>(segment.from.kind));
            hash.add(segment.from.stableId);
            hash.add(static_cast<std::uint64_t>(segment.to.kind));
            hash.add(segment.to.stableId);
            hash.add(segment.restLengthMeters);
            hash.add(segment.axialStiffnessNewtons);
            hash.add(segment.axialDampingNewtonSecondsPerMeter);
            hash.add(static_cast<std::uint64_t>(segment.role));
        }
        hash.add(static_cast<std::uint64_t>(suspension.solverIterations));
        hash.add(suspension.attachmentTolerance);
        hash.add(suspension.minimumLineLengthMeters);
        hash.add(suspension.maximumLineResidualMeters);
        hash.add(suspension.maximumControlWorkJoules);
    }
    return hash.value();
}

[[nodiscard]] std::string suspensionId(
    StructureSuspensionEndpointKind kind,
    std::uint64_t stableId) {
    switch (kind) {
    case StructureSuspensionEndpointKind::SurfaceAttachment:
        return "attachment:" + std::to_string(stableId);
    case StructureSuspensionEndpointKind::Junction:
        return "junction:" + std::to_string(stableId);
    case StructureSuspensionEndpointKind::PilotHarness:
        return "harness:" + std::to_string(stableId);
    }
    throw std::invalid_argument("Unknown structure suspension endpoint kind");
}

[[nodiscard]] softwing::SuspensionEndpointKind toSoftwing(
    StructureSuspensionEndpointKind kind) {
    switch (kind) {
    case StructureSuspensionEndpointKind::SurfaceAttachment:
        return softwing::SuspensionEndpointKind::Attachment;
    case StructureSuspensionEndpointKind::Junction:
        return softwing::SuspensionEndpointKind::Junction;
    case StructureSuspensionEndpointKind::PilotHarness:
        return softwing::SuspensionEndpointKind::HangPoint;
    }
    throw std::invalid_argument("Unknown structure suspension endpoint kind");
}

[[nodiscard]] softwing::SuspensionSystem buildSuspension(
    softwing::SoftBody& body,
    const StructureSuspensionDefinition& source) {
    constexpr const char* provenanceId = "simwing-structure-resolved";
    softwing::SuspensionDefinition definition;
    definition.identifier = "simwing-structure";
    definition.description =
        "Resolved scene-v2 suspension and rigid pilot";
    definition.unitsFrameTag = std::string(softwing::suspensionStage5FrameTag);
    definition.provenance = {
        {provenanceId, "authoritative scene-v2 stable topology"}};
    definition.attachments.reserve(source.attachments.size());
    std::vector<softwing::ResolvedSuspensionAttachment> resolved;
    resolved.reserve(source.attachments.size());
    for (const StructureSuspensionAttachmentDefinition& attachment :
         source.attachments) {
        const std::string id = suspensionId(
            StructureSuspensionEndpointKind::SurfaceAttachment,
            attachment.stableId);
        const std::string panelId =
            "resolved-node:" + std::to_string(attachment.node);
        definition.attachments.push_back(
            {id, panelId, {}, softwing::SuspensionSide::Centre, {},
             provenanceId});
        resolved.push_back(
            {id, panelId, {}, attachment.node,
             body.nodes()[attachment.node].position, provenanceId});
    }
    definition.junctions.reserve(source.junctions.size());
    std::vector<std::pair<std::string, std::size_t>> junctionNodes;
    junctionNodes.reserve(source.junctions.size());
    for (const StructureSuspensionJunctionDefinition& junction :
         source.junctions) {
        const std::string id = suspensionId(
            StructureSuspensionEndpointKind::Junction,
            junction.stableId);
        const softwing::Node& node = body.nodes()[junction.node];
        definition.junctions.push_back(
            {id, node.position, 1.0 / node.inverseMass,
             softwing::SuspensionSide::Centre, provenanceId});
        junctionNodes.emplace_back(id, junction.node);
    }

    definition.payload.mass = source.pilotMassKg;
    definition.payload.centreOfMassLocal = toSoftwing(
        source.pilotCenterOfMassLocalMeters);
    definition.payload.inertiaBody = {{
        source.pilotPrincipalInertiaKgSquareMeters.x, 0.0, 0.0,
        0.0, source.pilotPrincipalInertiaKgSquareMeters.y, 0.0,
        0.0, 0.0, source.pilotPrincipalInertiaKgSquareMeters.z}};
    definition.payload.initialState = {
        toSoftwing(source.pilotInitialCenterOfMassWorldMeters),
        toSoftwing(source.pilotInitialBodyToWorld),
        toSoftwing(source.pilotInitialLinearVelocityMetersPerSecond),
        toSoftwing(source.pilotInitialAngularVelocityRadiansPerSecond)};
    definition.payload.provenanceId = provenanceId;
    definition.payload.hangPoints.reserve(source.harnessPoints.size());
    for (const StructurePilotHarnessDefinition& harness :
         source.harnessPoints) {
        definition.payload.hangPoints.push_back(
            {suspensionId(
                 StructureSuspensionEndpointKind::PilotHarness,
                 harness.stableId),
             toSoftwing(harness.localPositionMeters),
             softwing::SuspensionSide::Centre,
             provenanceId});
    }
    definition.segments.reserve(source.segments.size());
    for (const StructureSuspensionSegmentDefinition& segment :
         source.segments) {
        definition.segments.push_back(
            {"segment:" + std::to_string(segment.stableId),
             {toSoftwing(segment.from.kind),
              suspensionId(segment.from.kind, segment.from.stableId)},
             {toSoftwing(segment.to.kind),
              suspensionId(segment.to.kind, segment.to.stableId)},
             segment.restLengthMeters,
             segment.axialStiffnessNewtons,
             segment.axialDampingNewtonSecondsPerMeter,
             softwing::SuspensionSide::Centre,
             {"role:" + std::to_string(segment.role)},
             provenanceId});
    }
    definition.solver.lineIterations = source.solverIterations;
    definition.solver.attachmentTolerance = source.attachmentTolerance;
    definition.solver.minimumLineLength = source.minimumLineLengthMeters;
    definition.solver.maximumLineResidual =
        source.maximumLineResidualMeters;
    definition.solver.maximumControlWork = source.maximumControlWorkJoules;
    definition.ground.mode = softwing::PayloadGroundMode::Free;
    return softwing::SuspensionSystem::buildResolved(
        body, definition, resolved, junctionNodes);
}

[[nodiscard]] softwing::SoftBody buildBody(
    const StructureDefinition& definition) {
    softwing::SoftBody body;
    for (const StructureNodeDefinition& node : definition.nodes) {
        if (node.fixed) {
            body.addFixedNode(toSoftwing(node.positionMeters));
        } else {
            body.addNode(toSoftwing(node.positionMeters), node.massKg);
        }
    }
    for (const StructureTriangleDefinition& triangle : definition.triangles) {
        body.addTriangle(
            triangle.nodes[0], triangle.nodes[1], triangle.nodes[2]);
    }
    for (const StructureConstraintDefinition& constraint :
         definition.constraints) {
        switch (constraint.kind) {
        case StructureConstraintKind::Distance:
            body.addDistanceConstraint(
                constraint.firstNode,
                constraint.secondNode,
                constraint.restLengthMeters,
                constraint.complianceMetersPerNewton);
            break;
        case StructureConstraintKind::Cable:
            body.addCableConstraint(
                constraint.firstNode,
                constraint.secondNode,
                constraint.restLengthMeters,
                constraint.complianceMetersPerNewton);
            break;
        case StructureConstraintKind::SuspensionTie:
            body.addSuspensionTieConstraint(
                constraint.firstNode,
                constraint.secondNode,
                constraint.restLengthMeters,
                constraint.complianceMetersPerNewton);
            break;
        case StructureConstraintKind::SeamStitch:
            body.addSeamStitchConstraint(
                constraint.firstNode,
                constraint.secondNode);
            break;
        }
    }
    if (!definition.membranes.empty()) {
        std::vector<softwing::MembraneElementDefinition> membranes;
        membranes.reserve(definition.membranes.size());
        for (const StructureMembraneDefinition& membrane :
             definition.membranes) {
            softwing::MembraneElementDefinition converted;
            converted.triangle = membrane.triangle;
            for (std::size_t corner = 0; corner < 3; ++corner) {
                converted.chart[corner] = {
                    membrane.materialCoordinates[corner].x,
                    membrane.materialCoordinates[corner].y};
            }
            converted.material = toSoftwing(membrane.material);
            converted.role = toSoftwing(membrane.role);
            membranes.push_back(converted);
        }
        static_cast<void>(body.addMembraneElements(membranes));
    }
    for (const StructureDihedralDefinition& dihedral : definition.dihedrals) {
        body.addDihedralBendingConstraint(
            dihedral.nodes[0],
            dihedral.nodes[1],
            dihedral.nodes[2],
            dihedral.nodes[3],
            dihedral.restAngleRadians,
            dihedral.complianceRadiansPerNewtonMeter);
    }
    if (definition.fabricSelfContact) {
        const StructureFabricContactDefinition& contact =
            *definition.fabricSelfContact;
        const softwing::ContactSurfaceHandle surface = body.addContactSurface(
            body.surfaceGroup(0, definition.triangles.size()),
            contact.halfThicknessMeters);
        static_cast<void>(body.addContactPair(
            surface.collider(), surface.collider(),
            {contact.normalComplianceMetersPerNewton,
             contact.staticFriction,
             contact.dynamicFriction}));
    }
    return body;
}

[[nodiscard]] softwing::StepSettings toSoftwing(
    const StructureStepSettings& settings) {
    softwing::StepSettings converted;
    converted.timeStep = settings.timeStepSeconds;
    converted.substeps = settings.substeps;
    converted.constraintIterations = settings.constraintIterations;
    converted.cableConstraintSweepPairs =
        settings.cableConstraintSweepPairs;
    converted.gravity = toSoftwing(settings.gravityMetersPerSecondSquared);
    converted.velocityDampingPerSecond =
        settings.velocityDampingPerSecond;
    converted.dampingReferenceVelocity = toSoftwing(
        settings.dampingReferenceVelocityMetersPerSecond);
    converted.workerThreads = settings.workerThreads;
    return converted;
}

} // namespace

struct Structure::Impl {
    explicit Impl(StructureDefinition value)
        : definition(std::move(value)),
          fingerprint(computeDefinitionFingerprint(definition)),
          body(buildBody(definition)),
          pendingForces(definition.nodes.size()) {
        if (definition.suspension) {
            suspension = std::make_unique<softwing::SuspensionSystem>(
                buildSuspension(body, *definition.suspension));
        }
    }

    StructureDefinition definition;
    std::uint64_t fingerprint = 0;
    softwing::SoftBody body;
    std::unique_ptr<softwing::SuspensionSystem> suspension;
    std::vector<StructureVector3> pendingForces;
    StructureVector3 lastAppliedForce;
    std::uint64_t acceptedSteps = 0;
    double simulationTime = 0.0;
};

Structure::Structure(StructureDefinition definition) {
    validateDefinition(definition);
    impl_ = std::make_unique<Impl>(std::move(definition));
}

Structure::~Structure() = default;
Structure::Structure(Structure&&) noexcept = default;
Structure& Structure::operator=(Structure&&) noexcept = default;

const StructureDefinition& Structure::definition() const noexcept {
    return impl_->definition;
}

std::uint64_t Structure::definitionFingerprint() const noexcept {
    return impl_->fingerprint;
}

std::uint64_t Structure::acceptedStepCount() const noexcept {
    return impl_->acceptedSteps;
}

double Structure::simulationTimeSeconds() const noexcept {
    return impl_->simulationTime;
}

std::vector<StructureNodeState> Structure::nodeStates() const {
    std::vector<StructureNodeState> result;
    result.reserve(impl_->body.nodes().size());
    for (const softwing::Node& node : impl_->body.nodes()) {
        result.push_back({fromSoftwing(node.position),
                          fromSoftwing(node.previousPosition),
                          fromSoftwing(node.velocity)});
    }
    return result;
}

std::optional<StructureSuspensionState> Structure::suspensionState() const {
    if (!impl_->suspension || !impl_->definition.suspension) {
        return std::nullopt;
    }
    const StructureSuspensionDefinition& definition =
        *impl_->definition.suspension;
    const softwing::SuspensionSystem& suspension = *impl_->suspension;
    const softwing::RigidPayloadState& payload = suspension.payloadState();
    StructureSuspensionState result;
    result.payload = {
        fromSoftwing(payload.centreOfMassWorld),
        fromSoftwing(payload.orientation),
        fromSoftwing(payload.linearVelocity),
        fromSoftwing(payload.angularVelocity)};
    result.harnessPositionsMeters.reserve(definition.harnessPoints.size());
    result.harnessVelocitiesMetersPerSecond.reserve(
        definition.harnessPoints.size());
    for (const StructurePilotHarnessDefinition& harness :
         definition.harnessPoints) {
        const std::string id = suspensionId(
            StructureSuspensionEndpointKind::PilotHarness,
            harness.stableId);
        const auto found = std::ranges::find_if(
            suspension.definition().payload.hangPoints,
            [&](const softwing::PayloadPointDefinition& point) {
                return point.id == id;
            });
        if (found == suspension.definition().payload.hangPoints.end()) {
            throw std::logic_error(
                "Structure suspension lost a pilot harness mapping");
        }
        const std::size_t index = static_cast<std::size_t>(
            found - suspension.definition().payload.hangPoints.begin());
        result.harnessPositionsMeters.push_back(
            fromSoftwing(suspension.hangPointPosition(index)));
        result.harnessVelocitiesMetersPerSecond.push_back(
            fromSoftwing(suspension.hangPointVelocity(index)));
    }
    result.segments.reserve(definition.segments.size());
    for (const StructureSuspensionSegmentDefinition& segment :
         definition.segments) {
        const std::string id = "segment:" + std::to_string(segment.stableId);
        const auto found = std::ranges::find_if(
            suspension.segmentDiagnostics(),
            [&](const softwing::SuspensionSegmentDiagnostics& diagnostic) {
                return diagnostic.id == id;
            });
        if (found == suspension.segmentDiagnostics().end()) {
            result.segments.push_back(
                {segment.stableId, segment.restLengthMeters,
                 segment.restLengthMeters, 0.0});
        } else {
            result.segments.push_back(
                {segment.stableId, found->length,
                 found->commandedRestLength, found->tension});
        }
    }
    return result;
}

void Structure::clearExternalForces() noexcept {
    std::ranges::fill(impl_->pendingForces, StructureVector3{});
}

void Structure::addExternalForce(std::size_t node,
                                 const StructureVector3& forceNewtons) {
    if (node >= impl_->pendingForces.size()) {
        throw std::out_of_range("Structure force node is out of range");
    }
    if (!finite(forceNewtons)) {
        throw std::invalid_argument("Structure force must be finite");
    }
    impl_->pendingForces[node] = add(impl_->pendingForces[node], forceNewtons);
}

void Structure::setExternalForces(
    std::span<const StructureVector3> forcesNewtons) {
    if (forcesNewtons.size() != impl_->pendingForces.size()) {
        throw std::invalid_argument(
            "Structure force count must equal the node count");
    }
    if (!std::ranges::all_of(
            forcesNewtons,
            [](const StructureVector3& value) { return finite(value); })) {
        throw std::invalid_argument("Structure forces must be finite");
    }
    std::ranges::copy(forcesNewtons, impl_->pendingForces.begin());
}

StructureDiagnostics Structure::step(const StructureStepSettings& settings) {
    const StructureCheckpoint before = checkpoint();
    StructureVector3 applied;
    try {
        impl_->body.clearExternalForces();
        for (std::size_t node = 0; node < impl_->pendingForces.size(); ++node) {
            impl_->body.addForce(node, toSoftwing(impl_->pendingForces[node]));
            applied = add(applied, impl_->pendingForces[node]);
        }
        if (impl_->suspension) {
            impl_->body.stepCoupled(
                toSoftwing(settings), *impl_->suspension);
        } else {
            impl_->body.step(toSoftwing(settings));
        }
        clearExternalForces();
        impl_->lastAppliedForce = applied;
        ++impl_->acceptedSteps;
        impl_->simulationTime += settings.timeStepSeconds;
    } catch (...) {
        restore(before);
        throw;
    }
    return diagnostics();
}

StructureDiagnostics Structure::diagnostics() const {
    StructureDiagnostics result;
    result.nodeCount = impl_->body.nodes().size();
    result.triangleCount = impl_->body.triangles().size();
    result.constraintCount = impl_->body.constraints().size();
    result.membraneCount = impl_->body.membraneElements().size();
    result.dihedralCount = impl_->body.dihedralConstraints().size();
    if (impl_->definition.fabricSelfContact) {
        result.contactPairCount = 1;
        const softwing::ContactDiagnostics& contact =
            impl_->body.contactDiagnostics();
        result.activeContactCount = contact.activeCount;
        result.maximumContactPenetrationMeters = contact.maximumPenetration;
        result.finite = result.finite && contact.solveSucceeded
            && !contact.hasFailure
            && std::isfinite(contact.maximumPenetration);
    }

    StructureVector3 weightedPosition;
    for (const softwing::Node& node : impl_->body.nodes()) {
        const StructureVector3 position = fromSoftwing(node.position);
        const StructureVector3 previous = fromSoftwing(node.previousPosition);
        const StructureVector3 velocity = fromSoftwing(node.velocity);
        result.finite = result.finite && finite(position) && finite(previous)
            && finite(velocity) && std::isfinite(node.inverseMass)
            && node.inverseMass >= 0.0;
        if (node.inverseMass > 0.0) {
            const double mass = 1.0 / node.inverseMass;
            ++result.dynamicNodeCount;
            result.totalDynamicMassKg += mass;
            weightedPosition = add(weightedPosition, scaled(position, mass));
            result.linearMomentumKgMetersPerSecond = add(
                result.linearMomentumKgMetersPerSecond,
                scaled(velocity, mass));
            result.kineticEnergyJoules +=
                0.5 * mass * squaredLength(velocity);
        }
    }
    if (result.totalDynamicMassKg > 0.0) {
        result.centerOfMassMeters = scaled(
            weightedPosition, 1.0 / result.totalDynamicMassKg);
    }

    if (impl_->suspension) {
        const softwing::SuspensionSystem& suspension = *impl_->suspension;
        const softwing::RigidPayloadDefinition& payloadDefinition =
            suspension.definition().payload;
        const softwing::RigidPayloadState& payload =
            suspension.payloadState();
        const double mass = payloadDefinition.mass;
        const StructureVector3 position = fromSoftwing(
            payload.centreOfMassWorld);
        const StructureVector3 momentum = fromSoftwing(
            softwing::payloadLinearMomentum(payloadDefinition, payload));
        weightedPosition = add(weightedPosition, scaled(position, mass));
        result.totalDynamicMassKg += mass;
        result.centerOfMassMeters = scaled(
            weightedPosition, 1.0 / result.totalDynamicMassKg);
        result.linearMomentumKgMetersPerSecond = add(
            result.linearMomentumKgMetersPerSecond, momentum);
        result.kineticEnergyJoules += softwing::payloadKineticEnergy(
            payloadDefinition, payload);
        const softwing::SuspensionDiagnostics& suspensionDiagnostics =
            suspension.diagnostics();
        result.suspensionSegmentCount =
            suspensionDiagnostics.segmentCount;
        result.maximumSuspensionResidualMeters =
            suspensionDiagnostics.maximumResidual;
        result.finite = result.finite
            && finite(position) && finite(momentum)
            && finite(fromSoftwing(payload.angularVelocity))
            && finite(fromSoftwing(payload.linearVelocity))
            && finite(fromSoftwing(payload.orientation))
            && std::isfinite(result.kineticEnergyJoules)
            && std::isfinite(result.maximumSuspensionResidualMeters);
    }

    for (const softwing::DistanceConstraint& constraint :
         impl_->body.constraints()) {
        const double currentLength = softwing::length(
            impl_->body.nodes()[constraint.b].position
            - impl_->body.nodes()[constraint.a].position);
        const double error = currentLength - constraint.restLength;
        if (constraint.kind == softwing::ConstraintKind::Cable) {
            result.maximumCableExtensionMeters = std::max(
                result.maximumCableExtensionMeters, std::max(0.0, error));
        } else {
            result.maximumDistanceErrorMeters = std::max(
                result.maximumDistanceErrorMeters, std::abs(error));
        }
    }
    for (std::size_t membrane = 0;
         membrane < impl_->body.membraneElements().size();
         ++membrane) {
        const softwing::MembraneElementDiagnostics value =
            impl_->body.membraneDiagnostics(membrane);
        result.maximumAbsoluteMembraneStrain = std::max(
            result.maximumAbsoluteMembraneStrain,
            std::max({std::abs(value.greenStrain.x),
                      std::abs(value.greenStrain.y),
                      std::abs(value.greenStrain.z)}));
        result.maximumMembraneResidual = std::max(
            result.maximumMembraneResidual,
            std::abs(value.normalizedResidual));
        result.finite = result.finite
            && std::isfinite(value.elasticEnergy)
            && std::isfinite(value.normalizedResidual);
    }
    for (const StructureVector3& force : impl_->pendingForces) {
        result.pendingExternalForceNewtons = add(
            result.pendingExternalForceNewtons, force);
        result.finite = result.finite && finite(force);
    }
    result.lastAppliedExternalForceNewtons = impl_->lastAppliedForce;
    result.finite = result.finite && finite(result.centerOfMassMeters)
        && finite(result.linearMomentumKgMetersPerSecond)
        && std::isfinite(result.kineticEnergyJoules)
        && std::isfinite(result.maximumDistanceErrorMeters)
        && std::isfinite(result.maximumCableExtensionMeters)
        && std::isfinite(result.maximumAbsoluteMembraneStrain)
        && std::isfinite(result.maximumMembraneResidual)
        && std::isfinite(result.maximumContactPenetrationMeters)
        && std::isfinite(result.maximumSuspensionResidualMeters)
        && finite(result.pendingExternalForceNewtons)
        && finite(result.lastAppliedExternalForceNewtons);
    return result;
}

StructureCheckpoint Structure::checkpoint() const {
    StructureCheckpoint result;
    result.definitionFingerprint = impl_->fingerprint;
    result.acceptedStepCount = impl_->acceptedSteps;
    result.simulationTimeSeconds = impl_->simulationTime;
    result.nodes = nodeStates();
    result.pendingExternalForcesNewtons = impl_->pendingForces;
    result.lastAppliedExternalForceNewtons = impl_->lastAppliedForce;
    auto detail = std::make_shared<StructureCheckpoint::Detail>();
    detail->body = impl_->body.checkpoint();
    if (impl_->suspension) {
        detail->suspension = impl_->suspension->checkpoint();
    }
    detail->publicNodes = result.nodes;
    result.detail = std::move(detail);
    return result;
}

void Structure::restore(const StructureCheckpoint& checkpointValue) {
    if (checkpointValue.version != structureCheckpointVersion) {
        throw std::invalid_argument("Unsupported structure checkpoint version");
    }
    if (checkpointValue.definitionFingerprint != impl_->fingerprint) {
        throw std::invalid_argument(
            "Structure checkpoint belongs to a different definition");
    }
    if (checkpointValue.nodes.size() != impl_->definition.nodes.size()
        || checkpointValue.pendingExternalForcesNewtons.size()
            != impl_->definition.nodes.size()
        || !std::isfinite(checkpointValue.simulationTimeSeconds)
        || checkpointValue.simulationTimeSeconds < 0.0
        || !finite(checkpointValue.lastAppliedExternalForceNewtons)
        || !checkpointValue.detail
        || !checkpointValue.detail->body.valid()
        || checkpointValue.detail->publicNodes != checkpointValue.nodes
        || checkpointValue.detail->suspension.has_value()
            != impl_->suspension.operator bool()
        || !std::ranges::all_of(
            checkpointValue.pendingExternalForcesNewtons,
            [](const StructureVector3& value) { return finite(value); })) {
        throw std::invalid_argument("Invalid structure checkpoint");
    }
    for (const StructureNodeState& node : checkpointValue.nodes) {
        if (!finite(node.positionMeters)
            || !finite(node.previousPositionMeters)
            || !finite(node.velocityMetersPerSecond)) {
            throw std::invalid_argument("Structure checkpoint is non-finite");
        }
    }

    // Rebuild and restore off to the side so a failed suspension or core
    // restore cannot expose a partially restored coupled adapter.
    auto restored = std::make_unique<Impl>(impl_->definition);
    restored->body.restore(checkpointValue.detail->body);
    if (restored->suspension) {
        restored->suspension->restore(
            *checkpointValue.detail->suspension);
    }
    restored->pendingForces = checkpointValue.pendingExternalForcesNewtons;
    restored->lastAppliedForce =
        checkpointValue.lastAppliedExternalForceNewtons;
    restored->acceptedSteps = checkpointValue.acceptedStepCount;
    restored->simulationTime = checkpointValue.simulationTimeSeconds;
    impl_.swap(restored);
}

} // namespace simwing::fsi

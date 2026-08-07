#include "structure.h"

#include <softwing/soft_body.h>

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

[[nodiscard]] softwing::Vec3 toSoftwing(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

[[nodiscard]] StructureVector3 fromSoftwing(const softwing::Vec3& value) {
    return {value.x, value.y, value.z};
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

    [[nodiscard]] std::uint64_t value() const { return value_; }

private:
    std::uint64_t value_ = 14695981039346656037ULL;
};

[[nodiscard]] std::uint64_t computeDefinitionFingerprint(
    const StructureDefinition& definition) {
    Fingerprint hash;
    hash.add(std::uint64_t{1});
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
    return hash.value();
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
          pendingForces(definition.nodes.size()) {}

    StructureDefinition definition;
    std::uint64_t fingerprint = 0;
    softwing::SoftBody body;
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
        impl_->body.step(toSoftwing(settings));
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

    // Reassembly is the public, topology-safe way to clear core-private
    // per-iteration multipliers. They are scratch for the wrapped primitives
    // and are reset before every subsequent substep.
    softwing::SoftBody restored = buildBody(impl_->definition);
    for (std::size_t index = 0; index < checkpointValue.nodes.size(); ++index) {
        softwing::Node& target = restored.nodes()[index];
        const StructureNodeState& source = checkpointValue.nodes[index];
        target.position = toSoftwing(source.positionMeters);
        target.previousPosition = toSoftwing(source.previousPositionMeters);
        target.velocity = toSoftwing(source.velocityMetersPerSecond);
        target.force = {};
    }
    impl_->body = std::move(restored);
    impl_->pendingForces = checkpointValue.pendingExternalForcesNewtons;
    impl_->lastAppliedForce = checkpointValue.lastAppliedExternalForceNewtons;
    impl_->acceptedSteps = checkpointValue.acceptedStepCount;
    impl_->simulationTime = checkpointValue.simulationTimeSeconds;
}

} // namespace simwing::fsi

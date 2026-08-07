#include "transfer.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

StructureVector3 add(const StructureVector3& first,
                     const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

StructureVector3 subtract(const StructureVector3& first,
                          const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

StructureVector3 scale(const StructureVector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double dot(const StructureVector3& first, const StructureVector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

StructureVector3 cross(const StructureVector3& first,
                       const StructureVector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double norm(const StructureVector3& value) {
    return std::hypot(value.x, value.y, value.z);
}

template<typename Value>
std::array<Value, 3> canonicalOrientedNodes(
    const std::array<Value, 3>& nodes) {
    const std::array<Value, 3> second{
        nodes[1], nodes[2], nodes[0]};
    const std::array<Value, 3> third{
        nodes[2], nodes[0], nodes[1]};
    return std::min({nodes, second, third});
}

class FingerprintBuilder final {
public:
    void add(const std::uint64_t value) noexcept {
        for (unsigned shift = 0; shift < 64; shift += 8) {
            value_ ^= static_cast<std::uint8_t>(value >> shift);
            value_ *= 1099511628211ULL;
        }
    }

    [[nodiscard]] std::uint64_t finish() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = 14695981039346656037ULL;
};

} // namespace

std::uint64_t ConservativeTransferResult::surfaceFingerprint() const noexcept {
    return surfaceFingerprint_;
}

std::uint64_t
ConservativeTransferResult::targetDefinitionFingerprint() const noexcept {
    return targetDefinitionFingerprint_;
}

std::span<const CouplingNodeLoad>
ConservativeTransferResult::nodeLoads() const noexcept {
    return nodeLoads_;
}

const ConservativeTransferDiagnostics&
ConservativeTransferResult::diagnostics() const noexcept {
    return diagnostics_;
}

ConservativeSurfaceTransfer::ConservativeSurfaceTransfer(
    const Structure& target,
    std::vector<CouplingSurfaceNodeDefinition> nodes,
    std::vector<CouplingSurfaceTriangleDefinition> triangles)
    : targetDefinitionFingerprint_(target.definitionFingerprint()),
      nodes_(std::move(nodes)),
      triangles_(std::move(triangles)) {
    if (nodes_.empty() || triangles_.empty()) {
        throw std::invalid_argument(
            "a coupling surface requires nodes and triangles");
    }
    std::sort(nodes_.begin(), nodes_.end(), [](const auto& first, const auto& second) {
        return first.stableId < second.stableId;
    });
    std::map<std::uint64_t, std::size_t> nodeIndices;
    std::set<std::size_t> targetNodes;
    const auto& targetDefinition = target.definition();
    for (std::size_t index = 0; index < nodes_.size(); ++index) {
        const auto& node = nodes_[index];
        if (node.stableId == 0) {
            throw std::invalid_argument(
                "coupling surface node stable IDs must be nonzero");
        }
        if (node.structureNode >= targetDefinition.nodes.size()) {
            throw std::invalid_argument(
                "coupling surface node target is out of range");
        }
        if (!nodeIndices.emplace(node.stableId, index).second) {
            throw std::invalid_argument(
                "coupling surface node stable IDs must be unique");
        }
        if (!targetNodes.insert(node.structureNode).second) {
            throw std::invalid_argument(
                "coupling surface nodes must map one-to-one to structural nodes");
        }
    }

    for (auto& triangle : triangles_) {
        triangle.nodeStableIds = canonicalOrientedNodes(
            triangle.nodeStableIds);
    }
    std::sort(triangles_.begin(), triangles_.end(),
              [](const auto& first, const auto& second) {
                  return first.stableId < second.stableId;
              });
    std::set<std::uint64_t> triangleStableIds;
    std::set<std::size_t> claimedTargetTriangles;
    std::map<std::array<std::size_t, 3>, std::size_t>
        targetTrianglesByNodes;
    for (std::size_t index = 0;
         index < targetDefinition.triangles.size(); ++index) {
        const auto key = canonicalOrientedNodes(
            targetDefinition.triangles[index].nodes);
        if (!targetTrianglesByNodes.emplace(key, index).second) {
            throw std::invalid_argument(
                "target structure contains duplicate oriented triangles");
        }
    }
    std::vector<bool> usedNodes(nodes_.size(), false);
    triangleNodeIndices_.reserve(triangles_.size());
    for (const auto& triangle : triangles_) {
        if (triangle.stableId == 0
            || !triangleStableIds.insert(triangle.stableId).second) {
            throw std::invalid_argument(
                "coupling surface triangle stable IDs must be nonzero and unique");
        }
        std::array<std::size_t, 3> couplingNodes{};
        std::array<std::size_t, 3> structureNodes{};
        for (std::size_t corner = 0; corner < 3; ++corner) {
            const auto found = nodeIndices.find(triangle.nodeStableIds[corner]);
            if (found == nodeIndices.end()) {
                throw std::invalid_argument(
                    "coupling triangle references an unknown surface node");
            }
            couplingNodes[corner] = found->second;
            structureNodes[corner] = nodes_[found->second].structureNode;
            usedNodes[found->second] = true;
        }
        if (couplingNodes[0] == couplingNodes[1]
            || couplingNodes[1] == couplingNodes[2]
            || couplingNodes[2] == couplingNodes[0]) {
            throw std::invalid_argument(
                "coupling triangle must reference three distinct nodes");
        }
        const auto matchingTargetTriangle = targetTrianglesByNodes.find(
            canonicalOrientedNodes(structureNodes));
        if (matchingTargetTriangle == targetTrianglesByNodes.end()) {
            throw std::invalid_argument(
                "coupling triangle does not match an oriented structural triangle");
        }
        if (!claimedTargetTriangles.insert(matchingTargetTriangle->second).second) {
            throw std::invalid_argument(
                "a structural triangle cannot appear twice in one coupling surface");
        }
        triangleNodeIndices_.push_back(couplingNodes);
    }
    if (std::ranges::find(usedNodes, false) != usedNodes.end()) {
        throw std::invalid_argument(
            "every coupling surface node must belong to a coupling triangle");
    }

    FingerprintBuilder fingerprint;
    fingerprint.add(couplingSurfaceTopologyVersion);
    fingerprint.add(targetDefinitionFingerprint_);
    fingerprint.add(static_cast<std::uint64_t>(nodes_.size()));
    for (const auto& node : nodes_) {
        fingerprint.add(node.stableId);
        fingerprint.add(static_cast<std::uint64_t>(node.structureNode));
    }
    fingerprint.add(static_cast<std::uint64_t>(triangles_.size()));
    for (const auto& triangle : triangles_) {
        fingerprint.add(triangle.stableId);
        for (const auto nodeStableId : triangle.nodeStableIds) {
            fingerprint.add(nodeStableId);
        }
    }
    fingerprint_ = fingerprint.finish();
}

std::uint64_t ConservativeSurfaceTransfer::fingerprint() const noexcept {
    return fingerprint_;
}

std::uint64_t
ConservativeSurfaceTransfer::targetDefinitionFingerprint() const noexcept {
    return targetDefinitionFingerprint_;
}

std::span<const CouplingSurfaceNodeDefinition>
ConservativeSurfaceTransfer::nodes() const noexcept {
    return nodes_;
}

std::span<const CouplingSurfaceTriangleDefinition>
ConservativeSurfaceTransfer::triangles() const noexcept {
    return triangles_;
}

std::vector<CouplingNodeKinematics>
ConservativeSurfaceTransfer::captureKinematics(const Structure& target) const {
    if (target.definitionFingerprint() != targetDefinitionFingerprint_) {
        throw std::invalid_argument(
            "coupling surface cannot capture a foreign structure definition");
    }
    const auto states = target.nodeStates();
    std::vector<CouplingNodeKinematics> result;
    result.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        const auto& state = states[node.structureNode];
        result.push_back(
            {node.stableId, state.positionMeters,
             state.velocityMetersPerSecond});
    }
    return result;
}

ConservativeTransferResult ConservativeSurfaceTransfer::evaluate(
    const std::span<const CouplingNodeKinematics> nodeKinematics,
    const std::span<const CouplingTriangleTraction> triangleTractions,
    const ConservativeTransferSettings& settings) const {
    if (nodeKinematics.size() != nodes_.size()
        || triangleTractions.size() != triangles_.size()) {
        throw std::invalid_argument(
            "coupling kinematics and traction counts must match the surface");
    }
    if (!finite(settings.momentReferenceMeters)
        || !std::isfinite(settings.minimumTriangleAreaSquareMeters)
        || !(settings.minimumTriangleAreaSquareMeters > 0.0)) {
        throw std::invalid_argument(
            "transfer settings must have finite reference and positive minimum area");
    }
    for (std::size_t index = 0; index < nodes_.size(); ++index) {
        const auto& kinematics = nodeKinematics[index];
        if (kinematics.stableId != nodes_[index].stableId
            || !finite(kinematics.positionMeters)
            || !finite(kinematics.velocityMetersPerSecond)) {
            throw std::invalid_argument(
                "coupling node kinematics are non-finite or out of canonical order");
        }
    }
    for (std::size_t index = 0; index < triangles_.size(); ++index) {
        if (triangleTractions[index].stableId != triangles_[index].stableId
            || !finite(triangleTractions[index].tractionPascals)) {
            throw std::invalid_argument(
                "triangle tractions are non-finite or out of canonical order");
        }
    }

    ConservativeTransferResult result;
    result.surfaceFingerprint_ = fingerprint_;
    result.targetDefinitionFingerprint_ = targetDefinitionFingerprint_;
    result.nodeLoads_.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        result.nodeLoads_.push_back({node.stableId, node.structureNode, {}});
    }
    auto& diagnostics = result.diagnostics_;
    diagnostics.nodeCount = nodes_.size();
    diagnostics.triangleCount = triangles_.size();
    diagnostics.momentReferenceMeters = settings.momentReferenceMeters;

    for (std::size_t triangleIndex = 0;
         triangleIndex < triangles_.size(); ++triangleIndex) {
        const auto indices = triangleNodeIndices_[triangleIndex];
        const auto& first = nodeKinematics[indices[0]];
        const auto& second = nodeKinematics[indices[1]];
        const auto& third = nodeKinematics[indices[2]];
        const auto twiceAreaVector = cross(
            subtract(second.positionMeters, first.positionMeters),
            subtract(third.positionMeters, first.positionMeters));
        const double area = 0.5 * norm(twiceAreaVector);
        if (!std::isfinite(area)
            || !(area >= settings.minimumTriangleAreaSquareMeters)) {
            throw std::invalid_argument(
                "current coupling triangle is degenerate or below minimum area");
        }
        const auto force = scale(
            triangleTractions[triangleIndex].tractionPascals, area);
        const auto centroid = scale(
            add(add(first.positionMeters, second.positionMeters),
                third.positionMeters),
            1.0 / 3.0);
        const auto averageVelocity = scale(
            add(add(first.velocityMetersPerSecond,
                    second.velocityMetersPerSecond),
                third.velocityMetersPerSecond),
            1.0 / 3.0);
        diagnostics.surfaceAreaSquareMeters += area;
        diagnostics.integratedSurfaceForceNewtons = add(
            diagnostics.integratedSurfaceForceNewtons, force);
        diagnostics.integratedSurfaceMomentNewtonMeters = add(
            diagnostics.integratedSurfaceMomentNewtonMeters,
            cross(subtract(centroid, settings.momentReferenceMeters), force));
        diagnostics.integratedSurfacePowerWatts +=
            dot(force, averageVelocity);
        const auto cornerForce = scale(force, 1.0 / 3.0);
        for (const auto nodeIndex : indices) {
            result.nodeLoads_[nodeIndex].forceNewtons = add(
                result.nodeLoads_[nodeIndex].forceNewtons, cornerForce);
        }
    }

    for (std::size_t index = 0; index < nodes_.size(); ++index) {
        const auto& force = result.nodeLoads_[index].forceNewtons;
        diagnostics.transferredNodalForceNewtons = add(
            diagnostics.transferredNodalForceNewtons, force);
        diagnostics.transferredNodalMomentNewtonMeters = add(
            diagnostics.transferredNodalMomentNewtonMeters,
            cross(subtract(nodeKinematics[index].positionMeters,
                           settings.momentReferenceMeters),
                  force));
        diagnostics.transferredNodalPowerWatts += dot(
            force, nodeKinematics[index].velocityMetersPerSecond);
    }
    diagnostics.forceResidualNewtons = subtract(
        diagnostics.transferredNodalForceNewtons,
        diagnostics.integratedSurfaceForceNewtons);
    diagnostics.forceResidualNormNewtons = norm(
        diagnostics.forceResidualNewtons);
    diagnostics.momentResidualNewtonMeters = subtract(
        diagnostics.transferredNodalMomentNewtonMeters,
        diagnostics.integratedSurfaceMomentNewtonMeters);
    diagnostics.momentResidualNormNewtonMeters = norm(
        diagnostics.momentResidualNewtonMeters);
    diagnostics.powerResidualWatts =
        diagnostics.transferredNodalPowerWatts
        - diagnostics.integratedSurfacePowerWatts;
    diagnostics.finite =
        std::isfinite(diagnostics.surfaceAreaSquareMeters)
        && finite(diagnostics.integratedSurfaceForceNewtons)
        && finite(diagnostics.transferredNodalForceNewtons)
        && finite(diagnostics.forceResidualNewtons)
        && std::isfinite(diagnostics.forceResidualNormNewtons)
        && finite(diagnostics.integratedSurfaceMomentNewtonMeters)
        && finite(diagnostics.transferredNodalMomentNewtonMeters)
        && finite(diagnostics.momentResidualNewtonMeters)
        && std::isfinite(diagnostics.momentResidualNormNewtonMeters)
        && std::isfinite(diagnostics.integratedSurfacePowerWatts)
        && std::isfinite(diagnostics.transferredNodalPowerWatts)
        && std::isfinite(diagnostics.powerResidualWatts);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "surface traction transfer produced non-finite ledgers");
    }
    return result;
}

void ConservativeSurfaceTransfer::addLoadsTo(
    Structure& target,
    const ConservativeTransferResult& result) const {
    if (target.definitionFingerprint() != targetDefinitionFingerprint_
        || result.targetDefinitionFingerprint_ != targetDefinitionFingerprint_
        || result.surfaceFingerprint_ != fingerprint_
        || result.nodeLoads_.size() != nodes_.size()) {
        throw std::invalid_argument(
            "transfer result does not belong to this coupling surface and structure");
    }
    for (std::size_t index = 0; index < nodes_.size(); ++index) {
        const auto& load = result.nodeLoads_[index];
        if (load.stableId != nodes_[index].stableId
            || load.structureNode != nodes_[index].structureNode
            || !finite(load.forceNewtons)) {
            throw std::invalid_argument(
                "transfer result nodal loads do not match the coupling surface");
        }
    }
    for (const auto& load : result.nodeLoads_) {
        target.addExternalForce(load.structureNode, load.forceNewtons);
    }
}

} // namespace simwing::fsi

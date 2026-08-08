#pragma once

#include "softwing/soft_body.h"

namespace softwing {

struct SoftBodyCheckpoint::State {
    std::uint64_t topologyFingerprint = 0;
    std::vector<Node> nodes;
    std::vector<Triangle> triangles;
    std::vector<DistanceConstraint> constraints;
    std::vector<MembraneElement> membranes;
    std::vector<DihedralBendingConstraint> dihedrals;
    std::vector<std::pair<ContactFeatureKey, double>> contactMultipliers;
    std::vector<ContactRecord> contactRecords;
    ContactDiagnostics contactDiagnostics;
    std::vector<ContactDiagnostics> contactPairDiagnostics;
    std::vector<ContactFeatureKey> iterationCandidateKeys;
    std::vector<ContactFeatureKey> iterationQueryKeys;
    std::vector<ContactFeatureKey> certificationCandidateKeys;
    std::vector<ContactFeatureKey> certificationQueryKeys;
};

} // namespace softwing

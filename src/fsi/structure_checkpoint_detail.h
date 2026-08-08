#pragma once

#include "structure.h"

#include <softwing/soft_body.h>
#include <softwing/suspension.h>

namespace simwing::fsi {

struct StructureCheckpoint::Detail {
    softwing::SoftBodyCheckpoint body;
    std::optional<softwing::SuspensionCheckpoint> suspension;
    std::vector<StructureNodeState> publicNodes;
};

} // namespace simwing::fsi

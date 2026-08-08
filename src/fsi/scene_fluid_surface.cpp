#include "scene_fluid_surface.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

template<typename Value>
std::optional<std::size_t> indexOf(
    const std::vector<Value>& values,
    const Value value) noexcept {
    const auto found = std::lower_bound(values.begin(), values.end(), value);
    if (found == values.end() || *found != value) {
        return std::nullopt;
    }
    return static_cast<std::size_t>(found - values.begin());
}

template<typename Entity>
std::vector<const Entity*> sortedPointers(const std::vector<Entity>& values) {
    std::vector<const Entity*> result;
    result.reserve(values.size());
    for (const Entity& value : values) {
        result.push_back(&value);
    }
    std::sort(result.begin(), result.end(), [](const Entity* first,
                                               const Entity* second) {
        return first->id < second->id;
    });
    return result;
}

template<typename Entity>
const Entity& entityWithId(
    const std::vector<const Entity*>& values,
    const StableId id) {
    const auto found = std::lower_bound(
        values.begin(), values.end(), id,
        [](const Entity* value, const StableId expected) {
            return value->id < expected;
        });
    if (found == values.end() || (*found)->id != id) {
        throw std::logic_error("validated scene entity lookup failed");
    }
    return **found;
}

void sortUnique(std::vector<StableId>& values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

bool checkedAdd(std::size_t& total, const std::size_t value) {
    if (value > std::numeric_limits<std::size_t>::max() - total) {
        return false;
    }
    total += value;
    return true;
}

bool checkedBytes(std::size_t& total,
                  const std::size_t count,
                  const std::size_t itemBytes) {
    if (count != 0
        && itemBytes > std::numeric_limits<std::size_t>::max() / count) {
        return false;
    }
    return checkedAdd(total, count * itemBytes);
}

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            byte(static_cast<std::uint8_t>(value & 0xffU));
            value >>= 8U;
        }
    }

    template<typename Enum>
    void enumeration(const Enum value) {
        using Underlying = std::underlying_type_t<Enum>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    void byte(const std::uint8_t value) noexcept {
        value_ ^= value;
        value_ *= fnvPrime;
    }

    std::uint64_t value_ = fnvOffsetBasis;
};

std::uint64_t definitionFingerprint(
    const SceneFluidSurfaceDefinition& definition) {
    Fingerprint fingerprint;
    fingerprint.integer(definition.version);
    fingerprint.integer(
        static_cast<std::uint64_t>(definition.regions.size()));
    for (const auto& region : definition.regions) {
        fingerprint.integer(region.id);
        fingerprint.enumeration(region.kind);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(definition.materials.size()));
    for (const auto& material : definition.materials) {
        fingerprint.integer(material.id);
        fingerprint.real(material.porosityFraction);
        fingerprint.real(material.permeabilitySquareMeters);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(definition.vertices.size()));
    for (const auto& vertex : definition.vertices) {
        fingerprint.integer(vertex.id);
        fingerprint.real(vertex.referencePositionMeters.x);
        fingerprint.real(vertex.referencePositionMeters.y);
        fingerprint.real(vertex.referencePositionMeters.z);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(definition.triangles.size()));
    for (const auto& triangle : definition.triangles) {
        fingerprint.integer(triangle.id);
        for (const std::size_t index : triangle.vertexIndices) {
            fingerprint.integer(static_cast<std::uint64_t>(index));
        }
        fingerprint.integer(static_cast<std::uint64_t>(
            triangle.negativeSideRegionIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            triangle.positiveSideRegionIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            triangle.materialIndex));
        fingerprint.integer(triangle.sheetId);
        fingerprint.enumeration(triangle.role);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(definition.openings.size()));
    for (const auto& opening : definition.openings) {
        fingerprint.integer(opening.id);
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.orderedVertexIndices.size()));
        for (const std::size_t index : opening.orderedVertexIndices) {
            fingerprint.integer(static_cast<std::uint64_t>(index));
        }
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.negativeSideRegionIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.positiveSideRegionIndex));
        fingerprint.enumeration(opening.role);
    }
    return fingerprint.value();
}

std::uint64_t stateFingerprint(const SceneFluidSurfaceState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.definitionFingerprint);
    fingerprint.integer(state.structureDefinitionFingerprint);
    fingerprint.integer(state.acceptedStepCount);
    fingerprint.real(state.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(state.vertices.size()));
    for (const auto& vertex : state.vertices) {
        fingerprint.integer(vertex.id);
        fingerprint.real(vertex.positionMeters.x);
        fingerprint.real(vertex.positionMeters.y);
        fingerprint.real(vertex.positionMeters.z);
        fingerprint.real(vertex.velocityMetersPerSecond.x);
        fingerprint.real(vertex.velocityMetersPerSecond.y);
        fingerprint.real(vertex.velocityMetersPerSecond.z);
    }
    return fingerprint.value();
}

template<typename Entity>
bool idsMatch(const std::vector<Entity>& entities,
              const std::vector<StableId>& ids) {
    if (entities.size() != ids.size()) {
        return false;
    }
    for (std::size_t index = 0; index < ids.size(); ++index) {
        if (entities[index].id != ids[index]
            || ids[index] == invalidStableId
            || (index != 0 && ids[index - 1] >= ids[index])) {
            return false;
        }
    }
    return true;
}

bool validStableIds(const std::vector<StableId>& ids) {
    for (std::size_t index = 0; index < ids.size(); ++index) {
        if (ids[index] == invalidStableId
            || (index != 0 && ids[index - 1] >= ids[index])) {
            return false;
        }
    }
    return true;
}

void validateDefinitionImpl(const SceneFluidSurfaceDefinition& definition) {
    if (definition.version != sceneFluidSurfaceDefinitionVersion
        || definition.fingerprint == 0
        || !idsMatch(definition.regions, definition.mappings.regionIds)
        || !idsMatch(definition.materials, definition.mappings.materialIds)
        || !idsMatch(definition.vertices, definition.mappings.vertexIds)
        || !idsMatch(definition.triangles, definition.mappings.triangleIds)
        || !idsMatch(definition.openings, definition.mappings.openingIds)) {
        throw std::invalid_argument(
            "scene fluid-surface definition identity is invalid");
    }
    for (const auto& material : definition.materials) {
        if (!std::isfinite(material.porosityFraction)
            || !std::isfinite(material.permeabilitySquareMeters)) {
            throw std::invalid_argument(
                "scene fluid-surface material is non-finite");
        }
    }
    for (const auto& vertex : definition.vertices) {
        if (!std::isfinite(vertex.referencePositionMeters.x)
            || !std::isfinite(vertex.referencePositionMeters.y)
            || !std::isfinite(vertex.referencePositionMeters.z)) {
            throw std::invalid_argument(
                "scene fluid-surface vertex is non-finite");
        }
    }
    for (const auto& triangle : definition.triangles) {
        if (triangle.vertexIndices[0] >= definition.vertices.size()
            || triangle.vertexIndices[1] >= definition.vertices.size()
            || triangle.vertexIndices[2] >= definition.vertices.size()
            || triangle.negativeSideRegionIndex >= definition.regions.size()
            || triangle.positiveSideRegionIndex >= definition.regions.size()
            || triangle.materialIndex >= definition.materials.size()) {
            throw std::invalid_argument(
                "scene fluid-surface triangle mapping is invalid");
        }
    }
    for (const auto& opening : definition.openings) {
        if (opening.negativeSideRegionIndex >= definition.regions.size()
            || opening.positiveSideRegionIndex >= definition.regions.size()
            || std::any_of(
                opening.orderedVertexIndices.begin(),
                opening.orderedVertexIndices.end(),
                [&](const std::size_t index) {
                    return index >= definition.vertices.size();
                })) {
            throw std::invalid_argument(
                "scene fluid-surface opening mapping is invalid");
        }
    }
    if (definition.fingerprint != definitionFingerprint(definition)) {
        throw std::invalid_argument(
            "scene fluid-surface fingerprint does not match its definition");
    }
}

SceneFluidSurfaceAssembly failedAssembly(
    SceneFluidSurfaceDiagnostic diagnostic) {
    SceneFluidSurfaceAssembly result;
    result.diagnostics.push_back(std::move(diagnostic));
    return result;
}

} // namespace

void validateSceneFluidSurfaceDefinition(
    const SceneFluidSurfaceDefinition& definition) {
    validateDefinitionImpl(definition);
}

void validateSceneFluidSurfaceState(const SceneFluidSurfaceState& state) {
    if (state.version != sceneFluidSurfaceStateVersion
        || state.fingerprint == 0
        || state.definitionFingerprint == 0
        || state.structureDefinitionFingerprint == 0
        || !std::isfinite(state.simulationTimeSeconds)
        || state.simulationTimeSeconds < 0.0) {
        throw std::invalid_argument(
            "scene fluid-surface state identity is invalid");
    }
    StableId previousId = invalidStableId;
    for (const auto& vertex : state.vertices) {
        if (vertex.id == invalidStableId || vertex.id <= previousId
            || !std::isfinite(vertex.positionMeters.x)
            || !std::isfinite(vertex.positionMeters.y)
            || !std::isfinite(vertex.positionMeters.z)
            || !std::isfinite(vertex.velocityMetersPerSecond.x)
            || !std::isfinite(vertex.velocityMetersPerSecond.y)
            || !std::isfinite(vertex.velocityMetersPerSecond.z)) {
            throw std::invalid_argument(
                "scene fluid-surface state vertices are invalid");
        }
        previousId = vertex.id;
    }
    if (state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "scene fluid-surface state fingerprint does not match its payload");
    }
}

void validateSceneFluidSurfaceState(
    const SceneFluidSurfaceDefinition& definition,
    const SceneFluidSurfaceState& state) {
    validateSceneFluidSurfaceDefinition(definition);
    validateSceneFluidSurfaceState(state);
    if (state.definitionFingerprint != definition.fingerprint
        || state.vertices.size() != definition.vertices.size()) {
        throw std::invalid_argument(
            "scene fluid-surface state does not match its definition");
    }
    for (std::size_t index = 0; index < state.vertices.size(); ++index) {
        if (state.vertices[index].id != definition.vertices[index].id) {
            throw std::invalid_argument(
                "scene fluid-surface state vertex order does not match its definition");
        }
    }
}

std::optional<std::size_t> SceneFluidSurfaceMappings::regionIndex(
    const StableId id) const noexcept {
    return indexOf(regionIds, id);
}

std::optional<std::size_t> SceneFluidSurfaceMappings::materialIndex(
    const StableId id) const noexcept {
    return indexOf(materialIds, id);
}

std::optional<std::size_t> SceneFluidSurfaceMappings::vertexIndex(
    const StableId id) const noexcept {
    return indexOf(vertexIds, id);
}

std::optional<std::size_t> SceneFluidSurfaceMappings::triangleIndex(
    const StableId id) const noexcept {
    return indexOf(triangleIds, id);
}

std::optional<std::size_t> SceneFluidSurfaceMappings::openingIndex(
    const StableId id) const noexcept {
    return indexOf(openingIds, id);
}

bool SceneFluidSurfaceAssembly::ok() const noexcept {
    return assembled && diagnostics.empty();
}

namespace {

SceneFluidSurfaceAssembly assembleSceneFluidSurfaceImpl(
    const Scene& scene,
    const SceneFluidSurfaceLimits& limits) {
    const ValidationReport validation = validateScene(scene);
    if (!validation.ok()) {
        SceneFluidSurfaceAssembly result;
        result.diagnostics.reserve(validation.diagnostics.size());
        for (const auto& diagnostic : validation.diagnostics) {
            result.diagnostics.push_back({
                SceneFluidSurfaceDiagnosticCode::InvalidScene,
                diagnostic.entityKind,
                diagnostic.entityId,
                diagnostic.code,
                diagnostic.message,
            });
        }
        return result;
    }

    if (scene.triangles.size() > limits.maximumTriangles
        || scene.openings.size() > limits.maximumOpenings) {
        return failedAssembly({
            SceneFluidSurfaceDiagnosticCode::LimitExceeded,
            EntityKind::Scene,
            invalidStableId,
            std::nullopt,
            "scene fluid-surface entity count exceeds its configured limit",
        });
    }

    std::size_t openingVertexCount = 0;
    const auto orderedOpenings = sortedPointers(scene.openings);
    for (const Opening* opening : orderedOpenings) {
        if (!checkedAdd(
                openingVertexCount,
                opening->orderedVertexIds.size())
            || openingVertexCount > limits.maximumOpeningVertices) {
            return failedAssembly({
                SceneFluidSurfaceDiagnosticCode::LimitExceeded,
                EntityKind::Opening,
                opening->id,
                std::nullopt,
                "scene fluid-surface opening vertex count exceeds its configured limit",
            });
        }
    }

    std::vector<StableId> vertexIds;
    std::vector<StableId> materialIds;
    std::vector<StableId> regionIds;
    vertexIds.reserve(scene.triangles.size());
    materialIds.reserve(scene.triangles.size());
    regionIds.reserve(scene.triangles.size());
    for (const Triangle& triangle : scene.triangles) {
        vertexIds.insert(
            vertexIds.end(), triangle.vertexIds.begin(),
            triangle.vertexIds.end());
        materialIds.push_back(triangle.materialId);
        regionIds.push_back(triangle.negativeSideRegionId);
        regionIds.push_back(triangle.positiveSideRegionId);
    }
    sortUnique(vertexIds);
    sortUnique(materialIds);
    for (const Opening& opening : scene.openings) {
        regionIds.push_back(opening.negativeSideRegionId);
        regionIds.push_back(opening.positiveSideRegionId);
    }
    sortUnique(regionIds);

    std::vector<SceneFluidSurfaceDiagnostic> openingDiagnostics;
    for (const Opening* opening : orderedOpenings) {
        for (const StableId vertexId : opening->orderedVertexIds) {
            if (!std::binary_search(
                    vertexIds.begin(), vertexIds.end(), vertexId)) {
                openingDiagnostics.push_back({
                    SceneFluidSurfaceDiagnosticCode::
                        OpeningVertexWithoutFabricMotion,
                    EntityKind::Opening,
                    opening->id,
                    std::nullopt,
                    "opening references a vertex without a fabric-triangle Structure degree of freedom",
                });
                break;
            }
        }
    }
    if (!openingDiagnostics.empty()) {
        SceneFluidSurfaceAssembly result;
        result.diagnostics = std::move(openingDiagnostics);
        return result;
    }

    if (regionIds.size() > limits.maximumRegions
        || materialIds.size() > limits.maximumMaterials
        || vertexIds.size() > limits.maximumVertices
        || openingVertexCount > limits.maximumOpeningVertices) {
        return failedAssembly({
            SceneFluidSurfaceDiagnosticCode::LimitExceeded,
            EntityKind::Scene,
            invalidStableId,
            std::nullopt,
            "scene fluid-surface entity count exceeds its configured limit",
        });
    }

    std::size_t mappingBytes = 0;
    std::size_t stableIdCount = 0;
    if (!checkedAdd(stableIdCount, regionIds.size())
        || !checkedAdd(stableIdCount, materialIds.size())
        || !checkedAdd(stableIdCount, vertexIds.size())
        || !checkedAdd(stableIdCount, scene.triangles.size())
        || !checkedAdd(stableIdCount, scene.openings.size())
        || !checkedBytes(mappingBytes, stableIdCount, sizeof(StableId))
        || !checkedBytes(
            mappingBytes, scene.triangles.size(),
            6 * sizeof(std::size_t))
        || !checkedBytes(
            mappingBytes, openingVertexCount, sizeof(std::size_t))
        || !checkedBytes(
            mappingBytes, scene.openings.size(),
            2 * sizeof(std::size_t))
        || mappingBytes > limits.maximumMappingBytes) {
        return failedAssembly({
            SceneFluidSurfaceDiagnosticCode::MappingOverflow,
            EntityKind::Scene,
            invalidStableId,
            std::nullopt,
            "scene fluid-surface mapping exceeds its byte limit",
        });
    }

    const auto sourceRegions = sortedPointers(scene.regions);
    const auto sourceMaterials = sortedPointers(scene.fabricMaterials);
    const auto sourceVertices = sortedPointers(scene.vertices);
    const auto sourceTriangles = sortedPointers(scene.triangles);

    SceneFluidSurfaceDefinition definition;
    definition.mappings.regionIds = regionIds;
    definition.mappings.materialIds = materialIds;
    definition.mappings.vertexIds = vertexIds;
    definition.mappings.triangleIds.reserve(sourceTriangles.size());
    definition.mappings.openingIds.reserve(orderedOpenings.size());

    definition.regions.reserve(regionIds.size());
    for (const StableId id : regionIds) {
        const FluidRegion& source = entityWithId(sourceRegions, id);
        definition.regions.push_back({source.id, source.kind});
    }
    definition.materials.reserve(materialIds.size());
    for (const StableId id : materialIds) {
        const FabricMaterial& source = entityWithId(sourceMaterials, id);
        definition.materials.push_back({
            source.id,
            source.porosityFraction,
            source.permeabilitySquareMeters,
        });
    }
    definition.vertices.reserve(vertexIds.size());
    for (const StableId id : vertexIds) {
        const Vertex& source = entityWithId(sourceVertices, id);
        definition.vertices.push_back({source.id, source.positionMeters});
    }
    definition.triangles.reserve(sourceTriangles.size());
    for (const Triangle* source : sourceTriangles) {
        SceneFluidSurfaceTriangle triangle;
        triangle.id = source->id;
        for (std::size_t corner = 0; corner < 3; ++corner) {
            triangle.vertexIndices[corner] =
                *definition.mappings.vertexIndex(source->vertexIds[corner]);
        }
        triangle.negativeSideRegionIndex =
            *definition.mappings.regionIndex(source->negativeSideRegionId);
        triangle.positiveSideRegionIndex =
            *definition.mappings.regionIndex(source->positiveSideRegionId);
        triangle.materialIndex =
            *definition.mappings.materialIndex(source->materialId);
        triangle.sheetId = source->sheetId;
        triangle.role = source->role;
        definition.triangles.push_back(std::move(triangle));
        definition.mappings.triangleIds.push_back(source->id);
    }
    definition.openings.reserve(orderedOpenings.size());
    for (const Opening* source : orderedOpenings) {
        SceneFluidSurfaceOpening opening;
        opening.id = source->id;
        opening.orderedVertexIndices.reserve(source->orderedVertexIds.size());
        for (const StableId id : source->orderedVertexIds) {
            opening.orderedVertexIndices.push_back(
                *definition.mappings.vertexIndex(id));
        }
        opening.negativeSideRegionIndex =
            *definition.mappings.regionIndex(source->negativeSideRegionId);
        opening.positiveSideRegionIndex =
            *definition.mappings.regionIndex(source->positiveSideRegionId);
        opening.role = source->role;
        definition.openings.push_back(std::move(opening));
        definition.mappings.openingIds.push_back(source->id);
    }
    definition.fingerprint = definitionFingerprint(definition);
    validateSceneFluidSurfaceDefinition(definition);

    SceneFluidSurfaceAssembly result;
    result.assembled = true;
    result.definition = std::move(definition);
    return result;
}

} // namespace

SceneFluidSurfaceAssembly assembleSceneFluidSurface(
    const Scene& scene,
    const SceneFluidSurfaceLimits& limits) {
    try {
        return assembleSceneFluidSurfaceImpl(scene, limits);
    } catch (const std::bad_alloc&) {
        return failedAssembly({
            SceneFluidSurfaceDiagnosticCode::MappingOverflow,
            EntityKind::Scene,
            invalidStableId,
            std::nullopt,
            "scene fluid-surface assembly exhausted its bounded allocation",
        });
    }
}

SceneFluidSurfaceState captureSceneFluidSurfaceState(
    const SceneFluidSurfaceDefinition& definition,
    const SceneStructureMappings& structureMappings,
    const Structure& structure) {
    validateSceneFluidSurfaceDefinition(definition);
    const StructureDefinition& structureDefinition = structure.definition();
    std::size_t mappedNodeCount = structureMappings.nodeVertexIds.size();
    if (!validStableIds(structureMappings.nodeVertexIds)
        || !validStableIds(structureMappings.nodeSuspensionJunctionIds)
        || !validStableIds(structureMappings.triangleIds)
        || !checkedAdd(
            mappedNodeCount,
            structureMappings.nodeSuspensionJunctionIds.size())
        || mappedNodeCount != structureDefinition.nodes.size()
        || structureMappings.triangleIds.size()
            != structureDefinition.triangles.size()
        || std::any_of(
            structureMappings.nodeSuspensionJunctionIds.begin(),
            structureMappings.nodeSuspensionJunctionIds.end(),
            [&](const StableId id) {
                return std::binary_search(
                    structureMappings.nodeVertexIds.begin(),
                    structureMappings.nodeVertexIds.end(), id);
            })) {
        throw std::invalid_argument(
            "scene fluid-surface Structure mapping is invalid");
    }

    const StructureCheckpoint checkpoint = structure.checkpoint();
    if (checkpoint.definitionFingerprint
            != structure.definitionFingerprint()
        || checkpoint.nodes.size() != structureDefinition.nodes.size()) {
        throw std::invalid_argument(
            "scene fluid-surface Structure checkpoint identity is invalid");
    }
    const auto& nodeStates = checkpoint.nodes;
    std::vector<std::size_t> surfaceNodeIndices;
    surfaceNodeIndices.reserve(definition.vertices.size());
    std::vector<bool> usedNodes(nodeStates.size(), false);
    for (const auto& vertex : definition.vertices) {
        const std::optional<std::size_t> nodeIndex =
            structureMappings.nodeIndex(vertex.id);
        if (!nodeIndex || *nodeIndex >= nodeStates.size()
            || usedNodes[*nodeIndex]) {
            throw std::invalid_argument(
                "scene fluid-surface vertex is not uniquely mapped to Structure");
        }
        const StructureVector3& reference =
            structureDefinition.nodes[*nodeIndex].positionMeters;
        if (vertex.referencePositionMeters.x != reference.x
            || vertex.referencePositionMeters.y != reference.y
            || vertex.referencePositionMeters.z != reference.z) {
            throw std::invalid_argument(
                "scene fluid-surface reference geometry does not match Structure");
        }
        usedNodes[*nodeIndex] = true;
        surfaceNodeIndices.push_back(*nodeIndex);
    }
    for (const auto& triangle : definition.triangles) {
        const std::optional<std::size_t> triangleIndex =
            structureMappings.triangleIndex(triangle.id);
        if (!triangleIndex
            || *triangleIndex >= structureDefinition.triangles.size()) {
            throw std::invalid_argument(
                "scene fluid-surface triangle is not mapped to Structure");
        }
        const auto& structureTriangle =
            structureDefinition.triangles[*triangleIndex];
        for (std::size_t corner = 0; corner < 3; ++corner) {
            if (structureTriangle.nodes[corner]
                != surfaceNodeIndices[triangle.vertexIndices[corner]]) {
                throw std::invalid_argument(
                    "scene fluid-surface winding does not match Structure");
            }
        }
    }

    SceneFluidSurfaceState result;
    result.definitionFingerprint = definition.fingerprint;
    result.structureDefinitionFingerprint = checkpoint.definitionFingerprint;
    result.acceptedStepCount = checkpoint.acceptedStepCount;
    result.simulationTimeSeconds = checkpoint.simulationTimeSeconds;
    result.vertices.reserve(definition.vertices.size());
    for (std::size_t index = 0; index < definition.vertices.size(); ++index) {
        const auto& vertex = definition.vertices[index];
        const StructureNodeState& node =
            nodeStates[surfaceNodeIndices[index]];
        result.vertices.push_back({
            vertex.id,
            {node.positionMeters.x,
             node.positionMeters.y,
             node.positionMeters.z},
            {node.velocityMetersPerSecond.x,
             node.velocityMetersPerSecond.y,
             node.velocityMetersPerSecond.z},
        });
    }
    result.fingerprint = stateFingerprint(result);
    validateSceneFluidSurfaceState(definition, result);
    return result;
}

} // namespace simwing::fsi

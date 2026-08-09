#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <set>
#include <span>
#include <vector>

namespace simwing::fsi::detail {

// Validates only the stable combinatorial contract shared by scene-v2 and the
// compact fluid adapter. Geometry and accepted-motion checks belong to the
// opening-cap owner.
template<typename Index>
bool validOrientedBoundaryDisk(
    const std::span<const Index> boundary,
    const std::span<const std::array<Index, 3>> triangles) {
    if (triangles.empty()) {
        return true;
    }
    if (boundary.size() < 3 || triangles.size() != boundary.size() - 2) {
        return false;
    }

    using DirectedEdge = std::array<Index, 2>;
    const auto edgeKey = [](const Index first, const Index second) {
        return std::array{std::min(first, second),
                          std::max(first, second)};
    };
    const std::set<Index> boundaryVertices(
        boundary.begin(), boundary.end());
    if (boundaryVertices.size() != boundary.size()) {
        return false;
    }
    std::map<std::array<Index, 2>, DirectedEdge> boundaryEdges;
    for (std::size_t edge = 0; edge < boundary.size(); ++edge) {
        const DirectedEdge directed{
            boundary[edge], boundary[(edge + 1) % boundary.size()]};
        if (!boundaryEdges.emplace(
                edgeKey(directed[0], directed[1]), directed).second) {
            return false;
        }
    }

    std::set<std::array<Index, 3>> triangleKeys;
    std::map<std::array<Index, 2>, std::vector<DirectedEdge>> incidences;
    for (const auto& triangle : triangles) {
        if (triangle[0] == triangle[1] || triangle[1] == triangle[2]
            || triangle[2] == triangle[0]
            || !boundaryVertices.contains(triangle[0])
            || !boundaryVertices.contains(triangle[1])
            || !boundaryVertices.contains(triangle[2])) {
            return false;
        }
        auto triangleKey = triangle;
        std::sort(triangleKey.begin(), triangleKey.end());
        if (!triangleKeys.insert(triangleKey).second) {
            return false;
        }
        for (std::size_t edge = 0; edge < 3; ++edge) {
            const DirectedEdge directed{
                triangle[edge], triangle[(edge + 1) % 3]};
            incidences[edgeKey(directed[0], directed[1])].push_back(
                directed);
        }
    }
    for (const auto& [edge, directed] : boundaryEdges) {
        const auto found = incidences.find(edge);
        if (found == incidences.end() || found->second.size() != 1
            || found->second.front() != directed) {
            return false;
        }
    }
    for (const auto& [edge, directed] : incidences) {
        if (boundaryEdges.contains(edge)) {
            continue;
        }
        if (directed.size() != 2
            || directed[0][0] != directed[1][1]
            || directed[0][1] != directed[1][0]) {
            return false;
        }
    }
    return true;
}

} // namespace simwing::fsi::detail

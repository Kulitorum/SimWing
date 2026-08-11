#pragma once

#include "viewer_protocol.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

namespace simwing::viewer {

inline constexpr std::string_view cfdGridVertexBeginFieldName =
    "viewer.cfd_grid.vertex_begin";
inline constexpr std::string_view cfdGridCellCountXFieldName =
    "viewer.cfd_grid.cells_x";
inline constexpr std::string_view cfdGridCellCountYFieldName =
    "viewer.cfd_grid.cells_y";
inline constexpr std::string_view cfdGridCellCountZFieldName =
    "viewer.cfd_grid.cells_z";

enum class CfdSliceAxis : std::uint8_t {
    X = 0,
    Y = 1,
    Z = 2,
};

struct CfdGridDescriptor {
    std::size_t vertexBegin = 0;
    std::array<std::size_t, 3> cellCounts{};
    std::array<std::vector<double>, 3> cellCentreCoordinatesMetres;
    std::array<double, 3> cellSpacingMetres{};

    [[nodiscard]] std::size_t cellCount() const noexcept;
};

struct CfdSliceQuad {
    std::uint32_t sourceVertexIndex = 0;
    std::array<Vec3d, 4> cornersMetres;

    [[nodiscard]] bool operator==(
        const CfdSliceQuad& other) const noexcept;
};

struct CfdSliceGeometry {
    CfdSliceAxis axis = CfdSliceAxis::Z;
    std::size_t sliceIndex = 0;
    double coordinateMetres = 0.0;
    std::vector<CfdSliceQuad> quads;

};

// Reads the explicit X-fast structured-grid descriptor carried by a frame.
// Malformed or incomplete metadata returns no descriptor; trace input remains
// an untrusted boundary and cannot make the renderer index outside a frame.
[[nodiscard]] std::optional<CfdGridDescriptor> describeCfdGrid(
    const DiagnosticFrame& frame);

// Builds one cell-centred quad per cell on the selected axis-aligned plane.
// Each quad retains its source vertex index so scalar colour and vector glyphs
// use the exact immutable diagnostic sample that generated the cell.
[[nodiscard]] std::optional<CfdSliceGeometry> buildCfdSliceGeometry(
    const DiagnosticFrame& frame,
    CfdSliceAxis axis,
    std::size_t sliceIndex);

} // namespace simwing::viewer

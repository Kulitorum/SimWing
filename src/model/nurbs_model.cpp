#include "nurbs_model.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BinXCAFDrivers.hxx>
#include <BSplCLib.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GCPnts_QuasiUniformDeflection.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomConvert.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_HArray1.hxx>
#include <NCollection_Sequence.hxx>
#include <OSD_Parallel.hxx>
#include <PCDM_StoreStatus.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Precision.hxx>
#include <Quantity_Color.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <Standard_ConstructionError.hxx>
#include <Standard_Failure.hxx>
#include <TCollection_ExtendedString.hxx>
#include <TDF_Label.hxx>
#include <TDataStd_Name.hxx>
#include <TDocStd_Application.hxx>
#include <TDocStd_Document.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_ColorType.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <gp_XY.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <exception>
#include <limits>
#include <list>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

// LEparagliding's engineering coordinates are centimetres. OCCT and the STEP
// files produced here use millimetres, which avoids an ambiguous unitless CAD
// model and matches the convention used by most STEP consumers.
constexpr double millimetresPerCentimetre = 10.0;
constexpr double pointToleranceMillimetres = 0.005;
constexpr double maximumSourceDeviationMillimetres = 0.01;
constexpr double maximumLegacyAgreementMillimetres = 0.00001;
// The coarse scene mesh welds each panel edge onto the independently captured
// rib station. Keep that explicit canonicalization well below the rib mesher's
// 1 mm deflection; larger disagreement indicates incompatible captures.
constexpr double maximumSceneRibWeldDeviationMillimetres = 0.25;
// A curve whose points all stay this close to the symmetry plane is on the
// wing centre; its mirror copy would coincide with it.
constexpr double symmetryPlaneToleranceMillimetres = 0.5;
// Face boundaries are already exact B-spline edges. These sample counts add
// only a clear interior surface wireframe, without duplicating boundaries or
// bloating the STEP presentation.
constexpr int maximumDisplayedSpanSamples = 5;
constexpr int maximumDisplayedChordSamples = 10;

// Stable part-tree vocabulary. The viewport recognises these exact names in
// the STEP assembly, so change them only together with the GUI.
constexpr const char *wingGroupName = "Wing";
constexpr const char *extradosGroupName = "Extrados";
constexpr const char *ventsGroupName = "Vents";
constexpr const char *intradosGroupName = "Intrados";
// The interior surface-wireframe curves live in sibling groups next to the
// skin groups so CAD tools can hide them independently of the surfaces.
constexpr const char *extradosCurvesGroupName = "Extrados curves";
constexpr const char *ventCurvesGroupName = "Vent curves";
constexpr const char *intradosCurvesGroupName = "Intrados curves";
constexpr const char *ribsGroupName = "Ribs";
constexpr const char *linesGroupName = "Lines";
constexpr const char *brakeGroupName = "Brake lines";
constexpr const char *diagonalsGroupName = "Diagonals";
constexpr const char *miniribsGroupName = "Mini-ribs";
constexpr const char *otherCurvesName = "Other curves";
constexpr const char *rightSideName = "Right";
constexpr const char *leftSideName = "Left";
constexpr const char *centerSideName = "Center";

enum class Region
{
    Extrados,
    Vent,
    Intrados,
};

const char *regionGroupName(Region region)
{
    switch (region) {
    case Region::Extrados:
        return extradosGroupName;
    case Region::Vent:
        return ventsGroupName;
    case Region::Intrados:
        return intradosGroupName;
    }
    return extradosGroupName;
}

const char *regionCurvesGroupName(Region region)
{
    switch (region) {
    case Region::Extrados:
        return extradosCurvesGroupName;
    case Region::Vent:
        return ventCurvesGroupName;
    case Region::Intrados:
        return intradosCurvesGroupName;
    }
    return extradosCurvesGroupName;
}

// Region names as used in diagnostics; these match the historical engine
// messages ("upper"/"vent"/"lower").
const char *regionDiagnosticName(Region region)
{
    switch (region) {
    case Region::Extrados:
        return "upper";
    case Region::Vent:
        return "vent";
    case Region::Intrados:
        return "lower";
    }
    return "upper";
}

struct PartColor
{
    double red = 0.0;
    double green = 0.0;
    double blue = 0.0;

    Quantity_Color color() const
    {
        return {red, green, blue, Quantity_TOC_RGB};
    }
};

// Default part colours embedded in the STEP file so external CAD tools show
// a structured wing. The Studio viewport applies its own user-configurable
// palette on top of these.
constexpr PartColor extradosColor{0.20, 0.57, 0.88};
constexpr PartColor intradosColor{0.45, 0.69, 0.90};
constexpr PartColor ventColor{0.93, 0.60, 0.23};
constexpr PartColor wireframeColor{0.37, 0.82, 1.0};
constexpr PartColor ribColor{0.78, 0.80, 0.84};
constexpr PartColor brakeColor{0.95, 0.83, 0.28};
constexpr PartColor otherCurveColor{0.62, 0.66, 0.72};
constexpr PartColor diagonalColor{0.83, 0.45, 0.74};
constexpr PartColor planColors[6] = {
    {0.89, 0.29, 0.29}, // Plan A
    {0.95, 0.62, 0.19}, // Plan B
    {0.35, 0.79, 0.42}, // Plan C
    {0.29, 0.74, 0.86}, // Plan D
    {0.72, 0.47, 0.90}, // Plan E
    {0.62, 0.66, 0.72}, // Plan F
};

PartColor regionColor(Region region)
{
    switch (region) {
    case Region::Extrados:
        return extradosColor;
    case Region::Vent:
        return ventColor;
    case Region::Intrados:
        return intradosColor;
    }
    return extradosColor;
}

struct CapturedLine
{
    gp_Pnt start;
    gp_Pnt end;
    int colorIndex = 0;
    int planIndex = 0;
    bool brake = false;
    std::string label;
    std::string typeName;
    double capturedDiameterMillimetres = 0.0;
};

// A diagonal-rib or mini-rib sheet, captured as its two boundary polylines.
// Matching sample counts are in exact rung correspondence: sample j of
// curveA pairs with sample j of curveB, and the writer spans an exact ruled
// surface through every rung. Mini-rib boundaries have independent counts
// and are paired proportionally instead.
struct CapturedStrip
{
    bool minirib = false;
    std::string label;
    std::vector<gp_Pnt> curveA;
    std::vector<gp_Pnt> curveB;
};

struct SourcePanel
{
    const double *u = nullptr;
    const double *v = nullptr;
    const double *w = nullptr;
    const double *shapingHeight = nullptr;
    const double *legacyTessellation = nullptr;
    int panelIndex = 0;
    int segmentCount = 0;
};

// One skin region of one panel, sampled coarsely on the exact ballooning
// law for the Playground simulation mesh: a row per kept chordwise
// station, spanColumnCount samples across the cell. Rows keep their
// legacy station index so rib loops and neighbouring panels weld on
// identical stations.
constexpr int simSpanColumnCount = 5;
struct SimRegionCapture
{
    int panelIndex = 0;
    // Which skin the rows came from, so the Playground can draw the
    // surfaces separately (hiding the extrados to look inside).
    Region surface = Region::Extrados;
    // Vent samples exist for the legacy Playground even when the authored
    // model leaves the intake open. Scene-v2 must never mistake that toy-only
    // closure for fabric.
    bool authoredSurface = true;
    bool singleSkin = false;
    std::vector<int> stations;
    std::vector<std::array<gp_Pnt, simSpanColumnCount>> rows;

    // True when the panel spans symmetrically across x=0, so mirroring it
    // would reproduce the panel itself.
    [[nodiscard]] bool selfMirrored() const
    {
        if (rows.empty()) {
            return false;
        }
        const auto &row = rows.front();
        return std::abs(row.front().X() + row.back().X()) < 0.5;
    }
};

struct PanelSurface
{
    Region region = Region::Extrados;
    int panelIndex = 0;
    // Chordwise airfoil point range covered by this surface; rib faces use
    // it to know which parts of a rib outline are panel boundaries.
    int firstPoint = 0;
    int lastPoint = 0;
    occ::handle<Geom_BSplineSurface> surface;
    TopoDS_Face rightFace;
    TopoDS_Face leftFace;
    std::vector<TopoDS_Edge> rightWireframe;
    std::vector<TopoDS_Edge> leftWireframe;
};

// One airfoil hole from the legacy hole table, scaled to model millimetres.
// The rotation stays in the legacy convention: radians for ellipses (type 1),
// degrees for triangles and rectangles (types 3 and 4).
struct RibHole
{
    int type = 0;
    double x = 0.0;
    double y = 0.0;
    double a = 0.0;
    double b = 0.0;
    double rotation = 0.0;
    double cornerRadius = 0.0;
};

// A rib station: the chord-scaled planar profile and its rigidly placed 3D
// contour, in exact point correspondence, plus the rib's hole table. The
// planar-to-spatial correspondence recovers the rib plane frame that places
// the hole outlines, which the legacy core only ever draws in 2D.
struct CapturedRib
{
    std::vector<gp_XY> planarPoints;
    std::vector<gp_Pnt> spatialPoints;
    std::vector<RibHole> holes;
};

// The orthonormal frame mapping the rib's planar coordinates into model
// space, fitted from the captured point correspondence and only trusted
// when the fit reproduces every contour point.
struct RibFrame
{
    gp_Pnt origin;
    gp_Vec axisX;
    gp_Vec axisY;

    gp_Pnt point(const gp_XY &planar) const
    {
        return origin.Translated(
            axisX.Multiplied(planar.X()) + axisY.Multiplied(planar.Y()));
    }

    gp_Vec direction(const gp_XY &planar) const
    {
        return axisX.Multiplied(planar.X()) + axisY.Multiplied(planar.Y());
    }
};

// One chordwise airfoil point range of a rib outline covered by a panel
// surface. The rib face splits its outline at these region boundaries.
struct RibBoundarySegment
{
    int firstPoint = 0;
    int lastPoint = 0;
};

struct QuantizedPoint
{
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::int64_t z = 0;

    bool operator==(const QuantizedPoint &) const = default;

    bool operator<(const QuantizedPoint &other) const
    {
        return std::tie(x, y, z) < std::tie(other.x, other.y, other.z);
    }
};

struct QuantizedPointHash
{
    std::size_t operator()(const QuantizedPoint &point) const noexcept
    {
        std::size_t result = std::hash<std::int64_t>{}(point.x);
        result ^= std::hash<std::int64_t>{}(point.y) + 0x9e3779b9U
                  + (result << 6U) + (result >> 2U);
        result ^= std::hash<std::int64_t>{}(point.z) + 0x9e3779b9U
                  + (result << 6U) + (result >> 2U);
        return result;
    }
};

struct QuantizedSegment
{
    QuantizedPoint first;
    QuantizedPoint second;

    bool operator==(const QuantizedSegment &) const = default;
};

struct QuantizedSegmentHash
{
    std::size_t operator()(const QuantizedSegment &segment) const noexcept
    {
        QuantizedPointHash hash;
        std::size_t result = hash(segment.first);
        result ^= hash(segment.second) + 0x9e3779b9U
                  + (result << 6U) + (result >> 2U);
        return result;
    }
};

// A coincident segment may intentionally carry a different authored row or
// brake role.  Geometry alone is therefore not enough to identify a duplicate
// in the Playground handoff: keep those semantic alternatives, while folding
// repeated captures (including reversed endpoints) of the same physical line.
struct QuantizedSemanticSegment
{
    QuantizedSegment geometry;
    int planIndex = 0;
    bool brake = false;

    bool operator==(const QuantizedSemanticSegment &) const = default;
};

struct QuantizedSemanticSegmentHash
{
    std::size_t operator()(
        const QuantizedSemanticSegment &segment) const noexcept
    {
        QuantizedSegmentHash hash;
        std::size_t result = hash(segment.geometry);
        result ^= std::hash<int>{}(segment.planIndex) + 0x9e3779b9U
                  + (result << 6U) + (result >> 2U);
        result ^= std::hash<bool>{}(segment.brake) + 0x9e3779b9U
                  + (result << 6U) + (result >> 2U);
        return result;
    }
};

struct ParametricOpeningVertex
{
    simwing::fsi::StableId id = simwing::fsi::invalidStableId;
    double u = 0.0;
    double v = 0.0;
};

std::vector<std::array<simwing::fsi::StableId, 3>>
triangulateParametricOpeningBoundary(
    const std::vector<ParametricOpeningVertex> &boundary)
{
    using simwing::fsi::StableId;
    if (boundary.size() < 3) {
        return {};
    }
    std::set<StableId> uniqueIds;
    for (const ParametricOpeningVertex &vertex : boundary) {
        if (vertex.id == simwing::fsi::invalidStableId
            || !uniqueIds.insert(vertex.id).second) {
            return {};
        }
    }

    const auto turn = [&boundary](std::size_t first,
                                  std::size_t second,
                                  std::size_t third) {
        const auto &a = boundary[first];
        const auto &b = boundary[second];
        const auto &c = boundary[third];
        return (b.u - a.u) * (c.v - b.v)
            - (b.v - a.v) * (c.u - b.u);
    };
    const auto insideTriangle = [&turn](std::size_t point,
                                        std::size_t first,
                                        std::size_t second,
                                        std::size_t third) {
        return turn(first, second, point) >= 0
            && turn(second, third, point) >= 0
            && turn(third, first, point) >= 0;
    };

    std::vector<std::size_t> active(boundary.size());
    for (std::size_t index = 0; index < active.size(); ++index) {
        active[index] = index;
    }
    std::vector<std::array<StableId, 3>> triangles;
    triangles.reserve(boundary.size() - 2);
    while (active.size() > 3) {
        bool clipped = false;
        for (std::size_t activeIndex = 0;
             activeIndex < active.size(); ++activeIndex) {
            const std::size_t previous = active[
                (activeIndex + active.size() - 1) % active.size()];
            const std::size_t current = active[activeIndex];
            const std::size_t next = active[
                (activeIndex + 1) % active.size()];
            if (turn(previous, current, next) <= 0) {
                continue;
            }
            bool containsVertex = false;
            for (const std::size_t candidate : active) {
                if (candidate == previous || candidate == current
                    || candidate == next) {
                    continue;
                }
                if (insideTriangle(
                        candidate, previous, current, next)) {
                    containsVertex = true;
                    break;
                }
            }
            if (containsVertex) {
                continue;
            }
            triangles.push_back({boundary[previous].id,
                                 boundary[current].id,
                                 boundary[next].id});
            active.erase(active.begin()
                         + static_cast<std::ptrdiff_t>(activeIndex));
            clipped = true;
            break;
        }
        if (!clipped) {
            return {};
        }
    }
    if (turn(active[0], active[1], active[2]) <= 0) {
        return {};
    }
    triangles.push_back({boundary[active[0]].id,
                         boundary[active[1]].id,
                         boundary[active[2]].id});
    return triangles;
}

QuantizedPoint quantize(const gp_Pnt &point)
{
    constexpr double scale = 1.0 / pointToleranceMillimetres;
    return {
        std::llround(point.X() * scale),
        std::llround(point.Y() * scale),
        std::llround(point.Z() * scale),
    };
}

bool pointLess(const QuantizedPoint &left, const QuantizedPoint &right)
{
    if (left.x != right.x) {
        return left.x < right.x;
    }
    if (left.y != right.y) {
        return left.y < right.y;
    }
    return left.z < right.z;
}

QuantizedSegment quantizeSegment(const gp_Pnt &start, const gp_Pnt &end)
{
    QuantizedPoint first = quantize(start);
    QuantizedPoint second = quantize(end);
    if (pointLess(second, first)) {
        std::swap(first, second);
    }
    return {first, second};
}

gp_Trsf mirrorTransform()
{
    gp_Trsf transform;
    transform.SetMirror(
        gp_Ax2(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(1.0, 0.0, 0.0)));
    return transform;
}

gp_Pnt modelPoint(double x, double y, double z, bool mirrored = false)
{
    // The calculation core uses Z down; CAD and the OCCT viewport use Z up.
    return {
        (mirrored ? -x : x) * millimetresPerCentimetre,
        y * millimetresPerCentimetre,
        -z * millimetresPerCentimetre,
    };
}

bool isFinite(const gp_Pnt &point)
{
    return std::isfinite(point.X())
           && std::isfinite(point.Y())
           && std::isfinite(point.Z());
}

std::string trimmedLabel(const char *label, int labelLength)
{
    std::string result;
    for (int index = 0; index < labelLength; ++index) {
        const char character = label[index];
        if (character != ' ' && character != '\0') {
            result.push_back(character);
        }
    }
    return result;
}

TopoDS_Shape mirrored(const TopoDS_Shape &shape)
{
    BRepBuilderAPI_Transform mirror(shape, mirrorTransform(), true);
    if (!mirror.IsDone()) {
        return {};
    }
    return mirror.Shape();
}

// Recovers the rigid planar-to-spatial frame of a rib by least squares over
// the captured contour correspondence. The legacy chain from the chord-scaled
// profile to the placed airfoil is a composition of rotations and
// translations, so the fit must come back orthonormal and reproduce every
// point; anything else means the correspondence was misread.
bool fitRibFrame(const std::vector<gp_XY> &planarPoints,
                 const std::vector<gp_Pnt> &spatialPoints,
                 RibFrame &frame)
{
    if (planarPoints.size() != spatialPoints.size()
        || planarPoints.size() < 3) {
        return false;
    }
    const double count = static_cast<double>(planarPoints.size());

    gp_XY planarCentroid(0.0, 0.0);
    gp_Vec spatialCentroid(0.0, 0.0, 0.0);
    for (std::size_t index = 0; index < planarPoints.size(); ++index) {
        planarCentroid += planarPoints[index];
        spatialCentroid += gp_Vec(spatialPoints[index].XYZ());
    }
    planarCentroid /= count;
    spatialCentroid /= count;

    double sumXX = 0.0;
    double sumXY = 0.0;
    double sumYY = 0.0;
    gp_Vec sumXP(0.0, 0.0, 0.0);
    gp_Vec sumYP(0.0, 0.0, 0.0);
    for (std::size_t index = 0; index < planarPoints.size(); ++index) {
        const gp_XY planar = planarPoints[index] - planarCentroid;
        const gp_Vec spatial =
            gp_Vec(spatialPoints[index].XYZ()) - spatialCentroid;
        sumXX += planar.X() * planar.X();
        sumXY += planar.X() * planar.Y();
        sumYY += planar.Y() * planar.Y();
        sumXP += spatial.Multiplied(planar.X());
        sumYP += spatial.Multiplied(planar.Y());
    }
    const double determinant = sumXX * sumYY - sumXY * sumXY;
    if (std::abs(determinant) <= Precision::Confusion()) {
        return false;
    }
    frame.axisX = (sumXP.Multiplied(sumYY) - sumYP.Multiplied(sumXY))
                      .Divided(determinant);
    frame.axisY = (sumYP.Multiplied(sumXX) - sumXP.Multiplied(sumXY))
                      .Divided(determinant);
    frame.origin = gp_Pnt(
        (spatialCentroid
         - frame.axisX.Multiplied(planarCentroid.X())
         - frame.axisY.Multiplied(planarCentroid.Y()))
            .XYZ());

    constexpr double rigidityTolerance = 1.0e-6;
    if (std::abs(frame.axisX.Magnitude() - 1.0) > rigidityTolerance
        || std::abs(frame.axisY.Magnitude() - 1.0) > rigidityTolerance
        || std::abs(frame.axisX.Dot(frame.axisY)) > rigidityTolerance) {
        return false;
    }
    for (std::size_t index = 0; index < planarPoints.size(); ++index) {
        if (frame.point(planarPoints[index])
                .Distance(spatialPoints[index])
            > pointToleranceMillimetres) {
            return false;
        }
    }
    return true;
}

// Winding of a closed planar wire about the given plane normal, measured as
// the signed enclosed area. Used only for a sign decision, so a modest
// sampling of each edge is exact enough.
double signedAreaAlongNormal(const TopoDS_Wire &wire,
                             const gp_Pnt &origin,
                             const gp_Dir &normal)
{
    gp_Vec area(0.0, 0.0, 0.0);
    for (BRepTools_WireExplorer explorer(wire);
         explorer.More();
         explorer.Next()) {
        const TopoDS_Edge &edge = explorer.Current();
        double first = 0.0;
        double last = 0.0;
        const occ::handle<Geom_Curve> curve =
            BRep_Tool::Curve(edge, first, last);
        if (curve.IsNull()) {
            continue;
        }
        if (edge.Orientation() == TopAbs_REVERSED) {
            std::swap(first, last);
        }
        constexpr int sampleCount = 32;
        gp_Pnt previous = curve->Value(first);
        for (int sample = 1; sample <= sampleCount; ++sample) {
            const gp_Pnt current = curve->Value(
                first
                + (last - first) * static_cast<double>(sample)
                      / static_cast<double>(sampleCount));
            area += gp_Vec(origin, previous)
                        .Crossed(gp_Vec(origin, current))
                        .Multiplied(0.5);
            previous = current;
        }
    }
    return area.Dot(gp_Vec(normal));
}

gp_XY projectOntoLine(const gp_XY &lineStart,
                      const gp_XY &lineEnd,
                      const gp_XY &point)
{
    gp_XY direction = lineEnd - lineStart;
    const double length = direction.Modulus();
    if (length <= Precision::Confusion()) {
        return lineStart;
    }
    direction /= length;
    return lineStart + direction * ((point - lineStart) * direction);
}

// One leaf part of the exported assembly: a named shape with default
// surface/curve display colours.
struct AssemblyPart
{
    std::string name;
    TopoDS_Shape shape;
    PartColor faceColor;
    PartColor curveColor;
    bool hasFaces = false;
};

// A named group of leaf parts (one level of the assembly tree). Children
// live in a std::list so references and pointers to a group stay valid
// while sibling groups are still being created.
struct AssemblyGroup
{
    std::string name;
    std::list<AssemblyGroup> groups;
    std::vector<AssemblyPart> parts;

    AssemblyGroup &group(const std::string &groupName)
    {
        for (AssemblyGroup &child : groups) {
            if (child.name == groupName) {
                return child;
            }
        }
        groups.push_back({groupName, {}, {}});
        return groups.back();
    }

    bool empty() const
    {
        return parts.empty()
               && std::all_of(
                   groups.begin(),
                   groups.end(),
                   [](const AssemblyGroup &child) { return child.empty(); });
    }
};

class NurbsModel
{
public:
    void reset()
    {
        panels_.clear();
        capturedRibs_.clear();
        capturedLines_.clear();
        capturedStrips_.clear();
        simRegions_.clear();
        errors_.clear();
        warnings_.clear();
        captureLines_ = false;
        currentLineTag_ = {};
        maximumSourceDeviation_ = 0.0;
        maximumLegacyAgreement_ = 0.0;
    }

    void capturePanel(const double *u,
                      const double *v,
                      const double *w,
                      const double *shapingHeight,
                      const double *legacyTessellation,
                      int panelIndex,
                      int totalPointCount,
                      int upperPointCount,
                      int ventPointCount,
                      int segmentCount,
                      bool includeVentSurface,
                      bool singleSkin)
    {
        if (u == nullptr
            || v == nullptr
            || w == nullptr
            || shapingHeight == nullptr
            || legacyTessellation == nullptr
            || panelIndex < 1
            || panelIndex > 100
            || totalPointCount < 2
            || totalPointCount > 500
            || upperPointCount < 2
            || upperPointCount > totalPointCount
            || ventPointCount < 2
            || segmentCount < 1
            || segmentCount > 98) {
            errors_.push_back(
                "Rejected invalid source-shape dimensions for panel "
                + std::to_string(panelIndex));
            return;
        }

        const int ventLast = upperPointCount + ventPointCount - 1;
        if (ventLast > totalPointCount) {
            errors_.push_back(
                "Vent point range exceeds the airfoil range for panel "
                + std::to_string(panelIndex));
            return;
        }

        const SourcePanel source{
            u,
            v,
            w,
            shapingHeight,
            legacyTessellation,
            panelIndex,
            segmentCount,
        };
        addRegion(source,
                  panelIndex,
                  1,
                  upperPointCount,
                  Region::Extrados);
        if (includeVentSurface) {
            addRegion(source,
                      panelIndex,
                      upperPointCount,
                      ventLast,
                      Region::Vent);
        }
        if (!singleSkin && ventLast < totalPointCount) {
            addRegion(source,
                      panelIndex,
                      ventLast,
                      totalPointCount,
                      Region::Intrados);
        }

        // Coarse skin samples for the Playground simulation mesh. Unlike
        // the STEP model, the vent wall is always captured: the toy
        // stamps a uniform pressure difference per face, so open intakes
        // would both leave flapping rims and give the closed-surface
        // pressure field a spurious net thrust through the hole.
        captureSimRegion(source,
                         panelIndex,
                         1,
                         upperPointCount,
                         16,
                         Region::Extrados,
                         true,
                         singleSkin);
        captureSimRegion(source,
                         panelIndex,
                         upperPointCount,
                         ventLast,
                         3,
                         Region::Vent,
                         includeVentSurface,
                         singleSkin);
        if (!singleSkin && ventLast < totalPointCount) {
            captureSimRegion(source,
                             panelIndex,
                             ventLast,
                             totalPointCount,
                             12,
                             Region::Intrados,
                             true,
                             singleSkin);
        }
    }

    void captureRib(const double *u,
                    const double *v,
                    const double *w,
                    const double *holes,
                    double chordCentimetres,
                    int ribIndex,
                    int totalPointCount)
    {
        if (u == nullptr
            || v == nullptr
            || w == nullptr
            || holes == nullptr
            || ribIndex < 0
            || ribIndex > 100
            || totalPointCount < 3
            || totalPointCount > 500) {
            errors_.push_back(
                "Rejected invalid source-shape dimensions for rib "
                + std::to_string(ribIndex));
            return;
        }

        CapturedRib rib;
        rib.planarPoints.reserve(
            static_cast<std::size_t>(totalPointCount));
        rib.spatialPoints.reserve(
            static_cast<std::size_t>(totalPointCount));
        for (int pointIndex = 1;
             pointIndex <= totalPointCount;
             ++pointIndex) {
            const std::size_t planarIndex =
                sourceFieldIndex(ribIndex, pointIndex, 3);
            rib.planarPoints.emplace_back(
                u[planarIndex] * millimetresPerCentimetre,
                v[planarIndex] * millimetresPerCentimetre);
            const std::size_t spatialIndex =
                sourceFieldIndex(ribIndex, pointIndex, 47);
            rib.spatialPoints.push_back(
                modelPoint(u[spatialIndex],
                           v[spatialIndex],
                           w[spatialIndex]));
            if (!isFinite(rib.spatialPoints.back())) {
                errors_.push_back(
                    "Non-finite contour point in rib "
                    + std::to_string(ribIndex));
                return;
            }
        }

        const double chordMillimetres =
            chordCentimetres * millimetresPerCentimetre;
        const int holeCount = std::clamp(
            static_cast<int>(holes[ribIndex]), 0, 20);
        for (int holeIndex = 1; holeIndex <= holeCount; ++holeIndex) {
            const auto field = [&](int fieldIndex) {
                return holes[holeFieldIndex(ribIndex, holeIndex, fieldIndex)];
            };
            RibHole hole;
            hole.type = static_cast<int>(field(9));
            // Types 1, 3 and 4 are the airfoil holes the legacy core draws
            // into rib patterns; type 11 parameterizes unloaded miniribs,
            // which are not rib stations, and other values are ignored by
            // the legacy drawing code as well.
            if (hole.type != 1 && hole.type != 3 && hole.type != 4) {
                continue;
            }
            hole.x = field(2) * chordMillimetres / 100.0;
            hole.y = field(3) * chordMillimetres / 100.0;
            hole.a = field(4) * chordMillimetres / 100.0;
            hole.b = field(5) * chordMillimetres / 100.0;
            hole.rotation = field(6);
            hole.cornerRadius = field(7) * chordMillimetres / 100.0;
            rib.holes.push_back(hole);
        }

        capturedRibs_[ribIndex] = std::move(rib);
    }

    void captureLine(CapturedLine line)
    {
        if (!captureLines_
            || !isFinite(line.start)
            || !isFinite(line.end)
            || line.start.Distance(line.end) <= Precision::Confusion()) {
            return;
        }
        line.planIndex = currentLineTag_.planIndex;
        line.brake = currentLineTag_.brake;
        line.label = currentLineTag_.label;
        line.typeName = currentLineTag_.typeName;
        line.capturedDiameterMillimetres =
            currentLineTag_.diameterMillimetres;
        capturedLines_.push_back(std::move(line));
    }

    void setLineCapture(bool enabled)
    {
        captureLines_ = enabled;
        if (!enabled) {
            currentLineTag_ = {};
        }
    }

    void setLineTag(const char *label,
                    int labelLength,
                    int planIndex,
                    bool brake,
                    const char *typeName,
                    int typeNameLength,
                    double diameterMm)
    {
        std::string title =
            label != nullptr && labelLength > 0
                ? trimmedLabel(label, labelLength)
                : std::string();
        const std::string type =
            typeName != nullptr && typeNameLength > 0
                ? trimmedLabel(typeName, typeNameLength)
                : std::string();
        if (!type.empty()
            && std::isfinite(diameterMm)
            && diameterMm > 0.0) {
            std::ostringstream suffix;
            suffix << type << " " << diameterMm << " mm";
            title = title.empty()
                        ? suffix.str()
                        : title + " (" + suffix.str() + ")";
        }
        currentLineTag_.label = std::move(title);
        currentLineTag_.planIndex = std::clamp(planIndex, 0, 6);
        currentLineTag_.brake = brake;
        currentLineTag_.typeName = type;
        currentLineTag_.diameterMillimetres =
            std::isfinite(diameterMm) && diameterMm > 0.0
                ? diameterMm
                : 0.0;
    }

    void captureDiagonalStrip(const char *kind,
                              int kindLength,
                              int index,
                              const double *xA,
                              const double *yA,
                              const double *zA,
                              const double *xB,
                              const double *yB,
                              const double *zB,
                              int pointCount,
                              int stride)
    {
        if (xA == nullptr || yA == nullptr || zA == nullptr
            || xB == nullptr || yB == nullptr || zB == nullptr
            || pointCount < 2 || stride < 1) {
            return;
        }
        CapturedStrip strip;
        strip.label = (kind != nullptr && kindLength > 0
                           ? trimmedLabel(kind, kindLength)
                           : std::string("Diagonal"))
                      + " " + std::to_string(index);
        strip.curveA.reserve(pointCount);
        strip.curveB.reserve(pointCount);
        for (int sample = 0; sample < pointCount; ++sample) {
            const int offset = sample * stride;
            const gp_Pnt a =
                modelPoint(xA[offset], yA[offset], zA[offset]);
            const gp_Pnt b =
                modelPoint(xB[offset], yB[offset], zB[offset]);
            if (!isFinite(a) || !isFinite(b)) {
                return;
            }
            strip.curveA.push_back(a);
            strip.curveB.push_back(b);
        }
        capturedStrips_.push_back(std::move(strip));
    }

    void captureMiniribs(const double *x,
                         const double *y,
                         const double *z,
                         const int *np,
                         const double *rib,
                         int ribCount,
                         bool singleSkin)
    {
        if (x == nullptr || y == nullptr || z == nullptr || np == nullptr
            || rib == nullptr || singleSkin) {
            return;
        }
        // The legacy arrays keep their f2c layouts: x/y/z and rib are
        // Fortran (0:100,500) matrices, np is (0:100,9), all column-major
        // with 101 rows. Station 0 already holds the mirrored first rib.
        const auto midPoint = [&](int station, int j) {
            const int offset = station + (j - 1) * 101;
            return modelPoint(0.5 * (x[offset - 1] + x[offset]),
                              0.5 * (y[offset - 1] + y[offset]),
                              0.5 * (z[offset - 1] + z[offset]));
        };
        const auto ribValue = [&](int station, int column) {
            return rib[station + (column - 1) * 101];
        };
        for (int i = 1; i <= std::min(ribCount, 100); ++i) {
            const double percent = ribValue(i, 56);
            if (percent <= 1.0) {
                continue;
            }
            const int totalPoints = np[i];
            const int upperPoints = np[i + 101];
            if (totalPoints < 4 || totalPoints > 500 || upperPoints < 2
                || upperPoints >= totalPoints) {
                continue;
            }
            // Mid-cell section between ribs i-1 and i, where the legacy
            // program sews the mini-rib but never draws it in 3D.
            std::vector<gp_Pnt> section(totalPoints + 1);
            bool finite = true;
            for (int j = 1; j <= totalPoints; ++j) {
                section[j] = midPoint(i, j);
                finite = finite && isFinite(section[j]);
            }
            const gp_Pnt trailing = section[1];
            const gp_Pnt leading = section[upperPoints];
            if (!finite
                || trailing.Distance(leading) <= Precision::Confusion()) {
                continue;
            }
            const gp_Dir chordDirection(gp_Vec(trailing, leading));
            const auto chordwise = [&](const gp_Pnt &point) {
                return gp_Vec(trailing, point).Dot(chordDirection);
            };
            // The cut runs rib(i-1,5)*pct/100 centimetres ahead of the
            // trailing edge, measured on the left rib chord exactly as the
            // legacy pattern code does; 100% keeps the whole section (a
            // complete unloaded middle rib).
            double chord = ribValue(i - 1, 5);
            if (chord <= 0.0) {
                chord = ribValue(i, 5);
            }
            const double cut = chord * millimetresPerCentimetre
                               * std::min(percent, 100.0) / 100.0;
            const auto appendCut = [&](std::vector<gp_Pnt> &curve,
                                       int fromIndex,
                                       int toIndex) {
                const double before = chordwise(section[fromIndex]);
                const double after = chordwise(section[toIndex]);
                if (after <= before + Precision::Confusion()) {
                    return;
                }
                const double t =
                    std::clamp((cut - before) / (after - before), 0.0, 1.0);
                if (t <= Precision::Confusion()) {
                    return;
                }
                const gp_Pnt &from = section[fromIndex];
                const gp_Pnt &to = section[toIndex];
                curve.emplace_back(from.X() + t * (to.X() - from.X()),
                                   from.Y() + t * (to.Y() - from.Y()),
                                   from.Z() + t * (to.Z() - from.Z()));
            };

            CapturedStrip strip;
            strip.minirib = true;
            strip.label = "Mini-rib " + std::to_string(i);
            // Upper boundary: trailing edge forward along the extrados.
            int j = 1;
            strip.curveA.push_back(section[1]);
            while (j + 1 <= upperPoints
                   && chordwise(section[j + 1]) < cut) {
                ++j;
                strip.curveA.push_back(section[j]);
            }
            if (j + 1 <= upperPoints) {
                appendCut(strip.curveA, j, j + 1);
            }
            // Lower boundary: trailing edge forward along the intrados,
            // continuing across the vent contour when a 100% rib reaches
            // the leading edge.
            j = totalPoints;
            strip.curveB.push_back(section[totalPoints]);
            while (j - 1 > upperPoints
                   && chordwise(section[j - 1]) < cut) {
                --j;
                strip.curveB.push_back(section[j]);
            }
            if (j - 1 > upperPoints) {
                appendCut(strip.curveB, j, j - 1);
            }
            if (strip.curveA.size() >= 2 && strip.curveB.size() >= 2) {
                capturedStrips_.push_back(std::move(strip));
            }
        }
    }

    lep::NurbsWriteResult writeStep(const std::filesystem::path &path,
                                    bool includeConstructionCurves)
    {
        warnings_.clear();
        lep::NurbsWriteResult result;
        result.maximumSourceDeviationMillimetres =
            maximumSourceDeviation_;
        result.maximumLegacyAgreementMillimetres =
            maximumLegacyAgreement_;

        if (!errors_.empty()) {
            result.error = errors_.front();
            return result;
        }
        if (panels_.empty()) {
            result.error = "The calculation produced no NURBS surfaces.";
            return result;
        }

        try {
            validateSewing(result);
            if (!result.error.empty()) {
                return result;
            }

            AssemblyGroup wing{wingGroupName, {}, {}};
            addPanelParts(wing, result, includeConstructionCurves);
            addRibParts(wing, result);
            if (!errors_.empty()) {
                result.error = errors_.front();
                return result;
            }
            addLineParts(wing, result);
            addStripParts(wing, result);

            writeAssembly(wing, path, result);
        } catch (const Standard_Failure &failure) {
            result.error =
                std::string("OCCT STEP export failed: ")
                + (failure.GetMessageString() != nullptr
                       ? failure.GetMessageString()
                       : "unknown OCCT error");
            return result;
        } catch (const std::exception &exception) {
            result.error =
                std::string("STEP export failed: ") + exception.what();
            return result;
        }

        result.warnings = warnings_;
        result.success = result.error.empty();
        return result;
    }

private:
    struct LineTag
    {
        std::string label;
        int planIndex = 0;
        bool brake = false;
        std::string typeName;
        double diameterMillimetres = 0.0;
    };

    static std::size_t sourceFieldIndex(int panelIndex,
                                        int pointIndex,
                                        int fieldIndex)
    {
        // f2c layout for real*8 u/v/w(0:100,500,99), using the
        // unadjusted arrays owned by MAIN__.
        return static_cast<std::size_t>(
            panelIndex + (pointIndex + fieldIndex * 500) * 101
            - 50601);
    }

    static std::size_t holeFieldIndex(int ribIndex,
                                      int holeIndex,
                                      int fieldIndex)
    {
        // f2c layout of the legacy hol array: hole property fieldIndex
        // (2 x%, 3 y%, 4 a%, 5 b%, 6 rotation, 7 corner radius %, 9 type)
        // of hole holeIndex, encoded as hol(rib, holeIndex + fieldIndex*200).
        // The hole count itself sits at the plain rib index.
        return static_cast<std::size_t>(
            ribIndex + (holeIndex + fieldIndex * 200) * 101 - 20301);
    }

    static std::size_t shapingHeightIndex(int panelIndex,
                                          int pointIndex)
    {
        // f2c layout for real*8 hautok(0:100,500).
        return static_cast<std::size_t>(
            panelIndex + pointIndex * 101 - 101);
    }

    static std::size_t tessellationIndex(int panelIndex, int pointIndex, int segmentIndex)
    {
        // f2c layout for real*8 tesse3d(3,0:100,500,99). This is the
        // unadjusted array passed by MAIN__; X/Y/Z occupy consecutive slots.
        return static_cast<std::size_t>(
            (panelIndex + (pointIndex + segmentIndex * 500) * 101) * 3
            - 151803);
    }

    gp_Pnt tessellationPoint(const double *tessellation,
                             int panelIndex,
                             int pointIndex,
                             int segmentIndex) const
    {
        const std::size_t index =
            tessellationIndex(panelIndex, pointIndex, segmentIndex);
        return modelPoint(tessellation[index],
                          tessellation[index + 1],
                          tessellation[index + 2]);
    }

    static gp_Pnt sourceControlPoint(const SourcePanel &source,
                                     int panelIndex,
                                     int pointIndex,
                                     int fieldIndex)
    {
        const std::size_t index =
            sourceFieldIndex(panelIndex, pointIndex, fieldIndex);
        return modelPoint(source.u[index],
                          source.v[index],
                          source.w[index]);
    }

    static double sourceShapingHeight(const SourcePanel &source,
                                      int pointIndex)
    {
        return source.shapingHeight[
                   shapingHeightIndex(source.panelIndex, pointIndex)]
               * millimetresPerCentimetre;
    }

    gp_Pnt sourceShapePoint(const SourcePanel &source,
                            int pointIndex,
                            double spanParameter) const
    {
        const gp_Pnt start =
            sourceControlPoint(
                source,
                source.panelIndex,
                pointIndex,
                47);
        const gp_Pnt end =
            sourceControlPoint(
                source,
                source.panelIndex - 1,
                pointIndex,
                47);
        // The circular expression below contains a large-radius subtraction
        // for lightly shaped stations. Its mathematical endpoints are these
        // captured rib points, so preserve them exactly: simulation skin,
        // rib, seam, and opening topology must share one vertex identity.
        if (spanParameter <= 0.0) {
            return start;
        }
        if (spanParameter >= 1.0) {
            return end;
        }
        const gp_Vec spanVector(start, end);
        const double spanLength = spanVector.Magnitude();
        if (spanLength <= Precision::Confusion()) {
            return start;
        }

        gp_Vec spanDirection = spanVector;
        spanDirection.Normalize();

        const double height = sourceShapingHeight(source, pointIndex);
        double alongSpan = spanParameter * spanLength;
        double normalOffset = 0.0;

        // This is the analytical circular ballooning law in tessella_. It is
        // evaluated directly from the transformed airfoil stations and
        // shaping height; no STL/DXF tessellation is involved.
        if (height >= 0.01 * millimetresPerCentimetre) {
            const double q =
                spanLength * spanLength / (8.0 * height)
                - 0.5 * height;
            const double radius = q + height;
            const double theta =
                std::atan(q / (0.5 * spanLength));
            const double omega =
                std::numbers::pi - 2.0 * theta;
            const double alpha =
                theta + omega * spanParameter;
            alongSpan =
                0.5 * spanLength - radius * std::cos(alpha);
            normalOffset =
                radius * std::sin(alpha) - q;
        }

        gp_Vec offset = spanDirection.Multiplied(alongSpan);
        if (normalOffset != 0.0) {
            const gp_Pnt normalStart =
                sourceControlPoint(
                    source,
                    source.panelIndex,
                    pointIndex,
                    48);
            const gp_Pnt normalEnd =
                sourceControlPoint(
                    source,
                    source.panelIndex,
                    pointIndex,
                    49);
            gp_Vec normalDirection(normalStart, normalEnd);
            if (normalDirection.Magnitude()
                <= Precision::Confusion()) {
                throw Standard_ConstructionError(
                    "Zero source shaping direction");
            }
            normalDirection.Normalize();
            offset += normalDirection.Multiplied(normalOffset);
        }
        return start.Translated(offset);
    }

    occ::handle<Geom_BSplineCurve> makeSourceSpanCurve(
        const SourcePanel &source,
        int pointIndex) const
    {
        const gp_Pnt start =
            sourceControlPoint(
                source,
                source.panelIndex,
                pointIndex,
                47);
        const gp_Pnt end =
            sourceControlPoint(
                source,
                source.panelIndex - 1,
                pointIndex,
                47);
        const gp_Vec spanVector(start, end);
        const double spanLength = spanVector.Magnitude();
        if (spanLength <= Precision::Confusion()) {
            return {};
        }

        const double height = sourceShapingHeight(source, pointIndex);
        if (height < 0.01 * millimetresPerCentimetre) {
            return makeLinearSpline(start, end);
        }

        gp_Vec spanDirection = spanVector;
        spanDirection.Normalize();

        const gp_Pnt normalStart =
            sourceControlPoint(
                source,
                source.panelIndex,
                pointIndex,
                48);
        const gp_Pnt normalEnd =
            sourceControlPoint(
                source,
                source.panelIndex,
                pointIndex,
                49);
        gp_Vec normalDirection(normalStart, normalEnd);
        if (normalDirection.Magnitude() <= Precision::Confusion()) {
            return {};
        }
        normalDirection.Normalize();

        // A circle converted to a rational B-spline is exact. Transforming
        // its poles by the source span/shaping basis also exactly preserves
        // the source law if that local basis is slightly skewed.
        GC_MakeArcOfCircle makeArc(
            gp_Pnt(0.0, 0.0, 0.0),
            gp_Pnt(0.5 * spanLength, 0.0, height),
            gp_Pnt(spanLength, 0.0, 0.0));
        if (!makeArc.IsDone()) {
            return {};
        }
        occ::handle<Geom_BSplineCurve> curve =
            GeomConvert::CurveToBSplineCurve(
                makeArc.Value(),
                Convert_QuasiAngular);
        if (curve.IsNull()) {
            return {};
        }

        // The conversion leaves an arbitrary uniform scale on the arc's
        // weights (scaling every weight by one constant is the identity on
        // a rational curve). Pin the end weights to exactly 1 so the arc
        // sections agree with the unit-weight straight sections at the
        // strip boundaries and the loft's chordwise weight function has no
        // step where straight and ballooned spans meet.
        const double endWeight = curve->Weight(1);
        if (endWeight > 0.0) {
            for (int poleIndex = 1;
                 poleIndex <= curve->NbPoles();
                 ++poleIndex) {
                curve->SetWeight(
                    poleIndex,
                    curve->Weight(poleIndex) / endWeight);
            }
        }

        for (int poleIndex = 1;
             poleIndex <= curve->NbPoles();
             ++poleIndex) {
            const gp_Pnt localPole = curve->Pole(poleIndex);
            gp_Vec offset =
                spanDirection.Multiplied(localPole.X());
            offset += normalDirection.Multiplied(localPole.Z());
            curve->SetPole(
                poleIndex,
                start.Translated(offset));
        }
        return curve;
    }

    // Skins the span sections into a tensor-product surface that contains
    // every section exactly: the sections are merged onto one shared
    // degree and knot structure (degree elevation and knot insertion are
    // exact), and their homogeneous poles are interpolated across the
    // chord at the given parameters. GeomFill_NSections used to build
    // this surface with a global least-squares fit in homogeneous space;
    // where degree-1 straight spans meet rational ballooning arcs (at the
    // trailing edge) its fitted weight function rang by about 1.5%,
    // growing spanwise fins that scale with the distance from the model
    // origin — worst at the wing tips and invisible to sampling at the
    // stations themselves. Interpolated poles stay bounded by the pole
    // data, so no such fins can appear.
    static occ::handle<Geom_BSplineSurface> skinSections(
        const std::vector<occ::handle<Geom_BSplineCurve>> &sections,
        const NCollection_Sequence<double> &sectionParameters)
    {
        const int sectionCount = static_cast<int>(sections.size());
        if (sectionCount < 2) {
            return {};
        }

        int uDegree = 1;
        bool rational = false;
        for (const occ::handle<Geom_BSplineCurve> &section : sections) {
            uDegree = std::max(uDegree, section->Degree());
            rational = rational || section->IsRational();
        }

        // Unify to the common degree and the union of the normalized
        // interior knots at each knot's maximum multiplicity.
        constexpr double knotTolerance = 1.0e-9;
        std::vector<occ::handle<Geom_BSplineCurve>> unified;
        unified.reserve(sections.size());
        std::vector<std::pair<double, int>> interiorKnots;
        for (const occ::handle<Geom_BSplineCurve> &section : sections) {
            occ::handle<Geom_BSplineCurve> copy =
                occ::handle<Geom_BSplineCurve>::DownCast(section->Copy());
            copy->IncreaseDegree(uDegree);
            const double first = copy->FirstParameter();
            const double range = copy->LastParameter() - first;
            for (int knotIndex = 2;
                 knotIndex < copy->NbKnots();
                 ++knotIndex) {
                const double knot =
                    (copy->Knot(knotIndex) - first) / range;
                const int multiplicity = copy->Multiplicity(knotIndex);
                const auto match = std::find_if(
                    interiorKnots.begin(),
                    interiorKnots.end(),
                    [&](const std::pair<double, int> &entry) {
                        return std::abs(entry.first - knot)
                               <= knotTolerance;
                    });
                if (match == interiorKnots.end()) {
                    interiorKnots.emplace_back(knot, multiplicity);
                } else {
                    match->second =
                        std::max(match->second, multiplicity);
                }
            }
            unified.push_back(copy);
        }
        std::sort(interiorKnots.begin(), interiorKnots.end());
        for (const occ::handle<Geom_BSplineCurve> &curve : unified) {
            const double first = curve->FirstParameter();
            const double range = curve->LastParameter() - first;
            for (const auto &[knot, multiplicity] : interiorKnots) {
                curve->InsertKnot(first + knot * range,
                                  multiplicity,
                                  knotTolerance * range);
            }
        }
        const int uPoleCount = unified.front()->NbPoles();
        for (const occ::handle<Geom_BSplineCurve> &curve : unified) {
            if (curve->NbPoles() != uPoleCount) {
                return {};
            }
        }

        // Chordwise interpolation: clamped knots averaged from the
        // parameters, satisfying Schoenberg-Whitney so the system below
        // is solvable.
        const int vDegree = std::min(3, sectionCount - 1);
        NCollection_Array1<double> parameters(1, sectionCount);
        for (int index = 1; index <= sectionCount; ++index) {
            parameters.SetValue(index, sectionParameters.Value(index));
        }
        NCollection_Array1<double> flatKnots(1, sectionCount + vDegree + 1);
        for (int index = 1; index <= vDegree + 1; ++index) {
            flatKnots.SetValue(index, parameters.Value(1));
            flatKnots.SetValue(sectionCount + index,
                               parameters.Value(sectionCount));
        }
        for (int index = 1;
             index <= sectionCount - vDegree - 1;
             ++index) {
            double sum = 0.0;
            for (int offset = 1; offset <= vDegree; ++offset) {
                sum += parameters.Value(index + offset);
            }
            flatKnots.SetValue(vDegree + 1 + index, sum / vDegree);
        }

        // One interpolation solve over all pole rows at once, on
        // homogeneous (weighted) coordinates so rational sections are
        // reproduced exactly.
        const int pointDimension = rational ? 4 : 3;
        const int rowDimension = uPoleCount * pointDimension;
        NCollection_Array1<double> data(1, sectionCount * rowDimension);
        for (int sectionIndex = 1;
             sectionIndex <= sectionCount;
             ++sectionIndex) {
            const occ::handle<Geom_BSplineCurve> &curve =
                unified[static_cast<std::size_t>(sectionIndex - 1)];
            int slot = (sectionIndex - 1) * rowDimension + 1;
            for (int poleIndex = 1;
                 poleIndex <= uPoleCount;
                 ++poleIndex) {
                const gp_Pnt pole = curve->Pole(poleIndex);
                const double weight =
                    curve->IsRational() ? curve->Weight(poleIndex) : 1.0;
                data.SetValue(slot++, pole.X() * weight);
                data.SetValue(slot++, pole.Y() * weight);
                data.SetValue(slot++, pole.Z() * weight);
                if (rational) {
                    data.SetValue(slot++, weight);
                }
            }
        }
        NCollection_Array1<int> contactOrder(1, sectionCount);
        contactOrder.Init(0);
        int inversionProblem = 0;
        BSplCLib::Interpolate(vDegree,
                              flatKnots,
                              parameters,
                              contactOrder,
                              rowDimension,
                              data.ChangeValue(1),
                              inversionProblem);
        if (inversionProblem != 0) {
            return {};
        }

        NCollection_Array2<gp_Pnt> poles(1, uPoleCount, 1, sectionCount);
        NCollection_Array2<double> weights(1, uPoleCount, 1, sectionCount);
        for (int sectionIndex = 1;
             sectionIndex <= sectionCount;
             ++sectionIndex) {
            int slot = (sectionIndex - 1) * rowDimension + 1;
            for (int poleIndex = 1;
                 poleIndex <= uPoleCount;
                 ++poleIndex) {
                const double x = data.Value(slot);
                const double y = data.Value(slot + 1);
                const double z = data.Value(slot + 2);
                const double weight =
                    rational ? data.Value(slot + 3) : 1.0;
                slot += pointDimension;
                if (weight <= Precision::Confusion()) {
                    return {};
                }
                poles.SetValue(poleIndex,
                               sectionIndex,
                               gp_Pnt(x / weight, y / weight, z / weight));
                weights.SetValue(poleIndex, sectionIndex, weight);
            }
        }

        NCollection_Array1<double> uKnots(
            1, static_cast<int>(interiorKnots.size()) + 2);
        NCollection_Array1<int> uMults(
            1, static_cast<int>(interiorKnots.size()) + 2);
        uKnots.SetValue(1, 0.0);
        uMults.SetValue(1, uDegree + 1);
        for (int index = 0;
             index < static_cast<int>(interiorKnots.size());
             ++index) {
            uKnots.SetValue(index + 2, interiorKnots[index].first);
            uMults.SetValue(index + 2, interiorKnots[index].second);
        }
        uKnots.SetValue(uKnots.Upper(), 1.0);
        uMults.SetValue(uMults.Upper(), uDegree + 1);

        const int interiorVCount = sectionCount - vDegree - 1;
        NCollection_Array1<double> vKnots(1, interiorVCount + 2);
        NCollection_Array1<int> vMults(1, interiorVCount + 2);
        vKnots.SetValue(1, parameters.Value(1));
        vMults.SetValue(1, vDegree + 1);
        for (int index = 1; index <= interiorVCount; ++index) {
            vKnots.SetValue(index + 1,
                            flatKnots.Value(vDegree + 1 + index));
            vMults.SetValue(index + 1, 1);
        }
        vKnots.SetValue(vKnots.Upper(), parameters.Value(sectionCount));
        vMults.SetValue(vMults.Upper(), vDegree + 1);

        if (rational) {
            return new Geom_BSplineSurface(poles,
                                           weights,
                                           uKnots,
                                           vKnots,
                                           uMults,
                                           vMults,
                                           uDegree,
                                           vDegree);
        }
        return new Geom_BSplineSurface(poles,
                                       uKnots,
                                       vKnots,
                                       uMults,
                                       vMults,
                                       uDegree,
                                       vDegree);
    }

    void addRegion(const SourcePanel &source,
                   int panelIndex,
                   int firstPoint,
                   int lastPoint,
                   Region region)
    {
        if (lastPoint - firstPoint + 1 < 2) {
            return;
        }

        addSurface(source,
                   panelIndex,
                   firstPoint,
                   lastPoint,
                   region);
    }

    void addSurface(const SourcePanel &source,
                    int panelIndex,
                    int firstPoint,
                    int lastPoint,
                    Region region)
    {
        const char *regionName = regionDiagnosticName(region);
        const int chordPointCount = lastPoint - firstPoint + 1;

        try {
            std::vector<occ::handle<Geom_BSplineCurve>> sections;
            sections.reserve(static_cast<std::size_t>(chordPointCount));
            NCollection_Sequence<double> sectionParameters;

            for (int pointIndex = firstPoint;
                 pointIndex <= lastPoint;
                 ++pointIndex) {
                occ::handle<Geom_BSplineCurve> curve =
                    makeSourceSpanCurve(source, pointIndex);
                if (curve.IsNull()) {
                    errors_.push_back(
                        "Could not create the analytical "
                        + std::string(regionName)
                        + " span curve for panel "
                        + std::to_string(panelIndex));
                    return;
                }

                const gp_Pnt midpoint =
                    sourceShapePoint(source, pointIndex, 0.5);
                if (!isFinite(midpoint)) {
                    errors_.push_back(
                        "Non-finite " + std::string(regionName)
                        + " source point in panel "
                        + std::to_string(panelIndex));
                    return;
                }
                sections.push_back(curve);
                // Point correspondence is established by the legacy remap
                // stage. Keeping this canonical parameter on every panel
                // makes adjacent lofts share the same airfoil boundary curve,
                // allowing OCCT to sew them into real shell topology.
                sectionParameters.Append(
                    static_cast<double>(pointIndex - firstPoint));
            }

            const occ::handle<Geom_BSplineSurface> surface =
                skinSections(sections, sectionParameters);
            if (!surface.IsNull()) {
                debugDumpPanel(source,
                               surface,
                               panelIndex,
                               firstPoint,
                               lastPoint,
                               regionName);
            }
            if (surface.IsNull()) {
                errors_.push_back(
                    "Could not loft the source "
                    + std::string(regionName)
                    + " NURBS surface for panel "
                    + std::to_string(panelIndex));
                return;
            }

            validateSourceSurface(
                source,
                surface,
                sectionParameters,
                firstPoint,
                lastPoint,
                panelIndex,
                regionName);
            if (!errors_.empty()) {
                return;
            }

            BRepBuilderAPI_MakeFace makeFace(surface, pointToleranceMillimetres);
            if (!makeFace.IsDone()) {
                errors_.push_back(
                    "Could not create the "
                    + std::string(regionName)
                    + " face for panel "
                    + std::to_string(panelIndex));
                return;
            }

            const TopoDS_Face face = makeFace.Face();
            if (!BRepCheck_Analyzer(face, true).IsValid()) {
                errors_.push_back(
                    "OCCT rejected the "
                    + std::string(regionName)
                    + " face for panel "
                    + std::to_string(panelIndex));
                return;
            }

            const TopoDS_Shape mirroredFace = mirrored(face);
            if (mirroredFace.IsNull()
                || mirroredFace.ShapeType() != TopAbs_FACE) {
                errors_.push_back(
                    "Could not mirror the "
                    + std::string(regionName)
                    + " face for panel "
                    + std::to_string(panelIndex));
                return;
            }

            PanelSurface panel;
            panel.region = region;
            panel.panelIndex = panelIndex;
            panel.firstPoint = firstPoint;
            panel.lastPoint = lastPoint;
            panel.surface = surface;
            panel.rightFace = face;
            panel.leftFace = TopoDS::Face(mirroredFace);
            addIsoparametricSplines(
                panel,
                maximumDisplayedSpanSamples,
                std::min(chordPointCount, maximumDisplayedChordSamples));
            panels_.push_back(std::move(panel));
        } catch (const Standard_Failure &failure) {
            std::ostringstream message;
            message << "OCCT failed to fit the " << regionName
                    << " NURBS surface for panel " << panelIndex;
            if (failure.GetMessageString() != nullptr) {
                message << ": " << failure.GetMessageString();
            }
            errors_.push_back(message.str());
        }
    }

    // Diagnostic: set LEP_DEBUG_PANEL to a comma-separated list of panel
    // indices to dump each region's section curves (span length, shaping
    // height, pole structure, weights) and the lofted surface's rib-side
    // boundary pole row to stderr. This is how the v0.2.2 trailing-edge
    // fin artifact was isolated.
    void debugDumpPanel(const SourcePanel &source,
                        const occ::handle<Geom_BSplineSurface> &surface,
                        int panelIndex,
                        int firstPoint,
                        int lastPoint,
                        const char *regionName) const
    {
        const char *debugPanels = std::getenv("LEP_DEBUG_PANEL");
        if (debugPanels == nullptr) {
            return;
        }
        const std::string padded = "," + std::string(debugPanels) + ",";
        if (padded.find("," + std::to_string(panelIndex) + ",")
            == std::string::npos) {
            return;
        }
        std::ostringstream dump;
        dump << "=== panel " << panelIndex << " region " << regionName
             << " points " << firstPoint << ".." << lastPoint << "\n";
        for (int pointIndex = firstPoint;
             pointIndex <= lastPoint;
             ++pointIndex) {
            const gp_Pnt start =
                sourceControlPoint(source, panelIndex, pointIndex, 47);
            const gp_Pnt end =
                sourceControlPoint(source, panelIndex - 1, pointIndex, 47);
            const double spanLength = start.Distance(end);
            const double height = sourceShapingHeight(source, pointIndex);
            const occ::handle<Geom_BSplineCurve> section =
                makeSourceSpanCurve(source, pointIndex);
            dump << "j=" << pointIndex
                 << " L=" << spanLength
                 << " h=" << height;
            if (!section.IsNull()) {
                dump << " poles=" << section->NbPoles()
                     << " deg=" << section->Degree()
                     << " rational=" << (section->IsRational() ? 1 : 0)
                     << " knots=" << section->NbKnots();
                if (pointIndex - firstPoint < 8) {
                    dump << " w=[";
                    for (int poleIndex = 1;
                         poleIndex <= section->NbPoles();
                         ++poleIndex) {
                        dump << (poleIndex > 1 ? " " : "")
                             << section->Weight(poleIndex);
                    }
                    dump << "]";
                }
            }
            dump << " start=(" << start.X() << ',' << start.Y() << ','
                 << start.Z() << ")\n";
        }
        dump << "surface uPoles=" << surface->NbUPoles()
             << " vPoles=" << surface->NbVPoles()
             << " uDeg=" << surface->UDegree()
             << " vDeg=" << surface->VDegree() << "\n";
        for (int vPole = 1; vPole <= surface->NbVPoles(); ++vPole) {
            const gp_Pnt pole = surface->Pole(1, vPole);
            dump << "  bpole " << vPole << " (" << pole.X() << ','
                 << pole.Y() << ',' << pole.Z() << ") w="
                 << surface->Weight(1, vPole) << "\n";
        }
        std::cerr << dump.str();
    }

    // Samples one skin region on the exact ballooning law at a decimated
    // set of chordwise stations for the Playground simulation mesh. The
    // station selection is deterministic from the region bounds, so
    // neighbouring panels with equal point tables produce bit-identical
    // rib columns and weld in the exporter.
    void captureSimRegion(const SourcePanel &source,
                          int panelIndex,
                          int firstPoint,
                          int lastPoint,
                          int targetRows,
                          Region surface,
                          bool authoredSurface,
                          bool singleSkin)
    {
        const int stationCount = lastPoint - firstPoint + 1;
        if (stationCount < 2) {
            return;
        }
        const int stride =
            std::max(1, stationCount / std::max(2, targetRows));

        SimRegionCapture capture;
        capture.panelIndex = panelIndex;
        capture.surface = surface;
        capture.authoredSurface = authoredSurface;
        capture.singleSkin = singleSkin;
        try {
            for (int station = firstPoint;
                 station <= lastPoint;
                 station += stride) {
                // Always land exactly on the region boundary so adjacent
                // regions share their junction station.
                const int kept =
                    station + stride > lastPoint ? lastPoint : station;
                std::array<gp_Pnt, simSpanColumnCount> row;
                for (int column = 0; column < simSpanColumnCount; ++column) {
                    const double spanParameter =
                        static_cast<double>(column)
                        / (simSpanColumnCount - 1);
                    row[column] =
                        sourceShapePoint(source, kept, spanParameter);
                    if (!isFinite(row[column])) {
                        return;
                    }
                }
                capture.stations.push_back(kept);
                capture.rows.push_back(row);
                if (kept == lastPoint) {
                    break;
                }
            }
        } catch (const Standard_Failure &) {
            return;
        }
        if (capture.rows.size() >= 2) {
            simRegions_.push_back(std::move(capture));
        }
    }

    void validateSourceSurface(
        const SourcePanel &source,
        const occ::handle<Geom_BSplineSurface> &surface,
        const NCollection_Sequence<double> &sectionParameters,
        int firstPoint,
        int lastPoint,
        int panelIndex,
        const char *regionName)
    {
        double uFirst = 0.0;
        double uLast = 0.0;
        double vFirst = 0.0;
        double vLast = 0.0;
        surface->Bounds(uFirst, uLast, vFirst, vLast);

        double regionSourceDeviation = 0.0;
        constexpr int validationIntervals = 8;
        for (int pointIndex = firstPoint;
             pointIndex <= lastPoint;
             ++pointIndex) {
            const std::size_t curveIndex =
                static_cast<std::size_t>(pointIndex - firstPoint);
            const double vParameter =
                sectionParameters.Value(
                    static_cast<int>(curveIndex) + 1);
            if (vParameter < vFirst - Precision::PConfusion()
                || vParameter > vLast + Precision::PConfusion()) {
                errors_.push_back(
                    "Invalid chord parameter in the "
                    + std::string(regionName)
                    + " source loft for panel "
                    + std::to_string(panelIndex));
                return;
            }

            for (int sample = 0;
                 sample <= validationIntervals;
                 ++sample) {
                const double spanParameter =
                    static_cast<double>(sample)
                    / static_cast<double>(validationIntervals);
                const gp_Pnt expected =
                    sourceShapePoint(
                        source,
                        pointIndex,
                        spanParameter);
                const gp_Pnt actual =
                    surface->Value(
                        uFirst
                            + (uLast - uFirst) * spanParameter,
                        vParameter);
                regionSourceDeviation =
                    std::max(
                        regionSourceDeviation,
                        expected.Distance(actual));
            }
        }
        maximumSourceDeviation_ =
            std::max(
                maximumSourceDeviation_,
                regionSourceDeviation);
        if (regionSourceDeviation
            > maximumSourceDeviationMillimetres) {
            std::ostringstream message;
            message << "The " << regionName
                    << " NURBS loft for panel " << panelIndex
                    << " deviates from its analytical source curves by "
                    << regionSourceDeviation << " mm";
            errors_.push_back(message.str());
            return;
        }

        // The loft is law-constrained only AT the chord stations; between
        // them it is free, which is exactly where an ill-behaved fit can
        // grow artifacts that station sampling cannot see. Bound the
        // mid-station surface against the midpoint of the neighbouring law
        // points: legitimate inter-station behaviour stays within one
        // station spacing plus the local ballooning-height crease (the law
        // may drop its shaping height discontinuously, e.g. at an air
        // intake lip, and a smooth loft must absorb that step locally).
        const auto effectiveHeight = [&](int pointIndex) {
            const double height = sourceShapingHeight(source, pointIndex);
            return height >= 0.01 * millimetresPerCentimetre ? height : 0.0;
        };
        for (int pointIndex = firstPoint;
             pointIndex < lastPoint;
             ++pointIndex) {
            const int curveIndex = pointIndex - firstPoint;
            const double vMid =
                0.5
                * (sectionParameters.Value(curveIndex + 1)
                   + sectionParameters.Value(curveIndex + 2));
            // A crease rings the interpolant over its neighbouring
            // intervals too, so measure the largest height step in a
            // one-interval neighbourhood.
            double heightStep = 0.0;
            for (int neighbour = std::max(firstPoint, pointIndex - 1);
                 neighbour <= std::min(lastPoint - 1, pointIndex + 1);
                 ++neighbour) {
                heightStep =
                    std::max(heightStep,
                             std::abs(effectiveHeight(neighbour)
                                      - effectiveHeight(neighbour + 1)));
            }
            for (int sample = 0; sample <= 4; ++sample) {
                const double spanParameter = sample / 4.0;
                const gp_Pnt nearLaw =
                    sourceShapePoint(source, pointIndex, spanParameter);
                const gp_Pnt farLaw =
                    sourceShapePoint(source, pointIndex + 1, spanParameter);
                const gp_Pnt midpoint(
                    (nearLaw.XYZ() + farLaw.XYZ()) * 0.5);
                const gp_Pnt actual =
                    surface->Value(
                        uFirst + (uLast - uFirst) * spanParameter,
                        vMid);
                const double allowance =
                    nearLaw.Distance(farLaw) + heightStep + 1.0;
                if (actual.Distance(midpoint) > allowance) {
                    std::ostringstream message;
                    message << "The " << regionName
                            << " NURBS loft for panel " << panelIndex
                            << " bulges "
                            << actual.Distance(midpoint)
                            << " mm between chord stations "
                            << pointIndex << " and " << pointIndex + 1;
                    errors_.push_back(message.str());
                    return;
                }
            }
        }

        double legacyDeviation = 0.0;
        for (int pointIndex = firstPoint;
             pointIndex <= lastPoint;
             ++pointIndex) {
            for (int segment = 0;
                 segment <= source.segmentCount;
                 ++segment) {
                const double spanParameter =
                    static_cast<double>(segment)
                    / static_cast<double>(source.segmentCount);
                const gp_Pnt expected =
                    sourceShapePoint(
                        source,
                        pointIndex,
                        spanParameter);
                const gp_Pnt legacy =
                    tessellationPoint(
                        source.legacyTessellation,
                        source.panelIndex,
                        pointIndex,
                        segment + 1);
                legacyDeviation =
                    std::max(
                        legacyDeviation,
                        expected.Distance(legacy));
            }
        }
        maximumLegacyAgreement_ =
            std::max(
                maximumLegacyAgreement_,
                legacyDeviation);
        if (legacyDeviation
            > maximumLegacyAgreementMillimetres) {
            std::ostringstream message;
            message << "The interpreted " << regionName
                    << " source shape for panel " << panelIndex
                    << " differs from the legacy validation grid by "
                    << legacyDeviation << " mm";
            errors_.push_back(message.str());
        }
    }

    void addIsoparametricSplines(
        PanelSurface &panel,
        int spanPointCount,
        int chordPointCount)
    {
        double uFirst = 0.0;
        double uLast = 0.0;
        double vFirst = 0.0;
        double vLast = 0.0;
        panel.surface->Bounds(uFirst, uLast, vFirst, vLast);

        auto addCurve = [&panel](const occ::handle<Geom_Curve> &curve) {
            if (curve.IsNull()) {
                return;
            }
            BRepBuilderAPI_MakeEdge makeEdge(curve);
            if (makeEdge.IsDone()) {
                const TopoDS_Edge edge = makeEdge.Edge();
                panel.rightWireframe.push_back(edge);
                const TopoDS_Shape mirroredEdge = mirrored(edge);
                if (!mirroredEdge.IsNull()
                    && mirroredEdge.ShapeType() == TopAbs_EDGE) {
                    panel.leftWireframe.push_back(
                        TopoDS::Edge(mirroredEdge));
                }
            }
        };

        // Index 0 and the final index coincide with the face boundaries.
        for (int index = 1; index + 1 < spanPointCount; ++index) {
            const double parameter =
                uFirst
                + (uLast - uFirst)
                      * static_cast<double>(index)
                      / static_cast<double>(spanPointCount - 1);
            addCurve(panel.surface->UIso(parameter));
        }
        for (int index = 1; index + 1 < chordPointCount; ++index) {
            const double parameter =
                vFirst
                + (vLast - vFirst)
                      * static_cast<double>(index)
                      / static_cast<double>(chordPointCount - 1);
            addCurve(panel.surface->VIso(parameter));
        }
    }

    static occ::handle<Geom_BSplineCurve> makeLinearSpline(
        const gp_Pnt &start,
        const gp_Pnt &end)
    {
        NCollection_Array1<gp_Pnt> poles(1, 2);
        poles.SetValue(1, start);
        poles.SetValue(2, end);

        NCollection_Array1<double> knots(1, 2);
        knots.SetValue(1, 0.0);
        knots.SetValue(2, 1.0);

        NCollection_Array1<int> multiplicities(1, 2);
        multiplicities.SetValue(1, 2);
        multiplicities.SetValue(2, 2);

        return new Geom_BSplineCurve(
            poles,
            knots,
            multiplicities,
            1,
            false);
    }

    // Sewing remains the topology-quality gate for the lofted skins even
    // though the export is now a structured assembly of individual panels:
    // adjacent panels must still meet on shared boundary curves.
    void validateSewing(lep::NurbsWriteResult &result) const
    {
        BRepBuilderAPI_Sewing sewing(
            pointToleranceMillimetres,
            true,
            true,
            true,
            false);
        for (const PanelSurface &panel : panels_) {
            sewing.Add(panel.rightFace);
            sewing.Add(panel.leftFace);
        }
        sewing.Perform();
        const TopoDS_Shape skins = sewing.SewedShape();
        if (skins.IsNull()
            || !BRepCheck_Analyzer(skins, true).IsValid()) {
            result.error =
                "OCCT could not sew the NURBS faces into valid skin topology.";
            return;
        }
        result.sewnEdgeCount = sewing.NbContigousEdges();
        result.freeEdgeCount = sewing.NbFreeEdges();
    }

    static TopoDS_Compound makeCompound(
        const std::vector<TopoDS_Shape> &shapes)
    {
        BRep_Builder builder;
        TopoDS_Compound compound;
        builder.MakeCompound(compound);
        for (const TopoDS_Shape &shape : shapes) {
            builder.Add(compound, shape);
        }
        return compound;
    }

    void addPanelParts(AssemblyGroup &wing,
                       lep::NurbsWriteResult &result,
                       bool includeConstructionCurves) const
    {
        for (const PanelSurface &panel : panels_) {
            AssemblyGroup &regionGroup =
                wing.group(regionGroupName(panel.region));
            const std::string partName =
                "Panel " + std::to_string(panel.panelIndex);
            const PartColor faceColor = regionColor(panel.region);

            regionGroup.group(rightSideName).parts.push_back(
                {partName,
                 panel.rightFace,
                 faceColor,
                 wireframeColor,
                 true});
            regionGroup.group(leftSideName).parts.push_back(
                {partName,
                 panel.leftFace,
                 faceColor,
                 wireframeColor,
                 true});
            result.surfaceCount += 2;

            if (!includeConstructionCurves) {
                continue;
            }
            AssemblyGroup &curvesGroup =
                wing.group(regionCurvesGroupName(panel.region));
            const auto addWireframePart =
                [&](const char *sideName,
                    const std::vector<TopoDS_Edge> &wireframe) {
                    if (wireframe.empty()) {
                        return;
                    }
                    const std::vector<TopoDS_Shape> shapes(
                        wireframe.begin(), wireframe.end());
                    curvesGroup.group(sideName).parts.push_back(
                        {partName,
                         makeCompound(shapes),
                         wireframeColor,
                         wireframeColor,
                         false});
                };
            addWireframePart(rightSideName, panel.rightWireframe);
            addWireframePart(leftSideName, panel.leftWireframe);
            result.splineCount += static_cast<int>(
                panel.rightWireframe.size() + panel.leftWireframe.size());
        }
    }

    // The exact chordwise rib contour over [firstPoint, lastPoint], as one
    // degree-1 B-spline through the captured station points. The legacy
    // shape definition is chordwise piecewise linear (the cut patterns and
    // the 3D reference output are polylines through these points), and a
    // smoothing interpolation can overshoot into the opposite side of a
    // thin profile, so the polyline is both the faithful and the robust
    // boundary.
    TopoDS_Edge contourEdge(const CapturedRib &rib,
                            int firstPoint,
                            int lastPoint) const
    {
        const int count = lastPoint - firstPoint + 1;
        NCollection_Array1<gp_Pnt> poles(1, count);
        NCollection_Array1<double> knots(1, count);
        NCollection_Array1<int> multiplicities(1, count);
        for (int index = 0; index < count; ++index) {
            poles.SetValue(
                index + 1,
                rib.spatialPoints[
                    static_cast<std::size_t>(firstPoint - 1 + index)]);
            multiplicities.SetValue(index + 1, 1);
        }
        multiplicities.SetValue(1, 2);
        multiplicities.SetValue(count, 2);
        knots.SetValue(1, 0.0);
        for (int index = 2; index <= count; ++index) {
            knots.SetValue(
                index,
                knots.Value(index - 1)
                    + std::max(
                        poles.Value(index - 1)
                            .Distance(poles.Value(index)),
                        Precision::Confusion()));
        }

        BRepBuilderAPI_MakeEdge makeEdge(
            occ::handle<Geom_BSplineCurve>(
                new Geom_BSplineCurve(
                    poles,
                    knots,
                    multiplicities,
                    1,
                    false)));
        if (!makeEdge.IsDone()) {
            return {};
        }
        return makeEdge.Edge();
    }

    // Builds one closed rounded-polygon hole wire (types 3 and 4) from the
    // legacy corner construction: circle centres sit on the corner
    // bisectors, so the straight sides are the common tangents between
    // consecutive corner circles.
    TopoDS_Wire roundedPolygonHoleWire(const std::vector<gp_XY> &corners,
                                       const std::vector<gp_XY> &centers,
                                       double radius,
                                       const RibFrame &frame,
                                       int &edgeCount) const
    {
        const std::size_t count = corners.size();
        BRepBuilderAPI_MakeWire makeWire;
        int added = 0;
        const auto addEdge = [&makeWire, &added](const TopoDS_Edge &edge) {
            if (edge.IsNull()) {
                return false;
            }
            makeWire.Add(edge);
            ++added;
            return makeWire.IsDone() == Standard_True;
        };
        const auto lineEdge = [](const gp_Pnt &from,
                                 const gp_Pnt &to) -> TopoDS_Edge {
            BRepBuilderAPI_MakeEdge makeEdge(makeLinearSpline(from, to));
            if (!makeEdge.IsDone()) {
                return {};
            }
            return makeEdge.Edge();
        };

        if (radius <= Precision::Confusion()) {
            for (std::size_t index = 0; index < count; ++index) {
                const gp_Pnt from = frame.point(corners[index]);
                const gp_Pnt to =
                    frame.point(corners[(index + 1) % count]);
                if (from.Distance(to) <= Precision::Confusion()
                    || !addEdge(lineEdge(from, to))) {
                    return {};
                }
            }
        } else {
            std::vector<gp_XY> tangentIn(count);
            std::vector<gp_XY> tangentOut(count);
            for (std::size_t index = 0; index < count; ++index) {
                const gp_XY &corner = corners[index];
                tangentIn[index] = projectOntoLine(
                    corners[(index + count - 1) % count],
                    corner,
                    centers[index]);
                tangentOut[index] = projectOntoLine(
                    corner,
                    corners[(index + 1) % count],
                    centers[index]);
                if (std::abs(
                        (tangentIn[index] - centers[index]).Modulus()
                        - radius) > pointToleranceMillimetres
                    || std::abs(
                        (tangentOut[index] - centers[index]).Modulus()
                        - radius) > pointToleranceMillimetres) {
                    return {};
                }
            }
            for (std::size_t index = 0; index < count; ++index) {
                if ((tangentOut[index] - tangentIn[index]).Modulus()
                    > Precision::Confusion()) {
                    const gp_XY bisector =
                        corners[index] - centers[index];
                    const double bisectorLength = bisector.Modulus();
                    if (bisectorLength <= Precision::Confusion()) {
                        return {};
                    }
                    const gp_XY arcMid =
                        centers[index]
                        + bisector * (radius / bisectorLength);
                    GC_MakeArcOfCircle makeArc(
                        frame.point(tangentIn[index]),
                        frame.point(arcMid),
                        frame.point(tangentOut[index]));
                    if (!makeArc.IsDone()) {
                        return {};
                    }
                    BRepBuilderAPI_MakeEdge makeEdge(makeArc.Value());
                    if (!makeEdge.IsDone()
                        || !addEdge(makeEdge.Edge())) {
                        return {};
                    }
                }
                const gp_Pnt from = frame.point(tangentOut[index]);
                const gp_Pnt to =
                    frame.point(tangentIn[(index + 1) % count]);
                if (from.Distance(to) > Precision::Confusion()
                    && !addEdge(lineEdge(from, to))) {
                    return {};
                }
            }
        }
        if (added < 2 || !makeWire.IsDone()) {
            return {};
        }
        const TopoDS_Wire wire = makeWire.Wire();
        if (!wire.Closed()) {
            return {};
        }
        edgeCount = added;
        return wire;
    }

    // Builds the closed wire of one airfoil hole in the rib plane. The
    // shapes reproduce the legacy 2D construction exactly, but as ellipse,
    // line and arc geometry instead of the drawing polylines.
    TopoDS_Wire ribHoleWire(const RibHole &hole,
                            const RibFrame &frame,
                            int ribIndex,
                            int holeNumber,
                            int &edgeCount,
                            std::vector<std::string> &warnings) const
    {
        // A hole that cannot be built is left uncut with a warning; the
        // legacy core also draws inconsistent hole definitions without
        // complaint.
        const auto fail =
            [&warnings, ribIndex, holeNumber](const std::string &reason) {
                warnings.push_back(
                    "Left hole " + std::to_string(holeNumber)
                    + " of rib " + std::to_string(ribIndex)
                    + " uncut: " + reason);
                return TopoDS_Wire();
            };

        if (hole.type == 1) {
            const double semiA = std::abs(hole.a);
            const double semiB = std::abs(hole.b);
            if (semiA <= Precision::Confusion()
                || semiB <= Precision::Confusion()) {
                return fail("Degenerate ellipse axes");
            }
            // The legacy core rotates ellipses in the y-flipped drawing
            // sheet, taking the file value directly as radians; in the rib
            // frame the axis directions therefore come out as below.
            const double rotation = hole.rotation;
            gp_Vec direction1 = frame.direction(
                {std::cos(rotation), -std::sin(rotation)});
            gp_Vec direction2 = frame.direction(
                {-std::sin(rotation), -std::cos(rotation)});
            double major = semiA;
            double minor = semiB;
            if (major < minor) {
                std::swap(major, minor);
                std::swap(direction1, direction2);
            }
            const gp_Ax2 axes(
                frame.point({hole.x, hole.y}),
                gp_Dir(direction1.Crossed(direction2)),
                gp_Dir(direction1));
            BRepBuilderAPI_MakeEdge makeEdge(gp_Elips(axes, major, minor));
            if (!makeEdge.IsDone()) {
                return fail("Could not build the ellipse");
            }
            BRepBuilderAPI_MakeWire makeWire(makeEdge.Edge());
            if (!makeWire.IsDone()) {
                return fail("Could not close the ellipse");
            }
            edgeCount = 1;
            return makeWire.Wire();
        }

        // Types 3 (rounded triangle) and 4 (rounded rectangle) rotate in
        // degrees and mirror through the sign of the first side length.
        const double alpha =
            hole.rotation * std::numbers::pi / 180.0;
        if (std::abs(hole.a) <= Precision::Confusion()
            || hole.b <= Precision::Confusion()
            || hole.cornerRadius < 0.0) {
            return fail("Degenerate hole sides");
        }
        const double sideSign = hole.a >= 0.0 ? 1.0 : -1.0;
        const double sideA = std::abs(hole.a);
        const double sideB = hole.b;
        const double radius = hole.cornerRadius;

        std::vector<gp_XY> corners;
        std::vector<gp_XY> centers;
        const gp_XY corner1(hole.x, hole.y);
        if (hole.type == 3) {
            const gp_XY corner2(
                hole.x + sideSign * sideA * std::cos(alpha),
                hole.y + sideA * std::sin(alpha));
            const gp_XY corner3(hole.x, hole.y + sideB);
            const double gammaC = 0.5 * std::numbers::pi - alpha;
            const double sideC = std::sqrt(std::max(
                0.0,
                sideA * sideA + sideB * sideB
                    - 2.0 * sideA * sideB * std::cos(gammaC)));
            if (sideC <= Precision::Confusion()) {
                return fail("Degenerate triangle");
            }
            const double angleA = std::acos(std::clamp(
                (sideC * sideC + sideB * sideB - sideA * sideA)
                    / (2.0 * sideB * sideC),
                -1.0,
                1.0));
            const double angleB = std::acos(std::clamp(
                (sideC * sideC + sideA * sideA - sideB * sideB)
                    / (2.0 * sideA * sideC),
                -1.0,
                1.0));
            double height1 = 0.0;
            double height2 = 0.0;
            double height3 = 0.0;
            if (radius > Precision::Confusion()) {
                const double sin1 = std::sin(0.5 * gammaC);
                const double sin2 = std::sin(0.5 * angleB);
                const double sin3 = std::sin(0.5 * angleA);
                if (sin1 <= Precision::Confusion()
                    || sin2 <= Precision::Confusion()
                    || sin3 <= Precision::Confusion()) {
                    return fail("Invalid corner rounding");
                }
                height1 = radius / sin1;
                height2 = radius / sin2;
                height3 = radius / sin3;
            }
            corners = {corner1, corner2, corner3};
            centers = {
                {hole.x
                     + sideSign * height1 * std::cos(alpha + 0.5 * gammaC),
                 hole.y + height1 * std::sin(alpha + 0.5 * gammaC)},
                {corner2.X()
                     - sideSign * height2 * std::cos(-alpha + 0.5 * angleB),
                 corner2.Y() + height2 * std::sin(-alpha + 0.5 * angleB)},
                {corner3.X() + sideSign * height3 * std::sin(0.5 * angleA),
                 corner3.Y() - height3 * std::cos(0.5 * angleA)},
            };
        } else {
            const gp_XY corner2(
                hole.x + sideSign * sideA * std::cos(alpha),
                hole.y + sideA * std::sin(alpha));
            const gp_XY corner3(
                corner2.X() - sideSign * sideB * std::sin(alpha),
                corner2.Y() + sideB * std::cos(alpha));
            const gp_XY corner4(
                hole.x - sideSign * sideB * std::sin(alpha),
                hole.y + sideB * std::cos(alpha));
            const double radiusCos = radius * std::cos(alpha);
            const double radiusSin = radius * std::sin(alpha);
            corners = {corner1, corner2, corner3, corner4};
            centers = {
                {corner1.X() + sideSign * (radiusCos - radiusSin),
                 corner1.Y() + radiusSin + radiusCos},
                {corner2.X() - sideSign * (radiusCos + radiusSin),
                 corner2.Y() - radiusSin + radiusCos},
                {corner3.X() - sideSign * (radiusCos - radiusSin),
                 corner3.Y() - radiusSin - radiusCos},
                {corner4.X() + sideSign * (radiusCos + radiusSin),
                 corner4.Y() + radiusSin - radiusCos},
            };
        }

        const TopoDS_Wire wire = roundedPolygonHoleWire(
            corners, centers, radius, frame, edgeCount);
        if (wire.IsNull()) {
            return fail("Could not build the outline");
        }
        return wire;
    }

    // Builds the closed planar face of one rib. The outer boundary reuses
    // the exact panel boundary edges wherever a panel exists, interpolates
    // the captured contour where none does, and closes the trailing edge
    // with a straight seam. Airfoil holes become inner wires.
    // Builds one rib face; runs on worker threads, so every warning and
    // error goes to the caller's sinks instead of the shared members.
    TopoDS_Face makeRibFace(int ribIndex,
                            std::vector<RibBoundarySegment> segments,
                            const CapturedRib &rib,
                            int &edgeCount,
                            TopoDS_Shape &curveFallback,
                            std::vector<std::string> &warnings,
                            std::vector<std::string> &errors,
                            const std::set<int> *sceneContourStations = nullptr) const
    {
        const auto fail = [&errors, ribIndex](const std::string &reason) {
            errors.push_back(
                reason + " for rib " + std::to_string(ribIndex));
            return TopoDS_Face();
        };

        const int totalPointCount =
            static_cast<int>(rib.spatialPoints.size());
        std::sort(segments.begin(),
                  segments.end(),
                  [](const RibBoundarySegment &left,
                     const RibBoundarySegment &right) {
                      return left.firstPoint < right.firstPoint;
                  });

        // The rib plane uses a Newell normal so the winding test below is
        // well conditioned. Detect a collapsed wingtip before choosing the
        // scene-only coarse edge topology: its retraced zero-area
        // contour is deliberately retained as the compact curve fallback.
        gp_XYZ centroid(0.0, 0.0, 0.0);
        for (const gp_Pnt &point : rib.spatialPoints) {
            centroid += point.XYZ();
        }
        centroid /= static_cast<double>(totalPointCount);
        gp_XYZ normalAccumulator(0.0, 0.0, 0.0);
        for (int index = 0; index < totalPointCount; ++index) {
            const gp_XYZ current =
                rib.spatialPoints[
                    static_cast<std::size_t>(index)].XYZ() - centroid;
            const gp_XYZ next =
                rib.spatialPoints[
                    static_cast<std::size_t>(
                        (index + 1) % totalPointCount)].XYZ() - centroid;
            normalAccumulator += current.Crossed(next);
        }
        double contourExtent = 0.0;
        for (const gp_Pnt &point : rib.spatialPoints) {
            contourExtent = std::max(
                contourExtent,
                (point.XYZ() - centroid).Modulus());
        }
        const bool collapsedContour = normalAccumulator.Modulus()
            <= 1.0e-6 * contourExtent * contourExtent;
        const bool useSceneContour =
            sceneContourStations != nullptr && !collapsedContour;

        // CAD output interpolates through the exact rib station points rather
        // than reusing panel boundary edges. The scene-only path uses straight
        // segments between retained coarse skin stations so its rib boundary
        // shares the same finite-element topology while keeping every endpoint
        // on the authoritative captured rib contour. Both remain split at the
        // panel region boundaries (vent corners).
        std::vector<std::pair<int, int>> ranges;
        int cursor = 1;
        for (const RibBoundarySegment &segment : segments) {
            if (segment.firstPoint < cursor
                || segment.firstPoint >= segment.lastPoint
                || segment.lastPoint > totalPointCount) {
                return fail("Inconsistent panel boundary ranges");
            }
            if (segment.firstPoint > cursor) {
                ranges.emplace_back(cursor, segment.firstPoint);
            }
            ranges.emplace_back(segment.firstPoint, segment.lastPoint);
            cursor = segment.lastPoint;
        }
        if (cursor < totalPointCount) {
            ranges.emplace_back(cursor, totalPointCount);
        }
        if (ranges.empty()) {
            return fail("No outline curves");
        }

        struct BoundaryEdge
        {
            TopoDS_Edge edge;
            gp_Pnt start;
            gp_Pnt end;
        };
        std::vector<BoundaryEdge> boundary;
        for (const auto &[firstPoint, lastPoint] : ranges) {
            std::vector<int> edgeStations{firstPoint};
            if (useSceneContour) {
                for (auto station = sceneContourStations->upper_bound(
                         firstPoint);
                     station != sceneContourStations->end()
                         && *station < lastPoint;
                     ++station) {
                    edgeStations.push_back(*station);
                }
            }
            edgeStations.push_back(lastPoint);
            for (std::size_t station = 0;
                 station + 1 < edgeStations.size(); ++station) {
                const int edgeFirst = edgeStations[station];
                const int edgeLast = edgeStations[station + 1];
                TopoDS_Edge edge;
                if (useSceneContour) {
                    BRepBuilderAPI_MakeEdge makeEdge(
                        makeLinearSpline(
                            rib.spatialPoints[static_cast<std::size_t>(
                                edgeFirst - 1)],
                            rib.spatialPoints[static_cast<std::size_t>(
                                edgeLast - 1)]));
                    if (makeEdge.IsDone()) {
                        edge = makeEdge.Edge();
                    }
                } else {
                    edge = contourEdge(rib, edgeFirst, edgeLast);
                }
                if (edge.IsNull()) {
                    return fail("Could not build the outline");
                }
                boundary.push_back(
                    {edge,
                     rib.spatialPoints[
                         static_cast<std::size_t>(edgeFirst - 1)],
                     rib.spatialPoints[
                         static_cast<std::size_t>(edgeLast - 1)]});
            }
        }
        // Straight trailing-edge seam whenever the airfoil is open there.
        if (boundary.back().end.Distance(boundary.front().start)
            > Precision::Confusion()) {
            BRepBuilderAPI_MakeEdge closing(
                makeLinearSpline(
                    boundary.back().end,
                    boundary.front().start));
            if (!closing.IsDone()) {
                return fail("Could not close the trailing edge");
            }
            boundary.push_back(
                {closing.Edge(),
                 boundary.back().end,
                 boundary.front().start});
        }

        BRepBuilderAPI_MakeWire makeWire;
        for (const BoundaryEdge &edge : boundary) {
            makeWire.Add(edge.edge);
            if (!makeWire.IsDone()) {
                return fail("Could not chain the outline");
            }
        }
        TopoDS_Wire outerWire = makeWire.Wire();
        if (!outerWire.Closed()) {
            return fail("The outline is not closed");
        }

        // A collapsed rib encloses no area (the wingtip typically closes
        // the wing with a zero-thickness profile whose lower side retraces
        // the upper one). No face exists there; export the outline curves.
        if (collapsedContour) {
            std::vector<TopoDS_Shape> outline;
            outline.reserve(boundary.size());
            for (const BoundaryEdge &edge : boundary) {
                outline.push_back(edge.edge);
            }
            curveFallback = makeCompound(outline);
            edgeCount = static_cast<int>(boundary.size());
            return {};
        }
        gp_Dir normal(normalAccumulator);
        for (const gp_Pnt &point : rib.spatialPoints) {
            if (std::abs((point.XYZ() - centroid).Dot(normal.XYZ()))
                > pointToleranceMillimetres) {
                return fail("The contour is not planar");
            }
        }

        const gp_Pnt planeOrigin(centroid);
        const double outerWinding =
            signedAreaAlongNormal(outerWire, planeOrigin, normal);
        if (std::abs(outerWinding) <= Precision::SquareConfusion()) {
            return fail("Could not orient the outline");
        }
        if (outerWinding < 0.0) {
            normal.Reverse();
        }

        BRepBuilderAPI_MakeFace makeFace(
            gp_Pln(planeOrigin, normal), outerWire);
        if (!makeFace.IsDone()) {
            return fail("Could not build the face");
        }
        TopoDS_Face face = makeFace.Face();
        edgeCount = static_cast<int>(boundary.size());

        // Every face curve is built from in-plane geometry, so any residual
        // out-of-plane deviation is numerical noise. Measure it, reject
        // anything conspicuous, and record the measured bound as the edge
        // tolerance where it exceeds the model precision.
        constexpr double maximumPlaneDeviationMillimetres = 0.05;
        BRep_Builder toleranceUpdater;
        for (TopExp_Explorer explorer(face, TopAbs_EDGE);
             explorer.More();
             explorer.Next()) {
            const TopoDS_Edge &edge = TopoDS::Edge(explorer.Current());
            double first = 0.0;
            double last = 0.0;
            const occ::handle<Geom_Curve> curve =
                BRep_Tool::Curve(edge, first, last);
            if (curve.IsNull()) {
                return fail("Could not evaluate an outline curve");
            }
            double deviation = 0.0;
            constexpr int sampleCount = 64;
            for (int sample = 0; sample <= sampleCount; ++sample) {
                const double parameter =
                    first
                    + (last - first) * static_cast<double>(sample)
                          / static_cast<double>(sampleCount);
                deviation = std::max(
                    deviation,
                    std::abs(
                        (curve->Value(parameter).XYZ() - centroid)
                            .Dot(normal.XYZ())));
            }
            if (deviation > maximumPlaneDeviationMillimetres) {
                std::ostringstream message;
                message << "A face curve leaves the rib plane by "
                        << deviation << " mm";
                return fail(message.str());
            }
            if (deviation > Precision::Confusion()) {
                toleranceUpdater.UpdateEdge(edge, deviation * 1.25);
            }
        }

        const auto invalidFaceDetail = [](const BRepCheck_Analyzer &analyzer,
                                          const TopoDS_Face &checked) {
            std::ostringstream detail;
            detail << '(';
            const auto report = [&](const TopoDS_Shape &shape,
                                    const char *kind) {
                const occ::handle<BRepCheck_Result> checkResult =
                    analyzer.Result(shape);
                if (checkResult.IsNull()) {
                    return;
                }
                for (const BRepCheck_Status status :
                     checkResult->Status()) {
                    if (status != BRepCheck_NoError) {
                        detail << ' ' << kind
                               << static_cast<int>(status);
                    }
                }
                if (!checkResult->IsStatusOnShape(checked)) {
                    return;
                }
                for (const BRepCheck_Status status :
                     checkResult->StatusOnShape(checked)) {
                    if (status != BRepCheck_NoError) {
                        detail << ' ' << kind << "ctx"
                               << static_cast<int>(status);
                    }
                }
            };
            report(checked, "face:");
            for (TopExp_Explorer explorer(checked, TopAbs_WIRE);
                 explorer.More();
                 explorer.Next()) {
                report(explorer.Current(), "wire:");
            }
            for (TopExp_Explorer explorer(checked, TopAbs_EDGE);
                 explorer.More();
                 explorer.Next()) {
                report(explorer.Current(), "edge:");
            }
            detail << " )";
            return detail.str();
        };

        const auto outlineFallback = [&](const BRepCheck_Analyzer &analyzer) {
            // Very thin profiles (single-skin ribs) can produce an
            // outline OCCT cannot classify as a face; keep the exact
            // outline curves instead of failing the whole model.
            warnings.push_back(
                "Exported rib " + std::to_string(ribIndex)
                + " as outline curves: OCCT rejected its face "
                + invalidFaceDetail(analyzer, face));
            std::vector<TopoDS_Shape> outline;
            outline.reserve(boundary.size());
            for (const BoundaryEdge &edge : boundary) {
                outline.push_back(edge.edge);
            }
            curveFallback = makeCompound(outline);
            edgeCount = static_cast<int>(boundary.size());
        };

        // Prepare every hole wire up front. Their warnings are committed
        // only when the rib keeps its face: a rib that falls back to
        // outline curves never reported hole problems.
        struct PreparedHole
        {
            std::size_t number = 0;
            TopoDS_Wire wire;
            int edgeCount = 0;
        };
        std::vector<PreparedHole> preparedHoles;
        std::vector<std::string> holeWarnings;
        if (!rib.holes.empty()) {
            RibFrame frame;
            if (!fitRibFrame(rib.planarPoints, rib.spatialPoints, frame)) {
                holeWarnings.push_back(
                    "Left the holes of rib " + std::to_string(ribIndex)
                    + " uncut: could not recover its rigid planar frame");
            } else {
                for (std::size_t holeIndex = 0;
                     holeIndex < rib.holes.size();
                     ++holeIndex) {
                    int holeEdgeCount = 0;
                    TopoDS_Wire holeWire = ribHoleWire(
                        rib.holes[holeIndex],
                        frame,
                        ribIndex,
                        static_cast<int>(holeIndex) + 1,
                        holeEdgeCount,
                        holeWarnings);
                    if (holeWire.IsNull()) {
                        // ribHoleWire already recorded why.
                        continue;
                    }
                    // Holes must wind against the outer boundary.
                    if (signedAreaAlongNormal(holeWire, planeOrigin, normal)
                        > 0.0) {
                        holeWire.Reverse();
                    }
                    preparedHoles.push_back(
                        {holeIndex + 1, holeWire, holeEdgeCount});
                }
            }
        }

        // Cut every hole, then validate the finished face once: the
        // analyzer dominates the whole model build, so the per-hole
        // diagnosis below runs only when this combined face fails.
        TopoDS_Face holedFace = face;
        int holedEdgeCount = 0;
        bool allHolesCut = true;
        for (const PreparedHole &hole : preparedHoles) {
            BRepBuilderAPI_MakeFace cut(holedFace, hole.wire);
            if (!cut.IsDone()) {
                allHolesCut = false;
                break;
            }
            holedFace = cut.Face();
            holedEdgeCount += hole.edgeCount;
        }
        if (allHolesCut) {
            const BRepCheck_Analyzer analyzer(holedFace, true);
            if (analyzer.IsValid()) {
                warnings.insert(warnings.end(),
                                holeWarnings.begin(),
                                holeWarnings.end());
                edgeCount += holedEdgeCount;
                return holedFace;
            }
        }

        // Something in the combined face is invalid. Check the bare face
        // first (it decides between the outline fallback and hole trouble),
        // then re-cut hole by hole to skip exactly the offenders.
        {
            const BRepCheck_Analyzer analyzer(face, true);
            if (!analyzer.IsValid()) {
                outlineFallback(analyzer);
                return {};
            }
        }
        warnings.insert(warnings.end(),
                        holeWarnings.begin(),
                        holeWarnings.end());
        for (const PreparedHole &hole : preparedHoles) {
            BRepBuilderAPI_MakeFace cut(face, hole.wire);
            const BRepCheck_Analyzer holedCheck(
                cut.IsDone() ? cut.Face() : face, true);
            if (!cut.IsDone() || !holedCheck.IsValid()) {
                // The legacy core draws such holes over the outline in
                // the 2D patterns; in the solid model they cannot be
                // cut, so leave them out and say so.
                warnings.push_back(
                    "Left hole " + std::to_string(hole.number)
                    + " of rib " + std::to_string(ribIndex)
                    + " uncut: it does not fit inside the rib outline");
                continue;
            }
            face = cut.Face();
            edgeCount += hole.edgeCount;
        }
        return face;
    }

    // Every captured panel spans rib i to rib i-1, so panel i declares the
    // panel-covered chord ranges of rib i, and of rib i-1 only when panel
    // i-1 was not built (which also yields the centre rib 0 from panel 1).
    // Each rib station becomes a closed planar face over its full captured
    // contour, with its airfoil holes cut out.
    void addRibParts(AssemblyGroup &wing,
                     lep::NurbsWriteResult &result) const
    {
        std::unordered_set<int> capturedPanels;
        for (const PanelSurface &panel : panels_) {
            capturedPanels.insert(panel.panelIndex);
        }

        std::map<int, std::vector<RibBoundarySegment>> ribSegments;
        for (const PanelSurface &panel : panels_) {
            ribSegments[panel.panelIndex].push_back(
                {panel.firstPoint, panel.lastPoint});
            if (!capturedPanels.contains(panel.panelIndex - 1)) {
                ribSegments[panel.panelIndex - 1].push_back(
                    {panel.firstPoint, panel.lastPoint});
            }
        }

        struct RibJob
        {
            int ribIndex = 0;
            const std::vector<RibBoundarySegment> *segments = nullptr;
            const CapturedRib *rib = nullptr;
        };
        std::vector<RibJob> jobs;
        jobs.reserve(ribSegments.size());
        for (const auto &[ribIndex, segments] : ribSegments) {
            const auto captured = capturedRibs_.find(ribIndex);
            if (captured == capturedRibs_.end()) {
                errors_.push_back(
                    "Missing captured rib station for rib "
                    + std::to_string(ribIndex));
                return;
            }
            jobs.push_back({ribIndex, &segments, &captured->second});
        }

        struct RibBuild
        {
            TopoDS_Shape shape;
            TopoDS_Shape mirroredShape;
            bool hasFace = false;
            bool onCenter = false;
            int edgeCount = 0;
            std::vector<std::string> warnings;
            std::vector<std::string> errors;
        };
        std::vector<RibBuild> builds(jobs.size());

        // Every rib face is geometrically independent, so they build on
        // OCCT's thread pool; the merge below runs in rib order to keep
        // the assembly, warnings, and errors deterministic.
        OSD_Parallel::For(
            0,
            static_cast<int>(jobs.size()),
            [this, &jobs, &builds](int jobIndex) {
                const RibJob &job =
                    jobs[static_cast<std::size_t>(jobIndex)];
                RibBuild &build =
                    builds[static_cast<std::size_t>(jobIndex)];
                try {
                    TopoDS_Shape outlineCurves;
                    const TopoDS_Face face = makeRibFace(
                        job.ribIndex,
                        *job.segments,
                        *job.rib,
                        build.edgeCount,
                        outlineCurves,
                        build.warnings,
                        build.errors);
                    if (!build.errors.empty()) {
                        return;
                    }
                    build.hasFace = !face.IsNull();
                    build.shape = build.hasFace ? TopoDS_Shape(face)
                                                : outlineCurves;
                    if (build.shape.IsNull()) {
                        build.errors.push_back(
                            "Could not build the shape of rib "
                            + std::to_string(job.ribIndex));
                        return;
                    }
                    build.onCenter = std::all_of(
                        job.rib->spatialPoints.begin(),
                        job.rib->spatialPoints.end(),
                        [](const gp_Pnt &point) {
                            return std::abs(point.X())
                                   <= symmetryPlaneToleranceMillimetres;
                        });
                    if (!build.onCenter) {
                        build.mirroredShape = mirrored(build.shape);
                        if (build.mirroredShape.IsNull()) {
                            build.errors.push_back(
                                "Could not mirror the shape of rib "
                                + std::to_string(job.ribIndex));
                        }
                    }
                } catch (const Standard_Failure &failure) {
                    build.errors.push_back(
                        "OCCT failed building rib "
                        + std::to_string(job.ribIndex) + ": "
                        + (failure.GetMessageString() != nullptr
                               ? failure.GetMessageString()
                               : "unknown OCCT error"));
                } catch (const std::exception &exception) {
                    build.errors.push_back(
                        "Failed building rib "
                        + std::to_string(job.ribIndex) + ": "
                        + exception.what());
                }
            });

        AssemblyGroup &ribs = wing.group(ribsGroupName);
        for (std::size_t index = 0; index < jobs.size(); ++index) {
            RibBuild &build = builds[index];
            warnings_.insert(warnings_.end(),
                             build.warnings.begin(),
                             build.warnings.end());
            if (!build.errors.empty()) {
                errors_.insert(errors_.end(),
                               build.errors.begin(),
                               build.errors.end());
                return;
            }

            const std::string partName =
                "Rib " + std::to_string(jobs[index].ribIndex);
            if (build.onCenter) {
                ribs.group(centerSideName).parts.push_back(
                    {partName, build.shape, ribColor, ribColor,
                     build.hasFace});
                result.surfaceCount += build.hasFace ? 1 : 0;
                result.splineCount += build.edgeCount;
            } else {
                ribs.group(rightSideName).parts.push_back(
                    {partName, build.shape, ribColor, ribColor,
                     build.hasFace});
                ribs.group(leftSideName).parts.push_back(
                    {partName, build.mirroredShape, ribColor, ribColor,
                     build.hasFace});
                result.surfaceCount += build.hasFace ? 2 : 0;
                result.splineCount += build.edgeCount * 2;
            }
            ++result.ribCount;
        }
    }

    void addLineParts(AssemblyGroup &wing,
                      lep::NurbsWriteResult &result) const
    {
        if (capturedLines_.empty()) {
            return;
        }

        AssemblyGroup &lines = wing.group(linesGroupName);
        std::unordered_set<QuantizedSegment, QuantizedSegmentHash> added;

        struct LabelPart
        {
            AssemblyGroup *group = nullptr;
            std::string label;
            PartColor color;
            std::vector<TopoDS_Shape> segments;
        };
        std::vector<LabelPart> labelParts;

        const auto labelPartFor = [&](const CapturedLine &line) -> LabelPart & {
            AssemblyGroup *group = nullptr;
            PartColor color = otherCurveColor;
            std::string label = line.label;
            if (line.brake) {
                group = &lines.group(brakeGroupName);
                color = brakeColor;
            } else if (line.planIndex >= 1 && line.planIndex <= 6) {
                group = &lines.group(
                    std::string("Plan ")
                    + static_cast<char>('A' + line.planIndex - 1));
                color = planColors[line.planIndex - 1];
            } else {
                group = &lines.group(otherCurvesName);
            }
            if (label.empty()) {
                label = "Curve";
            }
            for (LabelPart &part : labelParts) {
                if (part.group == group && part.label == label) {
                    return part;
                }
            }
            labelParts.push_back({group, label, color, {}});
            return labelParts.back();
        };

        for (const CapturedLine &line : capturedLines_) {
            const QuantizedSegment key =
                quantizeSegment(line.start, line.end);
            if (!added.insert(key).second) {
                continue;
            }
            BRepBuilderAPI_MakeEdge makeEdge(
                makeLinearSpline(line.start, line.end));
            if (!makeEdge.IsDone()) {
                continue;
            }
            labelPartFor(line).segments.push_back(makeEdge.Edge());
        }

        for (const LabelPart &part : labelParts) {
            if (part.segments.empty()) {
                continue;
            }
            part.group->parts.push_back(
                {part.label,
                 makeCompound(part.segments),
                 part.color,
                 part.color,
                 false});
            ++result.lineCount;
            result.splineCount += static_cast<int>(part.segments.size());
        }
    }

    // Interpolates one strip boundary through its samples at the given
    // parameters. Identical parameters across sibling boundaries yield
    // identical knot vectors, which makeStripFace relies on to pair poles.
    static occ::handle<Geom_BSplineCurve> interpolateBoundary(
        const std::vector<gp_Pnt> &points,
        const std::vector<double> &parameters)
    {
        const int count = static_cast<int>(points.size());
        occ::handle<NCollection_HArray1<gp_Pnt>> pointArray =
            new NCollection_HArray1<gp_Pnt>(1, count);
        occ::handle<NCollection_HArray1<double>> parameterArray =
            new NCollection_HArray1<double>(1, count);
        for (int index = 0; index < count; ++index) {
            pointArray->SetValue(index + 1, points[index]);
            parameterArray->SetValue(index + 1, parameters[index]);
        }
        GeomAPI_Interpolate interpolate(pointArray,
                                        parameterArray,
                                        false,
                                        Precision::Confusion());
        interpolate.Perform();
        if (!interpolate.IsDone()) {
            return {};
        }
        return interpolate.Curve();
    }

    static void dropCoincidentSamples(std::vector<gp_Pnt> &points)
    {
        std::vector<gp_Pnt> kept;
        kept.reserve(points.size());
        for (const gp_Pnt &point : points) {
            if (kept.empty()
                || point.Distance(kept.back()) > Precision::Confusion()) {
                kept.push_back(point);
            }
        }
        points = std::move(kept);
    }

    // Re-samples a boundary polyline to targetCount points spaced evenly
    // along its interpolated curve, so differently sampled boundaries can
    // be paired rung by rung.
    static bool resampleBoundary(std::vector<gp_Pnt> &points,
                                 int targetCount)
    {
        dropCoincidentSamples(points);
        if (points.size() < 2 || targetCount < 2) {
            return false;
        }
        std::vector<double> parameters(points.size(), 0.0);
        for (std::size_t index = 1; index < points.size(); ++index) {
            parameters[index] = parameters[index - 1]
                                + points[index].Distance(points[index - 1]);
        }
        const occ::handle<Geom_BSplineCurve> curve =
            interpolateBoundary(points, parameters);
        if (curve.IsNull()) {
            return false;
        }
        const double length = parameters.back();
        std::vector<gp_Pnt> resampled(targetCount);
        for (int index = 0; index < targetCount; ++index) {
            resampled[index] =
                curve->Value(length * index / (targetCount - 1));
        }
        points = std::move(resampled);
        return true;
    }

    // An exact ruled loft between the strip's two boundary polylines: both
    // interpolate their samples on one shared chord-length
    // parameterization, so each legacy rung (sample j to sample j) and the
    // skin-following curvature of both boundary edges lie exactly on the
    // face.
    static TopoDS_Face makeStripFace(const CapturedStrip &strip)
    {
        std::vector<gp_Pnt> a = strip.curveA;
        std::vector<gp_Pnt> b = strip.curveB;
        if (a.size() != b.size()) {
            const int targetCount =
                static_cast<int>(std::max(a.size(), b.size()));
            if (!resampleBoundary(a, targetCount)
                || !resampleBoundary(b, targetCount)) {
                return {};
            }
        }
        // Drop rungs that collapse on either boundary (piecewise legacy
        // curves repeat their junction samples).
        std::vector<gp_Pnt> firstBoundary{a.front()};
        std::vector<gp_Pnt> secondBoundary{b.front()};
        for (std::size_t index = 1; index < a.size(); ++index) {
            if (a[index].Distance(firstBoundary.back())
                    > Precision::Confusion()
                && b[index].Distance(secondBoundary.back())
                       > Precision::Confusion()) {
                firstBoundary.push_back(a[index]);
                secondBoundary.push_back(b[index]);
            }
        }
        if (firstBoundary.size() < 2) {
            return {};
        }
        std::vector<double> parameters(firstBoundary.size(), 0.0);
        for (std::size_t index = 1; index < firstBoundary.size();
             ++index) {
            parameters[index] =
                parameters[index - 1]
                + 0.5
                      * (firstBoundary[index].Distance(
                             firstBoundary[index - 1])
                         + secondBoundary[index].Distance(
                             secondBoundary[index - 1]));
        }
        const occ::handle<Geom_BSplineCurve> firstCurve =
            interpolateBoundary(firstBoundary, parameters);
        const occ::handle<Geom_BSplineCurve> secondCurve =
            interpolateBoundary(secondBoundary, parameters);
        if (firstCurve.IsNull() || secondCurve.IsNull()
            || firstCurve->NbPoles() != secondCurve->NbPoles()
            || firstCurve->NbKnots() != secondCurve->NbKnots()) {
            return {};
        }
        NCollection_Array2<gp_Pnt> poles(1, firstCurve->NbPoles(), 1, 2);
        for (int index = 1; index <= firstCurve->NbPoles(); ++index) {
            poles.SetValue(index, 1, firstCurve->Pole(index));
            poles.SetValue(index, 2, secondCurve->Pole(index));
        }
        NCollection_Array1<double> uKnots(1, firstCurve->NbKnots());
        firstCurve->Knots(uKnots);
        NCollection_Array1<int> uMultiplicities(1, firstCurve->NbKnots());
        firstCurve->Multiplicities(uMultiplicities);
        NCollection_Array1<double> vKnots(1, 2);
        vKnots.SetValue(1, 0.0);
        vKnots.SetValue(2, 1.0);
        NCollection_Array1<int> vMultiplicities(1, 2);
        vMultiplicities.SetValue(1, 2);
        vMultiplicities.SetValue(2, 2);
        const occ::handle<Geom_BSplineSurface> surface =
            new Geom_BSplineSurface(poles,
                                    uKnots,
                                    vKnots,
                                    uMultiplicities,
                                    vMultiplicities,
                                    firstCurve->Degree(),
                                    1);
        BRepBuilderAPI_MakeFace makeFace(surface, Precision::Confusion());
        if (!makeFace.IsDone()) {
            return {};
        }
        return makeFace.Face();
    }

    // The H/V/VH diagonal strips and the trailing-edge mini-ribs become
    // ruled faces; strips sharing a label (e.g. the three sheets of one
    // VH-rib row) merge into a single named part.
    void addStripParts(AssemblyGroup &wing,
                       lep::NurbsWriteResult &result) const
    {
        if (capturedStrips_.empty()) {
            return;
        }
        struct StripPart
        {
            AssemblyGroup *group = nullptr;
            std::string label;
            PartColor color;
            std::vector<TopoDS_Shape> shapes;
            bool hasFaces = false;
        };
        std::vector<StripPart> stripParts;
        for (const CapturedStrip &strip : capturedStrips_) {
            AssemblyGroup *group =
                strip.minirib
                    ? &wing.group(ribsGroupName).group(miniribsGroupName)
                    : &wing.group(diagonalsGroupName);
            StripPart *part = nullptr;
            for (StripPart &candidate : stripParts) {
                if (candidate.group == group
                    && candidate.label == strip.label) {
                    part = &candidate;
                    break;
                }
            }
            if (part == nullptr) {
                stripParts.push_back(
                    {group,
                     strip.label,
                     strip.minirib ? ribColor : diagonalColor,
                     {},
                     false});
                part = &stripParts.back();
            }
            TopoDS_Shape shape;
            try {
                shape = makeStripFace(strip);
            } catch (const Standard_Failure &) {
                shape = {};
            }
            const bool face = !shape.IsNull();
            if (!face) {
                // Keep at least the boundary outline if OCCT could not
                // rule a face between the two curves.
                std::vector<TopoDS_Shape> outline;
                for (const std::vector<gp_Pnt> *boundary :
                     {&strip.curveA, &strip.curveB}) {
                    for (std::size_t index = 1; index < boundary->size();
                         ++index) {
                        BRepBuilderAPI_MakeEdge makeEdge(
                            makeLinearSpline((*boundary)[index - 1],
                                             (*boundary)[index]));
                        if (makeEdge.IsDone()) {
                            outline.push_back(makeEdge.Edge());
                        }
                    }
                }
                if (outline.empty()) {
                    continue;
                }
                warnings_.push_back(
                    strip.label
                    + " could not be ruled into a surface and stays a "
                      "wireframe outline.");
                shape = makeCompound(outline);
            }
            // A strip whose samples all sit on the symmetry plane, or one
            // that spans the centre cell mapping onto itself, must not get
            // a mirror copy.
            const auto centered = [](const gp_Pnt &point) {
                return std::abs(point.X())
                       <= symmetryPlaneToleranceMillimetres;
            };
            const bool onCenter =
                std::all_of(strip.curveA.begin(), strip.curveA.end(),
                            centered)
                && std::all_of(strip.curveB.begin(), strip.curveB.end(),
                               centered);
            const auto mirrorsOntoItself = [&]() {
                if (strip.curveA.size() != strip.curveB.size()) {
                    return false;
                }
                for (std::size_t index = 0; index < strip.curveA.size();
                     ++index) {
                    const gp_Pnt &a = strip.curveA[index];
                    if (gp_Pnt(-a.X(), a.Y(), a.Z())
                            .Distance(strip.curveB[index])
                        > symmetryPlaneToleranceMillimetres) {
                        return false;
                    }
                }
                return true;
            };
            part->shapes.push_back(shape);
            part->hasFaces = part->hasFaces || face;
            result.surfaceCount += face ? 1 : 0;
            if (!onCenter && !mirrorsOntoItself()) {
                const TopoDS_Shape mirroredShape = mirrored(shape);
                if (!mirroredShape.IsNull()) {
                    part->shapes.push_back(mirroredShape);
                    result.surfaceCount += face ? 1 : 0;
                }
            }
        }
        for (const StripPart &part : stripParts) {
            if (part.shapes.empty()) {
                continue;
            }
            part.group->parts.push_back({part.label,
                                         makeCompound(part.shapes),
                                         part.color,
                                         part.color,
                                         part.hasFaces});
        }
    }

    void writeAssembly(const AssemblyGroup &wing,
                       const std::filesystem::path &path,
                       lep::NurbsWriteResult &result) const
    {
        // A .xbf target serializes the XCAF document in OCCT's binary
        // format: identical geometry and assembly names, but no STEP
        // entity translation, so previews write and load far faster.
        const bool binary = path.extension() == ".xbf";
        if (binary) {
            BinXCAFDrivers::DefineFormat(
                XCAFApp_Application::GetApplication());
        }
        const occ::handle<TDocStd_Document> document =
            newDocument(binary ? "BinXCAF" : "MDTV-XCAF");
        if (document.IsNull()) {
            result.error = "OCCT could not create the XCAF document.";
            return;
        }
        const occ::handle<XCAFDoc_ShapeTool> shapeTool =
            XCAFDoc_DocumentTool::ShapeTool(document->Main());
        const occ::handle<XCAFDoc_ColorTool> colorTool =
            XCAFDoc_DocumentTool::ColorTool(document->Main());
        XCAFDoc_ShapeTool::SetAutoNaming(false);

        const TDF_Label wingLabel = shapeTool->NewShape();
        setLabelName(wingLabel, wing.name);
        addGroupContents(shapeTool, colorTool, wingLabel, wing, result);
        shapeTool->UpdateAssemblies();

        if (binary) {
            const auto encoded = path.u8string();
            const std::string encodedPath{
                reinterpret_cast<const char *>(encoded.data()),
                encoded.size()};
            if (XCAFApp_Application::GetApplication()->SaveAs(
                    document,
                    TCollection_ExtendedString(encodedPath.c_str(), true))
                != PCDM_SS_OK) {
                result.error =
                    "OCCT could not write the binary XCAF model.";
            }
            return;
        }

        // AP242 preserves exact B-spline geometry, carries the assembly
        // names and colours, and is the current STEP schema. All model
        // coordinates are already represented in mm.
        STEPCAFControl_Writer writer;
        // The writer initializes the shared STEP parameters in its
        // constructor, so select the schema only after constructing it and
        // recreate the model with those parameters.
        Interface_Static::SetIVal("write.step.schema", 5); // AP242DIS
        Interface_Static::SetCVal("write.step.unit", "MM");
        Interface_Static::SetIVal("write.surfacecurve.mode", 1);
        writer.ChangeWriter().Model(true);
        writer.ChangeWriter().SetTolerance(pointToleranceMillimetres);
        writer.SetColorMode(true);
        writer.SetNameMode(true);
        if (!writer.Transfer(document, STEPControl_AsIs)) {
            result.error = "OCCT could not transfer the NURBS model to STEP.";
            return;
        }

        const auto encoded = path.u8string();
        const std::string encodedPath{
            reinterpret_cast<const char *>(encoded.data()),
            encoded.size()};
        if (writer.Write(encodedPath.c_str()) != IFSelect_RetDone) {
            result.error = "OCCT could not write the STEP file.";
        }
    }

    static occ::handle<TDocStd_Document> newDocument(const char *format)
    {
        occ::handle<TDocStd_Document> document;
        XCAFApp_Application::GetApplication()->NewDocument(
            format,
            document);
        return document;
    }

    static void setLabelName(const TDF_Label &label, const std::string &name)
    {
        TDataStd_Name::Set(
            label,
            TCollection_ExtendedString(name.c_str(), true));
    }

    void addGroupContents(const occ::handle<XCAFDoc_ShapeTool> &shapeTool,
                          const occ::handle<XCAFDoc_ColorTool> &colorTool,
                          const TDF_Label &groupLabel,
                          const AssemblyGroup &group,
                          lep::NurbsWriteResult &result) const
    {
        for (const AssemblyGroup &child : group.groups) {
            if (child.empty()) {
                continue;
            }
            const TDF_Label childLabel = shapeTool->NewShape();
            setLabelName(childLabel, child.name);
            const TDF_Label instance =
                shapeTool->AddComponent(
                    groupLabel,
                    childLabel,
                    TopLoc_Location());
            setLabelName(instance, child.name);
            addGroupContents(shapeTool, colorTool, childLabel, child, result);
        }
        for (const AssemblyPart &part : group.parts) {
            const TDF_Label partLabel =
                shapeTool->AddShape(part.shape, false);
            setLabelName(partLabel, part.name);
            const TDF_Label instance =
                shapeTool->AddComponent(
                    groupLabel,
                    partLabel,
                    TopLoc_Location());
            setLabelName(instance, part.name);
            colorTool->SetColor(
                partLabel,
                part.hasFaces ? part.faceColor.color()
                              : part.curveColor.color(),
                XCAFDoc_ColorGen);
            if (part.hasFaces) {
                colorTool->SetColor(
                    partLabel,
                    part.faceColor.color(),
                    XCAFDoc_ColorSurf);
            }
            colorTool->SetColor(
                partLabel,
                part.curveColor.color(),
                XCAFDoc_ColorCurv);
            ++result.partCount;
        }
    }

public:
    lep::SimWingSceneExportResult buildSimWingScene(
        const lep::SimWingSceneExportSettings &settings) const
    {
        using namespace simwing::fsi;
        lep::SimWingSceneExportResult result;
        Scene scene;
        scene.metadata.designChecksum = settings.designChecksum;
        scene.metadata.exporterVersion = settings.exporterVersion;
        scene.metadata.sourceLengthToMeters = 0.01;

        const auto addError = [&result](std::string message) {
            result.errors.push_back(std::move(message));
        };
        const auto finalizeFailure = [&result]() {
            std::sort(result.errors.begin(), result.errors.end());
            result.errors.erase(
                std::unique(result.errors.begin(), result.errors.end()),
                result.errors.end());
            std::sort(result.warnings.begin(), result.warnings.end());
            result.scene = {};
            result.success = false;
        };
        if (simRegions_.empty()) {
            addError("No analytical skin was captured");
        }
        if (settings.designChecksum.empty()) {
            addError("A design checksum is required for scene-v2 export");
        }
        if (settings.exporterVersion.empty()) {
            addError("An exporter version is required for scene-v2 export");
        }
        if (!std::isfinite(settings.surfaceEndpointMatchToleranceMeters)
            || !(settings.surfaceEndpointMatchToleranceMeters > 0.0)) {
            addError(
                "A finite positive surface endpoint match tolerance is required");
        }
        if (!std::isfinite(settings.pilot.endpointMatchToleranceMeters)
            || !(settings.pilot.endpointMatchToleranceMeters > 0.0)) {
            addError(
                "A finite positive pilot endpoint match tolerance is required");
        }
        if (!capturedLines_.empty()
            && (!std::isfinite(settings.suspensionJunctionMassKg)
                || !(settings.suspensionJunctionMassKg > 0.0))) {
            addError(
                "A finite positive suspension-junction mass is required when lines are captured");
        }
        for (const SimRegionCapture &capture : simRegions_) {
            if (capture.singleSkin) {
                addError(
                    "Single-skin captured geometry is not yet a closed scene-v2 fluid region");
                break;
            }
        }

        constexpr std::uint64_t ordinalMask =
            (std::uint64_t{1} << 56U) - 1U;
        const auto stableId = [&addError, ordinalMask](std::uint8_t category,
                                                       std::size_t ordinal) {
            const std::uint64_t oneBased =
                static_cast<std::uint64_t>(ordinal) + 1U;
            if (oneBased > ordinalMask) {
                addError("Scene entity count exceeds the stable-ID range");
                return invalidStableId;
            }
            return (static_cast<std::uint64_t>(category) << 56U)
                | oneBased;
        };

        std::map<SurfaceRole, const lep::SimWingFabricExportSettings *>
            fabricSettings;
        for (const auto &material : settings.fabricMaterials) {
            if (!fabricSettings.emplace(material.role, &material).second) {
                addError("Duplicate fabric assignment for surface role "
                         + std::to_string(static_cast<int>(material.role)));
            }
        }
        std::set<SurfaceRole> requiredRoles{SurfaceRole::Skin};
        if (!capturedRibs_.empty()) {
            requiredRoles.insert(SurfaceRole::Rib);
        }
        for (const CapturedStrip &strip : capturedStrips_) {
            requiredRoles.insert(strip.minirib ? SurfaceRole::MiniRib
                                               : SurfaceRole::Diagonal);
        }
        std::map<SurfaceRole, StableId> fabricMaterialIds;
        std::size_t fabricOrdinal = 0;
        for (const SurfaceRole role : requiredRoles) {
            const auto found = fabricSettings.find(role);
            if (found == fabricSettings.end()) {
                addError("Missing explicit fabric assignment for surface role "
                         + std::to_string(static_cast<int>(role)));
                continue;
            }
            const auto &source = *found->second;
            const StableId id = stableId(3, fabricOrdinal++);
            fabricMaterialIds.emplace(role, id);
            scene.fabricMaterials.push_back(
                {id,
                 source.name,
                 source.warpStiffnessNewtonsPerMeter,
                 source.weftStiffnessNewtonsPerMeter,
                 source.shearStiffnessNewtonsPerMeter,
                 source.bendingStiffnessNewtonMeters,
                 source.arealDensityKgPerSquareMeter,
                 source.dampingSeconds,
                 source.porosityFraction,
                 source.permeabilitySquareMeters});
        }

        std::map<std::string, const lep::SimWingLineExportSettings *>
            lineSettings;
        for (const auto &material : settings.lineMaterials) {
            if (!lineSettings.emplace(material.captureTypeName, &material)
                     .second) {
                addError("Duplicate line material assignment for captured type '"
                         + material.captureTypeName + "'");
            }
        }
        std::set<std::string> requiredLineTypes;
        for (const CapturedLine &line : capturedLines_) {
            requiredLineTypes.insert(line.typeName);
        }
        std::map<std::string, StableId> lineMaterialIds;
        std::size_t lineMaterialOrdinal = 0;
        for (const std::string &type : requiredLineTypes) {
            const auto found = lineSettings.find(type);
            if (found == lineSettings.end()) {
                addError("Missing explicit line material assignment for captured type '"
                         + type + "'");
                continue;
            }
            const auto &source = *found->second;
            const StableId id = stableId(6, lineMaterialOrdinal++);
            lineMaterialIds.emplace(type, id);
            scene.lineMaterials.push_back(
                {id,
                 source.name,
                 source.diameterMeters,
                 source.linearDensityKgPerMeter,
                 source.axialStiffnessNewtons,
                 source.dragCoefficient});
        }
        if (!result.errors.empty()) {
            finalizeFailure();
            return result;
        }

        const StableId outsideRegionId = stableId(1, 0);
        scene.regions.push_back(
            {outsideRegionId, RegionKind::Outside, "outside"});
        struct CellKey
        {
            int panel = 0;
            bool mirror = false;
            bool operator<(const CellKey &other) const
            {
                return std::tie(panel, mirror)
                    < std::tie(other.panel, other.mirror);
            }
        };
        struct CellAccumulator
        {
            gp_XYZ sum{0.0, 0.0, 0.0};
            std::size_t count = 0;
            StableId regionId = invalidStableId;
        };
        std::map<CellKey, CellAccumulator> cells;
        for (const SimRegionCapture &capture : simRegions_) {
            for (const bool mirror : {false, true}) {
                if (mirror && capture.selfMirrored()) {
                    continue;
                }
                CellAccumulator &cell = cells[{capture.panelIndex, mirror}];
                for (const auto &row : capture.rows) {
                    for (gp_Pnt point : row) {
                        if (mirror) {
                            point.SetX(-point.X());
                        }
                        cell.sum += point.XYZ();
                        ++cell.count;
                    }
                }
            }
        }
        std::size_t cellOrdinal = 1;
        for (auto &[key, cell] : cells) {
            cell.regionId = stableId(1, cellOrdinal++);
            scene.regions.push_back(
                {cell.regionId,
                 RegionKind::Cell,
                 std::string(key.mirror ? "left-cell-" : "right-cell-")
                     + std::to_string(key.panel)});
        }
        const auto cellCenter = [&cells](const CellKey &key) {
            const CellAccumulator &cell = cells.at(key);
            return gp_Pnt(cell.sum / static_cast<double>(cell.count));
        };

        std::unordered_map<QuantizedPoint, StableId, QuantizedPointHash>
            vertexIds;
        const auto place = [](gp_Pnt point, bool mirror) {
            const double snappedX =
                std::abs(point.X()) < symmetryPlaneToleranceMillimetres
                    ? 0.0
                    : point.X();
            point.SetX(mirror ? -snappedX : snappedX);
            return point;
        };
        const auto vertexId = [&](const gp_Pnt &point) {
            const QuantizedPoint key = quantize(point);
            const auto found = vertexIds.find(key);
            if (found != vertexIds.end()) {
                return found->second;
            }
            const StableId id = stableId(2, scene.vertices.size());
            vertexIds.emplace(key, id);
            scene.vertices.push_back(
                {id,
                 {point.X() * 0.001,
                  point.Y() * 0.001,
                  point.Z() * 0.001}});
            return id;
        };

        const auto pointForVertexId = [&scene](const StableId id) {
            const auto vertex = std::lower_bound(
                scene.vertices.begin(), scene.vertices.end(), id,
                [](const Vertex &candidate, const StableId value) {
                    return candidate.id < value;
                });
            if (vertex == scene.vertices.end() || vertex->id != id) {
                throw std::logic_error(
                    "Scene triangle references a missing vertex");
            }
            return gp_Pnt(vertex->positionMeters.x * 1000.0,
                          vertex->positionMeters.y * 1000.0,
                          vertex->positionMeters.z * 1000.0);
        };

        const auto chartFor = [](const std::array<gp_Pnt, 3> &points,
                                 gp_Vec preferredWarp) {
            std::array<Vec2, 3> chart{};
            gp_Vec first(points[0], points[1]);
            gp_Vec second(points[0], points[2]);
            gp_Vec normal = first.Crossed(second);
            if (normal.Magnitude() <= Precision::Confusion()) {
                return chart;
            }
            normal.Normalize();
            preferredWarp -= normal.Multiplied(preferredWarp.Dot(normal));
            if (preferredWarp.Magnitude() <= Precision::Confusion()) {
                preferredWarp = first;
            }
            preferredWarp.Normalize();
            gp_Vec weft = normal.Crossed(preferredWarp);
            weft.Normalize();
            for (std::size_t index = 1; index < 3; ++index) {
                const gp_Vec offset(points[0], points[index]);
                chart[index] = {offset.Dot(preferredWarp) * 0.001,
                                offset.Dot(weft) * 0.001};
            }
            const double determinant =
                chart[1].x * chart[2].y
                - chart[2].x * chart[1].y;
            if (determinant < 0.0) {
                for (Vec2 &coordinate : chart) {
                    coordinate.y = -coordinate.y;
                }
            }
            return chart;
        };

        std::size_t triangleOrdinal = 0;
        std::size_t sheetOrdinal = 0;
        std::map<std::array<StableId, 3>,
                 std::pair<StableId, StableId>> directedSheetEdges;
        const auto addTriangle = [&](std::array<gp_Pnt, 3> points,
                                      gp_Vec preferredWarp,
                                      StableId negativeRegion,
                                      StableId positiveRegion,
                                      StableId sheetId,
                                      SurfaceRole role) {
            gp_Vec normal(points[0], points[1]);
            normal.Cross(gp_Vec(points[0], points[2]));
            if (normal.Magnitude() <= Precision::Confusion()) {
                std::ostringstream message;
                message.precision(17);
                message
                    << "Captured surface role "
                    << static_cast<int>(role)
                    << ", sheet " << sheetId
                    << " produced a degenerate triangle at ("
                    << points[0].X() << ", " << points[0].Y() << ", "
                    << points[0].Z() << "), ("
                    << points[1].X() << ", " << points[1].Y() << ", "
                    << points[1].Z() << "), ("
                    << points[2].X() << ", " << points[2].Y() << ", "
                    << points[2].Z() << ')';
                addError(message.str());
                return;
            }
            std::array<StableId, 3> ids{
                vertexId(points[0]),
                vertexId(points[1]),
                vertexId(points[2])};
            for (std::size_t corner = 0; corner < 3; ++corner) {
                points[corner] = pointForVertexId(ids[corner]);
            }
            normal = gp_Vec(points[0], points[1]);
            normal.Cross(gp_Vec(points[0], points[2]));
            if (normal.Magnitude() <= Precision::Confusion()) {
                addError("Captured surface collapsed after canonical vertex mapping");
                return;
            }
            std::optional<bool> reverseWinding;
            for (std::size_t corner = 0; corner < 3; ++corner) {
                const StableId from = ids[corner];
                const StableId to = ids[(corner + 1) % 3];
                const std::array<StableId, 3> key{
                    sheetId, std::min(from, to), std::max(from, to)};
                const auto found = directedSheetEdges.find(key);
                if (found == directedSheetEdges.end()) {
                    continue;
                }
                const bool needsReverse = found->second
                    == std::pair{from, to};
                if (reverseWinding && *reverseWinding != needsReverse) {
                    addError("Captured sheet " + std::to_string(sheetId)
                             + " has contradictory triangle winding");
                    return;
                }
                reverseWinding = needsReverse;
            }
            if (reverseWinding.value_or(false)) {
                std::swap(points[1], points[2]);
                std::swap(ids[1], ids[2]);
            }
            std::array<Vec2, 3> chart = chartFor(points, preferredWarp);
            scene.triangles.push_back(
                {stableId(4, triangleOrdinal++),
                 ids,
                 chart,
                  negativeRegion,
                  positiveRegion,
                  fabricMaterialIds.at(role),
                  sheetId,
                  role});
            for (std::size_t corner = 0; corner < 3; ++corner) {
                const StableId from = ids[corner];
                const StableId to = ids[(corner + 1) % 3];
                const std::array<StableId, 3> key{
                    sheetId, std::min(from, to), std::max(from, to)};
                directedSheetEdges.try_emplace(key,
                                               std::pair{from, to});
            }
        };

        struct PendingIntake
        {
            StableId id = invalidStableId;
            StableId cellRegionId = invalidStableId;
            int firstRibIndex = 0;
            int lastRibIndex = 0;
            std::vector<StableId> frontLipVertexIds;
            std::vector<StableId> backLipVertexIds;
        };
        struct RibOuterBoundary
        {
            int ribIndex = 0;
            std::vector<StableId> vertexIds;
            std::vector<gp_Pnt> points;
        };
        std::size_t openingOrdinal = 0;
        std::vector<PendingIntake> pendingIntakes;
        std::set<std::pair<StableId, int>> pendingRibBoundaryVertices;
        std::vector<RibOuterBoundary> ribOuterBoundaries;
        double maximumRibWeldDeviationMillimetres = 0.0;
        for (const SimRegionCapture &capture : simRegions_) {
            for (const bool mirror : {false, true}) {
                if (mirror && capture.selfMirrored()) {
                    continue;
                }
                const CellKey cellKey{capture.panelIndex, mirror};
                const StableId cellRegion = cells.at(cellKey).regionId;
                std::vector<std::array<gp_Pnt, simSpanColumnCount>> grid;
                grid.reserve(capture.rows.size());
                for (const auto &row : capture.rows) {
                    std::array<gp_Pnt, simSpanColumnCount> placed{};
                    for (int column = 0; column < simSpanColumnCount;
                         ++column) {
                        placed[column] = place(row[column], mirror);
                    }
                    grid.push_back(placed);
                }
                bool ribWeldValid = true;
                for (std::size_t row = 0; row < grid.size(); ++row) {
                    const auto canonicalRibPoint =
                        [&](int ribIndex) -> std::optional<gp_Pnt> {
                        const auto rib = capturedRibs_.find(ribIndex);
                        const int station = capture.stations[row];
                        if (rib == capturedRibs_.end()
                            || station < 1
                            || static_cast<std::size_t>(station)
                                > rib->second.spatialPoints.size()) {
                            return std::nullopt;
                        }
                        return place(
                            rib->second.spatialPoints[
                                static_cast<std::size_t>(station - 1)],
                            mirror);
                    };
                    const auto first = canonicalRibPoint(
                        capture.panelIndex);
                    const auto last = canonicalRibPoint(
                        capture.panelIndex - 1);
                    if (!first || !last) {
                        addError(
                            "Captured scene skin boundary has no corresponding rib station");
                        ribWeldValid = false;
                        break;
                    }
                    for (const auto &[column, point] :
                         {std::pair{0, *first},
                          std::pair{simSpanColumnCount - 1, *last}}) {
                        const double deviation =
                            grid[row][column].Distance(point);
                        maximumRibWeldDeviationMillimetres = std::max(
                            maximumRibWeldDeviationMillimetres, deviation);
                        if (deviation
                            > maximumSceneRibWeldDeviationMillimetres) {
                            addError(
                                "Captured scene skin and rib stations disagree beyond the weld tolerance");
                            ribWeldValid = false;
                            break;
                        }
                        grid[row][column] = point;
                    }
                    if (!ribWeldValid) {
                        break;
                    }
                }
                if (!ribWeldValid) {
                    continue;
                }
                if (capture.authoredSurface) {
                    for (std::size_t row = 0; row < grid.size(); ++row) {
                        pendingRibBoundaryVertices.emplace(
                            vertexId(grid[row].front()),
                            capture.panelIndex);
                        pendingRibBoundaryVertices.emplace(
                            vertexId(grid[row].back()),
                            capture.panelIndex - 1);
                    }
                }
                if (!capture.authoredSurface) {
                    PendingIntake intake;
                    intake.id = stableId(5, openingOrdinal++);
                    intake.cellRegionId = cellRegion;
                    intake.firstRibIndex = capture.panelIndex;
                    intake.lastRibIndex = capture.panelIndex - 1;
                    intake.frontLipVertexIds.reserve(simSpanColumnCount);
                    intake.backLipVertexIds.reserve(simSpanColumnCount);
                    for (int column = 0; column < simSpanColumnCount;
                         ++column) {
                        intake.frontLipVertexIds.push_back(
                            vertexId(grid.front()[column]));
                        intake.backLipVertexIds.push_back(
                            vertexId(grid.back()[column]));
                    }
                    pendingIntakes.push_back(std::move(intake));
                    continue;
                }

                const StableId sheetId = stableId(11, sheetOrdinal++);
                const gp_Pnt center = cellCenter(cellKey);
                const auto oriented = [&center](std::array<gp_Pnt, 3> points) {
                    gp_Vec normal(points[0], points[1]);
                    normal.Cross(gp_Vec(points[0], points[2]));
                    const gp_Pnt centroid(
                        (points[0].XYZ() + points[1].XYZ()
                         + points[2].XYZ())
                        / 3.0);
                    if (normal.Dot(gp_Vec(centroid, center)) < 0.0) {
                        std::swap(points[1], points[2]);
                    }
                    return points;
                };
                for (std::size_t row = 0; row + 1 < grid.size(); ++row) {
                    for (int column = 0;
                         column + 1 < simSpanColumnCount;
                         ++column) {
                        const gp_Vec chordDirection(
                            grid[row][column], grid[row + 1][column]);
                        const gp_Vec spanDirection(
                            grid[row][column], grid[row][column + 1]);
                        const auto material = fabricSettings.at(
                            SurfaceRole::Skin);
                        const gp_Vec warpDirection =
                            material->warpAxis
                                    == lep::SimWingMaterialAxis::Primary
                                ? chordDirection
                                : spanDirection;
                        addTriangle(
                            oriented({grid[row][column],
                                      grid[row][column + 1],
                                      grid[row + 1][column + 1]}),
                            warpDirection,
                            outsideRegionId,
                            cellRegion,
                            sheetId,
                            SurfaceRole::Skin);
                        addTriangle(
                            oriented({grid[row][column],
                                      grid[row + 1][column + 1],
                                      grid[row + 1][column]}),
                            warpDirection,
                            outsideRegionId,
                            cellRegion,
                            sheetId,
                            SurfaceRole::Skin);
                    }
                }
            }
        }
        if (maximumRibWeldDeviationMillimetres
            > pointToleranceMillimetres) {
            std::ostringstream warning;
            warning << "Scene skin boundaries were canonicalized to captured rib stations by up to "
                    << maximumRibWeldDeviationMillimetres << " mm";
            result.warnings.push_back(warning.str());
        }

        // Build the same exact planar, holed rib faces used by the CAD
        // exporter, then triangulate them for scene-v2. A rib separates its
        // two adjacent cells; a boundary rib separates its sole cell from
        // the outside. The centre rib is emitted only once and separates the
        // mirrored centre cells.
        std::unordered_set<int> capturedPanelIndices;
        std::map<int, std::vector<RibBoundarySegment>> ribSegments;
        for (const PanelSurface &panel : panels_) {
            capturedPanelIndices.insert(panel.panelIndex);
        }
        for (const PanelSurface &panel : panels_) {
            ribSegments[panel.panelIndex].push_back(
                {panel.firstPoint, panel.lastPoint});
            if (!capturedPanelIndices.contains(panel.panelIndex - 1)) {
                ribSegments[panel.panelIndex - 1].push_back(
                    {panel.firstPoint, panel.lastPoint});
            }
        }
        std::map<int, std::set<int>> sceneRibContourStations;
        for (const SimRegionCapture &capture : simRegions_) {
            auto &first = sceneRibContourStations[capture.panelIndex];
            auto &last = sceneRibContourStations[capture.panelIndex - 1];
            first.insert(capture.stations.begin(), capture.stations.end());
            last.insert(capture.stations.begin(), capture.stations.end());
        }

        const auto ribCandidateCells = [&cells](int ribIndex,
                                                 bool mirror,
                                                 bool onCenter) {
            std::vector<CellKey> candidates;
            const auto add = [&cells, &candidates](const CellKey &key) {
                if (cells.contains(key)) {
                    candidates.push_back(key);
                }
            };
            if (onCenter) {
                add({ribIndex + 1, false});
                add({ribIndex + 1, true});
            } else {
                add({ribIndex, mirror});
                add({ribIndex + 1, mirror});
            }
            return candidates;
        };
        const auto sideRegionsFor = [&](const std::array<gp_Pnt, 3> &points,
                                        const std::vector<CellKey> &candidates,
                                        const std::string &label) {
            std::pair<StableId, StableId> sides{invalidStableId,
                                                invalidStableId};
            gp_Vec normal(points[0], points[1]);
            normal.Cross(gp_Vec(points[0], points[2]));
            if (normal.Magnitude() <= Precision::Confusion()) {
                addError(label + " produced a degenerate mesh triangle");
                return sides;
            }
            normal.Normalize();
            const gp_Pnt centroid(
                (points[0].XYZ() + points[1].XYZ() + points[2].XYZ())
                / 3.0);
            std::vector<std::pair<double, StableId>> distances;
            distances.reserve(candidates.size());
            for (const CellKey &candidate : candidates) {
                distances.emplace_back(
                    normal.Dot(gp_Vec(centroid, cellCenter(candidate))),
                    cells.at(candidate).regionId);
            }
            if (distances.empty()) {
                addError(label + " has no adjacent captured cell");
                return sides;
            }
            std::sort(distances.begin(), distances.end());
            constexpr double sideToleranceMillimetres = 0.01;
            if (distances.size() == 1) {
                if (std::abs(distances.front().first)
                    <= sideToleranceMillimetres) {
                    addError(label
                             + " has an indeterminate adjacent-cell side");
                    return sides;
                }
                if (distances.front().first > 0.0) {
                    sides = {outsideRegionId, distances.front().second};
                } else {
                    sides = {distances.front().second, outsideRegionId};
                }
                return sides;
            }
            if (distances.front().first >= -sideToleranceMillimetres
                || distances.back().first <= sideToleranceMillimetres) {
                addError(label
                         + " does not geometrically separate its adjacent cells");
                return sides;
            }
            sides = {distances.front().second, distances.back().second};
            return sides;
        };
        const auto meshedWireBoundary =
            [&](const TopoDS_Wire &wire,
                const occ::handle<Poly_Triangulation> &triangulation,
                const TopLoc_Location &location,
                const std::string &description)
                -> std::optional<RibOuterBoundary> {
            RibOuterBoundary boundary;
            for (BRepTools_WireExplorer edgeExplorer(wire);
                 edgeExplorer.More(); edgeExplorer.Next()) {
                const TopoDS_Edge edge = edgeExplorer.Current();
                const occ::handle<Poly_PolygonOnTriangulation> polygon =
                    BRep_Tool::PolygonOnTriangulation(
                        edge, triangulation, location);
                if (polygon.IsNull() || polygon->NbNodes() < 2) {
                    addError("Could not recover a meshed " + description);
                    return std::nullopt;
                }
                const bool reversed = edge.Orientation() == TopAbs_REVERSED;
                for (int point = 1; point <= polygon->NbNodes(); ++point) {
                    const int polygonIndex = reversed
                        ? polygon->NbNodes() - point + 1 : point;
                    const int node = polygon->Node(polygonIndex);
                    const gp_Pnt position =
                        triangulation->Node(node).Transformed(
                            location.Transformation());
                    const StableId id = vertexId(position);
                    if (!boundary.vertexIds.empty()
                        && boundary.vertexIds.back() == id) {
                        continue;
                    }
                    boundary.vertexIds.push_back(id);
                    boundary.points.push_back(position);
                }
            }
            if (boundary.vertexIds.size() >= 2
                && boundary.vertexIds.front()
                    == boundary.vertexIds.back()) {
                boundary.vertexIds.pop_back();
                boundary.points.pop_back();
            }
            const std::set<StableId> unique(
                boundary.vertexIds.begin(), boundary.vertexIds.end());
            if (boundary.vertexIds.size() < 3
                || unique.size() != boundary.vertexIds.size()) {
                addError("A meshed " + description
                         + " does not form a unique closed loop");
                return std::nullopt;
            }
            return boundary;
        };

        for (const auto &[ribIndex, segments] : ribSegments) {
            const auto captured = capturedRibs_.find(ribIndex);
            if (captured == capturedRibs_.end()) {
                addError("Missing captured rib station for rib "
                         + std::to_string(ribIndex));
                continue;
            }
            const bool onCenter = std::all_of(
                captured->second.spatialPoints.begin(),
                captured->second.spatialPoints.end(),
                [](const gp_Pnt &point) {
                    return std::abs(point.X())
                        <= symmetryPlaneToleranceMillimetres;
                });
            int ignoredEdgeCount = 0;
            TopoDS_Shape curveFallback;
            std::vector<std::string> ribWarnings;
            std::vector<std::string> ribErrors;
            TopoDS_Face rightFace;
            try {
                rightFace = makeRibFace(ribIndex,
                                        segments,
                                        captured->second,
                                        ignoredEdgeCount,
                                        curveFallback,
                                        ribWarnings,
                                        ribErrors,
                                        &sceneRibContourStations.at(
                                            ribIndex));
            } catch (const Standard_Failure &failure) {
                ribErrors.push_back(
                    "OCCT failed building rib " + std::to_string(ribIndex)
                    + ": "
                    + (failure.GetMessageString() != nullptr
                           ? failure.GetMessageString()
                           : "unknown OCCT error"));
            } catch (const std::exception &exception) {
                ribErrors.push_back(
                    "Failed building rib " + std::to_string(ribIndex)
                    + ": " + exception.what());
            }
            result.warnings.insert(result.warnings.end(),
                                   ribWarnings.begin(),
                                   ribWarnings.end());
            result.errors.insert(result.errors.end(),
                                 ribErrors.begin(),
                                 ribErrors.end());
            if (!ribErrors.empty()) {
                continue;
            }
            if (rightFace.IsNull()) {
                if (!captured->second.holes.empty()) {
                    addError("Holed rib " + std::to_string(ribIndex)
                             + " collapsed to an outline");
                } else {
                    result.warnings.push_back(
                        "Skipped collapsed zero-area rib "
                        + std::to_string(ribIndex));
                }
                continue;
            }

            RibFrame materialFrame;
            if (!fitRibFrame(captured->second.planarPoints,
                             captured->second.spatialPoints,
                             materialFrame)) {
                addError("Could not recover the material axes of rib "
                         + std::to_string(ribIndex));
                continue;
            }

            std::size_t innerWireCount = 0;
            const TopoDS_Wire outer = BRepTools::OuterWire(rightFace);
            for (TopExp_Explorer explorer(rightFace, TopAbs_WIRE);
                 explorer.More(); explorer.Next()) {
                if (!TopoDS::Wire(explorer.Current()).IsSame(outer)) {
                    ++innerWireCount;
                }
            }
            if (innerWireCount != captured->second.holes.size()) {
                addError("Rib " + std::to_string(ribIndex) + " captured "
                         + std::to_string(captured->second.holes.size())
                         + " crossports but its exact face contains "
                         + std::to_string(innerWireCount));
                continue;
            }

            for (const bool mirror : {false, true}) {
                if (mirror && onCenter) {
                    continue;
                }
                TopoDS_Face face = rightFace;
                if (mirror) {
                    const TopoDS_Shape reflected = mirrored(rightFace);
                    if (reflected.IsNull()
                        || reflected.ShapeType() != TopAbs_FACE) {
                        addError("Could not mirror rib "
                                 + std::to_string(ribIndex));
                        continue;
                    }
                    face = TopoDS::Face(reflected);
                }
                const std::string label = std::string(mirror ? "left rib "
                                                              : "right rib ")
                    + std::to_string(ribIndex);
                const StableId sheetId = stableId(11, sheetOrdinal++);
                const std::vector<CellKey> candidates =
                    ribCandidateCells(ribIndex, mirror, onCenter);
                BRepTools::Clean(face);
                BRepMesh_IncrementalMesh mesher(face, 1.0, false, 0.20, true);
                if (!mesher.IsDone()) {
                    addError("Could not triangulate " + label);
                    continue;
                }
                TopLoc_Location location;
                const occ::handle<Poly_Triangulation> triangulation =
                    BRep_Tool::Triangulation(face, location);
                if (triangulation.IsNull()
                    || triangulation->NbTriangles() == 0) {
                    addError("OCCT returned no triangles for " + label);
                    continue;
                }

                std::pair<StableId, StableId> faceSides{
                    invalidStableId, invalidStableId};
                for (int index = 1;
                     index <= triangulation->NbTriangles(); ++index) {
                    Standard_Integer first = 0;
                    Standard_Integer second = 0;
                    Standard_Integer third = 0;
                    triangulation->Triangle(index).Get(
                        first, second, third);
                    std::array<gp_Pnt, 3> points{
                        triangulation->Node(first).Transformed(
                            location.Transformation()),
                        triangulation->Node(second).Transformed(
                            location.Transformation()),
                        triangulation->Node(third).Transformed(
                            location.Transformation())};
                    const auto sides = sideRegionsFor(points,
                                                     candidates,
                                                     label);
                    if (sides.first == invalidStableId) {
                        continue;
                    }
                    if (faceSides.first == invalidStableId) {
                        faceSides = sides;
                    } else if (faceSides != sides) {
                        addError(label
                                 + " triangulation has inconsistent winding");
                        continue;
                    }
                    gp_Vec preferredWarp =
                        fabricSettings.at(SurfaceRole::Rib)->warpAxis
                                == lep::SimWingMaterialAxis::Primary
                            ? materialFrame.axisX
                            : materialFrame.axisY;
                    if (mirror) {
                        preferredWarp.SetX(-preferredWarp.X());
                    }
                    addTriangle(points,
                                preferredWarp,
                                sides.first,
                                sides.second,
                                sheetId,
                                SurfaceRole::Rib);
                }
                if (faceSides.first == invalidStableId) {
                    continue;
                }

                const TopoDS_Wire faceOuter = BRepTools::OuterWire(face);
                const auto outerBoundary = meshedWireBoundary(
                    faceOuter, triangulation, location,
                    "outer boundary in " + label);
                if (!outerBoundary) {
                    continue;
                }
                RibOuterBoundary taggedOuterBoundary = *outerBoundary;
                taggedOuterBoundary.ribIndex = ribIndex;
                ribOuterBoundaries.push_back(
                    std::move(taggedOuterBoundary));
                for (TopExp_Explorer explorer(face, TopAbs_WIRE);
                     explorer.More(); explorer.Next()) {
                    const TopoDS_Wire wire =
                        TopoDS::Wire(explorer.Current());
                    if (wire.IsSame(faceOuter)) {
                        continue;
                    }
                    auto crossport = meshedWireBoundary(
                        wire, triangulation, location,
                        "crossport boundary in " + label);
                    if (!crossport) {
                        continue;
                    }
                    scene.openings.push_back(
                        {stableId(5, openingOrdinal++),
                         std::move(crossport->vertexIds),
                         faceSides.first,
                         faceSides.second,
                         OpeningRole::Crossport});
                }
            }
        }

        const auto ribBoundaryPath =
            [&ribOuterBoundaries](StableId start,
                                  StableId end,
                                  int ribIndex)
                -> std::optional<RibOuterBoundary> {
            std::optional<std::pair<double, RibOuterBoundary>> best;
            for (const RibOuterBoundary &loop : ribOuterBoundaries) {
                if (loop.ribIndex != ribIndex) {
                    continue;
                }
                const auto startIterator = std::find(
                    loop.vertexIds.begin(), loop.vertexIds.end(), start);
                const auto endIterator = std::find(
                    loop.vertexIds.begin(), loop.vertexIds.end(), end);
                if (startIterator == loop.vertexIds.end()
                    || endIterator == loop.vertexIds.end()) {
                    continue;
                }
                const std::size_t startIndex = static_cast<std::size_t>(
                    startIterator - loop.vertexIds.begin());
                const std::size_t endIndex = static_cast<std::size_t>(
                    endIterator - loop.vertexIds.begin());
                for (const int direction : {1, -1}) {
                    RibOuterBoundary candidate;
                    std::size_t index = startIndex;
                    double length = 0.0;
                    candidate.vertexIds.push_back(loop.vertexIds[index]);
                    candidate.points.push_back(loop.points[index]);
                    while (index != endIndex) {
                        const std::size_t next = direction > 0
                            ? (index + 1) % loop.vertexIds.size()
                            : (index + loop.vertexIds.size() - 1)
                                  % loop.vertexIds.size();
                        length += loop.points[index].Distance(
                            loop.points[next]);
                        index = next;
                        candidate.vertexIds.push_back(
                            loop.vertexIds[index]);
                        candidate.points.push_back(loop.points[index]);
                    }
                    const bool preferred = !best
                        || length < best->first
                        || (length == best->first
                            && candidate.vertexIds
                                < best->second.vertexIds);
                    if (preferred) {
                        best = std::pair{
                            length, std::move(candidate)};
                    }
                }
            }
            if (!best) {
                return std::nullopt;
            }
            return std::move(best->second);
        };
        std::map<StableId, StableId> ribBoundaryAliases;
        std::set<StableId> retiredVertexIds;
        double maximumRibMeshSnapDeviationMillimetres = 0.0;
        const auto trianglePoints = [&](const Triangle &triangle) {
            return std::array<gp_Pnt, 3>{
                pointForVertexId(triangle.vertexIds[0]),
                pointForVertexId(triangle.vertexIds[1]),
                pointForVertexId(triangle.vertexIds[2]),
            };
        };
        const auto trianglePreferredWarp = [&trianglePoints](
                                               const Triangle &triangle) {
            const auto points = trianglePoints(triangle);
            const Vec2 firstChart{
                triangle.materialCoordinates[1].x
                    - triangle.materialCoordinates[0].x,
                triangle.materialCoordinates[1].y
                    - triangle.materialCoordinates[0].y};
            const Vec2 secondChart{
                triangle.materialCoordinates[2].x
                    - triangle.materialCoordinates[0].x,
                triangle.materialCoordinates[2].y
                    - triangle.materialCoordinates[0].y};
            const double determinant = firstChart.x * secondChart.y
                - secondChart.x * firstChart.y;
            const gp_Vec first(points[0], points[1]);
            const gp_Vec second(points[0], points[2]);
            if (!std::isfinite(determinant) || determinant == 0.0) {
                return first;
            }
            return first.Multiplied(secondChart.y / determinant)
                - second.Multiplied(firstChart.y / determinant);
        };
        const auto canonicalRibBoundaryVertex =
            [&](StableId id, int ribIndex) -> std::optional<StableId> {
            const auto alias = ribBoundaryAliases.find(id);
            if (alias != ribBoundaryAliases.end()) {
                return alias->second;
            }
            for (const RibOuterBoundary &boundary : ribOuterBoundaries) {
                if (boundary.ribIndex == ribIndex
                    && std::find(boundary.vertexIds.begin(),
                              boundary.vertexIds.end(), id)
                    != boundary.vertexIds.end()) {
                    ribBoundaryAliases.emplace(id, id);
                    return id;
                }
            }
            const auto vertex = std::find_if(
                scene.vertices.begin(), scene.vertices.end(),
                [id](const Vertex &candidate) {
                    return candidate.id == id;
                });
            if (vertex == scene.vertices.end()) {
                return std::nullopt;
            }
            const gp_Pnt point(vertex->positionMeters.x * 1000.0,
                               vertex->positionMeters.y * 1000.0,
                               vertex->positionMeters.z * 1000.0);

            StableId replacement = invalidStableId;
            double bestDistance = std::numeric_limits<double>::max();
            for (const RibOuterBoundary &boundary : ribOuterBoundaries) {
                if (boundary.ribIndex != ribIndex) {
                    continue;
                }
                for (std::size_t index = 0;
                     index < boundary.vertexIds.size(); ++index) {
                    const double distance =
                        point.Distance(boundary.points[index]);
                    if (distance < bestDistance) {
                        bestDistance = distance;
                        replacement = boundary.vertexIds[index];
                    }
                }
            }
            if (replacement == invalidStableId
                || bestDistance
                    > maximumSceneRibWeldDeviationMillimetres) {
                return std::nullopt;
            }
            maximumRibMeshSnapDeviationMillimetres = std::max(
                maximumRibMeshSnapDeviationMillimetres, bestDistance);
            for (const Triangle &triangle : scene.triangles) {
                if (std::find(triangle.vertexIds.begin(),
                              triangle.vertexIds.end(), id)
                        != triangle.vertexIds.end()
                    && std::find(triangle.vertexIds.begin(),
                                 triangle.vertexIds.end(), replacement)
                        != triangle.vertexIds.end()) {
                    return std::nullopt;
                }
            }
            for (Triangle &triangle : scene.triangles) {
                if (std::find(triangle.vertexIds.begin(),
                              triangle.vertexIds.end(), id)
                    == triangle.vertexIds.end()) {
                    continue;
                }
                const gp_Vec preferredWarp =
                    trianglePreferredWarp(triangle);
                std::replace(triangle.vertexIds.begin(),
                             triangle.vertexIds.end(), id, replacement);
                triangle.materialCoordinates = chartFor(
                    trianglePoints(triangle), preferredWarp);
            }
            for (Opening &opening : scene.openings) {
                std::replace(opening.orderedVertexIds.begin(),
                             opening.orderedVertexIds.end(), id,
                             replacement);
                for (auto &triangle : opening.capTriangleVertexIds) {
                    std::replace(triangle.begin(), triangle.end(), id,
                                 replacement);
                }
            }
            vertexIds[quantize(point)] = replacement;
            retiredVertexIds.insert(id);
            ribBoundaryAliases.emplace(id, replacement);
            return replacement;
        };

        for (const auto &[id, ribIndex] : pendingRibBoundaryVertices) {
            const bool hasFabricRib = std::ranges::any_of(
                ribOuterBoundaries,
                [ribIndex](const RibOuterBoundary &boundary) {
                    return boundary.ribIndex == ribIndex;
                });
            if (hasFabricRib
                && !canonicalRibBoundaryVertex(id, ribIndex)) {
                addError(
                    "Captured skin boundary could not be canonicalized to its rib-mesh station");
            }
        }

        // A fluid intake must follow actual Structure-owned fabric edges.
        // The two spanwise lips already come from the skin grid. Complete
        // their loops with the shorter exact rib-mesh paths between the same
        // captured lip corners, then author one deterministic virtual disk in
        // the capture's rectangular parameter chart. The disk is fluid
        // topology only; none of its triangles enters the fabric Structure.
        for (PendingIntake &intake : pendingIntakes) {
            bool boundaryOwned = true;
            for (const auto &[id, ribIndex] : {
                     std::pair{&intake.frontLipVertexIds.front(),
                               intake.firstRibIndex},
                     std::pair{&intake.frontLipVertexIds.back(),
                               intake.lastRibIndex},
                     std::pair{&intake.backLipVertexIds.front(),
                               intake.firstRibIndex},
                     std::pair{&intake.backLipVertexIds.back(),
                               intake.lastRibIndex}}) {
                const auto canonical =
                    canonicalRibBoundaryVertex(*id, ribIndex);
                boundaryOwned = canonical.has_value() && boundaryOwned;
                if (canonical) {
                    *id = *canonical;
                }
            }
            if (!boundaryOwned) {
                addError(
                    "Captured intake lip corner could not be canonicalized to its rib-mesh boundary");
                continue;
            }
            const auto rightRib = ribBoundaryPath(
                intake.frontLipVertexIds.back(),
                intake.backLipVertexIds.back(),
                intake.lastRibIndex);
            const auto leftRib = ribBoundaryPath(
                intake.backLipVertexIds.front(),
                intake.frontLipVertexIds.front(),
                intake.firstRibIndex);
            if (!rightRib || !leftRib) {
                addError(
                    "Captured intake lip corners do not share exact rib-mesh boundary paths");
                continue;
            }

            std::vector<ParametricOpeningVertex> boundary;
            boundary.reserve(
                intake.frontLipVertexIds.size()
                + rightRib->vertexIds.size()
                + intake.backLipVertexIds.size()
                + leftRib->vertexIds.size() - 4);
            const auto fraction = [](std::size_t index,
                                     std::size_t count) {
                return static_cast<double>(index)
                    / static_cast<double>(count - 1);
            };
            for (std::size_t index = 0;
                 index < intake.frontLipVertexIds.size(); ++index) {
                boundary.push_back({
                    intake.frontLipVertexIds[index],
                    fraction(index, intake.frontLipVertexIds.size()),
                    0.0});
            }
            for (std::size_t index = 1;
                 index < rightRib->vertexIds.size(); ++index) {
                boundary.push_back({
                    rightRib->vertexIds[index], 1.0,
                    fraction(index, rightRib->vertexIds.size())});
            }
            for (std::size_t index =
                     intake.backLipVertexIds.size() - 1;
                 index-- > 0;) {
                boundary.push_back({
                    intake.backLipVertexIds[index],
                    fraction(index, intake.backLipVertexIds.size()),
                    1.0});
            }
            for (std::size_t index = 1;
                 index + 1 < leftRib->vertexIds.size(); ++index) {
                boundary.push_back({
                    leftRib->vertexIds[index], 0.0,
                    1.0 - fraction(index, leftRib->vertexIds.size())});
            }

            std::vector<StableId> loop;
            loop.reserve(boundary.size());
            for (const ParametricOpeningVertex &vertex : boundary) {
                loop.push_back(vertex.id);
            }
            auto capTriangles =
                triangulateParametricOpeningBoundary(boundary);
            if (capTriangles.size() + 2 != loop.size()) {
                addError(
                    "Captured intake boundary could not be triangulated as an oriented disk");
                continue;
            }
            scene.openings.push_back({
                intake.id,
                std::move(loop),
                outsideRegionId,
                intake.cellRegionId,
                OpeningRole::Intake,
                std::move(capTriangles)});
        }
        if (maximumRibMeshSnapDeviationMillimetres
            > pointToleranceMillimetres) {
            std::ostringstream warning;
            warning << "Intake lip corners were canonicalized to the rib mesh by up to "
                    << maximumRibMeshSnapDeviationMillimetres << " mm";
            result.warnings.push_back(warning.str());
        }

        // Internal diagonal and mini-rib sheets do not divide the connected
        // cell volume. Their two oriented sides deliberately retain the same
        // containing-cell region.
        for (const CapturedStrip &strip : capturedStrips_) {
            if (strip.curveA.size() < 2 || strip.curveB.size() < 2) {
                continue;
            }
            const SurfaceRole role = strip.minirib ? SurfaceRole::MiniRib
                                                   : SurfaceRole::Diagonal;
            double minX = std::numeric_limits<double>::max();
            double maxX = std::numeric_limits<double>::lowest();
            for (const gp_Pnt &point : strip.curveA) {
                minX = std::min(minX, point.X());
                maxX = std::max(maxX, point.X());
            }
            const bool selfMirrored = std::abs(minX + maxX)
                < symmetryPlaneToleranceMillimetres;
            for (const bool mirror : {false, true}) {
                if (mirror && selfMirrored) {
                    continue;
                }
                const StableId sheetId = stableId(11, sheetOrdinal++);
                gp_XYZ centroid(0.0, 0.0, 0.0);
                for (const gp_Pnt &point : strip.curveA) {
                    centroid += place(point, mirror).XYZ();
                }
                for (const gp_Pnt &point : strip.curveB) {
                    centroid += place(point, mirror).XYZ();
                }
                const gp_Pnt stripCenter(
                    centroid
                    / static_cast<double>(strip.curveA.size()
                                          + strip.curveB.size()));
                auto nearest = cells.end();
                double nearestDistance = std::numeric_limits<double>::max();
                for (auto candidate = cells.begin(); candidate != cells.end();
                     ++candidate) {
                    if (candidate->first.mirror != mirror
                        && !selfMirrored) {
                        continue;
                    }
                    const double distance = cellCenter(candidate->first)
                                                .SquareDistance(stripCenter);
                    if (distance < nearestDistance) {
                        nearestDistance = distance;
                        nearest = candidate;
                    }
                }
                if (nearest == cells.end()) {
                    addError("Could not assign internal sheet '" + strip.label
                             + "' to a captured cell");
                    continue;
                }
                const StableId region = nearest->second.regionId;
                if (!strip.minirib) {
                    if (strip.curveA.size() != strip.curveB.size()) {
                        addError("Captured internal sheet '" + strip.label
                                 + "' lost its authored rung correspondence");
                        continue;
                    }
                    for (std::size_t index = 0;
                         index + 1 < strip.curveA.size(); ++index) {
                        const gp_Pnt a0 = place(strip.curveA[index], mirror);
                        const gp_Pnt a1 =
                            place(strip.curveA[index + 1], mirror);
                        const gp_Pnt b0 = place(strip.curveB[index], mirror);
                        const gp_Pnt b1 =
                            place(strip.curveB[index + 1], mirror);
                        const gp_Vec progression(a0, a1);
                        const gp_Vec rung(a0, b0);
                        const gp_Vec warp =
                            fabricSettings.at(role)->warpAxis
                                    == lep::SimWingMaterialAxis::Primary
                                ? progression
                                : rung;
                        addTriangle({a0, b0, b1}, warp,
                                    region, region, sheetId, role);
                        addTriangle({a0, b1, a1}, warp,
                                    region, region, sheetId, role);
                    }
                    continue;
                }

                // Mini-rib upper and lower boundaries have independently
                // authored profile samples and may meet at either profile
                // tip.
                // Zipper the two polylines by normalized arc length so every
                // boundary vertex survives without nearest-index repetition.
                // The one candidate at a coincident tip is a triangulation
                // artifact, not an authored zero-area fabric element.
                std::vector<gp_Pnt> curveA;
                std::vector<gp_Pnt> curveB;
                curveA.reserve(strip.curveA.size());
                curveB.reserve(strip.curveB.size());
                for (const gp_Pnt &point : strip.curveA) {
                    curveA.push_back(place(point, mirror));
                }
                for (const gp_Pnt &point : strip.curveB) {
                    curveB.push_back(place(point, mirror));
                }
                const auto progress = [](const std::vector<gp_Pnt> &curve) {
                    std::vector<double> values(curve.size(), 0.0);
                    for (std::size_t index = 1; index < curve.size(); ++index) {
                        values[index] = values[index - 1]
                            + curve[index - 1].Distance(curve[index]);
                    }
                    if (values.back() > Precision::Confusion()) {
                        for (double &value : values) {
                            value /= values.back();
                        }
                    }
                    return values;
                };
                const std::vector<double> progressA = progress(curveA);
                const std::vector<double> progressB = progress(curveB);
                if (progressA.back() == 0.0 || progressB.back() == 0.0) {
                    addError("Captured mini-rib '" + strip.label
                             + "' has a collapsed boundary");
                    continue;
                }
                const bool hasSharedLeadingTip =
                    curveA.front().Distance(curveB.front())
                    <= Precision::Confusion();
                const bool hasSharedTrailingTip =
                    curveA.back().Distance(curveB.back())
                    <= Precision::Confusion();
                const std::size_t trianglesBefore = scene.triangles.size();
                const auto addZipperTriangle = [&](std::array<gp_Pnt, 3> points,
                                                   const gp_Vec &progression,
                                                   const gp_Vec &rung) {
                    gp_Vec normal(points[0], points[1]);
                    normal.Cross(gp_Vec(points[0], points[2]));
                    if (normal.Magnitude() <= Precision::Confusion()) {
                        const auto countAt = [&](const gp_Pnt &tip) {
                            return static_cast<std::size_t>(std::count_if(
                                points.begin(), points.end(),
                                [&](const gp_Pnt &point) {
                                    return point.Distance(tip)
                                        <= Precision::Confusion();
                                }));
                        };
                        if ((hasSharedLeadingTip
                             && countAt(curveA.front()) >= 2)
                            || (hasSharedTrailingTip
                                && countAt(curveA.back()) >= 2)) {
                            return;
                        }
                        addError("Captured mini-rib '" + strip.label
                                 + "' produced a degenerate zipper triangle");
                        return;
                    }
                    // Independently sampled boundaries can leave a vanishing
                    // final wedge next to a shared profile tip. It carries no
                    // meaningful fabric area, while mass lumping would give
                    // its vertices extreme inverse masses and an unusable
                    // membrane system. Treat it like the coincident-tip
                    // zero-area candidate documented above.
                    constexpr double minimumTipSliverAreaSquareMillimetres =
                        0.05;
                    if (0.5 * normal.Magnitude()
                        < minimumTipSliverAreaSquareMillimetres) {
                        return;
                    }
                    const gp_Vec warp =
                        fabricSettings.at(role)->warpAxis
                                == lep::SimWingMaterialAxis::Primary
                            ? progression
                            : rung;
                    addTriangle(points, warp,
                                region, region, sheetId, role);
                };
                std::size_t indexA = 0;
                std::size_t indexB = 0;
                while (indexA + 1 < curveA.size()
                       || indexB + 1 < curveB.size()) {
                    const bool advanceA = indexB + 1 == curveB.size()
                        || (indexA + 1 < curveA.size()
                            && progressA[indexA + 1]
                                   < progressB[indexB + 1]);
                    const gp_Vec rung(curveA[indexA], curveB[indexB]);
                    if (advanceA) {
                        const gp_Vec progression(
                            curveA[indexA], curveA[indexA + 1]);
                        addZipperTriangle(
                            {curveA[indexA], curveB[indexB],
                             curveA[indexA + 1]},
                            progression, rung);
                        ++indexA;
                    } else {
                        const gp_Vec progression(
                            curveB[indexB], curveB[indexB + 1]);
                        addZipperTriangle(
                            {curveA[indexA], curveB[indexB],
                             curveB[indexB + 1]},
                            progression, rung);
                        ++indexB;
                    }
                }
                if (scene.triangles.size() == trianglesBefore) {
                    addError("Captured mini-rib '" + strip.label
                             + "' collapsed to zero fabric area");
                }
            }
        }

        // A zero-area outer rib is an authored sewn closure, not a free
        // material edge. The two skin paths generally have different source
        // samples (gnuC2 has 20 versus 15 vertices including endpoints).
        // Refine both boundary paths to the union of their monotone chordwise
        // samples, retaining every original vertex and splitting only its
        // incident boundary triangle. Corresponding refined vertices then
        // occupy exactly the same point while keeping distinct IDs; the two
        // already-welded path endpoints remain shared.
        struct BoundaryIncidence
        {
            std::size_t triangleIndex = 0;
            StableId from = invalidStableId;
            StableId to = invalidStableId;
            StableId negativeRegion = invalidStableId;
            StableId positiveRegion = invalidStableId;
        };
        using BoundaryEdge = std::array<StableId, 2>;
        const auto boundaryEdge = [](StableId first, StableId second) {
            return BoundaryEdge{std::min(first, second),
                                std::max(first, second)};
        };
        std::map<BoundaryEdge, std::vector<BoundaryIncidence>>
            separatingEdges;
        for (std::size_t triangleIndex = 0;
             triangleIndex < scene.triangles.size(); ++triangleIndex) {
            const Triangle &triangle = scene.triangles[triangleIndex];
            if (triangle.negativeSideRegionId
                == triangle.positiveSideRegionId) {
                continue;
            }
            for (std::size_t corner = 0; corner < 3; ++corner) {
                const StableId from = triangle.vertexIds[corner];
                const StableId to =
                    triangle.vertexIds[(corner + 1) % 3];
                separatingEdges[boundaryEdge(from, to)].push_back(
                    {triangleIndex, from, to,
                     triangle.negativeSideRegionId,
                     triangle.positiveSideRegionId});
            }
        }
        std::set<BoundaryEdge> openingEdges;
        for (const Opening &opening : scene.openings) {
            for (std::size_t index = 0;
                 index < opening.orderedVertexIds.size(); ++index) {
                openingEdges.insert(boundaryEdge(
                    opening.orderedVertexIds[index],
                    opening.orderedVertexIds[
                        (index + 1) % opening.orderedVertexIds.size()]));
            }
        }
        using BoundaryRegions = std::pair<StableId, StableId>;
        std::map<BoundaryRegions, std::vector<BoundaryIncidence>>
            unresolvedBoundaries;
        for (const auto &[edge, incidences] : separatingEdges) {
            if (incidences.size() == 1 && !openingEdges.contains(edge)) {
                const BoundaryIncidence &incidence = incidences.front();
                unresolvedBoundaries[
                    {incidence.negativeRegion,
                     incidence.positiveRegion}]
                    .push_back(incidence);
            }
        }

        struct SeamSample
        {
            double parameter = 0.0;
            Vec3 positionMeters;
            std::array<std::optional<StableId>, 2> existingIds;
        };
        struct EdgeInsertion
        {
            StableId id = invalidStableId;
            double fraction = 0.0;
        };
        const auto scenePosition = [&](StableId id) {
            const auto vertex = std::lower_bound(
                scene.vertices.begin(), scene.vertices.end(), id,
                [](const Vertex &candidate, StableId value) {
                    return candidate.id < value;
                });
            if (vertex == scene.vertices.end() || vertex->id != id) {
                throw std::logic_error(
                    "Collapsed seam references a missing scene vertex");
            }
            return vertex->positionMeters;
        };
        const auto squaredDistance = [](const Vec3 &first,
                                        const Vec3 &second) {
            const double x = first.x - second.x;
            const double y = first.y - second.y;
            const double z = first.z - second.z;
            return x * x + y * y + z * z;
        };
        const auto splitBoundaryTriangle =
            [&](StableId from,
                StableId to,
                const std::vector<EdgeInsertion> &insertions) {
            if (insertions.empty()) {
                return;
            }
            std::optional<std::size_t> foundTriangle;
            std::size_t foundCorner = 0;
            bool followsTriangle = false;
            for (std::size_t triangleIndex = 0;
                 triangleIndex < scene.triangles.size(); ++triangleIndex) {
                const Triangle &candidate = scene.triangles[triangleIndex];
                if (candidate.negativeSideRegionId
                    == candidate.positiveSideRegionId) {
                    continue;
                }
                for (std::size_t corner = 0; corner < 3; ++corner) {
                    const StableId first = candidate.vertexIds[corner];
                    const StableId second =
                        candidate.vertexIds[(corner + 1) % 3];
                    if ((first == from && second == to)
                        || (first == to && second == from)) {
                        if (foundTriangle) {
                            addError(
                                "Collapsed seam boundary edge has multiple material triangles");
                            return;
                        }
                        foundTriangle = triangleIndex;
                        foundCorner = corner;
                        followsTriangle = first == from;
                    }
                }
            }
            if (!foundTriangle) {
                addError(
                    "Collapsed seam boundary edge has no material triangle");
                return;
            }

            const Triangle original = scene.triangles[*foundTriangle];
            const std::size_t oppositeCorner = (foundCorner + 2) % 3;
            const StableId opposite = original.vertexIds[oppositeCorner];
            const Vec2 fromChart = followsTriangle
                ? original.materialCoordinates[foundCorner]
                : original.materialCoordinates[(foundCorner + 1) % 3];
            const Vec2 toChart = followsTriangle
                ? original.materialCoordinates[(foundCorner + 1) % 3]
                : original.materialCoordinates[foundCorner];
            const Vec2 oppositeChart =
                original.materialCoordinates[oppositeCorner];

            std::vector<StableId> edgeIds;
            std::vector<Vec2> edgeCharts;
            edgeIds.reserve(insertions.size() + 2);
            edgeCharts.reserve(insertions.size() + 2);
            edgeIds.push_back(from);
            edgeCharts.push_back(fromChart);
            for (const EdgeInsertion &insertion : insertions) {
                edgeIds.push_back(insertion.id);
                edgeCharts.push_back(
                    {fromChart.x
                         + insertion.fraction
                               * (toChart.x - fromChart.x),
                     fromChart.y
                         + insertion.fraction
                               * (toChart.y - fromChart.y)});
            }
            edgeIds.push_back(to);
            edgeCharts.push_back(toChart);

            for (std::size_t segment = 0;
                 segment + 1 < edgeIds.size(); ++segment) {
                Triangle piece = original;
                if (followsTriangle) {
                    piece.vertexIds = {
                        edgeIds[segment], edgeIds[segment + 1], opposite};
                    piece.materialCoordinates = {
                        edgeCharts[segment], edgeCharts[segment + 1],
                        oppositeChart};
                } else {
                    piece.vertexIds = {
                        edgeIds[segment + 1], edgeIds[segment], opposite};
                    piece.materialCoordinates = {
                        edgeCharts[segment + 1], edgeCharts[segment],
                        oppositeChart};
                }
                if (segment == 0) {
                    scene.triangles[*foundTriangle] = piece;
                } else {
                    piece.id = stableId(4, triangleOrdinal++);
                    scene.triangles.push_back(std::move(piece));
                }
            }
        };

        std::optional<StableId> collapsedSeamMaterialId;
        std::size_t seamOrdinal = 0;
        const auto requireCollapsedSeamMaterial = [&]() {
            if (collapsedSeamMaterialId) {
                return *collapsedSeamMaterialId;
            }
            const auto &material = settings.collapsedBoundarySeam;
            if (material.name.empty()
                || !std::isfinite(material.linearDensityKgPerMeter)
                || !(material.linearDensityKgPerMeter > 0.0)
                || !std::isfinite(material.axialStiffnessNewtons)
                || !(material.axialStiffnessNewtons > 0.0)) {
                addError(
                    "Collapsed sewn boundaries require explicit finite positive seam properties");
            }
            const StableId id = stableId(12, 0);
            scene.seamMaterials.push_back(
                {id, material.name,
                 material.linearDensityKgPerMeter,
                 material.axialStiffnessNewtons});
            collapsedSeamMaterialId = id;
            return id;
        };

        constexpr double collapsedBoundaryRelativeAreaTolerance = 1.0e-10;
        constexpr double seamParameterTolerance = 1.0e-12;
        for (const auto &[regions, incidences] : unresolvedBoundaries) {
            std::map<StableId, BoundaryIncidence> next;
            std::set<StableId> incoming;
            for (const BoundaryIncidence &incidence : incidences) {
                if (!next.emplace(incidence.from, incidence).second
                    || !incoming.insert(incidence.to).second) {
                    addError(
                        "Unresolved material boundary is branched or inconsistently wound");
                }
            }
            if (next.size() != incoming.size()) {
                addError("Unresolved material boundary is not a closed cycle");
                continue;
            }
            std::set<StableId> unvisited;
            for (const auto &[from, incidence] : next) {
                static_cast<void>(incidence);
                unvisited.insert(from);
            }
            while (!unvisited.empty()) {
                const StableId start = *unvisited.begin();
                std::vector<StableId> loop;
                StableId current = start;
                do {
                    if (!unvisited.erase(current)) {
                        addError(
                            "Unresolved material boundary cycles overlap");
                        break;
                    }
                    loop.push_back(current);
                    const auto found = next.find(current);
                    if (found == next.end()) {
                        addError(
                            "Unresolved material boundary cycle is incomplete");
                        break;
                    }
                    current = found->second.to;
                } while (current != start);
                if (current != start || loop.size() < 4) {
                    continue;
                }

                Vec3 minimum = scenePosition(loop.front());
                Vec3 maximum = minimum;
                Vec3 areaVector;
                for (std::size_t index = 0; index < loop.size(); ++index) {
                    const Vec3 first = scenePosition(loop[index]);
                    const Vec3 second =
                        scenePosition(loop[(index + 1) % loop.size()]);
                    minimum.x = std::min(minimum.x, first.x);
                    minimum.y = std::min(minimum.y, first.y);
                    minimum.z = std::min(minimum.z, first.z);
                    maximum.x = std::max(maximum.x, first.x);
                    maximum.y = std::max(maximum.y, first.y);
                    maximum.z = std::max(maximum.z, first.z);
                    areaVector.x +=
                        (first.y - second.y) * (first.z + second.z);
                    areaVector.y +=
                        (first.z - second.z) * (first.x + second.x);
                    areaVector.z +=
                        (first.x - second.x) * (first.y + second.y);
                }
                const double extentSquared = squaredDistance(minimum, maximum);
                const double relativeArea = 0.5 * std::hypot(
                    areaVector.x, areaVector.y, areaVector.z)
                    / extentSquared;
                if (!(extentSquared > 0.0)
                    || !std::isfinite(relativeArea)
                    || relativeArea
                        > collapsedBoundaryRelativeAreaTolerance) {
                    addError(
                        "Captured separating surface has an unauthored finite-area boundary");
                    continue;
                }

                std::size_t firstEnd = 0;
                std::size_t secondEnd = 1;
                double maximumDistanceSquared = -1.0;
                for (std::size_t first = 0; first < loop.size(); ++first) {
                    for (std::size_t second = first + 1;
                         second < loop.size(); ++second) {
                        const double distanceSquared = squaredDistance(
                            scenePosition(loop[first]),
                            scenePosition(loop[second]));
                        if (distanceSquared > maximumDistanceSquared) {
                            maximumDistanceSquared = distanceSquared;
                            firstEnd = first;
                            secondEnd = second;
                        }
                    }
                }
                std::vector<StableId> paths[2];
                for (std::size_t index = firstEnd;
                     index <= secondEnd; ++index) {
                    paths[0].push_back(loop[index]);
                }
                std::size_t index = firstEnd;
                paths[1].push_back(loop[index]);
                while (index != secondEnd) {
                    index = (index + loop.size() - 1) % loop.size();
                    paths[1].push_back(loop[index]);
                }
                if (paths[0].size() < 3 || paths[1].size() < 3) {
                    addError(
                        "Collapsed sewn boundary cannot be split into two fabric paths");
                    continue;
                }

                const Vec3 axisStart = scenePosition(paths[0].front());
                const Vec3 axisEnd = scenePosition(paths[0].back());
                const Vec3 axis{axisEnd.x - axisStart.x,
                                axisEnd.y - axisStart.y,
                                axisEnd.z - axisStart.z};
                const double axisSquared =
                    axis.x * axis.x + axis.y * axis.y + axis.z * axis.z;
                const auto parameter = [&](StableId id) {
                    const Vec3 point = scenePosition(id);
                    return ((point.x - axisStart.x) * axis.x
                            + (point.y - axisStart.y) * axis.y
                            + (point.z - axisStart.z) * axis.z)
                        / axisSquared;
                };
                bool monotone = axisSquared > 0.0;
                for (const auto &path : paths) {
                    double previous = -std::numeric_limits<double>::infinity();
                    for (const StableId id : path) {
                        const double value = parameter(id);
                        monotone = monotone && std::isfinite(value)
                            && value > previous + seamParameterTolerance;
                        previous = value;
                    }
                }
                if (!monotone) {
                    addError(
                        "Collapsed sewn boundary paths are not monotone between their farthest endpoints");
                    continue;
                }

                std::vector<SeamSample> samples;
                for (std::size_t path = 0; path < 2; ++path) {
                    for (const StableId id : paths[path]) {
                        SeamSample sample;
                        sample.parameter = parameter(id);
                        sample.positionMeters = scenePosition(id);
                        sample.existingIds[path] = id;
                        samples.push_back(sample);
                    }
                }
                std::ranges::sort(
                    samples,
                    [](const SeamSample &first, const SeamSample &second) {
                        if (first.parameter != second.parameter) {
                            return first.parameter < second.parameter;
                        }
                        return first.existingIds < second.existingIds;
                    });
                std::vector<SeamSample> mergedSamples;
                const auto sampleExistingId = [](const SeamSample &sample) {
                    return sample.existingIds[0]
                        ? *sample.existingIds[0]
                        : *sample.existingIds[1];
                };
                for (SeamSample sample : samples) {
                    if (!mergedSamples.empty()) {
                        SeamSample &previous = mergedSamples.back();
                        const StableId previousId =
                            sampleExistingId(previous);
                        const StableId sampleId =
                            sampleExistingId(sample);
                        if (previousId == sampleId) {
                            previous.existingIds[0] =
                                previous.existingIds[0]
                                    ? previous.existingIds[0]
                                    : sample.existingIds[0];
                            previous.existingIds[1] =
                                previous.existingIds[1]
                                    ? previous.existingIds[1]
                                    : sample.existingIds[1];
                            continue;
                        }
                        if (!(sample.parameter
                              > previous.parameter
                                  + seamParameterTolerance)) {
                            addError(
                                "Collapsed sewn boundary has ambiguous coincident samples");
                            break;
                        }
                    }
                    mergedSamples.push_back(std::move(sample));
                }
                if (mergedSamples.size() != loop.size()) {
                    continue;
                }

                const auto refinePath =
                    [&](std::size_t pathIndex) {
                    const auto &path = paths[pathIndex];
                    std::vector<double> pathParameters;
                    pathParameters.reserve(path.size());
                    for (const StableId id : path) {
                        pathParameters.push_back(parameter(id));
                    }
                    std::vector<std::vector<EdgeInsertion>> insertions(
                        path.size() - 1);
                    std::vector<StableId> refined;
                    refined.reserve(mergedSamples.size());
                    for (const SeamSample &sample : mergedSamples) {
                        if (sample.existingIds[pathIndex]) {
                            refined.push_back(*sample.existingIds[pathIndex]);
                            continue;
                        }
                        const auto upper = std::upper_bound(
                            pathParameters.begin(), pathParameters.end(),
                            sample.parameter);
                        if (upper == pathParameters.begin()
                            || upper == pathParameters.end()) {
                            addError(
                                "Collapsed seam sample lies outside its opposite path");
                            return std::vector<StableId>{};
                        }
                        const std::size_t segment =
                            static_cast<std::size_t>(
                                upper - pathParameters.begin() - 1);
                        const double fraction =
                            (sample.parameter - pathParameters[segment])
                            / (pathParameters[segment + 1]
                               - pathParameters[segment]);
                        const StableId id = stableId(2, scene.vertices.size());
                        scene.vertices.push_back(
                            {id, sample.positionMeters});
                        insertions[segment].push_back({id, fraction});
                        refined.push_back(id);
                    }
                    for (std::size_t segment = 0;
                         segment < insertions.size(); ++segment) {
                        splitBoundaryTriangle(
                            path[segment], path[segment + 1],
                            insertions[segment]);
                    }
                    return refined;
                };
                std::vector<StableId> firstChain = refinePath(0);
                std::vector<StableId> secondChain = refinePath(1);
                if (firstChain.size() != mergedSamples.size()
                    || secondChain.size() != mergedSamples.size()) {
                    continue;
                }
                scene.seams.push_back(
                    {stableId(13, seamOrdinal++),
                     requireCollapsedSeamMaterial(),
                     std::move(firstChain),
                     std::move(secondChain)});
            }
        }

        if (!result.errors.empty()) {
            finalizeFailure();
            return result;
        }

        // Lines are added after surface vertices so endpoints can be matched
        // without moving the authored line graph. Non-matches remain explicit
        // SuspensionJunction entities.
        std::set<StableId> fabricVertexIds;
        for (const Triangle &triangle : scene.triangles) {
            fabricVertexIds.insert(triangle.vertexIds.begin(),
                                   triangle.vertexIds.end());
        }
        const StableId pilotId = stableId(7, 0);
        scene.pilots.push_back(
            {pilotId,
             settings.pilot.name,
             settings.pilot.massKg,
             settings.pilot.centerOfMassPositionMeters,
             settings.pilot.linearVelocityMetersPerSecond,
             settings.pilot.bodyToWorld,
             settings.pilot.principalInertiaKgSquareMeters});

        std::map<QuantizedPoint, StableId, decltype(&pointLess)>
            junctionIds(&pointLess);
        std::map<std::tuple<int, StableId, std::size_t>, StableId>
            attachmentIds;
        std::map<StableId, Vec3> attachmentWorldPositions;
        const auto endpointAttachment = [&](const gp_Pnt &point) {
            const Vec3 metres{point.X() * 0.001,
                              point.Y() * 0.001,
                              point.Z() * 0.001};
            std::size_t harnessIndex = settings.pilot.harnessPoints.size();
            double harnessDistance =
                settings.pilot.endpointMatchToleranceMeters;
            for (std::size_t index = 0;
                 index < settings.pilot.harnessPoints.size();
                 ++index) {
                const Vec3 &candidate =
                    settings.pilot.harnessPoints[index].worldPositionMeters;
                const double distance = std::hypot(
                    metres.x - candidate.x,
                    metres.y - candidate.y,
                    metres.z - candidate.z);
                if (distance <= harnessDistance) {
                    harnessDistance = distance;
                    harnessIndex = index;
                }
            }
            if (harnessIndex < settings.pilot.harnessPoints.size()) {
                const auto key = std::make_tuple(2, pilotId, harnessIndex);
                const auto found = attachmentIds.find(key);
                if (found != attachmentIds.end()) {
                    return found->second;
                }
                const StableId id = stableId(9, scene.attachments.size());
                const auto &harness =
                    settings.pilot.harnessPoints[harnessIndex];
                scene.attachments.push_back(
                    {id,
                     AttachmentKind::PilotHarness,
                     invalidStableId,
                     pilotId,
                     harness.pilotLocalPositionMeters,
                     invalidStableId});
                attachmentIds.emplace(key, id);
                attachmentWorldPositions.emplace(
                    id, harness.worldPositionMeters);
                return id;
            }

            StableId surfaceId = invalidStableId;
            double surfaceDistance =
                settings.surfaceEndpointMatchToleranceMeters;
            for (const Vertex &vertex : scene.vertices) {
                if (!fabricVertexIds.contains(vertex.id)) {
                    continue;
                }
                const double distance = std::hypot(
                    metres.x - vertex.positionMeters.x,
                    metres.y - vertex.positionMeters.y,
                    metres.z - vertex.positionMeters.z);
                if (distance <= surfaceDistance) {
                    surfaceDistance = distance;
                    surfaceId = vertex.id;
                }
            }
            if (surfaceId != invalidStableId) {
                const auto key = std::make_tuple(1, surfaceId, 0U);
                const auto found = attachmentIds.find(key);
                if (found != attachmentIds.end()) {
                    return found->second;
                }
                const StableId id = stableId(9, scene.attachments.size());
                scene.attachments.push_back(
                    {id,
                     AttachmentKind::SurfaceVertex,
                     surfaceId,
                     invalidStableId,
                     {},
                     invalidStableId});
                attachmentIds.emplace(key, id);
                const auto vertex = std::lower_bound(
                    scene.vertices.begin(), scene.vertices.end(), surfaceId,
                    [](const Vertex &candidate, const StableId value) {
                        return candidate.id < value;
                    });
                attachmentWorldPositions.emplace(
                    id, vertex->positionMeters);
                return id;
            }

            const QuantizedPoint key = quantize(point);
            auto junction = junctionIds.find(key);
            StableId junctionId = invalidStableId;
            if (junction == junctionIds.end()) {
                junctionId = stableId(8, scene.suspensionJunctions.size());
                junctionIds.emplace(key, junctionId);
                scene.suspensionJunctions.push_back(
                    {junctionId,
                     metres,
                     settings.suspensionJunctionMassKg,
                     false});
            } else {
                junctionId = junction->second;
            }
            const auto attachmentKey =
                std::make_tuple(3, junctionId, 0U);
            const auto existing = attachmentIds.find(attachmentKey);
            if (existing != attachmentIds.end()) {
                return existing->second;
            }
            const StableId id = stableId(9, scene.attachments.size());
            scene.attachments.push_back(
                {id,
                 AttachmentKind::SuspensionJunction,
                 invalidStableId,
                 invalidStableId,
                 {},
                 junctionId});
            attachmentIds.emplace(attachmentKey, id);
            attachmentWorldPositions.emplace(id, metres);
            return id;
        };

        std::set<std::tuple<QuantizedPoint, QuantizedPoint, int, bool,
                            std::string>> emittedLines;
        std::size_t lineOrdinal = 0;
        double maximumLineRestLengthRebaseMeters = 0.0;
        for (const CapturedLine &line : capturedLines_) {
            const QuantizedSegment segment =
                quantizeSegment(line.start, line.end);
            const auto key = std::make_tuple(
                segment.first, segment.second, line.planIndex,
                line.brake, line.typeName);
            if (!emittedLines.insert(key).second) {
                continue;
            }
            const StableId firstAttachment = endpointAttachment(line.start);
            const StableId secondAttachment = endpointAttachment(line.end);
            if (firstAttachment == secondAttachment) {
                addError("Captured suspension segment collapses to one attachment: "
                         + line.label);
                continue;
            }
            const Vec3 &firstPosition =
                attachmentWorldPositions.at(firstAttachment);
            const Vec3 &secondPosition =
                attachmentWorldPositions.at(secondAttachment);
            const double restLengthMeters = std::hypot(
                secondPosition.x - firstPosition.x,
                secondPosition.y - firstPosition.y,
                secondPosition.z - firstPosition.z);
            maximumLineRestLengthRebaseMeters = std::max(
                maximumLineRestLengthRebaseMeters,
                std::abs(restLengthMeters
                         - line.start.Distance(line.end) * 0.001));
            scene.suspensionLines.push_back(
                {stableId(10, lineOrdinal++),
                 firstAttachment,
                 secondAttachment,
                 lineMaterialIds.at(line.typeName),
                 restLengthMeters,
                 line.brake ? SuspensionLineRole::Brake
                            : SuspensionLineRole::Suspension});
        }
        if (maximumLineRestLengthRebaseMeters > 1.0e-9) {
            std::ostringstream message;
            message.precision(6);
            message
                << "Discrete suspension endpoint matching rebased line rest"
                << " lengths by at most "
                << maximumLineRestLengthRebaseMeters * 1000.0 << " mm";
            result.warnings.push_back(message.str());
        }

        if (!retiredVertexIds.empty()) {
            std::erase_if(scene.vertices, [&](const Vertex &vertex) {
                return retiredVertexIds.contains(vertex.id);
            });
        }

        if (!result.errors.empty()) {
            finalizeFailure();
            return result;
        }
        result.validation = validateScene(scene);
        if (!result.validation.ok()) {
            for (const ValidationDiagnostic &diagnostic :
                 result.validation.diagnostics) {
                result.errors.push_back(diagnostic.message);
            }
            finalizeFailure();
            return result;
        }
        result.scene = std::move(scene);
        result.success = true;
        std::sort(result.warnings.begin(), result.warnings.end());
        return result;
    }

    // Writes the Playground simulation mesh: the coarse skin samples of
    // both wing halves welded into one node table, quad connectivity per
    // region, rib outline loops for internal webs, and the labelled
    // suspension lines. Everything is in millimetres, Z up, matching the
    // STEP model.
    bool writeSimMesh(const std::filesystem::path &path,
                      std::string &error) const
    {
        if (simRegions_.empty()) {
            error = "No simulation skin was captured";
            return false;
        }

        std::vector<gp_Pnt> nodes;
        std::unordered_map<QuantizedPoint, int, QuantizedPointHash> welded;
        const auto nodeIndex = [&](const gp_Pnt &point, bool mirror) {
            // The centre rib lands within float residue of the symmetry
            // plane (~0.1 mm); snap it on so the mirrored half welds to
            // the right half instead of floating 0.2 mm away from it.
            const double snappedX =
                std::abs(point.X()) < 0.5 ? 0.0 : point.X();
            const gp_Pnt placed(mirror ? -snappedX : snappedX,
                                point.Y(),
                                point.Z());
            const auto [entry, inserted] =
                welded.try_emplace(quantize(placed),
                                   static_cast<int>(nodes.size()));
            if (inserted) {
                nodes.push_back(placed);
            }
            return entry->second;
        };

        std::vector<std::array<int, 4>> quads;
        // Which skin each quad came from, parallel to `quads`, so the
        // Playground can draw the surfaces separately.
        std::vector<int> quadSurfaces;
        // Rib outline loops keyed by (ribIndex, mirrored side), assembled
        // from the regions' rib-side sample columns in station order.
        std::map<std::pair<int, bool>, std::vector<std::pair<int, int>>>
            ribColumns;

        for (const SimRegionCapture &region : simRegions_) {
            for (const bool mirror : {false, true}) {
                // A panel that straddles the centreline (the centre cell of
                // a wing whose innermost ribs sit either side of x=0) is its
                // own mirror image, and emitting the mirrored copy would lay
                // duplicate faces over it. Decided from the geometry rather
                // than from panelIndex == 1: wings whose innermost rib sits
                // *on* the centreline have a genuine second half-panel there
                // that must be mirrored.
                if (mirror && region.selfMirrored()) {
                    continue;
                }
                std::vector<std::array<int, simSpanColumnCount>> grid;
                grid.reserve(region.rows.size());
                for (const auto &row : region.rows) {
                    std::array<int, simSpanColumnCount> ids{};
                    for (int column = 0;
                         column < simSpanColumnCount;
                         ++column) {
                        ids[column] = nodeIndex(row[column], mirror);
                    }
                    grid.push_back(ids);
                }
                for (std::size_t rowIndex = 0;
                     rowIndex + 1 < grid.size();
                     ++rowIndex) {
                    for (int column = 0;
                         column + 1 < simSpanColumnCount;
                         ++column) {
                        quads.push_back({grid[rowIndex][column],
                                         grid[rowIndex][column + 1],
                                         grid[rowIndex + 1][column + 1],
                                         grid[rowIndex + 1][column]});
                        quadSurfaces.push_back(
                            static_cast<int>(region.surface));
                    }
                }
                // Span column 0 lies on rib panelIndex, the last column on
                // rib panelIndex-1. Both edges are registered: on wings
                // whose innermost rib sits on the centreline that rib is no
                // panel's column 0, and harvesting only column 0 dropped it
                // — leaving the two centre bays joined into one, which then
                // ballooned as a single double-width cell. Ribs shared by
                // neighbouring panels yield identical loops and are removed
                // by the uniqueLoops filter below.
                for (std::size_t rowIndex = 0;
                     rowIndex < grid.size();
                     ++rowIndex) {
                    ribColumns[{region.panelIndex, mirror}].emplace_back(
                        region.stations[rowIndex], grid[rowIndex][0]);
                    ribColumns[{region.panelIndex - 1, mirror}].emplace_back(
                        region.stations[rowIndex],
                        grid[rowIndex][simSpanColumnCount - 1]);
                }
            }
        }

        std::vector<std::vector<int>> ribLoops;
        // (rib index, mirrored) per emitted loop, so each one can be paired
        // with the hole outlines of the rib it actually is.
        std::vector<std::pair<int, bool>> ribLoopKeys;
        std::set<std::vector<int>> uniqueLoops;
        for (auto &[key, column] : ribColumns) {
            std::stable_sort(column.begin(),
                             column.end(),
                             [](const auto &a, const auto &b) {
                                 return a.first < b.first;
                             });
            std::vector<int> loop;
            for (const auto &[station, node] : column) {
                if (loop.empty() || loop.back() != node) {
                    loop.push_back(node);
                }
            }
            // The mirrored centre rib welds onto the unmirrored one; keep
            // one copy of each distinct loop.
            if (loop.size() >= 3 && uniqueLoops.insert(loop).second) {
                ribLoops.push_back(std::move(loop));
                ribLoopKeys.push_back(key);
            }
        }

        // Hole outlines for one rib, as polylines in model space. The wires
        // come from the same builder that cuts the STEP rib faces, so every
        // hole type — ellipse, triangle, rounded rectangle — is exact here
        // rather than reconstructed downstream. 1 mm of chord deflection is
        // far finer than the mesh that tests against them.
        const auto ribHoleOutlines =
            [this](int ribIndex,
                   bool mirror) -> std::vector<std::vector<gp_Pnt>> {
            std::vector<std::vector<gp_Pnt>> outlines;
            const auto captured = capturedRibs_.find(ribIndex);
            if (captured == capturedRibs_.end()
                || captured->second.holes.empty()) {
                return outlines;
            }
            RibFrame frame;
            if (!fitRibFrame(captured->second.planarPoints,
                             captured->second.spatialPoints,
                             frame)) {
                return outlines;
            }
            for (std::size_t holeIndex = 0;
                 holeIndex < captured->second.holes.size();
                 ++holeIndex) {
                int edgeCount = 0;
                std::vector<std::string> ignored;
                const TopoDS_Wire wire =
                    ribHoleWire(captured->second.holes[holeIndex],
                                frame,
                                ribIndex,
                                static_cast<int>(holeIndex) + 1,
                                edgeCount,
                                ignored);
                if (wire.IsNull()) {
                    continue;
                }
                std::vector<gp_Pnt> outline;
                for (BRepTools_WireExplorer explorer(wire);
                     explorer.More();
                     explorer.Next()) {
                    BRepAdaptor_Curve curve(explorer.Current());
                    GCPnts_QuasiUniformDeflection sampler(curve, 1.0);
                    if (!sampler.IsDone()) {
                        continue;
                    }
                    // The last point of one edge is the first of the next.
                    for (int point = 1; point < sampler.NbPoints(); ++point) {
                        gp_Pnt placed = sampler.Value(point);
                        if (mirror) {
                            placed.SetX(-placed.X());
                        }
                        outline.push_back(placed);
                    }
                }
                if (outline.size() >= 3) {
                    outlines.push_back(std::move(outline));
                }
            }
            return outlines;
        };

        const auto escaped = [](const std::string &text) {
            std::string result;
            for (const char character : text) {
                if (character == '"' || character == '\\') {
                    result.push_back('\\');
                }
                result.push_back(character);
            }
            return result;
        };

        std::ostringstream json;
        json.setf(std::ios::fixed);
        json.precision(3);
        json << "{\n\"version\": 1,\n\"unit\": \"mm\",\n\"nodes\": [";
        for (std::size_t index = 0; index < nodes.size(); ++index) {
            json << (index == 0 ? "\n" : ",\n")
                 << '[' << nodes[index].X() << ',' << nodes[index].Y()
                 << ',' << nodes[index].Z() << ']';
        }
        json << "\n],\n\"quads\": [";
        for (std::size_t index = 0; index < quads.size(); ++index) {
            const auto &quad = quads[index];
            json << (index == 0 ? "\n" : ",\n")
                 << '[' << quad[0] << ',' << quad[1] << ',' << quad[2]
                 << ',' << quad[3] << ']';
        }
        // Per-quad skin tag, indexing surfaceNames. Readers that predate
        // this field simply draw every quad alike.
        json << "\n],\n\"surfaceNames\": "
                "[\"extrados\",\"vent\",\"intrados\"],\n"
                "\"quadSurfaces\": [";
        for (std::size_t index = 0; index < quadSurfaces.size(); ++index) {
            json << (index == 0 ? "" : ",") << quadSurfaces[index];
        }
        json << "],\n\"ribLoops\": [";
        for (std::size_t index = 0; index < ribLoops.size(); ++index) {
            json << (index == 0 ? "\n" : ",\n") << '[';
            for (std::size_t node = 0;
                 node < ribLoops[index].size();
                 ++node) {
                json << (node == 0 ? "" : ",") << ribLoops[index][node];
            }
            json << ']';
        }
        // Hole outlines per rib loop, in the same order, so the Playground
        // can mesh the ribs with their holes cut out.
        json << "\n],\n\"ribHoles\": [";
        for (std::size_t index = 0; index < ribLoops.size(); ++index) {
            const std::vector<std::vector<gp_Pnt>> outlines =
                index < ribLoopKeys.size()
                    ? ribHoleOutlines(ribLoopKeys[index].first,
                                      ribLoopKeys[index].second)
                    : std::vector<std::vector<gp_Pnt>>{};
            json << (index == 0 ? "\n" : ",\n") << '[';
            for (std::size_t outline = 0; outline < outlines.size();
                 ++outline) {
                json << (outline == 0 ? "" : ",") << '[';
                for (std::size_t point = 0;
                     point < outlines[outline].size();
                     ++point) {
                    const gp_Pnt &placed = outlines[outline][point];
                    json << (point == 0 ? "" : ",") << '['
                         << placed.X() << ',' << placed.Y() << ','
                         << placed.Z() << ']';
                }
                json << ']';
            }
            json << ']';
        }

        // Internal diagonal/mini-rib sheets as paired sample rows: the
        // Playground ties each pair together so line load spreads across
        // ribs the way the real V-ribs make it.
        json << "\n],\n\"straps\": [";
        bool firstStrap = true;
        for (const CapturedStrip &strip : capturedStrips_) {
            if (strip.curveA.size() < 2 || strip.curveB.size() < 2) {
                continue;
            }
            constexpr int strapSamples = 6;
            const auto sampleAt = [](const std::vector<gp_Pnt> &curve,
                                     int index) {
                const std::size_t position =
                    static_cast<std::size_t>(
                        std::llround(static_cast<double>(index)
                                     * (curve.size() - 1)
                                     / (strapSamples - 1)));
                return curve[position];
            };
            // A strip lying symmetric on the centre plane is its own
            // mirror image.
            double minX = std::numeric_limits<double>::max();
            double maxX = std::numeric_limits<double>::lowest();
            for (const gp_Pnt &point : strip.curveA) {
                minX = std::min(minX, point.X());
                maxX = std::max(maxX, point.X());
            }
            const bool selfMirrored = std::abs(minX + maxX) < 0.5;
            for (const bool mirror : {false, true}) {
                if (mirror && selfMirrored) {
                    continue;
                }
                const double sign = mirror ? -1.0 : 1.0;
                json << (firstStrap ? "\n" : ",\n") << "{\"a\":[";
                for (int index = 0; index < strapSamples; ++index) {
                    const gp_Pnt point = sampleAt(strip.curveA, index);
                    json << (index == 0 ? "" : ",") << '['
                         << sign * point.X() << ',' << point.Y() << ','
                         << point.Z() << ']';
                }
                json << "],\"b\":[";
                for (int index = 0; index < strapSamples; ++index) {
                    const gp_Pnt point = sampleAt(strip.curveB, index);
                    json << (index == 0 ? "" : ",") << '['
                         << sign * point.X() << ',' << point.Y() << ','
                         << point.Z() << ']';
                }
                json << "]}";
                firstStrap = false;
            }
        }

        json << "\n],\n\"lines\": [";
        bool firstLine = true;
        std::unordered_set<QuantizedSemanticSegment,
                           QuantizedSemanticSegmentHash> emittedLines;
        for (const CapturedLine &line : capturedLines_) {
            const QuantizedSemanticSegment key{
                quantizeSegment(line.start, line.end),
                line.planIndex,
                line.brake,
            };
            if (!emittedLines.insert(key).second) {
                continue;
            }
            // The legacy drawing pass has already visited both sides of the
            // full wing. Mirroring here used to emit every physical segment
            // twice. Preserve first-seen order and metadata from that full
            // captured stream instead.
            json << (firstLine ? "\n" : ",\n")
                 << "{\"a\":[" << line.start.X() << ','
                 << line.start.Y() << ',' << line.start.Z()
                 << "],\"b\":[" << line.end.X() << ','
                 << line.end.Y() << ',' << line.end.Z()
                 << "],\"label\":\"" << escaped(line.label)
                 << "\",\"plan\":" << line.planIndex
                 << ",\"brake\":" << (line.brake ? 1 : 0) << '}';
            firstLine = false;
        }
        json << "\n]\n}\n";

        std::ofstream file(path, std::ios::binary);
        if (!file) {
            error = "Could not open " + path.string() + " for writing";
            return false;
        }
        file << json.str();
        if (!file) {
            error = "Could not write " + path.string();
            return false;
        }
        return true;
    }

private:
    std::vector<PanelSurface> panels_;
    std::map<int, CapturedRib> capturedRibs_;
    std::vector<CapturedLine> capturedLines_;
    std::vector<CapturedStrip> capturedStrips_;
    std::vector<SimRegionCapture> simRegions_;
    mutable std::vector<std::string> errors_;
    mutable std::vector<std::string> warnings_;
    bool captureLines_ = false;
    LineTag currentLineTag_;
    double maximumSourceDeviation_ = 0.0;
    double maximumLegacyAgreement_ = 0.0;
};

NurbsModel &model()
{
    static NurbsModel instance;
    return instance;
}

} // namespace

namespace lep {

void resetNurbsModel()
{
    model().reset();
}

NurbsWriteResult writeNurbsStep(const std::filesystem::path &path,
                                bool includeConstructionCurves)
{
    return model().writeStep(path, includeConstructionCurves);
}

bool writeSimMesh(const std::filesystem::path &path, std::string &error)
{
    return model().writeSimMesh(path, error);
}

SimWingSceneExportResult buildSimWingScene(
    const SimWingSceneExportSettings &settings)
{
    return model().buildSimWingScene(settings);
}

SimWingSceneExportResult writeSimWingScene(
    const std::filesystem::path &path,
    const SimWingSceneExportSettings &settings)
{
    SimWingSceneExportResult result = buildSimWingScene(settings);
    if (!result.success) {
        return result;
    }

    std::ostringstream encoded(std::ios::binary);
    std::string error;
    if (!simwing::fsi::writeScene(result.scene, encoded, &error)) {
        result.errors.push_back("Could not encode scene-v2: " + error);
    } else {
        std::ofstream file(path, std::ios::binary);
        if (!file) {
            result.errors.push_back(
                "Could not open " + path.string() + " for writing");
        } else {
            const std::string bytes = encoded.str();
            file.write(bytes.data(),
                       static_cast<std::streamsize>(bytes.size()));
            if (!file) {
                result.errors.push_back(
                    "Could not write " + path.string());
            }
        }
    }
    if (!result.errors.empty()) {
        std::sort(result.errors.begin(), result.errors.end());
        result.scene = {};
        result.success = false;
    }
    return result;
}

} // namespace lep

extern "C" void lep_nurbs_capture_panel(const double *u,
                                         const double *v,
                                         const double *w,
                                         const double *shapingHeight,
                                         const double *legacyTessellation,
                                         int panelIndex,
                                         int totalPointCount,
                                         int upperPointCount,
                                         int ventPointCount,
                                         int segmentCount,
                                         int includeVentSurface,
                                         int singleSkin)
{
    model().capturePanel(
        u,
        v,
        w,
        shapingHeight,
        legacyTessellation,
        panelIndex,
        totalPointCount,
        upperPointCount,
        ventPointCount,
        segmentCount,
        includeVentSurface != 0,
        singleSkin != 0);
}

extern "C" void lep_nurbs_capture_rib(const double *u,
                                       const double *v,
                                       const double *w,
                                       const double *holes,
                                       double chordCentimetres,
                                       int ribIndex,
                                       int totalPointCount)
{
    model().captureRib(
        u,
        v,
        w,
        holes,
        chordCentimetres,
        ribIndex,
        totalPointCount);
}

extern "C" void lep_nurbs_set_line_capture(int enabled)
{
    model().setLineCapture(enabled != 0);
}

extern "C" void lep_nurbs_set_line_tag(const char *label,
                                        int labelLength,
                                        int planIndex,
                                        int isBrake,
                                        const char *typeName,
                                        int typeNameLength,
                                        double diameterMm)
{
    model().setLineTag(
        label,
        labelLength,
        planIndex,
        isBrake != 0,
        typeName,
        typeNameLength,
        diameterMm);
}

extern "C" void lep_nurbs_capture_diagonal_strip(const char *kind,
                                                  int kindLength,
                                                  int index,
                                                  const double *xA,
                                                  const double *yA,
                                                  const double *zA,
                                                  const double *xB,
                                                  const double *yB,
                                                  const double *zB,
                                                  int pointCount,
                                                  int stride)
{
    model().captureDiagonalStrip(
        kind, kindLength, index, xA, yA, zA, xB, yB, zB, pointCount, stride);
}

extern "C" void lep_nurbs_capture_miniribs(const double *x,
                                            const double *y,
                                            const double *z,
                                            const int *np,
                                            const double *rib,
                                            int ribCount,
                                            int singleSkin)
{
    model().captureMiniribs(x, y, z, np, rib, ribCount, singleSkin != 0);
}

extern "C" void lep_nurbs_capture_line(double x1,
                                        double y1,
                                        double z1,
                                        double x2,
                                        double y2,
                                        double z2,
                                        int colorIndex)
{
    model().captureLine({
        modelPoint(x1, y1, z1),
        modelPoint(x2, y2, z2),
        colorIndex,
    });
}

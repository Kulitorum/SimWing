cmake_minimum_required(VERSION 3.21)

# Ship binaries that run several macOS releases back rather than only on the
# build host: without an explicit deployment target the minimum OS defaults to
# the CI runner's version (macOS 15), so the app refuses to launch on older
# systems. Must be set before project(). The CI release build additionally
# passes -DCMAKE_OSX_ARCHITECTURES="arm64;x86_64" for a universal binary.
if(APPLE AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0" CACHE STRING
        "Minimum macOS version")
endif()

project(
    LEparagliding
    VERSION 0.4.1
    DESCRIPTION "Qt C++ port of LEparagliding 3.28"
    LANGUAGES C CXX
)

set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(WIN32 AND "${CMAKE_PREFIX_PATH}" STREQUAL "")
    file(GLOB qt_candidates LIST_DIRECTORIES true "C:/Qt/6.*/msvc2022_64")
    list(SORT qt_candidates COMPARE NATURAL ORDER DESCENDING)
    foreach(qt_candidate IN LISTS qt_candidates)
        if(EXISTS "${qt_candidate}/lib/cmake/Qt6/Qt6Config.cmake")
            set(CMAKE_PREFIX_PATH "${qt_candidate}" CACHE PATH "Qt installation" FORCE)
            message(STATUS "Auto-detected Qt at ${CMAKE_PREFIX_PATH}")
            break()
        endif()
    endforeach()
endif()

find_package(Qt6 6.5 REQUIRED COMPONENTS Widgets OpenGL OpenGLWidgets Xml)
find_package(OpenGL REQUIRED)

# Open CASCADE Technology is the authoritative 3D kernel for model creation,
# STEP serialization, tessellation, and viewport rendering. The Windows
# default mirrors C:\CODE\cobod-slicer and can be overridden for another OCCT
# build with -DLEP_OCCT_ROOT=<path>.
set(LEP_OCCT_ROOT "" CACHE PATH "Open CASCADE build or installation root")
if(WIN32)
    if(NOT LEP_OCCT_ROOT)
        foreach(occt_candidate
                "C:/OpenCASCADE-8.0/build2"
                "C:/OpenCASCADE-8.0/build"
                "C:/OpenCASCADE-7.9.0-vc14-64/opencascade-7.9.0")
            if(EXISTS "${occt_candidate}/inc/Standard.hxx")
                set(LEP_OCCT_ROOT "${occt_candidate}" CACHE PATH
                    "Open CASCADE build or installation root" FORCE)
                break()
            endif()
        endforeach()
    endif()
    if(EXISTS "${LEP_OCCT_ROOT}/inc/Standard.hxx")
        # Classic OCCT development layout (local builds under C:/OpenCASCADE).
        set(LEP_OCCT_INCLUDE_DIR "${LEP_OCCT_ROOT}/inc")
        set(LEP_OCCT_LIBRARY_DIR "${LEP_OCCT_ROOT}/win64/vc14/lib")
        set(LEP_OCCT_BINARY_DIR "${LEP_OCCT_ROOT}/win64/vc14/bin")
    else()
        # OCCT installed with INSTALL_DIR_LAYOUT=Unix (the CI build).
        set(LEP_OCCT_INCLUDE_DIR "${LEP_OCCT_ROOT}/include/opencascade")
        set(LEP_OCCT_LIBRARY_DIR "${LEP_OCCT_ROOT}/lib")
        set(LEP_OCCT_BINARY_DIR "${LEP_OCCT_ROOT}/bin")
    endif()
else()
    find_package(OpenCASCADE REQUIRED)
    set(LEP_OCCT_INCLUDE_DIR "${OpenCASCADE_INCLUDE_DIR}")
    set(LEP_OCCT_LIBRARY_DIR "${OpenCASCADE_LIBRARY_DIR}")
    set(LEP_OCCT_BINARY_DIR "${OpenCASCADE_BINARY_DIR}")
endif()

if(NOT EXISTS "${LEP_OCCT_INCLUDE_DIR}/Standard.hxx")
    message(FATAL_ERROR
        "Open CASCADE headers were not found. Set LEP_OCCT_ROOT to an OCCT "
        "build/install root (current value: '${LEP_OCCT_ROOT}').")
endif()

set(LEP_OCCT_TOOLKITS
    TKernel
    TKMath
    TKG2d
    TKG3d
    TKGeomBase
    TKBRep
    TKGeomAlgo
    TKTopAlgo
    TKMesh
    TKService
    TKV3d
    TKOpenGl
    TKDE
    TKXSBase
    TKDESTEP
    TKCDF
    TKLCAF
    TKCAF
    TKXCAF
    TKBinL
    TKBin
    TKBinXCAF
)

add_library(lep_occt INTERFACE)
target_include_directories(lep_occt SYSTEM BEFORE INTERFACE
    "${LEP_OCCT_INCLUDE_DIR}"
)
target_compile_definitions(lep_occt INTERFACE
    OCCT_NO_DEBUG
    OCCT_NO_DEPRECATED
)
foreach(occt_toolkit IN LISTS LEP_OCCT_TOOLKITS)
    if(WIN32)
        set(occt_library
            "${LEP_OCCT_LIBRARY_DIR}/${occt_toolkit}.lib")
    else()
        # find_library skips the search when the result variable is already
        # set, so clear it or every toolkit after the first reuses TKernel.
        unset(occt_library)
        find_library(occt_library
            NAMES "${occt_toolkit}"
            HINTS "${LEP_OCCT_LIBRARY_DIR}"
            REQUIRED
            NO_CACHE
        )
    endif()
    if(NOT EXISTS "${occt_library}")
        message(FATAL_ERROR
            "Required Open CASCADE toolkit not found: ${occt_library}")
    endif()
    target_link_libraries(lep_occt INTERFACE "${occt_library}")
endforeach()
message(STATUS "Using Open CASCADE at ${LEP_OCCT_ROOT}")

if(NOT WIN32)
    # OCCT usually lives in a custom prefix (see README); bake its library
    # directory into the build-tree rpath so binaries and tests run directly.
    list(APPEND CMAKE_BUILD_RPATH "${LEP_OCCT_LIBRARY_DIR}")
endif()

file(GLOB F2C_RUNTIME_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/libf2c/*.c"
)

add_library(f2c_runtime STATIC ${F2C_RUNTIME_SOURCES})
target_include_directories(f2c_runtime PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/libf2c"
)
target_compile_definitions(f2c_runtime PRIVATE
    SIGNED_ZEROS
    WANT_LEAD_0
)
if(WIN32)
    target_compile_definitions(f2c_runtime PRIVATE
        USE_CLOCK
        MSDOS
        NO_ONEXIT
        NO_My_ctype
        NO_ISATTY
        _CRT_SECURE_NO_WARNINGS
    )
    target_compile_options(f2c_runtime PRIVATE /W0)
else()
    target_compile_options(f2c_runtime PRIVATE -w)
    if(UNIX)
        target_link_libraries(f2c_runtime PUBLIC m)
    endif()
endif()

# Vendored SoftWingLab XPBD physics core (src/softwing, see its README.md
# for provenance): dependency-free soft-body/membrane/pneumatics/suspension
# simulation used by the Playground tab's toy live-wing model.
add_library(softwing_core STATIC
    src/softwing/canopy.cpp
    src/softwing/cell_volume.cpp
    src/softwing/contact.cpp
    src/softwing/membrane.cpp
    src/softwing/parallel.cpp
    src/softwing/pneumatics.cpp
    src/softwing/soft_body.cpp
    src/softwing/soft_body_checkpoint_persistence.cpp
    src/softwing/suspension.cpp
    src/softwing/suspension_checkpoint_persistence.cpp
)
target_include_directories(softwing_core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/softwing/include"
)
if(MSVC)
    target_compile_options(softwing_core PRIVATE /W0 /bigobj)
else()
    target_compile_options(softwing_core PRIVATE -w)
endif()

# SimWing remake foundations. These targets deliberately do not link any
# inherited Playground simulation library: only the structural adapter uses
# the isolated XPBD core.
add_library(simwing_scene STATIC
    src/fsi/scene.cpp
    src/fsi/scene.h
)
target_include_directories(simwing_scene PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)

add_library(lep_nurbs_model STATIC
    src/model/nurbs_model.cpp
    src/model/nurbs_model.h
)
target_include_directories(lep_nurbs_model PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/model"
)
target_link_libraries(lep_nurbs_model PUBLIC
    lep_occt
    simwing_scene
)

add_library(simwing_structure STATIC
    src/fsi/structure.cpp
    src/fsi/structure.h
    src/fsi/structure_checkpoint_persistence.cpp
    src/fsi/structure_checkpoint_persistence.h
)
target_include_directories(simwing_structure PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_structure PUBLIC softwing_core)

add_library(simwing_transfer STATIC
    src/fsi/transfer.cpp
    src/fsi/transfer.h
)
target_include_directories(simwing_transfer PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_transfer PUBLIC simwing_structure)

add_library(simwing_coupling STATIC
    src/fsi/coupling.cpp
    src/fsi/coupling.h
)
target_include_directories(simwing_coupling PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_coupling PUBLIC simwing_transfer)

add_library(simwing_scene_structure STATIC
    src/fsi/scene_structure.cpp
    src/fsi/scene_structure.h
)
target_include_directories(simwing_scene_structure PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_scene_structure PUBLIC
    simwing_scene
    simwing_structure
)

add_library(simwing_viewer_protocol STATIC
    src/viewer/viewer_protocol.cpp
    src/viewer/viewer_protocol.h
)
target_include_directories(simwing_viewer_protocol PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/viewer"
)

add_library(simwing_structure_frame STATIC
    src/viewer/structure_frame.cpp
    src/viewer/structure_frame.h
)
target_include_directories(simwing_structure_frame PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/viewer"
)
target_link_libraries(simwing_structure_frame PUBLIC
    simwing_scene_structure
    simwing_structure
    simwing_viewer_protocol
)

add_library(simwing_viewer_geometry STATIC
    src/viewer/vector_glyphs.cpp
    src/viewer/vector_glyphs.h
)
target_include_directories(simwing_viewer_geometry PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/viewer"
)
target_link_libraries(simwing_viewer_geometry PUBLIC
    simwing_viewer_protocol
)

add_library(simwing_canonical_case STATIC
    src/fsi/canonical_case.cpp
    src/fsi/canonical_case.h
)
target_include_directories(simwing_canonical_case PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_canonical_case PUBLIC
    simwing_structure_frame
)

add_library(simwing_fluid STATIC
    src/fsi/fluid/advection.cpp
    src/fsi/fluid/advection.h
    src/fsi/fluid/checkpoint.cpp
    src/fsi/fluid/checkpoint.h
    src/fsi/fluid/checkpoint_persistence.cpp
    src/fsi/fluid/diffusion.cpp
    src/fsi/fluid/diffusion.h
    src/fsi/fluid/evolution.cpp
    src/fsi/fluid/evolution.h
    src/fsi/fluid/grid.cpp
    src/fsi/fluid/grid.h
    src/fsi/fluid/interface_jump.cpp
    src/fsi/fluid/interface_jump.h
    src/fsi/fluid/moving_control_volume.cpp
    src/fsi/fluid/moving_control_volume.h
    src/fsi/fluid/moving_interface.cpp
    src/fsi/fluid/moving_interface.h
    src/fsi/fluid/planar_cut_surface.cpp
    src/fsi/fluid/planar_cut_surface.h
    src/fsi/fluid/porous_interface.cpp
    src/fsi/fluid/porous_interface.h
    src/fsi/fluid/porous_flow.cpp
    src/fsi/fluid/porous_flow.h
    src/fsi/fluid/projection.cpp
    src/fsi/fluid/projection.h
    src/fsi/fluid/projected_advection.cpp
    src/fsi/fluid/projected_advection.h
)
target_include_directories(simwing_fluid PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)

add_library(simwing_fluid_frame STATIC
    src/viewer/fluid_frame.cpp
    src/viewer/fluid_frame.h
    src/viewer/pressure_jump_frame.cpp
    src/viewer/pressure_jump_frame.h
)
target_include_directories(simwing_fluid_frame PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/viewer"
)
target_link_libraries(simwing_fluid_frame PUBLIC
    simwing_fluid
    simwing_viewer_protocol
)

add_library(simwing_pressure_jump_case STATIC
    src/fsi/pressure_jump_case.cpp
    src/fsi/pressure_jump_case.h
)
target_include_directories(simwing_pressure_jump_case PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_pressure_jump_case PUBLIC
    simwing_fluid_frame
)

add_library(simwing_porous_flow_case STATIC
    src/fsi/porous_flow_case.cpp
    src/fsi/porous_flow_case.h
)
target_include_directories(simwing_porous_flow_case PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_porous_flow_case PUBLIC
    simwing_fluid_frame
)

add_library(simwing_periodic_flow_case STATIC
    src/fsi/periodic_flow_case.cpp
    src/fsi/periodic_flow_case.h
    src/fsi/periodic_flow_checkpoint.cpp
)
target_include_directories(simwing_periodic_flow_case PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_periodic_flow_case PUBLIC
    simwing_fluid_frame
)

add_library(simwing_worker_control_protocol STATIC
    src/fsi/worker_control_protocol.cpp
    src/fsi/worker_control_protocol.h
)
target_include_directories(simwing_worker_control_protocol PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)

add_library(simwing_worker_control_stream STATIC
    src/fsi/worker_control_stream.cpp
    src/fsi/worker_control_stream.h
)
target_include_directories(simwing_worker_control_stream PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_worker_control_stream PUBLIC
    simwing_worker_control_protocol
)

add_library(simwing_worker_control_session STATIC
    src/fsi/worker_control_session.cpp
    src/fsi/worker_control_session.h
)
target_include_directories(simwing_worker_control_session PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_worker_control_session PUBLIC
    simwing_viewer_protocol
    simwing_worker_control_protocol
)

add_library(simwing_periodic_flow_control STATIC
    src/fsi/periodic_flow_control.cpp
    src/fsi/periodic_flow_control.h
)
target_include_directories(simwing_periodic_flow_control PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_periodic_flow_control PUBLIC
    simwing_periodic_flow_case
    simwing_worker_control_session
)

add_library(simwing_fluid_structure_bridge STATIC
    src/fsi/face_resolved_bridge.cpp
    src/fsi/face_resolved_bridge.h
    src/fsi/fluid_structure_bridge.cpp
    src/fsi/fluid_structure_bridge.h
)
target_include_directories(simwing_fluid_structure_bridge PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_fluid_structure_bridge PUBLIC
    simwing_fluid
    simwing_transfer
)

add_library(simwing_piston_case STATIC
    src/fsi/piston_case.cpp
    src/fsi/piston_case.h
)
target_include_directories(simwing_piston_case PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_piston_case PUBLIC
    simwing_coupling
    simwing_fluid_structure_bridge
    simwing_structure_frame
)

add_library(simwing_porous_sheet_case STATIC
    src/fsi/porous_sheet_case.cpp
    src/fsi/porous_sheet_case.h
    src/fsi/porous_sheet_checkpoint_detail.h
    src/fsi/porous_sheet_checkpoint_persistence.cpp
    src/fsi/porous_sheet_checkpoint_persistence.h
)
target_include_directories(simwing_porous_sheet_case PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_porous_sheet_case PUBLIC
    simwing_coupling
    simwing_fluid_structure_bridge
    simwing_structure
    simwing_structure_frame
)

add_library(simwing_open_piston_case STATIC
    src/fsi/open_piston_case.cpp
    src/fsi/open_piston_case.h
    src/fsi/open_piston_checkpoint_persistence.cpp
    src/fsi/open_piston_checkpoint_persistence.h
)
target_include_directories(simwing_open_piston_case PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_open_piston_case PUBLIC
    simwing_coupling
    simwing_fluid_structure_bridge
    simwing_structure_frame
)

add_library(simwing_open_piston_control STATIC
    src/fsi/open_piston_control.cpp
    src/fsi/open_piston_control.h
)
target_include_directories(simwing_open_piston_control PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fsi"
)
target_link_libraries(simwing_open_piston_control PUBLIC
    simwing_open_piston_case
    simwing_worker_control_session
)

# The Playground's wing: mesh parsing, refinement, SoftBody assembly and the
# per-frame step, with no widget attached. A library rather than GUI sources
# so the headless solver benchmark drives the identical simulation.
add_library(playground_pressure STATIC
    src/gui/playground_pressure_solve.cpp
    src/gui/playground_pressure_solve.h
)
target_include_directories(playground_pressure PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gui)

add_library(playground_cell_air STATIC
    src/gui/playground_cell_air.cpp
    src/gui/playground_cell_air.h
)
target_link_libraries(playground_cell_air PUBLIC softwing_core)
target_include_directories(playground_cell_air PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gui)

add_library(playground_sim STATIC
    src/gui/playground_sim.cpp
    src/gui/playground_sim.h
    src/gui/playground_metrics.cpp
    src/gui/playground_metrics.h
)
add_library(playground_contact STATIC
    src/gui/playground_contact.cpp
    src/gui/playground_contact.h
)
target_link_libraries(playground_contact PUBLIC softwing_core)
target_include_directories(playground_contact PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/gui"
)
target_link_libraries(playground_sim PUBLIC
    playground_cell_air
    playground_pressure
    Qt6::Core
    softwing_core
    playground_contact
)
target_include_directories(playground_sim PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/gui"
)

# Vendored XFLR5 v6.62 (third_party/xflr5, see PROVENANCE.md): the XFoil
# solver library plus the complete xflr5v6 application compiled as a static
# library. The GUI hosts XFLR5's MainFrame in a tab; XFLR5's own
# QApplication bootstrap (globals/main.cpp, globals/xflr5app.*) is excluded.
add_library(xfoil STATIC
    third_party/xflr5/XFoil-lib/xfoil.cpp
)
target_include_directories(xfoil PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/xflr5/XFoil-lib"
)
target_compile_definitions(xfoil PUBLIC XFOILLIB_STATIC)
target_link_libraries(xfoil PUBLIC Qt6::Core)

# Explicit source list generated from the upstream qmake manifests: the
# v6.62 tree contains orphaned duplicate files that qmake never compiled,
# so a glob would sweep in code that does not build (see the header of the
# generated file for how to regenerate).
set(LEP_XFLR5_V6_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/xflr5/xflr5v6")
include(cmake/xflr5_sources.cmake)

add_library(xflr5core STATIC ${XFLR5_CORE_SOURCES})
set_target_properties(xflr5core PROPERTIES AUTOMOC ON)
target_include_directories(xflr5core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/xflr5/xflr5v6"
)
target_link_libraries(xflr5core PUBLIC
    xfoil
    Qt6::Widgets
    Qt6::OpenGL
    Qt6::OpenGLWidgets
    Qt6::Xml
)
if(WIN32)
    target_link_libraries(xflr5core PUBLIC opengl32 user32)
endif()
if(MSVC)
    target_compile_options(xfoil PRIVATE /W0 /bigobj)
    target_compile_options(xflr5core PRIVATE /W0 /bigobj)
    foreach(xflr5_target xfoil xflr5core)
        target_compile_definitions(${xflr5_target} PRIVATE
            _CRT_SECURE_NO_WARNINGS
            _USE_MATH_DEFINES
            NOMINMAX
        )
    endforeach()
else()
    target_compile_options(xfoil PRIVATE -w)
    target_compile_options(xflr5core PRIVATE -w)
endif()

# XFLR5's Qt resources (icons, dark stylesheet, GL shaders, script/text
# templates) must be compiled into the executable: rcc registration relies
# on static initializers that the linker may drop from a static library.
set(XFLR5_RESOURCES
    third_party/xflr5/xflr5v6/resources/images.qrc
    third_party/xflr5/xflr5v6/resources/qss.qrc
    third_party/xflr5/xflr5v6/xfl3d/resources/shaders.qrc
    third_party/xflr5/xflr5v6/xfl3d/resources/text_files.qrc
    third_party/xflr5/xflr5v6/xflscript/script.qrc
    third_party/xflr5/xflr5v6/xflwidgets/textfiles.qrc
)

add_library(leparagliding_engine_runtime STATIC
    src/engine/engine_paths.cpp
    src/engine/engine_paths.h
    src/engine/input_migration.cpp
    src/engine/input_migration.h
    src/flatparts/flat_capture.cpp
    src/flatparts/flat_capture.h
    src/legacy/leparagliding_core.cpp
)
target_include_directories(leparagliding_engine_runtime PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/engine"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/flatparts"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/model"
)
target_link_libraries(leparagliding_engine_runtime PUBLIC
    f2c_runtime
    lep_nurbs_model
)
add_executable(leparagliding-engine
    src/engine/main.cpp
)
target_link_libraries(leparagliding-engine PRIVATE
    leparagliding_engine_runtime
)
if(MSVC)
    set_source_files_properties(src/legacy/leparagliding_core.cpp PROPERTIES
        COMPILE_OPTIONS "/bigobj;/W0;/Zc:strictStrings-"
    )
else()
    # The translated core is C89 wrapped in extern "C" guards; compiling it
    # as C keeps GCC/Clang from rejecting f2c's literal-to-char* idioms.
    set_source_files_properties(src/legacy/leparagliding_core.cpp PROPERTIES
        LANGUAGE C
        COMPILE_OPTIONS "-w"
    )
endif()

# Flat-part model, nester and sheet writers, shared by the Print tab and the
# nesting bench.
add_library(flatparts STATIC
    src/flatparts/flat_parts.cpp
    src/flatparts/flat_parts.h
    src/flatparts/nesting.cpp
    src/flatparts/nesting.h
    src/flatparts/sheet_export.cpp
    src/flatparts/sheet_export.h
)
target_link_libraries(flatparts PUBLIC Qt6::Gui)
target_include_directories(flatparts PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/flatparts"
)

add_executable(nesting-bench tests/nesting_bench.cpp)
target_link_libraries(nesting-bench PRIVATE flatparts)

add_executable(reference-report-compare
    tests/reference_report_compare.cpp
)

add_executable(preprocessor-reference-compare
    tests/preprocessor_reference_compare.cpp
    src/gui/geometry_preprocessor.cpp
)
target_include_directories(preprocessor-reference-compare PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/gui"
)

add_executable(section1-curves-test
    tests/section1_curves_test.cpp
    src/gui/section1_curves.cpp
)
target_include_directories(section1-curves-test PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/gui"
)

add_executable(spline-fit-test
    tests/spline_fit_test.cpp
    src/gui/spline_fit.cpp
)
target_include_directories(spline-fit-test PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/gui"
)

add_executable(section-grid-test
    tests/section_grid_test.cpp
    src/gui/section_grid.cpp
)
target_include_directories(section-grid-test PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/gui"
)

add_executable(airfoil-file-test
    tests/airfoil_file_test.cpp
    src/gui/airfoil_file.cpp
    src/gui/spline_fit.cpp
)
target_include_directories(airfoil-file-test PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/gui"
)

add_executable(playground-metrics-test
    tests/playground_metrics_test.cpp
)
target_link_libraries(playground-metrics-test PRIVATE playground_sim)

add_executable(playground-pressure-solve-test
    tests/playground_pressure_solve_test.cpp
)
target_link_libraries(playground-pressure-solve-test PRIVATE playground_pressure)

add_executable(playground-cell-air-test
    tests/playground_cell_air_test.cpp
)
target_link_libraries(playground-cell-air-test PRIVATE playground_cell_air)

add_executable(softwing-material-test
    tests/softwing_material_test.cpp
)
target_link_libraries(softwing-material-test PRIVATE softwing_core)

add_executable(softwing-cell-volume-test
    tests/softwing_cell_volume_test.cpp
)
target_link_libraries(softwing-cell-volume-test PRIVATE softwing_core)

add_executable(softwing-suspension-checkpoint-test
    tests/softwing_suspension_checkpoint_test.cpp
)
target_link_libraries(softwing-suspension-checkpoint-test
    PRIVATE softwing_core)

add_executable(softwing-checkpoint-persistence-test
    tests/softwing_checkpoint_persistence_test.cpp
)
target_link_libraries(softwing-checkpoint-persistence-test
    PRIVATE softwing_core)

add_executable(simwing-scene-test
    tests/simwing_scene_test.cpp
)
target_link_libraries(simwing-scene-test PRIVATE simwing_scene)

add_executable(simwing-model-scene-export-test
    tests/simwing_model_scene_export_test.cpp
)
target_link_libraries(simwing-model-scene-export-test PRIVATE
    lep_nurbs_model
)

add_executable(simwing-model-scene-real-export-test
    tests/simwing_model_scene_real_export_test.cpp
)
target_link_libraries(simwing-model-scene-real-export-test PRIVATE
    leparagliding_engine_runtime
    simwing_structure_frame
)

add_executable(simwing-structure-test
    tests/simwing_structure_test.cpp
)
target_link_libraries(simwing-structure-test PRIVATE simwing_structure)

add_executable(simwing-scene-structure-test
    tests/simwing_scene_structure_test.cpp
)
target_link_libraries(simwing-scene-structure-test
    PRIVATE simwing_scene_structure)

add_executable(simwing-viewer-protocol-test
    tests/simwing_viewer_protocol_test.cpp
)
target_link_libraries(simwing-viewer-protocol-test
    PRIVATE simwing_viewer_protocol)

add_executable(simwing-vector-glyphs-test
    tests/simwing_vector_glyphs_test.cpp
)
target_link_libraries(simwing-vector-glyphs-test
    PRIVATE simwing_viewer_geometry)

add_executable(simwing-structure-frame-test
    tests/simwing_structure_frame_test.cpp
)
target_link_libraries(simwing-structure-frame-test
    PRIVATE simwing_structure_frame)

add_executable(simwing-fsi-worker-test
    tests/simwing_fsi_worker_test.cpp
)
target_link_libraries(simwing-fsi-worker-test
    PRIVATE simwing_canonical_case)

add_executable(simwing-periodic-flow-case-test
    tests/simwing_periodic_flow_case_test.cpp
)
target_link_libraries(simwing-periodic-flow-case-test
    PRIVATE simwing_periodic_flow_case)

add_executable(simwing-pressure-jump-frame-test
    tests/simwing_pressure_jump_frame_test.cpp
)
target_link_libraries(simwing-pressure-jump-frame-test
    PRIVATE simwing_fluid_frame)

add_executable(simwing-pressure-jump-case-test
    tests/simwing_pressure_jump_case_test.cpp
)
target_link_libraries(simwing-pressure-jump-case-test
    PRIVATE simwing_pressure_jump_case)

add_executable(simwing-fluid-porous-interface-test
    tests/simwing_fluid_porous_interface_test.cpp
)
target_link_libraries(simwing-fluid-porous-interface-test
    PRIVATE simwing_fluid)

add_executable(simwing-porous-flow-case-test
    tests/simwing_porous_flow_case_test.cpp
)
target_link_libraries(simwing-porous-flow-case-test
    PRIVATE simwing_porous_flow_case)

add_executable(simwing-worker-control-protocol-test
    tests/simwing_worker_control_protocol_test.cpp
)
target_link_libraries(simwing-worker-control-protocol-test
    PRIVATE simwing_worker_control_protocol)

add_executable(simwing-worker-control-stream-test
    tests/simwing_worker_control_stream_test.cpp
)
target_link_libraries(simwing-worker-control-stream-test
    PRIVATE simwing_worker_control_stream)

add_executable(simwing-control-stdio-fixture
    tests/simwing_control_stdio_fixture.cpp
)
target_link_libraries(simwing-control-stdio-fixture PRIVATE
    simwing_open_piston_case
    simwing_periodic_flow_case
    simwing_worker_control_stream
)

add_executable(simwing-periodic-flow-control-test
    tests/simwing_periodic_flow_control_test.cpp
)
target_link_libraries(simwing-periodic-flow-control-test
    PRIVATE simwing_periodic_flow_control)

add_executable(simwing-open-piston-control-test
    tests/simwing_open_piston_control_test.cpp
)
target_link_libraries(simwing-open-piston-control-test
    PRIVATE simwing_open_piston_control)

add_executable(simwing-fluid-projection-test
    tests/simwing_fluid_projection_test.cpp
)
target_link_libraries(simwing-fluid-projection-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-interface-jump-test
    tests/simwing_fluid_interface_jump_test.cpp
)
target_link_libraries(simwing-fluid-interface-jump-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-moving-interface-test
    tests/simwing_fluid_moving_interface_test.cpp
)
target_link_libraries(simwing-fluid-moving-interface-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-control-volume-test
    tests/simwing_fluid_control_volume_test.cpp
)
target_link_libraries(simwing-fluid-control-volume-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-cut-surface-test
    tests/simwing_fluid_cut_surface_test.cpp
)
target_link_libraries(simwing-fluid-cut-surface-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-checkpoint-test
    tests/simwing_fluid_checkpoint_test.cpp
)
target_link_libraries(simwing-fluid-checkpoint-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-diffusion-test
    tests/simwing_fluid_diffusion_test.cpp
)
target_link_libraries(simwing-fluid-diffusion-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-advection-test
    tests/simwing_fluid_advection_test.cpp
)
target_link_libraries(simwing-fluid-advection-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-variable-advection-test
    tests/simwing_fluid_variable_advection_test.cpp
)
target_link_libraries(simwing-fluid-variable-advection-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-projected-advection-test
    tests/simwing_fluid_projected_advection_test.cpp
)
target_link_libraries(simwing-fluid-projected-advection-test PRIVATE simwing_fluid)

add_executable(simwing-fluid-evolution-test
    tests/simwing_fluid_evolution_test.cpp
)
target_link_libraries(simwing-fluid-evolution-test PRIVATE simwing_fluid)

add_executable(simwing-transfer-test
    tests/simwing_transfer_test.cpp
)
target_link_libraries(simwing-transfer-test PRIVATE simwing_transfer)

add_executable(simwing-coupling-test
    tests/simwing_coupling_test.cpp
)
target_link_libraries(simwing-coupling-test PRIVATE simwing_coupling)

add_executable(simwing-fluid-structure-bridge-test
    tests/simwing_fluid_structure_bridge_test.cpp
)
target_link_libraries(simwing-fluid-structure-bridge-test PRIVATE
    simwing_coupling
    simwing_fluid_structure_bridge
)

add_executable(simwing-piston-case-test
    tests/simwing_piston_case_test.cpp
)
target_link_libraries(simwing-piston-case-test PRIVATE
    simwing_piston_case
)

add_executable(simwing-porous-sheet-case-test
    tests/simwing_porous_sheet_case_test.cpp
)
target_link_libraries(simwing-porous-sheet-case-test PRIVATE
    simwing_porous_sheet_case
)

add_executable(simwing-open-piston-case-test
    tests/simwing_open_piston_case_test.cpp
)
target_link_libraries(simwing-open-piston-case-test PRIVATE
    simwing_open_piston_case
)

add_executable(playground-cells-test
    tests/playground_cells_test.cpp
)
target_link_libraries(playground-cells-test PRIVATE playground_sim)

add_executable(playground-material-test
    tests/playground_material_test.cpp
)
target_link_libraries(playground-material-test PRIVATE playground_sim)

add_executable(playground-contact-test
    tests/playground_contact_test.cpp
)
target_link_libraries(playground-contact-test PRIVATE playground_contact)

add_executable(playground-contact-integration-test
    tests/playground_contact_integration_test.cpp
)
target_link_libraries(playground-contact-integration-test PRIVATE playground_sim)

add_executable(flatparts-export-test
    tests/flatparts_export_test.cpp
)
target_link_libraries(flatparts-export-test PRIVATE flatparts Qt6::Gui)

qt_add_executable(simwing-viewer
    tools/simwing_viewer_main.cpp
    src/viewer/viewer_window.cpp
    src/viewer/viewer_window.h
)
target_include_directories(simwing-viewer PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/viewer"
)
target_link_libraries(simwing-viewer PRIVATE
    simwing_viewer_geometry
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::OpenGL
    Qt6::OpenGLWidgets
    OpenGL::GL
)

add_executable(simwing-fsi
    tools/simwing_fsi_main.cpp
)
target_link_libraries(simwing-fsi PRIVATE
    simwing_canonical_case
    simwing_open_piston_control
    simwing_open_piston_case
    simwing_periodic_flow_control
    simwing_periodic_flow_case
    simwing_piston_case
    simwing_porous_flow_case
    simwing_porous_sheet_case
    simwing_pressure_jump_case
    simwing_worker_control_stream
)

qt_add_executable(LEparagliding WIN32
    src/gui/airfoil_file.cpp
    src/gui/airfoil_file.h
    src/gui/airfoil_panel.cpp
    src/gui/airfoil_panel.h
    src/gui/flat_parts_view.cpp
    src/gui/flat_parts_view.h
    src/gui/design_document.cpp
    src/gui/design_document.h
    src/gui/geometry_preprocessor.cpp
    src/gui/geometry_preprocessor.h
    src/gui/geometry_preprocessor_dialog.cpp
    src/gui/geometry_preprocessor_dialog.h
    src/gui/grid_curve_panel.cpp
    src/gui/grid_curve_panel.h
    src/gui/holes_panel.cpp
    src/gui/holes_panel.h
    src/gui/main.cpp
    src/gui/mainwindow.cpp
    src/gui/mainwindow.h
    src/gui/curve_editor.cpp
    src/gui/curve_editor.h
    src/gui/panel_undo.h
    src/gui/paraglider_view.cpp
    src/gui/paraglider_view.h
    src/gui/playground_analysis.cpp
    src/gui/playground_analysis.h
    src/gui/playground_page.cpp
    src/gui/playground_page.h
    src/gui/nest_worker.cpp
    src/gui/nest_worker.h
    src/gui/print_page.cpp
    src/gui/print_page.h
    src/gui/preset_catalog.cpp
    src/gui/preset_catalog.h
    src/gui/section1_curve_panel.cpp
    src/gui/section1_curve_panel.h
    src/gui/section1_curves.cpp
    src/gui/section1_curves.h
    src/gui/section_grid.cpp
    src/gui/section_grid.h
    src/gui/section_grid_panel.cpp
    src/gui/section_grid_panel.h
    src/gui/section_help.cpp
    src/gui/section_help.h
    src/gui/section_specs.cpp
    src/gui/section_specs.h
    src/gui/spline_fit.cpp
    src/gui/spline_fit.h
    resources/manual/manual.qrc
    ${XFLR5_RESOURCES}
)
set_target_properties(LEparagliding PROPERTIES AUTORCC ON AUTOMOC ON)
target_link_libraries(LEparagliding PRIVATE
    Qt6::Widgets Qt6::OpenGLWidgets lep_occt xflr5core softwing_core
    playground_sim flatparts)
add_dependencies(LEparagliding leparagliding-engine)

# Headless timing harness for the Playground's XPBD solve. Developer tool,
# not shipped: it is excluded from the installer and the release archives.
add_executable(softwing-bench
    tools/softwing_bench.cpp
    tools/softwing_gpu.cpp
    tools/softwing_gpu.h
)
target_link_libraries(softwing-bench PRIVATE playground_sim Qt6::OpenGL)

# Wing presets from laboratoridenvol.com: shipped with the executable so the
# Presets menu (and the installer wildcard) picks them up. On macOS they must
# land in Contents/Resources rather than next to the binary in Contents/MacOS:
# data files under MacOS/ make `codesign --deep` treat the folder as a nested
# code bundle and reject the app ("bundle format unrecognized, invalid, or
# unsuitable").
if(APPLE)
    set(lep_presets_deploy_dir
        "$<TARGET_BUNDLE_CONTENT_DIR:LEparagliding>/Resources/presets")
else()
    set(lep_presets_deploy_dir "$<TARGET_FILE_DIR:LEparagliding>/presets")
endif()
add_custom_command(TARGET LEparagliding POST_BUILD
    COMMAND "${CMAKE_COMMAND}" -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/presets"
        "${lep_presets_deploy_dir}"
    COMMENT "Deploying wing presets"
)

set_target_properties(
    LEparagliding
    leparagliding-engine
    reference-report-compare
    preprocessor-reference-compare
    section1-curves-test
    spline-fit-test
    section-grid-test
    airfoil-file-test
    playground-metrics-test
    playground-pressure-solve-test
    playground-cell-air-test
    softwing-material-test
    softwing-cell-volume-test
    softwing-suspension-checkpoint-test
    softwing-checkpoint-persistence-test
    simwing-scene-test
    simwing-model-scene-export-test
    simwing-model-scene-real-export-test
    simwing-structure-test
    simwing-scene-structure-test
    simwing-viewer-protocol-test
    simwing-vector-glyphs-test
    simwing-structure-frame-test
    simwing-fsi-worker-test
    simwing-periodic-flow-case-test
    simwing-pressure-jump-frame-test
    simwing-pressure-jump-case-test
    simwing-fluid-porous-interface-test
    simwing-porous-flow-case-test
    simwing-worker-control-protocol-test
    simwing-worker-control-stream-test
    simwing-control-stdio-fixture
    simwing-periodic-flow-control-test
    simwing-open-piston-control-test
    simwing-fluid-projection-test
    simwing-fluid-interface-jump-test
    simwing-fluid-moving-interface-test
    simwing-fluid-control-volume-test
    simwing-fluid-cut-surface-test
    simwing-fluid-checkpoint-test
    simwing-fluid-diffusion-test
    simwing-fluid-advection-test
    simwing-fluid-variable-advection-test
    simwing-fluid-projected-advection-test
    simwing-fluid-evolution-test
    simwing-transfer-test
    simwing-coupling-test
    simwing-fluid-structure-bridge-test
    simwing-piston-case-test
    simwing-open-piston-case-test
    simwing-viewer
    simwing-fsi
    playground-cells-test
    playground-material-test
    playground-contact-test
    playground-contact-integration-test
    flatparts-export-test
    nesting-bench
    softwing-bench
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
foreach(configuration Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER "${configuration}" configuration_upper)
    set_target_properties(
        LEparagliding
        leparagliding-engine
        reference-report-compare
        preprocessor-reference-compare
        section1-curves-test
        spline-fit-test
        section-grid-test
        airfoil-file-test
        playground-metrics-test
        playground-pressure-solve-test
        playground-cell-air-test
        softwing-material-test
        softwing-cell-volume-test
        softwing-suspension-checkpoint-test
        softwing-checkpoint-persistence-test
        simwing-scene-test
        simwing-model-scene-export-test
        simwing-model-scene-real-export-test
        simwing-structure-test
        simwing-scene-structure-test
        simwing-viewer-protocol-test
        simwing-vector-glyphs-test
        simwing-structure-frame-test
        simwing-fsi-worker-test
        simwing-periodic-flow-case-test
        simwing-pressure-jump-frame-test
        simwing-pressure-jump-case-test
        simwing-fluid-porous-interface-test
        simwing-porous-flow-case-test
        simwing-worker-control-protocol-test
        simwing-worker-control-stream-test
        simwing-control-stdio-fixture
        simwing-periodic-flow-control-test
        simwing-open-piston-control-test
        simwing-fluid-projection-test
        simwing-fluid-interface-jump-test
        simwing-fluid-moving-interface-test
        simwing-fluid-control-volume-test
        simwing-fluid-cut-surface-test
        simwing-fluid-checkpoint-test
        simwing-fluid-diffusion-test
        simwing-fluid-advection-test
        simwing-fluid-variable-advection-test
        simwing-fluid-projected-advection-test
        simwing-fluid-evolution-test
        simwing-transfer-test
        simwing-coupling-test
        simwing-fluid-structure-bridge-test
        simwing-piston-case-test
        simwing-open-piston-case-test
        simwing-viewer
        simwing-fsi
        playground-cells-test
        playground-material-test
        playground-contact-test
        playground-contact-integration-test
        flatparts-export-test
        nesting-bench
        softwing-bench
        PROPERTIES
        "RUNTIME_OUTPUT_DIRECTORY_${configuration_upper}" "${CMAKE_BINARY_DIR}/bin/${configuration}"
    )
endforeach()

if(APPLE)
    # Bundled so macdeployqt can produce a self-contained, DMG-able app.
    set_target_properties(LEparagliding PROPERTIES
        MACOSX_BUNDLE TRUE
        MACOSX_BUNDLE_BUNDLE_NAME "LEparagliding"
        MACOSX_BUNDLE_GUI_IDENTIFIER "com.kulitorum.leparagliding"
        MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
    )
endif()

if(WIN32)
    get_target_property(qmake_executable Qt6::qmake IMPORTED_LOCATION)
    get_filename_component(qt_bin_directory "${qmake_executable}" DIRECTORY)
    find_program(WINDEPLOYQT_EXECUTABLE windeployqt
        HINTS "${qt_bin_directory}"
        NO_CACHE
    )
    if(WINDEPLOYQT_EXECUTABLE)
        foreach(qt_runtime_target LEparagliding simwing-viewer)
            add_custom_command(TARGET ${qt_runtime_target} POST_BUILD
                COMMAND "${WINDEPLOYQT_EXECUTABLE}"
                    --no-translations
                    --no-system-d3d-compiler
                    "$<TARGET_FILE:${qt_runtime_target}>"
                COMMENT "Deploying Qt runtime for ${qt_runtime_target}"
            )
        endforeach()
    endif()

    # Runtime dependencies are broader than link-time imports: STEP and the
    # OpenGL viewer load several OCCT support toolkits transitively.
    set(LEP_OCCT_RUNTIME_TOOLKITS
        ${LEP_OCCT_TOOLKITS}
        TKBO
        TKCAF
        TKCDF
        TKHLR
        TKLCAF
        TKPrim
        TKShHealing
        TKVCAF
        TKXCAF
    )
    list(REMOVE_DUPLICATES LEP_OCCT_RUNTIME_TOOLKITS)
    set(lep_occt_runtime_files)
    foreach(occt_toolkit IN LISTS LEP_OCCT_RUNTIME_TOOLKITS)
        list(APPEND lep_occt_runtime_files
            "${LEP_OCCT_BINARY_DIR}/${occt_toolkit}.dll")
    endforeach()
    # Superset of the third-party runtimes across supported OCCT builds
    # (local build2, CI source build, official prebuilt); each is copied
    # only if the OCCT bin directory actually provides it.
    foreach(occt_dependency
            freetype.dll
            brotlicommon.dll
            brotlidec.dll
            bz2.dll
            libpng16.dll
            zlib1.dll
            FreeImage.dll
            avcodec-57.dll
            avformat-57.dll
            avutil-55.dll
            swscale-4.dll
            openvr_api.dll
            tbb12.dll
            tbbmalloc.dll
            jemalloc.dll)
        if(EXISTS "${LEP_OCCT_BINARY_DIR}/${occt_dependency}")
            list(APPEND lep_occt_runtime_files
                "${LEP_OCCT_BINARY_DIR}/${occt_dependency}")
        endif()
    endforeach()
    foreach(runtime_target LEparagliding leparagliding-engine)
        add_custom_command(TARGET ${runtime_target} POST_BUILD
            COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                ${lep_occt_runtime_files}
                "$<TARGET_FILE_DIR:${runtime_target}>"
            COMMENT "Deploying Open CASCADE runtime for ${runtime_target}"
            COMMAND_EXPAND_LISTS
        )
    endforeach()
endif()

enable_testing()
add_test(NAME engine_help COMMAND leparagliding-engine --help)
add_test(NAME section1_curves COMMAND section1-curves-test)
add_test(NAME spline_fit COMMAND spline-fit-test)
add_test(NAME section_grid COMMAND section-grid-test)
add_test(NAME airfoil_file COMMAND airfoil-file-test
    "${CMAKE_CURRENT_SOURCE_DIR}/Data/gnua.txt")
add_test(NAME playground_metrics COMMAND playground-metrics-test)
add_test(NAME playground_pressure_solve COMMAND playground-pressure-solve-test)
add_test(NAME playground_cell_air COMMAND playground-cell-air-test)
add_test(NAME softwing_material COMMAND softwing-material-test)
add_test(NAME softwing_cell_volume COMMAND softwing-cell-volume-test)
add_test(NAME softwing_suspension_checkpoint
    COMMAND softwing-suspension-checkpoint-test)
add_test(NAME softwing_checkpoint_persistence
    COMMAND softwing-checkpoint-persistence-test)
add_test(NAME simwing_scene COMMAND simwing-scene-test)
add_test(NAME simwing_model_scene_export
    COMMAND simwing-model-scene-export-test)
add_test(NAME simwing_model_scene_real_export
    COMMAND simwing-model-scene-real-export-test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/3.28/leparagliding.txt"
        "${CMAKE_BINARY_DIR}/simwing-model-scene-real-export-output")
add_test(NAME simwing_structure COMMAND simwing-structure-test)
add_test(NAME simwing_scene_structure COMMAND simwing-scene-structure-test)
add_test(NAME simwing_viewer_protocol COMMAND simwing-viewer-protocol-test)
add_test(NAME simwing_vector_glyphs COMMAND simwing-vector-glyphs-test)
add_test(NAME simwing_structure_frame COMMAND simwing-structure-frame-test)
add_test(NAME simwing_fsi_worker COMMAND simwing-fsi-worker-test)
add_test(NAME simwing_periodic_flow_case
    COMMAND simwing-periodic-flow-case-test)
add_test(NAME simwing_pressure_jump_frame
    COMMAND simwing-pressure-jump-frame-test)
add_test(NAME simwing_pressure_jump_case
    COMMAND simwing-pressure-jump-case-test)
add_test(NAME simwing_fluid_porous_interface
    COMMAND simwing-fluid-porous-interface-test)
add_test(NAME simwing_porous_flow_case
    COMMAND simwing-porous-flow-case-test)
add_test(NAME simwing_worker_control_protocol
    COMMAND simwing-worker-control-protocol-test)
add_test(NAME simwing_worker_control_stream
    COMMAND simwing-worker-control-stream-test)
add_test(NAME simwing_periodic_flow_control
    COMMAND simwing-periodic-flow-control-test)
add_test(NAME simwing_open_piston_control
    COMMAND simwing-open-piston-control-test)
add_test(NAME simwing_fluid_projection COMMAND simwing-fluid-projection-test)
add_test(NAME simwing_fluid_interface_jump
    COMMAND simwing-fluid-interface-jump-test)
add_test(NAME simwing_fluid_moving_interface
    COMMAND simwing-fluid-moving-interface-test)
add_test(NAME simwing_fluid_control_volume
    COMMAND simwing-fluid-control-volume-test)
add_test(NAME simwing_fluid_cut_surface
    COMMAND simwing-fluid-cut-surface-test)
add_test(NAME simwing_fluid_checkpoint
    COMMAND simwing-fluid-checkpoint-test)
add_test(NAME simwing_fluid_diffusion
    COMMAND simwing-fluid-diffusion-test)
add_test(NAME simwing_fluid_advection
    COMMAND simwing-fluid-advection-test)
add_test(NAME simwing_fluid_variable_advection
    COMMAND simwing-fluid-variable-advection-test)
add_test(NAME simwing_fluid_projected_advection
    COMMAND simwing-fluid-projected-advection-test)
add_test(NAME simwing_fluid_evolution
    COMMAND simwing-fluid-evolution-test)
add_test(NAME simwing_transfer COMMAND simwing-transfer-test)
add_test(NAME simwing_coupling COMMAND simwing-coupling-test)
add_test(NAME simwing_fluid_structure_bridge
    COMMAND simwing-fluid-structure-bridge-test)
add_test(NAME simwing_piston_case COMMAND simwing-piston-case-test)
add_test(NAME simwing_porous_sheet_case
    COMMAND simwing-porous-sheet-case-test)
add_test(NAME simwing_open_piston_case COMMAND simwing-open-piston-case-test)
add_test(NAME simwing_fsi_help COMMAND simwing-fsi --help)
add_test(NAME simwing_fsi_piston_headless
    COMMAND simwing-fsi
        --case piston
        --steps 4
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-piston-headless.swtrace")
add_test(NAME simwing_fsi_porous_sheet_headless
    COMMAND simwing-fsi
        --case porous-sheet
        --steps 120
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-porous-sheet-headless.swtrace")
set_tests_properties(simwing_fsi_porous_sheet_headless PROPERTIES
    PASS_REGULAR_EXPRESSION "completed 120 step\\(s\\)")
add_test(NAME simwing_fsi_periodic_flow_headless
    COMMAND simwing-fsi
        --case periodic-flow
        --steps 4
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-periodic-flow-headless.swtrace")
add_test(NAME simwing_fsi_pressure_jump_headless
    COMMAND simwing-fsi
        --case pressure-jump
        --steps 4
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-pressure-jump-headless.swtrace")
set_tests_properties(simwing_fsi_pressure_jump_headless PROPERTIES
    PASS_REGULAR_EXPRESSION "completed 4 pressure-jump step\\(s\\)")
add_test(NAME simwing_fsi_porous_flow_headless
    COMMAND simwing-fsi
        --case porous-flow
        --steps 120
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-porous-flow-headless.swtrace")
set_tests_properties(simwing_fsi_porous_flow_headless PROPERTIES
    PASS_REGULAR_EXPRESSION "completed 120 porous-flow step\\(s\\)")
add_test(NAME simwing_fsi_control_stdio
    COMMAND ${CMAKE_COMMAND}
        "-DSIMWING_FSI=$<TARGET_FILE:simwing-fsi>"
        "-DCONTROL_FIXTURE=$<TARGET_FILE:simwing-control-stdio-fixture>"
        "-DWORKER_CASE=periodic-flow"
        "-DVERIFY_MODE=verify"
        "-DRESUME_VERIFY_MODE=verify-resume"
        "-DCOMMAND_FILE=${CMAKE_BINARY_DIR}/simwing-control-commands.swcc"
        "-DRESPONSE_FILE=${CMAKE_BINARY_DIR}/simwing-control-responses.swcr"
        "-DCHECKPOINT_FILE=${CMAKE_BINARY_DIR}/simwing-control-checkpoint.swpc"
        "-DTRACE_FILE=${CMAKE_BINARY_DIR}/simwing-control.swtrace"
        "-DERROR_FILE=${CMAKE_BINARY_DIR}/simwing-control-stderr.txt"
        "-DRESUME_COMMAND_FILE=${CMAKE_BINARY_DIR}/simwing-control-resume-commands.swcc"
        "-DRESUME_RESPONSE_FILE=${CMAKE_BINARY_DIR}/simwing-control-resume-responses.swcr"
        "-DRESUME_TRACE_FILE=${CMAKE_BINARY_DIR}/simwing-control-resume.swtrace"
        "-DRESUME_ERROR_FILE=${CMAKE_BINARY_DIR}/simwing-control-resume-stderr.txt"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSimWingControlStdioTest.cmake")
set_tests_properties(simwing_fsi_control_stdio PROPERTIES TIMEOUT 45)
add_test(NAME simwing_fsi_open_piston_control_stdio
    COMMAND ${CMAKE_COMMAND}
        "-DSIMWING_FSI=$<TARGET_FILE:simwing-fsi>"
        "-DCONTROL_FIXTURE=$<TARGET_FILE:simwing-control-stdio-fixture>"
        "-DWORKER_CASE=open-piston"
        "-DVERIFY_MODE=verify-open"
        "-DRESUME_VERIFY_MODE=verify-open-resume"
        "-DCOMMAND_FILE=${CMAKE_BINARY_DIR}/simwing-open-control-commands.swcc"
        "-DRESPONSE_FILE=${CMAKE_BINARY_DIR}/simwing-open-control-responses.swcr"
        "-DCHECKPOINT_FILE=${CMAKE_BINARY_DIR}/simwing-open-control-checkpoint.swop"
        "-DTRACE_FILE=${CMAKE_BINARY_DIR}/simwing-open-control.swtrace"
        "-DERROR_FILE=${CMAKE_BINARY_DIR}/simwing-open-control-stderr.txt"
        "-DRESUME_COMMAND_FILE=${CMAKE_BINARY_DIR}/simwing-open-control-resume-commands.swcc"
        "-DRESUME_RESPONSE_FILE=${CMAKE_BINARY_DIR}/simwing-open-control-resume-responses.swcr"
        "-DRESUME_TRACE_FILE=${CMAKE_BINARY_DIR}/simwing-open-control-resume.swtrace"
        "-DRESUME_ERROR_FILE=${CMAKE_BINARY_DIR}/simwing-open-control-resume-stderr.txt"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSimWingControlStdioTest.cmake")
set_tests_properties(simwing_fsi_open_piston_control_stdio PROPERTIES
    TIMEOUT 45)
add_test(NAME simwing_fsi_control_rejects_unsupported_case
    COMMAND simwing-fsi --case structural --control-stdio --no-viewer)
set_tests_properties(simwing_fsi_control_rejects_unsupported_case PROPERTIES
    WILL_FAIL TRUE)
add_test(NAME simwing_fsi_periodic_flow_checkpoint_write
    COMMAND simwing-fsi
        --case periodic-flow
        --steps 4
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-periodic-flow-checkpoint-write.swtrace"
        --checkpoint-out "${CMAKE_BINARY_DIR}/simwing-periodic-flow.swpc"
        --checkpoint-every 2)
set_tests_properties(simwing_fsi_periodic_flow_checkpoint_write PROPERTIES
    PASS_REGULAR_EXPRESSION "completed 4 periodic-flow step\\(s\\).*checkpoint-writes=2")
add_test(NAME simwing_fsi_periodic_flow_checkpoint_resume
    COMMAND simwing-fsi
        --case periodic-flow
        --steps 3
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-periodic-flow-checkpoint-resume.swtrace"
        --checkpoint-in "${CMAKE_BINARY_DIR}/simwing-periodic-flow.swpc"
        --checkpoint-out "${CMAKE_BINARY_DIR}/simwing-periodic-flow.swpc"
        --checkpoint-every 2)
set_tests_properties(simwing_fsi_periodic_flow_checkpoint_resume PROPERTIES
    DEPENDS simwing_fsi_periodic_flow_checkpoint_write
    PASS_REGULAR_EXPRESSION "completed 7 periodic-flow step\\(s\\).*checkpoint-writes=2")
add_test(NAME simwing_fsi_periodic_flow_checkpoint_verify
    COMMAND simwing-fsi
        --case periodic-flow
        --steps 0
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-periodic-flow-checkpoint-verify.swtrace"
        --checkpoint-in "${CMAKE_BINARY_DIR}/simwing-periodic-flow.swpc")
set_tests_properties(simwing_fsi_periodic_flow_checkpoint_verify PROPERTIES
    DEPENDS simwing_fsi_periodic_flow_checkpoint_resume
    PASS_REGULAR_EXPRESSION "completed 7 periodic-flow step\\(s\\)")
add_test(NAME simwing_fsi_checkpoint_rejects_foreign_case
    COMMAND simwing-fsi
        --case structural
        --steps 0
        --no-viewer
        --checkpoint-out "${CMAKE_BINARY_DIR}/invalid-structural.swpc")
set_tests_properties(simwing_fsi_checkpoint_rejects_foreign_case PROPERTIES
    WILL_FAIL TRUE)
add_test(NAME simwing_fsi_checkpoint_interval_requires_output
    COMMAND simwing-fsi
        --case periodic-flow
        --steps 0
        --no-viewer
        --checkpoint-every 2)
set_tests_properties(simwing_fsi_checkpoint_interval_requires_output PROPERTIES
    WILL_FAIL TRUE)
add_test(NAME simwing_fsi_open_piston_headless
    COMMAND simwing-fsi
        --case open-piston
        --steps 4
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-open-piston-headless.swtrace")
add_test(NAME simwing_fsi_open_piston_rebase_headless
    COMMAND simwing-fsi
        --case open-piston
        --steps 2401
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-open-piston-rebase-headless.swtrace")
add_test(NAME simwing_fsi_open_piston_checkpoint_write
    COMMAND simwing-fsi
        --case open-piston
        --steps 1200
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-open-piston-checkpoint-write.swtrace"
        --checkpoint-out "${CMAKE_BINARY_DIR}/simwing-open-piston.swop"
        --checkpoint-every 600)
set_tests_properties(simwing_fsi_open_piston_checkpoint_write PROPERTIES
    PASS_REGULAR_EXPRESSION "completed 1200 open-piston step\\(s\\).*checkpoint-writes=2")
add_test(NAME simwing_fsi_open_piston_checkpoint_resume
    COMMAND simwing-fsi
        --case open-piston
        --steps 3
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-open-piston-checkpoint-resume.swtrace"
        --checkpoint-in "${CMAKE_BINARY_DIR}/simwing-open-piston.swop"
        --checkpoint-out "${CMAKE_BINARY_DIR}/simwing-open-piston.swop"
        --checkpoint-every 2)
set_tests_properties(simwing_fsi_open_piston_checkpoint_resume PROPERTIES
    DEPENDS simwing_fsi_open_piston_checkpoint_write
    PASS_REGULAR_EXPRESSION "completed 1203 open-piston step\\(s\\).*checkpoint-writes=2")
add_test(NAME simwing_fsi_open_piston_checkpoint_verify
    COMMAND simwing-fsi
        --case open-piston
        --steps 0
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-open-piston-checkpoint-verify.swtrace"
        --checkpoint-in "${CMAKE_BINARY_DIR}/simwing-open-piston.swop")
set_tests_properties(simwing_fsi_open_piston_checkpoint_verify PROPERTIES
    DEPENDS simwing_fsi_open_piston_checkpoint_resume
    PASS_REGULAR_EXPRESSION "completed 1203 open-piston step\\(s\\)")
add_test(NAME simwing_fsi_open_piston_rejects_foreign_checkpoint
    COMMAND simwing-fsi
        --case open-piston
        --steps 0
        --no-viewer
        --trace "${CMAKE_BINARY_DIR}/simwing-open-piston-foreign.swtrace"
        --checkpoint-in "${CMAKE_BINARY_DIR}/simwing-periodic-flow.swpc")
set_tests_properties(simwing_fsi_open_piston_rejects_foreign_checkpoint
    PROPERTIES
        DEPENDS simwing_fsi_periodic_flow_checkpoint_verify
        WILL_FAIL TRUE)
add_test(NAME simwing_viewer_help COMMAND simwing-viewer --help)
if(NOT WIN32)
    set_tests_properties(simwing_viewer_help PROPERTIES
        ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
    )
endif()
add_test(NAME playground_cells COMMAND playground-cells-test)
add_test(NAME playground_material COMMAND playground-material-test)
add_test(NAME playground_contact COMMAND playground-contact-test)
add_test(NAME playground_contact_integration
    COMMAND playground-contact-integration-test)
add_test(NAME flatparts_export COMMAND flatparts-export-test)
add_test(NAME engine_rejects_missing_input
    COMMAND leparagliding-engine
        "${CMAKE_BINARY_DIR}/does-not-exist.txt"
        "${CMAKE_BINARY_DIR}/missing-input-output"
)
set_tests_properties(engine_rejects_missing_input PROPERTIES WILL_FAIL TRUE)
add_test(NAME gui_smoke COMMAND LEparagliding --smoke-test)
set_tests_properties(gui_smoke PROPERTIES
    TIMEOUT 20
)
add_test(NAME studio_model_smoke
    COMMAND LEparagliding
        --studio-self-test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/3.28/leparagliding.txt"
)
set_tests_properties(studio_model_smoke PROPERTIES
    TIMEOUT 180
)
add_test(NAME gui_headless_help COMMAND LEparagliding --headless --help)
# One representative wing preset must calculate cleanly end to end; catches
# regressions in the packaged resources/presets tree.
add_test(NAME preset_gnua7_calculates
    COMMAND leparagliding-engine
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/presets/gnuA7/24/leparagliding.txt"
        "${CMAKE_BINARY_DIR}/preset-test-output"
)
set_tests_properties(preset_gnua7_calculates PROPERTIES
    TIMEOUT 120
)
# Every shipped preset must pass the Studio editor's stricter validation
# (unique section headers, no blank records) — the engine alone tolerates
# upstream quirks that would leave the GUI preview silently unavailable.
add_test(NAME presets_validate_in_studio
    COMMAND LEparagliding
        --validate-presets
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/presets"
)
set_tests_properties(presets_validate_in_studio PROPERTIES
    TIMEOUT 60
)
add_test(NAME engine_tolerates_blank_lines
    COMMAND leparagliding-engine
        --resource-dir "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/3.28"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/3.28/leparagliding-blank-lines.txt"
        "${CMAKE_BINARY_DIR}/blank-lines-test-output"
)
set_tests_properties(engine_tolerates_blank_lines PROPERTIES
    TIMEOUT 90
)
foreach(preprocessor_case
        basic
        vault1-cells1-odd
        vault2-cells2-even
        vault2-cells4)
    add_test(NAME "preprocessor_matches_fortran_${preprocessor_case}"
        COMMAND preprocessor-reference-compare
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/pre-1.6/${preprocessor_case}/pre-data.txt"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/pre-1.6/${preprocessor_case}/geometry-out.txt"
    )
endforeach()
set_tests_properties(gui_headless_help PROPERTIES
    TIMEOUT 20
)

set(reference_input
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/3.28/leparagliding.txt")
set(reference_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/reference/3.28")
if(WIN32
   AND EXISTS "${reference_input}"
   AND EXISTS "${reference_directory}/lep-out.txt"
   AND EXISTS "${reference_directory}/lines.txt"
   AND EXISTS "${reference_directory}/run-log.txt")
    add_test(NAME engine_matches_fortran_reference
        COMMAND "${CMAKE_COMMAND}"
            "-DENGINE=$<TARGET_FILE:leparagliding-engine>"
            "-DREPORT_COMPARATOR=$<TARGET_FILE:reference-report-compare>"
            "-DINPUT=${reference_input}"
            "-DOUTPUT_DIR=${CMAKE_BINARY_DIR}/reference-test-output"
            "-DREFERENCE_DIR=${reference_directory}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompareReferenceRun.cmake"
    )
    set_tests_properties(engine_matches_fortran_reference PROPERTIES
        TIMEOUT 90
    )
endif()

install(TARGETS LEparagliding leparagliding-engine simwing-fsi simwing-viewer
    RUNTIME DESTINATION .
    BUNDLE DESTINATION .
)
install(DIRECTORY resources/presets
    DESTINATION .
)
if(WIN32)
    install(FILES ${lep_occt_runtime_files}
        DESTINATION .
    )
endif()

# The Inno Setup script includes installer/version.iss, regenerated on every
# configure so the installer filename embeds the HEAD commit. Commit before
# building a release installer or the hash will be stale.
execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE LEP_GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT LEP_GIT_HASH)
    set(LEP_GIT_HASH "unknown")
endif()
set(LEP_INSTALLER_VERSION "${PROJECT_VERSION}-${LEP_GIT_HASH}")
configure_file(
    installer/version.iss.in
    "${CMAKE_CURRENT_SOURCE_DIR}/installer/version.iss"
    @ONLY
)

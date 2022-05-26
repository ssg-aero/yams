cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

project(gbs LANGUAGES CXX)

#
# c++ conf
#
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# #

include(GNUInstallDirs)

message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> binaries: " ${CMAKE_INSTALL_BINDIR})
message("  -> libs: " ${CMAKE_INSTALL_LIBDIR})
message("  -> includes: " ${CMAKE_INSTALL_INCLUDEDIR})
#
# set output dirs (vs is a pain)
#
function(SET_OUT_DIR REQUIRED_ARG)
    list(GET ARGV 0 TARGET_NAME)
    list(GET ARGV 1 OUT_DIR)
    message(${TARGET_NAME})
    message(${OUT_DIR})
    foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
        string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
        set_property(TARGET ${TARGET_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
        set_property(TARGET ${TARGET_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
        set_property(TARGET ${TARGET_NAME} PROPERTY ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
    endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
endfunction()

# Paths
file(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} lib_dir)
file(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} bin_dir)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/gbs)
link_directories(${lib_dir})
link_directories(${bin_dir})


# Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
# libIGES
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/iges/)
#
# xtensors
#
set(XTENSOR_USE_XSIMD BOOL:ON)
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xsimd REQUIRED)
add_definitions(-DXTENSOR_USE_XSIMD)
#
# VTK
#
find_package(
    VTK REQUIRED
    COMPONENTS
        CommonCore
        CommonColor
        CommonDataModel
        FiltersSources
        InteractionStyle
        RenderingCore
        RenderingAnnotation
        RenderingFreeType
        IOXML
        ChartsCore
        InteractionStyle
        RenderingContext2D
        RenderingContextOpenGL2
        RenderingCore
        RenderingFreeType
        ViewsContext2D
        RenderingOpenGL2 
        Python
        WrappingPythonCore
)
# NLopt
find_package(NLopt REQUIRED)

if(LINUX)
    # add_definitions(-D__TBB_show_deprecation_message_task_H))
    find_package(TBB REQUIRED)
endif()

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif()

if(${BUILD_PYTHON_BINDING})
    add_subdirectory(python)
endif()

add_executable(
    yams MACOSX_BUNDLE 
    yams.cpp
)
target_link_libraries(
    yams 
    NLopt::nlopt
    ${VTK_LIBRARIES}
    xtensor 
    xtensor::optimize 
    xtensor::use_xsimd
    gbs-render
    iges
    TBB::tbb
)

vtk_module_autoinit(
    TARGETS yams
    MODULES ${VTK_LIBRARIES}
)
install(TARGETS yams)
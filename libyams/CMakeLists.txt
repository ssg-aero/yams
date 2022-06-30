project(libyams)

# C++ conf
set(CMAKE_CXX_STANDARD 20)

# dependencies
## Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#include_directories(${EIGEN3_INCLUDE_DIR})

## libIGES
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/iges/)

## xtensor/xsimd
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xsimd REQUIRED)
add_definitions(-DXTENSOR_USE_XSIMD)

## VTK
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

## NLopt
find_package(NLopt REQUIRED)

# TBB
if(LINUX)
    # add_definitions(-D__TBB_show_deprecation_message_task_H))
    find_package(TBB REQUIRED)
endif()

set(YAMS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(yams_headers
    ${YAMS_INCLUDE_DIR}/yams/correlationscax.h
    ${YAMS_INCLUDE_DIR}/yams/curvaturesolver.h
    ${YAMS_INCLUDE_DIR}/yams/datastorage.h
    ${YAMS_INCLUDE_DIR}/yams/diffop.h
    ${YAMS_INCLUDE_DIR}/yams/eqcurvaturesolver.h
    ${YAMS_INCLUDE_DIR}/yams/gridmetrics.h
    ${YAMS_INCLUDE_DIR}/yams/gridreader.h
    ${YAMS_INCLUDE_DIR}/yams/gridrender.h
    ${YAMS_INCLUDE_DIR}/yams/meridionalsolvercase.h
    ${YAMS_INCLUDE_DIR}/yams/meshtools.h
    ${YAMS_INCLUDE_DIR}/yams/plots.h
    ${YAMS_INCLUDE_DIR}/yams/solvercasetools.h
)

add_library(libyams INTERFACE ${yams_headers})
target_compile_features(libyams INTERFACE cxx_std_20)

target_include_directories(
    libyams
    INTERFACE
        $<BUILD_INTERFACE:${YAMS_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(
    libyams INTERFACE
    NLopt::nlopt
    ${VTK_LIBRARIES}
    xtensor 
    xtensor::optimize 
    xtensor::use_xsimd
    Eigen3::Eigen
    gbs-render
    iges
    TBB::tbb
)

vtk_module_autoinit(
    TARGETS libyams
    MODULES ${VTK_LIBRARIES}
)

install(TARGETS libyams EXPORT ${PROJECT_NAME}-targets)
install(FILES ${yams_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/yams)

# CMake config
# Add helpers for cmake configuration
include(CMakePackageConfigHelpers)

# Definition
set(ConfigFileInstallDir ${CMAKECONFIG_INSTALL_DIR})
configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                            "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                            INSTALL_DESTINATION ${ConfigFileInstallDir})

# Makes the project importable from the build directory
export(EXPORT ${PROJECT_NAME}-targets
FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")

# Installation
install(EXPORT ${PROJECT_NAME}-targets 
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${ConfigFileInstallDir})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        DESTINATION ${ConfigFileInstallDir})

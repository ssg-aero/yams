cmake_minimum_required(VERSION 3.17.0)

project(pyams LANGUAGES CXX)

# set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Python 
    COMPONENTS 
        Interpreter 
        Development
)
message("Python3_SITELIB:" ${Python3_SITELIB})
find_package(pybind11 REQUIRED)

file(GLOB SRC_LIST "*.cpp")

pybind11_add_module(
    ${PROJECT_NAME}
    yamsbind.cpp
)

target_link_libraries(
        ${PROJECT_NAME} PRIVATE 
        ${VTK_LIBRARIES}
        TBB::tbb 
        gbs-render
        Eigen3::Eigen
        libIGES::iges
)

vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
)

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME "yams"
)

# The module keeps the compiler's precise floating point (no /fp:fast, no
# -ffast-math, which also made yams.mesh_channel throw an uncaught exception
# with clang on Windows): no reassociation, and the robust geometric
# predicates of gbs intact. The results are closer between OS, not identical:
# the compiler's default contraction (FMA on arm64) and the libm still differ.
# It is built for the packager's target CPU, not the build machine's:
# -march=native would make a conda package that dies on an illegal
# instruction on an older CPU.
# For a local build on this machine only. MSVC has no "native" target: the
# option asks for AVX2 there, which the CPU of this machine must support.
option(YAMS_PYTHON_NATIVE "Build the Python module for this machine's CPU (-march=native; MSVC: /arch:AVX2, the CPU must support it): not for packages" OFF)
if(YAMS_PYTHON_NATIVE)
    if(MSVC)
        message(STATUS "YAMS_PYTHON_NATIVE: MSVC has no native target, using /arch:AVX2")
        target_compile_options(${PROJECT_NAME} PRIVATE /arch:AVX2)
    else()
        target_compile_options(${PROJECT_NAME} PRIVATE -march=native)
    endif()
endif()


# install(
#     TARGETS ${PROJECT_NAME} 
#     DESTINATION ${Python3_SITELIB}/${PROJECT_NAME}
# )

file(GLOB PY_SRC_LIST "*.py")

if($ENV{CONDA_BUILD})
        install(TARGETS ${PROJECT_NAME} 
        DESTINATION $ENV{SP_DIR}/${PROJECT_NAME}
        )
        install(FILES 
                        ${PY_SRC_LIST}
                DESTINATION 
                        $ENV{SP_DIR}/${PROJECT_NAME}
                )
else($ENV{CONDA_BUILD})
        install(TARGETS ${PROJECT_NAME} 
                DESTINATION ${Python3_SITELIB}/${PROJECT_NAME}
        )
        install(FILES 
                        ${PY_SRC_LIST}
                DESTINATION 
                        ${Python3_SITELIB}/${PROJECT_NAME}
        )
endif()
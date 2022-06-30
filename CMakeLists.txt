cmake_minimum_required(VERSION 3.17.0)

project(yams)

# Build options
option(BUILD_LIB "Enable the build of header-only library" ON)
option(BUILD_CLI "Enable the build of the CLI" OFF)
option(BUILD_TESTS "Enable the build of the C++ tests" OFF)
option(BUILD_PYTHON_BINDINGS "Enable the build of Python bindings" OFF)

# Set installation directories
include(GNUInstallDirs)

message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> binaries: " ${CMAKE_INSTALL_BINDIR})
message("  -> libs: " ${CMAKE_INSTALL_LIBDIR})
message("  -> includes: " ${CMAKE_INSTALL_INCLUDEDIR})

set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for libyamsConfig.cmake")

# Paths
# TODO: remove this after GBS refactoring
file(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} lib_dir)
file(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} bin_dir)
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/gbs)
link_directories(${lib_dir})
link_directories(${bin_dir})

if (${BUILD_LIB})
    add_subdirectory(libyams)
endif()

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif()

if(${BUILD_PYTHON_BINDINGS})
    add_subdirectory(python)
endif()

if(${BUILD_CLI})
    add_subdirectory(yams-cli)
endif()

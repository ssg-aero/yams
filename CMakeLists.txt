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
# # CODAN ENV
# #
# if(EXISTS $ENV{LIBRARY_PREFIX})
#     file(TO_CMAKE_PATH $ENV{LIBRARY_BIN} LIBRARY_BIN)
#     file(TO_CMAKE_PATH $ENV{LIBRARY_INC} LIBRARY_INC)
#     file(TO_CMAKE_PATH $ENV{LIBRARY_LIB} LIBRARY_LIB)
# elseif(EXISTS $ENV{CONDA_PREFIX})
#     file(TO_CMAKE_PATH $ENV{CONDA_PREFIX} CONDA_PREFIX)
#     file(TO_CMAKE_PATH ${CONDA_PREFIX}/Library/bin LIBRARY_BIN)
#     file(TO_CMAKE_PATH ${CONDA_PREFIX}/Library/include LIBRARY_INC)
#     file(TO_CMAKE_PATH ${CONDA_PREFIX}/Library/lib LIBRARY_LIB)
# else()
#     message(FATAL_ERROR "A conda environement is required")
# endif()
# message(${LIBRARY_LIB} ${LIBRARY_BIN} ${LIBRARY_INC})
# c++ conf
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_compile_options(/bigobj)

# Set installation directories (CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_BINDIR, CMAKE_INSTALL_LIBDIR)
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

# include_directories(${LIBRARY_INC})
# link_directories(${LIBRARY_BIN})
# link_directories(${LIBRARY_LIB})
# Paths
file(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX}/lib lib_dir)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/gbs)
link_directories(${lib_dir})


# gbs
include_directories(C:/Users/sebastien/workspace/gbslib)
# Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
# libIGES
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/iges/)

add_subdirectory(tests)
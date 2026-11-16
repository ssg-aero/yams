cmake_minimum_required(VERSION 3.17.0)

project(yams_tests LANGUAGES CXX)

find_package(doctest REQUIRED)

option(YAMS_TESTS_PLOT "Open the interactive VTK plots of the C++ tests (blocks until each window is closed)" OFF)

if(MSVC)
    add_compile_options(/bigobj)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-nan-infinity-disabled ")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(VECTORIZATION_FLAGS /arch:AVX2 /fp:fast /Qvec-report:1)
endif()
# Tests built with VECTORIZATION_FLAGS (MSVC only)
set(VECTORIZED_TESTS tests_b2bsolver tests_vectorization tests_euler2D tests_biflux)

# One executable per tests_*.cpp; doctest provides main().
file(GLOB SRC_LIST "tests_*.cpp")
foreach(file_name ${SRC_LIST})
    get_filename_component(exe_name ${file_name} NAME_WE)
    add_executable(${exe_name} ${file_name})

    target_link_libraries(${exe_name}
        doctest::doctest
        ${VTK_LIBRARIES}
        gbs-render
        iges
        TBB::tbb
        Eigen3::Eigen
        libIGES::iges
    )

    target_compile_definitions(${exe_name} PRIVATE
        DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
        TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
        YAMS_TESTS_PLOT=$<BOOL:${YAMS_TESTS_PLOT}>
    )

    if(MSVC)
        SET_OUT_DIR(${exe_name} ${CMAKE_BINARY_DIR})
    endif()
    vtk_module_autoinit(
        TARGETS ${exe_name}
        MODULES ${VTK_LIBRARIES}
    )

    if(exe_name IN_LIST VECTORIZED_TESTS)
        target_compile_options(${exe_name} PRIVATE ${VECTORIZATION_FLAGS})
    endif()

    add_test(NAME ${exe_name} COMMAND ${exe_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    install(TARGETS ${exe_name})
endforeach()

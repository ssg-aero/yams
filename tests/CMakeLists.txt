cmake_minimum_required(VERSION 3.17.0)

find_package(GTest REQUIRED)

project(quiss_tests LANGUAGES CXX)
if(MSVC)
    add_compile_options(/bigobj)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-nan-infinity-disabled ")
endif()

message(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" )
    set(VECTORIZATION_FLAGS 
        /arch:AVX2 
        # /O2
        # /Qpar
        # /openmp
        /fp:fast
        # /fp:precise /fp:contract
        /Qvec-report:1
        # /Qpar-report:2
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
    set(VECTORIZATION_FLAGS 
        # -Rpass=loop-vectorize
        # -Wno-nan-infinity-disabled
        # -stdlib=libstdc++
    )
endif()

message(STATUS "${VTK_LIBRARIES} " ${VTK_LIBRARIES})

file(GLOB SRC_LIST "tests_*")
foreach(file_name ${SRC_LIST})
        get_filename_component(exe_name ${file_name} NAME_WE  )
        add_executable(${exe_name} ${file_name})

        target_link_libraries(${exe_name} 
                            GTest::gtest
                            GTest::gtest_main
                            # xtensor 
                            # xtensor::optimize 
                            # xtensor::use_xsimd
                            ${VTK_LIBRARIES}
                            NLopt::nlopt
                            gbs-render
                            iges
                            TBB::tbb
        )

        install(TARGETS ${exe_name})
        if ( MSVC )
        SET_OUT_DIR(${exe_name} ${CMAKE_BINARY_DIR})
        endif ( MSVC )
        vtk_module_autoinit(
                TARGETS ${exe_name}
                MODULES ${VTK_LIBRARIES}
                )
        
        if(exe_name STREQUAL "tests_b2bsolver")
            target_compile_options(${exe_name} PRIVATE
                # /arch:AVX2 
                # # /O2
                # # /Qpar
                # # /openmp
                # /fp:fast
                # # /fp:precise /fp:contract
                # /Qvec-report:1
                # # /Qpar-report:2
                ${VECTORIZATION_FLAGS}
            )
        endif()

        if(exe_name STREQUAL "tests_vectorization")
            target_compile_options(${exe_name} PRIVATE
                # /arch:AVX2
                # # /O2
                # # /Qpar
                # # /openmp
                # /fp:fast
                # # /fp:precise /fp:contract
                # # /Qvec-report:2
                # # /Qpar-report:2
                ${VECTORIZATION_FLAGS}
            )
        endif()

        if(exe_name STREQUAL "tests_euler2D")
            target_compile_options(${exe_name} PRIVATE
                # /arch:AVX2 
                # # /O2
                # # /Qpar
                # /openmp
                # /fp:fast
                # # /fp:precise /fp:contract
                # /Qvec-report:1
                # # /Qpar-report:1
                ${VECTORIZATION_FLAGS}
            )
        endif()

        if(exe_name STREQUAL "tests_biflux")
            target_compile_options(${exe_name} PRIVATE
                # /arch:AVX2 
                # # /O2
                # # /Qpar
                # /openmp
                # /fp:fast
                # # /fp:precise /fp:contract
                # # /Qvec-report:2
                # # /Qpar-report:1
                ${VECTORIZATION_FLAGS}
            )
        endif()
        
endforeach()
cmake_minimum_required(VERSION 3.17.0)

find_package(GTest REQUIRED)

project(quiss_tests LANGUAGES CXX)
if(MSVC)
    add_compile_options(/bigobj)
endif()

message(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

file(GLOB SRC_LIST "tests_*")
foreach(file_name ${SRC_LIST})
        get_filename_component(exe_name ${file_name} NAME_WE  )
        add_executable(${exe_name} ${file_name})

        target_link_libraries(${exe_name} 
                            GTest::gtest
                            GTest::gtest_main
                            xtensor 
                            xtensor::optimize 
                            xtensor::use_xsimd
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
                /arch:AVX2 
                # /O2
                # /Qpar
                # /openmp
                /fp:fast
                # /fp:precise /fp:contract
                # /Qvec-report:2
                # /Qpar-report:2
            )
        endif()

        if(exe_name STREQUAL "tests_vectorization")
            target_compile_options(${exe_name} PRIVATE
                /arch:AVX2 
                # /O2
                # /Qpar
                # /openmp
                /fp:fast
                # /fp:precise /fp:contract
                # /Qvec-report:2
                # /Qpar-report:2
            )
        endif()
        
endforeach()
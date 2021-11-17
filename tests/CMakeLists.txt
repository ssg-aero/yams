cmake_minimum_required(VERSION 3.17.0)

find_package(GTest REQUIRED)

project(quiss_tests LANGUAGES CXX)
if(MSVC)
    add_compile_options(/bigobj)
endif()

message(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

file(GLOB SRC_LIST "tests_*")
message(${SRC_LIST})
add_executable(${PROJECT_NAME} ${SRC_LIST})

if(MSVC) # xt
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
endif()

if ( MSVC )
    SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )

target_link_libraries(${PROJECT_NAME} 
                            GTest::gtest
                            GTest::gtest_main
                            xtensor 
                            xtensor::optimize 
                            xtensor::use_xsimd
                            ${VTK_LIBRARIES}
                            nlopt
                            gbs-render
                            iges
                            TBB::tbb
                    )

vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
)
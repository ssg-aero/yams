cmake_minimum_required(VERSION 3.17.0)

project(quiss_tests LANGUAGES CXX)
#
# quiss
#
message(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
#
# xtensors
#
# set(XTENSOR_USE_XSIMD BOOL:ON)
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)

file(GLOB SRC_LIST "tests_*")
message(${SRC_LIST})
add_executable(${PROJECT_NAME} ${SRC_LIST})

if(MSVC) # xt
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif()

if ( MSVC )
    SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )

target_link_libraries(${PROJECT_NAME} 
                            optimized gtest_dll
                            debug gtest_dlld
                            xtensor 
                            # xtensor::optimize 
                            # xtensor::use_xsimd
                    )
cmake_minimum_required(VERSION 3.17.0)

project(yams-cli)

if(NOT TARGET libyams)
    find_package(libyams REQUIRED)
endif ()

add_executable(
    yams MACOSX_BUNDLE 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/yams.cpp
)

target_link_libraries(
    yams PRIVATE
    libyams
)

install(TARGETS yams)
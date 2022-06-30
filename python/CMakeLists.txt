cmake_minimum_required(VERSION 3.17.0)

project(yams_bindings)

if(NOT TARGET libyams)
    find_package(libyams REQUIRED)
endif ()

find_package(pybind11 REQUIRED)

find_package(Python3 
    COMPONENTS 
        Interpreter 
        Developement
)
message("Python3_SITELIB:" ${Python3_SITELIB})

pybind11_add_module(yams_bindings
    src/yamsbind.cpp
)

target_link_libraries(
    ${PROJECT_NAME} PRIVATE 
    libyams
)

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME "yams_bindings"
)

install(TARGETS ${PROJECT_NAME} 
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pyams/)
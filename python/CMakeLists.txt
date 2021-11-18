cmake_minimum_required(VERSION 3.17.0)

project(pyams)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Python3 COMPONENTS Interpreter Developement)
find_package(pybind11 REQUIRED)

pybind11_add_module(
    ${PROJECT_NAME}
    yamsbind.cpp
)

target_link_libraries(
        ${PROJECT_NAME} PRIVATE 
        TBB::tbb 
)

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME "yams"
)

install(
    TARGETS ${PROJECT_NAME} 
    DESTINATION ${Python3_SITELIB}/${PROJECT_NAME}
)

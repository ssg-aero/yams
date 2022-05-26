cmake_minimum_required(VERSION 3.17.0)

project(pyams LANGUAGES CXX)

# set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Python3 
    COMPONENTS 
        Interpreter 
        Developement
)
message("Python3_SITELIB:" ${Python3_SITELIB})
find_package(pybind11 REQUIRED)

file(GLOB SRC_LIST "*.cpp")

pybind11_add_module(
    ${PROJECT_NAME}
    ${SRC_LIST}
)

target_link_libraries(
        ${PROJECT_NAME} PRIVATE 
        ${VTK_LIBRARIES}
        TBB::tbb 
        NLopt::nlopt
        gbs-render
)

vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
)

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME "yams"
)

install(
    TARGETS ${PROJECT_NAME} 
    DESTINATION ${Python3_SITELIB}/${PROJECT_NAME}
)

file(GLOB PY_SRC_LIST "*.py")

install(FILES ${PY_SRC_LIST}
        DESTINATION ${Python3_SITELIB}/${PROJECT_NAME}
)
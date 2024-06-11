cmake_minimum_required(VERSION 3.10)
project(simple_controller)
set(CMAKE_CXX_STANDARD 17)

# ###################
# cmake configures
# ###################
# add support for vscode clangd plugin
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ###################
# dependencies
# ###################
find_package(pinocchio CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -DEIGEN_DONT_PARALLELIZE -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -Ofast -flto -DEIGEN_NO_DEBUG")

set(RAISIM_OS linux)
find_package(raisim CONFIG REQUIRED)
set(Dependencies raisim::raisim)

# ########################
# src and linking
# ########################
include_directories(
    ${EIGEN3_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDE_DIRS}
    src
)

add_executable(main main.cpp src/controller.cpp src/env.cpp src/geometry.cpp src/kinematic.cpp)
target_link_libraries(main PRIVATE ${PINOCCHIO_LIBRARIES} ${Dependencies})

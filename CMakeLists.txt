cmake_minimum_required(VERSION 3.10)
project(controller)
set(CMAKE_CXX_STANDARD 17)

# ###################
# cmake configures
# ###################
# add support for vscode clangd plugin
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# colorful message
string(ASCII 27 ESC)
set(BOLD_RED "${ESC}[1;31m")
set(BOLD_GREEN "${ESC}[1;32m")
set(BOLD_YELLOW "${ESC}[1;33m")
set(BOLD_MAGENTA "${ESC}[1;35m")
set(BOLD_CYAN "${ESC}[1;36m")
set(CLR "${ESC}[m")

# ###################
# dependencies
# ###################
find_package(pinocchio CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DEIGEN_DONT_PARALLELIZE -march=native")
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

add_executable(main main.cpp src/controller.cpp src/env.cpp)
target_link_libraries(main PRIVATE ${PINOCCHIO_LIBRARIES} ${Dependencies})

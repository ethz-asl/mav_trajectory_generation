cmake_minimum_required(VERSION 2.8.3)
project(mav_trajectory_generation_example)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

set(CMAKE_MACOSX_RPATH 0)
add_definitions(-std=c++11)

############
# BINARIES #
############
cs_add_executable(trajectory_generation_example
        src/example_planner.cc
        src/example_planner_node.cc
)

cs_add_executable(trajectory_generation_6dof_example
        src/example_planner_6dof.cc
        src/example_planner_6dof_node.cc
)

##########
# EXPORT #
##########
cs_install()
cs_export()

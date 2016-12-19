cmake_minimum_required(VERSION 2.8.3)
project(mav_trajectory_generation)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

set(CMAKE_MACOSX_RPATH 0)
add_definitions(-std=c++11)

#############
# LIBRARIES #
#############
cs_add_library(${PROJECT_NAME}
  src/segment.cpp
)

############
# BINARIES #
############
#cs_add_executable(voxblox_rox
#  src/voxblox.cc
#)
#target_link_libraries(voxblox_rox ${PROJECT_NAME})

#########
# TESTS #
#########
#catkin_add_gtest(test_tsdf_map
#  test/test_tsdf_map.cc
#)
#target_link_libraries(test_tsdf_map ${PROJECT_NAME} ${catkin_LIBRARIES})

##########
# EXPORT #
##########
cs_install()
cs_export()

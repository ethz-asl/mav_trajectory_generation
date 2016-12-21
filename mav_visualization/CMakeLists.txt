cmake_minimum_required(VERSION 2.8.3)
project(mav_visualization)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

set(CMAKE_MACOSX_RPATH 0)
add_definitions(-std=c++11)

#############
# LIBRARIES #
#############
cs_add_library(${PROJECT_NAME}
  src/marker_group.cpp
  src/leica_marker.cpp
  src/hexacopter_marker.cpp
)

############
# BINARIES #
############
cs_add_executable(hexacopter_publisher
  src/hexacopter_publisher.cpp
)
target_link_libraries(hexacopter_publisher ${PROJECT_NAME})

cs_add_executable(leica_publisher
  src/leica_publisher.cpp
)
target_link_libraries(leica_publisher ${PROJECT_NAME})

##########
# EXPORT #
##########
cs_install()
cs_export()

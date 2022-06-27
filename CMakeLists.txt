cmake_minimum_required(VERSION 2.8)

project(calibration)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

add_executable(./build/calibration calibration.cc)
target_link_libraries(./build/calibration ${CERES_LIBRARIES})

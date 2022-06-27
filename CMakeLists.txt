cmake_minimum_required(VERSION 2.8)

project(calibration)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

add_executable(calibration calibration.cc)
target_link_libraries(calibration ${CERES_LIBRARIES})

add_executable(calibration_lamber calibration_lamber.cc)
target_link_libraries(calibration_lamber ${CERES_LIBRARIES})

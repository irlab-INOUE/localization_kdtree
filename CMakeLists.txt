cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
set (CMAKE_CXX_FLAGS "-O3 -g3 -std=c++17")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "")

project(kdTree)

include_directories(${PROJECT_SOURCE_DIR})
link_directories(${PROJECT_SOURCE_DIR})

find_package (OpenCV REQUIRED)

add_executable (global_localization global_localization.cpp map.cpp)
add_executable (slam slam.cpp map.cpp)
target_link_libraries(global_localization ${OpenCV_LIBRARIES})
target_link_libraries(slam ${OpenCV_LIBRARIES})

cmake_minimum_required(VERSION 3.1)
project(MSACONVERTER)


set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_BUILD_TYPE Debug)
set (CMAKE_CXX_STANDARD 11)

add_subdirectory(${PROJECT_SOURCE_DIR}/dependencies/libpll-2)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)



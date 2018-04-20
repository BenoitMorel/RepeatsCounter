

set(CMAKE_BUILD_TYPE Debug)

add_executable(convert  convert.cpp
  )

target_include_directories(convert PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/libpll-2/src/")
target_link_libraries(convert pll_static)

#target_link_libraries(multi-raxml ${MPI_LIBRARIES} dl gcc_s gcc)



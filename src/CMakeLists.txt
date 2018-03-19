
find_package(MPI REQUIRED)

set(CMAKE_BUILD_TYPE Debug)

add_executable(RepeatsCounter  RepeatsCounter.cpp
  )

#target_link_libraries(multi-raxml ${MPI_LIBRARIES} dl gcc_s gcc)



cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

project(fusion)


set(INCLUDE_DIR
    ../include
)

include_directories(${INCLUDE_DIR})

set(SOURCES
    ../src/main.cpp

    ../include/point.h
    ../src/point.cpp
 
    ../include/voxel.h
    ../src/voxel.cpp

    ../include/sdf.h
    ../src/sdf.cpp
)

add_executable(fusion ${SOURCES})

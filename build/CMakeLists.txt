cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

# choose compiler
set (CMAKE_CXX_COMPILER "/usr/bin/nvcc")
#set (CMAKE_CXX_COMPILER "/usr/bin/g++")

set (CMAKE_CXX_FLAGS "-pthread -lX11 -Dcimg_display=0")

project(fusion CXX)

set(INCLUDE_DIR
    ../include
)

include_directories(${INCLUDE_DIR})

set(SOURCES
    ../src/main.cpp
    ../src/point.cpp
    ../src/sdf.cpp
    ../src/canon_sdf.cpp
    ../src/fusion.cpp
    ../src/matrix.cpp
    ../src/gpu_sdf.cpp
)

add_executable(fusion ${SOURCES})

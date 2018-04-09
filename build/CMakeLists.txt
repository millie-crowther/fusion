cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-pthread -Dcimg_display=0")
project(fusion CXX)

set(INCLUDE_DIR
    ../include
)

include_directories(${INCLUDE_DIR})

set(SOURCES
    ../src/main.cpp
    ../src/point.cpp
    ../src/sdf.cpp
    ../src/canonical_sdf.cpp
    ../src/fusion.cpp
)

add_executable(fusion ${SOURCES})

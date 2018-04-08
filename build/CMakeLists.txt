cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

project(fusion)

set(INCLUDE_DIR
    ../include
    ../libpng
)

include_directories(${INCLUDE_DIR})

set(SOURCES
    ../src/point.cpp
    ../src/sdf.cpp
    ../src/fusion.cpp
    ../src/main.cpp
)

add_executable(fusion ${SOURCES})
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
target_link_libraries(fusion ${PNG_LIBRARY})

cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

project(fusion)


set(INCLUDE_DIR
    ../include
)

include_directories(${INCLUDE_DIR})

set(SOURCES
    ../src/main.cpp
)

add_executable(fusion ${SOURCES})

cmake_minimum_required(VERSION 3.6)
project(geometry)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp geometry.h geometry.h)
add_executable(geometry ${SOURCE_FILES})
cmake_minimum_required(VERSION 3.1)
project(hysteresis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    src/hysteresis.cpp
    src/PythonCopy.cpp
    #src/Iter.cpp
    #src/Iter.h
    src/PythonCopy.h
    )

add_executable(hysteresis ${SOURCE_FILES})
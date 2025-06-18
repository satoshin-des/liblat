cmake_minimum_required(VERSION 3.10)
project(LatLibrary)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR}/include)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB SRC_FILES
    ${CMAKE_SOURCE_DIR}/src/*.cpp
)

add_library(lat SHARED ${SRC_FILES})

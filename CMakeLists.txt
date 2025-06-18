cmake_minimum_required(VERSION 3.10)
project(LatLibrary)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR}/include)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB SRC_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)

add_library(lat SHARED ${SRC_FILES})

add_executable(test_exec ${CMAKE_SOURCE_DIR}/test.cpp)

target_link_libraries(test_exec lat)

target_include_directories(test_exec PRIVATE ${CMAKE_SOURCE_DIR}/include)


enable_testing()

add_test(NAME RunTestExec COMMAND test_exec)
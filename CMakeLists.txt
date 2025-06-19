cmake_minimum_required(VERSION 3.10)
project(LatLibrary)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_library(lat SHARED ${SRC_FILES})

target_include_directories(lat PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)

    add_executable(test_exec ${CMAKE_CURRENT_SOURCE_DIR}/test.cpp)

    target_link_libraries(test_exec lat)

    enable_testing()

    add_test(NAME RunTestExec COMMAND test_exec)

endif()

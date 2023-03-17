cmake_minimum_required(VERSION 3.24)
project(HwaUtil)

set(CMAKE_CXX_STANDARD 20)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/)

add_executable(try tryanything.cpp)

add_subdirectory(src)
add_subdirectory(test)


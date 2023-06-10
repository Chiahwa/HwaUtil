cmake_minimum_required(VERSION 3.20)
project(HwaUtil)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/)


if ("${MPI}" STREQUAL on)
    find_package(MPI REQUIRED)
    #find_package(SCALAPACK REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}" )
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")

    set(MPI_COMPILE_DEFINITIONS ${MPI_COMPILE_DEFINITIONS} __MPI__)

    message(STATUS "MPI_COMPILE_DEFINITIONS: ${MPI_COMPILE_DEFINITIONS}")
    message(STATUS "MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
    add_compile_definitions(${MPI_COMPILE_DEFINITIONS})
    message(STATUS "MPI is on")
else ()
    message(STATUS "MPI is off")
endif ()

add_executable(try tryanything.cpp)

target_link_libraries(try ${MPI_LIBRARIES})
target_include_directories(try PUBLIC ${MPI_INCLUDE_PATH})
message(STATUS ${SCALAPACK_LIBRARIES})
target_link_libraries(try ${SCALAPACK_LIBRARIES})

add_subdirectory(src)
add_subdirectory(exec)


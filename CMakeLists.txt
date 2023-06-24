cmake_minimum_required(VERSION 3.20)
project(HwaUtil)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/)

add_executable(try tryanything.cpp)


find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#find_package(GFortran REQUIRED)

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

    target_link_libraries(try ${MPI_LIBRARIES})

    if ("${SCALAPACK}" STREQUAL on)
        find_library(SCALAPACK_LIB scalapack REQUIRED)

        if (SCALAPACK_LIB)
            message(STATUS "SCALAPACK_LIB: ${SCALAPACK_LIB}")
            add_compile_definitions(__SCALAPACK__)
            message(STATUS "SCALAPACK is on")
            target_link_libraries(try ${SCALAPACK_LIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)
        else ()
            message(STATUS "SCALAPACK is off")
        endif ()
    endif ()
else ()
    message(STATUS "MPI is off")
endif ()

# OpenMP
if ("${OPENMP}" STREQUAL on)
    FIND_PACKAGE(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        message("OPENMP FOUND")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
    set(OPENMP_COMPILE_DEFINITIONS ${OPENMP_COMPILE_DEFINITIONS} __OPENMP__)
    add_compile_definitions(${OPENMP_COMPILE_DEFINITIONS})
endif ()


add_subdirectory(src)
add_subdirectory(exec)


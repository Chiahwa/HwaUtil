set(PROG_NAME ClassProject)
aux_source_directory(. SRC_LIST)
add_executable(${PROG_NAME} ${SRC_LIST})
target_link_libraries(${PROG_NAME} Mat_Demo ArgumentReader Timer)

include_directories(${CMAKE_SOURCE_DIR}/src)

#set(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/${PROG_NAME})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/${PROG_NAME})


# test file
execute_process(COMMAND ${CMAKE_COMMAND} -E
        make_directory "${EXECUTABLE_OUTPUT_PATH}/input"
        )
execute_process(COMMAND ${CMAKE_COMMAND} -E
        copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/test"
        "${EXECUTABLE_OUTPUT_PATH}/input"
        )

# MPI
if ("${MPI}" STREQUAL on)
    message(STATUS ${PROG_NAME} ": MPI is on")
    target_link_libraries(${PROG_NAME} ${MPI_LIBRARIES})
    # message(STATUS ${SCALAPACK_LIBRARIES})
    # target_link_libraries(${PROG_NAME} ${SCALAPACK_LIBRARIES})
endif ()

# OpenMP
if ("${OPENMP}" STREQUAL on)
    FIND_PACKAGE( OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message("OPENMP FOUND")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
    message(STATUS ${PROG_NAME} ": OpenMP is on")
    set(OPENMP_COMPILE_DEFINITIONS ${OPENMP_COMPILE_DEFINITIONS} __OPENMP__)
    add_compile_definitions(${OPENMP_COMPILE_DEFINITIONS})
    target_link_libraries(${PROG_NAME} ${OpenMP_CXX_LIBRARIES})
endif ()
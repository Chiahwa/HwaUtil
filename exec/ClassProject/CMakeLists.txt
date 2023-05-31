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
endif ()
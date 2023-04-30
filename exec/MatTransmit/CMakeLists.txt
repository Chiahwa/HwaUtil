set(PROG_NAME MatTransmit)
aux_source_directory(. DIR_LIB_SRCS)
add_executable(${PROG_NAME} ${DIR_LIB_SRCS})
target_link_libraries(${PROG_NAME} Mat_Demo ArgumentReader Timer)

include_directories(${CMAKE_SOURCE_DIR}/src)

#set(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/${PROG_NAME})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/${PROG_NAME})


# test file
configure_file(input.txt ${EXECUTABLE_OUTPUT_PATH}/input.txt COPYONLY)
configure_file(mat-32-16.txt ${EXECUTABLE_OUTPUT_PATH}/mat-32-16.txt COPYONLY)
configure_file(mat_1.txt ${EXECUTABLE_OUTPUT_PATH}/mat_1.txt COPYONLY)

if ("${MPI}" STREQUAL on)
    message(STATUS ${PROG_NAME} ": MPI is on")
    target_link_libraries(${PROG_NAME} ${MPI_LIBRARIES})
endif ()
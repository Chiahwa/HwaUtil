set(PROG_NAME MatManipTest)
aux_source_directory(. DIR_LIB_SRCS)
add_executable(${PROG_NAME} ${DIR_LIB_SRCS})
target_link_libraries(${PROG_NAME} Mat_Demo ArgumentReader Timer)

include_directories(${CMAKE_SOURCE_DIR}/src)


#set(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/${PROG_NAME})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/${PROG_NAME})


# test file



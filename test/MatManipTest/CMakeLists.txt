set(PROG_NAME MatManipTest)
aux_source_directory(. DIR_LIB_SRCS)
add_executable(${PROG_NAME} ${DIR_LIB_SRCS})
target_link_libraries(${PROG_NAME} Mat_Demo ArgumentReader Timer)

include_directories(${CMAKE_SOURCE_DIR}/src)


#set(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/${PROG_NAME})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/${PROG_NAME})


# test file
configure_file(args/matmul1.txt ${EXECUTABLE_OUTPUT_PATH}/args/matmul1.txt COPYONLY)
configure_file(args/matmul2.txt ${EXECUTABLE_OUTPUT_PATH}/args/matmul2.txt COPYONLY)
configure_file(args/matmul3.txt ${EXECUTABLE_OUTPUT_PATH}/args/matmul3.txt COPYONLY)
configure_file(args/matmul4.txt ${EXECUTABLE_OUTPUT_PATH}/args/matmul4.txt COPYONLY)
configure_file(args/matmul5.txt ${EXECUTABLE_OUTPUT_PATH}/args/matmul5.txt COPYONLY)
configure_file(args/diago1.txt ${EXECUTABLE_OUTPUT_PATH}/args/diago1.txt COPYONLY)
configure_file(args/diago2.txt ${EXECUTABLE_OUTPUT_PATH}/args/diago2.txt COPYONLY)
configure_file(args/diago3.txt ${EXECUTABLE_OUTPUT_PATH}/args/diago3.txt COPYONLY)

execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${EXECUTABLE_OUTPUT_PATH}/matmul")
execute_process(COMMAND ${CMAKE_COMMAND} -E
        copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/matmul"
        "${EXECUTABLE_OUTPUT_PATH}/matmul"
        )
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${EXECUTABLE_OUTPUT_PATH}/diago")
execute_process(COMMAND ${CMAKE_COMMAND} -E
        copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/diago"
        "${EXECUTABLE_OUTPUT_PATH}/diago"
        )
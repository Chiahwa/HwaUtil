set(PROG_NAME MatDemoTest)
add_executable(${PROG_NAME} MatDemoTest.cpp)
target_link_libraries(${PROG_NAME} Mat_Demo ArgumentReader)

include_directories(${CMAKE_SOURCE_DIR}/src)

#set(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/${PROG_NAME})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/${PROG_NAME})
# test file for MatDemoTest
configure_file(input.txt ${EXECUTABLE_OUTPUT_PATH}/input.txt COPYONLY)
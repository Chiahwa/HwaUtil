set(PROG_NAME ArgumentReaderTest)
add_executable(${PROG_NAME} ArgumentReaderTest.cpp)
target_link_libraries(${PROG_NAME} ArgumentReader)

include_directories(${CMAKE_SOURCE_DIR}/src)

#set(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/ArgumentReaderTest)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/${PROG_NAME})
# test file for ArgumentReaderTest
configure_file(input.txt ${EXECUTABLE_OUTPUT_PATH}/input.txt COPYONLY)
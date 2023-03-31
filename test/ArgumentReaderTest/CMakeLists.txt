add_executable(ArgumentReaderTest ArgumentReaderTest.cpp)
target_link_libraries(ArgumentReaderTest ArgumentReader)

include_directories(${CMAKE_SOURCE_DIR}/src)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/ArgumentReaderTest)

configure_file(input.txt ${EXECUTABLE_OUTPUT_PATH}/input.txt COPYONLY)
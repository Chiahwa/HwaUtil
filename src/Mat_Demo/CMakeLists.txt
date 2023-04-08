aux_source_directory(. DIR_LIB_SRCS)

set(LIB_NAME Mat_Demo)
add_library(${LIB_NAME} ${DIR_LIB_SRCS})
target_link_libraries(${PROG_NAME} Timer)
include_directories(${CMAKE_SOURCE_DIR}/src)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/${LIB_NAME})
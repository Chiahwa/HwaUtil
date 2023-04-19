aux_source_directory(. DIR_LIB_SRCS)
set(LIB_NAME ArgumentReader)

# 生成链接库
add_library(${LIB_NAME} STATIC  ${DIR_LIB_SRCS})
target_link_libraries(${LIB_NAME} Timer)
include_directories(${CMAKE_SOURCE_DIR}/src)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/${LIB_NAME})

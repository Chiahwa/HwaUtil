aux_source_directory(. DIR_LIB_SRCS)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/Mat_Demo)
add_library(Mat_Demo ${DIR_LIB_SRCS})

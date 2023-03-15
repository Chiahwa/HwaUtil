aux_source_directory(. DIR_LIB_SRCS)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/ArgumentReader)
add_library(ArgumentReader ${DIR_LIB_SRCS})

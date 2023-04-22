aux_source_directory(. DIR_LIB_SRCS)

set(LIB_NAME Mat_Demo)
add_library(${LIB_NAME} STATIC  ${DIR_LIB_SRCS})

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


target_link_libraries(${LIB_NAME}  Timer
        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
        )
include_directories(${CMAKE_SOURCE_DIR}/src )
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/${LIB_NAME})
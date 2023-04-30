aux_source_directory(. DIR_LIB_SRCS)

set(LIB_NAME Timer)
add_library(${LIB_NAME} STATIC ${DIR_LIB_SRCS})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/${LIB_NAME})
if ("${MPI}" STREQUAL on)
    message(STATUS ${LIB_NAME} ": MPI is on")
    target_link_libraries(${LIB_NAME} ${MPI_LIBRARIES})
endif ()


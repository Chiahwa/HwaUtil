aux_source_directory(. DIR_LIB_SRCS)

set(LIB_NAME Timer)
add_library(${LIB_NAME} STATIC ${DIR_LIB_SRCS})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/${LIB_NAME})

if ("${MPI}" STREQUAL on)
    target_link_libraries(${LIB_NAME} ${MPI_LIBRARIES})
    # target_compile_definitions(${LIB_NAME} PUBLIC ${MPI_COMPILE_DEFINITIONS})
    # target_include_directories(${LIB_NAME} PUBLIC ${MPI_INCLUDE_PATH})
    # message(${MPI_COMPILE_DEFINITIONS})
    # target_compile_definitions(${LIB_NAME} PUBLIC ${MPI_COMPILE_DEFINITIONS})
endif ()
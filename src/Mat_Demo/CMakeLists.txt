aux_source_directory(. DIR_LIB_SRCS)

set(LIB_NAME Mat_Demo)
add_library(${LIB_NAME} STATIC  ${DIR_LIB_SRCS})



target_link_libraries(${LIB_NAME}  Timer
        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
        )

# MPI
if ("${MPI}" STREQUAL on)
    message(STATUS ${LIB_NAME} ": MPI is on")
    target_link_libraries(${LIB_NAME} ${MPI_LIBRARIES})
    # Scalapack
    if ("${SCALAPACK}" STREQUAL on)
        message(STATUS ${LIB_NAME} ": Scalapack is on")
        target_link_libraries(${LIB_NAME} ${SCALAPACK_LIB} -lgfortran)
    endif ()
endif ()
include_directories(${CMAKE_SOURCE_DIR}/src )
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/${LIB_NAME})
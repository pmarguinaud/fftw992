cmake_minimum_required(VERSION 3.10)

project(dummy_mkl_user LANGUAGES Fortran)

# Use ecbuild's FindMKL.cmake
list(INSERT CMAKE_MODULE_PATH 0 /home/gmap/mrpm/marguina/ecbuild/cmake)

find_package(MKL REQUIRED)
find_package(BLAS REQUIRED)

message(STATUS "========================================")
message(STATUS "MKL_FOUND        = ${MKL_FOUND}")
message(STATUS "MKL_INCLUDE_DIRS = ${MKL_INCLUDE_DIRS}")
message(STATUS "MKL_LIBRARIES    = ${MKL_LIBRARIES}")
message(STATUS "BLAS_LIBRARIES   = ${BLAS_LIBRARIES}")
message(STATUS "========================================")

add_executable(dummy main.F90)
target_include_directories(dummy PRIVATE ${MKL_INCLUDE_DIRS})
target_link_libraries(dummy ${MKL_LIBRARIES} ${BLAS_LIBRARIES})

enable_testing()
add_test(NAME dummy_mkl_tests COMMAND dummy)

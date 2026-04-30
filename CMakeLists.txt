cmake_minimum_required(VERSION 3.10)

project(mkl LANGUAGES Fortran VERSION 1.0.0)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# Force lib/ instead of platform-dependent lib64/
set(CMAKE_INSTALL_LIBDIR lib)

#-----------------------------------------------------------------------
# Options
#-----------------------------------------------------------------------
option(BUILD_TESTING "Build tests" ON)

#-----------------------------------------------------------------------
# Find BLAS and create a modern imported target if missing
#-----------------------------------------------------------------------
find_package(BLAS REQUIRED)

if(NOT TARGET BLAS::BLAS)
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  if(BLAS_LIBRARIES)
    set_target_properties(BLAS::BLAS PROPERTIES
      INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}"
    )
  endif()
  if(BLAS_INCLUDE_DIRS)
    set_target_properties(BLAS::BLAS PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}"
    )
  endif()
endif()

#-----------------------------------------------------------------------
# Output directories inside the build tree
#-----------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Copy public header into build tree so the layout mirrors the install tree
configure_file(${CMAKE_SOURCE_DIR}/include/fftw3.f03
               ${CMAKE_BINARY_DIR}/include/fftw3.f03
               COPYONLY)

#-----------------------------------------------------------------------
# Library source files (everything except the test driver)
#-----------------------------------------------------------------------
set(LIB_SRCS
    src/set99_dp.F90
    src/fft992_dp.F90
    src/set99_sp.F90
    src/fft992_sp.F90
    src/fftw3.F90
)

#-----------------------------------------------------------------------
# Create both static and shared libraries
#-----------------------------------------------------------------------
add_library(mkl_static STATIC ${LIB_SRCS})
set_target_properties(mkl_static PROPERTIES
    OUTPUT_NAME mkl
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_link_libraries(mkl_static PUBLIC BLAS::BLAS)

add_library(mkl_shared SHARED ${LIB_SRCS})
set_target_properties(mkl_shared PROPERTIES
    OUTPUT_NAME mkl
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_link_libraries(mkl_shared PUBLIC BLAS::BLAS)

# Convenience alias target
add_library(mkl ALIAS mkl_shared)

#-----------------------------------------------------------------------
# Install targets + export
#-----------------------------------------------------------------------
install(TARGETS mkl_static mkl_shared
    EXPORT mkl-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES include/fftw3.f03 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.mod"
)

install(EXPORT mkl-targets
    FILE mkl-targets.cmake
    NAMESPACE mkl::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mkl
)

#-----------------------------------------------------------------------
# Package configuration file
#-----------------------------------------------------------------------
configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/mkl-config.cmake.in"
    "${CMAKE_BINARY_DIR}/mkl-config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mkl
)

write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/mkl-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_BINARY_DIR}/mkl-config.cmake"
    "${CMAKE_BINARY_DIR}/mkl-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mkl
)

#-----------------------------------------------------------------------
# Test executable (linked against the shared library)
#-----------------------------------------------------------------------
if(BUILD_TESTING)
    enable_testing()
    add_executable(test_fftw3 tests/test_fftw3.F90)
    target_link_libraries(test_fftw3 mkl_shared)
    add_test(NAME fftw3_wrapper COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_fftw3)
endif()

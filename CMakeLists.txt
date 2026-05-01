cmake_minimum_required(VERSION 3.10)

project(fftw3ial LANGUAGES Fortran VERSION 1.0.0)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# Force lib/ instead of platform-dependent lib64/
set(CMAKE_INSTALL_LIBDIR lib)

#-----------------------------------------------------------------------
# Options
#-----------------------------------------------------------------------
option(BUILD_TESTING "Build tests" ON)

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
# Library: libfftw3 (double precision)
#-----------------------------------------------------------------------
add_library(fftw3_static STATIC
    src/set99_dp.F90
    src/fft992_dp.F90
    src/fftw3_dp.F90
)
set_target_properties(fftw3_static PROPERTIES
    OUTPUT_NAME fftw3
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_include_directories(fftw3_static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(fftw3 SHARED
    src/set99_dp.F90
    src/fft992_dp.F90
    src/fftw3_dp.F90
)
set_target_properties(fftw3 PROPERTIES
    OUTPUT_NAME fftw3
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_include_directories(fftw3 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

#-----------------------------------------------------------------------
# Library: libfftw3f (single precision)
#-----------------------------------------------------------------------
add_library(fftw3f_static STATIC
    src/set99_sp.F90
    src/fft992_sp.F90
    src/fftw3_sp.F90
)
set_target_properties(fftw3f_static PROPERTIES
    OUTPUT_NAME fftw3f
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_include_directories(fftw3f_static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(fftw3f SHARED
    src/set99_sp.F90
    src/fft992_sp.F90
    src/fftw3_sp.F90
)
set_target_properties(fftw3f PROPERTIES
    OUTPUT_NAME fftw3f
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_include_directories(fftw3f PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

#-----------------------------------------------------------------------
# Install targets + export
#-----------------------------------------------------------------------
install(TARGETS fftw3 fftw3f fftw3_static fftw3f_static
    EXPORT fftw3ial-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES include/fftw3.f03 include/mkl.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.mod"
)

#-----------------------------------------------------------------------
# Create library symlinks so that ecbuild FindMKL.cmake discovers us
#-----------------------------------------------------------------------
install(CODE "
    set(_mkl_libdir \"\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")
    message(STATUS \"Creating MKL compatibility symlinks in \"\${_mkl_libdir})
    # .so symlinks (used by FindMKL.cmake on Linux)
    execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink libfftw3.so libmkl_intel_lp64.so
                    WORKING_DIRECTORY \${_mkl_libdir})
    execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink libfftw3.so libmkl_sequential.so
                    WORKING_DIRECTORY \${_mkl_libdir})
    execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink libfftw3.so libmkl_core.so
                    WORKING_DIRECTORY \${_mkl_libdir})
    # .a symlinks (fallback for static-only installs)
    execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink libfftw3.a libmkl_intel_lp64.a
                    WORKING_DIRECTORY \${_mkl_libdir})
    execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink libfftw3.a libmkl_sequential.a
                    WORKING_DIRECTORY \${_mkl_libdir})
    execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink libfftw3.a libmkl_core.a
                    WORKING_DIRECTORY \${_mkl_libdir})
")

install(EXPORT fftw3ial-targets
    FILE fftw3ial-targets.cmake
    NAMESPACE fftw3ial::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftw3ial
)

#-----------------------------------------------------------------------
# Package configuration file
#-----------------------------------------------------------------------
configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/fftw3ial-config.cmake.in"
    "${CMAKE_BINARY_DIR}/fftw3ial-config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftw3ial
)

write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/fftw3ial-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_BINARY_DIR}/fftw3ial-config.cmake"
    "${CMAKE_BINARY_DIR}/fftw3ial-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftw3ial
)

#-----------------------------------------------------------------------
# Test executable (linked against both shared libraries)
#-----------------------------------------------------------------------
if(BUILD_TESTING)
    enable_testing()
    add_executable(test_fftw3 tests/test_fftw3.F90)
    target_link_libraries(test_fftw3 fftw3 fftw3f)
    add_test(NAME fftw3_wrapper COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_fftw3)
endif()

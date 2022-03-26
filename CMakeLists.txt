project(SMPPHAT)
cmake_minimum_required(VERSION 2.4.6)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PC_FFTW3 REQUIRED fftw3f)

include_directories("${PROJECT_SOURCE_DIR}/include")

set(CMAKE_BUILD_TYPE "Release")

ADD_DEFINITIONS(-DUSE_SIMD)
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-msse3 -Wall -ffast-math")

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")

set(SRC_SMPPHAT
  src/signal.c
  src/system.c
)

add_library(smpphat SHARED
  ${SRC_SMPPHAT}
)

target_link_libraries(smpphat
  ${PC_FFTW3_LIBRARIES}
  m
)

add_executable(ssl
  demo/ssl.c
)
target_link_libraries(ssl
    smpphat
)

add_executable(test
  demo/test.c
)
target_link_libraries(test
    smpphat
)

cmake_minimum_required(VERSION 3.11)
project(QA_opt)
include(FetchContent)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

message("CMake install dir ${CMAKE_MODULE_PATH}")
message("CMake build dir ${CMAKE_CURRENT_BINARY_DIR}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# ========= Correct RPATHs for mac osx
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_MACOSX_RPATH 1)
endif()

#setup the main executable
set(APP_EXE main)
add_executable( ${APP_EXE} main.cpp)

target_compile_definitions( ${APP_EXE} PRIVATE UNIT_TESTING_MATERIALS=1)
#target_link_libraries( ${APP_EXE})

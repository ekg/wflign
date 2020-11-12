cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(edlign VERSION 0.0.1)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.                                                                                                      
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.                                                                                                             
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(deps/edlib EXCLUDE_FROM_ALL)
add_subdirectory(deps/wflambda EXCLUDE_FROM_ALL)

add_executable(edlign
  src/main.cpp
  src/edlign.cpp
  src/wfa_edit_callback.cpp)

target_link_libraries(edlign
  edlib
  wflambda
  z)

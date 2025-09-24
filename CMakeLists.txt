cmake_minimum_required(VERSION 3.16)
project(BitDraw)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/$<CONFIGURATION>")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/$<CONFIGURATION>")

# SDL as subproject
add_subdirectory(vendored/SDL EXCLUDE_FROM_ALL)

# Add your app
add_executable(bitdraw main.cpp)

# Add include path for headers in project root (so #include "OB.h" works)
target_include_directories(bitdraw PRIVATE ${CMAKE_SOURCE_DIR})

# Link SDL
target_link_libraries(bitdraw PRIVATE SDL3::SDL3)

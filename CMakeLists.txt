cmake_minimum_required(VERSION 3.16)
project(SDL3Skeleton)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/$<CONFIGURATION>")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/$<CONFIGURATION>")

# Add SDL3 as a subproject
add_subdirectory(vendored/SDL EXCLUDE_FROM_ALL)

# Build our app
add_executable(sdl3_app main.cpp)
target_link_libraries(sdl3_app PRIVATE SDL3::SDL3)

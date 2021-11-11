cmake_minimum_required(VERSION 3.5)

project(example0 LANGUAGES C)
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/")
file(GLOB SRC_FILES "${CMAKE_CURRENT_LIST_DIR}/src/*.c")
add_executable(${PROJECT_NAME} ${SRC_FILES})

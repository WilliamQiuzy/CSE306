cmake_minimum_required(VERSION 3.0.0)
project(lab2 VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

set(CMAKE_CXX_FLAGS_RELEASE " -O3 -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG " -O3 -fopenmp")


file(GLOB SOURCES
    "${PROJECT_SOURCE_DIR}/*.cpp"
    "${PROJECT_SOURCE_DIR}/*.c"
    "${PROJECT_SOURCE_DIR}/liblbfgs/*.cpp"
    "${PROJECT_SOURCE_DIR}/liblbfgs/*.c"
)

add_executable(lab2 ${SOURCES})

target_include_directories(lab2 PRIVATE "${PROJECT_SOURCE_DIR}/liblbfgs")



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


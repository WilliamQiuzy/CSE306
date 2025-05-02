cmake_minimum_required(VERSION 3.0.0)
project(cse306 VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

set(SOURCE_FILES 
    Assignment1/main.cpp
    # Assignment1/vector_class.cpp
)


add_executable(cse306 ${SOURCE_FILES})

file(COPY ${CMAKE_SOURCE_DIR}/Assignment1/model/ DESTINATION ${CMAKE_BINARY_DIR}/model/)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

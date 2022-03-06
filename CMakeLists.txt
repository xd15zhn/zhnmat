cmake_minimum_required(VERSION 3.12)
project(zhnmat)
set(CMAKE_BUILD_TYPE release)
add_library(${CMAKE_PROJECT_NAME} STATIC matrix.cpp matrixcontrib.cpp vector3d.cpp)

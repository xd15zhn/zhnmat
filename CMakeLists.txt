cmake_minimum_required(VERSION 3.12)
project(zhnmat VERSION 1.1.8)
set(CMAKE_BUILD_TYPE release)
add_library(${CMAKE_PROJECT_NAME} STATIC matrix.cpp matrixcontrib.cpp vector3d.cpp)
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

install(TARGETS zhnmat
    EXPORT zhnmatConfig
    ARCHIVE DESTINATION lib
)
install(EXPORT zhnmatConfig
    FILE zhnmatConfig.cmake
    NAMESPACE zhnmat::
    DESTINATION lib/cmake/zhnmat
)
INSTALL(FILES zhnmat.hpp
    DESTINATION ${CMAKE_INSTALL_PREFIX}
)

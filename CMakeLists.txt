cmake_minimum_required(VERSION 3.12)
project(zhnmat VERSION 1.3.7)
set(CMAKE_BUILD_TYPE release)
message(STATUS "project name: ${CMAKE_PROJECT_NAME}")
message(STATUS "project version: ${PROJECT_VERSION}")
message(STATUS "project build type: ${CMAKE_BUILD_TYPE}")
add_library(${CMAKE_PROJECT_NAME} STATIC
    matrix.cpp
    vector3d.cpp
    extra.cpp
)

if (WIN32)
    set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
    set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
elseif (UNIX)
    include(GNUInstallDirs)
endif ()
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "INSTALL_CONFIGDIR: ${INSTALL_CONFIGDIR}")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR}")

option(USE_TRACELOG "Dependent library tracelog, used to print logs." ON)
option(USE_EXTRA "Use extra functions." ON)

if(USE_TRACELOG)
    message(STATUS "Use dependent library tracelog.")
    find_package(tracelog REQUIRED)
    message(STATUS "tracelog_VERSION: ${tracelog_VERSION}")
    message(STATUS "tracelog_DIR: ${tracelog_DIR}")
    message(STATUS "tracelog_LIBS: ${tracelog_LIBS}")
    message(STATUS "tracelog_INCLUDE_DIRS: ${tracelog_INCLUDE_DIRS}")
    add_definitions(-DUSE_TRACELOG)
    target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${tracelog_LIBS})
    target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${tracelog_INCLUDE_DIRS})
endif()
if(USE_EXTRA)
    add_definitions(-DUSE_EXTRA)
endif ()

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)
install(TARGETS ${PROJECT_NAME}
    EXPORT install_targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(EXPORT install_targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${INSTALL_CONFIGDIR}
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/zhnmat.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

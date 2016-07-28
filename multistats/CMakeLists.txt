cmake_minimum_required (VERSION 3.1)

project(CollectStats CXX)
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../external/alglib/src/*.cpp)

# set our library and executable destination dirs
set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
  
add_executable(CollectStats main.cpp ${SOURCES})
   
set_property(TARGET CollectStats PROPERTY CXX_STANDARD 11)
set_property(TARGET CollectStats PROPERTY CXX_STANDARD_REQUIRED ON)

find_package(ZLIB REQUIRED)
include_directories(../external/boost_1_58_0)
include_directories(../external/alglib/src)
target_link_libraries( CollectStats boost_iostreams z)

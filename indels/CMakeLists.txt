cmake_minimum_required (VERSION 3.1)

project(Compare_VCFs_indels CXX)

# set our library and executable destination dirs
set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
  
add_executable(Compare_VCFs_indels main.cpp)
   
set_property(TARGET Compare_VCFs_indels PROPERTY CXX_STANDARD 11)
set_property(TARGET Compare_VCFs_indels PROPERTY CXX_STANDARD_REQUIRED ON)

find_package(ZLIB REQUIRED)
include_directories(../external/boost_1_58_0)
target_link_libraries( Compare_VCFs_indels boost_iostreams z)


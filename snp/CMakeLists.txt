cmake_minimum_required (VERSION 3.1)

project(Compare_VCFs CXX)
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

file(GLOB_RECURSE SOURCES ../external/alglib/src/*.cpp)

# set our library and executable destination dirs
set( EXECUTABLE_OUTPUT_PATH "../bin" )
set(CMAKE_PREFIX_PATH "/store/local/")
  
add_executable(Compare_VCFs main.cpp ${SOURCES})
   
set_property(TARGET Compare_VCFs PROPERTY CXX_STANDARD 11)
set_property(TARGET Compare_VCFs PROPERTY CXX_STANDARD_REQUIRED ON)

find_package(ZLIB REQUIRED)
# place the boost library into the external section
# download from http://www.boost.org/users/history/version_1_58_0.html
include_directories(../external/boost_1_58_0)
# we are expecting the ALGLIB CPP code being in the external directory
# to download see http://www.alglib.net/download.php
# usually the library extracts into a directory called "cpp"
# so do something like
# VCF_Comparisons/external $ tar xzvf alglib-3.10.0.cpp.gpl.tgz
# VCF_Comparisons/external $ mv cpp alglib
include_directories(../external/alglib/src)
target_link_libraries( Compare_VCFs boost_iostreams z)

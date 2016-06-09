cmake_minimum_required (VERSION 3.1)

project(Compare_VCFs CXX)
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

file(GLOB_RECURSE SOURCES /pica/h1/pelin/External/alglib/*.cpp)

 # set our library and executable destination dirs
  set( EXECUTABLE_OUTPUT_PATH "/pica/h1/pelin/bin" )
  
   add_executable(Compare_VCFs main.cpp ${SOURCES})
   
  set_property(TARGET Compare_VCFs PROPERTY CXX_STANDARD 11)
  set_property(TARGET Compare_VCFs PROPERTY CXX_STANDARD_REQUIRED ON)

  find_package(ZLIB REQUIRED)
  include_directories(/pica/h1/pelin/External/boost_1_58_0)
  include_directories(/pica/h1/pelin/External/alglib)
  target_link_libraries( Compare_VCFs boost_iostreams z)




# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(dataEval CXX C)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
# Put executables into a bin
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_MACOSX_RPATH 1)  # only needed to suppress cmake warning

# ROOT dict generation
add_custom_command (
  OUTPUT ${PROJECT_SOURCE_DIR}/src/SRawEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SRawEventDict.h
  COMMAND rootcint
  ARGS -f ${PROJECT_SOURCE_DIR}/src/SRawEventDict.cxx -c ${PROJECT_SOURCE_DIR}/inc/SRawEvent.h ${PROJECT_SOURCE_DIR}/inc/SRawEventLinkDef.h
  )

include_directories(${PROJECT_SOURCE_DIR}/inc)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx   ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h ${PROJECT_SOURCE_DIR}/include/*.hh)

# ROOT
find_program(ROOTCONF "root-config")
if(ROOTCONF)
  message("-- Detecting ROOT:    found at ${ROOTCONF}")
else()
  message(FATAL_ERROR "-- Detecting ROOT:    not found")
endif()
execute_process(COMMAND root-config --prefix OUTPUT_VARIABLE ROOT_PREFIX  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --libs   OUTPUT_VARIABLE ROOT_LINK    OUTPUT_STRIP_TRAILING_WHITESPACE)
set(ROOT_LINK "${ROOT_LINK}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x ${ROOT_CFLAGS}")
add_library(dataEval SHARED ${sources} ${headers} ${PROJECT_SOURCE_DIR}/src/SRawEventDict.cxx)
target_link_libraries(dataEval "${ROOT_LINK}")

file(GLOB execs ${PROJECT_SOURCE_DIR}/exe/*.cxx)
foreach(exec ${execs})
  string(REPLACE  "${PROJECT_SOURCE_DIR}/exe/" "" exec ${exec})
  string(REPLACE  ".cxx" "" exec ${exec})
  message(STATUS "Adding targets to the compile list: ${exec}")

  add_executable(${exec} ${PROJECT_SOURCE_DIR}/exe/${exec}.cxx)
  target_link_libraries(${exec} ${PROJECT_SOURCE_DIR}/lib/libdataEval.so)
  target_link_libraries(${exec} "${ROOT_LINK}")
  #  target_include_directories(${exec} PUBLIC ${PROJECT_SOURCE_DIR}/src)
  add_dependencies(${exec} dataEval)
endforeach(exec)

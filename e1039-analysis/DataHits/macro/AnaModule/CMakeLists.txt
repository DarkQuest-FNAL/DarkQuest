# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(anamodule CXX)

file(GLOB          sources ${PROJECT_SOURCE_DIR}/*.cxx)
file(GLOB     dist_headers ${PROJECT_SOURCE_DIR}/*.h)
file(GLOB non_dist_headers ${PROJECT_SOURCE_DIR}/*LinkDef.h)
list(REMOVE_ITEM dist_headers ${non_dist_headers})

# ROOT dict generation
add_custom_command (
  OUTPUT anamodule_Dict.cc
  COMMAND rootcint
  ARGS -f anamodule_Dict.cc -c
    -I$ENV{OFFLINE_MAIN}/include/ -I${ROOT_PREFIX}/include/
    ${dist_headers}
    ${PROJECT_SOURCE_DIR}/LinkDef.h
  DEPENDS ${dist_headers} ${PROJECT_SOURCE_DIR}/LinkDef.h
)

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

#include_directories("${PROJECT_SOURCE_DIR}/" "$ENV{MY_INSTALL}/include/" "$ENV{OFFLINE_MAIN}/include/")
#link_directories("./" "$ENV{MY_INSTALL}/lib/" "$ENV{OFFLINE_MAIN}/lib/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -I$ENV{OFFLINE_MAIN}/include/ ${ROOT_CFLAGS}")

add_library(anamodule SHARED ${sources} anamodule_Dict.cc)
target_link_libraries(anamodule  fun4all interface_main geom_svc  phg4hit g4detectors)
install(TARGETS anamodule DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

# Install the pcm files in case of ROOT 6.
execute_process(COMMAND root-config --version OUTPUT_VARIABLE ROOT_VER)
string(SUBSTRING ${ROOT_VER} 0 1 ROOT_VER)
if (ROOT_VER GREATER 5)
   add_custom_target(install_pcm ALL COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX}/lib COMMAND cp -up *_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)
   add_dependencies(install_pcm anamodule)
endif()

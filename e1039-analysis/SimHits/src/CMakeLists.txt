# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(sim_ana CXX)

file(GLOB          sources ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB     dist_headers ${PROJECT_SOURCE_DIR}/*.h)
file(GLOB non_dist_headers ${PROJECT_SOURCE_DIR}/*LinkDef.h)
list(REMOVE_ITEM dist_headers ${non_dist_headers})

# ROOT dict generation
add_custom_command (
  OUTPUT sim_ana_Dict.cc
  COMMAND rootcint
  ARGS -f sim_ana_Dict.cc -noIncludePaths -inlineInputHeader -c -p
    -I${PROJECT_SOURCE_DIR}  -I${ROOT_PREFIX}/include
    ${dist_headers}  ${PROJECT_SOURCE_DIR}/LinkDef.h
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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -I$ENV{OFFLINE_MAIN}/include/ ${ROOT_CFLAGS}")

add_library(sim_ana SHARED ${sources} sim_ana_Dict.cc)
target_link_libraries(sim_ana  fun4all interface_main geom_svc ktracker phg4hit PHPythia8 g4detectors UtilAna)
install(TARGETS sim_ana DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

install(FILES ${dist_headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

# Install the pcm files.
add_custom_target(install_pcm ALL COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX}/lib COMMAND cp -up *_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)
add_dependencies(install_pcm sim_ana)

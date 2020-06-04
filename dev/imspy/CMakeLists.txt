#compile this project with -A x64 flag on msvc
cmake_minimum_required(VERSION 3.5)

project(readPTU
	LANGUAGES CXX C
	VERSION 1.0)

######################################################
# export a .pyd python library
######################################################

set( SOURCES
	src/bindings.cpp
	src/ProcessPhotonStream.cpp
	)

#generate python module
add_subdirectory( thirdparty/pybind11 )
pybind11_add_module(imspy ${SOURCES})

target_include_directories(imspy PUBLIC
	thirdparty/
	include/
	)

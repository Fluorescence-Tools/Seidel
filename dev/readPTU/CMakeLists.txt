cmake_minimum_required(VERSION 3.5)

#find_package(Python 3.7 COMPONENTS Interpreter Development REQUIRED)

project(readPTU
	LANGUAGES CXX C
	VERSION 1.0)

#######################################################
# create a library
#######################################################

#generate a static library from the library sources
add_library(readPTU_c SHARED
	src/ProcessPhotonStream.cpp
	)

set_property(TARGET readPTU_c PROPERTY CXX_STANDARD 11)


target_include_directories(readPTU_c PUBLIC
	${PROJECT_SOURCE_DIR}/include
	${PROJECT_SOURCE_DIR}/thirdparty
	)

######################################################
# export a .pyd python library
######################################################

set( SOURCES
	src/ProcesSPhotonStream.cpp
	)

#generate python module
add_subdirectory( thirdparty/pybind11 )
pybind11_add_module(readPTU_py ${SOURCES})

target_include_directories(readPTU_py PUBLIC
	thirdparty/
	include/
	)

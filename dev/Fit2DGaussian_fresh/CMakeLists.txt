cmake_minimum_required(VERSION 3.5)

project(Fit2DGaussian
        LANGUAGES CXX C
        VERSION 1.0)

add_subdirectory( thirdparty/ap )


        
############################################################
# Create a library
############################################################

#Generate the static library from the library sources
add_library(Fit2DGaussian SHARED 
     src/fit2DGaussian.cpp
     src/i_lbfgs.cpp
     src/lbfgs.cpp
     src/twoIstar.cpp
)

target_include_directories(Fit2DGaussian
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/thirdparty
)

target_link_libraries( Fit2DGaussian PRIVATE LibAP )

############################################################
# export a pybind .pyd python library
############################################################

set( SOURCES 
	src/pywrap.cpp
	src/fit2DGaussian.cpp
	src/twoIstar.cpp
	src/i_lbfgs.cpp
	src/lbfgs.cpp
	src/bindings.cpp
)

#Generate python module
add_subdirectory( thirdparty/pybind11 )
pybind11_add_module(GaussFits ${SOURCES} )

target_include_directories(GaussFits PUBLIC 
	thirdparty/
	include/
)

target_link_libraries(GaussFits PRIVATE LibAP )

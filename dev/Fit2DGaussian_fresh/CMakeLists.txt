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
)

target_link_libraries( Fit2DGaussian LibAP )

cmake_minimum_required(VERSION 3.5)

project(LibAP
        LANGUAGES CXX C
        VERSION 1.0)


############################################################
# Create a library
############################################################

#Generate the static library from the library sources
add_library(LibAP STATIC 
    src/ap.cpp
)

target_include_directories(LibAP
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)

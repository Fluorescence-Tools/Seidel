cmake_minimum_required(VERSION 2.8.12)
project(Fit2DGaussian)

add_library(Fit2DGaussian SHARED "src/fit2DGaussian.cpp")


#look in folder src for c files and headers
#include_directories(src)

target_include_directories(Fit2DGaussian PUBLIC src)

#find_library(Istar_lib NAMES twoIstar.cpp PATHS "src")
#message(STATUS "library path Istar_lib is " ${Istar_lib})
#target_link_libraries(Fit2DGaussian Istar_lib)


set(SOURCES "src/twoIstar.cpp")
target_link_libraries(Fit2DGaussian src/twoIstar.h src/twoIstar.h)


set(LIBRARY_OUTPUT_PATH build)
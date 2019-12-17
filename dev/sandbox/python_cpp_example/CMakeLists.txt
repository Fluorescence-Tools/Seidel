cmake_minimum_required(VERSION 2.8.12)
project(python_cpp_example)

#set source dir
set(SOURCE_DIR "src/python_cpp_example")
#Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR} lib)
set(SOURCES "${SOURCE_DIR}/math.cpp")

add_library(math_lib SHARED "${SOURCE_DIR}/math.cpp")

#Generate python module
add_subdirectory(lib/pybind11)
pybind11_add_module(python_cpp_example ${SOURCES} "${SOURCE_DIR}/bindings.cpp")
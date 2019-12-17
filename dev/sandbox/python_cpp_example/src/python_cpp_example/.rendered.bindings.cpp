

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "math.hpp"

namespace py = pybind11;

PYBIND11_PLUGIN(python_cpp_example)
{
	py::module m("python_cpp_example");
	m.def("add", &add);
	m.def("subtract", &subtract);
	m.def("inv", &inv);
	return m.ptr();
}
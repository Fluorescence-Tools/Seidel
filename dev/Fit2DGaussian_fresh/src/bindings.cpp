//pybind11 wrapper calls wraps all function in eigen_wrap
//Author: Nicolaas van der Voort, AG Seidel, HHU Düsseldorf
//Created December 16, 2019

#include "eigen_wrap.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;

PYBIND11_PLUGIN(GaussFits)
{
	py::module m("GaussFits");
	m.def("inv", &inv);
	m.def("subtract", &subtract);
	return m.ptr();
}
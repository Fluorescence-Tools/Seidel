//pybind11 wrapper calls wraps all function in eigen_wrap
//Author: Nicolaas van der Voort, AG Seidel, HHU Düsseldorf
//Created December 16, 2019

#include "pywrap.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;

PYBIND11_PLUGIN(GaussFits)
{
	py::module m("GaussFits");
	m.def("inv", &inv);
	m.def("subtract", &subtract);
	m.def("Fit2DGauss", &Fit2DGaussian_pywrap);
	m.def("model2DGaussian", &model2DGaussian_pywrap);
	m.def("modelTwo2DGaussian", &modelTwo2DGaussian_pywrap);
	m.def("modelThree2DGaussian", &modelThree2DGaussian_pywrap);
	m.def("twoIstar", &W2DG_pywrap);
	return m.ptr();
}
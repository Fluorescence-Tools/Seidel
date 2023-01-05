//pybind11 wrapper calls wraps all function in eigen_wrap
//Author: Nicolaas van der Voort, AG Seidel, HHU Düsseldorf
//Created December 16, 2019

#include "ProcessPhotonStream.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> //needed for std::vector

namespace py = pybind11;

PYBIND11_PLUGIN(readPTUv2)
{
	py::module m("readPTUv2");
	m.def("inv", &inv);
	m.def("subtract", &subtract);
	
	py::class_<Pet>(m, "Pet")
		.def(py::init<const std::string &>())
		.def("setName", &Pet::setName)
		.def("getName", &Pet::getName)
		.def_readwrite("name", &Pet::name);
	//this way of binding passes the constructor, and allows for reading and writing the other variables.
	py::class_<imOpts>(m, "imOpts")
		.def(py::init<std::vector<int>, float, float, long long, float, float>())
		.def_readwrite("line_ids", &imOpts::line_ids)
		.def_readwrite("dwelltime", &imOpts::dwelltime)
		.def_readwrite("counttime", &imOpts::counttime)
		.def_readwrite("NumRecords", &imOpts::NumRecords)
		.def_readwrite("linestep", &imOpts::linestep, 
		"e.g. 1 FRET line, 1 PIE line, 10nm px: linstep = 5")
		.def_readwrite("pxsize", &imOpts::pxsize);
	
	py::class_<ph>(m, "ph")
		.def(py::init<int, long long, unsigned char, long long, float, float, int>())
		.def_readwrite("tac", &ph::tac)
		.def_readwrite("t", &ph::t)
		.def_readwrite("can", &ph::can)
		.def_readwrite("n", &ph::n)
		.def_readwrite("x", &ph::x)
		.def_readwrite("y", &ph::y)
		.def_readwrite("frame", &ph::frame);
		
	py::class_<imChannel>(m, "imChannel")
		.def(py::init<std::vector<unsigned char>, int, int, long long, long long, int, std::vector<ph>>())
		.def_readwrite("can", &imChannel::can)
		.def_readwrite("tacmin", &imChannel::tacmin)
		.def_readwrite("tacmax", &imChannel::tacmax)
		.def_readwrite("tmin", &imChannel::tmin)
		.def_readwrite("tmax", &imChannel::tmax)
		.def_readwrite("line_id", &imChannel::line_id)
		.def_readwrite("phstream", &imChannel::phstream);
	return m.ptr();
}
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
	//practice examples
	
	py::module m("readPTUv2");
	m.def("inv", &inv);
	m.def("subtract", &subtract);
	
	py::class_<Pet>(m, "Pet")
		.def(py::init<const std::string &>())
		.def("setName", &Pet::setName)
		.def("getName", &Pet::getName)
		.def_readwrite("name", &Pet::name);
	
	m.def("ProcessPhotonStream", &ProcessPhotonStream_pywrap);
	//this way of binding passes the constructor, and allows for reading and writing the other variables.
	py::class_<ImOpts>(m, "ImOpts")
		.def(py::init<std::vector<int>, float, float, long long, float, int, int, int, int>())
		.def_readwrite("line_ids", &ImOpts::line_ids)
		.def_readwrite("dwelltime", &ImOpts::dwelltime)
		.def_readwrite("counttime", &ImOpts::counttime)
		.def_readwrite("NumRecords", &ImOpts::NumRecords)
		.def_readwrite("pxsize", &ImOpts::pxsize)
		.def_readwrite("dimX", &ImOpts::dimX)
		.def_readwrite("dimY", &ImOpts::dimY)
		.def_readwrite("ntacs", &ImOpts::ntacs)
		.def_readwrite("TAC_range", &ImOpts::TAC_range);
	/*
	py::class_<ph>(m, "ph")
		.def(py::init<int, long long, unsigned char, long long, float, float, int>())
		.def_readwrite("tac", &ph::tac)
		.def_readwrite("t", &ph::t)
		.def_readwrite("can", &ph::can)
		.def_readwrite("n", &ph::n)
		.def_readwrite("x", &ph::x)
		.def_readwrite("y", &ph::y)
		.def_readwrite("frame", &ph::frame);
	*/
	py::class_<ImChannel>(m, "ImChannel", py::buffer_protocol())
		.def(py::init<std::string, std::vector<unsigned char>, int, int, long long, long long, int, ImOpts>())
		.def_readwrite("name", &ImChannel::name)
		.def_readwrite("can", &ImChannel::can)
		.def_readwrite("tacmin", &ImChannel::tacmin)
		.def_readwrite("tacmax", &ImChannel::tacmax)
		.def_readwrite("tmin", &ImChannel::tmin)
		.def_readwrite("tmax", &ImChannel::tmax)
		.def_readwrite("line_id", &ImChannel::line_id)
		//this seems to work now in maintaining the unit8 data type. The only remaining problem is that the data is copied when transfer from cpp to python.
		//unfortunately there seems to be no direct way to solve this other then using Eigen or resorting to black magic. The former not supporting 3D arrays.
		// another suggestion was to use the def_buffer function to also change the array shape. I can't get it to work though, so I will just do that in python. link: https://alexsm.com/pybind11-buffer-protocol-opencv-to-numpy/
		.def("get_ltImage", [](ImChannel &imCh) -> py::array {
			//auto func = &ImChannel::get_ltImage;
			return py::array(imCh.ltImage.size(), imCh.ltImage.data());
		});
		// next time do:
		// 1 clean upper -- done
		// 2 test ProcessPhotonStream function
		// 3 implement in Image Manipulation class
	return m.ptr();
}
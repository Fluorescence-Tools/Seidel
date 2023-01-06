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
	
	//practice example
	py::class_<Pet>(m, "Pet")
		.def(py::init<const std::string &>())
		.def("setName", &Pet::setName)
		.def("getName", &Pet::getName)
		.def_readwrite("name", &Pet::name);
		
	//this way of binding passes the constructor, and allows for reading and writing the other variables.
	py::class_<ImOpts>(m, "ImOpts")
		.def(py::init<std::vector<int>, float, float, long long, float, float, int, int, int, int>())
		.def_readwrite("line_ids", &ImOpts::line_ids)
		.def_readwrite("dwelltime", &ImOpts::dwelltime)
		.def_readwrite("counttime", &ImOpts::counttime)
		.def_readwrite("NumRecords", &ImOpts::NumRecords)
		.def_readwrite("linestep", &ImOpts::linestep, 
		"e.g. 1 FRET line, 1 PIE line, 10nm px: linstep = 5")
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
		//.def_property_readonly("imOpts", &ImChannel::imOpts) #some compiler problem here, stating that there is no matching function call, investigate later
		//.def("get_ltImage", &ImChannel::get_ltImage)
		
		//this seems to work now in maintaining the unit8 data type. The only remaining problem is that the data is copied when transfer from cpp to python.
		//unfortunately there seems to be no direct way to solve this other then using Eigen or resorting to black magic. The former not supporting 3D arrays.
		.def("get_ltImage", [](ImChannel &imCh) -> py::array {
			//auto func = &ImChannel::get_ltImage;
			return py::array(imCh.ltImage.size(), imCh.ltImage.data());
		});
		//from example https://alexsm.com/pybind11-buffer-protocol-opencv-to-numpy/
		/*.def_buffer([](ImChannel & imCh) -> py::buffer_info{
			return py::buffer_info(
				//pointer to buffer // is this right?
				&imCh.ltImage,
				//size of one scalar
				sizeof(unsigned char),
				// Python struct-style format descriptor
				py::format_descriptor<unsigned char>::format(),
				//number of dimensions
				3,
				// Buffer dimensions
				{imCh.imOpts.dimY, imCh.imOpts.dimX, imCh.imOpts.ntacs},
				// Strides in bytes for each index
				{
					sizeof(unsigned char) * imCh.imOpts.dimX * imCh.imOpts.ntacs,
					sizeof(unsigned char) * imCh.imOpts.ntacs,
					sizeof(unsigned char)
				}
			);
			
		});*/
	return m.ptr();
}
//Author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//June 3, 2020
//TODO add description
//TODO change to readonly when applicable
//TODO add setters and getters for all setting variables

#include "ProcessPhotonStream.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <vector>


namespace py = pybind11;

PYBIND11_MODULE(imspy, m) {

	py::class_<imOpts>(m, "imOpts")
		.def(py::init<>())
		.def_readwrite("line_ids", &imOpts::line_ids)
		.def_readwrite("dwelltime", &imOpts::dwelltime)
		.def_readwrite("counttime", &imOpts::counttime)
		.def_readwrite("NumRecords", &imOpts::NumRecords)
		.def_readwrite("linestep", &imOpts::linestep, 
		"e.g. 1 FRET line, 1 PIE line, 10nm px: linstep = 5")
		.def_readwrite("pxsize", &imOpts::pxsize);
	
/*
	py::class_<ph>(m, "ph")
		.def(py::init<>())
		.def_readwrite("tac", &ph::tac)
		.def_readwrite("t", &ph::t)
		.def_readwrite("can", &ph::can)
		.def_readwrite("n", &ph::n)
		.def_readwrite("x", &ph::x)
		.def_readwrite("y", &ph::y)
		.def_readwrite("frame", &ph::frame);
*/

	py::class_<imChannel>(m, "imChannel", py::dynamic_attr())
		.def(py::init<>())
		.def("gentacdecay", &imChannel::gentacdecay)
		.def_readwrite("can", &imChannel::can)
		.def_readwrite("tacmin", &imChannel::tacmin)
		.def_readwrite("tacmax", &imChannel::tacmax)
		.def_readwrite("tmin", &imChannel::tmin)
		.def_readwrite("tmax", &imChannel::tmax)
		.def_readwrite("line_id", &imChannel::line_id)
		;

	py::class_<imspy>(m, "imspy")
		.def(py::init<
			std::vector<int>, 
			int64vec, 
			std::vector<unsigned char>,
			imOpts
		>())
		.def("ProcessPhotonStream", &imspy::ProcessPhotonStream)
		//.def("gettac", (Eigen::ArrayXi (imspy::*)()) &imspy::gettac, "get all tacs")
		.def("gettac", &imspy::gettac)
//		.def_readonly("tac", &imspy::tac)
//		.def_readonly("t", &imspy::t)
//		.def_readonly("can", &imspy::can)
//		.def_readonly("ImOpts", &imspy::ImOpts)
		;
	//return m.ptr();
}

//Author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//June 3, 2020

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "ProcessPhotonStream.h"

namespace py = pybind11;

int add(int i, int j) {
	return i + j;
}

ph make_ph(int val) {
	ph Ph;
	Ph.tac = val;
	return Ph;
}

imOpts make_imOpts(){
	imOpts ImOpts;
	return ImOpts;
}

int gettac(ph Ph) {
	return Ph.tac;
}

PYBIND11_PLUGIN(imspy) {
	py::module m("imspy");

	m.def("add", &add);
	m.def("make_ph", &make_ph);
	m.def("make_imOpts", &make_imOpts);
	m.def("gettac", &gettac);

	//to enable dynamic attributes, replace
	//py::class_<imOpts>(m, "imOpts", py::dynamic_attr())
	py::class_<imOpts>(m, "imOpts")
		.def_readwrite("line_ids", &imOpts::line_ids)
		.def_readwrite("dwelltime", &imOpts::dwelltime)
		.def_readwrite("counttime", &imOpts::counttime)
		.def_readwrite("NumRecords", &imOpts::NumRecords)
		.def_readwrite("linestep", &imOpts::linestep, 
		"e.g. 1 FRET line, 1 PIE line, 10nm px: linstep = 5")
		.def_readwrite("pxsize", &imOpts::pxsize);
	

	py::class_<ph>(m, "ph")
		.def_readwrite("tac", &ph::tac)
		.def_readwrite("t", &ph::t)
		.def_readwrite("can", &ph::can)
		.def_readwrite("n", &ph::n)
		.def_readwrite("x", &ph::x)
		.def_readwrite("y", &ph::y)
		.def_readwrite("frame", &ph::frame);

	return m.ptr();
}

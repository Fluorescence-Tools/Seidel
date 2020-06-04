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

typedef struct {
	int tac;
	std::vector<int> line_ids;
} ph;

ph make_ph(int val) {
	ph Ph;
	Ph.tac = val;
	Ph.line_ids.push_back(10);
	return Ph;
}

int gettac(ph Ph) {
	return Ph.tac;
}

std::vector<int> getlineids(ph Ph) {
	return Ph.line_ids;
}

PYBIND11_PLUGIN(imspy) {
	py::module m("imspy");

	m.def("add", &add);
	m.def("make_ph", &make_ph);
	m.def("gettac", &gettac);
	m.def("getlineids", &getlineids);

	py::class_<ph>(m, "ph")
		.def_readwrite("tac", &ph::tac)
		.def_readwrite("line_ids", &ph::line_ids);
	return m.ptr();
}

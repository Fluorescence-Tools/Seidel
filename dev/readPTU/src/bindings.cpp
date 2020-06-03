//Author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//June 3, 2020

#include "ProcessPhotonStream.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_PLUGIN(ProcessPhotonStream)
{
	py::module m("ProcessPhotonStream");
	m.def("ProcessPhotonStream", &ProcessPhotonStream);

	py::class_<imOpts>(m, "imOpts")
		.def("line_ids", &imOpts.line_ids)
		.def("dwelltime", &imOpts.dwelltime)
		.def("counttime", &imOpts.counttime)
		.def("NumRecords", &imOpts.NumRecords)
		.def("linestep", &imOpts.linstep, 
		"e.g. 1 FRET line, 1 PIE line, 10nm px: linstep = 5")
		.def("pxsize", &imOpts.pxsize);
	
	py::class_<ph>(m, "ph")
		.def("tac", &ph.tac)
		.def("t", &ph.t)
		.def("can", &ph.can)
		.def("n", &ph.n)
		.def("x", &ph.x)
		.def("y", &ph.y)
		.def("frame", &ph.frame);

	py::class_<imChannel>(m, "imChannel")
		.def("can", &imChannel.can)
		.def("tacmin", &imChannel.tacmin)
		.def("tacmax", &imChannel.tacmax)
		.def("tmin", &imChannel.tmin)
		.def("tmax", &imChannel.tmax)
		.def("line_id", &imChannel.line_id)
		.def("phstream", &imChannel.phstream);

	return m.ptr
}

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <meshtools.h>
namespace py = pybind11;

PYBIND11_MODULE(yams, m)
{
    using namespace yams;
    using T = double;

    m.def( "mesh_channel",
        py::overload_cast<const crv_vector<T> &, const std::vector<T> &, size_t, size_t, size_t>(
            &mesh_channel<T>
        ),
        "Mesh channel defined by a set of stream line curves",
        py::arg("crv_iso_eth"), py::arg("u"), py::arg("nu"), py::arg("nv"), py::arg("max_deg") = 3
    );
}
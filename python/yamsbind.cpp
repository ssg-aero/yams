#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <meshtools.h>
#include <gridreader.h>
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
        py::arg("iso_eth"), py::arg("ksi_i"), py::arg("n_iso_ksi"), py::arg("n_iso_eth"), py::arg("max_deg") = 3
    );

    // m.def( "read_vtk_grid",
    //     &read_vtk_grid,
    //     "Read Channel Grid containing meridional mesh and blades informations",
    //     py("fname")
    // );
}
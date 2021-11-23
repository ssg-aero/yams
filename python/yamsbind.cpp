#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <meshtools.h>
#include <gridreader.h>
#include <meridionalsolvercase.h>
#include <vtk_bind.h>
namespace py = pybind11;

PYBIND11_MODULE(yams, m)
{
    using namespace yams;
    using T = double;

    py::class_< MeridionalGridPoint<T> >(m, "MeridionalGridPoint")
    .def(py::init<>())
    .def_readwrite("x", &MeridionalGridPoint<T>::x)
    .def_readwrite("y", &MeridionalGridPoint<T>::y)
    .def_readwrite("Vm", &MeridionalGridPoint<T>::Vm)
    .def_readwrite("Vu", &MeridionalGridPoint<T>::Vu)
    .def_readwrite("rho", &MeridionalGridPoint<T>::rho)
    .def_readwrite("Ps", &MeridionalGridPoint<T>::Ps)
    .def_readwrite("Ts", &MeridionalGridPoint<T>::Ts)
    .def_readwrite("Pt", &MeridionalGridPoint<T>::Pt)
    .def_readwrite("Tt", &MeridionalGridPoint<T>::Tt)
    .def_readwrite("s", &MeridionalGridPoint<T>::s)
    ;

    py::class_< MeridionalGrid<T>, std::shared_ptr<MeridionalGrid<T>> >(m, "MeridionalGrid")
    .def(py::init<>())
    .def(py::init<size_t, size_t>())
    // .def(
    //     "__call__", 
    //     py::overload_cast<size_t, size_t>(&Array2d<T>::operator(), py::const_), 
    //     py::arg("i"), py::arg("j") 
    // )
    .def(
        "__call__", 
        py::overload_cast<size_t, size_t>(&MeridionalGrid<T>::operator()), 
        py::arg("i"), py::arg("j") 
    )
    .def(
        "nCols",
        &MeridionalGrid<T>::nCols
    )
    .def(
        "nRows",
        &MeridionalGrid<T>::nRows
    )
    .def(
        "resize",
        &MeridionalGrid<T>::resize
    )
    .def(
        "init",
        &MeridionalGrid<T>::init
    )
    ;

    py::class_< Grid2dMetricsPoint<T> >(m, "Grid2dMetricsPoint")
    .def(py::init<>())
    .def_readwrite("x1_ksi", &Grid2dMetricsPoint<T>::x1_ksi)
    .def_readwrite("x1_eth", &Grid2dMetricsPoint<T>::x1_eth)
    .def_readwrite("x2_ksi", &Grid2dMetricsPoint<T>::x2_ksi)
    .def_readwrite("x2_eth", &Grid2dMetricsPoint<T>::x2_eth)
    .def_readwrite("J", &Grid2dMetricsPoint<T>::J)
    ;

    py::class_< Grid2dMetrics<T>, std::shared_ptr<Grid2dMetrics<T>> >(m, "Grid2dMetrics")
    .def(py::init<>())
    ;

    py::class_< GridInfo<T>, std::shared_ptr<GridInfo<T>> >( m, "GridInfo" )
    .def(py::init<>())
    .def_readwrite("g", &GridInfo<T>::g)
    // .def_readwrite("g_metrics", &GridInfo<T>::g_metrics)
    // .def_readwrite("d_ksi", &GridInfo<T>::d_ksi)
    ;


    m.def( "mesh_channel",
        py::overload_cast<const crv_vector<T> &, const std::vector<T> &, size_t, size_t, size_t>(
            &mesh_channel<T>
        ),
        "Mesh channel defined by a set of stream line curves",
        py::arg("iso_eth"), py::arg("ksi_i"), py::arg("n_iso_ksi"), py::arg("n_iso_eth"), py::arg("max_deg") = 3
    );

    m.def( "read_vtk_grid",
        py::overload_cast<const std::string&>( &read_vtk_grid<T> ),
        "Read Channel Grid containing meridional mesh and blades informations",
        py::arg("fname")
    );

    m.def( "read_vtk_grid",
        py::overload_cast<const vtkSmartPointer<vtkStructuredGrid> &>( &read_vtk_grid<T> ),
        "Read Channel Grid containing meridional mesh and blades informations",
        py::arg("fname")
    );
}
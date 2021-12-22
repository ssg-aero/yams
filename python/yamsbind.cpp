#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <meshtools.h>
#include <gridreader.h>
#include <meridionalsolvercase.h>
#include <gridmetrics.h>
#include <eqcurvaturesolver.h>
#include <gridrender.h>
#include <curvaturesolver.h>
#include <plots.h>
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
    .def_readwrite("m", &MeridionalGridPoint<T>::m)
    .def_readwrite("l", &MeridionalGridPoint<T>::l)
    .def_readwrite("cur", &MeridionalGridPoint<T>::cur)

    .def_readwrite("Vm", &MeridionalGridPoint<T>::Vm)
    .def_readwrite("Vu", &MeridionalGridPoint<T>::Vu)
    .def_readwrite("rho", &MeridionalGridPoint<T>::rho)
    .def_readwrite("Ps", &MeridionalGridPoint<T>::Ps)
    .def_readwrite("Ts", &MeridionalGridPoint<T>::Ts)
    .def_readwrite("Pt", &MeridionalGridPoint<T>::Pt)
    .def_readwrite("Tt", &MeridionalGridPoint<T>::Tt)
    .def_readwrite("s", &MeridionalGridPoint<T>::s)
    .def_readwrite("iB", &MeridionalGridPoint<T>::iB)
    .def_readwrite("k", &MeridionalGridPoint<T>::k)
    .def_readwrite("bet", &MeridionalGridPoint<T>::bet)
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
    .def(py::init<size_t, size_t>())
    ;

    py::class_< GridInfo<T>, std::shared_ptr<GridInfo<T>> >( m, "GridInfo" )
    .def(py::init<>())
    .def_readwrite("g", &GridInfo<T>::g)
    .def_readwrite("g_metrics", &GridInfo<T>::g_metrics)
    .def_readwrite("rho_cst", &GridInfo<T>::rho_cst)
    .def_readwrite("RF", &GridInfo<T>::RF)
    .def_readwrite("Pref", &GridInfo<T>::Pref)
    .def_readwrite("Tref", &GridInfo<T>::Tref)
    .def_readwrite("tol_newtow_mf_f", &GridInfo<T>::tol_newtow_mf_f)
    .def_readwrite("tol_newtow_mf_u", &GridInfo<T>::tol_newtow_mf_u)
    .def_readwrite("vm_distribution_max_count", &GridInfo<T>::vm_distribution_max_count)
    ;

    py::enum_<MeridionalBladeMode>(m, "MeridionalBladeMode", py::arithmetic())
    .value("DESIGN_BETA_OUT", MeridionalBladeMode::DESIGN_BETA_OUT)
    .value("DESIGN_PSI", MeridionalBladeMode::DESIGN_PSI)
    .value("DIRECT", MeridionalBladeMode::DIRECT)
    ;

    py::class_<BladeInfo<T>>(m,"BladeInfo")
    .def(py::init<>())
    .def_readwrite("name",&BladeInfo<T>::name)
    .def_readwrite("i1",&BladeInfo<T>::i1)
    .def_readwrite("i_s",&BladeInfo<T>::is)
    .def_readwrite("i2",&BladeInfo<T>::i2)
    .def_readwrite("omg",&BladeInfo<T>::omg)
    .def_readwrite("omg_",&BladeInfo<T>::omg_)
    .def_readwrite("mode",&BladeInfo<T>::mode)
    .def_readwrite("beta_out",&BladeInfo<T>::beta_out)
    .def_readwrite("psi",&BladeInfo<T>::psi)
    ;

    py::enum_<MeridionalBC>(m, "MeridionalBC", py::arithmetic())
    .value("INLET_Mf_Ts_Ps_Vu", MeridionalBC::INLET_Mf_Ts_Ps_Vu)
    .value("INLET_VmMoy_Ts_Ps_Vu", MeridionalBC::INLET_VmMoy_Ts_Ps_Vu)
    ;

    py::class_<InletBC<T>>(m,"InletBC")
    .def(py::init<>())
    .def_readwrite("mode",&InletBC<T>::mode)
    .def_readwrite("Ps",&InletBC<T>::Ps)
    .def_readwrite("Ts",&InletBC<T>::Ts)
    .def_readwrite("Vu",&InletBC<T>::Vu)
    .def_readwrite("Pt",&InletBC<T>::Pt)
    .def_readwrite("Tt",&InletBC<T>::Tt)
    .def_readwrite("Mf",&InletBC<T>::Mf)
    .def_readwrite("Vm_moy",&InletBC<T>::Vm_moy)
    ;

    py::class_<SolverLog<T>>(m,"SolverLog")
    .def(py::init<>())
    .def_readonly("delta_pos_max",&SolverLog<T>::delta_pos_max)
    .def_readonly("delta_pos_moy",&SolverLog<T>::delta_pos_moy)
    .def_readonly("delta_pos",&SolverLog<T>::delta_pos)
    .def("clear",&SolverLog<T>::clear)
    ;

    py::class_<SolverCase<T>>(m,"SolverCase")
    .def(py::init<>())
    .def_readwrite("gi",&SolverCase<T>::gi)
    .def_readwrite("bld_info_lst",&SolverCase<T>::bld_info_lst)
    .def_readwrite("inlet",&SolverCase<T>::inlet)
    .def_readwrite("mf",&SolverCase<T>::mf)
    .def_readwrite("log",&SolverCase<T>::log)
    .def_readwrite("max_geom",&SolverCase<T>::max_geom)
    .def_readwrite("eps",&SolverCase<T>::eps)
    .def_readwrite("tol_rel_mf",&SolverCase<T>::tol_rel_mf)
    .def_readwrite("tol_rel_pos",&SolverCase<T>::tol_rel_pos)
    .def_readwrite("relocate",&SolverCase<T>::relocate)
    .def("__copy__",  [](const  SolverCase<T> &self) {
        return  SolverCase<T>(self);
    })
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

    m.def( "make_vtk_grid",
        py::overload_cast<const MeridionalGrid<T> &>(&make_vtk_grid<T>),
        "Convert grid to vtk strutured grid",
        py::arg("g")
    );

    m.def( "make_grid_info",
        py::overload_cast<vtkStructuredGrid*>(&make_grid_info<T>),
        "Make and init GridInfo",
        py::arg( "sgrid" )
    );

    m.def( "make_grid_info",
        py::overload_cast<const std::string&>(&make_grid_info<T>),
        "Make and init GridInfo metrics",
        py::arg( "fname" )
    );

    m.def("make_solver_case",
        py::overload_cast<
            vtkStructuredGrid* , 
            const std::vector< std::tuple< BladeInfo<T> , gbs::BSSfunction<T> > > &
        >(&make_solver_case<T,gbs::BSSfunction<T>>),
        "Make solver case and init GridInfo metrics",
        py::arg("sgrid"), py::arg("bld_info_lst")
    );

    m.def("make_solver_case",
        py::overload_cast<
            vtkStructuredGrid* , 
            const std::vector< std::tuple< BladeInfo<T> , std::function<T(T,T)> > > &
        >(&make_solver_case<T,std::function<T(T,T)>>),
        "Make solver case and init GridInfo metrics",
        py::arg("sgrid"), py::arg("bld_info_lst")
    );



    m.def("curvature_solver",
        &curvature_solver<T>,
        "Solve case with curvature solver",
        py::arg( "solver_case" )
    );

    m.def( "plot",
        [](const MeridionalGrid<T> &g, const std::string &value , bool edges_on, bool countour_on)
        {
            plot_vtkStructuredGrid(make_vtk_grid<T>(g), value.c_str(), edges_on, countour_on);
        },
        "Plot grid results",
        py::arg("g"), py::arg("value") = "Vm", py::arg("edges_on") = false, py::arg("countour_on") = false
    );

    m.def( "plot",
        [](const SolverCase<T> &solver_case, const std::string &value , bool edges_on, bool countour_on)
        {
            plot_vtkStructuredGrid(solver_case, value.c_str(), edges_on, countour_on);
        },
        "Plot grid results",
        py::arg("solver_case"), py::arg("value") = "Vm", py::arg("edges_on") = false, py::arg("countour_on") = false
    );

    m.def( "plot",
        [](vtkStructuredGrid* sgrid, const std::string &value , bool edges_on, bool countour_on)
        {
            plot_vtkStructuredGrid(sgrid, value.c_str(), edges_on, countour_on);
        },
        "Plot grid results",
        py::arg("g"), py::arg("value") = "Vm", py::arg("edges_on") = false, py::arg("countour_on") = false
    );
    m.def( "plot_residual",&plot_residual<T>,py::arg("log"));

}
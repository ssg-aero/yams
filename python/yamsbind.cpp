#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;

#include <meshtools.h>
#include <gridreader.h>
#include <meridionalsolvercase.h>
#include <solvercasetools.h>
#include <gridmetrics.h>
#include <eqcurvaturesolver.h>
#include <gridrender.h>
#include <curvaturesolver.h>
#include <plots.h>
#include <vtk_bind.h>
#include <gbs-render/vtkgridrender.h>

#include <baldeToBlade/bladeToBladeCurvatureSolver.h>
#include <baldeToBlade/gridReader.h>
#include <baldeToBlade/bladeToBladeCurvatureSolverPost.h>

// #include <meridionalsolvercase_experiemental.h>
// #include <curvaturesolver_experimental.h>

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
    .def_readwrite("Cp", &MeridionalGridPoint<T>::Cp)
    .def_readwrite("ga", &MeridionalGridPoint<T>::ga)
    .def_readwrite("s", &MeridionalGridPoint<T>::s)
    .def_readwrite("iB", &MeridionalGridPoint<T>::iB)
    .def_readwrite("k", &MeridionalGridPoint<T>::k)
    .def_readwrite("th_", &MeridionalGridPoint<T>::th_)
    .def_readwrite("bet", &MeridionalGridPoint<T>::bet)
    .def_readwrite("omg", &MeridionalGridPoint<T>::omg)
    .def_readwrite("omg_", &MeridionalGridPoint<T>::omg_)
    .def_readwrite("gam", &MeridionalGridPoint<T>::gam)
    .def_readwrite("phi", &MeridionalGridPoint<T>::phi)
    .def_readwrite("H", &MeridionalGridPoint<T>::H)
    .def_readwrite("I", &MeridionalGridPoint<T>::I)
    .def_readwrite("s", &MeridionalGridPoint<T>::s)
    .def_readwrite("q", &MeridionalGridPoint<T>::q)
    ;

    py::class_< MeridionalGrid<T>, std::shared_ptr<MeridionalGrid<T>> >(m, "MeridionalGrid")
    .def(py::init<>())
    .def(py::init<size_t, size_t>())
    .def(
        "__call__", 
        py::overload_cast<size_t, size_t>(&MeridionalGrid<T>::operator()), 
        py::arg("i"), py::arg("j"),
        py::return_value_policy::reference
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
    .def_readwrite("R", &GridInfo<T>::R)
    .def_readwrite("RF", &GridInfo<T>::RF)
    .def_readwrite("Pref", &GridInfo<T>::Pref)
    .def_readwrite("Tref", &GridInfo<T>::Tref)
    .def_readwrite("tol_newtow_mf_f", &GridInfo<T>::tol_newtow_mf_f)
    .def_readwrite("tol_newtow_mf_u", &GridInfo<T>::tol_newtow_mf_u)
    .def_readwrite("vm_distribution_max_count", &GridInfo<T>::vm_distribution_max_count)
    .def_readwrite("ni", &GridInfo<T>::ni)
    .def_readwrite("nj", &GridInfo<T>::nj)
    .def_readwrite("d_ksi", &GridInfo<T>::d_ksi)
    .def_readwrite("d_eth", &GridInfo<T>::d_eth)
    .def_readwrite("j_0", &GridInfo<T>::j_0)
    .def_readwrite("KD",&GridInfo<T>::KD)
    ;

    py::enum_<MeridionalBladeMode>(m, "MeridionalBladeMode", py::arithmetic())
    .value("DESIGN_BETA_OUT", MeridionalBladeMode::DESIGN_BETA_OUT)
    .value("DESIGN_ALPHA_OUT", MeridionalBladeMode::DESIGN_ALPHA_OUT)
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
    .def_readwrite("Ksh",&BladeInfo<T>::Ksh)
    .def_readwrite("z_",&BladeInfo<T>::z_)
    .def_readwrite("omg_",&BladeInfo<T>::omg_)
    .def_readwrite("dev",&BladeInfo<T>::dev)
    .def_readwrite("mu",&BladeInfo<T>::mu)
    .def_readwrite("mode",&BladeInfo<T>::mode)
    .def_readwrite("beta_out",&BladeInfo<T>::beta_out)
    .def_readwrite("alpha_out",&BladeInfo<T>::alpha_out)
    .def_readwrite("psi",&BladeInfo<T>::psi)
    .def_readwrite("k",&BladeInfo<T>::k)
    .def_readwrite("tb",&BladeInfo<T>::tb)
    .def_readwrite("eps",&BladeInfo<T>::eps)
    .def_readwrite("gauge",&BladeInfo<T>::gauge)
    .def_readwrite("max_thickness",&BladeInfo<T>::max_thickness)
    .def_readwrite("compute_dev",&BladeInfo<T>::compute_dev)
    ;

    py::enum_<MeridionalBC>(m, "MeridionalBC", py::arithmetic())
    .value("INLET_Mf_Tt_Pt_Vu", MeridionalBC::INLET_Mf_Tt_Pt_Vu)
    .value("INLET_Mf_Ts_Ps_Vu", MeridionalBC::INLET_Mf_Ts_Ps_Vu)
    .value("INLET_VmMoy_Ts_Ps_Vu", MeridionalBC::INLET_VmMoy_Ts_Ps_Vu)
    .value("INLET_Vm_Ts_Ps_Vu", MeridionalBC::INLET_Vm_Ts_Ps_Vu)
    ;

    py::class_<InletBC<T>>(m,"InletBC")
    .def(py::init<>())
    .def_readwrite("mode",&InletBC<T>::mode)
    .def_readwrite("Ps",&InletBC<T>::Ps)
    .def_readwrite("Ts",&InletBC<T>::Ts)
    .def_readwrite("Vm",&InletBC<T>::Vm)
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
    .def_readwrite("mf_ref_span",&SolverCase<T>::mf_ref_span)
    .def_readwrite("mf_uniform",&SolverCase<T>::mf_uniform)
    .def_readwrite("use_meridional_grad",&SolverCase<T>::use_meridional_grad)
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

    m.def( "mesh_channel",
        py::overload_cast<
            const crv_vector<T> &, 
            const std::vector<gbs::points_vector<T,2>> &,
            // const std::vector<std::vector<std::array<T,2>>> &,
            size_t, 
            size_t, 
            T,
            size_t>
        (&mesh_channel<T>),
        "Mesh channel defined by a set of stream line curves and hard points",
        py::arg("master_streams"), 
        py::arg("hard_points"), 
        py::arg("n_streams"), 
        py::arg("n_spans"), 
        py::arg("tol"), 
        py::arg("max_deg") = 3
    );

    m.def("smooth_mesh",
        py::overload_cast<
            gbs::points_vector<T,2> &,
            size_t,
            const std::vector<size_t> &,
            size_t,
            T>(&smooth_mesh<T>),
            "Smooth mesh blocs",
            py::arg("pts"),
            py::arg("n_streams"),
            py::arg("n_span_per_bloc"),
            py::arg("max_it") = 100,
            py::arg("tol") = 1.e-5
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
            const std::vector< BladeInfo<T> > &
        >(&make_solver_case<T>),
        "Make solver case and init GridInfo metrics",
        py::arg("sgrid"), py::arg("bld_info_lst")
    );

    m.def("make_solver_case",
        [](gbs::points_vector<T,2> &pts, size_t  n_span,size_t n_stream, const std::vector< BladeInfo<T> > &bld_info_lst)
        {
            auto sgrid = gbs::make_structuredgrid(pts, n_span, n_stream);
            return make_solver_case(sgrid, bld_info_lst);
        },
        "Make solver case and init GridInfo metrics",
        py::arg("pts"), py::arg("n_span"), py::arg("n_stream"), py::arg("bld_info_lst")
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
        py::overload_cast<SolverCase<T> &>( &curvature_solver<T> ),
        "Solve case with curvature solver",
        py::arg( "solver_case" )
    );
    // m.def("curvature_solver",
    //     py::overload_cast<SolverCaseSet<T> &>( &curvature_solver<T> ),
    //     "Solve case with curvature solver",
    //     py::arg( "set" )
    // );
    m.def("apply_bc", py::overload_cast<SolverCase<T> &>( &apply_bc<T> ), "Apply the Boundary conditions", py::arg( "solver_case" ));
    m.def("apply_mf", py::overload_cast<SolverCase<T> &>( &apply_mf<T> ), "Apply MassFlow", py::arg( "solver_case" ));
    m.def("init_values",
        py::overload_cast<SolverCase<T> &, T, T>(&init_values<T>),
        "Mean Line Computation",
        py::arg( "solver_case" ), py::arg("tol_rel_mf"), py::arg("eps")
    );

    m.def("eq_vu",
        [](const MeridionalGrid<T> &g, const Grid2dMetrics<T> &g_metrics, size_t i, size_t j, T d_ksi, T d_eth)
        {
            return eq_vu(g,g_metrics, i,j,d_ksi,d_eth);
        }
    );
    m.def("G",
        [](const MeridionalGridPoint<T> &gp)
        {
            return G(gp);
        }
    );
    m.def("eq_bet",
        [](const MeridionalGrid<T> &g, const Grid2dMetrics<T> &g_metrics, size_t i, size_t j, T d_ksi, T d_eth)
        {
            return eq_bet(g,g_metrics, i,j,d_ksi,d_eth);
        }
    );
    m.def("D",
        [](const MeridionalGrid<T> &g, const Grid2dMetrics<T> &g_metrics, size_t i, size_t j, T d_ksi, T d_eth)
        {
            return D(g,g_metrics, i,j,d_ksi,d_eth);
        }
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

    m.def( "switch_to_direct",&switch_to_direct<T>,py::arg("solver_case"));

    py::class_< BladeToBladeCurvatureSolverData<T>, std::shared_ptr<BladeToBladeCurvatureSolverData<T>> >(m, "BladeToBladeCurvatureSolverData")
    .def(py::init<>())
    .def(py::init<size_t, size_t>())
    .def_readonly("ni", &BladeToBladeCurvatureSolverData<T>::ni,"Stream line number")
    .def_readonly("nj", &BladeToBladeCurvatureSolverData<T>::nj,"Computation plane number")
    .def_readonly("W", &BladeToBladeCurvatureSolverData<T>::W,"Relative velocity")
    .def_readonly("V", &BladeToBladeCurvatureSolverData<T>::V,"Absolute velocity")
    .def_readonly("Vm", &BladeToBladeCurvatureSolverData<T>::Vm,"Meridional velocity")
    .def_readonly("Vu", &BladeToBladeCurvatureSolverData<T>::Vu,"Absolute tangential velocity")
    .def_readonly("Wu", &BladeToBladeCurvatureSolverData<T>::Wu,"Relative tangential velocity")
    .def_readonly("M", &BladeToBladeCurvatureSolverData<T>::M,"Absolute Mach number")
    .def_readonly("MW", &BladeToBladeCurvatureSolverData<T>::MW,"Relative Mach number")
    .def_readonly("R", &BladeToBladeCurvatureSolverData<T>::R,"Density")
    .def_readonly("PS", &BladeToBladeCurvatureSolverData<T>::PS,"Static Pressure")
    .def_readonly("TS", &BladeToBladeCurvatureSolverData<T>::TS,"Static Temperature")
    .def_readonly("PT", &BladeToBladeCurvatureSolverData<T>::PT,"Total Pressure")
    .def_readonly("TT", &BladeToBladeCurvatureSolverData<T>::TT,"Total Temperature")
    .def_readonly("H", &BladeToBladeCurvatureSolverData<T>::H,"Enthalpy")
    .def_readonly("S", &BladeToBladeCurvatureSolverData<T>::S,"Entropy")
    .def_readonly("Q", &BladeToBladeCurvatureSolverData<T>::Q,"Accumulated mass flow")
    ;

    py::class_< BladeToBladeCurvatureSolverMesh<T>, std::shared_ptr<BladeToBladeCurvatureSolverMesh<T>> >(m, "BladeToBladeCurvatureSolverMesh")
    .def(py::init<>())
    .def(py::init<size_t, size_t>())
    .def_readonly("ni", &BladeToBladeCurvatureSolverMesh<T>::ni,"Stream line number")
    .def_readonly("nj", &BladeToBladeCurvatureSolverMesh<T>::nj,"Computation plane number")
    .def_readonly("M", &BladeToBladeCurvatureSolverMesh<T>::M,"Curvilinear abscissa on stream line")
    .def_readonly("U", &BladeToBladeCurvatureSolverMesh<T>::U,"Stream line parameter")
    .def_readonly("TH", &BladeToBladeCurvatureSolverMesh<T>::TH,"Tangential coordinate")
    .def_readonly("PH", &BladeToBladeCurvatureSolverMesh<T>::PH,"Stream slope")
    .def_readonly("R", &BladeToBladeCurvatureSolverMesh<T>::R,"Radial position on stream line")
    .def_readonly("BT", &BladeToBladeCurvatureSolverMesh<T>::BT,"Relative speed flow angle")
    .def_readonly("TAU", &BladeToBladeCurvatureSolverMesh<T>::TAU,"Stream tube thickness")
    ;

    py::class_< BladeToBladeCurvatureSolver<T>, std::shared_ptr<BladeToBladeCurvatureSolver<T>> >(m, "BladeToBladeCurvatureSolver")
    .def(py::init<>())
    .def(
        py::init<
            const gbs::points_vector<T,2> &,
            size_t,
            const std::shared_ptr<gbs::Curve<T, 2>> &,
            size_t,
            size_t,
            size_t
        >(),
        py::arg("points"),
        py::arg("n_computation_planes"),
        py::arg("stream_line"),
        py::arg("n_blades"),
        py::arg("j_le"),
        py::arg("j_te")
    )
    .def_property(
        "compressibility_rel_tol",
        &BladeToBladeCurvatureSolver<T>::compressibilityRelTol,
        &BladeToBladeCurvatureSolver<T>::setCompressibilityRelTol
    )
    .def_property(
        "geom_rel_tol",
        &BladeToBladeCurvatureSolver<T>::geomRelTol,
        &BladeToBladeCurvatureSolver<T>::setGeomRelTol
    )
    .def_property(
        "mass_flow_rel_tol",
        &BladeToBladeCurvatureSolver<T>::massFlowRelTol,
        &BladeToBladeCurvatureSolver<T>::setMassFlowRelTol
    )
    .def_property(
        "relax_factor_geom",
        &BladeToBladeCurvatureSolver<T>::relaxFactorGeom,
        &BladeToBladeCurvatureSolver<T>::setRelaxFactorGeom
    )
    .def_property(
        "relax_factor_per",
        &BladeToBladeCurvatureSolver<T>::relaxFactorPeriodic,
        &BladeToBladeCurvatureSolver<T>::setRelaxFactorPeriodic
    )
    .def_property(
        "relax_factor_compress",
        &BladeToBladeCurvatureSolver<T>::relaxFactorCompressibility,
        &BladeToBladeCurvatureSolver<T>::setRelaxFactorCompressibility
    )
    .def_property(
        "periodicity",
        &BladeToBladeCurvatureSolver<T>::periodicity,
        &BladeToBladeCurvatureSolver<T>::setPeriodicity
    )
    .def_property(
        "full_passage_mass_flow",
        &BladeToBladeCurvatureSolver<T>::fullBladePassageMassFlow,
        &BladeToBladeCurvatureSolver<T>::setFullBladePassageMassFlow
    )
    .def_property(
        "rotation_speed",
        &BladeToBladeCurvatureSolver<T>::rotationSpeed,
        &BladeToBladeCurvatureSolver<T>::setRotationSpeed
    )
    // .def_property(
    //     "jLe",
    //     &BladeToBladeCurvatureSolver<T>::leadingEdgeIndex,
    //     &BladeToBladeCurvatureSolver<T>::setLeadingEdgeIndex
    // )
    // .def_property(
    //     "jTe",
    //     &BladeToBladeCurvatureSolver<T>::trailingEdgeIndex,
    //     &BladeToBladeCurvatureSolver<T>::setTrailingEdgeIndex
    // )
    .def_property(
        "fix_upstream_flow_periodicity",
        &BladeToBladeCurvatureSolver<T>::fixUpStreamFlowPeriodicity,
        &BladeToBladeCurvatureSolver<T>::setFixUpStreamFlowPeriodicity
    )
    .def_property(
        "fix_downstream_flow_periodicity",
        &BladeToBladeCurvatureSolver<T>::fixDownStreamFlowPeriodicity,
        &BladeToBladeCurvatureSolver<T>::setFixDownStreamFlowPeriodicity
    )
    .def_property(
        "max_convergence_iterations",
        &BladeToBladeCurvatureSolver<T>::maxConvergenceIterations,
        &BladeToBladeCurvatureSolver<T>::setMaxConvergenceIterations
    )
    .def_property_readonly(
        "jLe", &BladeToBladeCurvatureSolver<T>::leadingEdgeIndex )
    .def_property_readonly(
        "jTe", &BladeToBladeCurvatureSolver<T>::trailingEdgeIndex )
    .def_property_readonly(
        "dimensions",&BladeToBladeCurvatureSolver<T>::dimensions
    )
    .def_property_readonly(
        "mesh",&BladeToBladeCurvatureSolver<T>::mesh
    )
    .def_property_readonly(
        "data",&BladeToBladeCurvatureSolver<T>::data
    )
    .def_property_readonly(
        "meridionalStreamLine", &BladeToBladeCurvatureSolver<T>::meridionalStreamLine
    )
    .def_property_readonly(
        "compressibilityResidual",&BladeToBladeCurvatureSolver<T>::compressibilityResidualAverage
    )
    .def_property_readonly(
        "periodicityUpStreamResidual",&BladeToBladeCurvatureSolver<T>::periodicityUpStreamResidual
    )
    .def_property_readonly(
        "periodicityDownStreamResidual",&BladeToBladeCurvatureSolver<T>::periodicityDownStreamResidual
    )
    .def_property_readonly(
        "geomResidualMax",&BladeToBladeCurvatureSolver<T>::geomResidualMax
    )
    .def("setPtIn",&BladeToBladeCurvatureSolver<T>::setPtIn)
    .def("setTtIn",&BladeToBladeCurvatureSolver<T>::setTtIn)
    .def("setDr",&BladeToBladeCurvatureSolver<T>::setDr)
    .def(
        "computeW", 
        py::overload_cast<>(&BladeToBladeCurvatureSolver<T>::computeW)
    )
    ;

    m.def(
        "plot",
        py::overload_cast<
            const BladeToBladeCurvatureSolverMesh<T> &,
            const vector<T> &,
            const char *,
            bool,
            bool,
            T
        > ( &plot<T> ),
        "Plot Blade to blade flow",
        py::arg("msh"),
        py::arg("value"),
        py::arg("name"),
        py::arg("edges_on")=true,
        py::arg("contour_on")=true,
        py::arg("scale")=1.
    );

    m.def(
        "getPsSide1",&getPsSide1<T>,
        "Get static pressure on blade side 1",
        py::arg("solver")
    );

    m.def(
        "getPsSide2",&getPsSide2<T>,
        "Get static pressure on blade side 2",
        py::arg("solver")
    );

    m.def(
        "getValuesSide1",&getValuesSide1<T>,
        "Get values on blade side 1",
        py::arg("solver"),
        py::arg("val")
    );

    m.def(
        "getValuesSide2",&getValuesSide2<T>,
        "Get values on blade side 2",
        py::arg("solver"),
        py::arg("val")
    );

    m.def(
        "getXYZSide1",&getXYZSide1<T>,
        "Get blade side 1 cut 3d coordinates",
        py::arg("solver")
    );

    m.def(
        "getXYZSide2",&getXYZSide2<T>,
        "Get blade side 2 cut 3d coordinates",
        py::arg("solver")
    );
}


#pragma once
#include <datastorage.h>
#include <meridionalsolvercase.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridReader.h>
namespace yams
{
    template <typename T>
    inline auto distance(const MeridionalGridPoint<T> &gp1, const MeridionalGridPoint<T> &gp2) -> T
    {
        return sqrt((gp1.y - gp2.y) * (gp1.y - gp2.y) + (gp1.x - gp2.x) * (gp1.x - gp2.x));
    }

    template <typename T>
    inline auto compute_abscissas(MeridionalGrid<T> &g)
    {
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        for (auto i = 0; i < ni; i++)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).m = i == 0 ? 0. : g(i - 1, j).m + distance(g(i, j), g(i - 1, j));
                g(i, j).l = j == 0 ? 0. : g(i, j - 1).l + distance(g(i, j), g(i, j - 1));
            }
        }
    }


    inline auto fz = [](const auto &gp) { return gp.x; };

    inline auto fr = [](const auto &gp) { return gp.y; };

    inline auto fm = [](const auto &gp) { return gp.m; };
    
    inline auto fl = [](const auto &gp) { return gp.l; };
    
    inline auto fphi = [](const auto &gp) { return gp.phi; };

    template <typename T>
    inline auto compute_angles(MeridionalGrid<T> &g)
    {
        size_t ni = g.nRows();
        size_t nj = g.nCols();

        T drqdm, dzqdm, drqdl, dzqdl;
        for (auto i = 0; i < ni; i++)
        {
            for (auto j = 0; j < nj; j++)
            {
                // auto gam_ = -PI / 2 + PI * j / (nj - 1.);
                drqdm = D1_O2_i(g, i, j, fr, fm);
                dzqdm = D1_O2_i(g, i, j, fz, fm);
                g(i,j).phi = atan2(drqdm, dzqdm); // Stream line angle
                drqdl = D1_O2_j(g, i, j, fr, fl);
                dzqdl = D1_O2_j(g, i, j, fz, fl);
                g(i,j).gam = atan2(dzqdl, drqdl); // Span line angle
            }
        }
    }
    
    template <typename T>
    inline auto compute_curvature(MeridionalGrid<T> &g, bool interpolate = false)
    {
        size_t ni = g.nRows()-1;
        size_t nj = g.nCols();

        for (auto i = 1; i < ni; i++)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i,j).cur = D1_O2_i(g, i, j, fphi, fm);// TODO check why Aungier put -DphiDm
            }
        }
        for (auto j = 0; j < nj; j++)
        {
            if( interpolate )
            {
                g(0, j).cur = 2. * g(1, j).cur - g(2, j).cur;
                g(ni, j).cur = 2. * g(ni - 1, j).cur - g(ni - 2, j).cur;
            }
            else
            {
                g(0, j).cur  =  D1_O1_i_fw(g, 0 , j, fphi, fm);
                g(ni, j).cur =  D1_O1_i_bw(g, ni, j, fphi, fm);
            }

        }
    }

    template <typename T>
    inline auto compute_geom_values(MeridionalGrid<T> &g)
    {
        compute_abscissas(g);
        compute_angles(g);
        compute_curvature(g);
    }

    // template <typename T>
    // inline auto compute_angles(MeridionalGrid<T> &g,const Array2d<Grid2dMetricsPoint<T>> &g_metrics)
    // {
    //     size_t ni = g.nRows();
    //     size_t nj = g.nCols();
    //     T d_ksi = 1. / (ni - 1.);
    //     T d_eth = 1. / (nj - 1.);

    //     T drqdm, dzqdm, drqdl, dzqdl;
    //     for (auto i = 0; i < ni; i++)
    //     {
    //         for (auto j = 0; j < nj; j++)
    //         {
    //             // auto gam_ = -PI / 2 + PI * j / (nj - 1.);
    //             // drqdm = D1_O2_i(g, i, j, fr, fm);
    //             // dzqdm = D1_O2_i(g, i, j, fz, fm);
    //             drqdm = D1_O2_dx1(g,g_metrics, i, j, d_ksi, d_eth, fr);
    //             dzqdm = D1_O2_dx1(g,g_metrics, i, j, d_ksi, d_eth, fz);
    //             g(i,j).phi = atan2(drqdm, dzqdm); // Stream line angle
    //             // drqdl = D1_O2_j(g, i, j, fr, fl);
    //             // dzqdl = D1_O2_j(g, i, j, fz, fl);
    //             drqdl = D1_O2_dx2(g,g_metrics, i, j, d_ksi, d_eth, fr);
    //             dzqdl = D1_O2_dx2(g,g_metrics, i, j, d_ksi, d_eth, fz);
    //             g(i,j).gam = atan2(dzqdl, drqdl); // Span line angle
    //             g(i,j).cgp = std::cos(g(i,j).gam+g(i,j).phi);
    //             g(i,j).sgp = std::sin(g(i,j).gam+g(i,j).phi);
    //         }
    //     }
    // }
    
    // template <typename T>
    // inline auto compute_curvature(MeridionalGrid<T> &g,const Array2d<Grid2dMetricsPoint<T>> &g_metrics)
    // {
    //     size_t ni = g.nRows();
    //     size_t nj = g.nCols();
    //     T d_ksi = 1. / (ni - 1.);
    //     T d_eth = 1. / (nj - 1.);

    //     for (auto i = 1; i < ni-1; i++)
    //     {
    //         for (auto j = 0; j < nj; j++)
    //         {
    //             // g(i,j).cur = D1_O2_i(g, i, j, fphi, fm);// TODO check why Aungier put -DphiDm
    //             g(i,j).cur = D1_O2_dx1(g,g_metrics, i, j, d_ksi, d_eth, fphi);// TODO check why Aungier put -DphiDm
    //         }
    //     }
    //     for (auto j = 0; j < nj; j++) // extrapolate on bounds
    //     {
    //         g(0, j).cur = 2. * g(1, j).cur - g(2, j).cur;
    //         g(ni-1, j).cur = 2. * g(ni - 2, j).cur - g(ni - 3, j).cur;
    //     }
    // }

    template <typename T>
    inline auto compute_grid_metrics(MeridionalGrid<T> &g, Array2d<Grid2dMetricsPoint<T>> &g_metrics, const auto &f_m, const auto &f_l)
    {
        compute_abscissas(g);
        compute_metrics(g,f_m,f_l,g_metrics);
        compute_angles(g);
        compute_curvature(g);
        // compute_angles(g,g_metrics); <- problematics Novak 1977 indicate derivatte has to be made along stream lines hence a priori more correct without metrics
        // compute_curvature(g,g_metrics);
    }

    template <typename T>
    auto make_grid_info(vtkStructuredGrid* sgrid)
    {
        auto dims =sgrid->GetDimensions();
        size_t ni = dims[0];
        size_t nj = dims[1];
        double ksi = 1. / (ni-1.);
        double eth = 1. / (nj-1.);

        auto g = read_vtk_grid<T>(sgrid);
        auto g_metrics = Grid2dMetrics<T>{ni,nj};
        compute_grid_metrics(g,g_metrics,fm,fl);

        auto gi = GridInfo<T>{
                .g = std::make_shared<  MeridionalGrid<T> >(g),
                .g_metrics = std::make_shared< Grid2dMetrics<T> >( g_metrics ),
                .d_ksi = ksi,
                .d_eth = eth,
                .ni = ni,
                .nj = nj,
            };

        return std::make_shared<GridInfo<T>>( gi );
    }

    template <typename T>
    auto make_grid_info(const std::string &fname)
    {
        vtkNew<vtkXMLStructuredGridReader> reader;
        reader->SetFileName(fname.c_str());
        reader->Update();
        vtkSmartPointer<vtkStructuredGrid> sgrid {reader->GetOutput()};
        return make_grid_info<T>( sgrid );
    }


    template <typename T, typename F>
    auto make_solver_case( vtkStructuredGrid* sgrid,  const std::vector< std::tuple< BladeInfo<T> , F > > &bld_info_lst )
    {
        SolverCase<T> solver_case{};
        solver_case.gi = make_grid_info<T>(sgrid);
        size_t iB{};
        auto &g = *(solver_case.gi->g);
        for(const auto &[ bld_info, f_k] : bld_info_lst)
        {
            solver_case.bld_info_lst.push_back(bld_info);

            for( auto i{bld_info.i1} ; i <= bld_info.i2; i++  )
            {
                auto m = ( i - bld_info.i1 ) / T( bld_info.i2 - bld_info.i1 );
                for( size_t j{}; j < solver_case.gi->nj; j++ )
                {
                    auto l = j / T(solver_case.gi->nj -1 );
                    g(i,j).k = f_k(m,l);
                    g(i,j).iB = iB;
                }
            }
            iB++;
        }
        return solver_case;
    }
} // namespace yams
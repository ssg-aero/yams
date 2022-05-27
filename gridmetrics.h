#pragma once
#include <datastorage.h>
#include <meridionalsolvercase.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridReader.h>
#include <gridreader.h>
#include <diffop.h>
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
                drqdm = D1_O2_i(g, i, j, fr, fm);
                dzqdm = D1_O2_i(g, i, j, fz, fm);
                g(i,j).phi = atan2(drqdm, dzqdm); // Stream line angle
                drqdl = D1_O2_j(g, i, j, fr, fl);
                dzqdl = D1_O2_j(g, i, j, fz, fl);
                g(i,j).gam = atan2(dzqdl, drqdl); // Span line angle
                // if(i==0)
                // {
                //     g(i,j).phi = atan2( (g(i+1,j).y-g(i,j).y) , (g(i+1,j).x-g(i,j).x) );
                // }
                // else if(i==ni-1)
                // {
                //     g(i,j).phi = atan2( (g(i,j).y-g(i-1,j).y) , (g(i,j).x-g(i-1,j).x) );
                // }
                // else
                // {
                //     g(i,j).phi = 0.5 * ( atan2( (g(i+1,j).y-g(i,j).y) , (g(i+1,j).x-g(i,j).x) ) +
                //                          atan2( (g(i,j).y-g(i-1,j).y) , (g(i,j).x-g(i-1,j).x) ));
                // }
                // if(j==0)
                // {
                //     g(i,j).gam = atan2( (g(i,j+1).x-g(i,j).x) , (g(i,j+1).y-g(i,j).y) );
                // }
                // else if(j==nj-1)
                // {
                //     g(i,j).gam = atan2( (g(i,j).x-g(i,j-1).x) , (g(i,j).y-g(i,j-1).y) );
                // }
                // else
                // {
                //     g(i,j).gam = 0.5 * ( atan2( (g(i,j+1).x-g(i,j).x) , (g(i,j+1).y-g(i,j).y) ) +
                //                          atan2( (g(i,j).x-g(i,j-1).x) , (g(i,j).y-g(i,j-1).y) ));
                // }

                g(i,j).cgp = std::cos( g(i,j).gam + g(i,j).phi);
                g(i,j).sgp = std::sin( g(i,j).gam + g(i,j).phi);
            }
        }
    }
    
    template <typename T>
    inline auto compute_curvature(MeridionalGrid<T> &g, bool interpolate = true)
    {
        size_t nim = g.nRows()-1;
        size_t nj = g.nCols();

        for (auto i = 1; i < nim; i++)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i,j).cur = D1_O2_i(g, i, j, fphi, fm);// TODO check why Aungier put -DphiDm
                // g(i,j).cur = 2 * ( atan2( (g(i+1,j).y-g(i,j).y) , (g(i+1,j).x-g(i,j).x) ) -
                //                    atan2( (g(i,j).y-g(i-1,j).y) , (g(i,j).x-g(i-1,j).x) ) ) /
                //                    (g(i+1,j).m - g(i-1,j).m);

            }
        }
        for (auto j = 0; j < nj; j++)
        {
            if( interpolate )
            {
                g(0, j).cur = 2. * g(1, j).cur - g(2, j).cur;
                g(nim, j).cur = 2. * g(nim - 1, j).cur - g(nim - 2, j).cur;
            }
            else
            {
                // g(0, j).cur  =  D1_O1_i_fw(g, 0 , j, fphi, fm);
                // g(nim, j).cur =  D1_O1_i_bw(g, nim, j, fphi, fm);
                g(0, j).cur  =  0.;
                g(nim, j).cur =  0.;
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
    
    template <typename T>
    inline auto compute_curvature(MeridionalGrid<T> &g,const Array2d<Grid2dMetricsPoint<T>> &g_metrics)
    {
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        T d_ksi = 1. / (ni - 1.);
        T d_eth = 1. / (nj - 1.);

        // for (auto i = 1; i < ni-1; i++)
        // {
        //     for (auto j = 0; j < nj; j++)
        //     {
        //         // g(i,j).cur = D1_O2_i(g, i, j, fphi, fm);// TODO check why Aungier put -DphiDm
        //         g(i,j).cur = D1_O2_dx1(g,g_metrics, i, j, d_ksi, d_eth, fphi);// TODO check why Aungier put -DphiDm
        //     }
        // }
        // for (auto j = 0; j < nj; j++) // extrapolate on bounds
        // {
        //     g(0, j).cur = 2. * g(1, j).cur - g(2, j).cur;
        //     g(ni-1, j).cur = 2. * g(ni - 2, j).cur - g(ni - 3, j).cur;
        // }
        for (auto i = 0; i < ni; i++)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i,j).cur = D1_O2_dx1(g,g_metrics, i, j, d_ksi, d_eth, fphi);// TODO check why Aungier put -DphiDm
            }
        }
    }

    template <typename T>
    inline auto compute_grid_metrics(MeridionalGrid<T> &g, Array2d<Grid2dMetricsPoint<T>> &g_metrics, const auto &f_m, const auto &f_l)
    {
        compute_abscissas(g);
        compute_metrics(g,f_m,f_l,g_metrics);
        compute_angles(g);
        // compute_curvature(g,g_metrics);
        compute_curvature(g);
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

    template <typename T>
    auto apply_blade_info(SolverCase<T> &solver_case)
    {
        auto nj = solver_case.gi->nj;
        auto &g = *(solver_case.gi->g);
        size_t iB{};
        for(const auto &bld_info : solver_case.bld_info_lst)
        {
            for( auto i{bld_info.i1} ; i <= bld_info.i2; i++  )
            {
                for( size_t j{}; j < nj; j++ )
                {
                    auto u = ( g(i,j).m - g(bld_info.i1,j).m ) / ( g(bld_info.i2,j).m - g(bld_info.i1,j).m );
                    auto v = ( g(i,j).l - g(i,0).l ) / ( g(i,nj-1).l - g(i,0).l );
                    if(bld_info.k)
                        g(i,j).k = bld_info.k(u,v);
                    if(bld_info.tb)
                    {
                        auto r  = g(i, j).y;
                        auto tb_= bld_info.z_ * bld_info.tb(u, v) / std::cos( g(i,j).k ); // projected total thickness
                        g(i, j).th_ = tb_ / r; // effective tangential span
                    }
                    if(bld_info.eps)
                    {
                        g(i, j).eps = bld_info.eps(u, v);
                    }
                    g(i,j).iB = iB;
                }
            }
            iB++;
        }
    }

    template <typename T>
    auto make_solver_case( vtkStructuredGrid* sgrid,  const std::vector< BladeInfo<T> > &bld_info_lst )
    {
        SolverCase<T> solver_case{};
        solver_case.gi = make_grid_info<T>(sgrid);
        size_t iB{};
        auto &g = *(solver_case.gi->g);
        for(const auto &bld_info : bld_info_lst)
        {
            solver_case.bld_info_lst.push_back(bld_info);
        }
        apply_blade_info(solver_case);
        return solver_case;
    }
} // namespace yams
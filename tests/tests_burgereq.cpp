#include <gtest/gtest.h>
// #include <tbb/parallel_for.h>
#include <gbs/maths.h>


#include <vtkChartXY.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
// #include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPen.h>
#include <vtkPlot.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>
template <typename T>
void plotXY(const std::vector<T> &X, const std::vector<std::vector<T>> &Y)
{
    // vtkNew<vtkNamedColors> colors;

    // Create a table with some points in it.
    vtkNew<vtkTable> table;

    vtkNew<vtkFloatArray> arrX;
    arrX->SetName("X Axis");
    table->AddColumn(arrX);
    // Set up the view
    vtkNew<vtkContextView> view;
    view->GetRenderWindow()->SetWindowName("LinePlot");
    // Add multiple line plots, setting the colors etc.
    vtkNew<vtkChartXY> chart;

    for (size_t j{}; j < Y.size(); j++)
    {
        vtkNew<vtkFloatArray> arrU;
        arrU->SetName("U");
        table->AddColumn(arrU);

        // Fill in the table with some example values.
        int numPoints = X.size();
        table->SetNumberOfRows(numPoints);
        for (int i = 0; i < numPoints; ++i)
        {
            table->SetValue(i, 0, X[i]);
            table->SetValue(i, j + 1, Y[j][i]);
        }

        view->GetScene()->AddItem(chart);
        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetInputData(table, 0, j + 1);
        line->SetColor(0, 255, 0, 255);
        line->SetWidth(1.0);
    }

    // Start interactor
    view->GetRenderWindow()->Render();
    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();
}

TEST(tests_burgereq, lax)
{
    auto f_diff_ct_O2 = [](size_t i, const auto &E, const auto &XM)
    {
        return (E[i + 1] - E[i - 1]) / XM[i] /2.;
    };

    auto f_diff_fw_O1 = [](size_t i, const auto &E, const auto &XM)
    {
        return (E[i + 1] - E[i]) / XM[i];
    };

    auto f_diff_bw_O1 = [](size_t i, const auto &E, const auto &XM)
    {
        return (E[i] - E[i - 1]) / XM[i];
    };

    auto f_E = [](const auto &Q)
    {
        return Q * Q / 2.;
    };

    // auto f_lax_scheme = [&f_diff_ct_O2](size_t i, const auto &E, const auto &U, const auto &XM, auto dt)
    // {
    //     return 0.5 * ( U[i+1] + U[i-1]) -  dt * f_diff_ct_O2(i,E,XM);
    // };


    auto f_lax_scheme = [&f_diff_ct_O2, &f_E]
    (
        const auto &U, 
        auto &E, 
        auto &U_, 
        const auto &X, 
        const auto  &Xm, 
        auto dt
    )
    {
        // tbb::parallel_for(tbb::blocked_range<size_t>(1, U.size() - 1),
        //                   [&](tbb::blocked_range<size_t> r)
        //                   {
        //                       for (int i = r.begin(); i < r.end(); ++i)
        //                       {
        //                           E[i] = f_E(U[i]);
        //                       }
        //                   });
        std::transform( // <- seems faster
            // std::execution::par, // even faster if not //
            U.begin(),U.end(),E.begin(),f_E
        );

        auto residual{0.};
        // tbb::parallel_for(tbb::blocked_range<size_t>(1, U.size() - 1),
        //                   [&](tbb::blocked_range<size_t> r)
        //                   {
        //                       for (int i = r.begin(); i < r.end(); ++i)
        //                       {
        //                           U_[i] = f_scheme(i, E, U, Xm, dt);
        //                           residual += fabs(U_[i] = U[i]);
        //                       }
        //                   });
        for (int i = 1; i < U.size() -1; ++i)
        {
            U_[i] =0.5 * ( U[i+1] + U[i-1]) -  dt * f_diff_ct_O2(i,E,Xm);;
            residual += fabs(U_[i] - U[i]);
        }
        return log10(residual / (U.size() - 2));
    };

        auto f_mac_cormack = [&f_diff_fw_O1,&f_diff_bw_O1, &f_E]
    (
        const auto &U, 
        auto &E, 
        auto &U_, 
        const auto &X, 
        const auto  &Xm, 
        auto dt
    )
    {
        std::transform( 
            // std::execution::par, // <- seems faster
            U.begin(),U.end(),E.begin(),f_E
        );

        auto residual{0.};

        for (int i = 1; i < U.size() -1; ++i)
        {
            U_[i] = U[i] - dt  * f_diff_fw_O1(i,E,Xm);
        }
        std::transform( 
            // std::execution::par, // <- seems faster
            U_.begin(),U_.end(),E.begin(),f_E
        );

        for (int i = 1; i < U.size() -1; ++i)
        {
            U_[i] = 0.5 *( U_[i] + U[i] - dt  * f_diff_bw_O1(i,E,Xm) );
            residual += fabs(U_[i] - U[i]);
        }
        return log10(residual / (U.size() - 2));
    };

    auto f_init = [](auto &U,const auto &X)
    {
        std::transform(
            std::execution::par,
            X.begin(),X.end(),U.begin(),[](auto x){return x < 2.0 ? 1. : 0.;}
        );
    };

    auto f_metrics = [](const auto &X,auto &Xm)
    {
        Xm.front() = X[1] - X [0];
        Xm.back() = X[X.size()-1] - X [X.size()-2];
        for( auto i {1} ; i < X.size() -1 ; i++)
        {
            Xm[i] = (X[i+1]-X[i-1])/2.;
        }
    };

    size_t n  = 41;
    size_t nit= 18;
    
    auto X = gbs::make_range(0.,4.,n);

    auto solve = [&f_metrics,&f_init](const auto &X,size_t n_it,const auto &f_scheme)
    {

        std::vector<double> U(X.size());
        std::vector<double> U_(X.size());
        std::vector<double> E(X.size());
        std::vector<double> Xm(X.size());
        std::vector<double> residual;
        std::vector<double> It;

        f_init(U,X);
        f_init(U_,X);
        f_metrics(X,Xm);
        size_t it =0;
        double dt = Xm[0] ;
        while (it <n_it)
        {
            residual.push_back(f_scheme(U, E, U_, X, Xm,  dt));
            It.push_back(it);
            std::swap(U_,U);
            // apply_bc
            it++;
        }
        return std::make_tuple(It,residual,U);
    };

    // auto [It,residual,U] = solve(X,nit,f_lax_scheme);
    auto [It,residual,U] = solve(X,nit,f_mac_cormack);

    plotXY(It,{residual});
    plotXY(X,{U});

}
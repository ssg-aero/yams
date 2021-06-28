#include <gtest/gtest.h>
#include <gbs/curves>


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
void plotXY(const std::vector<T> &X, const std::vector<T> &Y)
{
    // vtkNew<vtkNamedColors> colors;

    // Create a table with some points in it.
    vtkNew<vtkTable> table;

    vtkNew<vtkFloatArray> arrX;
    arrX->SetName("X Axis");
    table->AddColumn(arrX);

    vtkNew<vtkFloatArray> arrU;
    arrU->SetName("U");
    table->AddColumn(arrU);

    // Fill in the table with some example values.
    int numPoints = X.size();
    table->SetNumberOfRows(numPoints);
    for (int i = 0; i < numPoints; ++i)
    {
        table->SetValue(i, 0, X[i]);
        table->SetValue(i, 1, Y[i]);
    }

    // Set up the view
    vtkNew<vtkContextView> view;
    view->GetRenderWindow()->SetWindowName("LinePlot");
    // view->GetRenderer()->SetBackground(colors->GetColor3d("SlateGray").GetData());

    // Add multiple line plots, setting the colors etc.
    vtkNew<vtkChartXY> chart;
    view->GetScene()->AddItem(chart);
    vtkPlot *line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 1);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(1.0);

    // For dotted line, the line type can be from 2 to 5 for different dash/dot
    // patterns (see enum in vtkPen containing DASH_LINE, value 2):
    // #ifndef WIN32
    //   line->GetPen()->SetLineType(vtkPen::DASH_LINE);
    // #endif
    // (ifdef-ed out on Windows because DASH_LINE does not work on Windows
    //  machines with built-in Intel HD graphics card...)

    // view->GetRenderWindow()->SetMultiSamples(0);

    // Start interactor
    view->GetRenderWindow()->Render();
    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();
}

template <typename T>
struct fl_data_inv_1d
{
    T rho;
    T u;
    T et;
};

using gbs::operator*;
using gbs::operator/;
using gbs::operator-;
using gbs::operator+;

auto assert_nan_container(const auto array){std::for_each(array.begin(),array.end(),[](const auto &v_){assert(v_==v_);});}

TEST(tests_euler1D, explicit_rk4)
{
    using T  = double;
    T R = 1716.;
    T gamma = 1.4;
    size_t n = 50;
    T l = 10.;
    auto dx = l / (n-1.);
    T t1 = 520.;
    T p1 = 2000.;
    T rho1=p1/R/t1;
    T a1 = std::sqrt(gamma*R*t1);
    T M1 = 0.5;
    T u1 = M1 * a1;
    T dt = 0.00005; 
    size_t n_it = 2500000;
    auto cfl_ = (u1+a1) * dt / dx;
    

    std::cout << "Computation made at CFL: " << cfl_ << std::endl;

    auto f_S = [](auto x)
    { return 1.398 + 0.347 * std::tanh(0.8 * x - 4); };
    // {return 1.398 - x * 0.02;};
    // {return 1.398;};

    auto f_Sx = [](auto x)
    { return 0.347 * 0.8 * (1 - std::tanh(0.8 * x - 4) * std::tanh(0.8 * x - 4)); };
    // {return -0.02;};
    // {return 0.;};

    auto Q_to_p = [gamma](const std::array<T, 3> &Q)
    {
        return (gamma - 1) * (Q[2] - 0.5 * Q[1] * Q[1] / Q[0]);
    };

    auto Q_to_t = [gamma,R,&Q_to_p](const std::array<T, 3> &Q)
    {
        auto p = Q_to_p(Q);
        return p / R / Q[0];
    };

    auto Q_to_u = [](const std::array<T, 3> &Q)
    {
        return Q[1] / Q[0];
    };

    auto Q_to_M = [&](const std::array<T, 3> &Q)
    {
        auto u = Q_to_u(Q);
        auto a = std::sqrt( gamma * R * Q_to_t(Q) );
        return u / a;
    };

    auto ptu_to_Q = [gamma,R](const std::array<T, 3> &ptu)
    {
        auto rho   = ptu[0]/R/ptu[1];
        auto et = ptu[0] / rho / (gamma-1.) + 0.5 * ptu[2] * ptu[2];
        return std::array<T,3>{
            rho,
            rho*ptu[2],
            rho*et
        };
    };

    auto Q_to_E = [gamma,&Q_to_u,&Q_to_p](const std::array<T, 3> &Q, T S)
    {
        auto rho = Q[0];
        auto u = Q_to_u(Q);
        auto p = Q_to_p(Q);
        std::array<T, 3> E{
            S*Q[1],
            S*(rho * u*u +p),
            S*(Q[2]+p)*u
        };
        // std::array<T, 3> E{
        //     Q[1],
        //     Q[1] * Q[1] / Q[0] + (gamma - 1) * (Q[2] - 0.5 * Q[1] * Q[1] / Q[0]),
        //     gamma * Q[2] - (gamma - 1) * 0.5 * (Q[1] * Q[1] / Q[0]) * Q[1] / Q[0]};
        assert(E[0]==E[0]);
        assert(E[1]==E[1]);
        assert(E[2]==E[2]);
        return E;
        // return S*E;
    };

    auto Q_to_H = [&] (const std::array<T, 3> &Q,T Sx)
    {
        return std::array<T,3>{
            0.,
            Sx*Q_to_p(Q),
            0.
        };
    };


    auto inlet = [u1,t1,R,gamma,&Q_to_p](auto &Q){
        auto p   = 2 * Q_to_p(Q[1]) - Q_to_p(Q[2]);
        // auto p   = Q_to_p(Q[1]) ;
        auto rho = p / R / t1;
        auto et  = p / (gamma-1) / rho + 0.5 * u1*u1;
        return std::array<T,3>{
            rho,
            rho*u1,
            rho * et
        };
    };

    auto outlet = [p1,R,gamma,&Q_to_t,&Q_to_u](auto &Q){
        auto Q_m1 = *std::next(Q.end(),-2);
        auto Q_m2 = *std::next(Q.end(),-3);
        auto u      = 2 * Q_m1[1] / Q_m1[0] - Q_m2[1] / Q_m2[0];
        auto t      = 2* Q_to_t(Q_m1)-Q_to_t(Q_m2);
        // auto u = Q_to_u(*std::next(Q.end(),-2));
        // auto t = Q_to_t(*std::next(Q.end(),-2));
        auto rho    = p1 / R / t;
        auto et  = p1 / (gamma-1) / rho + 0.5 * u*u;
        return std::array<T,3>{
            rho,
            rho*u,
            rho * et
        }; 
    };

    auto RK4_mod = [dt,&inlet,&outlet](auto &Q, auto &Q_next, const auto &F_Q, const auto &S, const auto &Sx)
    {
        auto count = Q.size() - 1;
        std::array<T,3> residual {};
        // Q_next.front() = inlet(Q);
        for (auto i{1}; i < count; i++)
        {
            Q_next[i] = Q[i] - dt / 4. / S[i] * F_Q(i, Q,S, Sx[i]);
            assert_nan_container(Q_next);
        }
        // Q_next.back() = outlet(Q);
        // Q_next.front() = inlet(Q_next);
        for (auto i{1}; i < count; i++)
        {
            Q_next[i] = Q[i] - dt / 3. / S[i] * F_Q(i, Q_next,S, Sx[i]);
            assert_nan_container(Q_next);
        }
        // Q_next.back() = outlet(Q_next);
        // Q_next.front() = inlet(Q_next);
        for (auto i{1}; i < count; i++)
        {
            Q_next[i] = Q[i] - dt / 2. / S[i] * F_Q(i, Q_next,S, Sx[i]);
            assert_nan_container(Q_next);
        }
        // Q_next.back() = outlet(Q_next);
        // Q_next.front() = inlet(Q_next);
        for (auto i{1}; i < count; i++)
        {
            Q_next[i] = Q[i] - dt / S[i] * F_Q(i, Q_next,S, Sx[i]);
            assert_nan_container(Q_next);
            for(size_t j {}; j < 3 ; j++){ residual[j]+=fabs(Q[i][j]-Q_next[i][j])/dt;} 
        }
        return residual / (count -1.);
        // Q_next.back() = outlet(Q_next);
    };

    auto cfl = 0.1;
    auto McCormack = [cfl,gamma,R, &inlet, &outlet,&Q_to_E,&Q_to_H,dx,&Q_to_t](auto &Q, auto &Q_next, const auto &S, const auto &Sx)
    {
        auto count = Q.size() - 1;

        // auto ua_max = 0.;
        // for (auto i{1}; i < count; i++)
        // {
        //     auto t = Q_to_t(Q[i]);
        //     auto a = std::sqrt(gamma*R*t);
        //     ua_max = fmax(Q[i][1]/Q[i][0]+a,ua_max);
        // }
        // auto dt = cfl / ua_max * dx;

        std::array<T,3> residual {};
        for (auto i{1}; i < count; i++)
        {
            auto t = Q_to_t(Q[i]);
            auto a = std::sqrt(gamma*R*t);
            auto dt = cfl / (Q[i][1]/Q[i][0]+a) * dx;

            auto Ep = Q_to_E(Q[i + 1], S[i + 1]);
            auto Em = Q_to_E(Q[i], S[i]);
            Q_next[i] = Q[i] - dt / S[i] * ((Ep - Em) / dx - Q_to_H(Q[i], Sx[i]));
        }
        Q_next.front() = inlet(Q);
        Q_next.back() = outlet(Q);
        for (auto i{1}; i < count; i++)
        {
            auto t = Q_to_t(Q[i]);
            auto a = std::sqrt(gamma*R*t);
            auto dt = cfl / (Q[i][1]/Q[i][0]+a) * dx;
            
            auto Ep = Q_to_E(Q_next[i], S[i]);
            auto Em = Q_to_E(Q_next[i - 1], S[i - 1]);
            Q_next[i] = 0.5 * (Q[i] + Q_next[i] - dt / S[i] * ((Ep - Em) / dx + Q_to_H(Q_next[i], Sx[i])));
            for(size_t j {}; j < 3 ; j++){ residual[j]+=fabs(Q[i][j]-Q_next[i][j])/dt;} 
        }
        return residual / (count -1.);
    };

    auto F_Q = [&Q_to_E,&Q_to_H,dx](size_t i,const auto &Q,const auto &S,const T &Sx)
    {
        auto Ep = Q_to_E(Q[i+1],S[i+1]);
        auto Em = Q_to_E(Q[i-1],S[i-1]);
        auto res = (Ep-Em)/(2.*dx)-Q_to_H(Q[i],Sx);
        assert(res==res); 
        return res;
    };



    std::vector<std::array<T,3>> Q(n,ptu_to_Q({p1,t1,u1*0.95}));
    std::vector<std::array<T,3>> Q_next(n,ptu_to_Q({p1,t1,u1*0.95}));
    std::vector<T> S(n),Sx(n),X(n);
    X = gbs::make_range<T>(0.,l,n);
    std::transform(X.begin(),X.end(),S.begin(),f_S);
    std::transform(X.begin(),X.end(),Sx.begin(),f_Sx);
    auto count = n-1;
    T t = 0.;
    // while (t < 0.01)
    size_t it =0;
    std::vector<T> It,Res;
    while (it < n_it)
    {
        // auto residual = RK4_mod(Q,Q_next, F_Q, S, Sx);
        auto residual = McCormack(Q,Q_next, S, Sx);
        Q_next.front() = inlet(Q);
        Q_next.back() = outlet(Q);

        std::swap(Q, Q_next);
        It.push_back(it);
        Res.push_back(log10(residual[0]/Q[0][0]));
        t += dt;
        it++;
        // if(fmod(it,n_it/10)==0)
        // {
        //     std::vector<T> u(n);
        //     std::transform(Q.begin(), Q.end(), u.begin(), Q_to_M);
        //     plotXY(X, u);
        // }
    }

    std::vector<T> u(n);
    std::transform(Q.begin(),Q.end(),u.begin(),Q_to_M);
    plotXY(X,u);
    plotXY(It,Res);
    // plotXY(X,S);
    // plotXY(X,Sx);
}
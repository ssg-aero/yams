#include <gtest/gtest.h>
#include <chrono>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <execution>
// #include <xtensor/xarray.hpp>
TEST(vectorization, 1d)
{

    // using namespace std;
    // using T = double;
    // using TT = chrono::nanoseconds;
    // T factor_time = 1e6;
    // size_t n = 1e2;
    // xt::xarray<T>::shape_type shape = {n};
    // xt::xarray<T> X_xt(shape);
    // xt::xarray<T> Y_xt(shape);

    // vector<T> X(n);
    // vector<T> Y(n);

    // auto start = chrono::steady_clock::now();

    // auto f_sin = [](T x){return sin(x);};

    // // for(size_t i{}; i < n; i++)
    // // {
    // //     Y_xt[i] = sin(X_xt[i]);
    // // }
    // Y_xt = xt::sin(X_xt);

    // auto elapsed = chrono::duration_cast<TT>(chrono::steady_clock::now() - start).count() / factor_time;
    // cout << elapsed << endl;

    // start = chrono::steady_clock::now(); 

    // for(size_t i{}; i < n; i++)
    // {
    //     // Y[i] = f_sin(X[i]);
    //     Y[i] = std::sin(X[i]);
    // }

    // elapsed = chrono::duration_cast<TT>(chrono::steady_clock::now() - start).count() / factor_time;
    // cout << elapsed << endl;

    // for(size_t i{}; i < n; i++)
    // {
    //     ASSERT_NEAR(Y_xt[i], Y[i], 1e-6);
    // }

    // start = chrono::steady_clock::now();

    // // #pragma omp parallel for 
    // #pragma loop( no_vector )
    // // #pragma loop( hint_parallel( 0 ) )
    // for(int i{}; i < n; i++)
    // {
    //     Y[i] = f_sin(X[i]);
    // }

    // elapsed = chrono::duration_cast<TT>(chrono::steady_clock::now() - start).count() / factor_time;
    // cout << elapsed << endl;

    // start = chrono::steady_clock::now();

    // std::transform(
    //     X.begin(), X.end(), Y.begin(), f_sin
    // );

    // elapsed = chrono::duration_cast<TT>(chrono::steady_clock::now() - start).count() / factor_time;
    // cout << elapsed << endl;

    // start = chrono::steady_clock::now();

    // std::transform(
    //     std::execution::par,
    //     X.begin(), X.end(), Y.begin(), f_sin
    // );

    // elapsed = chrono::duration_cast<TT>(chrono::steady_clock::now() - start).count() / factor_time;
    // cout << elapsed << endl;
}

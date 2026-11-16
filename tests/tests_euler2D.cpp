#include <gtest/gtest.h>
#include <gbs/curves>
#include <gbs-mesh/tfi.h>
#include <gbs-mesh/smoothing.h>
#include <gbs-render/vtkgridrender.h>

#include <gas/gasModel.h>
#include <euler2d/euler2dSolver.h>

#include <chrono>
#include <limits>

using namespace gbs;
using namespace yams;

#include <euler2d/solverData.h>

TEST(euler2D, SolverData4Init)
{
    using T = double;
    size_t ni{11};
    size_t nj{5};
    auto n = ni * nj;
    SolverData4<T> U{ni, nj};
    ASSERT_EQ(U.ni,ni);
    ASSERT_EQ(U.nj,nj);
    T u1{1.}, u2{2.}, u3{3.}, u4{4.};
    U.init(u1, u2, u3, u4);
    auto & [ U1, U2, U3, U4] = U.d;
    ASSERT_EQ(U1.size(),n);
    ASSERT_EQ(U2.size(),n);
    ASSERT_EQ(U3.size(),n);
    ASSERT_EQ(U4.size(),n);
    std::for_each(U1.begin(), U1.end(),
        [u1](auto u){
            ASSERT_LT( std::abs(u-u1), std::numeric_limits<T>::epsilon());
        }
    );
    std::for_each(U2.begin(), U2.end(),
        [u2](auto u){
            ASSERT_LT( std::abs(u-u2), std::numeric_limits<T>::epsilon());
        }
    );
    std::for_each(U3.begin(), U3.end(),
        [u3](auto u){
            ASSERT_LT( std::abs(u-u3), std::numeric_limits<T>::epsilon());
        }
    );
    std::for_each(U4.begin(), U4.end(),
        [u4](auto u){
            ASSERT_LT( std::abs(u-u4), std::numeric_limits<T>::epsilon());
        }
    );
}

TEST(euler2D, SolverData4Perf)
{
    using T = double;
    const T half = 0.5;
    const size_t ni{100}; // n cells
    const size_t nj{50};  // 
    // const size_t ni{10}; // n cells
    // const size_t nj{5};  // 
    const auto nig = ni + 4; // n cells with ghosts
    const auto njg = nj + 4; //
    const auto nip = ni + 1; // n grid pts
    const auto njp = nj + 1; // n grid pts
    SolverData4<T> U{nig, njg,yams::make_tuple(1.,10.,20.,1e5/0.4+0.5*(10*10+20+20))};
    SolverData4<T> E_EW{nip, nj,yams::make_tuple(0.,0.,0.,0.)};
    SolverData4<T> F_EW{nip, nj,yams::make_tuple(0.,0.,0.,0.)};
    vector<T> Sx_EW(nip*nj, 0.02), Sy_EW(nip*nj, 0.0 );
    vector<T> Sx_SN(ni*njp, 0.0 ), Sy_SN(ni*njp, 0.02);
    SolverData4<T> E_SN{ni, njp,yams::make_tuple(0.,0.,0.,0.)};
    SolverData4<T> F_SN{ni, njp,yams::make_tuple(0.,0.,0.,0.)};
    std::vector<T> gamma(nig* njg, 1.4);
    {
        size_t n =101;
        std::vector<double> A(n),B(n+1);
        auto pA  = A.data();
        auto p1B = B.data();
        auto p2B = B.data()+1;
        for(size_t i{}; i < n; i++)
        {
            pA[i] = ( p1B[i] + p2B[i] ) * half;
        }

        for(size_t i{}; i < n; i++)
        {
            A[i] = ( B[i] + B[i+1] ) * half;
        }
    }
    {
        vector<T> g(nip), u1(nip), u2(nip), u3(nip), u4(nip), p(nip);
        auto start = std::chrono::steady_clock::now();
        auto &[U1, U2, U3, U4] = U.d;
        auto &[E_EW1, E_EW2, E_EW3, E_EW4] = E_EW.d;
        auto &[F_EW1, F_EW2, F_EW3, F_EW4] = F_EW.d;

        auto p1gamma =gamma.data()+1;
        auto p2gamma =gamma.data()+2;
        auto p1U1 =U1.data()+1;
        auto p2U1 =U1.data()+2;
        auto p1U2 =U2.data()+1;
        auto p2U2 =U2.data()+2;
        auto p1U3 =U3.data()+1;
        auto p2U3 =U3.data()+2;
        auto p1U4 =U4.data()+1;
        auto p2U4 =U4.data()+2;
        auto pE_EW1 =E_EW1.data();
        auto pE_EW2 =E_EW2.data();
        auto pE_EW3 =E_EW3.data();
        auto pE_EW4 =E_EW4.data();
        auto pF_EW1 =F_EW1.data();
        auto pF_EW2 =F_EW2.data();
        auto pF_EW3 =F_EW3.data();
        auto pF_EW4 =F_EW4.data();

        // #pragma omp for
        // for (long j{}; j < nj; j++)
        for (size_t j{}; j < nj; j++)
        {

            // #pragma loop( no_vector )
            for (size_t i{}; i < nip; i++)
            {
                g[i]  = ( p1gamma[i] + p2gamma[i] ) * half;
                u1[i] = ( p1U1[i] + p2U1[i] ) * half;
                u2[i] = ( p1U2[i] + p2U2[i] ) * half;
                u3[i] = ( p1U3[i] + p2U3[i] ) * half;
                u4[i] = ( p1U4[i] + p2U4[i] ) * half;
                p[i]  = (g[i]-1)*(u4[i]-(u2[i]*u2[i]/u1[i]+u3[i]*u3[i]/u1[i]) * half);
            }
            for (size_t i{}; i < nip; i++)
            {
                pE_EW1[i] = u2[i];
                pE_EW2[i] = u2[i]*u2[i]/u1[i] + p[i];
                pE_EW3[i] = u2[i]*u3[i]/u1[i];
                pE_EW4[i] = (u4[i]+p[i])*u2[i]/u1[i];
            // }
            // for (size_t i{}; i < nip; i++)
            // {
                pF_EW1[i] = u3[i];
                pF_EW2[i] = u2[i]*u3[i]/u1[i];
                pF_EW3[i] = u3[i]*u3[i]/u1[i] + p[i];
                pF_EW4[i] = (u4[i]+p[i])*u3[i]/u1[i];

            }
            p1gamma+=nig;
            p2gamma+=nig;
            p1U1+=nig;
            p2U1+=nig;
            p1U2+=nig;
            p2U2+=nig;
            p1U3+=nig;
            p2U3+=nig;
            p1U4+=nig;
            p2U4+=nig;
            pE_EW1+=nip;
            pE_EW2+=nip;
            pE_EW3+=nip;
            pE_EW4+=nip;
            pF_EW1+=nip;
            pF_EW2+=nip;
            pF_EW3+=nip;
            pF_EW4+=nip;
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Cells volume and center single loop:" << elapsed << "mcs" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        auto &[U1, U2, U3, U4] = U.d;
        auto &[E_EW1, E_EW2, E_EW3, E_EW4] = E_EW.d;
        auto &[F_EW1, F_EW2, F_EW3, F_EW4] = F_EW.d;

        auto offset1 = 1;
        auto offset2 = 2;

        auto p1gamma =gamma.data()+offset1;
        auto p2gamma =gamma.data()+2;
        auto p1U1 =U1.data()+offset1;
        auto p2U1 =U1.data()+offset2;
        auto p1U2 =U2.data()+offset1;
        auto p2U2 =U2.data()+offset2;
        auto p1U3 =U3.data()+offset1;
        auto p2U3 =U3.data()+offset2;
        auto p1U4 =U4.data()+offset1;
        auto p2U4 =U4.data()+offset2;
        auto pE_EW1 =E_EW1.data();
        auto pE_EW2 =E_EW2.data();
        auto pE_EW3 =E_EW3.data();
        auto pE_EW4 =E_EW4.data();
        auto pF_EW1 =F_EW1.data();
        auto pF_EW2 =F_EW2.data();
        auto pF_EW3 =F_EW3.data();
        auto pF_EW4 =F_EW4.data();
        auto pSx_EW =Sx_EW.data();
        auto pSy_EW =Sy_EW.data();


        
        // #pragma omp for
        // for (long j{}; j < nj; j++)
        for (size_t j{}; j < nj; j++)
        {

            // #pragma loop( no_vector )
            for (size_t i{}; i < nip; i++)
            {
                auto g  = ( p1gamma[i] + p2gamma[i] ) * half;
                auto u1 = ( p1U1[i] + p2U1[i] ) * half;
                auto u2 = ( p1U2[i] + p2U2[i] ) * half;
                auto u3 = ( p1U3[i] + p2U3[i] ) * half;
                auto u4 = ( p1U4[i] + p2U4[i] ) * half;
                auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                auto S = pSx_EW[i];
                pE_EW1[i] = u2 * S;
                pE_EW2[i] = (u2*u2/u1 + p)*S;
                pE_EW3[i] = u2*u3/u1*S;
                pE_EW4[i] = (u4+p)*u2/u1*S;
            }
            for (size_t i{}; i < nip; i++)
            {
                auto g  = ( p1gamma[i] + p2gamma[i] ) * half ;
                auto u1 = ( p1U1[i] + p2U1[i] ) * half;
                auto u2 = ( p1U2[i] + p2U2[i] ) * half;
                auto u3 = ( p1U3[i] + p2U3[i] ) * half;
                auto u4 = ( p1U4[i] + p2U4[i] ) * half;
                auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                auto S = pSy_EW[i];
                pF_EW1[i] = u3*S;
                pF_EW2[i] = u2*u3/u1*S;
                pF_EW3[i] = (u3*u3/u1 + p)*S;
                pF_EW4[i] = (u4+p)*u3/u1*S;

            }
            p1gamma+=nig;
            p2gamma+=nig;
            p1U1+=nig;
            p2U1+=nig;
            p1U2+=nig;
            p2U2+=nig;
            p1U3+=nig;
            p2U3+=nig;
            p1U4+=nig;
            p2U4+=nig;
            pE_EW1+=nip;
            pE_EW2+=nip;
            pE_EW3+=nip;
            pE_EW4+=nip;
            pF_EW1+=nip;
            pF_EW2+=nip;
            pF_EW3+=nip;
            pF_EW4+=nip;
            pSx_EW+=nip;
            pSy_EW+=nip;
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Faces Flux EW vectorized:" << elapsed << "mcs" << std::endl;
    }

    {
        auto start = std::chrono::steady_clock::now();
        auto &[U1, U2, U3, U4] = U.d;
        auto &[E_EW1, E_EW2, E_EW3, E_EW4] = E_EW.d;
        auto &[F_EW1, F_EW2, F_EW3, F_EW4] = F_EW.d;

        auto offset1 = 1; // ngvtx
        auto offset2 = 2; // nghost
       
        // #pragma omp for
        // for (long j{}; j < nj; j++)
        for (size_t j{}; j < nj; j++)
        // std::vector<size_t> idx_j(nj);
        // std::generate(idx_j.begin(), idx_j.end(), [n=0]() mutable{return n++;});
        // std::for_each(
        //     std::execution::par,
        //     idx_j.begin(),idx_j.end(),[&](auto j)
        {
            auto p1gamma =gamma.data()+offset1;
            auto p2gamma =gamma.data()+offset2;
            auto p1U1 =U1.data()+offset1 + nig*j;
            auto p2U1 =U1.data()+offset2 + nig*j;
            auto p1U2 =U2.data()+offset1 + nig*j;
            auto p2U2 =U2.data()+offset2 + nig*j;
            auto p1U3 =U3.data()+offset1 + nig*j;
            auto p2U3 =U3.data()+offset2 + nig*j;
            auto p1U4 =U4.data()+offset1 + nig*j;
            auto p2U4 =U4.data()+offset2 + nig*j;
            auto pE_EW1 =E_EW1.data() + nip*j;
            auto pE_EW2 =E_EW2.data() + nip*j;
            auto pE_EW3 =E_EW3.data() + nip*j;
            auto pE_EW4 =E_EW4.data() + nip*j;
            auto pF_EW1 =F_EW1.data() + nip*j;
            auto pF_EW2 =F_EW2.data() + nip*j;
            auto pF_EW3 =F_EW3.data() + nip*j;
            auto pF_EW4 =F_EW4.data() + nip*j;
            auto pSx_EW =Sx_EW.data() + nip*j;
            auto pSy_EW =Sy_EW.data() + nip*j;
            // #pragma loop( no_vector )
            for (size_t i{}; i < nip; i++)
            {
                auto g  = ( p1gamma[i] + p2gamma[i] ) * half;
                auto u1 = ( p1U1[i] + p2U1[i] ) * half;
                auto u2 = ( p1U2[i] + p2U2[i] ) * half;
                auto u3 = ( p1U3[i] + p2U3[i] ) * half;
                auto u4 = ( p1U4[i] + p2U4[i] ) * half;
                auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                auto S = pSx_EW[i];
                pE_EW1[i] = u2 * S;
                pE_EW2[i] = (u2*u2/u1 + p)*S;
                pE_EW3[i] = u2*u3/u1*S;
                pE_EW4[i] = (u4+p)*u2/u1*S;
            }
            for (size_t i{}; i < nip; i++)
            {
                auto g  = ( p1gamma[i] + p2gamma[i] ) * half;
                auto u1 = ( p1U1[i] + p2U1[i] ) * half;
                auto u2 = ( p1U2[i] + p2U2[i] ) * half;
                auto u3 = ( p1U3[i] + p2U3[i] ) * half;
                auto u4 = ( p1U4[i] + p2U4[i] ) * half;
                auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                auto S = pSy_EW[i];
                pF_EW1[i] = u3*S;
                pF_EW2[i] = u2*u3/u1*S;
                pF_EW3[i] = (u3*u3/u1 + p)*S;
                pF_EW4[i] = (u4+p)*u3/u1*S;

            }

        }
        // );
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Faces Flux EW vectorized, pralelalizable:" << elapsed << "mcs" << std::endl;
    }

    {
        std::vector<size_t> idx_j(nj);
        std::generate(idx_j.begin(), idx_j.end(), [n=0]() mutable{return n++;});
        auto start = std::chrono::steady_clock::now();
        auto &[U1, U2, U3, U4] = U.d;
        auto &[E_EW1, E_EW2, E_EW3, E_EW4] = E_EW.d;
        auto &[F_EW1, F_EW2, F_EW3, F_EW4] = F_EW.d;

        auto p1gamma =gamma.data()+1;
        auto p2gamma =gamma.data()+2;
        auto p1U1 =U1.data()+1;
        auto p2U1 =U1.data()+2;
        auto p1U2 =U2.data()+1;
        auto p2U2 =U2.data()+2;
        auto p1U3 =U3.data()+1;
        auto p2U3 =U3.data()+2;
        auto p1U4 =U4.data()+1;
        auto p2U4 =U4.data()+2;
        auto pE_EW1 =E_EW1.data();
        auto pE_EW2 =E_EW2.data();
        auto pE_EW3 =E_EW3.data();
        auto pE_EW4 =E_EW4.data();
        auto pF_EW1 =F_EW1.data();
        auto pF_EW2 =F_EW2.data();
        auto pF_EW3 =F_EW3.data();
        auto pF_EW4 =F_EW4.data();

        vector<T> g(nip), u1(nip), u2(nip), u3(nip), u4(nip), p(nip);
        std::for_each(
            // std::execution::par,
            idx_j.begin(),idx_j.end(),[&](auto j)
        {

            // #pragma loop( no_vector )
            for (size_t i{}; i < nip; i++)
            {
                g[i]  = ( p1gamma[i] + p2gamma[i] ) * half;
                u1[i] = ( p1U1[i] + p2U1[i] ) * half;
                u2[i] = ( p1U2[i] + p2U2[i] ) * half;
                u3[i] = ( p1U3[i] + p2U3[i] ) * half;
                u4[i] = ( p1U4[i] + p2U4[i] ) * half;
                p[i]  = (g[i]-1)*(u4[i]-(u2[i]*u2[i]/u1[i]+u3[i]*u3[i]/u1[i]) * half);
            }
            for (size_t i{}; i < nip; i++)
            {
                pE_EW1[i] = u2[i];
                pE_EW2[i] = u2[i]*u2[i]/u1[i] + p[i];
                pE_EW3[i] = u2[i]*u3[i]/u1[i];
                pE_EW4[i] = (u4[i]+p[i])*u2[i]/u1[i];
            }
            for (size_t i{}; i < nip; i++)
            {
                pF_EW1[i] = u3[i];
                pF_EW2[i] = u2[i]*u3[i]/u1[i];
                pF_EW3[i] = u3[i]*u3[i]/u1[i] + p[i];
                pF_EW4[i] = (u4[i]+p[i])*u3[i]/u1[i];

            }
            p1gamma+=nig;
            p2gamma+=nig;
            p1U1+=nig;
            p2U1+=nig;
            p1U2+=nig;
            p2U2+=nig;
            p1U3+=nig;
            p2U3+=nig;
            p1U4+=nig;
            p2U4+=nig;
            pE_EW1+=nip;
            pE_EW2+=nip;
            pE_EW3+=nip;
            pE_EW4+=nip;
            pF_EW1+=nip;
            pF_EW2+=nip;
            pF_EW3+=nip;
            pF_EW4+=nip;
        });
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Cells volume and center single loop:" << elapsed << "mcs" << std::endl;
    }
    {
        std::vector<size_t> idx_j(nj);
        std::generate(idx_j.begin(), idx_j.end(), [n=0]() mutable{return n++;});
        auto start = std::chrono::steady_clock::now();
        auto &[U1, U2, U3, U4] = U.d;
        auto &[E_EW1, E_EW2, E_EW3, E_EW4] = E_EW.d;
        auto &[F_EW1, F_EW2, F_EW3, F_EW4] = F_EW.d;

        auto p1gamma =gamma.data()+1;
        auto p2gamma =gamma.data()+2;
        auto p1U1 =U1.data()+1;
        auto p2U1 =U1.data()+2;
        auto p1U2 =U2.data()+1;
        auto p2U2 =U2.data()+2;
        auto p1U3 =U3.data()+1;
        auto p2U3 =U3.data()+2;
        auto p1U4 =U4.data()+1;
        auto p2U4 =U4.data()+2;
        auto pE_EW1 =E_EW1.data();
        auto pE_EW2 =E_EW2.data();
        auto pE_EW3 =E_EW3.data();
        auto pE_EW4 =E_EW4.data();
        auto pF_EW1 =F_EW1.data();
        auto pF_EW2 =F_EW2.data();
        auto pF_EW3 =F_EW3.data();
        auto pF_EW4 =F_EW4.data();

        
        std::for_each(
            // std::execution::par,
            idx_j.begin(),idx_j.end(),[&](auto j)
        {

            // #pragma loop( no_vector )
            for (size_t i{}; i < nip; i++)
            {
                auto g  = ( p1gamma[i] + p2gamma[i] ) * half;
                auto u1 = ( p1U1[i] + p2U1[i] ) * half;
                auto u2 = ( p1U2[i] + p2U2[i] ) * half;
                auto u3 = ( p1U3[i] + p2U3[i] ) * half;
                auto u4 = ( p1U4[i] + p2U4[i] ) * half;
                auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                pE_EW1[i] = u2;
                pE_EW2[i] = u2*u2/u1 + p;
                pE_EW3[i] = u2*u3/u1;
                pE_EW4[i] = (u4+p)*u2/u1;
            }
            for (size_t i{}; i < nip; i++)
            {
                auto g  = ( p1gamma[i] + p2gamma[i] ) * half;
                auto u1 = ( p1U1[i] + p2U1[i] ) * half;
                auto u2 = ( p1U2[i] + p2U2[i] ) * half;
                auto u3 = ( p1U3[i] + p2U3[i] ) * half;
                auto u4 = ( p1U4[i] + p2U4[i] ) * half;
                auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                pF_EW1[i] = u3;
                pF_EW2[i] = u2*u3/u1;
                pF_EW3[i] = u3*u3/u1 + p;
                pF_EW4[i] = (u4+p)*u3/u1;

            }
            p1gamma+=nig;
            p2gamma+=nig;
            p1U1+=nig;
            p2U1+=nig;
            p1U2+=nig;
            p2U2+=nig;
            p1U3+=nig;
            p2U3+=nig;
            p1U4+=nig;
            p2U4+=nig;
            pE_EW1+=nip;
            pE_EW2+=nip;
            pE_EW3+=nip;
            pE_EW4+=nip;
            pF_EW1+=nip;
            pF_EW2+=nip;
            pF_EW3+=nip;
            pF_EW4+=nip;
        });
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Cells volume and center single loop:" << elapsed << "mcs" << std::endl;
    }
    {
        std::vector<std::array<T,4>> E_EW(nip*nj),F_EW(nip*nj), U(nig*njg);
        auto start = std::chrono::steady_clock::now();
        auto E_EW_start = E_EW.begin();
        auto E_EW_end   = std::next(E_EW.begin(),nip);
        auto F_EW_start = F_EW.begin();
        auto F_EW_end   = std::next(F_EW.begin(),nip);
        auto Ul = std::next(U.begin(),1);
        auto Ur = std::next(U.begin(),2);
        for (size_t j{}; j < nj; j++)
        {
            std::transform(
                E_EW_start, E_EW_end,
                Ul,Ur,
                [g=1.4, half](const auto &ul, const auto &ur )
                {
                    auto u1 = 0.5 * (ul[0]+ur[0] );
                    auto u2 = 0.5 * (ul[1]+ur[1] );
                    auto u3 = 0.5 * (ul[2]+ur[2] );
                    auto u4 = 0.5 * (ul[3]+ur[3] );
                    auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                    return std::array<T,4>{
                        u2, u2*u2/u1 + p, u2*u3/u1, (u4+p)*u2/u1
                    };
                }
            );
            if(j<nj-1)
            {
                E_EW_end = std::next(E_EW_end, nip);
                E_EW_start = std::next(E_EW_start, nip);
            }
            std::transform(
                F_EW_start, F_EW_end,
                Ul,Ur,
                [g=1.4, half](const auto &ul, const auto &ur )
                {
                    auto u1 = 0.5 * (ul[0]+ur[0] );
                    auto u2 = 0.5 * (ul[1]+ur[1] );
                    auto u3 = 0.5 * (ul[2]+ur[2] );
                    auto u4 = 0.5 * (ul[3]+ur[3] );
                    auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                    return std::array<T,4>{
                        u3, u2*u3/u1, u3*u3/u1 + p, (u4+p)*u3/u1
                    };
                }
            );
            if(j<nj-1)
            {
                F_EW_end   = std::next(F_EW_end,nip);
                F_EW_start = std::next(F_EW_start,nip);
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Cells volume and center single loop:" << elapsed << "mcs" << std::endl;
    }
    
    {
        auto start = std::chrono::steady_clock::now();
        auto &[U1, U2, U3, U4] = U.d;
        auto &[E_SN1, E_SN2, E_SN3, E_SN4] = E_SN.d;
        auto &[F_SN1, F_SN2, F_SN3, F_SN4] = F_SN.d;

        auto offset1 = 1*nig+2;
        auto offset2 = 2*nig+2;

        auto p1gamma =gamma.data()+offset1;
        auto p2gamma =gamma.data()+offset1;
        auto p1U1 =U1.data()+offset1;
        auto p2U1 =U1.data()+offset2;
        auto p1U2 =U2.data()+offset1;
        auto p2U2 =U2.data()+offset2;
        auto p1U3 =U3.data()+offset1;
        auto p2U3 =U3.data()+offset2;
        auto p1U4 =U4.data()+offset1;
        auto p2U4 =U4.data()+offset2;
        auto pE_SN1 =E_SN1.data();
        auto pE_SN2 =E_SN2.data();
        auto pE_SN3 =E_SN3.data();
        auto pE_SN4 =E_SN4.data();
        auto pF_SN1 =F_SN1.data();
        auto pF_SN2 =F_SN2.data();
        auto pF_SN3 =F_SN3.data();
        auto pF_SN4 =F_SN4.data();
        auto pSx_SN =Sx_SN.data();
        auto pSy_SN =Sy_SN.data();

        
        // #pragma omp for
        // for (long j{}; j < njp; j++)
        for (size_t j{}; j < njp; j++)
        {

            // #pragma loop( no_vector )
            for (size_t i{}; i < ni; i++)
            {
                auto g  = ( p1gamma[i] + p2gamma[i] ) * half;
                auto u1 = ( p1U1[i] + p2U1[i] ) * half;
                auto u2 = ( p1U2[i] + p2U2[i] ) * half;
                auto u3 = ( p1U3[i] + p2U3[i] ) * half;
                auto u4 = ( p1U4[i] + p2U4[i] ) * half;
                auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                auto S = pSx_SN[i];
                pE_SN1[i] = u2*S;
                pE_SN2[i] = (u2*u2/u1 + p)*S;
                pE_SN3[i] = u2*u3/u1*S;
                pE_SN4[i] = (u4+p)*u2/u1*S;
            }
            for (size_t i{}; i < ni; i++)
            {
                auto g  = ( p1gamma[i] + p2gamma[i] ) * half;
                auto u1 = ( p1U1[i] + p2U1[i] ) * half;
                auto u2 = ( p1U2[i] + p2U2[i] ) * half;
                auto u3 = ( p1U3[i] + p2U3[i] ) * half;
                auto u4 = ( p1U4[i] + p2U4[i] ) * half;
                auto p  = (g-1)*(u4-(u2*u2/u1+u3*u3/u1) * half);
                auto S = pSy_SN[i];
                pF_SN1[i] = u3*S;
                pF_SN2[i] = u2*u3/u1*S;
                pF_SN3[i] = (u3*u3/u1 + p)*S;
                pF_SN4[i] = (u4+p)*u3/u1*S;

            }
            p1gamma+=nig;
            p2gamma+=nig;
            p1U1+=nig;
            p2U1+=nig;
            p1U2+=nig;
            p2U2+=nig;
            p1U3+=nig;
            p2U3+=nig;
            p1U4+=nig;
            p2U4+=nig;
            pE_SN1+=ni;
            pE_SN2+=ni;
            pE_SN3+=ni;
            pE_SN4+=ni;
            pF_SN1+=ni;
            pF_SN2+=ni;
            pF_SN3+=ni;
            pF_SN4+=ni;
            pSx_SN+=ni;
            pSy_SN+=ni;
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Faces Flux SN vectorized: " << elapsed << "mcs" << std::endl;
    }
}

// using T = double;
// using BSCurve2d = BSCurve<T,2>;
// template <typename T>
// using vector = std::vector<T>;
// template <typename... _Types>
// using tuple = std::tuple<_Types...>;



// namespace yams{
//     template <typename... _Types>
//     inline auto make_tuple(_Types&&... _Args) -> tuple<_Types...>
//     {
//         return std::make_tuple(_Args...);
//     }
//     template <size_t i>
//     auto get(auto &tuple)
//     {
//         return std::get<i>(tuple);
//     }
//     template <typename Tuple, typename Functor, size_t Index = 0>
//     auto tuple_for_each(const Tuple &tpl, const Functor &f) -> void
//     {
//         constexpr auto tuple_size = std::tuple_size_v<Tuple>;
//         if constexpr (Index < tuple_size)
//         {
//             f(std::get<Index>(tpl));
//             tuple_for_each<Tuple, Functor, Index + 1>(tpl, f);
//         }
//     }
// }

template<typename T>
inline auto build_channel_mesh(size_t ni_ = 30,size_t nj_ = 30, T bump = 0.1, T l = 2., T h = 1.)
{
    using BSCurve2d = BSCurve<T,2>;
    auto south_crv = std::make_shared<BSCurve2d>(BSCurve2d{
        {
            {0.0,0.0},
            {0.2 * l,0.0},
            {0.4 * l,bump},
            {0.6 * l,bump},
            {0.8 * l,0.0},
            {1.0 * l,0.0}
        },
        {0.,0.25, 0.75, 1.},
        {4,1,1, 4},
        3
    });
    auto north_crv = std::make_shared<BSCurve2d>(BSCurve2d{
        {
            {0.0, h},
            {l, h}
        },
        {0.,1.},
        {2, 2},
        1
    });
    auto west_crv = std::make_shared<BSCurve2d>(BSCurve2d{
        {
            {0.0,0.0},
            {0.0,h}
        },
        {0.,1.},
        {2, 2},
        1
    });
    auto east_crv = std::make_shared<BSCurve2d>(BSCurve2d{
        {
            {l,0.0},
            {l,h}
        },
        {0.,1.},
        {2, 2},
        1
    });

    return tfi_mesh_2d<T,2,1,1>({south_crv, north_crv},{west_crv, east_crv},{0., 1.},{0., 1.}, ni_, nj_);
}


// struct FlowData
// {
//     tuple<vector<T>, vector<T>, vector<T>, vector<T>> d; // tuple can be used thanks to thrust in cuda
//     FlowData(size_t n) {
//         d = yams::make_tuple(
//             vector<T>(n),
//             vector<T>(n),
//             vector<T>(n),
//             vector<T>(n)
//         );
//     }
//     void init(T v1, T v2, T v3, T v4)
//     {
//         std::fill(std::get<0>(d).begin(), std::get<0>(d).end(), v1);
//         std::fill(std::get<1>(d).begin(), std::get<1>(d).end(), v2);
//         std::fill(std::get<2>(d).begin(), std::get<2>(d).end(), v3);
//         std::fill(std::get<3>(d).begin(), std::get<3>(d).end(), v4);
//     }
// };

// template <typename T>
// class rk_integrator
// {
// private:
//     T
//     const std::array<T,5> rk_coefs;
// public:
//     rk_integrator(const std::array<T,5> &rk_coefs = {0.0533,0.1263,0.2375, 0.4414,1.0}) rk_coefs{rk_coefs} {}
//     template <typename U, typename _Func>
//     U integrate(U W, T dt, _Func R)
//     {
//         for(auto alf : rk_coefs)
//         {

//         }
//     }
// };

// TEST(tests_euler2D, rk_integrator)
// {

// }



TEST(tests_euler2D, gris_access_speeds)
{
    using T = double;
    const T half = 0.5;
    size_t ni_{51}, nj_{25};
    // size_t ni_{5}, nj_{3};
    T bump = 0.0;
    auto [ pts, ni, nj, n_iso_ksi, n_iso_eth] = build_channel_mesh(ni_, nj_, bump);
    elliptic_structured_smoothing(pts, ni, 0, nj-1, 0, ni-1, 300);
    auto njm=nj-1;
    auto nim=ni-1;
    auto njp=nj+1;
    auto nip=ni+1;
    auto nm = nim * njm; // N cells without ghosts
    auto np = nip * njp; // N cells with ghosts
    auto n = ni * nj;    // N vertices

    // Define flow description
    SolverData4<T> U{nip, njp}; // U = (rho, rho*u, rho*v, rho*et)
    SolverData4<T> R{nip, njp};
    T Pt_inlet = 1.05e5;
    T Tt_inlet = 288;
    T Ps_outlet = 1e5;
    // Init flow
    yams::GasPerfectGammaConstant<T> air{287.04, 1.4};
    T g = air.g(Tt_inlet, Pt_inlet);
    T M_outlet = sqrt( 2/(g-1)*(pow(Pt_inlet/Ps_outlet,(g-1)/g)-1));
    T Ts_outlet= Tt_inlet * (1+(g-1) * half*M_outlet*M_outlet);
    T u_outlet = M_outlet * air.a(Ts_outlet, Ps_outlet);
    T rho_outlet = air.rho(Ts_outlet, Ps_outlet);
    T et_outlet  = air.et(Ts_outlet, Ps_outlet, u_outlet);
    U.init(
        rho_outlet,
        rho_outlet*u_outlet,
        0.,
        et_outlet
    );
    // compute gog and volume

    vector<T> V(nm),X(n),Y(n),Xc(nm),Yc(nm);
    std::transform(std::execution::par, pts.begin(), pts.end(),X.begin(),[](const auto &pt){return pt[0];});
    std::transform(std::execution::par, pts.begin(), pts.end(),Y.begin(),[](const auto &pt){return pt[1];});
    
    {    
        auto start = std::chrono::steady_clock::now();
        for (long j{}; j < njm; j++)
        {
            for (long i{}; i < nim; i++)
            {
                auto x1 = X[i + ni * j];
                auto y1 = Y[i + ni * j];
                auto x2 = X[i + 1 + ni * j];
                auto y2 = Y[i + 1 + ni * j];
                auto x3 = X[i + 1 + ni * (j + 1)];
                auto y3 = Y[i + 1 + ni * (j + 1)];
                auto x4 = X[i + ni * (j + 1)];
                auto y4 = Y[i + ni * (j + 1)];
                V[i + nim * j]  = T(0.5) * ((x1 - x3) * (y2 - y4) + (x4 - x2) * (y1 - y3));
                Xc[i + nim * j] = T(0.25)* (x1 + x2 + x3 + x4);
                Yc[i + nim * j] = T(0.25)* (y1 + y2 + y3 + y4); 
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout<< "Cells volume and center, double loop: " << elapsed << "mcs" << std::endl;
    }
    
    {
        auto start = std::chrono::steady_clock::now();
        // #pragma omp for
        for( int id{}; id< nm; id++)
        {
            size_t i = id % nim;
            size_t j = (id-i) / nim;
            auto id1 = i + ni * j;
            auto id2 = i + 1 + ni * j;
            auto id3 = i + 1 + ni * (j + 1);
            auto id4 = i + ni * (j + 1);

            auto x1 = X[id1];
            auto y1 = Y[id1];
            auto x2 = X[id2];
            auto y2 = Y[id2];
            auto x3 = X[id3];
            auto y3 = Y[id3];
            auto x4 = X[id4];
            auto y4 = Y[id4];
            V[ id]  = T(0.5) * ((x1 - x3) * (y2 - y4) + (x4 - x2) * (y1 - y3));
            Xc[id] = T(0.25)* (x1 + x2 + x3 + x4);
            Yc[id] = T(0.25)* (y1 + y2 + y3 + y4); 
            
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Cells volume and center single loop:" << elapsed << "mcs" << std::endl;
    }

    gbs::points_vector<T,2> pts_c(nm);
    for( int id{}; id< nm; id++)
        pts_c[id] = {Xc[id],Yc[id]};
    // Compute faces areas

    vector<T> S_EW_x((ni)*(njm)),S_EW_y((ni)*(njm)),S_SN_x((nim)*(nj)),S_SN_y((nim)*(nj));
    {
        auto start = std::chrono::steady_clock::now();
        for(size_t j{}; j < njm; j++)
        {
            for(size_t i{}; i < ni; i++)
            {
                auto id1 = i + ni * j;
                auto id2 = i + 1 + ni * j;
                S_EW_x[id1] = Y[id2] - Y[id2];
                S_EW_y[id1] = X[id1] - X[id2];
            }
        }
        for(size_t i{}; i < nim; i++)
        {
            for(size_t j{}; j < nj; j++)
            {
                auto id1 = j + nj * i;
                auto id2 = j + 1 + nj * i;
                S_SN_x[id1] = Y[id1] - Y[id2];
                S_SN_y[id1] = X[id2] - X[id2];
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Cells Surf/Normals:" << elapsed << "mcs" << std::endl;
    }

    SolverData4<T> E{nip, njp};
    SolverData4<T> F{nip, njp};
    SolverData4<T> E_EW{ni, njm}, E_SN{nim, nj};
    SolverData4<T> F_EW{ni, njm}, F_SN{nim, nj};

    auto & [U1, U2, U3, U4] = U.d;
    auto & [E1, E2, E3, E4] = E.d;
    auto & [F1, F2, F3, F4] = F.d;
    vector<T> gamma(np, 1.4);
    {
        auto start = std::chrono::steady_clock::now();
// #pragma omp for
        for( int id{}; id<np ; id++)
        {
            auto u = U2[id]/U1[id];
            auto v = U3[id]/U1[id];
            auto p = (gamma[id]-1) * ( U4[id] - 0.5 * (u*u+v*v));

            E1[id] = U2[id];
            E2[id] = U2[id]/U1[id] + p;
            E3[id] = U2[id]*U3[id]/U1[id];
            E4[id] = (U4[id] + p) * u; 

            F1[id] = U3[id];
            F2[id] = U2[id]*U3[id]/U1[id];
            F3[id] = U3[id]/U1[id] + p;
            F4[id] = (U4[id] + p) * v; 
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "E/F cells centers single loop:" << elapsed << "mcs" << std::endl;
    }

    {
        auto start = std::chrono::steady_clock::now();
        for(size_t j{}; j < njp; j++)
        {
            for( size_t id{nip*j}; id < nip * (j+1); id++)
            {
                // auto id = i + nip * j;
                auto u = U2[id]/U1[id];
                auto v = U3[id]/U1[id];
                auto p = (gamma[id]-1) * ( U4[id] - 0.5 * (u*u+v*v));

                E1[id] = U2[id];
                E2[id] = U2[id]/U1[id] + p;
                E3[id] = U2[id]*U3[id]/U1[id];
                E4[id] = (U4[id] + p) * u; 

                F1[id] = U3[id];
                F2[id] = U2[id]*U3[id]/U1[id];
                F3[id] = U3[id]/U1[id] + p;
                F4[id] = (U4[id] + p) * v; 
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "E/F cells centers double loop:" << elapsed << "mcs" << std::endl;
    }

    // vector<size_t> ids_jp(njp);
    // std::generate(ids_jp.begin(), ids_jp.end(), [j=0]() mutable {return j++;});
    // {
    //     auto start = std::chrono::steady_clock::now();
    //     // size_t id{};
    //     std::for_each(
    //         // std::execution::par,
    //         ids_jp.begin(), ids_jp.end(),
    //         [&](auto j)
    //         {
    //             for( size_t id{nip * j}; id < nip * (j+1); id++)
    //             {
    //                 // auto id = i + nip * j;
    //                 auto u = U2[id]/U1[id];
    //                 auto v = U3[id]/U1[id];
    //                 auto p = (gamma[id]-1) * ( U4[id] - 0.5 * (u*u+v*v));

    //                 E1[id] = U2[id];
    //                 E2[id] = U2[id]/U1[id] + p;
    //                 E3[id] = U2[id]*U3[id]/U1[id];
    //                 E4[id] = (U4[id] + p) * u; 

    //                 F1[id] = U3[id];
    //                 F2[id] = U2[id]*U3[id]/U1[id];
    //                 F3[id] = U3[id]/U1[id] + p;
    //                 F4[id] = (U4[id] + p) * v; 
    //             }
    //         }
    //     );
    //     auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
    //     std::cout << "E/F cells centers for_each vectorized loop:" << elapsed << "mcs" << std::endl;
    // }

    auto & [E1_EW, E2_EW, E3_EW, E4_EW] = E_EW.d;
    auto & [F1_EW, F2_EW, F3_EW, F4_EW] = F_EW.d;

    {
        auto start = std::chrono::steady_clock::now();
        for(size_t j{}; j < njm; j++)
        {
            for(size_t i{}; i < ni; i++)
            {
                auto id  = i + ni * j;
                auto id1p= i + nip * j;
                auto id2p= i+1 + nip * j;
                E1_EW[id] = 0.5 * ( E1[id1p] + E1[id2p]);
                F1_EW[id] = 0.5 * ( F1[id1p] + F1[id2p]);
                E2_EW[id] = 0.5 * ( E2[id1p] + E2[id2p]);
                F2_EW[id] = 0.5 * ( F2[id1p] + F2[id2p]);
                E3_EW[id] = 0.5 * ( E3[id1p] + E3[id2p]);
                F3_EW[id] = 0.5 * ( F3[id1p] + F3[id2p]);
                E4_EW[id] = 0.5 * ( E4[id1p] + E4[id2p]);
                F4_EW[id] = 0.5 * ( F4[id1p] + F4[id2p]);
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "E/F_EW cells face double loop:" << elapsed << "mcs" << std::endl;
    }

    auto & [E1_SN, E2_SN, E3_SN, E4_SN] = E_SN.d;
    auto & [F1_SN, F2_SN, F3_SN, F4_SN] = F_SN.d;

    {
        auto start = std::chrono::steady_clock::now();
        for(size_t i{}; i < nim; i++)
        {
            for(size_t j{}; j < nj; j++)
            {
                auto id  = j + nj * i;
                auto id1p= i + nip * j;
                auto id2p= i + nip * (j+1);
                E1_SN[id] = 0.5 * ( E1[id1p] + E1[id2p]);
                F1_SN[id] = 0.5 * ( F1[id1p] + F1[id2p]);
                E2_SN[id] = 0.5 * ( E2[id1p] + E2[id2p]);
                F2_SN[id] = 0.5 * ( F2[id1p] + F2[id2p]);
                E3_SN[id] = 0.5 * ( E3[id1p] + E3[id2p]);
                F3_SN[id] = 0.5 * ( F3[id1p] + F3[id2p]);
                E4_SN[id] = 0.5 * ( E4[id1p] + E4[id2p]);
                F4_SN[id] = 0.5 * ( F4[id1p] + F4[id2p]);
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "E/F_SN cells faces double loop:" << elapsed << "mcs" << std::endl;
    }

    auto & [R1, R2, R3, R4] = R.d;
    {
        auto start = std::chrono::steady_clock::now();
        for(size_t j{}; j < njm; j++)
        {
            for(size_t i{}; i < nim; i++)
            {
                auto id1 = i + ni*j;
                auto id2 = j + nj*i;
                auto id3 = i+1 + ni*j;
                auto id4 = j+1 + nj *i;
                auto S_EW_x_id1=S_EW_x[id1];
                auto S_SN_x_id2=S_SN_x[id2];
                auto S_EW_x_id3=S_EW_x[id3];
                auto S_SN_x_id4=S_SN_x[id4];
                auto S_EW_y_id1=S_EW_y[id1];
                auto S_SN_y_id2=S_SN_y[id2];
                auto S_EW_y_id3=S_EW_y[id3];
                auto S_SN_y_id4=S_SN_y[id4];
                R1[i + j *ni] = -(
                    -S_EW_x_id1 * E1_EW[id1] - S_EW_y_id1 * F1_EW[id1]
                    -S_SN_x_id2 * E1_SN[id2] - S_SN_y_id2 * F1_SN[id2]
                    +S_EW_x_id3 * E1_EW[id3] + S_EW_y_id3 * F1_EW[id3]
                    +S_SN_x_id4 * E1_SN[id4] + S_SN_y_id4 * F1_SN[id4]
                );
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Residual cells  double loop:" << elapsed << "mcs" << std::endl;
    }

    const std::array<T,5> rk_coefs{0.0533,0.1263,0.2375, 0.4414,1.0};



    plot(
        make_structuredgrid_actor(pts, ni, nj)
        , pts_c
    );
}

TEST(tests_euler2D, base_solver)
{
    using T = double;
    size_t ni_{51}, nj_{25};
    // size_t ni_{11}, nj_{5};
    T bump = 0.;
    auto [ pts, ni, nj, n_iso_ksi, n_iso_eth] = build_channel_mesh(ni_, nj_, bump);

    yams::Euler2DSolver<T> solver{pts, ni, nj}; 
    // yams::Euler2DSolver solver{pts, ni, nj}; 

    {
        T u = 123;
        T v = 65;
        T p = 1.53e5;
        T t = 364;
        auto U = solver.U({u, v, t, p}, 1.39);
        auto P = solver.P(U, 1.39);
        ASSERT_NEAR( std::get<0>(P), u, 1e-7);
        ASSERT_NEAR( std::get<1>(P), v, 1e-7);
        ASSERT_NEAR( std::get<2>(P), t, 1e-7);
        ASSERT_NEAR( std::get<3>(P), p, 1e-7);
    }
    {
        auto start = std::chrono::steady_clock::now();
        solver.computeCellsCenterAndVolume();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "computeCellsCenterAndVolume:" << elapsed << "mcs" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        solver.computeCellsNormalsEW();
        solver.computeCellsNormalsSN();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "computeCellsNormals:" << elapsed << "mcs" << std::endl;
    }
    solver.init();
    {
        auto start = std::chrono::steady_clock::now();
        solver.computeCellsFlowInfo();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "computeCellsFlowInfo:" << elapsed << "mcs" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        solver.computeCellFacesFlowInfo();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "computeCellFacesFlowInfo:" << elapsed << "mcs" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        solver.computeCellsFlowInfoEW();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "computeCellsFlowInfoEW:" << elapsed << "mcs" << std::endl;
    }
        solver.evalResidual();
}


#pragma once
#include <gbs/curves>
namespace yams
{

    template <typename T>
    using vector = std::vector<T>;

    template <typename T>
    struct BladeToBladeCurvatureSolverData
    {
        size_t ni,nj;
        vector<T> W;  // Relative velocity
        vector<T> V;  // Absolute velocity
        vector<T> U;  // Blade tangential velocity
        vector<T> M;  // Absolute Mach number
        vector<T> MW; // Relative Mach number
        vector<T> Vm; // Meridional velocity
        vector<T> Vu; // Tangential absolute velocity
        vector<T> Wu; // Tangential relative velocity
        vector<T> TS; // Static temp
        vector<T> TT; // Total temp
        vector<T> PS; // Total press
        vector<T> PT; // Total press
        vector<T> H;  // Enthalpy
        vector<T> S;  // Entropy
        vector<T> R;  // Density
        vector<T> Q;  // Accumulated mass flow
        vector<T> W1; // W + d_theta * d_W / d_th
        vector<T> W2; // W + d_theta * d_W1 / d_th
        vector<T> G1; // d_Vm / d_m
        vector<T> G2; // d_S_ / d_m
        vector<T> G3; // d_r^2Vm / d_m
        BladeToBladeCurvatureSolverData() = default;
        BladeToBladeCurvatureSolverData(size_t ni, size_t nj) : ni{ni}, nj{nj}
        {
            auto n = ni * nj;
            M.resize(n);
            MW.resize(n);
            W.resize(n);
            V.resize(n);
            U.resize(n);
            Vm.resize(n);
            Vu.resize(n);
            Wu.resize(n);
            TS.resize(n);
            TT.resize(n);
            PS.resize(n);
            PT.resize(n);
            H.resize(n);
            S.resize(n);
            R.resize(n);
            Q.resize(n);
            W1.resize(n);
            W2.resize(n);
            G1.resize(n);
            G2.resize(n);
            G3.resize(n);
        }
    };

    template <typename T>
    struct BladeToBladeCurvatureSolverMesh
    {
        size_t ni, nj;
        vector<T> U;
        vector<T> M;
        vector<T> TH;
        vector<T> DTH;
        vector<T> R;
        vector<T> Z;
        vector<T> TAU; // delta r, stream sheet thickness
        vector<T> BT;  // beta = atan()
        vector<T> RTB; // r * tan( beta )
        vector<T> CB;  // Cos( beta )
        vector<T> SB;  // Sin( beta )
        vector<T> PH;  // Phi d_r / d_z
        vector<T> SP;  // Sin( phi )
        vector<T> CP;  // Cos( phi )
        vector<T> G1;  // d_th / d_m
        vector<T> G2;  // d_(r*tan(beta)) / d_m
        // vector<T> G2; // d_r / d_m
        // vector<T> G3; // d2_th / d_m2
        vector<T> G4; // dr / dz
        BladeToBladeCurvatureSolverMesh() = default;
        BladeToBladeCurvatureSolverMesh(size_t ni, size_t nj) : ni{ni}, nj{nj}
        {
            auto n = ni * nj;
            R.resize(n);
            Z.resize(n);
            TH.resize(n);
            DTH.resize(n);
            TAU.resize(n);
            M.resize(n);
            U.resize(n);
            BT.resize(n);
            RTB.resize(n);
            CB.resize(n);
            SB.resize(n);
            PH.resize(n);
            SP.resize(n);
            CP.resize(n);
            G1.resize(n);
            G2.resize(n);
            // G3.resize(n);
            G4.resize(n);
        }
    };

    template <typename T>
    class BladeToBladeCurvatureSolver
    {
    private:
        // Geom info
        BladeToBladeCurvatureSolverData<T> dat;
        BladeToBladeCurvatureSolverMesh<T> msh;
        size_t ni{}, nj{}, jLe{}, jTe{}, i0{};
        vector<size_t> computational_planes_offsets;
        vector<size_t> stream_lines_indices;
        std::shared_ptr<gbs::Curve<T, 2>> stream_line;
        // Stage info
        size_t z_;
        T Mf;
        T omg;
        // Gas prop
        T gamma{1.4};
        T R{287.04};
        T Cp{1004.};
        // Blockage
        size_t max_compressibility_iterations{100};
        T compressibility_rel_tol{1e-3};
        T blockage_mass_flow_reduction_factor{0.999};
        vector<T> compressibility_residual;
        // W newton solve on computation plane
        T eps_newton = 1e-4;
        T tol_rel_mf = 1e-3;
        size_t max_iter_newton{100};
        //**************
        void f_grad(vector<T> &G, const vector<T> &Y, const vector<T> &X);
        // void f_grad2(vector<T> &G, const vector<T> &Y, const vector<T> &X);
        void computeMeshData();
        void computeFlowData();
        void computeVm();
        void computeW(size_t id_start);
        void computeW(size_t id_start, T Wi);
        T evalMassFlow(size_t id_start);

    public:
        std::array<size_t,2> dimensions() const {return {ni,nj};}
        const auto &mesh() const { return msh; }
        const auto &data() const { return dat; }
        size_t periodicity() const { return z_; }
        T fullBladePassageMassFlow() const { return Mf * z_;}
        T rotationSpeed() { return omg;}
        size_t leadingEdgeIndex() const {return jLe;}
        size_t trailingEdgeIndex() const {return jTe;}
        const auto & meridionalStreamLine() const {return *stream_line;}
        void setFullBladePassageMassFlow(T mf) { Mf = mf / z_; }
        void setPeriodicity(size_t n)
        {
            auto mf =fullBladePassageMassFlow();
            z_ = n;
            setFullBladePassageMassFlow(mf);
        }
        void setRotationSpeed(T rdPerSec) { omg = rdPerSec; }
        void setCp(T cp){ Cp = cp;}
        void setGasConstant(T RG){ R = RG;}
        void setGamma(T Gamma){gamma=Gamma;}
        void setPtIn(T Pt){for( size_t i {}; i <ni; i++) dat.PT[i]=Pt;}
        void setTtIn(T Tt){for( size_t i {}; i <ni; i++) dat.TT[i]=Tt;}
        // void setLeadingEdgeIndex(size_t j){jLe = j;}
        // void setTrailingEdgeIndex(size_t j){jTe = j;}

        void computeW();
        
        BladeToBladeCurvatureSolver() = default;
        BladeToBladeCurvatureSolver(
            const gbs::points_vector<T,2> &pts, size_t n_computation_planes,
            const std::shared_ptr<gbs::Curve<T, 2>> &stream_line,
            size_t n_blades,
            size_t j_le,
            size_t j_te
            );
    };

    template <typename T>
    BladeToBladeCurvatureSolver<T>::BladeToBladeCurvatureSolver(
        const gbs::points_vector<T,2> &pts, size_t n_computation_planes,
        const std::shared_ptr<gbs::Curve<T, 2>> &stream_line,
        size_t n_blades, size_t j_le, size_t j_te) : stream_line{stream_line}, z_{n_blades}, jLe{j_le}, jTe{j_te}
    {
        auto n = pts.size();
        if (n % n_computation_planes != 0)
        {
            throw std::invalid_argument("Wrong array size.");
        }
        nj = n_computation_planes;
        ni = n / n_computation_planes;
        i0 = (ni-1)/2+1;
        dat = BladeToBladeCurvatureSolverData<T>{ni, nj};
        std::fill(dat.W.begin(), dat.W.end(), 10.);
        msh = BladeToBladeCurvatureSolverMesh<T>{ni, nj};
        std::fill(msh.TAU.begin(), msh.TAU.end(), .1);
        std::fill(dat.R.begin(), dat.R.end(), 1.255);
        std::fill(dat.TS.begin(), dat.TS.end(), 288.15);
        std::fill(dat.TT.begin(), dat.TT.end(), 288.15);
        std::fill(dat.PT.begin(), dat.PT.end(), 1e5);
        std::fill(dat.H.begin(), dat.H.end(), 288.15 * 1004);

        int im = n - 1;
        auto [u1, u2] = stream_line->bounds();
        // #pragma omp parallel for
        for (int i{}; i <= im; i++)
        {
            auto [u, th] = pts[i];
            auto [z, r] = stream_line->value(u);
            auto m = gbs::length(*stream_line, u1, u);
            msh.U[i] = u;
            msh.M[i] = m;
            msh.R[i] = r;
            msh.Z[i] = z;
            msh.TH[i] = th;
        }


        for (size_t j{}; j < nj; j++)
        {
            computational_planes_offsets.push_back(j * ni);
        }
        for(size_t i{}; i < ni; i++)
        {
            stream_lines_indices.push_back(i);
        }
        computeMeshData();
    }

    template <typename T>
    void BladeToBladeCurvatureSolver<T>::f_grad(vector<T> & G, const vector<T> &Y, const vector<T> &X)
    {
        int nm = Y.size();

        for (int j{}; j < nj; j++)
        {
            int jOffset = j * ni;
            int jOffset1 = std::max<int>(0, jOffset - ni);
            int jOffset2 = std::min<int>(nm - ni, jOffset + ni);
            // #pragma omp parallel for
            for (int i{}; i < ni; i++)
            {
                G[i + jOffset] = (Y[i + jOffset2] - Y[i + jOffset1]) / (X[i + jOffset2] - X[i + jOffset1]);
            }
        }
    }

    // template <typename T>
    // void BladeToBladeCurvatureSolver<T>::f_grad2(vector<T> & G, const vector<T> &Y, const vector<T> &X)
    // {
    //     int nm = Y.size() - 1;

    //     for (int j{1}; j < nj - 1; j++)
    //     {
    //         int jOffset = j * ni;
    //         int jOffset1 = std::max<int>(0, jOffset - ni);
    //         int jOffset2 = std::min<int>(nm, jOffset + ni);
    //         // #pragma omp parallel for
    //         for (int i{}; i < ni; i++)
    //         {
    //             auto dX = (X[i + jOffset2] - X[i + jOffset1]);
    //             G[i + jOffset] = (Y[i + jOffset2] - 2 * Y[i + jOffset] + Y[i + jOffset1]) / dX / dX;
    //         }
    //     }

    //     {
    //         int jOffset = 0;
    //         int jOffset1 = ni;
    //         int jOffset2 = 2 * ni;
    //         for (int i{}; i < ni; i++)
    //         {
    //             G[i + jOffset] = 2 * G[i + jOffset1] - G[i + jOffset2];
    //         }
    //     }

    //     {
    //         int jOffset = ni * (nj - 1);
    //         int jOffset1 = ni * (nj - 2);
    //         int jOffset2 = ni * (nj - 3);
    //         for (int i{}; i < ni; i++)
    //         {
    //             G[i + jOffset] = 2 * G[i + jOffset1] - G[i + jOffset2];
    //         }
    //     }
    // }

    template <typename T>
    void BladeToBladeCurvatureSolver<T>::computeMeshData()
    {
        // compute beta
        f_grad(msh.G1, msh.TH, msh.M);
        transform(
            std::execution::par,
            msh.G1.begin(), msh.G1.end(),
            msh.R.begin(),
            msh.BT.begin(),
            [](auto g, auto r)
            { return atan(r * g); });
        transform(
            std::execution::par,
            msh.BT.begin(), msh.BT.end(),
            msh.R.begin(),
            msh.RTB.begin(),
            [](auto beta, auto r)
            { return r * tan(beta); });
        // eval geom grads
        // f_grad(msh.G2, msh.R , msh.M);
        // f_grad2(msh.G3, msh.TH , msh.M);
        f_grad(msh.G2, msh.RTB, msh.M);
        f_grad(msh.G4, msh.R, msh.Z);
        // angles trigonometry
        transform(
            std::execution::par,
            msh.G4.begin(), msh.G4.end(),
            msh.PH.begin(),
            [](auto d_r_d_z)
            { return atan(d_r_d_z); });
        transform(
            std::execution::par,
            msh.PH.begin(), msh.PH.end(),
            msh.SP.begin(),
            [](auto phi)
            { return sin(phi); });
        transform(
            std::execution::par,
            msh.PH.begin(), msh.PH.end(),
            msh.CP.begin(),
            [](auto phi)
            { return cos(phi); });
        transform(
            std::execution::par,
            msh.BT.begin(), msh.BT.end(),
            msh.CB.begin(),
            [](auto beta)
            { return cos(beta); });
        transform(
            std::execution::par,
            msh.BT.begin(), msh.BT.end(),
            msh.SB.begin(),
            [](auto beta)
            { return sin(beta); });
        std::for_each(
            std::execution::par,
            computational_planes_offsets.begin(),
            computational_planes_offsets.end(),
            [ni=ni, &DTH=msh.DTH, &TH=msh.TH](auto id_start)
                { 
                    auto id_end = id_start + ni - 1;
                    for(size_t id{id_start+1}; id <= id_end; id++)
                        DTH[id] = TH[id] - TH[id-1];
                }
        );
    }

    template <typename T>
    void BladeToBladeCurvatureSolver<T>::computeFlowData()
    {
        computeVm();
        f_grad(dat.G1, dat.Vm, msh.M);
        vector<T> Vec(ni*nj);
        std::transform(
            msh.R.begin(), msh.R.end(),
            dat.Vm.begin(),
            Vec.begin(),
            [](T r, T Vm){return r*r*Vm;}
        );
        f_grad(dat.G3, Vec, msh.M);
        auto ga_ = gamma / (gamma - 1);
        for (int j{0}; j < nj; j++)
        {
            int jOffset = j * ni;
            int jOffset2 = jOffset + ni;
            for (int id{jOffset}; id < jOffset2; id++)
            {
                auto U = msh.R[id] * omg;
                auto Wu = msh.SB[id] * dat.W[id];
                auto Vu = Wu + U;
                dat.V[id] = std::sqrt(Vu * Vu + dat.Vm[id] * dat.Vm[id]);
                dat.Vu[id] = Vu;
                dat.Wu[id] = Wu;
                dat.U[id] = U;
            }
        }
        for (int j{1}; j < nj; j++)
        {
            int jOffset = j * ni;
            int jOffset2 = jOffset + ni;
            for (int id{jOffset}; id < jOffset2; id++)
            {
                auto DH = dat.U[id - ni] * dat.Vu[id - ni] - dat.U[id] * dat.Vu[id];
                dat.H[id] = dat.H[id - ni] + DH;
                dat.TT[id] = dat.TT[id - ni] - DH / Cp;
                dat.PT[id] = dat.PT[id - ni] * std::pow(dat.TT[id] / dat.TT[id - ni], ga_); // when applying losses use leading edge
            }
        }
        for (int j{0}; j < nj; j++)
        {
            int jOffset = j * ni;
            int jOffset2 = jOffset + ni;
            for (int id{jOffset}; id < jOffset2; id++)
            {
                dat.M[id] = dat.V[id] / std::sqrt(gamma * R * dat.TS[id]);
                dat.MW[id] = dat.W[id] / std::sqrt(gamma * R * dat.TS[id]);
                dat.TS[id] = dat.TT[id] - dat.V[id] * dat.V[id] / 2 / Cp;
                dat.PS[id] = dat.PT[id] * std::pow(dat.TS[id] / dat.TT[id], ga_);
                dat.R[id] = dat.PS[id] / R / dat.TS[id];
                dat.S[id] = std::log(pow(dat.TS[id] / 288.15, Cp) / std::pow(dat.PS[id] / 1e5, R));
            }
        }
        f_grad(dat.G2, dat.S, msh.M);
    }

    template <typename T>
    void BladeToBladeCurvatureSolver<T>::computeVm()
    {
        std::transform(
            std::execution::par,
            dat.W.begin(), dat.W.end(),
            msh.BT.begin(),
            dat.Vm.begin(),
            [](auto W, auto beta)
            { return W * cos(beta); });
    }

    template <typename T>
    void BladeToBladeCurvatureSolver<T>::computeW()
    {
        for (int i{}; i < max_compressibility_iterations; i++) // TODO find a convergence criteria
        {
            std::for_each(
                std::execution::par,
                computational_planes_offsets.begin(),
                computational_planes_offsets.end(),
                [this](auto id_start)
                { this->computeW(id_start); });
            // Update energy and compressibility
            compressibility_residual = {dat.R};
            computeFlowData();
            // Eval flow convergence
            std::transform(
                std::execution::par,
                compressibility_residual.begin(), compressibility_residual.end(),
                dat.R.begin(),
                compressibility_residual.begin(),
                [](auto rho_prev, auto rho)
                { return std::abs(rho_prev - rho) / rho; });
            auto residual = std::reduce(
                                std::execution::par, compressibility_residual.begin(), compressibility_residual.end()) /
                            compressibility_residual.size();
            // std::cout << residual << std::endl;

            if (std::isnan(residual))
            {
                i = 0;
                std::fill(dat.W.begin(), dat.W.end(), 10.);
                std::fill(dat.G1.begin(), dat.G1.end(), 0.);
                std::fill(dat.G2.begin(), dat.G2.end(), 0.);
                std::fill(dat.R.begin(), dat.R.end(), dat.R.front()); // TODO improve reset
                std::cout << "Reducing Mass-Flow from " << Mf * z_ << " to ";
                Mf *= blockage_mass_flow_reduction_factor;
                std::cout << Mf * z_ << std::endl;
            }
            else
            {
                if (residual < compressibility_rel_tol)
                    break;
            }
        }
    }

    template <typename T>
    T BladeToBladeCurvatureSolver<T>::evalMassFlow(size_t id_start)
    {
        auto id_end = id_start + ni - 1;
        dat.Q[0] = 0.;
        vector<T> q(ni);
        for (size_t id{id_start}; id <= id_end; id++)
        {
            q[id - id_start] = dat.R[id] * dat.W[id] * msh.CB[id] * msh.CP[id] * msh.TAU[id];
        }
        
        for (size_t i{1}; i < ni; i++)
        {
            dat.Q[id_start + i] = dat.Q[id_start + i -1] + (msh.TH[id_start + i] - msh.TH[id_start + i - 1]) * (q[i] + q[i - 1]) / 2;
        }
        return dat.Q[id_end];
    }

    template <typename T>
    void BladeToBladeCurvatureSolver<T>::computeW(size_t id_start)
    {
        auto err_mf = tol_rel_mf * 10.;
        size_t id0_integration = id_start;
        T Wi = dat.W[id0_integration];
        size_t count{};
        while (err_mf > tol_rel_mf && count < max_iter_newton)
        {
            computeW(id_start, Wi - eps_newton);
            auto mf_pre = evalMassFlow(id_start);
            computeW(id_start, Wi);
            auto mf_ = evalMassFlow(id_start);
            Wi = Wi - eps_newton * (mf_ - Mf) / (mf_ - mf_pre);
            err_mf = fabs(mf_ - Mf) / Mf;
            count++;
        }
    }

    template <typename T>
    void BladeToBladeCurvatureSolver<T>::computeW(size_t id_start, T Wi)
    {

        auto id_end = id_start + ni - 1;
        dat.W[id_start + i0] = Wi;
//*************** unsuccessful aggressive vectorization *****************//
        // const auto &W   = dat.W;
        // const auto &CB  = msh.CB;
        // const auto &SB  = msh.SB;
        // const auto &SP  = msh.SP;
        // const auto &MG2 = msh.G2;
        // const auto &DG1 = dat.G1;
        // const auto &DG2 = dat.G2;
        // const auto &TS  = dat.TS;
        // const auto &R   = msh.R;

        // for (size_t id{id_start + 1}; id <= id_end; id++)
        // {
        //     auto dth = msh.DTH[id];

        //     auto W = dat.W[id - 1];
        //     auto P_ = CB[id - 1] * CB[id - 1] * MG2[id - 1];
        //     auto Q_ = SB[id - 1] * DG1[id - 1] * R[id - 1] + 2 * omg * R[id - 1] * SP[id - 1] * CB[id - 1];
        //     auto R_ = TS[id - 1] * R[id - 1] * SB[id - 1] * CB[id - 1] * DG2[id - 1];

        //     auto R = W * P_ + Q_ + R_ / W;
        //     dat.W1[id] = W + dth *R;

        // }
        // for (size_t id{id_start + 1}; id <= id_end; id++)
        // {
        //     auto dth = msh.DTH[id];
            
        //     auto W1 = dat.W1[id];
        //     auto W = dat.W[id - 1];
        //     auto P_ = CB[id] * CB[id] * MG2[id];
        //     auto Q_ = SB[id] * DG1[id] * R[id] + 2 * omg * R[id] * SP[id] * CB[id];
        //     auto R_ = TS[id] * R[id] * SB[id] * CB[id] * DG2[id];

        //     auto R = W1 * P_ + Q_ + R_ / W1;

        //     dat.W2[id] = W + dth *R;
        // }

        // for (size_t id{id_start + 1}; id <= id_end; id++)
        // {
        //     dat.W[id] = ( dat.W1[id] + dat.W2[id] ) /2;
        // }
//**********************************************//
        auto f = [omg = omg,
                  &CB = msh.CB,
                  &SB = msh.SB,
                  &SP = msh.SP,
                  &MG2 = msh.G2,
                  &DG1 = dat.G1,
                  &DG2 = dat.G2,
                  &TS = dat.TS,
                  &R = msh.R
                ](auto W, auto id) // <- lambda to enable vectorization
        {
            auto P_ = CB[id] * CB[id] * MG2[id];
            auto Q_ = SB[id] * DG1[id] * R[id] + 2 * omg * R[id] * SP[id] * CB[id];
            auto R_ = TS[id] * R[id] * SB[id] * CB[id] * DG2[id];

            return W * P_ + Q_ + R_ / W;
        };

        for (size_t id{id_start + i0 + 1}; id <= id_end; id++)
        {
            auto dth = msh.DTH[id];
            dat.W1[id] = dat.W[id - 1] + dth * f(dat.W[id - 1], id - 1);
            dat.W[id] = (dat.W1[id] + (dat.W[id - 1] + dth * f(dat.W1[id], id))) / 2;
        }

        if(i0>0)
        {
            int id_min = id_start;
            for (int id = id_start + i0 - 1; id >= id_min; id--)
            {
                auto dth = msh.DTH[id];
                dat.W1[id] = dat.W[id + 1] - dth * f(dat.W[id + 1], id + 1);
                dat.W[id] = (dat.W1[id] + (dat.W[id + 1] - dth * f(dat.W1[id], id))) / 2;
            }
        }

        // std::transform(
        //     std::next(stream_lines_indices.begin()), stream_lines_indices.end(),
        //     std::next(dat.W.begin(),id_start + 1),
        //     [Wmin = T{1.},id_start, &msh = this->msh, &dat = this->dat, &f](size_t i)
        //     {
        //         auto id = id_start + i;
        //         auto dth = msh.DTH[id];
        //         dat.W1[id] = dat.W[id - 1] + dth * f(dat.W[id - 1], id - 1);
        //         return  std::max(
        //             (dat.W1[id] + (dat.W[id - 1] + dth * f(dat.W1[id], id))) / 2,
        //             Wmin
        //         );
        //     }
        // );

    }

    }
}
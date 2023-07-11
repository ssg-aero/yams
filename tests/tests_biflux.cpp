#include <vector>
#include <gtest/gtest.h>
#include <gbs/curves>
#include <gbs-render/vtkGbsRender.h>
#include <gbs-mesh/tfi.h>
#include <gridrender.h>
#include <gridreader.h>
#include <memory>
#include <ranges>
#include <boost/math/tools/roots.hpp>

template<typename T>
using vector = std::vector<T>;
template<typename T, size_t d>
using array = std::array<T,d>;

template<std::floating_point T>
struct MeridianMesh
{
    // const size_t dim{9};
    array<vector<T>, 10> vectors;
    vector<int> iB;
    vector<T> omg;
    size_t ni{};
    size_t nj{};
    void init(const array<T, 10 - 2> &values = array<T, 10 - 2>{})
    {
        assert(std::size(values) == std::size(vectors) - 2);
        auto n{vectors.size()};
        for (size_t i{2}; i < n; i++)
        {
            std::ranges::fill(vectors[i], values[i - 2]);
        }
    }
    /// @brief Radius
    /// @return 
    auto & r(){return vectors[0];}
    const auto & r() const {return vectors[0];}
    /// @brief Axial position
    /// @return 
    auto & z(){return vectors[1];}
    const auto & z() const {return vectors[1];}
    /// @brief Curvilinear abscissa, stream direction
    /// @return 
    auto & m(){return vectors[2];}
    const auto & m() const {return vectors[2];}
    /// @brief Curvilinear abscissa, spans direction
    /// @return 
    auto & l(){return vectors[3];}
    const auto & l() const{return vectors[3];}
    /// @brief Computational plane angle, atan2(dz, dr)
    /// @return 
    auto & g(){return vectors[4];}
    const auto & g() const{return vectors[4];}
    /// @brief Stream line angle, atan2(dr, dz)
    /// @return 
    auto & p(){return vectors[5];}
    const auto & p() const{return vectors[5];}
    /// @brief stream sheet lean angle, atan2(rdd, dq)
    /// @return 
    auto & e(){return vectors[6];}
    /// @brief Stream line curvature
    /// @return 
    auto & c(){return vectors[7];}
    const auto & c() const {return vectors[7];}
    auto & ksi(){return vectors[8];}
    const auto & ksi() const {return vectors[8];}
    auto & cgp(){return vectors[9];}
    const auto & cgp() const {return vectors[9];}
    /// @brief Allocation
    /// @param ni_ 
    /// @param nj_ 
    void resize(size_t ni_, size_t nj_);
    size_t index(size_t i, size_t j){return j + nj * i;}
    /// @brief Write point set to mesh
    /// @param pts 
    /// @param nj 
    void copy_pts(const vector<array<T,2>> &pts, size_t nj);
    void compute_metrics();
};

template <std::floating_point T>
void MeridianMesh<T>::resize(size_t ni_, size_t nj_)
{
    ni = ni_;
    nj = nj_;
    std::ranges::for_each(vectors, [n = ni * nj](auto &V)
                          { V.resize(n); });
    iB.resize(ni,-1);
    omg.resize(ni,0.);
}

template <std::floating_point T>
void MeridianMesh<T>::copy_pts(const vector<array<T, 2>> &pts, size_t nj)
{
    size_t ni = pts.size() / nj;
    assert(ni!=0);
    resize(ni, nj);
    std::ranges::transform(pts, r().begin(), [](const auto &pt)
                           { return pt[1]; });
    std::ranges::transform(pts, z().begin(), [](const auto &pt)
                           { return pt[0]; });
}

template <std::floating_point T>
void MeridianMesh<T>::compute_metrics()
{
    auto &r_ = r();
    auto &z_ = z();

    auto &m_ = m();
    for (size_t i = 1; i < ni; i++)
    {
        for (size_t j = 0; j < nj; j++)
        {
            auto id = index(i, j);
            auto idm = index(i - 1, j);
            m_[id] = hypot(z_[id] - z_[idm], r_[id] - r_[idm]) + m_[idm];
        }
    }

    auto &l_ = l();
    for (size_t i = 0; i < ni; i++)
    {
        for (size_t j = 1; j < nj; j++)
        {
            auto id = index(i, j);
            auto idm = index(i, j - 1);
            l_[id] = hypot(z_[id] - z_[idm], r_[id] - r_[idm]) + l_[idm];
        }
    }

    size_t i_max = ni - 1;
    size_t j_max = nj - 1;

    auto &p_ = p();
    for (size_t i = 0; i < ni; i++)
    {
        for (size_t j = 0; j < nj; j++)
        {
            auto idp = index(std::min<int64_t>(i + 1, i_max), j);
            auto id = index(i, j);
            auto idm = index(std::max<int64_t>(i - 1, 0), j);
            auto dz = z_[idp] - z_[idm];
            auto dr = r_[idp] - r_[idm];
            p_[id] = atan2(dr, dz);
        }
    }

    auto &g_ = g();
    for (size_t i = 0; i < ni; i++)
    {
        for (size_t j = 0; j < nj; j++)
        {
            auto idp = index(i, std::min<int64_t>(j + 1, j_max));
            auto id = index(i, j);
            auto idm = index(i, std::max<int64_t>(j - 1, 0));
            auto dz = z_[idp] - z_[idm];
            auto dr = r_[idp] - r_[idm];
            g_[id] = atan2(dz, dr);
        }
    }

    std::ranges::transform(g(), p(), cgp().begin(), [](auto g_, auto p_){return std::cos(g_+p_);});

    auto &c_ = c();
    for (size_t i = 1; i < i_max; i++)
    {
        for (size_t j = 0; j < nj; j++)
        {
            auto idp = index(i + 1, j);
            auto id = index(i, j);
            auto idm = index(i - 1, j);
            auto dzp = z_[idp] - z_[id];
            auto drp = r_[idp] - r_[id];
            auto dzm = z_[id] - z_[idm];
            auto drm = r_[id] - r_[idm];
            auto dm  = m_[idp] - m_[idm] ;
            c_[id] = 2 * (atan2(drp, dzp)-atan2(drm, dzm)) / dm;
        }
    }
    // Extrapolate bounds values
    for (size_t j = 0; j < nj; j++)
    {
        c_[j] = 2 * c_[j + nj] - c_[j + 2 * nj];
        c_[j + i_max * nj] = 2 * c_[j + nj * (i_max-1)] - c_[j + nj * (i_max-2)];
    }

}

template<std::floating_point T>
struct MeridianData
{
    size_t ni, nj;
    vector<T> mf;
    vector<T> geom_residual;
    T compressibility_residual{};
    array<vector<T>, 13> vectors;
    auto &Vm(){return vectors[0];}
    const auto &Vm() const{return vectors[0];}
    auto &Vu(){return vectors[1];}
    const auto &Vu() const{return vectors[1];}
    auto &G() {return vectors[2];}
    auto &J() {return vectors[3];}
    auto &K() {return vectors[4];}
    auto &Tt(){return vectors[5];}
    auto &Pt(){return vectors[6];}
    auto &Ts(){return vectors[7];}
    auto &Ps(){return vectors[8];}
    auto &V() {return vectors[9];}
    auto &Rho() {return vectors[10];}
    auto &q() {return vectors[11];}
    const auto &q() const{return vectors[11];}
    auto &M() {return vectors[12];}
    const auto &M() const{return vectors[12];}
    void resize(size_t ni_, size_t nj_);
    void init(T Vm, T Vu, T Tt = 288., T Pt = 1.01325e5);
};

template <std::floating_point T>
void MeridianData<T>::resize(size_t ni_, size_t nj_)
{
    ni = ni_;
    nj = nj_;
    std::ranges::for_each(vectors, [n = ni * nj](auto &vec)
                          { vec.resize(n, T{}); });
    mf.resize(ni);
    geom_residual.resize(ni);
}

template <std::floating_point T>
void MeridianData<T>::init(T Vm_, T Vu_, T Tt_, T Pt_)
{
    std::ranges::fill(Vm(), Vm_);
    std::ranges::fill(Vu(), Vu_);
    std::ranges::fill(Tt(), Tt_);
    std::ranges::fill(Pt(), Pt_);
}

template<std::floating_point T>
class gasModel
{
    public:
    const T r;
    gasModel(T r_) : r{r_} {}
    virtual T  H(T t) const = 0;
    virtual T cp(T t) const = 0;
    T cv(T t) const{
        return cp(t) - r;
    }
    T cp(T Tt, T V, T factor = 2., size_t max_iter = 100) const{
        // boost::uintmax_t max_iter{100};
        auto guess = cp(Tt);
        bool rising{true};
        const int digits = std::numeric_limits<T>::digits; 
        int get_digits = digits - 3; 
        boost::math::tools::eps_tolerance<T> tol(get_digits);
        auto result = boost::math::tools::bracket_and_solve_root(
            [this, Tt, V](T cp_) {
                auto Ts = Tt - V*V / (2*cp_);
                return cp_ - this->cp(Ts);
            }, guess, factor, rising, tol, max_iter
        );
        return result.first + (result.second - result.first)/2;
    }
    T gamma(T t) const{
        return cp(t) / cv(t);
    }
    T rho(T t, T p) const{
        return p / ( r * t );
    }
    T a(T t) const { 
        return std::sqrt( sqa(t));
    }
    T sqa(T t) const { 
        return r * gamma(t) * t;
    }
    T t_from_H(T H)
    {
        auto t_guess = H / cp(288.);
        auto t_min = 0.1*t_guess;
        auto t_max = 10*t_guess;
        const int digits = std::numeric_limits<T>::digits;  
        return boost::math::tools::newton_raphson_iterate(
            [this, H](auto t){return std::make_tuple(H - this->H(t),this->cp(t));},
            t_guess, t_min, t_max, digits
        );
    }
};

template<std::floating_point T>
class airGp : public gasModel<T>
{
    public:
    airGp() : gasModel<T>{287.04} {}
    T H(T t) const override{
        auto Tr = 3090 / t;
        return  this->r * (3.5 * t - 1.4e-5 * t*t + 7.464e-9 * t*t*t + 3090 / ( std::exp(Tr) - 1 ));
    }
    T cp(T t) const override{
        auto Tr = 3090 / t;
        auto a  = exp(Tr);
        return this->r * ( 3.5 + t * (-2.8e-5 + 2.24e-8 * t) + Tr*Tr * a / ((a - 1)*(a - 1)));
    }
};

template<std::floating_point T>
class MeridianSolver
{
private:

public:
    void integrate_Vm(size_t i);
    void compute_gradients();
    void compute_tangential_flow();
    void compute_fluid_prop();
    T compute_mass_flow(size_t i, T Vmj0);
    void update_streams(size_t i, size_t i0);
    void update_streams(size_t i0);
    void solve_Vmj0(size_t i, size_t max_iter = 100, T factor = 1.);
    void solve_Vm(size_t max_iter = 100, T factor = 1.);
    void solve();
    size_t ni, nj;
    size_t j0;
    MeridianMesh<T> mesh;
    MeridianData<T> data;
    std::unique_ptr<gasModel<T>> gm;
    /////////// public
    auto n_spans() const {return ni;}
    auto n_streams() const {return nj;}
    const auto& msh() const { return mesh;}
    const auto& dat() const { return data;}
    MeridianSolver(const vector<array<T,2>> &pts, size_t nj, T Vm0=100., T Vu0=0., T Tt_in=288.15, T Pt_in = 1.01325e5);
    void integrate_Vm();
    void setMassFlow(T mf_){std::ranges::fill(data.mf, mf_);}
    ~MeridianSolver() = default;
};

template<std::floating_point T>
MeridianSolver<T>::MeridianSolver(const vector<array<T,2>> &pts, size_t nj, T Vm0, T Vu0, T Tt_in, T Pt_in) : gm{std::make_unique<airGp<T>>()}
{
    mesh.copy_pts(pts, nj);
    mesh.init();
    mesh.compute_metrics();

    data.resize(mesh.ni, mesh.nj);

    this->ni = mesh.ni;
    this->nj = mesh.nj;

    T s = 0.;
    this->j0 = (this->nj-1) * s;

    data.init(Vm0, Vu0, Tt_in, Pt_in);

    compute_tangential_flow();
    compute_fluid_prop();
    compute_gradients();

}

template<std::floating_point T>
void MeridianSolver<T>::integrate_Vm(size_t i)
{
    const vector<T> &G  = data.G();   
    const vector<T> &J  = data.J();
    const vector<T> &K  = data.K();
    const vector<T> &l  = mesh.l();
    vector<T> &Vm = data.Vm();

    auto f = [](T Vm, T G, T J, T K)
    {
        return (G * Vm * Vm + J * Vm + K) / (2 * Vm);
    };

    auto rk4_step = [&f](T Vm, T G, T J, T K, T dl)
    {
        T k1 = dl * f(Vm, G, J, K);
        T k2 = dl * f(Vm + 0.5 * k1, G, J, K);
        T k3 = dl * f(Vm + 0.5 * k2, G, J, K);
        T k4 = dl * f(Vm + k3, G, J, K);

        return Vm + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    };

    for (int j = j0 + 1; j < nj; j++)
    {
        int id = i * nj + j;
        int idm= i * nj + j - 1;
        T dl = l[id] - l[idm];
        Vm[id] = rk4_step(Vm[idm], G[idm], J[idm], K[idm], dl);
    }

    for(int j = j0 - 1 ; j >= 0; j--)
    {
        int id = i * nj + j;
        int idp= i * nj + j + 1;
        T dl = l[idp] - l[id];
        Vm[id] = rk4_step(Vm[idp], G[idp], J[idp], K[idp], dl);
    }
}

template<std::floating_point T>
void MeridianSolver<T>::integrate_Vm()
{
    #ifndef _DEBUG
        #pragma omp parallel for
    #endif
    for(int i = 0 ; i < ni ; i++) //omp needs simple type
    {
        integrate_Vm(i);
    }

}

template<std::floating_point T>
void MeridianSolver<T>::compute_tangential_flow()
{
    auto &Vu = data.Vu();
    auto &r  = mesh.r();
    vector<T> rVu(nj);
    std::transform(r.begin(), std::next(r.begin(), nj), Vu.begin(), rVu.begin(), std::multiplies{});
    #ifndef _DEBUG
        #pragma omp parallel for
    #endif
    for(int i = 1 ; i < ni ; i++)//omp needs simple type
    {
        std::transform(
            rVu.begin(), rVu.end(),
            std::next(r.begin(), i * nj),
            std::next(Vu.begin(), i * nj),
            std::divides{}
        );
    }
}

template<std::floating_point T>
void MeridianSolver<T>::compute_fluid_prop()
{
    auto &Tt       = data.Tt();
    auto &Pt       = data.Pt();
    auto &V        = data.V();
    const auto &Vu = data.Vu();
    const auto &Vm = data.Vm();
    const auto &r  = mesh.r();
    const auto & ksi= mesh.ksi();
    const auto &iB = mesh.iB;
    const auto &omg= mesh.omg;
    // Compute V
    std::ranges::transform(Vu, Vm, V.begin(),[](auto a, auto b) { return std::hypot(a, b); });

    // Compute Tt
    int i_le{-1};
    #ifndef _DEBUG
        #pragma omp parallel for
    #endif
    for (int i = 1; i < ni; i++) // omp needs simple type
    {
        if (iB[i] >= 0)
        {
            if(i_le<0) i_le = i-1;
            auto omg_ = omg[i];
            for (int j = 0; j < nj; j++)
            {
                auto id = j + i * nj;
                auto idm = j + (i - 1) * nj;
                auto dH = omg_ * (r[id] * Vu[id] - r[idm] * Vu[idm]);
                auto H = gm->H(Tt[idm]);
                Tt[id] = gm->t_from_H(H);
                auto id_le = j + i_le * nj;
                auto V_ = V[id_le];
                Pt[id] = Pt[id_le] - ksi[id] * V_ * V_;
            }
        }
        else
        {
            if(i_le>=0) i_le = -1;
            for (int j = 0; j < nj; j++)
            {
                auto id = j + i * nj;
                auto idm = j + (i - 1) * nj;
                Tt[id] = Tt[idm];
                auto V_ = V[idm];
                Pt[i] = Pt[idm] - ksi[id] * V_ * V_;
            }
        }
    }

    auto &Ts       = data.Ts();
    auto &Ps       = data.Ps();
    auto &Rho      = data.Rho();
    auto &M        = data.M();
    auto n = ni * nj;
    T compress_residual{};
    #ifndef _DEBUG
        #pragma omp parallel for
    #endif
    for (int id = 0; id < n; id++) // omp needs simple type
    {
        auto V_ = V[id];
        auto Ts_ = Ts[id] = Tt[id] - V_ * V_ / (2 * gm->cp(Tt[id], V_));
        Ts[id] = Ts_;
        auto gamma = gm->gamma(Ts_);
        auto a = gm->a(Ts_);
        auto M_ = V_ / a;
        auto Ps_ = Ps[id] = Pt[id]*std::pow(1 + 0.5 * (gamma - 1) * M_ * M_, -gamma / (gamma - 1));
        auto rho_= gm->rho(Ts_, Ps_);
        auto res_=std::abs(rho_-Rho[id]);
        compress_residual = std::max(compress_residual, res_);
        Rho[id] = rho_;
        M[id]   = M_;
    }
    data.compressibility_residual = compress_residual;
}

template<std::floating_point T>
void MeridianSolver<T>::compute_gradients()
{
    auto &G = data.G();
    const auto &cgp = mesh.cgp();
    const auto &c = mesh.c();
    #ifndef _DEBUG
        #pragma omp parallel for
    #endif
    for(int i = 0 ; i < ni ; i++)//omp needs simple type
    {
        for(int j = 0 ; j < nj ; j++)
        {
            auto id = j + i * nj;
            G[id] = cgp[id] * c[id];
        }
    }
}

template<std::floating_point T>
T MeridianSolver<T>::compute_mass_flow(size_t i, T Vmj0)
{
    auto &Vm   = data.Vm();
    auto &q    = data.q();
    const auto &Rho   = data.Rho();
    const auto &r = mesh.r();
    const auto &cgp = mesh.cgp();
    const auto &l = mesh.l();

    T mf_{};

    Vm[j0+nj*i] = Vmj0;
    integrate_Vm(i);

    for(size_t j{1}; j < nj ; j++)
    {
        auto id = j + i * nj;
        auto idm= j-1 + i * nj;
        auto q_ = std::numbers::pi_v<T> * ( Rho[id] * Vm[id] * cgp[id] * r[id] + Rho[idm] * Vm[idm] * cgp[idm] * r[idm] ) * (l[id]-l[idm]);
        mf_ += q_;
        q[id] = mf_;
    }

    return mf_;

}

template<std::floating_point T>
void MeridianSolver<T>::solve_Vmj0(size_t i, size_t max_iter, T factor)
{
    const auto &Vm   = data.Vm();
    const auto &Ts   = data.Ts();
    auto guess =  Vm[j0 + nj * i];
    auto min = 0.;
    auto max = gm->a(Ts[j0 + nj * i]);
    const int digits = std::numeric_limits<T>::digits;  
    auto Vmj0 = boost::math::tools::newton_raphson_iterate(
        [this, mf_=data.mf[i], i, eps=1e-5](T Vmj0) {
            auto mf1 = this->compute_mass_flow(i, Vmj0-0.5*eps);
            auto mf2 = this->compute_mass_flow(i, Vmj0+0.5*eps);
            auto res = mf1 - mf_;
            auto dmf_dVm = (mf2-mf1)/eps;

            return std::make_tuple( res , dmf_dVm ) ;},
        guess, min, max, digits
    );
    compute_mass_flow(i, Vmj0);
}

template<std::floating_point T>
void MeridianSolver<T>::solve_Vm(size_t max_iter, T factor)
{
    #ifndef _DEBUG
        #pragma omp parallel for
    #endif
    for(int i = 0 ; i < ni ; i++)//omp needs simple type
    {
        solve_Vmj0(i, max_iter, factor);
    }
}

template<std::floating_point T>
void MeridianSolver<T>::update_streams(size_t i, size_t i0)
{
    auto &r      = mesh.r();
    auto &z      = mesh.z();
    const auto &l= mesh.l();
    const auto &m= mesh.m();
    const auto &c= mesh.c();
    const auto &q= data.q();
    auto q_begin = std::next(q.begin(), i * nj);
    auto q_end   = std::next(q.begin(), (i+1) * nj);
    T RF = 0.1;
    auto njm = nj-1;
    T geom_residual_{};
    vector<T> delta(nj-1);
    for(size_t j{1}; j < njm ; j++)
    {
        auto q_    = q[i0*nj + j];
        auto low   = std::lower_bound(q_begin, q_end, q_);
        auto id    = std::distance(q.begin(), low);
        auto delta_= (q_ - *low) / q_;
        delta[j-1] = delta_; 
        geom_residual_ = std::max(std::abs(delta_), geom_residual_);
    }

    auto dm_min = std::numeric_limits<T>::max(); 
    for(size_t j{1}; j < njm ; j++)
    {
        auto id = j + i * nj;
        if(i>0)
        {
            auto idm= j + (i-1) * nj;
            dm_min = std::min(
                dm_min,
                (m[id]-m[idm])
            );
        }
        if(i<ni-1)
        {
            auto idp= j + (i+1) * nj;
            dm_min = std::min(
                dm_min,
                (m[idp]-m[id])
            );
        }
    }

    auto kmin = - std::numbers::pi_v<T>*std::numbers::pi_v<T>;
    const auto &M= data.M();
    T A{}; 
    T Mm{};
    auto W = l[(i+1)*nj-1];
    for(size_t j{1}; j < nj ; j++)
    {
        auto id = j + i * nj;
        Mm += 0.5 * ( M[id] + M[id-1] ) * ( l[id] - l[id-1] ); // Integrate
        A = std::max(A, W/( l[id] - l[id-1] ) );               // Mean
    }
    Mm /= W;
    Mm = std::min(0.6, Mm);

    T cst = 8; // 8 -> 4
    auto f = 1 / ( 1 + (1-Mm*Mm) / cst * (W/dm_min) * (W/dm_min));
    for(size_t j{1}; j < njm ; j++)
    {
        auto id = j + i * nj;
        auto delta_ = delta[j-1];
        // r[id] = r[id] + delta_ * f * (r[id+1] - r[id]);
        // z[id] = z[id] + delta_ * f * (z[id+1] - z[id]);
        // r[id] = r[id] + delta_ * f * (r[id+1] - r[id]) *0.5;
        // z[id] = z[id] + delta_ * f * (z[id+1] - z[id])*0.5;
        r[id] = r[id] + delta_ * RF * (r[id+1] - r[id]);
        z[id] = z[id] + delta_ * RF * (z[id+1] - z[id]);
    }
    data.geom_residual[i] = geom_residual_;
    // std::cout<<"i: "<<i << " residual: "<< geom_residual_ << " dm_min: " << dm_min  << " Mm: " << Mm  << " f: " << f << std::endl;
}

template<std::floating_point T>
void MeridianSolver<T>::update_streams(size_t i0)
{
    #ifndef _DEBUG
        #pragma omp parallel for
    #endif
    for(int i = 0 ; i < i0 ; i++)//omp needs simple type
    {
        update_streams(i, i0);
    }
    #ifndef _DEBUG
        #pragma omp parallel for
    #endif
    for(int i = i0+1 ; i < ni ; i++)//omp needs simple type
    {
        update_streams(i, i0);
    }
    // #ifndef _DEBUG
    //     #pragma omp parallel for
    // #endif
    // for(int i = 0 ; i < ni ; i++)//omp needs simple type
    // {
    //     update_streams(i, i0);
    // }
}

template <std::floating_point T>
void MeridianSolver<T>::solve()
{
    size_t nit = 1000;
    size_t nit_compressibility = 20;
    for (size_t it{}; it < nit; it++)
    {
        for (size_t it_sub{}; it_sub < nit_compressibility; it_sub++)
        {
            solve_Vm();
            compute_tangential_flow();
            compute_fluid_prop();
            compute_gradients();
            std::cout<< "\tCompressibility_residual: " << data.compressibility_residual << std::endl;
            if(data.compressibility_residual < 1e-5)
            {
                break;
            }
        }
        update_streams(0);
        mesh.compute_metrics();
        auto max_it = std::ranges::max_element(data.geom_residual);
        std::cout<< "Geometrical residual max: " << *max_it << " on span: " << std::distance(data.geom_residual.begin(), max_it) << std::endl;
        if(*max_it < 1e-5)
        {
            break;
        }
    }
}

using namespace gbs;
using T =double;
using p_curve2d = std::shared_ptr<Curve<T,2>>;
using p_bsc2d   = std::shared_ptr<BSCurve<T,2>>;
using bsc2d     = BSCurve<T,2>;
auto mesh_biflux()
{
    BSCurve<T,2> hub_main{
        { {0.,0.} , {0.4, 0.2}, {0.6, 0.2}, {0.7, 0.2} },
        {0., 1.},
        {4,4},
        3
    };

    BSCurve<T,2> bpr_main{
        { {0.,0.4} , {0.4, 0.4}, {0.6, 0.4}, {0.7, 0.4} },
        {0., 1.},
        {4,4},
        3
    };

    BSCurve<T,2> shr_main{
        { {0.,0.6} , {0.4, 0.6}, {0.6, 0.6}, {0.7, 0.6} },
        {0., 1.},
        {4,4},
        3
    };

    BSCurve<T,2> hub_sec{
        { {0.7,0.4} , {0.9, 0.4}, {1.1, 0.5}, {1.4, 0.5} },
        {0., 1.},
        {4,4},
        3
    };

    BSCurve<T,2> shr_sec{
        { {0.7,0.6} , {0.9, 0.6}, {1.3, 0.7}, {1.4, 0.7} },
        {0., 1.},
        {4,4},
        3
    };

    BSCurve<T,2> hub_pri{
        { {0.7,0.2} , {0.9, 0.2}, {1.2, 0.1}, {1.2, 0.1} },
        {0., 1.},
        {4,4},
        3
    };

    BSCurve<T,2> shr_pri{
        { {0.7,0.4} , {0.9, 0.35}, {1.1, 0.3}, {1.2, 0.3} },
        {0., 1.},
        {4,4},
        3
    };

    BSCurve<T,2> in_main{{ hub_main.begin(), bpr_main.begin(), shr_main.begin()},{0.,0.5,1.},{2,1,2},1};
    BSCurve<T,2> out_main{{ hub_main.end(), bpr_main.end(), shr_main.end()},{0.,0.5,1.},{2,1,2},1};

    BSCurve<T,2> in_main_pri{{ hub_main.begin(), bpr_main.begin()},{0.,1.},{2,2},1};
    BSCurve<T,2> in_main_sec{{ bpr_main.begin(), shr_main.begin()},{0.,1.},{2,2},1};

    BSCurve<T,2> in_pri{{ hub_pri.begin(), shr_pri.begin()},{0.,1.},{2,2},1};
    BSCurve<T,2> out_pri{{ hub_pri.end(), shr_pri.end()},{0.,1.},{2,2},1};

    BSCurve<T,2> in_sec{{ hub_sec.begin(), shr_sec.begin()},{0.,1.},{2,2},1};
    BSCurve<T,2> out_sec{{ hub_sec.end(), shr_sec.end()},{0.,1.},{2,2},1};


    size_t ni_main{13}, ni_pri{10}, ni_sec{13};
    size_t nj_main{11};

    vector<p_curve2d> iso_eth_main{
            std::make_shared<bsc2d>(hub_main),
            std::make_shared<bsc2d>(bpr_main),
            std::make_shared<bsc2d>(shr_main)
        };
    vector<p_curve2d> iso_ksi_main{
            std::make_shared<bsc2d>(in_main),
            std::make_shared<bsc2d>(out_main),
        };
    auto [pts_main, n_main, m_main, n_iso_ksi_main, n_iso_eth_main] = tfi_mesh_2d<T,2,1,1>(
        iso_ksi_main,
        iso_eth_main,
        nj_main,
        ni_main,
        1e-6
    );
    vector<p_curve2d> iso_eth_pri{
            std::make_shared<bsc2d>(hub_pri),
            std::make_shared<bsc2d>(shr_pri)
        };
    vector<p_curve2d> iso_ksi_pri{
            std::make_shared<bsc2d>(in_pri),
            std::make_shared<bsc2d>(out_pri),
        };
    auto [pts_pri, n_pri, m_pri, n_iso_ksi_pri, n_iso_eth_pri] = tfi_mesh_2d<T,2,1,1>(
        iso_ksi_pri,
        iso_eth_pri,
        n_iso_ksi_main[0],
        ni_pri,
        1e-6
    );
    vector<p_curve2d> iso_eth_sec{
            std::make_shared<bsc2d>(hub_sec),
            std::make_shared<bsc2d>(shr_sec)
        };
    vector<p_curve2d> iso_ksi_sec{
            std::make_shared<bsc2d>(in_sec),
            std::make_shared<bsc2d>(out_sec),
        };
    auto [pts_sec, n_sec, m_sec, n_iso_ksi_sec, n_iso_eth_sec] = tfi_mesh_2d<T,2,1,1>(
        iso_ksi_sec,
        iso_eth_sec,
        n_iso_ksi_main[1],
        ni_sec,
        1e-6
    );

    vector<point<T,2>> pts_flux_1, pts_flux_2;
    size_t nj = std::reduce(n_iso_ksi_main.begin(), n_iso_ksi_main.end()) - (n_iso_ksi_main.size() - 1); // Actual number of points
    size_t ni = std::reduce(n_iso_eth_main.begin(), n_iso_eth_main.end()) - (n_iso_eth_main.size() - 1); // Actual number of points
    auto nj_flux_1 = n_iso_ksi_main[0];
    auto nj_flux_2 = n_iso_ksi_main[1];
    auto end = std::next(pts_main.end(),-nj);
 
    for( auto start = pts_main.begin(); start != end ; std::advance(start, nj))
    {
        std::copy(
            start, 
            std::next(start, n_iso_ksi_main[0]), 
            std::back_inserter(pts_flux_1)
        );
        std::copy(
            std::next(start, n_iso_ksi_main[0]-1), 
            std::next(start, nj), 
            std::back_inserter(pts_flux_2)
        );
    }
    std::ranges::copy(pts_pri, std::back_inserter(pts_flux_1));
    std::ranges::copy(pts_sec, std::back_inserter(pts_flux_2));
 

    // plot(
    //     vector<BSCurve<T,2>>{
    //         hub_main,
    //         bpr_main,
    //         shr_main,
    //         hub_sec,
    //         shr_sec,
    //         hub_pri,
    //         shr_pri,
    //         in_main,
    //         out_main,
    //         in_pri,
    //         out_pri,
    //         in_sec,
    //         out_sec
    //     },
    //     // pts_main,
    //     // pts_pri,
    //     // pts_sec
    //     // pts_flux_1, 
    //     pts_flux_2
    // );

    return std::make_tuple(
        pts_main,
        pts_flux_1,
        pts_flux_2,
        nj,
        nj_flux_1,
        nj_flux_2
    );
}

template <std::floating_point T>
inline auto make_vtkStructuredGrid(const vector<T> &z, const vector<T> &r, size_t nj) -> vtkSmartPointer<vtkStructuredGrid>
{
    assert(z.size() == r.size());
    vtkSmartPointer<vtkStructuredGrid> structuredGrid =
        vtkSmartPointer<vtkStructuredGrid>::New();

    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    vector<array<T,2>> g(z.size());
    std::ranges::transform(
        z, r,
        g.begin(),
        [](const auto &x1, const auto &x2){return array<T,2>{x1, x2};}
    );

    std::ranges::for_each( g, [&](const auto &gp) {
            points->InsertNextPoint(gp[0], gp[1], 0);
    });
    auto ni = z.size() / nj;
    structuredGrid->SetDimensions(nj, ni, 1); // i j inverted Grid inner loop is i
    structuredGrid->SetPoints(points);
    return structuredGrid;
}

template <std::floating_point T>
inline auto add_value(const vector<T> &v, vtkStructuredGrid *structuredGrid, const char *name)
{
    vtkSmartPointer<vtkDoubleArray> value_array =
        vtkSmartPointer<vtkDoubleArray>::New();

    value_array->SetName(name);

    std::ranges::for_each(
        v,
        [&](T value) {
            value_array->InsertNextValue(value);
    });

    structuredGrid->GetPointData()->AddArray(value_array);
}

template <std::floating_point T>
void main_to_partial_flow(const vector<T> &main_value, vector<T> &partial_flow_value, size_t i1, size_t i2, size_t j1, size_t j2, size_t nj)
{
    auto end = main_value.end();
    auto nj_loc = j2-j1;
    auto start_main = main_value.begin();
    auto start_partial = partial_flow_value.begin();
    for(size_t i{i1}; i <= i2; i++ )
    {
        auto index_offset = i * nj;
        auto start_loc = std::next(start_main, j1 + index_offset);
        auto end_loc = std::next(start_main, nj_loc);
        std::copy(start_loc,end_loc, std::next(start_partial, index_offset) );
    }
}

TEST(solver, design_novak)
{
 auto [
        pts_main,
        pts_flux_1,
        pts_flux_2,
        nj,
        nj_flux_1,
        nj_flux_2
    ] = mesh_biflux();
    
    // MeridianMesh<T> mesh_main, mesh_flux_1, mesh_flux_2;
    // mesh_main.copy_pts(pts_main, nj);
    // mesh_flux_1.copy_pts(pts_flux_1, nj_flux_1);
    // mesh_flux_2.copy_pts(pts_flux_2, nj_flux_2);

    // mesh_main.init(array<T,6>{});
    // mesh_flux_1.init(array<T,6>{});
    // mesh_flux_2.init(array<T,6>{});

    // mesh_main.compute_metrics();
    // mesh_flux_1.compute_metrics();
    // mesh_flux_2.compute_metrics();

    // MeridianSolver<T> solver_main(pts_main, nj, 120, 30);
    // MeridianSolver<T> solver_main(pts_flux_2,nj_flux_2, 120, 30);
    MeridianSolver<T> solver_main(pts_flux_1,nj_flux_1, 30, 0);
    solver_main.compute_tangential_flow();
    solver_main.compute_fluid_prop();
    solver_main.compute_gradients();
    auto mf = solver_main.compute_mass_flow(0,30);
    solver_main.setMassFlow(mf);
    // solver_main.solve_Vm();
    // solver_main.update_streams(0);
    solver_main.solve();

    auto sgrid = make_vtkStructuredGrid(solver_main.msh().z(), solver_main.msh().r(), solver_main.n_streams());
    add_value(solver_main.msh().c(), sgrid, "Curvature");
    add_value(solver_main.msh().g(), sgrid, "Gamma");
    add_value(solver_main.msh().p(), sgrid, "Phi");
    add_value(solver_main.msh().m(), sgrid, "m");
    add_value(solver_main.msh().l(), sgrid, "l");
    add_value(solver_main.dat().Vm(), sgrid, "Vm");
    add_value(solver_main.dat().Vu(), sgrid, "Vu");
    add_value(solver_main.dat().q(), sgrid, "q");
    yams::plot_vtkStructuredGrid(sgrid, "Vm", true);
    

    auto f_name = "novak_design_biflux.vts";

    vtkNew<vtkXMLStructuredGridWriter> writer;
    writer->SetFileName(f_name);
    writer->SetInputData(sgrid);
    writer->Write();

}

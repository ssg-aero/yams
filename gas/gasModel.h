#pragma once
#include <cmath>
#include <stdexcept>

namespace yams
{
    /// Perfect gas of the meridional solver (curvaturesolver.h): gas constant R and
    /// constant cp, gamma = cp / (cp - R) follows, entropy from (Tref, Pref) like
    /// GridInfo. The blade-to-blade solver uses the same model, so a stream line
    /// analysed after the meridional run sees the same fluid.
    template <typename T>
    struct PerfectGas
    {
        T R{287.04};      ///< gas constant, J/(kg K)
        T cp{1004.5};     ///< specific heat at constant pressure, J/(kg K)
        T Tref{288.};     ///< entropy origin, K
        T Pref{1.01325e5};///< entropy origin, Pa

        /// @throw std::invalid_argument unless R > 0, cp > R (gamma > 1), Tref > 0
        /// and Pref > 0
        void validate() const
        {
            if (!(R > 0 && cp > R && Tref > 0 && Pref > 0))
                throw std::invalid_argument("yams: a perfect gas needs R > 0, cp > R (gamma > 1), Tref > 0 and Pref > 0");
        }

        T gamma() const { return cp / (cp - R); }
        T rho(T Ts, T Ps) const { return Ps / (R * Ts); }
        T sound_speed(T Ts) const { return std::sqrt(gamma() * R * Ts); }
        T entropy(T Ts, T Ps) const { return cp * std::log(Ts / Tref) - R * std::log(Ps / Pref); }
    };

    template <typename T>
    class GasBase
    {
        virtual T cp(T t, T p) = 0;
        virtual T g(T t, T p) = 0;
        virtual T a(T t, T p) = 0;
        virtual T rho(T t, T p) = 0;
        virtual T et(T t, T p, T V) = 0;
    };

    template <typename T>
    class GasPerfectGammaConstant : public GasBase<T>
    {
        T r_;
        T g_;
        public:
        GasPerfectGammaConstant(T r, T g) : r_{r}, g_{g} {}
        T cp(T t, T p) override 
        {
            g_ = g(t, p); // used storage for gas real
            return g_ * r_ / ( g_ - 1);
        }
        T g(T t, T p) override 
        {
            return g_;
        }
        T a(T t, T p) override
        {
            return sqrt(g(t, p) * r_ * t);
        }
        T rho(T t, T p) override
        {
            return p / (r_ * t);
        }
        T et(T t, T p , T V) override
        {
            auto a_ = a(t, p);
            g_ = g(t, p); // used storage for gas real
            return a_ * a_ / (g_ * (g_ - 1)) + 0.5 * V * V;
        }
        T r(){return r_;}
        // T MachFromPr()
    };
    template <typename T>
    class AirGasReal : public GasPerfectGammaConstant<T>
    {

        public:
        AirGasReal(T r) : GasPerfectGammaConstant<T>(r, 1.4) {} // gamma value has no importance, it will be computed later
        T cp(T t, T p) override 
        {
            auto Tr = 3090 / t;
            auto expTr = exp(Tr);
            return this->r() * (
                3.5 + t * (-2.8e-5 + 2.24e-8 * t)
                + Tr*Tr * expTr / ( (expTr - 1)*(expTr - 1) )
            );
        }
        T g(T t, T p) override 
        {
            auto cp_ = cp(t, p);
            return cp_ / ( cp_ - this->r());
        }
        // T MachFromPr()
    };
}
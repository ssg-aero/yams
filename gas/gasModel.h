#pragma once


namespace yams
{
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
        T a(T t, T p)
        {
            return sqrt(g(t, p) * r_ * t);
        }
        T rho(T t, T p)
        {
            return sqrt(p / (r_ * t));
        }
        T et(T t, T p , T V)
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
#pragma once

template <typename U>
class BoundaryCondition
{
public:
    virtual U eval(U & Fc, const U & U_extropolated) const  = 0;
    virtual U eval(U & Fc) const  = 0;
};

template <typename U>
class BoundaryConditionDirichlet : public BoundaryCondition<U>
{
    public:
    U value;
    BoundaryConditionDirichlet(const U& val) : value{val} {}
};

template< typename T >
class ConstantTemperature : public BoundaryConditionDirichlet<T>
{
    public:
    ConstantTemperature(T val) : BoundaryConditionDirichlet<T>{val} {}
    T eval(T & Fc) const  override
    {
        return value*Fc;
    }
    T eval(T & Fc,const T & T_extropolated) const  override
    {
        return value*Fc;
    }
};

template< typename T >
class ThermalWall : public BoundaryConditionDirichlet<T>
{
    public:
    ThermalWall() : BoundaryConditionDirichlet<T>{0.} {}
    T eval(T & Fc) const  override
    {
        Fc = T{};
        return value*Fc;
    }
    T eval(T & Fc,const T & T_extropolated) const  override
    {
        Fc = T{};
        return value*Fc;
    }
};

template< typename T >
class VelocityWall : public BoundaryConditionDirichlet<T>
{
    public:
    VelocityWall() : BoundaryConditionDirichlet<T>{T{}} {}
    T eval(T & Fc) const  override
    {
        // Fc = T{};
        // return value*Fc;
        return T{};
    }
    T eval(T & Fc,const T & T_extropolated) const  override
    {
        // Fc = T{};
        // return value*Fc;
        return T{};
    }
};

template <typename U>
class BoundaryConditionNeumann : public BoundaryCondition<U>
{
    
};
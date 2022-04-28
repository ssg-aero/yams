#pragma once
#include <datastorage.h>
#include <optional>
#include <functional>
#include <gbs/curves>
namespace yams
{
    template <typename T>
    struct GridInfo{
        std::shared_ptr< MeridionalGrid<T> > g;
        std::shared_ptr< Array2d<Grid2dMetricsPoint<T>> >   g_metrics;
        T d_ksi;
        T d_eth;
        size_t ni;
        size_t nj;
        T j_0 {0.5};
        const char* gas_name;
        bool rho_cst = true;
        T R = 287.04; // perfect gas constant
        T Pref = 1.01325e5;
        T Tref = 288.;
        T RF = 0.01;
        T tol_newtow_mf_f = 1e-5;
        T tol_newtow_mf_u = 1e-6;
        size_t vm_distribution_max_count =200;
    };

    enum class MeridionalBladeMode
    {
        DESIGN_BETA_OUT,
        DESIGN_ALPHA_OUT,
        DESIGN_PSI,
        DIRECT
    };

    template <typename T>
    struct BladeInfo{
        std::string name;
        int i1 = -1;
        int is = -1; 
        int i2 = -1;
        T omg     = 0.;
        size_t z_{2};
        std::function<T(T)> omg_;
        MeridionalBladeMode mode = MeridionalBladeMode::DIRECT;
        std::function<T(T)> beta_out;
        std::function<T(T)> alpha_out;
        std::function<T(T)> psi;
        std::function<T(T,T)> k;
        std::function<T(T,T)> tb;
        std::function<T(T,T)> eps;
    };

    enum class MeridionalBC
    {
        INLET_Mf_Ts_Ps_Vu,
        INLET_VmMoy_Ts_Ps_Vu,
        INLET_Vm_Ts_Ps_Vu,
        CON,
    };

    template <typename T>
    struct InletBC{
        MeridionalBC mode = MeridionalBC::INLET_VmMoy_Ts_Ps_Vu;
        std::function<T(T)> Ps = [](auto l_rel){return 1.01325e5;};
        std::function<T(T)> Ts = [](auto l_rel){return 300.;};
        std::function<T(T)> Vm = [](auto l_rel){return 30.;};
        std::function<T(T)> Vu = [](auto l_rel){return 0.;};
        std::function<T(T)> Pt = [](auto l_rel){return 1.01325e5;};
        std::function<T(T)> Tt = [](auto l_rel){return 300.;};
        T Mf =1.;
        T Vm_moy = 30.;
    };

    template <typename T>
    struct SolverLog
    {
        std::vector<T> delta_pos_max;
        std::vector<T> delta_pos_moy;
        std::vector<std::vector<T>> delta_pos;
        void clear(){delta_pos.clear(),delta_pos_max.clear();delta_pos_moy.clear();}
    };

    // enum class MassFlowDistrib
    template <typename T>
    struct 
    SolverCase
    {
        std::shared_ptr< GridInfo<T> > gi;
        std::vector<BladeInfo<T>> bld_info_lst;
        InletBC<T> inlet;
        std::vector<T> mf;
        SolverLog<T> log;
        size_t max_geom = 200;
        T eps = 0.001;
        T tol_rel_mf = 1e-4;
        T tol_rel_pos = 1e-5;
        bool relocate = true;
        size_t mf_ref_span = 0;
        bool mf_uniform = false;
        bool use_meridional_grad=false;
    };

    template <typename T>
    class MeridionalGridConnection
    {
        std::shared_ptr< SolverCase<T> > m_left;
        std::shared_ptr< SolverCase<T> > m_right;
        size_t j1_left;
        size_t j2_left;
        size_t j1_right;
        size_t j2_right;
        public:
        MeridionalGridConnection( const std::shared_ptr< SolverCase<T> > &left, const std::shared_ptr< SolverCase<T> > &right, const std::array<size_t, 4> &connection_spans) : 
            m_left{left}, 
            m_right{right}, 
            j1_left{connection_spans[0]}, 
            j2_left{connection_spans[1]},
            j1_right{connection_spans[2]}, 
            j2_right{connection_spans[3]}
        {
            if( (j1_left - j2_left) != (j1_right - j2_right) )
            {
                throw std::invalid_argument("Default connection requires same stream lines number.");
            }
        }
        MeridionalGridConnection( const std::shared_ptr< SolverCase<T> > &left, const std::shared_ptr< SolverCase<T> > &right) : 
            MeridionalGridConnection<T>{left, right, {0, left->nRows()-1, 0, right->nRows()-1}} {}

        auto isLeft(const std::shared_ptr< SolverCase<T> > &solver_case)
        {
            return m_left == solver_case;
        }
        auto isRight(const std::shared_ptr< SolverCase<T> > &solver_case)
        {
            return m_right == solver_case;
        }

        // auto copyLeftToRightGhost() -> void
        // {
        //     std::copy(
        //         std::next( m_left->begin(m_left->nRows()-1), j1_left),
        //         std::next( m_left->end(m_left->nRows()-1), j2_left),
        //         std::next( m_right->begin(0), j1_right)
        //     );
        //     std::copy(
        //         std::next( m_left->begin(m_left->nRows()-2), j1_left),
        //         std::next( m_left->end(m_left->nRows()-2), j2_left),
        //         std::next( m_right->begin(-1), j1_right)
        //     );
        //     std::copy(
        //         std::next( m_right->begin(1), j1_right),
        //         std::next( m_right->begin(1), j2_right),
        //         std::next( m_left->begin(m_left->nRows()), j1_left)
        //     );
        // }
        auto getIterators(size_t i_left, size_t i_right)
        {
            auto g_left  = m_left->gi->g ;
            auto g_right = m_right->gi->g;

            auto to_begin = std::next( g_right->begin(i_right), j1_right);
            auto to_end   = std::next( g_right->begin(i_right), j2_right+1);
            auto from_begin = std::next( g_left->begin(i_left), j1_left);
            auto from_end   = std::next( g_left->begin(i_left), j2_left+1);
            return std::make_tuple(to_begin, to_end, from_begin, from_end);
        }
        
        auto copyLeftToRight() -> void
        {
            auto [ to_begin, to_end, from_begin, from_end ] = getIterators(m_left->gi->ni - 1, 0);

            std::copy(
                from_begin,
                from_end,
                to_begin
            );
        }

        auto interpolateCurvature() -> void
        {
            auto [ to_begin, to_end, from_begin, from_end ] = getIterators(m_left->gi->ni - 2, 1);
            // std::vector<T> cur(size());
            // std::transform(
            //     from_begin, from_end,
            //     to_begin,
            //     cur.begin(),
            //     [](const auto &gp_left, const auto &gp_right)
            //     {
            //         return (gp_left.cur+gp_right.cur) / 2;
            //     }
            // );

            std::vector<MeridionalGridPoint<T>> cur(size());
            std::transform(
                from_begin, from_end,
                to_begin,
                cur.begin(),
                [](const auto &gp_left, const auto &gp_right)
                {
                    auto gp{gp_left };
                    gp.cur= (gp_left.cur+gp_right.cur) / 2;
                    gp.gam= (gp_left.gam+gp_right.gam) / 2;
                    gp.phi= (gp_left.phi+gp_right.phi) / 2;
                    return gp;
                }
            );

            auto [ to_right_begin, to_right_end, to_left_begin, to_left_end ] = getIterators(m_left->gi->ni - 1, 0);

            std::transform(
                cur.begin(),
                cur.end(),
                to_right_begin,
                to_right_begin,
                [](const auto &gp_, const auto &gp_orig )
                {
                    auto gp{gp_orig };
                    gp.cur = gp_.cur;
                    gp.gam = gp_.gam;
                    gp.phi = gp_.phi;
                    return gp;
                }
            );

            std::transform(
                cur.begin(),
                cur.end(),
                to_left_begin,
                to_left_begin,
                [](const auto &gp_, const auto &gp_orig )
                {
                    auto gp{gp_orig };
                    gp.cur = gp_.cur;
                    gp.gam = gp_.gam;
                    gp.phi = gp_.phi;
                    return gp;
                }
            );

            // std::transform(
            //     cur.begin(),
            //     cur.end(),
            //     to_right_begin,
            //     to_right_begin,
            //     [](auto cur, const auto &gp )
            //     {
            //         auto gp_{gp};
            //         gp_.cur = cur;
            //         return gp_;
            //     }
            // );

            // std::transform(
            //     cur.begin(),
            //     cur.end(),
            //     to_left_begin,
            //     to_left_begin,
            //     [](auto cur, const auto &gp )
            //     {
            //         auto gp_{gp};
            //         gp_.cur = cur;
            //         return gp_;
            //     }
            // );

        }

        auto setBoth(size_t i,T v, const auto &f)
        {
            // f()
        }

        auto size() const { return j2_left - j1_left + 1;}
    };

    template <typename T>
    class SolverCaseSet
    {
        std::list< std::shared_ptr< SolverCase<T> > > m_cases;
        std::list< MeridionalGridConnection<T> > m_connections;
        auto containsCase(const std::shared_ptr< SolverCase<T> > &solver_case) -> bool
        {
            return std::end(m_cases) != std::find(m_cases.begin(), m_cases.end(), solver_case);
        }
        public:
        // SolverCaseSet() = default;
        SolverCaseSet() = delete;
        SolverCaseSet(const std::shared_ptr< SolverCase<T> > & inlet) : m_cases{inlet} {}

        auto addCases(const std::shared_ptr< SolverCase<T> > &left, const std::shared_ptr< SolverCase<T> > &right, const std::array<size_t, 4> &connection_spans)
        {
            left->mf_ref_span = left->gi->ni - 1;
            right->mf_ref_span = 0;
            auto has_left = containsCase(left);
            auto has_right= containsCase(right);

            if( has_right && has_left)
            {
                throw std::invalid_argument("Channels already connected.");
            }
            MeridionalGridConnection<T> connection{left, right, connection_spans};
            m_connections.push_back( connection );

            if(!has_left)
            {
                m_cases.push_back(left);
            }

            if(!has_right)
            {
                m_cases.push_back(right);
            }
        }

        auto addCases(const std::shared_ptr< SolverCase<T> > &left, const std::shared_ptr< SolverCase<T> > &right )// -> bool
        {
            if( left->gi->nj != right->gi->nj)
            {
                throw std::invalid_argument("Default connection requires same stream lines number.");
            }
            addCases(left, right, {0, left->gi->nj-1, 0, right->gi->nj-1});
        }

        auto & cases() {return m_cases;}
        auto & connections() {return m_connections;}
        const auto & cases() const {return m_cases;}
        const auto & connections() const {return m_connections;}
        
    };


    template <typename T>
    class SolverCaseBiPass
    {
        std::shared_ptr<SolverCase<T> > m_prim;
        std::shared_ptr<SolverCase<T> > m_sec;
        SolverCase<T> m_inlet;
        size_t m_i_bip;
        void computeInletGrid(){
            auto ni = m_i_bip + 1;
            auto nj_prim =  m_prim->gi->nj;
            auto nj_sec  =  m_sec->gi->nj;
            auto nj = nj_prim + nj_prim - 1;
            gbs::points_vector<T,2> pts(ni*nj);
            for(size_t i{}; i <= m_i_bip; i++)
            {
                
            }
        }
        public:
        // InletBC<T> inlet;
        auto & primary(){return *m_prim;}
        auto & secondary(){return *m_sec;}
        auto & inlet(){return m_inlet;}
        auto iBip(){return m_i_bip;}
        T BPR = 2.;
        SolverCaseBiPass(
            const std::shared_ptr< SolverCase<T> > &prim, 
            const std::shared_ptr< SolverCase<T> > &sec, size_t iBip
        ) : m_prim{prim}, m_sec{sec}, m_i_bip{iBip} 
        {
            computeInletGrid();
        }

    };

}
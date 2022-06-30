#include "meridionalsolvercase.h"

namespace yams
{
    /**
     * @brief Set all blades in direct mode, copy relative fluid angle to metal angle (deviation model is not taken into account yet)
     * 
     * @tparam T 
     * @param sover_case 
     */
    template <typename T>
    void switch_to_direct(SolverCase<T> &solver_case)
    {
        auto &g = *(solver_case.gi->g);
        for(auto & bld_info: solver_case.bld_info_lst)
        {
            if( bld_info.mode != MeridionalBladeMode::DIRECT)
            {
                bld_info.mode = MeridionalBladeMode::DIRECT;
                auto i1 = bld_info.i1;
                auto i2 = bld_info.i2;
                for( auto i = i1; i <= i2; i++ )
                {
                    std::for_each(
                        g.begin(i), g.end(i),
                        [](auto &gp)
                        {
                            gp.k = gp.bet; // TODO: add deviation model here
                        }
                    );
                }
            }
        }
    }
}
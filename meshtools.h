#pragma once
#include <vector>
#include <memory>
#include <gbs/curves>
#include <gbs/bscinterp.h>
#include <gbs-mesh/tfi.h>
namespace yams
{
    template <typename T>
    using crv_vector = std::vector<std::shared_ptr<gbs::Curve<T, 2>>>;

    template <typename T>
    auto build_iso_ksi_curves(const crv_vector<T> &crv_iso_eth, const std::vector<T> &u, size_t max_deg = 3)
    {
        auto n_eth = crv_iso_eth.size();
        if (n_eth < 2)
        {
            throw std::length_error("yams::mesh_channel, 2 Curves at least needed.");
        }
        auto n_ksi = u.size();
        crv_vector<T> crv_iso_ksi(n_ksi);

        auto v1 = static_cast<T>(0.);
        auto v2 = static_cast<T>(1.);
        auto v = gbs::make_range(v1, v2, n_eth);

        std::transform(
            u.begin(), u.end(), crv_iso_ksi.begin(),
            [&crv_iso_eth, &v, max_deg, n_eth](T u_)
            {
                gbs::points_vector<T, 2> pts(n_eth);
                std::transform(
                    crv_iso_eth.begin(), crv_iso_eth.end(), pts.begin(),
                    [u_](const auto &crv)
                    {
                        return crv->value(u_);
                    });
                auto deg = std::min(max_deg, n_eth-1);
                return std::make_shared<gbs::BSCurve<T, 2>>(gbs::interpolate(pts, v, deg));
            });
        return std::make_pair(
            crv_iso_ksi,
            v
        ); 
    }

    template <typename T>
    auto mesh_channel(const crv_vector<T> &crv_iso_eth, const crv_vector<T> &crv_iso_ksi, const std::vector<T> &u, const std::vector<T> &v, size_t nu, size_t nv)
    {
        const size_t dim = 2;
        const size_t P = 2;
        const size_t Q = 1;

        return gbs::tfi_mesh_2d<T, 2, P, Q, true>(crv_iso_ksi, crv_iso_eth, u, v, nu, nv);
    }
    /**
     * @brief 
     * 
     * @tparam T 
     * @param crv_iso_eth 
     * @param u 
     * @param max_deg 
     * @return auto 
     */
    template <typename T>
    auto mesh_channel(const crv_vector<T> &crv_iso_eth, const std::vector<T> &u, size_t nu, size_t nv, size_t max_deg = 3)
    {

        auto [ crv_iso_ksi, v ] = build_iso_ksi_curves<T>(crv_iso_eth,u, max_deg);

        return mesh_channel(crv_iso_eth, crv_iso_ksi, u, v, nu, nv);
    }



}
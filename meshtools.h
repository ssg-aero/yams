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
    auto build_iso_ksi_curves(const crv_vector<T> &iso_eth, const std::vector<T> &ksi_i, size_t max_deg = 3)
    {
        auto n_eth = iso_eth.size();
        if (n_eth < 2)
        {
            throw std::length_error("yams::mesh_channel, 2 Curves at least needed.");
        }
        auto n_ksi = ksi_i.size();
        crv_vector<T> iso_ksi(n_ksi);

        auto eth_1 = static_cast<T>(0.);
        auto eth_2 = static_cast<T>(1.);
        auto eth_j = gbs::make_range(eth_1, eth_2, n_eth);

        std::transform(
            ksi_i.begin(), ksi_i.end(), iso_ksi.begin(),
            [&iso_eth, &eth_j, max_deg, n_eth](T u_)
            {
                gbs::points_vector<T, 2> pts(n_eth);
                std::transform(
                    iso_eth.begin(), iso_eth.end(), pts.begin(),
                    [u_](const auto &crv)
                    {
                        return crv->value(u_);
                    });
                auto deg = std::min(max_deg, n_eth-1);
                return std::make_shared<gbs::BSCurve<T, 2>>(gbs::interpolate(pts, eth_j, deg));
            });
        return std::make_pair(
            iso_ksi,
            eth_j
        ); 
    }

    template <typename T>
    auto mesh_channel(const crv_vector<T> &iso_ksi, const crv_vector<T> &iso_eth, const std::vector<T> &ksi_i, const std::vector<T> &eth_j, size_t n_iso_ksi, size_t n_iso_eth)
    {
        const size_t dim = 2;
        const size_t P = 1;
        const size_t Q = 2;

        return gbs::tfi_mesh_2d<T, 2, P, Q, true>(iso_ksi, iso_eth, ksi_i, eth_j, n_iso_ksi, n_iso_eth);
    }
    /**
     * @brief 
     * 
     * @tparam T 
     * @param iso_eth 
     * @param ksi_i 
     * @param n_ksi 
     * @param n_eth 
     * @param max_deg 
     * @return auto [pts, nj, ni, n_iso_eth, n_iso_ksi]
     */
    template <typename T>
    auto mesh_channel(const crv_vector<T> &iso_eth, const std::vector<T> &ksi_i, size_t n_iso_ksi, size_t n_iso_eth, size_t max_deg = 3)
    {

        auto [ iso_ksi, eth_j ] = build_iso_ksi_curves<T>(iso_eth,ksi_i, max_deg);

        return mesh_channel(iso_eth, iso_ksi, eth_j, ksi_i, n_iso_eth, n_iso_ksi);
    }

    template <typename T>
    auto mesh_channel(const crv_vector<T> &iso_eth, const std::vector<T> &ksi_i, const std::vector<T> &eth_j, size_t n_iso_ksi, size_t n_iso_eth, size_t max_deg = 3)
    {

        auto [ iso_ksi, eth_j_ ] = build_iso_ksi_curves<T>(iso_eth,ksi_i, max_deg);

        return mesh_channel(iso_eth, iso_ksi, eth_j, ksi_i, n_iso_eth, n_iso_ksi);
    }

}
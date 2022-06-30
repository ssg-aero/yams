#pragma once

#include <gbs/curves>
#include <gbs/bscinterp.h>
#include <gbs-mesh/tfi.h>
#include <gbs-mesh/smoothing.h>

#include <vector>
#include <memory>

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
        const size_t Q = 1;

        return gbs::tfi_mesh_2d<T, 2, P, Q, true>(iso_ksi, iso_eth, ksi_i, eth_j, n_iso_ksi, n_iso_eth);
    }
    template <typename T>
    auto mesh_channel(const crv_vector<T> &iso_ksi, const crv_vector<T> &iso_eth, const std::vector<std::vector<T>> &ksi_i, const std::vector<std::vector<T>> &eth_j, size_t n_iso_ksi, size_t n_iso_eth)
    {
        const size_t dim = 2;
        const size_t P = 1;
        const size_t Q = 1;

        return gbs::tfi_mesh_2d<T, dim, P, Q, true>(iso_ksi, iso_eth, ksi_i, eth_j, n_iso_ksi, n_iso_eth);
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
/**
 * @brief Compute channel mesh 
 * 
 * @tparam T 
 * @param master_streams 
 * @param hard_points    : blades points on master stream
 * @param n_streams      : stream lines number
 * @param n_spans        : indicative span number
 * @param tol            : tolerance
 * @param max_deg        : degree for interpolation between hard points in span direction.
 * @return auto 
 */
    template <typename T>
    auto mesh_channel(const crv_vector<T> &master_streams, const std::vector<gbs::points_vector<T,2>> &hard_points, size_t n_streams, size_t n_spans, T tol, size_t max_deg = 3)
    // auto mesh_channel(const crv_vector<T> &master_streams, const std::vector<std::vector<std::array<T,2>>> &hard_points, size_t n_streams, size_t n_spans, T tol, size_t max_deg = 3)
    {

        auto nc_iso_eth = master_streams.size();
        if(nc_iso_eth != hard_points.size())
        {
            throw std::invalid_argument("incompatible size.");
        }
        // Project hard points on curve
        std::vector<std::vector<T>> eth_j(nc_iso_eth);
        auto hard_points_it = hard_points.begin();
        auto eth_j_it       = eth_j.begin();
        for(const auto &eth_crv : master_streams)
        {
            for(const auto &pt : *hard_points_it)
            {
                auto [u,d] = gbs::extrema_curve_point(*eth_crv, pt, tol);
                if(d>tol)
                {
                    throw std::invalid_argument("Hard point too far from curve.");
                }
                eth_j_it->push_back(u);
            }
            std::advance(hard_points_it,1);
            std::advance(eth_j_it,1);
        }
        // Build transverse curve set
        auto nc_iso_ksi = hard_points.front().size();
        crv_vector<T> iso_ksi(nc_iso_ksi);
        std::vector<std::vector<T>> ksi_i(nc_iso_ksi);
        for(size_t i{} ; i < nc_iso_ksi; i++)
        {
            gbs::points_vector<T,2> pts(nc_iso_eth);
            for(size_t j{} ; j < nc_iso_eth; j++)
            {
                if(nc_iso_ksi != hard_points[j].size())
                {
                    throw std::invalid_argument("incompatible size.");
                }
                pts[j] = hard_points[j][i];
            }
            auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH);
            auto deg = std::min(max_deg,nc_iso_eth-1);
            iso_ksi[i] = std::make_shared<gbs::BSCurve<T,2>>( gbs::interpolate(pts, u ,deg) );
            ksi_i[i]   = u;
        }

        return mesh_channel(master_streams, iso_ksi, ksi_i, eth_j, n_spans, n_streams);
    }

    template <typename T>
    auto mesh_channel(const crv_vector<T> &iso_eth, const std::vector<T> &ksi_i, const std::vector<T> &eth_j, size_t n_iso_ksi, size_t n_iso_eth, size_t max_deg = 3)
    {

        auto [ iso_ksi, eth_j_ ] = build_iso_ksi_curves<T>(iso_eth,ksi_i, max_deg);

        return mesh_channel(iso_eth, iso_ksi, eth_j, ksi_i, n_iso_eth, n_iso_ksi);
    }
/**
 * @brief Smooth mesh blocs
 * 
 * @tparam T 
 * @param pts             : mesh points
 * @param n_streams       : number of stream
 * @param n_span_per_bloc : blocs span count
 * @param max_it          : max smoothing iteration
 * @param tol             : tolerance eon grid
 * @return auto 
 */
    template <typename T>
    auto smooth_mesh(gbs::points_vector<T,2> &pts, size_t n_streams, const std::vector<size_t> &n_span_per_bloc, size_t max_it=100, T tol =1e-5)
    {
        std::vector<size_t> block_indices(n_span_per_bloc.size()+1);
        block_indices[0] = 0;
        std:transform(
            n_span_per_bloc.begin(), n_span_per_bloc.end(),
            block_indices.begin(),
            std::next(block_indices.begin()),
            []( auto nvi_, auto i_prev ){return i_prev + nvi_-1;}
        );
        auto n_spans = pts.size() / n_streams;
        for( size_t i{1}; i < block_indices.size(); i++)
        {
            auto [it, err_max] = gbs::elliptic_structured_smoothing(pts,n_spans,0, n_streams-1,block_indices[i-1],block_indices[i],max_it, tol);
            printf("iterations : %i, error max: %.3e\n", int(it), err_max);
        }
        return pts;
    }


}
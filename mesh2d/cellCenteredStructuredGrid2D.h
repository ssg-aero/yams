#pragma once
#include "geom.h"
// #include <thrust/host_vector.h>
#include "boundaryCondition.h"


template <typename T>
struct QuadCellStructured
{
    point2d<T> Xc;
    T V;
};

enum class QuadCellStructuredFaces {W = 0, S = 1, E = 2, N = 3};
enum class StructuredGrid2DFaces {W = 0, S = 1, E = 2, N = 3};
template <typename T, typename U>
class CellCenteredStructuredGrid2DConnection;
template <typename T, typename U>
class CellCenteredStructuredGrid2DSet;

template <typename T, typename U>
class CellCenteredStructuredGrid2D
{
    friend class CellCenteredStructuredGrid2DConnection<T,U>;
    friend class CellCenteredStructuredGrid2DSet<T,U>;
    template <typename V>
    // using host_vector = thrust::host_vector<V>;
    using host_vector = std::vector<V>;
    using shared_connection = std::shared_ptr<CellCenteredStructuredGrid2DConnection<T,U>>;
    using shared_boco   = std::shared_ptr<BoundaryCondition<U>>;

    host_vector<point2d<T>> m_grid_points;
    host_vector<QuadCellStructured<T>> m_cells;
    host_vector<point2d<T>> m_normals_WE;
    host_vector<point2d<T>> m_normals_SN;
    host_vector<U> m_U;
    host_vector<U> m_Q;
    size_t m_ni;
    size_t m_nj;

    shared_boco bc_W{};
    shared_boco bc_S{};
    shared_boco bc_E{};
    shared_boco bc_N{};

    shared_connection con_W;
    shared_connection con_S;
    shared_connection con_E;
    shared_connection con_N;

    public:

    /**
     * @brief Construct a new Cell Centered Structured Grid 2D object and takes ownership of points
     * 
     * @tparam _V 
     * @param ni 
     * @param grid_points 
     */

    template<typename _V>
    CellCenteredStructuredGrid2D(size_t ni, _V&grid_points)
    {
        // m_grid_points = std::move(grid_points);
        m_grid_points = grid_points;
        size_t nj = m_grid_points.size() / ni;
        assert( ni*nj == m_grid_points.size());

        resize(ni,nj);

        auto id_grid_points = [ni](auto i, auto j){return i+ni*j;};

        // --- Build cells centroids and volumes
        for (size_t j{}; j < nj - 1; j++)
        {
            for (size_t i{}; i < ni - 1; i++)
            {
                const auto &p1 = m_grid_points[id_grid_points(i, j)];
                const auto &p2 = m_grid_points[id_grid_points(i + 1, j)];
                const auto &p3 = m_grid_points[id_grid_points(i + 1, j + 1)];
                const auto &p4 = m_grid_points[id_grid_points(i, j + 1)];
                m_cells[id(i, j)].Xc = T(0.25) * (p1 + p2 + p3 + p4);
                m_cells[id(i, j)].V = T(0.50) * ((p1.x - p3.x) * (p2.y - p4.y) + (p4.x - p2.x) * (p1.y - p3.y));
            }
        }
        // --- Build cells normals
        size_t index{};
        for (size_t j{}; j < nj; j++)
        {
            for(size_t i{}; i < ni-1 ; i++)
            {
                m_normals_SN[index] = rot90(m_grid_points[id_grid_points(i+1,j)]-m_grid_points[id_grid_points(i,j)]);
                index++;
            }
        }
        index = 0;
        for (size_t j{}; j < nj-1; j++)
        {
            for(size_t i{}; i < ni ; i++)
            {
                m_normals_WE[index] = rot90(m_grid_points[id_grid_points(i,j)]-m_grid_points[id_grid_points(i,j+1)]);
                index++;
            }
        }
    };

/**
 * @brief Allocate host grid storage cells, points must be allocated elsewhere
 * 
 * @param ni : number of points in i direction
 * @param nj : number of points in j direction
 */
    void resize(size_t ni, size_t nj)
    {
        auto nc = ( ni - 1 ) * ( nj - 1 );
        host_vector<QuadCellStructured<T>> cells(nc);
        host_vector<point2d<T>> normals_WE(ni*(nj-1));
        host_vector<point2d<T>> normals_SN(nj*(ni-1));
        host_vector<U> U_(nc,U{});
        host_vector<U> Q_(nc,U{});

        m_cells      = std::move(cells);
        m_normals_WE = std::move(normals_WE);
        m_normals_SN = std::move(normals_SN);
        m_U          = std::move(U_);
        m_Q          = std::move(Q_);
        m_ni         = ni -1;
        m_nj         = nj -1;
    }

    auto ni() const {return m_ni;}
    auto nj() const {return m_nj;}

    __host__ __device__
    size_t id(size_t i, size_t j) const
    {
        return i + m_ni * j;
    }

    __host__ __device__
    size_t id_WE(size_t i, size_t j) const
    {
        return i + (m_ni+1) * j;
    }

    __host__ __device__
    size_t id_SN(size_t i, size_t j) const
    {
        return i + m_ni * j;
    }

    size_t size(){return m_cells.size();}
    auto & cells() {return m_cells;}
    auto & values() {return m_U;}
    auto & sources() {return m_Q;}
    auto & normals_WE() { return m_normals_WE;}
    auto & normals_SN() { return m_normals_SN;}

    const auto & cells() const {return m_cells;}
    const auto & points() const {return m_grid_points;}
    const auto & values() const {return m_U;}
    const auto & sources() const {return m_Q;}
    const auto & normals_WE() const { return m_normals_WE;}
    const auto & normals_SN() const { return m_normals_SN;}

    auto & cell(size_t i, size_t j)
    {
        return m_cells[id(i,j)];
    }

    const auto & cell(size_t i, size_t j) const
    {
        return m_cells[id(i,j)];
    }

    auto & cell(size_t id)
    {
        return m_cells[id];
    }

    const auto & cell(size_t id) const
    {
        return m_cells[id];
    }

    auto & value(size_t id)
    {
        return m_U[id];
    }

    const auto & value(size_t id) const
    {
        return m_U[id];
    }

    const auto & point(size_t i, size_t j) const
    {
        return m_grid_points[i + (m_ni+1)*j];
    }

    auto cell_W(size_t i, size_t j)
    {
        if(i==0) return m_cells.end();
        else     return m_cells.begin()+id(i-1,j);
    }

    const auto& connection_W() const {return con_W;}
    const auto& connection_S() const {return con_S;}
    const auto& connection_E() const {return con_E;}
    const auto& connection_N() const {return con_N;}
    const auto& boco_W() const {return bc_W;}
    const auto& boco_S() const {return bc_S;}
    const auto& boco_E() const {return bc_E;}
    const auto& boco_N() const {return bc_N;}
    auto& boco(StructuredGrid2DFaces face)
    {
        switch (face)
        {
        case StructuredGrid2DFaces::W:
            return bc_W;
            break;
        case StructuredGrid2DFaces::S:
            return bc_S;
            break;
        case StructuredGrid2DFaces::E:
            return bc_E;
            break;
        default:// case StructuredGrid2DFaces::N:
            return bc_N;
            break;
        }
    }
    const auto& boco(StructuredGrid2DFaces face) const
    {
        return boc(face);
    }

};
#pragma once

#include <sparse>
#include "cellCenteredStructuredGrid2D.h"
#include "cellCenteredStructuredGrid2DConnection.h"
#include "boundaryCondition.h"
template <typename T, typename U>
class CellCenteredStructuredGrid2DSet
{
    template <typename V>
    // using host_vector = thrust::host_vector<V>;
    using host_vector = std::vector<V>;
    using grid2d      = CellCenteredStructuredGrid2D<T,U>;
    using shared_grid2d = std::shared_ptr<grid2d>;
    using shared_boco   = std::shared_ptr<BoundaryCondition<U>>;
    using shared_connection = std::shared_ptr<CellCenteredStructuredGrid2DConnection<T,U>>;

    T m_tol_connection{1e-6};
    host_vector<CellCenteredStructuredGrid2DConnection<T,U>> m_connections;
    host_vector<shared_boco> m_boco;
    host_vector<shared_grid2d> m_grids;
    host_vector<size_t> m_bloc_sizes{0}; //add 0 ad init


    public:

    auto & grids()
    {
        return m_grids;
    }

    const auto & grids() const
    {
        return m_grids;
    }

    auto & cell(size_t id)
    {
        auto it = std::lower_bound(m_bloc_sizes.begin(), m_bloc_sizes.end(), id+1) - 1;
        auto block_offset = *it;
        auto iB = std::distance(m_bloc_sizes.begin(),it);
        return m_grids[iB]->cell(id-block_offset);
    }

    auto & value(size_t id)
    {
        auto it = std::lower_bound(m_bloc_sizes.begin(), m_bloc_sizes.end(), id+1) - 1;
        auto block_offset = *it;
        auto iB = std::distance(m_bloc_sizes.begin(),it);
        return m_grids[iB]->value(id-block_offset);
    }

    auto & grid(size_t i) // TODO set const
    {
        return *m_grids[i];
    }

    const  auto & connections() const
    {
        return m_connections;
    }

    const  auto & connection(size_t id) const
    {
        return m_connections[id];
    }

    size_t blocId(const shared_grid2d &g)
    {
        auto it_g = std::find(m_grids.begin(), m_grids.end(), g);
        return std::distance(m_grids.begin(), it_g);
    }

    size_t blocOffset(size_t iB)
    {
        return m_bloc_sizes[iB];
    }

    const auto & blocOffsets() const
    {
        return m_bloc_sizes;
    }

    size_t idNeighbor(const shared_grid2d &g, const shared_connection &con_neighbor, size_t i_j)
    {
        auto iB_neighbor = idGrid(con_neighbor->g_neighbor(g));
        auto offset = m_bloc_sizes[iB_neighbor];
        return con_neighbor->getIdNeighbor(g, i_j) + offset;   
    }

    int id(size_t iB, size_t i, size_t j){

        const auto &g = grids()[iB];
        auto ni = g->ni();
        auto nj = g->nj();

        if( i == -1 || j == -1 || i == ni || j == nj)
        {
            if( i <= -1)
            {
                if(auto con_neighbor = g->connection_W())
                {
                    return idNeighbor(g, con_neighbor, j);
                }
                else
                {
                    return -idBoco(g->boco_W())-1;
                }
            }
            else if (j <= -1)
            {
                if(auto con_neighbor = g->connection_S())
                {
                    return idNeighbor(g, con_neighbor, i);
                }
                else
                {
                    return -idBoco(g->boco_S())-1;
                }
            }
            else if (i >= ni)
            {
                if(auto con_neighbor = g->connection_E())
                {
                    return idNeighbor(g, con_neighbor, j);
                }
                else
                {
                    return -idBoco(g->boco_E())-1;
                }
            }
            else // j == nj
            {
                if(auto con_neighbor = g->connection_N())
                {
                    return idNeighbor(g, con_neighbor, i);
                }
                else
                {
                    return -idBoco(g->boco_N())-1;
                }
            }
        }
        else
        {
            return m_bloc_sizes[iB] + grid(iB).id(i,j);
        }
    }

    size_t idGrid(const shared_grid2d &g)
    {
        return std::distance(m_grids.begin(),std::find(m_grids.begin(),m_grids.end(),g));
    }

    size_t idBoco(const shared_boco &bc)
    {
        return std::distance(m_boco.begin(),std::find(m_boco.begin(),m_boco.end(),bc));
    }

    void addBoco(const shared_boco &boco)
    {
        m_boco.push_back(boco);
    }

    void setBoco(size_t iB, StructuredGrid2DFaces face, size_t iBoco)
    {
        m_grids[iB]->boco(face) = m_boco[iBoco];
    }

    bool addToGrid(const CellCenteredStructuredGrid2D<T,U> &g, bool connect=true)
    {

        m_grids.push_back(std::make_shared<CellCenteredStructuredGrid2D<T,U>>(g));
        m_bloc_sizes.push_back(m_bloc_sizes.back() + m_grids.back()->cells().size());

        bool connected{false};
        if(connect)
        {
            auto it = m_grids.begin();
            auto it_end = std::next(it, m_grids.size() - 1);
            while (it != it_end && !connected)
            {
                connected = connectGrids(*it, *it_end);
                it++;
            }
        }
        return connected;
    }

    template <typename I>
    void buildSolverMatrix(SparseMatrixSolverBase<T,I> &spm_solver)
    {
        size_t iB{};
        size_t bloc_offset{};
        for(auto &msh: m_grids)
        {
            bloc_offset = m_bloc_sizes[iB];
            for (size_t j{}; j < msh->nj(); j++)
            {
                for(size_t i{}; i < msh->ni() ; i++)
                {
                    std::vector<int> C_;
                    std::vector<int> V_;

                    auto id_S = id(iB,i,j-1);
                    auto id_W = id(iB,i-1,j);
                    auto id_C = id(iB,i,j);
                    auto id_E = id(iB,i+1,j);
                    auto id_N = id(iB,i,j+1);

                    if(id_S>=0) add_to_col(V_, C_, T(1.), id_S);
                    if(id_W>=0) add_to_col(V_, C_, T(1.), id_W);
                    add_to_col(V_, C_, T(5.), id_C);// For diag dominant
                    if(id_E>=0) add_to_col(V_, C_, T(1.), id_E);
                    if(id_N>=0) add_to_col(V_, C_, T(1.), id_N);

                    spm_solver.add_Row(C_.begin(), C_.end(), V_.begin());
                }
            }
            iB++;
        }
    }

    auto connectGrids( shared_grid2d& g1,  shared_grid2d&g2) -> bool
    {
        auto ni1m = g1->ni();
        auto ni2m = g2->ni();
        auto nj1m = g1->nj();
        auto nj2m = g2->nj();
        if( g1->nj() == g2->nj() )
        {
            if( 
                sq_norm( g1->point(0,0)   , g2->point(ni2m,0)   ) < m_tol_connection &&
                sq_norm( g1->point(0,nj1m), g2->point(ni2m,nj2m)) < m_tol_connection
            )
            {
                m_connections.push_back(CellCenteredStructuredGrid2DConnection<T,U>{g1,g2,StructuredGrid2DConnectionType::WE});
                return true;
            }
            if( 
                sq_norm( g1->point(ni1m,0)   , g2->point(0,0)   ) < m_tol_connection &&
                sq_norm( g1->point(ni1m,nj1m), g2->point(0,nj2m)) < m_tol_connection
            )
            {
                m_connections.push_back(CellCenteredStructuredGrid2DConnection<T,U>{g1,g2,StructuredGrid2DConnectionType::EW});
                return true;
            }
        }
        if( g1->ni() == g2->ni() )
        {
            if(
                sq_norm( g1->point(0,0)   , g2->point(0,nj2m)    ) < m_tol_connection &&
                sq_norm( g1->point(ni1m,0), g2->point(ni2m,nj2m) ) < m_tol_connection
            )
            {
                m_connections.push_back(CellCenteredStructuredGrid2DConnection<T,U>{g1,g2,StructuredGrid2DConnectionType::SN});
                return true;
            }
            if(
                sq_norm( g1->point(0,nj1m)   , g2->point(0,0)    ) < m_tol_connection &&
                sq_norm( g1->point(ni1m,nj1m), g2->point(ni2m,0) ) < m_tol_connection
            )
            {
                m_connections.push_back(CellCenteredStructuredGrid2DConnection<T,U>{g1,g2,StructuredGrid2DConnectionType::NS});
                return true;
            }
        }
        return false;
    }
};
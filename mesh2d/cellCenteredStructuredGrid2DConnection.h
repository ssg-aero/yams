#pragma once
#include "cellCenteredStructuredGrid2D.h"

enum class StructuredGrid2DConnectionType {EW, WE, NS, SN};

template <typename T, typename U>
class CellCenteredStructuredGrid2DConnection
{
    using shared_grid2d = std::shared_ptr<CellCenteredStructuredGrid2D<T,U>>;
    // friend class CellCenteredStructuredGrid2DSet<T,U>;
    StructuredGrid2DConnectionType m_type;

    const shared_grid2d m_g1;
    const shared_grid2d m_g2;

public:
    CellCenteredStructuredGrid2DConnection(const shared_grid2d &g1, const shared_grid2d &g2, StructuredGrid2DConnectionType type)
    : m_g1{g1}, m_g2{g2}, m_type{type}
    {
        switch (m_type)
        {
        case StructuredGrid2DConnectionType::EW:
            {
                m_g1->con_E = std::make_shared<CellCenteredStructuredGrid2DConnection<T,U>>(*this);
                m_g2->con_W = std::make_shared<CellCenteredStructuredGrid2DConnection<T,U>>(*this);
            }
            break;
        case StructuredGrid2DConnectionType::WE:
            {
                m_g1->con_W = std::make_shared<CellCenteredStructuredGrid2DConnection<T,U>>(*this);
                m_g2->con_E = std::make_shared<CellCenteredStructuredGrid2DConnection<T,U>>(*this);
            }
            break;
        case StructuredGrid2DConnectionType::NS:
            {
                m_g1->con_N = std::make_shared<CellCenteredStructuredGrid2DConnection<T,U>>(*this);
                m_g2->con_S = std::make_shared<CellCenteredStructuredGrid2DConnection<T,U>>(*this);
            }
            break;
        case StructuredGrid2DConnectionType::SN:
            {
                m_g1->con_S = std::make_shared<CellCenteredStructuredGrid2DConnection<T,U>>(*this);
                m_g2->con_N = std::make_shared<CellCenteredStructuredGrid2DConnection<T,U>>(*this);
            }
            break;
        }
    }

    auto getId2(size_t i_j1, size_t offset = 0) const
    {
        switch (m_type)
        {
        case StructuredGrid2DConnectionType::EW:
            return m_g2->id(offset,i_j1);
            break;
        case StructuredGrid2DConnectionType::WE:
            return m_g2->id(m_g2->ni()-1-offset,i_j1);
            break;
        case StructuredGrid2DConnectionType::NS:
            return m_g2->id(i_j1,offset);
            break;
        // case StructuredGrid2DConnectionType::SN:
        default:
            return m_g2->id(i_j1,m_g2->nj()-1-offset);
            break;
        }
    }

    auto getId1(size_t i_j2, size_t offset = 0) const
    {
        switch (m_type)
        {
        case StructuredGrid2DConnectionType::EW:
            return m_g1->id(m_g1->ni()-1-offset,i_j2);
            break;
        case StructuredGrid2DConnectionType::WE:
            return m_g1->id(offset,i_j2);
            break;
        case StructuredGrid2DConnectionType::NS:
            return m_g1->id(i_j2,m_g1->nj()-1-offset);
            break;
        // case StructuredGrid2DConnectionType::SN:
        default:
            return m_g2->id(i_j2,offset);
            break;
        }
    }

    auto getIdNeighbor(const shared_grid2d &g, size_t i_j, size_t offset = 0) const
    {
        return m_g1 == g ? getId2(i_j, offset) : getId1(i_j, offset);
    }

    const auto &g1()
    {
        return m_g1;  
    }

    const auto &g2()
    {
        return m_g2;
    }

    const auto &g_neighbor(const shared_grid2d &g)
    {
        return m_g1 == g ? m_g2 : m_g1;
    }
    const auto &cell_neighbor(const shared_grid2d &g, size_t i_j)
    {
        return  g_neighbor(g)->cell(getIdNeighbor(g,i_j));
    }
};

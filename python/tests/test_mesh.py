import pygbs.gbs as gbs
from pygbs import vtkplot as vbs
import pyams.yams as yams
from pytest import approx
import pytest
from math import cos, sin, radians

plot_on = True

@pytest.fixture
def channel1():
    """Fixture returning a simple channel"""
    knots = [0.,0.33,0.66,1.0]
    mult  = [3,1,1,3]

    crv1 = gbs.BSCurve2d(
        [
            [0.0,0.0],
            [0.2,0.0],
            [0.5,0.25],
            [0.8,0.3],
            [1.0,0.3],
        ],
        knots,
        mult,
        2
    )
    crv2 = gbs.BSCurve2d(
        [
            [0.1,0.25],
            [0.2,0.25],
            [0.50,0.55],
            [0.8,0.65],
            [0.9,0.65],
        ],
        knots,
        mult,
        2
    )
    crv3 = gbs.BSCurve2d(
        [
            [0.0,0.5],
            [0.2,0.5],
            [0.5,0.85],
            [0.8,1.0],
            [1.0,1.0],
        ],
        knots,
        mult,
        2
    )

    return [ crv1, crv2, crv3]

@pytest.fixture
def channel2():
    """Fixture returning a simple channel"""
    knots = [0.,0.33,0.66,1.0]
    mult  = [3,1,1,3]

    crv1 = gbs.BSCurve2d(
        [
            [0.0,0.0],
            [0.2,0.0],
            [0.5,0.25],
            [0.8,0.3],
            [1.0,0.3],
        ],
        knots,
        mult,
        2
    )

    crv2 = gbs.BSCurve2d(
        [
            [0.0,0.5],
            [0.2,0.5],
            [0.5,0.85],
            [0.8,1.0],
            [1.0,1.0],
        ],
        knots,
        mult,
        2
    )

    return [ crv1, crv2]  

@pytest.fixture
def channel3():
    """Fixture returning a simple channel"""
    hub_start = (0.0,[0.0,0.1])
    hub_end   = (1.0,[0.3,0.2])
    hub_le    = (0.3,[0.1,0.15],0)
    hub_te    = (0.6,[0.2,0.2],0)
    hub_te_tg = (0.2,[0.01,0.0],1)
    hub_end_tg= (0.3,[0.001,0.0],1)

    shr_start = (0.0,[0.0,0.3])
    shr_end   = (1.0,[0.3,0.3])
    shr_le    = (0.3,[0.12,0.3],0)
    shr_te    = (0.6,[0.17,0.3],0)
    shr_te_tg = (0.2,[0.01,0.0],1)
    shr_end_tg= (0.3,[0.001,0.0],1)

    return [ 
        gbs.interpolate(hub_start, hub_end, 
            [
                hub_le, 
                hub_te, 
                # hub_te_tg,
                # hub_end_tg,
            ], 
        2), 
        gbs.interpolate(shr_start, shr_end, 
            [
                shr_le, 
                shr_te, 
                # shr_te_tg,
                # shr_end_tg,
            ], 
        2)    
    ]  

@pytest.fixture
def channel4():
    hub_constains = [
        [
            [0.,0.1],
            [cos(radians(30)),sin(radians(30))],
        ],
        [
            [0.1,0.15],
            [cos(radians(30)),sin(radians(30))],
        ],
        [
            [0.2,0.2],
            [cos(radians(0)),sin(radians(0))],
        ],
        [
            [0.3,0.2],
            [cos(radians(0)),sin(radians(0))],
        ],
    ]
    shr_constains = [
        [
            [0.,0.3],
            [cos(radians(0)),sin(radians(0))],
        ],
        [
            [0.11,0.3],
            [cos(radians(0)),sin(radians(0))],
        ],
        [
            [0.18,0.28],
            [cos(radians(0)),sin(radians(0))],
        ],
        [
            [0.3,0.28],
            [cos(radians(0)),sin(radians(0))],
        ],
    ]

    return [
        gbs.interpolate_c1(hub_constains, [0.,0.1,0.2,0.3]),
        gbs.interpolate_c1(shr_constains, [0.,0.1,0.2,0.3]),
    ]


@pytest.mark.parametrize("state", [
    {
        "nu": 30,
        "nv": 15,
    },
])

def test_mesh_base_channel(channel1, state):

    crv_lst = channel1
    crv1, crv2, crv3 = crv_lst
    knots = crv_lst[1].knots()

    nu = state['nu']
    nv = state['nv']

    pts, ni, nj, n_iso_ksi, n_iso_eth = yams.mesh_channel(crv_lst, knots, nv, nu)

    assert gbs.dist(pts[0],crv1.begin()) == approx(0.)
    assert gbs.dist(pts[ni-1],crv1.end()) == approx(0.)
    assert gbs.dist(pts[ni * (nj - 1)],crv3.begin()) == approx(0.)
    assert gbs.dist(pts[ni-1 + ni * (nj - 1)],crv3.end()) == approx(0.)

    if plot_on:
        sgrid = gbs.make_structuredgrid(pts, ni, nj)
        sgridActor = vbs.make_sgrid_actor(sgrid, color='Black' ,grid_only=True)
        vbs.render_actors([
            gbs.make_actor(crv1), 
            gbs.make_actor(crv2), 
            gbs.make_actor(crv3), 
            gbs.make_actor(pts,render_as_sphere=True), 
            sgridActor,
        ])

@pytest.mark.parametrize("state", [
    {
        "nu": 30,
        "nv": 15,
    },
])

def test_grid_base_channel(channel1, state):
    crv_lst = channel1
    knots = crv_lst[1].knots()

    nu = state['nu']
    nv = state['nv']

    pts, ni, nj, n_iso_ksi, n_iso_eth = yams.mesh_channel(crv_lst, knots, nv, nu)
    sgrid = gbs.make_structuredgrid(pts, ni, nj)
    g = yams.read_vtk_grid(sgrid)

    assert g(0,0).x == approx( crv_lst[0].begin()[0] )
    assert g(0,0).y == approx( crv_lst[0].begin()[1] )
    assert g(0,0).Vm == approx( 0.0 )
    assert g(0,0).Vu == approx( 0.0 )

    gp = yams.MeridionalGridPoint()
    gp.Vm = 120.0
    gp.Vu =  30.0
    gp.Pt = 2e5
    gp.Tt = 600.0

    g.init(gp)
    assert g(ni-1,nj-1).Pt == approx(gp.Pt)
    assert g(ni-1,nj-1).Tt == approx(gp.Tt)
    assert g(ni-1,nj-1).Vm == approx(gp.Vm)
    assert g(ni-1,nj-1).Vu == approx(gp.Vu)


    try:
        gi = yams.GridInfo()
        gi.g = g
    except:
        assert False

@pytest.mark.parametrize("state", [
    {
        "nu": 30,
        "nv": 15,
    },
])

def test_grid_metrics_base_channel(channel1, state):
    crv_lst = channel1
    knots = crv_lst[1].knots()

    nu = state['nu']
    nv = state['nv']

    pts, ni, nj, n_iso_ksi, n_iso_eth = yams.mesh_channel(crv_lst, knots, nv, nu)
    sgrid = gbs.make_structuredgrid(pts, ni, nj)
    gi =yams. make_grid_info(sgrid)

    assert gi.g(0,0).cur    != 0.
    if plot_on: yams.plot(gi.g,"cur",False)

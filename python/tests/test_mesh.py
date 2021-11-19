import pygbs.gbs as gbs
from pygbs import vtkplot as vbs
import pyams.yams as yams
from pytest import approx


plot_on = True

def test_mesh_base_channel():
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

    crv_lst = [ crv1, crv2, crv3]

    nu = 30
    nv = 15

    pts, ni, nj, n_iso_ksi, n_iso_eth = yams.mesh_channel(crv_lst, knots, nv, nu)

    assert gbs.dist(pts[0],crv1.begin()) == approx(0.)
    assert gbs.dist(pts[ni-1],crv3.begin()) == approx(0.)
    assert gbs.dist(pts[ni * (nj - 1)],crv1.end()) == approx(0.)
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
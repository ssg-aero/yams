import pygbs.gbs as gbs
from pygbs import vtkplot as vbs
import pyams.yams as yams
from pytest import approx


plot_on = True

def test_mesh_base_channel():
    crv1 = gbs.BSCurve2d(
        [
            [0.0,0.0],
            [0.2,0.0],
            [0.8,0.3],
            [1.0,0.3],
        ],
        [0.,0.5,1.],
        [3,1,3],
        2
    )
    crv2 = gbs.BSCurve2d(
        [
            [0.0,0.25],
            [0.35,0.25],
            [0.8,0.65],
            [1.0,0.65],
        ],
        [0.,0.5,1.],
        [3,1,3],
        2
    )
    crv3 = gbs.BSCurve2d(
        [
            [0.0,0.5],
            [0.4,0.5],
            [0.8,1.0],
            [1.0,1.0],
        ],
        [0.,0.5,1.],
        [3,1,3],
        2
    )

    knots   = crv1.knots()
    crv_lst = [ crv1, crv2, crv3]

    nu = 30
    nv = 15

    pts, ni, nj, n_iso_eth, n_iso_ksi = yams.mesh_channel(crv_lst, knots, nu, nv)

    assert gbs.dist(pts[0],crv1.begin()) == approx(0.)
    assert gbs.dist(pts[ni-1],crv3.begin()) == approx(0.)
    assert gbs.dist(pts[ni * (nj - 1)],crv1.end()) == approx(0.)
    assert gbs.dist(pts[ni-1 + ni * (nj - 1)],crv3.end()) == approx(0.)

    if plot_on:
        sgrid = gbs.make_structuredgrid(pts, ni, nj)
        sgridActor = vbs.make_sgrid_actor(sgrid, color='Black' ,grid_only=True)
        vbs.render_actors([
            gbs.make_actor(crv1), 
            gbs.make_actor(crv3), 
            gbs.make_actor(pts,render_as_sphere=True), 
            sgridActor,
        ])
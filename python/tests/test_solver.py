from test_mesh import channel1, channel2
import pygbs.gbs as gbs
import pyams.yams as yams
from pytest import approx
import pytest
import vtk
from math import pi, sin
plot_on = True

@pytest.mark.parametrize("state", [
    {
        "nu": 30,
        "nv": 21,
    },
])

def test_grid_metrics_base_channel(channel2, state):
    crv_lst = channel2
    knots = crv_lst[1].knots()

    nu = state['nu']
    nv = state['nv']

    pts, ni, nj, n_iso_ksi, n_iso_eth = yams.mesh_channel(crv_lst, knots, nv, nu)
    sgrid = gbs.make_structuredgrid(pts, ni, nj)

    solver_case = yams.SolverCase()
    solver_case.gi = yams. make_grid_info(sgrid)
    solver_case.max_geom = 500
    solver_case.gi.RF = 0.05
    solver_case.inlet.Ts = lambda l_rel : 300. + 30 * sin(l_rel * pi)

    yams.curvature_solver(solver_case)

    if plot_on: 
        yams.plot(solver_case.gi.g,"Vm",False)
        yams.plot(solver_case.gi.g,"Ts",False)
        yams.plot(solver_case.gi.g,"cur",True)
        yams.plot_residual(solver_case.log)


def test_grid_metrics_mixed_channel( ):
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001.vts")
    reader.Update() 
    sgrid = reader.GetOutput()
    solver_case = yams.SolverCase()
    solver_case.gi = yams.make_grid_info(sgrid)

    yams.curvature_solver(solver_case)

    if plot_on: 
        yams.plot(solver_case.gi.g,"Vm",False)
        yams.plot_residual(solver_case.log)
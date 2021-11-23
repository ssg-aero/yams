from test_mesh import channel1
import pygbs.gbs as gbs
import pyams.yams as yams
from pytest import approx
import pytest
import vtk
plot_on = True

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
    solver_case = yams.SolverCase()
    solver_case.gi = yams. make_grid_info(sgrid)
    solver_case.inlet.Mf = 1.
    solver_case.inlet.mode

    yams.curvature_solver(solver_case)

    if plot_on: 
        yams.plot(solver_case.gi.g,"Vm",False)
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
from test_mesh import channel1, channel2
import pygbs.gbs as gbs
import pyams.yams as yams
from pytest import approx
import pytest
import vtk
from math import pi, radians, sin
plot_on = True

def config_solver(state, solver_case):
    solver_case.max_geom = state['max_geom']
    solver_case.gi.RF = state['relax_factor']
    solver_case.inlet.Ts = lambda l_rel : state['Ts_min'] + state['Ts_delta'] * sin(l_rel * pi)
    solver_case.inlet.Vm_moy = state['Vm_moy']


@pytest.mark.parametrize("state, expected", [
    ({
        "nu": 30,
        "nv": 21,
        "Vm_moy":30.,
        "Ts_min":300.,
        "Ts_delta":30,
        "max_geom":500,
        "relax_factor":0.05,
    },
    {
        "ral":1.4
    })
])


def test_solve_base_channel(channel2, state, expected):
    crv_lst = channel2
    knots = crv_lst[1].knots()

    nu = state['nu']
    nv = state['nv']

    pts, ni, nj, n_iso_ksi, n_iso_eth = yams.mesh_channel(crv_lst, knots, nv, nu)
    sgrid = gbs.make_structuredgrid(pts, ni, nj)

    solver_case = yams.SolverCase()
    solver_case.gi = yams. make_grid_info(sgrid)
    config_solver(state, solver_case)

    yams.curvature_solver(solver_case)

    assert solver_case.gi.g(0,0).Tt == approx( solver_case.gi.g(ni-1,0).Tt )
    assert solver_case.gi.g(0,int(nj/2)).Vm > solver_case.gi.g(ni-1,int(nj/2)).Vm *expected["ral"]
    assert solver_case.gi.g(0,int(nj/2)).Tt == approx( solver_case.gi.g(ni-1,int(nj/2)).Tt )
    if plot_on: 
        yams.plot(solver_case.gi.g,"Vm",False)
        yams.plot(solver_case.gi.g,"Ts",False)
        yams.plot(solver_case.gi.g,"cur",True)
        yams.plot_residual(solver_case.log)


def test_solve_mixed_channel( ):
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

@pytest.mark.parametrize("state, expected", [
    ({
        "nu": 15,
        "nv": 11,
        "Vm_moy":30.,
        "Ts_min":300.,
        "Ts_delta":15,
        "max_geom":500,
        "relax_factor":0.05,
        "k1":radians(0.),
        "k2":radians(30.),
    },
    {
        "beta1":radians(0.),
        "beta2":radians(30.),
    })
])

def test_solve_base_stator(channel2,state, expected):
    crv_lst = channel2
    knots = crv_lst[1].knots()

    nu = state['nu']
    nv = state['nv']

    pts, ni, nj, n_iso_ksi, n_iso_eth = yams.mesh_channel(crv_lst, knots, nv, nu)
    sgrid = gbs.make_structuredgrid(pts, ni, nj)

    blade_info = yams.BladeInfo()
    blade_info.i1 = int(ni/3.0)
    blade_info.i_s = int(ni/2.0)
    blade_info.i2 = int(2*ni/3.0)
    blade_info.mode = yams.MeridionalBladeMode.DIRECT

    solver_case = yams.make_solver_case(sgrid, [blade_info], lambda m, l : (1-m) * state['k1'] + m *state['k2'] )
    config_solver(state, solver_case)

    yams.curvature_solver(solver_case)

    for j in range(nj):
        assert solver_case.gi.g( blade_info.i2,j).bet == approx(expected['beta2'])
        assert solver_case.gi.g( blade_info.i1,j).bet == approx(expected['beta1'])

    if plot_on: 
        yams.plot(solver_case.gi.g,"Vu",True)
        yams.plot(solver_case.gi.g,"Vm",True)
        yams.plot(solver_case.gi.g,"Ts",True)
        yams.plot_residual(solver_case.log)


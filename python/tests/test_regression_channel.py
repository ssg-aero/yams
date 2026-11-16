"""Non-regression of the meridional curvature solver on an annular channel without blades.

The case and the sequence of calls reproduce what cosapp-turbo does
(ChannelSetAeroMeridianPreprocessor + ChannelAeroSolverMeridian.compute): inlet
VmMoy/Ts/Ps/Vu, compressible, stream lines relocation, uniform mass flow.

The reference values are the results of the solver as of 2026-11-16 (branch
refactor/core-for-tbf, before the solver bug fixes). A fix that changes the
solution must update them on purpose and say why in its commit.
"""
import math

import numpy as np
import pygbs.gbs as gbs
import pyams.yams as yams
import pytest
from pytest import approx

# Relative tolerance on the reference values: tight, to detect any change of the solution.
REL = 1e-5

# The solver stops at max_geom without reaching tol_rel_pos (known issue, see the review).
REF_GEOM_ITERATIONS = 200
REF_INLET_VM = [29.592071702021883, 27.098181787303805, 28.781955057304785, 29.115141168773057,
                28.919474240941945, 28.481496385001726, 27.92142017790463]
REF_OUTLET_VM = [39.132031220899655, 38.32818496755033, 37.908954190564074, 37.69977685020357,
                 37.62643170780902, 37.65081686945557, 37.75047065115762]
REF_OUTLET_PS = [99077.34885029127, 99114.74646516851, 99133.94643420634, 99143.44837394622,
                 99146.76780160255, 99145.66489538027, 99141.15036095279]
REF_MASS_FLOW_IN = 428.815882047719


@pytest.fixture
def channel_case():
    hub = gbs.BSCurve2d([[0., 0.], [0.9, 1.], [3., 1.]], [0., 1.], [3, 3], 2)
    shroud = gbs.BSCurve2d(
        [[0., 2.], [0.45, 2.], [0.9, 1.9], [1.5, 1.9], [2.1, 1.9], [3., 2.]],
        [0., 0.2, 0.4, 0.6, 1.], [3, 1, 1, 1, 3], 2,
    )
    # 31 stations from inlet to outlet, 7 stream lines from hub to shroud
    pts, ni, nj, _, _ = yams.mesh_channel([hub, shroud], [0., 1.], 7, 31)

    solver_case = yams.SolverCase()
    solver_case.gi = yams.make_grid_info(gbs.make_structuredgrid(pts, ni, nj))

    solver_case.max_geom = 200
    solver_case.tol_rel_pos = 1e-5
    solver_case.relocate = True
    solver_case.mf_uniform = True
    solver_case.gi.rho_cst = False
    solver_case.gi.RF = 0.05
    solver_case.gi.j_0 = 0.
    solver_case.gi.R = 287.04

    solver_case.inlet.mode = yams.MeridionalBC.INLET_VmMoy_Ts_Ps_Vu
    solver_case.inlet.Vm_moy = 30.
    solver_case.inlet.Vu = lambda l_rel: 0.
    solver_case.inlet.Ps = lambda l_rel: 1.e5
    solver_case.inlet.Ts = lambda l_rel: 288.15

    g = solver_case.gi.g
    for i, j in np.ndindex(g.nRows(), g.nCols()):
        g(i, j).Cp = 1004.5
    return solver_case


def mass_flow(g, i):
    """Mass flow through station i, trapezoidal rule along the span."""
    nj = g.nCols()
    return 2 * math.pi * sum(
        0.5 * (g(i, j).rho * g(i, j).Vm * g(i, j).y + g(i, j - 1).rho * g(i, j - 1).Vm * g(i, j - 1).y)
        * (g(i, j).l - g(i, j - 1).l)
        for j in range(1, nj)
    )


def test_channel_without_blades(channel_case):
    solver_case = channel_case
    g = solver_case.gi.g

    yams.curvature_solver(solver_case)

    ni, nj = g.nRows(), g.nCols()
    assert (ni, nj) == (31, 7)
    values = np.array([[g(i, j).Vm, g(i, j).Ps, g(i, j).Ts, g(i, j).rho] for i in range(ni) for j in range(nj)])
    assert np.all(np.isfinite(values))

    # Physics: no work and no loss without blades, mass conserved along the channel
    assert all(g(i, j).Tt == approx(288.15) for i in range(ni) for j in range(nj))
    mf = [mass_flow(g, i) for i in range(ni)]
    assert max(mf) - min(mf) < 0.03 * mf[0]

    # Non-regression (see module docstring)
    assert len(solver_case.log.delta_pos_max) == REF_GEOM_ITERATIONS
    assert [g(0, j).Vm for j in range(nj)] == approx(REF_INLET_VM, rel=REL)
    assert [g(ni - 1, j).Vm for j in range(nj)] == approx(REF_OUTLET_VM, rel=REL)
    assert [g(ni - 1, j).Ps for j in range(nj)] == approx(REF_OUTLET_PS, rel=REL)
    assert mf[0] == approx(REF_MASS_FLOW_IN, rel=REL)

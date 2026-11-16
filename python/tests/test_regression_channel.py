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


# Rotor row in DIRECT mode (metal angle imposed) with a non-uniform inlet total
# enthalpy, so that the span gradient of rothalpy dI_dl matters. Reference values
# computed after the dI_dl fix (it was only computed at the hub before): the
# solution was oscillating, with negative work at the hub, and did not converge.
REF_ROTOR_GEOM_ITERATIONS = 188
REF_ROTOR_OUTLET_VM = [30.957613653401285, 33.064532764177976, 34.462343434363035, 35.7156156089264,
                       36.84608407629142, 37.87890105620643, 38.83257679600115, 39.72084960376482,
                       40.55416125853153, 41.34070082387587, 42.087332724218015]
REF_ROTOR_OUTLET_VU = [7.637182116993376, 11.539949536233609, 14.659411791537616, 17.28025792141593,
                       19.564018082068486, 21.6029568262659, 23.45519096351598, 25.16024168923725,
                       26.746731576453413, 28.23660112281076, 29.648052634776512]
REF_ROTOR_WORK = [143.95749711757526, 266.570056863653, 388.9752652159659, 509.4677502489649,
                  627.8917660520528, 744.3535636184388, 859.0268806785461, 972.1070721958531,
                  1083.8019635217497, 1194.3354902406572, 1303.9894387023523]


@pytest.fixture
def rotor_direct_case():
    hub = gbs.BSCurve2d([[0., 0.5], [1.5, 0.6], [3., 0.6]], [0., 1.], [3, 3], 2)
    shroud = gbs.BSCurve2d([[0., 1.5], [1.5, 1.4], [3., 1.4]], [0., 1.], [3, 3], 2)
    # 41 stations from inlet to outlet, 11 stream lines from hub to shroud
    pts, ni, nj, _, _ = yams.mesh_channel([hub, shroud], [0., 1.], 11, 41)

    blade = yams.BladeInfo()
    blade.i1, blade.i2 = 15, 25
    blade.mode = yams.MeridionalBladeMode.DIRECT
    blade.omg = 300. * math.pi / 30  # 300 rpm
    k1, k2 = math.radians(-50.), math.radians(-20.)
    solver_case = yams.make_solver_case(
        gbs.make_structuredgrid(pts, ni, nj),
        [(blade, lambda m, l: (1 - m) * k1 + m * k2)],
    )

    solver_case.max_geom = 300
    solver_case.tol_rel_pos = 1e-5
    solver_case.relocate = True
    solver_case.mf_uniform = True
    solver_case.gi.rho_cst = False
    solver_case.gi.RF = 0.05
    solver_case.gi.R = 287.04

    solver_case.inlet.mode = yams.MeridionalBC.INLET_VmMoy_Ts_Ps_Vu
    solver_case.inlet.Vm_moy = 30.
    solver_case.inlet.Vu = lambda l_rel: 0.
    solver_case.inlet.Ps = lambda l_rel: 1.e5
    solver_case.inlet.Ts = lambda l_rel: 288.15 + 20. * l_rel  # non-uniform total enthalpy

    g = solver_case.gi.g
    for i, j in np.ndindex(g.nRows(), g.nCols()):
        g(i, j).Cp = 1004.5
    return solver_case, blade


def test_rotor_direct(rotor_direct_case):
    solver_case, blade = rotor_direct_case
    g = solver_case.gi.g

    yams.curvature_solver(solver_case)

    ni, nj = g.nRows(), g.nCols()
    assert np.all(np.isfinite([g(i, j).Vm for i in range(ni) for j in range(nj)]))
    # Converged before max_geom
    assert len(solver_case.log.delta_pos_max) < solver_case.max_geom

    # Physics: the flow leaves at the metal angle, rothalpy is conserved through the
    # rotor, the rotor gives work on every stream line, increasing with radius
    assert all(g(blade.i2, j).bet == approx(math.radians(-20.)) for j in range(nj))
    i_in = blade.i1 - 1
    assert all(g(blade.i2, j).I == approx(g(i_in, j).I, rel=1e-12) for j in range(nj))
    work = [g(blade.i2, j).H - g(i_in, j).H for j in range(nj)]
    assert all(w > 0 for w in work)
    assert all(b > a for a, b in zip(work, work[1:]))

    # Non-regression
    assert len(solver_case.log.delta_pos_max) == REF_ROTOR_GEOM_ITERATIONS
    assert [g(ni - 1, j).Vm for j in range(nj)] == approx(REF_ROTOR_OUTLET_VM, rel=REL)
    assert [g(ni - 1, j).Vu for j in range(nj)] == approx(REF_ROTOR_OUTLET_VU, rel=REL)
    assert work == approx(REF_ROTOR_WORK, rel=REL)

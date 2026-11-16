from pygbs import gbs
from pyams import yams

def test_issue_18():
    # Params
    points = [[0.0, -1.6771781385954805], [0.8087768025549643, 1.6207941488717532]]
    n_computation_planes = 2
    n_blades = 2
    j_le = 2
    j_te = 2

    # Curve
    poles = [[5.2e-18, 0.033], [1.25, 0.244]]
    knots = [0.0, 0.4, 0.60, 0.809, 0.958, 1.258]
    degree = 3

    stream = gbs.BSCurve2d(poles, knots, degree)

    # Yams Solver
    yams.BladeToBladeCurvatureSolver(
        points,
        n_computation_planes, 
        stream,
        n_blades,
        j_le,
        j_te,
    )
# -*- coding: utf-8 -*-

import sympy
from ..salcs import (
    calc_salcs_projection,
    calc_salcs_func,
    _expand_irreducible,
    _angles_to_vectors,

)
from ..tables import (
    tables,
    headers,
    mulliken,
    rot_trans_modes,
    IR_active,
    Raman_active,
    masks,
    atom_contribution
)


def test_calc_salcs_projection():
    # ammonia hydrogens
    a, b, c = sympy.symbols('a b c')
    assert (calc_salcs_projection([a, b, c, a, b, c], 'c3v') ==
            [2*a + 2*b + 2*c, 0, 2*a - b - c])

    # trigonal bipyramidal
    a1, a2, e1, e2, e3 = sympy.symbols('a1, a2, e1, e2, e3')
    assert (calc_salcs_projection([e1, e2, e3, -e1, -e2, -e3, -e1,
                                   -e2, -e3, e1, e2, e3], 'd3h') ==
            [0, 0, 0, 0, 4*e1 + 4*e2 + 4*e3, 4*e1 - 2*e2 - 2*e3])
    assert (calc_salcs_projection([a1, a1, a1, -a2, -a2, -a2, -a2,
                                   -a2, -a2, a1, a1, a1], 'd3h') ==
            [6*a1 - 6*a2, 0, 0, 0, 6*a1 + 6*a2, 0])


def test_calc_salcs_func():
    a, b, c, d = sympy.symbols('a b c d')

    # salc_true = [a + b + c, 0, [a - 0.5*b - 0.5*c, b - c, a -
    #              0.5*b - 0.5*c, -b + c], 0, 0, 0]
    # assert (calc_salcs_func([[0, 0], [120, 0], [240, 0]], [a, b, c],
    #                         'd3h', mode='angle') == salc_true)

    # square planar
    salc_true = [a + b + c + d, 0, a - b + c - d, 0, 0, 0, 0, 0, 0,
                 [a - c, b - d]]
    assert (calc_salcs_func([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
                            [a, b, c, d], 'd4h', mode='vector') == salc_true)

    # trigonal bipyramidal
    a1, a2, e1, e2, e3 = sympy.symbols('a1, a2, e1, e2, e3')
    salc_true = [[1.0*e1 + 1.0*e2 + 1.0*e3, 1.0*a1 + 1.0*a2], 0,
                 [1.0*e1 - 0.5*e2 - 0.5*e3, 1.0*e2 - 1.0*e3,
                  1.0*e1 - 0.5*e2 - 0.5*e3, -1.0*e2 + 1.0*e3], 0,
                 1.0*a1 - 1.0*a2, 0]
    angles = [[0, 0], [120, 0], [240, 0], [0, 90], [0, -90]]
    assert(calc_salcs_func(angles, [e1, e2, e3, a1, a2], 'd3h', mode='angle')
           == salc_true)

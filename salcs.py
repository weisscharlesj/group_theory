# -*- coding: utf-8 -*-

"""
Contains chemical group theory functions for calculating symmetry adapted
linear combinations (SALCs) or group orbitals using either the projection
operator method or using the symmetry functions in the character tables.
"""

from math import cos, sin, radians
import numpy as np

from .tables import (
    tables,
    symmetry_func_dict,
    table_coeff
)


# PROJECTION OPERATOR METHOD

def _expand_irreducible(irred, group):
    """
    Return expanded irreducible.

    Expands irreducible from the common condensed form to the full form. This
    entails repeating characters for every equivalent operation. Foe example,
    C3v retruns the rotation twice and reflection thrice.

    Parameters
    ----------
    irred : tuple
        Condensed irreducible representation from character table.
    group : str
        Point group (e.g., 'C2v').

    Returns
    -------
    Expanded irreducible representation as a list.

    Example
    ------
    >>> _expand_irreducible([2, -1, 0], 'c3v')
    >>> [2, -1, -1, 0, 0, 0]

    """
    expanded_irred = []
    for i in range(len(irred)):
        expanded_irred.extend([irred[i]] * table_coeff[group.lower()][i])

    return expanded_irred


def calc_salcs_projection(projection, group):
    """
    Return SALCs using projection operator method.

    Given the projections of orbitals as a result of a point group symmetry
    operations, returns the SALCs. This is a two-step process.
    1. Provide all ligands or outer atoms a variable name (e.g., a, b, c, etc.)
    2. Track one orbital to see how it transforms after each symmeter operation
    3. Provide a list of strings of the results

    Parameters
    ----------
    projection : List of SymPy symbols
        Results of projection operations for symmetry operations of group.
    group : str
        Point group (e.g., 'C2v').

    Returns
    -------
    Nested list of strings of the SALCS for each irreducible representation.
    Returns None for irreducibles with no SALC.

    Example
    -------
    >>> import sympy
    >>> a, b, c = sympy.symbols('a b c')
    >>> projection_operator([a, b, c, a, b, c], 'C3v')
    >>> [2*a + 2*b + 2*c, , 0, 2*a - b - c]

    """
    salcs = []

    for irred in tables[group.lower()]:
        product = np.array(_expand_irreducible(irred, group.lower()) *
                           np.array(projection))
        salcs.append(np.sum(product))

    return salcs


# USING SYMMETRY FUNCTIONS

def _angles_to_vectors(ligand_angles):
    """
    Calculate xyz vectors from angles around central atom.

    Given angles for outer ligands/atoms with respect to x-axis and elevation,
    this functiuon returns a list of [x, y, z] unit vectors.

    Parameters
    ----------
    ligand_angles : nested list or tuple of numbers
        Angle of each outer atom/ligand in degrees with respect
        azimuth/yaw(phi) and pitch/elevation(theta)

    Returns
    -------
    Nested Numpy array containing [x,y,z] vectors.

    Example
    -------
    >>> _angles_to_vectors([[0, 0], [90, 0], [180, 0], [-90, 0]])
    >>> array([[ 1.0,  0.0, 0.0],
               [ 0.0,  1.0,  0.0],
               [-1.0,  0.0,  0.0],
               [ 0.0, -1.0,  0.0]])

    """
    ligand_vectors = []
    for orbital in ligand_angles:
        x = cos(radians(orbital[0]))*cos(radians(orbital[1]))
        y = sin(radians(orbital[0]))*cos(radians(orbital[1]))
        z = sin(radians(orbital[1]))
        ligand_vectors.append([x, y, z])

    return ligand_vectors


def _eval_sym_func(coords, funcs):
    """
    Evaluate symmetry function.

    Evalutes symmetr functions using the supplied series of xyz coordinates.
    If all values evaluate as zeros, the irreducible has no SALC, and 0 is
    returned.

    Parameters
    ----------
    coords : array-like nested containing values in threes
        xyz coordinates of ligand unit vectors.
    funcs : str
        The symmetry function supplied as a string or tuple of strings
        (e.g., 'x**2-y**2' or ('z**2', 'x**2+y**2')).

    Returns
    -------
    List or 0.

    """
    salcs = []

    for func in funcs:
        ligand_contribs = []
        for unit_vector in coords:
            x, y, z = unit_vector[0], unit_vector[1], unit_vector[2]
            ligand_contrib = eval(func, {'x': x, 'y': y, 'z': z})
            ligand_contribs.append(round(ligand_contrib, 2))

        if np.any(ligand_contribs):
            salcs.append(ligand_contribs)

    if not salcs:
        return 0
    elif len(salcs) == 1:
        return salcs[0]
    else:
        return salcs


def _normalize_salcs(salcs):
    """
    Normalize SALC.

    Normalizes SALC by dividing each SALC through by largest value in
    that SALC.

    Parameters
    ----------
    salcs : List or nested list
        Nested list of SALCS.

    Returns
    -------
    List or nested list.

    """
    normalized_values = []
    for value in salcs:
        if isinstance(value, list):
            normalized_values.append(_normalize_salcs(value))
        elif isinstance(value, int):
            normalized_values.append(value)
        else:
            normalized_values.append(round(value / max(salcs), 2))

    return normalized_values


def calc_salcs_func(ligands, group, mode='vector'):
    """
    Return SALCs from symmetry functions.

    Returns SALCs (symmetry adapted linear combinations) using the character
    functions for a given point group. This requires the user to give the
    position of each ligand atomic orbital with respect to the x-axis, y-axis,
    and z-axis (mode='vector') or .

    Parameters
    ----------
    ligand : list or nested list
        Nested list of SALCs.
    group : str
        Point group (e.g., 'C2v').
    mode : 'vector' or 'angle'
        Whether the position of ligands or outer atoms are provided in 3D
        coordinates ('vector') or [angle from x-axis, elevation angle] pair
        ('angle').


    Returns
    -------
    Array of values indicating the weight and sign of each atomic orbital
    contribution to the SALC

    Examples
    --------
    >>> calc_salcs_func([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
                            'd4h', mode='vector')
    >>> np.array([[1, 1, 1, 1],
                  0,
                 [1, -1, 1, -1],
                 0, 0, 0, 0, 0, 0,
                 [[1, 0, -1, 0],
                  [0, 1, 0, -1]]])
    >>> calc_salcs_func([[0, 0], [120, 0], [240, 0]], 'd3d', mode='angle')
    >>> np.array([[1.0, 1.0, 1.0],
                  0,
                  [[1.0, -0.5, -0.5],
                   [0.0, 1.0, -1.0],
                   [1.0, -0.5, -0.5],
                   [0.0, -1.0, 1.0]]])
    """
    if mode == 'angle':
        ligand_vectors = _angles_to_vectors(ligands)
    elif mode == 'vector':
        ligand_vectors = ligands
    else:
        raise Exception("Invalide mode input: must be 'angle' or 'vector'")

    salcs = []
    for basis_func in symmetry_func_dict[group]:
        if basis_func == 0 or basis_func == 1:
            salcs.append(basis_func)
        else:
            salcs.append(_eval_sym_func(ligand_vectors, basis_func))

    return _normalize_salcs(salcs)

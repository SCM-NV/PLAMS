import numpy as np


try:
    from scipy.spatial.distance import cdist

    scipy_present = True
except ImportError:
    scipy_present = False

from scm.plams.tools.units import Units

__all__ = [
    "rotation_matrix",
    "axis_rotation_matrix",
    "distance_array",
    "angle",
    "dihedral",
    "cell_shape",
    "cellvectors_from_shape",
]

HALF_PI = np.pi / 2


def rotation_matrix(vec1, vec2):
    """
    Calculate the rotation matrix rotating *vec1* to *vec2*. Vectors can be any containers with 3 numerical values. They don't need to be normalized. Returns 3x3 numpy array.
    """
    vec1, vec2 = np.array(vec1), np.array(vec2)
    a = vec1 / np.linalg.norm(vec1)
    b = vec2 / np.linalg.norm(vec2)

    # avoid division by zero in case of antiparallel vectors
    if abs(1 + np.dot(a, b)) < 1e-8:
        return -np.eye(3)

    v1, v2, v3 = np.cross(a, b)
    M = np.array([[0, -v3, v2], [v3, 0, -v1], [-v2, v1, 0]])
    return np.identity(3) + M + np.dot(M, M) / (1 + np.dot(a, b))


def axis_rotation_matrix(vector, angle, unit="radian"):
    """
    Calculate the rotation matrix rotating along the *vector* by *angle* expressed in *unit*.

    *vector* can be any container with 3 numerical values. They don't need to be normalized. A positive angle denotes counterclockwise rotation, when looking along *vector*. Returns 3x3 numpy array.
    """

    vector /= np.linalg.norm(vector)
    v0, v1, v2 = vector

    W = np.array([[0, -v2, v1], [v2, 0, -v0], [-v1, v0, 0]])

    angle = Units.convert(angle, unit, "radian")
    a1 = np.sin(angle)
    a2 = 1.0 - np.cos(angle)

    return np.identity(3) + a1 * W + a2 * W @ W


def distance_array(array1, array2):
    """Calculates distance between each pair of points in *array1* and *array2*. Returns 2D ``numpy`` array.

    Uses fast ``cdist`` function if ``scipy`` is present, otherwise falls back to slightly slower ``numpy`` loop. Arguments should be 2-dimensional ``numpy`` arrays with the same second dimension. If *array1* is A x N and *array2* is B x N, the returned array is A x B.
    """
    return cdist(array1, array2) if scipy_present else np.array([np.linalg.norm(i - array2, axis=1) for i in array1])


def angle(vec1, vec2, result_unit="radian"):
    """Calculate an angle between vectors *vec1* and *vec2*.

    *vec1* and *vec2* should be iterable containers of length 3 (for example: tuple, list, numpy array). Values stored in them are expressed in Angstrom. Returned value is expressed in *result_unit*.

    This method requires all atomic coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
    """
    vec1 = np.array([*vec1], dtype=float)
    vec2 = np.array([*vec2], dtype=float)

    num = np.dot(vec1, vec2)
    den = np.sqrt(((vec1) ** 2).sum()) * np.sqrt(((vec2) ** 2).sum())
    return Units.convert(np.arccos(num / den), "radian", result_unit)


def dihedral(p1, p2, p3, p4, unit="radian"):
    """Calculate the value of diherdal angle formed by points *p1*, *p2*, *p3* and *p4* in a 3D space. Arguments can be any containers with 3 numerical values, also instances of |Atom|. Returned value is always non-negative, measures the angle clockwise (looking along *p2-p3* vector) and is expressed in *unit*."""
    p1 = np.array([*p1], dtype=float)
    p2 = np.array([*p2], dtype=float)
    p3 = np.array([*p3], dtype=float)
    p4 = np.array([*p4], dtype=float)

    b0 = p1 - p2
    b1 = p3 - p2
    b2 = p4 - p3

    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    ret = np.arctan2(y, x)
    ret = 2 * np.pi + ret if ret < 0 else ret
    return Units.convert(ret, "radian", unit)


def cell_shape(lattice):
    """
    Converts lattice vectors to lengths and angles (in radians)
    Sets internal cell size data, based on set of cell vectors.

    *lattice* is list containing three cell vectors (a 3x3 matrix)
    """
    lattice = np.asarray(lattice)
    a, b, c = np.sqrt((lattice**2).sum(axis=1))

    if a == 0.0 and b == 0.0 and c == 0.0:
        return

    alpha, beta, gamma = (90.0, 90.0, 90.0)

    if c != 0:
        alpha = angle(lattice[1], lattice[2])
        beta = angle(lattice[0], lattice[2])
    if b != 0:
        gamma = angle(lattice[0], lattice[1])

    return [a, b, c, alpha, beta, gamma]


def cell_lengths(lattice, unit="angstrom"):
    """Return the lengths of the lattice vector. Returns a list with the same length as the number of lattice vector."""

    if lattice is None or len(lattice) == 0:
        raise ValueError("Cannot calculate cell_lengths for nonperiodic system")
    lattice = np.asarray(lattice)
    ret = np.sqrt((lattice**2).sum(axis=1)) * Units.conversion_ratio("angstrom", unit)
    return ret.tolist()


def cell_angles(lattice, unit="degree"):
    """Return the angles between lattice vectors.

    unit : str
        output unit

    For 2D systems, returns a list [gamma]

    For 3D systems, returns a list [alpha, beta, gamma]
    """
    ndim = len(lattice)

    if ndim < 2:
        raise ValueError(
            "Cannot calculate cell_angles for fewer than 2 lattice vectors. Tried with {} lattice vectors".format(ndim)
        )

    gamma = angle(lattice[0], lattice[1], result_unit=unit)

    if ndim == 2:
        return [gamma]

    if ndim >= 3:
        alpha = angle(lattice[1], lattice[2], result_unit=unit)
        beta = angle(lattice[0], lattice[2], result_unit=unit)
        return [alpha, beta, gamma]


def cellvectors_from_shape(box):
    """
    Converts lengths and angles (in radians) of lattice vectors to the lattice vectors.

    *box* should be an iterable of length 3, containing the lengths a, b, c, or an iterable of length 6 additionally containing the angles alpha, beta, gamma
    """
    a = box[0]
    b = box[1]
    c = box[2]
    alpha, beta, gamma = HALF_PI, HALF_PI, HALF_PI
    if len(box) == 6:
        alpha = box[3]
        beta = box[4]
        gamma = box[5]

    va = [a, 0.0, 0.0]
    vb = [b * np.cos(gamma), b * np.sin(gamma), 0.0]

    cx = c * np.cos(beta)
    cy = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) * c / np.sin(gamma)
    volume = 1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2
    volume += 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
    volume = np.sqrt(volume)
    cz = c * volume / np.sin(gamma)
    vc = [cx, cy, cz]

    lattice = [va, vb, vc]

    return lattice

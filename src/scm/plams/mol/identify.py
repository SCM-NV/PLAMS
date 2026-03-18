from collections import OrderedDict
from itertools import combinations
from typing import Optional, Dict, Tuple, TYPE_CHECKING, Sequence, List, Any, Union

import numpy as np
from numpy.typing import ArrayLike

from scm.plams.core.functions import requires_optional_package
from scm.plams.core.private import sha256
from scm.plams.tools.units import Units

if TYPE_CHECKING:
    from scm.plams.mol.atom import Atom
    from scm.plams.mol.bond import Bond
    from scm.plams.mol.molecule import Molecule
    from networkx import Graph

__all__ = ["label_atoms"]


possible_flags = ["BO", "RS", "EZ", "DH", "CO", "H2"]


def twist(v1: ArrayLike, v2: ArrayLike, v3: ArrayLike, tolerance: Optional[float] = None) -> Tuple[int, Optional[int]]:
    """
    Given 3 vectors in 3D space measure their "chirality" with *tolerance*.

    Returns a pair. The first element is an integer number measuring the orientation (clockwise vs counterclockwise) of *v1* and *v3* while looking along *v2*. Values 1 and -1 indicate this case and the second element of returned pair is ``None``. Value 0 indicates that *v1*, *v2*, and *v3* are coplanar, and the second element of the returned pair is indicating if two turns made by going *v1*->*v2*->*v3* are the same (left-left, right-right) or the opposite (left-right, right-left).
    """
    tolerance = 1e-2 if tolerance is None else tolerance
    v1 = np.asarray(v1)
    v2 = np.asarray(v2)
    v3 = np.asarray(v3)
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    v3 /= np.linalg.norm(v3)
    x = np.dot(v3, np.cross(v1, v2))
    if abs(x) <= tolerance:  # v1, v2, v3 are coplanar
        return 0, int(np.sign(np.dot(np.cross(v1, v2), np.cross(v2, v3))))
    return int(np.sign(x)), None


def bend(v1: ArrayLike, v2: ArrayLike, tolerance: Optional[float] = None) -> int:
    """Check if two vectors in 3D space are parallel or perpendicular, with *tolerance* (in degrees).

    Returns 1 if *v1* and *v2* are collinear, 2 if they are perpendicular, 0 otherwise."""
    tolerance = 7.5 if tolerance is None else tolerance
    v1 = np.asarray(v1)
    v2 = np.asarray(v2)
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    angle = Units.convert(abs(np.arccos(np.dot(v1, v2))), "rad", "deg")
    if angle <= tolerance:
        return 1
    if abs(angle - 90.0) <= tolerance:
        return 2
    return 0


def unique_atoms(atomlist: Union[Sequence["Atom"], "Molecule"]) -> List["Atom"]:
    """Filter *atomlist* (list or |Molecule|) for atoms with unique ``IDname``."""
    d = {}
    for atom in atomlist:
        if atom.IDname not in d:  # type: ignore[attr-defined]
            d[atom.IDname] = 0  # type: ignore[attr-defined]
        d[atom.IDname] += 1  # type: ignore[attr-defined]
    return [atom for atom in atomlist if d[atom.IDname] == 1]  # type: ignore[attr-defined]


def initialize(molecule: "Molecule") -> None:
    """Initialize atom labeling algorithm by setting ``IDname`` and ``IDdone`` attributes for all atoms in *molecule*."""
    for at in molecule:
        at.IDname = at.symbol  # type: ignore[attr-defined]
        at.IDdone = False  # type: ignore[attr-defined]


def clear(molecule: "Molecule") -> None:
    """Remove ``IDname`` and ``IDdone`` attributes from all atoms in *molecule*."""
    for at in molecule:
        if hasattr(at, "IDname"):
            del at.IDname  # type: ignore[attr-defined]
        if hasattr(at, "IDdone"):
            del at.IDdone  # type: ignore[attr-defined]


def iterate(molecule: "Molecule", flags: Dict[str, Any]) -> bool:
    """Perform one iteration of atom labeling algorithm.

    First, mark all atoms that are unique and have only unique neighbors as "done". Then calculate new label for each atom that is not done. Return True if the number of different atom labels increased during this iteration.
    """
    names = len(set(at.IDname for at in molecule))  # type: ignore[attr-defined]
    unique = set(unique_atoms(molecule))

    for atom in molecule:
        if atom in unique and all(N in unique for N in atom.neighbors()):
            atom.IDdone = True  # type: ignore[attr-defined]
        if not atom.IDdone:  # type: ignore[attr-defined]
            atom.IDnew = new_name(atom, flags)  # type: ignore[attr-defined]

    for atom in molecule:
        if not atom.IDdone:  # type: ignore[attr-defined]
            atom.IDname = atom.IDnew  # type: ignore[attr-defined]

    new_names = len(set(atom.IDname for atom in molecule))  # type: ignore[attr-defined]
    return new_names > names  # True means this iteration increased the number of distinct names


def new_name(atom: "Atom", flags: Dict[str, Any]) -> str:
    """Compute new label for *atom*.

    The new label is based on the existing label of *atom*, labels of all its neighbors and (possibly) some additional conformational information. The labels of neighbors are not obtained directly by reading neighbor's ``IDname`` but rather by a process called "knocking". The *atom* knocks all its bonds. Each knocked bond returns an identifier describing the atom on the other end of the bond. The identifier is composed of knocked atom's ``IDname`` together with some additional information desribing the character of the bond and knocked atom's spatial environment. The exact behavior of this mechanism is adjusted by the contents of *flags* dictionary (see :func:`label_atoms` for details).
    """

    knocks = [knock(atom, bond, flags) for bond in atom.bonds]
    knocks.sort(key=lambda x: x[0])
    labels = set(i[0] for i in knocks)  # types of neighbors

    more = []
    if flags["RS"] and len(knocks) == 4 and len(labels) == 4:  # 4 different neighbors
        v1 = atom.vector_to(knocks[0][1])
        v2 = atom.vector_to(knocks[1][1])
        v3 = atom.vector_to(knocks[2][1])
        more.append("RS" + str(twist(v1, v2, v3, flags.get("twist_tol"))))

    if flags["CO"] and len(knocks) >= 4:
        d: OrderedDict[str, List[Atom]] = OrderedDict()
        for label, at in knocks:
            if label not in d:
                d[label] = []
            d[label].append(at)
        for label in d:
            if len(d[label]) > 1:
                angles = [
                    bend(atom.vector_to(a), atom.vector_to(b), flags.get("bend_tol"))
                    for a, b in combinations(d[label], 2)
                ]
                angles.sort()
            else:  # len(d[label]) == 1
                v1 = atom.vector_to(d[label][0])
                angles = []
                for k in d:
                    if k != label:
                        angles.append(
                            sorted(bend(v1, atom.vector_to(a), flags.get("bend_tol")) for a in d[k])  # type: ignore[arg-type]
                        )
            more.append("CO" + str(angles))

    return sha256("|".join([atom.IDname] + [i[0] for i in knocks] + more))  # type: ignore[attr-defined]


def knock(A: "Atom", bond: "Bond", flags: Dict[str, Any]) -> Tuple[str, "Atom"]:
    """Atom *A* knocks one of its bonds.

    *bond* has to be a bond formed by atom *A*. The other end of this bond (atom S) returns its description, consisting of its ``IDname`` plus, possibly, some additional information. If *BO* flag is set, the description includes the bond order of *bond*. If *EZ* flag is set, the description includes additional bit of information whenever E/Z isomerism is possible. If *DH* flag is set, the description includes additional information for all dihedrals A-S-N-F such that A is a unique neighbor of S and F is a unique neighbor of N.
    """

    S = bond.other_end(A)
    ret: str = S.IDname  # type: ignore[attr-defined]

    if flags["BO"] and bond.order != 1:
        ret += "BO" + str(bond.order)

    if flags["H2"]:
        S_nbors = S.neighbors()
        S_nbors.remove(A)
        if len(S_nbors) == 3 and len(unique_atoms(S_nbors)) == 3:
            S_nbors.sort(key=lambda x: x.IDname)  # type: ignore[attr-defined]
            v1, v2, v3 = [S.vector_to(i) for i in S_nbors]
            t = twist(v1, v2, v3, flags.get("twist_tol"))
            ret += "H2" + str(t)

    if flags["EZ"] or flags["DH"]:
        S_unique = unique_atoms(S.neighbors())
        if A in S_unique:  # *A* is a unique neighbor of *S*
            S_unique.remove(A)
            S_unique.sort(key=lambda x: x.IDname)  # type: ignore[attr-defined]
            for b in S.bonds:
                N = b.other_end(S)
                if N in S_unique:
                    N_unique = unique_atoms(N.neighbors())
                    if S in N_unique:
                        N_unique.remove(S)
                    N_unique.sort(key=lambda x: x.IDname)  # type: ignore[attr-defined]
                    if N_unique:
                        F = N_unique[0]
                        v1 = A.vector_to(S)
                        v2 = S.vector_to(N)
                        v3 = N.vector_to(F)
                        t = twist(v1, v2, v3, flags.get("twist_tol"))
                        if flags["DH"]:
                            ret += "DH" + str(t)
                        elif flags["EZ"] and b.order == 2 and t[0] == 0:
                            # A-S=N-F are coplanar
                            ret += "EZ" + str(t[1])
                        break

    return (ret, S)


def label_atoms(molecule: "Molecule", **kwargs: Any) -> "Molecule":
    """Label atoms in *molecule*.

    Boolean keyword arguments:

    *   *BO* -- include bond orders
    *   *RS* -- include R/S stereoisomerism
    *   *EZ* -- include E/Z stereoisomerism
    *   *DH* -- include some dihedrals to detect alkane rotamers and syn-anti conformers in cycloalkanes
    *   *CO* -- include more spatial info to detect different conformation of coordination complexes (flat square, octaedr etc.)

    Numerical keyword arguments:

    *   *twist_tol* -- tolerance for :func:`twist` function
    *   *bend_tol* -- tolerance for :func:`bend` function

    Diherdals considered with *DH* are all the dihedrals A-B-C-D such that A is a unique neighbor or B and D is a unique neighbor of C.

    For atoms with 4 or more neighbors, the *CO* flag includes information about relative positions of equivalent/non-equivalent neighbors by checking if vectors from the central atom to the neighbors form 90 or 180 degrees angles.
    """
    initialize(molecule)
    while iterate(molecule, kwargs):
        pass
    return molecule


def molecule_name(molecule: "Molecule") -> str:
    """Compute the label of the whole *molecule* based on ``IDname`` attributes of all the atoms."""
    names = [atom.IDname for atom in molecule]  # type: ignore[attr-defined]
    names.sort()
    return sha256(" ".join(names))


@requires_optional_package("networkx")
def get_graph(mol: "Molecule", dic: Dict[str, Any], level: int = 1) -> Optional["Graph"]:
    """
    Create a networkx graph for this molecule that can be used to compare (all info is in the edge.weight attribute)
    """
    import networkx

    nats = len(mol)
    mol = mol.copy()
    if len(mol.bonds) == 0:
        mol.guess_bonds()
    if not hasattr(mol.atoms[0], "IDname"):
        mol.label(level=level, keep_labels=True)

    # Get the connectivity matrix (remove bond orders)
    matrix = mol.bond_matrix()
    matrix[matrix > 0] = 1
    matrix = matrix.astype(np.int32)

    # Multiply the graph entries with the unique labels for each atom
    identifiers = np.array([dic[at.IDname] if at.IDname in dic.keys() else None for at in mol.atoms])  # type: ignore[attr-defined]
    if None in identifiers:
        return None
    identifiers = identifiers.astype(np.int32)
    matrix *= identifiers.reshape((1, nats))
    matrix *= identifiers.reshape((nats, 1))
    np.fill_diagonal(matrix, identifiers)

    # Create the graph
    graph = networkx.from_numpy_array(matrix)

    return graph

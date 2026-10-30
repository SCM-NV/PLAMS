"""Public reaction-network behavior, data preservation, and optional dependencies."""

import importlib
import sys
from typing import Any, Dict, Iterator, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

pytest.importorskip("networkx")

from scm.plams import AMSResults, Atom, Molecule, plot_reaction_network
from scm.plams.core.errors import MissingOptionalPackageError
from scm.plams.tools.plot_reaction_network_utils import _build_reaction_graph
from scm.plams.tools.units import Units


def add_state(
    landscape: AMSResults.EnergyLandscape,
    molecule: Molecule,
    energy: float,
    is_ts: bool = False,
    reactant: Optional[int] = None,
    product: Optional[int] = None,
) -> None:
    landscape._states.append(
        AMSResults.EnergyLandscape.State(
            landscape,
            "state",
            energy,
            molecule.copy(),
            1,
            is_ts,
            reactant,
            product,
            prefactorsFromReactant=10.0 if is_ts else None,
            prefactorsFromProduct=20.0 if is_ts else None,
        )
    )


@pytest.fixture
def landscape() -> AMSResults.EnergyLandscape:
    result = AMSResults.EnergyLandscape(None)
    molecule = Molecule()
    molecule.add_atom(Atom(symbol="H", coords=(0, 0, 0)))
    molecule.add_atom(Atom(symbol="H", coords=(0, 0, 0.74)))
    molecule.add_bond(molecule[1], molecule[2], 1)
    add_state(result, molecule, -0.10)
    molecule[2].coords = (0, 0, 4.0)
    add_state(result, molecule, -0.08)
    add_state(result, molecule, -0.05, True, 1, 2)
    return result


@pytest.fixture(autouse=True)
def close_figures() -> Iterator[None]:
    yield
    plt.close("all")


def test_graph_preserves_stationary_points_and_direction(landscape: AMSResults.EnergyLandscape) -> None:
    before = [(s.energy, len(s.molecule.bonds), s.molecule.as_array().copy()) for s in landscape]
    graph = _build_reaction_graph(landscape)
    assert set(graph) == {1, 2, 3}
    assert set(graph.edges()) == {(1, 3), (3, 2)}
    assert graph.nodes[1]["coarse_label"] == "H2"
    assert graph.nodes[2]["coarse_label"] == "H+H"
    assert graph.nodes[3]["prefactor_from_reactant"] == 10
    assert graph.nodes[3]["prefactor_from_product"] == 20
    for state, (energy, bonds, coords) in zip(landscape, before):
        assert state.energy == energy and len(state.molecule.bonds) == bonds
        np.testing.assert_array_equal(state.molecule.as_array(), coords)


def test_multiple_ts_and_missing_endpoints(landscape: AMSResults.EnergyLandscape) -> None:
    add_state(landscape, landscape[1].molecule, -0.04, True, 1, 2)
    add_state(landscape, landscape[1].molecule, -0.03, True, 1, None)
    graph = _build_reaction_graph(landscape)
    assert set(graph.edges()) == {(1, 3), (3, 2), (1, 4), (4, 2), (1, 5)}
    assert graph.out_degree(5) == 0
    ax = plot_reaction_network(landscape)
    assert len(ax.patches) == 5


def test_original_ids_and_supplied_axes(landscape: AMSResults.EnergyLandscape) -> None:
    selected = landscape.select_states([3, 1, 2], keep_original_ids=True)
    fig, ax = plt.subplots(figsize=(5, 4))
    before = fig.get_size_inches().copy()
    assert plot_reaction_network(selected, ax=ax) is ax
    np.testing.assert_array_equal(fig.get_size_inches(), before)
    assert {t.get_text() for t in ax.texts} == {"1", "2", "3"}
    assert _build_reaction_graph(selected).nodes[1]["label"] == "3"


@pytest.mark.parametrize("color_by", ["formula", "energy"])
def test_empty_landscape(color_by: str) -> None:
    ax = plot_reaction_network(AMSResults.EnergyLandscape(None), color_by=color_by)
    assert not ax.collections and not ax.patches


def test_recoloring_preserves_geometry(landscape: AMSResults.EnergyLandscape) -> None:
    positions = {1: (0, 0), 2: (2, 0), 3: (1, 0)}
    axes = [plot_reaction_network(landscape, positions=positions, color_by=mode) for mode in ("formula", "energy")]
    for ax in axes:
        ax.figure.canvas.draw()
    np.testing.assert_array_equal(axes[0].figure.get_size_inches(), axes[1].figure.get_size_inches())
    np.testing.assert_array_equal(axes[0].get_position().bounds, axes[1].get_position().bounds)
    for first, second in zip(axes[0].collections, axes[1].collections):
        np.testing.assert_array_equal(first.get_offsets(), second.get_offsets())
        np.testing.assert_array_equal(first.get_sizes(), second.get_sizes())
    for first, second in zip(axes[0].patches, axes[1].patches):
        np.testing.assert_array_equal(first.get_path().vertices, second.get_path().vertices)
    assert len(axes[1].child_axes) == 1


def test_energy_scale_and_midpoint_arrows(landscape: AMSResults.EnergyLandscape) -> None:
    positions = {1: (0, 0), 2: (2, 0), 3: (1, 0)}
    ax = plot_reaction_network(
        landscape, positions=positions, color_by="energy", cmap="viridis", energy_reference=-0.1, energy_limits=(0, 2)
    )
    ax.figure.canvas.draw()
    factor = Units.convert(1.0, "hartree", "eV")
    np.testing.assert_allclose(ax.collections[0].get_facecolors(), plt.get_cmap("viridis")([0, 0.02 * factor / 2]))
    np.testing.assert_allclose(ax.collections[1].get_facecolors(), plt.get_cmap("viridis")([0.05 * factor / 2]))
    for (left, right), patch in zip([(1, 3), (3, 2)], ax.patches):
        vertices = patch.get_path().vertices
        np.testing.assert_allclose(vertices[0], positions[left], atol=1e-12)
        np.testing.assert_allclose(vertices[2], positions[right], atol=1e-12)
        np.testing.assert_allclose(vertices[3], (np.array(positions[left]) + positions[right]) / 2, atol=1e-12)
        notch = vertices[5]
        assert np.dot(vertices[3] - notch, np.array(positions[right]) - positions[left]) > 0


def test_constant_energies_and_units(landscape: AMSResults.EnergyLandscape) -> None:
    for state in landscape:
        state.energy = -1.0
    ax = plot_reaction_network(landscape, color_by="energy", unit="kcal/mol")
    ax.figure.canvas.draw()
    assert np.isfinite(ax.collections[0].get_facecolors()).all()
    np.testing.assert_allclose(ax.collections[0].get_facecolors()[0], ax.collections[1].get_facecolors()[0])
    assert ax.child_axes[0].get_xlabel() == "Relative energy (kcal/mol)"


@pytest.mark.parametrize("dpi", [72, 100, 200])
@pytest.mark.parametrize("node", [1, 3])
@pytest.mark.parametrize("font_family", ["sans-serif", "serif"])
def test_long_labels_fit(landscape: AMSResults.EnergyLandscape, dpi: int, node: int, font_family: str) -> None:
    landscape[node].originalID = 123456789
    with plt.rc_context({"figure.dpi": dpi, "font.family": font_family}):
        ax = plot_reaction_network(landscape)
    ax.figure.canvas.draw()
    text = next(t for t in ax.texts if t.get_text() == "123456789")
    bbox = text.get_window_extent(ax.figure.canvas.get_renderer())
    diameter = np.sqrt(ax.collections[1 if node == 3 else 0].get_sizes()[0]) * ax.figure.dpi / 72
    padding = 3.6 * ax.figure.dpi / 72
    required = (
        max(bbox.width, bbox.height) + padding if node == 3 else np.hypot(bbox.width + padding, bbox.height + padding)
    )
    assert diameter >= required


def test_native_fallback_and_periodic_requirement(
    landscape: AMSResults.EnergyLandscape, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(sys.modules, "scm.base", None)
    graph = _build_reaction_graph(landscape)
    assert graph.nodes[1]["coarse_label"] == "H2"
    assert graph.nodes[2]["coarse_label"] == "H+H"
    landscape[1].molecule.lattice = [(10, 0, 0), (0, 10, 0)]
    with pytest.raises(MissingOptionalPackageError, match="scm.base"):
        plot_reaction_network(landscape)


def test_surface_selection() -> None:
    result = AMSResults.EnergyLandscape(None)
    mol = Molecule()
    mol.add_atom(Atom(symbol="Cu", coords=(0, 0, 0)))
    mol[1].properties.region = "surface"
    mol.add_atom(Atom(symbol="H", coords=(0, 0, 1.5)))
    for energy in (-1.0, -0.9):
        add_state(result, mol, energy)
    add_state(result, mol, -0.8, True, 1, 2)
    mol[2].coords = (0, 0, 8)
    add_state(result, mol, -0.7)
    add_state(result, mol, -0.6, True, 1, 4)
    graph = _build_reaction_graph(result, surface_only=True)
    assert set(graph) == {1, 2, 3}
    assert graph.nodes[1]["coarse_label"] == "H*"
    assert set(graph.edges()) == {(1, 3), (3, 2)}


@pytest.mark.parametrize("package", ["matplotlib", "networkx", "scipy"])
def test_optional_dependencies(
    landscape: AMSResults.EnergyLandscape, monkeypatch: pytest.MonkeyPatch, package: str
) -> None:
    functions = importlib.import_module("scm.plams.core.functions")
    original = functions.find_spec
    monkeypatch.setattr(functions, "find_spec", lambda name: None if name == package else original(name))
    with pytest.raises(MissingOptionalPackageError, match=package):
        plot_reaction_network(landscape, layout="tree")
    if package == "scipy":
        ax = plot_reaction_network(landscape, layout="tree", positions={1: (0, 0), 2: (2, 0), 3: (1, 0)})
        assert len(ax.patches) == 2


@pytest.mark.parametrize(
    "kwargs",
    [
        {"size_scale": 0},
        {"label_padding": -1},
        {"node_spacing": -1},
        {"layout": "unknown"},
        {"layout_iterations": 0},
        {"within_class_spacing": 0},
        {"between_class_spacing": 0},
        {"guess_bonds_dmax": 0},
        {"surface_only": True, "surface_region": None},
        {"edge_gray": -1},
        {"edge_alpha": 2},
        {"color_by": "unknown"},
        {"positions": {1: (0, 0)}},
        {"positions": {1: (0, 0), 2: (0, float("nan")), 3: (1, 1)}},
        {"color_by": "energy", "energy_reference": float("nan")},
        {"color_by": "energy", "energy_limits": (2, 1)},
    ],
)
def test_invalid_options(landscape: AMSResults.EnergyLandscape, kwargs: Dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        plot_reaction_network(landscape, **kwargs)


def test_tree_layout(landscape: AMSResults.EnergyLandscape) -> None:
    ax = plot_reaction_network(landscape, layout="tree", layout_iterations=3)
    assert len(ax.patches) == 2
    assert all(np.isfinite(collection.get_offsets()).all() for collection in ax.collections)

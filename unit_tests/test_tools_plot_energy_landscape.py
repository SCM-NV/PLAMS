"""Behavioral checks for the public energy-landscape plotter and its helpers."""

from typing import List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from PIL import Image

from scm.plams import AMSResults, Atom, Molecule
from scm.plams.tools.plot import plot_energy_landscape
from scm.plams.tools.units import Units


@pytest.fixture
def landscape() -> AMSResults.EnergyLandscape:
    result = AMSResults.EnergyLandscape(None)
    molecule = Molecule()
    molecule.add_atom(Atom(symbol="H"))
    states: List[Tuple[float, bool, Optional[int], Optional[int]]] = [
        (-1.0, False, None, None),
        (-0.8, False, None, None),
        (-0.7, True, 1, 2),
    ]
    for energy, is_ts, reactant, product in states:
        result._states.append(
            AMSResults.EnergyLandscape.State(result, "state", energy, molecule.copy(), 1, is_ts, reactant, product)
        )
    return result


@pytest.mark.parametrize("layout", ["auto", "dfs", "bfs", "longest_path", "force", "crossings"])
def test_energy_landscape_layouts_and_highlighting(landscape: AMSResults.EnergyLandscape, layout: str) -> None:
    fig, ax = plt.subplots()
    try:
        returned = plot_energy_landscape(landscape, ax=ax, layout=layout, force_iterations=3, highlight_states=[1, 3])
        assert returned is ax
        assert {text.get_text() for text in ax.texts} == {"1", "2", "3"}
        assert len(ax.collections) == 5  # Three levels plus two highlighted underlays.
        assert len(ax.lines) == 3  # Two connectors plus one highlighted underlay.
        assert ax.get_ylabel() == "Relative Energy (eV)"
        assert ax.get_ylim()[0] < 0
        assert ax.get_ylim()[1] > Units.convert(0.3, "hartree", "eV")
    finally:
        plt.close(fig)


def test_energy_landscape_disconnected_spacing(landscape: AMSResults.EnergyLandscape) -> None:
    landscape._states.append(AMSResults.EnergyLandscape.State(landscape, "isolated", -0.9, Molecule(), 1, False))
    fig, ax = plt.subplots()
    try:
        plot_energy_landscape(landscape, ax=ax, layout="bfs", spacing=2.0, landscape_width=0.6)
        positions = {int(text.get_text()): text.get_position()[0] for text in ax.texts}
        assert positions == {1: 0.0, 3: 2.0, 2: 4.0, 4: 7.0}
    finally:
        plt.close(fig)


def test_energy_landscape_empty() -> None:
    fig, ax = plt.subplots()
    try:
        returned = plot_energy_landscape(AMSResults.EnergyLandscape(None), ax=ax, unit="kcal/mol")
        assert returned is ax
        assert ax.get_ylabel() == "Relative Energy (kcal/mol)"
        assert not ax.collections and not ax.lines and not ax.texts
    finally:
        plt.close(fig)


def test_energy_landscape_invalid_layout(landscape: AMSResults.EnergyLandscape) -> None:
    fig, ax = plt.subplots()
    try:
        with pytest.raises(ValueError, match="Unsupported layout"):
            plot_energy_landscape(landscape, ax=ax, layout="invalid")  # type: ignore[arg-type]
    finally:
        plt.close(fig)


@pytest.mark.parametrize("mode", ["RGB", "RGBA", "L"])
def test_energy_landscape_view_image_crop(
    landscape: AMSResults.EnergyLandscape, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    from scm.plams.tools import view as view_module

    image = Image.new(mode, (12, 10), (0, 0, 0, 0) if mode == "RGBA" else "white")
    image.paste("black", (3, 2, 9, 6))
    calls = []

    def fake_view(molecule: Molecule, **kwargs: object) -> Image.Image:
        calls.append(kwargs)
        return image.copy()

    monkeypatch.setattr(view_module, "view", fake_view)
    selected = landscape.select_states([2], keep_original_ids=True)
    fig, ax = plt.subplots()
    try:
        plot_energy_landscape(
            selected,
            ax=ax,
            show_molecules=True,
            molecule_plot_kwargs={"width": 128, "guess_bonds": False},
            molecule_plot_kwargs_by_state={2: {"width": 96}},
        )
        assert calls == [{"width": 96, "guess_bonds": False}]
        assert [text.get_text() for text in ax.texts] == ["2"]
        assert len(ax.child_axes) == 1
        assert ax.child_axes[0].images[0].get_array().shape[:2] == (4, 6)
        assert np.asarray(image).shape[:2] == (10, 12)
    finally:
        plt.close(fig)

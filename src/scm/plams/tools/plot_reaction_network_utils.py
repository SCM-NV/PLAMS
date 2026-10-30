"""Private graph, layout and rendering helpers for ``plot_reaction_network``.

Optional plotting/graph dependencies are imported only when the public function
is called. No molecules, landscapes or user-supplied positions are modified.
"""

from collections import defaultdict
import hashlib
import math
from typing import (
    Any,
    DefaultDict,
    Dict,
    List,
    Literal,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
    TYPE_CHECKING,
    cast,
)

import networkx as nx
import numpy as np

from scm.plams.core.errors import MissingOptionalPackageError
from scm.plams.core.functions import requires_optional_package
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.tools.periodic_table import PT

if TYPE_CHECKING:
    import matplotlib
    from matplotlib.backend_bases import RendererBase
    from scm.plams.interfaces.adfsuite.ams import AMSResults
    from scm.plams.mol.molecule import Molecule

__all__: List[str] = []
Position = Tuple[float, float]


def _bond_adjacency(
    molecule: "Molecule", surface_region: Optional[str], dmax: float
) -> Tuple[Set[int], Dict[int, Set[int]]]:
    """Guess bonds on a separate chemical system, or on a nonperiodic PLAMS copy."""
    try:
        from scm.base import ChemicalSystem
    except ImportError:
        if molecule.lattice:
            raise MissingOptionalPackageError("scm.base")
        clone = molecule.copy()
        if clone.atoms:
            clone.guess_bonds(dmax=dmax)
        atom_index = {id(atom): index for index, atom in enumerate(clone.atoms)}
        neighbors: Dict[int, Set[int]] = {i: set() for i in range(len(clone))}
        surface = set()
        for i, atom in enumerate(clone.atoms):
            regions = atom.properties.get("region", [])
            if isinstance(regions, str):
                regions = regions.replace(",", " ").split()
            if surface_region is not None and regions and surface_region in regions:
                surface.add(i)
        for bond in clone.bonds:
            left, right = atom_index[id(bond.atom1)], atom_index[id(bond.atom2)]
            neighbors[left].add(right)
            neighbors[right].add(left)
        return surface, neighbors
    # Same conversion as scm.utils, without depending on that separate package.
    system = ChemicalSystem(AMSJob(molecule=molecule).get_input())
    surface = set(int(i) for i in system.get_atoms_in_region(surface_region)) if surface_region else set()
    system.guess_bonds(dmax=dmax)
    return surface, {i: set(int(j) for j in system.bonds.get_bonded_atoms(i)) for i in range(len(system))}


def _coarse_info(
    state: "AMSResults.EnergyLandscape.State", surface_region: Optional[str], dmax: float
) -> Dict[str, Any]:
    """Classify connected nonsurface fragments by formula and surface attachment."""
    molecule = state.molecule
    surface: Set[int]
    neighbors: Dict[int, Set[int]]
    if not molecule.atoms:
        surface, neighbors = set(), {}
    else:
        surface, neighbors = _bond_adjacency(molecule, surface_region, dmax)
    nonsurface = [i for i in range(len(molecule)) if i not in surface]
    visited: Set[int] = set()
    fragments: List[Dict[str, Any]] = []
    for start in nonsurface:
        if start in visited:
            continue
        stack = [start]
        component = []
        touches_surface = False
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            touches_surface = touches_surface or bool(neighbors[current] & surface)
            stack.extend(neighbors[current] - surface - visited)
        symbols: DefaultDict[str, int] = defaultdict(int)
        for i in component:
            symbols[molecule.atoms[i].symbol] += 1
        order = (
            ([s for s in ("C", "H") if s in symbols] + sorted(s for s in symbols if s not in {"C", "H"}))
            if "C" in symbols
            else sorted(symbols)
        )
        formula = "".join(symbol + (str(symbols[symbol]) if symbols[symbol] != 1 else "") for symbol in order)
        if touches_surface:
            formula += "*"
        fragments.append(
            dict(
                label=formula,
                total_mass=sum(PT.get_mass(molecule.atoms[i].atnum) for i in component),
                touches_surface=touches_surface,
            )
        )
    fragments.sort(key=lambda f: (-f["total_mass"], f["label"]))
    return dict(
        label="+".join(f["label"] for f in fragments) if fragments else "∅",
        fragment_count=len(fragments),
        fragment_masses=tuple(f["total_mass"] for f in fragments),
        has_surface_fragment=any(f["touches_surface"] for f in fragments),
        all_fragments_surface=bool(fragments) and all(f["touches_surface"] for f in fragments),
    )


def _build_reaction_graph(
    energy_landscape: "AMSResults.EnergyLandscape",
    surface_region: Optional[str] = "surface",
    surface_only: bool = False,
    guess_bonds_dmax: float = 1.35,
) -> nx.DiGraph:
    """Represent every stationary point explicitly, retaining directional TS data."""
    if not math.isfinite(guess_bonds_dmax) or guess_bonds_dmax <= 0:
        raise ValueError("guess_bonds_dmax must be finite and positive")
    if surface_only and not surface_region:
        raise ValueError("surface_only requires a surface_region")
    states = list(energy_landscape)
    state_by_id = {state.id: state for state in states}
    coarse = {state.id: _coarse_info(state, surface_region, guess_bonds_dmax) for state in states if not state.isTS}
    graph = nx.DiGraph()
    for node, state in state_by_id.items():
        graph.add_node(
            node,
            kind="ts" if state.isTS else "minimum",
            energy=state.energy,
            formula=state.molecule.get_formula(False),
            count=state.count,
            label=str(getattr(state, "display_id", node)),
        )
        if not state.isTS:
            info = coarse[node]
            graph.nodes[node].update(
                coarse_label=info["label"],
                coarse_fragment_count=info["fragment_count"],
                coarse_fragment_masses=info["fragment_masses"],
                has_surface_fragment=info["has_surface_fragment"],
                all_fragments_surface=info["all_fragments_surface"],
            )
            continue
        for endpoint in (state.reactantsID, state.productsID):
            if endpoint is not None and endpoint not in coarse:
                raise ValueError(f"Transition state {node} refers to an absent minimum {endpoint}")
        left = coarse.get(state.reactantsID, {}) if state.reactantsID is not None else {}
        right = coarse.get(state.productsID, {}) if state.productsID is not None else {}
        graph.nodes[node].update(
            reactant_coarse_label=left.get("label"),
            product_coarse_label=right.get("label"),
            coarse_label=f"{left['label']} -> {right['label']}" if left and right else None,
            reactant_all_fragments_surface=left.get("all_fragments_surface", False),
            product_all_fragments_surface=right.get("all_fragments_surface", False),
            prefactor_from_reactant=state.prefactorsFromReactant,
            prefactor_from_product=state.prefactorsFromProduct,
        )
    for node, state in state_by_id.items():
        if state.isTS:
            if state.reactantsID is not None:
                graph.add_edge(state.reactantsID, node)
            if state.productsID is not None:
                graph.add_edge(node, state.productsID)
    if surface_only:
        kept = {
            node
            for node, data in graph.nodes(data=True)
            if (data["kind"] == "minimum" and data["all_fragments_surface"])
            or (
                data["kind"] == "ts"
                and data["reactant_all_fragments_surface"]
                and data["product_all_fragments_surface"]
            )
        }
        graph = graph.subgraph(kept).copy()
        graph.remove_nodes_from(list(nx.isolates(graph)))
    return graph


def _stable_angle_and_radius(label: str, node_id: int) -> Tuple[float, float]:
    digest = hashlib.sha256(f"{label}:{node_id}".encode()).digest()
    angle = 2.0 * math.pi * int.from_bytes(digest[:8], "big") / 2**64
    radius = 0.03 + 0.12 * int.from_bytes(digest[8:16], "big") / 2**64
    return angle, radius


def _build_coarse_layout_graph(graph: nx.DiGraph) -> Tuple[Dict[str, np.ndarray], Dict[str, List[int]]]:
    coarse_graph = nx.Graph()
    minimum_nodes_by_label = defaultdict(list)

    for node, data in graph.nodes(data=True):
        if data["kind"] == "minimum":
            label = data["coarse_label"]
            minimum_nodes_by_label[label].append(node)
            coarse_graph.add_node(label)

    for _, data in graph.nodes(data=True):
        if data.get("kind") != "ts":
            continue
        reactant_label = data.get("reactant_coarse_label")
        product_label = data.get("product_coarse_label")
        if reactant_label is None or product_label is None:
            continue
        coarse_graph.add_node(reactant_label)
        coarse_graph.add_node(product_label)
        if reactant_label == product_label:
            coarse_graph.add_edge(reactant_label, product_label, weight=1.0)
        elif coarse_graph.has_edge(reactant_label, product_label):
            coarse_graph[reactant_label][product_label]["weight"] += 1.0
        else:
            coarse_graph.add_edge(reactant_label, product_label, weight=1.0)

    if coarse_graph.number_of_nodes() == 0:
        return {}, minimum_nodes_by_label

    coarse_positions = nx.spring_layout(
        coarse_graph,
        seed=7,
        iterations=200,
        weight="weight",
    )
    return coarse_positions, minimum_nodes_by_label


def _build_detailed_layout_graph(graph: nx.DiGraph, minimum_nodes_by_label: Mapping[str, List[int]]) -> nx.Graph:
    layout_graph = nx.Graph()
    layout_graph.add_nodes_from(graph.nodes(data=True))

    for source, target in graph.edges():
        layout_graph.add_edge(source, target, weight=1.0)

    for node_ids in minimum_nodes_by_label.values():
        if len(node_ids) < 2:
            continue
        anchor = min(node_ids)
        for node_id in node_ids:
            if node_id != anchor:
                layout_graph.add_edge(anchor, node_id, weight=8.0)

    for node, data in graph.nodes(data=True):
        if data["kind"] != "ts":
            continue
        for label in {data.get("reactant_coarse_label"), data.get("product_coarse_label")}:
            if label is None:
                continue
            node_ids = minimum_nodes_by_label.get(label, [])
            if node_ids:
                layout_graph.add_edge(node, min(node_ids), weight=4.0)

    return layout_graph


def _compute_initial_positions(
    graph: nx.DiGraph, coarse_positions: Mapping[str, np.ndarray], minimum_nodes_by_label: Mapping[str, List[int]]
) -> Dict[int, Position]:
    initial_positions = {}

    for label, node_ids in minimum_nodes_by_label.items():
        center_x, center_y = coarse_positions[label]
        for node_id in node_ids:
            angle, radius = _stable_angle_and_radius(label, node_id)
            initial_positions[node_id] = (
                center_x + radius * math.cos(angle),
                center_y + radius * math.sin(angle),
            )

    for node, data in graph.nodes(data=True):
        if data["kind"] != "ts":
            continue
        reactant_pos = coarse_positions.get(data.get("reactant_coarse_label"))
        product_pos = coarse_positions.get(data.get("product_coarse_label"))
        if reactant_pos is not None and product_pos is not None:
            center_x = 0.5 * (reactant_pos[0] + product_pos[0])
            center_y = 0.5 * (reactant_pos[1] + product_pos[1])
        elif reactant_pos is not None:
            center_x, center_y = reactant_pos
        elif product_pos is not None:
            center_x, center_y = product_pos
        else:
            center_x, center_y = 0.0, 0.0

        angle, radius = _stable_angle_and_radius("ts", node)
        initial_positions[node] = (
            center_x + 0.6 * radius * math.cos(angle),
            center_y + 0.6 * radius * math.sin(angle),
        )

    return initial_positions


@requires_optional_package("scipy")
def _compute_layout_positions(graph: nx.DiGraph, seed: int = 11, iterations: int = 500) -> Dict[int, Position]:
    coarse_positions, minimum_nodes_by_label = _build_coarse_layout_graph(graph)
    if not coarse_positions:
        return nx.spring_layout(graph.to_undirected(), seed=seed)

    initial_positions = _compute_initial_positions(graph, coarse_positions, minimum_nodes_by_label)
    layout_graph = _build_detailed_layout_graph(graph, minimum_nodes_by_label)
    return nx.spring_layout(
        layout_graph,
        k=3.0 / math.sqrt(float(graph.number_of_nodes())),
        pos=initial_positions,
        seed=seed,
        iterations=iterations,
        weight="weight",
    )


def _make_coarse_label_colors(graph: nx.DiGraph) -> Dict[str, str]:
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_hex

    cycle = [to_hex(plt.get_cmap("tab20")(i)) for i in range(20)]
    label_info = {}
    for _, data in graph.nodes(data=True):
        if data["kind"] == "minimum":
            label = data["coarse_label"]
            label_info[label] = {
                "label": label,
                "fragment_count": data["coarse_fragment_count"],
                "fragment_masses": data["coarse_fragment_masses"],
            }

    minimum_labels = sorted(label_info, key=lambda label: _coarse_label_sort_key(label_info[label]))
    return {label: cycle[index % len(cycle)] for index, label in enumerate(minimum_labels)}


def _coarse_label_sort_key(label_info: Mapping[str, Any]) -> Tuple[Any, ...]:
    return (
        tuple(-mass for mass in label_info["fragment_masses"]),
        label_info["fragment_count"],
        label_info["label"],
    )


@requires_optional_package("scipy")
def _tree_preserving_layout(
    graph: nx.DiGraph,
    reference_positions: Mapping[int, Sequence[float]],
    within_class_spacing: float = 0.016,
    between_class_spacing: float = 1.8,
    iterations: int = 150,
) -> Dict[int, Position]:
    """Spread nodes locally while retaining the reference layout's branching shape.

    Class centers keep their exact reference directions and relative distances.
    ``between_class_spacing`` uniformly scales those centers about the origin;
    it is a positive multiplier, unlike the gap argument of ``grouped_layout``.
    ``within_class_spacing`` is a target minimum separation in coordinate units,
    rather than a guaranteed separation in the rendered image. Local repulsion
    relaxes minima without replacing their geometry with a regular spiral.
    Same-class TS nodes keep their reference angles around their class, while
    inter-class TS nodes retain their reference locations along the branches.

    No graph attributes or edges are modified. This is a network visualization,
    not a conversion to a spanning tree; cycles and directional links are retained.
    Requires scipy for local nearest-neighbor queries.
    """
    from scipy.spatial import cKDTree

    if not math.isfinite(within_class_spacing) or within_class_spacing <= 0:
        raise ValueError("within_class_spacing must be finite and positive")
    if not math.isfinite(between_class_spacing) or between_class_spacing <= 0:
        raise ValueError("between_class_spacing must be finite and positive")
    if iterations < 0:
        raise ValueError("iterations must be non-negative")
    reference = {n: np.asarray(reference_positions[n], dtype=float) for n in graph}
    if any(x.shape != (2,) or not np.isfinite(x).all() for x in reference.values()):
        raise ValueError("Reference coordinates must be finite 2D positions")
    groups = defaultdict(list)
    for node, data in graph.nodes(data=True):
        if data["kind"] == "minimum":
            groups[data["coarse_label"]].append(node)
    centers = {label: np.mean([reference[n] for n in nodes], axis=0) for label, nodes in groups.items()}
    positions = {node: xy * between_class_spacing for node, xy in reference.items()}
    radii = {}
    for label, nodes in groups.items():
        nodes = sorted(nodes)
        local = np.array([reference[n] - centers[label] for n in nodes])
        old_radius = float(np.linalg.norm(local, axis=1).max())
        radius = max(old_radius, within_class_spacing * math.sqrt(len(nodes)) / 1.75)
        if old_radius:
            local *= radius / old_radius
        # Preserve the initial angular structure; distribute radial crowding.
        lengths = np.linalg.norm(local, axis=1)
        if len(nodes) > 1 and old_radius:
            order = np.argsort(lengths, kind="stable")
            ranks = np.empty(len(nodes), dtype=int)
            ranks[order] = np.arange(len(nodes))
            targets = radius * np.sqrt((ranks + 0.5) / len(nodes))
            nonzero = lengths > 1e-12
            local[nonzero] *= (targets[nonzero] / lengths[nonzero])[:, None]
        local -= local.mean(axis=0)
        # Only close pairs interact, avoiding a full-network spring calculation.
        for _ in range(iterations):
            pairs = np.array(sorted(cKDTree(local).query_pairs(within_class_spacing)))
            if not len(pairs):
                break
            delta = local[pairs[:, 1]] - local[pairs[:, 0]]
            distance = np.linalg.norm(delta, axis=1)
            coincident = distance < 1e-12
            delta[coincident] = [1e-8, 0]
            distance[coincident] = 1e-8
            push = delta * (0.3 * (within_class_spacing - distance) / distance)[:, None]
            shift = np.zeros_like(local)
            np.add.at(shift, pairs[:, 0], -push)
            np.add.at(shift, pairs[:, 1], push)
            local += shift
            local -= local.mean(axis=0)
        radii[label] = float(np.linalg.norm(local, axis=1).max())
        for node, xy in zip(nodes, local):
            positions[node] = centers[label] * between_class_spacing + xy
    for node, data in graph.nodes(data=True):
        if data["kind"] != "ts":
            continue
        label = data.get("reactant_coarse_label")
        if label != data.get("product_coarse_label") or label not in centers:
            continue
        offset = reference[node] - centers[label]
        length = float(np.linalg.norm(offset))
        direction = offset / length if length else np.array([1.0, 0.0])
        # Keep the original TS angular pattern in a thin halo outside the minima.
        halo_width = within_class_spacing * 2.5
        fraction = 0.5 + math.atan(length * 10) / math.pi
        positions[node] = centers[label] * between_class_spacing + direction * (
            radii[label] + within_class_spacing + fraction * halo_width
        )
    return {node: tuple(xy) for node, xy in positions.items()}


def _adaptive_node_areas(node_count: int, size_scale: float = 1.0) -> Tuple[float, float]:
    """Return minimum/TS marker areas in pt², bounded for sparse/dense networks.

    A bounded power law decreases marker area as network density increases.
    Small networks have a larger TS/minimum area ratio to accommodate state IDs.
    size_scale multiplies marker area, rather than radius or diameter.
    These are base areas; visible labels may require larger individual markers.
    """
    if node_count < 0:
        raise ValueError("node_count must be non-negative")
    if not math.isfinite(size_scale) or size_scale <= 0:
        raise ValueError("size_scale must be finite and positive")
    base = min(1300.0, max(100.0, 4400.0 / max(node_count, 1) ** 0.4))
    minimum = min(2400.0, max(50.0, base * size_scale))
    ratio = 1 / 6 + (0.45 - 1 / 6) * max(0.0, 1 - node_count / 40)
    return minimum, minimum * ratio


def _label_fitted_area(
    label: str,
    base_area: float,
    shape: str,
    font_size: float,
    padding: float = 2.0,
    renderer: Optional["RendererBase"] = None,
) -> float:
    """Enlarge a marker to enclose its text with padding, in physical points.

    Text metrics use the same bold sans-serif font as the rendered labels.
    Raster metrics account for pixel hinting when a renderer is available;
    vector metrics remain a lower bound and support canvases without a renderer.
    Squares must contain the padded text rectangle; circles must contain its diagonal.
    The label-fit minimum takes precedence over the density-based size cap.
    """
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextToPath

    if shape not in {"o", "s"}:
        raise ValueError("shape must be 'o' or 's'")
    if not all(math.isfinite(x) for x in (base_area, font_size, padding)):
        raise ValueError("Area, font size and padding must be finite")
    if base_area <= 0 or font_size <= 0 or padding < 0:
        raise ValueError("Area and font size must be positive; padding must be non-negative")
    if not label:
        return base_area
    font = FontProperties(family="sans-serif", size=font_size, weight="bold")
    width, height, _ = TextToPath().get_text_width_height_descent(label, font, ismath=False)
    if renderer is not None:
        rendered_width, rendered_height, _ = renderer.get_text_width_height_descent(label, font, ismath=False)
        pixels_per_point = float(np.asarray(renderer.points_to_pixels(1.0)).item())
        width = max(width, rendered_width / pixels_per_point)
        height = max(height, rendered_height / pixels_per_point)
    width += 2 * padding + 0.5
    # Full line height also accommodates font ascenders/descenders and export backends.
    height = max(height, font_size) + 2 * padding + 0.5
    diameter = max(width, height) if shape == "s" else math.hypot(width, height)
    return max(base_area, diameter**2)


def _connectivity_layout(graph: nx.DiGraph, seed: int = 11) -> Dict[int, Tuple[float, float]]:
    """Layer acyclic networks by connectivity; use spring layout for cycles.

    Every directed edge is retained. Disconnected nodes share appropriate layers.
    This does not infer reaction direction from energies or remove cycles.
    """
    if not graph:
        return {}
    if not nx.is_directed_acyclic_graph(graph):
        return {n: tuple(p) for n, p in nx.spring_layout(graph, seed=seed, iterations=200).items()}
    positions: Dict[int, Position] = {}
    for level, nodes in enumerate(nx.topological_generations(graph)):
        nodes = sorted(nodes)
        nodes.sort(
            key=lambda n: np.mean([positions[p][1] for p in graph.predecessors(n)]) if graph.in_degree(n) else 0.0
        )
        for index, node in enumerate(nodes):
            positions[node] = (float(level), (len(nodes) - 1) / 2 - index)
    return positions


def _midpoint_arrow(
    path: "matplotlib.path.Path", mutation_size: float, linewidth: float, aspect_ratio: float = 1.0
) -> Tuple[List["matplotlib.path.Path"], List[bool]]:
    """Retain the shaft and add a slender, notched head halfway along its curve.

    FancyArrowPatch calls this in display coordinates, so arrowhead dimensions
    remain in physical points through zooming, resizing and DPI changes.
    """
    from matplotlib.path import Path

    vertices = np.asarray(path.vertices, dtype=float)
    codes = np.asarray(path.codes) if path.codes is not None else None
    start, end = vertices[0], vertices[-1]
    if len(vertices) == 3 and codes is not None and codes[1] == Path.CURVE3:
        midpoint = (start + 2 * vertices[1] + end) / 4
        tangent = end - start
    elif len(vertices) == 4 and codes is not None and codes[1] == Path.CURVE4:
        midpoint = (start + 3 * vertices[1] + 3 * vertices[2] + end) / 8
        tangent = -start - vertices[1] + vertices[2] + end
    else:
        midpoint = (start + end) / 2
        tangent = end - start
    length = float(np.linalg.norm(tangent))
    if length < 1e-12:
        return [path], [False]
    direction = tangent / length
    normal = np.array([-direction[1], direction[0]])
    base = midpoint - 0.65 * mutation_size * direction
    half_width = 0.16 * mutation_size * normal
    notch = midpoint - 0.44 * mutation_size * direction
    head = Path(
        [midpoint, base + half_width, notch, base - half_width, midpoint],
        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY],
    )
    return [path, head], [False, True]


def _plot_reaction_graph(
    graph: nx.DiGraph,
    positions: Optional[Mapping[int, Sequence[float]]] = None,
    ax: Optional["matplotlib.axes.Axes"] = None,
    size_scale: float = 1.0,
    label_states: Optional[bool] = None,
    edge_gray: float = 0.11666666666666667,
    edge_alpha: float = 0.7066666666666667,
    title: str = "Reaction network",
    label_padding: float = 2.0,
    node_spacing: float = 14.0,
    color_by: Literal["formula", "energy"] = "formula",
    cmap: str = "jet",
    layout: str = "auto",
    seed: int = 11,
    layout_iterations: int = 500,
    within_class_spacing: float = 0.014,
    between_class_spacing: float = 1.8,
    energy_unit: str = "eV",
    energy_reference: Optional[float] = None,
    energy_limits: Optional[Tuple[float, float]] = None,
) -> "matplotlib.axes.Axes":
    """Draw the network with marker sizes and canvas adapted to its node count.

    Formula classes set minimum colors by default; squares denote transition states.
    Plain state IDs are shown by default for at most 40 nodes; marker shape
    distinguishes minima from transition states. Position coordinates have no
    energetic meaning. An existing ax is reused without changing its figure size.
    Saving is left to the caller. Large networks without supplied positions
    use a formula-grouped, deterministic spring-layout backbone.
    When labels are visible, each marker grows as needed to enclose its label
    with ``label_padding`` points of clearance. This applies at any network size.
    For automatically created figures with at most 40 nodes, ``node_spacing``
    specifies a minimum target gap between marker boundaries in physical points.
    The canvas expands as needed; node sizes, labels and positions stay fixed.
    A supplied ax keeps its existing figure size. Unlabeled large-network figures
    retain their compact tree presentation. Edge shafts connect node centers
    beneath opaque marker fills; one arrowhead marks the midpoint of each edge.

    Set ``color_by='energy'`` to color both minima and TS nodes using the same
    colormap and normalization. Graph ``energy`` values must be finite hartree
    energies. Relative energies use the lowest plotted state as reference unless
    ``energy_reference`` is supplied in hartree. ``energy_limits`` are relative
    values in ``energy_unit`` and allow a shared scale across different plots.
    A colorbar is added without changing the main axes or layout. Edges retain
    their grayscale style. Energy-colored nodes are opaque to match the colorbar;
    label colors adapt to the palette for contrast.
    """
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    if not 0 <= edge_gray <= 1 or not 0 <= edge_alpha <= 1:
        raise ValueError("edge_gray and edge_alpha must be between 0 and 1")
    if not math.isfinite(label_padding) or label_padding < 0:
        raise ValueError("label_padding must be finite and non-negative")
    if not math.isfinite(node_spacing) or node_spacing < 0:
        raise ValueError("node_spacing must be finite and non-negative")
    if color_by not in {"formula", "energy"}:
        raise ValueError("color_by must be 'formula' or 'energy'")
    if layout not in {"auto", "connectivity", "tree"}:
        raise ValueError("layout must be 'auto', 'connectivity' or 'tree'")
    if not isinstance(layout_iterations, int) or layout_iterations < 1:
        raise ValueError("layout_iterations must be a positive integer")
    if not math.isfinite(within_class_spacing) or within_class_spacing <= 0:
        raise ValueError("within_class_spacing must be finite and positive")
    if not math.isfinite(between_class_spacing) or between_class_spacing <= 0:
        raise ValueError("between_class_spacing must be finite and positive")
    _adaptive_node_areas(len(graph), size_scale)
    if not graph:
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 3.8))
        ax.set_axis_off()
        ax.set_title(title, loc="left")
        return ax
    if positions is not None:
        if not set(graph).issubset(positions):
            raise ValueError("positions must include every plotted state ID")
        position_arrays = {n: np.asarray(positions[n], dtype=float) for n in graph}
        if any(xy.shape != (2,) or not np.isfinite(xy).all() for xy in position_arrays.values()):
            raise ValueError("positions must contain finite 2D coordinates")
        positions = {n: (float(xy[0]), float(xy[1])) for n, xy in position_arrays.items()}
    energy_mappable = None
    energy_colors: Dict[int, np.ndarray] = {}
    if color_by == "energy":
        from scm.plams.tools.units import Units

        energies = {n: float(data["energy"]) for n, data in graph.nodes(data=True)}
        if not energies or not all(math.isfinite(e) for e in energies.values()):
            raise ValueError("Energy coloring requires finite energies for every node")
        reference = min(energies.values()) if energy_reference is None else float(energy_reference)
        if not math.isfinite(reference):
            raise ValueError("energy_reference must be finite")
        factor = Units.convert(1.0, "hartree", energy_unit)
        relative = {n: (e - reference) * factor for n, e in energies.items()}
        if energy_limits is None:
            low, high = min(relative.values()), max(relative.values())
            if low == high:
                low, high = low - 0.5, high + 0.5
        else:
            low, high = energy_limits
            if not all(math.isfinite(e) for e in (low, high)) or low >= high:
                raise ValueError("energy_limits must contain two finite increasing values")
        energy_mappable = ScalarMappable(norm=Normalize(low, high), cmap=plt.get_cmap(cmap))
        energy_colors = {n: energy_mappable.to_rgba(np.asarray([e], dtype=float))[0] for n, e in relative.items()}
    count = len(graph)
    min_area, ts_area = _adaptive_node_areas(count, size_scale)
    small = count <= 40
    owns_figure = ax is None
    if ax is None:
        figsize = (10, 3.8) if count <= 15 else (10, 7) if small else (18, 16)
        _, ax = plt.subplots(figsize=figsize)
        if small:
            ax.figure.subplots_adjust(left=0.04, right=0.97, top=0.78, bottom=0.22)
        else:
            ax.figure.subplots_adjust(left=0.02, right=0.81, bottom=0.04, top=0.94)
    if positions is None:
        if layout == "connectivity" or (layout == "auto" and small):
            positions = _connectivity_layout(graph, seed=seed)
        else:
            positions = _tree_preserving_layout(
                graph,
                _compute_layout_positions(graph, seed, layout_iterations),
                within_class_spacing=within_class_spacing,
                between_class_spacing=between_class_spacing,
            )
    if label_states is None:
        label_states = small
    minima = [n for n, d in graph.nodes(data=True) if d["kind"] == "minimum"]
    ts = [n for n, d in graph.nodes(data=True) if d["kind"] == "ts"]
    colors = _make_coarse_label_colors(graph)
    minimum_labels = {n: str(graph.nodes[n].get("label", n)) for n in minima}
    ts_labels = {n: str(graph.nodes[n].get("label", n)) for n in ts}
    get_renderer = getattr(ax.figure.canvas, "get_renderer", None)
    renderer = get_renderer() if label_states and get_renderer is not None else None
    minimum_areas = [
        _label_fitted_area(minimum_labels[n], min_area, "o", 10, label_padding, renderer) if label_states else min_area
        for n in minima
    ]
    ts_areas = [
        _label_fitted_area(ts_labels[n], ts_area, "s", 8, label_padding, renderer) if label_states else ts_area
        for n in ts
    ]
    area_by_node = dict(zip(minima + ts, minimum_areas + ts_areas))
    if positions:
        xy = np.array(list(positions.values()))
        lo, hi = xy.min(axis=0), xy.max(axis=0)
        pad = np.maximum((hi - lo) * (0.09 if small else 0.04), 0.7 if small else 0.05)
        ax.set_xlim(lo[0] - pad[0], hi[0] + pad[0])
        ax.set_ylim(lo[1] - pad[1], hi[1] + pad[1])
    ax.set_aspect("equal")
    if owns_figure and small and node_spacing > 0:
        ax.apply_aspect()
        nodes = list(graph)
        centers = ax.transData.transform([positions[n] for n in nodes]) * 72 / ax.figure.dpi
        expansion = 1.0
        for i, left in enumerate(nodes):
            for j in range(i):
                right = nodes[j]
                delta = centers[i] - centers[j]
                distance = float(np.linalg.norm(delta))
                if distance < 1e-12:
                    continue  # Resizing cannot separate coincident supplied positions.
                direction = np.abs(delta / distance)
                supports = []
                for node in (left, right):
                    radius = math.sqrt(area_by_node[node]) / 2
                    supports.append(radius * direction.sum() if graph.nodes[node]["kind"] == "ts" else radius)
                expansion = max(expansion, (sum(supports) + node_spacing) / distance)
        if expansion > 1:
            figure = ax.figure
            assert isinstance(figure, Figure)
            width, height = figure.get_size_inches()
            figure.set_size_inches(float(width * expansion), float(height * expansion), forward=True)
    for source, target in graph.edges():
        ax.add_patch(
            FancyArrowPatch(
                (float(positions[source][0]), float(positions[source][1])),
                (float(positions[target][0]), float(positions[target][1])),
                # Matplotlib accepts callables, but its stubs only list ArrowStyle.
                arrowstyle=cast("matplotlib.patches.ArrowStyle", _midpoint_arrow),
                shrinkA=0,
                shrinkB=0,
                mutation_scale=15 if small else 4,
                connectionstyle="arc3,rad=0.0" if small else "arc3,rad=0.035",
                linewidth=1.3 if small else 0.35,
                color=str(edge_gray),
                alpha=edge_alpha,
                zorder=1,
                joinstyle="miter",
                capstyle="round",
            )
        )
    nx.draw_networkx_nodes(
        graph,
        positions,
        ax=ax,
        nodelist=minima,
        node_size=minimum_areas,
        node_color=[
            energy_colors[n] if color_by == "energy" else colors[graph.nodes[n]["coarse_label"]] for n in minima
        ],
        edgecolors="white",
        linewidths=0.8 if small else 0.25,
    )
    nx.draw_networkx_nodes(
        graph,
        positions,
        ax=ax,
        nodelist=ts,
        node_size=ts_areas,
        node_color=(
            [energy_colors[n] for n in ts]
            if color_by == "energy"
            else [0.7 * np.array([0.6] * 3) + 0.3 * np.array(ax.get_facecolor()[:3])] * len(ts)
        ),
        node_shape="s",
        edgecolors="white",
        linewidths=0.6 if small else 0.1,
    )
    if label_states:
        for labels, font_size, default_color in ((minimum_labels, 10, "white"), (ts_labels, 8, "0.15")):
            label_groups: Dict[str, Dict[int, str]] = {}
            for node, text in labels.items():
                color = default_color
                if color_by == "energy":
                    rgb = np.array(energy_colors[node][:3])
                    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
                    luminance = float(linear @ np.array([0.2126, 0.7152, 0.0722]))
                    color = "black" if luminance > 0.179 else "white"
                label_groups.setdefault(color, {})[node] = text
            for color, subset in label_groups.items():
                nx.draw_networkx_labels(
                    graph,
                    positions,
                    ax=ax,
                    labels=subset,
                    font_color=color,
                    font_size=font_size,
                    font_family="sans-serif",
                    font_weight="bold",
                )
    ax.set_axis_off()
    ax.set_title(title, loc="left", fontsize=13 if small else 16, pad=20 if small else 18)
    if color_by == "energy":
        assert energy_mappable is not None
        handles: List[Union[Line2D, Patch]] = [
            Line2D([], [], marker="o", linestyle="", color="0.5", markersize=8, label="Minima"),
            Line2D([], [], marker="s", linestyle="", color="0.5", markersize=7, label="Transition states"),
        ]
        if small:
            ax.legend(handles=handles, loc="center", bbox_to_anchor=(0.22, -0.15), ncol=2, frameon=False, fontsize=9)
            colorbar_ax = ax.inset_axes((0.55, -0.20, 0.4, 0.04))
            colorbar = ax.figure.colorbar(energy_mappable, cax=colorbar_ax, orientation="horizontal")
        else:
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=9)
            colorbar_ax = ax.inset_axes((1.04, 0.18, 0.03, 0.6))
            colorbar = ax.figure.colorbar(energy_mappable, cax=colorbar_ax)
        colorbar.set_label(f"Relative energy ({energy_unit})", fontsize=9)
        colorbar.ax.tick_params(labelsize=8)
    elif small:
        handles = [
            Line2D([], [], marker="o", linestyle="", color=c, markersize=8, label=label) for label, c in colors.items()
        ]
        handles.append(Line2D([], [], marker="s", linestyle="", color="0.6", markersize=7, label="Transition states"))
        ax.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.05),
            ncol=min(4, len(handles)),
            frameon=False,
            fontsize=9,
        )
    else:
        handles = [
            Patch(facecolor=c, label=f"{label} ({sum(graph.nodes[n]['coarse_label'] == label for n in minima)})")
            for label, c in colors.items()
        ]
        handles.append(Patch(facecolor="0.6", label=f"Transition states ({len(ts)})"))
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            frameon=False,
            fontsize=9,
            title="Formula classes (minima)",
        )
    ax.figure.text(
        0.04,
        0.06 if small else 0.015,
        f"{len(minima)} minima · {len(ts)} transition states · " f"{graph.number_of_edges()} directed edges",
        fontsize=9,
        color="0.35",
    )
    return ax

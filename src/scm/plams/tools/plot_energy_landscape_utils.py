"""Private layout and image helpers for :func:`plot_energy_landscape`."""

from collections import deque
import itertools
from typing import Any, Callable, Dict, Iterable, List, Optional, Set, Tuple, TYPE_CHECKING

import numpy as np

from scm.plams.tools.units import Units

if TYPE_CHECKING:
    from scm.plams.interfaces.adfsuite.ams import AMSResults

__all__: List[str] = []


class _EnergyLandscapeLayout:
    """Keep connectivity and ordering state shared by the layout algorithms."""

    def __init__(
        self,
        states: List["AMSResults.EnergyLandscape.State"],
        spacing: float,
        unit: str,
        force_iterations: int,
    ) -> None:
        self.states = states
        self.spacing = spacing
        self.unit = unit
        self.force_iterations = force_iterations
        self.state_map = {state.id: state for state in self.states}
        # Build a graph representation of the landscape from the TS reactant/product links.
        self.adjacency: Dict[int, Set[int]] = {state.id: set() for state in self.states}
        self.edge_pairs: Set[Tuple[int, int]] = set()
        for state in self.states:
            if state.reactants is not None:
                self.adjacency[state.id].add(state.reactants.id)
                self.adjacency[state.reactants.id].add(state.id)
                left_id, right_id = sorted((state.id, state.reactants.id))
                self.edge_pairs.add((left_id, right_id))
            if state.products is not None:
                self.adjacency[state.id].add(state.products.id)
                self.adjacency[state.products.id].add(state.id)
                left_id, right_id = sorted((state.id, state.products.id))
                self.edge_pairs.add((left_id, right_id))

    def _state_sort_key(self, state_id: int) -> Tuple[bool, float, int]:
        state = self.state_map[state_id]
        return (state.isTS, state.energy, state.id)

    def _component(self, start_id: int, remaining_ids: Set[int]) -> Set[int]:
        component: Set[int] = set()
        stack = [start_id]
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend(self.adjacency[node] & remaining_ids)
        return component

    def _component_start(self, component: Set[int]) -> int:
        endpoints = [node for node in component if len(self.adjacency[node] & component) <= 1]
        candidates = endpoints if endpoints else list(component)
        return min(candidates, key=self._state_sort_key)

    def _connected_components(self) -> List[Set[int]]:
        # Layout disconnected reaction networks independently before placing them side by side.
        remaining_ids = {state.id for state in self.states}
        components: List[Set[int]] = []
        while remaining_ids:
            start_id = min(remaining_ids, key=self._state_sort_key)
            component = self._component(start_id, remaining_ids)
            components.append(component)
            remaining_ids -= component
        components.sort(key=lambda component: self._state_sort_key(self._component_start(component)))
        return components

    def _ordered_components(self, orderer: Callable[[Set[int]], List[int]]) -> List[int]:
        ordered_ids: List[int] = []
        for component in self._connected_components():
            ordered_ids.extend(orderer(component))
        return ordered_ids

    def _dfs_order(self, component: Set[int]) -> List[int]:
        # Follow one branch deeply before backtracking, which often matches a reaction path view.
        start_id = self._component_start(component)
        visited: Set[int] = set()
        ordered_ids: List[int] = []

        def _visit(node: int) -> None:
            visited.add(node)
            ordered_ids.append(node)
            neighbors = sorted(self.adjacency[node] & component, key=self._state_sort_key)
            for neighbor in neighbors:
                if neighbor not in visited:
                    _visit(neighbor)

        _visit(start_id)
        return ordered_ids

    def _bfs_order(self, component: Set[int]) -> List[int]:
        # Expand level by level from one endpoint to keep nearby states grouped together.
        start_id = self._component_start(component)
        queue = deque([start_id])
        visited: Set[int] = {start_id}
        ordered_ids: List[int] = []
        while queue:
            node = queue.popleft()
            ordered_ids.append(node)
            neighbors = sorted(self.adjacency[node] & component, key=self._state_sort_key)
            for neighbor in neighbors:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        return ordered_ids

    def _farthest(self, start_id: int, component: Set[int]) -> Tuple[int, Dict[int, int], Dict[int, Optional[int]]]:
        distances: Dict[int, int] = {start_id: 0}
        parents: Dict[int, Optional[int]] = {start_id: None}
        queue = deque([start_id])
        while queue:
            node = queue.popleft()
            for neighbor in sorted(self.adjacency[node] & component, key=self._state_sort_key):
                if neighbor not in distances:
                    distances[neighbor] = distances[node] + 1
                    parents[neighbor] = node
                    queue.append(neighbor)
        farthest_id = max(
            distances,
            key=lambda node: (distances[node], -int(self.state_map[node].isTS), -self.state_map[node].energy, -node),
        )
        return farthest_id, distances, parents

    def _longest_path_order(self, component: Set[int]) -> List[int]:
        # Use the graph backbone first, then attach side branches around that main path.
        start_id = self._component_start(component)
        end_a, _, _ = self._farthest(start_id, component)
        end_b, _, parents = self._farthest(end_a, component)

        backbone: List[int] = []
        node: Optional[int] = end_b
        while node is not None:
            backbone.append(node)
            node = parents[node]
        backbone.reverse()

        ordered_ids: List[int] = []
        placed: Set[int] = set()

        def _add_side_branch(root_id: int, blocked: Set[int]) -> None:
            neighbors = sorted((self.adjacency[root_id] & component) - blocked - placed, key=self._state_sort_key)
            for neighbor in neighbors:
                ordered_ids.append(neighbor)
                placed.add(neighbor)
                _add_side_branch(neighbor, blocked | {root_id})

        for node in backbone:
            if node not in placed:
                ordered_ids.append(node)
                placed.add(node)
            _add_side_branch(node, set(backbone))

        leftovers = sorted(component - placed, key=self._state_sort_key)
        ordered_ids.extend(leftovers)
        return ordered_ids

    def _crossings_for_order(self, order: Iterable[int]) -> Tuple[int, int, Tuple[int, ...]]:
        # Score a 1D ordering by counting connector crossings and total connector length.
        positions = {state_id: idx for idx, state_id in enumerate(order)}
        edge_list = sorted(self.edge_pairs)
        crossings = 0
        edge_length = 0
        for idx, (a1, b1) in enumerate(edge_list):
            x1, x2 = sorted((positions[a1], positions[b1]))
            edge_length += x2 - x1
            for a2, b2 in edge_list[idx + 1 :]:
                if len({a1, b1, a2, b2}) < 4:
                    continue
                y1, y2 = sorted((positions[a2], positions[b2]))
                if (x1 < y1 < x2 < y2) or (y1 < x1 < y2 < x2):
                    crossings += 1
        return crossings, edge_length, tuple(order)

    def _crossings_order(self, component: Set[int]) -> List[int]:
        # For small graphs try all permutations; otherwise improve a good initial guess locally.
        component_list = sorted(component, key=self._state_sort_key)
        if len(component_list) <= 8:
            return list(min(itertools.permutations(component_list), key=lambda perm: self._crossings_for_order(perm)))

        best = min(
            (self._dfs_order(component), self._bfs_order(component), self._longest_path_order(component)),
            key=self._crossings_for_order,
        )
        improved = list(best)
        improved_flag = True
        while improved_flag:
            improved_flag = False
            best_score = self._crossings_for_order(improved)
            for i in range(len(improved) - 1):
                candidate = improved.copy()
                candidate[i], candidate[i + 1] = candidate[i + 1], candidate[i]
                candidate_score = self._crossings_for_order(candidate)
                if candidate_score < best_score:
                    improved = candidate
                    best_score = candidate_score
                    improved_flag = True
        return improved

    def _force_order(self, component: Set[int]) -> List[int]:
        # Relax 1D positions with attractive edges and repulsive nodes, then sort by the relaxed positions.
        component_ids = sorted(component, key=self._state_sort_key)
        initial = self._longest_path_order(component)
        x_pos: Dict[int, float] = {state_id: float(index) for index, state_id in enumerate(initial)}
        ideal_gap = max(self.spacing, 1.0)
        for _ in range(max(self.force_iterations, 1)):
            delta: Dict[int, float] = {state_id: 0.0 for state_id in component_ids}
            for i, left in enumerate(component_ids):
                for right in component_ids[i + 1 :]:
                    distance = x_pos[right] - x_pos[left]
                    if abs(distance) < 1e-8:
                        distance = 1e-8
                    repulsion = 0.02 / abs(distance)
                    delta[left] -= repulsion
                    delta[right] += repulsion
            for left, right in self.edge_pairs:
                if left in component and right in component:
                    distance = x_pos[right] - x_pos[left]
                    attraction = 0.08 * (distance - ideal_gap)
                    delta[left] += attraction
                    delta[right] -= attraction
            for state_id in component_ids:
                state = self.state_map[state_id]
                delta[state_id] += 0.03 * Units.convert(
                    state.energy - min(s.energy for s in self.states), "hartree", self.unit
                )
            for state_id in component_ids:
                x_pos[state_id] += delta[state_id]
            centered = float(np.mean([x_pos[state_id] for state_id in component_ids]))
            for state_id in component_ids:
                x_pos[state_id] -= centered
        return sorted(component_ids, key=lambda state_id: (x_pos[state_id],) + self._state_sort_key(state_id))

    def ordered_ids(self, layout: str) -> List[int]:
        layout_orderers: Dict[str, Callable[[], List[int]]] = {
            "dfs": lambda: self._ordered_components(self._dfs_order),
            "bfs": lambda: self._ordered_components(self._bfs_order),
            "longest_path": lambda: self._ordered_components(self._longest_path_order),
            "force": lambda: self._ordered_components(self._force_order),
            "crossings": lambda: self._ordered_components(self._crossings_order),
        }

        if layout == "auto":
            # Compare all available strategies and keep the one with the cleanest connector pattern.
            candidate_orders: Dict[str, List[int]] = {name: orderer() for name, orderer in layout_orderers.items()}
            ordered_ids: List[int] = min(
                candidate_orders.values(),
                key=lambda order: self._crossings_for_order(order)
                + (sum(abs(i - order.index(state.id)) for i, state in enumerate(self.states)),),
            )
        else:
            if layout not in layout_orderers:
                raise ValueError(
                    "Unsupported layout '{}'. Choose from 'auto', 'dfs', 'bfs', "
                    "'longest_path', 'force', or 'crossings'.".format(layout)
                )
            ordered_ids = layout_orderers[layout]()
        return ordered_ids

    def x_positions(self, ordered_ids: List[int], landscape_width: float) -> Dict[int, float]:
        # Convert the chosen ordering into actual x-coordinates, leaving a gap between components.
        component_gap = max(1.5 * self.spacing, self.spacing + landscape_width)
        x_map: Dict[int, float] = {}
        current_x = 0.0
        for component in self._connected_components():
            component_order = [state_id for state_id in ordered_ids if state_id in component]
            for offset, state_id in enumerate(component_order):
                x_map[state_id] = current_x + offset * self.spacing
            current_x = x_map[component_order[-1]] + component_gap
        return x_map


def _crop_view_image(image: Any) -> Any:
    """Crop empty borders from molecule images returned by ``view()``."""
    image_array = np.asarray(image)
    if image_array.ndim < 2:
        return image

    if image_array.ndim == 3 and image_array.shape[2] == 4:
        mask = image_array[:, :, 3] > 0
    elif image_array.ndim == 3:
        mask = np.any(image_array[:, :, :3] < 250, axis=2)
    else:
        mask = image_array < 250

    nonempty = np.argwhere(mask)
    if nonempty.size == 0:
        return image

    y0, x0 = nonempty.min(axis=0)
    y1, x1 = nonempty.max(axis=0) + 1
    if x0 == 0 and y0 == 0 and y1 == image_array.shape[0] and x1 == image_array.shape[1]:
        return image

    return image.crop((int(x0), int(y0), int(x1), int(y1)))

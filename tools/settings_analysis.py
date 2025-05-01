import ast
import warnings
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

from scm.plams import Settings

__all__ = ["compare_settings", "CategorizeSettings"]


def compare_settings(
    settings_list: List[Settings],
    analyze_blocks: bool = True,
    analyze_keys: bool = False,
    default_settings: Optional[Settings] = None,
    flatten_list: bool = True,
    remove_unimportant_columns: bool = True,
    unimportant_variation_threshold: int = 1,
    none_is_unimportant: bool = False,
) -> Dict[Tuple, List[Any]]:
    """Deprecated use plams.JobAnalysis instead"""
    comparison_summary = defaultdict(list)

    if default_settings is None:
        default_settings = Settings()
    else:
        check_each_flatten = [isinstance(x[0], tuple) for x in list(default_settings.flatten())]
        is_already_flatten = all(check_each_flatten)
        if not is_already_flatten:
            default_settings = default_settings.flatten()

    all_s_not_flatten = Settings()
    settings_list_flatten = []
    for s in settings_list:
        all_s_not_flatten += s
        settings_list_flatten.append(s.flatten(flatten_list=flatten_list))

    all_blocks = all_s_not_flatten.get_blocks_path(flatten_list=flatten_list, add_empty_settings_as_block=True)
    if analyze_blocks:
        for s in settings_list:
            for k in all_blocks:
                value = bool(s.get_nested(k, default=False))
                comparison_summary[k].append(value)

    if analyze_keys:
        for s in settings_list_flatten:
            keys_paths = set(all_s_not_flatten.flatten().keys()) - set(all_blocks)
            for k in keys_paths:
                value = s.get(k, default=default_settings.get(k, default=None))
                try:
                    # trying to convert 'True' to True and other values
                    value = ast.literal_eval(value)
                except (ValueError, SyntaxError):
                    pass
                comparison_summary[k].append(value)

    if remove_unimportant_columns:
        for k in list(comparison_summary.keys()):
            if isinstance(comparison_summary[k][0], list):
                warnings.warn(
                    "comparison_summary[k] is a list is better to set flatten_list=True, for now list converted to tuple"
                )
                comparison_summary[k] = [tuple(x) for x in comparison_summary[k]]
            set_options = set(comparison_summary[k])
            if none_is_unimportant:
                set_options = set_options - set([None])
            if len(set_options) == unimportant_variation_threshold:
                comparison_summary.pop(k)
    return dict(comparison_summary)


def compare_all_subtuples(tuple1, tuple_long, key_sensitive=False):
    if key_sensitive:
        return [tuple1 == tuple_long[i:] for i in range(len(tuple_long) + 1 - len(tuple1))]
    return [str(tuple1).lower() == str(tuple_long[i:]).lower() for i in range(len(tuple_long) + 1 - len(tuple1))]


def is_contained(tuple1, tuple_long):
    return any(compare_all_subtuples(tuple1, tuple_long))


def equal(a, b):
    return a == b


def compare_all_susettings(s_in: Settings, s_large: Settings, comparison_method=equal, flatten_list=True):
    flatten_dict = s_large.copy().flatten(flatten_list=flatten_list).as_dict()
    tagged_dict = s_in.copy().flatten(flatten_list=flatten_list).as_dict()
    matches = {}
    for k_tag, v_tag in tagged_dict.items():
        for k, v in flatten_dict.items():
            if is_contained(k_tag, k):
                matches[k_tag] = comparison_method(v_tag, v)
        if k_tag not in matches:
            matches[k_tag] = False
    return matches


def is_contained_susettings(s: Settings, s_large: Settings, comparison_method=equal, flatten_list=True):
    matches = compare_all_susettings(s, s_large, comparison_method=comparison_method, flatten_list=flatten_list)
    return all(matches.values())


DictLike = Union[Settings, Dict]


@dataclass
class CategorizeSettings:
    collection_categories: Dict[str, DictLike]
    default: Optional[str] = None
    comparison_method: Callable[[Any, Any], bool] = equal
    collision: Literal["raise", "overwrite", "concatenate"] = "concatenate"
    flatten_list: bool = True

    def run(self, s_collection: List[DictLike]):
        ret = []
        for s in s_collection:
            ret.append(self.run_one(s))
        return ret

    def run_one(self, s: DictLike):
        ret = self.default
        for k, set_in in self.collection_categories.items():
            if is_contained_susettings(
                Settings(set_in).copy(),
                Settings(s),
                comparison_method=self.comparison_method,
                flatten_list=self.flatten_list,
            ):
                if ret == self.default:
                    ret = k
                else:
                    ret = self._handle_collision(ret, s, k)
        return ret

    def _handle_collision(self, ret, s, k):
        if self.collision == "raise":
            raise ValueError(
                f"The tag was already assigned with value {ret}, but wanted to assign to {k} for \n{Settings(s)}"
            )
        elif self.collision == "concatenate":
            ret = f"{ret}-{k}"
        elif self.collision == "overwrite":
            ret = k
        else:
            raise ValueError(f"Unknown collision method {self.collision}")
        return ret

import ast
import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from scm.plams import Settings


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

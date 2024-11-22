import warnings
from collections import defaultdict
from typing import List, Optional

import numpy as np
from scm.plams import Settings


def print_in_table(data):
    # Extract keys and values
    keys = list(data.keys())  # Keys will form the table header
    values = np.array(list(data.values()), dtype=object).T  # Transpose values to make rows the data points
    # Prepare table headers
    header = [str(key) for key in keys]
    # Calculate column widths dynamically
    col_widths = [max(len(str(key)), max(len(str(row[i])) for row in values)) for i, key in enumerate(keys)]
    # Create the separator line
    separator = "-".join("-" * (width + 2) for width in col_widths)
    # Print the table
    header_row = " | ".join(f"{key:<{col_widths[i]}}" for i, key in enumerate(header))
    print(header_row)
    print(separator)
    # Print each row of values
    for row in values:
        row_text = " | ".join(f"{str(value):<{col_widths[i]}}" for i, value in enumerate(row))
        print(row_text)


def compare_settings(
    settings_list: List[Settings],
    blocks_analysis: bool = True,
    keys_analysis: bool = False,
    default_settings: Optional[Settings] = None,
    flatten_list: bool = True,
    clean_singular_value: bool = True,
):
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

    if blocks_analysis:
        all_blocks = all_s_not_flatten.get_blocks_path(flatten_list=flatten_list, add_empty_settings_as_block=True)
        for s in settings_list:
            for k in all_blocks:
                value = bool(s.get_nested(k, default=False))
                comparison_summary[k].append(value)

    if keys_analysis:
        for s in settings_list_flatten:
            for k in all_s_not_flatten.flatten().keys():
                value = s.get(k, default=default_settings.get(k, default=None))
                comparison_summary[k].append(value)

    if clean_singular_value:
        for k in list(comparison_summary.keys()):
            if isinstance(comparison_summary[k][0], list):
                warnings.warn(
                    "comparison_summary[k] is a list is better to set flatten_list=True, for now list converted to tuple"
                )
                comparison_summary[k] = [tuple(x) for x in comparison_summary[k]]
            if len(set(comparison_summary[k])) == 1:
                comparison_summary.pop(k)
    return dict(comparison_summary)

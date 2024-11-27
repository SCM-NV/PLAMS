from typing import Dict, List

import numpy as np


def format_in_table(
    data: Dict[str, List],
    ret_str: bool = False,
    print_on: bool = True,
    max_col_length: int = -1,
    max_rows_displayed: int = 30,
):
    """_summary_

    :param data: _description_
    :type data: Dict[str, List]
    :param ret_str: _description_, defaults to False
    :type ret_str: bool, optional
    :param print_on: _description_, defaults to True
    :type print_on: bool, optional
    :param max_col_length: can be integer positive value or -1, defaults to -1
    :type max_col_length: int, optional
    :param max_rows_displayed: can be integer positive value or -1, defaults to 10
    :type max_rows_displayed: int, optional
    :return: _description_
    :rtype: _type_
    """

    # Function to truncate strings
    def truncate(s, max_len):
        s = str(s)
        if max_len > 0 and len(s) > max_len:
            return s[:max_len] + "..."
        else:
            return s

    # Extract keys and values
    keys = list(data.keys())  # Keys will form the table header

    # Truncate headers
    truncated_header = [truncate(str(key), max_col_length) for key in keys]

    # Truncate values
    columns = [[truncate(value, max_col_length) for value in col_values] for col_values in data.values()]

    # Transpose columns to get rows
    values = np.array(columns, dtype=object).T  # Transpose to get rows

    # Calculate column widths dynamically based on truncated headers and values
    col_widths = [max(len(truncated_header[i]), max(len(str(row[i])) for row in values)) for i in range(len(keys))]

    # Prepare the table
    strs_to_print = []
    header_row = " | ".join(f"{truncated_header[i]:<{col_widths[i]}}" for i in range(len(keys)))
    strs_to_print.append(header_row)
    # Create the separator line
    separator = "-+-".join("-" * width for width in col_widths)
    strs_to_print.append(separator)

    num_rows = len(values)
    if num_rows <= max_rows_displayed or max_rows_displayed == -1:
        # Print all rows
        for row in values:
            row_text = " | ".join(f"{str(row[i]):<{col_widths[i]}}" for i in range(len(keys)))
            strs_to_print.append(row_text)
    else:
        # Determine how many rows to show at the start and end
        rows_to_show = max_rows_displayed - 1  # Subtract 1 for the separator row
        rows_at_start = rows_to_show // 2
        rows_at_end = rows_to_show - rows_at_start

        # First part of the rows
        for row in values[:rows_at_start]:
            row_text = " | ".join(f"{str(row[i]):<{col_widths[i]}}" for i in range(len(keys)))
            strs_to_print.append(row_text)

        # Separator row with '...'
        dots_row = " | ".join(f"{'...'[:col_widths[i]]:<{col_widths[i]}}" for i in range(len(keys)))
        strs_to_print.append(dots_row)

        # Last part of the rows
        for row in values[-rows_at_end:]:
            row_text = " | ".join(f"{str(row[i]):<{col_widths[i]}}" for i in range(len(keys)))
            strs_to_print.append(row_text)
    str_to_print = "\n".join(strs_to_print)
    if print_on:
        print(str_to_print)
    if ret_str:
        return str_to_print

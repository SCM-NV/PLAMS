from typing import Dict, List

import numpy as np


__all__ = ["format_in_table"]


def format_in_table(
    data: Dict[str, List],
    max_col_width: int = -1,
    max_rows: int = 30,
) -> str:
    """
    Create a table from a dictionary of data, with the keys as column headers and the values as rows.

    This will appear in a Markdown compatible format, as follows:

    .. code-block:: python

        | A | B    | C        |
        |---|------|----------|
        | 1 | one  | row o... |
        | . | ...  | ...      |
        | 5 | five | row f... |

    :param data: data to format in the table
    :param max_col_width: can be integer positive value or -1, defaults to -1 (no maximum width)
    :param max_rows: can be integer positive value or -1, defaults to 30
    :return: formatted string table
    """
    elip = "..."

    def truncate(s, max_len):
        s = str(s)
        return f"{s[:max_len]}{elip}" if 0 < max_len < len(s) else s

    # Extract keys for headers
    keys = list(data.keys())

    # Truncate headers and values based on maximum column width
    truncated_header = [truncate(str(key), max_col_width) for key in keys]
    columns = [[truncate(value, max_col_width) for value in col_values] for col_values in data.values()]

    # Transpose columns to get rows
    values = np.array(columns, dtype=object).T

    # Calculate column widths dynamically based on truncated headers and values and values of displayed rows
    num_rows = len(values)
    write_all_rows = num_rows <= max_rows or max_rows == -1
    rows_at_start = num_rows if write_all_rows else max_rows // 2
    rows_at_end = 0 if write_all_rows else max_rows - rows_at_start
    width_values = [v for i, v in enumerate(values) if i < rows_at_start or i >= (num_rows - rows_at_end)]
    col_widths = [
        max(len(truncated_header[i]), max(len(str(row[i])) for row in width_values)) for i in range(len(keys))
    ]

    # Create the separator line
    separator = f"|-{'-|-'.join('-' * width for width in col_widths)}-|"

    # Prepare the table
    table_rows = []
    header_row = f"| {' | '.join(f'{truncated_header[i].ljust(col_widths[i])}' for i in range(len(keys)))} |"
    table_rows.append(header_row)
    table_rows.append(separator)

    if write_all_rows:
        # Make all rows
        for row in values:
            row_text = f"| {' | '.join(f'{str(row[i]).ljust(col_widths[i])}' for i in range(len(keys)))} |"
            table_rows.append(row_text)
    else:
        # First part of the rows
        for row in values[:rows_at_start]:
            row_text = f"| {' | '.join(f'{str(row[i]).ljust(col_widths[i])}' for i in range(len(keys)))} |"
            table_rows.append(row_text)

        # Separator row with '...'
        dots_row = f"| {' | '.join(f'{elip[:col_widths[i]].ljust(col_widths[i])}' for i in range(len(keys)))} |"
        table_rows.append(dots_row)

        # Last part of the rows
        for row in values[-rows_at_end:]:
            row_text = f"| {' | '.join(f'{str(row[i]).ljust(col_widths[i])}' for i in range(len(keys)))} |"
            table_rows.append(row_text)

    table = "\n".join(table_rows)

    return table

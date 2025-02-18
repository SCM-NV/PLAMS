from typing import Dict, List, Literal, Union

import numpy as np


__all__ = ["format_in_table"]


def format_in_table(
    data: Dict[str, List],
    max_col_width: int = -1,
    max_rows: int = 30,
    fmt: Union[Literal["markdown", "html", "rst"]] = "markdown",
    monospace=False,
) -> str:
    """
    Create a table from a dictionary of data, with the keys as column headers and the values as rows.

    This will appear either in a Markdown as follows:

    .. code-block:: python

        | A | B    | C        |
        |---|------|----------|
        | 1 | one  | row o... |
        | . | ...  | ...      |
        | 5 | five | row f... |

    In rst as follows:

    .. code-block:: python

        +---+------+----------+
        | A | B    | C        |
        +===+======+==========+
        | 1 | one  | row o... |
        +---+------+----------+
        | . | ...  | ...      |
        +---+------+----------+
        | 5 | five | row f... |
        +---+------+----------+

    or HTML format.

    :param data: data to format in the table
    :param max_col_width: can be integer positive value or -1, defaults to -1 (no maximum width)
    :param max_rows: can be integer positive value or -1, defaults to 30
    :param fmt: format of the table, either markdown (default) or html
    :param monospace: whether to use monospace text, defaults to False
    :return: formatted string table
    """
    elip = "..."

    def truncate(s, max_len):
        s = str(s)
        return f"{s[:max_len]}{elip}" if 0 < max_len < len(s) else s

    if fmt not in ["markdown", "html", "rst"]:
        raise ValueError(f"'fmt' was '{fmt}' but must be one of: 'markdown', 'html'")

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
        max(len(truncated_header[i]), max([len(str(row[i])) for row in width_values], default=0))
        for i in range(len(keys))
    ]

    # Prepare the table
    table_rows = []
    if fmt == "html":

        def make_html_row(row=None, skip=False):
            return f"<tr>{''.join(f'<td>{(str(row[i]) if not skip else elip[:col_widths[i]]).ljust(col_widths[i])}</td>' for i in range(len(keys)))}</tr>"

        table_rows = [
            '<div style="max-width: 100%; overflow-x: auto;">',
            f'<table border="1" style="border-collapse: collapse; width: auto; {("font-family: monospace; " if monospace else "")}">',
        ]
        header_row = f"<thead><tr>{''.join(f'<th>{truncated_header[i].ljust(col_widths[i])}' for i in range(len(keys)))}</th></tr></thead>"

        separator = "<tbody>"

        make_row = make_html_row
    elif fmt == "rst":

        def make_rst_row(row=None, skip=False):
            row = f"| {' | '.join(f'{(str(row[i]) if not skip else elip[:col_widths[i]]).ljust(col_widths[i])}' for i in range(len(keys)))} |"
            sep = f"+{'+'.join(['-' * (w + 2) for w in col_widths])}+"
            return f"{row}\n{sep}"

        header_row = f"+{'+'.join(['-' * (w + 2) for w in col_widths])}+\n"
        header_row += (
            f"| {' | '.join(f'{(str(truncated_header[i])).ljust(col_widths[i])}' for i in range(len(keys)))} |"
        )

        separator = f"+{'+'.join(['=' * (w + 2) for w in col_widths])}+"

        make_row = make_rst_row
    else:

        def make_md_row(row=None, skip=False):
            return f"| {' | '.join(f'{(str(row[i]) if not skip else elip[:col_widths[i]]).ljust(col_widths[i])}' for i in range(len(keys)))} |"

        if monospace:
            table_rows.append("<pre>")

        separator = f"|-{'-|-'.join('-' * width for width in col_widths)}-|"
        header_row = make_md_row(truncated_header)

        make_row = make_md_row
    table_rows.append(header_row)
    table_rows.append(separator)

    if write_all_rows:
        # Make all rows
        for row in values:
            row_text = make_row(row)
            table_rows.append(row_text)
    else:
        # First part of the rows
        for row in values[:rows_at_start]:
            row_text = make_row(row)
            table_rows.append(row_text)

        # Separator row with '...'
        dots_row = make_row(skip=True)
        table_rows.append(dots_row)

        # Last part of the rows
        for row in values[-rows_at_end:]:
            row_text = make_row(row)
            table_rows.append(row_text)

    if fmt == "html":
        table_rows.append("</tbody>\n</table>\n</div>")
    elif fmt == "rst":
        pass
    else:
        if monospace:
            table_rows.append("</pre>")
    table = "\n".join(table_rows)

    return table

import pytest

from scm.plams.tools.table_formatter import format_in_table


class TestFormatInTable:

    @pytest.fixture
    def data(self):
        return {
            "A": [1, 22, -3, 44, -55],
            "B": ["one", "two", "three", "four", "five"],
            "CCCCC": ["max", "col", "width", "is", "five"],
            "D": ["evenbigger", "than", "maximum", "column", "width"],
            "EEEEEEEE": ["header", "also", "evenbiggerrrrrr", "than", "max"],
        }

    def test_format_in_table_default(self, data):
        t = format_in_table(data)
        assert (
            t
            == """\
| A   | B     | CCCCC | D          | EEEEEEEE        |
|-----|-------|-------|------------|-----------------|
| 1   | one   | max   | evenbigger | header          |
| 22  | two   | col   | than       | also            |
| -3  | three | width | maximum    | evenbiggerrrrrr |
| 44  | four  | is    | column     | than            |
| -55 | five  | five  | width      | max             |"""
        )

    def test_format_in_table_with_max_column_width_and_max_rows(self, data):
        t = format_in_table(data, max_col_width=5, max_rows=3)
        assert (
            t
            == """\
| A   | B    | CCCCC | D        | EEEEE... |
|-----|------|-------|----------|----------|
| 1   | one  | max   | evenb... | heade... |
| ... | ...  | ...   | ...      | ...      |
| 44  | four | is    | colum... | than     |
| -55 | five | five  | width    | max      |"""
        )

        t = format_in_table(data, max_col_width=3, max_rows=2)
        assert (
            t
            == """\
| A   | B      | CCC... | D      | EEE... |
|-----|--------|--------|--------|--------|
| 1   | one    | max    | eve... | hea... |
| ... | ...    | ...    | ...    | ...    |
| -55 | fiv... | fiv... | wid... | max    |"""
        )

        t = format_in_table(data, max_col_width=100, max_rows=1)
        assert (
            t
            == """\
| A   | B    | CCCCC | D     | EEEEEEEE |
|-----|------|-------|-------|----------|
| ... | ...  | ...   | ...   | ...      |
| -55 | five | five  | width | max      |"""
        )

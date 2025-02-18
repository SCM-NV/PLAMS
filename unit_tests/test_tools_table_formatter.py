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

    def test_format_in_table_default_md(self, data):
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

    def test_format_in_table_with_max_column_width_and_max_rows_md(self, data):
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

    def test_format_in_table_default_html(self, data):
        t = format_in_table(data, fmt="html")
        assert (
            t
            == """\
<div style="max-width: 100%; overflow-x: auto;">
<table border="1" style="border-collapse: collapse; width: auto; ">
<thead><tr><th>A  <th>B    <th>CCCCC<th>D         <th>EEEEEEEE       </th></tr></thead>
<tbody>
<tr><td>1  </td><td>one  </td><td>max  </td><td>evenbigger</td><td>header         </td></tr>
<tr><td>22 </td><td>two  </td><td>col  </td><td>than      </td><td>also           </td></tr>
<tr><td>-3 </td><td>three</td><td>width</td><td>maximum   </td><td>evenbiggerrrrrr</td></tr>
<tr><td>44 </td><td>four </td><td>is   </td><td>column    </td><td>than           </td></tr>
<tr><td>-55</td><td>five </td><td>five </td><td>width     </td><td>max            </td></tr>
</tbody>
</table>
</div>"""
        )

    def test_format_in_table_with_max_column_width_and_max_rows_html(self, data):
        t = format_in_table(data, max_col_width=5, max_rows=3, fmt="html")
        assert (
            t
            == """\
<div style="max-width: 100%; overflow-x: auto;">
<table border="1" style="border-collapse: collapse; width: auto; ">
<thead><tr><th>A  <th>B   <th>CCCCC<th>D       <th>EEEEE...</th></tr></thead>
<tbody>
<tr><td>1  </td><td>one </td><td>max  </td><td>evenb...</td><td>heade...</td></tr>
<tr><td>...</td><td>... </td><td>...  </td><td>...     </td><td>...     </td></tr>
<tr><td>44 </td><td>four</td><td>is   </td><td>colum...</td><td>than    </td></tr>
<tr><td>-55</td><td>five</td><td>five </td><td>width   </td><td>max     </td></tr>
</tbody>
</table>
</div>"""
        )

        t = format_in_table(data, max_col_width=3, max_rows=2, fmt="html")
        assert (
            t
            == """\
<div style="max-width: 100%; overflow-x: auto;">
<table border="1" style="border-collapse: collapse; width: auto; ">
<thead><tr><th>A  <th>B     <th>CCC...<th>D     <th>EEE...</th></tr></thead>
<tbody>
<tr><td>1  </td><td>one   </td><td>max   </td><td>eve...</td><td>hea...</td></tr>
<tr><td>...</td><td>...   </td><td>...   </td><td>...   </td><td>...   </td></tr>
<tr><td>-55</td><td>fiv...</td><td>fiv...</td><td>wid...</td><td>max   </td></tr>
</tbody>
</table>
</div>"""
        )

        t = format_in_table(data, max_col_width=100, max_rows=1, fmt="html")
        assert (
            t
            == """\
<div style="max-width: 100%; overflow-x: auto;">
<table border="1" style="border-collapse: collapse; width: auto; ">
<thead><tr><th>A  <th>B   <th>CCCCC<th>D    <th>EEEEEEEE</th></tr></thead>
<tbody>
<tr><td>...</td><td>... </td><td>...  </td><td>...  </td><td>...     </td></tr>
<tr><td>-55</td><td>five</td><td>five </td><td>width</td><td>max     </td></tr>
</tbody>
</table>
</div>"""
        )

    def test_format_in_table_default_rst(self, data):
        t = format_in_table(data, fmt="rst")
        assert (
            t
            == """\
+-----+-------+-------+------------+-----------------+
| A   | B     | CCCCC | D          | EEEEEEEE        |
+=====+=======+=======+============+=================+
| 1   | one   | max   | evenbigger | header          |
+-----+-------+-------+------------+-----------------+
| 22  | two   | col   | than       | also            |
+-----+-------+-------+------------+-----------------+
| -3  | three | width | maximum    | evenbiggerrrrrr |
+-----+-------+-------+------------+-----------------+
| 44  | four  | is    | column     | than            |
+-----+-------+-------+------------+-----------------+
| -55 | five  | five  | width      | max             |
+-----+-------+-------+------------+-----------------+"""
        )

    def test_format_in_table_with_max_column_width_and_max_rows_rst(self, data):
        t = format_in_table(data, max_col_width=5, max_rows=3, fmt="rst")
        assert (
            t
            == """\
+-----+------+-------+----------+----------+
| A   | B    | CCCCC | D        | EEEEE... |
+=====+======+=======+==========+==========+
| 1   | one  | max   | evenb... | heade... |
+-----+------+-------+----------+----------+
| ... | ...  | ...   | ...      | ...      |
+-----+------+-------+----------+----------+
| 44  | four | is    | colum... | than     |
+-----+------+-------+----------+----------+
| -55 | five | five  | width    | max      |
+-----+------+-------+----------+----------+"""
        )

        t = format_in_table(data, max_col_width=3, max_rows=2, fmt="rst")
        assert (
            t
            == """\
+-----+--------+--------+--------+--------+
| A   | B      | CCC... | D      | EEE... |
+=====+========+========+========+========+
| 1   | one    | max    | eve... | hea... |
+-----+--------+--------+--------+--------+
| ... | ...    | ...    | ...    | ...    |
+-----+--------+--------+--------+--------+
| -55 | fiv... | fiv... | wid... | max    |
+-----+--------+--------+--------+--------+"""
        )

        t = format_in_table(data, max_col_width=100, max_rows=1, fmt="rst")
        assert (
            t
            == """\
+-----+------+-------+-------+----------+
| A   | B    | CCCCC | D     | EEEEEEEE |
+=====+======+=======+=======+==========+
| ... | ...  | ...   | ...   | ...      |
+-----+------+-------+-------+----------+
| -55 | five | five  | width | max      |
+-----+------+-------+-------+----------+"""
        )

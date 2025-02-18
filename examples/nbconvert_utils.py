from nbconvert.preprocessors import Preprocessor
import re
from scm.plams import JobStatus


class PlamsPreprocessor(Preprocessor):
    """
    Preprocessor to remove empty cells and truncate long plams output cells in Jupyter notebooks.
    """

    max_log_lines = 10
    log_line_pattern = f"(.*JOB.*({'|'.join(s.upper() for s in JobStatus)}).*|.*Renaming job.*)"

    def preprocess(self, nb, resources):
        """Process the notebook to remove empty cells and truncate long outputs."""
        new_cells = []
        for cell in nb.cells:
            if self._is_empty(cell):
                continue

            if cell.cell_type == "code":
                self._truncate_outputs(cell)

            new_cells.append(cell)

        nb.cells = new_cells
        return nb, resources

    def _is_empty(self, cell):
        """Check if the cell is completely empty (no source, no output, no metadata)."""
        return not cell.source.strip() and cell.cell_type == "code" and not cell.outputs

    def _truncate_outputs(self, cell):
        """Truncate long log lines."""
        for output in cell.outputs:
            if "text" in output and isinstance(output["text"], str):
                lines = output["text"].splitlines()
                if len(lines) > self.max_log_lines:
                    log_line_count = 0
                    filtered_lines = []
                    for line in lines:
                        if re.fullmatch(self.log_line_pattern, line):
                            log_line_count += 1
                            if log_line_count <= self.max_log_lines:
                                filtered_lines.append(line)
                            elif log_line_count == self.max_log_lines + 1:
                                filtered_lines.append("... (PLAMS log lines truncated) ...")
                        else:
                            filtered_lines.append(line)

                    output["text"] = "\n".join(filtered_lines)

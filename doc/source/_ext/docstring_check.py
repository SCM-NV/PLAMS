from __future__ import annotations
from typing import Any, Callable, Iterable, List, Optional, Tuple
from sphinx.application import Sphinx
from sphinx.util import logging
import inspect

logger = logging.getLogger(__name__)

# -----------------------
# Example check functions
# -----------------------

def get_params(what, obj):
    if what not in {"function", "method"}:
        return None
    try:
        sig = inspect.signature(obj)
    except Exception:
        return None

    params = [p for p in sig.parameters.values()
              if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)
              and p.name not in {"self", "cls"}]
    return params

def check_params_format(what: str, name: str, obj: Any, lines: List[str]) -> Optional[str]:
    """If the object has user params, check if :param: sections exist for each."""
    params = get_params(what, obj)
    if not params:
        return None
    docstring = '\n'.join(lines)
    def param_missing_in_doc(param):
        return f':param {param.name}:' not in docstring

    missing_params = [param for param in params if param_missing_in_doc(param)]

    if not missing_params:
        return None
    return f'missing formatted documentation for parameters {[param.name for param in missing_params]}'

def check_params_mentioned(what: str, name: str, obj: Any, lines: List[str]) -> Optional[str]:
    """If the object has user params, check for each parameter if they are mentioned."""
    params = get_params(what, obj)
    if not params:
        return None
    docstring = '\n'.join(lines)
    def param_missing_in_doc(param):
        return param.name not in docstring

    missing_params = [param for param in params if param_missing_in_doc(param)]

    if not missing_params:
        return None
    return f'The parameters {[param.name for param in missing_params]} are not mentioned in the docstring'



# Register the checks you want to run (order matters)
CHECKS: Tuple[Callable[[str, str, Any, List[str]], Optional[str]], ...] = (
    check_params_mentioned,
    check_params_format,
)

# ---------------------------------
# Sphinx integration / event hookup
# ---------------------------------

def _format_problems(problems: Iterable[str]) -> str:
    return "; ".join(problems)

def _process_docstring(app: Sphinx, what: str, name: str, obj: Any, options: Any, lines: List[str]) -> None:
    problems: list[str] = []
    for check in CHECKS:
        try:
            msg = check(what, name, obj, lines)
        except Exception as exc:
            msg = f"docstring check crashed: {exc!r}"
        if msg:
            problems.append(msg)
    if problems:
        # location=(docname, lineno) is optional; name is a good fallback
        logger.info(f"[docstring-check] {what} '{name}': {_format_problems(problems)}", location=(name, 1))

def setup(app: Sphinx):
    app.connect("autodoc-process-docstring", _process_docstring)
    return {"version": "1.0", "parallel_read_safe": True}

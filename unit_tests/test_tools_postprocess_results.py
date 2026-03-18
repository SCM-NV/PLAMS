import numpy as np
import pytest

from scm.plams.tools.postprocess_results import broaden_results


@pytest.mark.parametrize("broadening_type", ["gaussian_height", "lorentzian_height", "gaussian"])
def test_broaden_results_height_methods_have_peak_height_at_least_max_input_area(broadening_type):
    centers = np.array([1.0, 3.0])
    areas = np.array([1.5, 4.0])

    _, y_result = broaden_results(
        centers=centers,
        areas=areas,
        broadening_width=0.01,
        broadening_type=broadening_type,
    )

    assert np.max(y_result) >= np.max(areas)


@pytest.mark.parametrize("broadening_type", ["gaussian_area", "lorentzian_area", "lorentzian"])
def test_broaden_results_area_methods_preserve_total_area(broadening_type):
    centers = np.array([10.0, 30.0])
    areas = np.array([1.5, 4.0])

    x_result, y_result = broaden_results(
        centers=centers,
        areas=areas,
        broadening_width=1.0,
        broadening_type=broadening_type,
    )
    if hasattr(np, "trapezoid"):
        # avoid warnings
        integrated_area = np.trapezoid(y_result, x_result)
    else:
        # retrocompatibility
        integrated_area = np.trapz(y_result, x_result)
    assert integrated_area <= np.sum(areas)

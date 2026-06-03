import numpy as np
import pytest

from scm.plams.tools.postprocess_results import moving_average, broaden_results


def test_moving_average():
    avg_x, avg_y = moving_average([1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0], window=2)

    assert avg_x.tolist() == pytest.approx([1.5, 2.5, 3.5])
    assert avg_y.tolist() == pytest.approx([3.0, 5.0, 7.0])

    avg_x, avg_y = moving_average([1.0, 2.0], [3.0, 5.0], window=0)

    assert avg_x.tolist() == pytest.approx([1.0, 2.0])
    assert avg_y.tolist() == pytest.approx([3.0, 5.0])


@pytest.mark.parametrize("broadening_type", ["gaussian_height", "lorentzian_height", "gaussian"])
def test_broaden_results_height_methods_have_peak_height_at_least_max_input_area(
    broadening_type,
):
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

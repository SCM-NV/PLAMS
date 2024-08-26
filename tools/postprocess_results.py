import numpy as np
from typing import Union, Literal, Tuple


def _gaussian(x, A, x0, sigma):
    return A * np.exp(-((x - x0) ** 2) / (2 * sigma**2))


def _lorentzian(x, A, x0, sigma):
    return (A / np.pi) * (0.5 * sigma) / ((x - x0) ** 2 + (0.5 * sigma) ** 2)


def _generate_broadening_widths(x_data, broadening_width: Union[float, np.ndarray], centers):
    if not isinstance(broadening_width, np.ndarray):
        sigmas = x_data * 0 + broadening_width
    else:
        if len(broadening_width) == len(centers):
            sigmas = broadening_width
        else:
            raise ValueError(f"{len(broadening_width)=} != {len(centers)=} they should be equal")
    return sigmas


def broaden_results(
    centers: np.ndarray,
    areas: np.ndarray,
    broadening_width: Union[float, np.ndarray] = 40,
    broadening_type: Literal["gaussian", "lorentzian"] = "gaussian",
    x_data: Union[np.ndarray, Tuple[float, float, float]] = (0, 4000, 0.5),
    post_process: Literal["max_to_1"] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """convenient function for create xy curve with gaussian or lorentzian peaks. Used to obtain IR or Raman spectra for example.

    :param centers: x positions of the centers of the peaks
    :type centers: np.ndarray
    :param areas: areas of the peaks
    :type areas: np.ndarray
    :param broadening_width: if is np.ndarray each peak broadening is assigned otherwise apply the same value for all the peaks, defaults to 40
    :type broadening_width: Union[float, np.ndarray], optional
    :param broadening_type: the line shape of the peak, defaults to "gaussian"
    :type broadening_type: Literal['gaussian', 'lorentzian'], optional
    :param x_data: it can be a np.ndarray or you can set a tuple with (min,max,step_size), defaults to (0, 4000, 0.5)
    :type x_data: Union[np.ndarray, Tuple[float, float, float]], optional
    :param post_process: if max_to_1 the resulted spectrum has the max peak hight =1, defaults to None
    :type post_process: Literal['max_to_1'], optional
    :return: the x and y arrays
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    Function_Broaden = {
        "gaussian": _gaussian,
        "lorentzian": _lorentzian,
    }
    if isinstance(x_data, tuple) and len(x_data) == 3:
        x_data = np.arange(*x_data)

    y_data: np.ndarray = x_data * 0
    sigmas = _generate_broadening_widths(x_data, broadening_width, centers)

    for freq_i, Abs_i, sigma_i in zip(centers, areas, sigmas):
        y_data += Function_Broaden[broadening_type](x_data, A=Abs_i, x0=freq_i, sigma=sigma_i)

    if post_process == "max_to_1":
        y_data /= np.max(y_data)
    return x_data, y_data

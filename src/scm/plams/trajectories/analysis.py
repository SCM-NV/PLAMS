import numpy as np

__all__ = ["autocorrelation", "power_spectrum"]


def autocorrelation(np_data, max_dt=None, normalize=False):
    """
    Calculate an autocorrelation function

    np_data : numpy array of shape (n_time_steps, n_data_per_time_step)
        Data for which to calculate the autocorrelation function
    max_dt : int
        Smallest integer delta_time for which the correlation function is not calculated. If max_dt is greater than n_time_steps, n_time_steps is used instead. If None, the largest possible value is used.
    normalize : bool
        Whether to return a normalized autocorrelation function (where the first value is 1)

    Returns : A numpy array of shape (max_dt,) or (max_dt, n_data_per_time_step) (depending on ``average``).
    """
    num_timesteps = np_data.shape[0]
    if max_dt is None or max_dt > num_timesteps:
        max_dt = num_timesteps
    C = np.zeros(max_dt)
    C[0] = 1
    for dt in range(max_dt):
        try:
            Va = np_data[: num_timesteps - dt, :]  # the first num_timesteps-dt rows
            Vb = np_data[dt:, :]  # ignore the first dt rows
            C[dt] = np.mean(Va * Vb)
        except KeyboardInterrupt:
            print("Stopping at dt = {} of max_dt = {}".format(dt, max_dt))
            break
    if normalize:
        C = C / C[0]
    return C


def power_spectrum(times, C, max_freq=None, number_of_points=None):
    """
    Calculate power spectrum from an autocorrelation function
    times : 1D numpy array
        times (in femtoseconds), same length as C. Only the delta between the first two points is needed.
    C : numpy array
        Contains the autocorrelation function
    max_freq: float
        the maximum frequency (in cm-1) to be output
    number_of_points: int
        number of points to include in the FFT.

    Returns: a 2-tuple
        an (x, y) tuple where x = frequency in cm-1, and y = power spectrum as 1D numpy arrays only the real part of the power spectrum is returned
    """
    speedoflight = 29979245800
    # in cm/s
    dt = (times[1] - times[0]) * 1e-15  # in seconds
    max_freq = max_freq or 5000
    if number_of_points is None:
        number_of_points = int(
            1 / (dt * speedoflight)
        )  # the resulting spectrum_x will have points spaced by about 1 cm-1
        if number_of_points < C.size:
            number_of_points = C.size
    spectrum_y = np.fft.fft(C, number_of_points)
    spectrum_x = np.arange(number_of_points) / (dt * speedoflight * number_of_points)
    spectrum_y = np.real(spectrum_y[spectrum_x < max_freq])
    spectrum_x = spectrum_x[spectrum_x < max_freq]
    return spectrum_x.ravel(), spectrum_y.ravel()

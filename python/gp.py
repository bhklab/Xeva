from GPy.models import GPRegression
from GPy.kern import RBF
import numpy as np

from scipy.integrate import quad
from scipy.stats import norm


def normalize_data(X, Y):
    """
    Normalizes all growths using normalize_first_day_and_log_transform helper function.

    :param X: Pandas dataframe of times
    :param Y: Pandas dataframe of volumes
    :return:
    """
    X = X[0]
    x = np.asarray([[day] for day in X])
    y = np.asarray([[size for size in Y[replicate]] for replicate in Y])
    y_norm = __normalize_first_day_and_log_transform(y)

    return x, y, y_norm

def __normalize_first_day_and_log_transform(y):
    """
    Normalize by dividing every y element-wise by the first day's median
    and then taking the log.

    :param y:
    :return:
    """
    if y.ndim == 1:
        return np.log(y / np.median(y[0]))
    else:
        return np.log(y / np.median(y.T[0]))


def fit_gaussian_process(x, y_norm, num_replicates):
    """

    :param x: Numpy array of times
    :param y_norm: Numpy array of normalized tumor volumes
    :param num_replicates: i.e. Number of columns in DF
    :return:
    """

    kernel = RBF(input_dim=1, variance=1., lengthscale=10.)

    x = np.tile(x, (num_replicates, 1))
    y = np.resize(y_norm, (y_norm.shape[0] * y_norm.shape[1], 1))

    gp = GPRegression(x, y, kernel)
    gp.optimize_restarts(num_restarts=9, messages=False)

    return gp

def calculate_kl_divergence(control_x, case_x, gp_control, gp_case):
    """
    Calculates the KL divergence between the GPs fit for both the
    batched controls and batched cases.

    :param control_x:
    :param case_x:
    :param gp_control:
    :param gp_case:
    :return: kl_divergence
    """

    kl_divergence = None

    def kl_integrand(t):
        mean_control, var_control = gp_control.predict(np.asarray([[t]]))
        mean_case, var_case = gp_case.predict(np.asarray([[t]]))

        return np.log10(var_case / var_control) + ((var_control + (mean_control - mean_case) ** 2) / (2 * var_case))

    # TODO: Need to replace zero with drug start day - this data needs to come in from Xeva
    if len(control_x) > len(case_x):
        kl_divergence = abs(quad(kl_integrand, 0, max(case_x))[0]
                                 - max(case_x) / 2)[0]
    else:
        kl_divergence = abs(quad(kl_integrand, 0, max(control_x))[0]
                                 - max(control_x) / 2)[0]

    return kl_divergence

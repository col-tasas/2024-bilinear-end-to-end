"""
Functions to compute the error bounds of the bilinear setting
based on the theretical derivations

"""

__author__ = "Nicolas Chatzikiriakos"
__contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
__date__ = "2024/08/23"


import numpy as np


def boundsAffine(T: int, delta: float, n_x: int, stdNoise: float) -> float:
    """
    Bound for id of the affine linear system resulting from the application
    of u = 1 to the bilin system.
    In accordance to the derivation yields confidence delta/n_u
    Inputs:
        - T: Number of timesteps
        - delta: confidence
        - n_x: state space dimension
        - n_u: inout dimension of the bilin system
        - stdNoise: Standard deviation of the noise
    Output:
        - bound: Error Bound
    """

    bound = (
        4
        / 3
        * np.sqrt(10)
        * stdNoise
        * (np.sqrt(2 * T * np.log(2 * 9 ** (n_x) / delta)))
        / (T / 2 - 4 / 3 * np.sqrt(2 * T * np.log(2 * 9 ** (n_x) / delta)))
    )

    return bound


def boundsAffineData(
    Y: np.array, T: int, delta: float, n_x: int, stdNoise: float
) -> float:
    """
    A postiori (data-dependent) bound for id of the affine linear system resulting
    from the application of u = 1 to the bilin system.
    In accordance to the derivation yields confidence delta/n_u
    Inputs:
        - Z: np.array
        - T: Number of timesteps
        - delta: confidence
        - n_x: state space dimension
        - n_u: inout dimension of the bilin system
        - stdNoise: Standard deviation of the noise
    Output:
        - bound: errorBound
    """
    Z = np.ones((n_x + 1, int(T)))
    Z[:n_x, :] = Y
    M = Z @ Z.transpose()
    EVs, _ = np.linalg.eig(M)
    lambda_min = np.min(np.abs(EVs))

    bound = (
        4
        / 3
        * np.sqrt(10)
        * stdNoise
        * (np.sqrt(2 * T * np.log(9 ** (n_x) / delta)))
        / lambda_min
    )

    return bound

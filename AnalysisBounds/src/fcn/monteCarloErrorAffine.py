"""
Function to compute the approximation error
of Bilin System trough Monte Carlo simulations
"""

__author__ = "Nicolas Chatzikiriakos"
__contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
__date__ = "2024/08/23"


import numpy as np
from numpy.linalg import norm as norm

from fcn.affineFcn import affineFcn


def monteCarloError(
    nMonteCarlo: int,
    fcn: affineFcn,
    T_max: int,
    T_min: int,
    stdSampling: float,
    evals: float,
):
    """
    Inputs:
        - nMonteCarlo: int, number of Monte Carlo rollouts
        - sys: biLinDiscrete, bilinear system
        - u: np.array, dim: (n_u, 1), chosen input -> chooses the scenarios
        - T_max: int, maximum rollout time
        - T_min: int, burn-in time
        - stdSampling: float, sampling std
        - evals: float, number of evals of the error along sim
    """

    n_x = fcn.n_x
    data_x = np.zeros((nMonteCarlo, T_max, 2 * n_x))
    data_z = np.zeros((nMonteCarlo, T_max, n_x + 1))
    t_eval = np.zeros(evals)
    errorA = np.zeros((evals, nMonteCarlo))
    errorB = np.zeros((evals, nMonteCarlo))

    stateSampler = np.random.default_rng()

    t_eval = np.arange(T_min, T_max, np.floor((T_max - T_min) / evals).astype(int))[
        : int(evals)
    ]

    for nMC in range(nMonteCarlo):
        # Simulation
        for kk in range(T_max):
            x_sample = stateSampler.normal(0, stdSampling, (n_x, 1))
            x_p = fcn.sample(x_sample)
            data_x[nMC, kk, :n_x] = x_sample.flatten()
            data_x[nMC, kk, n_x:] = x_p.flatten()

            data_z[nMC, kk, :] = np.append(x_sample.flatten(), 1)

        for index_kk, kk in enumerate(t_eval):

            kk = int(kk)

            # Compute closed-form sol of OLS
            theta_hat = np.linalg.inv(
                data_z[nMC, : int(kk), :].transpose() @ data_z[nMC, : int(kk), :]
            ) @ (data_z[nMC, : int(kk), :].transpose() @ data_x[nMC, :kk, n_x:])

            A_hat = theta_hat[:n_x, :n_x].transpose()
            B_hat = theta_hat[n_x:, :].transpose()

            # Compute the error
            errorB[index_kk, nMC] = norm(B_hat - fcn.B, 2)
            errorA[index_kk, nMC] = norm(A_hat - fcn.A, 2)

    return t_eval, data_x, errorA, errorB

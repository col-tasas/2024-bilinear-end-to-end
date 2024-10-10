""" Computes the Monte Carlo error as well as error bounds for random affine function of specified size"""

__author__ = "Nicolas Chatzikiriakos"
__contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
__date__ = "24/08/23"


import numpy as np
from scipy.io import savemat
import os

from fcn.affineFcn import affineFcn
from fcn.bounds import boundsAffine, boundsAffineData
from fcn.monteCarloErrorAffine import monteCarloError


if __name__ == "__main__":
    # Define parameters
    delta = 0.05
    T_max = 250000
    nMonteCarlo = 100  # number of monte carlo simulations
    total_eval = 250

    # For reproducability
    np.random.seed(42)

    # Define Affine Function
    n_x = 2  # 25

    # Random
    scale1 = 1
    scale0 = 1

    B_i = scale1 * (np.random.rand(n_x, n_x) - 0.5)
    B_0 = scale0 * (np.random.rand(n_x, 1))

    # Noise variance
    stdNoise = 0.5

    fcn = affineFcn(B_i, B_0, stdNoise)

    # Sampling
    stdSampling = 1

    # Burn-in Times
    T_burnIn = np.ceil(64 * (3 + np.sqrt(2) * 2) * np.log(2 * 9 ** (n_x) / delta))
    T_burnInData = np.ceil(0.5 * np.log(9 ** (2 * n_x) / delta))

    if T_max <= T_burnIn:
        raise NameError(
            "T_max smaller than T_burnIn, increase T_max to be larger than "
            + str(T_burnIn)
        )

    # Data Matrices
    errorBound = np.zeros(total_eval)
    errorBoundData = np.zeros((nMonteCarlo, total_eval))
    error = np.zeros((total_eval, nMonteCarlo))
    t_eval = np.zeros(total_eval)

    t_eval, data_x, errorA, errorB = monteCarloError(
        nMonteCarlo, fcn, T_max, T_burnIn, stdSampling, total_eval
    )

    # Compute theoretical error bounds
    for index, t in enumerate(t_eval):

        errorBound[index] = 1 / stdSampling * boundsAffine(t, delta, n_x, stdNoise)

        for nMC in range(nMonteCarlo):
            Y = 1 / stdSampling * data_x[nMC, :, :n_x].transpose()
            errorBoundData[nMC, index] = (
                1
                / stdSampling
                * boundsAffineData(Y[:, : int(t)], t, delta, n_x, stdNoise)
            )

    # save data
    if not (os.path.isdir("../data")):
        os.mkdir("../data")

    np.savez(
        "data/dataError"
        + "nx"
        + str(n_x)
        + "Tmax"
        + str(T_max)
        + "nMC"
        + str(nMonteCarlo)
        + "delta"
        + str(int(100 * delta)),
        t_eval=t_eval,
        errorData=errorA,
        errorBound=errorBound,
        errorBoundData=errorBoundData,
    )

    data_dict = {
        "t_eval": t_eval,
        "errorData": errorA,
        "errorBound": errorBound,
        "errorBoundData": errorBoundData,
    }

    # Save to a .mat file
    savemat(
        "data/dataError"
        + "nx"
        + str(n_x)
        + "Tmax"
        + str(T_max)
        + "nMC"
        + str(nMonteCarlo)
        + "delta"
        + str(int(100 * delta))
        + ".mat",
        data_dict,
    )

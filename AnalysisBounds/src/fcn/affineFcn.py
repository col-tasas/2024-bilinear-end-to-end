"""
Builds affine Function
    x+ = Ax + B + w
"""

__author__ = "Nicolas Chatzikiriakos"
__contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
__date__ = "24/08/08"

import numpy as np


class affineFcn:

    def __init__(self, A: np.array, B: np.array, stdNoise) -> None:
        """Init for a affine function with noise
                x^+ = Ax + B + w
        Inputs:
            - A: A matrix
            - B: B vector
            - stdNoise: std of Gaussian noise
        """
        self.A = A
        self.B = B
        self.n_x, _ = np.shape(A)
        self.stdNoise = stdNoise
        self.noiseGen = np.random.default_rng()

    def sample(self, x: np.array) -> np.array:
        """Evaluetes the affine function at given query point
        Input:
            - x: query point
        Output:
            - x_p: function value at query point
        """
        w = self.noiseGen.normal(0, self.stdNoise, (self.n_x, 1))
        x_p = self.A @ x + self.B + w

        return x_p

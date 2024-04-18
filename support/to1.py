"""Find tower shading efficiency"""

import numpy as np

def calc(alpha_s): # Measured in radians
    """Calculate tower shading efficiency"""

    a = 1 - (7 * ((88 / np.tan(alpha_s)) - 100)) / (np.pi * (250) * (450))
    eta = np.where(a > 1, 1, a)
    return eta

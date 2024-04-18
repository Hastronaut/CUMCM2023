"""Find cosine efficiency"""

import numpy as np

def calc(alpha_o, alpha_r, alpha_s, gamma_s): # Measured in radians
    """Calculate cosine efficiency"""
    return np.sin((alpha_s + alpha_r) / 2)
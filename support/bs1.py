"""Find blocking and shading efficiency"""

import numpy as np

def calc(i, alpha_r, alpha_s): # Greek letters measured in radians
    """Calculate bs efficiency (excluding tower shading)"""

    from nc import mirror, my_mirror

    j = my_mirror[i]
    if (j == -1):
        return 1
    hi = mirror[i, 3]
    hj = mirror[j, 3]
    di = mirror[i, 4]
    dj = mirror[j, 4]
    a = np.minimum(di, dj) * ((hi + hj) * 0.5 - np.linalg.norm(mirror[i, 0:3] - mirror[j, 0:3]) * np.cos((alpha_r - alpha_s) * 0.5))
    return np.where(a <= 0, 1, 1 - a / (hi * hj))
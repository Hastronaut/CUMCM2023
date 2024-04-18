"""Find truncation efficiency"""

import numpy as np

def calc(i):
    """Calculate truncation efficiency"""

    from nc import mirror
    import ef1

    bound = 4.65 / 1000 # 4.65mrad
    h1 = 84 + mirror[i, 3] / (2 * np.cos((ef1.alpha_r - ef1.alpha_s) / 2 + bound))
    h2 = 84 - mirror[i, 3] / (2 * np.cos((ef1.alpha_r - ef1.alpha_s) / 2 - bound))
    a = 56 / np.square(h1 - h2)
    eta_tr = np.where(h1 >= 88 and h2 <= 80, a, 1)
    return eta_tr
"""Find DNI"""

import numpy as np

def Dni(G_0, H, alpha_s): # alpha_s measured in radians
    """Calculate DNI"""

    a = 0.4237 - 0.00821 * np.square(6 - H)
    b = 0.5055 + 0.00595 * np.square(6.5 - H)
    c = 0.2711 + 0.01858 * np.square(2.5 - H)
    dni = G_0 * (a + b * np.exp(-c / np.sin(alpha_s)))
    return dni
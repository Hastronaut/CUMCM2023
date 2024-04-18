"""Find the efficiency of a mirror at a certain time"""

alpha_o, alpha_r, omega, delta, alpha_s, gamma_s, eta_sb, eta_cos, eta_at, eta_trunc, eta_ref, eta = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0


import numpy as np

def calc(i, month, hour):
    """Calculate the said efficiency
    
    This function has side effects."""

    import nc, bs1, to1, co, at, tr

    global alpha_o, alpha_r, omega, delta, alpha_s, gamma_s, eta_sb, eta_cos, eta_at, eta_trunc, eta_ref, eta
    
    global eta_bs1, eta_to1


    na = np.array
    alpha_o = np.arctan2(nc.mirror[i, 1], nc.mirror[i, 0])
    x = 3.5 * np.cos(alpha_o)
    y = 3.5 * np.sin(alpha_o)
    alpha_r = nc.Angle(nc.mirror[i, 0:3], na([x, y, nc.mirror[i, 2]]), na([x, y, 84]))
    phi = nc.phi
    omega = nc.Omega(nc.ST[hour])
    delta = nc.Delta(nc.D[month])
    alpha_s = nc.Alpha_s(delta, phi, omega)
    gamma_s = nc.Gamma_s(delta, phi, alpha_s)
    eta_bs1 = bs1.calc(i, alpha_r, alpha_s)
    eta_to1 = to1.calc(alpha_s)
    eta_sb = eta_bs1 * eta_to1
    eta_cos = co.calc(alpha_o, alpha_r, alpha_s, gamma_s)
    eta_at = at.calc(nc.mirror[i][5])
    eta_trunc = tr.calc(i)
    eta_ref = 0.92
    eta = eta_sb * eta_cos * eta_at * eta_trunc * eta_ref
    return eta
"""Find the efficiency of a mirror at a certain time"""

alpha_o, alpha_r, omega, delta, alpha_s, gamma_s, eta_sb, eta_cos, eta_at, eta_trunc, eta_ref, eta = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0


import numpy as np

def calc(i, month, hour):
    """Calculate the said efficiency
    
    This function has side effects."""

    import nc, co, at, tr
    import ps, p2

    global alpha_o, alpha_r, omega, delta, alpha_s, gamma_s, eta_sb, eta_cos, eta_at, eta_trunc, eta_ref, eta
    
    na = np.array
    x, y = ps.pol2cart(p2.particle[i, 0], p2.particle[i, 1])
    alpha_o = np.arctan2(y, x)
    alpha_r = nc.Angle(na([x, y, p2.particle[i, 2]]), na([0, 0, p2.particle[i, 2]]), na([0, 0, 84]))
    phi = nc.phi
    omega = nc.Omega(nc.ST[hour])
    delta = nc.Delta(nc.D[month])
    alpha_s = nc.Alpha_s(delta, phi, omega)
    gamma_s = nc.Gamma_s(delta, phi, alpha_s)
    eta_sb = 0.999
    eta_cos = co.calc(alpha_o, alpha_r, alpha_s, gamma_s)
    eta_at = at.calc(at.D_hr(i))
    eta_trunc = tr.calc(i)
    eta_ref = 0.92
    eta = eta_sb * eta_cos * eta_at * eta_trunc * eta_ref
    return eta
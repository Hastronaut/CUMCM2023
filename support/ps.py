"""Incorporate global and local search"""

import numpy as np
import p2

shuffle_count = 0
shuffle_result = np.array([])
adjust_ret = None

def place(rou, theta, z, length, width):
    """Five parameters are set: coordinates(rou, theta, z), length, and width."""

    newone = np.array([[rou, theta, z, length, width]])
    if (len(p2.particle) <= p2.num):
        p2.particle = np.append(p2.particle, newone, axis=0)
    else:
        p2.particle[p2.num] = newone[0]
    p2.num += 1

def remove():
    """Remove (pop) the last particle"""

    p2.num -= 1

def gridpos(i):
    """Find the grid containing the particle"""

    rou, theta = p2.particle[i, 0:2]
    nt = int(theta // (2 * np.pi / p2.maxt))
    nr = int((rou - 100) // (600 / p2.maxr))
    return (nr, nt)

def nbhd(i):
    """Find the four neighbors of a particle
    
    Return in shape=(4, 2), invalid result set to [-1, -1]"""

    nr, nt = gridpos(i)
    a = np.array([[nr - 1, nt - 1], [nr - 1, nt + 1], [nr + 1, nt - 1], [nr + 1, nt + 1]], dtype=int)
    for i in range(4):
        if (a[i, 0] < 0) or (a[i, 0] >= p2.maxr):
            a[i] = np.array([-1, -1], dtype=int)
        if (a[i, 1] == -1):
            a[i, 1] = p2.maxt - 1
        if (a[i, 1] == p2.maxt):
            a[i, 1] = 0
    return a

def cart2pol(x, y):
    """Cart to pol, return (rho, theta)
    
    This function takes two scalars as input."""

    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return (rho, theta)

def pol2cart(rho, theta):
    """Pol to cart, return (x, y)
    
    This function takes two scalars as input."""

    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return (x, y)

def dist2p(i, j):
    """Calculate the distance between two particles in 3d-space"""

    if (i == -1) or (j == -1):
        return 100
    x1, y1 = pol2cart(p2.particle[i, 0], p2.particle[i, 1])
    x2, y2 = pol2cart(p2.particle[j, 0], p2.particle[j, 1])
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)

def dkp(i):
    """Return True if proper Distance is Kept between Particle i and its neighbors"""

    neii = nbhd(i)
    for i in range(4):
        if (neii[i, 0] == -1):
            continue
        rou, theta = neii[i]
        if (p2.occupied[rou, theta] == -1):
            continue
        j = p2.occupied[rou, theta]
        dij = dist2p(i, j)
        width = np.maximum(p2.particle[i, 4], p2.particle[j, 4])
        if (dij < width + 5):
            return False
    return True

def infield(rou, theta, opgt, oprt):
    """Return True if (rou, theta) is in the heliostat field"""

    x1, y1 = pol2cart(rou, theta)
    x2, y2 = pol2cart(oprt, opgt - np.pi)
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2) < 350

def run(opgt, oprt, opt, opr, opl, opw, opz, maxnum, testround):
    """Run with gamma_t, rou_t ..."""

    p2.maxnum = maxnum
    p2.rou_t = oprt
    preset(opgt, oprt, opt, opr, opl, opw, opz)
    for i in range(testround):
        shuffle()
    print(shuffle_result)

def preset(opgt, oprt, opt, opr, opl, opw, opz):
    """Preset initial conditions"""

    lt = p2.maxt
    lr = p2.maxr
    for i in range(lr):
        for j in range(lt):
            if (p2.otrange[j] <= i):
                continue
            if (not infield(i, j, opgt, oprt)):
                p2.otrange[j] = i
                continue
    for i in range(lr - 1):
        for j in range(lt - 1):
            if (p2.num >= p2.maxnum):
                break
            if (np.minimum(p2.otrange[j], p2.otrange[j + 1]) <= lr + 1):
                continue
            rou = (250 + oprt) / p2.maxr * (i * 2 + 1) / 2
            theta = 2 * np.pi / p2.maxt * (j * 2 + 1) / 2
            place(rou, theta, opz, opl, opw)
            if (not dkp(p2.num - 1)):
                remove()
                continue
            p2.occupied[i, j] = int(p2.num - 1)

def adjust(i):
    """Adjust parameters of a mirror locally
    
    Maximize DNI * eta.
    This function has side effects."""

    from scipy.optimize import minimize

    power = 0

    def target(x):
        """The target function to be minimized"""

        import ef2
        import DNI, nc

        nonlocal i, power

        power = 0
        for month in range(12):
            for hour in range(5):
                place(x[0], x[1], p2.particle[i, 2], p2.particle[i, 3], p2.particle[i, 4])
                eta = ef2.calc(p2.num - 1, month, hour)
                dni = DNI.Dni(nc.G_0, nc.H, ef2.alpha_s)
                remove()
                power += dni * eta / 60
        return -power
    
    def consdkp(x):
        """dkp as a contraint function"""

        nonlocal i

        place(x[0], x[1], p2.particle[i, 2], p2.particle[i, 3], p2.particle[i, 4])
        a = np.where(dkp(p2.num - 1), 1, -1)
        remove()
        return a
    
    global adjust_ret
    
    nr, nt = gridpos(i)
    sr = (250 + p2.rou_t) / p2.maxr
    st = 2 * np.pi / p2.maxt
    bnds = ((nr * sr, (nr + 1) * sr), (nt * st, (nt + 1) * st))
    cons = ({'type':'ineq', 'fun':consdkp})
    res = minimize(fun=target, x0=p2.particle[i, 0:2], method='SLSQP', bounds=bnds, constraints=cons)
    adjust_ret = res
    return power

def shuffle():
    """Shuffle mirrors globally, one at a time
    
    This function has side effects."""

    global adjust_ret, shuffle_count, shuffle_result

    power = 0
    for i in range(p2.num):
        a = adjust(i)
        b = adjust_ret.x
        power += a / p2.num
        p2.particle[i, 0:2] = b
    shuffle_count += 1
    shuffle_result = np.append(shuffle_result, power)
    print(power * p2.num * 35)
    return power
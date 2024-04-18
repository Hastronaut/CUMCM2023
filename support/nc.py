"""Number crunching"""

import numpy as np

G_0 = 1.366 # kW/m2
H = 3 # km
phi = np.radians(39.4)
D = np.array([-59, -28, 0, 31, 61, 92, 122, 153, 184, 214, 245, 275]) # Year: 2023; index starts from 0
ST = np.array([9, 10.5, 12, 13.5, 15])

num = 0
area = 0
mirror = np.ndarray(shape=(0, 6))
my_mirror = np.ndarray(shape=(0))

def Delta(D):
    """Calculate solar declination
    
    This function takes a value from the array D as input."""

    a = np.sin(2 * np.pi * D / 365) * np.sin(np.pi * 23.45 / 180)
    a = np.where(a > 1, 1, a)
    a = np.where(a < -1, -1, a)
    return np.arcsin(a)

def Omega(hour):
    """Calculate solar hour angle
    
    This function takes a value from the array ST as input."""

    return np.pi / 12 * (hour - 12)

def Alpha_s(delta, phi, omega):
    """Calculate solar altitude angle"""

    a = (np.cos(delta) * np.cos(phi) * np.cos(omega)) + (np.sin(delta) * np.sin(phi))
    a = np.where(a > 1, 1, a)
    a = np.where(a < -1, -1, a)
    return np.arcsin(a)

def Gamma_s(delta, phi, alpha_s):
    """Calculate solar azimuth angle"""
    a = (np.sin(delta) - np.sin(alpha_s) * np.sin(phi)) / (np.cos(alpha_s) * np.cos(phi))
    a = np.where(a > 1, 1, a)
    a = np.where(a < -1, -1, a)
    ans = np.arccos(a)
    return ans

def Angle(a, b, c): # A is vertex
    """Calculate angle BAC, return value in radian
    
    a, b, c should be numpy arrays containing coordinates (x, y, z).
    """

    ab = b - a
    ac = c - a
    cosine_angle = np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac))
    cosine_angle = np.where(cosine_angle > 1, 1, cosine_angle)
    cosine_angle = np.where(cosine_angle < -1, -1, cosine_angle)
    return np.arccos(cosine_angle)

def find_my_mirror(mirror):
    """For each mirror find the closest mirror that has a larger d_hr"""

    global my_mirror
    my_mirror = np.ndarray(shape=(num), dtype=int)
    meerror = mirror[:, 5].argsort()
    for i in range(num):
        myne = -1
        my_dist = 100000
        for j in range(i + 1, num):
            new_dist = np.linalg.norm(mirror[meerror[i], 0:3] - mirror[meerror[j], 0:3])
            if (new_dist < my_dist):
                myne = j
                my_dist = new_dist
        my_mirror[meerror[i]] = meerror[myne]
        
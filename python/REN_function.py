import numpy as np

def ren(theta, k=8, lam=1.6, q=2.0):
    """Base REN formula: a normalised radial profile with k lobes.

    r(theta) = exp(-lam * |sin(k * theta)|^q)
    """
    return np.exp(-lam * np.abs(np.sin(k * theta))**q)

def ren_warp(theta, k=8, lam=1.6, q=2.0, theta0=np.pi/2, mu=0.5):
    """REN with a simple angular warp centred at theta0.

    The warp shifts angles toward theta0 without preserving angular
    measure, producing a local deformation of the base profile.
    """
    phi = theta + mu * np.sin(theta - theta0)
    return np.exp(-lam * np.abs(np.sin(k * phi))**q)

def ren_cluster(theta, k=8, lam=1.6, q=2.0,
                theta0=np.pi/2, mu=0.2, kappa=8.0):
    """REN with a proportional (CDF-based) angular clustering around theta0.

    Unlike ren_warp, the redistribution of angular space is derived from a
    normalised density, which preserves the total angular measure.
    """
    dtheta = theta[1] - theta[0]

    bump = np.exp(kappa * np.cos(theta - theta0))
    bump /= bump.mean()

    w = (1.0 - mu) + mu * bump
    F = np.cumsum(w) * dtheta

    phi = 2*np.pi * (F - F.min()) / (F.max() - F.min())

    return np.exp(-lam * np.abs(np.sin(k * phi))**q)

import logging
import numpy as np
from phoebe import dynamics

from scipy.optimize import newton

logger = logging.getLogger("ETVS")

def crossing(b, time, cind1=0, cind2=1, dynamics_method='keplerian', ltte=True, tol=1e-4, maxiter=1000):
    """
    Compute a crossing ~ time of eclipse ~ light curve minimum.

    Args:
        b: (Bundle) bundle
        time: (float) time close to eclipse
        cind1: (int) index of primary
        cind2: (int) index of secondary
        dynamics_method: (string) dynamics method
        ltte: (bool) compute light-time effects
        tol: (float) tolerance in days
        maxiter: (int) maximum number of iterations

    Returns:
        time_ecl: (float) time of eclipse

    """

    return newton(_projected_separation_sq, x0=time, args=(b, cind1, cind2, dynamics_method, ltte), tol=tol, maxiter=maxiter)


def _projected_separation_sq(time, b, cind1, cind2, dynamics_method, ltte=True):
    """
    Projected separation (^2) to minimize.

    """
    times = np.array([time])

    if dynamics_method in ['nbody', 'rebound']:
        ts, xs, ys, zs, vxs, vys, vzs = dynamics.nbody.dynamics_from_bundle(b, times, compute=None, ltte=ltte, return_roche_euler=False)

    elif dynamics_method=='xyz':
        ts, xs, ys, zs, vxs, vys, vzs = dynamics.xyz.dynamics_from_bundle(b, times, compute=None, ltte=ltte, return_roche_euler=False)

    elif dynamics_method=='bs':
        ts, xs, ys, zs, vxs, vys, vzs = dynamics.nbody.dynamics_from_bundle_bs(b, times, compute=None, ltte=ltte, return_roche_euler=False)

    elif dynamics_method=='keplerian':
        ts, xs, ys, zs, vxs, vys, vzs = dynamics.keplerian.dynamics_from_bundle(b, times, compute=None, ltte=ltte, return_euler=False)

    else:
        raise NotImplementedError

    return (xs[cind2][0]-xs[cind1][0])**2 + (ys[cind2][0]-ys[cind1][0])**2



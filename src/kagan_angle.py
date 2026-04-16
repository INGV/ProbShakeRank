"""
This module includes functions taken from the repository
'pyPTF_data_update' by Louise Cordrie:
https://github.com/louisecordrie/pyPTF_data_update

References
----------
1. Cordrie, L., Selva, J., Bernardi, F., Tonini, R., Romano, F., Volpe, M.,
   & Lorito, S. (2025). Dynamic management of uncertainty in rapid tsunami
   forecasting. Communications Earth & Environment, 6(1), 637.

2. Selva, J., Lorito, S., Volpe, M., Romano, F., Tonini, R., Perfetti, P.,
   Bernardi, F., Taroni, M., Scala, A., Babeyko, A., Løvholt, F., Gibbons, S. J.,
   Macías, J., Castro, M. J., González-Vida, J. M., Sánchez-Linares, C.,
   Bayraktar, H. B., Basili, R., Maesano, F. E., Tiberti, M. M., Mele, F.,
   Piatanesi, A., & Amato, A. (2021). Probabilistic tsunami forecasting for
   early warning. Nature Communications, 12, 5677.
"""

import numpy as np
from copy import deepcopy

def get_kagan_angle(strike1, dip1, rake1, strike2, dip2, rake2):
    """Calculate the Kagan angle between two moment tensors defined by strike,dip and
    rake.
    Kagan, Y. "Simplified algorithms for calculating double-couple rotation",
    Geophysical Journal, Volume 171, Issue 1, pp. 411-418.
    Args:
        strike1 (float): strike of slab or moment tensor
        dip1 (float): dip of slab or moment tensor
        rake1 (float): rake of slab or moment tensor
        strike2 (float): strike of slab or moment tensor
        dip2 (float): dip of slab or moment tensor
        rake2 (float): rake of slab or moment tensor
    Returns:
        float: Kagan angle between two moment tensors
    """
    # convert from strike, dip , rake to moment tensor
    tensor1 = plane_to_tensor(strike1, dip1, rake1)
    tensor2 = plane_to_tensor(strike2, dip2, rake2)

    kagan = calc_theta(tensor1, tensor2)

    return kagan

def calc_theta(vm1, vm2):
    """Calculate angle between two moment tensor matrices.
    Args:
        vm1 (ndarray): Moment Tensor matrix (see plane_to_tensor).
        vm2 (ndarray): Moment Tensor matrix (see plane_to_tensor).
    Returns:
        float: Kagan angle (degrees) between input moment tensors.
    """
    # calculate the eigenvectors of either moment tensor
    V1 = calc_eigenvec(vm1)
    V2 = calc_eigenvec(vm2)

    # find angle between rakes
    th = ang_from_R1R2(V1, V2)

    # calculate kagan angle and return
    for j in range(3):
        k = (j + 1) % 3
        V3 = deepcopy(V2)
        V3[:, j] = -V3[:, j]
        V3[:, k] = -V3[:, k]
        x = ang_from_R1R2(V1, V3)
        if x < th:
            th = x
    return th * 180.0 / np.pi


def calc_eigenvec(TM):
    """Calculate eigenvector of moment tensor matrix.
    Args:
        ndarray: moment tensor matrix (see plane_to_tensor)
    Returns:
        ndarray: eigenvector representation of input moment tensor.
    """

    # calculate eigenvector
    V, S = np.linalg.eigh(TM)
    inds = np.argsort(V)
    S = S[:, inds]
    S[:, 2] = np.cross(S[:, 0], S[:, 1])
    return S

def ang_from_R1R2(R1, R2):
    """Calculate angle between two eigenvectors.
    Args:
        R1 (ndarray): eigenvector of first moment tensor
        R2 (ndarray): eigenvector of second moment tensor
    Returns:
        float: angle between eigenvectors
    """

    #    return np.arccos((np.trace(np.dot(R1, R2.transpose())) - 1.) / 2.)
    return np.arccos(np.clip((np.trace(np.dot(R1, R2.transpose())) - 1.0) / 2.0, -1, 1))

def plane_to_tensor(strike, dip, rake, mag=6.0):
    """Convert strike,dip,rake values to moment tensor parameters.
    Args:
        strike (float): Strike from (assumed) first nodal plane (degrees).
        dip (float): Dip from (assumed) first nodal plane (degrees).
        rake (float): Rake from (assumed) first nodal plane (degrees).
        magnitude (float): Magnitude for moment tensor
            (not required if using moment tensor for angular comparisons.)
    Returns:
        nparray: Tensor representation as 3x3 numpy matrix:
            [[mrr, mrt, mrp]
            [mrt, mtt, mtp]
            [mrp, mtp, mpp]]
    """
    # define degree-radian conversions
    d2r = np.pi / 180.0
    # get exponent and moment magnitude
    magpow = mag * 1.5 + 16.1
    mom = np.power(10, magpow)

    # get tensor components
    mrr = mom * np.sin(2 * dip * d2r) * np.sin(rake * d2r)
    mtt = -mom * ((np.sin(dip * d2r) * np.cos(rake * d2r) * np.sin(2 * strike * d2r)) +
                  (np.sin(2 * dip * d2r) * np.sin(rake * d2r) *
                  (np.sin(strike * d2r) * np.sin(strike * d2r))))
    mpp = mom * ((np.sin(dip * d2r) * np.cos(rake * d2r) * np.sin(2 * strike * d2r)) -
                 (np.sin(2 * dip * d2r) * np.sin(rake * d2r) *
                 (np.cos(strike * d2r) * np.cos(strike * d2r))))
    mrt = -mom * ((np.cos(dip * d2r) * np.cos(rake * d2r) * np.cos(strike * d2r)) +
                  (np.cos(2 * dip * d2r) * np.sin(rake * d2r) * np.sin(strike * d2r)))
    mrp = mom * ((np.cos(dip * d2r) * np.cos(rake * d2r) * np.sin(strike * d2r)) -
                 (np.cos(2 * dip * d2r) * np.sin(rake * d2r) * np.cos(strike * d2r)))
    mtp = -mom * ((np.sin(dip * d2r) * np.cos(rake * d2r) * np.cos(2 * strike * d2r)) +
                  (0.5 * np.sin(2 * dip * d2r) * np.sin(rake * d2r) *
                  np.sin(2 * strike * d2r)))

    mt_matrix = np.array([[mrr, mrt, mrp],
                          [mrt, mtt, mtp],
                          [mrp, mtp, mpp]])
    mt_matrix = mt_matrix * 1e-7  # convert from dyne-cm to N-m
    return mt_matrix



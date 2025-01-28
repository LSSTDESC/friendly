import numpy as np


def ab2cov(infos: dict):
    """
    Computes the covariance matrix and its inverse for a two dimensional gaussian distribution,
    characterize by its semi-major axis a, semi-minor axis b, and rotation angle theta.

    Args:
        infos (dict): x, y, a, b, theta information for one specific ellipse

    Returns:
        ndarray, ndarray: Covariance matrix and inverse covariance matrix.
    """

    x, y, a, b, theta = infos
    
    sin = np.sin(np.radians(theta))
    cos = np.cos(np.radians(theta))

    #covariance matrix
    C11 = (a*cos)**2 + (b*sin)**2
    C22 = (a*sin)**2 + (b*cos)**2
    C12 = (a**2 - b**2)*sin*cos
    C21 = C12

    cov = np.array( [[C11, C12], [C21, C22]]) / (2*np.log(2))
    
    #inverse covariance matrix
    D11 = (a*sin)**2 + (b*cos)**2
    D22 = (a*cos)**2 + (b*sin)**2
    D12 = (b**2 - a**2)*sin*cos
    D21 = D12

    cov_inv = 2*np.log(2)* np.array( [[D11, D12], [D21, D22]]) / (a*b)**2
    
    return cov, cov_inv  


def gaussian_square_int(infos: dict, F: float=1) -> float:
    """
    Calculates the integral of a squared Gaussian profile specified by the x, y, a, b, theta values from
    the parameter infos. Represents the total brightness of an ellipse-like structure.

    Args:
        infos (dict): x, y, a, b, theta information for one specific ellipse
        F (float, optional): Scaling parameter for adjusting the amplitude of the Gaussian. Defaults to 1.

    Returns:
        float: Integral of the squared Gaussian profile.
    """

    x, y, a, b, theta = infos
    cov, cov_inv = ab2cov(infos)
    
    res = F / (4*np.pi*np.sqrt(np.linalg.det(cov)))
    
    return res


def gaussian_overlap(infos1: dict, infos2: dict, F1: float=1, F2: float=1) -> float:
    """
    Calculates the overlap between two 2D Gaussian profiles, specified by their x, y, a, b, theta information. 
    The function computes the overlap of these two Gaussians and returns a value that quantifies the degree 
    to which the two Gaussians intersect.

    Args:
        infos1 (dict): x1, y1, a1, b1, theta1 information for ellipse 1
        infos2 (dict): x2, y2, a2, b2, theta2 information for ellipse 2
        F1 (float, optional): Scaling parameter for adjusting the amplitude of the Gaussian 1. Defaults to 1.
        F2 (float, optional): Scaling parameter for adjusting the amplitude of the Gaussian 2. Defaults to 1.

    Returns:
        float: overlap between two 2D Gaussian profiles.
    """
    
    x1, y1, a1, b1, theta1 = infos1     
    x2, y2, a2, b2, theta2 = infos2
    
    mu = np.array([x2-x1, y2-y1]) #difference in positions between the two Gaussian
    
    cov1, cov_inv1 = ab2cov(infos1)
    cov2, cov_inv2 = ab2cov(infos2)
    
    det1 = np.linalg.det(cov1)
    det2 = np.linalg.det(cov2)
    det_inv_sum = np.linalg.det(cov_inv1 + cov_inv2)
    
    part1 = F1 * F2 / (2*np.pi * np.sqrt(det1 * det2 * det_inv_sum)) #scaling factor
    part2 = np.exp(-mu @ cov_inv2 @ (cov_inv1 + cov_inv2) @ cov_inv1 @ mu / 2 ) #exponential term

    return part1 * part2 #degree of overlap
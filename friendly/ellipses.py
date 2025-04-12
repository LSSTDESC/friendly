import pandas as pd
import numpy as np
from utils import Group, FCatalog

def ellipse_equation(x: float, y: float, param: dict) -> float:
    """
    Evaluates the equation of an ellipse at a given point (x, y) with specified A, B, C, D, E, F ellipse 
    parameters.

    Args:
        x (float): x-axis position of the ellipse center (in pixels)
        y (float): y-axis position of the ellipse center (in pixels)
        param (dict): A, B, C, D, E, F ellipse parameters.

    Returns:
        float: Value of the ellipse equation at the given point (x, y).
    """
    
    A, B, C, D, E, F = param
    
    return A*x**2 + B*y**2 + 2*C*x*y + 2*D*x + 2*E*y + F


def ellipse_lambda(param: dict) -> float:
    """
    Returns a lambda function to evaluates the equation of an ellipse with specified A, B, C, D, E, F ellipse 
    parameters.

    Args:
        param (dict): A, B, C, D, E, F ellipse parameters.

    Returns:
        lambda: f(x,y) to evaluate the ellipse equation
    """
    
    A, B, C, D, E, F = param

    fn = lambda x,y : A*x**2 + B*y**2 + 2*C*x*y + 2*D*x + 2*E*y + F
    
    return fn

def moments2ab(Ixx: float, Iyy: float, Ixy: float) -> list:
    """
    Convert second moments Ixx, Iyy, Ixy of luminous flux in a, b and theta ellipse parameters.
    This is using a Gaussian flux distribution.

    Args:
        Ixx (float): Ixx second moment of luminous flux
        Iyy (float): Iyy second moment of luminous flux
        Ixy (float): Ixy second moment of luminous flux

    Returns:
        list: list of a, b and theta ellipse parameters.
    """
    theta2 = np.arctan2(2*Ixy, Ixx-Iyy)
    theta2[theta2<0] += 2*np.pi
    theta = theta2 / 2
    
    sin = np.sin(theta)
    cos = np.cos(theta)
    cos2t = cos**2 - sin**2
    
    a_squared= 2*np.log(2)*(Ixx*cos**2 - Iyy*sin**2) / cos2t
    b_squared = 2*np.log(2)*(Iyy*cos**2 - Ixx*sin**2) / cos2t
    
    return [np.sqrt(a_squared), np.sqrt(b_squared), np.degrees(theta)]


def ellipse_infos(group: list, truth_cat: dict, obj_cat: dict, dc2_type: str="object") -> dict:
    """
    Recover the pixel position centre (x, y), the major axis (a), minor axis (b) and orientation
    angle (theta) of the galaxy or object ellipses.

    Args:
        group (list): Friends-of-Friends group
        truth_cat (dict): Truth (galaxy) catalog
        obj_cat (dict): Object catalog
        dc2_type (str, optional): "galaxy" or "object". Defaults to "object".

    Returns:
        dict: Pandas Dataframe containing the x, y, a, b (in pixels) and theta (in degrees) ellipse parameters.
    """
    keys = ['x', 'y', 'a', 'b', 'theta']
    infos = {k:[] for k in keys}
    
    if dc2_type == 'galaxy':
        idx = group[0]
        psf = PSF(group, obj_cat)
        out = deg2pix(group, truth_cat, obj_cat, psf)
        infos = {k:j for (k,j) in zip(keys, out)}

    if dc2_type == 'object':
        idx = group[1]
        infos['x'] = obj_cat['x'][idx]
        infos['y'] = obj_cat['y'][idx]
        infos['a'], infos['b'], infos['theta'] = moments2ab(obj_cat['Ixx_pixel_i'][idx], 
                                                            obj_cat['Iyy_pixel_i'][idx], obj_cat['Ixy_pixel_i'][idx])

    return pd.DataFrame(infos, index=idx)



def f_ellipse_infos(group: list, cat1: FCatalog, cat2: FCatalog,  naming: dict) -> dict:
    """
    Recover the sky position centre (x, y), the major axis (a), minor axis (b) and orientation
    angle (theta) of the galaxy or object ellipses.

    Args:
        group (list): Friends-of-Friends group
        truth_cat (dict): Truth (galaxy) catalog
        obj_cat (dict): Object catalog
        dc2_type (str, optional): "galaxy" or "object". Defaults to "object".

    Returns:
        dict: Pandas Dataframe containing the x, y, a, b (in arcseconds) and theta (in degrees) ellipse parameters.
    """
    keys = ['X', 'Y', 'A', 'B', 'THETA']
    infos = {k:[] for k in keys}

    for idx1 in group.idx1:


    

    

    
    



def ab2AB(x: float, y: float, a: float, b: float, theta: float, sky=False) -> list:
    """
    Convert the major axis (a), minor axis (b) and orientation angle (theta) into A, B, C, D, E, F ellipse parameters.

    Args:
        x (float): x-axis position of the ellipse center (in pixels)
        y (float): y-axis position of the ellipse center (in pixels)
        a (float): major axis (in pixels)
        b (float): minor axis (in pixels)
        theta (float): orientation angle (in degrees)

    Returns:
        list: A, B, C, D, E, F ellipse parameters.
    """
    sin = np.sin(np.radians(theta))
    cos = np.cos(np.radians(theta))
    
    A = (a*sin)**2 + (b*cos)**2
    B = (a*cos)**2 + (b*sin)**2
    C = 2*(b**2 - a**2)*sin*cos
    if sky:
        C *= -1
    D = -2*A*x - C*y
    E = -C*x - 2*B*y
    F = A*x**2 + C*x*y + B*y**2 - (a*b)**2
    C, D, E = C/2, D/2, E/2
    
    return [A, B, C, D, E, F]

def ab2AB_np(x: float, y: float, a: float, b: float, theta: float, sky=False):
    """
    Convert the major axis (a), minor axis (b) and orientation angle (theta) into A, B, C, D, E, F ellipse parameters.

    Args:
        x (float): x-axis position of the ellipse center (in pixels)
        y (float): y-axis position of the ellipse center (in pixels)
        a (float): major axis (in pixels)
        b (float): minor axis (in pixels)
        theta (float): orientation angle (in degrees)

    Returns:
        numpy array: A, B, C, D, E, F ellipse parameters.
    """
    sin = np.sin(np.radians(theta))
    cos = np.cos(np.radians(theta))
    
    A = (a*sin)**2 + (b*cos)**2
    B = (a*cos)**2 + (b*sin)**2
    C = 2*(b**2 - a**2)*sin*cos
    if sky:
        C *= -1
    D = -2*A*x - C*y
    E = -C*x - 2*B*y
    F = A*x**2 + C*x*y + B*y**2 - (a*b)**2
    C, D, E = C/2, D/2, E/2
    
    return np.array([A, B, C, D, E, F])


def ellipse_parameters(infos: dict) -> dict:
    """
    Return the A, B, C, D, E, F ellipse parameters, based on the x, y, a, b, theta ellipse provided information.

    Args:
        infos (dict): Pandas Dataframe containing the x, y, a, b and theta ellipse parameters

    Returns:
        dict: Pandas Dataframe containing A, B, C, D, E, F ellipse parameters.
    """
    out = pd.DataFrame(ab2AB(infos['x'], infos['y'], infos['a'], infos['b'], infos['theta']))
    out.index = ['A', 'B', 'C', 'D', 'E', 'F']
    
    return out.T


def PSF(group: list, cat: dict) -> float:
    """
    Determine the mean PSF of the tract or the PSF around a specific group.
    To be applied to the galaxy shapes.

    Args:
        group (list): Friends-of-Friends group
        cat (dict): Object catalog.

    Returns:
        float: Measured PSF for each group in degrees.
    """

    idx = group[1]

    if len(idx) == 0:
        psf = np.sqrt(2)*np.mean(cat['psf_fwhm_i'])/(2*np.sqrt(2*np.log(2)))/3600

    else:
        psf = np.sqrt(2)*cat['psf_fwhm_i'][idx[0]]/(2*np.sqrt(2*np.log(2)))/3600
        
    return psf


def deg2pix(group, truth_cat, obj_cat, psf):
    """
    Converts major axis (a) and minor axis (b) in degrees into pixels and convolves them with the psf.
    Converts the oriantation angle (theta) from the degree plane into pixel plane.
    Computes the center coordinates (x, y) of galaxies in pixels, using center coordinates of one detected object.

    Args:
        group (list): Friends-of-Friends group
        truth_cat (dict): Truth (galaxy) catalog
        obj_cat (dict): Object catalog
        psf (float): PSF measured around the group

    Returns:
        list: List of [x, y, a, b, theta] parameters.
              x, y, a ,b in pixels
              theta in degrees
    """

    idx1 = group[0]
    idx2 = group[1]

    a_deg = truth_cat['size_true'][idx1]/3600 #in degrees
    b_deg = truth_cat['size_minor_true'][idx1]/3600 #in degrees
    a_psf = np.sqrt(a_deg**2 + psf**2) #in degrees (with the effect of the PSF)
    b_psf = np.sqrt(b_deg**2 + psf**2) #in degrees (with the effect of the PSF)
    
    a = a_psf*3600/0.2 #in pixels (pixel ratio is 0.2"/pix)
    b = b_psf*3600/0.2 #in pixels (pixel ratio is 0.2"/pix)
    theta = 270 - truth_cat['position_angle_true'][idx1]
    theta[theta>180] -= 180

    if len(idx2) == 0:

        x = (((truth_cat['ra'][idx1[0]] - truth_cat['ra'][idx1]) * np.cos(np.radians((truth_cat['dec'][idx1]+truth_cat['dec'][idx1[0]])/2)) * 3600 / 0.2)) 
        y = ((truth_cat['dec'][idx1[0]] - truth_cat['dec'][idx1]) * 3600 / 0.2)

    else:
        
        #pixel displacement with respect to the position of the object #1
        dx = (((obj_cat['ra'][idx2[0]] - truth_cat['ra'][idx1]) * np.cos(np.radians((truth_cat['dec'][idx1]+obj_cat['dec'][idx2[0]])/2)) * 3600 / 0.2)) 
        dy = ((obj_cat['dec'][idx2[0]] - truth_cat['dec'][idx1]) * 3600 / 0.2)
        
        x = obj_cat['x'][idx2[0]] + dx
        y = obj_cat['y'][idx2[0]] - dy

    return [x, y, a, b, theta]


def is_overlapping(p1: dict, p2: dict) -> bool:
    """
    Determines if two ellipses overlap, based on their A, B, C, D, E, F ellipse parameters.
    Returns "True" if they overlap, "False" otherwise.

    Args:
        p1 (dict): Pandas Dataframe containing the A, B, C, D, E, F ellipse parameters of the first ellipse
        p2 (dict): Pandas Dataframe containing the A, B, C, D, E, F ellipse parameters of the second ellipse

    Returns:
        bool: Returns "True" if the two ellipses overlap, "False" otherwise.
    """
    
    if (np.isnan(p1)).any() or (np.isnan(p2)).any():
        return False

    A1, B1, C1, D1, E1, F1 = p1
    A2, B2, C2, D2, E2, F2 = p2
        
    M1 = [[A1, C1, D1], 
          [C1, B1, E1], 
          [D1, E1, F1]]
    
    M2 = [[A2, C2, D2], 
          [C2, B2, E2], 
          [D2, E2, F2]]
    
    l0 = np.linalg.det(M2)
    
    l1 = (np.linalg.det([[A1,C2,D2],
                         [C1,B2,E2],
                         [D1,E2,F2]])
        + np.linalg.det([[A2,C1,D2],
                         [C2,B1,E2],
                         [D2,E1,F2]]) 
        + np.linalg.det([[A2,C2,D1],
                         [C2,B2,E1],
                         [D2,E2,F1]])) / 3
    
    l2 = (np.linalg.det([[A2,C1,D1],
                         [C2,B1,E1],
                         [D2,E1,F1]]) 
         + np.linalg.det([[A1,C2,D1],
                          [C1,B2,E1],
                          [D1,E2,F1]]) 
         + np.linalg.det([[A1,C1,D2],
                          [C1,B1,E2],
                          [D1,E1,F2]])) / 3
    
    l3 = np.linalg.det(M1)
    
    delta1 = np.linalg.det([[l3, l2],
                            [l2, l1]])
    delta2 = np.linalg.det([[l3, l1],
                            [l2, l0]])
    delta3 = np.linalg.det([[l2, l1],
                            [l1, l0]])
    
    discriminant_P = np.linalg.det([[2*delta1, delta2], 
                                    [delta2, 2*delta3]])
    
    if (discriminant_P>=0) and ((l1>0) or (l2>0)):
        return False
    else:
        return True 
    

def overlap_area_MC(param1: dict, param2: dict, xlim: list=[0,1], ylim: list=[0,1], Nx: int=1000, Ny: int=1000):
    """
    Calculates the total area and the area of overlap between two 2D ellipses using Monte Carlo (MC) integration.
    The 2D ellipses are defined by their A, B, C, D, E, F ellipse parameters.

    Args:
        param1 (dict): A1, B1, C1, D1, E1, F1 ellipse parameters of ellipse 1
        param2 (dict): A2, B2, C2, D2, E2, F2 ellipse parameters of ellipse 2
        xlim (list, optional): limits of the x-axis. Defaults to [0,1].
        ylim (list, optional): limits of the y-axis. Defaults to [0,1].
        Nx (int, optional): number of points in the x dimension. Defaults to 1000.
        Ny (int, optional): number of points in the y dimension. Defaults to 1000.

    Returns:
        (float, float): (total area of the two ellipses, overlap area between the two ellipses).
    """
    
    N_total = Nx * Ny
    A_total = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0])
    
    x_grid, y_grid = np.meshgrid(np.linspace(xlim[0], xlim[1], Nx), np.linspace(ylim[0], ylim[1], Ny))
    
    Z1 = ellipse_equation(x_grid, y_grid, param1)
    Z2 = ellipse_equation(x_grid, y_grid, param2)
    
    N_1 = np.sum(Z1 < 0)
    N_2 = np.sum(Z2 < 0)
    
    N_overlap = np.sum((Z1 < 0) & (Z2 < 0))
    
    A_total_ellipses = (N_1 + N_2 - N_overlap) * A_total / N_total
    A_overlap = N_overlap * A_total / N_total
    
    return A_total_ellipses, A_overlap


def plot_shape(x0: float, y0: float, rough_size: float, params: dict, c: str='b', ls: str='-', ax=None, linewidth: float=2):
    """
    Plot an ellipse contour based on the ellipse equation defined by the A,B,C,D,E,F parameters.

    Parameters
    ----------
    x0 : float
        X-coordinate of the center of the ellipse.
    y0 : float
        Y-coordinate of the center of the ellipse.
    rough_size : float
        Approximate half-width of the plotting region.
    params : dict or array-like
        A,B,C,D,E,F parameters defining the ellipse.
    c : str, optional
        Color of the contour line (default is 'b' for blue).
    ls : str, optional
        Line style of the contour (default is '-' for solid line).
    ax : matplotlib.axes.Axes, optional
        Axes object on which to plot the contour. If None, the function assumes an existing plot.
    linewidth : float, optional
        Line width of the contour (default is 2).
    """
    
    x = np.linspace(x0 - rough_size, x0 + rough_size, 100)
    y = np.linspace(y0 - rough_size, y0 + rough_size, 100)
    x, y = np.meshgrid(x,y)
 
    z = ellipse_equation(x, y, params)
    
    ax.contour(x, y, z, [0], colors=c, linestyles=ls, linewidths=linewidth) 

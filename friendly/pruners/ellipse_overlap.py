from .. import Pruner

import numpy as np

class EllipseOverlap(Pruner):
    def __init__(self, tune_params):
        self.fudge = tune_params.get('fudge', 1.)
    
    
    
    def __precompute_one__(self, cat1, cat2, group):
        assert(["a", "b", "theta", "x0", "y0"] in cat1.columns)
        assert(["a", "b", "theta", "x0", "y0"] in cat2.columns)
        p1 = ab2AB(a, b, theta, x0, y0)
    
    def __call__(self, cat1, cat2, groups)->List[Group]:
        # for now, assuming ra, dec rather than pixel coordinates
        
        for group in groups:
            for i1 in group.idx1:
                abt = moments2ab(*[cat1.get_quantity(i1, k) for k in ['Ixx', 'Iyy', 'Ixy']])
                for i2 in group.idx2:
                    p2 = moments2ab(*[cat2.get_quantity(i2, k) for k in ['Ixx', 'Iyy', 'Ixy']])
                    is_overlapping(p1, p2)
        
        

def ellipse_equation(x, y, params):
    
    A, B, C, D, E, F = params
    
    return A*x**2 + B*y**2 + 2*C*x*y + 2*D*x + 2*E*y + F

def ab2AB(a, b, theta, x0, y0):
    
    ### a, b, in any length units; theta in degrees ###
    sin = np.sin(np.radians(theta))
    cos = np.cos(np.radians(theta))
    
    A = (a*sin)**2 + (b*cos)**2
    B = (a*cos)**2 + (b*sin)**2
    C = 2*(b**2 - a**2)*sin*cos
    D = -2*A*x0 - C*y0
    E = -C*x0 - 2*B*y0
    F = A*x0**2 + C*x0*y0 + B*y0**2 - (a*b)**2
    C, D, E = C/2, D/2, E/2
    
    return [A, B, C, D, E, F]

def is_overlapping(p1, p2):
    
    A1, B1, C1, D1, E1, F1 = p1
    A2, B2, C2, D2, E2, F2 = p2
        
    M1 = [[A1, C1, D1], [C1, B1, E1], [D1, E1, F1]]
    M2 = [[A2, C2, D2], [C2, B2, E2], [D2, E2, F2]]
    
    l0 = np.linalg.det(M2)
    l1 = (np.linalg.det([[A1,C2,D2],[C1,B2,E2],[D1,E2,F2]]) +
            np.linalg.det([[A2,C1,D2],[C2,B1,E2],[D2,E1,F2]]) +
            np.linalg.det([[A2,C2,D1],[C2,B2,E1],[D2,E2,F1]])) / 3
    l2 = (np.linalg.det([[A2,C1,D1],[C2,B1,E1],[D2,E1,F1]]) +
            np.linalg.det([[A1,C2,D1],[C1,B2,E1],[D1,E2,F1]]) +
            np.linalg.det([[A1,C1,D2],[C1,B1,E2],[D1,E1,F2]])) / 3
    l3 = np.linalg.det(M1)

    delta1 = np.linalg.det([[l3, l2],[l2, l1]])
    delta2 = np.linalg.det([[l3, l1],[l2, l0]])
    delta3 = np.linalg.det([[l2, l1],[l1, l0]])
    discriminant_P = np.linalg.det([[2*delta1, delta2], [delta2, 2*delta3]])
    
    if (discriminant_P >=0) and ((l1>0) or (l2>0)):
        return False  ## not over-lapping
    else:
        return True
    
# def moments2AB(Ixx, Iyy, Ixy, x0, y0):
    
#     ### Ixx, Iyy, Ixy in same unit as (x0, y0) ###
   
#     pi_ab = 2 * np.sqrt(np.pi) * np.power(Ixx*Iyy - Ixy**2, 1/4)

#     A = 4*Iyy/pi_ab
#     B = 4*Ixx/pi_ab
#     C = -8*Ixy/pi_ab

#     D = -2*A*x0 - C*y0
#     E = -C*x0 - 2*B*y0

#     F = A*x0**2 + C*x0*y0 + B*y0**2 - (pi_ab/np.pi)**2
#     C, D, E = C/2, D/2, E/2
    
#     return [A, B, C, D, E, F]


def moments2ab(Ixx, Iyy, Ixy):
    
    ### Ixx, Iyy, Ixy in same unit as (x0, y0) ###
   
    pi_ab = 2 * np.sqrt(np.pi) * np.power(Ixx*Iyy - Ixy**2, 1/4)
    
    theta2 = np.arctan2(2*Ixy, Ixx-Iyy)
    theta2[theta2<0] += 2*np.pi
    
    theta = theta2 / 2
    
    sin = np.sin(theta)
    cos = np.cos(theta)
    cos2t = cos**2 - sin**2
    
    a_square = 4*(Ixx*cos**2 - Iyy*sin**2) / (pi_ab * cos2t)
    b_square = 4*(Iyy*cos**2 - Ixx*sin**2) / (pi_ab * cos2t)
    
    return [np.sqrt(a_square), np.sqrt(b_square), np.degrees(theta)]

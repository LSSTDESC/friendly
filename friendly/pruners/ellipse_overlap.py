from .. import Pruner

import numpy as np

class EllipseOverlap(Pruner):
    def __init__(self, tune_params):
        self.fudge = tune_params.get('fudge', 1.)
        self.pixel_scale = tune_params.get('pixel_scale', 0.2)
    
    def __precompute_one__(self, cat1, cat2, group):
        assert(["size_semi_major", "size_semi_minor", "theta_orientation", "ra", "dec"] in cat1.columns)
        assert(["size_semi_major", "size_semi_minor", "theta_orientation", "ra", "dec"] in cat2.columns)
    
    def __call__(self, cat1, cat2, groups)->List[Group]:
        # check for ellipse overlapping in pixel coordinates.
	# put ellipse 1 in the origin of pixel coordinate, and put ellipse 2 according to its relative position wrt. ellipse 1
        
        for group in groups:
            for i1 in group.idx1:
                a1, b1, theta1, ra1, dec1 = [cat1.get_quantity(i1, k) for k in ['semi_major', 'semi_minor', 'theta_orientation', 'ra', 'dec']]
		p1 = ab2AB(a1, b1, theta1, 0, 0)
                for i2 in group.idx2:
                	a2, b2, theta2, ra2, dec2 = [cat2.get_quantity(i2, k) for k in ['semi_major', 'semi_minor', 'theta_orientation', 'ra', 'dec']]
			dx = (ra2 - ra1) * np.cos(np.radians( (dec1+dec2) / 2 )) * 3600 / self.pixel_scale
			dy = (dec2 - dec1) * 3600 / self.pixel_scale
			p2 = ab2AB(a2, b2, theta2, dx, dy)
                    	is_overlapping(p1, p2)
        
        

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
    
def moments2ab(Ixx, Iyy, Ixy):
    
    ### Ixx, Iyy, Ixy in same unit as (x0, y0) ###
   
    theta2 = np.arctan2(2*Ixy, Ixx-Iyy)
    theta2[theta2<0] += 2*np.pi
    
    theta = theta2 / 2
    
    sin = np.sin(theta)
    cos = np.cos(theta)
    cos2t = cos**2 - sin**2
    
    a_square = 2*np.log(2)*(Ixx*cos**2 - Iyy*sin**2) / cos2t
    b_square = 2*np.log(2)*(Iyy*cos**2 - Ixx*sin**2) / cos2t
    
    return [np.sqrt(a_square), np.sqrt(b_square), np.degrees(theta)]


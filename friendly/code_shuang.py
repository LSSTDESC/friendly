import numpy as np
import matplotlib.pyplot as plt

###################################################################################################################################################################
def ellipse_equation(x, y, params): #ellipse equation
    
    A, B, C, D, E, F = params
    
    return A*x**2 + B*y**2 + 2*C*x*y + 2*D*x + 2*E*y + F

###################################################################################################################################################################
def ab2AB(x0, y0, a, b, theta): #converts semi-axis a and b and theta angle in A,B,C,D,E,F parameters
    
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

###################################################################################################################################################################
def moments2ab(Ixx, Iyy, Ixy): #converts second moments in a,b and theta parameters !!! Using a Gaussian flux distribution !!!
    
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

###################################################################################################################################################################
def plot_shape(x0, y0, rough_size, params, c='b', ls='-', ax=None, linewidth=2): #plot function
    
    x = np.linspace(x0 - rough_size, x0 + rough_size, 100)
    y = np.linspace(y0 - rough_size, y0 + rough_size, 100)
    x, y = np.meshgrid(x,y)
 
    z = ellipse_equation(x, y, params)
    
    plt.contour(x, y, z, [0], colors=c, linestyles=ls, linewidths=linewidth) 

###################################################################################################################################################################
def is_overlapping(p1, p2):
    
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
    
    if (discriminant_P >=0) and ((l1>0) or (l2>0)):
        return False  ## not over-lapping
    else:
        return True
    
###################################################################################################################################################################
def ab2cov(p_ab):
    
    x, y, a, b, theta = p_ab
    
    sin = np.sin(np.radians(theta))
    cos = np.cos(np.radians(theta))

    C11 = (a*cos)**2 + (b*sin)**2
    C22 = (a*sin)**2 + (b*cos)**2
    C12 = (a**2 - b**2)*sin*cos
    C21 = C12

    cov = np.array( [[C11, C12], [C21, C22]]) / (2*np.log(2))
    
    ## the inverse matrix
    
    D11 = (a*sin)**2 + (b*cos)**2
    D22 = (a*cos)**2 + (b*sin)**2
    D12 = (b**2 - a**2)*sin*cos
    D21 = D12
    
    cov_inv = 2*np.log(2)* np.array( [[D11, D12], [D21, D22]]) / (a*b)**2
    
    return cov, cov_inv   

###################################################################################################################################################################
def gaussian_overlap(p1_ab, p2_ab, F1=1, F2=1):
    
    ## continuous version of eq (9.2) of https://arxiv.org/abs/2103.02078, but only for i=1, j=[1,2]
    
    x1, y1, a1, b1, theta1 = p1_ab     
    x2, y2, a2, b2, theta2 = p2_ab
    
    mu = np.array([x2-x1, y2-y1])
    
    cov1, cov_inv1 = ab2cov(p1_ab)
    cov2, cov_inv2 = ab2cov(p2_ab)
    
    det1 = np.linalg.det(cov1)
    det2 = np.linalg.det(cov2)
    det_inv_sum = np.linalg.det( cov_inv1 + cov_inv2 )
    
    part1 = F1*F2 / (2*np.pi * np.sqrt( det1 * det2 * det_inv_sum ) )
    part2 = np.exp( - mu@cov_inv2@(cov_inv1+cov_inv2)@cov_inv1@mu / 2 )

    return part1 * part2

###################################################################################################################################################################
def gaussian_square_int(p_ab, F=1):
    
    ## calculate the integral of a sqaured gaussian profile
    ## equivalent to gaussian_overlap(p_ab, p_ab), but faster.
    x, y, a, b, theta = p_ab
    cov, cov_inv = ab2cov(p_ab)
    
    res = F / (4*np.pi*np.sqrt(np.linalg.det(cov)))
    
    return res

###################################################################################################################################################################
def overlap_area_MC(p1_AB, p2_AB, xlim=[0,1], ylim=[0,1], Nx=1000, Ny=1000):
    
    N_total = Nx * Ny
    A_total = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0])
    
    x_grid, y_grid = np.meshgrid(np.linspace(xlim[0], xlim[1], Nx), np.linspace(ylim[0], ylim[1], Ny))
    
    Z1 = ellipse_equation(x_grid, y_grid, p1_AB)
    Z2 = ellipse_equation(x_grid, y_grid, p2_AB)
    
    N_1 = np.sum(Z1<0)
    N_2 = np.sum(Z2<0)
    
    N_overlap = np.sum((Z1<0) & (Z2<0))
    
    A_ellipses = (N_1 + N_2 - N_overlap) * A_total / N_total
    A_overlap = N_overlap * A_total / N_total
    
    return A_ellipses, A_overlap

###################################################################################################################################################################
def quadroot(p1_ab, p2_AB):
    
    ## solve interception in ellipse 1 frame: [x, y, theta] = 0
    x, y, a, b, theta = p1_ab
    A, B, C, D, E, F = p2_AB
    
    e = b/a
    dA = A - B*e**2
    dB = B*b**2 + F
    
    p4 = -4*C**2*e**2 - dA**2
    p3 = -8*C*E*e**2 - 4*dA*D
    p2 = -4*E**2*e**2 + 4*C**2*b**2 - 4*D**2 - 2*dA*dB
    p1 = 8*C*E*b**2 - 4*D*dB
    p0 = 4*E**2*b**2 - dB**2
    
    ## solve for x; accept only real roots
    xroots_prop = np.roots([p4, p3, p2, p1, p0])
    xroots_prop = np.real(xroots_prop[np.isreal(xroots_prop)])

    ## substitute for y
    xroots, yroots = [], []
    
    eps_root = 1e-9
    eps_inter = 1e-3
    
    for x in xroots_prop:
        
        Delta = b**2 - e**2*x**2
        if np.abs(Delta)<eps_root:
            y = 0.0
        else:
            y = np.sqrt(Delta)

        # z1 = ellipse_equation_general(p2_AB, x, y)
        # z2 = ellipse_equation_general(p2_AB, x, -y)

        z1 = ellipse_equation(x, y, p2_AB)
        z2 = ellipse_equation(x, -y, p2_AB)
        
        if (np.abs(z1) <= np.abs(z2)) & (np.abs(z1) < eps_inter):
            xroots.append(x)
            yroots.append(y)
        elif (np.abs(z2) < np.abs(z1)) & (np.abs(z2) < eps_inter):
            xroots.append(x)
            yroots.append(-y)

    ### remove duplicated roots 
    if len(xroots) == 0:
        xroots_unique, yroots_unique = xroots, yroots
        
    else:
        xroots_unique = [xroots[0]]
        yroots_unique = [yroots[0]]

        x_dup = []
        eps_dup = 1e-6

        for x, y in zip(xroots[1:], yroots[1:]):
            dist = (np.array(xroots_unique) - x)**2 + (np.array(yroots_unique) - y)**2
            n_dup = np.sum( (dist<eps_dup) )
            if n_dup == 0:
                xroots_unique.append(x)
                yroots_unique.append(y)
            else:
                x_dup.append(x)
                
    ### if there are 3 roots left:
    ### ( [2 unique roots + 2 double root] -> [2 unique roots + 1 double root] ),
    ### further remove the remaining double root.
    ### (it is a tangential interception which should be ignored).
    
    if len(xroots_unique) == 3:
        for x, y in zip(xroots_unique, yroots_unique):
            if x in x_dup:
                xroots_unique.remove(x)
                yroots_unique.remove(y)
            
    return np.array(xroots_unique), np.array(yroots_unique)
    
###################################################################################################################################################################
def map_onto_ell(p_ab, x_global, y_global):
    
    # map a point (x,y) in global coords onto the local coords defined by ellipse:
    # x axis along a, y axis along b, centered on ellipse origin.
    
    x0, y0, a, b, theta = p_ab
    
    sin = np.sin(np.radians(theta))
    cos = np.cos(np.radians(theta))
    
    x_ell = cos* (x_global-x0) + sin*(y_global-y0)
    y_ell = -sin*(x_global-x0) + cos*(y_global-y0)
    
    return x_ell, y_ell

###################################################################################################################################################################
def map_from_ell(p_ab, x_ell, y_ell):
    
    # map a point (x,y) from the ellipse coordinates back to the global coordinates.
    
    x0, y0, a, b, theta = p_ab
    
    sin = np.sin(np.radians(theta))
    cos = np.cos(np.radians(theta))
    
    x_global = cos*x_ell - sin*y_ell + x0
    y_global = sin*x_ell + cos*y_ell + y0
    
    return x_global, y_global

###################################################################################################################################################################
def segment_area2(p1_ab, p2_ab, x_inter, y_inter):
    
    x0, y0, a, b, theta = p1_ab
    p2_AB = ab2AB(p2_ab)
    
    x1_ell, y1_ell = map_onto_ell(p1_ab, x_inter[0], y_inter[0])
    x2_ell, y2_ell = map_onto_ell(p1_ab, x_inter[1], y_inter[1])
    
    t1_ell = np.arctan2(y1_ell*a/b, x1_ell)
    t2_ell = np.arctan2(y2_ell*a/b, x2_ell)
    t_mid = (t1_ell + t2_ell) / 2

    xT_mid_ell, yT_mid_ell = a*np.cos(t_mid), b*np.sin(t_mid)              ## midpoint on This side
    xO_mid_ell, yO_mid_ell = a*np.cos(t_mid+np.pi), b*np.sin(t_mid+np.pi)  ## midpoint on Other side
    
    xT_mid, yT_mid = map_from_ell(p1_ab, xT_mid_ell, yT_mid_ell)               
    xO_mid, yO_mid = map_from_ell(p1_ab, xO_mid_ell, yO_mid_ell)
    
    # zT = ellipse_equation_general(p2_AB, xT_mid, yT_mid)
    # zO = ellipse_equation_general(p2_AB, xO_mid, yO_mid)

    zT = ellipse_equation(xT_mid, yT_mid, p2_AB)
    zO = ellipse_equation(xO_mid, yO_mid, p2_AB)
    
    if zT < zO:   ## this side of arc closer to center of ellipse 2
        theta_sector = np.abs(t1_ell - t2_ell)
        x_mid, y_mid = xT_mid, yT_mid
    else:         ## other side closer; flip to other side
        theta_sector = 2*np.pi - np.abs(t1_ell - t2_ell)
        x_mid, y_mid = xO_mid, yO_mid
   
    area_triangle = np.abs( x1_ell*y2_ell - x2_ell*y1_ell ) / 2            
    area_segment = theta_sector * a * b / 2 + np.sign(theta_sector-np.pi)*area_triangle
        
    return area_segment

###################################################################################################################################################################
def segment_area4(p1_ab, p2_ab, x_inter, y_inter):

    x0, y0, a, b, theta = p1_ab
    p2_AB = ab2AB(p2_ab)
    
    N = len(x_inter)
    
    x_ell, y_ell = map_onto_ell(p1_ab, x_inter, y_inter)
    t_ell = np.arctan2(y_ell*a/b, x_ell)

    x_ell = [x for t, x in sorted(zip(t_ell, x_ell))]
    y_ell = [y for t, y in sorted(zip(t_ell, y_ell))]
    t_ell = sorted(t_ell)
    
    x_mid, y_mid, z_mid = np.empty(N), np.empty(N), np.empty(N)
    t_mid, t_diff = np.empty(N), np.empty(N)
    area_triangle = np.empty(N)
   
    for i in range(N):
        j = (i+1)%N     ## periodic boundary N->0
        t_mid[i] = (t_ell[i] + t_ell[j])/2
        t_diff[i] = np.abs(t_ell[i] - t_ell[j])
        
        ## the last two points jump from +pi to -pi due to arctan2; fixing it
        if j == 0:
            t_mid[i] += np.pi
            t_diff[-1] = 2*np.pi - t_diff[-1]

        x_mid[i], y_mid[i] = map_from_ell(p1_ab, a*np.cos(t_mid[i]), b*np.sin(t_mid[i]))
        #z_mid[i] = ellipse_equation_general(p2_AB, x_mid[i], y_mid[i])

        z_mid[i] = ellipse_equation(x_mid[i], y_mid[i], p2_AB)

        area_triangle[i] = np.abs( x_ell[i]*y_ell[j] - x_ell[j]*y_ell[i] ) / 2 

    t_sign = np.sign( t_diff - np.pi )
    
    if (z_mid[0] < 0) & (z_mid[2] < 0):

        area_segment_01 = t_diff[0] * a * b / 2 + t_sign[0]*area_triangle[0]
        area_segment_23 = t_diff[2] * a * b / 2 + t_sign[2]*area_triangle[2]
        
        area_segments = area_segment_01 + area_segment_23
        return area_segments, x_ell, y_ell
    
    elif (z_mid[1] < 0) & (z_mid[3] < 0):
        
        area_segment_12 = t_diff[1] * a * b / 2 + t_sign[1]*area_triangle[1]
        area_segment_30 = t_diff[3] * a * b / 2 + t_sign[3]*area_triangle[3]
        
        area_segments = area_segment_12 + area_segment_30
        return area_segments, x_ell, y_ell

###################################################################################################################################################################
def overlap_area_HC(p1_ab, p2_ab, x_inter=[], y_inter=[]):
    
    x1, y1, a1, b1, theta1 = p1_ab
    x2, y2, a2, b2, theta2 = p2_ab
    
    area = -1
    
    area_1 = (np.pi*a1*b1)
    area_2 = (np.pi*a2*b2)

    num_inter = len(x_inter)
        
    if num_inter == 0:  ## 0 interception points

#         p1_AB = ab2AB(p1_ab)
#         p2_AB = ab2AB(p2_ab)

        p1_AB =ab2AB(x1, y1, a1, b1, theta1)
        p2_AB =ab2AB(x2, y2, a2, b2, theta2)
        
#         z1 = ellipse_equation_general(p1_AB, x2, y2)
#         z2 = ellipse_equation_general(p2_AB, x1, y1)

        z1 = ellipse_equation(x2, y2, p1_AB)
        z2 = ellipse_equation(x1, y1, p2_AB)
    
        if (z1<0) | (z2<0):  ## one inside another
            area = min(np.pi*a1*b1, np.pi*a2*b2)
        else:                ## separated
            area = 0.0

    elif num_inter == 1:  ## 1 interception points; same procedure as 0

#         p1_AB = ab2AB(p1_ab)
#         p2_AB = ab2AB(p2_ab)

        p1_AB =ab2AB(x1, y1, a1, b1, theta1)
        p2_AB =ab2AB(x2, y2, a2, b2, theta2)

#         z1 = ellipse_equation_general(p1_AB, x2, y2)
#         z2 = ellipse_equation_general(p2_AB, x1, y1)

        z1 = ellipse_equation(x2, y2, p1_AB)
        z2 = ellipse_equation(x1, y1, p2_AB)

        if (z1<0) | (z2<0):  ## one inside another
            area = min(np.pi*a1*b1, np.pi*a2*b2)
        else:                ## separated
            area = 0.0

    elif num_inter == 2: ## 2 interception points

        A1 = segment_area2(p1_ab, p2_ab, x_inter, y_inter)
        A2 = segment_area2(p2_ab, p1_ab, x_inter, y_inter)

        area = A1 + A2

    elif num_inter == 4: ## 4 interception points

        A1, x_ell, y_ell = segment_area4(p1_ab, p2_ab, x_inter, y_inter)
        A2, _____, _____ = segment_area4(p2_ab, p1_ab, x_inter, y_inter)

        A_quad = np.abs( (x_ell[0] - x_ell[2]) * (y_ell[1] - y_ell[3]) \
                        - (x_ell[1] - x_ell[3]) * (y_ell[0] - y_ell[2]) ) /2

        area = A1 + A2 + A_quad
        
    area_ellipses = (area_1 + area_2) - area 

    return area_ellipses, area #area = area overlap

###################################################################################################################################################################

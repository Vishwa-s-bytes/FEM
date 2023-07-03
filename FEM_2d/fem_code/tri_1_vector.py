import numpy as np
import numpy.polynomial.legendre as g

def gauss_quadrature(f,points,limits):
    a = limits[0]
    b = limits[1]
    J = (b - a)/2
    gauss = g.leggauss(points)
    integral = np.zeros_like(f(1))
    for i in range(points):
        integral+=gauss[1][i]*f((a+b)/2+(b-a)*gauss[0][i]/2)
    return (J)*integral

def Area(xy_coord):
    m = np.transpose([[1, 1, 1], xy_coord[:, 0], xy_coord[:, 1]])
    return np.linalg.det(m)/2

def shape(x,y,xy_coord):
    A = Area(xy_coord)
    x1, x2, x3 = xy_coord[:, 0]
    y1, y2, y3 = xy_coord[:, 1]
    N1 = (x2*y3-x3*y2+(y2-y3)*x+(x3-x2)*y)
    N2 = (x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y)
    N3 = (x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y)
    return np.array([N1,N2,N3])/(2*A)

def B_vector(xy_coord):
    A = Area(xy_coord)
    x1, x2, x3 = xy_coord[:, 0]
    y1, y2, y3 = xy_coord[:, 1]
    B_m = np.array([[(y2 - y3), (y3 - y1), (y1 - y2)],
                         [(x3 - x2), (x1 - x3), (x2 - x1)]]) / (2 * A)
    B = np.zeros([3, 6])
    B[0, [0, 2, 4]] = B_m[0, :]
    B[1, [1, 3, 5]] = B_m[1, :]
    B[2, [1, 3, 5]] = B_m[0, :]
    B[2, [0, 2, 4]] = B_m[1, :]
    return B

def element_stiffness_matrix(D,xy_coord):
    B_matrix = B_vector(xy_coord)
    K = Area(xy_coord)*np.dot(np.transpose(B_matrix), np.dot(D,B_matrix))
    return K

def traction_force(nodes, traction_elem, T, side):
    xy_coord = nodes[traction_elem]
    if side == "a":
        x1, y1 = xy_coord[0]
        x2, y2 = xy_coord[1]
    elif side == "b":
        x1, y1 = xy_coord[1]
        x2, y2 = xy_coord[2]
    elif side == "c":
        x1, y1 = xy_coord[2]
        x2, y2 = xy_coord[0]
    else:
        print("unknown element no")
        quit()

    if (x2 != x1):
        limits = [x2, x1]
        def inside_integral(x):
            y = ((y2 - y1) / (x2 - x1)) * (x - x1) + y1
            return np.transpose(shape(x, y, xy_coord))
    else:
        limits = [y2, y1]
        def inside_integral(y):
            return np.transpose(shape(x2, y, xy_coord))

    limits = np.array(limits)
    limits.sort()
    I = gauss_quadrature(inside_integral, 1, limits)
    fx = T[0] * I
    fy = T[1] * I
    f = np.zeros([6, 1])
    f[[0, 2, 4], 0] = fx
    f[[1, 3, 5], 0] = fy
    return f

#usage
# xy_coord = np.array([[2,0.5],[2,1],[0,1]])
# E = 3
# v = 0.3
# print(element_stiffness_matrix(plane_stress(E,v),xy_coord))
#
# def traction_eqn(x):
#     y = 1
#     return y
# print(traction_force(2,0, [0,2], xy_coord,y=traction_eqn))
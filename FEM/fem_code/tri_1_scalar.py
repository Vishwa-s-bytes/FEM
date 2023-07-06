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
    return J*integral

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
    return np.array([[N1,N2,N3]])/(2*A)

def B_scalar(xy_coord):
    A = Area(xy_coord)
    x1, x2, x3 = xy_coord[:, 0]
    y1, y2, y3 = xy_coord[:, 1]
    B_matrix = np.array([[(y2-y3),(y3-y1),(y1-y2)],
                [(x3-x2),(x1-x3),(x2-x1)]])/(2*A)
    return B_matrix

def element_stiffness_matrix(k,xy_coord):
    B_matrix = B_scalar(xy_coord)
    K = k*Area(xy_coord)*np.dot(np.transpose(B_matrix),B_matrix)
    return K

def traction(limits,q,xy_coord,y=None, x=None,):
    if y is not None:
        def inside_integral(x):
            return np.transpose(shape(x,y(x),xy_coord))
    else:
        def inside_integral(y):
            return np.transpose(shape(x(y), y, xy_coord))
    f = -q * gauss_quadrature(inside_integral,1,limits)
    return f

#usage
# k = 5
# xy_coord = np.array([[2,0.5],[2,1],[0,1]])
#
# print(element_stiffness_matrix(k,xy_coord))
#
# def traction_line(x):
#     y = 1
#     return y
# print(traction([0,2],20,xy_coord,y = traction_line))


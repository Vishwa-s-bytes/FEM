import numpy as np
import numpy.polynomial.legendre as g

def shape(eta,psi):
    N1 = (1 - psi) * (1 - eta) / 4
    N2 = (1 + psi) * (1 - eta) / 4
    N3 = (1 + psi) * (1 + eta) / 4
    N4 = (1 - psi) * (1 + eta) / 4
    N = [[N1,N2,N3,N4]]
    return N

def GN(eta,psi):
    return (1/4)*np.array([[eta - 1, 1 - eta, 1 + eta, -1 * eta - 1],
                           [psi - 1, -1 - psi, 1 + psi, 1 - psi]])
def B(eta,psi,xy_coord):
    dN = np.dot(np.linalg.inv(np.dot(GN(eta, psi), xy_coord)), GN(eta, psi))
    return dN

def gauss_quadrature(f,points,limits):
    a = limits[0]
    b = limits[1]
    J = (b - a)/2
    gauss = g.leggauss(points)
    integral = np.zeros_like(f(1))
    for i in range(points):
        integral+=gauss[1][i]*f((a+b)/2+(b-a)*gauss[0][i]/2)
    return J*integral

def quad_gauss_quadrature(f,points,xy_coord):
    gauss = g.leggauss(points)
    integral = np.zeros_like(f(1,1))
    for i in range(points):
        for j in range(points):
            psi = gauss[0][i]
            eta = gauss[0][j]
            integral+=gauss[1][i]*gauss[1][j]*np.linalg.det(np.dot(GN(eta,psi),xy_coord))*f(eta,psi)
    return integral

def element_stiffness_matrix(k,xy_coord):
    def inside_integral(eta,psi):
        B_matrix = B(eta,psi,xy_coord)
        return k*np.dot(np.transpose(B_matrix),B_matrix)
    K = quad_gauss_quadrature(inside_integral, 2, xy_coord)
    return K

def traction_force(q,psi=None,eta=None):
    if psi is not None:
        def inside_integral(eta):
            return np.transpose(shape(eta, psi))
    else:
        def inside_integral(psi):
            return np.transpose(shape(eta, psi))
    f = -q * gauss_quadrature(inside_integral, 1, [-1, 1])
    return f


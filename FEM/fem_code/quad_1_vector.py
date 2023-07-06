import numpy as np
import numpy.polynomial.legendre as g

def shape(eta,psi):
    N1 = (1 - psi) * (1 - eta) / 4
    N2 = (1 + psi) * (1 - eta) / 4
    N3 = (1 + psi) * (1 + eta) / 4
    N4 = (1 - psi) * (1 + eta) / 4
    N = [N1,N2,N3,N4]
    return N

def GN(eta,psi):
    return (1/4)*np.array([[eta - 1, 1 - eta, 1 + eta, -1 * eta - 1],
                           [psi - 1, -1 - psi, 1 + psi, 1 - psi]])

def B_vector(eta,psi,xy_coord):
    dN = np.dot(np.linalg.inv(np.dot(GN(eta, psi), xy_coord)), GN(eta, psi))
    B = np.zeros([3,8])
    B[0, [0, 2, 4, 6]] = dN[0, :]
    B[1, [1, 3, 5, 7]] = dN[1, :]
    B[2, [1, 3, 5, 7]] = dN[0, :]
    B[2, [0, 2, 4, 6]] = dN[1, :]
    return B

def gauss_quadrature(f,points,limits):
    a = limits[0]
    b = limits[1]
    J = (b - a)/2
    gauss = g.leggauss(points)
    integral = np.zeros_like(f(1))
    for i in range(points):
        integral+=gauss[1][i]*f((a+b)/2+(b-a)*gauss[0][i]/2)
    return J*integral

def quad_gauss_quadrature(f,p,xy_coord):
    # gauss = g.leggauss(points)
    gauss = [[p,-p],[1,1]]
    integral = np.zeros_like(f(1,1))
    for i in range(2):
        for j in range(2):
            psi = gauss[0][i]
            eta = gauss[0][j]
            integral+=gauss[1][i]*gauss[1][j]*np.linalg.det(np.dot(GN(eta,psi),xy_coord))*f(eta,psi)
    return integral

def element_stiffness_matrix(D,xy_coord,point = 1/np.sqrt(3)):
    def inside_integral(eta, psi):
        B_matrix = B_vector(eta, psi, xy_coord)
        return np.dot(np.transpose(B_matrix), np.dot(D,B_matrix))
    K = quad_gauss_quadrature(inside_integral, point, xy_coord)
    return K

def traction_force(nodes,traction_elem,T,psi=None,eta=None):
    xy = nodes[traction_elem]
    if psi is not None:
        if psi==1:
            dis = np.linalg.norm(xy[2] - xy[1])
        else:
            dis = np.linalg.norm(xy[3] - xy[0])

        def inside_integral(eta):
            return np.transpose(shape(eta, psi))
    else:
        if eta==1:
            dis = np.linalg.norm(xy[2] - xy[3])
        else:
            dis = np.linalg.norm(xy[1] - xy[0])

        def inside_integral(psi):
            return np.transpose(shape(eta, psi))

    I = abs(dis) * gauss_quadrature(inside_integral, 1, [-1,1])/2
    fx = T[0] * I
    fy = T[1] * I
    f = np.zeros([8, 1])
    f[[0, 2, 4, 6], 0] = fx
    f[[1, 3, 5, 7], 0] = fy
    return f


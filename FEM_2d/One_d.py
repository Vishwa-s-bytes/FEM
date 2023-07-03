import numpy as np
import numpy.polynomial.legendre as g

def gauss_quadrature(f,points,x_coord):
    a = x_coord[0]
    b = x_coord[1]
    J = (b - a)/2
    gauss = g.leggauss(points)
    integral = np.zeros_like(f(1))
    for i in range(points):
        integral+=gauss[1][i]*f((a+b)/2+(b-a)*gauss[0][i]/2)
    return J*integral

def d1_quad_shape(nodes,x):
    x1, x2, x3 = nodes[0], nodes[1], nodes[2]
    N1 = ((x - x2) * (x - x3)) / ((x1 - x2) * (x1 - x3))
    N2 = ((x - x1) * (x - x3)) / ((x2 - x1) * (x2 - x3))
    N3 = ((x - x1) * (x - x2)) / ((x3 - x1) * (x3 - x2))
    N = [[N1,N2,N3]]
    return N

def d1_quad_shape_grad(nodes,x):
    x1, x2, x3 = nodes[0], nodes[1], nodes[2]
    dN = [[(2 * x - x2 - x3) / ((x1 - x2) * (x1 - x3)), (2 * x - x1 - x3) / ((x2 - x1) * (x2 - x3)),
          (2 * x - x1 - x2) / ((x3 - x1) * (x3 - x2))]]
    return dN

def d1_linear_shape(nodes,x):
    x1, x2 = nodes[0], nodes[1]
    N1 = (x - x2) / (x1 - x2)
    N2 = (x - x1) / (x2 - x1)
    N = [[N1,N2]]
    return N

def d1_linear_shape_grad(nodes,x):
    x1, x2 = nodes[0], nodes[1]
    dN = np.array([[-1, 1]]) / (x2 - x1)
    return dN

def element_stiffness_matrix(E, Area, order_of_area_polynomial, nodes):
    if len(nodes)==3:
        def inside_integral(x):
            B = d1_quad_shape_grad(nodes,x)
            return Area(x)*E*np.dot(np.transpose(B),B)
    elif len(nodes)==2:
        def inside_integral(x):
            B = d1_linear_shape_grad(nodes,x)
            return Area(x)*E*np.dot(np.transpose(B),B)
    else:
        print("Incorrect no.of nodes per element")
        quit()

    points = int(np.ceil((2*(len(nodes)-2)+order_of_area_polynomial+1)/2))
    k_matrix = (gauss_quadrature(inside_integral,points,[nodes[0],nodes[-1]]))
    return k_matrix

#usage
E = 1
def A(x):
    return 2
order_of_area = 0 #2x
print(element_stiffness_matrix(E,A,order_of_area,[2,3]))




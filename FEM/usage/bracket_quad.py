import fem_code as fc
import matplotlib.pyplot as plt
import numpy as np

E = 1e7
v = 0.3
fx = 0
fy = -1e5

def magnitude(value):
    v1 = value[:,0,:]
    v2 = value[:, 1, :]
    v3 = value[:, 2, :]
    return np.sqrt(v1**2+v2**2+v3**2)

nodes,quad_elem = fc.read_mesh("bracket_quad.msh")

quad_k = fc.k_global(nodes, quad_elem, fc.plane_stress(E,v))
dispv = np.array([0,1,53,54,55,56,57,58,59,60,61,62])
dispv = np.concatenate((dispv*2,dispv*2+1))
dispv = dispv.astype("int32")
displacement = fc.disp_boundary(nodes,dispv,0)

force1 = fc.quad_traction_force(nodes,quad_elem[81],[fx,fy],eta=-1)

global_force = fc.global_force(nodes,quad_elem,[81],[force1])
r,u = fc.solve(quad_k, displacement, global_force)
strain = fc.strain(nodes,quad_elem,u)
stress = fc.stress(strain,fc.plane_stress(E,v))

ax = fc.set_axis()
plt.title('nodes and elements')
fc.show(nodes,quad_elem,ax)

ax1 = fc.set_axis()
plt.title('Displaced nodes and elements')
fc.show(fc.displaced_nodes(nodes,u),quad_elem,ax1,color="red")
ax2 = fc.set_axis()
plt.title('Total strain')
fc.quad_display_value(nodes,quad_elem,ax2,magnitude(strain),resolution=4)
ax3 = fc.set_axis()
plt.title('Total stress')
fc.quad_display_value(nodes,quad_elem,ax3,magnitude(stress),resolution=4)

fc.display_all()
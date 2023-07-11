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

nodes,tri_elem = fc.read_mesh("bracket_tri.msh")
tri_k = fc.k_global(nodes, tri_elem, fc.plane_stress(E,v))

dispv = np.array([0,1,53,54,55,56,57,58,59,60,61,62])
dispv = np.concatenate((dispv*2,dispv*2+1))
dispv = dispv.astype("int32")
displacement = fc.disp_boundary(nodes,dispv,0)
force1 = fc.tri_traction_force(nodes, tri_elem[186], [fx,fy], "c")
global_force = fc.global_force(nodes,tri_elem,[186],[force1])
r,u = fc.solve(tri_k, displacement, global_force)

strain = fc.strain(nodes,tri_elem,u)
stress = fc.stress(strain,fc.plane_stress(E,v))

ax1 = fc.set_axis()
plt.title('nodes and elements')
fc.show(nodes,tri_elem,ax1)

ax1 = fc.set_axis()
plt.title('Displaced nodes and elements')
fc.show(nodes,tri_elem,ax1)
fc.show(fc.displaced_nodes(nodes,u),tri_elem,ax1,color="red")
ax2 = fc.set_axis()
plt.title('Total strain')
fc.tri_display_value(nodes,tri_elem,ax2,magnitude(strain))
ax3 = fc.set_axis()
plt.title('Total stress')
fc.tri_display_value(nodes,tri_elem,ax3,magnitude(stress))


fc.display_all()
import fem_code as fc
import matplotlib.pyplot as plt

E = 1e7
v = 0.3
fx = 0
fy = -1e5

nodes,quad_elem = fc.read_mesh("quad.msh")
quad_k = fc.k_global(nodes, quad_elem, fc.plane_stress(E,v))
#For higher meshes use multiprocessing
# quad_k = fc.k_global_threading(nodes, quad_elem, gk.plane_stress(E,v))
displacement = fc.disp_boundary(nodes,[0,1,6,7,22,23],0)
force1 = fc.quad_traction_force(nodes,quad_elem[3],[fx,fy],psi=1)
force2 = fc.quad_traction_force(nodes,quad_elem[4],[fx,fy],psi=1)
global_force = fc.global_force(nodes,quad_elem,[3,4],[force1,force2])
r,u = fc.solve(quad_k, displacement, global_force)
strain = fc.strain(nodes,quad_elem,u)
stress = fc.stress(strain,fc.plane_stress(E,v))
ax1 = fc.set_axis()
plt.title('displaced nodes and elements')
fc.show(nodes,quad_elem,ax1)
fc.show(fc.displaced_nodes(nodes,u),quad_elem,ax1,color="red")
ax2 = fc.set_axis()
plt.title('x strain')
fc.quad_display_value(nodes,quad_elem,ax2,strain[:,0,:])
ax3 = fc.set_axis()
plt.title('y stress')
fc.quad_display_value(nodes,quad_elem,ax3,stress[:,1,:])

fc.display_all()

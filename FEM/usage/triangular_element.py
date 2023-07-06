import fem_code as fc
import matplotlib.pyplot as plt

E = 1e7
v = 0.3
fx = 0
fy = -1e5

nodes,tri_elem = fc.read_mesh("tri.msh")
tri_k = fc.k_global_threading(nodes, tri_elem, fc.plane_stress(E,v))
displacement = fc.disp_boundary(nodes,[0,1,6,7,22,23],0)
force1 = fc.tri_traction_force(nodes, tri_elem[10], [fx,fy], "c")
force2 = fc.tri_traction_force(nodes, tri_elem[16], [fx,fy], "c")
global_force = fc.global_force(nodes,tri_elem,[10,16],[force1,force2])
r,u = fc.solve(tri_k, displacement, global_force)
strain = fc.strain(nodes,tri_elem,u)
stress = fc.stress(strain,fc.plane_stress(E,v))
ax = fc.set_axis()
plt.title('triangular nodes and elements')
fc.show(nodes,tri_elem,ax,1)
ax1 = fc.set_axis()
plt.title('displaced nodes and elements')
fc.show(nodes,tri_elem,ax1)
fc.show(fc.displaced_nodes(nodes,u),tri_elem,ax1,color="red")
ax2 = fc.set_axis()
plt.title('x stress')
fc.tri_display_value(nodes,tri_elem,ax2,stress[:,0,:])
ax3 = fc.set_axis()
plt.title('x strain')
fc.tri_display_value(nodes,tri_elem,ax3,strain[:,0,:])

fc.display_all()

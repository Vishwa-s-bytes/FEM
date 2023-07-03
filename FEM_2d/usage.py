from fem_code import globalise as gk
from fem_code import tri_1_vector as tv
from fem_code import quad_1_vector as qv
from fem_code import display as dp
from fem_code import getvalue as gv

E = 1e7
v = 0.3
fx = 0
fy = -1e5

#triangular_element

nodes,tri_elem = gv.read_mesh("tri.msh","triangle")
tri_k = gk.k_global_threading(nodes, tri_elem, gk.plane_strain(E,v))
displacement = gk.disp_boundary(nodes,[0,1,6,7,22,23],0)
force1 = tv.traction_force(nodes, tri_elem[10], [fx,fy], "c")
force2 = tv.traction_force(nodes, tri_elem[16], [fx,fy], "c")
global_force = gk.global_force(nodes,tri_elem,[10,16],[force1,force2])
r,u = gk.solve(tri_k, displacement, global_force)
strain = gk.strain(nodes,tri_elem,u)
stress = gk.stress(strain,gk.plane_stress(E,v))
ax1 = dp.set_axis()
dp.show(nodes,tri_elem,ax1,1)
dp.show(dp.displaced_nodes(nodes,u),tri_elem,ax1,color="red")
ax2 = dp.set_axis()
dp.tri_display_value(nodes,tri_elem,ax2,strain[:,0,:])


#quadrilateral_element

# nodes,quad_elem = gv.read_mesh("quad.msh","quad")
# quad_k = gk.k_global(nodes, quad_elem, gk.plane_stress(E,v))
# # using threading for more no.of elements
# # quad_k = gk.k_global_threading(nodes, quad_elem, gk.plane_stress(E,v))
# displacement = gk.disp_boundary(nodes,[0,1,6,7,22,23],0)
# force1 = qv.traction_force(nodes,quad_elem[3],[fx,fy],psi=1)
# force2 = qv.traction_force(nodes,quad_elem[4],[fx,fy],psi=1)
# global_force = gk.global_force(nodes,quad_elem,[3,4],[force1,force2])
# r,u = gk.solve(quad_k, displacement, global_force)
# strain = gk.strain(nodes,quad_elem,u)
# stress = gk.stress(strain,gk.plane_stress(E,v))
# ax1 = dp.set_axis()
# dp.show(nodes,quad_elem,ax1,1)
# dp.show(dp.displaced_nodes(nodes,u),quad_elem,ax1,color="red")
# ax2 = dp.set_axis()
# dp.quad_display_value(nodes,quad_elem,ax2,strain[:,0,:])


dp.display_all()
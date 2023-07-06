from fem_code.getvalue import read_mesh
from fem_code.display import show,tri_display_value,quad_display_value,displaced_nodes,set_axis,display_all
from fem_code.tri_1_vector import traction_force as tri_traction_force
from fem_code.quad_1_vector import traction_force as quad_traction_force
from fem_code.globalise import plane_stress, plane_strain, k_global, k_global_threading, disp_boundary, global_force, solve, strain, stress

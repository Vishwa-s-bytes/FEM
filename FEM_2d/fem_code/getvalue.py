import meshio
# import matplotlib.pyplot as plt
# from . import display as dp



def read_mesh(file_path,element_type):
    mesh = meshio.read(file_path)
    nodes = mesh.points[:,[0,1]]
    elements = mesh.cells_dict[element_type]
    return nodes,elements

#usage
# nodes,elements = read_mesh("t1.msh","triangle")
# _,ax = plt.subplots()
# dp.show(nodes,elements,ax,1)
# plt.show()

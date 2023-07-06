import meshio
# import matplotlib.pyplot as plt
# from . import display as dp



def read_mesh(file_path):
    mesh = meshio.read(file_path)
    nodes = mesh.points[:,[0,1]]
    elem_dict = mesh.cells_dict
    if 'triangle' in elem_dict:
        elements = elem_dict["triangle"]
    else:
        elements = elem_dict["quad"]
    return nodes,elements

#usage
# nodes,elements = read_mesh("t1.msh","triangle")
# _,ax = plt.subplots()
# dp.show(nodes,elements,ax,1)
# plt.show()

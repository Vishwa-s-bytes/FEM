import meshio



def read_mesh(file_path):
    mesh = meshio.read(file_path)
    nodes = mesh.points[:,[0,1]]
    elem_dict = mesh.cells_dict
    if 'triangle' in elem_dict:
        elements = elem_dict["triangle"]
    else:
        elements = elem_dict["quad"]
    return nodes,elements

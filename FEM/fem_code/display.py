import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def show(nodes,elem,ax,numbering=0,color="black"):
    if numbering == 1:
        for i in range(len(nodes)):
            ax.text(nodes[i][0],nodes[i][1],str(i),color="black")
        for i in range(len(elem)):
            points = nodes[elem[i]]
            x1 = np.mean(points[:, 0])
            y1 = np.mean(points[:, 1])
            ax.text(x1,y1,str(i),color="red")
            x2 = (x1 + points[0][0]) / 2
            y2 = (y1 + points[0][1]) / 2
            ax.text(x2,y2,"*",color="red")
    else:
        pass
    ax.set_aspect('equal')
    ax.scatter(nodes[:,0],nodes[:,1])
    for i in range(len(elem)):
        coord = nodes[elem[i]]
        x = coord[:, 0]
        y = coord[:, 1]
        ax.fill(x, y, edgecolor = color, fill = False)

def displaced_nodes(nodes,disp):
    nodes = np.array(nodes)
    disp1 = np.reshape(disp,(len(disp)//2,2))
    return nodes+disp1

def tri_display_value(nodes,elem,ax,value):
    ax.set_aspect('equal')
    ax.scatter(nodes[:, 0], nodes[:, 1])
    cmap = cm.get_cmap('jet')
    norm = plt.Normalize(value.min(), value.max())
    colors = cmap(norm(value))
    for i in range(len(elem)):
        coord = nodes[elem[i]]
        x = coord[:, 0]
        y = coord[:, 1]
        ax.fill(x, y, color = colors[i])
    scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    plt.colorbar(scalar_map, ax=plt.gca())
    return None

def quad_display_value(nodes,elem,ax,value,resolution = 20):
    ax.set_aspect('equal')
    ax.scatter(nodes[:, 0], nodes[:, 1])
    vmin = min(value.flatten())
    vmax = max(value.flatten())
    for i in range(len(elem)):
        co = nodes[elem[i]]
        corner = value[i]
        a = np.repeat([np.linspace(0, 1, resolution)], resolution, axis=0)
        b = np.repeat(np.linspace(0, 1, resolution), resolution)
        a = np.hstack(a)
        n = len(a)
        inter = (1 - a) * (1 - b) * corner[0] + a * (1 - b) * corner[1] + b * (1 - a) * corner[3] + a * b * corner[2]
        ref = np.transpose(np.concatenate(([a], [b]), axis=0))
        psi = np.reshape(ref[:, 0] * 2 - 1, (n, 1))
        eta = np.reshape(ref[:, 1] * 2 - 1, (n, 1))
        N1 = (1 - psi) * (1 - eta) / 4
        N2 = (1 + psi) * (1 - eta) / 4
        N3 = (1 + psi) * (1 + eta) / 4
        N4 = (1 - psi) * (1 + eta) / 4
        N = np.concatenate((N1, N2, N3, N4), axis=1)
        xy = np.dot(N, co)
        plt.scatter(xy[:, 0], xy[:, 1], c=inter, cmap="jet",vmin=vmin,vmax=vmax)
    plt.colorbar()
    return None

def set_axis():
    _,ax = plt.subplots()
    return ax

def display_all():
    plt.show()
    return None

import numpy as np
from . import quad_1_vector as qv
from . import tri_1_vector as tv
import concurrent.futures

def plane_stress(E,v):
    D = np.array([[1,v,0],[v,1,0],[0,0,(1-v)/2]])*E/(1-v**2)
    return D

def plane_strain(E,v):
    D = np.array([[1-v,v,0],[v,1-v,0],[0,0,(1-2*v)/2]])*E/((1+v)*(1-2*v))
    return D

def k_global(nodes,elem,D):
    s = len(elem[0])
    e = len(elem)
    n = len(nodes)
    k = np.zeros([2*n,2*n])
    if s == 4:
        for i in range(e):
            ele_node = np.array(elem[i])
            knode = np.zeros([2*s])
            knode[np.arange(s)*2] = 2 * ele_node
            knode[np.arange(s)*2+1] = 2 * ele_node + 1
            knode = knode.astype("int32")
            xy_coord = nodes[elem[i]]
            update = qv.element_stiffness_matrix(D,xy_coord)
            k[knode[:, None], knode] += update
    elif s == 3:
        for i in range(e):
            ele_node = np.array(elem[i])
            knode = np.zeros([2 * s])
            knode[np.arange(s)*2] = 2 * ele_node
            knode[np.arange(s)*2+1] = 2 * ele_node + 1
            knode = knode.astype("int32")
            xy_coord = nodes[elem[i]]
            update = tv.element_stiffness_matrix(D, xy_coord)
            k[knode[:, None], knode] += update
    return k

def thread_K_compute(D, nodes, elem):
    l = len(elem[0]) * 2
    k_collection = np.zeros([len(elem), l, l])
    if len(elem[0]) == 3:
        for i in range(len(elem)):
            xy_coord = nodes[elem[i]]
            k_collection[i] = tv.element_stiffness_matrix(D, xy_coord)
    elif len(elem[0]) == 4:
        for i in range(len(elem)):
            xy_coord = nodes[elem[i]]
            k_collection[i] = qv.element_stiffness_matrix(D, xy_coord)
    return k_collection

def k_global_threading(nodes, elem, D, num_processes = 8):
    nodes = np.array(nodes)
    elem = np.array(elem)
    s = len(elem[0])
    n = len(nodes)
    k = np.zeros([2*n,2*n])
    arr_process = np.zeros(num_processes)
    arr_process[:] = len(elem) // (num_processes)
    arr_process[:len(elem) % (num_processes)] += 1
    arr_process = arr_process.astype("int32")
    processess = []
    slicing = [0]
    slice = 0
    for i in range(num_processes):
        slice += arr_process[i]
        slicing.append(slice)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in range(num_processes):
            f = executor.submit(thread_K_compute,D,nodes,elem[slicing[i]:slicing[i+1],:])
            processess.append(f)
    slice = 0
    for i in range(num_processes):
        k_collection = processess[i].result()
        for j in range(len(k_collection)):
            ele_node = np.array(elem[int(slice + j)])
            knode = np.zeros([2 * s])
            knode[np.arange(s)*2] = 2 * ele_node
            knode[np.arange(s)*2+1] = 2 * ele_node + 1
            knode = knode.astype("int32")
            k[knode[:, None], knode] += k_collection[j]
        slice += len(k_collection)
    return k

def disp_boundary(nodes,boundary,value):
    displacement = np.full((2 * len(nodes), 1), np.nan)
    displacement[boundary, :] = value
    return displacement

def assemble_force(global_force,force,elem):
    node = vector(elem)
    global_force[node,:] += force
    return global_force

def global_force(nodes,elements,elem_no,forces):
    gforce = np.zeros([2*len(nodes),1])
    for i in range(len(elem_no)):
        gforce = assemble_force(gforce,forces[i],elements[elem_no[i]])
    return gforce

def bounded_k(stiff_k, displacement):
    free = np.argwhere(np.isnan(displacement))[:, 0]
    return stiff_k[free[:, None], free]

def solve(stiff_k, displacement, global_force):
    free = np.argwhere(np.isnan(displacement))[:,0]
    displacement[free,:] = np.dot(np.linalg.inv(bounded_k(stiff_k, displacement)),global_force[free,:])
    total_force = np.dot(stiff_k,displacement)
    reaction = total_force - global_force
    return reaction, displacement

def vector(arr):
    arr = np.array(arr)
    s = len(arr)
    newarr = np.zeros([2*s])
    newarr[np.array(range(0,2*s,2))] = 2 * arr
    newarr[np.array(range(1,2*s,2))] = 2 * arr + 1
    newarr = newarr.astype("int32")
    return newarr

def strain(nodes,elem,disp):
    res = []
    for j in range(len(elem)):
        dispv = disp[vector(elem[j]), :]
        if len(elem[j]) == 4:
            s = np.zeros([3, 4])
            l = np.array([[-1,-1],[1,-1],[1,1],[-1,1]])#actual
            # l = np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]])/(3**0.5)
            for i in range(4):
                s[:,i] = np.dot(qv.B_vector(l[i][1], l[i][0], nodes[elem[j]]), dispv)[:,0]
        elif len(elem[j]) == 3:
            s = np.dot(tv.B_vector(nodes[elem[j]]), dispv)
        else:
            quit()
        res.append(s)
    return np.array(res)

def stress(strain, D):
    stress = np.zeros_like(strain)
    for i in range(len(stress)):
        stress[i] = np.dot(D,strain[i])
    return stress


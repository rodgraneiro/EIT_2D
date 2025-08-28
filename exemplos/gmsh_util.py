# Exemplos de uso da API do gmsh em Python
#
# Autor: Erick Dario León Bueno de Camargo, 27/ago/2025

import numpy as np
import gmsh # para instalar: pip install gmsh


# Abre a GUI do gmsh e mostra malha
def view(name_msh):
    gmsh.initialize()
    gmsh.open(name_msh)
    gmsh.fltk.run()
    gmsh.finalize()


# Cria malha 2D circular com eletrodos pontuais,
# formada por arcos entre os eletrodos.
# Nó central (terra) é identificado pelo número 10000
# Cada eletrodo é identificado por seu número, de 10001 a 10000+nelectrodes
# O domínio (sem objetos) é identificado por '1000'
#
# Autor: Erick Dario León Bueno de Camargo, 27/ago/2025
def create_2D_circular_mesh_point_electrodes(nelectrodes, view=False, oldversion=True, name='', radius=0.1, lc=1e-2):
    if name=='':
        ver = '2.2' if oldversion else '4.1'
        name = f'circular_el{nelectrodes}_r{radius:.2f}_lc{lc}_v{ver}'
    name_msh = f'malhas/{name}.msh'
    print(f"Creting mesh {name_msh}")

    gmsh.initialize()
    gmsh.model.add(name)
    p_center = gmsh.model.geo.addPoint(0, 0, 0, lc) # centro em (0,0,0)

    # cria pontos dos eletrodos:
    p_elect = {}
    for idx in range(nelectrodes):
        x = radius*np.cos(2.0*np.pi*idx/nelectrodes)
        y = radius*np.sin(2.0*np.pi*idx/nelectrodes)
        p_elect[idx] = gmsh.model.geo.addPoint(x, y, 0, lc)
    
    # cria arcos entre os eletrodos:
    l_arcs = {}
    arcs_vec = []
    for idx in range(nelectrodes):
        el_ini = idx
        el_fim = (idx+1)%nelectrodes # no último arco, termina com o primeiro eletrodo
        l_arcs[idx] = gmsh.model.geo.addCircleArc(p_elect[el_ini], p_center, p_elect[el_fim])
        arcs_vec.append(l_arcs[idx])

    # cria loop com todos os arcos:
    cl1 = gmsh.model.geo.addCurveLoop(arcs_vec)
    
    # definindo uma superfície plana
    ps1 = gmsh.model.geo.addPlaneSurface([cl1])

    # criando as estruturas de dados do gmsh
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0, [p_center], 2, ps1) # força ponto central (dim=0) a fazer parte da superfície (dim=2)

    # criando entidades físicas...
    # nó central (terra) é identificado pelo número 10000
    # cada eletrodo é identificado por seu número, de 10001 a 10000+nelectrodes
    # o domínio (sem objetos) é identificado por '1000'
    gmsh.model.addPhysicalGroup(0, [p_center], 10000, name='central_node')
    for idx in range(nelectrodes):
        gmsh.model.addPhysicalGroup(0, [p_elect[idx]], 10000+idx+1, name=f'electrode_{idx+1}')
    gmsh.model.addPhysicalGroup(2, [ps1], 1000, name='body')

    # gerando malha 2D
    gmsh.model.mesh.generate(2)
    # gravando malha
    if oldversion:
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
    gmsh.write(name_msh)
    if view: # se quiser visualizar no gmsh
        gmsh.fltk.run()
    # Ao final, finalizar API do gmsh
    gmsh.finalize()
    return name_msh

# retorna nós de um determinado elemento da malha
# ex: elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim,tag_body)
def get_nodes_tri(elemNodeTags, idx):
    nodes_per_element = 3 # triângulo
    idx_noh1 = nodes_per_element*idx
    idx_noh2 = nodes_per_element*idx + 1
    idx_noh3 = nodes_per_element*idx + 2
    noh1 = elemNodeTags[0][idx_noh1]
    noh2 = elemNodeTags[0][idx_noh2]
    noh3 = elemNodeTags[0][idx_noh3]
    return (noh1, noh2, noh3)

# retorna coordenadas de um determinado nó da malha
# ex:
# nodeTags, nodeCoords = gmsh.model.mesh.getNodesForPhysicalGroup(dim, 1000)
def get_coords_node(nodeCoords, node):
    coords_per_node = 3 # gmsh sempre usa 3 coordenadas (x,y,z), mesmo para malhas 2D
    idx_node = node-1 # numeração dos nós começa em 1
    idx_x = coords_per_node*idx_node
    idx_y = coords_per_node*idx_node + 1
    idx_z = coords_per_node*idx_node + 2
    x = nodeCoords[idx_x]
    y = nodeCoords[idx_y]
    z = nodeCoords[idx_z]
    return (x,y,z)

# Mostra  dados da malha.... funciona para qualquer versão de msh (2.2 ou 4.1)
def show_mesh_data(name_msh, mostra=2):
    gmsh.initialize()
    gmsh.open(name_msh)

    entities = gmsh.model.getEntities()

    # Mostra info de todos os elementos encontrados no msh
    for e in entities:
        dim = e[0]
        tag = e[1]
        physicalTags = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
        # Get the mesh nodes for the entity (dim, tag):
        nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, tag)
        # Get the mesh elements for the entity (dim, tag):
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
        # Caso necessário, pega tipo do elemento...
        mytype = gmsh.model.getType(dim, tag)

        if len(physicalTags): # Pertence a um grupo físico
            if physicalTags[0]==10000:
                print(f"Nó central encontrado (nó {nodeTags}). Coords: {nodeCoords}")
            elif physicalTags[0]>10000:
                el = physicalTags[0]-10000
                print(f"Eletrodo {el} encontrado (nó {nodeTags}). Coords: {nodeCoords}")
            elif physicalTags[0]==1000:
                numElem = sum(len(i) for i in elemTags)
                numNodes = len(nodeTags)
                print(f'Corpo encontrado. numElem {numElem} e numNodes {numNodes} coords {nodeCoords.shape}')
            else:
                print(f'physicalTag desconhecida: {physicalTags}')

    # pega só elementos do corpo:
    print('***** Elementos do corpo *****')
    minha_tag = 1000
    dim=2
    tags_body = gmsh.model.getEntitiesForPhysicalGroup(dim,minha_tag) # é só uma entidade, então pega só a 1a
    nodes_per_element = 3
    coords_per_node = 3

    for tag_body in tags_body:
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim,tag_body)
        nodeTags, nodeCoords = gmsh.model.mesh.getNodesForPhysicalGroup(dim, minha_tag)
        numElemTotal = sum(len(i) for i in elemTags)
        print(f'tag_body {tag_body} com {numElemTotal} elementos, {len(nodeTags)} nós e {len(nodeCoords)} coordenadas')

        for elemTagsvec in elemTags:
            numElem = len(elemTagsvec)
            print(f'  Bloco com {numElem} elementos: {min(elemTagsvec)} a {max(elemTagsvec)}')
            for idx,elemTag in enumerate(elemTagsvec):
                (noh1, noh2, noh3) = get_nodes_tri(elemNodeTags, idx)
                if idx < mostra or idx > numElem-mostra-1: # mostrando só primeiros e últimos
                    print(f'    Elemento {elemTag}, nós {noh1} {noh2} {noh3}:')
                    for node in (noh1, noh2, noh3):
                        (x,y,z) = get_coords_node(nodeCoords, node)
                        print(f'     Nó {node}: coordenadas [{x} {y} {z}]')
    
    gmsh.finalize()


# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 20:35:29 2025

@author: rodgr
"""

import numpy as np
import meshio
import elements

class MyMesh:
    def __init__(self, nome_msh=None):
        if not(nome_msh is None):
            self.MshFileName = nome_msh
        else:
            self.MshFileName = ""
        
        self.NumberOfElements = 0
        self.NumberOfElectrodes = 0
        self.NumberOfNodes = 0

        self.Elements = None
        self.dim = 0

        self.GndNode = 0
        self.FlagRhoBased = False
        self.KGlobal = None

        self.element_type = ''

        self.Coordinates = None
        self.ElectrodeNodes = None
        
        self.msh_physical_groups = None


    # Ajusta um valor por elemento
    # Cada elemento tem um valor distinto
    # - array: vetor de dimensão NumberOfElements
    def SetRhoElements(self, array):
        for idx in range(self.NumberOfElements):
            self.Elements[idx].SetRho(array[idx])

    # Ajusta um valor por elemento
    # Cada elemento tem um valor distinto
    # - array: vetor de dimensão NumberOfElements
    def SetSigmaElements(self, array):
        for idx in range(self.NumberOfElements):
            self.Elements[idx].SetSigma(array[idx])

    # Ajusta um valor igual para todos os elementos
    # - array: valor escalar (float)
    def SetRhoValue(self, value):
        for idx in range(self.NumberOfElements):
            self.Elements[idx].SetRho(value)
    
    # Ajusta um valor igual para todos os elementos
    # - array: valor escalar (float)
    def SetSigmaValue(self, value):
        for idx in range(self.NumberOfElements):
            self.Elements[idx].SetSigma(value)

    # Ajusta valor por physical entity:
    # dic: dicionário
    # ex:
    # dic = {}
    # dic[1] = 5 
    # dic[2] = 10
    def SetRhoPhysicaEntity(self, dic):
        for idx in range(self.NumberOfElements):
            tag = self.Elements[idx].PhysicalEntity
            self.Elements[idx].SetRho(dic[tag])

    # Ajusta valor por physical entity:
    # dic: dicionário
    # ex:
    # dic = {}
    # dic[1] = 0.2 
    # dic[2] = 0.1
    def SetSigmaPhysicaEntity(self, dic):
        print("CHAVES recebidas em dic:", sorted(dic.keys()))
        for idx in range(self.NumberOfElements):
            tag = self.Elements[idx].PhysicalEntity
            self.Elements[idx].SetSigma(dic[tag])


    def CalcKGlobal(self):
        if self.Elements[0].Rho == 0.0:
            raise Exception("MyMesh:CalcKGlobal(): Valor de Rho/Sigma nao definido.")
            
        if self.KGlobal is None:
            self.KGlobal = np.zeros((self.NumberOfNodes, self.NumberOfNodes), dtype=float)
        else:
            self.KGlobal.fill(0)
            
        for elem in range(self.NumberOfElements): # para cada elemento:

            for i in range(len(self.Elements[elem].Topology)): # para cada i (noh local):
                no_i = self.Elements[elem].Topology[i] # pega noh_i (noh global)

                for j in range(len(self.Elements[elem].Topology)): # para cada j (noh local):
                    no_j = self.Elements[elem].Topology[j] # pega noh_j (noh global)

                    if self.FlagRhoBased:
                        valor = self.Elements[elem].KGeo[i, j] / self.Elements[elem].Rho
                    else:
                        valor = self.Elements[elem].KGeo[i, j] * self.Elements[elem].Sigma
                    self.KGlobal[no_i, no_j] += valor


    def ReadMesh(self):
        raise NotImplementedError("A função ReadMesh() tem que ser implementada na subclasse.")
    '''
    ###############################################################################
    # Esta função aplica condições de contorno conforme exemplo abaixo:
    #
    # [ 1   0          0          0          0        ] [u1]   [u1]    
    # [ 0  k1+k2+k3    -k3        0         -k2       ] [u2] = [F2 +k1*u1] 
    # [ 0   -k3      k3+k5+k4     -k5        0        ] [u3]   [F3 +k4*u4] 
    # [ 0   0          -k5      k2+k6       -k6       ] [u4]   [F4 +k6*u5] 
    # [ 0   0           0          0         1        ] [u5]   [u5]
    ###############################################################################  
    def aplica_cond_contorno(I_CC,
                            Y,
                            n_nodes,
                            V_imposto):           # Aplica condições de contorno

    I_CC_cond_contorno = I_CC[:].copy()

    for [noh_cond_contorno,valor_cond_contorno] in V_imposto:
        for i in range(0, n_nodes):                    # corrige matriz de corrente
        I_CC_cond_contorno[i] = (
            I_CC_cond_contorno[i] - \
            Y[i][noh_cond_contorno]*valor_cond_contorno
        )

    for [noh_cond_contorno,valor_cond_contorno] in V_imposto:
        I_CC_cond_contorno[noh_cond_contorno] = (
            valor_cond_contorno )            # coloca valor de contorno conhecido

    Y_cond_contorno = Y[:].copy()                        # Criar matriz solução
    for [noh_cond_contorno,valor_cond_contorno] in V_imposto:
        for k in range(0, n_nodes):                # laço para zerar linha e coluna
            Y_cond_contorno[noh_cond_contorno][k] = 0
            Y_cond_contorno[k][noh_cond_contorno] = 0

        Y_cond_contorno[noh_cond_contorno][noh_cond_contorno] = 1
    return Y_cond_contorno, I_CC_cond_contorno
    ###############################################################################
    '''


'''
    Malhas 2D do Edson onde os eletrodos são os primeiros x pontos do msh
    Os elementos são sempre triangulares.
    Não possui modelo de eletrodo.
    O corpo (background) é physical_group 1.
    O índice do GND é o ponto = (número de eletrodos) + 1
    Os objetos são physical_group 2, 3, 4 etc.
'''
class PointElectrodes2DMeshEdson(MyMesh):
    def __init__(self, NumberOfEletrodes, nome_msh=None, altura2D = 0.1, useEdson=False):
        super().__init__(nome_msh)

        if type(NumberOfEletrodes) == int:
            self.NumberOfElectrodes = NumberOfEletrodes
        else:
            raise Exception("PointElectrodes2DMeshEdson(): Invalid NumberOfEletrodes.")

        self.altura2D = altura2D
        self.useEdson = useEdson
    

    '''
    - O corpo (background) é physical_group 1.
    - Os objetos são physical_group 2, 3, 4 etc.
    - os eletrodos são os primeiros (n_elect) pontos do msh
    - O índice do GND é o ponto = (n_elect) + 1
    '''
    def ReadMesh(self):
        if self.MshFileName == "":
            raise Exception("PointElectrodes2DMeshEdson(): MshFileName not defined.")

        print(f"Reading {self.MshFileName}.")
        self.__mshdata = meshio.read(self.MshFileName)

        # Check msh dimension
        if 'gmsh:physical' in self.__mshdata.cell_data_dict.keys():
            if 'tetra' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                raise Exception("PointElectrodes2DMeshEdson(): dimension identification error. Should be 2D mesh.")
            elif 'triangle' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                self.dim = 2 # 2D mesh
                self.element_type = 'triangle'
            else:
                raise Exception("PointElectrodes2DMeshEdson(): dimension identification error.")
        else:
            raise Exception("PointElectrodes2DMeshEdson(): invalid mesh (no physical entities found).")        

        self.Coordinates = self.__mshdata.points
        msh_topology = self.__mshdata.cells_dict[self.element_type]
        msh_physical_groups = self.__mshdata.cell_data_dict["gmsh:physical"][self.element_type]
        self.physical_tags = np.unique(msh_physical_groups)
        print(f"Physical tags found: {self.physical_tags}.")

        self.NumberOfNodes = self.Coordinates.shape[0]
        self.NumberOfElements = msh_topology.shape[0]

        print(f"{self.NumberOfElements} Elements and {self.NumberOfNodes} Nodes found.")

        self.ElectrodeNodes = np.arange(self.NumberOfElectrodes, dtype=int)
        self.GndNode = self.NumberOfElectrodes

        print(f"ElectrodeNodes: {self.ElectrodeNodes}")
        print(f"GndNode: {self.GndNode}")

        self.Elements = [None] * self.NumberOfElements # alocando vetor de elementos

        if self.useEdson:
            elements.LinearTriangleEdson.Coordinates = self.Coordinates
            elements.LinearTriangleEdson.Altura2D = self.altura2D # define a altura padrão como 1cm
        else:
            elements.LinearTriangle.Coordinates = self.Coordinates
            elements.LinearTriangle.Altura2D = self.altura2D # define a altura padrão como 1cm
        for idx in range(self.NumberOfElements):
            if self.useEdson:
                self.Elements[idx] = elements.LinearTriangleEdson()
            else:
                self.Elements[idx] = elements.LinearTriangle()
            self.Elements[idx].Topology = msh_topology[idx]
            self.Elements[idx].PhysicalEntity = msh_physical_groups[idx]
        

            self.Elements[idx].CalcCentroid()
            self.Elements[idx].CalcKgeo()



'''
    Malhas 2D do Edson onde os eletrodos são definidos por linhas para implementação do modelo Hua.
    - Os elementos são sempre triangulares.
    - Usa modelo de eletrodo Hua.
    - O corpo (background) é physical_group 1000.
    - Os objetos são physical_group 1001, 1002, 1003...
    - Os eletrodos são physical_group 5001, 5002, 5003 etc...
    - O gnd está o physical_group 10000.
'''
class HuaElectrodes2DMeshEdson(MyMesh):
    def __init__(self, NumberOfEletrodes, nome_msh=None, altura2D = 0.1):
        super().__init__(nome_msh)

        if type(NumberOfEletrodes) == int:
            self.NumberOfElectrodes = NumberOfEletrodes
        else:
            raise Exception("HuaElectrodes2DMeshEdson(): Invalid NumberOfEletrodes.")

        self.altura2D = altura2D
    
    '''
    - O corpo (background) é physical_group 1000 formado por triângulos.
    - Os objetos são physical_group 1001, 1002, 1003... formados por triângulos.
    - Os eletrodos linhas identificadas por physical_group 5001, 5002, 5003 etc...
    - O gnd está o physical_group 10000.
    '''
    def ReadMesh(self):
        if self.MshFileName == "":
            raise Exception("HuaElectrodes2DMeshEdson(): MshFileName not defined.")

        print(f"Reading {self.MshFileName}.")
        self.__mshdata = meshio.read(self.MshFileName)

        # Check msh dimension
        if 'gmsh:physical' in self.__mshdata.cell_data_dict.keys():
            if 'tetra' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                raise Exception("PointElectrodes2DMeshEdson(): dimension identification error. Should be 2D mesh.")
            elif 'triangle' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                self.dim = 2 # 2D mesh
                self.element_type = 'triangle'
            else:
                raise Exception("PointElectrodes2DMeshEdson(): dimension identification error.")
        else:
            raise Exception("PointElectrodes2DMeshEdson(): invalid mesh (no physical entities found).")  
        
        # Verifica se tem as linhas dos eletrodos
        if not ('line' in self.__mshdata.cell_data_dict["gmsh:physical"].keys()):
            raise Exception("HuaElectrodes2DMeshEdson(): Não encontrei as linhas.")    
        
        self.Coordinates = self.__mshdata.points
        msh_topology = self.__mshdata.cells_dict[self.element_type]

        self.msh_physical_groups = self.__mshdata.cell_data_dict["gmsh:physical"][self.element_type]
        self.physical_tags = np.unique(self.msh_physical_groups)
        print(f"msh_physical_groups found (type {self.element_type}): {self.msh_physical_groups}.")
        print(f"Physical tags found (type {self.element_type}): {self.physical_tags}.")


        physical_groups_lines = self.__mshdata.cell_data_dict["gmsh:physical"]['line']

        physical_tags_lines = np.unique(physical_groups_lines)
        physical_tags_points = np.unique(self.__mshdata.cell_data_dict["gmsh:physical"]['vertex'])

        print(f"Physical tags: lines: {physical_tags_lines}; points: {physical_tags_points}")

        # Verifica se o physical do GND está no arquivo msh
        if not (10000 in physical_tags_points):
            raise Exception("HuaElectrodes2DMeshEdson(): GND vertex not found.")  

        n_electrodes = len(physical_tags_lines)
        print(f"{n_electrodes} electrodes found.")

        n_nohs_msh = self.Coordinates.shape[0]
        n_elementos_msh = msh_topology.shape[0] 
        print(f"MSH file with {n_elementos_msh} elements and {n_nohs_msh} nodes.")

        self.NumberOfNodes = n_nohs_msh + n_electrodes # incluindo os nós virtuais
        self.ElectrodeNodes = np.arange(n_nohs_msh, n_nohs_msh + n_electrodes, dtype=int) # nessa malha os eletrodos são os últimos nós (virtuais))

        print(f"ElectrodeNodes: {self.ElectrodeNodes}")

        electrodes_topology = self.__mshdata.cells_dict['line'] # só linhas dos eletrodos
        points_topology = self.__mshdata.cells_dict['vertex']   # só pontos (deveria ser só o GND)

        n_elementos_eletrodos = electrodes_topology.shape[0] # só linhas 

        self.NumberOfElements = n_elementos_msh + n_elementos_eletrodos # inclui triângulos  e as linhas dos eletrodos

        print(f"{self.NumberOfElements} Elements and {self.NumberOfNodes} Nodes found on model.")

        self.GndNode = points_topology[0][0] # o primeiro noh do primeiro physical vertex
        print(f'GndNode: {self.GndNode}')

        self.Elements = [None] * self.NumberOfElements # alocando vetor de elementos (triângulos do meio + elementos dos eletrodos)


        elements.LinearTriangle.Coordinates = self.Coordinates
        elements.LinearTriangle.Altura2D = self.altura2D # define a altura padrão como 1cm

        elements.LinearLineHua.Coordinates = self.Coordinates
        elements.LinearLineHua.Altura2D = self.altura2D # define a altura padrão como 1cm

        # Pegando elementos triangulares:
        for idx in range(n_elementos_msh):
            self.Elements[idx] = elements.LinearTriangle()
            self.Elements[idx].Topology = msh_topology[idx]
            self.Elements[idx].PhysicalEntity = self.msh_physical_groups[idx]
            self.Elements[idx].CalcCentroid()
            self.Elements[idx].CalcKgeo()

        # Pegando elementos dos eletrodos:
        for idy in range(n_elementos_eletrodos): # idy começa em zero
            idx = idy + n_elementos_msh          # idx continua a partir de n_elementos_msh
            self.Elements[idx] = elements.LinearLineHua()
            self.Elements[idx].Topology = electrodes_topology[idy]
            self.Elements[idx].PhysicalEntity = physical_groups_lines[idy]
            self.Elements[idx].CalcCentroid()
            self.Elements[idx].CalcKgeo()

    
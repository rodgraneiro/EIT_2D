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
        self.msh_topology = None


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
            #print(f' elemento_{elem} {self.Elements[elem].Topology}')
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
        self.msh_topology = self.__mshdata.cells_dict[self.element_type]   # só dos triângulos

        self.msh_physical_groups = self.__mshdata.cell_data_dict["gmsh:physical"][self.element_type]   # só dos triângulos
        self.physical_tags = np.unique(self.msh_physical_groups)
        print(f"msh_physical_groups found (type {self.element_type}): {self.msh_physical_groups}.")
        print(f"Physical tags found (type {self.element_type}): {self.physical_tags}.")

        # pegando informações dos eletrodos (elementos linha)
        physical_groups_lines = self.__mshdata.cell_data_dict["gmsh:physical"]['line']                 # só das linhas (eletrodos = 500X)
        physical_tags_lines = np.unique(physical_groups_lines)
        physical_tags_points = np.unique(self.__mshdata.cell_data_dict["gmsh:physical"]['vertex'])     # só dos pontos (GND = 10000)
        print(f"Physical tags: lines: {physical_tags_lines}; points: {physical_tags_points}")

        # Verifica se o physical do GND está no arquivo msh
        if not (10000 in physical_tags_points):
            raise Exception("HuaElectrodes2DMeshEdson(): GND vertex not found.")  

        n_electrodes = len(physical_tags_lines)
        print(f"{n_electrodes} electrodes found.")

        n_nohs_msh = self.Coordinates.shape[0]
        n_elementos_msh = self.msh_topology.shape[0] 
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

        # Setando variáveis das classes elemento (a mesma variável para toda a classe)
        elements.LinearTriangle.Coordinates = self.Coordinates
        elements.LinearTriangle.Altura2D = self.altura2D # define a altura padrão como 1cm

        elements.LinearLineHua.Coordinates = self.Coordinates
        elements.LinearLineHua.Altura2D = self.altura2D # define a altura padrão como 1cm

        # Pegando elementos triangulares:
        for idx in range(n_elementos_msh):
            self.Elements[idx] = elements.LinearTriangle()
            self.Elements[idx].Topology = self.msh_topology[idx]                           # só dos triângulos
            self.Elements[idx].PhysicalEntity = self.msh_physical_groups[idx]         # só dos triângulos
            self.Elements[idx].CalcCentroid()
            self.Elements[idx].CalcKgeo()

        # Pegando elementos dos eletrodos:
        for idy in range(n_elementos_eletrodos): # idy começa em zero
            idx = idy + n_elementos_msh          # idx continua a partir de n_elementos_msh
            self.Elements[idx] = elements.LinearLineHua()
            self.Elements[idx].PhysicalEntity = physical_groups_lines[idy]
            noh_virtual = self.ElectrodeNodes[self.Elements[idx].PhysicalEntity - 5001]
            self.Elements[idx].Topology = np.append(electrodes_topology[idy], noh_virtual)
            self.Elements[idx].CalcCentroid()
            self.Elements[idx].CalcKgeo()

####################################################################
# A seguinte classe PointElectrodes1DMeshEdson() lê uma malha unidimensional
# - O corpo possui XXX physicals (1, 2, 3, 4 etc)
# - Pode ter YYY eletrodos (nós de medida ou injeção)
# - O gnd é sempre o primeiro noh do primeiro physical vertex
###################################################################

class PointElectrodes1DMeshEdson(MyMesh):
    
    def __init__(self, ElectrodeNodes, nome_msh=None, altura1D = 0.1):
    
        super().__init__(nome_msh)

        if ( type(ElectrodeNodes) == list ) or ( type(ElectrodeNodes) == np.ndarray ):
            self.NumberOfElectrodes = len(ElectrodeNodes)
            self.ElectrodeNodes = ElectrodeNodes
        else:
            raise Exception("PointElectrodes1DMeshEdson(): Invalid ElectrodeNodes.")
        
        print(f'NumberOfElectrodes: {self.NumberOfElectrodes}')

        self.altura1D = altura1D


    '''
    - O corpo possui XXX physicals (1, 2, 3, 4 etc)
    - Pode ter YYY eletrodos (nós de medida ou injeção)
    - O gnd é sempre o primeiro noh do primeiro physical vertex
    '''
    def ReadMesh(self):
        if self.MshFileName == "":
            raise Exception("PointElectrodes1DMeshEdson(): MshFileName not defined.")

        print(f"Reading {self.MshFileName}.")
        self.__mshdata = meshio.read(self.MshFileName)
        
        # Check msh dimension
        if 'gmsh:physical' in self.__mshdata.cell_data_dict.keys():
            if 'tetra' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                raise Exception("PointElectrodes1DMeshEdson(): dimension identification error. Should be 1D mesh.")
            elif 'triangle' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                raise Exception("PointElectrodes1DMeshEdson(): dimension identification error. Should be 1D mesh.")
            elif 'line' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                self.dim = 1 # 1D mesh
                self.element_type = 'line'
            else:
                raise Exception("PointElectrodes1DMeshEdson(): dimension identification error.")
        else:
            raise Exception("PointElectrodes1DMeshEdson(): invalid mesh (no physical entities found).")  
        
        
        self.Coordinates = self.__mshdata.points
        self.msh_topology = self.__mshdata.cells_dict[self.element_type]
        
        self.msh_physical_groups = self.__mshdata.cell_data_dict["gmsh:physical"][self.element_type]
        self.physical_tags = np.unique(self.msh_physical_groups)
        print(f"msh_physical_groups found (type {self.element_type}): {self.msh_physical_groups}.")
        print(f"Physical tags found (type {self.element_type}): {self.physical_tags}.")
        
        n_nohs_msh = self.Coordinates.shape[0]
        n_elementos_msh = self.msh_topology.shape[0] 
        print(f"MSH file with {n_elementos_msh} elements and {n_nohs_msh} nodes.")
        
        self.GndNode = self.msh_topology[0][0] # o primeiro noh do primeiro physical vertex
        print(f'GndNode: {self.GndNode}')
        
        self.NumberOfNodes = n_nohs_msh
        self.NumberOfElements = n_elementos_msh
        self.Elements = [None] * self.NumberOfElements # alocando vetor de elementos
        
        elements.LinearLineEdson.Coordinates = self.Coordinates
        elements.LinearLineEdson.Altura1D = self.altura1D # define a altura padrão como 1cm
        print(f'first five coordinates: {elements.LinearLineEdson.Coordinates[:5]}')
        print(f'Altura1D: {elements.LinearLineEdson.Altura1D}')
        
        # Pegando elementos lineares 'lines':
        for idx in range(n_elementos_msh):
            self.Elements[idx] = elements.LinearLineEdson()
            self.Elements[idx].Topology = self.msh_topology[idx]
            #print(f'ElementsTopo1 {self.Elements[idx].Topology}')
            self.Elements[idx].PhysicalEntity = self.msh_physical_groups[idx]
            self.Elements[idx].CalcCentroid()
            self.Elements[idx].CalcKgeo()
        
        
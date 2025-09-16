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
        self.MeshType = 0
        self.dim = 0

        self.GndNode = 0
        self.FlagRhoBased = False
        self.KGlobal = None

        self.element_type = ''

        self.Coordinates = None
        self.ElectrodeNodes = None

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
        for idx in range(self.NumberOfElements):
            tag = self.Elements[idx].PhysicalEntity
            self.Elements[idx].SetSigma(dic[tag])


    def CalcKGlobal(self):
        raise NotImplementedError("A função CalcKGlobal() tem que ser implementada na subclasse.")

    def ReadMesh(self):
        raise NotImplementedError("A função ReadMesh() tem que ser implementada na subclasse.")


'''
Malhas tipo 1:
    Malhas 2D do Edson onde os eletrodos são os primeiros x pontos do msh
    Os elementos são sempre triangulares.
    Não possui modelo de eletrodo.
    O corpo (background) é physical_group 1.
    Os objetos são physical_group 2, 3, 4 etc.
'''
class PointElectrodes2DMeshEdson(MyMesh):
    def __init__(self, NumberOfEletrodes, nome_msh=None, altura2D = 0.1):
        super().__init__(nome_msh)

        self.MeshType = 1
        if type(NumberOfEletrodes) == int:
            self.NumberOfElectrodes = NumberOfEletrodes
        else:
            raise Exception("PointElectrodes2DMeshEdson(): Invalid NumberOfEletrodes.")

        self.altura2D = altura2D
    

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
        elements.LinearTriangle.Coordinates = self.Coordinates
        elements.LinearTriangle.Altura2D = self.altura2D # define a altura padrão como 1cm
        for idx in range(self.NumberOfElements):
            self.Elements[idx] = elements.LinearTriangle()
            self.Elements[idx].Topology = msh_topology[idx]
            self.Elements[idx].PhysicalEntity = msh_physical_groups[idx]
        

            self.Elements[idx].CalcCentroid()
            self.Elements[idx].CalcKgeo()


    def CalcKGlobal(self):
        if self.Elements[0].Rho == 0.0:
            raise Exception("PointElectrodes2DMeshEdson:CalcKGlobal(): Valor de Rho/Sigma nao definido.")
            
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




'''
Malhas tipo 2:
    Malhas 2D do Edson onde os eletrodos são definidos pelos physical groups xxxx
    Os elementos são sempre triangulares.
    Usa modelo de eletrodo Hua.
    O corpo (background) é physical_group 1.
    Os objetos são physical_group 2, 3, 4 etc.
    Os eletrodos são physical_group 1001, 1002, 1003 etc...
'''
class HuaElectrodes2DMeshEdson(MyMesh):
    def __init__(self, NumberOfEletrodes, nome_msh=None, altura2D = 0.1):
        super().__init__(nome_msh)

        self.MeshType = 2
        if type(NumberOfEletrodes) == int:
            self.NumberOfElectrodes = NumberOfEletrodes
        else:
            raise Exception("PointElectrodes2DMeshEdson(): Invalid NumberOfEletrodes.")

        self.altura2D = altura2D

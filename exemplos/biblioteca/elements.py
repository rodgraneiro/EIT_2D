import numpy as np


class MyElement:
    Coordinates = None

    def __init__(self):
        self.Centroid = None
        self.ElementType = 0
        self.PhysicalEntity = None
        self.FlagIsElectrode = False

        self.Topology = None
        self.KGeo = None
        self.Rho = 0.0
        self.Sigma = 0.0
    
    def SetRho(self, value):
        if value == 0:
            raise Exception("MyElement(): SetRho(): value não pode ser 0.")
        self.Rho = value
        self.Sigma = 1.0/value

    def SetSigma(self, value):
        if value == 0:
            raise Exception("MyElement(): SetRho(): value não pode ser 0.")
        self.Sigma = value
        self.Rho = 1.0/value


    def CalcCentroid(self):
        if self.Topology is None:
            raise Exception("MyElement(): CalcCentroid(): Topology not defined.")
        
        self.Centroid = np.zeros(3, dtype=float)

        for i in range(self.Topology.size):
            self.Centroid[0] += self.Coordinates[int(self.Topology[i])][0]
            self.Centroid[1] += self.Coordinates[int(self.Topology[i])][1]
            self.Centroid[2] += self.Coordinates[int(self.Topology[i])][2]
        
        self.Centroid[0] /= self.Topology.size
        self.Centroid[1] /= self.Topology.size
        self.Centroid[2] /= self.Topology.size


    def CalcKgeo(self):
        raise NotImplementedError("A função CalcKgeo() tem que ser implementada na subclasse.")


class LinearTriangle(MyElement):
    Altura2D = 1.0

    def __init__(self):
        super().__init__()

    # Gera a matriz local dos elementos triangulares lineares com 3 nohs
    # Equacionamento deduzido na tese do Fernando Moura, Apendice A.1.1
    def CalcKgeo(self):
        area = 0.0
        mat_area = np.zeros((3, 3), dtype=float)
        coeficientes = np.zeros((2, 3), dtype=float)
        
        for i in range(3):
            mat_area[i][0] = 1.0
            mat_area[i][1] = self.Coordinates[self.Topology[i]][0]
            mat_area[i][2] = self.Coordinates[self.Topology[i]][1]
        
        area = np.linalg.det(mat_area) / 2.0
        
        if area < 0:
            area *= -1
        
        coeficientes[0][0] = self.Coordinates[self.Topology[1]][1] - self.Coordinates[self.Topology[2]][1]
        coeficientes[0][1] = self.Coordinates[self.Topology[2]][1] - self.Coordinates[self.Topology[0]][1]
        coeficientes[0][2] = self.Coordinates[self.Topology[0]][1] - self.Coordinates[self.Topology[1]][1]
        coeficientes[1][0] = self.Coordinates[self.Topology[2]][0] - self.Coordinates[self.Topology[1]][0]
        coeficientes[1][1] = self.Coordinates[self.Topology[0]][0] - self.Coordinates[self.Topology[2]][0]
        coeficientes[1][2] = self.Coordinates[self.Topology[1]][0] - self.Coordinates[self.Topology[0]][0]
        
        # MATRIZ DE RIGIDEZ DO ELEMENTO TRIANGULAR
        self.KGeo = np.zeros((3, 3), dtype=float)
        
        for i in range(3):
            for j in range(i, 3):
                value = (self.Altura2D / (4.0 * area)) * (coeficientes[0][i] * coeficientes[0][j] + coeficientes[1][i] * coeficientes[1][j])
                self.KGeo[j][i] = value
                self.KGeo[i][j] = value



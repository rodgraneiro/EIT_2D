# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 20:33:28 2025

@author: rodgr
"""

import numpy as np


class MyElement:
    Coordinates = None
    Altura2D = None
    Altura1D = None
    thetaAngle = None
    sigmaX = None
    sigmaY = None

    def __init__(self):
        self.Centroid = None
        self.ElementType = 0
        self.PhysicalEntity = None
        self.FlagIsElectrode = False

        self.Topology = None
        self.KGeo = None
        self.KHua = None ########################################################
        self.Rho = 0.0
        self.Sigma = 0.0
        
        self.msh_physical_groups = None  ########################################################
    
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
    
    def CalcKHua(self):
        raise NotImplementedError("A função CalcKgeo() tem que ser implementada na subclasse.")

class LinearTriangle(MyElement):
    #Altura2D = 1.0

    def __init__(self):
        super().__init__()

    # Gera a matriz local dos elementos triangulares lineares com 3 nohs
    # Equacionamento deduzido na tese do Fernando Moura, Apendice A.1.1
    # OBS: só parte geométrica (que nunca muda), sem multiplicar por sigma (que muda)
    #def CalcKgeo(self, elem):

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



class LinearTriangleEdson(MyElement):
    #Altura2D = 1.0

    def __init__(self):
        super().__init__()

    # Gera a matriz local dos elementos triangulares lineares com 3 nohs
    # Equacionamento deduzido no livro "Numerical techniques in electromagnetics-Sadiku(2000)"
    # Seção 6.2.2. Veja equacionamento  em  "Método Elementos Finitos 2D - Triângulos.pdf"
    ###############################################################################
    # Esta função calcula a matriz local de condutividade da malha de elementos
    # finitos
    #
    #            sigma_i                                      
    # Y_local = ------------- * [matriz 2d ***completar***]             
    #             4*A_i                                             
    ###############################################################################
    # OBS: só parte geométrica (que nunca muda), sem multiplicar por sigma (que muda)
    #def CalcKgeo(self, elem):
    
    def CalcKgeo(self):
        noh1 = int(self.Topology[0])
        noh2 = int(self.Topology[1])
        noh3 = int(self.Topology[2])
        #print(noh1, noh2, noh3)
        
        x = [self.Coordinates[noh1][0], self.Coordinates[noh2][0], self.Coordinates[noh3][0]]
        y = [self.Coordinates[noh1][1], self.Coordinates[noh2][1], self.Coordinates[noh3][1]]

        triangulo = np.array([ [1, x[0],  y[0]], [1, x[1],  y[1]], [1, x[2],  y[2]]], dtype=np.float64)
        area_triangulo = abs((np.linalg.det(triangulo))/2)
    
        Y_11 =  ( y[1]-y[2])**2 + (x[2]-x[1])**2
        Y_12 =  ( y[1]-y[2])*(y[2]-y[0])+(x[2]-x[1])*(x[0]-x[2])          #Y_21 = Y_12
        Y_13 =  ( y[1]-y[2])*(y[0]-y[1])+(x[2]-x[1])*(x[1]-x[0])          #Y_31 = Y_13
        Y_22 =  (y[2]-y[0])**2 + (x[0]-x[2])**2                           #Y_32 = Y_23
        #Y_23 =  (y[2]-y[0])*(y[1]- y[1])+(x[0]-x[2])*(x[1]-x[0])          # O ERRO ESTAVA AQUI. Onde lê-se "(y[1]- y[1])",
                                                                           # leia-se (y[0]- y[1])
        Y_23 =  (y[2]-y[0])*(y[0]- y[1])+(x[0]-x[2])*(x[1]-x[0])          # 
        Y_33 =  (y[0]- y[1])**2 + (x[1]-x[0])**2  

        '''
        self.KGeo = (self.Altura2D /(4.0*area_triangulo))*np.array([[Y_11, Y_12, Y_13], 
                                        [Y_12, Y_22, Y_23],      
                                        [Y_13, Y_23, Y_33]       
                                        ])
        '''
        self.KGeo = (1 /(4.0*area_triangulo*self.Altura2D))*np.array([[Y_11, Y_12, Y_13], 
                                        [Y_12, Y_22, Y_23],      
                                        [Y_13, Y_23, Y_33]       
                                        ])
        print(noh1,noh2,noh3)
        print('x',x)
        print('y',y)
        print('KGeo1 zz \n', self.KGeo)
    

# Equacionamento deduzido na tese do Fernando Moura, Apendice A.2.1
class LinearLineHua(MyElement):
    

    #sigma_linha = 0.1
    
    def __init__(self):
        super().__init__()
        self.FlagIsElectrode = True
        
    def CalcKgeo(self):    
        mtrz_lenth_a = np.zeros((2, 2), dtype=float)
        coeficientes = np.zeros((2,2), dtype=float)
        #self.Topology = np.append(self.Topology, (17))
        #self.msh_physical_groups = self.__mshdata.cell_data_dict["gmsh:physical"][self.element_type]
        for i in range(2):
            mtrz_lenth_a[i][0] =self.Coordinates[self.Topology[i]][0]
            mtrz_lenth_a[i][1] =self.Coordinates[self.Topology[i]][1]
            node1 = self.Topology[i]
            #self.msh_physical_groups = self.__mshdata.cell_data_dict["gmsh:physical"][self.element_type]
            #print('self.msh_physical_groups ???', self.msh_physical_groups)
            #print('nodes', node1)
            
        # comprimento 'a' sqtr( (x2-x1)^2 + (y2-y1)^2) )

        lenth_a = np.sqrt((mtrz_lenth_a[1][0] - mtrz_lenth_a[0][0])**2 +(mtrz_lenth_a[1][1] - mtrz_lenth_a[0][1])**2)            
        
        self.KGeo = np.zeros((3, 3), dtype=float)
        mHua = np.array([[2.0,1.0,-3.0],[1.0,2.0,-3.0],[-3.0,-3.0,6.0]])
        
        # MATRIZ DE RIGIDEZ DO ELEMENTO Hua
        self.KGeo = ((self.Altura2D*lenth_a)/6)*mHua
        #print('KGeo_Hua', self.KGeo)

    def CalcCentroid(self):
        if self.Topology is None:
            raise Exception("MyElement(): CalcCentroid(): Topology not defined.")
        
        self.Centroid = np.zeros(3, dtype=float)

        for i in range(self.Topology.size - 1): # o último nó é o nó virtual, não tem coordenada
            self.Centroid[0] += self.Coordinates[int(self.Topology[i])][0]
            self.Centroid[1] += self.Coordinates[int(self.Topology[i])][1]
            self.Centroid[2] += self.Coordinates[int(self.Topology[i])][2]
        
        self.Centroid[0] /= self.Topology.size
        self.Centroid[1] /= self.Topology.size
        self.Centroid[2] /= self.Topology.size

# Elemento linear unidimensional
# Equacionamento deduzido em ..

class LinearLineEdson(MyElement):
    
    def __init__(self):
        super().__init__()
        
        
    def CalcKgeo(self):
        mtrz_lenth_a = np.zeros((2, 2), dtype=float)
        coeficientes = np.zeros((2,2), dtype=float)
                    
        for i in range(2):
            mtrz_lenth_a[i][0] =self.Coordinates[self.Topology[i]][0]
            #print(f'mtrz_lenth_a {mtrz_lenth_a[i][0]}')
            mtrz_lenth_a[i][1] =self.Coordinates[self.Topology[i]][1]
            #print(f'mtrz_lenth_a {mtrz_lenth_a[i][1]}')
            node1 = self.Topology
            
            #print(f'node1 {node1}')
            #print('self.msh_physical_groups ???', self.msh_physical_groups)
        #print('nodes', node1)
            #print(f'groups_lines {self.PhysicalEntity[i]}')
            #n_nohs_msh = self.Coordinates.shape[0] 
        #print(f'mtrz_lenth_a {mtrz_lenth_a}')
        
        # comprimento 'a' sqtr( (x2-x1)^2 + (y2-y1)^2) )
    
        lenth_a = np.sqrt((mtrz_lenth_a[1][0] - mtrz_lenth_a[0][0])**2 +(mtrz_lenth_a[1][1] - mtrz_lenth_a[0][1])**2)            
        #print(f'lenth_a {lenth_a}')
###############################################################################
# Esta função calcula a matriz local de condutividade da malha de elementos
# finitos
#
#               A_i                                       
# Y_local = ------------- * [[1, -1], [-1, 1]]             
#                L                                             
###############################################################################        
        self.KGeo = np.zeros((2, 2), dtype=float)
        mtz = np.array([[1.0,-1.0],[-1.0,1.0]])
              
        
        # MATRIZ DE RIGIDEZ DO ELEMENTO Hua
        self.KGeo = ((self.Altura2D/lenth_a))*mtz
        #print('KGeo \n', self.KGeo)


################################################################################
######################## ANISOTROPICO ##########################################
################################################################################




class LinearTriangleAnisotropic(MyElement):
    #Altura2D = 1.0

    def __init__(self):
        super().__init__()

    # Gera a matriz local dos elementos triangulares lineares com 3 nohs
    # Equacionamento deduzido no livro "Numerical techniques in electromagnetics-Sadiku(2000)"
    # Seção 6.2.2. Veja equacionamento  em  "Método Elementos Finitos 2D - Triângulos.pdf"
    ###############################################################################
    # Esta função calcula a matriz local de condutividade da malha de elementos
    # finitos
    #
    #            sigma_i                                      
    # Y_local = ------------- * [matriz 2d ***completar***]             
    #             4*A_i                                             
    ###############################################################################
    # OBS: só parte geométrica (que nunca muda), sem multiplicar por sigma (que muda)
    #def CalcKgeo(self, elem):
    
    def CalcKgeo(self):
        noh1 = int(self.Topology[0])
        noh2 = int(self.Topology[1])
        noh3 = int(self.Topology[2])
        #print(noh1, noh2, noh3)
        
        x = [self.Coordinates[noh1][0], self.Coordinates[noh2][0], self.Coordinates[noh3][0]]
        y = [self.Coordinates[noh1][1], self.Coordinates[noh2][1], self.Coordinates[noh3][1]]

        triangulo = np.array([ [1, x[0],  y[0]], [1, x[1],  y[1]], [1, x[2],  y[2]]], dtype=np.float64)
        area_triangulo = abs((np.linalg.det(triangulo))/2)
        #print('area_triangulo',area_triangulo)
        B_l = (y[1]-y[2])
        B_m = (y[2]-y[0])
        B_n = (y[0]-y[1])
        G_l = (x[2]-x[1])
        G_m = (x[0]-x[2])
        G_n = (x[1]-x[0])
        Sx = self.sigmaX
        Sy = self.sigmaY
        atheta_deg = self.thetaAngle
        atheta =  np.deg2rad(atheta_deg)
        Sxx =  Sx*np.cos(atheta)**2 + Sy*np.sin(atheta)**2
        Sxy = Sx*np.sin(atheta)*np.cos(atheta) - Sy*np.sin(atheta)*np.cos(atheta)
        Syy = Sx*np.sin(atheta)**2 + Sy*np.cos(atheta)**2
        
        #print('Xs',x[0], x[1], x[2] )
        #print('Ys',y[0], y[1], y[2] )
        #print(Sxx,Sxy,Syy)
        '''
        C_11 = Sxx*B_l**2 + 2*B_l*G_l*Sxy + Syy*G_l**2
        C_12 = B_l*(Sxx*B_m + Sxy*G_m) + G_l*(Sxy*B_m + Syy*G_m)        #C_21 = C_12
        C_13 = B_l*(Sxx*B_n + Sxy*G_n) + G_l*(Sxy*B_n + Syy*G_n)        #C_31 = C_13
        C_21 = B_m*(Sxx*B_l + Sxy*G_l) + G_m*(Sxy*B_l + Syy*G_l)
        C_22 = Sxx*B_m**2 + 2*B_m*G_m*Sxy + Syy*G_m**2
        C_23 = B_m*(Sxx*B_n + Sxy*G_n) + G_m*(Sxy*B_n + Syy*G_n)        #C_32 = C_23
        C_31 = B_n*(Sxx*B_l + Sxy*G_l) + G_n*(Sxy*B_l + Syy*G_l)       
        C_32 = B_n*(Sxx*B_m + Sxy*G_m) + G_n*(Sxy*B_m + Syy*G_m)
        C_33 = Sxx*B_n**2 + 2*B_n*G_n*Sxy + Syy*G_n**2
        '''
        C_11 = Sxx*B_l**2 + 2.0*Sxy*B_l*G_l + Syy*G_l**2
        C_22 = Sxx*B_m**2 + 2.0*Sxy*B_m*G_m + Syy*G_m**2
        C_33 = Sxx*B_n**2 + 2.0*Sxy*B_n*G_n + Syy*G_n**2

        C_12 = Sxx*B_l*B_m + Sxy*(B_l*G_m + G_l*B_m) + Syy*G_l*G_m
        C_13 = Sxx*B_l*B_n + Sxy*(B_l*G_n + G_l*B_n) + Syy*G_l*G_n
        C_23 = Sxx*B_m*B_n + Sxy*(B_m*G_n + G_m*B_n) + Syy*G_m*G_n

        print('C_xx',C_11,C_12,C_13,C_22,C_23,C_33)
        
        self.KGeo = (self.Altura2D /(4.0*area_triangulo))*np.array([[C_11, C_12, C_13], 
                                [C_12, C_22, C_23],      
                                [C_13, C_23, C_33]       
                                ])

        print('KGeo1', self.KGeo)

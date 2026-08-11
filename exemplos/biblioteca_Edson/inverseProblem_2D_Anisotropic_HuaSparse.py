# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 13:56:12 2025

@author: rodgr
"""
# este comentário está no branch test


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
#from mpl_toolkits.mplot3d import art3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


import mesh
import elements
import gmsh
import sys
import plotly.graph_objects as go
import os
import glob
from collections import defaultdict
#import matplotlib.cm as cm


import pandas as pd

from scipy.sparse import lil_matrix, csr_matrix, eye, kron
from scipy.sparse.linalg import factorized, spsolve
from scipy.spatial import cKDTree

from datetime import datetime

timestamp = datetime.now().strftime("%Y%m%d_%H%M")


#import subprocess
#import os


class inverse_problem: 
    def __init__(self, mymesh, V_imposto=None, Pcorrente=None, formatImage = None, SkipPattern=None, VirtualNode = False, I =1.0e-3):
        if not hasattr(mymesh, "KGlobal"): # verifica se o objeto mymesh tem um atributo chamado KGlobal.
            raise TypeError("Parâmetro incorreto: mymesh.")

        self.mymesh = mymesh
        self.Vmedido = None
        self.TempJ = None
        self.KGlobalTemp = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes), dtype=float)
        self.imgFotmat = formatImage

        if V_imposto is None:
            V_imposto = [[mymesh.GndNode, 0.0]]
        self.V_imposto = V_imposto

        if Pcorrente is None:
            # ToDo: implementar a matriz de medidas com padrão pula informado (SkipPattern)
            self.corrente = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfElectrodes))
            
            if VirtualNode == True:
                ajuste = self.mymesh.NumberOfNodes-self.mymesh.NumberOfElectrodes
            else: ajuste = 0
            
            k = np.arange(self.mymesh.NumberOfElectrodes)  
            self.corrente[k+ajuste, k] = I
            self.corrente[((k + SkipPattern+1) % self.mymesh.NumberOfElectrodes)+ajuste, k] = -I
            #raise Exception("forward_problem: Pcorrente padrão ainda não implementada.")
        else:
            self.corrente = Pcorrente
        
        matriz_coordenadas = self.mymesh.Coordinates
        matriz_topologia = self.mymesh.msh_topology
        self.Y_jacobian = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))
        self.Y_temp = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))
        self.Y_Vcalc = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))
        self.jacob_1 = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfElements))   # inicia a variável 1 do jacobiano
        self.jacob_2 = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes))          # inicia a variável 2 do jacobiano
        self.length = np.zeros(self.mymesh.NumberOfElements)
        self.sigmaStar = np.full(self.mymesh.NumberOfElements, 1.0)                           # Monta vetor chute
        
        #global comprimento
        #comprimento = np.zeros(nro_elementos_J)
        
    #def calc_Y_local_1D(self, xl, xm, Area, sigma): # função p/ calc a matriz local p/ elemento 1D
    #  return ((Area*sigma)/(xm-xl))*np.array([[1, -1], [-1, 1]])
    ###############################################################################
 
    ###############################################################################
    # Essa função calcula as derivadas paraciais para montar a matriz Jacobiana p/ 
    # método direto de soma e sobreposição das matrizes locais.
    # A função pega as coordenadas do elemento atual e calcula a derivada parcial  
    # da matriz local atual  por meio da equação:
    #
    #        ∂                                                              
    # ---------- [ Y_local ] = ( A_i / L_i ) * [ [ 1  -1 ], [ -1  1 ] ]    
    #      ∂σ_i                                                            
    ###############################################################################
    #################
    def area_triangulo(x1, y1, x2, y2, x3, y3):                      # calcula a área do triângulo com 3 pontos
        matriz = [ [1, x1, y1], [1, x2, y2], [1, x3, y3]]

        area_triangulo = (np.linalg.det(matriz))/2
        return abs(area_triangulo)
    ################
            

    def calc_Y_jacobian(self):   # função para calcular a matriz global
        self.Y_jacobian = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))
        self.listJacobian = []
        #for i in range(0, self.mymesh.NumberOfElements):        # laço para montar a matriz global
        for i in range(0, self.mymesh.msh_topology.shape[0]):
        #for i in range(0, 1):        # laço para montar a matriz global
            #print('nro_elementos_J',i)
            node_l = int(self.mymesh.msh_topology[i][0]) +1       # pegar dados da matriz elementos
            node_m = int(self.mymesh.msh_topology[i][1]) +1
            node_n = int(self.mymesh.msh_topology[i][2]) +1
            #print(f'node_l {node_l}, node_m {node_m}, node_n {node_n}')
            x1, y1, z1 = self.mymesh.Coordinates[node_l-1]    # coordenada x do ponto 1
            x2, y2, z2 = self.mymesh.Coordinates[node_m-1]    # coordenada x do ponto 2
            x3, y3, z3 = self.mymesh.Coordinates[node_n-1]    # coordenada x do ponto 2

            tri_coord = [ [1, x1, y1], [1, x2, y2], [1, x3, y3]]
            
            area_triangle = (np.linalg.det(tri_coord))/2
            #print(area_triangle)
            #print(f'area_triangle {i} {area_triangle}')
            C_11 =  (y2-y3)**2 + (x3-x2)**2
            C_12 =  (y2-y3)*(y3-y1)+(x3-x2)*(x1-x3)               #C_21 = C_12
            C_13 =  (y2-y3)*(y1-y2)+(x3-x2)*(x2-x1)               #C_31 = C_13
            C_22 =  (y3-y1)**2 + (x1-x3)**2                       #C_32 = C_23
            C_23 =  (y3-y1)*(y1-y2)+(x1-x3)*(x2-x1)                   # Na tese está: (Bt_m*Bt_o)+(Gm_m*Gm_o)+(Dt_m*Dt_o)
            C_33 =  (y1-y2)**2 + (x2-x1)**2  

            
            Y_local = ( self.mymesh.altura2D /(4*area_triangle))*np.array([[C_11, C_12, C_13], 
                                            [C_12, C_22, C_23],      
                                            [C_13, C_23, C_33]]      
                                            )
            #print('Y_local_c',Y_local)
            self.Y_jacobian[node_l-1, node_l-1] += Y_local[0, 0]
            self.Y_jacobian[node_m-1, node_l-1] += Y_local[1, 0]
            self.Y_jacobian[node_n-1, node_l-1] += Y_local[2, 0]
           
            self.Y_jacobian[node_l-1, node_m-1] += Y_local[0, 1]
            self.Y_jacobian[node_m-1, node_m-1] += Y_local[1, 1]
            self.Y_jacobian[node_n-1, node_m-1] += Y_local[2, 1]

            
            self.Y_jacobian[node_l-1, node_n-1] += Y_local[0, 2]
            self.Y_jacobian[node_m-1, node_n-1] += Y_local[1, 2]
            self.Y_jacobian[node_n-1, node_n-1] += Y_local[2, 2]
            #temp = self.apply_boundary_conditions(self.Y_jacobian)
            temp = self.apply_boundary_conditions_sparse(self.Y_jacobian)
            self.listJacobian.append(temp)
            
            #print('Y_jacobian calc_Y_jacobian \n',self.Y_jacobian)
            
            self.Y_jacobian = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))
        #print('self.listJacobian.append \n',self.listJacobian)
        #return Y_jacobiano
    ###############################################################################
    '''
    def calc_Y_jacobian_anisotropic_hua(self):
        """
        Monta uma lista de derivadas parciais globais:
        uma para cada parâmetro reconstruído.
        
        Ordem sugerida das colunas:
        [elem0_sxx, elem0_sxy, elem0_syy,
        elem1_sxx, elem1_sxy, elem1_syy, ...]
        """

        self.listJacobian = []

        for elem_idx in range(self.mymesh.NumberOfElements):
            elem = self.mymesh.Elements[elem_idx]

            # eletrodo Hua não é parâmetro do inverso
            if elem.FlagIsElectrode:
                continue

            topo = elem.Topology

            # três derivadas locais
            dKe_list = [elem.Kxx, elem.Kxy, elem.Kyy]

            for dKe in dKe_list:
                dK_global = np.zeros(
                    (self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes),
                    dtype=float
                )

                for i in range(len(topo)):
                    no_i = int(topo[i])
                    for j in range(len(topo)):
                        no_j = int(topo[j])
                        dK_global[no_i, no_j] += dKe[i, j]

                #dK_global_cc = self.apply_boundary_conditions(dK_global.copy())
                #dK_global_cc = self.apply_boundary_conditions_sparse(dK_global.copy())
                dK_global_cc = lil_matrix((self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes), dtype=float)
                self.listJacobian.append(dK_global_cc)

        #print("Número de derivadas parciais =", len(self.listJacobian))
    '''
    def calc_Y_jacobian_anisotropic_hua(self):
    
        self.listJacobian = []
    
        for elem_idx in range(self.mymesh.NumberOfElements):
            elem = self.mymesh.Elements[elem_idx]
    
            if elem.FlagIsElectrode:
                continue
    
            topo = elem.Topology
            dKe_list = [elem.Kxx, elem.Kxy, elem.Kyy]
    
            for dKe in dKe_list:
    
                dK_global = lil_matrix(
                    (self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes),
                    dtype=float
                )
    
                for i in range(len(topo)):
                    no_i = int(topo[i])
                    for j in range(len(topo)):
                        no_j = int(topo[j])
                        dK_global[no_i, no_j] += dKe[i, j]
    
                dK_global_cc = self.apply_boundary_conditions_sparse(dK_global)
                self.listJacobian.append(dK_global_cc)






    ###############################################################################
    # Esta função aplica condições de contorno conforme exemplo abaixo:
    #
    # [ 1   0          0          0          0        ] [u1]   [u1]    
    # [ 0  k1+k2+k3    -k3        0         -k2       ] [u2] = [F2 +k1*u1] 
    # [ 0   -k3      k3+k5+k4     -k5        0        ] [u3]   [F3 +k4*u4] 
    # [ 0   0          -k5      k2+k6       -k6       ] [u4]   [F4 +k6*u5] 
    # [ 0   0           0          0         1        ] [u5]   [u5]
    ###############################################################################  
    def apply_boundary_conditions(self, Kmatrix):
        self.vetor_corrente_cond_contorno = self.corrente[:].copy()
        #print(f'Vetor de corrente: \n {self.vetor_corrente_cond_contorno}')

        #KJacobian = KJacobian      # Criar matriz solução

        # atualiza vetor de correntes:
        for [noh_cond_contorno,valor_cond_contorno] in self.V_imposto:   # Necessário quando valor imposto é diferente de zero
            for i in range(0, self.mymesh.NumberOfNodes):                    # corrige matriz de corrente
                self.vetor_corrente_cond_contorno[i] = (
                    self.vetor_corrente_cond_contorno[i] - \
                    Kmatrix[i][noh_cond_contorno]*valor_cond_contorno
                )

        for [noh_cond_contorno,valor_cond_contorno] in self.V_imposto:
            self.vetor_corrente_cond_contorno[noh_cond_contorno] = valor_cond_contorno  # coloca valor de contorno conhecido
        
        # atualiza matriz KGlobal:
        for [noh_cond_contorno,valor_cond_contorno] in self.V_imposto:
          for k in range(0, self.mymesh.NumberOfNodes):                # laço para zerar linha e coluna
              Kmatrix[noh_cond_contorno][k] = 0
              Kmatrix[k][noh_cond_contorno] = 0
              Kmatrix[noh_cond_contorno][noh_cond_contorno] = 1
        return Kmatrix
    ###############################################################################
    def apply_boundary_conditions_sparse(self, Kmatrix):

        K = Kmatrix.tolil()
        self.vetor_corrente_cond_contorno = self.corrente.copy()
    
        for noh, valor in self.V_imposto:
            noh = int(noh)
    
            if valor != 0:
                col = K[:, noh].toarray()
                self.vetor_corrente_cond_contorno -= col * valor
    
            self.vetor_corrente_cond_contorno[noh, :] = valor
    
        for noh, valor in self.V_imposto:
            noh = int(noh)
    
            K[noh, :] = 0.0
            K[:, noh] = 0.0
            K[noh, noh] = 1.0
    
        return K.tocsr()
    
    ###############################################################################
    # Essa função monta a matriz do Filtro Passa Alta
    # FPA = M - I
    # Onde I é a matriz Identidade e M a matriz gaussiana.
    ###############################################################################
    def calc_L2_gauss_1D(self, centroids_1D, std=0.12, tol = 1e-9):#, covariance_vector):
        """
        Calcula L2 Gaussiana para malha 1D de EIT em uma barra de 1 metro.

        Parâmetros:
        - std: desvio padrão do kernel Gaussiano
        - centroids_1D: vetor (nelements,) com posições dos centroides
        - covariance_vector: vetor (nelements,) com covariância anatômica

        Retorna:
        - F: matriz Filtro Passa Alta
        """
        #print('centroids_1D',centroids_1D)
        nelements = len(centroids_1D)
        #tol = 1e-9

        F = np.zeros((nelements, nelements), dtype=np.float32)

        for i in range(nelements):
            soma = 0.0
            ci = centroids_1D[i]

            for j in range(nelements):
                cj = centroids_1D[j]
                dist = np.abs(cj - ci)

                if dist <= 1.5 * std:
                    fator = 1.0 if i == j else np.exp(-dist**2 / (2 * std**2))
                    soma += fator
                    F[i, j] = fator

            for j in range(nelements):
                if i == j:
                    aux = (1.0 - F[i, j] / soma)
                    F[i, j] = aux if np.abs(aux) > tol else 0.0
                else:
                    aux = (-F[i, j] / soma)
                    F[i, j] = aux if np.abs(aux) > tol else 0.0

        return F
    ###############################################################################
    ###############################################################################
    # Essa função calcula a variação de sigma estimado para a próxima iteração
    #
    # Δ = α_k * [ (J_kᵗ W₁ J_k + λ² L₂ᵗ L₂ )⁻¹ ] *                                             
    #         [ J_kᵗ W₁ ( z - h(θ̂_k) )                                                        
    #           - λ² L₂ᵗ L₂ ( θ̂_k - θ* ) ]                                                    
    ###############################################################################
    def calc_delta_sigma(self, J_b, residue, L2, sigma_inicial, alpha=0.01, Lambda=0.006):
        zW1=np.eye(J_b.shape[0])
        JT = J_b.T
        zJTW = JT @ zW1
        zJTWJ = zJTW @ J_b
        zLTL = L2.T @ L2
        termo_L = (Lambda**2) * zLTL
        primeiroTermo = zJTWJ + termo_L
        inv_primeiroTermo = np.linalg.inv(primeiroTermo)

        JTW_zh = zJTW @ residue
        ztermo_reg = (sigma_inicial)# - self.sigmaStar)
        zregularizacao = (Lambda**2)*zLTL @ ztermo_reg
        segundoTermo = JTW_zh - zregularizacao
        return -alpha*(inv_primeiroTermo @ segundoTermo)
    ###############################################################################


    def CalcTempKGlobal(self, SigmaTemp):
        self.KGlobalTemp = np.zeros((self.mymesh.msh_topology.shape[0], self.mymesh.msh_topology.shape[0]), dtype=float)
        
        #print(f'(self.mymesh.Elements.KGeo {self.Elements.KGeo})')
        for elem in range(self.mymesh.msh_topology.shape[0]): # para cada elemento:
        #for elem in range(self.mymesh.NumberOfElements): # para cada elemento:
            #print(f' SigmaTemp {SigmaTemp}')
            for i in range(len(self.mymesh.Elements[elem].Topology)): # para cada i (noh local):
                no_i = self.mymesh.Elements[elem].Topology[i] # pega noh_i (noh global)
                #print(f' no_i = {no_i}')
                for j in range(len(self.mymesh.Elements[elem].Topology)): # para cada j (noh local):
                    no_j = self.mymesh.Elements[elem].Topology[j] # pega noh_j (noh global)
                    valorTemp = self.mymesh.Elements[elem].KGeo[i, j] * SigmaTemp[elem]
                    #print(f' self.Elements.KGeo {self.mymesh.Elements[elem].KGeo}')
                    #print(f' valorTemp \n {valorTemp}')
                    self.KGlobalTemp[no_i, no_j]  += valorTemp
        print(f' self.KGlobalTemp \n {self.KGlobalTemp.shape}')
        np.savetxt("CalcTempKGlobal.txt", self.KGlobalTemp, fmt="%e")
        return self.KGlobalTemp
    ###############################################################################
    
    
    ###############################################################################
    def CalcTempKGlobalAnisotropicHua(self, SigmaTemp):
        """
        Monta a matriz global completa:
        - triângulos físicos: dependem de SigmaTemp
        - eletrodos Hua: entram com KGeo constante
        """

        #print("SigmaTemp.shape =", SigmaTemp.shape)

        self.KGlobalTemp = np.zeros(
            (self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes),
            dtype=float
        )

        k_phys = 0

        for elem_idx in range(self.mymesh.NumberOfElements):
            elem = self.mymesh.Elements[elem_idx]

            # =========================
            # Elemento Hua
            # =========================
            if elem.FlagIsElectrode:
                Ke = elem.KGeo

            # =========================
            # Triângulo físico anisotrópico
            # =========================
            else:
                sxx, sxy, syy = SigmaTemp[k_phys]

                Kxx = elem.Kxx
                Kxy = elem.Kxy
                Kyy = elem.Kyy

                Ke = sxx * Kxx + sxy * Kxy + syy * Kyy

                k_phys += 1

            # =========================
            # Montagem global
            # =========================
            topo = elem.Topology
            for i in range(len(topo)):
                no_i = int(topo[i])

                for j in range(len(topo)):
                    no_j = int(topo[j])
                    self.KGlobalTemp[no_i, no_j] += Ke[i, j]

        linhas_zeradas = np.where(~self.KGlobalTemp.any(axis=1))[0]
        colunas_zeradas = np.where(~self.KGlobalTemp.any(axis=0))[0]

        #print("KGlobalTemp.shape =", self.KGlobalTemp.shape)
        #print("Linhas zeradas:", linhas_zeradas)
        #print("Colunas zeradas:", colunas_zeradas)

        #np.savetxt("CalcTempKGlobal.txt", self.KGlobalTemp, fmt="%e")

        return self.KGlobalTemp

    ###############################################################################
    def CalcTempKGlobalAnisotropicHua_sparse(self, SigmaTemp):

        N = self.mymesh.NumberOfNodes
        K = lil_matrix((N, N), dtype=float)
    
        k_phys = 0
    
        for elem_idx in range(self.mymesh.NumberOfElements):
            elem = self.mymesh.Elements[elem_idx]
    
            if elem.FlagIsElectrode:
                Ke = elem.KGeo
            else:
                sxx, sxy, syy = SigmaTemp[k_phys]
                Ke = sxx * elem.Kxx + sxy * elem.Kxy + syy * elem.Kyy
                k_phys += 1
    
            topo = elem.Topology
    
            for i, no_i in enumerate(topo):
                no_i = int(no_i)
                for j, no_j in enumerate(topo):
                    no_j = int(no_j)
                    K[no_i, no_j] += Ke[i, j]
    
        return K.tocsr()
    ###############################################################################
    def Calc_J(self, invVtemp):
        listTempJ = []

        for idx in range(len(self.listJacobian)):
            termo1 = np.dot(invVtemp, self.vetor_corrente_cond_contorno)
            termo2 = np.dot(self.listJacobian[idx], termo1)
            termo3 = -np.dot(invVtemp, termo2)

            # mantém apenas os nós de eletrodo medidos
            termo3 = termo3[self.mymesh.ElectrodeNodes]
            termo3 = termo3.reshape(-1, 1, order='F')

            listTempJ.append(termo3)

        self.TempJ = np.concatenate(listTempJ, axis=1)
        self.JTJ = np.dot(self.TempJ.T, self.TempJ)

        #print("TempJ.shape =", self.TempJ.shape)
        #print("JTJ.shape =", self.JTJ.shape)
        '''
        J = self.TempJ

        corr_xx_yy = np.corrcoef(J[:, 0::3].ravel(), J[:, 2::3].ravel())[0, 1]
        print("corr σxx-σyy =", corr_xx_yy)
        '''
        J = self.TempJ

        a = J[:, 0::3].ravel()
        b = J[:, 2::3].ravel()
        
        if np.std(a) > 1e-14 and np.std(b) > 1e-14:
            corr_xx_yy = np.corrcoef(a, b)[0, 1]
        else:
            corr_xx_yy = np.nan
        
        print("corr σxx-σyy =", corr_xx_yy)
        
    ###############################################################################    
        
    def Calc_J_sparse(self, solve_K):
    
        listTempJ = []
    
        termo1 = solve_K(self.vetor_corrente_cond_contorno)
    
        for dK in self.listJacobian:
    
            dK = dK.tocsr()
    
            termo2 = dK @ termo1
            termo3 = -solve_K(termo2)
    
            termo3 = termo3[self.mymesh.ElectrodeNodes]
            termo3 = termo3.reshape(-1, 1, order='F')
    
            listTempJ.append(termo3)
    
        self.TempJ = np.concatenate(listTempJ, axis=1)
        self.JTJ = self.TempJ.T @ self.TempJ
    
        #print("norm J =", np.linalg.norm(self.TempJ))
        #print("min J =", np.min(self.TempJ))
        #print("max J =", np.max(self.TempJ))

    ###############################################################################
    ###############################################################################
    # Essa função calcula FPA com distância de cada elemento
    ############################################################################### 
    
    def calc_L2_gauss_2D(self, centroids_2D, std=0.007, tol=1e-9):
        
        #nelements = centroids_2D.shape[0]
        L2 = np.zeros((self.mymesh.msh_topology.shape[0], self.mymesh.msh_topology.shape[0]), dtype=np.float32)
    
        for i in range(self.mymesh.msh_topology.shape[0]):
            ci = centroids_2D[i]
            soma = 0.0
    
            # --- primeira passada: calcula fatores gaussianos ---
            for j in range(self.mymesh.msh_topology.shape[0]):
                cj = centroids_2D[j]
                #print(f'cj = {cj}')
                dist = np.linalg.norm(ci - cj)  # distância Euclidiana
                #print(f'dist = {dist}')
                if dist <= 5.0 * std:
                    #print(f'dist if = {2.0 * std}')
                    fator = 1.0 if i == j else np.exp(-dist**2 / (2 * std**2))
                    soma += fator
                    L2[i, j] = fator
    
            # --- segunda passada: normaliza e transforma em passa-alta ---
            for j in range(self.mymesh.msh_topology.shape[0]):
                if soma > 0:
                    if i == j:
                        aux = 1.0 - L2[i, j] / soma
                    else:
                        aux = -L2[i, j] / soma
                    L2[i, j] = aux if np.abs(aux) > tol else 0.0
       
        # plot  matrix sparsity 
        plt.figure(figsize=(6, 5))
        plt.spy(L2, markersize=1)
        plt.title('HPFilter matrix sparsity pattern', fontsize=15)
        plt.xlabel('Colun', fontsize=12)
        plt.ylabel('Line', fontsize=12)
        #plt.tight_layout()
        plt.show(block=False)   # mostra sem travar
        plt.pause(3)            # mantém aberto por 3 segundos
        plt.close('all')        # fecha automaticamente


        N = L2.shape[0]
        X, Y = np.meshgrid(np.arange(N), np.arange(N))
        '''
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        surf = ax.plot_surface( X, Y, L2, cmap='viridis',  linewidth=0, antialiased=True )
        
        ax.set_xlabel('Coluna')
        ax.set_ylabel('Linha')
        #ax.set_zlabel('L[i,j]')
        ax.set_title('HPFilter – Superfície 3D')
        
        fig.colorbar(surf, shrink=0.5)
        #plt.pause(0.1)
        #plt.show()
        #plt.pause(0.1)
        #plt.show(block = False)
        #plt.pause(0.1)
        plt.show(block=False)   # mostra sem travar
        plt.pause(3)            # mantém aberto por 3 segundos
        plt.close('all')        # fecha automaticamente
        '''
        return L2

    
    ###############################################################################
    # Essa função calcula FPA com distância média entre os elementoss
    ###############################################################################    

    def calc_L2_gauss_2D_only_domain(self, centroids_2D, std=None, tol=1e-9, raio=1.5):
        """
        Filtro passa-alta gaussiano apenas no domínio físico.
        Exclui elementos de eletrodo Hua.
        """

        # índices dos elementos que realmente pertencem ao domínio da imagem
        idx_dom = np.array(
            [i for i, elem in enumerate(self.mymesh.Elements)
            if not elem.FlagIsElectrode],
            dtype=int
        )

        cent_dom = centroids_2D[idx_dom, :2]   # usa só x,y
        n_dom = len(idx_dom)

        # ------------------------------------------------------------
        # std automático baseado na distância ao vizinho mais próximo
        # ------------------------------------------------------------
        if std is None:
            dmins = []
            for i in range(n_dom):
                ci = cent_dom[i]
                dmin = np.inf
                for j in range(n_dom):
                    if i == j:
                        continue
                    cj = cent_dom[j]
                    dij = np.linalg.norm(ci - cj)
                    if dij < dmin:
                        dmin = dij
                if np.isfinite(dmin):
                    dmins.append(dmin)

            dmed = np.mean(dmins)
            std = 4.5 * dmed
            print('WARNING: Std', std)
        else: 
            print('Std', std)

        # ------------------------------------------------------------
        # monta gaussiana normalizada
        # ------------------------------------------------------------
        W = np.zeros((n_dom, n_dom), dtype=np.float64)

        for i in range(n_dom):
            ci = cent_dom[i]
            soma = 0.0

            for j in range(n_dom):
                cj = cent_dom[j]
                dist = np.linalg.norm(ci - cj)

                if dist <= raio * std:
                    g = np.exp(-(dist**2) / (2.0 * std**2))
                    W[i, j] = g
                    soma += g

            if soma > 0.0:
                W[i, :] /= soma

        W[np.abs(W) < tol] = 0.0

        # filtro passa-alta
        L2 = np.eye(n_dom) - W


        '''
        # plot  matrix sparsity 
        plt.figure(figsize=(6, 5))
        plt.spy(L2, markersize=1)
        plt.title('HPFilter matrix sparsity pattern', fontsize=15)
        plt.xlabel('Colun', fontsize=12)
        plt.ylabel('Line', fontsize=12)
        #plt.tight_layout()
        plt.show(block=False)   # mostra sem travar
        plt.pause(3)            # mantém aberto por 3 segundos
        plt.close('all')        # fecha automaticamente


        N = L2.shape[0]
        X, Y = np.meshgrid(np.arange(N), np.arange(N))
        '''
        
        '''
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        surf = ax.plot_surface( X, Y, L2, cmap='viridis',  linewidth=0, antialiased=True )
        
        ax.set_xlabel('Coluna')
        ax.set_ylabel('Linha')
        #ax.set_zlabel('L[i,j]')
        ax.set_title('HPFilter – Superfície 3D')
        
        fig.colorbar(surf, shrink=0.5)
        #plt.pause(0.1)
        #plt.show()
        #plt.pause(0.1)
        #plt.show(block = False)
        #plt.pause(0.1)
        plt.show(block=False)   # mostra sem travar
        plt.pause(3)            # mantém aberto por 3 segundos
        plt.close('all')        # fecha automaticamente
        '''
        return L2, idx_dom, std
    ###############################################################################
    # Essa função plota o gráfico convergência das iterações
    ###############################################################################

    def plotar_iteracoes(self, lista_indice, lista_valor, nome_arquivo="Optimization"):
        
        imgFormat = 'webp'

        fig, ax = plt.subplots(figsize=(6, 5))
    
        # NÃO usar aspect='equal' para gráfico de linhas
        # ax.set_aspect('equal', adjustable='box')
    
        ax.plot(
            lista_indice,
            lista_valor,
            marker='.',
            linestyle='-',
            color='b',
            label=r'Norm $\Delta\sigma$'
        )
    
        ax.set_xlabel("Iteration", fontsize=12)
        ax.set_ylabel(r"Norm [$\Delta\sigma$]", fontsize=12)
        ax.set_title("Optimization", fontsize=15)
    
        ax.legend()
        ax.grid(True)
    
        fig.tight_layout()
    
        ax.ticklabel_format(style='plain')
    
        #fig.savefig(f"{nome_arquivo}.svg", format="svg", dpi=200) 
        fig.savefig(f"{nome_arquivo}.{imgFormat}", format=imgFormat, dpi=200)
    
        plt.show(block=False)
        plt.pause(0.1)
    
        plt.close(fig)
    ###############################################################################
    ###############################################################################
    # Essa função plota o gráfico da condutividade da malha
    ###############################################################################

    def plotMSH(self, sigma, Lambda=None, iteration=None, save=False,   SigmaXXXYYY=None, DifAniso=None, nome_arquivo=None):


        x = self.mymesh.Coordinates[:, 0]
        y = self.mymesh.Coordinates[:, 1]
        topo = self.mymesh.msh_topology
        imgFormat = self.imgFotmat
        
        elems_2D = np.array([el for el in topo if len(el) == 3])
        elems_1D = np.array([el for el in topo if len(el) == 2])
    
        if len(elems_2D) == 0:
            print("Nenhum elemento triangular encontrado.")
            return
    
        triang = tri.Triangulation(x, y, elems_2D)
        ntri = triang.triangles.shape[0]
        fc = sigma.ravel()[:ntri]
    
        # ============================================================
        # Função interna para desenhar e salvar cada figura
        # ============================================================
        def desenhar_figura(tipo_escala):
    
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.set_aspect('equal', adjustable='box')
    
            # -------------------------------
            # ESCALA AUTOMÁTICA
            # -------------------------------
            if tipo_escala == "AutoScale":
    
                if SigmaXXXYYY in ('xy', 'θ°'):
                    tpc = ax.tripcolor(
                        triang,
                        facecolors=fc,
                        edgecolors='k',
                        linewidths=0.1,
                        cmap='RdBu_r'
                    )
                else:
                    tpc = ax.tripcolor(
                        triang,
                        facecolors=fc,
                        edgecolors='k',
                        linewidths=0.1,
                        cmap='rainbow'
                    )
    
            # -------------------------------
            # ESCALA FIXA / SAME SCALE
            # -------------------------------
            elif tipo_escala == "SameScale":
    
                if SigmaXXXYYY in ('xy', 'θ°'):
    
                    lim = np.max(np.abs(fc))
    
                    if lim == 0:
                        lim = 1e-12
    
                    norm = TwoSlopeNorm(
                        vmin=-lim,
                        vcenter=0,
                        vmax=lim
                    )
    
                    tpc = ax.tripcolor(
                        triang,
                        facecolors=fc,
                        edgecolors='k',
                        linewidths=0.1,
                        cmap='RdBu_r',
                        norm=norm
                    )
    
                else:
                    tpc = ax.tripcolor(
                        triang,
                        facecolors=fc,
                        edgecolors='k',
                        linewidths=0.1,
                        cmap='rainbow',
                        vmin=0.0,
                        vmax=4.0
                    )
    
            # -------------------------------
            # Plotar elementos 1D / eletrodos
            # -------------------------------
            if len(elems_1D) > 0:
                for n1, n2 in elems_1D:
                    x_coords = [x[n1], x[n2]]
                    y_coords = [y[n1], y[n2]]
                    ax.plot(x_coords, y_coords, color='red', linewidth=2)
    
            # -------------------------------
            # Colorbar
            # -------------------------------
            
            '''
            if SigmaXXXYYY != 'θ°':
                fig.colorbar(
                    tpc,
                    ax=ax,
                    shrink=0.70,
                    label='Conductivity σ [S/m]'
                )
            else:
                fig.colorbar(
                    tpc,
                    ax=ax,
                    shrink=0.70,
                    label='Angle θ°'
                )
            '''
            if SigmaXXXYYY != 'θ°':
                cbar = fig.colorbar(tpc, ax=ax, shrink=0.7)
                cbar.set_label('σ (Conductivity) (S/m)', fontsize=15)
            else:
                cbar = fig.colorbar(tpc, ax=ax, shrink=0.7)
                cbar.set_label('Angle θ°', fontsize=15)
            # -------------------------------
            # Título
            # -------------------------------
            if save:
                ax.set_title(
                    #f"σ{SigmaXXXYYY} - λ_{Lambda:.2e}-it_{iteration}",
                    f"Conductivity σ{SigmaXXXYYY}",
                    fontsize=15
                )
            else:
                ax.set_title("Conductivity Real (σ)", fontsize=15)
    
            ax.set_xlabel("Lenght [m]", fontsize=15)
            ax.set_ylabel("Lenght [m]", fontsize=15)
    
            plt.tight_layout()
            plt.ticklabel_format(style='plain')
    
            # -------------------------------
            # Salvar
            # -------------------------------
            if save:
    
                if nome_arquivo is None:
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
                    #nome_saida = f"Conductivity_{tipo_escala}_{timestamp}.svg" 
                    nome_saida = f"Conductivity_{tipo_escala}_{timestamp}.{imgFormat}" 
                    
                else:
                    #nome_saida = f"{nome_arquivo}_{tipo_escala}.svg"
                    nome_saida = f"{nome_arquivo}_{tipo_escala}.{imgFormat}"
                #plt.savefig(nome_saida, format="svg", dpi=200)
                plt.savefig(nome_saida, format=imgFormat, dpi=200)
    
            plt.show(block=False)
            plt.pause(0.1)
            plt.close(fig)
    
        # ============================================================
        # Gera as duas imagens
        # ============================================================
        desenhar_figura("AutoScale")
        desenhar_figura("SameScale")
        ###############################################################################
    
    ###############################################################################   
    def plotElipse(self, sigma, sigma_L, sigma_T, theta_deg, Lambda = None, iteration = None, save = False, DifAniso=None, nome_arquivo= None):
        dados_elipses_phantom = np.loadtxt("dados_elipses_phantom.txt")
        #print(f'dados_elipses_phantom\n {dados_elipses_phantom}')
        
        x, y = self.mymesh.Coordinates[:, 0], self.mymesh.Coordinates[:, 1]
        topo = self.mymesh.msh_topology
        ajuste  = self.mymesh.NumberOfElectrodes*4 +1
        imgFormat = self.imgFotmat
    
        # --- Separar elementos 2D (triangulares) e 1D (linhas) ---
        elems_2D = np.array([el for el in topo if len(el) == 3])
        elems_1D = np.array([el for el in topo if len(el) == 2])

        # parâmetros
        R_circ = 0.15
        R_temp = 0.145      # raio do círculo
        dx = 0.01    # espaçamento
        
        x_pts = []
        y_pts = []
        
        # número de passos
        N = int(np.floor(R_temp/dx))
        
        # gera pontos internos
        for i in range(-N, N+1):
            for j in range(-N, N+1):
        
                x_temp = i * dx
                y_temp = j * dx
        
                # mantém apenas pontos internos
                if x_temp**2 +  y_temp**2 <= R_temp**2:
                    x_pts.append(x_temp)
                    y_pts.append(y_temp)
        #if len(elems_2D) > 0:
        triang = tri.Triangulation(x, y, elems_2D)
        #tpc = ax.tripcolor(triang,facecolors=sigma[:len(elems_2D)],edgecolors='k', cmap='Blues')#,vmin=1.0 )
        ntri = triang.triangles.shape[0]
        fc = sigma.ravel()[:ntri]
        
        finder = triang.get_trifinder()
        pontos = np.column_stack((x_pts, y_pts))
        print("pontos:", pontos.shape)

        ex = pontos[:,0]
        ey = pontos[:,1]
        #ex, ey = 0.02, 0.07
    
        idx_elem = finder(ex, ey) 
        idx_elem_global = finder(ex, ey) + ajuste
        #print("Elemento:", idx_elem)
        
        mask = idx_elem  >= 0

        pontos_validos = pontos[mask]
        idx_validos = idx_elem[mask]
        
        
        # pega os valores correspondentes aos elementos onde caíram os pontos
        sigma_L_pts = sigma_L[idx_validos]
        sigma_T_pts = sigma_T[idx_validos]
        theta_pts   = theta_deg[idx_validos]
        
        dados_elipses = np.column_stack([
                                        idx_validos,
                                        pontos_validos[:, 0],
                                        pontos_validos[:, 1],
                                        sigma_L_pts,
                                        sigma_T_pts,
                                        theta_pts
                                    ])

        self.quantitative_elipses = np.column_stack((dados_elipses_phantom, dados_elipses))
        np.savetxt("quantitative_elipses_artigo_2b.txt",  self.quantitative_elipses)
        fig, ax = plt.subplots(figsize=(7, 7))

        #sigmaL_max = max(np.max(np.abs(sigma_L_pts)), 1e-12)
        #sigmaT_max = max(np.max(np.abs(sigma_T_pts)), 1e-12)
        sigma_max = max(np.max(np.abs(sigma_L_pts)), np.max(np.abs(sigma_T_pts)), 1e-12)
        
        escala = 0.005
        
        
        
        # ===== índice de anisotropia para todos os pontos =====
        #AI_all = sigma_T_pts / sigma_L_pts

        ##sigma_max = max(np.max(np.abs(AI_all)), 1e-12)
        
        #AI_min = np.min(AI_all)
        #AI_max = np.max(AI_all)
        
        # evita erro se todos os valores forem iguais
        #if abs(AI_max - AI_min) < 1e-12:
        #    AI_min = AI_min - 1.0
        #    AI_max = AI_max + 1.0
        
        #norm = colors.Normalize(
        #    vmin=AI_min,
        #    vmax=AI_max
        #)
        
        
        
        norm = colors.Normalize(vmin=0.0, vmax=1.0)
        cmap = cm.jet
        
        for linha in dados_elipses:
        
            _, x0, y0, sL, sT, theta = linha
        
            AI = sT / sL
            if AI > 1.0:
                AI = 1 / AI
        
            a = escala * abs(sL) / sigma_max
            b = escala * abs(sT) / sigma_max
        
            cor = cmap(norm(AI))
        
            elipse = Ellipse(
                xy=(x0, y0),
                width=2 * a,
                height=2 * b,
                angle=theta,
                fill=True,
                facecolor=cor,
                edgecolor='black',
                linewidth=0.3,
                alpha=0.8
            )
        
            ax.add_patch(elipse)
        
        # ===== colorbar =====
        sm = cm.ScalarMappable(
            cmap=cmap,
            norm=norm
        )
        
        sm.set_array([])
        
        cbar = plt.colorbar(sm, ax=ax, shrink=0.7)
        cbar.set_label(r"Anisotropy $AI = \sigma_{min} / \sigma_{max}$", fontsize=15)
        
        theta_circ = np.linspace(0, 2*np.pi, 400)
        
        ax.plot(
            R_circ*np.cos(theta_circ),
            R_circ*np.sin(theta_circ),
            'k-',
            linewidth=1.0
        )
        
        ax.set_xlim(-R_circ - 0.02, R_circ + 0.02)
        ax.set_ylim(-R_circ - 0.02, R_circ + 0.02)
        
        ax.set_title(
            #f"σ Anisotropy - λ_{Lambda:.2e}-it_{iteration} - Aniso_{DifAniso:.1f}",
            f"σ Anisotropy",
            fontsize=20
        )
        
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('Lenght [m]', fontsize = 15)
        ax.set_ylabel('Lenght [m]', fontsize = 15)
        ax.grid(False)
        
        #plt.show()
        
                
        if save:
            #plt.savefig(nome_arquivo, dpi=100)
            #plt.savefig(nome_arquivo, dpi=200,bbox_inches="tight")
            #plt.savefig(f'{nome_arquivo}.svg',   format="svg", dpi=300) 
            plt.savefig(f'{nome_arquivo}.{imgFormat}',   format=imgFormat, dpi=300) 
            #plt.show()
            plt.show(block=False)
            plt.pause(0.1)  
            plt.close()   # importante
        else:
            plt.show()    
        
    '''
    ###############################################################################   
    def plotElipse(self, sigma, sigma_L, sigma_T, theta_deg, Lambda = None, iteration = None, save = False, DifAniso=None, nome_arquivo= None):
        dados_elipses_phantom = np.loadtxt("dados_elipses_phantom.txt")
        #print(f'dados_elipses_phantom\n {dados_elipses_phantom}')
        
        x, y = self.mymesh.Coordinates[:, 0], self.mymesh.Coordinates[:, 1]
        topo = self.mymesh.msh_topology
        ajuste  = self.mymesh.NumberOfElectrodes*4 +1
        imgFormat = self.imgFotmat
    
        # --- Separar elementos 2D (triangulares) e 1D (linhas) ---
        elems_2D = np.array([el for el in topo if len(el) == 3])
        elems_1D = np.array([el for el in topo if len(el) == 2])

        # parâmetros
        R_circ = 0.15
        R_temp = 0.145      # raio do círculo
        dx = 0.01    # espaçamento
        
        x_pts = []
        y_pts = []
        
        # número de passos
        N = int(np.floor(R_temp/dx))
        
        # gera pontos internos
        for i in range(-N, N+1):
            for j in range(-N, N+1):
        
                x_temp = i * dx
                y_temp = j * dx
        
                # mantém apenas pontos internos
                if x_temp**2 +  y_temp**2 <= R_temp**2:
                    x_pts.append(x_temp)
                    y_pts.append(y_temp)
        #if len(elems_2D) > 0:
        triang = tri.Triangulation(x, y, elems_2D)
        #tpc = ax.tripcolor(triang,facecolors=sigma[:len(elems_2D)],edgecolors='k', cmap='Blues')#,vmin=1.0 )
        ntri = triang.triangles.shape[0]
        fc = sigma.ravel()[:ntri]
        
        finder = triang.get_trifinder()
        pontos = np.column_stack((x_pts, y_pts))
        print("pontos:", pontos.shape)

        ex = pontos[:,0]
        ey = pontos[:,1]
        #ex, ey = 0.02, 0.07
    
        idx_elem = finder(ex, ey) 
        idx_elem_global = finder(ex, ey) + ajuste
        #print("Elemento:", idx_elem)
        
        mask = idx_elem  >= 0

        pontos_validos = pontos[mask]
        idx_validos = idx_elem[mask]
        
        
        # pega os valores correspondentes aos elementos onde caíram os pontos
        sigma_L_pts = sigma_L[idx_validos]
        sigma_T_pts = sigma_T[idx_validos]
        theta_pts   = theta_deg[idx_validos]
        #HPF_pts = HPF[idx_validos]
        HPF_pts = (sigma_T[idx_validos] / sigma_L[idx_validos])/20
        
        dados_elipses = np.column_stack([
                                        idx_validos,
                                        pontos_validos[:, 0],
                                        pontos_validos[:, 1],
                                        sigma_L_pts,
                                        sigma_T_pts,
                                        theta_pts,
                                        HPF_pts
                                    ])

        self.quantitative_elipses = np.column_stack((dados_elipses_phantom, dados_elipses))
        np.savetxt("quantitative_elipses_artigo_2b.txt",  self.quantitative_elipses)
        
        #fig, ax = plt.subplots(figsize=(7, 7))
        fig = plt.figure(figsize=(10, 9))
        ax = fig.add_subplot(111, projection='3d')

        #sigmaL_max = max(np.max(np.abs(sigma_L_pts)), 1e-12)
        #sigmaT_max = max(np.max(np.abs(sigma_T_pts)), 1e-12)
        sigma_max = max(np.max(np.abs(sigma_L_pts)), np.max(np.abs(sigma_T_pts)), 1e-12)
        
        escala = 0.005
    '''
    '''
        # ===== índice de anisotropia para todos os pontos =====
        AI_all = sigma_T_pts / sigma_L_pts

        #sigma_max = max(np.max(np.abs(AI_all)), 1e-12)
        
        AI_min = np.min(AI_all)
        AI_max = np.max(AI_all)
        
        # evita erro se todos os valores forem iguais
        if abs(AI_max - AI_min) < 1e-12:
            AI_min = AI_min - 1.0
            AI_max = AI_max + 1.0
        
        norm = colors.Normalize(
            vmin=AI_min,
            vmax=AI_max
        )
    '''
    '''
        norm = colors.Normalize(vmin=0.0, vmax=1.0)
        cmap = cm.jet
    '''
    '''
        for linha in dados_elipses:
            
            #_, x0, y0, sL, sT, theta = linha
            _, x0, y0, sL, sT, theta, z0 = linha
        
            AI = sT / sL
            if AI > 1.0:
                AI = 1 / AI
        
            a = escala * abs(sL) / sigma_max
            b = escala * abs(sT) / sigma_max
        
            cor = cmap(norm(AI))
    '''
    '''
        for linha in dados_elipses:
        
            _, x0, y0, sL, sT, theta, z0 = linha
        
            # ========================================================
            # Índice de anisotropia
            # ========================================================
        
            AI = sT / sL
        
            if AI > 1.0:
                AI = 1 / AI
        
            # ========================================================
            # Dimensões da elipse
            # ========================================================
        
            a = escala * abs(sL) / sigma_max
            b = escala * abs(sT) / sigma_max
        
            # ========================================================
            # Cor
            # ========================================================
        
            cor = cmap(norm(AI))
        
            # ========================================================
            # Coordenadas da elipse
            # ========================================================
        
            phi = np.linspace(
                0,
                2*np.pi,
                50
            )
        
            theta_rad = np.deg2rad(theta)
        
            # elipse sem rotação
            xe = a * np.cos(phi)
            ye = b * np.sin(phi)
        
            # rotação
            xr = (
                xe * np.cos(theta_rad)
                - ye * np.sin(theta_rad)
            )
        
            yr = (
                xe * np.sin(theta_rad)
                + ye * np.cos(theta_rad)
            )
        
            # translação para x0,y0
            Xelipse = x0 + xr
            Yelipse = y0 + yr
        
            # altura = efeito do filtro
            Zelipse = np.full_like(
                Xelipse,
                z0
            )
        
            # ========================================================
            # Cria polígono 3D
            # ========================================================
        
            vertices = [
                list(
                    zip(
                        Xelipse,
                        Yelipse,
                        Zelipse
                    )
                )
            ]
        
            poly = Poly3DCollection(
                vertices,
                facecolor=cor,
                edgecolor='black',
                linewidth=0.3,
                alpha=0.8
            )
        
            ax.add_collection3d(poly)
        # ===== colorbar =====
        sm = cm.ScalarMappable(
            cmap=cmap,
            norm=norm
        )
        
        #sm.set_array([])
        
        #cbar = plt.colorbar(sm, ax=ax, shrink=0.7, pad=0.10)
        #cbar.set_label(r"Anisotropy $AI = \sigma_{min} / \sigma_{max}$", fontsize=15)
        
        theta_circ = np.linspace(0, 2*np.pi, 400)
        
        ax.plot(
            R_circ*np.cos(theta_circ),
            R_circ*np.sin(theta_circ),
            np.zeros_like(theta_circ),
            'k-',
            linewidth=1.0
        )
        
        ax.set_xlim(-R_circ - 0.02, R_circ + 0.02)
        ax.set_ylim(-R_circ - 0.02, R_circ + 0.02)
        ax.set_zlim(0,0.1)
        ax.set_title(
            #f"σ Anisotropy - λ_{Lambda:.2e}-it_{iteration} - Aniso_{DifAniso:.1f}",
            f"σ Anisotropy",
            fontsize=20
        )
        
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('Lenght [m]', fontsize = 15)
        ax.set_ylabel('Lenght [m]', fontsize = 15)
        # ============================================================
        # Ângulo da câmera
        # ============================================================

        ax.view_init(elev=30, azim=-60 )
        #ax.set_box_aspect((1, 1, 0.4), zoom=2.0)
        ax.grid(False)
        
        #plt.show()
        
                
        if save:
            #plt.savefig(nome_arquivo, dpi=100)
            #plt.savefig(nome_arquivo, dpi=200,bbox_inches="tight")
            #plt.savefig(f'{nome_arquivo}.svg',   format="svg", dpi=300) 
            plt.savefig(f'{nome_arquivo}.{imgFormat}',   format=imgFormat, dpi=300) 
            #plt.show()
            plt.tight_layout()
            plt.show(block=False)
            plt.pause(0.1)  
            plt.close()   # importante
        else:
            plt.show()
    '''
    ###############################################################################

    
    ###############################################################################
    #          Plotar máscara de anisotropia
    
    ###############################################################################


    '''
    def plotMask(self, sigma, Lambda = None, iteration = None, save = False, SigmaXXXYYY = None, DifAniso = None, nome_arquivo= None):

        x, y = self.mymesh.Coordinates[:, 0], self.mymesh.Coordinates[:, 1]
        topo = self.mymesh.msh_topology
    
        # --- Separar elementos 2D (triangulares) e 1D (linhas) ---
        elems_2D = np.array([el for el in topo if len(el) == 3])
        elems_1D = np.array([el for el in topo if len(el) == 2])
    
        fig, ax = plt.subplots(figsize=(6, 6.2))
        #fig, ax = plt.subplots(figsize=(5,5))
        #plt.figure(figsize=(6, 6))
        #ax.set_aspect('equal')
        #ax.set_aspect('equal', adjustable='box')

        # ============================================================
        #   PLOTAR ELEMENTOS 2D (TRIANGULARES)
        # ============================================================
        if len(elems_2D) > 0:
            triang = tri.Triangulation(x, y, elems_2D)
            #tpc = ax.tripcolor(triang,facecolors=sigma[:len(elems_2D)],edgecolors='k', cmap='Blues')#,vmin=1.0 )
            ntri = triang.triangles.shape[0]
            fc = sigma.ravel()[:ntri]
            #level = np.max(np.abs(fc))
            #threshold = 0.7*level
            threshold = np.percentile(np.abs(fc), 80)

            # máscara binária
            fc_bin = np.where(fc >= threshold, 1, 0)

            
            # mapa de cores fixo
            cmap = ListedColormap(['white', 'black'])

            # força apenas 2 níveis
            normMask = BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)
            
            tpc = ax.tripcolor(triang,facecolors=fc_bin, edgecolors='k',cmap=cmap, norm=normMask )
            
            #cbar = plt.colorbar(tpc, ticks=[0,1])
            #cbar.set_label('Máscara binária')
            #cbar.set_ticklabels(['< 1.0', '≥ 1.0'])
            #fig.colorbar(tpc, ax=ax, label='Conductivity σ [S/m]')
            #plt.colorbar(tpc, label='Máscara binária')
            #ax.set_aspect('equal')
            if save == True:
                timestamp = datetime.now().strftime("%m%d_%H%M")
                ax.set_title(f"σ{SigmaXXXYYY} - λ_{Lambda:.2e}-it_{iteration} - Aniso_{DifAniso:.1f} - th_{threshold:.3f}", fontsize=11)
            if save == False:
                ax.set_title(f"Conductivity Real (σ) ", fontsize=15)

        if save == True:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M")
            #plt.savefig(f"Conductivity_itr_{iteration}.png", dpi=150, bbox_inches='tight')
            plt.savefig(f'{nome_arquivo}.svg', format="svg", dpi=200)
            #plt.savefig(f'{nome_arquivo}', dpi=300, bbox_inches='tight')
            #plt.savefig(f'{nome_arquivo}',  dpi=150, bbox_inches='tight', pil_kwargs={"quality": 70})
        plt.show(block=False)   # mostra sem travar
        plt.pause(0.1)            # mantém aberto por 3 segundos
        plt.close('all')        # fecha automaticamente
    ###############################################################################
    '''





    ###############################################################################
    # Essa função plota o resultado do sigma calculado no problema inverso
    ############################################################################### 

    def plot_sigma(self, sigmaResult, ref_sL = 3.0, ref_sT = 1.0, titulo="Results of  σxx, σxy, σyy, σL and σT", salvar=False, nome_arquivo="plot_sigma.png"):
    
        data = sigmaResult
        imgFormat = self.imgFotmat
    
        sigma_xx = data[:, 0]
        sigma_xy = data[:, 1]
        sigma_yy = data[:, 2]
        #sigma_Dif = sigma_xx - sigma_yy

        
        Smed = 0.5 * (sigma_xx + sigma_yy)
        D = np.sqrt(((sigma_xx - sigma_yy)/2)**2 + sigma_xy**2)

        sigma_L = Smed + D
        sigma_T = Smed - D

        sigma_Dif = sigma_L - sigma_T

        x = np.arange(len(sigma_xx))
        plot_ref_sL = np.ones(len(sigma_xx))*ref_sL
        plot_ref_sT = np.ones(len(sigma_xx))*ref_sT
        fig, ax = plt.subplots(figsize=(6.0, 5.0))

        ax.plot(x, sigma_xx, label='σxx', linewidth=1.0)
        ax.plot(x, sigma_xy, label='σxy', linewidth=1.0)
        ax.plot(x, sigma_yy, label='σyy', linewidth=1.0)
        ax.plot(x, sigma_L, label='σL', linewidth=2.0, linestyle=':')
        ax.plot(x, sigma_T, label='σT', linewidth=2.0, linestyle=':')
        ax.plot(x, plot_ref_sL, label='Ref_σL', linewidth=3.0, linestyle='--')
        ax.plot(x, plot_ref_sT, label='Ref_σT', linewidth=3.0, linestyle='--')
        
        ax.set_title(titulo, fontsize=14)
        ax.set_xlabel('Element', fontsize=10)
        ax.set_ylabel('Conductivity [S/m]', fontsize=10)
        
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=4, fontsize=10)
        
        fig.tight_layout()

        
        if salvar:
            #plt.savefig(nome_arquivo, dpi=100)
            #plt.savefig(nome_arquivo, dpi=200,bbox_inches="tight")
            #plt.savefig(f'{nome_arquivo}.svg',  format="svg", dpi=200)
            plt.savefig(f'{nome_arquivo}.{imgFormat}',  format= imgFormat, dpi=200)
            #plt.show()
            plt.show(block=False)
            plt.pause(0.1)  
            plt.close()   # importante
        else:
            plt.show()    
    ############################################################################### 





    '''
    ###############################################################################
    # Essa função plota o resultado do sigma calculado no problema inverso
    ############################################################################### 

    def plot_Sorting(self, sigma_Dif, titulo="Results (σx-σy)", salvar=False, nome_arquivo="plot_Sorting.svg"):
        
        #sigma_dif = sigma_x - sigma_y
        
        # ordenar do maior para o menor
        sigma_dif_ord = np.sort(sigma_Dif)[::-1]
        
        # eixo x
        x = np.arange(1, len(sigma_dif_ord) + 1)
        
        # plot
        plt.figure(figsize=(6, 6))
        plt.plot(x, sigma_dif_ord, linewidth=1.5)
        
        plt.title(r'($\sigma_x - \sigma_y$) Sorting')
        plt.xlabel('Sorting position')
        plt.ylabel(r'$\sigma_x - \sigma_y$')
        
        plt.grid(True, linestyle='--', alpha=0.5)
    
    
        plt.tight_layout()
 
        
        if salvar:
            #plt.savefig(nome_arquivo, dpi=100)
            #plt.savefig(nome_arquivo, dpi=200,bbox_inches="tight")
            plt.savefig(f'{nome_arquivo}.svg', format="svg", dpi=200)
            #plt.show()
            plt.show(block=False)
            plt.pause(0.1)  
            plt.close()   # importante
        else:
            plt.show()    
    ###############################################################################    
    '''







       
            
    def plot_theta_deg(self, theta_deg, titulo="Results of the calculated θ angle ", salvar=False, nome_arquivo="plot_theta.png"):
    
        data_theta = theta_deg
        imgFormat = self.imgFotmat
    
        x = np.arange(len(data_theta))
    
        plt.figure(figsize=(6, 6))
    

        plt.plot(x, data_theta, label='θ', linewidth=1.5)
    
        plt.title(titulo, fontsize=14)
        plt.xlabel('Element', fontsize=12)
        plt.ylabel('Angle θ [°]', fontsize=12)
        plt.ylim(-90, 90)   # exemplo: de -65° a 66°
        
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend()
    
        plt.tight_layout()
        
        
        if salvar:
            #plt.savefig(nome_arquivo, dpi=150)
            #plt.savefig(f'{nome_arquivo}.svg', format="svg", dpi=200)
            plt.savefig(f'{nome_arquivo}.{imgFormat}', format=imgFormat, dpi=200)
            #plt.show()
            plt.show(block=False)
            plt.pause(0.1)  
            plt.close()   # importante
        else:
            plt.show()    
    ###############################################################################
    '''    
    def plot_espectro(self,x, titulo="Espectro (FFT)"):
        X = np.fft.fft(x)
        w = np.fft.fftfreq(len(x))
        
        X_shifted = np.fft.fftshift(X)
        w_shifted = np.fft.fftshift(w)
        
        



        plt.figure(figsize=(8,4))
        #plt.stem(w, np.abs(X))
        #plt.plot(w, np.abs(X))
        plt.plot(w_shifted, np.abs(X_shifted))
        plt.yscale("log")
        plt.xlim(left=0)
        plt.ylim(top=1e4) 
        #plt.ylim(bottom=1e-3) 
        plt.title(titulo)
        plt.xlabel("Frequência (w)")
        plt.ylabel("|X(w)|")
        plt.grid(True)
        plt.tight_layout()
        #plt.show()
        plt.show(block=False)
        plt.pause(0.1)      
    '''

        
    def configurar_celula(self,
            celula,
            fontsize=10,
            bold=False,
            linewidth=0.8,
            facecolor="white"
    ):
        """Aplica a formatação básica a uma célula da tabela."""
    
        celula.set_edgecolor("black")
        celula.set_linewidth(linewidth)
        celula.set_facecolor(facecolor)
    
        celula.get_text().set_fontsize(fontsize)
        celula.get_text().set_weight("bold" if bold else "normal")
        celula.get_text().set_ha("center")
        celula.get_text().set_va("center")
        
    
    def plotar_tabela_condutividade(self, ax, background, object_1, object_2, object_3):
        """Desenha a TABLE I."""
        #print('background', background)
        Z11 = f"{np.mean(background[:, 3]):.2e}"
        Z12 = f"{np.mean(background[:, 4]):.2e}"
        Z13 = f"{np.mean(background[:,10]):.2e}"
        Z14 = f"{np.mean(background[:,11]):.2e}"
        Z15 = f"{(np.mean(background[:,3]) - np.mean(background[:,10])):.2e}"
        Z16 = f"{(np.mean(background[:,4]) - np.mean(background[:,11])):.2e}"
        
        Z21 = f"{np.mean(object_1[:, 3]):.2e}"
        Z22 = f"{np.mean(object_1[:, 4]):.2e}"
        Z23 = f"{np.mean(object_1[:,10]):.2e}"
        Z24 = f"{np.mean(object_1[:,11]):.2e}"
        Z25 = f"{(np.mean(object_1[:,3]) - np.mean(object_1[:,10])):.2e}"
        Z26 = f"{(np.mean(object_1[:,4]) - np.mean(object_1[:,11])):.2e}"
        
        Z31 = f"{np.mean(object_2[:, 3]):.2e}"
        Z32 = f"{np.mean(object_2[:, 4]):.2e}"
        Z33 = f"{np.mean(object_2[:,10]):.2e}"
        Z34 = f"{np.mean(object_2[:,11]):.2e}"
        Z35 = f"{(np.mean(object_2[:,3]) - np.mean(object_2[:,10])):.2e}"
        Z36 = f"{(np.mean(object_2[:,4]) - np.mean(object_2[:,11])):.2e}"
        
        Z41 = f"{np.mean(object_3[:, 3]):.2e}"
        Z42 = f"{np.mean(object_3[:, 4]):.2e}"
        Z43 = f"{np.mean(object_3[:,10]):.2e}"
        Z44 = f"{np.mean(object_3[:,11]):.2e}"
        Z45 = f"{(np.mean(object_3[:,3]) - np.mean(object_3[:,10])):.2e}"
        Z46 = f"{(np.mean(object_3[:,4]) - np.mean(object_3[:,11])):.2e}"
        
        Z51 = f"{np.mean(self.quantitative_elipses[:, 3]):.2e}"
        Z52 = f"{np.mean(self.quantitative_elipses[:, 4]):.2e}"
        Z53 = f"{np.mean(self.quantitative_elipses[:,10]):.2e}"
        Z54 = f"{np.mean(self.quantitative_elipses[:,11]):.2e}"
        Z55 = f"{(np.mean(self.quantitative_elipses[:,3]) - np.mean(self.quantitative_elipses[:,10])):.2e}"
        Z56 = f"{(np.mean(self.quantitative_elipses[:,4]) - np.mean(self.quantitative_elipses[:,11])):.2e}"
        
        
        
        ax.axis("off")
    
        # Títulos
        ax.text(
            0.5,
            0.97,
            "TABLE I",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=11,
            fontfamily="serif"
        )
    
        ax.text(
            0.5,
            0.89,
            "REFERENCE AND RECONSTRUCTED CONDUCTIVITY PARAMETERS",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=10,
            fontfamily="serif",
            fontvariant="small-caps"
        )
    
        # Duas linhas de cabeçalho
        cabecalho_1 = [
            "",
            "Reference", "",
            "Results", "",
            "Mean Error", ""
        ]
    
        cabecalho_2 = [
            "Region",
            r"${\overline{\sigma L}}_{ref}$",
            r"$\overline{\sigma T}_{ref}$",
            r"$\overline{\sigma L}_{calc}$",
            r"$\overline{\sigma T}_{calc}$",
            r"$\overline{\sigma L}_{error}$",
            r"$\overline{\sigma T}_{error}$"
        ]
    
        dados = [
            ["Bkgd", Z11, Z12, Z13, Z14, Z15, Z16],
            ["Obj 1", Z21, Z22, Z23, Z24, Z25, Z26],
            ["Obj 2", Z31, Z32, Z33, Z34, Z35, Z36],
            ["Obj 3", Z41, Z42, Z43, Z44, Z45, Z46],
            ["All",  Z51, Z52, Z53, Z54, Z55, Z56]
        ]
    
        conteudo = [cabecalho_1, cabecalho_2] + dados
    
        tabela = ax.table(
            cellText=conteudo,
            cellLoc="center",
            #bbox=[0.16, 0.20, 0.68, 0.56],
            #colWidths=[0.16, 0.14, 0.14, 0.14, 0.14, 0.16, 0.16]
            bbox=[0.09, 0.08, 0.82, 0.61],
            colWidths=[0.14, 0.14, 0.14, 0.15, 0.15, 0.16, 0.16]
        )
    
        tabela.auto_set_font_size(False)
    
        numero_linhas = len(conteudo)
        numero_colunas = len(conteudo[0])
    
        for linha in range(numero_linhas):
            for coluna in range(numero_colunas):
    
                celula = tabela[(linha, coluna)]
    
                self.configurar_celula(
                    celula,
                    fontsize=9,
                    bold=(linha <= 1),
                    linewidth=0.8
                )
    
                celula.get_text().set_fontfamily("serif")
    
        # Ajusta a altura das linhas
        for coluna in range(numero_colunas):
            tabela[(0, coluna)].set_height(0.12)
            tabela[(1, coluna)].set_height(0.14)
    
        for linha in range(2, numero_linhas):
            for coluna in range(numero_colunas):
                tabela[(linha, coluna)].set_height(0.105)
    
        # Simulação das células mescladas
        # Reference: colunas 1 e 2
        tabela[(0, 1)].visible_edges = "LTB"
        tabela[(0, 2)].visible_edges = "RTB"
    
        # Results: colunas 3 e 4
        tabela[(0, 3)].visible_edges = "LTB"
        tabela[(0, 4)].visible_edges = "RTB"
    
        # Mean Error: colunas 5 e 6
        tabela[(0, 5)].visible_edges = "LTB"
        tabela[(0, 6)].visible_edges = "RTB"
    
        # Célula superior de Region
        tabela[(0, 0)].visible_edges = "LTR"
        tabela[(1, 0)].visible_edges = "LBR"
    
        # Remove o texto da primeira célula
        tabela[(0, 0)].get_text().set_text("")
    
        # Nota
        ax.text(
            0.18,
            0.04,
            r"Note: Conductivities $\sigma$ in S/m.",
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=9,
            fontfamily="serif"
        )
    
    
    def plotar_tabela_anisotropia(self, ax, background, object_1, object_2, object_3):
        """Desenha a TABLE II."""

        #Acalc = (1 - background[:, 11]/background[:, 10])
        #S2 = 1 - Acalc
       
        #S2a = np.sqrt((1 - background[:, 11]/background[:, 10])*(background[:, 12] - background[:, 5])**2)  
        #print('S2a',S2a)
        
        ZZ11 = f"{np.mean(background[:, 4]/background[:, 3]):.2e}"
        ZZ12 = f"{np.mean(background[:, 5]):.2e}°"
        ZZ13 =  f"{np.mean(background[:, 11]/background[:, 10]):.2e}"
        ZZ14 = f"{np.mean(background[:,12]):.2e}°"
        ZZ15 = f"{((np.mean(background[:, 11]/background[:, 10]))-(np.mean(background[:, 4]/background[:, 3]))):.2e}"
        ZZ16 = f"{np.sqrt(np.mean(((1 - background[:, 11] / background[:, 10])*(background[:, 12] - background[:, 5]))**2)):.2e}°"
        
        ZZ21 = f"{np.mean(object_1[:, 4]/object_1[:, 3]):.2e}"
        ZZ22 = f"{np.mean(object_1[:, 5]):.2e}°"
        ZZ23 =  f"{np.mean(object_1[:, 11]/object_1[:, 10]):.2e}"
        ZZ24 = f"{np.mean(object_1[:,12]):.2e}°"
        ZZ25 =f"{((np.mean(object_1[:, 11]/object_1[:, 10]))-(np.mean(object_1[:, 4]/object_1[:, 3]))):.2e}"
        ZZ26 = f"{np.sqrt(np.mean(((1 - object_1[:, 11] / object_1[:, 10])*(object_1[:, 12] - object_1[:, 5]))**2)):.2e}°"
        
        ZZ31 = f"{np.mean(object_2[:, 4]/object_2[:, 3]):.2e}"
        ZZ32 = f"{np.mean(object_2[:, 5]):.2e}°"
        ZZ33 =  f"{np.mean(object_2[:, 11]/object_2[:, 10]):.2e}"
        ZZ34 = f"{np.mean(object_2[:,12]):.2e}°"
        ZZ35 = f"{((np.mean(object_2[:, 11]/object_2[:, 10]))-(np.mean(object_2[:, 4]/object_2[:, 3]))):.2e}"
        ZZ36 = f"{np.sqrt(np.mean(((1 - object_2[:, 11] / object_2[:, 10])*(object_2[:, 12] - object_2[:, 5]))**2)):.2e}°"
        
        ZZ41 = f"{np.mean(object_3[:, 4]/object_3[:, 3]):.2e}"
        ZZ42 = f"{np.mean(object_3[:, 5]):.2e}°"
        ZZ43 =  f"{np.mean(object_3[:, 11]/object_3[:, 10]):.2e}"
        ZZ44 = f"{np.mean(object_3[:,12]):.2e}°"
        ZZ45 =f"{((np.mean(object_3[:, 11]/object_3[:, 10]))-(np.mean(object_3[:, 4]/object_3[:, 3]))):.2e}"
        ZZ46 = f"{np.sqrt(np.mean(((1 - object_3[:, 11] / object_3[:, 10])*(object_3[:, 12] - object_3[:, 5]))**2)):.2e}°"
        
        ZZ51 = f"{np.mean(self.quantitative_elipses[:, 4]/self.quantitative_elipses[:, 3]):.2e}"
        ZZ52 = f"{np.mean(self.quantitative_elipses[:, 5]):.2e}°"
        ZZ53 =  f"{np.mean(self.quantitative_elipses[:, 11]/self.quantitative_elipses[:, 10]):.2e}"
        ZZ54 = f"{np.mean(self.quantitative_elipses[:,12]):.2e}°"
        ZZ55 = f"{((np.mean(self.quantitative_elipses[:, 11]/self.quantitative_elipses[:, 10]))-(np.mean(self.quantitative_elipses[:, 4]/self.quantitative_elipses[:, 3]))):.2e}"
        ZZ56 = f"{np.sqrt(np.mean(((1 - self.quantitative_elipses[:, 11] / self.quantitative_elipses[:, 10])*(self.quantitative_elipses[:, 12] - self.quantitative_elipses[:, 5]))**2)):.2e}°"
        
        
    
        ax.axis("off")
    
        # Títulos
        ax.text(
            0.5,
            0.97,
            "TABLE II",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=11,
            fontfamily="serif"
        )
    
        ax.text(
            0.5,
            0.89,
            "REFERENCE AND RECONSTRUCTED ANISOTROPY PARAMETERS",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=10,
            fontfamily="serif",
            fontvariant="small-caps"
        )
    
        cabecalho_1 = [
            "",
            "Reference", "",
            "Results", "",
            "Mean Error", ""
        ]
    
        cabecalho_2 = [
            "Region",
            r"$\overline{AI}_{ref}$",
            r"$\overline{\theta}_{ref}$",
            r"$\overline{AI}_{calc}$",
            r"$\overline{\theta}_{calc}$",
            r"$\overline{AI}_{error}$",
            r"$\overline{\theta}_{error}$"
        ]
    
        dados = [
            ["Bkgd", ZZ11, ZZ12, ZZ13, ZZ14, ZZ15, ZZ16],
            ["Obj 1", ZZ21, ZZ22, ZZ23, ZZ24, ZZ25, ZZ26],
            ["Obj 2", ZZ31, ZZ32, ZZ33, ZZ34, ZZ35, ZZ36],
            ["Obj 3", ZZ41, ZZ42, ZZ43, ZZ44, ZZ45, ZZ46],
            ["All",  ZZ51, ZZ52, ZZ53, ZZ54, ZZ55, ZZ56]
        ]
    
        conteudo = [cabecalho_1, cabecalho_2] + dados
    
        tabela = ax.table(
            cellText=conteudo,
            cellLoc="center",
            bbox=[0.09, 0.08, 0.82, 0.61],
            colWidths=[0.14, 0.14, 0.14, 0.15, 0.15, 0.16, 0.16]
        )
    
        tabela.auto_set_font_size(False)
    
        numero_linhas = len(conteudo)
        numero_colunas = len(conteudo[0])
    
        for linha in range(numero_linhas):
            for coluna in range(numero_colunas):
    
                celula = tabela[(linha, coluna)]
    
                self.configurar_celula(
                    celula,
                    fontsize=9,
                    bold=(linha <= 1),
                    linewidth=0.8
                )
    
                celula.get_text().set_fontfamily("serif")
    
        # Altura das linhas
        for coluna in range(numero_colunas):
            tabela[(0, coluna)].set_height(0.12)
            tabela[(1, coluna)].set_height(0.15)
    
        for linha in range(2, numero_linhas):
            for coluna in range(numero_colunas):
                tabela[(linha, coluna)].set_height(0.105)
    
        # Simula células agrupadas
        tabela[(0, 1)].visible_edges = "LTB"
        tabela[(0, 2)].visible_edges = "RTB"
    
        tabela[(0, 3)].visible_edges = "LTB"
        tabela[(0, 4)].visible_edges = "RTB"
    
        tabela[(0, 5)].visible_edges = "LTB"
        tabela[(0, 6)].visible_edges = "RTB"
    
        # Region ocupando as duas linhas de cabeçalho
        tabela[(0, 0)].visible_edges = "LTR"
        tabela[(1, 0)].visible_edges = "LBR"
        tabela[(0, 0)].get_text().set_text("")


    ###############################################################################
    # Essa função plota o resultado do sigma calculado no problema inverso
    ############################################################################### 
    
    def plot_Table(self, nome_arquivo="plot_theta.png"):
        
        background = self.quantitative_elipses[self.quantitative_elipses[:, 6] == 1000]
        object_1 = self.quantitative_elipses[self.quantitative_elipses[:, 6] == 1001]
        object_2 = self.quantitative_elipses[self.quantitative_elipses[:, 6] == 1002]
        object_3 = self.quantitative_elipses[self.quantitative_elipses[:, 6] == 1003]
        
        imgFormat = self.imgFotmat
   

    # ============================================================
    # Figura principal
    # ============================================================
    
        fig, eixos = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(7.3, 6.2),
            gridspec_kw={
                "height_ratios": [1.0, 1.08],
                "hspace": 0.05
            }
        )
    
        self.plotar_tabela_condutividade(eixos[0], background, object_1, object_2, object_3)
        self.plotar_tabela_anisotropia(eixos[1], background, object_1, object_2, object_3)
    
        # Ajustes das margens
        plt.subplots_adjust(
            left=0.03,
            right=0.97,
            top=0.98,
            bottom=0.03,
            hspace=0.02
        )
    
        # Salvar a imagem
        #plt.savefig("Tabelas_IEEE.png", dpi=300, bbox_inches="tight", pad_inches=0.05)
        plt.savefig(f'{nome_arquivo}.{imgFormat}',   format = imgFormat, dpi=300)
        # Também pode salvar em formato vetorial:
        # plt.savefig("Tabelas_IEEE.svg", bbox_inches="tight")
        # plt.savefig("Tabelas_IEEE.pdf", bbox_inches="tight")
    
        plt.show()
    ###############################################################################
    '''
    # Essa função salva arqui html do resultado doproblema inverso
    ###############################################################################         
    @staticmethod
    def salvar_html(self, lista_imagens, nome_html="resultado"):
        with open(nome_html, "w", encoding="utf-8") as f:
            f.write("""
            <html>
            <head>
            <title>Resultados TIE</title>
            </head>
            <body>
            """)
    
            f.write("<h1>Results of the Anisotropic Inverse Problem</h1>\n")
    
            f.write('<div style="display:flex; gap:10px;">\n')
            f.write("<h1>nome_html</h1>\n")
            for img in lista_imagens:
                #f.write(f'<img src="{img}" width="150">\n')
                f.write(f'<img src="{img}" style="width:150px; height: auto;">\n')
                
                
    
            f.write('</div>\n')
            

            f.write("</body>\n</html>")
    '''
    ###############################################################################
    # Essa função salva arqui html do resultado doproblema inverso
    ###############################################################################         
        
    def salvar_html_todos_lambdas(self, pasta, html_name="resultado_completo"):
    
        print('pasta_salvar_html_todos_lambdas', pasta)
        
        imgFormat = self.imgFotmat
    
        #arquivos = glob.glob(os.path.join(pasta, f"{html_name}*.svg"))
        arquivos = glob.glob(os.path.join(pasta, f"{html_name}*.{imgFormat}"))
    
        grupos = defaultdict(lambda: {
            "AutoScale": {},
            "SameScale": {}
        })
    
        for arq in arquivos:
    
            nome = os.path.basename(arq)
            partes = nome.replace(f".{imgFormat}", "").split("_")
    
            if partes[-1] in ("AutoScale", "SameScale"):
                tipo_escala = partes[-1]
                lambda_str = partes[-2]
            else:
                tipo_escala = "AutoScale"
                lambda_str = partes[-1]
    
            if "sigma_xx" in nome:
                grupos[lambda_str][tipo_escala]["sigma_xx"] = nome
            elif "sigma_xy" in nome:
                grupos[lambda_str][tipo_escala]["sigma_xy"] = nome
            elif "sigma_yy" in nome:
                grupos[lambda_str][tipo_escala]["sigma_yy"] = nome
            elif "sigma_L" in nome:
                grupos[lambda_str][tipo_escala]["sigma_L"] = nome
            elif "sigma_T" in nome:
                grupos[lambda_str][tipo_escala]["sigma_T"] = nome
            #elif "sigma_Dif" in nome:
                #grupos[lambda_str][tipo_escala]["sigma_Dif"] = nome
            elif "theta" in nome:
                grupos[lambda_str][tipo_escala]["theta"] = nome
            elif "Elipses" in nome:
                grupos[lambda_str][tipo_escala]["Elipses"] = nome
            elif "sigma_linhas" in nome:
                grupos[lambda_str][tipo_escala]["sigma_linhas"] = nome
            elif "iterations" in nome:
                grupos[lambda_str][tipo_escala]["iterations"] = nome
            elif "plot_theta_deg" in nome:
                grupos[lambda_str][tipo_escala]["Sorting"] = nome
            elif "Table" in nome:
                 grupos[lambda_str][tipo_escala]["Table"] = nome
    
        colunas = [
            ("sigma_xx", "σxx"),
            ("sigma_xy", "σxy"),
            ("sigma_yy", "σyy"),
            ("sigma_L", "σ_L"),
            ("sigma_T", "σ_T"),
            #("sigma_Dif", "(σ_L-σ_T)"),
            ("theta", "Angle [°]"),
            ("Elipses", "Anisotropy"),
            ("sigma_linhas", "Summary"),
            ("iterations", "Optimization"),
            ("Table", "Results"),
        ]
    
        html_path = os.path.join(pasta, f"{html_name}.html")
    
        with open(html_path, "w", encoding="utf-8") as f:
    
            f.write("""
    <html>
    <head>
    <meta charset="utf-8">
    <title>Resultados TIE</title>
    
    <style>
    body {
        background-color: black;
        color: white;
        font-family: Arial, sans-serif;
    }
    
    .info {
        color: lime;
        font-weight: bold;
        margin-bottom: 6px;
    }
    
    table {
        border-collapse: collapse;
        width: 100%;
        margin-bottom: 35px;
    }
    
    th {
        font-size: 24px;
        padding: 10px;
        color: white;
    }
    
    td {
        padding: 2px;
        margin: 0px;
        text-align: center;
        vertical-align: top;
    }
    
    img {
        width: 200px;
        height: auto;
        background-color: white;
    }
    
    .lambda {
        color: yellow;
        font-size: 12px;
        margin-bottom: 5px;
    }
    
    .titulo {
        color: white;
        font-size: 18px;
        font-weight: bold;
        margin-bottom: 3px;
    }
    
    .escala {
        color: cyan;
        font-size: 13px;
        font-weight: bold;
        margin-bottom: 4px;
    }
    
    .lambda_bloco {
        color: orange;
        font-size: 22px;
        font-weight: bold;
        padding-top: 18px;
        padding-bottom: 8px;
        text-align: left;
    }
    
    .separador_lambda td {
        border-top: 3px solid white;
        height: 18px;
        padding: 0px;
    }
    </style>
    </head>
    
    <body>
    """)
    
            f.write("""
    <p class="info">
    Results obtained for a maximum of 25 iterations.
    &nbsp;&nbsp;&nbsp;
    First row: AutoScale.
    &nbsp;&nbsp;&nbsp;
    Second row: SameScale.
    &nbsp;&nbsp;&nbsp;
    Lambda values are obtained by: lambdas = np.logspace(-6, 1, 12).
    </p>
    <hr>
    """)
    
            f.write("<table>\n")
    
            f.write("<tr>\n")
            for _, titulo in colunas:
                f.write(f"<th>{titulo}</th>\n")
            f.write("</tr>\n")
    
            lambdas_ordenados = sorted(grupos.keys(), key=lambda x: float(x))
    
            for lambda_str in lambdas_ordenados:
    
                # separador ANTES de cada lambda
                f.write(f"""
    <tr class="separador_lambda">
        <td colspan="{len(colunas)}"></td>
    </tr>
    """)
    
                f.write(f"""
    <tr>
        <td colspan="{len(colunas)}" class="lambda_bloco">
            λ = {float(lambda_str):.2e}
        </td>
    </tr>
    """)
    
                # AutoScale
                f.write("<tr>\n")
    
                for chave, _ in colunas:
    
                    img = grupos[lambda_str]["AutoScale"].get(chave)
    
                    if img is not None:
                        titulo = dict(colunas)[chave]
    
                        f.write(f"""
    <td>
        <div class="titulo">{titulo}</div>
        <div class="escala">AutoScale</div>
        <div class="lambda">λ = {float(lambda_str):.2e}</div>
        <img src="{img}">
    </td>
    """)
                    else:
                        f.write("<td></td>\n")
    
                f.write("</tr>\n")
    
                # SameScale
                f.write("<tr>\n")
    
                for chave, _ in colunas:
    
                    img = grupos[lambda_str]["SameScale"].get(chave)
    
                    if img is not None:
                        titulo = dict(colunas)[chave]
    
                        f.write(f"""
    <td>
        <div class="titulo">{titulo}</div>
        <div class="escala">SameScale</div>
        <div class="lambda">λ = {float(lambda_str):.2e}</div>
        <img src="{img}">
    </td>
    """)
                    else:
                        f.write("<td></td>\n")
    
                f.write("</tr>\n")
    
            f.write("""
    </table>
    </body>
    </html>
    """)

    #print("HTML salvo em:", html_path)
    ###############################################################################
    # Essa função calcula o problema inverso
    ###############################################################################
    def solve(self, V_measured,initialEstimate=1.0, alpha =1.0,  Lambda = 0.50, max_iter=500, Tol=1.0e-4, iteration=0, html_name = None):
        print('html_name_solve',html_name)
        itr_start = int(iteration)
        ultimos10 = []
        ultimaNorma =[99,99,99]
        lastResidue = [99,99,99]
        listXplot = []                                                         # Lista armazenar eixo x das iterações para plotar
        listaItrPlot = []                                                      # Lista Valores da normaSigma para plotar
        #centroids_2D = np.array([elem.Centroid for elem in self.mymesh.Elements])
        
        #self.plotMSH(self.mymesh.sigma_vec, save = False)
        
        #L2 = self.calc_L2_gauss_2D(centroids_2D)

        centroids_2D = np.array([elem.Centroid for elem in self.mymesh.Elements])

        L2, idx_dom, std_auto = self.calc_L2_gauss_2D_only_domain(centroids_2D, std=0.01)

        #print(f"Filtro gaussiano: std = {std_auto:.6e}")
        #print(f"Número de elementos no domínio físico = {len(idx_dom)}")

        #L2 = self.calc_L2_gauss_mean_2D(centroids_2D)
        difResidue = 0
        normaDeltaTemp = 0
        fatorAlpha = 0.99
        contNorma = 0
        contNeg = 0
        contItr = 0
        V_measured = V_measured.T
        V_measured = V_measured.reshape(-1, 1)
        #print('V_measured', V_measured.shape)
        


        initialEstimate = np.array(initialEstimate, dtype=float)



        n_elem_phys = self.mymesh.NumberOfPhysicalElements
        sigmaInicial = np.tile(initialEstimate, (n_elem_phys, 1))
        #sigmaInicial = initialEstimate
        self.sigmaStar = np.tile(initialEstimate, (n_elem_phys, 1))


        #sigmaInicial_vec = sigmaInicial.reshape(-1, 1, order='C')
        #sigmaStar_vec    = self.sigmaStar.reshape(-1, 1, order='C')



        #######################################################################
        #  Regularização Filtro Gaussiano FP_Alta
        # fonte: Erick equação C.37
        #######################################################################
        """
            Atualização Gauss-Newton regularizada:

            θ̂(k+1) = θ̂(k) + αk(JkᵀW1Jk + λ²L2ᵀL2)⁻¹ [ JkᵀW1(z-h(θ̂k)) - λ²L2ᵀL2(θ̂k-θ*) ]

            Termos:

            θ̂  -> condutividade reconstruída
            J   -> Jacobiano
            W1  -> pesos das medidas
            L2  -> regularização espacial
            λ   -> parâmetro de regularização
            α   -> passo iterativo
            z   -> tensão medida
            h() -> problema direto
            θ*  -> prior
        """
        #######################################################################
        ###################        MAIN LOOP   ################################
        #######################################################################
        for itr in range(itr_start,max_iter):
            #np.savetxt("lastIteration.txt", np.array([itr]), fmt="%d") # Main Loop
            contItr = contItr + 1
            #Vtemp = self.CalcTempKGlobalAnisotropicHua(sigmaInicial)
            Vtemp = self.CalcTempKGlobalAnisotropicHua_sparse(sigmaInicial)
            #Vtemp = self.apply_boundary_conditions(Vtemp.copy())
            Vtemp = self.apply_boundary_conditions_sparse(Vtemp.copy())

            #print("posto Vtemp =", np.linalg.matrix_rank(Vtemp))
            #print("shape Vtemp =", Vtemp.shape)

            #invVtemp = np.linalg.inv(Vtemp)

            
            solve_K = factorized(Vtemp.tocsc())

            V_calc = solve_K(self.vetor_corrente_cond_contorno)

            #V_calc = np.dot(invVtemp, self.vetor_corrente_cond_contorno)       # Calcula Valor estimado
            V_calc_noh = V_calc[self.mymesh.ElectrodeNodes]                    # pega somente valores dos eletrodos
            
            V_calc_noh = V_calc_noh.T
            V_calc_noh = V_calc_noh.reshape(-1, 1)



                        # ***** Determinação do resíduo *****
            residue = V_calc_noh - V_measured                                  # Calcula resíduo matriz Nele X Nele
            normaResidue = np.linalg.norm(residue)
            listaItrPlot.append(normaResidue)
            #print('residue',residue.shape)
            #print('normaResidue',normaResidue)
            #self.calc_Y_jacobian()      # Calcula (dY/dσ_k) do Jacobiano
            self.calc_Y_jacobian_anisotropic_hua()      # Calcula (dY/dσ_k) do Jacobiano
            
    
            
            # ***** Cálculo J = - Y_inv * (dY/ds_k) * (Y_inv * C) *****
            #self.Calc_J(invVtempJ)
            #self.Calc_J(invVtemp)
            self.Calc_J_sparse(solve_K)
            
            # ***** Cálculo do termo 1a JT_W1_J *****

            #L2_aniso = np.kron(np.eye(3), L2)
            ###L2_aniso = np.kron(L2, np.eye(3))
            # ***** Cálculo do termo 1b Lambda^2 * LT_L *****
            #LTL =np.dot(L2.T, L2)
            ####LTL = L2_aniso.T @ L2_aniso
            
            L2_aniso = kron(csr_matrix(L2), eye(3, format="csr"), format="csr")
            LTL = L2_aniso.T @ L2_aniso
            termo_L = (Lambda**2) * LTL
            
            termo_L = (Lambda**2) * LTL
            #firstTerm = self.JTJ + termo_L
            #ztermo_reg = (sigma_inicial - self.sigmaStar)
            
            sigmaInicial_vec = sigmaInicial.reshape(-1, 1, order='C')
            sigmaStar_vec = self.sigmaStar.reshape(-1, 1, order='C')

            #zregularizacao = (Lambda**2)*zLTL @ ztermo_reg            
            

            #print('termo_L',termo_L.shape)
            #print('self.JTJ',self.JTJ.shape)
            
            '''
            # ***** Cálcula e inverte termo 1 -> (JTWJ + Lambda^2*LTL)^-1 *****
            firstTerm = self.JTJ + termo_L
            
            #inv_firstTerm = np.linalg.inv(firstTerm)
            
            firstTerm = csr_matrix(self.JTJ) + termo_L
            
            
            # ***** Cálculo do termo 2a (JT_W1_residue) *****
    
            JTW_H = np.dot(self.TempJ.T,residue)      
            
            # ***** Cálculo do termo 2b (Lambda^2 * LTL*residue) *****
            #regTerm = (sigmaInicial)# - sigmaStar)* (Lambda**2)
            regTerm = (sigmaInicial_vec)# - sigmaInicial_vec)
            

            

            regularization = np.dot(termo_L, regTerm)
            

            
            # ***** Cálculo final do termo 2 *****
            secondTerm = -JTW_H - regularization
            
            
            # ***** Produto entre termo 1 e termo 2 *****
            #deltaSigma = np.dot(inv_firstTerm, secondTerm)
            
            deltaSigma = spsolve(firstTerm.tocsc(), secondTerm).reshape(-1, 1)
            '''
            JTJ_sparse = csr_matrix(np.asarray(self.JTJ, dtype=np.float64))

            firstTerm = JTJ_sparse + termo_L
            firstTerm = firstTerm.astype(np.float64).tocsc()
            
            JTW_H = self.TempJ.T @ residue
            JTW_H = np.asarray(JTW_H, dtype=np.float64)
            
            regTerm = np.asarray(sigmaInicial_vec, dtype=np.float64)
            
            regularization = termo_L @ regTerm
            regularization = np.asarray(regularization, dtype=np.float64)
            
            secondTerm = -JTW_H - regularization
            secondTerm = np.asarray(secondTerm, dtype=np.float64).ravel()
            
            #print("firstTerm dtype:", firstTerm.dtype)
            #print("secondTerm dtype:", secondTerm.dtype)
            #print("firstTerm shape:", firstTerm.shape)
            #print("secondTerm shape:", secondTerm.shape)
            
            deltaSigma = spsolve(firstTerm, secondTerm).reshape(-1, 1)
            
            ######## deltaSigma[1::3] = 0.0 # zerar Δσxy   # tentar usa sigma XY
            
            #print('deltaSigma', deltaSigma)#[:1])
            #deltaSigma = deltaSigma[:, 0]*alpha
            
            alphaDeltaSigma = alpha*deltaSigma
            
            ###### sigmaInicial[:,1] = 0  # tentar usa sigma XY
            
            #print('alphaDeltaSigma',alphaDeltaSigma[:1])
            normaDelta = np.linalg.norm(alphaDeltaSigma)                       # Calcula norma  delta sigma
            plotItr = np.linalg.norm(alphaDeltaSigma)                          # Armazena delta sigma para plot
            #print('normaDelta',normaDelta)
            listXplot.append(itr)                                              # Armazena o índice da iteração
                                             # Armazena o valor a ser plotado
            
            ultimaNorma.append(normaDelta)
            if len(ultimaNorma) > 2:
                ultimaNorma.pop(0)
                
            lastResidue.append(normaResidue)                                   # Calcula norma  resíduo
            
            if len(lastResidue) > 2:
                lastResidue.pop(0)
                difResidue =  lastResidue[2]- lastResidue[1]                                           # Armazena 3 últimos valores da norma lastResidue
            #print(f'{itr} - nDelta = {normaDelta}, nResidue = {normaResidue}, alfa = {alpha} ')
            print(f'{normaResidue} - {normaDelta} - {itr}')
            
            if normaDelta < Tol:    # Convergência atingida se a norma de # delta_sigam < que  1e-6
              print(f'Convergência atingida após {itr} iterações.')
              convergencia = True
              break
                                                                        # interrompe o processo de iteração
            if normaResidue < Tol:    # Convergência atingida se a norma de # delta_sigam < que  1e-6
              print(f'Convergência atingida após {itr} iterações.')
              convergencia = True
              break
            
            #print('sigmaInicial_antes', sigmaInicial[:1])
            sigmaPlusOne = (sigmaInicial_vec + alphaDeltaSigma)
            #sigmaPlusOne = np.clip(sigmaPlusOne, 1.99, 3.01)
            sigmaInicial_vec = sigmaPlusOne
            #print('sigmaInicial_depois', sigmaInicial_vec[:3])
            sigmaInicial = sigmaInicial_vec.reshape(-1, 3, order='C')
            
            ########### sigmaInicial[:,1] = 0 # zerar σxy  # tentar usa sigma XY

            # Impõe condutividade mínima
            sigmaInicial[:, 0] = np.clip(sigmaInicial[:, 0], 0.1, 6.0)  # σxx           
            sigmaInicial[:, 2] = np.clip(sigmaInicial[:, 2], 0.1, 6.0)  # σyy

            alphaFator = 0.5
            limite = alphaFator * np.sqrt(sigmaInicial[:, 0] * sigmaInicial[:, 2])        
            sigmaInicial[:, 1] = np.clip(sigmaInicial[:, 1], -limite, limite)
            

            
            ultimos10.append(sigmaPlusOne)                                     # Armazena 10 últimos valores de sigmaPlusOne
            if len(ultimos10) > 5:
                ultimos10.pop(0)
            
            
            if ultimaNorma[2] > ultimaNorma[1]:
               print(f'Encontrou norma normaDelta maior  que a anterior.')
               break
                
            
            if lastResidue[2] > lastResidue[1]:
               print(f'Encontrou norma lastResidue maior  que a anterior.')
               break
               #self.plotMSH(sigmaInicial,itr, save = True)
               #alpha = alpha*fatorAlpha
               #break

            #if contItr ==50:
            #    np.savetxt('sigma_inicial_bk_{itr}.txt', sigmaInicial, fmt="%.8f")
            #    contItr = 0
            
            #if itr % 10 == 0:   # salva de 1000 em 1000 ...
            #    self.plotMSH(sigmaInicial[:, 0],Lambda, itr, save = True)
            #    self.plotMSH(sigmaInicial[:, 2],Lambda, itr, save = True)
            
        ##################################################################################
        # Deteminação das condutividas longitudinal e transversal do material anisotrópico
        # e do ângulo theta.
        ##################################################################################
        """
        Condutividades principais:

                σxx + σyy
        Sx = ---------------- + √( ((σxx-σyy)/2)^2 + σxy² )
                    2
                σxx + σyy
        Sy = ---------------- - √( ((σxx-σyy)/2)^2 + σxy² )
                    2
        Ângulo:
        θ = 0.5 atan2(2σxy , σxx-σyy)
        Tensor:
                [σxx  σxy]
        σ =     [        ]
                [σxy  σyy]
        """
        ########################################################################################
        
        Smed = 0.5 * (sigmaInicial[:, 0] + sigmaInicial[:, 2])
        D = np.sqrt(((sigmaInicial[:, 0] - sigmaInicial[:, 2])/2)**2 + sigmaInicial[:, 1]**2)

        sigma_L = Smed + D
        sigma_T = Smed - D
        
        
        DifAnisotropia = sigma_L - sigma_T
        DifAnisotropia_Med =  np.abs(np.mean(DifAnisotropia))
        
        # para evitar saturação para  +/- 90 devido a arctan

        #anisotropia = np.abs(sigma_L - sigma_T)       # verifica se a diferença Sxx - Syy é muito pequena
        #np.savetxt("sigmaInicialTheta.txt", sigmaInicial)
        
        #anisotropia = np.abs(sigmaInicial[:, 0] - sigmaInicial[:, 2])       # verifica se a diferença Sxx - Syy é muito pequena
        #np.savetxt("anisotropia.txt", anisotropia)  # formato binário
        
        sigma_theta = sigmaInicial[:, 0] - sigmaInicial[:, 2]
        #np.savetxt("sigma_theta.txt", sigma_theta)
        theta_rad = 0.5 * np.arctan2(2.0 * sigmaInicial[:, 1], sigma_theta)
        #np.savetxt("theta_rad.txt", theta_rad)
        theta_deg = np.rad2deg(theta_rad)
        #np.savetxt("theta_deg.txt", theta_deg)
        #tolerance = 0.01* np.max(anisotropia)

        #theta_deg[anisotropia < tolerance] = 0.0
        #np.savetxt("theta_final.txt", theta_deg)
        #print('tolerance',tolerance)
        #print('anisotropia',anisotropia)
        tol_aniso = 1e-4

        quase_isotropico = np.abs(sigma_L - sigma_T) < tol_aniso

        theta_deg[quase_isotropico] = 0.0
        
     
          # formato binário
        #theta_rad = 0.5 * np.arctan2(2*sigma_xy, sigma_xx - sigma_xy)
        #theta_deg = np.rad2deg(theta_rad)
        
        
        
        lista_imgs = []
        pasta_base = "../../docs/"
        pasta_teste = os.path.join(pasta_base, html_name)
        print('pasta_teste', pasta_teste)
        os.makedirs(pasta_teste, exist_ok=True)
        # σxx, σxy, σyy (MSH) f"{Lambda:.6f}"
        
        
        nome1 =  f'{pasta_teste}/{html_name}_sigma_xx_{Lambda:.10f}' 
        self.plotMSH(sigmaInicial[:,0], Lambda, itr, save=True, SigmaXXXYYY='xx', DifAniso = DifAnisotropia_Med, nome_arquivo=nome1)
        lista_imgs.append(nome1)
        
        nome2 =  f'{pasta_teste}/{html_name}_sigma_xy_{Lambda:.10f}' 
        self.plotMSH(sigmaInicial[:,1], Lambda, itr, save=True, SigmaXXXYYY='xy', DifAniso = DifAnisotropia_Med, nome_arquivo=nome2)
        lista_imgs.append(nome2)
        
        nome3 = f'{pasta_teste}/{html_name}_sigma_yy_{Lambda:.10f}'  
        self.plotMSH(sigmaInicial[:,2], Lambda, itr, save=True, SigmaXXXYYY='yy', DifAniso = DifAnisotropia_Med, nome_arquivo=nome3)
        lista_imgs.append(nome3)
             

        nome4 = f'{pasta_teste}/{html_name}_sigma_L_{Lambda:.10f}'          
        self.plotMSH(sigma_L, Lambda, itr, save=True, SigmaXXXYYY='L', DifAniso = DifAnisotropia_Med, nome_arquivo=nome4)
        lista_imgs.append(nome4)
        
        nome5 = f'{pasta_teste}/{html_name}_sigma_T_{Lambda:.10f}'  
        self.plotMSH(sigma_T, Lambda, itr, save=True, SigmaXXXYYY='T', DifAniso = DifAnisotropia_Med, nome_arquivo=nome5)
        lista_imgs.append(nome5)
        
        # Diferença anisotrópica
        #nome6 =  f'{pasta_teste}/{html_name}_sigma_Dif_{Lambda:.6f}' 
        #self.plotMSH(DifAnisotropia, Lambda, itr, save=True, SigmaXXXYYY='L-σT', DifAniso = DifAnisotropia_Med, nome_arquivo=nome6)
        #lista_imgs.append(nome6)
        
        # Theta amgle
        nome7 =  f'{pasta_teste}/{html_name}_theta_deg_{Lambda:.10f}' 
        self.plotMSH(theta_deg, Lambda, itr, save=True, SigmaXXXYYY='θ°', DifAniso = DifAnisotropia_Med, nome_arquivo=nome7)
        lista_imgs.append(nome7)  
        print('theta_deg', theta_deg)
        # Gráfico tipo linha (o que você mandou)
        nome8 = f'{pasta_teste}/{html_name}_sigma_linhas_{Lambda:.10f}'
        self.plot_sigma(sigmaInicial,  ref_sL = 3.0, ref_sT = 1.0, salvar=True, nome_arquivo=nome8)
        lista_imgs.append(nome8)
        #print('plot_sigma_banana_solver', theta_deg)
        
        nome9 = f'{pasta_teste}/{html_name}_iterations_{Lambda:.10f}'
        self.plotar_iteracoes(listXplot, listaItrPlot, nome_arquivo = nome9)
        lista_imgs.append(nome9)  
        
        
        '''
        # Gráfico tipo linha (o que você mandou)
        nome10 = f'{pasta_teste}/{html_name}_theta_{Lambda:.6f}'
        self.plot_theta_deg(theta_deg, salvar=True, nome_arquivo=nome10)
        lista_imgs.append(nome10)
        

        
        nome11 = f'{pasta_teste}/{html_name}_Mask_{Lambda:.6f}'  
        self.plotMask(DifAnisotropia, Lambda, itr, save=True, SigmaXXXYYY='y', DifAniso = DifAnisotropia_Med, nome_arquivo=nome11)
        lista_imgs.append(nome11)

        nome12 = f'{pasta_teste}/{html_name}_Sorting_{Lambda:.6f}'
        self.plot_Sorting(DifAnisotropia, titulo="Results (σx-σy)", salvar=True, nome_arquivo=nome12)
        lista_imgs.append(nome12)
        '''
        nome13 = f'{pasta_teste}/{html_name}_Elipses_{Lambda:.10f}'     
        self.plotElipse(sigmaInicial, sigma_L, sigma_T, theta_deg, Lambda, itr, save=True, DifAniso = DifAnisotropia_Med, nome_arquivo=nome13)
        lista_imgs.append(nome13)
        # Theta angle

        nome14 = f'{pasta_teste}/{html_name}_Table_{Lambda:.10f}'     
        self.plot_Table(nome_arquivo=nome14)
        lista_imgs.append(nome14)
        
        
        # Criar HTML
        #nome_html = '../docs/figureTemp/{html_name}_{Lambda}.html'
        #nome_html = f'../docs/figureTemp/{html_name}_{Lambda:.6f}.html'
        #self.salvar_html(lista_imgs, nome_html)
        
        self.salvar_html_todos_lambdas(pasta_teste, html_name)
        
    

        
        #self.plot_espectro(sigmaInicial)
        #print('sigmaInicial',sigmaInicial_vec)
        np.savetxt("sigmaInicial_result.txt", sigmaInicial)  # formato binário
        s = np.linalg.svd(self.TempJ, compute_uv = False)
        '''
        # --- 2) Gráfico tipo semilogy
        plt.figure()
        plt.semilogy(s, 'o-')
        plt.grid(True, which='both')
        plt.xlabel('índice')
        plt.ylabel('valor singular (log)')
        plt.title('Espectro singular do Jacobiano')
        plt.show()
        #plt.show(block=False)
        #plt.pause(0.01)  
        # --- 3) Rank efetivo (mesma ideia do seu código)
        '''
        tol = 1e-6 * s[0]          # ou outro fator
        rank_eff = np.sum(s > tol)
        print(f'rank efetivo ~ {rank_eff} (tol={tol:g})')
     
 
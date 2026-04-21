# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 13:56:12 2025

@author: rodgr
"""



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import mesh
import elements
import gmsh
import sys
from datetime import datetime
from matplotlib import cm
timestamp = datetime.now().strftime("%Y%m%d_%H%M")


#import subprocess
#import os


class inverse_problem: 
    def __init__(self, mymesh, V_imposto=None, Pcorrente=None, SkipPattern=None, VirtualNode = False, I =1.0e-3, name = None, pLambda = None, saveVideo = False):
        if not hasattr(mymesh, "KGlobal"): # verifica se o objeto mymesh tem um atributo chamado KGlobal.
            raise TypeError("Parâmetro incorreto: mymesh.")

        self.mymesh = mymesh
        self.KGlobalBoundary = self.CalcKGlobalBoundary()
        self.Vmedido = None
        self.TempJ = None
        self.KGlobalTemp = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes), dtype=float)
        self.name = name
        self.pLambda = pLambda
        self.saveVideo = saveVideo
        self.contsaveVideo = 0
        self.KGlobal= self.mymesh.KGlobal
        print('self.KGlobal_tst',self.KGlobal)
        





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
    

    ###############################################################################
    def CalcKGlobalBoundary(self):
        KGlobalBoundary = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes), dtype=float)

        for elem in range(self.mymesh.NumberOfElements):
            if not self.mymesh.Elements[elem].FlagIsElectrode:
                continue

            topo = self.mymesh.Elements[elem].Topology
            Ke = self.mymesh.Elements[elem].KGeo

            for i in range(len(topo)):
                no_i = topo[i]
                for j in range(len(topo)):
                    no_j = topo[j]
                    KGlobalBoundary[no_i, no_j] += Ke[i, j]
        #np.savetxt("matrizKGlobalBoundary.txt", KGlobalBoundary, fmt="%.6f")
        return KGlobalBoundary
    ###############################################################################

    def CalcKGlobalTriangle(self, sigma_tri):
        KGlobalTriangle = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes), dtype=float)

        idx_sigma = 0
        for elem in range(self.mymesh.NumberOfElements):
            if self.mymesh.Elements[elem].FlagIsElectrode:
                continue

            topo = self.mymesh.Elements[elem].Topology
            Ke = self.mymesh.Elements[elem].KGeo * sigma_tri[idx_sigma]

            for i in range(len(topo)):
                no_i = topo[i]
                for j in range(len(topo)):
                    no_j = topo[j]
                    KGlobalTriangle[no_i, no_j] += Ke[i, j]

            idx_sigma += 1
        #np.savetxt("matrizKGlobalTriangle.txt", KGlobalTriangle, fmt="%.6f")
        return KGlobalTriangle

    def area_triangulo(x1, y1, x2, y2, x3, y3):                      # calcula a área do triângulo com 3 pontos
        matriz = [ [1, x1, y1], [1, x2, y2], [1, x3, y3]]

        area_triangulo = (np.linalg.det(matriz))/2
        return abs(area_triangulo)
    ################
            

    def calc_Y_jacobian(self):   # função para calcular a matriz global
        self.Y_jacobian = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))
        self.listJacobian = []
        for i in range(0, self.mymesh.NumberOfElements):        # laço para montar a matriz global
        #for i in range(0, 1):        # laço para montar a matriz global
            if self.mymesh.Elements[i].FlagIsElectrode:
                continue   # pula elementos Hua/eletrodo
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

            
            Y_local = (1/(4*area_triangle))*np.array([[C_11, C_12, C_13], 
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
            temp = self.apply_boundary_conditions(self.Y_jacobian)
            self.listJacobian.append(temp)
            
            #print('Y_jacobian calc_Y_jacobian \n',self.Y_jacobian)
            
            self.Y_jacobian = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))
        #print('self.listJacobian.append \n',self.listJacobian)
        #return Y_jacobiano
    ###############################################################################
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
        #ztermo_reg = (sigma_inicial - self.sigmaStar)
        ztermo_reg = (sigma_inicial)# - self.sigmaStar)
        zregularizacao = (Lambda**2)*zLTL @ ztermo_reg
        segundoTermo = JTW_zh - zregularizacao
        return -alpha*(inv_primeiroTermo @ segundoTermo)
    ###############################################################################


    def CalcTempKGlobal(self, SigmaTemp):
        self.KGlobalTemp = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes), dtype=float)
        #print(f'(self.mymesh.Elements.KGeo {self.Elements.KGeo})')
        for elem in range(self.mymesh.NumberOfElements): # para cada elemento:
            #if self.mymesh.Elements[elem].FlagIsElectrode:
            #    continue   # pula elementos Hua/eletrodo
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
        #print(f' self.KGlobalTemp \n {self.KGlobalTemp}')
        return self.KGlobalTemp
    ###############################################################################

    def Calc_J(self, invVtemp):
        listTempJ = []
        termo1 = np.dot(invVtemp, self.vetor_corrente_cond_contorno)

        for idx in range(len(self.listJacobian)):
            termo2 = np.dot(self.listJacobian[idx], termo1)
            termo3 = -np.dot(invVtemp, termo2)
            termo3 = termo3[self.mymesh.ElectrodeNodes]
            termo3 = termo3.reshape(-1, 1, order='F')
            listTempJ.append(termo3)

        self.TempJ = np.concatenate(listTempJ, axis=1)
        self.JTJ = np.dot(self.TempJ.T, self.TempJ)

    ###############################################################################
    # Essa função calcula FPA com distância de cada elemento
    ############################################################################### 
    
    def calc_L2_gauss_2D(self, centroids_2D, std=0.0055000, tol=1e-9):
        ne = len(centroids_2D)
        L2 = np.zeros((ne, ne), dtype=np.float32)       
        #L2 = np.zeros((self.mymesh.NumberOfElements-self.mymesh.NumberOfElectrodes, self.mymesh.NumberOfElements-self.mymesh.NumberOfElectrodes), dtype=np.float32)   
        #for i in range(self.mymesh.NumberOfElements-self.mymesh.NumberOfElectrodes):
        for i in range(ne):
            ci = centroids_2D[i]
            soma = 0.0   
            # --- primeira passada: calcula fatores gaussianos ---
            #for j in range(self.mymesh.NumberOfElements-self.mymesh.NumberOfElectrodes):
            for j in range(ne):
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
            #for j in range(self.mymesh.NumberOfElements-self.mymesh.NumberOfElectrodes):
            for j in range(ne):
                if soma > 0:
                    if i == j:
                        aux = 1.0 - L2[i, j] / soma
                    else:
                        aux = -L2[i, j] / soma
                    L2[i, j] = aux if np.abs(aux) > tol else 0.0
        '''
        # plot  matrix sparsity 
        plt.figure(figsize=(6, 5))
        plt.spy(L2, markersize=1)
        plt.title('HPFilter matrix sparsity pattern', fontsize=15)
        plt.xlabel('Colun', fontsize=12)
        plt.ylabel('Line', fontsize=12)
        #plt.tight_layout()
        N = L2.shape[0]
        X, Y = np.meshgrid(np.arange(N), np.arange(N))      
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection='3d')   
        surf = ax.plot_surface( X, Y, L2, cmap='viridis',  linewidth=0, antialiased=True)      
        ax.set_xlabel('Coluna')
        ax.set_ylabel('Linha')
        #ax.set_zlabel('L[i,j]')
        ax.set_title('HPFilter – Superfície 3D')       
        fig.colorbar(surf, shrink=0.5)
        plt.show()    #plt.show(block = false)     #plt.pause(0.1)
        plt.show(block=False)   # NÃO bloqueia
        plt.pause(3)            # tempo que o gráfico fica na tela (segundos)
        plt.close('all')        # fecha tudo automaticamente
        '''
        return L2

    
    ###############################################################################
    # Essa função calcula FPA com distância média entre os elementoss
    ###############################################################################    

    def calc_L2_gauss_mean_2D(self, centroids_2D,  tol=1e-9):
   
        d_media = np.mean(np.linalg.norm(centroids_2D[1:] - centroids_2D[:-1], axis=1))
        std = 0.05*d_media
        print(f'std = {std}')
        ne = centroids_2D.shape[0]
        L = np.zeros((ne, ne), dtype=np.float64)
    
        for j in range(ne):                    # coluna 
            cj = centroids_2D[j]
            soma = 0.0
    
            for i in range(ne):                # linha 
                ci = centroids_2D[i]
                dist = np.linalg.norm(ci - cj)
                if dist <= 5 * std:
                    g = np.exp(-dist**2 / (2 * std**2))
                else:
                    g = 0.0   
                L[i, j] = g
                soma += g
    
            if soma > 0:
                L[:, j] /= soma
    
        L[np.abs(L) < tol] = 0.0 
        I = np.eye(ne)
    
        # filtro passa–alta
        L2 = I - L
        
        i, j =  np.nonzero(L)
        values = L2[i, j]


        # plot  matrix sparsity 
        plt.figure(figsize=(6, 5))
        plt.spy(L2, markersize=1)
        plt.title('HPFilter matrix sparsity pattern', fontsize=15)
        plt.xlabel('Colun', fontsize=12)
        plt.ylabel('Line', fontsize=12)
        #plt.tight_layout()
        #plt.show()

        #plt.show(block = false)
        #plt.pause(0.1)

        return L2

    ###############################################################################
    # Essa função plota o gráfico convergência das iterações
    ###############################################################################

    def plotar_iteracoes(self,lista_indice, lista_valor):
        plt.figure(figsize=(6, 5))
        plt.plot(lista_indice,
                lista_valor,
                marker='.',
                linestyle='-',
                color='b',
                label='Norm $\Delta\sigma$')        # plota gráfico das iterações,
        plt.xlabel("Iteration", fontsize=12)
        plt.ylabel("Norm [$\Delta\sigma$]", fontsize=12)
        plt.title("Optimization", fontsize=15)
        plt.legend()
        plt.grid(True)
        #plt.tight_layout()
        #plt.show(block=False)   # NÃO bloqueia
        #plt.pause(3)            # tempo que o gráfico fica na tela (segundos)
        #plt.close('all')        # fecha tudo automaticamente
    ###############################################################################
    ###############################################################################
    # Essa função plota o gráfico da condutividade da malha
    ###############################################################################

    
    def plotMSH(self, sigma, iteration = None, save = False, saveVideo = False):
        self.saveVideo = saveVideo
        
        x, y = self.mymesh.Coordinates[:, 0], self.mymesh.Coordinates[:, 1]
        topo = self.mymesh.msh_topology
    
        # --- Separar elementos 2D (triangulares) e 1D (linhas) ---
        elems_2D = np.array([el for el in topo if len(el) == 3])
        elems_1D = np.array([el for el in topo if len(el) == 2])
    
        fig, ax = plt.subplots(figsize=(6, 5))
    
        # ============================================================
        #  1) PLOTAR ELEMENTOS 2D (TRIANGULARES)
        # ============================================================
        if len(elems_2D) > 0:
            triang = tri.Triangulation(x, y, elems_2D)
            #tpc = ax.tripcolor(triang,facecolors=sigma[:len(elems_2D)],edgecolors='k', cmap='Blues'#,vmin=1.0 )
            ntri = triang.triangles.shape[0]
            fc = sigma.ravel()[:ntri]

            tpc = ax.tripcolor(triang,facecolors = fc,edgecolors='k', cmap='Blues',vmin=1.75,vmax=3.55 )

            fig.colorbar(tpc, ax=ax, label='σ (Conductivity)')
            if save == True:
                timestamp = datetime.now().strftime("%m%d_%H%M")
                ax.set_title(f"Conductivity Real (σ) - Lambda_{self.pLambda:.2e} - itr_{iteration}", fontsize=12)
            if save == False:
                ax.set_title(f"Conductivity Real (σ) ", fontsize=15)
        # ============================================================
        #  2) PLOTAR ELEMENTOS 1D (SEGMENTOS)
        # ============================================================
        if len(elems_1D) > 0:
            for (n1, n2) in elems_1D:
                x_coords = [x[n1], x[n2]]
                y_coords = [y[n1], y[n2]]
                ax.plot(x_coords, y_coords, color='red', linewidth=2)
    
        # ------------------------------------------------------------
        ax.set_xlabel("[m]", fontsize=12)
        ax.set_ylabel("[m]", fontsize=12)
        plt.tight_layout()
        #if save == True:
        #    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
        if self.saveVideo == False:
            nome_arquivo = f"../../docs/figures/{self.name}_{self.pLambda:.4f}.png"
        if self.saveVideo == True:
            self.contsaveVideo = self.contsaveVideo + 1
            nome_arquivo = f"../../docs/figures/{self.name}_{self.contsaveVideo:03d}.png"
        plt.savefig(nome_arquivo, dpi=100, bbox_inches='tight')

        #    plt.savefig(f"cond3_obj_skip2_v2_{timestamp}.png",     dpi=300, bbox_inches='tight')
        #plt.show()
        plt.show(block=False)   # NÃO bloqueia
        plt.pause(0.1)            # tempo que o gráfico fica na tela (segundos)
        plt.close('all')        # fecha tudo automaticamente
    ###############################################################################    
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
        #plt.show(block=False)   # NÃO bloqueia
        #plt.pause(3)            # tempo que o gráfico fica na tela (segundos)
        #plt.close('all')        # fecha tudo automaticamente
    ###############################################################################
    # Essa função calcula o problema inverso
    ###############################################################################
    #def solve(self, V_measured,initialEstimate=1.0, alpha =1.0,  Lambda = 0.50, max_iter=500, Tol=1.0e-6, iteration=0):
    def solve(self, V_measured, initialEstimate=1.0, alpha=1.0,Lambda=0.50, max_iter=500,Tol=1.0e-6, iteration=0, name=None, pLambda=None, saveVideo = False):
        self.name = name
        self.pLambda = pLambda
        self.saveVideo = saveVideo

        itr_start = int(iteration)
        ultimos10 = []
        ultimaNorma =[99,99,99]
        lastResidue = [99,99,99]
        listXplot = []                                                         # Lista armazenar eixo x das iterações para plotar
        listaItrPlot = []                                                      # Lista Valores da normaSigma para plotar
        centroids_2D = np.array([elem.Centroid for elem in self.mymesh.Elements])
        

        tri_elements = [elem for elem in self.mymesh.Elements if not elem.FlagIsElectrode]
        centroids_2D = np.array([elem.Centroid for elem in tri_elements])
        L2 = self.calc_L2_gauss_2D(centroids_2D)

        difResidue = 0
        normaDeltaTemp = 0
        fatorAlpha = 0.99
        contNorma = 0
        contNeg = 0
        contItr = 0
        V_measured = V_measured.T
        V_measured = V_measured.reshape(-1, 1)

        sigmaInicial = np.ones(len(centroids_2D))*initialEstimate
        
        sigmaInicial = sigmaInicial.reshape(-1,1)


        sigmaOne = np.ones(self.mymesh.NumberOfElements)
        '''
        ###
        VtempJ = self.CalcTempKGlobal(sigmaOne)                            # calcula derivadas parciais da matriz jacobiana
        
        # ***** Determinação do Valor calculado *****
        VtempJ = self.apply_boundary_conditions(VtempJ)                        # aplica cond contorno na matriz jacobiana
        invVtempJ = np.linalg.inv(VtempJ)                                      # inverte matriz TempKGobal para jacobiana
        '''
        #######################################################################
        ###################        MAIN LOOP   ################################
        #######################################################################
        for itr in range(itr_start,max_iter):
            #if itr == 0:
            #    self.plotMSH(sigmaInicial, itr, save = True, saveVideo = self.saveVideo)
            contItr = contItr + 1
            
            K_cuba = self.CalcKGlobalTriangle(sigmaInicial)
            Vtemp = K_cuba + self.KGlobalBoundary
            Vtemp = self.apply_boundary_conditions(Vtemp)
            
            # ***** Determinação do Valor calculado *****

            invVtemp = np.linalg.inv(Vtemp)                                    # inverte matriz TempKGobal para jacobiana                    
            
            V_calc = np.dot(invVtemp, self.vetor_corrente_cond_contorno)       # Calcula Valor estimado
            V_calc_noh = V_calc[self.mymesh.ElectrodeNodes]                    # pega somente valores dos eletrodos

            V_calc_noh = V_calc_noh.T
            V_calc_noh = V_calc_noh.reshape(-1, 1)




                        # ***** Determinação do resíduo *****
            residue = V_calc_noh - V_measured                                  # Calcula resíduo matriz Nele X Nele

            normaResidue = np.linalg.norm(residue)
            listaItrPlot.append(normaResidue)

            
            self.calc_Y_jacobian()      # Calcula (dY/dσ_k) do Jacobiano
    
    
            # ***** Cálculo J = - Y_inv * (dY/ds_k) * (Y_inv * C) *****
            #self.Calc_J(invVtempJ)
            self.Calc_J(invVtemp)
            
            # ***** Cálculo do termo 1a JT_W1_J *****
            
            # ***** Cálculo do termo 1b Lambda^2 * LT_L *****
            LTL =np.dot(L2.T, L2)
            termo_L = (Lambda**2)*LTL
            
            # ***** Cálcula e inverte termo 1 -> (JTWJ + Lambda^2*LTL)^-1 *****

            firstTerm = self.JTJ + termo_L
            inv_firstTerm = np.linalg.inv(firstTerm)
            
            # ***** Cálculo do termo 2a (JT_W1_residue) *****
                 
            JTW_H = np.dot(self.TempJ.T,residue)      
            
            # ***** Cálculo do termo 2b (Lambda^2 * LTL*residue) *****
            regTerm = (sigmaInicial)# - sigmaStar)* (Lambda**2)

            regTerm = regTerm.reshape(-1, 1)
            regularization = np.dot(termo_L, regTerm)
                        
            # ***** Cálculo final do termo 2 *****
            secondTerm = -JTW_H - regularization
            
            
            # ***** Produto entre termo 1 e termo 2 *****
            deltaSigma = np.dot(inv_firstTerm, secondTerm)

            
            alphaDeltaSigma = alpha*deltaSigma
            normaDelta = np.linalg.norm(alphaDeltaSigma)                       # Calcula norma  delta sigma
            plotItr = np.linalg.norm(alphaDeltaSigma)                          # Armazena delta sigma para plot
            listXplot.append(itr)                                              # Armazena o índice da iteração
                                             # Armazena o valor a ser plotado

            ultimaNorma.append(normaDelta)
            if len(ultimaNorma) > 2:
                ultimaNorma.pop(0)
                
            lastResidue.append(normaResidue)                                   # Calcula norma  resíduo
            
            if len(lastResidue) > 2:
                lastResidue.pop(0)
                difResidue =  lastResidue[2]- lastResidue[1]                                           # Armazena 3 últimos valores da norma lastResidue

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

            sigmaPlusOne = (sigmaInicial + alphaDeltaSigma)
            sigmaPlusOne = np.clip(sigmaPlusOne, 0.1, 5.01)
            sigmaInicial = sigmaPlusOne
            
            
            ultimos10.append(sigmaPlusOne)                                     # Armazena 10 últimos valores de sigmaPlusOne
            if len(ultimos10) > 5:
                ultimos10.pop(0)

            
                
           
            if lastResidue[2] > lastResidue[1]:

               print(f'Encontrou norma lastResidue maior  que a anterior.')
               np.savetxt('sigma_inicial_cont.txt', sigmaInicial, fmt="%.8f")
               self.plotMSH(sigmaInicial, itr, save = True, saveVideo = self.saveVideo)

               break
               

            if contItr ==50:
                np.savetxt('sigma_inicial_cont.txt', sigmaInicial, fmt="%.8f")
                contItr = 0
            
            if self.saveVideo == True and  itr % 4 == 0:   # salva de 1000 em 1000 ...
                self.plotMSH(sigmaInicial, itr, save = True, saveVideo = self.saveVideo)

        #print('sigmaInicial \n', sigmaInicial) 
        np.savetxt('sigma_inicial_cont.txt', sigmaInicial, fmt="%.8f")
        #self.plotar_iteracoes(listXplot, listaItrPlot)
        self.plotMSH(sigmaInicial, itr, save = True, saveVideo = self.saveVideo)
        #self.plot_espectro(sigmaInicial)

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
        
        
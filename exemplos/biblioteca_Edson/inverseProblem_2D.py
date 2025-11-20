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
#import subprocess
#import os


class inverse_problem: 
    def __init__(self, mymesh, V_imposto=None, Pcorrente=None, SkipPattern=None, VirtualNode = False, I =1.0e-3):
        if not hasattr(mymesh, "KGlobal"): # verifica se o objeto mymesh tem um atributo chamado KGlobal.
            raise TypeError("Parâmetro incorreto: mymesh.")

        self.mymesh = mymesh
        self.Vmedido = None
        self.TempJ = None
        self.KGlobalTemp = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfNodes), dtype=float)

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
        self.chute = np.full(self.mymesh.NumberOfElements, 1.0)                           # Monta vetor chute
        
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
        for i in range(0, self.mymesh.NumberOfElements):        # laço para montar a matriz global
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
            
            print('Y_jacobian calc_Y_jacobian \n',self.Y_jacobian)
            
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
    def calc_L2_gauss_1D(self, centroids_1D, std=0.1):#, covariance_vector):
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
        tol = 1e-9

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
        ztermo_reg = (sigma_inicial - self.chute)
        zregularizacao = (Lambda**2)*zLTL @ ztermo_reg
        segundoTermo = JTW_zh - zregularizacao
        return -alpha*(inv_primeiroTermo @ segundoTermo)
    ###############################################################################
    ###############################################################################
    # Essa função plota o gráfico convergência das iterações
    ###############################################################################

    def plotar_iteracoes(self, lista_indice, lista_valor):
        plt.plot(lista_indice,
                lista_valor,
                marker='.',
                linestyle='-',
                color='b',
                label='Norma $\Delta\sigma$')        # plota gráfico das iterações,
        plt.xlabel("Iteração", fontsize=15)
        plt.ylabel("Norma Delta_Sigma", fontsize=15)
        plt.title("Otimização da Barra de Cobre", fontsize=15)
        plt.legend()
        plt.grid(True)
        plt.show()
    ###############################################################################

    def CalcTempKGlobal(self, SigmaTemp):
        #print(f'(self.mymesh.Elements.KGeo {self.Elements.KGeo})')
        for elem in range(self.mymesh.NumberOfElements): # para cada elemento:
            #print(f' SigmaTemp {SigmaTemp}')
            for i in range(len(self.mymesh.Elements[elem].Topology)): # para cada i (noh local):
                no_i = self.mymesh.Elements[elem].Topology[i] # pega noh_i (noh global)

                for j in range(len(self.mymesh.Elements[elem].Topology)): # para cada j (noh local):
                    no_j = self.mymesh.Elements[elem].Topology[j] # pega noh_j (noh global)
                    valorTemp = self.mymesh.Elements[elem].KGeo[i, j] * SigmaTemp[elem]
                    #print(f' self.Elements.KGeo {self.mymesh.Elements[elem].KGeo}')
                    #print(f' valorTemp \n {valorTemp}')
                    self.KGlobalTemp[no_i, no_j]+= valorTemp
        #print(f' self.KGlobalTemp \n {self.KGlobalTemp}')
        return self.KGlobalTemp
    ###############################################################################
    def CalcJTJ(self, invKGlobalGeo):
        listTempJ=[]
        for idx in range(self.mymesh.NumberOfElements):
            termo1 = np.dot(invKGlobalGeo, self.vetor_corrente_cond_contorno) 
            #print('termo1 \n', termo1)
            termo2 = np.dot(self.listJacobian[idx], termo1)
            #print('termo2 \n', termo2)
            termo3 = -np.dot(invKGlobalGeo, termo2)
            termo3 = termo3.reshape(-1,1,  order='F')
            #print('termo3 \n', termo3)
            listTempJ.append(termo3)
            #print('listTempJTJ \n', listTempJ)
        self.TempJ = np.concatenate(listTempJ, axis=1)

        JTJ = np.dot(self.TempJ.T, self.TempJ)
        print('JTJ \n',JTJ)
    def Calc_J(self, invVtemp):
        listTempJ=[]
        for idx in range(self.mymesh.NumberOfElements):
            termo1 = np.dot(invVtemp, self.vetor_corrente_cond_contorno) 
            #print('termo1 \n', termo1)
            termo2 = np.dot(self.listJacobian[idx], termo1)
            #print('termo2 \n', termo2)
            termo3 = -np.dot(invVtemp, termo2)
            termo3 = termo3[self.mymesh.ElectrodeNodes]
            #print('termo3 \n', termo3.shape)
            termo3 = termo3.reshape(-1,1,  order='F')
            #print('termo3b \n', termo3.shape)
            listTempJ.append(termo3)
            listTempJa = np.array(listTempJ)
            #print('listTempJxxxxxxxxxxx \n', listTempJa.shape)
        self.TempJ = np.concatenate(listTempJ, axis=1)
        #print('self.TempJ \n',self.TempJ.shape)
        JTJ = np.dot(self.TempJ.T, self.TempJ)
        #print('JTJ \n',JTJ.shape)
    ###############################################################################
    def calc_L2_gauss_2D(self, centroids_2D, std=0.1, tol=1e-9):
        
        #nelements = centroids_2D.shape[0]
        F = np.zeros((self.mymesh.NumberOfElements, self.mymesh.NumberOfElements), dtype=np.float32)
    
        for i in range(self.mymesh.NumberOfElements):
            ci = centroids_2D[i]
            soma = 0.0
    
            # --- primeira passada: calcula fatores gaussianos ---
            for j in range(self.mymesh.NumberOfElements):
                cj = centroids_2D[j]
                dist = np.linalg.norm(ci - cj)  # distância Euclidiana
    
                if dist <= 1.50 * std:
                    fator = 1.0 if i == j else np.exp(-dist**2 / (2 * std**2))
                    soma += fator
                    F[i, j] = fator
    
            # --- segunda passada: normaliza e transforma em passa-alta ---
            for j in range(self.mymesh.NumberOfElements):
                if soma > 0:
                    if i == j:
                        aux = 1.0 - F[i, j] / soma
                    else:
                        aux = -F[i, j] / soma
                    F[i, j] = aux if np.abs(aux) > tol else 0.0
        return F
    def plotMSH(self,sigma):
        x, y = self.mymesh.Coordinates[:, 0], self.mymesh.Coordinates[:, 1]
        triang = tri.Triangulation(x, y, self.mymesh.msh_topology)
        fig, ax = plt.subplots(figsize=(6, 5))
        tpc = ax.tripcolor(triang, facecolors=sigma, edgecolors='k', cmap= 'Blues', vmin=0, vmax=6)#'Greys')
        fig.colorbar(tpc, ax=ax, label='σ (condutividade)')
        ax.set_title("Condutividade Real (σ)")
        plt.xlabel("[m]")
        plt.ylabel("[m]")
        plt.tight_layout()
        plt.show()
    def solve(self, V_measured,initialEstimate=1.0, max_iter=1):
        alpha =0.0100
        lista_i = []                                        # Lista armazenar iterações
        lista_plotar = []                                      # Lista Valores de sigma
        centroids_2D = np.array([elem.Centroid for elem in self.mymesh.Elements])
        #print(centroids_2D.shape[0])  # (n_elements, 2)
        
        self.plotMSH(self.mymesh.sigma_vec)

        Lambda = 0.06
        #print(f'self.mymesh.Elements.Centroid \n {centroids_2D}')
        L2 = self.calc_L2_gauss_2D(centroids_2D)
        
        #L2 = np.eye(self.mymesh.NumberOfElements)
        print(f'L2 \n {L2.shape}')
        #plt.figure(figsize=(8, 8))
        #plt.spy(L2, markersize=1)
        #plt.title('Padrão de esparsidade da matriz FtF')
        #plt.xlabel('Coluna')
        #plt.ylabel('Linha')
        #plt.show()
        #print(f'L2 \n {L2.shape}')
        
        
        sigmaInicial = np.ones(self.mymesh.NumberOfElements)*initialEstimate
        sigmaStar = sigmaInicial*0 #np.ones(self.mymesh.NumberOfElements)*0
        sigmaOne = np.ones(self.mymesh.NumberOfElements)
        ###
        VtempJ = self.CalcTempKGlobal(sigmaInicial)                            # calcula derivadas parciais da matriz jacobiana
        
        # ***** Determinação do Valor calculado *****
        #print(f'(Vtemp {Vtemp})')
        VtempJ = self.apply_boundary_conditions(VtempJ)                          # aplica cond contorno na matriz jacobiana
        #print('Vtemp \n', Vtemp)
        invVtempJ = np.linalg.inv(VtempJ)                                        # inverte matriz TempKGobal para jacobiana
        #print('invVtemp \n', invVtemp)
        ###
        for itr in range(2):
            #print('sigmaInicial \n', sigmaInicial)
            Vtemp = self.CalcTempKGlobal(sigmaInicial)                            # calcula derivadas parciais da matriz jacobiana
            
            # ***** Determinação do Valor calculado *****
            #print(f'(Vtemp {Vtemp})')
            Vtemp = self.apply_boundary_conditions(Vtemp)                          # aplica cond contorno na matriz jacobiana
            #print('Vtemp \n', Vtemp)
            invVtemp = np.linalg.inv(Vtemp)                                        # inverte matriz TempKGobal para jacobiana
            #print('invVtemp \n', invVtemp)
            
            
            
            V_calc = np.dot(invVtemp, self.vetor_corrente_cond_contorno)           # Calcula Valor estimado
            V_calc_noh = V_calc[self.mymesh.ElectrodeNodes]                        # pega somente valores dos eletrodos
            #print('V_calc \n', V_calc.shape)
            #print('V_calc_noh \n', V_calc_noh)
            #print('self.vetor_corrente_cond_contorno \n', self.vetor_corrente_cond_contorno.shape)
            
            #print(f'V_measured \n {V_measured[0]} \n Valc \n {V_calc_noh[0]}')
            # ***** Determinação do resíduo *****
            residue = V_calc_noh - V_measured                                      # Calcula resíduo matriz Nele X Nele
            normaResidue = np.linalg.norm(residue)
            residue = np.repeat(residue, self.mymesh.NumberOfElectrodes, axis=0)       # rearranja vetor resíduo para dim do jacobiano
            #residue = np.vstack([residue] * self.mymesh.NumberOfElectrodes) 
            #print('residue \n', residue.shape)
            
            
            
            #KGlobalGeo = self.CalcTempKGlobal(sigmaOne)
            #KGlobalGeo = self.apply_boundary_conditions(KGlobalGeo)
            #invKGlobalGeo = np.linalg.inv(KGlobalGeo)
            
                                                 # calc dif entre Vmedido e VCalculado
            #residue = np.vstack([residue] * 4)
            #residue = np.repeat(residue, 4, axis=0)
            #print('bbb \n',residue)
            #print('ccc \n',residue.shape)
    
            #print('residue \n', residue.shape)
            self.calc_Y_jacobian()      # Calcula (dY/dσ_k) do Jacobiano
    
    
            # ***** Cálculo J = - Y_inv * (dY/ds_k) * (Y_inv * C) *****
            self.Calc_J(invVtempJ)
            #print('TempJ \n', self.TempJ.shape)
            
            # ***** Cálculo do termo 1a JT_W1_J *****
            W1=np.eye(self.TempJ.shape[0])
            #W1=np.eye(self.mymesh.NumberOfElements)
            JTW = np.dot(self.TempJ.T, W1)
            #JTW = JTW[self.mymesh.ElectrodeNodes]                        # pega somente valores dos eletrodos
            #print('W1 \n', W1.shape)
            #print('self.TempJ \n', self.TempJ.shape)
            #print('JTW \n', JTW)
            JTWJ = np.dot(JTW, self.TempJ)
            
            # ***** Cálculo do termo 1b Lambda^2 * LT_L *****
            #print('JTWJ \n', JTWJ)
            #print('L2 \n', L2)
            LTL =np.dot(L2.T, L2)
            #LTL = np.eye(4)
            #print('LTL \n', LTL.shape)
            termo_L = (Lambda**2)*LTL
            #print('termo_L \n', termo_L)
            
            # ***** Cálcula e inverte termo 1 -> (JTWJ + Lambda^2*LTL)^-1 *****
            firstTerm = JTWJ + termo_L
            #print('firstTerm \n', firstTerm)
            inv_firstTerm = np.linalg.inv(firstTerm)
            #print('inv_firstTerm \n', inv_firstTerm)
            
            # ***** Cálculo do termo 2a (JT_W1_residue) *****
            
            JTW_H = np.dot(JTW,residue)
            #print('JTW1 \n', JTW)
            #print('JTW_H \n', JTW_H)
            #print('residue \n', residue)
            
            
            # ***** Cálculo do termo 2b (Lambda^2 * LTL*residue) *****
            regTerm = (sigmaInicial - sigmaStar)* (Lambda**2)
            #print('regTerm \n', regTerm.shape)
            #print('sigmaInicial \n', sigmaInicial.shape)
            #print('sigmaStar \n', sigmaStar.shape)
            regTerm = regTerm.reshape(-1, 1)
            regTermC = np.repeat(regTerm, self.mymesh.NumberOfElectrodes, axis=1)
            
            #print('sigmaStar \n', sigmaStar.shape)
            #print('regTermC \n', regTermC.shape)
            regularization = np.dot(termo_L, regTermC)
            #regularization = termo_L
            
            
            # ***** Cálculo final do termo 2 *****
            secondTerm = JTW_H - regularization
            #print('secondTerm \n', secondTerm)
            
            
            # ***** Produto entre termo 1 e termo 2 *****
            deltaSigma = np.dot(inv_firstTerm, secondTerm)
            deltaSigma = deltaSigma[:, 0]
            #print('deltaSigma \n', deltaSigma)
            alphaDeltaSigama = -alpha*deltaSigma
            normaDelta = np.linalg.norm(alphaDeltaSigama)
            
            print(f'{itr} - normaDelta = {normaDelta}, normaResidue = {normaResidue} ')
            if normaDelta < 1e-9:    # Convergência atingida se a norma de
                                                    # delta_sigam < que  1e-6
              print(f'Convergência atingida após {itr} iterações.')
              #print('Vmedido_b \n', Vmedido_b)
              #print('Valor calculado \n',V_calc_b)
              print('Sigma k+1 \n', sigmaInicial)
              convergencia = True
              break                                   # interrompe o processo de iteração
            
            sigmaPlusOne = sigmaInicial + alphaDeltaSigama
            
            #print('sigmaPlusOne \n', sigmaPlusOne[:10])
            
            sigmaInicial = sigmaPlusOne
            #print(f'iteração = {itr}')
            #print('sigmaInicial \n', sigmaInicial[:5])
            #print(f'{itr} {sigmaInicial[:5]} ')
        print('sigmaInicial \n', sigmaInicial) 
        
        self.plotMSH(sigmaInicial)

        '''
        # 4. Gráfico tipo steam
        eixoX = np.arange(1, self.mymesh.NumberOfElements+1)   # 1 até 256
        plt.figure()
        #plt.stem(centros, valores)
        plt.plot(eixoX, sigmaInicial ,
                marker='None',
                label='$\sigma$ calculado')
        
        #plt.xlim(0, 1.01)
        plt.ylim(0.9, 1.5)
        #plt.xlabel('Posição [m]')
        #plt.ylabel('Condutividade σ')
        #plt.title('Distribuição de  σ nos elementos 1D')
        #plt.text(0.5, 0.87,f'Para α = {alpha_b} , λ = {lambda_b} e std={std}', ha='center', va='center',    plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        '''

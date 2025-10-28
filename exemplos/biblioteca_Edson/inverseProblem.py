# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 13:56:12 2025

@author: rodgr
"""

import numpy as np
import matplotlib.pyplot as plt
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
    # Esta função calcula a matriz global de condutividade da malha de elementos
    # finitos
      
    ###############################################################################
    def calc_Y_global_1D(self, sigma):               # função para calcular a matriz global
      self.Y_temp = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))
      for i in range(0, self.mymesh.NumberOfElements):                 #  laço para montar a matriz global

            node_l = int(self.mymesh.msh_topology[i][0]) +1        # pegar dados da matriz elementos
            node_m = int(self.mymesh.msh_topology[i][1]) +1

            coord_1 = [self.mymesh.Coordinates[node_l-1][0]]                # coordenada x do ponto 1
            coord_2 = [self.mymesh.Coordinates[node_m-1][0]]                # coordenada x do ponto 2
        #######################################################################
        # Calcula a matriz local de condutividade para V_calc
        #        A_i * sigma_i                                      
        # Y_local = ------------- * [[1, -1], [-1, 1]]             
        #             L                                             
        #######################################################################
            Y_local = ((self.mymesh.altura1D*sigma[i])/(coord_2[0]-coord_1[0]))*np.array([[1, -1], [-1, 1]])
            #Y_local = (self.calc_Y_local_1D(coord_1[0], coord_2[0], self.mymesh.altura1D, sigma[i]) )                                                                             # calcula a matriz local
            print('Y_local \n',Y_local)

            self.Y_temp[node_l - 1, node_l - 1] += Y_local[0, 0] # monta a matriz global
            self.Y_temp[node_m - 1, node_l - 1] += Y_local[0, 1]
            self.Y_temp[node_l - 1, node_m - 1] += Y_local[1, 0]
            self.Y_temp[node_m - 1, node_m - 1] += Y_local[1, 1]
            #print('self.Y_temp \n',self.Y_temp[:5])
      return self.Y_temp
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

    def calc_Y_jacobian(self):   # função para calcular a matriz global
        self.Y_jacobian = np.zeros((int(self.mymesh.NumberOfNodes),int(self.mymesh.NumberOfNodes)))

        for i in range(0, self.mymesh.NumberOfElements):        # laço para montar a matriz global
            #print('nro_elementos_J',i)
            node_l = int(self.mymesh.msh_topology[i][0]) +1       # pegar dados da matriz elementos
            node_m = int(self.mymesh.msh_topology[i][1]) +1
            print(node_l, node_m)
            coord_1 = [self.mymesh.Coordinates[node_l-1][0]]    # coordenada x do ponto 1
            coord_2 = [self.mymesh.Coordinates[node_m-1][0]]    # coordenada x do ponto 2
            
            
            Y_local_c =((self.mymesh.altura1D/(coord_2[0]-coord_1[0]))*np.array([[1, -1], [-1, 1]]))  # derivada parcial da matriz local
            #print('Y_local_c',Y_local_c)
            
            self.Y_jacobian[node_l - 1, node_l - 1] += Y_local_c[0, 0] # monta a
            self.Y_jacobian[node_m - 1, node_l - 1] += Y_local_c[0, 1] # M global
            self.Y_jacobian[node_l - 1, node_m - 1] += Y_local_c[1, 0]
            self.Y_jacobian[node_m - 1, node_m - 1] += Y_local_c[1, 1]
            self.length[i] = coord_2[0] - coord_1[0]  # calc vetor comprimento
            #print('Y_jacobian',self.Y_jacobian)
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
    # Essa função monta a matriz Jacobiana
    ###############################################################################
    def calc_Jacobian(self, inv_KJacob ):  # Calc jacobiano em função de sigma


      limite = 1e-12                               # limite inferior para jacobiano
      YI_corrente = inv_KJacob @ self.corrente
      for i in range(self.mymesh.NumberOfElements):         # loop para montar a matriz jacobiana 1D
        termo_2 = (self.mymesh.altura1D/self.length[i])       # calc termo da matriz da
                                                    # derivada parcial do jacobiano

        self.jacob_2[i, i] = termo_2
        self.jacob_2[i, (i+1)] = -termo_2
        self.jacob_2[(i+1), i] = -termo_2
        self.jacob_2[(i+1), (i+1)] = termo_2

        self.jacob_1[:, i] = inv_KJacob @ (self.jacob_2.T @ YI_corrente) # calc matriz jacobiana x
                                                          # vetor de corrente
        self.jacob_1[np.abs(self.jacob_1) < limite] = 0      # se valor < 1e-12 forçar zero
                                                     # na matriz jacobiana


        self.jacob_2[:] = 0.0                                                # zerar matriz
      return self.jacob_1[self.mymesh.ElectrodeNodes, :]     # retorna matriz jacobiana x vetor de corrente

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
        print('centroids_1D',centroids_1D)
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
    ###############################################################################
    # Essa função plota a curva de valores reais de sigma versus valores de sigma
    # calculados nas iterações
    ###############################################################################
    def plotar_grafico(self, sigma_inicial_b,sigma_real_b, noh_medidos, ylim=(0.0, 0.60), figsize=(16, 4) ):
        # Coordenadas x dos nós
        x_coords = self.mymesh.Coordinates.flatten()
        topologia_b = self.mymesh.msh_topology

        centros = np.mean(x_coords[topologia_b], axis=1)
        #centros = (x_coords[:-1] + x_coords[1:]) / 2
        valores = sigma_inicial_b  # ou delta_sig_b
        valores_real = sigma_real_b  # ou delta_sig_b
        noh_medidos[-1] -= 1
        pos_medidas = [centros[i] for i in noh_medidos]
        pos_valores_med = [valores[i] for i in noh_medidos]
        #pos_valores_real = [valores_real[i] for i in noh_medidos]
        # 4. Gráfico tipo steam
        #plt.figure(figsize=figsize)
        #plt.stem(centros, valores)
        #plt.plot(centros, valores , marker='None', label='$\sigma$ calculado')
        #plt.plot(centros, sigma_real_b, marker='None', linestyle=':', label='$\sigma$ real' )
        #plt.plot(pos_medidas, pos_valores_real, marker='x', linestyle='None', markersize=10, color='red',        label='Ptos medidos')
        plt.xlim(0, 1.01)
        plt.ylim(ylim)
        plt.xlabel('Posição [m]')
        plt.ylabel('Condutividade σ')
        plt.title('Distribuição de  σ nos elementos 1D')
        #plt.text(0.5, 0.87, f'Para α = {alpha_b} , λ = {lambda_b} e std={std}', ha='center', va='center',transform=plt.gcf().transFigure)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    ###############################################################################
    
    def solve(self, sigma_inicial, V_measured, meus_sigmas, max_iter=200):
        lista_i = []                                        # Lista armazenar iterações
        lista_plotar = []                                      # Lista Valores de sigma
        self.calc_Y_jacobian()
        print('Y_jacobiano \n', self.Y_jacobian)
        KJacobian = self.apply_boundary_conditions(self.Y_jacobian)
        print('KJacobian \n', KJacobian)
        #print('self.KJacobian', self.KJacobian)
        for i in range(max_iter):
            self.Y_Vcalc = self.calc_Y_global_1D(sigma_inicial)
            print('self.Y_Vcalc \n', self.Y_Vcalc)
            KVcalc = self.apply_boundary_conditions(self.Y_Vcalc)
            print('KVcalc \n', KVcalc)
            
            invKVcalc = np.linalg.inv(KVcalc)
    
            V_calc = np.dot(invKVcalc, self.vetor_corrente_cond_contorno)
            V_calc_noh = V_calc[self.mymesh.ElectrodeNodes]
            print('V_calc_noh \n', V_calc_noh)
            residue = V_measured - V_calc_noh        # calc dif entre Vmedido e VCalculado
            print('residue \n', residue)
            
            inv_KJacobian = np.linalg.inv(KJacobian)
            print('inv_KJacobian \n', inv_KJacobian)
            Jacobian = self.calc_Jacobian(inv_KJacobian)
            print('Jacobian', Jacobian)
            #Centroid=self.Elements.CalcCentroid()
            #print('centroid', Centroid)
            x_coords_b = self.mymesh.Coordinates
            topologia_bc = self.mymesh.msh_topology
            centroids_1D = np.mean(x_coords_b[topologia_bc], axis=1)
            centroids_1D = centroids_1D[:,0]
            
            print('centroids_1D', centroids_1D)
            L2 = self.calc_L2_gauss_1D(centroids_1D)
            print('L2',L2)
            
            delta_sig = self.calc_delta_sigma(Jacobian, residue, L2, sigma_inicial, alpha=0.01, Lambda=0.006)
            print('delta_sig \n', delta_sig)
            
            plotar = np.linalg.norm(delta_sig)  # calc norma vetor delta_sigama para plot
            lista_i.append(i)                             # Armazena o índice da iteração
            lista_plotar.append(plotar)                  # Armazena o valor a ser plotado
    
            if np.linalg.norm(delta_sig) < 1e-3:    # Convergência atingida se a norma de
                                                    # delta_sigam < que  1e-6
              print(f'Convergência atingida após {i+1} iterações.')
              #print('Vmedido_b \n', Vmedido_b)
              #print('Valor calculado \n',V_calc_b)
              print('Sigma k+1 \n', sigma_inicial)
              convergencia = True
              break                                   # interrompe o processo de iteração
    
            Sig_kMais1 = sigma_inicial + delta_sig           # ajusta vetor sigma k+1
            sigma_inicial = Sig_kMais1            # armazena vetor sigma k+1 anterior
        self.plotar_iteracoes(lista_i, lista_plotar)
        
        
        self.plotar_grafico(sigma_inicial, meus_sigmas, self.mymesh.ElectrodeNodes)
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 14:17:52 2025

@author: rodgr
"""
# criar pull request

import numpy as np
import mesh
import elements
import gmsh
import sys
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import datetime
#import subprocess
#import os


class forward_problem: 
    def __init__(self, mymesh, V_imposto=None, Pcorrente=None, SkipPattern=None, VirtualNode = False, I =1.0e-3):
        if not hasattr(mymesh, "KGlobal"): # verifica se o objeto mymesh tem um atributo chamado KGlobal.
            raise TypeError("Parâmetro incorreto: mymesh.")

        self.mymesh = mymesh
        self.Vmedido = None
        self.Yinversa = None

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


    ###############################################################################
    # Esta função aplica condições de contorno conforme exemplo abaixo:
    #
    # [ 1   0          0          0          0        ] [u1]   [u1]    
    # [ 0  k1+k2+k3    -k3        0         -k2       ] [u2] = [F2 +k1*u1] 
    # [ 0   -k3      k3+k5+k4     -k5        0        ] [u3]   [F3 +k4*u4] 
    # [ 0   0          -k5      k2+k6       -k6       ] [u4]   [F4 +k6*u5] 
    # [ 0   0           0          0         1        ] [u5]   [u5]
    ###############################################################################  
    def apply_boundary_conditions(self):
        self.vetor_corrente_cond_contorno = self.corrente[:].copy()
        #print(f'Vetor de corrente: \n {self.vetor_corrente_cond_contorno}')

        self.KGlobal = self.mymesh.KGlobal.copy()       # Criar matriz solução
        print(f' AppCC self.KGlobal: \n {self.KGlobal}')
        # atualiza vetor de correntes:
        for [noh_cond_contorno,valor_cond_contorno] in self.V_imposto:   # Necessário quando valor imposto é diferente de zero
            for i in range(0, self.mymesh.NumberOfNodes):                    # corrige matriz de corrente
                self.vetor_corrente_cond_contorno[i] = (
                    self.vetor_corrente_cond_contorno[i] - \
                    self.KGlobal[i][noh_cond_contorno]*valor_cond_contorno
                )

        for [noh_cond_contorno_b,valor_cond_contorno] in self.V_imposto:
            self.vetor_corrente_cond_contorno[noh_cond_contorno] = valor_cond_contorno  # coloca valor de contorno conhecido
        
        # atualiza matriz KGlobal:
        for [noh_cond_contorno,valor_cond_contorno] in self.V_imposto:
          for k in range(0, self.mymesh.NumberOfNodes):                # laço para zerar linha e coluna
              self.KGlobal[noh_cond_contorno][k] = 0
              self.KGlobal[k][noh_cond_contorno] = 0
        self.KGlobal[noh_cond_contorno][noh_cond_contorno] = 1
    ###############################################################################
    def plotMSH(self, sigma, iteration = None, save = False):

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
            tpc = ax.tripcolor(
                triang,
                facecolors=sigma[:len(elems_2D)],
                edgecolors='k',
                cmap='Blues'
            )
            fig.colorbar(tpc, ax=ax, label='σ (Conductivity)')
            if save == True:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M")
                ax.set_title(f"Conductivity Real (σ) ", fontsize=12)
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
        
        plt.show()    

    def Solve(self, forceKGolbalCalc=False):
        #self.plotMSH(self.mymesh.sigma_vec)
        #print('self.mymesh.sigma_vec',self.mymesh.sigma_vec)
        if (self.mymesh.KGlobal is None) or (forceKGolbalCalc):
            self.mymesh.CalcKGlobal()
        
        self.apply_boundary_conditions()
        print('solve self.KGlobal \n',self.KGlobal)
        self.Yinversa = np.linalg.inv(self.KGlobal)

        self.Vmedido = np.dot(self.Yinversa, self.vetor_corrente_cond_contorno)
        print(f' Tensões medidas em todos os nós \n {self.Vmedido})')
        
        print('solve vetor_corrente_cond_contorno \n', self.vetor_corrente_cond_contorno)


        self.Vmedido_eletrodos = self.Vmedido[self.mymesh.ElectrodeNodes]
        #print(f' Tensões nos eletrodos \n {self.Vmedido_eletrodos})')


    ###############################################################################
    # Esta função cria arquivos .pos (Post-Processing) apara vizulização no Gmsh.
    #  opções: 'rho', 'sigma', 'V', 'I'
    ###############################################################################
    #def criar_arquivo_pos_2D(self, matriz_coordenadas, matriz_topologia, N_padraoCC, Post_Processing, nome_arquivo):
    def criar_arquivo_pos_2D(self, Post_Processing, nome_arquivo):

        matriz_coordenadas = self.mymesh.Coordinates
        matriz_topologia = self.mymesh.msh_topology
        dim = len(matriz_topologia[1])

        if dim == 2:
            self.n_PCurrent = 1
            header = "V(x) 1D"
            Sx = "SL"
            vetor_Post_Processing = Post_Processing
            arquivo = open('../../malhasPOS/'+ nome_arquivo + str(0) + '.pos', 'w+')
            arquivo.writelines('View "' + header + '" { \n')
            for i in range(0, len(matriz_topologia)):
            
                node_l = int(matriz_topologia[i][0]) # pegar dados do dataframe
                node_m = int(matriz_topologia[i][1])
                
                arquivo.writelines(Sx + '(')
                
                
                arquivo.writelines([str(matriz_coordenadas[node_l][0]),',', str(matriz_coordenadas[node_l][1]),',0,'])
                arquivo.writelines([str(matriz_coordenadas[node_m][0]),',', str(matriz_coordenadas[node_m][1]),',0'])
                
                arquivo.writelines(u')')
                arquivo.writelines(u'{')
                arquivo.writelines([str(Post_Processing[node_l]),',', str(Post_Processing[node_m])])
                arquivo.writelines(u'};')
                arquivo.writelines(u'\n')
            arquivo.writelines(u'};')
            
        else:
            header = "A list-based view"
            Sx = "ST"
            self.n_PCurrent = self.mymesh.NumberOfElectrodes
        

            for n_ele in range(self.n_PCurrent):
                vetor_Post_Processing = Post_Processing[:, n_ele]
                    
                arquivo = open('../../malhasPOS/'+ nome_arquivo + str(n_ele) + '.pos', 'w+')
                arquivo.writelines('View "' + header + '" { \n')
                
                for i in range(0, len(matriz_topologia)):
                
                    node_l = int(matriz_topologia[i][0]) # pegar dados do dataframe
                    node_m = int(matriz_topologia[i][1])
                    node_n = int(matriz_topologia[i][2]) 
                    arquivo.writelines(Sx + '(')
                    
                    
                    arquivo.writelines([str(matriz_coordenadas[node_l][0]),',', 
                                        str(matriz_coordenadas[node_l][1]),',0,'])
                    arquivo.writelines([str(matriz_coordenadas[node_m][0]),',', 
                                        str(matriz_coordenadas[node_m][1]),',0,'])
                    arquivo.writelines([str(matriz_coordenadas[node_n][0]),',', 
                                        str(matriz_coordenadas[node_n][1]),',0'])
                    
                    arquivo.writelines(u')')
                    arquivo.writelines(u'{')
    
                    arquivo.writelines([str(vetor_Post_Processing[node_l]),',', 
                                        str(vetor_Post_Processing[node_m]),',', 
                                        str(vetor_Post_Processing[node_n])])
            
                    arquivo.writelines(u'};')
                    arquivo.writelines(u'\n')
                arquivo.writelines(u'};')
            
        arquivo.close()  
##############################################################################    



###############################################################################
# Esta função abre arquivos .pos (Post-Processing), criados pela função
# 'criar_arquivo_pos'  para vizulização no Gmsh.
#  
###############################################################################
    def abrir_Gmsh_pos(self,nome_arquivo, runGmsh = False):
        gmsh.initialize()
        for pos in range(self.n_PCurrent):
            # Inicialize o Gmsh
            #gmsh.initialize()
            
            # Carregue o arquivo .geo
            #gmsh.open('D:/GIT_EIT_2D/EIT_2D/malhasPOS/' + nome_arquivo + str(pos) + '.pos') 
            gmsh.open('../../malhasPOS/' + nome_arquivo + str(pos) + '.pos')
            #gmsh.open( nome_arquivo + str(pos) + '.pos')
            #gmsh.open('../../malhasPOS/' + nome_arquivo + '.pos') 
            #gmsh.View[0].IntervalsType = 2;
            #gmsh.option.setNumber("View["+str(pos)+"].IntervalsType", 1)
            #gmsh.option.setNumber("View["+str(pos)+"].NbIso", 500)
           
            
        if runGmsh and '-nopopup' not in sys.argv:
            gmsh.fltk.run()
        gmsh.finalize()
###############################################################################
        

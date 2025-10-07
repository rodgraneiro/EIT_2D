# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 14:17:52 2025

@author: rodgr
"""

import numpy as np
import mesh
import elements
import gmsh
import sys


class forward_problem: 
    def __init__(self, mymesh: "mesh.MyMesh", Pcorrente=None):
        if not hasattr(mymesh, "KGlobal"): # verifica se o objeto mymesh tem um atributo chamado KGlobal.
            raise TypeError("Parâmetro incorreto: mymesh.")
        if mymesh.KGlobal is None:
            raise TypeError("mymesh não tem atributo KGlobal (calcule CalcKGlobal()).")
        if Pcorrente is None:
            raise TypeError("não tem atributo 'Pcorrente' (defina um vetor de correntes).")

        self.mymesh = mymesh
        self.KGlobal = np.asarray(mymesh.KGlobal, dtype=float)
        self.corrente = Pcorrente
        self.Vmedido = None


###############################################################################
# Esta função aplica condições de contorno conforme exemplo abaixo:
#
# [ 1   0          0          0          0        ] [u1]   [u1]    
# [ 0  k1+k2+k3    -k3        0         -k2       ] [u2] = [F2 +k1*u1] 
# [ 0   -k3      k3+k5+k4     -k5        0        ] [u3]   [F3 +k4*u4] 
# [ 0   0          -k5      k2+k6       -k6       ] [u4]   [F4 +k6*u5] 
# [ 0   0           0          0         1        ] [u5]   [u5]
###############################################################################  
    def apply_boundary_conditions(self, V_imposto ):
       
        print(f' KGlobal: \n {self.KGlobal}')
        vetor_corrente_cond_contorno = self.corrente[:].copy()
        print(f'Vetor de corrente: \n {vetor_corrente_cond_contorno}')
          
        n_nodes = vetor_corrente_cond_contorno.shape[0]
        for [noh_cond_contorno,valor_cond_contorno] in V_imposto:
          for i in range(0, n_nodes):                    # corrige matriz de corrente
            vetor_corrente_cond_contorno[i] = (
                vetor_corrente_cond_contorno[i] - \
                self.KGlobal[i][noh_cond_contorno]*valor_cond_contorno
            )
        for [noh_cond_contorno_b,valor_cond_contorno] in V_imposto:
          vetor_corrente_cond_contorno[noh_cond_contorno] = (
              valor_cond_contorno )            # coloca valor de contorno conhecido
    
        self.KGlobal_cond_contorno = self.KGlobal[:].copy()                        # Criar matriz solução
        for [noh_cond_contorno,valor_cond_contorno] in V_imposto:
          for k in range(0, n_nodes):                # laço para zerar linha e coluna
              self.KGlobal_cond_contorno[noh_cond_contorno][k] = 0
              self.KGlobal_cond_contorno[k][noh_cond_contorno] = 0
    
        self.KGlobal_cond_contorno[noh_cond_contorno][noh_cond_contorno] = 1
        print(f'KGlobal_cond_contorno \n {self.KGlobal_cond_contorno})')
        return vetor_corrente_cond_contorno, self.KGlobal_cond_contorno

###############################################################################

###############################################################################
# Essa função calcula o valor medido conforme equação abaixo
#
# { V_medido } = [ Y_cond_contorno ]^(-1) * { C_cond_contorno }          
###############################################################################
    def solucao_Vmedido(self, vetor_correnteVM,
                          Y_cond_contorno_VM):        # calc valor observado/medido
      Yinversa = np.linalg.inv(Y_cond_contorno_VM)
      Vmedido = np.dot(Yinversa,
                         vetor_correnteVM)
      return Vmedido                             # retorna valor observado/medido
###############################################################################

###############################################################################
# Esta função cria arquivos .pos (Post-Processing) apara vizulização no Gmsh.
#  
###############################################################################
    def criar_arquivo_pos_2D(self, matriz_coordenadas, matriz_topologia, 
                          n_eletrodos, V_sol, nome_arquivo):
        #for kk in range(n_eletrodos):
        for kk in range(1):

            
            # 9) Gerar arquivo POS
            
            #nome_arquivo = 'frame_1_5_11.pos'
            #arquivo = open('D:/GIT_EIT_2D/EIT_2D/malhasPOS/'+ nome_arquivo + str(kk) + '.pos', 'w+')
            arquivo = open('../../malhasPOS/'+ nome_arquivo + str(kk) + '.pos', 'w+')
            arquivo.writelines(u'View "A list-based view" { \n')
            
            for i in range(0, len(matriz_topologia)):
            
                node_l = int(matriz_topologia[i][0]) # pegar dados do dataframe
                node_m = int(matriz_topologia[i][1])
                node_n = int(matriz_topologia[i][2])
            
                arquivo.writelines(u'ST(')
              
                arquivo.writelines([str(matriz_coordenadas[node_l][0]),',', 
                                    str(matriz_coordenadas[node_l][1]),',0,'])
                arquivo.writelines([str(matriz_coordenadas[node_m][0]),',', 
                                    str(matriz_coordenadas[node_m][1]),',0,'])
                arquivo.writelines([str(matriz_coordenadas[node_n][0]),',', 
                                    str(matriz_coordenadas[node_n][1]),',0'])
                
                arquivo.writelines(u')')
                arquivo.writelines(u'{')
                arquivo.writelines([str(V_sol[node_l]),',', 
                                    str(V_sol[node_m]),',', 
                                    str(V_sol[node_n])])
               
            
                arquivo.writelines(u'};')
                arquivo.writelines(u'\n')
            arquivo.writelines(u'};')
            
            arquivo.close()  
###############################################################################    

###############################################################################
# Esta função abre arquivos .pos (Post-Processing), criados pela função
# 'criar_arquivo_pos'  para vizulização no Gmsh.
#  
###############################################################################
    def abrir_Gmsh_pos(self,nome_arquivo, n_eletrodos):
        for pos in range(n_eletrodos):
            # Inicialize o Gmsh
            gmsh.initialize()
            
            # Carregue o arquivo .geo
            #gmsh.open('D:/GIT_EIT_2D/EIT_2D/malhasPOS/' + nome_arquivo + str(pos) + '.pos') #open('TesteBanana0.pos')
            gmsh.open('../../malhasPOS/' + nome_arquivo + str(pos) + '.pos') #open('TesteBanana0.pos')
            #gmsh.View[0].IntervalsType = 2;
            #gmsh.option.setNumber("View["+str(pos)+"].IntervalsType", 1)
            #gmsh.option.setNumber("View["+str(pos)+"].NbIso", 500)
            # Obtenha as entidades do modelo
            #entities = gmsh.model.getEntities()
            #gmsh.model.geo.synchronize()
            
            #gmsh.model.mesh.generate(3)
            #gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
            #gmsh.write("poligono_2D_Rho_py.msh")
            #gmsh.write("poligono_2D_Rho.geo_unrolled")
            
            
            
            #gmsh.model.mesh.generate(3)
            # Feche o arquivo .geo
            
        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()
        gmsh.finalize()
###############################################################################      
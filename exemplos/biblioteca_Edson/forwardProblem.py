# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 14:17:52 2025

@author: rodgr
"""

import numpy as np
import mesh
import elements


class forward_problem: 
    def __init__(self, mymesh: "mesh.MyMesh"):
        if not hasattr(mymesh, "KGlobal"): # verifica se o objeto mymesh tem um atributo chamado KGlobal.
            raise TypeError("mymesh não tem atributo KGlobal (instancie MyMesh e calcule CalcKGlobal()).")
        if not hasattr(mymesh, "corrente"): # verifica se o objeto mymesh tem um atributo chamado corrente.
            raise TypeError("mymesh não tem atributo 'corrente' (defina um vetor de correntes).")

        self.mymesh = mymesh
        self.KGlobal = np.asarray(mymesh.KGlobal, dtype=float)
        self.corrente = np.asarray(mymesh.corrente, dtype=float)
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
  
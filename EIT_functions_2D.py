###############################################################################
# Esse arquivo tem as funções necessáras para rodar o progrma EIT_2D_main.py
# 
###############################################################################


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
#import meshio
from matplotlib.collections import LineCollection
import time




###############################################################################
# Essa função plota  a distribuição de condutividade da malha de Elementos
# Finitos
###############################################################################
def plot_mesh(matriz_coordenadas, matriz_topologia, sigma_real):
    
    # Plot Grafico da condutividade originl do fantoma            
    x, y = matriz_coordenadas[:, 0], matriz_coordenadas[:, 1]                                                 
    triang = tri.Triangulation(x, y, matriz_topologia)
    fig, ax = plt.subplots(figsize=(6, 5))
    tpc = ax.tripcolor(triang, facecolors=sigma_real, edgecolors='k',cmap='Blues' )
    fig.colorbar(tpc, ax=ax, label='sigma (condutividade)')
    ax.set_title("Condutividade Real (sigma)")
    plt.xlabel("[m]")
    plt.ylabel("[m]")
    plt.tight_layout()
    plt.show()

###############################################################################
# Esta função aplica condições de contorno conforme exemplo abaixo:
#
# [ 1   0          0          0          0        ] [u1]   [u1]    
# [ 0  k1+k2+k3    -k3        0         -k2       ] [u2] = [F2 +k1*u1] 
# [ 0   -k3      k3+k5+k4     -k5        0        ] [u3]   [F3 +k4*u4] 
# [ 0   0          -k5      k2+k6       -k6       ] [u4]   [F4 +k6*u5] 
# [ 0   0           0          0         1        ] [u5]   [u5]
###############################################################################  
def aplica_cond_contorno(I_CC,
                         Y,
                         n_nodes,
                         V_imposto):           # Aplica condições de contorno

  I_CC_cond_contorno = I_CC[:].copy()

  for [noh_cond_contorno,valor_cond_contorno] in V_imposto:
    for i in range(0, n_nodes):                    # corrige matriz de corrente
      I_CC_cond_contorno[i] = (
          I_CC_cond_contorno[i] - \
          Y[i][noh_cond_contorno]*valor_cond_contorno
      )

  for [noh_cond_contorno,valor_cond_contorno] in V_imposto:
    I_CC_cond_contorno[noh_cond_contorno] = (
        valor_cond_contorno )            # coloca valor de contorno conhecido

  Y_cond_contorno = Y[:].copy()                        # Criar matriz solução
  for [noh_cond_contorno,valor_cond_contorno] in V_imposto:
    for k in range(0, n_nodes):                # laço para zerar linha e coluna
        Y_cond_contorno[noh_cond_contorno][k] = 0
        Y_cond_contorno[k][noh_cond_contorno] = 0

    Y_cond_contorno[noh_cond_contorno][noh_cond_contorno] = 1
  return Y_cond_contorno, I_CC_cond_contorno
###############################################################################

###############################################################################
# Esta função calcula a matriz local de condutividade da malha de elementos
# finitos
#
#            sigma_i                                      
# Y_local = ------------- * [matriz 2d ***completar***]             
#             4*A_i                                             
###############################################################################
def calc_Y_local_triangle(coords, sigma):
    x = coords[:, 0]
    y = coords[:, 1]

    triangulo = np.array([ [1, x[0],  y[0]], [1, x[1],  y[1]], [1, x[2],  y[2]]], dtype=np.float64)
    area_triangulo = abs((np.linalg.det(triangulo))/2)
 
    C_11 =  ( y[1]-y[2])**2 + (x[2]-x[1])**2
    C_12 =  ( y[1]-y[2])*(y[2]-y[0])+(x[2]-x[1])*(x[0]-x[2])               #A_21 = A_12
    C_13 =  ( y[1]-y[2])*(y[0]-y[1])+(x[2]-x[1])*(x[1]-x[0])                 #A_31 = A_13
    C_22 =  (y[2]-y[0])**2 + (x[0]-x[2])**2                       #A_32 = A_23
    C_23 =  (y[2]-y[0])*(y[1]- y[1])+(x[0]-x[2])*(x[1]-x[0])                   # Na tese está: (Bt_m*Bt_o)+(Gm_m*Gm_o)+(Dt_m*Dt_o)
    C_33 =  (y[0]- y[1])**2 + (x[1]-x[0])**2  

    
    matriz_local = (sigma/(4*area_triangulo))*np.array([[C_11, C_12, C_13], 
                                    [C_12, C_22, C_23],      
                                    [C_13, C_23, C_33]       
                                    ])

    return matriz_local
###############################################################################




###############################################################################
# Esta função calcula a matriz global de condutividade da malha de elementos
# finitos
  
###############################################################################
def montar_Y_global(matriz_coordenadas, matriz_topologia, n_nodes, sigma_elem):
    Y = np.zeros((n_nodes, n_nodes))
    for e, conn in enumerate(matriz_topologia):
        coords_elem = matriz_coordenadas[conn]
        Y_loc = calc_Y_local_triangle(coords_elem, sigma_elem[e])
        for i in range(3):
            for j in range(3):
                Y[conn[i], conn[j]] += Y_loc[i, j]
    return Y
###############################################################################
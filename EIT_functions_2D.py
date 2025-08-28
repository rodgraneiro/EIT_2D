###############################################################################
# Esse arquivo tem as funções necessáras para rodar o progrma EIT_2D_main.py
# 
###############################################################################


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
#import meshio
import gmsh
import sys
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

###############################################################################
# Esta função cria arquivos .pos (Post-Processing) apara vizulização no Gmsh.
#  
###############################################################################
def criar_arquivo_pos(matriz_coordenadas, matriz_topologia, 
                      n_eletrodos, V_sol, nome_arquivo):
    for kk in range(n_eletrodos):
        #nome_arquivo = 'TesteBanana'
        
        
        # 9) Gerar arquivo POS
        
        #nome_arquivo = 'frame_1_5_11.pos'
        arquivo = open('D:/GIT_EIT_2D/EIT_2D/malhasPOS/'+ nome_arquivo + str(kk) + '.pos', 'w+')
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
            arquivo.writelines([str(V_sol[node_l][kk]),',', 
                                str(V_sol[node_m][kk]),',', 
                                str(V_sol[node_n][kk])])
           
        
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
def abrir_Gmsh_pos(nome_arquivo, n_eletrodos):
    for pos in range(n_eletrodos):
        # Inicialize o Gmsh
        gmsh.initialize()
        
        # Carregue o arquivo .geo
        gmsh.open('D:/GIT_EIT_2D/EIT_2D/malhasPOS/' + nome_arquivo + str(pos) + '.pos') #open('TesteBanana0.pos')
        #gmsh.View[0].IntervalsType = 2;
        gmsh.option.setNumber("View["+str(pos)+"].IntervalsType", 1)
        gmsh.option.setNumber("View["+str(pos)+"].NbIso", 500)
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


###############################################################################
# Esta função cria a matriz com o padrão de correntes
# 'criar_arquivo_pos'  para vizulização no Gmsh.
#  
###############################################################################
def matriz_corrente(linhas,colunas, s, I=1e-3):
    MtrCC = np.zeros((linhas, colunas))
    k = np.arange(colunas)
    MtrCC[k, k] = I
    MtrCC[(k + s) % colunas, k] = -I
    return MtrCC
###############################################################################
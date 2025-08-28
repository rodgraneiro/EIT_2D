import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.tri as tri
from scipy.sparse import lil_matrix
from matplotlib.collections import LineCollection
import meshio
import gmsh
#import shutil
import time
import importlib



import time
import sys

#from EIT_mesh_2D import files_mesh
from EIT_functions_2D import plot_mesh, matriz_corrente
from EIT_functions_2D import plot_mesh, aplica_cond_contorno
from EIT_functions_2D import calc_Y_local_triangle, montar_Y_global 
from EIT_functions_2D import criar_arquivo_pos, abrir_Gmsh_pos



###############################################################################
########################## CARREGA MALHA 1D ###################################
###############################################################################

print('\n ##### NOME DA MALHA #####')
nome_malha = input("Digite o nome da malha Gmsh (sem .msh): ")
malha_msh = meshio.read('D:/GIT_EIT_2D/EIT_2D/malhasMSH/' + nome_malha + '.msh')  
  
arquivoConfig = importlib.import_module('Configs.'+ nome_malha)


#nome_arquivo = arquivoConfig.nome_arquivo
n_eletrodos = arquivoConfig.n_eletrodos
pto_central= n_eletrodos + 1 
V_imposto = [[n_eletrodos, 0.0]]

puloMaximo = int(n_eletrodos // 2)


print('\n ##### VALOR DA CONDUTIVIDADE DA BASE (CUBA) - type float #####')
sigma_phyfical_1 = float(input("Digite a condutividade da base (ex: 0.1)): "))

print('\n ##### VALOR DA CONDUTIVIDADE DA anomalia - type float #####')
sigma_phyfical_2 = float(input("Digite a condutividade da anomalia ( ex: 0.75)): "))

print("\n ##### O VALOR MÁXIMO DO PULO DO PADRÃO DE CORRENTE É :", puloMaximo, '- type int #####')
puloPadraoCC = int(input("Digite o pulo do padrão de corrente (type int ex: 2)): "))



matriz_coordenadas = malha_msh.points[:, :2]       # Monta matriz coordenadas
matriz_topologia = malha_msh.cells_dict["triangle"]      # Monta Matriz topologia
physical = np.array(malha_msh.cell_data["gmsh:physical"])
n_nodes = matriz_coordenadas.shape[0]                             # Nr de nós
n_elements = matriz_topologia.shape[0]                      # Nr de elementos

I_matrix = matriz_corrente(n_nodes, n_eletrodos, puloPadraoCC, I=1e-3)




sigma_real = np.zeros(n_elements)

for phy in range(n_elements):
    if physical[0, phy] == 1:
        sigma_real[phy] = sigma_phyfical_1
    elif physical[0, phy] == 2:
        sigma_real[phy] = sigma_phyfical_2
    else:
        sigma_real[phy] = 0.0  # valor padrão, se aparecer outro grupo

plot_mesh(matriz_coordenadas, matriz_topologia, sigma_real)

# === Resolução do problema direto ===
nos_contorno = np.arange(n_eletrodos)
V_medido = []

#for k in range(n_eletrodos):
for k in range(n_eletrodos):
    I_vec = np.zeros(n_nodes)
    I_vec = I_matrix[k]
    print(k, I_matrix[k])
    Y_global = montar_Y_global(matriz_coordenadas, matriz_topologia, n_nodes, sigma_real)
    Y_mod, I_mod = aplica_cond_contorno(I_matrix, Y_global, n_nodes, V_imposto)
    #V_sol = np.linalg.solve(Y_mod, I_mod)
    inv_Y_mod = np.linalg.inv(Y_mod)
    V_sol = np.dot(inv_Y_mod, I_mod)
    V_medido.append(V_sol)
    #print('matriz_coordenadas', matriz_coordenadas[k])
    

V_medido = np.array(V_medido[0])
V_medido = np.array(V_medido[0:n_eletrodos, 0:n_eletrodos])



np.set_printoptions(formatter={'float': '{:0.2e}'.format}, linewidth=300)

print('V_medido \n', V_medido)
print('V_sol', V_sol)

nome_arquivo = 'ParaVernoGmsh'

criar_arquivo_pos(matriz_coordenadas, matriz_topologia, n_eletrodos, V_sol, nome_arquivo)

abrir_Gmsh_pos(nome_arquivo, n_eletrodos)



###############################################################################
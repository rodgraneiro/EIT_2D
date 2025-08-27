###############################################################################
# Programa para reconstrução de Tomografia por Impedância Elétrica 
# (Anisotrópica)
###############################################################################



import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.tri as tri
from scipy.sparse import lil_matrix
from matplotlib.collections import LineCollection
import meshio

import time




import time
import sys

from EIT_mesh_2D import files_mesh
from EIT_functions_2D import plot_mesh, aplica_cond_contorno
from EIT_functions_2D import calc_Y_local_triangle, montar_Y_global
from EIT_functions_2D import criar_arquivo_pos, abrir_Gmsh_pos

###############################################################################
########################## CARREGA MALHA 1D ###################################
###############################################################################

#opcao = input('Escolha a malha(1, 2): ')
opcao = '1'




if opcao not in ['1', '2']:             # verifica se existe a opçaõ
    raise ValueError("Opção inválida.")

caminho, sigma_real, V_imposto, I_matrix, n_eletrodos = files_mesh(opcao)

malha_msh = meshio.read(caminho)                            # Lê o arquivo .msh

matriz_coordenadas = malha_msh.points[:, :2]       # Monta matriz coordenadas
matriz_topologia = malha_msh.cells_dict["triangle"]      # Monta Matriz topologia


n_nodes = matriz_coordenadas.shape[0]                             # Nr de nós
n_elements = matriz_topologia.shape[0]                      # Nr de elementos


print('n_nodes', n_nodes)
print('n_elements', n_elements)


plot_mesh(matriz_coordenadas, matriz_topologia, sigma_real)


#print('I_matrix', I_matrix)



# === Resolução do problema direto ===
nos_contorno = np.arange(n_eletrodos)
V_medido = []

#for k in range(n_eletrodos):
for k in range(n_eletrodos):
    I_vec = np.zeros(n_nodes)
    I_vec = I_matrix[k]
    Y_global = montar_Y_global(matriz_coordenadas, matriz_topologia, n_nodes, sigma_real)
    Y_mod, I_mod = aplica_cond_contorno(I_matrix, Y_global, n_nodes, V_imposto)
    #V_sol = np.linalg.solve(Y_mod, I_mod)
    inv_Y_mod = np.linalg.inv(Y_mod)
    V_sol = np.dot(inv_Y_mod, I_mod)
    V_medido.append(V_sol)

V_medido = np.array(V_medido[0])
V_medido = np.array(V_medido[0:n_eletrodos, 0:n_eletrodos])

#V_medido.flatten().tofile("Novo_Vmedido_gerado_por_sigma.txt", sep="\n")
#print("Novo vetor V_medido salvo como 'Novo_Vmedido_gerado_por_sigma.txt'")

#print('V_medido', V_medido)
print('V_medido', V_medido.shape)
#print(np.array2string(V_medido, formatter={'float': '{:0.2e}'.format}))

np.set_printoptions(formatter={'float': '{:0.2e}'.format}, linewidth=300)

print('V_medido \n', V_medido)

nome_arquivo = 'Pos_To_Gmsh'

criar_arquivo_pos(matriz_coordenadas, matriz_topologia, n_eletrodos, V_sol, nome_arquivo)

abrir_Gmsh_pos(nome_arquivo, n_eletrodos)
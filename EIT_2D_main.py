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
from EIT_functions_2D import plot_mesh

###############################################################################
########################## CARREGA MALHA 1D ###################################
###############################################################################

#opcao = input('Escolha a malha(71, 100, 200, 300, 1000): ')
opcao = '1'

#fundo = input('Escolha a condutividade da base (Ex: 0.5; 1.0, etc...): ')

#anomalia = input('Escolha a condutividade da anomalia (Ex: 0.5; 1.0, etc...): ')


try:
    anomalia = float(input("Escolha a condutividade da anomalia (Ex: 0.5; 1.0, etc...): "))
    print(f"Você digitou: {anomalia}")
except ValueError:
    print("Valor inválido! Digite um número válido, usando ponto decimal.")
    sys.exit() 
    
try:
    fundo = float(input("Escolha a condutividade da base da 'cuba' (Ex: 0.5; 1.0, etc...): "))
    print(f"Você digitou: {fundo}")
except ValueError:
    print("Valor inválido! Digite um número válido, usando ponto decimal.")
    sys.exit() 
    
#fundo = 0.3
#anomalia = 0.2


if opcao not in ['1']: #,'100', '100', '200', '300', '1000']:
    raise ValueError("Opção inválida.")

caminho = files_mesh(opcao)

malha_msh = meshio.read(caminho)                            # Lê o arquivo .msh

matriz_coordenadas = malha_msh.points[:, :2]       # Monta matriz coordenadas
#matriz_topologia = malha_msh.cells_dict["line"]      # Monta Matriz topologia
matriz_topologia = malha_msh.cells_dict["triangle"]      # Monta Matriz topologia
#matriz_topologia = matriz_topologia +1

n_nodes = matriz_coordenadas.shape[0]                             # Nr de nós
n_elements = matriz_topologia.shape[0]                      # Nr de elementos

print('n_nodes', n_nodes)
print('n_elements', n_elements)
#print('matriz_coordenadas', matriz_coordenadas)
#print('matriz_topologia', matriz_topologia)


# Excolher sigma da base e sigma corpo interno respectivamente 
sigma_real = np.concatenate([np.full(128, fundo), np.full(18, anomalia)])

plot_mesh(matriz_coordenadas, matriz_topologia, sigma_real)


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


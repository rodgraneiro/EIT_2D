###############################################################################
# Esse arquivo tem as malhas de Elementos Finitos necessáras para rodar o 
# progrma EIT_2D_main.py
###############################################################################


import numpy as np

def files_mesh(opcao):
    caminho_malhas = {
    '1' : r'octogono_2208.msh',     # malha elem de tamanhos diferentes    
    # '100' : r'd:/ParaGitHub/Estudo-EIT/unidimensional_100e_py.msh',
    }
    return caminho_malhas[opcao]

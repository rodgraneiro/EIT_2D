###############################################################################
# Esse arquivo tem as malhas de Elementos Finitos necessáras para rodar o 
# progrma EIT_2D_main.py
###############################################################################


import numpy as np

def files_mesh(opcao):
    caminho_malhas = {
    '1' : r'octogono_2208_mod.msh',     # malha elem de tamanhos diferentes
    '2' : r'DoisTriangulos_4elem.msh', # DoisTriangulos  4 elementos    
    # '100' : r'd:/ParaGitHub/Estudo-EIT/unidimensional_100e_py.msh',
    }
    sigma_reais = {
        '1': np.concatenate([np.full(128, 0.1), np.full(18, 0.75)]),
        '2': np.concatenate([np.full(1, 0.1), np.full(2, 1.0), np.full(1, 0.1)]),
    }
    V_imposto = {
        '1':  [[13, 0.0]]  ,
        '2': [[4, 0.0]],
    }
    I_matrix = {
        '1': np.loadtxt("PadraoCC_pula_1a.txt")  ,
        '2': np.loadtxt("PadraoCC_pula_1_4elem.txt"),
    }
    n_eletrodos = {
        '1': 8  ,
        '2': 4,
    }
    return caminho_malhas[opcao], sigma_reais[opcao], V_imposto[opcao], I_matrix[opcao], n_eletrodos[opcao]


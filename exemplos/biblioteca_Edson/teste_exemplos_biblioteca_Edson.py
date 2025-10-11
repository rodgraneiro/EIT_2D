import numpy as np
import mesh
import forwardProblem
import matplotlib.pyplot as plt


'''
# TESTE TAREFA UNIDIMENSIONAL
###############################################################################
nome = 'D:/Git_EIT_2D/EIT_2D/malhasMSH/unidimensional_100e_py.msh'

noh_eletrodos = [0,  10, 20, 30, 40, 50, 60, 70, 80, 90,  100]

MinhaMalha1D = mesh.PointElectrodes1DMeshEdson(noh_eletrodos, nome_msh=nome, altura1D = 0.02)
MinhaMalha1D.ReadMesh()


meus_sigmas = {}
meus_sigmas[1] = 0.25
meus_sigmas[2] = 0.5
MinhaMalha1D.SetSigmaPhysicaEntity(meus_sigmas) # Informando sigma (e já calculando o rho de cada elemento)

amplitude_corrente = 0.001
Pcorrente = np.zeros(MinhaMalha1D.NumberOfNodes)                    # Monta vetor de corrente
Pcorrente[0] = -amplitude_corrente                            # Nó de saída de corrente
Pcorrente[MinhaMalha1D.NumberOfElements] = amplitude_corrente                  # Nó de entrada de corrente    

fwd = forwardProblem.forward_problem(MinhaMalha1D, Pcorrente=Pcorrente)   # __init__ roda aqui

fwd.Solve()

plt.plot(fwd.Vmedido)

nome_arquivo = 'ParaVernoGmsh1D'

N_padraoCC = 1

fwd.criar_arquivo_pos_2D(N_padraoCC, fwd.Vmedido, nome_arquivo)

fwd.abrir_Gmsh_pos(nome_arquivo, N_padraoCC)
###############################################################################
'''




# TESTE COM HUA

###############################################################################
#nome = 'D:/Git_EIT_2D/EIT_2D/malhasMSH/circ2_tst_Hua_v2_2_lc_especial.msh'
nome = '../../malhasMSH/circ2_tst_Hua_v2_2_lc_especial.msh'
#nome = '../malhasMSH/circ8_anom4_tst_Hua_v2_2_lc_especial.msh'

MinhaMalha = mesh.HuaElectrodes2DMeshEdson(2, nome_msh=nome, altura2D = 0.02)
MinhaMalha.ReadMesh() 

print(MinhaMalha.Elements[2])
print(f"Centroid: {MinhaMalha.Elements[2].Centroid}")
print(f"KGeo: \n{MinhaMalha.Elements[2].KGeo}")

meus_sigmas = {}
meus_sigmas[1000] = 0.1
meus_sigmas[5001] = 0.2
meus_sigmas[5002] = 0.2
MinhaMalha.SetSigmaPhysicaEntity(meus_sigmas) # Informando sigma (e já calculando o rho de cada elemento)

#MinhaMalha.CalcKGlobal() # calculando KGlobal usando Sigmas

#coordenadas = MinhaMalha.Coordinates
#topologia = MinhaMalha.msh_topology

#MinhaMalha.KGlobal

#KGlobal =  MinhaMalha.KGlobal


#print(f'n_nodes = {MinhaMalha.NumberOfNodes}')

Pcorrente = np.zeros(MinhaMalha.NumberOfNodes)                    # Monta vetor de corrente
Pcorrente[17] = -0.001                            # Nó de saída de corrente
Pcorrente[18] = 0.001                             # Nó de entrada de corrente

V_imposto = [[11, 0.0]]  

noh_eletrodos = [17,  18]

fwd = forwardProblem.forward_problem(MinhaMalha, Pcorrente=Pcorrente)   # __init__ roda aqui

fwd.Solve()

#fwd = forwardProblem.forward_problem(MinhaMalha, V_imposto, noh_eletrodos, Pcorrente)   # __init__ roda aqui
####

print(f'Vmedido \n {fwd.Vmedido}')


nome_arquivo = 'ParaVernoGmshzz2'

N_padraoCC = 1

fwd.criar_arquivo_pos_2D(N_padraoCC, fwd.Vmedido, nome_arquivo)

fwd.abrir_Gmsh_pos(nome_arquivo, N_padraoCC)

###############################################################################

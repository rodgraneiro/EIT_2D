# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:11:19 2026

@author: rodgr
"""

import numpy as np
import mesh
import forwardProblem
import inverseProblem
import inverseProblem_2D
import matplotlib.pyplot as plt

#nome = '../../malhasMSH/circ2_tst_Hua_v2_2_lc_especial.msh'
#nome = '../../malhasMSH/circ8_anom4_tst_Hua_v4_1_lc_0_01.msh'
#nome = '../../malhasMSH/circ8_anom4_tst_Hua_v2_2_lc_especial.msh'
#nome = '../../malhasMSH/circ8_anomZero_Hua.msh'
#nome = '../../malhasMSH/circ16_anomZero_Hua.msh'
#nome = '../../malhasMSH/circ16_anom1_Hua.msh'

#nome = '../../malhasMSH/Hua_quadrado_2e.msh'
#nome = '../../malhasMSH/circ4_anomZero_Hua.msh'
#nome = '../../malhasMSH/Hua_quadrado_4eletrodos_new.msh'
#nome = '../../malhasMSH/circ16_anom1_Square_Hua.msh'

nome = '../../malhasMSH/circ4_objetoUm_Hua.msh'


MinhaMalha = mesh.HuaElectrodes2DAnisotropic(4, nome_msh=nome, altura2D = 0.02, thetaAngle = 0.0, sigmaX = 1.00, sigmaY = 0.1000)

#MinhaMalha = mesh.HuaElectrodes2DAnisotropic(8, nome_msh=nome, altura2D = 0.02, thetaAngle = -45.0, sigmaX = 1000.00, sigmaY = 1.0)

#MinhaMalha = mesh.HuaElectrodes2DAnisotropic(2, nome_msh=nome, altura2D = 0.02)
#MinhaMalha = mesh.HuaElectrodes2DMeshEdson(16, nome_msh=nome, altura2D = 0.02)
MinhaMalha.ReadMesh() 

print('MinhaMalha.Elements[2]',MinhaMalha.Elements[2])
print(f"Centroid: {MinhaMalha.Elements[2].Centroid}")
#print(f"KGeo: \n{MinhaMalha.Elements[2].KGeo}")


meus_sigmas = {
    1000: [1.0, 0.0, 1.0],
    1001: [0.010, 0.0, 0.001],
    5001: [1000000.0, 0.0, 1000000.0],
    5002: [1000000.0, 0.0, 1000000.0],
    5003: [1000000.0, 0.0, 1000000.0],
    5004: [1000000.0, 0.0, 1000000.0]
}

MinhaMalha.SetSigmaAnisotropicElements(meus_sigmas)

for idx in range(MinhaMalha.NumberOfElements):
    if not MinhaMalha.Elements[idx].FlagIsElectrode:
        MinhaMalha.Elements[idx].CalcKgeo()


#meus_sigmas= np.loadtxt("circ16_anom1_Hua_Sigma_v6.txt")
#meus_sigmas= np.loadtxt("circ8_anom4_tst_Hua_Sigma.txt")
#meus_sigmas= np.loadtxt("circ16_anom1_Square_Hua_Sigma.txt")


#MinhaMalha.SetSigmaPhysicaEntity(meus_sigmas) # Informando sigma (e já calculando o rho de cada elemento)

#MinhaMalha.SetSigmaElements(meus_sigmas)


#MinhaMalha.SetSigmaAnisotropicElements(meus_sigmas)

MinhaMalha.CalcKGlobal() # calculando KGlobal usando Sigmas

#coordenadas = MinhaMalha.Coordinates
#topologia = MinhaMalha.msh_topology

#MinhaMalha.KGlobal

#KGlobal =  MinhaMalha.KGlobal


#print(f'n_nodes = {MinhaMalha.NumberOfNodes}')


fwd = forwardProblem.forward_problem(MinhaMalha, Pcorrente=None, SkipPattern=1, VirtualNode = True, I =1.0e-3)   # __init__ roda aqui

#fwd = forwardProblem.forward_problem(MinhaMalha, Pcorrente=None, SkipPattern=3, VirtualNode = True)   # __init__ roda aqui


#print(f'Pcorrente \n {fwd.corrente[MinhaMalha.NumberOfNodes-MinhaMalha.NumberOfElectrodes: MinhaMalha.NumberOfNodes]}')

#print(f'Pcorrente \n {fwd.corrente.shape}')

mtz_Vmedido = fwd.Solve()
#print(f'Vmedido \n {fwd.Vmedido[:,0]}')

nome_arquivo = 'ParaVernoGmshPto'
fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)

fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)

V_measured = fwd.Vmedido_eletrodos

print(f'V_mesured \n {V_measured}')


#invProblem_2D = inverseProblem_2D.inverse_problem(MinhaMalha, Pcorrente=fwd.corrente)
#invProblem_2D.solve(V_measured, initialEstimate=2.9,alpha =2.5,  Lambda = 0.50, max_iter=1,Tol=5.0e-4)
#print('Y_jacobian',invProblem.Y_jacobian)

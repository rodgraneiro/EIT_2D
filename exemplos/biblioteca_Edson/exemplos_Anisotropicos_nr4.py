# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 21:59:28 2026

@author: rodgr
"""

import numpy as np
import mesh
import forwardProblem
import inverseProblem
import inverseProblem_2D
import matplotlib.pyplot as plt


#######################################################
#nome = '../../malhasMSH/circ16_base_coarse.msh'
#nome = '../../malhasMSH/quatro_base_22jan25.msh'
#nome = '../../malhasMSH/circ16_base_v2208.msh'
#nome = '../../malhasMSH/circ16_3_anomalia6_coarse_v2208.msh'

#nome = '../../malhasMSH/circ16_anomalia6B.msh'
nome = '../../malhasMSH/fantoma_base.msh'

MinhaMalhaPto2 = mesh.PointElectrodes2DMeshAnisotropic(4, nome_msh=nome, altura2D=0.002, thetaAngle= 0.0 )

MinhaMalhaPto2.ReadMesh()

meus_sigmas = {
    1000: [1.0, 0.0, 1.0],
    #1001: [1.0, 0.0, 0.00001]
}

MinhaMalhaPto2.SetSigmaPhysicaEntity(meus_sigmas)

for idx in range(MinhaMalhaPto2.NumberOfElements):
    elem = MinhaMalhaPto2.Elements[idx]
    if not elem.FlagIsElectrode:
        elem.CalcKgeo()

fwd = forwardProblem.forward_problem(MinhaMalhaPto2, Pcorrente=None, SkipPattern=1, I =1.0e-3)   # __init__ roda aqui

mtz_Vmedido = fwd.Solve()
nome_arquivo = 'ParaVernoGmshPto'
fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)

fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)

V_measured = fwd.Vmedido_eletrodos


print(f'V_mesured \n {V_measured}')






'''
MinhaMalhaBase16 = mesh.PointElectrodes2DMeshEdson(16, nome_msh=nome, altura2D = 0.02)
#MinhaMalhaPto2 = mesh.PointElectrodes2DMeshEdson(4, nome_msh=nome, altura2D = 0.02)
MinhaMalhaBase16.ReadMesh() 

#meus_sigmas = {1000 : 1.0}
meus_sigmas = {
1000 : 1.0,   
1001 : 1.0}#,
#1002 : 1.0,
#1003 : 1.0}

MinhaMalhaBase16.SetSigmaPhysicaEntity(meus_sigmas)


MinhaMalhaBase16.CalcKGlobal() # calculando KGlobal usando Sigmas


#print(f'MinhaMalhaPto2.KGlobal =  {MinhaMalhaBase16.KGlobal.shape}')
#PcorrenteBase = np.loadtxt("padrao128CC_392.txt")
fwd = forwardProblem.forward_problem(MinhaMalhaBase16, Pcorrente=None, SkipPattern=3, I =1.0e-3)   # __init__ roda aqui


#print(MinhaMalhaBase16.Elements[2])
#print(f"Centroid: {MinhaMalhaBase16.Elements[2].Centroid}")
#print(f"KGeo: \n{MinhaMalhaBase16.Elements[2].KGeo}")


#iteration= np.loadtxt("lastIteration.txt")
iteration=0

#sigma_inicial_rnd = np.random.uniform(2.6, 2.9, MinhaMalhaBase16.NumberOfElements)
#print('x0', sigma_inicial_rnd)


#sigma_inicial_cont = np.loadtxt("sigma_inicial_cont.txt")
invProblem_2D = inverseProblem_2D.inverse_problem(MinhaMalhaBase16, Pcorrente=fwd.corrente)
#invProblem_2D = inverseProblem_2D.inverse_problem(MinhaMalhaBase16, Pcorrente=PcorrenteBase)


#invProblem_2D.solve(V_measured,initialEstimate=sigma_inicial_rnd, alpha =0.1,  Lambda = 0.50, max_iter=3, Tol=1.0e-20, iteration=iteration)
invProblem_2D.solve(V_measured, initialEstimate=3.5, alpha =0.2500,  Lambda = 1.0e-3, max_iter=1000, Tol=1.0e-2, iteration=iteration)
#print('Y_jacobian',invProblem.Y_jacobian)
'''
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:11:19 2026

@author: rodgr
"""

import numpy as np
import mesh
import forwardProblem
import inverseProblem_2D_Anisotropic_Hua
import matplotlib.pyplot as plt

def runFWD_InverseProblemAnisotropicHua():

    #nome = '../../malhasMSH/circ4_objetoUm_Hua_coarse.msh'
    #nome = '../../malhasMSH/Hua_cuba16eletrodos_3objetos.msh'


    #nome = '../../malhasMSH/Hua_cuba4eletrodos_1objetoDireita.msh'
    nome = '../../malhasMSH/Hua_cuba16eletrodos_1objeto_denso.msh'

    MinhaMalha = mesh.HuaElectrodes2DAnisotropic(16, nome_msh=nome, altura2D = 0.02, thetaAngle = -30.0)#, sigmaX = 1.00, sigmaY = 1.0000)
    #MinhaMalha = mesh.HuaElectrodes2DAnisotropic(8, nome_msh=nome, altura2D = 0.02, thetaAngle = -45.0, sigmaX = 1000.00, sigmaY = 1.0)

    MinhaMalha.ReadMesh() 

    print('MinhaMalha.Elements[2]',MinhaMalha.Elements[2])
    print(f"Centroid: {MinhaMalha.Elements[2].Centroid}")
    #print(f"KGeo: \n{MinhaMalha.Elements[2].KGeo}")


    meus_sigmas = {
        1000: [3.0, 0.0, 1.0],
        1001: [3.0, 0.0, 1.0],
        1002: [1.0, 0.0, 1.0],
        1003: [1.0, 0.0, 1.0],
        5001: [1.0, 0.0, 1.0],
        5002: [1.0, 0.0, 1.0],
        5003: [1.0, 0.0, 1.0],
        5004: [1.0, 0.0, 1.0],
        5005: [1.0, 0.0, 1.0],
        5006: [1.0, 0.0, 1.0],
        5007: [1.0, 0.0, 1.0],
        5008: [1.0, 0.0, 1.0],
        5009: [1.0, 0.0, 1.0],
        5010: [1.0, 0.0, 1.0],
        5011: [1.0, 0.0, 1.0],
        5012: [1.0, 0.0, 1.0],
        5013: [1.0, 0.0, 1.0],
        5014: [1.0, 0.0, 1.0],
        5015: [1.0, 0.0, 1.0],
        5016: [1.0, 0.0, 1.0]
    }

    MinhaMalha.SetSigmaAnisotropicElementsHua(meus_sigmas)


    
    for idx in range(MinhaMalha.NumberOfElements):
        if not MinhaMalha.Elements[idx].FlagIsElectrode:
            MinhaMalha.Elements[idx].CalcKgeo()
    '''
    for idx in range(MinhaMalha.NumberOfElements):
        if MinhaMalha.Elements[idx].FlagIsElectrode:
            MinhaMalha.Elements[idx].CalcKgeo()
    '''



    MinhaMalha.CalcKGlobal() # calculando KGlobal usando Sigmas

    fwd = forwardProblem.forward_problem(MinhaMalha, Pcorrente=None, SkipPattern=3, VirtualNode = True, I =1.0e-3, name = None, imageSave = True)   # __init__ roda aqui

    mtz_Vmedido = fwd.Solve()
    print(f'Vmedido \n {fwd.Vmedido[:10]}')

    nome_arquivo = 'ParaVernoGmshPto'
    fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)
    fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)

    V_measured_phaton = fwd.Vmedido_eletrodos
    np.save("V_measured_phaton.npy", V_measured_phaton)  # formato binário
    print(f'V_mesured\n {V_measured_phaton}')


#############################################################################################
#############################################################################################
#############################################################################################

def runInverseProblemAnisotropicHua():

    #nome = '../../malhasMSH/circ4_objetoUm_Hua_coarse.msh'
    #nome = '../../malhasMSH/Hua_cuba16eletrodos_3objetos.msh'


    #nome = '../../malhasMSH/Hua_cuba4eletrodos_1objetoDireita.msh'
    nome = '../../malhasMSH/Hua_cuba16eletrodos_base.msh'

    MinhaMalha_base = mesh.HuaElectrodes2DAnisotropic(16, nome_msh=nome, altura2D = 0.02, thetaAngle = 0.0)#, sigmaX = 1.00, sigmaY = 1.0000)
    #MinhaMalha = mesh.HuaElectrodes2DAnisotropic(8, nome_msh=nome, altura2D = 0.02, thetaAngle = -45.0, sigmaX = 1000.00, sigmaY = 1.0)

    MinhaMalha_base.ReadMesh() 

    print('MinhaMalha.Elements[2]',MinhaMalha_base.Elements[2])
    print(f"Centroid: {MinhaMalha_base.Elements[2].Centroid}")
    #print(f"KGeo: \n{MinhaMalha.Elements[2].KGeo}")


    meus_sigmas = {
        1000: [1.0, 0.0, 1.0],
        1001: [1.0, 0.0, 1.0],
        1002: [1.0, 0.0, 1.0],
        1003: [1.0, 0.0, 1.0],
        5001: [1.0, 0.0, 1.0],
        5002: [1.0, 0.0, 1.0],
        5003: [1.0, 0.0, 1.0],
        5004: [1.0, 0.0, 1.0],
        5005: [1.0, 0.0, 1.0],
        5006: [1.0, 0.0, 1.0],
        5007: [1.0, 0.0, 1.0],
        5008: [1.0, 0.0, 1.0],
        5009: [1.0, 0.0, 1.0],
        5010: [1.0, 0.0, 1.0],
        5011: [1.0, 0.0, 1.0],
        5012: [1.0, 0.0, 1.0],
        5013: [1.0, 0.0, 1.0],
        5014: [1.0, 0.0, 1.0],
        5015: [1.0, 0.0, 1.0],
        5016: [1.0, 0.0, 1.0]
    }

    MinhaMalha_base.SetSigmaAnisotropicElementsHua(meus_sigmas)



    for idx in range(MinhaMalha_base.NumberOfElements):
        if not MinhaMalha_base.Elements[idx].FlagIsElectrode:
            MinhaMalha_base.Elements[idx].CalcKgeo()
    '''
    for idx in range(MinhaMalha.NumberOfElements):
        if MinhaMalha.Elements[idx].FlagIsElectrode:
            MinhaMalha.Elements[idx].CalcKgeo()
    '''

    start = [3.0, 0.0, 3.0]

    MinhaMalha_base.CalcKGlobal() # calculando KGlobal usando Sigmas

    fwd = forwardProblem.forward_problem(MinhaMalha_base, Pcorrente=None, SkipPattern=0, VirtualNode = True, I =1.0e-3)   # __init__ roda aqui

    mtz_Vmedido = fwd.Solve()
    #print(f'Vmedido \n {fwd.Vmedido[:10]}')

    nome_arquivo = 'ParaVernoGmshPto'
    #fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)
    #fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)
    V_measured_phaton = fwd.Vmedido_eletrodos
    print(f'V_measured_phaton\n {V_measured_phaton.shape}')
    

    invProblem_2D = inverseProblem_2D_Anisotropic_Hua.inverse_problem(MinhaMalha_base, Pcorrente=fwd.corrente)
    invProblem_2D.solve(V_measured_phaton, initialEstimate=start,alpha =0.01,  Lambda = 0.001, max_iter=10,Tol=1.0e-6)
    #print('Y_jacobian',invProblem.Y_jacobian)




runFWD_InverseProblemAnisotropicHua()

#runInverseProblemAnisotropicHua()




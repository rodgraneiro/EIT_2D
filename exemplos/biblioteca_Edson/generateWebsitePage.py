# %%time
###############################################################################
###############################################################################
####################### Problema Direto ######################################
###############################################################################

import numpy as np
import mesh
import forwardProblem
import inverseProblem_2D_Anisotropic_Hua
import matplotlib.pyplot as plt
import os



#############################################################################################
#############################################################################################
#############################################################################################


def runFWD_InverseProblemAnisotropicHua(NrElec, MSH_name, name = "batataBanana", hight2D=0.02, thetaAngle = 0.0, imageSave = True):

    bName = name

    MinhaMalha = mesh.HuaElectrodes2DAnisotropic(NrElec, MSH_name, hight2D, thetaAngle)#, sigmaX = 1.00, sigmaY = 1.0000)

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

    fwd = forwardProblem.forward_problem(MinhaMalha, Pcorrente=None, SkipPattern=3, VirtualNode = True, I =1.0e-3, name = bName, imageSave = True)   # __init__ roda aqui

    mtz_Vmedido = fwd.Solve()
    print(f'Vmedido \n {fwd.Vmedido[:10]}')

    nome_arquivo = 'ParaVernoGmshPto'
    fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)
    fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)

    V_measured_phaton = fwd.Vmedido_eletrodos
    np.save("V_measured_phaton.npy", V_measured_phaton)  # formato binário
    print(f'V_mesured\n {V_measured_phaton}')
    print(f'meus_sigmas\n {meus_sigmas}')

    
#############################################################################################
#############################################################################################
#############################################################################################

NrElec =16
mshName = '../../malhasMSH/test_Olavo_Hua.msh'
    
runFWD_InverseProblemAnisotropicHua(NrElec,mshName, thetaAngle = 0.0, imageSave = True)
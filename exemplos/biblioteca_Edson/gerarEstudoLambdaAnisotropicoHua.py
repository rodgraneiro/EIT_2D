# %%time
###############################################################################



###############################################################################
####################### Problema Direto ######################################
###############################################################################

import numpy as np
import mesh
import forwardProblem
#import inverseProblem
#import inverseProblem_2D
#import inverseProblem_2D_Hua
import inverseProblem_2D_Anisotropic_Hua
import matplotlib.pyplot as plt
import os

#############################################################################################
#############################################################################################
#############################################################################################


def rodar_simulacao(lambda_val, html_name="resultado"):
    #nome = '../../malhasMSH/Hua_cuba16eletrodos_3objetos.msh'
    #nome = '../../malhasMSH/circ4_objetoUm_Hua_coarse.msh'
    #nome = '../../malhasMSH/Hua_cuba16eletrodos_3objetos.msh'


    #nome = '../../malhasMSH/Hua_cuba4eletrodos_1objetoDireita.msh'
    #nome = '../../malhasMSH/Hua_cuba16eletrodos_base.msh'
    #nome = '../../malhasMSH/test_Olavo_baseZeroGrau.msh'
    #nome = '../../malhasMSH/test_Olavo_30grausNeg.msh'
    #nome = '../../malhasMSH/test_Olavo_60graus.msh'
    nome = '../../malhasMSH/Hua_domain16e_base_v2D.msh'
    #nome = '../../malhasMSH/Domain32_Hua_Base.msh'

    


    MinhaMalha_base = mesh.HuaElectrodes2DAnisotropic(16, nome_msh=nome, altura2D = 0.02)#, thetaAngle = 0.0)#, sigmaX = 1.00, sigmaY = 1.0000)
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
        5016: [1.0, 0.0, 1.0],
        5017: [1.0, 0.0, 1.0],
        5018: [1.0, 0.0, 1.0],
        5019: [1.0, 0.0, 1.0],
        5020: [1.0, 0.0, 1.0],
        5021: [1.0, 0.0, 1.0],
        5022: [1.0, 0.0, 1.0],
        5023: [1.0, 0.0, 1.0],
        5024: [1.0, 0.0, 1.0],
        5025: [1.0, 0.0, 1.0],
        5026: [1.0, 0.0, 1.0],
        5027: [1.0, 0.0, 1.0],
        5028: [1.0, 0.0, 1.0],
        5029: [1.0, 0.0, 1.0],
        5030: [1.0, 0.0, 1.0],
        5031: [1.0, 0.0, 1.0],
        5032: [1.0, 0.0, 1.0]        
        
        
    }

    MinhaMalha_base.SetSigmaAnisotropicElementsHua(meus_sigmas)



    for idx in range(MinhaMalha_base.NumberOfElements):
        if not MinhaMalha_base.Elements[idx].FlagIsElectrode:
            MinhaMalha_base.Elements[idx].CalcKgeo()


    start = [5.0, 0.0, 5.0]
    
    #start = np.loadtxt('sigmaInicial_result.txt')

    MinhaMalha_base.CalcKGlobal() # calculando KGlobal usando Sigmas

    fwd = forwardProblem.forward_problem(MinhaMalha_base, Pcorrente=None, SkipPattern=3, VirtualNode = True, I =1.0e-3)   # __init__ roda aqui

    mtz_Vmedido = fwd.Solve()
    #print(f'Vmedido \n {fwd.Vmedido[:10]}')

    nome_arquivo = 'ParaVernoGmshPto'
    #fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)
    #fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)
    #V_measured_phaton = fwd.Vmedido_eletrodos
    V_measured_phaton = np.load("V_measured_phaton.npy")
    print(f'V_measured_phaton\n {V_measured_phaton.shape}')
    #htmlName = 'XXXrectangularHomogeneousAnisotropy30Neg'
    htmlName = nome_html
    invProblem_2D = inverseProblem_2D_Anisotropic_Hua.inverse_problem(MinhaMalha_base, Pcorrente=fwd.corrente)
    invProblem_2D.solve(V_measured_phaton, initialEstimate=start,alpha =0.075,  Lambda = lambda_val, max_iter= 250,Tol=1.0e-9, html_name = htmlName)
    #print('Y_jacobian',invProblem.Y_jacobian)

#sigma_inicial_cont = np.loadtxt("sigma_inicial_cont.txt")
#print(sigma_inicial_cont)


#lambdas = np.logspace(-5, -4, 10)

lambdas =[7.02e-06]
resultados = {}

nome_html="CorrecaoArtigo"
#pasta="../../docs/figureTemp"
#pasta2="../../docs"

for lam in lambdas:
    print(f"\nRodando lambda = {lam:.5f}")
    
    resultados[lam] = rodar_simulacao(lam, html_name=nome_html)

#inverseProblem_2D_Anisotropic_Hua.inverse_problem.salvar_html_todos_lambdas(pasta2, nome_html)

#rodar_simulacao(6.57933225e-03, None)


# %%time
###############################################################################
###############################################################################
# Arquivo de teste para desenvolvimento do problema inverso 
# para malha isotrópica 
###############################################################################
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



def rodar_simulacao(lambda_val, sigma_saved):
    
    #nome = '../../malhasMSH/Hua_cuba16eletrodos_3objetos.msh'
    #nome = '../../malhasMSH/circ4_objetoUm_Hua_coarse.msh'
    #nome = '../../malhasMSH/Hua_cuba16eletrodos_3objetos.msh'


    #nome = '../../malhasMSH/Hua_cuba4eletrodos_1objetoDireita.msh'
    #nome = '../../malhasMSH/Hua_cuba16eletrodos_base.msh'
    #nome = '../../malhasMSH/test_Olavo_baseZeroGrau.msh'
    #nome = '../../malhasMSH/test_Olavo_30grausNeg.msh'
    nome = '../../malhasMSH/test_Olavo_60grausNeg.msh'

    


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


    start = [5.0, 0.0, 5.0]

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
    
    htmlName = 'test_Olavo_Hua_Ani301_Scale_zc_1_ang_60Neg'
    invProblem_2D = inverseProblem_2D_Anisotropic_Hua.inverse_problem(MinhaMalha_base, Pcorrente=fwd.corrente)
    invProblem_2D.solve(V_measured_phaton, initialEstimate=start,alpha =0.1,  Lambda = lambda_val, max_iter=5,Tol=1.0e-6, html_name = htmlName)
    #print('Y_jacobian',invProblem.Y_jacobian)

#sigma_inicial_cont = np.loadtxt("sigma_inicial_cont.txt")
#print(sigma_inicial_cont)
#lambdas = np.logspace(-4, 2, 10)
#lambdas [1.00000000e-04 4.64158883e-04 2.15443469e-03 1.00000000e-02
#4.64158883e-02 2.15443469e-01 1.00000000e+00 4.64158883e+00
# 2.15443469e+01 1.00000000e+02]


#lambdas = np.logspace(-4, 1, 10)
#lambdas [1.00000000e-04 3.59381366e-04 1.29154967e-03 4.64158883e-03
# 1.66810054e-02 5.99484250e-02 2.15443469e-01 7.74263683e-01
# 2.78255940e+00 1.00000000e+01]



#lambdas = np.logspace(-5, 0, 10)
#e-5
#lambdas [1.00000000e-05 3.59381366e-05 1.29154967e-04 4.64158883e-04
# 1.66810054e-03 5.99484250e-03 2.15443469e-02 7.74263683e-02
# 2.78255940e-01 1.00000000e+00]

#lambdas = np.logspace(-6, 0, 10)
#e-6
#lambdas [1.00000000e-06 4.64158883e-06 2.15443469e-05 1.00000000e-04
# 4.64158883e-04 2.15443469e-03 1.00000000e-02 4.64158883e-02
# 2.15443469e-01 1.00000000e+00]

#lambdas = np.logspace(-7, 0, 10)
#e-7
#lambdas [1.00000000e-07 5.99484250e-07 3.59381366e-06 2.15443469e-05
# 1.29154967e-04 7.74263683e-04 4.64158883e-03 2.78255940e-02
# 1.66810054e-01 1.00000000e+00]

lambdas = np.logspace(-9, 1, 20)
#lambdas = [1.00000000e-09, 3.35981829e-09]#, 1.12883789e-08, 3.79269019e-08,
# 1.27427499e-07, 4.28133240e-07, 1.43844989e-06, 4.83293024e-06,
# 1.62377674e-05, 5.45559478e-05]# 
#lambdas = [1.83298071e-04, 6.15848211e-04,
# 2.06913808e-03, 6.95192796e-03, 2.33572147e-02, 7.84759970e-02,
# 2.63665090e-01, 8.85866790e-01, 2.97635144e+00, 1.00000000e+01]

#lambdas = np.logspace(-9, 1, 10)
#lambdas [1.00000000e-09 1.29154967e-08 1.66810054e-07 2.15443469e-06
# 2.78255940e-05 3.59381366e-04 4.64158883e-03 5.99484250e-02
# 7.74263683e-01 1.00000000e+01]


#lambdas = np.logspace(-9, 0, 10)
#lambdas [1.e-09 1.e-08 1.e-07 1.e-06 1.e-05 1.e-04 1.e-03 1.e-02 1.e-01 1.e+00]



resultados = {}

for lam in lambdas:
    print(f"\nRodando lambda = {lam:.5f}")
    
    resultados[lam] = rodar_simulacao(lam, None)

#rodar_simulacao(1.00000000e-02, None)


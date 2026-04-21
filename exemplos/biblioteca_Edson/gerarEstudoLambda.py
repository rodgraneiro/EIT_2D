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
import inverseProblem_2D_Hua
import matplotlib.pyplot as plt



def rodar_simulacao(lambda_val, sigma_saved):
    
    # ======= SEU CÓDIGO =======
    # leitura da malha
    # montagem KGlobal
    # forward
    # inverse
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################
    nome = '../../malhasMSH/Hua_cuba16eletrodos_base.msh'
    nome_malha = 'Hua_cuba16eletrodos_3sigmas_video'

    V_measured_phaton = np.load('V_measured_phaton.npy', allow_pickle=True)
    print('dados:\n', V_measured_phaton)

    MinhaMalha_base = mesh.HuaElectrodes2DMeshEdson(16, nome_msh=nome, altura2D = 0.02)
    MinhaMalha_base.ReadMesh() 


    meus_sigmas = {
    1000 : 1.0,    
    1001 : 1.0,
    1002 : 1.0,
    1003 : 1.0,
    5001 : 1.0, 
    5002 : 1.0, 
    5003 : 1.0, 
    5004 : 1.0, 
    5005 : 1.0, 
    5006 : 1.0, 
    5007 : 1.0, 
    5008 : 1.0,
    5009 : 1.0, 
    5010 : 1.0, 
    5011 : 1.0, 
    5012 : 1.0, 
    5013 : 1.0, 
    5014 : 1.0, 
    5015 : 1.0, 
    5016 : 1.0
    }

    MinhaMalha_base.SetSigmaPhysicaEntity(meus_sigmas) # Informando sigma (e já calculando o rho de cada elemento)

    fwd = forwardProblem.forward_problem(MinhaMalha_base, Pcorrente=None, SkipPattern=3, VirtualNode = True, name = nome_malha, imageSave = False )   # __init__ roda aqui


    mtz_Vmedido_base = fwd.Solve()


    nome_arquivo = 'ParaVernoGmshPto'
    #fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)
    #fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)
    #nome_malha = 'Hua_cuba16eletrodos_base'

    invProblem_2D = inverseProblem_2D_Hua.inverse_problem(MinhaMalha_base, Pcorrente=fwd.corrente)
    invProblem_2D.solve(V_measured_phaton, initialEstimate=4.5,alpha =1.0,  Lambda = lambda_val, name = nome_malha, pLambda = lambda_val, max_iter=101,Tol=1.0e-6, iteration=0, saveVideo = True)
    #invProblem_2D.solve(V_measured_phaton, initialEstimate=sigma_saved,alpha =1.0,  Lambda = lambda_val, name = nome_malha, pLambda = lambda_val, max_iter=501,Tol=1.0e-9, iteration=201, saveVideo = True)
    return invProblem_2D  # ou o que quiser salvar

sigma_inicial_cont = np.loadtxt("sigma_inicial_cont.txt")
#print(sigma_inicial_cont)
lambdas = np.logspace(-4, 2, 10)
#lambdas [1.00000000e-04 4.64158883e-04 2.15443469e-03 1.00000000e-02
# 4.64158883e-02 2.15443469e-01 1.00000000e+00 4.64158883e+00
# 2.15443469e+01 1.00000000e+02]

resultados = {}
'''
for lam in lambdas:
    print(f"\nRodando lambda = {lam:.5f}")
    
    resultados[lam] = rodar_simulacao(lam, None)
'''
rodar_simulacao(2.15443469e-03, None)


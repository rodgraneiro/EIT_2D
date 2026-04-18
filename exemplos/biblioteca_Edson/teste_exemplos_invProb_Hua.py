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

#nome = '../../malhasMSH/circ2_tst_Hua_v2_2_lc_especial.msh'
#nome = '../../malhasMSH/circ8_anom4_tst_Hua_v4_1_lc_0_01.msh'
#nome = '../../malhasMSH/circ4_objetoUm_Hua.msh'
#nome = '../../malhasMSH/circ4_Um_objetoGrande_Hua.msh'
#nome = '../../malhasMSH/Hua_cruz_editado.msh'
#nome = '../../malhasMSH/circ4_objetoUm_Hua_coarse.msh'
#nome = '../../malhasMSH/Hua_4e_coarse_test.msh'

#nome = '../../malhasMSH/Hua_cuba16eletrodos_3objetos.msh'
nome = '../../malhasMSH/Hua_cuba16eletrodos_1objeto_denso.msh'

#nome = '../../malhasMSH/Hua_cuba16eletrodos_base.msh'




#MinhaMalha = mesh.HuaElectrodes2DMeshEdson(8, nome_msh=nome, altura2D = 0.02)
MinhaMalha = mesh.HuaElectrodes2DMeshEdson(16, nome_msh=nome, altura2D = 0.02)
MinhaMalha.ReadMesh() 

print(MinhaMalha.Elements[2])
print(f"Centroid: {MinhaMalha.Elements[2].Centroid}")
#print(f"KGeo: \n{MinhaMalha.Elements[2].KGeo}")


meus_sigmas = {
1000 : 3.0,    
1001 : 0.0001,
1002 : 2.0,
1003 : 2.0,
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

MinhaMalha.SetSigmaPhysicaEntity(meus_sigmas) # Informando sigma (e já calculando o rho de cada elemento)

#MinhaMalha.CalcKGlobal() # calculando KGlobal usando Sigmas

#coordenadas = MinhaMalha.Coordinates
#topologia = MinhaMalha.msh_topology

#MinhaMalha.KGlobal

#KGlobal =  MinhaMalha.KGlobal


#print(f'n_nodes = {MinhaMalha.NumberOfNodes}')


fwd = forwardProblem.forward_problem(MinhaMalha, Pcorrente=None, SkipPattern=3, VirtualNode = True)   # __init__ roda aqui


#print(f'Pcorrente \n {fwd.corrente[MinhaMalha.NumberOfNodes-MinhaMalha.NumberOfElectrodes: MinhaMalha.NumberOfNodes]}')

#print(f'Pcorrente \n {fwd.corrente.shape}')

mtz_Vmedido = fwd.Solve()
#print(f'Vmedido \n {fwd.Vmedido[:,0]}')

nome_arquivo = 'ParaVernoGmshPto'
#fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)

#fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)

V_measured_phaton = fwd.Vmedido_eletrodos
#np.save("V_measured_phaton.npy", V_measured_phaton)  # formato binário
#print(f'V_mesured\n {V_measured_phaton}')



########################################################################################################
########################################################################################################
########################################################################################################
nome = '../../malhasMSH/Hua_cuba16eletrodos_base.msh'

#dados = np.load('V_measured_phaton.npy', allow_pickle=True)
#print('dados:\n', dados)

#MinhaMalha = mesh.HuaElectrodes2DMeshEdson(8, nome_msh=nome, altura2D = 0.02)
MinhaMalha_base = mesh.HuaElectrodes2DMeshEdson(16, nome_msh=nome, altura2D = 0.02)
MinhaMalha_base.ReadMesh() 

print(MinhaMalha_base.Elements[2])
print(f"Centroid: {MinhaMalha_base.Elements[2].Centroid}")
#print(f"KGeo: \n{MinhaMalha.Elements[2].KGeo}")


meus_sigmas = {
1000 : 1.0,    
#1001 : 2.0,
1002 : 2.0,
1003 : 2.0,
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

fwd = forwardProblem.forward_problem(MinhaMalha_base, Pcorrente=None, SkipPattern=3, VirtualNode = True)   # __init__ roda aqui


mtz_Vmedido_base = fwd.Solve()
#print(f'Vmedido \n {fwd.Vmedido[:,0]}')

nome_arquivo = 'ParaVernoGmshPto'
#fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)

#fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)


invProblem_2D = inverseProblem_2D_Hua.inverse_problem(MinhaMalha_base, Pcorrente=fwd.corrente)
invProblem_2D.solve(V_measured_phaton, initialEstimate=3.5,alpha =1.0,  Lambda = 10.001, max_iter=50,Tol=1.0e-9)



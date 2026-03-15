# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 21:57:38 2026

@author: rodgr
"""


#import numpy as np
#import mesh
import forwardProblem
#import inverseProblem
import inverseProblem_2D
import matplotlib.pyplot as plt
'''
# TESTE COM eletrodo pontual mesh 4 elementos

###############################################################################
#nome = '../../malhasMSH/quatro_triangulos_03nov2025.msh'
#nome = '../../malhasMSH/dezesseis_triangulos_22jan25.msh'
#nome = '../../malhasMSH/circ16_base_plus.msh'
#nome = '../../malhasMSH/fantomaAnisotropia01.msh'
#nome = '../../malhasMSH/quadrado_16e.msh'
#nome = '../../malhasMSH/dezesseis_triangulos_1body.msh'
#nome = '../../malhasMSH/quadrado_16e_plus.msh'
#nome = '../../malhasMSH/quatro_4e_em_XY.msh'

#MinhaMalhaPto2 = mesh.PointElectrodes2DMeshEdson(4, nome_msh=nome, altura2D = 1.00)
MinhaMalhaPto2 = mesh.PointElectrodes2DMeshAnisotropic(4, nome_msh=nome, altura2D = 0.02, thetaAngle = 0.0, sigmaX = 1.0, sigmaY = 1.00)
#MinhaMalhaPto2 = mesh.PointElectrodes2DMeshEdson(4, nome_msh=nome, altura2D = 0.02, useEdson=True)

#LinearTriangleAnisotropic
MinhaMalhaPto2.ReadMesh() 
'''
#print(MinhaMalhaPto2.Elements[2])
#print(f"Centroid: {MinhaMalhaPto2.Elements[2].Centroid}")
#print(f"KGeo: \n{MinhaMalhaPto2.Elements[2].KGeo}")
#sigma_inicial = np.full(MinhaMalhaPto2.NumberOfElements, 1.0)          # Monta vetor sigma inicial


'''
meus_sigmas = {
1000 : 1.000}#,
#1001 : 10.00}#,
#1002 : 2.00,    
#1003 : 2.00}

meus_sigmas = {
1000 : 1.0,
1001 : 0.10,
10000 : 1.0e0,
10001 : 1.0e0,
10002 : 1.0e0,
10003 : 1.0e0, 
10004 : 1.0e0,     
}
'''

#meus_sigmas= np.loadtxt("fantomaAnisotropia01_SigmaAni.txt")

'''
#MinhaMalhaPto2.SetSigmaPhysicaEntity(meus_sigmas)
MinhaMalhaPto2.SetSigmaElements(meus_sigmas)

MinhaMalhaPto2.CalcKGlobal() # calculando KGlobal usando Sigmas

print('meus_sigmas',meus_sigmas)
#print(f'MinhaMalhaPto2.KGlobal =  {MinhaMalhaPto2.KGlobal}')

fwd = forwardProblem.forward_problem(MinhaMalhaPto2, Pcorrente=None, SkipPattern=1, I =1.0e-3)   # __init__ roda aqui

#print(f'Pcorrente \n {fwd.corrente[MinhaMalhaPto.NumberOfNodes-MinhaMalhaPto.NumberOfElectrodes: MinhaMalhaPto.NumberOfNodes]}')
#print(f'Pcorrente \n {fwd.corrente[:16]}')

#print(f'Pcorrente \n {fwd.corrente.shape}')

mtz_Vmedido = fwd.Solve()
#print(f'Vmedido \n {fwd.Vmedido[:,0]}')

nome_arquivo = 'ParaVernoGmshPto'
fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)

fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)

#print(f'self.Yinversa \n {fwd.Yinversa}')

V_measured = fwd.Vmedido_eletrodos
#print('V_measured ',V_measured)
print('meus_sigmas ',meus_sigmas)
'''
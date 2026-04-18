# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 21:54:37 2026

@author: rodgr
"""

import numpy as np
import mesh
import forwardProblem
import inverseProblem
import inverseProblem_2D
import matplotlib.pyplot as plt

# TESTE COM eletrodo pontual mesh 4 elementos

###############################################################################
nome = '../../malhasMSH/quatro_triangulos_03nov2025.msh'
#nome = '../../malhasMSH/dezesseis_triangulos_22jan25.msh'
#nome = '../../malhasMSH/circ16_base_plus.msh'
#nome = '../../malhasMSH/fantomaAnisotropia01.msh'
#nome = '../../malhasMSH/quadrado_16e.msh'
#nome = '../../malhasMSH/dezesseis_triangulos_1body.msh'
#nome = '../../malhasMSH/quadrado_16e_plus.msh'
#nome = '../../malhasMSH/quatro_4e_em_XY.msh'


MinhaMalhaPto2 = mesh.PointElectrodes2DMeshAnisotropic(4, nome_msh=nome, altura2D=0.002, thetaAngle=0.0 )

MinhaMalhaPto2.ReadMesh()

meus_sigmas = {
    1000: [1.0, 0.0, 0.10],
    1001: [3.0, 0.0, 10.0]
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
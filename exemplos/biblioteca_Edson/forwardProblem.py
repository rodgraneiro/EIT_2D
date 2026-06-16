# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 14:17:52 2025

@author: rodgr
"""

import numpy as np
import mesh
import elements
import gmsh
import sys
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import datetime
from matplotlib.colors import TwoSlopeNorm
#import subprocess
#import os


import plotly.graph_objects as go
import os
import glob
from collections import defaultdict
#import matplotlib.cm as cm
import matplotlib.colors as colors


from datetime import datetime
from matplotlib import cm
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Ellipse








class forward_problem: 
    def __init__(self, mymesh, V_imposto=None, Pcorrente=None, SkipPattern=None, VirtualNode = False, I =1.0e-3, name = None,  imageSave = False):
        if not hasattr(mymesh, "KGlobal"): # verifica se o objeto mymesh tem um atributo chamado KGlobal.
            raise TypeError("Parâmetro incorreto: mymesh.")

        self.mymesh = mymesh
        self.Vmedido = None
        self.Yinversa = None
        self.name = name
        self.imageSave = imageSave

        if V_imposto is None:
            V_imposto = [[mymesh.GndNode, 0.0]]
        self.V_imposto = V_imposto

        if Pcorrente is None:
            # ToDo: implementar a matriz de medidas com padrão pula informado (SkipPattern)
            self.corrente = np.zeros((self.mymesh.NumberOfNodes, self.mymesh.NumberOfElectrodes))
            
            if VirtualNode == True:
                ajuste = self.mymesh.NumberOfNodes-self.mymesh.NumberOfElectrodes
            else: ajuste = 0
            
            k = np.arange(self.mymesh.NumberOfElectrodes)  
            self.corrente[k+ajuste, k] = I
            self.corrente[((k + SkipPattern+1) % self.mymesh.NumberOfElectrodes)+ajuste, k] = -I
            #raise Exception("forward_problem: Pcorrente padrão ainda não implementada.")
        else:
            self.corrente = Pcorrente


    ###############################################################################
    # Esta função aplica condições de contorno conforme exemplo abaixo:
    #
    # [ 1   0          0          0          0        ] [u1]   [u1]    
    # [ 0  k1+k2+k3    -k3        0         -k2       ] [u2] = [F2 +k1*u1] 
    # [ 0   -k3      k3+k5+k4     -k5        0        ] [u3]   [F3 +k4*u4] 
    # [ 0   0          -k5      k2+k6       -k6       ] [u4]   [F4 +k6*u5] 
    # [ 0   0           0          0         1        ] [u5]   [u5]
    ###############################################################################  
    def apply_boundary_conditions(self):
        self.vetor_corrente_cond_contorno = self.corrente[:].copy()
        #print(f'Vetor de corrente: \n {self.vetor_corrente_cond_contorno}')

        self.KGlobal = self.mymesh.KGlobal.copy()       # Criar matriz solução
        #np.savetxt("FWD_KGlobal.txt", self.KGlobal, fmt="%e")
        #print(f' AppCC self.KGlobal: \n {self.KGlobal}')
        # atualiza vetor de correntes:
        for [noh_cond_contorno,valor_cond_contorno] in self.V_imposto:   # Necessário quando valor imposto é diferente de zero
            for i in range(0, self.mymesh.NumberOfNodes):                    # corrige matriz de corrente
                self.vetor_corrente_cond_contorno[i] = (
                    self.vetor_corrente_cond_contorno[i] - \
                    self.KGlobal[i][noh_cond_contorno]*valor_cond_contorno
                )

        for [noh_cond_contorno_b,valor_cond_contorno] in self.V_imposto:
            self.vetor_corrente_cond_contorno[noh_cond_contorno] = valor_cond_contorno  # coloca valor de contorno conhecido
        
        # atualiza matriz KGlobal:
        for [noh_cond_contorno,valor_cond_contorno] in self.V_imposto:
          for k in range(0, self.mymesh.NumberOfNodes):                # laço para zerar linha e coluna
              self.KGlobal[noh_cond_contorno][k] = 0
              self.KGlobal[k][noh_cond_contorno] = 0
        self.KGlobal[noh_cond_contorno][noh_cond_contorno] = 1
    ###############################################################################
    def plotMSH(self, sigma, iteration = None, save = False, SigmaXXXYYY = None, nome_arquivo= None):

        x, y = self.mymesh.Coordinates[:, 0], self.mymesh.Coordinates[:, 1]
        topo = self.mymesh.msh_topology
    
        # --- Separar elementos 2D (triangulares) e 1D (linhas) ---
        elems_2D = np.array([el for el in topo if len(el) == 3])
        elems_1D = np.array([el for el in topo if len(el) == 2])
    
        #fig, ax = plt.subplots(figsize=(6, 5))
        fig, ax = plt.subplots(figsize=(8,5))
        ax.set_aspect('equal')
    
        # ============================================================
        #  1) PLOTAR ELEMENTOS 2D (TRIANGULARES)
        # ============================================================
        if len(elems_2D) > 0:
            triang = tri.Triangulation(x, y, elems_2D)
            #tpc = ax.tripcolor(triang,facecolors=sigma[:len(elems_2D)],edgecolors='k', cmap='Blues')#,vmin=1.0 )
            ntri = triang.triangles.shape[0]
            fc = sigma.ravel()[:ntri]
            lim = np.max(np.abs(fc))
            if lim == 0:
                lim = 1e-12
            norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
            
            if SigmaXXXYYY == 'xy':
                
                tpc = ax.tripcolor(triang,facecolors = fc,edgecolors='k', cmap='RdBu_r', norm=norm )
            if not SigmaXXXYYY == 'xy':
                #tpc = ax.tripcolor(triang,facecolors = fc,edgecolors='k', cmap='Blues', vmin=-5.0, vmax=5.0 )
                #tpc = ax.tripcolor(triang,facecolors = fc,edgecolors='k', cmap='RdBu_r', vmin=0.0, vmax=4.0 )
                tpc = ax.tripcolor(triang,facecolors = fc,edgecolors='k', cmap='rainbow', vmin=0.0, vmax=4.0)
            #triang = tri.Triangulation(x, y, elems_2D)
            #tpc = ax.tripcolor(triang, facecolors=sigma[:len(elems_2D)], edgecolors='k', cmap='Blues',min=0.0, vmax=5.0 )
            #tpc = ax.tripcolor(triang,facecolors = fc,edgecolors='k', cmap='Blues', vmin=0.0, vmax=5.0 )
            fig.colorbar(tpc, ax=ax, label='σ (Conductivity)')
            if save == True:
                #timestamp = datetime.now().strftime("%Y%m%d_%H%M")
                #ax.set_title(f"Conductivity Real (σ) ", fontsize=12)
                ax.set_title(f"Conductivity (σ{SigmaXXXYYY})", fontsize=11)
            if save == False:
                ax.set_title(f"Conductivity Real (σ) ", fontsize=15)
        # ============================================================
        #  2) PLOTAR ELEMENTOS 1D (SEGMENTOS)
        # ============================================================
        if len(elems_1D) > 0:
            for (n1, n2) in elems_1D:
                x_coords = [x[n1], x[n2]]
                y_coords = [y[n1], y[n2]]
                ax.plot(x_coords, y_coords, color='red', linewidth=2)
    
        # ------------------------------------------------------------
        ax.set_xlabel("[m]", fontsize=12)
        ax.set_ylabel("[m]", fontsize=12)
        plt.tight_layout()
        nome_arquivo = f"../../docs/{self.name}.webp"
        plt.savefig(f'{nome_arquivo}',  dpi=200, pil_kwargs={"quality": 70})
        plt.savefig(nome_arquivo, dpi=300, bbox_inches='tight') 
        plt.show() 
         
        
        plt.show(block=False)

        plt.pause(0.1)  


    ###############################################################################   
    def plotElipse(self, sigma, sigma_L, sigma_T, theta_deg, save = False,  nome_arquivo= None):

        x, y = self.mymesh.Coordinates[:, 0], self.mymesh.Coordinates[:, 1]
        topo = self.mymesh.msh_topology
        ajuste  = self.mymesh.NumberOfElectrodes*4 +1
    
        # --- Separar elementos 2D (triangulares) e 1D (linhas) ---
        elems_2D = np.array([el for el in topo if len(el) == 3])
        elems_1D = np.array([el for el in topo if len(el) == 2])

        # parâmetros
        R_circ = 0.15
        R_temp = 0.145      # raio do círculo
        dx = 0.01    # espaçamento
        
        x_pts = []
        y_pts = []
        
        # número de passos
        N = int(np.floor(R_temp/dx))
        
        # gera pontos internos
        for i in range(-N, N+1):
            for j in range(-N, N+1):
        
                x_temp = i * dx
                y_temp = j * dx
        
                # mantém apenas pontos internos
                if x_temp**2 +  y_temp**2 <= R_temp**2:
                    x_pts.append(x_temp)
                    y_pts.append(y_temp)
        #if len(elems_2D) > 0:
        triang = tri.Triangulation(x, y, elems_2D)
        #tpc = ax.tripcolor(triang,facecolors=sigma[:len(elems_2D)],edgecolors='k', cmap='Blues')#,vmin=1.0 )
        ntri = triang.triangles.shape[0]
        fc = sigma.ravel()[:ntri]
        
        finder = triang.get_trifinder()
        pontos = np.column_stack((x_pts, y_pts))
        print("pontos:", pontos.shape)

        ex = pontos[:,0]
        ey = pontos[:,1]
        #ex, ey = 0.02, 0.07
    
        idx_elem = finder(ex, ey) 
        idx_elem_global = finder(ex, ey) + ajuste
        print("Elemento:", idx_elem)
        
        mask = idx_elem  >= 0

        pontos_validos = pontos[mask]
        idx_validos = idx_elem[mask]
        
        
        # pega os valores correspondentes aos elementos onde caíram os pontos
        sigma_L_pts = sigma_L[idx_validos]
        sigma_T_pts = sigma_T[idx_validos]
        theta_pts   = theta_deg[idx_validos]
        
        dados_elipses = np.column_stack([
                                        idx_validos,
                                        pontos_validos[:, 0],
                                        pontos_validos[:, 1],
                                        sigma_L_pts,
                                        sigma_T_pts,
                                        theta_pts
                                    ])
        fig, ax = plt.subplots(figsize=(7, 7))

        #sigmaL_max = max(np.max(np.abs(sigma_L_pts)), 1e-12)
        #sigmaT_max = max(np.max(np.abs(sigma_T_pts)), 1e-12)
        sigma_max = max(np.max(np.abs(sigma_L_pts)), np.max(np.abs(sigma_T_pts)), 1e-12)
        sigma_min = min(np.min(np.abs(sigma_L_pts)), np.min(np.abs(sigma_T_pts)), 1e-12)
        escala = 0.005
        
        # ===== índice de anisotropia para todos os pontos =====
        #AI_all =  sigma_T_pts / sigma_L_pts
        
        #sigma_max = max(np.max(np.abs(AI_all)), 1e-12)
        
        #AI_min = np.min(sigma_max)
        #AI_max = np.max(sigma_max)
        '''
        # evita erro se todos os valores forem iguais
        if abs(AI_max - AI_min) < 1e-12:
            AI_min = AI_min - 1.0
            AI_max = AI_max + 1.0
        
        norm = colors.Normalize(
            vmin=sigma_min,
            vmax=sigma_max
        )
        '''
        norm = colors.Normalize(vmin=0.0, vmax=1.0)
        cmap = cm.jet
        
        for linha in dados_elipses:
        
            _, x0, y0, sL, sT, theta = linha
        
            AI = sT / sL
        
            a = escala * abs(sL) / sigma_max
            b = escala * abs(sT) / sigma_max
            #print(a,b)
            
        
            cor = cmap(norm(AI))
        
            elipse = Ellipse(
                xy=(x0, y0),
                width=2 * a,
                height=2 * b,
                angle=theta,
                fill=True,
                facecolor=cor,
                edgecolor='black',
                linewidth=0.3,
                alpha=0.8
            )
        
            ax.add_patch(elipse)
        
        # ===== colorbar =====
        sm = cm.ScalarMappable(
            cmap=cmap,
            norm=norm
            
        )
        
        sm.set_array([])
        
        cbar = plt.colorbar(sm, ax=ax, shrink=0.7)
        cbar.set_label(r"$AI = \sigma_L - \sigma_T$")
        
        theta_circ = np.linspace(0, 2*np.pi, 400)
        
        ax.plot(
            R_circ*np.cos(theta_circ),
            R_circ*np.sin(theta_circ),
            'k-',
            linewidth=1.0
        )
        
        ax.set_xlim(-R_circ - 0.02, R_circ + 0.02)
        ax.set_ylim(-R_circ - 0.02, R_circ + 0.02)
        
        ax.set_title(
            f"Phantom  Anisotropy ",
            fontsize=11
        )
        
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('[m]')
        ax.set_ylabel('[m]')
        ax.grid(False)
        
        #plt.show()
        
                
        if save:
            #plt.savefig(nome_arquivo, dpi=100)
            #plt.savefig(nome_arquivo, dpi=200,bbox_inches="tight")
            plt.savefig(f'{nome_arquivo}',  dpi=150, bbox_inches='tight', pil_kwargs={"quality": 70})
            plt.show() 
            plt.close()   # importante
        else:
            plt.show()    

    def Solve(self, forceKGolbalCalc=False):
        if self.imageSave == True:
            #self.plotMSH(self.mymesh.sigma_vec)
            #self.plotMSH(self.mymesh.sigma_vec[:, 0],  save = True, SigmaXXXYYY='xx')
            #self.plotMSH(self.mymesh.sigma_vec[:, 1],  save = True, SigmaXXXYYY='xy')
            #self.plotMSH(self.mymesh.sigma_vec[:, 2],  save = True, SigmaXXXYYY='yy')
            #### plotMSH(self, sigma, iteration = None, save = False, SigmaXXXYYY = None):
                
        
            Smed = 0.5 * (self.mymesh.sigma_vec[:, 0] + self.mymesh.sigma_vec[:, 2])
            D = np.sqrt(((self.mymesh.sigma_vec[:, 0] - self.mymesh.sigma_vec[:, 2])/2)**2 + self.mymesh.sigma_vec[:, 1]**2)
    
            sigma_L = Smed + D
            sigma_T = Smed - D
            
            
            DifAnisotropia = sigma_L - sigma_T
            DifAnisotropia_Med =  np.abs(np.mean(DifAnisotropia))
            
            # para evitar saturação para  +/- 90 devido a arctan
    
            #anisotropia = np.abs(sigma_L - sigma_T)       # verifica se a diferença Sxx - Syy é muito pequena
            #np.savetxt("sigmaInicialTheta.txt", sigmaInicial)
            
            #anisotropia = np.abs(sigmaInicial[:, 0] - sigmaInicial[:, 2])       # verifica se a diferença Sxx - Syy é muito pequena
            #np.savetxt("anisotropia.txt", anisotropia)  # formato binário
            
            sigma_theta = self.mymesh.sigma_vec[:, 0] - self.mymesh.sigma_vec[:, 2]
            #np.savetxt("sigma_theta.txt", sigma_theta)
            theta_rad = 0.5 * np.arctan2(2.0 * self.mymesh.sigma_vec[:, 1], sigma_theta)
            #np.savetxt("theta_rad.txt", theta_rad)
            theta_deg = np.rad2deg(theta_rad)                
                
            
            #nome1 =  f'{pasta_teste}/{html_name}_sigma_xx_{Lambda:.6f}.webp'
            nome1 = f"../../docs/{self.name}_sigma_xx.webp" 
            self.plotMSH(self.mymesh.sigma_vec[:, 0],  save = True, SigmaXXXYYY='xx', nome_arquivo=nome1)
            #lista_imgs.append(nome1)
            
            #nome2 =  f'{pasta_teste}/{html_name}_sigma_xy_{Lambda:.6f}.webp'
            nome2 =  f"../../docs/{self.name}_sigma_xy.webp"
            self.plotMSH(self.mymesh.sigma_vec[:, 1],  save = True, SigmaXXXYYY='xy', nome_arquivo=nome2)
            #lista_imgs.append(nome2)
            
            #nome3 = f'{pasta_teste}/{html_name}_sigma_yy_{Lambda:.6f}.webp'
            nome3 = f"../../docs/{self.name}_sigma_yy.webp"
            self.plotMSH(self.mymesh.sigma_vec[:, 2],  save = True, SigmaXXXYYY='yy', nome_arquivo=nome3)
            #lista_imgs.append(nome3)


            nome4 =  f"../../docs/{self.name}_Elipse.webp"     
            #self.plotElipse(self.mymesh.sigma_vec, sigma_L, sigma_T, theta_deg, Lambda, itr, save=True, DifAniso = DifAnisotropia_Med, nome_arquivo=nome4)
            self.plotElipse(self.mymesh.sigma_vec, sigma_L, sigma_T, theta_deg,  save=True, nome_arquivo=nome4)
            #lista_imgs.append(nome13)            
            
        #print('self.mymesh.sigma_vec',self.mymesh.sigma_vec)
        if (self.mymesh.KGlobal is None) or (forceKGolbalCalc):
            self.mymesh.CalcKGlobal()
        
        self.apply_boundary_conditions()
        #print('solve self.KGlobal \n',self.KGlobal)
        self.Yinversa = np.linalg.inv(self.KGlobal)
        #np.savetxt("self.Yinversa_fwd.txt", self.Yinversa, fmt="%.6f")
        self.Vmedido = np.dot(self.Yinversa, self.vetor_corrente_cond_contorno)
        #print(f' Tensões medidas em todos os nós \n {self.Vmedido})')
        
        #print('solve vetor_corrente_cond_contorno \n', self.vetor_corrente_cond_contorno)


        self.Vmedido_eletrodos = self.Vmedido[self.mymesh.ElectrodeNodes]
        #print(f' Tensões nos eletrodos \n {self.Vmedido_eletrodos})')


    ###############################################################################
    # Esta função cria arquivos .pos (Post-Processing) apara vizulização no Gmsh.
    #  opções: 'rho', 'sigma', 'V', 'I'
    ###############################################################################
    #def criar_arquivo_pos_2D(self, matriz_coordenadas, matriz_topologia, N_padraoCC, Post_Processing, nome_arquivo):
    def criar_arquivo_pos_2D(self, Post_Processing, nome_arquivo):

        matriz_coordenadas = self.mymesh.Coordinates
        matriz_topologia = self.mymesh.msh_topology
        dim = len(matriz_topologia[1])

        if dim == 2:
            self.n_PCurrent = 1
            header = "V(x) 1D"
            Sx = "SL"
            vetor_Post_Processing = Post_Processing
            arquivo = open('../../malhasPOS/'+ nome_arquivo + str(0) + '.pos', 'w+')
            arquivo.writelines('View "' + header + '" { \n')
            for i in range(0, len(matriz_topologia)):
            
                node_l = int(matriz_topologia[i][0]) # pegar dados do dataframe
                node_m = int(matriz_topologia[i][1])
                
                arquivo.writelines(Sx + '(')
                
                
                arquivo.writelines([str(matriz_coordenadas[node_l][0]),',', str(matriz_coordenadas[node_l][1]),',0,'])
                arquivo.writelines([str(matriz_coordenadas[node_m][0]),',', str(matriz_coordenadas[node_m][1]),',0'])
                
                arquivo.writelines(u')')
                arquivo.writelines(u'{')
                arquivo.writelines([str(Post_Processing[node_l]),',', str(Post_Processing[node_m])])
                arquivo.writelines(u'};')
                arquivo.writelines(u'\n')
            arquivo.writelines(u'};')
            
        else:
            header = "A list-based view"
            Sx = "ST"
            self.n_PCurrent = self.mymesh.NumberOfElectrodes
        

            for n_ele in range(self.n_PCurrent):
                vetor_Post_Processing = Post_Processing[:, n_ele]
                    
                arquivo = open('../../malhasPOS/'+ nome_arquivo + str(n_ele) + '.pos', 'w+')
                arquivo.writelines('View "' + header + '" { \n')
                
                for i in range(0, len(matriz_topologia)):
                
                    node_l = int(matriz_topologia[i][0]) # pegar dados do dataframe
                    node_m = int(matriz_topologia[i][1])
                    node_n = int(matriz_topologia[i][2]) 
                    arquivo.writelines(Sx + '(')
                    
                    
                    arquivo.writelines([str(matriz_coordenadas[node_l][0]),',', 
                                        str(matriz_coordenadas[node_l][1]),',0,'])
                    arquivo.writelines([str(matriz_coordenadas[node_m][0]),',', 
                                        str(matriz_coordenadas[node_m][1]),',0,'])
                    arquivo.writelines([str(matriz_coordenadas[node_n][0]),',', 
                                        str(matriz_coordenadas[node_n][1]),',0'])
                    
                    arquivo.writelines(u')')
                    arquivo.writelines(u'{')
    
                    arquivo.writelines([str(vetor_Post_Processing[node_l]),',', 
                                        str(vetor_Post_Processing[node_m]),',', 
                                        str(vetor_Post_Processing[node_n])])
            
                    arquivo.writelines(u'};')
                    arquivo.writelines(u'\n')
                arquivo.writelines(u'};')
            
        arquivo.close()  
##############################################################################    



###############################################################################
# Esta função abre arquivos .pos (Post-Processing), criados pela função
# 'criar_arquivo_pos'  para vizulização no Gmsh.
#  
###############################################################################
    def abrir_Gmsh_pos(self,nome_arquivo, runGmsh = False):
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        for pos in range(self.n_PCurrent):
            # Inicialize o Gmsh
            #gmsh.initialize()
            
            # Carregue o arquivo .geo
            #gmsh.open('D:/GIT_EIT_2D/EIT_2D/malhasPOS/' + nome_arquivo + str(pos) + '.pos') 
            gmsh.open('../../malhasPOS/' + nome_arquivo + str(pos) + '.pos')
            #gmsh.open( nome_arquivo + str(pos) + '.pos')
            #gmsh.open('../../malhasPOS/' + nome_arquivo + '.pos') 
            #gmsh.View[0].IntervalsType = 2;
            #gmsh.option.setNumber("View["+str(pos)+"].IntervalsType", 1)
            #gmsh.option.setNumber("View["+str(pos)+"].NbIso", 500)
           
            
        if runGmsh and '-nopopup' not in sys.argv:
            gmsh.fltk.run()
        gmsh.finalize()
###############################################################################
        

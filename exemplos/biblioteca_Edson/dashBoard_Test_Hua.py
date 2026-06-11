import math
import numpy as np
import gmsh
import sys
import time
import os
import shutil
import importlib
import os
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import mesh
import forwardProblem
import inverseProblem_2D_Anisotropic_Hua
import matplotlib.pyplot as plt



def safe_float(s, default=0.0):
    try:
        return float(s)
    except Exception:
        return default


def safe_int(s, default=0):
    try:
        return int(float(s))
    except Exception:
        return default


class ConfigDashboard(tk.Tk):
    def _row_tensor(self, parent, r, titulo, sxx, sxy, syy, theta):
        ttk.Label(parent, text=titulo).grid(
            row=r, column=0, sticky="w", pady=3
        )
    
        linha = ttk.Frame(parent)
        linha.grid(row=r, column=1, sticky="w", pady=3)
    
        campos = [
            ("σxx", sxx),
            ("σxy", sxy),
            ("σyy", syy),
            ("θ°", theta),
        ]
    
        for texto, var in campos:
            ttk.Label(linha, text=texto).pack(side="left", padx=(0, 2))
            ttk.Entry(linha, textvariable=var, width=5).pack(side="left", padx=(0, 8))
    
        return r + 1
    def __init__(self):
        super().__init__()
        self.title("EIT TEST ANISOTROPIC CONFIGURATION")
        self.geometry("520x940")

        # =========================================================
        # Variáveis - Forward problem / Mesh
        # =========================================================
        self.nome_arquivo = tk.StringVar(value="test_Hua_01")
        #self.raio = tk.StringVar(value="0.15")
        self.phaton_msh = tk.StringVar(value="circ16_anomZero_Hua")
        self.n_eletrodos = tk.StringVar(value="16")
        self.SkipPattern = tk.StringVar(value="3")
        #self.lc1 = tk.StringVar(value="2e-2")
        #self.lenth_e = tk.StringVar(value="0.02")
        #self.out_dir = tk.StringVar(value=str(Path.cwd() / "malhasMSH"))
        self.domain_sxx = tk.DoubleVar(value=1.0)
        self.domain_sxy = tk.DoubleVar(value=0.0)
        self.domain_syy = tk.DoubleVar(value=1.0)
        self.domain_theta = tk.DoubleVar(value=0.0)
        
        self.obj1_sxx = tk.DoubleVar(value=1.0)
        self.obj1_sxy = tk.DoubleVar(value=0.0)
        self.obj1_syy = tk.DoubleVar(value=1.0)
        self.obj1_theta = tk.DoubleVar(value=0.0)
        
        self.obj2_sxx = tk.DoubleVar(value=1.0)
        self.obj2_sxy = tk.DoubleVar(value=0.0)
        self.obj2_syy = tk.DoubleVar(value=1.0)
        self.obj2_theta = tk.DoubleVar(value=0.0)
        
        self.obj3_sxx = tk.DoubleVar(value=1.0)
        self.obj3_sxy = tk.DoubleVar(value=0.0)
        self.obj3_syy = tk.DoubleVar(value=1.0)
        self.obj3_theta = tk.DoubleVar(value=0.0)

        # =========================================================
        # Variáveis - Objects / Phantom
        # =========================================================
        self.NrAnomalias = tk.IntVar(value=0)

        self.anomalia1_raio = tk.StringVar(value="0.04")
        self.anomalia1_lados = tk.StringVar(value="16")
        self.anomalia1_rotacao = tk.StringVar(value="0.0")
        self.x1_ptoCentral = tk.StringVar(value="0.075")
        self.y1_ptoCentral = tk.StringVar(value="0.05")
        self.lc2 = tk.StringVar(value="2e-2")

        self.anomalia2_lados = tk.StringVar(value="4")
        self.anomalia2_raio = tk.StringVar(value="0.04")
        self.anomalia2_rotacao = tk.StringVar(value="45")
        self.x2_ptoCentral = tk.StringVar(value="-0.075")
        self.y2_ptoCentral = tk.StringVar(value="0.05")
        self.lc3 = tk.StringVar(value="2e-2")

        self.anomalia3_lados = tk.StringVar(value="5")
        self.anomalia3_raio = tk.StringVar(value="0.03")
        self.anomalia3_rotacao = tk.StringVar(value="0")
        self.x3_ptoCentral = tk.StringVar(value="0.0")
        self.y3_ptoCentral = tk.StringVar(value="-0.04")
        self.lc4 = tk.StringVar(value="2e-2")

        self.superficie_PCentral = tk.StringVar(value="Domain")

        # =========================================================
        # Variáveis - Inverse problem / Optimization
        # =========================================================
        self.lambda_reg = tk.StringVar(value="0.001")
        self.alpha = tk.StringVar(value="0.01")
        self.max_iter = tk.StringVar(value="20")
        self.tolerance = tk.StringVar(value="1e-6")
        self.initial_estimate = tk.StringVar(value="1.0")

        self.clip_sxx_min = tk.StringVar(value="0.1")
        self.clip_sxx_max = tk.StringVar(value="5.0")
        self.clip_sxy_min = tk.StringVar(value="-5.0")
        self.clip_sxy_max = tk.StringVar(value="5.0")
        self.clip_syy_min = tk.StringVar(value="0.1")
        self.clip_syy_max = tk.StringVar(value="5.0")

        self.use_gauss_newton = tk.IntVar(value=1)
        self.use_autoscale = tk.IntVar(value=1)
        self.plot_ellipses = tk.IntVar(value=0)

        # =========================================================
        # Layout principal
        # =========================================================
        main = ttk.Frame(self, padding=12)
        main.pack(fill="both", expand=True)
        main.columnconfigure(0, weight=1)

        # =========================================================
        # Frame 1 - Forward problem parameters
        # =========================================================
        forward_frame = ttk.LabelFrame(
            main,
            text="Forward problem parameters",
            padding=10
        )
        forward_frame.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        forward_frame.columnconfigure(1, weight=1)

        r = 0
        r = self._row(forward_frame, r, "File name:", self.nome_arquivo)
        r = self._row(forward_frame, r, "Phantom Mesh:", self.phaton_msh)
        r = self._row(forward_frame, r, "Nº electrodes:", self.n_eletrodos)
        r = self._row(forward_frame, r, "SkipPattern:", self.SkipPattern)
        #r = self._row(forward_frame, r, "lc (domain):", self.lc1)
        
        r = self._row_tensor(forward_frame, r,
            "Domain Anisotropy",
            self.domain_sxx,
            self.domain_sxy,
            self.domain_syy,
            self.domain_theta)
        r = self._row_tensor(forward_frame, r,
            "1st Object Anisotropy",
            self.obj1_sxx,
            self.obj1_sxy,
            self.obj1_syy,
            self.obj1_theta)
        
        r = self._row_tensor(forward_frame, r,
            "2nd Object Anisotropy",
            self.obj2_sxx,
            self.obj2_sxy,
            self.obj2_syy,
            self.obj2_theta)
        
        r = self._row_tensor(forward_frame, r,
            "3rd Object Anisotropy",
            self.obj3_sxx,
            self.obj3_sxy,
            self.obj3_syy,
            self.obj3_theta)

        out_frame = ttk.Frame(forward_frame)
        out_frame.grid(row=r, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        out_frame.columnconfigure(1, weight=1)

        #ttk.Label(out_frame, text="Output folder:").grid(row=0, column=0, sticky="w")
        #self.out_entry = ttk.Entry(out_frame, textvariable=self.out_dir)
        #self.out_entry.grid(row=0, column=1, sticky="ew", padx=(6, 6))
        #ttk.Button(out_frame, text="Choose...", command=self.choose_dir).grid(row=0, column=2, sticky="e")
        ''' 
        central_frame = ttk.Frame(forward_frame)
        central_frame.grid(row=r + 1, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        ttk.Label(central_frame, text="Central Point in:").pack(side="left")
        self.cmb_reg = ttk.Combobox(
            central_frame,
            textvariable=self.superficie_PCentral,
            values=["Domain", "Object 1", "Object 2", "Object 3"],
            state="readonly",
            width=14
        )
        self.cmb_reg.pack(side="left", padx=(8, 0))
        '''
        # =========================================================
        # Frame 2 - Objects / Phantom parameters
        # =========================================================
        objects_frame = ttk.LabelFrame(
            main,
            text="Objects / Phantom parameters",
            padding=10
        )
        objects_frame.grid(row=1, column=0, sticky="ew", pady=8)
        objects_frame.columnconfigure(1, weight=1)

        r = 0
        r = self._row(objects_frame, r, "Nr objects:", self.NrAnomalias)

        ttk.Separator(objects_frame, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8
        )
        r += 1
        '''
        r = self._row(objects_frame, r, "Object 1 - nr of sides:", self.anomalia1_lados)

        '''
        # =========================================================
        # Frame 3 - Inverse problem parameters
        # =========================================================
        inverse_frame = ttk.LabelFrame(
            main,
            text="Inverse problem parameters",
            padding=10
        )
        inverse_frame.grid(row=2, column=0, sticky="ew", pady=8)
        inverse_frame.columnconfigure(1, weight=1)

        r = 0
        r = self._row(inverse_frame, r, "Lambda:", self.lambda_reg)
        r = self._row(inverse_frame, r, "Alpha:", self.alpha)
        r = self._row(inverse_frame, r, "Max iterations:", self.max_iter)
        r = self._row(inverse_frame, r, "Tolerance:", self.tolerance)
        r = self._row(inverse_frame, r, "Initial estimate:", self.initial_estimate)

        ttk.Separator(inverse_frame, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8
        )
        r += 1

        r = self._row(inverse_frame, r, "Clip σxx min:", self.clip_sxx_min)
        r = self._row(inverse_frame, r, "Clip σxx max:", self.clip_sxx_max)
        r = self._row(inverse_frame, r, "Clip σxy min:", self.clip_sxy_min)
        r = self._row(inverse_frame, r, "Clip σxy max:", self.clip_sxy_max)
        r = self._row(inverse_frame, r, "Clip σyy min:", self.clip_syy_min)
        r = self._row(inverse_frame, r, "Clip σyy max:", self.clip_syy_max)

        check_frame = ttk.Frame(inverse_frame)
        check_frame.grid(row=r, column=0, columnspan=2, sticky="w", pady=(8, 0))

        ttk.Checkbutton(
            check_frame,
            text="Gauss-Newton",
            variable=self.use_gauss_newton
        ).pack(side="left", padx=(0, 10))

        ttk.Checkbutton(
            check_frame,
            text="AutoScale plots",
            variable=self.use_autoscale
        ).pack(side="left", padx=(0, 10))

        ttk.Checkbutton(
            check_frame,
            text="Plot ellipses",
            variable=self.plot_ellipses
        ).pack(side="left")

        # =========================================================
        # Botões principais
        # =========================================================
        btns = ttk.Frame(main)
        btns.grid(row=3, column=0, sticky="ew", pady=(12, 0))

        ttk.Button(
            btns,
            text="Update preview",
            command=self.update_preview
        ).pack(side="left")

        ttk.Button(
            btns,
            text="Save mesh",
            command=self.save_file
        ).pack(side="left", padx=(10, 0))


    def _row(self, parent, r, label, var):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=3)
        e = ttk.Entry(parent, textvariable=var, width=22)
        e.grid(row=r, column=1, sticky="w", padx=(8, 0), pady=3)
        return r + 1

    #def choose_dir(self):
    #    d = filedialog.askdirectory(title="Choose output folder")
    #    if d:
    #        self.out_dir.set(d)
   
    def update_preview(self):
    
        #nome = '../../malhasMSH/circ4_objetoUm_Hua_coarse.msh'
        #nome = '../../malhasMSH/Hua_cuba16eletrodos_3objetos.msh'
    
    
        #nome = '../../malhasMSH/Hua_cuba4eletrodos_1objetoDireita.msh'
        #nome = '../../malhasMSH/Hua_cuba16eletrodos_1objeto_denso.msh'
        #nome = '../../malhasMSH/test_Olavo_Hua.msh'
        #nome = '../../malhasMSH/circ16_anom1_Square_Hua_a_esquerda_denso.msh'
        nome = '../../malhasMSH/circ16_1object_Square_left_v1C.msh'
        #nome = '../../malhasMSH/circ16_2object_SqrCirc_Hua_v1.msh'
        
        
        nomePhanton = 'circ16_1object_Square_left__Elipse_ZeroDegree'
        MinhaMalha = mesh.HuaElectrodes2DAnisotropic(16, nome_msh=nome, altura2D = 0.02, thetaAngle = 0.0)#, sigmaX = 1.00, sigmaY = 1.0000)
        #MinhaMalha = mesh.HuaElectrodes2DAnisotropic(8, nome_msh=nome, altura2D = 0.02, thetaAngle = -45.0, sigmaX = 1000.00, sigmaY = 1.0)
    
        MinhaMalha.ReadMesh() 
    
        print('MinhaMalha.Elements[2]',MinhaMalha.Elements[2])
        print(f"Centroid: {MinhaMalha.Elements[2].Centroid}")
        #print(f"KGeo: \n{MinhaMalha.Elements[2].KGeo}")
    
    
        meus_sigmas = {
            1000: [3.0, 0.0, 3.0],
            1001: [3.0, 0.0, 1.0],
            1002: [1.0, 0.0, 3.0],
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
    
        fwd = forwardProblem.forward_problem(MinhaMalha, Pcorrente=None, SkipPattern=3, VirtualNode = True, I =1.0e-3, name = nomePhanton, imageSave = True)   # __init__ roda aqui
    
        mtz_Vmedido = fwd.Solve()
        print(f'Vmedido \n {fwd.Vmedido[:10]}')
    
        nome_arquivo = 'ParaVernoGmshPto'
        fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)
        fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)
    
        V_measured_phaton = fwd.Vmedido_eletrodos
        np.save("V_measured_phaton.npy", V_measured_phaton)  # formato binário
        print(f'V_mesured\n {V_measured_phaton}')
        print(f'meus_sigmas\n {meus_sigmas}')
    
    
    
    #runFWD_InverseProblemAnisotropicHua()
    

        
        '''
        gmsh.model.add("Domain_Hua")

        pts_circulo = []
        
        self.dtheta_e = float(self.lenth_e.get()) / float(self.raio.get())
       
        print('dtheta_e', self.dtheta_e)
        for k in range(int(self.n_eletrodos.get())):
            #passo      = 2*math.pi / n_eletrodos
            #dtheta_e   = lenth_e / raio
            theta_c = k * (2*math.pi /int(self.n_eletrodos.get()))
            # lista para guardar os pontos do eletrodo k
            pts_eletrodo = []
            
            
            # ponto A
            x_circ = float(self.raio.get())*math.cos(theta_c - self.dtheta_e/2)
            y_circ = float(self.raio.get())*math.sin(theta_c - self.dtheta_e/2)
            pA = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pA)
            pts_eletrodo.append(pA)
            
            # ponto B
            x_circ = float(self.raio.get())*math.cos(theta_c - self.dtheta_e/4)
            y_circ = float(self.raio.get())*math.sin(theta_c - self.dtheta_e/4)
            pB = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pB)
            pts_eletrodo.append(pB)
            
            
            
            # ponto C
            x_circ = float(self.raio.get())*math.cos(theta_c)
            y_circ = float(self.raio.get())*math.sin(theta_c)
            pC = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pC)
            pts_eletrodo.append(pC)
            print('pts_eletrodo', pts_eletrodo)
            
            # ponto D
            x_circ = float(self.raio.get())*math.cos(theta_c + self.dtheta_e/4)
            y_circ = float(self.raio.get())*math.sin(theta_c + self.dtheta_e/4)
            pD = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pD)
            pts_eletrodo.append(pD)
            
            # ponto E
            x_circ = float(self.raio.get())*math.cos(theta_c + self.dtheta_e/2)
            y_circ = float(self.raio.get())*math.sin(theta_c + self.dtheta_e/2)
            pE = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pE)
            pts_eletrodo.append(pE)
            
        ptoCentral = gmsh.model.geo.addPoint(0,0, 0, float(self.lc1.get())) 
        #gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='ptoCentral')    

        print('pts_circulo',pts_circulo )
        
        ###############################################################################
        # criando os pontos da "anomalia"
        ###############################################################################
        if int(self.NrAnomalias.get()) > 3: # verifica nr máximo de  objetos.
            raise TypeError("The maximum number of anomalies is 3.")
        if int(self.NrAnomalias.get()) == 1: 
            pts_anomalia1 = []
            for anom in range(int(self.anomalia1_lados.get())):
                x_anomalia1 = (float(self.x1_ptoCentral.get()) 
                              + float(self.anomalia1_raio.get())*math.cos((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                              + math.radians(float(self.anomalia1_rotacao.get()))))
                y_anomalia1 = (float(self.y1_ptoCentral.get())  
                              + float(self.anomalia1_raio.get())*math.sin((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                                                       + math.radians(float(self.anomalia1_rotacao.get()))))
                pts_anomalia1.append(gmsh.model.geo.addPoint(x_anomalia1,y_anomalia1,0.0,float(self.lc2.get())))
        if int(self.NrAnomalias.get()) == 2: 
            pts_anomalia1 = []
            for anom in range(int(self.anomalia1_lados.get())):
                x_anomalia1 = (float(self.x1_ptoCentral.get()) 
                              + float(self.anomalia1_raio.get())*math.cos((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                              + math.radians(float(self.anomalia1_rotacao.get()))))
                y_anomalia1 = (float(self.y1_ptoCentral.get())  
                              + float(self.anomalia1_raio.get())*math.sin((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                                                       + math.radians(float(self.anomalia1_rotacao.get()))))
                pts_anomalia1.append(gmsh.model.geo.addPoint(x_anomalia1,y_anomalia1,0.0,float(self.lc2.get())))
                
            pts_anomalia2 = []
            for anom2 in range(int(self.anomalia2_lados.get())):
                x_anomalia2 = (float(self.x2_ptoCentral.get()) 
                              + float(self.anomalia2_raio.get())*math.cos((2*np.pi*anom2/(int(self.anomalia2_lados.get()))) 
                              + math.radians(float(self.anomalia2_rotacao.get()))))
                y_anomalia2 = (float(self.y2_ptoCentral.get())  
                              + float(self.anomalia2_raio.get())*math.sin((2*np.pi*anom2/(int(self.anomalia2_lados.get()))) 
                                                       + math.radians(float(self.anomalia2_rotacao.get()))))
                pts_anomalia2.append(gmsh.model.geo.addPoint(x_anomalia2,y_anomalia2,0.0,float(self.lc3.get())))
            
        if int(self.NrAnomalias.get()) == 3:
            pts_anomalia1 = []
            for anom in range(int(self.anomalia1_lados.get())):
                x_anomalia1 = (float(self.x1_ptoCentral.get()) 
                              + float(self.anomalia1_raio.get())*math.cos((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                              + math.radians(float(self.anomalia1_rotacao.get()))))
                y_anomalia1 = (float(self.y1_ptoCentral.get())  
                              + float(self.anomalia1_raio.get())*math.sin((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                                                       + math.radians(float(self.anomalia1_rotacao.get()))))
                pts_anomalia1.append(gmsh.model.geo.addPoint(x_anomalia1,y_anomalia1,0.0,float(self.lc2.get())))
                
            pts_anomalia2 = []
            for anom2 in range(int(self.anomalia2_lados.get())):
                x_anomalia2 = (float(self.x2_ptoCentral.get()) 
                              + float(self.anomalia2_raio.get())*math.cos((2*np.pi*anom2/(int(self.anomalia2_lados.get()))) 
                              + math.radians(float(self.anomalia2_rotacao.get()))))
                y_anomalia2 = (float(self.y2_ptoCentral.get())  
                              + float(self.anomalia2_raio.get())*math.sin((2*np.pi*anom2/(int(self.anomalia2_lados.get()))) 
                                                       + math.radians(float(self.anomalia2_rotacao.get()))))
                pts_anomalia2.append(gmsh.model.geo.addPoint(x_anomalia2,y_anomalia2,0.0,float(self.lc3.get())))
            pts_anomalia3 = []
            for anom3 in range(int(self.anomalia3_lados.get())):
                x_anomalia3 = (float(self.x3_ptoCentral.get()) 
                              + float(self.anomalia3_raio.get())*math.cos((2*np.pi*anom3/(int(self.anomalia3_lados.get()))) 
                              + math.radians(float(self.anomalia3_rotacao.get()))))
                y_anomalia3 = (float(self.y3_ptoCentral.get())  
                              + float(self.anomalia3_raio.get())*math.sin((2*np.pi*anom3/(int(self.anomalia3_lados.get()))) 
                                                       + math.radians(float(self.anomalia3_rotacao.get()))))
                pts_anomalia3.append(gmsh.model.geo.addPoint(x_anomalia3,y_anomalia3,0.0,float(self.lc4.get())))
            
        else: pass  
        
        
        
        
        
        

        n_eletrodosHua = int(self.n_eletrodos.get())*5
        print('n_eletrodosHua', n_eletrodosHua)
        # ------------------------------------------------------------
        # 1) Curvas do contorno externo (apenas arcos do círculo)
        # ------------------------------------------------------------
        linhas_circulo = []
        for Lc in range(n_eletrodosHua):
            L_circulo = gmsh.model.geo.addCircleArc(
                pts_circulo[Lc], ptoCentral, pts_circulo[(Lc+1) % n_eletrodosHua]
            )
            linhas_circulo.append(L_circulo)
        loop_circulo = gmsh.model.geo.addCurveLoop(linhas_circulo)
        

        #if int(self.NrAnomalias.get()) == 0:
        #    surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo])
            
        if int(self.NrAnomalias.get()) == 1: 
            linhas_anomalia1 = []
            for A in range(int(self.anomalia1_lados.get())):
                L_anomalia1 = gmsh.model.geo.addLine(
                    pts_anomalia1[A], pts_anomalia1[(A+1) % int(self.anomalia1_lados.get())]
                )
                linhas_anomalia1.append(L_anomalia1)

            loop_anomalia1 = gmsh.model.geo.addCurveLoop(linhas_anomalia1)
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo, loop_anomalia1])
        if int(self.NrAnomalias.get()) == 2:
            linhas_anomalia1 = []
            for A in range(int(self.anomalia1_lados.get())):
                L_anomalia1 = gmsh.model.geo.addLine(
                    pts_anomalia1[A], pts_anomalia1[(A+1) % int(self.anomalia1_lados.get())]
                )
                linhas_anomalia1.append(L_anomalia1)

            loop_anomalia1 = gmsh.model.geo.addCurveLoop(linhas_anomalia1)
            linhas_anomalia2 = []
            for B in range(int(self.anomalia2_lados.get())):
                L_anomalia2 = gmsh.model.geo.addLine(
                    pts_anomalia2[B], pts_anomalia2[(B+1) % int(self.anomalia2_lados.get())]
                )
                linhas_anomalia2.append(L_anomalia2)

            loop_anomalia2 = gmsh.model.geo.addCurveLoop(linhas_anomalia2)
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo, loop_anomalia1, loop_anomalia2])
        if int(self.NrAnomalias.get()) == 3:
            linhas_anomalia1 = []
            for A in range(int(self.anomalia1_lados.get())):
                L_anomalia1 = gmsh.model.geo.addLine(
                    pts_anomalia1[A], pts_anomalia1[(A+1) % int(self.anomalia1_lados.get())]
                )
                linhas_anomalia1.append(L_anomalia1)

            loop_anomalia1 = gmsh.model.geo.addCurveLoop(linhas_anomalia1)
            linhas_anomalia2 = []
            for B in range(int(self.anomalia2_lados.get())):
                L_anomalia2 = gmsh.model.geo.addLine(
                    pts_anomalia2[B], pts_anomalia2[(B+1) % int(self.anomalia2_lados.get())]
                )
                linhas_anomalia2.append(L_anomalia2)

            loop_anomalia2 = gmsh.model.geo.addCurveLoop(linhas_anomalia2)
            # cuba COM FURO (interface conformal)
            linhas_anomalia3 = []
            for C in range(int(self.anomalia3_lados.get())):
                L_anomalia3 = gmsh.model.geo.addLine(
                    pts_anomalia3[C], pts_anomalia3[(C+1) % int(self.anomalia3_lados.get())]
                )
                linhas_anomalia3.append(L_anomalia3)

            loop_anomalia3 = gmsh.model.geo.addCurveLoop(linhas_anomalia3)
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo, loop_anomalia1, loop_anomalia2,  loop_anomalia3])
            
            # anomalia preenchendo o buraco
        if int(self.NrAnomalias.get()) == 0:
           surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo]) 
        if int(self.NrAnomalias.get()) == 1:
            surface_anomalia1 = gmsh.model.geo.addPlaneSurface([loop_anomalia1])
            
        if int(self.NrAnomalias.get()) == 2:
            surface_anomalia1 = gmsh.model.geo.addPlaneSurface([loop_anomalia1])
            surface_anomalia2 = gmsh.model.geo.addPlaneSurface([loop_anomalia2])
            
        if int(self.NrAnomalias.get()) == 3:
            surface_anomalia1 = gmsh.model.geo.addPlaneSurface([loop_anomalia1])
            surface_anomalia2 = gmsh.model.geo.addPlaneSurface([loop_anomalia2])
            surface_anomalia3 = gmsh.model.geo.addPlaneSurface([loop_anomalia3])
            


        
        # cuba sem anomalia
        #else:
        #    surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo])
        
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        # criando physicalGroup da superfície
        physycal_circulo = gmsh.model.addPhysicalGroup(2, [surface_circulo], 1000)
        gmsh.model.setPhysicalName(2, physycal_circulo, "Domain")
        
        if int(self.NrAnomalias.get()) == 1:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "Object_1")
            
            
        if int(self.NrAnomalias.get()) == 2:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "Object_1")            
            physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
            gmsh.model.setPhysicalName(2, physycal_anomalia2, "Object_2")
                
                
        if int(self.NrAnomalias.get()) == 3:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "Object_1")            
            physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
            gmsh.model.setPhysicalName(2, physycal_anomalia2, "Object_2")          
            physycal_anomalia3 = gmsh.model.addPhysicalGroup(2, [surface_anomalia3], 1003)
            gmsh.model.setPhysicalName(2, physycal_anomalia3, "Object_3")
        
        
        tipo = self.superficie_PCentral.get()
        print("Superfície escolhida:", tipo)
        
        match tipo:
            case "Domain":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_circulo)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Object 1":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia1)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Object 2":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia2)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Object 3":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia3)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
                
                
        ###############################################################################
        # criando o physycalGroup das curvas (1D) dos eletrodos
        # 
        ###############################################################################
        linhas_ele=[]
        count = 0
        for L_e in range(0, n_eletrodosHua,5):
            count = count+1
            print('count', count)
            linhaA = linhas_circulo[L_e]
            linhaB = linhas_circulo[L_e+1]
            linhaC = linhas_circulo[L_e+2]
            linhaD = linhas_circulo[L_e+3]
            gmsh.model.addPhysicalGroup(1, [linhaA, linhaB, linhaC, linhaD], 5000+count , 
                                        name=f"curvaElectrode_{count}")
        ###############################################################################   
        
        gmsh.model.geo.synchronize()
        #gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_circulo)
        gmsh.model.mesh.generate(2)
        #gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
        

        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()

        gmsh.finalize()
        '''    
###############################################################################
# Botão "Salvar"
# 
###############################################################################
    def save_file(self):
        gmsh.initialize()
    '''
    def save_file(self):


        gmsh.initialize()

        gmsh.model.add("Domain_Hua")

        pts_circulo = []
        
        self.dtheta_e = float(self.lenth_e.get()) / float(self.raio.get())
       
        print('dtheta_e', self.dtheta_e)
        for k in range(int(self.n_eletrodos.get())):
            #passo      = 2*math.pi / n_eletrodos
            #dtheta_e   = lenth_e / raio
            theta_c = k * (2*math.pi /int(self.n_eletrodos.get()))
            # lista para guardar os pontos do eletrodo k
            pts_eletrodo = []
            
            
            # ponto A
            x_circ = float(self.raio.get())*math.cos(theta_c - self.dtheta_e/2)
            y_circ = float(self.raio.get())*math.sin(theta_c - self.dtheta_e/2)
            pA = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pA)
            pts_eletrodo.append(pA)
            
            # ponto B
            x_circ = float(self.raio.get())*math.cos(theta_c - self.dtheta_e/4)
            y_circ = float(self.raio.get())*math.sin(theta_c - self.dtheta_e/4)
            pB = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pB)
            pts_eletrodo.append(pB)
            
            
            
            # ponto C
            x_circ = float(self.raio.get())*math.cos(theta_c)
            y_circ = float(self.raio.get())*math.sin(theta_c)
            pC = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pC)
            pts_eletrodo.append(pC)
            print('pts_eletrodo', pts_eletrodo)
            
            # ponto D
            x_circ = float(self.raio.get())*math.cos(theta_c + self.dtheta_e/4)
            y_circ = float(self.raio.get())*math.sin(theta_c + self.dtheta_e/4)
            pD = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pD)
            pts_eletrodo.append(pD)
            
            # ponto E
            x_circ = float(self.raio.get())*math.cos(theta_c + self.dtheta_e/2)
            y_circ = float(self.raio.get())*math.sin(theta_c + self.dtheta_e/2)
            pE = gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get()))
            pts_circulo.append(pE)
            pts_eletrodo.append(pE)
            
        ptoCentral = gmsh.model.geo.addPoint(0,0, 0, float(self.lc1.get())) 
        #gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='ptoCentral')    

        print('pts_circulo',pts_circulo )
        
        ###############################################################################
        # criando os pontos da "anomalia"
        ###############################################################################
        if int(self.NrAnomalias.get()) > 3: # verifica nr máximo de  objetos.
            raise TypeError("The maximum number of anomalies is 3.")
        if int(self.NrAnomalias.get()) == 1: 
            pts_anomalia1 = []
            for anom in range(int(self.anomalia1_lados.get())):
                x_anomalia1 = (float(self.x1_ptoCentral.get()) 
                              + float(self.anomalia1_raio.get())*math.cos((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                              + math.radians(float(self.anomalia1_rotacao.get()))))
                y_anomalia1 = (float(self.y1_ptoCentral.get())  
                              + float(self.anomalia1_raio.get())*math.sin((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                                                       + math.radians(float(self.anomalia1_rotacao.get()))))
                pts_anomalia1.append(gmsh.model.geo.addPoint(x_anomalia1,y_anomalia1,0.0,float(self.lc2.get())))
        if int(self.NrAnomalias.get()) == 2: 
            pts_anomalia1 = []
            for anom in range(int(self.anomalia1_lados.get())):
                x_anomalia1 = (float(self.x1_ptoCentral.get()) 
                              + float(self.anomalia1_raio.get())*math.cos((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                              + math.radians(float(self.anomalia1_rotacao.get()))))
                y_anomalia1 = (float(self.y1_ptoCentral.get())  
                              + float(self.anomalia1_raio.get())*math.sin((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                                                       + math.radians(float(self.anomalia1_rotacao.get()))))
                pts_anomalia1.append(gmsh.model.geo.addPoint(x_anomalia1,y_anomalia1,0.0,float(self.lc2.get())))
                
            pts_anomalia2 = []
            for anom2 in range(int(self.anomalia2_lados.get())):
                x_anomalia2 = (float(self.x2_ptoCentral.get()) 
                              + float(self.anomalia2_raio.get())*math.cos((2*np.pi*anom2/(int(self.anomalia2_lados.get()))) 
                              + math.radians(float(self.anomalia2_rotacao.get()))))
                y_anomalia2 = (float(self.y2_ptoCentral.get())  
                              + float(self.anomalia2_raio.get())*math.sin((2*np.pi*anom2/(int(self.anomalia2_lados.get()))) 
                                                       + math.radians(float(self.anomalia2_rotacao.get()))))
                pts_anomalia2.append(gmsh.model.geo.addPoint(x_anomalia2,y_anomalia2,0.0,float(self.lc3.get())))
            
        if int(self.NrAnomalias.get()) == 3:
            pts_anomalia1 = []
            for anom in range(int(self.anomalia1_lados.get())):
                x_anomalia1 = (float(self.x1_ptoCentral.get()) 
                              + float(self.anomalia1_raio.get())*math.cos((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                              + math.radians(float(self.anomalia1_rotacao.get()))))
                y_anomalia1 = (float(self.y1_ptoCentral.get())  
                              + float(self.anomalia1_raio.get())*math.sin((2*np.pi*anom/(int(self.anomalia1_lados.get()))) 
                                                       + math.radians(float(self.anomalia1_rotacao.get()))))
                pts_anomalia1.append(gmsh.model.geo.addPoint(x_anomalia1,y_anomalia1,0.0,float(self.lc2.get())))
                
            pts_anomalia2 = []
            for anom2 in range(int(self.anomalia2_lados.get())):
                x_anomalia2 = (float(self.x2_ptoCentral.get()) 
                              + float(self.anomalia2_raio.get())*math.cos((2*np.pi*anom2/(int(self.anomalia2_lados.get()))) 
                              + math.radians(float(self.anomalia2_rotacao.get()))))
                y_anomalia2 = (float(self.y2_ptoCentral.get())  
                              + float(self.anomalia2_raio.get())*math.sin((2*np.pi*anom2/(int(self.anomalia2_lados.get()))) 
                                                       + math.radians(float(self.anomalia2_rotacao.get()))))
                pts_anomalia2.append(gmsh.model.geo.addPoint(x_anomalia2,y_anomalia2,0.0,float(self.lc3.get())))
            pts_anomalia3 = []
            for anom3 in range(int(self.anomalia3_lados.get())):
                x_anomalia3 = (float(self.x3_ptoCentral.get()) 
                              + float(self.anomalia3_raio.get())*math.cos((2*np.pi*anom3/(int(self.anomalia3_lados.get()))) 
                              + math.radians(float(self.anomalia3_rotacao.get()))))
                y_anomalia3 = (float(self.y3_ptoCentral.get())  
                              + float(self.anomalia3_raio.get())*math.sin((2*np.pi*anom3/(int(self.anomalia3_lados.get()))) 
                                                       + math.radians(float(self.anomalia3_rotacao.get()))))
                pts_anomalia3.append(gmsh.model.geo.addPoint(x_anomalia3,y_anomalia3,0.0,float(self.lc4.get())))
            
        else: pass  
        
        
        
        n_eletrodosHua = int(self.n_eletrodos.get())*5
        print('n_eletrodosHua', n_eletrodosHua)
        # ------------------------------------------------------------
        # 1) Curvas do contorno externo (apenas arcos do círculo)
        # ------------------------------------------------------------
        linhas_circulo = []
        for Lc in range(n_eletrodosHua):
            L_circulo = gmsh.model.geo.addCircleArc(
                pts_circulo[Lc], ptoCentral, pts_circulo[(Lc+1) % n_eletrodosHua]
            )
            linhas_circulo.append(L_circulo)
        loop_circulo = gmsh.model.geo.addCurveLoop(linhas_circulo)
        
        linhas_anomalia = None
        loop_anomalia = None
        surface_anomalia = None
        #if int(self.NrAnomalias.get()) == 0:
        #    surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo])
            
        if int(self.NrAnomalias.get()) == 1: 
            linhas_anomalia1 = []
            for A in range(int(self.anomalia1_lados.get())):
                L_anomalia1 = gmsh.model.geo.addLine(
                    pts_anomalia1[A], pts_anomalia1[(A+1) % int(self.anomalia1_lados.get())]
                )
                linhas_anomalia1.append(L_anomalia1)

            loop_anomalia1 = gmsh.model.geo.addCurveLoop(linhas_anomalia1)
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo, loop_anomalia1])
        if int(self.NrAnomalias.get()) == 2:
            linhas_anomalia1 = []
            for A in range(int(self.anomalia1_lados.get())):
                L_anomalia1 = gmsh.model.geo.addLine(
                    pts_anomalia1[A], pts_anomalia1[(A+1) % int(self.anomalia1_lados.get())]
                )
                linhas_anomalia1.append(L_anomalia1)

            loop_anomalia1 = gmsh.model.geo.addCurveLoop(linhas_anomalia1)
            linhas_anomalia2 = []
            for B in range(int(self.anomalia2_lados.get())):
                L_anomalia2 = gmsh.model.geo.addLine(
                    pts_anomalia2[B], pts_anomalia2[(B+1) % int(self.anomalia2_lados.get())]
                )
                linhas_anomalia2.append(L_anomalia2)

            loop_anomalia2 = gmsh.model.geo.addCurveLoop(linhas_anomalia2)
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo, loop_anomalia1, loop_anomalia2])
        if int(self.NrAnomalias.get()) == 3:
            linhas_anomalia1 = []
            for A in range(int(self.anomalia1_lados.get())):
                L_anomalia1 = gmsh.model.geo.addLine(
                    pts_anomalia1[A], pts_anomalia1[(A+1) % int(self.anomalia1_lados.get())]
                )
                linhas_anomalia1.append(L_anomalia1)

            loop_anomalia1 = gmsh.model.geo.addCurveLoop(linhas_anomalia1)
            linhas_anomalia2 = []
            for B in range(int(self.anomalia2_lados.get())):
                L_anomalia2 = gmsh.model.geo.addLine(
                    pts_anomalia2[B], pts_anomalia2[(B+1) % int(self.anomalia2_lados.get())]
                )
                linhas_anomalia2.append(L_anomalia2)

            loop_anomalia2 = gmsh.model.geo.addCurveLoop(linhas_anomalia2)
            # cuba COM FURO (interface conformal)
            linhas_anomalia3 = []
            for C in range(int(self.anomalia3_lados.get())):
                L_anomalia3 = gmsh.model.geo.addLine(
                    pts_anomalia3[C], pts_anomalia3[(C+1) % int(self.anomalia3_lados.get())]
                )
                linhas_anomalia3.append(L_anomalia3)

            loop_anomalia3 = gmsh.model.geo.addCurveLoop(linhas_anomalia3)
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo, loop_anomalia1, loop_anomalia2,  loop_anomalia3])
            
            # anomalia preenchendo o buraco
        if int(self.NrAnomalias.get()) == 0:
           surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo]) 
        if int(self.NrAnomalias.get()) == 1:
            surface_anomalia1 = gmsh.model.geo.addPlaneSurface([loop_anomalia1])
            
        if int(self.NrAnomalias.get()) == 2:
            surface_anomalia1 = gmsh.model.geo.addPlaneSurface([loop_anomalia1])
            surface_anomalia2 = gmsh.model.geo.addPlaneSurface([loop_anomalia2])
            
        if int(self.NrAnomalias.get()) == 3:
            surface_anomalia1 = gmsh.model.geo.addPlaneSurface([loop_anomalia1])
            surface_anomalia2 = gmsh.model.geo.addPlaneSurface([loop_anomalia2])
            surface_anomalia3 = gmsh.model.geo.addPlaneSurface([loop_anomalia3])
            


        
        # cuba sem anomalia
        #else:
        #    surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo])
        
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        # criando physicalGroup da superfície
        physycal_circulo = gmsh.model.addPhysicalGroup(2, [surface_circulo], 1000)
        gmsh.model.setPhysicalName(2, physycal_circulo, "domain")


        if int(self.NrAnomalias.get()) == 1:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "Object_1")
            
            
        if int(self.NrAnomalias.get()) == 2:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "Object_1")            
            physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
            gmsh.model.setPhysicalName(2, physycal_anomalia2, "Object_2")
                
                
        if int(self.NrAnomalias.get()) == 3:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "Object_1")            
            physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
            gmsh.model.setPhysicalName(2, physycal_anomalia2, "Object_2")          
            physycal_anomalia3 = gmsh.model.addPhysicalGroup(2, [surface_anomalia3], 1003)
            gmsh.model.setPhysicalName(2, physycal_anomalia3, "Object_3")
        
        
        tipo = self.superficie_PCentral.get()
        print("Superfície escolhida:", tipo)
        
        match tipo:
            case "Domain":
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Object 1":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia1)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Object 2":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia2)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Object 3":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia3)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')


                
        ###############################################################################
        # criando o physycalGroup das curvas (1D) dos eletrodos
        # 
        ###############################################################################
        linhas_ele=[]
        count = 0
        for L_e in range(0, n_eletrodosHua,5):
            count = count+1
            print('count', count)
            linhaA = linhas_circulo[L_e]
            linhaB = linhas_circulo[L_e+1]
            linhaC = linhas_circulo[L_e+2]
            linhaD = linhas_circulo[L_e+3]
            gmsh.model.addPhysicalGroup(1, [linhaA, linhaB, linhaC, linhaD], 5000+count , 
                                        name=f"curvaElectrode_{count}")
        ###############################################################################           
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_circulo)
        gmsh.model.mesh.generate(2)
        #gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
        
        caminho = os.path.join(self.out_dir.get(), self.nome_arquivo.get() + ".msh")
        gmsh.write(caminho)
        
        gmsh.write(self.nome_arquivo.get()  + '.msh')
        gmsh.write(self.nome_arquivo.get()  + ".geo_unrolled")
        os.rename(self.nome_arquivo.get()  + ".geo_unrolled", self.nome_arquivo.get()  + ".geo_unrolled")
        shutil.move(self.nome_arquivo.get()  + ".geo_unrolled",self.nome_arquivo.get()  + ".geo")
        # adicionando linha no final do .geo
        geo_path = self.nome_arquivo.get() + ".geo"

        with open(geo_path, "a", encoding="utf-8") as f:
            f.write("\nMesh.MshFileVersion = 2.2;\n")
                #if '-nopopup' not in sys.argv:
                #    gmsh.fltk.run()

        gmsh.finalize()
    '''
if __name__ == "__main__":
    # no Windows, Tkinter já vem junto com Python padrão
    app = ConfigDashboard()
    app.mainloop()

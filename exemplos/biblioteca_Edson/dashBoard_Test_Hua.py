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
    def _row(self, parent, r, titulo, var, width=50):
        ttk.Label(parent, text=titulo).grid(
            row=r, column=0, sticky="w", pady=3
        )
    
        ttk.Entry(parent, textvariable=var, width=width).grid(
            row=r, column=1, sticky="w", pady=3
        )
    
        return r + 1
    def __init__(self):
        super().__init__()
        self.title("EIT TEST ANISOTROPIC CONFIGURATION")
        self.geometry("820x740")

        # =========================================================
        # Variáveis - Forward problem / Mesh
        # =========================================================
        self.nome_arquivo = tk.StringVar(value="test_Hua_01")
        #self.raio = tk.StringVar(value="0.15")
        self.phaton_msh = tk.StringVar(value="../../malhasMSH/circ16_1object_Square_left_v1C.msh")
        self.n_eletrodos = tk.IntVar(value=16)
        self.SkipPattern = tk.IntVar(value=3)
        #self.lc1 = tk.StringVar(value="2e-2")
        #self.lenth_e = tk.StringVar(value="0.02")
        #self.out_dir = tk.StringVar(value=str(Path.cwd() / "malhasMSH"))
        self.domain_sxx = tk.DoubleVar(value=3.0)
        self.domain_sxy = tk.DoubleVar(value=0.0)
        self.domain_syy = tk.DoubleVar(value=3.0)
        self.domain_theta = tk.DoubleVar(value=0.0)
        
        self.obj1_sxx = tk.DoubleVar(value=3.0)
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
        r = self._row(forward_frame, r, "File name:", self.nome_arquivo, width=50)
        r = self._row(forward_frame, r, "Phantom Mesh:", self.phaton_msh, width=50)
        r = self._row(forward_frame, r, "Nº electrodes:", self.n_eletrodos, width=10)
        r = self._row(forward_frame, r, "SkipPattern:", self.SkipPattern, width=10)
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
        '''
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
        '''
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

    '''
    def _row(self, parent, r, label, var):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=3)
        e = ttk.Entry(parent, textvariable=var, width=22)
        e.grid(row=r, column=1, sticky="w", padx=(8, 0), pady=3)
        return r + 1
    
    def _row(self, parent, r, label, var, width=20):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=3)
    
        e = ttk.Entry(parent, textvariable=var, width=width)
        e.grid(row=r, column=1, sticky="ew", padx=(8, 0), pady=3)
    
        return r + 1
    '''
    #def choose_dir(self):
    #    d = filedialog.askdirectory(title="Choose output folder")
    #    if d:
    #        self.out_dir.set(d)
   
    def update_preview(self):
    
        #nome = '../../malhasMSH/circ16_1object_Square_left_v1C.msh'
        
        
        nomePhanton = 'LIXOOOOOO'
        MinhaMalha = mesh.HuaElectrodes2DAnisotropic(self.n_eletrodos.get(),
                                                     nome_msh=self.phaton_msh.get(), 
                                                     altura2D = 0.02, 
                                                     thetaAngle = 0.0)#, sigmaX = 1.00, sigmaY = 1.0000)

        MinhaMalha.ReadMesh() 
    
        print('MinhaMalha.Elements[2]',MinhaMalha.Elements[2])
        print(f"Centroid: {MinhaMalha.Elements[2].Centroid}")
        #print(f"KGeo: \n{MinhaMalha.Elements[2].KGeo}")
    
    
        meus_sigmas = {
            1000: [self.domain_sxx.get(), self.domain_sxy.get(), self.domain_syy.get()],
            1001: [self.obj1_sxx.get(), self.obj1_sxy.get(), self.obj1_syy.get()],
            1002: [self.obj2_sxx.get(), self.obj2_sxy.get(), self.obj2_syy.get()],
            1003: [self.obj3_sxx.get(), self.obj3_sxy.get(), self.obj3_syy.get()],
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
    
        fwd = forwardProblem.forward_problem(MinhaMalha, Pcorrente=None, SkipPattern=self.SkipPattern.get(), VirtualNode = True, I =1.0e-3, name = nomePhanton, imageSave = True)   # __init__ roda aqui
    
        mtz_Vmedido = fwd.Solve()
        print(f'Vmedido \n {fwd.Vmedido[:10]}')
    
        nome_arquivo = 'ParaVernoGmshPto'
        fwd.criar_arquivo_pos_2D( fwd.Vmedido, nome_arquivo)
        fwd.abrir_Gmsh_pos(nome_arquivo, runGmsh=True)
    
        V_measured_phaton = fwd.Vmedido_eletrodos
        np.save("V_measured_phaton.npy", V_measured_phaton)  # formato binário
        print(f'V_mesured\n {V_measured_phaton}')
        print(f'meus_sigmas\n {meus_sigmas}')
    

###############################################################################
# Botão "Salvar"
# 
###############################################################################
    def save_file(self):
        gmsh.initialize()
    
if __name__ == "__main__":
    # no Windows, Tkinter já vem junto com Python padrão
    app = ConfigDashboard()
    app.mainloop()

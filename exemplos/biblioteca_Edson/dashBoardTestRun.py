# -*- coding: utf-8 -*-
"""
Created on Sun May 24 20:45:33 2026

@author: rodgr
"""

import math
import numpy as np
import gmsh
import sys
import time
import os
import shutil
import importlib
import os
import mesh
import forwardProblem
import inverseProblem_2D_Anisotropic_Hua
import matplotlib.pyplot as plt

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
    def __init__(self):
        super().__init__()
        self.title("EIT TEST ANISOTROPIC CONFIGURATION")
        self.geometry("600x820")

        # ====== variáveis ======
        # ====== variáveis ======

        

        self.project_name = tk.StringVar(value="project01")
        self.phaton_msh = tk.StringVar(value="circ16_anomZero_Hua")
        self.n_eletrodos = tk.StringVar(value="16")
        self.SkipPattern = tk.StringVar(value="3")          
        #self.NrObjects = tk.IntVar(value=0)
        #self.comAnomalia = tk.StringVar(value='True') 
        
        self.domain_sxx = tk.DoubleVar(value=1.0)
        self.domain_sxy = tk.DoubleVar(value=0.0)
        self.domain_syy = tk.DoubleVar(value=1.0)
        self.domain_theta = tk.DoubleVar(value=0.0)
        
        self.Object_1_sxx = tk.DoubleVar(value=1.0)
        self.Object_1_sxy = tk.DoubleVar(value=0.0)
        self.Object_1_syy = tk.DoubleVar(value=1.0)
        self.Object_1_theta = tk.DoubleVar(value=0.0)
        
        self.Object_2_sxx = tk.DoubleVar(value=1.0)
        self.Object_2_sxy = tk.DoubleVar(value=0.0)
        self.Object_2_syy = tk.DoubleVar(value=1.0)
        self.Object_2_theta = tk.DoubleVar(value=0.0)
        
        self.Object_3_sxx = tk.DoubleVar(value=1.0)
        self.Object_3_sxy = tk.DoubleVar(value=0.0)
        self.Object_3_syy = tk.DoubleVar(value=1.0)
        self.Object_3_theta = tk.DoubleVar(value=0.0)
        
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

        self.out_dir = tk.StringVar(value=str(Path.cwd() / "malhasMSH"))

        # ====== layout ======
        main = ttk.Frame(self, padding=12)
        main.pack(fill="both", expand=True)

        form = ttk.LabelFrame(main, text="Forward problem parameters", padding=10)
        form.grid(row=0, column=0, sticky="nsew")



        main.columnconfigure(0, weight=0)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        # ====== campos ======
        
        r = 0

        
        r += 1
        r = self._row(form, r, "Project Name :", self.project_name)
        r = self._row(form, r, "Phantom Mesh:", self.phaton_msh)
        r = self._row(form, r, "Nº Electrodes:", self.n_eletrodos)
        r = self._row(form, r, "Skip Pattern:", self.SkipPattern)
        #r = self._row(form, r, "Nr Objects from 0 to 3", self.NrObjects)
        #ttk.Label(form, text="Nr Objects from 0 to 3").grid(row=r, column=0, sticky="w", padx=5, pady=3)


        










        
        # linha
        ttk.Label(form, text="Domain Anisotropy").grid(row=r, column=0, sticky="w", padx=0, pady=3)
        
        linha_domain = ttk.Frame(form)
        linha_domain.grid(row=r, column=1, sticky="w")
        
        campos = [
            ("σxx", self.domain_sxx),
            ("σxy", self.domain_sxy),
            ("σyy", self.domain_syy),
            ("θ°",  self.domain_theta),
        ]
        
        for texto, var in campos:
            ttk.Label(linha_domain, text=texto).pack(side="left", padx=(0, 2))
            ttk.Entry(linha_domain, textvariable=var, width=4).pack(side="left", padx=(0, 8))
            
        #ttk.Separator(form, orient="horizontal").grid(row=r, column=0, columnspan=2, sticky="ew", pady=8)
        r += 1
        # linha
        ttk.Label(form, text="1st Object Anisotropy").grid(row=r, column=0, sticky="w", padx=0, pady=3)
        
        linha_domain = ttk.Frame(form)
        linha_domain.grid(row=r, column=1, sticky="w")
        
        campos = [
            ("σxx", self.Object_1_sxx),
            ("σxy", self.Object_1_sxy),
            ("σyy", self.Object_1_syy),
            ("θ°",  self.Object_1_theta),
        ]
        
        for texto, var in campos:
            ttk.Label(linha_domain, text=texto).pack(side="left", padx=(0, 2))
            ttk.Entry(linha_domain, textvariable=var, width=4).pack(side="left", padx=(0, 8))

        r += 1

        # linha
        ttk.Label(form, text="2nd Object Anisotropy").grid(row=r, column=0, sticky="w", padx=0, pady=3)
        
        linha_domain = ttk.Frame(form)
        linha_domain.grid(row=r, column=1, sticky="w")
        
        campos = [
            ("σxx", self.Object_2_sxx),
            ("σxy", self.Object_2_sxy),
            ("σyy", self.Object_2_syy),
            ("θ°",  self.Object_2_theta),
        ]
        
        for texto, var in campos:
            ttk.Label(linha_domain, text=texto).pack(side="left", padx=(0, 2))
            ttk.Entry(linha_domain, textvariable=var, width=4).pack(side="left", padx=(0, 8))

        r += 1

        # linha
        ttk.Label(form, text="2nd Object Anisotropy").grid(row=r, column=0, sticky="w", padx=0, pady=3)
        
        linha_domain = ttk.Frame(form)
        linha_domain.grid(row=r, column=1, sticky="w")
        
        campos = [
            ("σxx", self.Object_3_sxx),
            ("σxy", self.Object_3_sxy),
            ("σyy", self.Object_3_syy),
            ("θ°",  self.Object_3_theta),
        ]
        
        for texto, var in campos:
            ttk.Label(linha_domain, text=texto).pack(side="left", padx=(0, 2))
            ttk.Entry(linha_domain, textvariable=var, width=4).pack(side="left", padx=(0, 8))

        r += 1



        

    def _row(self, parent, r, label, var):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=3)
        e = ttk.Entry(parent, textvariable=var, width=22)
        e.grid(row=r, column=1, sticky="w", padx=(8, 0), pady=3)
        return r + 1

    def choose_dir(self):
        d = filedialog.askdirectory(title="Escolher pasta de saída")
        if d:
            self.out_dir.set(d)





    #############################################################################################
    

    '''
    def update_preview(self):
        NrElec =16
        mshName = '../../malhasMSH/test_Olavo_Hua.msh'
            
        self.runFWD_InverseProblemAnisotropicHua(NrElec,mshName, thetaAngle = 0.0, imageSave = True)

    def save_file(self):
    '''
if __name__ == "__main__":
    # no Windows, Tkinter já vem junto com Python padrão
    app = ConfigDashboard()
    app.mainloop()
    
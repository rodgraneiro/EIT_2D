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
        self.title("Gerador de Malhas (.msh) - EIT")
        self.geometry("420x860")

        # ====== variáveis ======
        self.nome_arquivo = tk.StringVar(value="testeHua_01")
        self.raio = tk.StringVar(value="0.15")
        self.n_eletrodos = tk.StringVar(value="16")
        self.lc1 = tk.StringVar(value="2e-2")
        self.lenth_e = tk.StringVar(value="0.02")            
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

        self.out_dir = tk.StringVar(value=str(Path.cwd() / "malhasMSH"))

        # ====== layout ======
        main = ttk.Frame(self, padding=12)
        main.pack(fill="both", expand=True)

        form = ttk.LabelFrame(main, text="Parâmetros", padding=10)
        form.grid(row=0, column=0, sticky="nsew")



        main.columnconfigure(0, weight=0)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        # ====== campos ======
        r = 0
        r = self._row(form, r, "Nome do arquivo :", self.nome_arquivo)
        r = self._row(form, r, "Raio:", self.raio)
        r = self._row(form, r, "Nº eletrodos:", self.n_eletrodos)
        r = self._row(form, r, "lc (cuba):", self.lc1)
        r = self._row(form, r, "Nr anomalias", self.NrAnomalias)
        r = self._row(form, r, "Comprimento eletrodo", self.lenth_e)
        ttk.Separator(form, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8
        )
        

        r += 1

        r = self._row(form, r, "Anomalia1 - nº de lados:", self.anomalia1_lados)
        r = self._row(form, r, "Raio:", self.anomalia1_raio)
        r = self._row(form, r, "lc:", self.lc2)
        r = self._row(form, r, "Rotação (graus):", self.anomalia1_rotacao)
        r = self._row(form, r, "x_ptoCentral:", self.x1_ptoCentral)
        r = self._row(form, r, "y_ptoCentral:", self.y1_ptoCentral)
        ttk.Separator(form, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8)
        
        r += 1
        r = self._row(form, r, "Anomalia2 - nº de lados:", self.anomalia2_lados)
        r = self._row(form, r, "Raio:", self.anomalia2_raio)
        r = self._row(form, r, "lc:", self.lc3)
        r = self._row(form, r, "Rotação (graus):", self.anomalia2_rotacao)
        r = self._row(form, r, "x_ptoCentral:", self.x2_ptoCentral)
        r = self._row(form, r, "y_ptoCentral:", self.y2_ptoCentral)
        ttk.Separator(form, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=16)
        r += 1
        r = self._row(form, r, "Anomalia3 - nº de lados:", self.anomalia3_lados)
        r = self._row(form, r, "Raio:", self.anomalia3_raio)
        r = self._row(form, r, "lc :", self.lc4)
        r = self._row(form, r, "Rotação (graus):", self.anomalia3_rotacao)
        r = self._row(form, r, "x_ptoCentral:", self.x3_ptoCentral)
        r = self._row(form, r, "y_ptoCentral:", self.y3_ptoCentral)
        
        
        # saída
        out_frame = ttk.Frame(form)
        out_frame.grid(row=r, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        out_frame.columnconfigure(1, weight=1)

        ttk.Label(out_frame, text="Pasta de saída:").grid(row=0, column=0, sticky="w")
        self.out_entry = ttk.Entry(out_frame, textvariable=self.out_dir)
        self.out_entry.grid(row=0, column=1, sticky="ew", padx=(6, 6))
        ttk.Button(out_frame, text="Escolher...", command=self.choose_dir).grid(row=0, column=2, sticky="e")

        # botões
        btns = ttk.Frame(form)
        btns.grid(row=r + 1, column=0, columnspan=2, sticky="ew", pady=(12, 0))
        ttk.Button(btns, text="Atualizar prévia", command=self.update_preview).pack(side="left")
        ttk.Button(btns, text="Salvar malha", command=self.save_file).pack(side="left", padx=(10, 0))
        ########opcao_var = tk.StringVar(value="Laplace")

        
        
        # combobox 
        self.superficie_PCentral = tk.StringVar(value="Cuba")
        ttk.Label(btns, text="Pto Central em:").pack(side="left", padx=(6, 4))
        self.cmb_reg = ttk.Combobox(
            btns,  # <- pai correto
            textvariable=self.superficie_PCentral,
            values=["Cuba", "Anomalia 1", "Anomalia 2", "Anomalia 3"],
            state="readonly",
            width=14
        )
        self.cmb_reg.pack(side="left")

        

    def _row(self, parent, r, label, var):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=3)
        e = ttk.Entry(parent, textvariable=var, width=22)
        e.grid(row=r, column=1, sticky="w", padx=(8, 0), pady=3)
        return r + 1

    def choose_dir(self):
        d = filedialog.askdirectory(title="Escolher pasta de saída")
        if d:
            self.out_dir.set(d)
   
    def update_preview(self):

        gmsh.initialize()

        gmsh.model.add("Cuba_Hua")

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
            raise TypeError("O número máximo de anomalias é 3")
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

        if int(self.NrAnomalias.get()) == 1:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "anomalia1")
            
            
        if int(self.NrAnomalias.get()) == 2:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "anomalia1")            
            physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
            gmsh.model.setPhysicalName(2, physycal_anomalia2, "anomalia2")
                
                
        if int(self.NrAnomalias.get()) == 3:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "anomalia1")            
            physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
            gmsh.model.setPhysicalName(2, physycal_anomalia2, "anomalia2")          
            physycal_anomalia3 = gmsh.model.addPhysicalGroup(2, [surface_anomalia3], 1003)
            gmsh.model.setPhysicalName(2, physycal_anomalia3, "anomalia3")
        
        
        tipo = self.superficie_PCentral.get()
        print("Superfície escolhida:", tipo)
        
        match tipo:
            case "Cuba":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_circulo)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Anomalia 1":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia1)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Anomalia 2":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia2)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Anomalia 3":
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
###############################################################################
# Botão "Salvar"
# 
###############################################################################


    def save_file(self):


        gmsh.initialize()

        gmsh.model.add("Cuba_Hua")

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
            raise TypeError("O número máximo de anomalias é 3")
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


        if int(self.NrAnomalias.get()) == 1:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "anomalia1")
            
            
        if int(self.NrAnomalias.get()) == 2:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "anomalia1")            
            physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
            gmsh.model.setPhysicalName(2, physycal_anomalia2, "anomalia2")
                
                
        if int(self.NrAnomalias.get()) == 3:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "anomalia1")            
            physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
            gmsh.model.setPhysicalName(2, physycal_anomalia2, "anomalia2")          
            physycal_anomalia3 = gmsh.model.addPhysicalGroup(2, [surface_anomalia3], 1003)
            gmsh.model.setPhysicalName(2, physycal_anomalia3, "anomalia3")
        
        
        tipo = self.superficie_PCentral.get()
        print("Superfície escolhida:", tipo)
        
        match tipo:
            case "Cuba":
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Anomalia 1":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia1)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Anomalia 2":
                gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia2)
                gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
            case "Anomalia 3":
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
        
        #gmsh.write(self.nome_arquivo.get()  + '.msh')
        #gmsh.write(self.nome_arquivo.get()  + ".geo_unrolled")
        #os.rename(self.nome_arquivo.get()  + ".geo_unrolled", self.nome_arquivo.get()  + ".geo_unrolled")
        #shutil.move(self.nome_arquivo.get()  + ".geo_unrolled",self.nome_arquivo.get()  + ".geo")

        #if '-nopopup' not in sys.argv:
        #    gmsh.fltk.run()

        gmsh.finalize()
    
if __name__ == "__main__":
    # no Windows, Tkinter já vem junto com Python padrão
    app = ConfigDashboard()
    app.mainloop()

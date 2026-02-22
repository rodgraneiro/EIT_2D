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
        self.title("Gerador de Malha (.msh) - EIT")
        self.geometry("420x820")

        # ====== variáveis ======
        self.nome_arquivo = tk.StringVar(value="teste01")
        self.raio = tk.StringVar(value="0.15")
        self.n_eletrodos = tk.StringVar(value="16")
        self.lc1 = tk.StringVar(value="2e-2")
        self.comAnomalia = tk.StringVar(value='True')            
        self.NrAnomalias = tk.IntVar(value=0)
        
        self.anomalia1_raio = tk.StringVar(value="0.05")
        self.anomalia1_lados = tk.StringVar(value="16")
        self.anomalia1_rotacao = tk.StringVar(value="22.5")
        self.x1_ptoCentral = tk.StringVar(value="0.06")
        self.y1_ptoCentral = tk.StringVar(value="0.05")
        self.lc2 = tk.StringVar(value="2e-2")

        self.anomalia2_lados = tk.StringVar(value="0")
        self.anomalia2_raio = tk.StringVar(value="0.0")
        self.anomalia2_rotacao = tk.StringVar(value="0")
        self.x2_ptoCentral = tk.StringVar(value="0.0")
        self.y2_ptoCentral = tk.StringVar(value="0.0")
        self.lc3 = tk.StringVar(value="2e-2")
        
        self.anomalia3_lados = tk.StringVar(value="0")
        self.anomalia3_raio = tk.StringVar(value="0.0")
        self.anomalia3_rotacao = tk.StringVar(value="0")
        self.x3_ptoCentral = tk.StringVar(value="0.0")
        self.y3_ptoCentral = tk.StringVar(value="0.0")
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
        opcao_var = tk.StringVar(value="Laplace")

        
        
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

        gmsh.model.add("circulo")

        pts_circulo = []
        for k in range(int(self.n_eletrodos.get())):
            x_circ = float(self.raio.get())*math.cos(2*math.pi*k/int(self.n_eletrodos.get()))
            y_circ = float(self.raio.get())*math.sin(2*math.pi*k/int(self.n_eletrodos.get()))
            pts_circulo.append(gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get())))

        ptoCentral = gmsh.model.geo.addPoint(0,0, 0, float(self.lc1.get())) 
        #gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='ptoCentral') 
        
        if int(self.NrAnomalias.get()) > 0:
            print('passei')
            pts_anomalia = []
            for anom in range(int(self.anomalia1_lados.get())):
                x_anomalia = float(self.x1_ptoCentral.get()) + float(self.anomalia1_raio.get())*math.cos((2*np.pi*anom/int(self.anomalia1_lados.get()))+ math.radians(float(self.anomalia1_rotacao.get())))
                y_anomalia = float(self.y1_ptoCentral.get()) + float(self.anomalia1_raio.get())*math.sin((2*np.pi*anom/int(self.anomalia1_lados.get()))+ math.radians(float(self.anomalia1_rotacao.get())))
                pts_anomalia.append(gmsh.model.geo.addPoint(x_anomalia, y_anomalia, 0.0,  float(self.lc2.get())))
        
        
        
            pts_anomalia2 = []
            for anom2 in range(int(self.anomalia2_lados.get())):
                x_anomalia2 = float(self.x2_ptoCentral.get()) + float(self.anomalia2_raio.get())*math.cos((2*np.pi*anom2/int(self.anomalia2_lados.get()))+ math.radians(float(self.anomalia2_rotacao.get())))
                y_anomalia2 = float(self.y2_ptoCentral.get()) + float(self.anomalia2_raio.get())*math.sin((2*np.pi*anom2/int(self.anomalia2_lados.get()))+ math.radians(float(self.anomalia2_rotacao.get())))
                pts_anomalia2.append(gmsh.model.geo.addPoint(x_anomalia2, y_anomalia2, 0.0,  float(self.lc3.get())))
            
  
     
            pts_anomalia3 = []
            for anom3 in range(int(self.anomalia3_lados.get())):
                x_anomalia3 = float(self.x3_ptoCentral.get()) + float(self.anomalia3_raio.get())*math.cos((2*np.pi*anom3/int(self.anomalia3_lados.get()))+ math.radians(float(self.anomalia3_rotacao.get())))
                y_anomalia3 = float(self.y3_ptoCentral.get()) + float(self.anomalia3_raio.get())*math.sin((2*np.pi*anom3/int(self.anomalia3_lados.get()))+ math.radians(float(self.anomalia3_rotacao.get())))
                pts_anomalia3.append(gmsh.model.geo.addPoint(x_anomalia3, y_anomalia3, 0.0,  float(self.lc4.get())))        
  
    
    
        linhas_circulo = []
        for L in range(int(self.n_eletrodos.get())):
            L_circulo = gmsh.model.geo.addCircleArc(pts_circulo[L], ptoCentral, pts_circulo[(L+1) % int(self.n_eletrodos.get())])
            linhas_circulo.append(L_circulo)
 
    
        loop_circulo = gmsh.model.geo.addCurveLoop(linhas_circulo)
    
        if int(self.NrAnomalias.get()) > 0:
            
            linhas_anomalia1 = []
            for A in range(int(self.anomalia1_lados.get())):
                L_anomalia1 = gmsh.model.geo.addLine(pts_anomalia[A], pts_anomalia[(A+1) % int(self.anomalia1_lados.get())])
                linhas_anomalia1.append(L_anomalia1)
                #print('pts_anomalia',pts_anomalia[A])
            loop_anomalia1 = gmsh.model.geo.addCurveLoop(linhas_anomalia1)
            
            linhas_anomalia2 = []
            for B in range(int(self.anomalia2_lados.get())):
                L_anomalia2 = gmsh.model.geo.addLine(pts_anomalia2[B], pts_anomalia2[(B+1) % int(self.anomalia2_lados.get())])
                linhas_anomalia2.append(L_anomalia2)
                #print('pts_anomalia',pts_anomalia[B])
            loop_anomalia2 = gmsh.model.geo.addCurveLoop(linhas_anomalia2)
            # cuba COM FURO (interface conformal)
            
                        
            linhas_anomalia3 = []
            for C in range(int(self.anomalia3_lados.get())):
                L_anomalia3 = gmsh.model.geo.addLine(pts_anomalia3[C], pts_anomalia3[(C+1) % int(self.anomalia3_lados.get())])
                linhas_anomalia3.append(L_anomalia3)
                #print('pts_anomalia',pts_anomalia[C])
            loop_anomalia3 = gmsh.model.geo.addCurveLoop(linhas_anomalia3)
            # cuba COM FURO (interface conformal)
            
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo, loop_anomalia1, loop_anomalia2, loop_anomalia3])

            # anomalia preenchendo o buraco
            surface_anomalia1 = gmsh.model.geo.addPlaneSurface([loop_anomalia1])
            surface_anomalia2 = gmsh.model.geo.addPlaneSurface([loop_anomalia2])
            surface_anomalia3 = gmsh.model.geo.addPlaneSurface([loop_anomalia3])



        else:
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo])


        # colocando ptoCentral na superfície da cuba
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        # criando physicalGroup da superfície
        physycal_circulo = gmsh.model.addPhysicalGroup(2, [surface_circulo], 1000)

        
        if int(self.NrAnomalias.get()) > 0 and surface_anomalia1 is not None:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "anomalia1")
            if int(self.NrAnomalias.get()) == 2 or int(self.NrAnomalias.get()) == 3:
                physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
                gmsh.model.setPhysicalName(2, physycal_anomalia2, "anomalia2")
            if int(self.NrAnomalias.get()) == 3:
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
                
        for idx in range(int(self.n_eletrodos.get())):
            gmsh.model.addPhysicalGroup(0, [pts_circulo[idx]], 10000+idx+1, name=f'electrode_{idx+1}')




        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_circulo)
        gmsh.model.mesh.generate(2)
        #gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
        

        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()

        gmsh.finalize()
    
    def save_file(self):
        #nome = self.nome_arquivo.get().strip()
        '''
        if not nome:
            messagebox.showerror("Erro", "Digite um nome de arquivo .")
            return

        out_dir = Path(self.out_dir.get().strip() or (Path.cwd() / "MalhasMSH"))
        try:
            os.makedirs(out_dir, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Erro", f"Não consegui criar a pasta:\n{out_dir}\n\n{e}")
            return

        path_out = out_dir / f"{nome}.msh"
        try:
            path_out.write_text(self.build_content(), encoding="utf-8")
        except Exception as e:
            messagebox.showerror("Erro", f"Não consegui salvar o arquivo:\n{path_out}\n\n{e}")
            return

        messagebox.showinfo("OK", f"Arquivo salvo:\n{path_out}")
        '''
        # atualiza dica do limite raio/2
        #raio = safe_float(self.raio.get(), 0.1)
        #self.info_lbl.config(text=f"Dica: no seu script, o limite do ponto central era 0 até raio/2 = {raio/2:g}")
        #print('banana')
        
        #content = self.build_content()
        #self.txt.delete("1.0", "end")
        #self.txt.insert("1.0", content)
        gmsh.initialize()

        gmsh.model.add("circulo")

        pts_circulo = []
        for k in range(int(self.n_eletrodos.get())):
            x_circ = float(self.raio.get())*math.cos(2*math.pi*k/int(self.n_eletrodos.get()))
            y_circ = float(self.raio.get())*math.sin(2*math.pi*k/int(self.n_eletrodos.get()))
            pts_circulo.append(gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, float(self.lc1.get())))

        ptoCentral = gmsh.model.geo.addPoint(0,0, 0, float(self.lc1.get())) 
        #gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='ptoCentral') 
        
        if int(self.NrAnomalias.get()) > 0:
            print('passei')
            pts_anomalia = []
            for anom in range(int(self.anomalia1_lados.get())):
                x_anomalia = float(self.x1_ptoCentral.get()) + float(self.anomalia1_raio.get())*math.cos((2*np.pi*anom/int(self.anomalia1_lados.get()))+ math.radians(float(self.anomalia1_rotacao.get())))
                y_anomalia = float(self.y1_ptoCentral.get()) + float(self.anomalia1_raio.get())*math.sin((2*np.pi*anom/int(self.anomalia1_lados.get()))+ math.radians(float(self.anomalia1_rotacao.get())))
                pts_anomalia.append(gmsh.model.geo.addPoint(x_anomalia, y_anomalia, 0.0,  float(self.lc2.get())))
        
        
        
            pts_anomalia2 = []
            for anom2 in range(int(self.anomalia2_lados.get())):
                x_anomalia2 = float(self.x2_ptoCentral.get()) + float(self.anomalia2_raio.get())*math.cos((2*np.pi*anom2/int(self.anomalia2_lados.get()))+ math.radians(float(self.anomalia2_rotacao.get())))
                y_anomalia2 = float(self.y2_ptoCentral.get()) + float(self.anomalia2_raio.get())*math.sin((2*np.pi*anom2/int(self.anomalia2_lados.get()))+ math.radians(float(self.anomalia2_rotacao.get())))
                pts_anomalia2.append(gmsh.model.geo.addPoint(x_anomalia2, y_anomalia2, 0.0,  float(self.lc3.get())))
            
  
     
            pts_anomalia3 = []
            for anom3 in range(int(self.anomalia3_lados.get())):
                x_anomalia3 = float(self.x3_ptoCentral.get()) + float(self.anomalia3_raio.get())*math.cos((2*np.pi*anom3/int(self.anomalia3_lados.get()))+ math.radians(float(self.anomalia3_rotacao.get())))
                y_anomalia3 = float(self.y3_ptoCentral.get()) + float(self.anomalia3_raio.get())*math.sin((2*np.pi*anom3/int(self.anomalia3_lados.get()))+ math.radians(float(self.anomalia3_rotacao.get())))
                pts_anomalia3.append(gmsh.model.geo.addPoint(x_anomalia3, y_anomalia3, 0.0,  float(self.lc4.get())))        
  
    
    
        linhas_circulo = []
        for L in range(int(self.n_eletrodos.get())):
            L_circulo = gmsh.model.geo.addCircleArc(pts_circulo[L], ptoCentral, pts_circulo[(L+1) % int(self.n_eletrodos.get())])
            linhas_circulo.append(L_circulo)
 
    
        loop_circulo = gmsh.model.geo.addCurveLoop(linhas_circulo)
    
        if int(self.NrAnomalias.get()) > 0:
            
            linhas_anomalia1 = []
            for A in range(int(self.anomalia1_lados.get())):
                L_anomalia1 = gmsh.model.geo.addLine(pts_anomalia[A], pts_anomalia[(A+1) % int(self.anomalia1_lados.get())])
                linhas_anomalia1.append(L_anomalia1)
                #print('pts_anomalia',pts_anomalia[A])
            loop_anomalia1 = gmsh.model.geo.addCurveLoop(linhas_anomalia1)
            
            linhas_anomalia2 = []
            for B in range(int(self.anomalia2_lados.get())):
                L_anomalia2 = gmsh.model.geo.addLine(pts_anomalia2[B], pts_anomalia2[(B+1) % int(self.anomalia2_lados.get())])
                linhas_anomalia2.append(L_anomalia2)
                #print('pts_anomalia',pts_anomalia[B])
            loop_anomalia2 = gmsh.model.geo.addCurveLoop(linhas_anomalia2)
            # cuba COM FURO (interface conformal)
            
                        
            linhas_anomalia3 = []
            for C in range(int(self.anomalia3_lados.get())):
                L_anomalia3 = gmsh.model.geo.addLine(pts_anomalia3[C], pts_anomalia3[(C+1) % int(self.anomalia3_lados.get())])
                linhas_anomalia3.append(L_anomalia3)
                #print('pts_anomalia',pts_anomalia[C])
            loop_anomalia3 = gmsh.model.geo.addCurveLoop(linhas_anomalia3)
            # cuba COM FURO (interface conformal)
            
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo, loop_anomalia1, loop_anomalia2, loop_anomalia3])

            # anomalia preenchendo o buraco
            surface_anomalia1 = gmsh.model.geo.addPlaneSurface([loop_anomalia1])
            surface_anomalia2 = gmsh.model.geo.addPlaneSurface([loop_anomalia2])
            surface_anomalia3 = gmsh.model.geo.addPlaneSurface([loop_anomalia3])


        #for L in range(int(self.anomalia2_lados.get())):
        #    L_circulo = gmsh.model.geo.addLine(pts_anomalia2[L], pts_anomalia2[(L+1) % int(self.anomalia2_lados.get())])
        #    linhas_circulo.append(L_circulo)

        #curva_circulo = gmsh.model.geo.addCurveLoop(linhas_circulo)
        else:
            surface_circulo = gmsh.model.geo.addPlaneSurface([loop_circulo])

        #linhas_anomalia1 = []
        #for A in range(int(self.anomalia1_lados.get())):
        #    L_anomalia1 = gmsh.model.geo.addLine(pts_anomalia[A], pts_anomalia[(A+1) % int(self.anomalia1_lados.get())])
        #    linhas_anomalia1.append(L_anomalia1)
    

        # colocando ptoCentral na superfície da cuba
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()
        #gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
        
        # criando physicalGroup da superfície
        physycal_circulo = gmsh.model.addPhysicalGroup(2, [surface_circulo], 1000)
        #gmsh.model.setPhysicalName(2, physycal_circulo, "cuba")
        
        if int(self.NrAnomalias.get()) > 0 and surface_anomalia1 is not None:
            physycal_anomalia1 = gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001)
            gmsh.model.setPhysicalName(2, physycal_anomalia1, "anomalia1")
            if int(self.NrAnomalias.get()) == 2 or int(self.NrAnomalias.get()) == 3:
                physycal_anomalia2 = gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002)
                gmsh.model.setPhysicalName(2, physycal_anomalia2, "anomalia2")
            if int(self.NrAnomalias.get()) == 3:
                physycal_anomalia3 = gmsh.model.addPhysicalGroup(2, [surface_anomalia3], 1003)
                gmsh.model.setPhysicalName(2, physycal_anomalia3, "anomalia3")
        
        #gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_anomalia1)
        #gmsh.model.addPhysicalGroup(0, [ptoCentral], 10000, name='central_node')
        
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
                
        for idx in range(int(self.n_eletrodos.get())):
            gmsh.model.addPhysicalGroup(0, [pts_circulo[idx]], 10000+idx+1, name=f'electrode_{idx+1}')
        #gmsh.model.addPhysicalGroup(2, [surface_circulo], 1000, name='body')
        #gmsh.model.addPhysicalGroup(2, [surface_anomalia1], 1001, name='anomalia1')
        #gmsh.model.addPhysicalGroup(2, [surface_anomalia2], 1002, name='anomalia1')



        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_circulo)
        gmsh.model.mesh.generate(2)
        #gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
        
        caminho = os.path.join(self.out_dir.get(), self.nome_arquivo.get() + ".msh")
        gmsh.write(caminho)
        #gmsh.write(self.out_dir.get()  + self.nome_arquivo.get()  + '.msh')
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

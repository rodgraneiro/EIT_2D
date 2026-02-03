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
        self.title("Gerador de Config (.py) - EIT")
        self.geometry("860x620")

        # ====== variáveis ======
        self.nome_arquivo = tk.StringVar(value="teste01")
        self.raio = tk.StringVar(value="0.15")
        self.n_eletrodos = tk.StringVar(value="32")
        self.lc1 = tk.StringVar(value="1e-1")
        
        
        self.NrAnomalias = tk.StringVar(value="0") 
        self.anomalia1_raio = tk.StringVar(value="0.01")
        self.anomalia1_lados = tk.StringVar(value="16")
        self.anomalia1_rotacao = tk.StringVar(value="22.5")
        self.x1_ptoCentral = tk.StringVar(value="0.01")
        self.y1_ptoCentral = tk.StringVar(value="0.01")
        self.lc2 = tk.StringVar(value="1e-1")

        self.anomalia2_lados = tk.StringVar(value="0")
        self.anomalia2_raio = tk.StringVar(value="0")
        self.anomalia2_rotacao = tk.StringVar(value="0")
        self.x2_ptoCentral = tk.StringVar(value="0")
        self.y2_ptoCentral = tk.StringVar(value="0")
        self.lc3 = tk.StringVar(value="0")

        self.out_dir = tk.StringVar(value=str(Path.cwd() / "Configs"))

        # ====== layout ======
        main = ttk.Frame(self, padding=12)
        main.pack(fill="both", expand=True)

        form = ttk.LabelFrame(main, text="Parâmetros", padding=10)
        form.grid(row=0, column=0, sticky="nsew")

        #preview = ttk.LabelFrame(main, text="Prévia do arquivo", padding=10)
        #preview.grid(row=0, column=1, sticky="nsew", padx=(12, 0))

        main.columnconfigure(0, weight=0)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        # ====== campos ======
        r = 0
        r = self._row(form, r, "Nome do arquivo (sem .py):", self.nome_arquivo)
        r = self._row(form, r, "Raio:", self.raio)
        r = self._row(form, r, "Nº eletrodos:", self.n_eletrodos)
        r = self._row(form, r, "lc1 (cuba):", self.lc1)
        ttk.Separator(form, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8
        )
        r += 1

        r = self._row(form, r, "Anomalia1 - nº de lados:", self.anomalia1_lados)
        r = self._row(form, r, "Raio:", self.anomalia1_raio)
        r = self._row(form, r, "lc2 (anomalia1):", self.lc2)
        r = self._row(form, r, "Anomalia1 - rotação (graus):", self.anomalia1_rotacao)
        r = self._row(form, r, "x_ptoCentral1:", self.x1_ptoCentral)
        r = self._row(form, r, "y_ptoCentral1:", self.y1_ptoCentral)
        ttk.Separator(form, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8
        )
        r += 1
        r = self._row(form, r, "Anomalia2 - nº de lados:", self.anomalia2_lados)
        r = self._row(form, r, "Raio:", self.anomalia2_raio)
        r = self._row(form, r, "lc3 (anomalia2):", self.lc3)
        r = self._row(form, r, "Anomalia2 - rotação (graus):", self.anomalia2_rotacao)
        r = self._row(form, r, "x_ptoCentral2:", self.x2_ptoCentral)
        r = self._row(form, r, "y_ptoCentral2:", self.y2_ptoCentral)
        
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
        #ttk.Button(btns, text="Salvar .py", command=self.save_file).pack(side="left", padx=(10, 0))

        # dica do limite (igual seu script: max = raio/2)
        #self.info_lbl = ttk.Label(form, text="", foreground="#444")
        #self.info_lbl.grid(row=r + 2, column=0, columnspan=2, sticky="w", pady=(10, 0))

        # ====== prévia ======
        #self.txt = tk.Text(preview, wrap="none")
        #self.txt.pack(fill="both", expand=True)
        #self.txt.configure(font=("Consolas", 10))
        '''
        # eventos: atualizar prévia ao digitar
        for var in [
            self.nome_arquivo, self.raio, self.n_eletrodos, self.lc1, self.lc2,
            self.anomalia_lados, self.anomalia_rotacao, self.x_ptoCentral, self.y_ptoCentral, self.out_dir
        ]:
            var.trace_add("write", lambda *_: self.update_preview())

        self.update_preview()
        '''
    def _row(self, parent, r, label, var):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=3)
        e = ttk.Entry(parent, textvariable=var, width=22)
        e.grid(row=r, column=1, sticky="w", padx=(8, 0), pady=3)
        return r + 1

    def choose_dir(self):
        d = filedialog.askdirectory(title="Escolher pasta de saída")
        if d:
            self.out_dir.set(d)
    '''
    def build_content(self) -> str:
        # strings/valores (mantendo o estilo do seu gerador)
        nome = self.nome_arquivo.get().strip() or "sem_nome"
        raio = safe_float(self.raio.get(), 0.1)
        ne = safe_int(self.n_eletrodos.get(), 16)
        lc1 = safe_float(self.lc1.get(), 1e-1)
        lc2 = safe_float(self.lc2.get(), 1e-1)
        lados = safe_int(self.anomalia_lados.get(), 4)
        rot = safe_float(self.anomalia_rotacao.get(), 0.0)
        x0 = safe_float(self.x_ptoCentral.get(), 0.0)
        y0 = safe_float(self.y_ptoCentral.get(), 0.0)

        # no seu script ponto_medio_x/y eram "0"
        ponto_medio_x = "0"
        ponto_medio_y = "0"
        
        # monta o arquivo .py
        lines = []
        lines.append(f'nome_arquivo = "{nome}"\n')
        lines.append(f"raio = {raio}\n")
        lines.append(f"n_eletrodos = {ne}\n")
        lines.append(f"ponto_medio_x = {ponto_medio_x}\n")
        lines.append(f"ponto_medio_y = {ponto_medio_y}\n")
        lines.append(f"lc1 = {lc1}\n")
        lines.append(f"lc2 = {lc2}\n")
        lines.append(f"anomalia_lados = {lados}\n")
        lines.append(f"anomalia_rotacao = {rot}\n")
        lines.append(f"x_ptoCentral = {x0}\n")
        lines.append(f"y_ptoCentral = {y0}\n")
        return "".join(lines)
    '''    
    def update_preview(self):
        # atualiza dica do limite raio/2
        #raio = safe_float(self.raio.get(), 0.1)
        #self.info_lbl.config(text=f"Dica: no seu script, o limite do ponto central era 0 até raio/2 = {raio/2:g}")
        print('banana')
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
        
        gmsh.model.geo.synchronize()
        #gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_circulo)
        gmsh.model.mesh.generate(2)
        #gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
        '''
        gmsh.write(self.nome_arquivo  + '.msh')
        gmsh.write(self.nome_arquivo  + ".geo_unrolled")
        os.rename(self.nome_arquivo  + ".geo_unrolled", self.nome_arquivo  + ".geo_unrolled")
        shutil.move(self.nome_arquivo  + ".geo_unrolled",self.nome_arquivo  + ".geo")
        
        gmsh.write(os.path.join("malhasMSH", nome_arquivo + '.msh'))
        gmsh.write(os.path.join("malhasGEO", nome_arquivo + ".geo_unrolled"))

        shutil.move(os.path.join("malhasGEO", nome_arquivo + ".geo_unrolled"), 
                    os.path.join("malhasGEO", nome_arquivo + ".geo"))
        # Launch the GUI to see the results:
        '''
        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()

        gmsh.finalize()
    '''
    def save_file(self):
        nome = self.nome_arquivo.get().strip()
        if not nome:
            messagebox.showerror("Erro", "Digite um nome de arquivo (sem .py).")
            return

        out_dir = Path(self.out_dir.get().strip() or (Path.cwd() / "Configs"))
        try:
            os.makedirs(out_dir, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Erro", f"Não consegui criar a pasta:\n{out_dir}\n\n{e}")
            return

        path_out = out_dir / f"{nome}.py"
        try:
            path_out.write_text(self.build_content(), encoding="utf-8")
        except Exception as e:
            messagebox.showerror("Erro", f"Não consegui salvar o arquivo:\n{path_out}\n\n{e}")
            return

        messagebox.showinfo("OK", f"Arquivo salvo:\n{path_out}")
    '''

if __name__ == "__main__":
    # no Windows, Tkinter já vem junto com Python padrão
    app = ConfigDashboard()
    app.mainloop()

import math
import numpy as np
import gmsh
import sys
import time
import os
import shutil
import importlib

nome_modulo = input("Digite o nome do arquivo de configuração (sem .py): ")
# adiciona a pasta "Configs" no caminho de importação
sys.path.append("Configs")
arquivoConfig = importlib.import_module(nome_modulo)


inicio = time.time()


# Dados iniciais
nome_arquivo = arquivoConfig.nome_arquivo
raio = arquivoConfig.raio
n_eletrodos = arquivoConfig.n_eletrodos
ponto_medio_x = arquivoConfig.ponto_medio_x
ponto_medio_y = arquivoConfig.ponto_medio_y
lc1 = arquivoConfig.lc1



x_arr = np.zeros([n_eletrodos,1])
y_arr = np.zeros([n_eletrodos,1])
z_1 = 0

for k in range(0,n_eletrodos):
    x_arr[k] = raio*math.cos(2*np.pi*k/n_eletrodos)
    y_arr[k] = raio*math.sin(2*np.pi*k/n_eletrodos)
    #print(x_arr[k], y_arr[k])




gmsh.initialize()

gmsh.model.add("poligono")



# Criar pontos z = 0
for j in range(0,n_eletrodos):
    gmsh.model.geo.addPoint(x_arr[j],y_arr[j], z_1, lc1, (j+1))
    #print('tag J= ', j)
ptoCentral = gmsh.model.geo.addPoint(0,0, 0, lc1)    
for L in range(1,n_eletrodos):
    gmsh.model.geo.addLine(L, (L+1), L)
    #print('tag= L', L)    

gmsh.model.geo.addLine(n_eletrodos, 1, n_eletrodos)

# tag das linhas
tag_linhas = []
for i in range(0,n_eletrodos):
    tag_linhas.append(i+1)


print('tag_linhas',tag_linhas)
Loop = gmsh.model.geo.addCurveLoop(tag_linhas)
base_superficie = gmsh.model.geo.addPlaneSurface([Loop])




gmsh.model.geo.synchronize()
gmsh.model.mesh.embed(0, [ptoCentral], 2, base_superficie)

base_Physical = gmsh.model.addPhysicalGroup(2, [base_superficie])



gmsh.model.mesh.generate(2)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   


gmsh.write(os.path.join(".\\malhasMSH\\", nome_arquivo + '.msh'))
gmsh.write(os.path.join(".\\malhasGEO\\", nome_arquivo + ".geo_unrolled"))

shutil.move(os.path.join(".\\malhasGEO\\", nome_arquivo + ".geo_unrolled"), 
            os.path.join(".\\malhasGEO\\", nome_arquivo + ".geo"))

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()


fim = time.time()
print(fim - inicio)
import math
import numpy as np
import gmsh
import sys
import time
import os
import shutil
import importlib

nome_modulo = input("Digite o nome do arquivo de configuração (sem .py): ")

arquivoConfig = importlib.import_module('Configs.'+ nome_modulo )


inicio = time.time()


# Dados iniciais
# Dados iniciais
nome_arquivo = arquivoConfig.nome_arquivo
raio = arquivoConfig.raio
n_eletrodos = arquivoConfig.n_eletrodos
ponto_medio_x = arquivoConfig.ponto_medio_x
ponto_medio_y = arquivoConfig.ponto_medio_y
lc1 = arquivoConfig.lc1
lc2 = arquivoConfig.lc2
anomalia_lados = arquivoConfig.anomalia_lados
anomalia_rotacao = arquivoConfig.anomalia_rotacao
x_ptoCentral = arquivoConfig.x_ptoCentral
y_ptoCentral = arquivoConfig.y_ptoCentral
anomalia_raio = raio/4


gmsh.initialize()

gmsh.model.add("circulo")

pts_circulo = []
for k in range(n_eletrodos):
    x_circ = raio*math.cos(2*math.pi*k/n_eletrodos)
    y_circ = raio*math.sin(2*math.pi*k/n_eletrodos)
    pts_circulo.append(gmsh.model.geo.addPoint(x_circ, y_circ, 0.0, lc1))

ptoCentral = gmsh.model.geo.addPoint(0,0, 0, lc1) 
    
pts_anomalia = []
for anom in range(anomalia_lados):
    x_anomalia = x_ptoCentral + anomalia_raio*math.cos((2*np.pi*anom/anomalia_lados)+ math.radians(anomalia_rotacao))
    y_anomalia = y_ptoCentral + anomalia_raio*math.sin((2*np.pi*anom/anomalia_lados)+ math.radians(anomalia_rotacao))
    pts_anomalia.append(gmsh.model.geo.addPoint(x_anomalia, y_anomalia, 0.0, lc2))



print('pts_circulo',pts_circulo )
print('pts_anomalia',pts_anomalia )


linhas_circulo = []
for L in range(n_eletrodos):
    L_circulo = gmsh.model.geo.addCircleArc(pts_circulo[L], ptoCentral, pts_circulo[(L+1) % n_eletrodos])
    linhas_circulo.append(L_circulo)

for L in range(anomalia_lados):
    L_circulo = gmsh.model.geo.addLine(pts_anomalia[L], pts_anomalia[(L+1) % anomalia_lados])
    linhas_circulo.append(L_circulo)

print('linhas_circulo',linhas_circulo )


curva_circulo = gmsh.model.geo.addCurveLoop(linhas_circulo)
surface_circulo = gmsh.model.geo.addPlaneSurface([curva_circulo])


physycal_circulo = gmsh.model.addPhysicalGroup(2, [surface_circulo])



linhas_anomalia = []
for A in range(anomalia_lados):
    L_anomalia = gmsh.model.geo.addLine(pts_anomalia[A], pts_anomalia[(A+1) % anomalia_lados])
    linhas_anomalia.append(L_anomalia)
    

curva_anomalia = gmsh.model.geo.addCurveLoop(linhas_anomalia)
surface_anomalia = gmsh.model.geo.addPlaneSurface([curva_anomalia])


physycal_anomalia = gmsh.model.addPhysicalGroup(2, [surface_anomalia])

#gmsh.model.mesh.generate(2)
#gmsh.option.setNumber("Mesh.SaveAll", 0)
#gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   





gmsh.model.geo.synchronize()
gmsh.model.mesh.embed(0, [ptoCentral], 2, surface_circulo)
gmsh.model.mesh.generate(2)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   

#gmsh.write(nome_arquivo  + '.msh')
#gmsh.write(nome_arquivo  + ".geo_unrolled")
#os.rename(nome_arquivo  + ".geo_unrolled", nome_arquivo  + ".geo_unrolled")
#shutil.move(nome_arquivo  + ".geo_unrolled",nome_arquivo  + ".geo")

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
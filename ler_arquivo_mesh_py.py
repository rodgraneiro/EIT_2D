#LER ARQUIVO MSH

import math
import numpy as np
import gmsh
import sys
import time

inicio = time.time()

#import gmsh

# Inicialize o Gmsh
gmsh.initialize()

# Carregue o arquivo .geo
gmsh.open('poligono_2D_Rho.msh')

# Obtenha as entidades do modelo
entities = gmsh.model.getEntities()
gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
#gmsh.write("poligono_2D_Rho_py.msh")
#gmsh.write("poligono_2D_Rho.geo_unrolled")



#gmsh.model.mesh.generate(3)
# Feche o arquivo .geo

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()

fim = time.time()
print(fim - inicio)
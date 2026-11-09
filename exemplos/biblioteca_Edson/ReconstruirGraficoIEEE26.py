# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:25:00 2026

@author: rodgr
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import meshio
from pathlib import Path

# ============================================================
# 1. ARQUIVOS
# ============================================================
'''
arquivo_sigma = (
    "sigma_reconstruidoIEEEstd027_300itr_"
    "1.2742749900e-03.npy"
)
'''
'''
sg_recIEEEstd027_300itr_Tol03lixo_1.0000000000e-06
sg_recIEEEstd027_300itr_Tol03lixo_1.0000000000e-06.npy

sg_recIEEEstd027_300itr_Tol03lixo_1.1288378917e-01.npy

sg_recIEEEstd027_300itr_Tol03lixo_1.2742749857e-02.npy

sg_recIEEEstd027_300itr_Tol03lixo_1.4384498883e-03.npy

sg_recIEEEstd027_300itr_Tol03lixo_1.6237767392e-04.npy

sg_recIEEEstd027_300itr_Tol03lixo_1.8329807108e-05.npy

sg_recIEEEstd027_300itr_Tol03lixo_2.0691380811e-06.npy

sg_recIEEEstd027_300itr_Tol03lixo_2.3357214691e-01.npy

sg_recIEEEstd027_300itr_Tol03lixo_2.6366508987e-02.npy

sg_recIEEEstd027_300itr_Tol03lixo_2.9763514416e-03.npy

sg_recIEEEstd027_300itr_Tol03lixo_3.3598182863e-04.npy

sg_recIEEEstd027_300itr_Tol03lixo_3.7926901907e-05.npy

sg_recIEEEstd027_300itr_Tol03lixo_4.2813323987e-06.npy

sg_recIEEEstd027_300itr_Tol03lixo_5.4555947812e-02.npy

sg_recIEEEstd027_300itr_Tol03lixo_6.1584821107e-03.npy

sg_recIEEEstd027_300itr_Tol03lixo_6.9519279618e-04.npy

sg_recIEEEstd027_300itr_Tol03lixo_7.8475997035e-05.npy

sg_recIEEEstd027_300itr_Tol03lixo_8.8586679041e-06.npy
'''


arquivo_sigma = ("sg_recIEEEstd027_300itr_af1_Tol1e3_6.1584821107e-03.npy"
)


arquivo_malha = "Hua_domain16e_base_v2D.msh"  # Substituir

lambda_reg = 6.1584821107e-03

# ============================================================
# 2. CARREGAR CONDUTIVIDADE RECONSTRUÍDA
# ============================================================

sigma = np.load(arquivo_sigma, allow_pickle=False)

print("Shape original:", sigma.shape)
print("Tipo:", sigma.dtype)

sigma = np.asarray(sigma).ravel()

print("Número de valores:", len(sigma))
print("Mínimo:", np.min(sigma))
print("Máximo:", np.max(sigma))

# ============================================================
# 3. CARREGAR MALHA
# ============================================================

mesh = meshio.read(arquivo_malha)

x = mesh.points[:, 0]
y = mesh.points[:, 1]

# Preservar a ordem dos blocos de triângulos
triangulos = np.vstack([
    bloco.data
    for bloco in mesh.cells
    if bloco.type == "triangle"
])

segmentos = [
    bloco.data
    for bloco in mesh.cells
    if bloco.type == "line"
]

print("Número de triângulos:", len(triangulos))

# ============================================================
# 4. VERIFICAR COMPATIBILIDADE
# ============================================================

ntri = len(triangulos)

if len(sigma) != ntri:
    raise ValueError(
        f"Sigma possui {len(sigma)} valores, "
        f"mas a malha possui {ntri} triângulos. "
        "Verifique a malha e a ordem dos elementos."
    )

# ============================================================
# 5. RECONSTRUIR GRÁFICO
# ============================================================

fig, ax = plt.subplots(figsize=(6, 5))

triang = tri.Triangulation(x, y, triangulos)

tpc = ax.tripcolor(
    triang,
    facecolors=sigma,
    edgecolors="k",
    linewidths=0.1,
    cmap="rainbow",
    vmin=0.0,
    vmax=4.0
)

cbar = fig.colorbar(
    tpc,
    ax=ax
)

cbar.set_label(
    "σ (Conductivity) (S/m)",
    fontsize=15
)



'''
fig.colorbar(
    tpc,
    ax=ax,
    label="σ (Conductivity) (S/m)"
)
'''
'''
# Elementos 1D
for bloco in segmentos:
    for n1, n2 in bloco:
        ax.plot(
            [x[n1], x[n2]],
            [y[n1], y[n2]],
            color="red",
            linewidth=2
        )
'''
#ax.set_title(f"Conductivity Isotropic - Lambda_{lambda_reg:.2e}",fontsize=12)

expoente = int(np.floor(np.log10(abs(lambda_reg))))
mantissa = lambda_reg / 10**expoente

ax.set_title(
    rf"Conductivity Isotropic",# $\rightarrow$ $\lambda = {mantissa:.2f}\times 10^{{{expoente}}}$",
    fontsize=15
)

ax.set_xlabel("Length [m]", fontsize=15)
ax.set_ylabel("Length [m]", fontsize=15)

ax.set_aspect("equal")

plt.tight_layout()

# ============================================================
# 6. SALVAR E EXIBIR
# ============================================================

saida = Path(f"sg_recIEEE_Isotropic_6e03b_{lambda_reg:.2e}.png")

plt.savefig(
    saida,
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Figura salva em:", saida.resolve())

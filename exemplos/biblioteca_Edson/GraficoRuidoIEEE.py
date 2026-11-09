# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 21:21:04 2026

@author: rodgr
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# ============================================================
# 1. Dados experimentais
# ============================================================

# Primeira linha: caso sem ruido (0%)
# Demais linhas: niveis de ruido em fracao do RMS das tensoes
'''
noise = np.array([
    0.0,
    1e-7,
    1e-6,
    2e-6,
    3e-6,
    4e-6,
    5e-6,
    6e-6,
    7e-6,
    8e-6,
    9e-6,
    1e-5
])
'''

noise = np.array([
    0.0,
    0.00001,
    0.0001,
    0.0002,
    0.0003,
    0.0004,
    0.0005,
    0.0006,
    0.0007,
    0.0008,
    0.0009,
    0.001,
])
error_background = np.array([
    11.76, 11.60, 12.23, 14.29,
    16.85, 20.53, 35.49, 36.25,
    37.39, 39.99, 42.52, 44.39
])

error_obj1 = np.array([
    19.92, 19.97, 21.22, 22.74,
    24.02, 27.05, 47.72, 52.37,
    54.73, 62.59, 74.84, 81.07
])

error_obj2 = np.array([
    21.48, 21.15, 19.36, 19.68,
    21.04, 24.24, 55.86, 64.61,
    70.69, 75.12, 75.21, 77.35
])

error_obj3 = np.array([
    20.74, 19.27, 26.84, 36.88,
    42.83, 48.30, 57.57, 55.48,
    56.72, 58.36, 60.46, 62.28
])

error_global = np.array([
    13.05, 12.82, 13.79, 16.29,
    18.92, 22.59, 37.96, 39.19,
    40.68, 43.56, 46.42, 48.53
])

# ============================================================
# 2. Configuracao da figura
# ============================================================

fig, ax = plt.subplots(figsize=(7.5, 5.0))

ax.plot(noise, error_background, '-o',
        linewidth=1.6, markersize=4,
        label='Background')

ax.plot(noise, error_obj1, '-s',
        linewidth=1.6, markersize=4,
        label='Object 1')

ax.plot(noise, error_obj2, '-^',
        linewidth=1.6, markersize=4,
        label='Object 2')

ax.plot(noise, error_obj3, '-D',
        linewidth=1.6, markersize=4,
        label='Object 3')

ax.plot(noise, error_global, '-o',
        linewidth=2.2, markersize=5,
        label='All')

# ============================================================
# 3. Eixos e formatacao
# ============================================================

ax.set_xlabel('Measurement noise level (%)', fontsize=15)


#ax.set_xlabel('Measurement noise level (fraction of RMS)', fontsize=15)

ax.set_ylabel('Error (%)',
              fontsize=15)
'''
ax.set_xlim(0, 1e-5)
ax.set_ylim(0, 90)

ax.xaxis.set_major_locator(MultipleLocator(2e-6))
ax.yaxis.set_major_locator(MultipleLocator(10))
'''

ax.set_xlim(0, 0.001)
ax.set_ylim(0, 90)
ax.xaxis.set_major_locator(MultipleLocator(0.0002))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))

# Exibir o eixo X como 0, 2, 4, 6, 8, 10
# com multiplicador 10^-6 no rotulo do eixo

#ax.xaxis.set_major_formatter(
#    FuncFormatter(lambda x, pos: f'{x / 1e-6:g}')
#)



ax.tick_params(axis='both', labelsize=10)

ax.grid(True, linestyle='--', alpha=0.35)

ax.legend(
    loc='upper left',
    fontsize=12,
    frameon=True
)
titulo="Sensitivity to measurement noise"
ax.set_title(titulo, fontsize=20)
# ============================================================
# 4. Salvar figura para o artigo IEEE
# ============================================================

plt.tight_layout()

plt.savefig(
    'IEEE_noise_Frobenius30out26.png',
    dpi=500,
    bbox_inches='tight'
)

plt.savefig(
    'IEEE_noise_Frobenius30out26.pdf',
    bbox_inches='tight'
)

plt.show()

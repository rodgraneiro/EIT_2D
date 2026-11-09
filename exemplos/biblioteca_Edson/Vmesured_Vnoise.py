# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 21:11:58 2026

@author: rodgr
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1. Carregar as tensões originais
# ============================================================

V_measured_phaton = np.load("V_measured_phatom_anisoIEEE.npy")

# ============================================================
# 2. Adicionar ruído gaussiano
# ============================================================

pct = '0.0001%'

rng = np.random.default_rng(42)

nivel_ruido = 0.00001  # 1%

V_rms = np.sqrt(np.mean(V_measured_phaton**2))

ruido = rng.normal(
    loc=0.0,
    scale=nivel_ruido * V_rms,
    size=V_measured_phaton.shape
)

V_measured_com_ruido = V_measured_phaton + ruido

ruido_rms = np.sqrt(np.mean(ruido**2))

# ============================================================
# 3. Preparar os dados para visualização
# ============================================================

V_original = V_measured_phaton.flatten(order='F')
V_ruidoso = V_measured_com_ruido.flatten(order='F')

medidas = np.arange(1, len(V_original) + 1)

# ============================================================
# 4. Plotar as curvas sobrepostas semi-log
# ============================================================

fig, ax = plt.subplots(figsize=(12, 5))

ax.plot(
    medidas,
    V_original,
    color='blue',
    linewidth=1.5,
    label='Sem ruído'
)

ax.plot(
    medidas,
    V_ruidoso,
    color='red',
    linewidth=1.0,
    alpha=0.7,
    label='Com ruído (1%)'
)

ax.set_xlabel('Número da medição')
ax.set_ylabel('Tensão (V)')
ax.set_yscale('symlog', linthresh=1e-6)
ax.set_title(f'Tensões medidas - Comparação com e sem ruído - {pct}')

ax.grid(True, linestyle='--', alpha=0.4)
ax.legend()

plt.tight_layout()
plt.show()

# ============================================================
# 5. Informações estatísticas
# ============================================================

snr_db = 20 * np.log10(V_rms / ruido_rms)

print(f"V RMS original : {V_rms:.8e} V")
print(f"Ruído RMS      : {ruido_rms:.8e} V")
print(f"Ruído aplicado : {nivel_ruido*100:.2f}%")
print(f"SNR            : {snr_db:.2f} dB")

# ============================================================
# 6. Plotar as curvas sobrepostas
# ============================================================

fig, ax = plt.subplots(figsize=(12, 5))

ax.plot(
    medidas,
    V_original,
    color='blue',
    linewidth=1.5,
    label='Sem ruído'
)

ax.plot(
    medidas,
    V_ruidoso,
    color='red',
    linewidth=1.0,
    alpha=0.7,
    label='Com ruído (1%)'
)

ax.set_xlabel('Número da medição')
ax.set_ylabel('Tensão (V)')
#ax.set_yscale('symlog', linthresh=1e-6)
ax.set_title(f'Tensões medidas - Comparação com e sem ruído -{pct}')

ax.grid(True, linestyle='--', alpha=0.4)
ax.legend()

plt.tight_layout()
plt.show()


'''
rng = np.random.default_rng(42)
nivel_ruido = 0.01

V_rms = np.sqrt(np.mean(V_measured_phaton**2))

ruido = rng.normal(
    loc=0.0,
    scale=nivel_ruido * V_rms,
    size=V_measured_phaton.shape
)

V_measured_com_ruido = V_measured_phaton + ruido
'''
'''
# Primeira injeção: eletrodos 0 e 4
V_original = V_measured_phaton[:, 0]
V_ruidoso = V_measured_com_ruido[:, 0]

eletrodos = np.arange(16)

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(eletrodos, V_original, 'o-',
        color='blue', label='Sem ruído')

ax.plot(eletrodos, V_ruidoso, 's--',
        color='red', label='Com ruído (1%)')

ax.axvline(0, color='green', linestyle=':',
           label='Injeção (+I)')

ax.axvline(4, color='purple', linestyle=':',
           label='Retorno (-I)')

ax.set_yscale('symlog', linthresh=1e-6)

ax.set_xticks(eletrodos)
ax.set_xlabel('Eletrodo')
ax.set_ylabel('Potencial (V)')
ax.set_title('Potenciais dos eletrodos - Skip-3 (0,4)')

ax.grid(True, linestyle='--', alpha=0.4)
ax.legend()

plt.tight_layout()
plt.show()
'''
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1. Carregar tensões originais
# ============================================================

V_measured_phaton = np.load(
    "V_measured_phatom_anisoIEEE.npy"
)

assert V_measured_phaton.shape == (16, 16)

# ============================================================
# 2. Adicionar ruído gaussiano de 1%
# ============================================================

rng = np.random.default_rng(42)
'''
nivel_ruido = 0.02

V_rms = np.sqrt(np.mean(V_measured_phaton**2))

ruido = rng.normal(
    loc=0.0,
    scale=nivel_ruido * V_rms,
    size=V_measured_phaton.shape
)

V_measured_com_ruido = V_measured_phaton + ruido

ruido_rms = np.sqrt(np.mean(ruido**2))
'''
SNR_dB = 20 * np.log10(V_rms / ruido_rms)

print(f"V RMS global   = {V_rms:.6e} V")
print(f"Ruído RMS      = {ruido_rms:.6e} V")
print(f"Ruído (%)      = {100*ruido_rms/V_rms:.4f}%")
print(f"SNR global     = {SNR_dB:.2f} dB")

# ============================================================
# 3. Organizar as 256 medições por padrão skip-3
# ============================================================

V_original = V_measured_phaton.flatten(order='F')
V_ruidoso = V_measured_com_ruido.flatten(order='F')

ruido_plot = V_ruidoso - V_original

medidas = np.arange(1, 257)

# ============================================================
# 4. Gráfico exclusivo do ruído
# ============================================================

fig, ax = plt.subplots(figsize=(12, 5))

ax.plot(
    medidas,
    ruido_plot,
    color='red',
    linewidth=1.0,
    label=(f'Ruído gaussiano - {pct}')
)

ax.axhline(0, color='black', linewidth=0.8)

# Limites de ±1 desvio-padrão teórico
sigma_ruido = nivel_ruido * V_rms

ax.axhline(
    sigma_ruido,
    color='blue',
    linestyle='--',
    label=r'$+\sigma_\eta$'
)

ax.axhline(
    -sigma_ruido,
    color='blue',
    linestyle='--',
    label=r'$-\sigma_\eta$'
)

# Separação dos 16 padrões de injeção
for k in range(1, 16):
    ax.axvline(
        16*k + 0.5,
        color='gray',
        linestyle=':',
        alpha=0.4
    )

ax.set_xlabel('Número da medição')
ax.set_ylabel('Ruído (V)')
ax.set_title(f'Ruído gaussiano de {pct} - 256 medições')

ax.grid(True, linestyle='--', alpha=0.3)
ax.legend()

plt.tight_layout()
plt.show()

# ============================================================
# 5. Calcular SNR para cada padrão skip-3
# ============================================================

print("\nSNR por padrão de injeção skip-3")
print("-" * 55)

snr_por_padrao = []

for k in range(16):

    V_k = V_measured_phaton[:, k]
    ruido_k = ruido[:, k]

    Vrms_k = np.sqrt(np.mean(V_k**2))
    Nrms_k = np.sqrt(np.mean(ruido_k**2))

    snr_k = 20 * np.log10(Vrms_k / Nrms_k)

    snr_por_padrao.append(snr_k)

    inj = k
    ret = (k + 4) % 16

    print(
        f"Padrão {k:02d} ({inj:02d},{ret:02d}): "
        f"SNR = {snr_k:7.2f} dB"
    )

# ============================================================
# 6. Gráfico da SNR por padrão de injeção
# ============================================================

fig, ax = plt.subplots(figsize=(11, 5))

padroes = np.arange(16)

ax.plot(
    padroes,
    snr_por_padrao,
    'o-',
    color='blue',
    linewidth=1.5
)

ax.axhline(
    SNR_dB,
    color='red',
    linestyle='--',
    label=f'SNR global = {SNR_dB:.2f} dB - {pct}'
)

ax.set_xticks(padroes)

ax.set_xlabel('Padrão de injeção skip-3')
ax.set_ylabel('SNR (dB)')
ax.set_title(f'Relação sinal-ruído por padrão de injeção - {pct}')

ax.grid(True, linestyle='--', alpha=0.4)
ax.legend()

plt.tight_layout()
plt.show()
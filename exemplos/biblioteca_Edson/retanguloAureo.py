# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:08:48 2026

@author: rodgr
"""

import matplotlib.pyplot as plt

phi = (1 + 5**0.5) / 2

a = 8
b = a * phi
print('b=', b)

fig, ax = plt.subplots()
rect = plt.Rectangle((0, 0), a, b, fill=False)
ax.add_patch(rect)

ax.set_xlim(0, b)
ax.set_ylim(0, b)
ax.set_aspect('equal')

plt.title(f"Retângulo Áureo ({a} x {b:.2f})")
plt.show()
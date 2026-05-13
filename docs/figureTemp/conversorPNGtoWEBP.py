# -*- coding: utf-8 -*-
"""
Created on Wed May 13 11:48:26 2026

@author: rodgr
"""

from PIL import Image

img = Image.open("batata_sigma_xy_0.00657933225.png")

img.save(
    "batata_sigma_xy_0.00657933225.webp",
    format="WEBP",
    quality=85
)
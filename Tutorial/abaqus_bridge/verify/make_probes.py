#!/usr/bin/env python3
"""Generate probes.txt: (x, y, z, t, znorm) points exercising every
branch of the flux model.

Columns are read by both driver.f and compare_flux.py.
"""
import numpy as np

tp = np.loadtxt("toolpath.crs")


def ty(t):
    return float(np.interp(t, tp[:, 0], tp[:, 2]))


rows = []

# --- radial sweep across both flux branches, on the layer-1 tool plane.
# t=60 s is mid layer-1 deposition (state=1, tool z=3).
t = 60.0
y0 = ty(t)
for dx in (0.0, 2.0, 4.999, 5.001, 8.0, 14.999, 15.001, 20.0):
    rows.append((dx, y0, 3.0, t, 1.0))

# --- orientation gating (production mode only)
rows.append((2.0, y0, 3.0, t, -1.0))   # downward-facing face
rows.append((2.0, y0, 3.0, t, 0.0))    # vertical / side face

# --- the substrate-underside bleed: 8 mm below a tool at z=3, so the
# 3-D radius is sqrt(2^2+8^2)=8.25 < r_shoulder and the reference model
# deposits q2 here. Production mode must return exactly zero.
rows.append((2.0, y0, -5.0, t, 1.0))

# --- deposition state gating
rows.append((2.0, ty(150.0), 6.0, 150.0, 1.0))    # inter-layer dwell
rows.append((2.0, ty(900.0), 15.0, 900.0, 1.0))   # past last waypoint
rows.append((2.0, ty(0.05), 0.0, 0.05, 1.0))      # before first waypoint

# --- a later layer, to confirm the toolpath lookup tracks z
rows.append((2.0, ty(300.0), 6.0, 300.0, 1.0))

np.savetxt("probes.txt", np.array(rows), fmt="%.10f")
print(f"wrote probes.txt ({len(rows)} probes)")

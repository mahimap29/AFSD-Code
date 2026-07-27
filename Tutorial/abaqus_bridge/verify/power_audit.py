#!/usr/bin/env python3
"""Closed-form power audit for the AFSD surface heat source.

Two purposes:

1. Give the total tool power in closed form, so an Abaqus run can be
   checked against it.  Requesting the integrated surface flux over the
   DFLUX surface and comparing to this number is the cheapest possible
   test that SCALE is right -- a units error shows up as a factor of
   1000, not as a subtle discrepancy.

2. Quantify the substrate-underside bleed caused by the reference
   implementation's 3-D radius (gamma_afsd.py:738).  This is the
   justification for the in-plane radius used in production mode.

Integrals (P in MPa, lengths in mm, omega in rad/s, result in W after
the reference /1000):

  K       = delta*mu + (1-delta)/sqrt(3)          (shoulder coefficient)
  Q_pin   = (2*pi/1000) * [ mu*P*omega*r_p^3/3 + mu*P*v*r_p^2/2 ]
  Q_shld  = (2*pi/1000) * K*P*omega*(r_s^3 - r_p^3)/3
  Q_under = (2*pi/1000) * K*P*omega * [ r_s^3 - dz^3 ] / 3
            over the disc where sqrt(rho^2+dz^2) < r_s, i.e.
            rho < sqrt(r_s^2 - dz^2); substituting u = rho^2+dz^2 gives
            the bracket above directly.
"""
import math

import numpy as np


def load_params(path="afsd_params.inp"):
    with open(path) as f:
        rows = [ln.split() for ln in f if ln.strip()]
    P, rs, rp, w, v, mu, slip = (float(x) for x in rows[0][:7])
    return P, rs, rp, w, v, mu, slip


def main():
    P, rs, rp, w, v, mu, slip = load_params()
    K = slip * mu + (1.0 - slip) / math.sqrt(3.0)
    c = 2.0 * math.pi / 1000.0  # 2*pi and the reference /1000, -> W

    q_pin = c * (mu * P * w * rp ** 3 / 3.0 + mu * P * v * rp ** 2 / 2.0)
    q_shld = c * K * P * w * (rs ** 3 - rp ** 3) / 3.0
    q_top = q_pin + q_shld

    print(f"P={P} MPa  r_s={rs} mm  r_p={rp} mm  omega={w} rad/s  "
          f"v={v} mm/s")
    print(f"mu={mu}  delta={slip}  ->  K = delta*mu + (1-delta)/sqrt(3) "
          f"= {K:.7f}")
    print()
    print("Intended tool power (top surface):")
    print(f"  pin branch      {q_pin:10.3f} W")
    print(f"  shoulder branch {q_shld:10.3f} W")
    print(f"  TOTAL           {q_top:10.3f} W   "
          f"= {q_top * 1000:.1f} mW  <-- compare against Abaqus's")
    print(f"{'':34}integrated DFLUX surface flux")

    # --- the underside bleed, for the tutorial geometry --------------
    # substrate underside at z=-5, tool plane at z=3 during layer 1
    dz = 8.0
    if dz < rs:
        rho_max = math.sqrt(rs ** 2 - dz ** 2)
        q_under = c * K * P * w * (rs ** 3 - dz ** 3) / 3.0
        print()
        print(f"Reference 3-D-radius bleed onto a surface {dz:.0f} mm "
              f"below the tool plane:")
        print(f"  in-plane extent   rho < {rho_max:.2f} mm")
        print(f"  spurious power    {q_under:10.3f} W")
        print(f"  as a fraction of intended power  "
              f"{100 * q_under / q_top:.1f}%")

        # numerical cross-check of the closed form
        rho = np.linspace(0.0, rho_max, 200001)
        r3d = np.sqrt(rho ** 2 + dz ** 2)
        q = K * P * w * r3d / 1000.0
        num = np.trapezoid(q * 2.0 * math.pi * rho, rho)
        print(f"  numerical check   {num:10.3f} W  "
              f"(rel diff {abs(num - q_under) / q_under:.2e})")
        print()
        print("  This is why production mode uses the in-plane radius:")
        print("  a flat shoulder occupies one plane and cannot heat a")
        print("  face 8 mm below it. Note mu=0.2/delta=0.90 were fitted")
        print("  with this bleed present, so they are not transferable")
        print("  to a corrected model without recalibration.")


if __name__ == "__main__":
    main()

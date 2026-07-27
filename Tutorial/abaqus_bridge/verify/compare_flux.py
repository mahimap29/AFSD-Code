#!/usr/bin/env python3
"""Prove the Fortran AFSD_QVAL reproduces gamma_afsd.py's flux.

This is the Phase 1 gate: nothing downstream (mu calibration, ring case)
is trustworthy until this passes.

The reference expression is transcribed directly from
Tutorial/gamma_afsd.py:741-744, using numpy instead of cupy so it runs
without a GPU.  The Fortran side is exercised by driver.f, which loads
the same afsd_params.inp and toolpath.crs through UEXTERNALDB.

The two differ by the unit factor: gamma_afsd.py hardcodes /1000 to land
on W/mm^2 (mm-gram-s), while the Fortran applies SCALE from
afsd_params.inp -- 1.0 for a mm-tonne-s Abaqus model (mW/mm^2).  So the
comparison is

    q_fortran  ==  q_python * (SCALE / 1e-3)

Usage:
    python3 compare_flux.py            # after running ./driver
"""

import sys
import numpy as np

TOL = 1e-10


def load_params(path="afsd_params.inp"):
    with open(path) as f:
        rows = [ln.split() for ln in f if ln.strip()]
    p = [float(v) for v in rows[0]]
    pressure, r_shoulder, r_pin, w, v, mu, slip = p[:7]
    scale = float(rows[1][0])
    irmode = int(float(rows[1][1]))
    ztol = float(rows[1][2])
    return dict(pressure=pressure, r_shoulder=r_shoulder, r_pin=r_pin,
                w=w, v=v, mu=mu, slip=slip, scale=scale, irmode=irmode,
                ztol=ztol)


def load_toolpath(path="toolpath.crs"):
    return np.loadtxt(path)


def tool_at(tp, t):
    """Piecewise-linear tool centre, plus the BINARY deposition state.

    gamma_afsd.py:156-158 interpolates position with np.interp, and
    line 160 interpolates the state column -- but lines 161-167 then
    re-binarise it:

        if laser_state[i+1] > laser_state[i] or laser_state[i] == 1:
            l[i] = 1
        else:
            l[i] = 0

    which reduces to "the segment is ON iff the state at its END
    waypoint is 1".  Line 168 additionally forces it off past the final
    waypoint.  Reproduced here so the reference and the Fortran agree.
    """
    x = np.interp(t, tp[:, 0], tp[:, 1])
    y = np.interp(t, tp[:, 0], tp[:, 2])
    z = np.interp(t, tp[:, 0], tp[:, 3])

    if t <= tp[0, 0] or t >= tp[-1, 0]:
        return x, y, z, 0.0
    i = int(np.searchsorted(tp[:, 0], t, side="right")) - 1
    i = min(max(i, 0), len(tp) - 2)
    s = 1.0 if tp[i + 1, 4] >= 1 else 0.0
    return x, y, z, s


def flux_python(pt, prm, tp):
    """gamma_afsd.py:741-744, verbatim in structure.

    Returns W/mm^2, i.e. with the reference /1000 applied.
    """
    x, y, z, t, znorm = pt
    xt, yt, zt, state = tool_at(tp, t)

    gated = False
    if prm["irmode"] == 1:
        # equivalence mode: 3-D distance, no orientation/plane filtering,
        # matching gamma_afsd.py's surface_xy (any z=const face, either
        # facing direction)
        r = np.sqrt((x - xt) ** 2 + (y - yt) ** 2 + (z - zt) ** 2)
    else:
        # production mode: in-plane radius, upward faces at the tool
        # plane only
        r = np.sqrt((x - xt) ** 2 + (y - yt) ** 2)
        if znorm <= 0 or abs(z - zt) > prm["ztol"]:
            gated = True

    P, mu, w, v = prm["pressure"], prm["mu"], prm["w"], prm["v"]
    slip = prm["slip"]

    q1 = (mu * P * w * r + mu * P * v) / 1000
    q2 = (slip * mu * P * w * r) / 1000 \
        + ((1 - slip) * P * w * r / np.sqrt(3)) / 1000

    q = np.where(r < prm["r_pin"], q1,
                 np.where(r < prm["r_shoulder"], q2, 0.0))

    if gated:
        return 0.0, float(r)
    return float(q * state), float(r)


def main():
    prm = load_params()
    tp = load_toolpath()
    probes = np.atleast_2d(np.loadtxt("probes.txt"))
    fort = np.atleast_2d(np.loadtxt("flux_fortran.txt"))

    if len(fort) != len(probes):
        sys.exit(f"FAIL: {len(probes)} probes but {len(fort)} Fortran rows")

    unit = prm["scale"] / 1e-3  # W/mm^2 -> SCALE-consistent units

    print(f"SCALE={prm['scale']}  IRMODE={prm['irmode']}  "
          f"unit factor={unit:g}")
    print(f"{'r [mm]':>10} {'q_python':>16} {'q_fortran':>16} "
          f"{'abs err':>11} {'rel err':>10}")

    worst_rel = 0.0
    worst_abs = 0.0
    for probe, frow in zip(probes, fort):
        q_py, r_py = flux_python(probe, prm, tp)
        r_f, q_f = frow[4], frow[5]
        expect = q_py * unit

        aerr = abs(q_f - expect)
        rerr = aerr / abs(expect) if expect != 0 else aerr
        worst_abs = max(worst_abs, aerr)
        worst_rel = max(worst_rel, rerr)

        if abs(r_f - r_py) > 1e-8:
            print(f"  radius mismatch: fortran {r_f} vs python {r_py}")

        print(f"{r_f:10.4f} {expect:16.8e} {q_f:16.8e} "
              f"{aerr:11.3e} {rerr:10.3e}")

    print(f"\nworst abs err {worst_abs:.3e}   worst rel err {worst_rel:.3e}")

    # documented sanity value from the plan / subroutine header
    q1_ref = (0.2 * 5.0 * 41.888 * 5.0 + 0.2 * 5.0 * 1.67) / 1000
    print(f"reference q1 at r=r_pin=5mm: {q1_ref:.4f} W/mm^2 "
          f"= {q1_ref * 1000:.1f} mW/mm^2")

    if worst_rel > TOL:
        sys.exit(f"\nFAIL: worst relative error {worst_rel:.3e} > {TOL:g}")
    print(f"\nPASS: all probes agree within {TOL:g}")


if __name__ == "__main__":
    main()

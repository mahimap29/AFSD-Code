#------------------------------------------------------------------------------
# generate_toolpath_ring.py
#
# Multi-pass concentric circular toolpath generator, for AFSD deposition of a
# ring stiffener onto a plate (the Agarwal et al. AIAA 2026-1670 benchmark
# geometry: 500x500x20 mm SS304L substrate, 3-layer Ti-6Al-4V ring, inner
# radius 100 mm, outer radius 200 mm).
#
# Emits the same 5-column toolpath.crs format the rest of AFSD-Code uses
# (time[s] x[mm] y[mm] z[mm] state[0/1]), so one file drives both the Python
# solver and the Abaqus DFLUX bridge in abaqus_bridge/.
#
# This is a separate module rather than an addition to generate_toolpath.py
# because that file calls its generator at import time, so importing it would
# overwrite toolpath.crs as a side effect.
#
# Two constraints govern whether a configuration is buildable at all. Both
# are computed and reported by check_feasibility() rather than assumed away:
#
#   1. COVERAGE. The ring is 100 mm wide. A 30 mm AFSD tool cannot cover that
#      in one pass, so multiple concentric passes are required. The pass
#      radii here are inset by the shoulder radius, so coverage is exact.
#
#   2. MASS BALANCE.  h_layer * s_stepover * v_travel = A_rod * v_feed
#      The feed rate is not free: it is fixed by the layer thickness, the
#      step-over, the travel speed and the feedstock cross-section. A 30 mm
#      tool fed by a 10 mm round rod needs ~198 mm/min of feed to build a
#      6.67 mm layer at 100 mm/min travel, against ~26 mm/min for a 50 mm
#      tool fed by 1-inch square stock. Both are reported so the trade is
#      explicit.
#
#      Note the VOLUMETRIC rate (~16,700 mm^3/min ~ 1000 cm^3/h ~ 4.4 kg/h
#      of Ti-6Al-4V) is essentially independent of the tool/rod choice --
#      it is set by the ring volume and the schedule. That figure is at the
#      high end of AFSD practice and should be checked against the rated
#      capability of the specific machine. It is reported, not asserted
#      against a hard bound, because published AFSD deposition rates span
#      more than an order of magnitude between bench and large-format
#      systems.
#
# Note also that Agarwal's stated 12 s/layer implies a traverse of ~18,850
# mm/min, roughly 190x a realistic AFSD speed, and a deposition rate ~800x
# physical. Preserving it does not give the same answer faster -- with a
# physics-based flux it gives a different wrong answer (a near-cold part),
# because line energy collapses. This generator computes the schedule from
# the travel speed instead.
#
# Licensed under the MIT License; see LICENSE file for details
#------------------------------------------------------------------------------

import math

# Feed stock cross-sections [mm^2]
ROD_10MM_ROUND = math.pi * 5.0 ** 2      #  78.54, as used by Hu et al.
ROD_1IN_SQUARE = 25.4 ** 2               # 645.16, typical MELD feedstock


def pass_radii(r_inner, r_outer, r_shoulder, n_passes):
    """Centre radii of n_passes concentric passes covering [r_inner, r_outer].

    The outermost pass centres are inset by the shoulder radius so the tool
    footprint reaches, but does not badly overrun, both ring edges.
    """
    if n_passes < 1:
        raise ValueError("n_passes must be >= 1")
    c0 = r_inner + r_shoulder
    c1 = r_outer - r_shoulder
    if n_passes == 1:
        return [0.5 * (c0 + c1)]
    if c1 < c0:
        raise ValueError(
            f"tool too wide: shoulder r={r_shoulder} mm cannot fit inside a "
            f"{r_outer - r_inner:.1f} mm wide ring")
    step = (c1 - c0) / (n_passes - 1)
    return [c0 + i * step for i in range(n_passes)]


def check_feasibility(r_inner, r_outer, r_shoulder, n_passes, n_layers,
                      layer_thickness, travel_speed_mm_min, rod_area,
                      verbose=True):
    """Report coverage, overlap and the feed rate mass balance demands.

    Returns a dict; raises nothing, so a marginal configuration can still be
    generated deliberately. Read the warnings.
    """
    radii = pass_radii(r_inner, r_outer, r_shoulder, n_passes)
    step = (radii[1] - radii[0]) if n_passes > 1 else 2 * r_shoulder
    tool_w = 2 * r_shoulder
    overlap = (tool_w - step) / tool_w if n_passes > 1 else 1.0

    covered_lo = radii[0] - r_shoulder
    covered_hi = radii[-1] + r_shoulder

    v = travel_speed_mm_min / 60.0                  # mm/s
    path_per_layer = sum(2 * math.pi * r for r in radii)
    t_per_layer = path_per_layer / v

    # mass balance: deposited volume rate must equal fed volume rate
    vol_rate = layer_thickness * step * v           # mm^3/s
    v_feed = vol_rate / rod_area                    # mm/s
    v_feed_mm_min = v_feed * 60.0

    info = dict(radii=radii, step_over=step, overlap=overlap,
                covered=(covered_lo, covered_hi),
                path_per_layer=path_per_layer, t_per_layer=t_per_layer,
                total_deposition=t_per_layer * n_layers,
                vol_rate=vol_rate, feed_rate_mm_min=v_feed_mm_min)

    if verbose:
        print(f"  pass centre radii : "
              f"{', '.join(f'{r:.1f}' for r in radii)} mm")
        print(f"  step-over         : {step:.1f} mm  "
              f"(overlap {100 * overlap:.0f}% of a {tool_w:.0f} mm tool)")
        print(f"  radial coverage   : {covered_lo:.1f} - {covered_hi:.1f} mm "
              f"(ring is {r_inner:.0f} - {r_outer:.0f})")
        print(f"  path per layer    : {path_per_layer:.0f} mm")
        print(f"  time per layer    : {t_per_layer:.0f} s "
              f"({t_per_layer / 60:.1f} min)")
        print(f"  total deposition  : {t_per_layer * n_layers:.0f} s "
              f"({t_per_layer * n_layers / 3600:.2f} h)")
        # Report the volumetric rate in the units machine specs use, since
        # this -- not the linear feed rate -- is what a given AFSD system
        # is actually limited by. Deliberately NOT flagged against a hard
        # bound: published AFSD rates vary by more than an order of
        # magnitude between bench and large-format systems, so this is a
        # number to check against your machine, not one to assert.
        cm3_hr = vol_rate * 3.6
        print(f"  volumetric rate   : {vol_rate:.1f} mm^3/s "
              f"({cm3_hr:.0f} cm^3/h; "
              f"{cm3_hr * 4.43e-3:.1f} kg/h of Ti-6Al-4V, "
              f"{cm3_hr * 2.78e-3:.1f} kg/h of Al)")
        print(f"                      ^ verify against your machine's rated "
              f"deposition rate")
        print(f"  required feed rate: {v_feed_mm_min:.1f} mm/min "
              f"(rod area {rod_area:.1f} mm^2)")

        if covered_lo > r_inner + 1e-6 or covered_hi < r_outer - 1e-6:
            print("  WARNING: passes do not cover the full ring width; "
                  "increase n_passes or the shoulder radius")
        if n_passes > 1 and overlap < 0.15:
            print(f"  WARNING: overlap {100 * overlap:.0f}% is below the "
                  "~20-50% typical of AFSD multi-pass practice")
        if not 10.0 <= v_feed_mm_min <= 500.0:
            print(f"  WARNING: feed rate {v_feed_mm_min:.0f} mm/min is "
                  "outside the ~10-500 mm/min AFSD window -- the tool/rod "
                  "geometry and layer thickness are not mutually consistent")

    return info


def generate_toolpath_ring(r_inner=100.0, r_outer=200.0, r_shoulder=25.0,
                           n_passes=3, n_layers=3, layer_thickness=20.0 / 3.0,
                           travel_speed_mm_min=100.0, dtheta_deg=1.0,
                           dwell_time_s=30.0, reposition_time_s=5.0,
                           rod_area=ROD_1IN_SQUARE,
                           filename="toolpath_ring.crs", verbose=True):
    """Write a multi-pass concentric ring toolpath.

    Conventions, all chosen to match the rest of AFSD-Code:

    * ``z`` is the TOP of the layer currently being deposited, measured from
      the substrate top at z=0. ``assign_birth_time`` in preprocessor.py
      births elements whose ``ele_topz <= toolpath[i,3]``, and the solver
      uses the toolpath z directly as the shoulder plane.

    * A segment is a DEPOSITION segment iff the state at its END waypoint is
      1. This is what ``get_toolpath`` in gamma_afsd.py:161-167 reduces to
      after it re-binarises the interpolated state column, and what
      AFSD_TOOLPOS in abaqus_bridge/afsd_dflux.f implements. Consequently
      the waypoint at the START of an arc carries state 0 (so the approach
      is off) while every subsequent arc waypoint carries state 1.

    * Pass direction alternates, and the radial ordering reverses on
      alternate layers (inner->outer, then outer->inner). Both are standard
      practice; the second is what keeps the multi-layer result from being
      an artefact of a single fixed path pattern.
    """
    if verbose:
        print(f"Ring toolpath: {n_layers} layers x {n_passes} passes, "
              f"shoulder r={r_shoulder} mm")
        info = check_feasibility(r_inner, r_outer, r_shoulder, n_passes,
                                 n_layers, layer_thickness,
                                 travel_speed_mm_min, rod_area)
    else:
        info = check_feasibility(r_inner, r_outer, r_shoulder, n_passes,
                                 n_layers, layer_thickness,
                                 travel_speed_mm_min, rod_area,
                                 verbose=False)

    radii = info["radii"]
    v = travel_speed_mm_min / 60.0
    dtheta = math.radians(dtheta_deg)
    n_steps = max(1, int(round(2 * math.pi / dtheta)))

    rows = []
    t = 0.0

    # start parked at the origin, nothing deposited
    rows.append((t, 0.0, 0.0, 0.0, 0))

    for layer in range(n_layers):
        z = (layer + 1) * layer_thickness
        order = radii if layer % 2 == 0 else list(reversed(radii))

        for ip, R in enumerate(order):
            # +1 / -1: alternate direction pass to pass and layer to layer
            sgn = 1.0 if (ip + layer) % 2 == 0 else -1.0

            # reposition to the start of this arc, deposition OFF.
            # state 0 here makes the approach segment off, and makes the
            # first arc chord (which ends at a state=1 waypoint) on.
            t += reposition_time_s
            rows.append((t, R, 0.0, z, 0))

            # the arc itself
            seg_len = R * dtheta
            for k in range(1, n_steps + 1):
                th = sgn * k * dtheta
                t += seg_len / v
                rows.append((t, R * math.cos(th), R * math.sin(th), z, 1))

            # brief dwell after each pass, deposition OFF
            if not (layer == n_layers - 1 and ip == len(order) - 1):
                t += dwell_time_s
                rows.append((t, R * math.cos(sgn * 2 * math.pi),
                             R * math.sin(sgn * 2 * math.pi), z, 0))

    with open(filename, "w") as f:
        for r in rows:
            f.write(f"{r[0]:<14.4f}{r[1]:<12.4f}{r[2]:<12.4f}"
                    f"{r[3]:<10.4f}{r[4]}\n")

    if verbose:
        sag = radii[-1] * dtheta ** 2 / 8.0
        print(f"  waypoints         : {len(rows)} "
              f"(dtheta={dtheta_deg} deg, chord sagitta "
              f"{1000 * sag:.1f} um at r={radii[-1]:.0f} mm)")
        print(f"  total path time   : {t:.0f} s ({t / 3600:.2f} h)")
        print(f"  written to        : {filename}")
        if len(rows) > 40000:
            print("  WARNING: exceeds MAXTP=40000 in "
                  "abaqus_bridge/afsd_common.inc -- raise it or coarsen "
                  "dtheta_deg")

    return rows, info


if __name__ == "__main__":
    # Recommended configuration. A 50 mm tool (r_shoulder=25) fed by 1-inch
    # square stock is the only combination found that is simultaneously
    # process-realistic and preserves Agarwal's 3 x 6.67 mm layering, which
    # is what makes a like-for-like baseline comparison possible.
    print("=" * 70)
    print("RECOMMENDED: 50 mm tool, 1in square rod, 3 concentric passes")
    print("=" * 70)
    generate_toolpath_ring(filename="toolpath_ring.crs")

    # For contrast: the calibrated 30 mm tool geometry from thinwall.k with
    # Hu's 10 mm round rod. Preserves the validated tool geometry but, as
    # the feed rate shows, demands ~4x the feedstock throughput.
    print()
    print("=" * 70)
    print("CONTRAST: calibrated 30 mm tool + 10 mm round rod (4x the feed rate)")
    print("=" * 70)
    generate_toolpath_ring(r_shoulder=15.0, n_passes=4,
                           rod_area=ROD_10MM_ROUND,
                           filename="toolpath_ring_30mm.crs")

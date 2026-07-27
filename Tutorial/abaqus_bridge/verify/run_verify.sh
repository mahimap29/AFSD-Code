#!/bin/sh
# Phase 1 gate: prove the Fortran AFSD_QVAL reproduces gamma_afsd.py's flux.
#
# Runs the check in both modes:
#   IRMODE=1  full 3-D radius, no filtering -- bit-comparable with
#             gamma_afsd.py, used to prove the translation is faithful
#   IRMODE=0  in-plane radius + tool-plane + upward-face filtering --
#             the production formulation (fixes the underside bleed)
#
# Requires gfortran and python3+numpy. No Abaqus needed: ABA_PARAM.INC
# here is a local stub and driver.f supplies GETOUTDIR/XIT.
set -e
cd "$(dirname "$0")"

echo "== building =="
gfortran -std=legacy -I. -o driver driver.f ../afsd_dflux.f

cp ../afsd_params_example.inp afsd_params.inp
[ -f toolpath.crs ] || cp ../../toolpath.crs .
python3 make_probes.py

for MODE in 1 0; do
    echo
    echo "== IRMODE=$MODE =="
    sed -i "s|^\([0-9.eE+-]*\)  *[01]  *\([0-9.eE+-]*\)$|\1  $MODE  \2|" \
        afsd_params.inp
    ./driver > /dev/null
    python3 compare_flux.py
done

echo
echo "== closed-form power audit =="
python3 power_audit.py

echo
echo "== both modes agree =="

## Abaqus bridge for the Ruvo et al. surface heat source

This folder ports the moving surface heat flux from `update_fluxes()` in
`../gamma_afsd.py` (Eqs. 13–17 of Ruvo et al., *"A surface-based heat source
model for thermomechanical simulation of AFSD process"*, IJHMT 257:128204)
into an Abaqus `DFLUX` user subroutine, so it can replace a prescribed
constant-temperature deposit boundary condition — for example the 1465 °C
placeholder in Agarwal et al. (AIAA 2026-1670), which those authors flag as
future work and want replaced by "a realistic heat flux".

### Files

| File | Purpose |
|---|---|
| `afsd_dflux.f` | The subroutine: `UEXTERNALDB` (loads inputs once), `AFSD_TOOLPOS` (toolpath interpolation), `AFSD_QVAL` (the flux physics), `DFLUX` (Abaqus entry point). |
| `afsd_common.inc` | Shared `COMMON` layout. **Must sit beside `afsd_dflux.f`.** Include it — never hand-copy the declarations (see the warning inside). |
| `afsd_params_example.inp` | Process parameters. Rename to `afsd_params.inp` in the job directory. |
| `verify/` | Standalone equivalence test — runs without Abaqus. See below. |

At runtime the job directory needs `afsd_params.inp` and a `toolpath.crs`
in the same 5-column format AFSD-Code already uses
(`time[s] x[mm] y[mm] z[mm] state[0/1]`, whitespace separated).

### Parameter file

```
5.0  15.0  5.0  41.888  1.67  0.2  0.90     <- P, r_s, r_p, omega, v, mu, slip
1.0  0  0.5                                 <- SCALE, IRMODE, ZTOL
```

Line 1 is the same seven values, in the same order, as the
`*P_Rs_Rp_w_v_mu_slip` block of an AFSD-Code `.k` deck. The values shipped
here are the **actual calibrated values** from `../thinwall.k:12`
(5 MPa, r_s=15 mm, r_p=5 mm, 400 rpm, 100 mm/min, μ=0.2, δ=0.90).

> These appear **nowhere in the paper.** Ruvo et al. write μ(T) and δ(r) but
> publish neither a value, a table, nor a functional form. The released deck
> is the only source. Treat μ and δ as lumped fitted parameters, not
> material properties — see *Calibration* below.

Line 2:
- **`SCALE`** — unit-consistency factor. Read the units section below before
  changing it.
- **`IRMODE`** — `0` production (in-plane radius, tool-plane and
  upward-face filtering); `1` equivalence mode (full 3-D radius, no
  filtering) for cross-checking against the Python solver.
- **`ZTOL`** — tool-plane tolerance in mm, `IRMODE=0` only. Set to a
  fraction of the deposited layer thickness (0.5 mm for a 3 mm layer).

### Units — the most error-prone part of this port

The reference solver and a standard Abaqus mm model do **not** share a unit
system. They differ by 1000 in power, and getting it wrong gives a result
1000× too hot or too cold with no error message.

| Quantity | AFSD-Code (mm-**gram**-s, J/W) | Abaqus mm (mm-**tonne**-s, mJ/mW) |
|---|---|---|
| Heat flux | W/mm² | **mW/mm²** |
| Conductivity | 0.164 W/(mm·K) | 164 mW/(mm·K) |
| Density | 0.00280 g/mm³ | 2.78e-9 tonne/mm³ |
| Specific heat | 0.881 J/(g·K) | 8.81e8 mJ/(tonne·K) |
| Film coefficient | 2e-5 / 1e-3 W/(mm²·K) | 0.02 / 1.0 mW/(mm²·K) |
| Stefan–Boltzmann | 5.6704e-14 | 5.6704e-11 mW/(mm²·K⁴) |
| Stress / modulus | MPa | MPa (unchanged) |

The decisive argument is dimensional, not a matter of copying the reference:

```
[mu][P][omega][r] = [-]·[N/mm²]·[1/s]·[mm] = N/(mm·s) = mJ/(mm²·s) = mW/mm²
```

`μ·P·ω·r` **already is** the flux in a mm-N-s system. AFSD-Code's `/1000`
exists only because that code works in J/W. So:

| Your model | `SCALE` |
|---|---|
| mm-tonne-s (mJ/mW) | **1.0** ← default |
| mm-gram-s (J/W), matching the Python | 1.0e-3 |
| SI (m-kg-s, W/m²) | 1.0e+3, and convert lengths too |

Apply the same conversion to your own Ti-6Al-4V or SS304L property data.
Agarwal's Table 3 covers only 304L — deposit properties are new input.

### Verifying before you trust it

```sh
cd verify && ./run_verify.sh          # needs gfortran + python3/numpy only
```

This is the gate: nothing downstream (μ calibration, ring case) is
meaningful until it passes. It compiles `afsd_dflux.f` against a local
`ABA_PARAM.INC` stub, evaluates `AFSD_QVAL` at probe points covering every
branch, and compares against the same expressions transcribed from
`gamma_afsd.py`. Current status: **both modes agree to 2.7e-16.**

It also prints a closed-form power audit:

```
Intended tool power = 350.016 W  (pin 11.097 + shoulder 338.918)
```

Request history output of the integrated flux over your `*DFLUX` surface and
compare. It must equal this in **mW**. Off by 1000 ⇒ `SCALE` is wrong; short
by some fraction ⇒ your surface does not cover the tool footprint. Both are
worth catching in five minutes rather than after a calibration study.

### Wiring it into a model

1. Put `afsd_params.inp` and `toolpath.crs` in the job's working directory.
2. Apply `*DFLUX` to the **upward-facing** faces of the deposit, magnitude
   `1.0`, and attach **no amplitude** — the subroutine does its own on/off
   gating and an amplitude would multiply it silently. Confirm the data-line
   syntax against the Keywords Reference for your Abaqus version.
3. Convection and radiation are constant-coefficient, so use native `*FILM`
   (sink 298 K; ~20 W/m²K air, ~1000 W/m²K backplate contact, converted per
   the table above) and `*RADIATE` (emissivity 0.7). They are deliberately
   not in the subroutine.
4. Compile: `abaqus job=yourjob user=afsd_dflux.f`.

Cap the increment so the tool moves well under a shoulder radius per
increment — `DFLUX` receives end-of-increment time, so coarse increments let
the tool position lead the thermal solution.

### Corrections made to the reference implementation

Reading `gamma_afsd.py` closely turned up four behaviours that are defects
rather than intended physics. Production mode (`IRMODE=0`) fixes them;
`IRMODE=1` reproduces them for verification.

1. **3-D radius → in-plane radius.** `gamma_afsd.py:738` computes `r` as a
   full 3-D distance from the tool centre, so flux lands on faces at other
   z-levels within `r_shoulder`. In the tutorial geometry the substrate
   underside sits 8 mm below the layer-1 tool plane, well inside
   `r_shoulder = 15 mm`, and receives a spurious **298.6 W against an
   intended 350.0 W — 85% extra heat into the bottom of the plate.**
   (`verify/power_audit.py` derives this in closed form and cross-checks it
   numerically.) A flat shoulder occupies one plane and cannot heat a face
   8 mm below it.
2. **Orientation.** `gamma_afsd.py:133` sets `surface_xy` for any face lying
   in a z=const plane, not only upward-facing ones, so downward-facing
   horizontal faces also collect tool flux. Production mode requires an
   upward normal and `|z − z_tool| ≤ ZTOL`. Together those two tests also
   subsume the reference's surface-death bookkeeping: with one element
   through each deposited layer's thickness, a top face leaves the tool
   plane exactly when the next layer begins, so **no surface reassignment is
   ever needed** — one `*SURFACE`, defined once.
3. **Ramped deposition state.** `gamma_afsd.py:160` interpolates the state
   column, but lines 161–167 then **re-binarise** it, so the reference flux
   is strictly on/off. An earlier version of this subroutine interpolated
   the state and applied a ramping partial flux across every lead-in
   segment. Fixed: a segment is ON iff the state at its *end* waypoint is 1,
   which is provably equivalent to the reference loop (derivation in the
   `AFSD_TOOLPOS` header).
4. **Double-counted bottom cooling.** The reference `surface_flux` mask is
   identically 1 — its test `min(face z) >= global_min_z` is trivially true
   — so convection *and* radiation are applied to every exposed face
   including the underside, which *separately* also receives `h_conv_2` at
   its nodes (`gamma_afsd.py:777-782`). Not reproduced. Apply exactly one
   condition to the underside: backplate contact conductance.

Also worth knowing: `gamma_afsd.py:498` casts nodal coordinates to `int64`
before the `surface_xy` planarity test. It survives the tutorial mesh by
luck and will misclassify faces on any mesh with sub-1 mm elements.

### On the paper's Eq. 17

As printed, Eq. 17 reads `δ·P·μ·ω·r + (1−δ)·P·μ·ω·r`, in which **δ cancels
identically** and the slip factor has no effect whatsoever. The released
code (`gamma_afsd.py:742`) uses `P/√3` — the von Mises shear limit — for the
stick term, so δ does matter. **The code is right and the paper has a
typographical error.** This subroutine follows the code. Any write-up
transplanting this model should print the code's form, cite Eq. 17 as
published, and note the discrepancy rather than propagating it.

Separately: the stick term uses the *applied pressure* `P/√3` where the FSW
literature (Schmidt & Hattel) uses the *material's* shear yield `σ_y(T)/√3`.
At P = 5 MPa that is ~2.9 MPa against ~30 MPa for AA2024 at 400 °C. This is
one reason μ here is a lumped fitting parameter rather than a friction
coefficient.

### Calibration

μ and δ are **not both identifiable** from thermocouple data. Total power
depends on them essentially only through

```
K = δ·μ + (1 − δ)/√3
```

so lines of constant `K` form a one-parameter family and the objective is
rank-deficient in (μ, δ). The only term breaking the degeneracy is the pin
branch, which is ~3% of total power (11.1 W of 350.0 W). **Fix δ and
calibrate μ alone**, then report `K*` and the iso-K curve so others can
transfer the result under a different slip assumption.

Note also that only the *product* `P·K` is identifiable, not μ separately —
any error in a reported axial force is absorbed directly into μ.

**Force → pressure.** `P` is applied over the whole shoulder disc: the code
uses the same `pressure` in both the pin branch (`r < r_p`) and the shoulder
branch (`r_p ≤ r < r_s`). Self-consistency with the total axial force
therefore requires

```
P = F_axial / (pi * r_s^2)
```

For Hu et al.'s 30 mm tool that gives 0.42 MPa (0.3 kN preheat) and
0.71 MPa (0.5 kN deposition). Flag honestly: real AFSD axial forces on a
30 mm tool are 20–90 kN, so 0.3–0.5 kN is either an unusually light lab
configuration or a differently-defined quantity — verify against the primary
source before calibrating, and record the substrate thickness at the same
time, since it determines whether defect 1's bleed was active when μ = 0.2
was fitted.

**μ does not transfer across material pairs.** A μ fitted for
AA2024-on-AA2024 has no defensible claim on Ti-6Al-4V-on-SS304L, and the
transfer is doubly ill-posed because Agarwal specify no axial force at all
for the ring. Anchor on absolute power instead (check the total against the
machine envelope), use the solid-state peak-temperature window as an
independent admissibility constraint, and report a sensitivity envelope
rather than a point value. The measurement that would close the gap is
spindle **torque**, which constrains power directly via `Q = T·ω`.

### What this does not do

`AFSD_QVAL` computes *what* flux applies at a point and time. It does not
track element activation. Drive that with the AM Modeler plug-in's
progressive activation (or `*MODEL CHANGE`), ideally from a birth table
generated by the same `assign_birth_time` logic in `../preprocessor.py` so
both solvers share one activation history — but emit it in **original mesh
numbering**, because `gamma_afsd.py:431-444` renumbers into birth order and
Abaqus must not be renumbered.

**Do not route this model through `UMDFlux`.** That interface returns a
*volumetric* flux at material points inside a moving box, which would smear
this surface flux over an element thickness and reintroduce exactly the
volumetric approximation the paper set out to replace — making layer
thickness an unphysical parameter of the heat source. Use `DFLUX` for the
flux and the plug-in for activation only; the two are orthogonal.

Finally, Abaqus permits only **one** `UEXTERNALDB` per job. If you add a
DAMASK/UMAT layer that also wants it, merge both into a single dispatcher
rather than defining a second one.

## Abaqus DFLUX bridge for the Ruvo et al. surface heat source

This folder ports the moving surface heat flux from `update_fluxes()` in
`../gamma_afsd.py` (implementing Eqs. 13-17 of Ruvo et al., "A surface-based
heat source model for thermomechanical simulation of AFSD process") into an
Abaqus `DFLUX` user subroutine, so it can be used in place of a simplified
constant-temperature deposit boundary condition in an Abaqus/Standard model.

### Files

- `afsd_dflux.f` - the subroutine itself: `UEXTERNALDB` (reads inputs once),
  `AFSD_TOOLPOS` (toolpath interpolation), `AFSD_QVAL` (the flux physics),
  `DFLUX` (the entry point Abaqus calls).
- `afsd_params_example.inp` - example process-parameter file, same 7 values
  and order as the `*P_Rs_Rp_w_v_mu_slip` block in an AFSD-Code `.k` deck:
  `pressure[MPa] r_shoulder[mm] r_pin[mm] omega[rad/s] v[mm/s] mu[-] slip[-]`.
  The numbers here (5 MPa, rs=15mm, rp=5mm, 400 rpm -> 41.888 rad/s,
  100 mm/min -> 1.667 mm/s) match the Tutorial case in `../Readme.md`; `mu`
  and `slip` are placeholders and must be replaced with calibrated values.

At runtime you also need a `toolpath.crs` in the job directory - the exact
same file format AFSD-Code already uses (`time[s] x[mm] y[mm] z[mm]
state[0/1]`, whitespace separated), e.g. generated with
`../generate_toolpath.py`.

### Wiring it into your model

1. Place `afsd_params_example.inp` (renamed `afsd_params.inp`) and
   `toolpath.crs` in the Abaqus job's working/output directory.
2. Reference the top surface of the deposit/wall with a `*DFLUX` load of
   type `S` (or `Sn` per element face), e.g.:
   ```
   *DFLUX
   WallTopSurf, S1, 1.
   ```
   Consult the *DFLUX entry in the Abaqus Keywords Reference for your
   version to confirm the exact data-line syntax - the key requirement is
   that this surface load references a heat-transfer (or coupled
   temp-disp) surface for which the user subroutine `DFLUX` is picked up
   automatically at every increment/integration point.
3. Convection and radiation (Eqs. 11-12) are constant-coefficient, so use
   Abaqus's native surface conditions directly instead of a subroutine:
   `*FILM` (sink temperature 298 K; film coefficient ~20 W/m^2K on
   air-exposed faces, ~1000 W/m^2K on the backplate-contact face) and
   `*RADIATE` (emissivity 0.7).
4. Compile/link `afsd_dflux.f` with your job, e.g.
   `abaqus job=yourjob user=afsd_dflux.f`.

### What this does *not* replicate

`gamma_afsd.py` tracks which surface is currently the exposed top of the
deposit via element birth/death and an `active_surface`/`surface_xy` mask -
this subroutine only computes the flux *value* at a given point and time,
not which surface it should be applied to. You still need to reproduce the
dynamic surface exposure yourself, either by:

- driving `*MODEL CHANGE` element activation in lockstep with the toolpath
  (the manual analogue of `element_birth` in AFSD-Code) and only applying
  `*DFLUX` to the current top surface, or
- using the Abaqus AM Modeler plug-in's `UEPACTIVATIONVOL` + `UMDFLUX` pair,
  which automates that tracking. If you go this route, drop the body of
  `AFSD_QVAL` into your `UMDFLUX` implementation - the physics is identical,
  only the calling signature differs, and that signature varies across
  Abaqus releases, so check the User Subroutines Reference Guide for your
  version rather than assuming the `DFLUX` argument list shown here applies.

Also note the friction coefficient `mu` and slip factor are treated as
scalar constants here, exactly as in the AFSD-Code `.k` deck. The paper's
more general `mu(T)` and `delta(r)` are not implemented; if you want
temperature- or radius-dependent friction, extend `AFSD_QVAL` using the
`TEMP` argument passed into `DFLUX` and/or the computed `RR`.

Finally, double check the unit-consistency factor `SCALE` in
`UEXTERNALDB` (currently `1/1000`, matching AFSD-Code's mm-N-s unit
system) against whatever consistent unit system your own Abaqus model
uses.

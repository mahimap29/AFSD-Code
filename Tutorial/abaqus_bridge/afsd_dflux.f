C-----------------------------------------------------------------------------
C afsd_dflux.f
C
C Abaqus user subroutine DFLUX implementing the Ruvo et al. surface-based
C moving heat source model ("A surface-based heat source model for
C thermomechanical simulation of AFSD process", Eqs. 13-17) as a
C term-by-term Fortran translation of update_fluxes() in
C AFSD-Code/Tutorial/gamma_afsd.py (lines 711-757).
C
C Intended use: replace a prescribed constant-temperature deposit BC
C (e.g. the 1465 degC placeholder in Agarwal et al., AIAA 2026-1670)
C with the physically-derived friction heat flux at the feed-rod and
C tool-shoulder interfaces, as a function of process parameters, radial
C position, and time.
C
C Companion routines:
C   UEXTERNALDB  - reads process parameters and the toolpath once, at
C                  the start of the analysis, into COMMON blocks.
C   AFSD_TOOLPOS - piecewise-linear interpolation of tool centre
C                  position + deposition on/off state at a given time,
C                  mirroring get_toolpath()/np.interp() in gamma_afsd.py.
C   AFSD_QVAL    - the flux-magnitude physics (Eqs. 14 and 17), factored
C                  out so it can also be called from UMDFLUX if you use
C                  the Abaqus AM Modeler plug-in (see README.md).
C
C-----------------------------------------------------------------------------
C UNITS  --  READ THIS BEFORE RUNNING ANYTHING
C-----------------------------------------------------------------------------
C The reference Python solver and a standard Abaqus mm model do NOT use
C the same unit system. They differ by 1000 in power, so the wrong SCALE
C silently gives a result that is 1000x too hot or too cold.
C
C   AFSD-Code (gamma_afsd.py) is  mm - GRAM - s - K - J/W:
C     flux W/mm^2 | k 0.164 W/(mm K) | rho 0.00280 g/mm^3
C     Cp 0.881 J/(g K) | film 2e-5 / 1e-3 W/(mm^2 K)
C     Stefan-Boltzmann 5.6704e-14 W/(mm^2 K^4)
C   Its "/1000" in q1,q2 converts mW/mm^2 -> W/mm^2.
C
C   A standard Abaqus mm model is mm - TONNE - s - K - mJ/mW:
C     flux mW/mm^2 | k 164 mW/(mm K) | rho 2.78e-9 tonne/mm^3
C     Cp 8.81e8 mJ/(tonne K) | film 0.02 / 1.0 mW/(mm^2 K)
C     Stefan-Boltzmann 5.6704e-11 mW/(mm^2 K^4)
C
C Because P[MPa] * omega*r[mm/s] is ALREADY mW/mm^2, the correct factor
C for a mm-tonne-s Abaqus model is SCALE = 1.0 -- the Python "/1000"
C must be DROPPED. SCALE is read from afsd_params.inp so you can change
C it without recompiling. Guidance:
C     mm-tonne-s (mJ/mW)     -> SCALE = 1.0        <-- default
C     mm-gram-s  (J/W)       -> SCALE = 1.0E-3     (matches Python)
C     SI (m-kg-s, W/m^2)     -> SCALE = 1.0E+3     (and convert lengths)
C
C Sanity check: with P=5 MPa, mu=0.2, omega=41.888 rad/s, v=1.67 mm/s,
C at r = r_pin = 5 mm, q1 = 211.1 mW/mm^2 (= 0.2111 W/mm^2 = 211 kW/m^2).
C-----------------------------------------------------------------------------
C DELIBERATE DEVIATIONS FROM THE REFERENCE PYTHON IMPLEMENTATION
C-----------------------------------------------------------------------------
C The reference update_fluxes() has three behaviours that are defects
C rather than intended physics. This file fixes them by default and can
C reproduce them on request (IRMODE / afsd_params.inp) for verification.
C
C 1. RADIUS. gamma_afsd.py:738 computes r as a full 3-D distance from the
C    tool centre. Consequently flux bleeds onto faces at other z levels
C    within r_shoulder -- including the substrate UNDERSIDE, which picks
C    up q2 at r = 8 mm while the tool is at z = 3. Default here is the
C    in-plane radius sqrt((x-xt)^2 + (y-yt)^2). Set IRMODE=1 for the 3-D
C    form when checking equivalence against the Python solver.
C
C 2. ORIENTATION. gamma_afsd.py:133 sets surface_xy for any face lying in
C    a z=const plane, NOT only upward-facing ones, so downward-facing
C    horizontal faces also receive tool flux. Here the caller supplies
C    the face normal and flux is applied only when it points up (see the
C    ZNORM argument of AFSD_QVAL and the JLTYP/SNAME note in DFLUX).
C
C 3. COOLING. gamma_afsd.py's surface_flux mask is identically 1 (its
C    test min(face z) >= global_min_z is trivially true), so convection
C    AND radiation are applied to every exposed face including the
C    underside, which separately also receives h_conv_2 at its nodes
C    (gamma_afsd.py:777-782) -- i.e. double-counted bottom cooling. Not
C    reproduced here: use native *FILM / *RADIATE and apply the
C    backplate coefficient to the bottom face only once.
C
C Convection/radiation are intentionally NOT implemented in this file.
C They are constant-coefficient, so use Abaqus *FILM (sink 298 K, ~20
C W/m^2K air / ~1000 W/m^2K backplate, converted per the table above)
C and *RADIATE (emissivity 0.7).
C-----------------------------------------------------------------------------
C REMAINING CAVEATS
C-----------------------------------------------------------------------------
C a) mu and the slip factor are scalar constants here, matching the
C    released deck (thinwall.k:12 -> mu=0.2, slip=0.90). Ruvo et al.
C    write mu(T) and delta(r) but publish neither a value nor a
C    functional form. To add temperature dependence, use the TEMP
C    argument of DFLUX; for radial dependence use RR in AFSD_QVAL.
C b) This routine computes WHAT flux applies at a point and time. It
C    does NOT track WHICH surface is currently the exposed top of the
C    deposit -- in gamma_afsd.py that is done by element birth/death
C    plus the half-open active_surface window
C    [surface_birth(:,0), surface_birth(:,1)). Reproduce it with
C    *MODEL CHANGE, or delegate to the AM plug-in's UEPACTIVATIONVOL.
C c) The printed Eq. 17 of the paper is degenerate: as published it
C    reads delta*P*mu*w*r + (1-delta)*P*mu*w*r, in which delta cancels
C    identically. The released CODE uses P/sqrt(3) (von Mises shear)
C    for the stick term, so delta does matter. This file follows the
C    CODE. See AFSD_QVAL.
C-----------------------------------------------------------------------------

      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
C     Reads the AFSD process parameters and the toolpath once, at the
C     very start of the analysis, into COMMON blocks shared with
C     AFSD_TOOLPOS and AFSD_QVAL.
C
C     NOTE: Abaqus permits only ONE UEXTERNALDB per job. If you later
C     add the DAMASK/UMAT coupling layer (which also wants UEXTERNALDB
C     for its I/O), merge that logic into this routine rather than
C     defining a second one.
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2)
C
      INCLUDE 'afsd_common.inc'
C
      CHARACTER*256 OUTDIR,PARMFILE,PATHFILE
C
      IF (LOP.EQ.0) THEN
C
        CALL GETOUTDIR(OUTDIR,LENOUTDIR)
C
C       ---- process parameters ---------------------------------------
C       afsd_params.inp, two lines:
C         line 1: the same 7 values, in the same order, as the
C                 *P_Rs_Rp_w_v_mu_slip block of an AFSD-Code .k deck:
C                 pressure[MPa] r_shoulder[mm] r_pin[mm] omega[rad/s]
C                 v[mm/s] mu[-] slip[-]
C         line 2: SCALE  IRMODE  ZTOL
C                 SCALE  = unit factor (see header; 1.0 for mm-tonne-s)
C                 IRMODE = 0 in-plane radius + tool-plane + upward-face
C                            filtering (default, physical)
C                          1 full 3-D radius, no filtering: reproduces
C                            the Python reference for equivalence
C                            checking only
C                 ZTOL   = tool-plane tolerance [mm], IRMODE=0 only.
C                          Set to a fraction of the deposited layer
C                          thickness (e.g. 0.5 mm for a 3 mm layer).
        PARMFILE = OUTDIR(1:LENOUTDIR)//'/afsd_params.inp'
        OPEN(UNIT=91,FILE=PARMFILE,STATUS='OLD')
        READ(91,*) PRES,RSHLD,RPIN,OMEGA,VTRAV,FRIC,SLIP
        READ(91,*) SCALE,IRMODE,ZTOL
        CLOSE(91)
C
C       ---- toolpath -------------------------------------------------
C       toolpath.crs: same 5-column layout AFSD-Code already uses:
C         time[s]  x[mm]  y[mm]  z[mm]  state[0 or 1]
C       Circular paths must arrive pre-discretised into short chords
C       (see generate_toolpath_ring.py) because the interpolation below
C       is piecewise LINEAR -- hence the large MAXTP.
        PATHFILE = OUTDIR(1:LENOUTDIR)//'/toolpath.crs'
        OPEN(UNIT=92,FILE=PATHFILE,STATUS='OLD')
        NTP = 0
   10   CONTINUE
        NTP = NTP+1
        IF (NTP.GT.MAXTP) THEN
          WRITE(7,*) 'AFSD_DFLUX: MAXTP exceeded, increase MAXTP'
          CALL XIT
        END IF
        READ(92,*,END=20) TPT(NTP),TPX(NTP),TPY(NTP),TPZ(NTP),
     1                     TPS(NTP)
        GO TO 10
   20   CONTINUE
        NTP = NTP-1
        CLOSE(92)
C
        IF (NTP.LT.2) THEN
          WRITE(7,*) 'AFSD_DFLUX: toolpath.crs needs >= 2 waypoints'
          CALL XIT
        END IF
C
C       ---- echo, so the .msg file records what was actually used ----
        WRITE(7,*) 'AFSD_DFLUX: P,Rs,Rp,w,v,mu,slip = ',
     1              PRES,RSHLD,RPIN,OMEGA,VTRAV,FRIC,SLIP
        WRITE(7,*) 'AFSD_DFLUX: SCALE,IRMODE,ZTOL,NTP = ',
     1              SCALE,IRMODE,ZTOL,NTP
C
      END IF
C
      RETURN
      END


      SUBROUTINE AFSD_TOOLPOS(TTIME,XT,YT,ZT,STATE)
C
C     Piecewise-linear interpolation of tool centre position and
C     deposition on/off state at total time TTIME. Direct analogue of
C     get_toolpath()'s np.interp() calls in gamma_afsd.py (lines
C     154-168), including the "hold flux off past the last waypoint"
C     behaviour (laser_state * (time <= toolpath_raw[-1,0])).
C
C     STATE is BINARY, not interpolated. gamma_afsd.py:160 does call
C     np.interp on the state column, but lines 161-167 then re-binarise
C     the result:
C         if laser_state[i+1] > laser_state[i] or laser_state[i] == 1:
C             l[i] = 1
C         else:
C             l[i] = 0
C     so the reference flux is strictly full-on or full-off. Working
C     through the four segment cases, that loop is exactly equivalent to
C     "the segment is ON iff the state at its END waypoint is 1":
C         TPS(I)=0, TPS(I+1)=1  (rising)  -> l=1  -> ON
C         TPS(I)=1, TPS(I+1)=1  (holding) -> l=1  -> ON
C         TPS(I)=1, TPS(I+1)=0  (falling) -> l=0  -> OFF
C         TPS(I)=0, TPS(I+1)=0  (idle)    -> l=0  -> OFF
C     Hence STATE = 1 iff TPS(ILO+1) >= 1 below. An earlier version of
C     this routine linearly interpolated STATE; that was a bug -- it
C     applied a ramping partial flux across every lead-in segment.
C
C     Uses a bisection search so cost is O(log NTP) per integration
C     point per increment. With ~15000 waypoints for a multi-pass ring
C     a linear scan here would dominate the analysis.
C
      INCLUDE 'ABA_PARAM.INC'
C
      INCLUDE 'afsd_common.inc'
C
      IF (TTIME.LE.TPT(1)) THEN
        XT = TPX(1)
        YT = TPY(1)
        ZT = TPZ(1)
C       flux is off before the first waypoint (never "on" at t<=t_0)
        STATE = 0.D0
        RETURN
      END IF
C
      IF (TTIME.GE.TPT(NTP)) THEN
        XT = TPX(NTP)
        YT = TPY(NTP)
        ZT = TPZ(NTP)
        STATE = 0.D0
        RETURN
      END IF
C
C     bisection: find ILO with TPT(ILO) <= TTIME < TPT(ILO+1)
      ILO = 1
      IHI = NTP
   30 CONTINUE
      IF (IHI-ILO.GT.1) THEN
        IMID = (ILO+IHI)/2
        IF (TPT(IMID).GT.TTIME) THEN
          IHI = IMID
        ELSE
          ILO = IMID
        END IF
        GO TO 30
      END IF
C
      DT = TPT(ILO+1)-TPT(ILO)
      IF (DT.LE.0.D0) THEN
        FRAC = 0.D0
      ELSE
        FRAC = (TTIME-TPT(ILO))/DT
      END IF
      XT = TPX(ILO)+FRAC*(TPX(ILO+1)-TPX(ILO))
      YT = TPY(ILO)+FRAC*(TPY(ILO+1)-TPY(ILO))
      ZT = TPZ(ILO)+FRAC*(TPZ(ILO+1)-TPZ(ILO))
C
C     binary segment state -- see the derivation in the header comment
      IF (TPS(ILO+1).GE.1.D0) THEN
        STATE = 1.D0
      ELSE
        STATE = 0.D0
      END IF
C
      RETURN
      END


      SUBROUTINE AFSD_QVAL(X,Y,Z,ZNORM,TTIME,Q)
C
C     Flux-magnitude physics, translated term-by-term from
C     update_fluxes() in gamma_afsd.py (lines 737-744):
C
C       r    = |ip_pos - tool_center|
C       q1   = (mu*P*w*r + mu*P*v) / 1000                     [Eq. 14]
C       q2   = (slip*mu*P*w*r)/1000
C              + ((1-slip)*P*w*r/sqrt(3))/1000                [Eq. 17]
C       qmov = where(r<r_pin, q1, where(r<r_shoulder, q2, 0))
C              * surface_xy * laser_state
C
C     Arguments:
C       X,Y,Z  integration point coordinates [mm]
C       ZNORM  z component of the OUTWARD face normal. Flux is applied
C              only when ZNORM > 0 (upward-facing), which is the fix for
C              defect 2 in the header. Pass +1.0 to bypass the check.
C       TTIME  total analysis time [s]
C       Q      returned flux, SCALE-consistent units (mW/mm^2 by default)
C
C     On Eq. 17: the stick term uses P/sqrt(3), the von Mises shear
C     limit, NOT P*mu. This follows the released code. The paper's
C     printed Eq. 17 uses P*mu for both terms, which makes delta cancel
C     out identically and renders the slip factor inoperative -- a typo.
C
C     Note the area weighting is NOT applied here. Q is a true flux;
C     Abaqus supplies the face area, exactly as surf_detJac does in the
C     Python reference (gamma_afsd.py:755-757).
C
      INCLUDE 'ABA_PARAM.INC'
C
      INCLUDE 'afsd_common.inc'
C
      Q = 0.D0
C
      CALL AFSD_TOOLPOS(TTIME,XT,YT,ZT,STATE)
C
      IF (STATE.LE.0.D0) RETURN
C
      IF (IRMODE.EQ.1) THEN
C       ---- equivalence mode: reproduce gamma_afsd.py exactly --------
C       3-D distance, and NO orientation or tool-plane filtering, since
C       the reference surface_xy mask accepts any z=const face
C       regardless of which way it faces. Use only for cross-checking
C       against the Python solver; not physically defensible (see
C       defect 1 in the header).
        RR = SQRT((X-XT)**2+(Y-YT)**2+(Z-ZT)**2)
      ELSE
C       ---- production mode -----------------------------------------
C       Fixes defects 1 and 2 together. The shoulder is a flat disc
C       occupying ONE plane, so a face only sees it if it lies in that
C       plane (|z - z_tool| <= ZTOL) and faces upward (ZNORM > 0).
C       Together these two tests also subsume the reference's
C       surface-death bookkeeping: with one element through the
C       thickness of each deposited layer, a top face stops being at
C       the tool plane exactly when the next layer starts, so no
C       surface reassignment is ever needed.
        IF (ZNORM.LE.0.D0) RETURN
        IF (ABS(Z-ZT).GT.ZTOL) RETURN
        RR = SQRT((X-XT)**2+(Y-YT)**2)
      END IF
C
      IF (RR.LT.RPIN) THEN
C       Eq. 14 - feed rod / pin-deposit interface (0 <= r < r_pin)
        Q = (FRIC*PRES*OMEGA*RR + FRIC*PRES*VTRAV)*SCALE
      ELSE IF (RR.LT.RSHLD) THEN
C       Eq. 17 - tool shoulder-deposit interface (r_pin <= r < r_shoulder)
        Q = (SLIP*FRIC*PRES*OMEGA*RR)*SCALE
     1    + ((1.D0-SLIP)*PRES*OMEGA*RR/SQRT(3.D0))*SCALE
      ELSE
        Q = 0.D0
        RETURN
      END IF
C
      Q = Q*STATE
C
      RETURN
      END


      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     1 JLTYP,TEMP,PRESS,SNAME)
C
C     Abaqus entry point. Verify this argument list against the User
C     Subroutines Reference for YOUR Abaqus release before relying on
C     it -- the signature has varied across versions.
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION FLUX(2),TIME(2),COORDS(3)
      CHARACTER*80 SNAME
C
C     ZNORM: z component of the outward normal of the loaded face.
C     DFLUX is not passed the face normal, so it cannot be computed
C     here. Restrict the *DFLUX surface definition to upward-facing
C     faces of the deposit (which is what you want anyway -- see defect
C     2 in the header) and leave this at +1.0. If you instead apply
C     *DFLUX to a surface containing faces of mixed orientation, this
C     routine cannot filter them and the flux WILL land on undersides.
      ZNORM = 1.D0
C
C     TIME(2) is total accumulated analysis time across steps, so one
C     continuous toolpath can span multiple *STEPs -- matching the
C     domain.current_time accumulation in AFSD-Code's time_integration().
      CALL AFSD_QVAL(COORDS(1),COORDS(2),COORDS(3),ZNORM,TIME(2),QQ)
C
      FLUX(1) = QQ
C     No temperature dependence in this translation (mu, P, omega, v and
C     slip are all constants, matching gamma_afsd.py), so the flux
C     Jacobian dq/dT is exactly zero. Leave it at zero rather than let
C     Abaqus estimate it numerically. If you add mu(T), you must supply
C     d(flux)/dT here.
      FLUX(2) = 0.D0
C
      RETURN
      END

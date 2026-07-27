C-----------------------------------------------------------------------------
C afsd_dflux.f
C
C Abaqus user subroutine DFLUX implementing the Ruvo et al. surface-based
C moving heat source model ("A surface-based heat source model for
C thermomechanical simulation of AFSD process", Eqs. 13-17) as a direct,
C term-by-term Fortran translation of update_fluxes() in
C AFSD-Code/Tutorial/gamma_afsd.py.
C
C This subroutine is a drop-in replacement for a constant-temperature
C deposit BC (e.g. Agarwal et al.'s current heat-source simplification):
C instead of prescribing a fixed temperature on newly activated elements,
C it applies the physically-derived friction heat flux at the feed-rod
C and tool-shoulder interfaces as a function of process parameters,
C radial position, and time.
C
C Companion routines:
C   UEXTERNALDB  - reads process parameters and the toolpath (once, at
C                  the start of the analysis) into COMMON blocks.
C   AFSD_TOOLPOS - piecewise-linear interpolation of tool centre
C                  position + deposition on/off state at a given time,
C                  mirroring get_toolpath()/np.interp() in gamma_afsd.py.
C   AFSD_QVAL    - the actual flux-magnitude physics (Eqs. 14 and 17),
C                  factored out so it can also be called from UMDFLUX
C                  if you are using the Abaqus AM Modeler plug-in
C                  instead of a plain *DFLUX (see README.md).
C
C IMPORTANT - things you must verify/adjust for your model:
C   1. Units. AFSD-Code uses a mm-N-s(-tonne) consistent unit system
C      (pressure in MPa, lengths in mm, velocities in mm/s, angular
C      velocity in rad/s) and divides q1,q2 by 1000 (variable SCALE
C      below) to land on the flux units consistent with that system.
C      If your Abaqus model uses a different consistent unit system,
C      recompute SCALE accordingly - do not assume 1/1000 is universal.
C   2. Friction coefficient mu and slip factor are treated as scalar
C      constants, exactly as domain.mu / domain.slip are in the
C      AFSD-Code .k deck (*P_Rs_Rp_w_v_mu_slip). The paper's mu(T) and
C      delta(r) are more general; if you want temperature- or radius-
C      dependent friction/slip, use the TEMP argument passed into
C      DFLUX and/or RR computed in AFSD_QVAL to look those up instead.
C   3. Surface selection. This subroutine only computes WHAT flux to
C      apply at a given point and time - it does not track WHICH
C      surface is currently exposed as the tool passes over it. In
C      gamma_afsd.py that bookkeeping is done by the element birth/
C      death + active_surface + surface_xy masks. In Abaqus you must
C      reproduce that by either (a) driving element activation with
C      *MODEL CHANGE in lockstep with the toolpath (the manual
C      equivalent of AFSD-Code's element birth times) and applying
C      *DFLUX only to the current top surface of the active deposit,
C      or (b) using the AM Modeler plug-in's UEPACTIVATIONVOL +
C      UMDFLUX pair, which automates that tracking for you.
C   4. Convection/radiation (Eqs. 11-12) are NOT reproduced here on
C      purpose - they are constant-coefficient, so just use Abaqus's
C      native *FILM (sink temp 298 K, film coeff ~20 W/m^2K on
C      air-exposed faces, ~1000 W/m^2K on the backplate-contact face)
C      and *RADIATE (emissivity 0.7) surface conditions directly; no
C      subroutine needed for those two terms.
C-----------------------------------------------------------------------------

      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
C     Reads the AFSD process parameters and the toolpath once, at the
C     very start of the analysis, into COMMON blocks shared with
C     AFSD_TOOLPOS and AFSD_QVAL.
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2)
C
      PARAMETER (MAXTP=2000)
      COMMON /AFSD_PATH/ NTP,TPT(MAXTP),TPX(MAXTP),TPY(MAXTP),
     1                    TPZ(MAXTP),TPS(MAXTP)
      COMMON /AFSD_PARM/ PRES,RSHLD,RPIN,OMEGA,VTRAV,FRIC,SLIP,SCALE
C
      CHARACTER*256 OUTDIR,PARMFILE,PATHFILE
C
      IF (LOP.EQ.0) THEN
C
        CALL GETOUTDIR(OUTDIR,LENOUTDIR)
C
C       ---- process parameters ---------------------------------------
C       afsd_params.inp: one line, whitespace separated, same 7 values
C       and order as *P_Rs_Rp_w_v_mu_slip in the AFSD-Code .k deck:
C         pressure[MPa]  r_shoulder[mm]  r_pin[mm]  omega[rad/s]
C         v[mm/s]  mu[-]  slip[-]
        PARMFILE = OUTDIR(1:LENOUTDIR)//'/afsd_params.inp'
        OPEN(UNIT=91,FILE=PARMFILE,STATUS='OLD')
        READ(91,*) PRES,RSHLD,RPIN,OMEGA,VTRAV,FRIC,SLIP
        CLOSE(91)
C
C       Unit-consistency scale factor - see note (1) above. Default
C       reproduces AFSD-Code's /1000 in q1, q2 (gamma_afsd.py L741-742).
        SCALE = 1.D0/1000.D0
C
C       ---- toolpath ----------------------------------------------------
C       toolpath.crs: same 5-column layout AFSD-Code already uses:
C         time[s]  x[mm]  y[mm]  z[mm]  state[0 or 1]
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
      END IF
C
      RETURN
      END


      SUBROUTINE AFSD_TOOLPOS(TTIME,XT,YT,ZT,STATE)
C
C     Piecewise-linear interpolation of tool centre position and
C     deposition on/off state at total time TTIME. Direct analogue of
C     get_toolpath()'s np.interp() calls in gamma_afsd.py (L154-168),
C     including the "hold flux off past the last waypoint" behaviour
C     (laser_state * (time <= toolpath_raw[-1,0])).
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (MAXTP=2000)
      COMMON /AFSD_PATH/ NTP,TPT(MAXTP),TPX(MAXTP),TPY(MAXTP),
     1                    TPZ(MAXTP),TPS(MAXTP)
C
      IF (TTIME.LE.TPT(1)) THEN
        XT = TPX(1)
        YT = TPY(1)
        ZT = TPZ(1)
        STATE = TPS(1)
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
      DO 30 I=1,NTP-1
        IF (TTIME.GE.TPT(I) .AND. TTIME.LE.TPT(I+1)) THEN
          FRAC = (TTIME-TPT(I))/(TPT(I+1)-TPT(I))
          XT = TPX(I)+FRAC*(TPX(I+1)-TPX(I))
          YT = TPY(I)+FRAC*(TPY(I+1)-TPY(I))
          ZT = TPZ(I)+FRAC*(TPZ(I+1)-TPZ(I))
          STATE = TPS(I)+FRAC*(TPS(I+1)-TPS(I))
          RETURN
        END IF
   30 CONTINUE
C
      RETURN
      END


      SUBROUTINE AFSD_QVAL(X,Y,Z,TTIME,Q)
C
C     Flux-magnitude physics, translated term-by-term from
C     update_fluxes() in gamma_afsd.py (L737-744):
C
C       r      = |ip_pos - tool_center|
C       q1     = (mu*P*w*r + mu*P*v) / 1000              [Eq. 14]
C       q2     = (slip*mu*P*w*r)/1000
C                + ((1-slip)*P*w*r/sqrt(3))/1000          [Eq. 17]
C       qmov   = where(r<r_pin, q1, where(r<r_shoulder, q2, 0))
C                * laser_state
C
C     (surface_xy is NOT reproduced - see note (3) at the top of this
C     file; it is handled by which surface you attach *DFLUX to.)
C
      INCLUDE 'ABA_PARAM.INC'
C
      COMMON /AFSD_PARM/ PRES,RSHLD,RPIN,OMEGA,VTRAV,FRIC,SLIP,SCALE
C
      CALL AFSD_TOOLPOS(TTIME,XT,YT,ZT,STATE)
C
      RR = SQRT((X-XT)**2+(Y-YT)**2+(Z-ZT)**2)
C
C     Eq. 14 - feed rod / pin-deposit interface (0 <= r < r_pin)
      QQ1 = (FRIC*PRES*OMEGA*RR + FRIC*PRES*VTRAV)*SCALE
C
C     Eq. 17 - tool shoulder-deposit interface (r_pin <= r < r_shoulder)
      QQ2 = (SLIP*FRIC*PRES*OMEGA*RR)*SCALE
     1     + ((1.D0-SLIP)*PRES*OMEGA*RR/SQRT(3.D0))*SCALE
C
      IF (RR.LT.RPIN) THEN
        Q = QQ1
      ELSE IF (RR.LT.RSHLD) THEN
        Q = QQ2
      ELSE
        Q = 0.D0
      END IF
C
      Q = Q*STATE
C
      RETURN
      END


      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     1 JLTYP,TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION FLUX(2),TIME(2),COORDS(3)
      CHARACTER*80 SNAME
C
C     TIME(2) is total (accumulated, across steps) analysis time, used
C     so one continuous toolpath can span multiple *STEPs, matching
C     domain.current_time accumulation in AFSD-Code's time_integration().
      CALL AFSD_QVAL(COORDS(1),COORDS(2),COORDS(3),TIME(2),QQ)
C
      FLUX(1) = QQ
C     No temperature dependence in this translation (mu, P, w, v, slip
C     are constants, matching gamma_afsd.py) so leave the Jacobian term
C     at zero rather than let Abaqus perturb/estimate it numerically.
      FLUX(2) = 0.D0
C
      RETURN
      END

C-----------------------------------------------------------------------------
C afsd_common.inc -- shared COMMON layout for the AFSD surface heat source.
C
C INCLUDE this in every routine that touches /AFSD_PATH/ or /AFSD_PARM/.
C Do NOT hand-copy the declarations: they were duplicated in three places
C once, drifted apart when a variable was added, and a REAL*8 slot was
C then silently reinterpreted as an INTEGER in one of the copies. A
C mismatched COMMON declaration does not fail to compile -- it corrupts
C data at runtime. Keep the single source of truth here.
C
C Must sit in the same directory as afsd_dflux.f so the compiler finds it.
C
C Ordering rule: all REAL*8 first, INTEGERs last. An integer ahead of
C REAL*8 arrays forces alignment padding, which gfortran warns about
C (-Walign-commons) and which is a portability hazard across compilers.
C
C /AFSD_PATH/  toolpath waypoints, loaded once by UEXTERNALDB
C   TPT  time [s]        TPX,TPY,TPZ  tool centre [mm]
C   TPS  deposition state at the waypoint (0 or 1)
C   NTP  number of waypoints actually read
C
C /AFSD_PARM/  process parameters, loaded once by UEXTERNALDB
C   PRES   normal pressure on the tool [MPa]
C   RSHLD  shoulder radius [mm]      RPIN  feed rod / pin radius [mm]
C   OMEGA  spindle speed [rad/s]     VTRAV traverse speed [mm/s]
C   FRIC   friction coefficient mu [-]
C   SLIP   slip factor delta [-]
C   SCALE  unit-consistency factor (1.0 for mm-tonne-s; see afsd_dflux.f)
C   ZTOL   tool-plane tolerance [mm], used when IRMODE=0
C   IRMODE 0 = in-plane radius + tool-plane + upward-face filtering
C          1 = full 3-D radius, no filtering (gamma_afsd.py equivalence)
C-----------------------------------------------------------------------------
      PARAMETER (MAXTP=40000)
      COMMON /AFSD_PATH/ TPT(MAXTP),TPX(MAXTP),TPY(MAXTP),
     1                    TPZ(MAXTP),TPS(MAXTP),NTP
      COMMON /AFSD_PARM/ PRES,RSHLD,RPIN,OMEGA,VTRAV,FRIC,SLIP,SCALE,
     1                    ZTOL,IRMODE

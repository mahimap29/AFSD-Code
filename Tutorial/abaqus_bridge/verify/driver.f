C-----------------------------------------------------------------------------
C driver.f -- standalone exerciser for AFSD_QVAL in ../afsd_dflux.f
C
C Compiled against the local ABA_PARAM.INC stub in this directory, so it
C runs without an Abaqus installation. It reads afsd_params.inp and
C toolpath.crs from the current working directory (exactly as
C UEXTERNALDB does), then prints the flux at a list of (x,y,z,t) probe
C points read from probes.txt.
C
C Output is one line per probe:  x y z t r q
C compare_flux.py reads it and checks it against the same expression
C evaluated from gamma_afsd.py's formulation.
C
C GETOUTDIR and XIT are Abaqus-supplied; minimal stubs are at the bottom.
C-----------------------------------------------------------------------------
      PROGRAM AFSD_DRIVER
      INCLUDE 'ABA_PARAM.INC'
C
      INCLUDE '../afsd_common.inc'
C
      DIMENSION TIME(2)
C
C     drive UEXTERNALDB's start-of-analysis branch (LOP=0) to load the
C     parameter file and toolpath into the COMMON blocks
      TIME(1) = 0.D0
      TIME(2) = 0.D0
      CALL UEXTERNALDB(0,0,TIME,0.D0,1,1)
C
      OPEN(UNIT=81,FILE='probes.txt',STATUS='OLD')
      OPEN(UNIT=82,FILE='flux_fortran.txt',STATUS='UNKNOWN')
C
   10 CONTINUE
      READ(81,*,END=20) XP,YP,ZP,TP,ZN
      CALL AFSD_QVAL(XP,YP,ZP,ZN,TP,QQ)
C     recompute r purely for reporting, same branch as AFSD_QVAL
      CALL AFSD_TOOLPOS(TP,XT,YT,ZT,ST)
      IF (IRMODE.EQ.1) THEN
        RR = SQRT((XP-XT)**2+(YP-YT)**2+(ZP-ZT)**2)
      ELSE
        RR = SQRT((XP-XT)**2+(YP-YT)**2)
      END IF
      WRITE(82,'(6(1X,E24.16))') XP,YP,ZP,TP,RR,QQ
      GO TO 10
   20 CONTINUE
C
      CLOSE(81)
      CLOSE(82)
      WRITE(6,*) 'driver: wrote flux_fortran.txt'
      END


C     ---- minimal stubs for Abaqus utility routines ----
      SUBROUTINE GETOUTDIR(OUTDIR,LENOUTDIR)
      CHARACTER*256 OUTDIR
      OUTDIR = '.'
      LENOUTDIR = 1
      RETURN
      END

      SUBROUTINE XIT
      STOP 'XIT'
      END

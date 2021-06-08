!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**************************************************************************%%
!%%** PROGRAM     PRESSURE FOR SQUARE-WELL FLUID FROM g(r)                 **%%
!%%**                                                                      **%%
!%%** AUTHOR      ALEXIS TORRES CARBAJAL                                   **%%
!%%** DATE        JUNE 4, 2021                                             **%%
!%%** VERSION     V00                                                      **%%
!%%** LICENSE     LGPL-V3                                                  **%%
!%%**                                                                      **%%
!%%** OBS         THE COMPUTATION DEPENDS ON THE FILE OF g(r)              **%%
!%%**             THE APPROACH IS BASED IN CAN. J. PHYS. 53, 5 (1975)      **%%
!%%**************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE PRESSVAR
 INTEGER, PARAMETER:: D     = KIND(1.0D0)   !PRECISION
 INTEGER, PARAMETER:: NDATX = 5000          !MAXIMUM NUMBER OF DATA

 REAL(D), PARAMETER:: ESTAR = 1.0D0         !ATTRACTIVE STRENGHT
 REAL(D), PARAMETER:: PI    = DACOS(-1.0D0) !PI NUMBER
 REAL(D), PARAMETER:: DBIN  = 0.01D0        !BIN WIDTH
 REAL(D), PARAMETER:: SIGMA = 1.0D0         !MOLECULES DIAMETER

 INTEGER:: I,NPART,NBIN,NSIG,NSIP,NSIM
 REAL(D):: RI(NDATX),GRI(NDATX)
 REAL(D):: RHO,BOXX,LAMBDA,LAMBDA3
 REAL(D):: L0,L1,L2,L3
 REAL(D):: TSTAR,PRESS
 REAL(D):: GSIG,GSIP,GSIM,PRESSC1
 REAL(D):: PRESSC2,Z1,Z2

 REAL(D):: FX,INTETW,XI0,XI1,XI2,Z
END MODULE PRESSVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM PRESSGR
 USE PRESSVAR

 OPEN(UNIT=10,FILE='Lambda2000Rho090T2000MDGr.dat',STATUS='OLD')
 REWIND(10) 

 NPART=1372                                  !NUMBER OF PARTICLES
 TSTAR=2.000D0                               !SYSTEM TEMPERATURE
 LAMBDA=2.00D0                               !PONTENTIAL RANGE
 RHO=0.90D0                                  !SYSTEM DENSITY
 BOXX=(DBLE(NPART)/RHO)**(1.0D0/3.0D0)       !SIMULATION BOX LENGHT
 NBIN=INT(0.5D0*BOXX/DBIN) - 1               !NUMBER OF DATA IN g(r) FILE
 NDAT=INT(LAMBDA/DBIN)                       !NUMBER OF DATA FOR INTEGRAL

! READ DATA FROM FILE
 DO I=1,NBIN
    READ(10,*)RI(I),GRI(I)
 ENDDO
 
 NSIG=101                                   !CONTACT VALUE OF g(r)
!SIGMA REPRESENTS THE VALUE AT EXTRAPOLATION
 L0=(SIGMA-RI(NSIG+1))*(SIGMA-RI(NSIG+2))*(SIGMA-RI(NSIG+3))
 L0=L0/((RI(NSIG) - RI(NSIG+1))*(RI(NSIG) - RI(NSIG+2))*(RI(NSIG) - RI(NSIG+3)))

 L1=(SIGMA-RI(NSIG))*(SIGMA-RI(NSIG+2))*(SIGMA-RI(NSIG+3))
 L1=L1/((RI(NSIG+1) - RI(NSIG))*(RI(NSIG+1) - RI(NSIG+2))*(RI(NSIG+1) - RI(NSIG+3)))

 L2=(SIGMA-RI(NSIG))*(SIGMA-RI(NSIG+1))*(SIGMA-RI(NSIG+3))
 L2=L2/((RI(NSIG+2) - RI(NSIG))*(RI(NSIG+2) - RI(NSIG+1))*(RI(NSIG+2) - RI(NSIG+3)))

 L3=(SIGMA-RI(NSIG))*(SIGMA-RI(NSIG+1))*(SIGMA-RI(NSIG+2))
 L3=L3/((RI(NSIG+3) - RI(NSIG))*(RI(NSIG+3) - RI(NSIG+1))*(RI(NSIG+3) - RI(NSIG+2)))
!EXTRAPOLATED g(r) AT CONTACT VALUE 
 GSIG=GRI(NSIG)*L0 + GRI(NSIG+1)*L1 + GRI(NSIG+2)*L2 + GRI(NSIG+3)*L3

!INTEGRATION SCHEME WITH THE COMPOSSED SIMPSON RULE
 FX=GSIG*(ESTAR/(LAMBDA - 1.0D0))*(1.0D0**3)
 XI0=FX
 FX=GRI(NDAT)*(ESTAR/(LAMBDA - 1.0D0))*RI(NDAT)**3
 XI0=XI0 + FX

 XI1=0.0D0
 XI2=0.0D0

 DO I=NSIG,NDAT-1
    FX=GRI(I)*(ESTAR/(LAMBDA-1.0D0))*RI(I)**3
    IF(MOD(I,2) .EQ. 0)THEN
      XI2=XI2 + FX
    ELSE
      XI1=XI1 + FX
    ENDIF
 ENDDO

!INTEGRAL
 INTETW=DBIN*( XI0 + 2.0D0*XI2 + 4.0D0*XI1)/3.0D0
!PRESSURE
 PRESS=RHO*TSTAR + (2.0D0*PI/3.0D0)*RHO*RHO*(TSTAR*GSIG - INTETW)
!COMPRESSIBILTY 
 Z=PRESS/(RHO*TSTAR)

!SHOW DATA ON DISPLAY 
 WRITE(6,101)RHO,PRESS,Z,GSIG,NDAT
 
 CLOSE(UNIT=10)
 101 FORMAT(4(F16.8),I10)

 STOP
END PROGRAM PRESSGR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


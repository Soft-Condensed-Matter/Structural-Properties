!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**************************************************************************%%
!%%** PROGRAM     <n> COORDINATION NUMBER FROM g(r)                        **%%
!%%**                                                                      **%%
!%%** AUTHOR      ALEXIS TORRES CARBAJAL                                   **%%
!%%** DATE        AGUST 2, 2021                                            **%%
!%%** VERSION     V00                                                      **%%
!%%** LICENSE     LGPL-V3                                                  **%%
!%%**                                                                      **%%
!%%** OBS         THE COMPUTATION DEPENDS ON THE FILE OF g(r)              **%%
!%%**             PROGRAM DEMANDS TEMPERATURE (TSTAR), NUMBER DENSITY (RHO)**%%
!%%**             AND PARTICLE NUMBER (NPART) AT WHICH g(r) WAS DETERMINED **%%
!%%**             SINCE FOR r < SIGMA g(r) IS EQUAL TO ZERO THE FIRST      **%%
!%%**             MINIMUM IS LOCATED AFTER THE FIRST MAXIMUM, THUS, THIS   **%%
!%%**             DETERMINATED FIRST AND TAKEN AS REFERENCE                **%%
!%%**************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE CORVAR
 INTEGER, PARAMETER:: D     = KIND(1.0D0)   !PRECISION
 INTEGER, PARAMETER:: NDATX = 5000          !MAXIMUM NUMBER OF DATA

 REAL(D), PARAMETER:: ESTAR = 1.0D0         !ATTRACTIVE STRENGHT
 REAL(D), PARAMETER:: PI    = DACOS(-1.0D0) !PI NUMBER
 REAL(D), PARAMETER:: DBIN  = 0.01D0        !BIN WIDTH
 REAL(D), PARAMETER:: SIGMA = 1.0D0         !MOLECULES DIAMETER

 INTEGER:: I,NPART,NBIN,NST
 REAL(D):: RI(NDATX),GRI(NDATX)
 REAL(D):: RHO,TSTAR,BOXX,NCOR
 REAL(D):: PHI,GMIN,GMAX

 REAL(D):: FX,INTE,XI0,XI1,XI2
END MODULE CORVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM NCOORDINATION
 USE CORVAR

 OPEN(UNIT=10,FILE='PhixxxTyyy/MDGr.dat',STATUS='OLD')
 REWIND(10) 

 NPART=1372                                  !NUMBER OF PARTICLES
 TSTAR=uuu                               !SYSTEM TEMPERATURE
 PHI=vvv                                     !
 RHO=6.0D0*PHI/PI                            !SYSTEM DENSITY
 BOXX=(DBLE(NPART)/RHO)**(1.0D0/3.0D0)       !SIMULATION BOX LENGHT
 NBIN=INT(0.5D0*BOXX/DBIN) - 1               !NUMBER OF DATA IN g(r) FILE

! READ DATA FROM FILE
 DO I=1,NBIN
    READ(10,*)RI(I),GRI(I)
 ENDDO
 
!THE FIRST MAXIMUM IS IDENTIFIED
 GMAX=0.0D0
 DO I=1,NBIN
    GMAX=MAX(GMAX,GRI(I))
    IF((GMAX .GT. GRI(I-1)) .AND. (GMAX .GT. GRI(I+1)) .AND. (RI(I) .GE. SIGMA))THEN
      NST=I
      EXIT
    ENDIF
 ENDDO

!NOW THE FIRS MIMIMUN IS IDENTIFIED TO SET INTEGRAL RANGE
 GMIN=100.0D0
 DO I=NST,NBIN
    GMIN=MIN(GMIN,GRI(I))
    IF((GMIN .LT. GRI(I-1)) .AND. (GMIN .LT. GRI(I+1)))THEN
       NDAT=I
       EXIT
    ENDIF
 ENDDO

!INTEGRATION SCHEME WITH THE COMPOSSED SIMPSON RULE
 FX=GRI(1)*RI(1)*RI(1)
 XI0=FX
 FX=GRI(NDAT)*RI(NDAT)*RI(NDAT)
 XI0=XI0 + FX

 XI1=0.0D0
 XI2=0.0D0

 DO I=2,NDAT-1
    FX=GRI(I)*RI(I)*RI(I)
    IF(MOD(I,2) .EQ. 0)THEN
      XI2=XI2 + FX  
    ELSE
      XI1=XI1 + FX  
    ENDIF
 ENDDO

!COORDINATION NUMBER INTEGRAL
 INTE=DBIN*( XI0 + 2.0D0*XI2 + 4.0D0*XI1)/3.0D0
 NCOR=4.0D0*PI*RHO*INTE

 WRITE(6,*)TSTAR,RHO,NCOR,GMAX,GMIN,NST,NDAT
 STOP
END PROGRAM NCOORDINATION

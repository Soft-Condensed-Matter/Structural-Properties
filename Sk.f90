!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**************************************************************************%%
!%%**  PROGRAM    STRUCTURE FACTOR                                         **%%
!%%**  AUTHOR     ALEXIS TORRES CARBAJAL                                   **%%
!%%**  LICENSE    LGPL-V3                                                  **%%
!%%**  DATE       OCTOBER 27, 2020                                         **%%
!%%**                                                                      **%%
!%%**  OBS        THIS PROGRAM PERFORMS THE FOURIER TRANSFORM IN 3D OF     **%%
!%%**             THE RADIAL DISTRIBUTION FUNCTION WHICH IS PROVIDED IN    **%%
!%%**             THE EXTERNAL FILE Gr.dat                                 **%%
!%%**                                                                      **%%
!%%**             THE NUMBER OF DATA IN Gr.dat FILE IS COMPUTE USING THE   **%%
!%%**             SYSTEM FRACTION VOLUME AND THE NUMBER OF PARTICLES IN    **%%
!%%**             THE SIMULATION. HOWEVER, SUCH DATA CAN BE OMMITED AND    **%%
!%%**             DIRECTLY SPECIFIES THE NDAT VALUE                        **%%
!%%**                                                                      **%%
!%%**             THE MINIMUM WAVE NUMBER IS SET TAKEN INTO ACCOUNT PBC    **%%
!%%**             IN THE SIMULATION AND THE INTEGRAL IS COMPUTED WITH THE  **%%
!%%**             SIMPLEST TRAPEZOIDAL RULE                                **%%
!%%**                                                                      **%%
!%%**             AT LOW K VALUES THE SIN FUNCTION CAN BE EXPANDED IN A    **%%
!%%**             TAYLOR SERIES TO GAIN SOME PRECISION. FEATURE NOT        **%%
!%%**             IMPLEMENTED                                              **%%
!%%**                                                                      **%%
!%%**************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE SKVAR
 INTEGER, PARAMETER:: D    = KIND(1.0D0)    !PRECISION
 INTEGER, PARAMETER:: NX   = 5000           !MAXIMUM NUMBER OF DATA
 REAL(D), PARAMETER:: PI   = DACOS(-1.0D0)  !PI NUMBER
 REAL(D), PARAMETER:: DBIN = 0.01D0         !BIN WIDTH
 REAL(D), PARAMETER:: KMAX = 35.0D0         !MAXIMUM WAVE NUMBER

 INTEGER:: NDAT,I,J,NK,NPART
 REAL(D):: R(NX),GR(NX),PHI
 REAL(D):: SUMMX,FX,K,DK,X
 REAL(D):: BOXX,FXA,FXB
 REAL(D):: INTE,RHO,SK
END MODULE SKVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM STRUCTURE_FACTOR
 USE SKVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='Gr.dat')                !g(r) FILE
 OPEN(UNIT=11,FILE='Sk.dat')                !S(k) FILE

 REWIND(10)

 PHI=0.47123890                             !PACKING FRACTION
 RHO=6.0D0*PHI/PI                           !NUMBER DENSITY
 NPART=1372                                 !NUMBER OF PARTICLES IN SIMULATION
 BOXX=(DBLE(NPART)/RHO)**(1.0D0/3.0D0)      !SIMULATION BOX LENGHT

 NDAT=INT(0.5D0*BOXX/DBIN)-2                !NUMBER OF DATA IN g(r) FILE
 DK=PI/BOXX                                 !WAVE NUMBER INCREMENT
 K=DK                                       !INITIAL WAVE NUMBER              
 NK=INT(KMAX/DK) + 1                        !NUMBER OF WAVE NUMBERS

 DO I=1,NDAT                                !g(r) DATA FROM FILE
    READ(10,*)R(I),GR(I)
 ENDDO

 DO J=1,NK
    SUMMX=0.0D0
    DO I=2,NDAT-1
       X=R(I)*K
       FX=(GR(I) - 1.0D0)*R(I)*DSIN(X)
       SUMMX=SUMMX + FX
    ENDDO

    X=R(1)*K
    FXA=(GR(1) - 1.0D0)*R(1)*DSIN(X)
    X=R(NDAT)*K
    FXB=(GR(NDAT) - 1.0D0)*R(NDAT)*DSIN(X)

    INTE=DBIN*0.5D0*(FXA + 2.0D0*SUMMX + FXB)
    SK=1.0D0 + 4.0D0*PI*RHO*INTE/K
    WRITE(11,*)K,SK
    K=K + DK
 ENDDO
 
 CLOSE(UNIT=10)
 CLOSE(UNIT=11)

 STOP
END PROGRAM STRUCTURE_FACTOR

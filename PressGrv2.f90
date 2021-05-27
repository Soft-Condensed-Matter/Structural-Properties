!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**************************************************************************%%
!%%** PROGRAM     PRESSURE FOR SQUARE-WELL FLUID FROM g(r)                 **%%
!%%**                                                                      **%%
!%%** AUTHOR      ALEXIS TORRES CARBAJAL                                   **%%
!%%** DATE        FEBRUARY 14, 2019                                        **%%
!%%** VERSION     V00                                                      **%%
!%%** LICENSE     LGPL-V3                                                  **%%
!%%**                                                                      **%%
!%%** OBS         THE COMPUTATION DEPENDS ON THE FILE OF g(r)              **%%
!%%**             THE APPROACH IS BASED IN J. CHEM. PHYS. 67, 5308 (1977)  **%%
!%%**************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE PRESSVAR
 INTEGER, PARAMETER:: D     = KIND(1.0D0)
 INTEGER, PARAMETER:: NDATX = 4000

 REAL(D), PARAMETER:: PI    = DACOS(-1.0D0)
 REAL(D), PARAMETER:: DBIN  = 0.01D0
 
 INTEGER:: I,NPART,NBIN,NSIG,NSIP,NSIM
 REAL(D):: RI(NDATX),GRI(NDATX)
 REAL(D):: RHO,BOXX,LAMBDA,LAMBDA3
 REAL(D):: L0,L1,L2,L3
 REAL(D):: TSTAR,PRESS
 REAL(D):: GSIG,GSIP,GSIM,PRESSC1
 REAL(D):: PRESSC2,Z1,Z2
END MODULE PRESSVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM PRESSGR
 USE PRESSVAR

 OPEN(UNIT=10,FILE='GrMD.dat',STATUS='OLD')
 REWIND(10) 

 NPART=2048                                  !NUMBER OF PARTICLES
 TSTAR=0.8370D0                              !SYSTEM TEMPERATURE
 LAMBDA=1.50D0                               !PONTENTIAL RANGE
 RHO=0.yyD0                                  !SYSTEM DENSITY
 BOXX=(DBLE(NPART)/RHO)**(1.0D0/3.0D0)      !SIMULATION BOX LENGHT
 NBIN=INT(0.5D0*BOXX/DBIN) - 1              !NBUMBER OF DATA IN g(r) FILE
 LAMBDA3=LAMBDA**3

 DO I=1,NBIN
    READ(10,*)RI(I),GRI(I)
 ENDDO
 
 NSIG=101                                   !CONTACT VALUE OF g(r)
 NSIM=124                                   !DISCONTINUITY VALUE - OF g(r)
 NSIP=127                                   !DISCONTINUITY VALUE + OF g(r)

!1.0DO REPRESENTS THE VALUE AT EXTRAPOLATION
 L0=(1.0D0-RI(NSIG+1))*(1.0D0-RI(NSIG+2))*(1.0D0-RI(NSIG+3))
 L0=L0/((RI(NSIG) - RI(NSIG+1))*(RI(NSIG) - RI(NSIG+2))*(RI(NSIG) - RI(NSIG+3)))

 L1=(1.0D0-RI(NSIG))*(1.0D0-RI(NSIG+2))*(1.0D0-RI(NSIG+3))
 L1=L1/((RI(NSIG+1) - RI(NSIG))*(RI(NSIG+1) - RI(NSIG+2))*(RI(NSIG+1) - RI(NSIG+3)))

 L2=(1.0D0-RI(NSIG))*(1.0D0-RI(NSIG+1))*(1.0D0-RI(NSIG+3))
 L2=L2/((RI(NSIG+2) - RI(NSIG))*(RI(NSIG+2) - RI(NSIG+1))*(RI(NSIG+2) - RI(NSIG+3)))

 L3=(1.0D0-RI(NSIG))*(1.0D0-RI(NSIG+1))*(1.0D0-RI(NSIG+2))
 L3=L3/((RI(NSIG+3) - RI(NSIG))*(RI(NSIG+3) - RI(NSIG+1))*(RI(NSIG+3) - RI(NSIG+2)))

 GSIG=GRI(NSIG)*L0 + GRI(NSIG+1)*L1 + GRI(NSIG+2)*L2 + GRI(NSIG+3)*L3


 L0=(LAMBDA-RI(NSIP+1))*(LAMBDA-RI(NSIP+2))*(LAMBDA-RI(NSIP+3))
 L0=L0/((RI(NSIP) - RI(NSIP+1))*(RI(NSIP) - RI(NSIP+2))*(RI(NSIP) - RI(NSIP+3)))

 L1=(LAMBDA-RI(NSIP))*(LAMBDA-RI(NSIP+2))*(LAMBDA-RI(NSIP+3))
 L1=L1/((RI(NSIP+1) - RI(NSIP))*(RI(NSIP+1) - RI(NSIP+2))*(RI(NSIP+1) - RI(NSIP+3)))

 L2=(LAMBDA-RI(NSIP))*(LAMBDA-RI(NSIP+1))*(LAMBDA-RI(NSIP+3))
 L2=L2/((RI(NSIP+2) - RI(NSIP))*(RI(NSIP+2) - RI(NSIP+1))*(RI(NSIP+2) - RI(NSIP+3)))

 L3=(LAMBDA-RI(NSIP))*(LAMBDA-RI(NSIP+1))*(LAMBDA-RI(NSIP+2))
 L3=L3/((RI(NSIP+3) - RI(NSIP))*(RI(NSIP+3) - RI(NSIP+1))*(RI(NSIP+3) - RI(NSIP+2)))

 GSIP=GRI(NSIP)*L0 + GRI(NSIP+1)*L1 + GRI(NSIP+2)*L2 + GRI(NSIP+3)*L3


 L0=(LAMBDA-RI(NSIM-1))*(LAMBDA-RI(NSIM-2))*(LAMBDA-RI(NSIM-3))
 L0=L0/((RI(NSIM) - RI(NSIM-1))*(RI(NSIM) - RI(NSIM-2))*(RI(NSIM) - RI(NSIM-3)))

 L1=(LAMBDA-RI(NSIM))*(LAMBDA-RI(NSIM-2))*(LAMBDA-RI(NSIM-3))
 L1=L1/((RI(NSIM-1) - RI(NSIM))*(RI(NSIM-1) - RI(NSIM-2))*(RI(NSIM-1) - RI(NSIM-3)))

 L2=(LAMBDA-RI(NSIM))*(LAMBDA-RI(NSIM-1))*(LAMBDA-RI(NSIM-3))
 L2=L2/((RI(NSIM-2) - RI(NSIM))*(RI(NSIM-2) - RI(NSIM-1))*(RI(NSIM-2) - RI(NSIM-3)))

 L3=(LAMBDA-RI(NSIM))*(LAMBDA-RI(NSIM-1))*(LAMBDA-RI(NSIM-2))
 L3=L3/((RI(NSIM-3) - RI(NSIM))*(RI(NSIM-3) - RI(NSIM-1))*(RI(NSIM-3) - RI(NSIM-2)))

 GSIM=GRI(NSIM)*L0 + GRI(NSIM-1)*L1 + GRI(NSIM-2)*L2 + GRI(NSIM-3)*L3

 PRESS=RHO*TSTAR + 2.0D0*PI*RHO*RHO*TSTAR*GRI(NSIG)/3.0D0
 PRESS=PRESS - PI*RHO*RHO*(GRI(NSIM) + GRI(NSIP))*LAMBDA3/3.0D0

 !HENDERSON
 PRESSC1=RHO*TSTAR + 2.0D0*PI*RHO*RHO*TSTAR*GSIG/3.0D0
 PRESSC1=PRESSC1 - 2.0D0*PI*RHO*RHO*TSTAR*(GSIM - GSIP)*LAMBDA3/3.0D0
 Z1=PRESSC1/(RHO*TSTAR)

 !TAGO
 PRESSC2=RHO*TSTAR + 2.0D0*PI*RHO*RHO*TSTAR*GSIG/3.0D0
 PRESSC2=PRESSC2 - RHO*RHO*PI*LAMBDA3*(GSIM + GSIP)/3.0D0
 Z2=PRESSC2/(RHO*TSTAR)
 
 PRINT*,RHO,PRESSC1,PRESSC2,Z1,Z2,GSIG,GSIM,GSIP

! PRINT*,'DENSITY     =',RHO
! PRINT*,'PARTICLES   =',NPART
! PRINT*,'LAMBDA      =',LAMBDA
! PRINT*,'TEMPERATURE =',TSTAR
! PRINT*,'PRESSURE    =',PRESS,PRESSC

 CLOSE(UNIT=10)

 STOP
END PROGRAM PRESSGR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


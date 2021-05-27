!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**************************************************************************%%
!%%**  PROGRAM     ANALYTICAL g(r) FOR LJ 12-6                             **%%
!%%**  AUTHOR      ALEXIS TORRES                                           **%%
!%%**  DATE        FEBREAURY 25, 2019                                      **%%
!%%**  VERSION     BETA 00                                                 **%%
!%%**  LICENSE     LGPL-V3                                                 **%%
!%%**                                                                      **%%
!%%**  OBS         THIS PROGRAM COMPUTES THE RADIAL DISTRIBUTION FUNCTION  **%%
!%%**              FOR THE LENNARD-JONES 12-6 POTENTIAL FOR THERMODYNAMIC  **%%
!%%**              STATES IN 0.35 < \rho < 1.1 AND 0.5 < T < 5.1           **%%
!%%**              THE PROGRAM IS BASED IN CHEM. PHYS. 310, 11 (2015)      **%%
!%%**                                                                      **%%
!%%**************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE LJGRVAR
 INTEGER, PARAMETER:: DP   = KIND(1.0D0) 
 REAL(DP), PARAMETER:: RI   = 0.0D0
 REAL(DP), PARAMETER:: RF   = 4.5D0
 INTEGER, PARAMETER:: NIT  = 100

 CHARACTER*50:: NAMEFILE 
 INTEGER:: I,J,O,IJ
 REAL(DP):: RHO,TEMP,R,GR,DR,RN
 REAL(DP):: IRHO,ITEMP
 REAL(DP):: AA,BB,CC,DD,GG,HH,KK,LL,MM,NN,SS
END MODULE LJGRVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM LJGR
 USE LJGRVAR

 OPEN(UNIT=11,FILE='IDfile.dat')
 
 IRHO=0.35D0
 ITEMP=0.50D0
 IJ=0
 
 DO I=1,101

    DO J=1,101
       IJ=IJ+1
       13 CONTINUE
       CALL RANDOM_NUMBER(RN)
       RHO=IRHO + RN
       IF(RHO .GT. 1.1D0)GOTO 13

       14 CONTINUE 
       CALL RANDOM_NUMBER(RN)
       TEMP=ITEMP + RN*10
       IF(TEMP .GT. 4.5D0)GOTO 14

       WRITE(NAMEFILE,'("Rand",I5.5,"T",I3.3,"Rho",I3.3,".dat")')IJ,J,I
       WRITE(11,*)J,TEMP,I,RHO,NAMEFILE
       OPEN(UNIT=10,FILE=NAMEFILE)

       CALL CONSTANTS                             !ADJUST CONSTANTS
       DR=(RF-RI)/DBLE(NIT)                       !BIN WIDTH
       R=RI                                       !INITIAL VALUE

       DO O=1,NIT
          IF(R .LE. 1.0D0)THEN
            GR=SS*EXP(-(MM*R + NN)**4)
          ELSE
            GR=1.0D0 + EXP(-(AA*R + BB))*DSIN(CC*R + DD)/R**2
            GR=GR + EXP (-(GG*R + HH))*DCOS(KK*R + LL)/R**2
         ENDIF

           WRITE(10,*)R,GR
           R=R + DR
       ENDDO
 
      CLOSE(UNIT=10)
    ENDDO
 ENDDO

 CLOSE(UNIT=11)

 STOP
END PROGRAM LJGR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SUBROUTINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE CONSTANTS 
 USE LJGRVAR

 CALL A
 CALL B
 CALL C
 CALL D
 CALL G
 CALL H
 CALL K
 CALL L
 CALL M
 CALL N
 CALL S

 RETURN
END SUBROUTINE CONSTANTS 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE A
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1A = 9.24792
 REAL(DP), PARAMETER:: Q2A = -2.64281
 REAL(DP), PARAMETER:: Q3A = 0.133386
 REAL(DP), PARAMETER:: Q4A = -1.35932
 REAL(DP), PARAMETER:: Q5A = 1.25338
 REAL(DP), PARAMETER:: Q6A = 0.45602
 REAL(DP), PARAMETER:: Q7A = -0.326422
 REAL(DP), PARAMETER:: Q8A = 0.045708
 REAL(DP), PARAMETER:: Q9A = -0.0287681

 AA=Q1A + Q2A*EXP(-Q3A*TEMP) + Q4A*EXP(-Q5A*TEMP)
 AA=AA + Q6A/RHO + Q7A/RHO**2 + Q8A*EXP(-Q3A*TEMP)/RHO**3
 AA=AA + Q9A*EXP(-Q5A*TEMP)/RHO**4
 
 RETURN
END SUBROUTINE A
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE B
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1B = -8.33289
 REAL(DP), PARAMETER:: Q2B = 2.1714
 REAL(DP), PARAMETER:: Q3B = 1.00063
 
 BB=Q1B + Q2B*EXP(-Q3B*RHO)

 RETURN
END SUBROUTINE B
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE C
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1C = -0.0677912
 REAL(DP), PARAMETER:: Q2C = -1.39505
 REAL(DP), PARAMETER:: Q3C = 0.512625
 REAL(DP), PARAMETER:: Q4C = 36.9323
 REAL(DP), PARAMETER:: Q5C = -36.8061
 REAL(DP), PARAMETER:: Q6C = 21.7353
 REAL(DP), PARAMETER:: Q7C = -7.76671
 REAL(DP), PARAMETER:: Q8C = 1.36342
 
 CC=Q1C + Q2C*EXP(-Q3C*TEMP) + Q4C*RHO + Q5C*RHO**2
 CC=CC + Q6C*RHO**3 + Q7C*RHO**4 + Q8C*RHO**5

 RETURN
END SUBROUTINE C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE D
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1D = -26.1615
 REAL(DP), PARAMETER:: Q2D = 27.4846
 REAL(DP), PARAMETER:: Q3D = 1.68124
 REAL(DP), PARAMETER:: Q4D = 6.74296

 DD=Q1D + Q2D*EXP(-Q3D*RHO) + Q4D*RHO
 
 RETURN
END SUBROUTINE D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE G
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1G = 0.663161
 REAL(DP), PARAMETER:: Q2G = -0.243089
 REAL(DP), PARAMETER:: Q3G = 1.24749
 REAL(DP), PARAMETER:: Q4G = -2.059
 REAL(DP), PARAMETER:: Q5G = 0.04261
 REAL(DP), PARAMETER:: Q6G = 1.65041
 REAL(DP), PARAMETER:: Q7G = -0.343652
 REAL(DP), PARAMETER:: Q8G = -0.037698
 REAL(DP), PARAMETER:: Q9G = 0.008899

 GG=Q1G + Q2G*EXP(-Q3G*TEMP) + Q4G*EXP(-Q5G*TEMP)
 GG=GG + Q6G/RHO + Q7G/RHO**2 + Q8G*EXP(-Q3G*TEMP)/RHO**3
 GG=GG + Q9G*EXP(-Q5G*TEMP)/RHO**4

 
 RETURN
END SUBROUTINE G
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE H
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1H = 0.0325239
 REAL(DP), PARAMETER:: Q2H = -1.28792
 REAL(DP), PARAMETER:: Q3H = 2.5487
 
 HH=Q1H + Q2H*EXP(-Q3H*RHO)

 RETURN
END SUBROUTINE H
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE K
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1K = 16.4821
 REAL(DP), PARAMETER:: Q2K = -0.300612
 REAL(DP), PARAMETER:: Q3K = 0.0937844
 REAL(DP), PARAMETER:: Q4K = -61.744
 REAL(DP), PARAMETER:: Q5K = 145.285
 REAL(DP), PARAMETER:: Q6K = -168.087
 REAL(DP), PARAMETER:: Q7K = 98.2181
 REAL(DP), PARAMETER:: Q8K = -23.0583
 
 KK=Q1K + Q2K*EXP(-Q3K*TEMP) + Q4K*RHO + Q5K*RHO**2
 KK=KK + Q6K*RHO**3 + Q7K*RHO**4 + Q8K*RHO**5
 
 RETURN
END SUBROUTINE K
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE L
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1L = -6.7293
 REAL(DP), PARAMETER:: Q2L = -59.5002
 REAL(DP), PARAMETER:: Q3L = 10.2466
 REAL(DP), PARAMETER:: Q4L = -0.43596

 LL=Q1L + Q2L*EXP(-Q3L*RHO) + Q4L*RHO
 
 RETURN
END SUBROUTINE L
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE M
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1M = -5.668
 REAL(DP), PARAMETER:: Q2M = -3.62671
 REAL(DP), PARAMETER:: Q3M = 0.680654
 REAL(DP), PARAMETER:: Q4M = 0.294481
 REAL(DP), PARAMETER:: Q5M = 0.186395
 REAL(DP), PARAMETER:: Q6M = -0.286954

 MM=Q1M + Q2M*EXP(-Q3M*TEMP) + Q4M/TEMP 
 MM=MM + Q5M*RHO + Q6M*RHO**2

 RETURN
END SUBROUTINE M
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE N
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1N = 6.01325
 REAL(DP), PARAMETER:: Q2N = 3.84098
 REAL(DP), PARAMETER:: Q3N = 0.60793

 NN=Q1N + Q2N*EXP(-Q3N*TEMP)
 
 RETURN
END SUBROUTINE N
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE S
 USE LJGRVAR
 REAL(DP), PARAMETER:: Q1S = 1.25225
 REAL(DP), PARAMETER:: Q2S = -1.0179
 REAL(DP), PARAMETER:: Q3S = 0.358564
 REAL(DP), PARAMETER:: Q4S = -0.18533
 REAL(DP), PARAMETER:: Q5S = 0.0482119
 REAL(DP), PARAMETER:: Q6S = 1.27592
 REAL(DP), PARAMETER:: Q7S = -1.78785
 REAL(DP), PARAMETER:: Q8S = 0.634741
 
 SS=Q1S + Q2S*RHO + Q3S/TEMP + Q4S/TEMP**2 + Q5S/TEMP**3
 SS=SS/(Q6S + Q7S*RHO + Q8S*RHO**2)

 RETURN
END SUBROUTINE S






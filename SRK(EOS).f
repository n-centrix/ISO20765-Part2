 
	  SUBROUTINE SRK(PRESS, TEMP, CFI1, CFI2, XI) 
!*******************************************************************************
!*******************SOAVE-REDLICH-KWONG SRK EOS PROCEDURE***********************
!******************************************************************************* 
!REFERENCE: API TECHNICAL DATA BOOK 6th Ed. - PROC. 8D1.1 
!CALCULATES ROOTS OF THE SRK CUBIC EOS WHILE SOVLING FOR COMPRESSIBILITY
!IDENTIFY THE TYPE OF ROOTS (IMAGINARY, REAL OR TRIPLE ROOT) PASSING TO: MTYPE
!RETURNS THE CORRESPONDING REDUCED DENSITIES (COMMON)
!RETURNS THE FUGACITY COEFFICIENTS OF THE OUTER ROOTS
!OUTPUT: CFI1, CFI2
!OUTPUT: DLV, DLL, DL99 (AS COMMON)
!OUTPUT: MTYPE (AS COMMON)	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, XI(21), TEMP, RT1,RT2,RT3
      REAL*8 ACT(21), PRESSC(21), CFI1(21), CFI2(21)
      REAL*8 RHOC(21),TEMPC(21), SRKALPH(21)
      REAL*8 BIP(21,21), PCP1(21),PCP2(21)	
      REAL*8 PRESSREF, TEMPREF, R, RSTAR, TA2, DLV, DLL
      REAL*8 ACFT, BCFT, CCFT, DCFT, TA, TB, TA1, TB1
      REAL*8 INVRHOR, TEMPR, TAU, Z1, Z2, Z99, DL99
                
	  INTEGER i,j,MTYPE
      COMMON /PARAM5/ RHOC, TEMPC
      COMMON /PARAM4/ ACT, PRESSC	
      COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
      COMMON /PARAM8/ BIP 
      COMMON /PARAM9/ DLV, DLL, DL99, MTYPE

      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)

       DLL=0D0
       DLV=0D0
       DL99=0D0
      
       TA=0D0      !TA = Alpha x a
       TB=0D0      !TB = b
       TA1=0D0     !TA1= A =alpha x a x P/ (RT)^2
       TB1=0D0     !TB1= B =b x P/ (RT)
       TA2=0D0
                  
       DO i=1,21
        PCP1(i) = 0D0
        PCP2(i) = 0D0
        SRKALPH(i) = 0D0
       ENDDO    
       
       DO i = 1,21
        DO j= 1,21
         BIP(j,i)=BIP(i,j)
        ENDDO 
       ENDDO  
            
!ALPHA FUNCTION    
       DO i=1,21
        PCP1(i) = 0.48508D0+1.55171D0*ACT(i)-0.15613D0*ACT(i)**2D0
       ENDDO
       
       PCP1(18)=  1.243997D0
       PCP2(18)= -0.201789D0  !USE ONLY IN LOW H2O CONCENTRATION IN HC's

       DO i=1,21      
       SRKALPH(i)=(1D0+PCP1(i)* (1D0- DSQRT(TEMP/TEMPC(i))) + PCP2(i)*  
     &  (1D0- DSQRT(TEMP/TEMPC(i))) / DSQRT(TEMP/TEMPC(i))  ) **2D0     
       ENDDO 
!--------------------           
       DO i=1,21 
        DO j=1,21      
        TA = TA+XI(i)*XI(j)*(1D0-BIP(i,j))* DSQRT( 
     &  SRKALPH(i)*SRKALPH(j)*
     &  (  0.42747D0*(R*TEMPC(i))**2D0/PRESSC(i)  ) * 
     &  (  0.42747D0*(R*TEMPC(j))**2D0/PRESSC(j)  )   ) 
        ENDDO   
       ENDDO             
       TA1 = TA * PRESS / (R*TEMP)**2D0     
              
       DO i=1,21       
        TB = TB+XI(i)*0.08664D0*R*TEMPC(i) / PRESSC(i)      
       ENDDO                      
       TB1 = TB * PRESS / (R*TEMP)    
         
!RE-USE VECTOR PCP1(i) FOR DIFFERENT PURPOSE
       DO i=1,21
        PCP1(i) = 0D0
       ENDDO  
                               
       DO i=1,21
        DO j=1,21       
        PCP1(i)=PCP1(i)+2D0*XI(j)*(1D0-BIP(i,j))* DSQRT( 
     &  SRKALPH(i)*SRKALPH(j)*
     &  (  0.42747D0*(R*TEMPC(i))**2D0/PRESSC(i)  ) * 
     &  (  0.42747D0*(R*TEMPC(j))**2D0/PRESSC(j)  )   )        
        ENDDO     
       ENDDO 
       
       DO i=1,21     
        TA2 = TA2 + XI(i)*PCP1(i)*PRESS / (R*TEMP)**2D0     
       ENDDO        
                                                           
!CUBIC EQUATION COEFFICIENTS AxZ^3+ BxZ^2+CxZ+D=0
       ACFT= 1D0
       BCFT= -1D0
       CCFT= TA1-TB1-TB1**2D0
       DCFT= -TA1*TB1
       
       CALL CUBICSOL(ACFT, BCFT, CCFT, DCFT ,RT1,RT2,RT3, MTYPE) 
!ASSIGN OUTER ROOTS
       Z1=MAX(RT1, RT2, RT3)
       Z2=MIN(RT1, RT2, RT3)
       DLV=(PRESS/1D-3)/( Z1* R*TEMPR/(TAU*INVRHOR))      
       DLL=(PRESS/1D-3)/( Z2* R*TEMPR/(TAU*INVRHOR))  
       
!ASSIGN INTERMEDIATE ROOT - TYPICALLY NOT USED
       Z99=RT1+RT2+RT3-MIN(RT1, RT2, RT3)-MAX(RT1, RT2, RT3)       
       DL99=(PRESS/1D-3)/(Z99* R*TEMPR/(TAU*INVRHOR))  
      
!CFI1 fugacity coefficient Fi for 1st outer root
!CFI2 fugacity coefficient Fi for 2nd outer root 
!Note: for fugacity divide CFI1 by (XI(i) * PRESS)  
     
       DO i=1,21
         CFI1(i)=0D0
         CFI2(i)=0D0
       ENDDO
       
      DO i=1,21
       CFI1(i)=DEXP( (Z1-1D0)*
     &  (0.08664D0*R*TEMPC(i)/PRESSC(i))/TB-DLOG(Z1-TB1)-
     &  TA1*( (PCP1(i)/TA)-(0.08664D0*R*TEMPC(i)/PRESSC(i))/TB )
     &  *DLOG(1D0+TB1/Z1)/TB1 ) 
       ENDDO       

      DO i=2,21
       CFI2(i)=DEXP( (Z2-1D0)*
     &  (0.08664D0*R*TEMPC(i)/PRESSC(i))/TB-DLOG(Z2-TB1)-
     &  TA1*( (PCP1(i)/TA)-(0.08664D0*R*TEMPC(i)/PRESSC(i))/TB )
     &  *DLOG(1D0+TB1/Z2)/TB1 )
       ENDDO         
                    
       END
!******************************************************************************* 	

!*******************************************************************************
       SUBROUTINE CUBICSOL(ACFT, BCFT, CCFT, DCFT , RT1,RT2,RT3, MTYPE) 
!GENERAL PROCEDURE FOR CUBIC EQUATION OF THE FORM A x Z^3+ B x Z^2 + C x Z + D=0 
!RETURNS THE ROOTS AND TYPE (IMAGINERY, REAL OR TRIPLE ROOT)
!OUTPUT: RT1, RT2, RT3
!OUTPUT: MTYPE     
       IMPLICIT DOUBLE PRECISION (A-Z)
       REAL*8 B(3)
       INTEGER MTYPE
         
       B(1)=BCFT/ACFT
       B10V3=B(1)/3.0D0
       B(2)=CCFT/ACFT
       B(3)=DCFT/ACFT
       ALF=B(2)-B(1)*B10V3
       BBT=2.0D0*B10V3**3D0-B(2)*B10V3+B(3)
       BETOV=BBT/2.0D0
       ALFOV=ALF/3.0D0
       CUAOV=ALFOV**3D0
       SQBOV=BETOV**2D0
       DEL=SQBOV+CUAOV
       IF (DEL) 90,10,40
10     MTYPE = 0
! THREE EQUAL ROOTS -> CRITICAL ISOTHERMS  (T = Tc FOR P = Pc)
       GAM=SQRT(-ALFOV)
       IF (BBT) 30,30,20
20     RT1 = -2.0D0*GAM-B10V3
       RT2 = GAM-B10V3
       RT3 = RT2
       GOTO 130
30     RT1 = 2.0D0*GAM-B10V3
       RT2 = -GAM-B10V3
       RT3 = RT2
       GOTO 130
40     MTYPE = 1
! ONE REAL ROOT & TWO IMAGINARY CONJUGATE ROOTS -> SUPERCRITICAL ISOTHERMS (T>Tc) 
       EPS=DSQRT(DEL)
       TAU1=-BETOV
       RCU=TAU1+EPS
       SCU=TAU1-EPS
       SIR=1.0D0
       SIS=1.0D0
       IF (RCU) 50,60,60
50     SIR=-1.0D0
60     IF (SCU) 70,80,80
70     SIS=-1.0D0
80     R1=SIR*(SIR*RCU)**(1D0/3D0)
       S1=SIS*(SIS*SCU)**(1D0/3D0)
       RT1=R1+S1-B10V3
       RT2=-(R1+S1)/2.0D0-B10V3
       RT3=3D0**0.5D0*0.5D0*(R1-S1)
       GOTO 130
90     MTYPE = -1
! THREE DISSIMILAR AND REAL ROOTS  ->  SUBCRITICAL ISOTHERMS (T<Tc) 
       QUOT=SQBOV/CUAOV
       RCOT=SQRT(-QUOT)
       IF (BBT) 110,100,100
100    PEI=(1.570796326791001D0+DATAN(RCOT/SQRT(1.0D0-RCOT**2D0)))/3.0D0
       GOTO 120
110    PEI=DATAN(SQRT(1.0D0-RCOT**2D0)/RCOT)/3.0D0
120    FACT=2.0*DSQRT(-ALFOV)
       RT1=FACT*DCOS(PEI)-B10V3
       PEI=PEI+ ACOS(-0.5D0)
       RT2=FACT*DCOS(PEI)-B10V3
       PEI=PEI+ ACOS(-0.5D0)
       RT3=FACT*DCOS(PEI)-B10V3
130    CONTINUE
       IF (MTYPE .EQ. 1) THEN
       RT2 = -99.99D0
       RT3 = -99.99D0
       ENDIF
       RETURN
       END 
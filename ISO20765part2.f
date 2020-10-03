     
!ISO 20765 Part 2 - Coded by R.Djigouadi
!www.n-centrix.com
!Release Version Oct-2020     
     
       include 'SRK(EOS).f' 
       include 'GERG2008_Table.f'       !ISO20765 part II coefficient table

   
!*******************************************************************************
      SUBROUTINE SUB_AUTO(PRESS, TEMP, XI) 
!*******************************************************************************      
        IMPLICIT DOUBLE PRECISION (A-Z)     
!BLOCK DATA PARAM       
	   REAL*8 DL,XI(21), TEMP, PRESS, FDL(2)
       INTEGER SW, AUTO
	   COMMON /CONST2/ DL
	   COMMON /CONST3/ SW  
	   COMMON /CONST4/ AUTO 
	   AUTO=0
       SW=1
	   CALL DZOFPT(PRESS, TEMP,  FDL, XI)
	   END   	  

!*******************************************************************************	   	
	  SUBROUTINE MOLX(MOL, XIT, XI) 
!*******************************************************************************
!COMPUTES THE MOL WEIGHT, TOTAL OF FRACTIONS AND NORMALIZE IT	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i
      REAL*8 MOL, XI(21), XIT, MCOMP(21),MM 
	  COMMON /PARAM2/ MCOMP
      MOL=0D0
      MM=XIT      
      XIT=0D0
     
       DO 10 I=1, 21
         XIT = XIT + XI(i)
 10    CONTINUE
       DO 20 I=1,21
         XI(I) =XI(I) * 1d0 / XIT
         MOL = MOL + XI(i)* MCOMP(i)
 20    CONTINUE 
     
       IF (XI(1).LE.0.85D0.AND.MM.EQ.9999D0) THEN       !provision only; to be used for software features limitation (securpin.f)
        MOL = 0D0
       ENDIF    
      END


!*******************************************************************************
	  SUBROUTINE AOTPX(PRESS, TEMP, A_SUM, XI) 
!*******************************************************************************
!COMPUTES THE HELMHOLTZ FREE ENERGY kJ/kg	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR              
      REAL*8 MOL
      CALL MOLX(MOL, XIT, XI) 
      CALL AROTPX(PRESS, TEMP, ALPHA_SUM, XI)         
      A_SUM= ALPHA_SUM*R*TEMP/MOL
      END


!*******************************************************************************
	  SUBROUTINE SOTPX(PRESS, TEMP, ENTRO, XI) 
!*******************************************************************************
!COMPUTES THE ENTROPY	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR              
      REAL*8 MOL, ENTRO
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI) 
      CALL AROTPX(PRESS, TEMP, ALPHA_SUM, XI)        
      CALL D4OTPX(PRESS, TEMP, DALPHA_T, XI)
      ENTRO= (TAU * DALPHA_T - ALPHA_SUM)*R/MOL                     
      END


!*******************************************************************************
	  SUBROUTINE CPOTPX(PRESS, TEMP, CP, XI) 
!*******************************************************************************
!COMPUTES CP	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR              
      REAL*8 MOL, CP
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI) 
      CALL D5OTPX(PRESS, TEMP, DALPHA_TT, XI)       
      CALL D6OTPX(PRESS, TEMP, DALPHA_1, XI)
      CALL D7OTPX(PRESS, TEMP, DALPHA_2, XI)    
      CP= (-TAU*TAU*DALPHA_TT+DALPHA_2*DALPHA_2/DALPHA_1)*R/MOL
      END


!*******************************************************************************
	  SUBROUTINE CVOTPX(PRESS, TEMP, CV, XI) 
!*******************************************************************************
!COMPUTES CV	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR              
      REAL*8 MOL, CV
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI) 
      CALL D5OTPX(PRESS, TEMP, DALPHA_TT, XI)         
      CV= (-TAU*TAU*DALPHA_TT)*R/MOL                   
      END

      
!*******************************************************************************
	  SUBROUTINE WOTPX(PRESS, TEMP, SOUND, XI) 
!*******************************************************************************
!COMPUTES SPEED OF SOUND	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR              
      REAL*8 MOL, SOUND, CP, CV
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI)      
      CALL D6OTPX(PRESS, TEMP, DALPHA_1, XI)
      CALL CVOTPX(PRESS, TEMP, CV, XI)    
      CALL CPOTPX(PRESS, TEMP, CP, XI)                          
      SOUND=   DSQRT( (DALPHA_1*CP*1000D0/CV)*R*TEMP/MOL  )       
      END


!*******************************************************************************
	  SUBROUTINE CAPVOTPX(PRESS, TEMP, ISENCOEFV, XI)
!*******************************************************************************
!COMPUTES ISENTROPIC VOLUME COEFFICIENT	   
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU     
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF, COMPR
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR              
      REAL*8 MOL, ISENCOEFV, CP, CV
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI)      
      CALL D6OTPX(PRESS, TEMP, DALPHA_1, XI)
      CALL CVOTPX(PRESS, TEMP, CV, XI)    
      CALL CPOTPX(PRESS, TEMP, CP, XI)
      CALL ZOTPX(PRESS, TEMP, COMPR, XI)                          
      ISENCOEFV= DALPHA_1 * CP/(CV*COMPR) 
      END


!*******************************************************************************
	  SUBROUTINE CAPTOTPX(PRESS, TEMP, ISENCOEFT, XI) 
!*******************************************************************************
!COMPUTES ISENTROPIC TEMPERATURE COEFFICIENT	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU , INVRHOR      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF, DL
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL              
      REAL*8 MOL,  CP
      REAL*8 ISENCOEFT , DALPHA_1, DALPHA_2    
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI)         
      CALL D6OTPX(PRESS, TEMP, DALPHA_1, XI)
      CALL D7OTPX(PRESS, TEMP, DALPHA_2, XI)
      CALL CPOTPX(PRESS, TEMP, CP, XI)
        ISENCOEFT= 1D0 / (1D0 - (PRESS*1000D0 * DALPHA_2 
     &   *INVRHOR / (MOL*DL*CP* TEMP*DALPHA_1) ) ) 
      END


!*******************************************************************************
	  SUBROUTINE RJTOTPX(PRESS, TEMP, JTCOEF, XI) 
!*******************************************************************************
!COMPUTES JOULE-THOMSON COEFFICIENT	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU , INVRHOR      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF, DL
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL              
      REAL*8 MOL, JTCOEF, CP
      REAL*8 DALPHA_1, DALPHA_2    
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI)      
      CALL D6OTPX(PRESS, TEMP, DALPHA_1, XI)
      CALL D7OTPX(PRESS, TEMP, DALPHA_2, XI)
      CALL CPOTPX(PRESS, TEMP, CP, XI)
      JTCOEF = 1000D0*( (DALPHA_2/ DALPHA_1)-1D0 )*INVRHOR /(MOL*CP*DL)
      END


!*******************************************************************************
	  SUBROUTINE UOTPX(PRESS, TEMP, INTEN, XI) 
!*******************************************************************************
!COMPUTES THE INTERNAL ENERGY	(PHYSICAL DIMENSION)	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR             
      REAL*8 MOL, INTEN
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI)         
      CALL D4OTPX(PRESS, TEMP, DALPHA_T, XI) 
      INTEN = TAU * DALPHA_T*R*TEMP/MOL                         
      END


!*******************************************************************************
	  SUBROUTINE HOTPX(PRESS, TEMP, ENTHA, XI)
!*******************************************************************************
!COMPUTES THE ENTHALPY 	(PHYSICAL DIMENSION)	   
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21),TAU, COMPR      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR              
      REAL*8 MOL, ENTHA
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)       
      CALL MOLX(MOL, XIT, XI)         
      CALL D4OTPX(PRESS, TEMP, DALPHA_T, XI) 
      CALL ZOTPX(PRESS, TEMP, COMPR, XI)       
      ENTHA = (TAU * DALPHA_T+COMPR)*R*TEMP/MOL                         
      END


!*******************************************************************************
	  SUBROUTINE GOTPX(PRESS, TEMP, GIBBS, XI) 
!*******************************************************************************
!COMPUTES THE GIBBS FREE ENERGY  (PHYSICAL DIMENSION)	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)      
	  REAL*8 R, RSTAR, PRESSREF, TEMPREF
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR               
      REAL*8 GIBBS, MOL, COMPR
      CALL MOLX(MOL, XIT, XI) 
      CALL AROTPX(PRESS, TEMP, ALPHA_SUM, XI)   
      CALL ZOTPX(PRESS, TEMP, COMPR, XI)     
      GIBBS = (ALPHA_SUM+COMPR)*R*TEMP/MOL         
      END


!*******************************************************************************    	
	  SUBROUTINE AROTPX(PRESS, TEMP, ALPHA_SUM, XI) 
!*******************************************************************************
!COMPUTES THE REDUCED HELMHOLTZ FREE ENERGY (PHYSICAL DIMENSION)	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k, MINTH(21), MAXTH(21)
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM3/ MINTH, MAXTH
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL            
      REAL*8 ALPHA0_SUM, ALPHAR, ALPHAR0, ALPHA_SUM
      ALPHA_SUM=0.D0
      ALPHA0_SUM=0.D0
      ALPHAR0 = 0.D0
      ALPHAR = 0.D0
      TAU = 0.D0      
        
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
                              
      DO 10 i=1,19
        IF (XI(i).NE.0D0) THEN
      ALPHA0_SUM = ALPHA0_SUM  + XI(i)* (  DLOG(XI(i))+ 
     &   DLOG(  DL / (INVRHOR*RHOC(i)) ) )  + (RSTAR/R)* XI(i) * (
     &   N0_0(i,1) + N0_0(i,2) * (TEMPC(i)/ TEMP) +
     &   N0_0(i,3) * DLOG (TEMPC(i)/ TEMP)   )   
        ENDIF         
 10   CONTINUE
 
       DO 15 i=1,19
         DO 15 k=MINTH(i)+1, MAXTH(i),2       
        IF (XI(i).NE.0D0) THEN
      ALPHA0_SUM = ALPHA0_SUM  + (RSTAR/R)* XI(i) * (
     &   -N0_0(i,k)*DLOG(DABS(DCOSH(THETA_0(i,k)*TEMPC(i)/TEMP))) ) 
        ENDIF  
      
 15   CONTINUE
 
       DO 16 i=1,19
        DO 16 k=MINTH(i), MAXTH(i),2       
        IF (XI(i).NE.0D0) THEN
      ALPHA0_SUM = ALPHA0_SUM  + (RSTAR/R)* XI(i) * (
     &   N0_0(i,k)* 
     &   DLOG(DABS(DSINH(THETA_0(i,k)*TEMPC(i)/TEMP))) ) 
        ENDIF                  
 16   CONTINUE
 
       DO 18 i=20,21
        IF (XI(i).NE.0D0) THEN
      ALPHA0_SUM = ALPHA0_SUM  + XI(i)* (  DLOG(XI(i))+ 
     &   DLOG(  DL / (INVRHOR*RHOC(i)) ) )  + (RSTAR/R)* XI(i) * (
     &   N0_0(i,1) + N0_0(i,2) * TEMPC(i)/ TEMP +
     &   N0_0(i,3) * DLOG (TEMPC(i)/ TEMP) )  
        ENDIF
           
 18   CONTINUE     
 
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        ALPHAR0 = ALPHAR0 + XI(i) * N_0(i,k) * DL ** D_0(i,k) 
     &   * TAU ** t_0(i,k)       
20     CONTINUE      
       
       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        ALPHAR0 = ALPHAR0 + XI(i) * N_0(i,k) * (DL**D_0(i,k))
     &   * (TAU ** t_0(i,k)) * DEXP(-DL**C_0(i,k))        
30     CONTINUE  
      
      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN
         ALPHAR=ALPHAR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * (DL**D(i,j,k)) * (TAU**T(i,j,k))
         ENDIF
40     CONTINUE
   
       DO 50  i=1,20
        DO 50  j=i+1,21
          DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
          IF ( XI(i)*XI(j).NE.0D0 ) THEN
         ALPHAR = ALPHAR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   *(DL**D(i,j,k))* (TAU**T(i,j,k))
     &   * DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0
     &  -BT(i,j,k)*(DL-GM(i,j,k)))
          ENDIF
50     CONTINUE

      ALPHA_SUM= (ALPHAR0 + ALPHAR + ALPHA0_SUM)
      END


!*******************************************************************************
	  SUBROUTINE AR0OTPX(PRESS, TEMP, ALPHA0_SUM, XI) 
!*******************************************************************************
!COMPUTES THE REDUCED HELMHOLTZ FREE ENERGY - IDEAL PART	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,k, MINTH(21), MAXTH(21)
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM3/ MINTH, MAXTH
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL            
      REAL*8 ALPHA0_SUM
      ALPHA0_SUM=0.D0
      TAU = 0.D0      
        
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
                             
      DO 10 i=1,19
        IF (XI(i).NE.0D0) THEN
      ALPHA0_SUM = ALPHA0_SUM  + XI(i)* (  DLOG(XI(i))+ 
     &   DLOG(  DL / (INVRHOR*RHOC(i)) ) )  + (RSTAR/R)* XI(i) * (
     &   N0_0(i,1) + N0_0(i,2) * (TEMPC(i)/ TEMP) +
     &   N0_0(i,3) * DLOG (TEMPC(i)/ TEMP)   )   
        ENDIF         
 10   CONTINUE
 
       DO 15 i=1,19
         DO 15 k=MINTH(i)+1, MAXTH(i),2       
        IF (XI(i).NE.0D0) THEN
      ALPHA0_SUM = ALPHA0_SUM  + (RSTAR/R)* XI(i) * (
     &   -N0_0(i,k)*DLOG(DABS(DCOSH(THETA_0(i,k)*TEMPC(i)/TEMP))) ) 
        ENDIF       
 15   CONTINUE
 
       DO 16 i=1,19
        DO 16 k=MINTH(i), MAXTH(i),2       
        IF (XI(i).NE.0D0) THEN
      ALPHA0_SUM = ALPHA0_SUM  + (RSTAR/R)* XI(i) * (
     &   N0_0(i,k)* 
     &   DLOG(DABS(DSINH(THETA_0(i,k)*TEMPC(i)/TEMP))) ) 
        ENDIF                  
 16   CONTINUE
 
       DO 18 i=20,21
        IF (XI(i).NE.0D0) THEN
      ALPHA0_SUM = ALPHA0_SUM  + XI(i)* (  DLOG(XI(i))+ 
     &   DLOG(  DL / (INVRHOR*RHOC(i)) ) )  + (RSTAR/R)* XI(i) * (
     &   N0_0(i,1) + N0_0(i,2) * TEMPC(i)/ TEMP +
     &   N0_0(i,3) * DLOG (TEMPC(i)/ TEMP) )  
        ENDIF           
 18   CONTINUE     
 
      END


!*******************************************************************************
	  SUBROUTINE ARROTPX(PRESS, TEMP, ALPHAR_SUM, XI) 
!*******************************************************************************
!COMPUTES THE REDUCED HELMHOLTZ FREE ENERGY - RESIDUAL PART	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k, MINTH(21), MAXTH(21)
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM3/ MINTH, MAXTH
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL            
      REAL*8 ALPHAR, ALPHAR0, ALPHAR_SUM
      ALPHAR_SUM=0.D0
      ALPHAR0 = 0.D0
      ALPHAR = 0.D0
      TAU = 0.D0      
        
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
                              
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        ALPHAR0 = ALPHAR0 + XI(i) * N_0(i,k) * DL ** D_0(i,k) 
     &   * TAU ** t_0(i,k)       
20     CONTINUE      
       
       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        ALPHAR0 = ALPHAR0 + XI(i) * N_0(i,k) * (DL**D_0(i,k))
     &   * (TAU ** t_0(i,k)) * DEXP(-DL**C_0(i,k))        
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN
         ALPHAR=ALPHAR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * (DL**D(i,j,k)) * (TAU**T(i,j,k))
         ENDIF
40     CONTINUE
    
       DO 50  i=1,20
        DO 50  j=i+1,21
          DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
          IF ( XI(i)*XI(j).NE.0D0 ) THEN
         ALPHAR = ALPHAR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   *(DL**D(i,j,k))* (TAU**T(i,j,k))
     &   * DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0
     &  -BT(i,j,k)*(DL-GM(i,j,k)))
          ENDIF
50     CONTINUE

      ALPHAR_SUM= ALPHAR0 + ALPHAR
      END


!*******************************************************************************
	  SUBROUTINE D3OTPX(PRESS, TEMP, DALPHA_DT, XI) 
!*******************************************************************************
!COMPUTES DALPHA_DT	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP    
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL              
      REAL*8 DALPHA_DT, DALPHA_DT0, DALPHA_DTR

      DALPHA_DT=0.D0
      DALPHA_DT0=0.D0
      DALPHA_DTR = 0.D0
      TAU = 0.D0      
           
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
     
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_DTR = DALPHA_DTR + XI(i) * N_0(i,k) * DL** (D_0(i,k)-1D0)
     &   * TAU ** (t_0(i,k)-1D0) * D_0(i,k) * t_0(i,k)
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DTR = DALPHA_DTR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-1D0)
     &   * TAU**(t_0(i,k)-1D0) * (D_0(i,k)- C_0(i,k) *DL**C_0(i,k))
     &   * DEXP(-DL**C_0(i,k))  * t_0(i,k)      
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DTR=DALPHA_DTR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-1D0) *D(i,j,k)* TAU**(T(i,j,k)-1D0) *T(i,j,k)
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DTR=DALPHA_DTR+XI(i)*XI(j)*F(i,j)*N(i,j,k)
     &  *DL**(D(i,j,k)-1D0) * TAU ** (T(i,j,k)-1D0) *T(i,j,k)
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ( D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k) )
50     CONTINUE

      DALPHA_DT= DALPHA_DTR+DALPHA_DT0
      END

      
!*******************************************************************************
	  SUBROUTINE D7OTPX(PRESS, TEMP, DALPHA_2, XI)
!*******************************************************************************
!COMPUTES DALPHA_2	   
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, COMPR, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP    
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL              
      REAL*8 DALPHA_2, DALPHA_DTR

      DALPHA_2 = 0.D0
      DALPHA_DTR = 0.D0
      TAU = 0.D0      
           
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
      CALL ZOTPX(PRESS, TEMP, COMPR, XI)
                       
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_DTR = DALPHA_DTR + XI(i) * N_0(i,k) * DL** (D_0(i,k)-1D0)
     &   * TAU ** (t_0(i,k)-1D0) * D_0(i,k) * t_0(i,k)
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DTR = DALPHA_DTR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-1D0)
     &   * TAU**(t_0(i,k)-1D0) * (D_0(i,k)- C_0(i,k) *DL**C_0(i,k))
     &   * DEXP(-DL**C_0(i,k))  * t_0(i,k)      
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DTR=DALPHA_DTR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-1D0) *D(i,j,k)* TAU**(T(i,j,k)-1D0) *T(i,j,k)
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DTR=DALPHA_DTR+XI(i)*XI(j)*F(i,j)*N(i,j,k)
     &  *DL**(D(i,j,k)-1D0) * TAU ** (T(i,j,k)-1D0) *T(i,j,k)
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ( D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k) )
50     CONTINUE

      DALPHA_2= COMPR-TAU*DL*DALPHA_DTR
      END


!*******************************************************************************
	  SUBROUTINE D6OTPX(PRESS, TEMP, DALPHA_1, XI)
!*******************************************************************************
!COMPUTES DALPHA_1	   
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL         
      REAL*8 DALPHA_1, DALPHA_DDR, DALPHA_DR
      
      DALPHA_1=0.D0
      DALPHA_DDR = 0.D0
      DALPHA_DR = 0.D0
      TAU = 0.D0
      TEMPR=0D0      
           
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
         
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_DDR = DALPHA_DDR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-2D0) 
     &   * TAU ** t_0(i,k) * D_0(i,k) * (D_0(i,k) -1D0)    
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DDR = DALPHA_DDR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-2D0)
     &   * TAU**t_0(i,k) * (  (D_0(i,k)- C_0(i,k) *DL**C_0(i,k)) * 
     &   (D_0(i,k)-1D0- C_0(i,k) *DL**C_0(i,k)) 
     &   -C_0(i,k)**2D0*DL**C_0(i,k)  )
     &   * DEXP(-DL**C_0(i,k))        
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DDR=DALPHA_DDR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-2D0) *D(i,j,k)* 
     &   (D(i,j,k)-1D0)* (TAU**T(i,j,k))
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DDR=DALPHA_DDR+XI(i)*XI(j)*F(i,j)*
     &   N(i,j,k)*DL**(D(i,j,k)-2D0) * TAU ** T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ( (  D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k))**2D0-D(i,j,k)
     &   - 2D0*ETA(i,j,k)*DL**2D0   )
50     CONTINUE

      DO 60 i=1,21
        DO 60 k = 1, K1POL(i)             
        DALPHA_DR = DALPHA_DR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-1D0) 
     &   * TAU ** t_0(i,k) * D_0(i,k)     
60     CONTINUE      

       DO 70 i=1,21
        DO 70 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DR = DALPHA_DR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-1D0)
     &   * TAU**t_0(i,k) * (D_0(i,k)- C_0(i,k) *DL**C_0(i,k))
     &   * DEXP(-DL**C_0(i,k))        
70     CONTINUE  

       DO 80  i=1,20
        DO 80  j=i+1,21   
         DO 80 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DR=DALPHA_DR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-1D0) *D(i,j,k)* (TAU**T(i,j,k))
         ENDIF
80     CONTINUE
    
       DO 90  i=1,21-1
        DO 90  j=i+1,21
         DO 90 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DR=DALPHA_DR+XI(i)*XI(j)*F(i,j)*N(i,j,k)*
     &   DL**(D(i,j,k)-1D0)* TAU ** T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ( D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k) )
90     CONTINUE
      
      DALPHA_1= 1D0+2D0*DL*DALPHA_DR+DL*DL* DALPHA_DDR     
      END


!*******************************************************************************
	  SUBROUTINE D1OTPX(PRESS, TEMP, DALPHA_D, XI) 
!*******************************************************************************
!COMPUTES DALPHA_D	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP    
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL              
      REAL*8 DALPHA_D, DALPHA_D0, DALPHA_DR

      DALPHA_D=0.D0
      DALPHA_D0=0.D0
      DALPHA_DR = 0.D0
      TAU = 0.D0      
          
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 

      DALPHA_D0 = 1D0/DL
      
      DO 10 i=1,21
       DALPHA_DD0= DALPHA_DD0 + 1D0/DL * (RSTAR/R)* XI(i)
10    CONTINUE      
      
      
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_DR = DALPHA_DR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-1D0) 
     &   * TAU ** t_0(i,k) * D_0(i,k)     
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DR = DALPHA_DR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-1D0)
     &   * TAU**t_0(i,k) * (D_0(i,k)- C_0(i,k) *DL**C_0(i,k))
     &   * DEXP(-DL**C_0(i,k))        
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DR=DALPHA_DR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-1D0) *D(i,j,k)* (TAU**T(i,j,k))
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DR=DALPHA_DR+XI(i)*XI(j)*F(i,j)*N(i,j,k)*
     &   DL**(D(i,j,k)-1D0)* TAU ** T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ( D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k) )
50     CONTINUE

      DALPHA_D= DALPHA_DR+DALPHA_D0
      END
      
      
!*******************************************************************************
	  SUBROUTINE D5OTPX(PRESS, TEMP, DALPHA_TT, XI) 
!*******************************************************************************
!COMPUTES DALPHA_TT	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k, MINTH(21), MAXTH(21)
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM3/ MINTH, MAXTH	  
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL      
      REAL*8 DALPHA_TT, DALPHA_TT0, DALPHA_TTR
      
      DALPHA_TT=0.D0
      DALPHA_TT0=0.D0
      DALPHA_TTR = 0.D0
      TAU = 0.D0
      TEMPR=0D0      
            
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
     
      DO 10 i=1,19
        IF (XI(i).NE.0D0) THEN
      DALPHA_TT0 = DALPHA_TT0+(RSTAR/R)* XI(i)*
     &   (TEMPC(i)/TEMPR)**2D0 *(-
     &   N0_0(i,3) * (TEMP/TEMPC(i))**2D0)
       ENDIF    
 10   CONTINUE

      DO 14 i=1,19
       DO 14 k=MINTH(i),MAXTH(i),2
        IF (XI(i).NE.0D0) THEN
      DALPHA_TT0 = DALPHA_TT0+(RSTAR/R)* XI(i)*
     &   (TEMPC(i)/TEMPR)**2D0 *(-N0_0(i,k)* 
     &   THETA_0(i,k)**2D0/(DSINH(THETA_0(i,k)*TEMPC(i)/TEMP))**2D0)
       ENDIF    
 14   CONTINUE 
 
 
      DO 16 i=1,19
       DO 16 k=MINTH(i)+1,MAXTH(i),2      
        IF (XI(i).NE.0D0) THEN
      DALPHA_TT0 = DALPHA_TT0+(RSTAR/R)* XI(i)*
     &   (TEMPC(i)/TEMPR)**2D0 *(-N0_0(i,k)* 
     &   THETA_0(i,k)**2D0/(DCOSH(THETA_0(i,k)*TEMPC(i)/TEMP))**2D0)
       ENDIF    
 16   CONTINUE
    
      DO 18 i=20,21
        IF (XI(i).NE.0D0) THEN
      DALPHA_TT0 = DALPHA_TT0+(RSTAR/R)* XI(i)*
     &   (TEMPC(i)/TEMPR)**2D0 *(-
     &   N0_0(i,3) * (TEMP/TEMPC(i))**2D0)
       ENDIF    
 18   CONTINUE 
       
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_TTR = DALPHA_TTR + XI(i) * N_0(i,k)*DL**D_0(i,k)
     &   * TAU ** (t_0(i,k)-2D0) * t_0(i,k)* (t_0(i,k)-1D0)     
20     CONTINUE      
     
       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_TTR = DALPHA_TTR + XI(i) * N_0(i,k)* 
     &   t_0(i,k) *  (t_0(i,k)-1D0) 
     &   * DL**D_0(i,k) * TAU**(t_0(i,k)-2D0) * DEXP(-DL**C_0(i,k))
30     CONTINUE
            
      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_TTR=DALPHA_TTR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**D(i,j,k) *T(i,j,k)*(T(i,j,k)-1D0)* TAU**(T(i,j,k)-2D0)
         ENDIF
40     CONTINUE

       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_TTR=DALPHA_TTR+XI(i)*XI(j)*F(i,j)*N(i,j,k)*DL**D(i,j,k)
     &   * TAU ** (T(i,j,k)-2D0) * T(i,j,k)* (T(i,j,k)-1D0) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
50     CONTINUE

      DALPHA_TT=DALPHA_TT0+DALPHA_TTR
      END


!*******************************************************************************
	  SUBROUTINE D11OTPX(PRESS, TEMP, DALPHA_DR, XI) 
!*******************************************************************************
!COMPUTES DALPHA_D residual part	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP    
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL              
      REAL*8 DALPHA_DR

      DALPHA_DR = 0.D0
      TAU = 0.D0      
          
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
    
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_DR = DALPHA_DR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-1D0) 
     &   * TAU ** t_0(i,k) * D_0(i,k)     
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DR = DALPHA_DR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-1D0)
     &   * TAU**t_0(i,k) * (D_0(i,k)- C_0(i,k) *DL**C_0(i,k))
     &   * DEXP(-DL**C_0(i,k))        
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DR=DALPHA_DR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-1D0) *D(i,j,k)* (TAU**T(i,j,k))
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DR=DALPHA_DR+XI(i)*XI(j)*F(i,j)*N(i,j,k)*
     &   DL**(D(i,j,k)-1D0)* TAU ** T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ( D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k) )
50     CONTINUE

      END
      
      
!*******************************************************************************
	  SUBROUTINE D9OTPX(PRESS, TEMP, DALPHA_DDDR, XI) 
!*******************************************************************************
!COMPUTES DALPHA_DDD residual part	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL         
      REAL*8 DALPHA_DDDR
      
      DALPHA_DDDR = 0.D0
      TAU = 0.D0
      TEMPR=0D0      
           
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 

       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_DDDR = DALPHA_DDDR + XI(i)*N_0(i,k)* DL**(D_0(i,k)-3D0) 
     &   * TAU ** t_0(i,k) * D_0(i,k) * (D_0(i,k) -1D0)*(D_0(i,k) -2D0)
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DDDR = DALPHA_DDDR + XI(i)* N_0(i,k)* DL**(D_0(i,k)-3D0)
     &   * TAU**t_0(i,k) * (   ( (D_0(i,k)- C_0(i,k) *DL**C_0(i,k)) * 
     &   (D_0(i,k)-2D0- C_0(i,k) *DL**C_0(i,k)) 
     &   -3D0*C_0(i,k)**2D0*DL**C_0(i,k) ) * (D_0(i,k)- 1D0-C_0(i,k) 
     &   *DL**C_0(i,k)) - C_0(i,k)**3D0*DL**C_0(i,k) )
     &   * DEXP(-DL**C_0(i,k))        
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DDDR=DALPHA_DDDR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-3D0) *D(i,j,k)* 
     &   (D(i,j,k)-1D0)*(D(i,j,k)-2D0)* (TAU**T(i,j,k))
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DDDR=DALPHA_DDDR+XI(i)*XI(j)*F(i,j)*
     &   N(i,j,k)*DL**(D(i,j,k)) * TAU ** T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*
     &  (DL-GM(i,j,k)) ) * ( ( (D(i,j,k)/DL)-2D0*ETA(i,j,k)*
     &  (DL-EPS(i,j,k))-BT(i,j,k) )**3D0-3D0*
     &  ( (D(i,j,k)/DL)-2D0*ETA(i,j,k)*(DL-EPS(i,j,k))-BT(i,j,k) )* 
     &  (  D(i,j,k)/(DL*DL)+ 2D0*ETA(i,j,k) ) + 
     &  2D0* ( D(i,j,k)/(DL*DL*DL) )  )

50     CONTINUE

      END
      
      
!*******************************************************************************
	  SUBROUTINE D10OTPX(PRESS, TEMP, DALPHA_DDR, XI) 
!*******************************************************************************
!COMPUTES DALPHA_DD residual part	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL         
      REAL*8 DALPHA_DDR
      
      DALPHA_DDR = 0.D0
      TAU = 0.D0
      TEMPR=0D0      
           
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
     
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_DDR = DALPHA_DDR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-2D0) 
     &   * TAU ** t_0(i,k) * D_0(i,k) * (D_0(i,k) -1D0)    
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DDR = DALPHA_DDR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-2D0)
     &   * TAU**t_0(i,k) * (  (D_0(i,k)- C_0(i,k) *DL**C_0(i,k)) * 
     &   (D_0(i,k)-1D0- C_0(i,k) *DL**C_0(i,k)) 
     &   -C_0(i,k)**2D0*DL**C_0(i,k)  )
     &   * DEXP(-DL**C_0(i,k))        
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DDR=DALPHA_DDR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-2D0) *D(i,j,k)* 
     &   (D(i,j,k)-1D0)* (TAU**T(i,j,k))
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DDR=DALPHA_DDR+XI(i)*XI(j)*F(i,j)*
     &   N(i,j,k)*DL**(D(i,j,k)-2D0) * TAU ** T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ((   D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k)   )**2D0-D(i,j,k)
     &   - 2D0*ETA(i,j,k)*DL**2D0   )
50     CONTINUE
      END
      
      
!*******************************************************************************
	  SUBROUTINE D2OTPX(PRESS, TEMP, DALPHA_DD, XI) 
!*******************************************************************************
!COMPUTES DALPHA_DD	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL         
      REAL*8 DALPHA_DD, DALPHA_DD0, DALPHA_DDR
      
      DALPHA_DD=0.D0
      DALPHA_DD0=0.D0
      DALPHA_DDR = 0.D0
      TAU = 0.D0
      TEMPR=0D0      
           
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
      DO 10 i=1,21
       DALPHA_DD0= DALPHA_DD0 -(1D0/DL)**2D0 * (RSTAR/R)* XI(i)
10    CONTINUE
      
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_DDR = DALPHA_DDR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-2D0) 
     &   * TAU ** t_0(i,k) * D_0(i,k) * (D_0(i,k) -1D0)    
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_DDR = DALPHA_DDR + XI(i) * N_0(i,k) * DL**(D_0(i,k)-2D0)
     &   * TAU**t_0(i,k) * (  (D_0(i,k)- C_0(i,k) *DL**C_0(i,k)) * 
     &   (D_0(i,k)-1D0- C_0(i,k) *DL**C_0(i,k)) 
     &   -C_0(i,k)**2D0*DL**C_0(i,k)  )
     &   * DEXP(-DL**C_0(i,k))        
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_DDR=DALPHA_DDR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**(D(i,j,k)-2D0) *D(i,j,k)* 
     &   (D(i,j,k)-1D0)* (TAU**T(i,j,k))
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_DDR=DALPHA_DDR+XI(i)*XI(j)*F(i,j)*
     &   N(i,j,k)*DL**(D(i,j,k)-2D0) * TAU ** T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ((   D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k)   )**2D0-D(i,j,k)
     &   - 2D0*ETA(i,j,k)*DL**2D0   )
50     CONTINUE

      DALPHA_DD= DALPHA_DDR+DALPHA_DD0   
      END
      
      
!*******************************************************************************
	  SUBROUTINE D8OTPX(PRESS, TEMP, DALPHA_TR, XI) 
!*******************************************************************************
!COMPUTES DALPHA_T Residual Part	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k, MINTH(21), MAXTH(21)
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM3/ MINTH, MAXTH	  
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL        
      REAL*8 DALPHA_TR
      
      DALPHA_TR = 0.D0
      TAU = 0.D0
      TEMPR=0D0      
            
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
    
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_TR = DALPHA_TR + XI(i) * N_0(i,k)*DL**D_0(i,k)
     &   * TAU ** (t_0(i,k)-1D0) * t_0(i,k)     
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_TR = DALPHA_TR + XI(i) * N_0(i,k)* t_0(i,k)
     &   * DL**D_0(i,k) * TAU**(t_0(i,k)-1D0) * DEXP(-DL**C_0(i,k))
30     CONTINUE  

      DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_TR=DALPHA_TR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**D(i,j,k) *T(i,j,k)* TAU**(T(i,j,k)-1D0)
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_TR=DALPHA_TR+XI(i)*XI(j)*F(i,j)*N(i,j,k)*DL**D(i,j,k)
     &   * TAU ** (T(i,j,k)-1D0) * T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
50     CONTINUE
      
      END
      
      
!*******************************************************************************
	  SUBROUTINE D4OTPX(PRESS, TEMP, DALPHA_T, XI)
!*******************************************************************************
!COMPUTES DALPHA_T	   
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k, MINTH(21), MAXTH(21)
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM3/ MINTH, MAXTH	  
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL        
      REAL*8 DALPHA_T, DALPHA_T0, DALPHA_TR
      
      DALPHA_T=0.D0
      DALPHA_T0=0.D0
      DALPHA_TR = 0.D0
      TAU = 0.D0
      TEMPR=0D0      
            
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 

      DO 14 i=1,19
        IF (XI(i).NE.0D0) THEN
      DALPHA_T0 = DALPHA_T0+(RSTAR/R)* XI(i)*
     &   (TEMPC(i)/TEMPR)*(N0_0(i,2) + N0_0(i,3) * (TEMP/TEMPC(i)))
       ENDIF 
 14   CONTINUE

      DO 16 i=1,19
       DO 16 k=MINTH(i), MAXTH(i),2
        IF (XI(i).NE.0D0) THEN
      DALPHA_T0 = DALPHA_T0+(RSTAR/R)* XI(i)*(TEMPC(i)/TEMPR)*
     &   (+ N0_0(i,k)* 
     &   THETA_0(i,k)/(DTANH (THETA_0(i,k)*TEMPC(i)/TEMP)) )
       ENDIF      
 16   CONTINUE

       DO 17 i=1,19
         DO 17 k=MINTH(i)+1, MAXTH(i),2
        IF (XI(i).NE.0D0) THEN
      DALPHA_T0 = DALPHA_T0+(RSTAR/R)* XI(i)*(TEMPC(i)/TEMPR)*
     &   (- N0_0(i,k)* 
     &   THETA_0(i,k)*(DTANH (THETA_0(i,k)*TEMPC(i)/TEMP) ))
       ENDIF     
 17    CONTINUE
  
       DO 18 i=20,21
        IF (XI(i).NE.0D0) THEN
      DALPHA_T0 = DALPHA_T0+(RSTAR/R)* XI(i)*
     &   (TEMPC(i)/TEMPR)*(N0_0(i,2) +N0_0(i,3) * (TEMP/TEMPC(i)) )
       ENDIF 
 18    CONTINUE 
      
       DO 20 i=1,21
        DO 20 k = 1, K1POL(i)             
        DALPHA_TR = DALPHA_TR + XI(i) * N_0(i,k)*DL**D_0(i,k)
     &   * TAU ** (t_0(i,k)-1D0) * t_0(i,k)     
20     CONTINUE      

       DO 30 i=1,21
        DO 30 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        DALPHA_TR = DALPHA_TR + XI(i) * N_0(i,k)* t_0(i,k)
     &   * DL**D_0(i,k) * TAU**(t_0(i,k)-1D0) * DEXP(-DL**C_0(i,k))
30     CONTINUE  

       DO 40  i=1,20
        DO 40  j=i+1,21   
         DO 40 k=1,K2POL(i,j)
         IF ( XI(i)*XI(j).NE.0D0 ) THEN          
         DALPHA_TR=DALPHA_TR + XI(i)*XI(j)*F(i,j)* N(i,j,k) 
     &   * DL**D(i,j,k) *T(i,j,k)* TAU**(T(i,j,k)-1D0)
         ENDIF
40     CONTINUE
    
       DO 50  i=1,21-1
        DO 50  j=i+1,21
         DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        DALPHA_TR=DALPHA_TR+XI(i)*XI(j)*F(i,j)*N(i,j,k)*DL**D(i,j,k)
     &   * TAU ** (T(i,j,k)-1D0) * T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
50     CONTINUE
      
       DALPHA_T=DALPHA_T0+DALPHA_TR
      END
      
      
!*******************************************************************************
	  SUBROUTINE ZOTPX(PRESS, TEMP, COMPR, XI) 
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: ZOTPX
!*******************************************************************************
!COMPUTES THE COMPRESSIBILITY
!USES DL AS COMMON TO CALCULATE COMPRESSIBILITY	  
!INPUT: P,T, XI
!OUTPUT: COMPRESSIBILITY	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, DL, COMPR, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP
      COMMON /CONST2/ DL
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
      CALL ZOTDLX(COMPR, TEMP, DL, XI)
      END
      
      
!******************************************************************************* 	  	    
	  SUBROUTINE DOTPX(PRESS, TEMP, DENS, XI) 
!*******************************************************************************
!COMPUTES THE DENSITY
!INPUT: P, T, XI
!OUTPUT: DENSITY (DENS) (IN PHYSICAL DIMENSION)	  
      IMPLICIT DOUBLE PRECISION (L-Z)	  
      REAL*8 PRESS, DL, TAU, XI(21), MOL, XIT
      REAL*8 TEMPR, INVRHOR, TEMP, DENS, MCOMP(21)
	  COMMON /PARAM2/ MCOMP
	  COMMON /CONST2/ DL
	  MOL=0D0
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
 	  CALL MOLX(MOL, XIT, XI)           
      DENS=  DL*MOL/INVRHOR 
      END
      
      
!*******************************************************************************
	  SUBROUTINE REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
!*******************************************************************************
!COMPUTES REDUCED MIXTURE FUNCTIONS FOR DENSITY AND TEMPERATURE	  
!INPUT : TEMP, XI
!OUTPUT: DENSITY REDUCED FUNCTION (INVRHOR=1/RHO_r)
!OUTPUT: TEMPERATURE REDUCED FUNCTION (TEMPR)
!OUTPUT: REDUCED TEMPERATURE TAU = TEMP_r/TEMP
       IMPLICIT DOUBLE PRECISION (A-Z)
       INTEGER i,j
       REAL*8 RHOC(21),TEMPC(21)
	   REAL*8 BTV(21,21), GMV(21,21), XI(21), INVRHOR
	   REAL*8 BTT(21,21), GMT(21,21), TAU, TEMP, TEMPR	
	   COMMON /PARAM1/ BTV, GMV, BTT, GMT
	   COMMON /PARAM5/RHOC, TEMPC
      
         INVRHOR  = 0D0
         TEMPR = 0D0 
       
! OVERWRITE ALL THE BTV, GMV, BTT, GMT ZERO VALUES WITH UNITY             
         DO 10 I=1,21
         DO 10 J=1,21
        IF  (BTV(I,J).EQ.0D0) THEN
        BTV(I,J) =1D0
        ENDIF        
        IF  (GMV(I,J).EQ.0D0) THEN
        GMV(I,J) =1D0
        ENDIF        
        IF  (BTT(I,J).EQ.0D0) THEN
        BTT(I,J) =1D0
        ENDIF        
        IF  (GMT(I,J).EQ.0D0) THEN 
        GMT(I,J) =1D0
        ENDIF                            
 10      CONTINUE  

        DO 20 I=1,21
         DO 20 J=I+1,21
        BTV(J,I) = 1D0/BTV(I,J)
        BTT(J,I) = 1D0/BTT(I,J)       
 20      CONTINUE 
  
        DO 30 I=1,21
         DO 30 J=i+1,21
        GMV(J,I) = GMV(I,J)
        GMT(J,I) = GMT(I,J)     
 30      CONTINUE  
           
       DO 40 i=1,21     
        DO 40 j=1,21
        IF (XI(i)*XI(j).NE.0) Then        
        INVRHOR = INVRHOR + XI(i) * XI(j) * BTV(i,j) * GMV(i,j)
     &   * ((XI(i) + XI(j)) / (BTV(i,j)**2D0*XI(i) + XI(j))) * 0.125D0
     &   * (1D0/(RHOC(i)**(1D0/3D0)) + 1D0/(RHOC(j)**(1D0/3D0)))**3D0
        TEMPR = TEMPR + XI(i) * XI(j) * BTT(i,j) * GMT(i,j)
     &   * ((XI(i) + XI(j)) / (BTT(i,j)**2D0*XI(i) + XI(j)))
     &   *  (TEMPC(i) * TEMPC(j))**0.5D0           
        ENDIF      	
40      CONTINUE   
       
       TAU=TEMPR/TEMP              
      END
      
      
!*******************************************************************************	  
	  SUBROUTINE POTDLX(DL, PRESS, TEMP, XI) 
!*******************************************************************************
!COMPUTES PRESSURE BASED ON DL, TEMPERATURE AND COMPOSITION
!INPUT : REDUCED DENSITY DL, T, XI
!OUTPUT: PRESSURE	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, DL, COMPR, TAU,XI(21),TEMPR, INVRHOR, TEMP,R
      COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
      CALL ZOTDLX(COMPR, TEMP, DL, XI)
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)   
      PRESS=  COMPR * R*TEMPR/(TAU*INVRHOR) * DL 
      PRESS = PRESS / 1000D0   !CONVERT kPa to MPa 
      END
      
      
!*******************************************************************************
      SUBROUTINE ZOTDLX(COMPR, TEMP, DL, XI)
!*******************************************************************************
!COMPUTES COMPRESSIBILITY BASED ON DL
!OUTPUT: REDUCED DENSITY (DL) 
!REFER TO EQ. No. 34 - ISO20765-Part II 
!INPUT : REDUCED DENSITY DL, T, XI
!OUTPUT: COMPRESSIBILITY

      IMPLICIT DOUBLE PRECISION (A-Z)
!BLOCK DATA PARAM       
       REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	   REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	   REAL*8 BTV(21,21), GMV(21,21),PRESSREF, TEMPREF, TEMPR 
	   REAL*8 BTT(21,21), GMT(21,21),	R, RSTAR, INVRHOR
	   REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	   REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)

	   INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	   COMMON /PARAM/ N0_0, THETA_0, F
	   COMMON /PARAM2/ MCOMP	   
	   COMMON /PARAM/ N_0, C_0, D_0, T_0
	   COMMON /PARAM1/ BTV, GMV, BTT, GMT
	   COMMON /PARAM5/ RHOC, TEMPC
	   COMMON /PARAM/ D, T, N, ETA 
	   COMMON /PARAM/ EPS, BT, GM
       COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	   COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
!----       
       INTEGER i,j,k
       REAL*8 XI(21), TAU, DL
	   REAL*8 COMPR
	   
       COMPR= 0D0
	   CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI)  	   
       
       DO 10  i=1,21
        DO 10 k=1,K1POL(i)
        COMPR=COMPR + XI(i) * N_0(i,k) * D_0(i,k) * DL ** D_0(i,k) 
     &   * TAU ** t_0(i,k)
10     CONTINUE        

       DO 20  i=1,21
        DO 20 k=K1POL(i)+1,K1POL(i)+K1EXP(i)
        COMPR=COMPR + XI(i) * N_0(i,k) * (D_0(i,k) - C_0(i,k) 
     &   * DL ** C_0(i,k)) * DL ** D_0(i,k) 
     &   * TAU ** t_0(i,k) * DEXP(-DL**C_0(i,k))
20     CONTINUE   
             
       DO 30  i=1,21-1
        DO 30  j=i+1,21
         DO 30 k=1,K2POL(i,j)
         COMPR=COMPR + XI(i)*XI(j)*F(i,j)* N(i,j,k) * D(i,j,k) 
     &   * DL ** D(i,j,k) * TAU ** T(i,j,k)
30     CONTINUE
    
       DO 40  i=1,21-1
        DO 40  j=i+1,21
         DO 40 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
         COMPR=COMPR + XI(i)*XI(j)*F(i,j)* N(i,j,k) * DL ** D(i,j,k)
     &   * TAU ** T(i,j,k) 
     &  *DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0-BT(i,j,k)*(DL-GM(i,j,k)))
     &   * ( D(i,j,k)-2D0*ETA(i,j,k)*DL
     &   * (DL-EPS(i,j,k))-DL*BT(i,j,k) )
40     CONTINUE
       COMPR = COMPR+1D0                     
	  END


!      SUBROUTINE SWITCH
!      IMPLICIT DOUBLE PRECISION (A-Z)	
!	  INTEGER SW, M
!	  COMMON /CONST3/ SW
	  
!	  CALL D2PDT(PRESS, TEMP, M, XI)  
	  
!	  SW=2
!      END

!*******************************************************************************
	  SUBROUTINE DZOFPT(PRESS, TEMP,  FDL, XI)
!*******************************************************************************
!MULTIPLE ROOT DENSITY SOLVER - SOLVE FOR REDUCED DENSITY
!INTERVALS ARE FIRST ESTIMATED USING SRK EOS 
!REGULA FALSI METHOD IS APPLIED WITHIN THE DEFINED INTERVALS
!CALCULATES VAPOR ROOT AND RETURNS AS COMMON DATA (DEFAULT)
!ATTEMPT TO CALCULATE LIQUID ROOT (FOR PROVISIONAL USE)
!OUTPUT: FDL ARRAY
!FDL(1) = VAPOR ROOT
!FDL(2) = LIQUID ROOT (WHEN AVAILABLE)
!DL REDUCED DENSITY AS COMMON	  
      IMPLICIT DOUBLE PRECISION (A-Z)	
	  INTEGER MTYPE, m, SW, AUTO
      REAL*8 PRESS, XI(21), CFI1(21), CFI2(21)
	  REAL*8 FDL(2)
	  REAL*8 SH, UPPER, LOWER, LB
	  REAL*8 DLV, DLL, DL99, DL
	  COMMON /PARAM9/ DLV, DLL, DL99, MTYPE
	  COMMON /CONST2/ DL
	  COMMON /CONST3/ SW
	  COMMON /CONST4/ AUTO
!	  auto=1
!	  RT=1 ! VAPOR ROOT (DEFAULT)	  
      !INITIAL ROOT ESTIMATION - USING SRK EOS
	  CALL SRK(PRESS, TEMP, CFI1, CFI2, XI) 

!VAPOR ROOT
      m=0
	  LB=1D0  
4	  SH=0D0
5	  SH=SH+0.02D0
	  UPPER=1D0+SH
	  LOWER=1D0-SH
	  m=m+1
     
	  CALL RGFALSI(PRESS, TEMP, FDL(1), XI, LB*DLV*LOWER, DLV*UPPER)

	  IF (m.GE.20) THEN
	  FDL(1)=0D0
	  GOTO 7
	  ENDIF	 
	   
	  IF (FDL(1).EQ.0D0) THEN

	   GOTO 5
	  ENDIF
7     CONTINUE	  
	       
	  IF (FDL(1).LT.0D0) THEN
	  LB=0D0
	  GOTO 4
	  ENDIF	  

 28	  IF (AUTO.EQ.0) THEN
	  DL= FDL(1)
	  ENDIF
	  
!	  PRINT *, "DENSITY SOLVER - ", "ROOT TYPE:", SW, "DL=", DL
       
	  RETURN 
	  END

	  	  
!*******************************************************************************
	  SUBROUTINE RGFALSI(PRESS, TEMP,  DL, XI, X1, X2)
!*******************************************************************************
!"REGULA FALSI" METHOD SEE ISO20765-PART I OR AGA8DC-92 FORMULATION
!*******************************************************************************	  
      IMPLICIT DOUBLE PRECISION (A-Z)	
	  INTEGER i
      REAL*8 PRESS, DL, XI(21)
	  REAL*8 X1, X2, X3, FF, F1, F2, F3, TOL, TEMP

	  TOL = 0.5D-09
	  DL = 0D0
	  CALL POTDLX(X1, F1, TEMP, XI)
	  CALL POTDLX(X2, F2, TEMP, XI) 
	  F1 = F1 - PRESS
	  F2 = F2 - PRESS

	  IF (F1*F2.GE.0) RETURN
!----------------------------------------------------------------------
!BEGIN ITERATING
!----------------------------------------------------------------------
	  DO 60 I = 1, 50
!...Use False Position to get point 3.
	  X3 = X1 - F1*(X2 - X1)/(F2 - F1)
	  CALL POTDLX(X3, F3, TEMP, XI) 
	  F3 = F3 - PRESS
!...Use points 1, 2, and 3 to estimate the root using Chamber's
!...method (quadratic solution).
	  DL = X1*F2*F3/((F1 - F2)*(F1 - F3)) 
     & + X2*F1*F3/((F2 - F1)*(F2 - F3))
     & + X3*F1*F2/((F3 - F1)*(F3 - F2))
	  IF ((DL - X1)*(DL - X2).GE.0) DL = (X1 + X2)/2.D0
	  CALL POTDLX(DL, FF, TEMP, XI) 
	  FF = FF - PRESS
	  IF (DABS(FF).LE.TOL) RETURN
!...Discard quadratic solution if false position root is closer.
	  IF (DABS(F3).LT.DABS(FF) .AND. FF*F3.GT.0) THEN
	  IF (F3*F1.GT.0) THEN
	  X1 = X3
	  F1 = F3
	  ELSE
	  X2 = X3
	  F2 = F3
	  ENDIF
	  ELSE
!...Swap in new value from quadratic solution
	  IF (FF*F3.LT.0) THEN
	  X1 = DL
	  F1 = FF
	  X2 = X3
	  F2 = F3
	  ELSEIF (F3*F1.GT.0) THEN
	  X1 = DL
	  F1 = FF
	  ELSE
	  X2 = DL
	  F2 = FF
	  ENDIF
	  ENDIF
60	  CONTINUE
      DL=0D0
      RETURN
	  END
	  
	  
!*******************************************************************************
	  SUBROUTINE MUOTPX(PRESS, TEMP, CHEM, XI) 
!*******************************************************************************
!COMPUTES CHEMICAL POTENTIAL OF COMPONENT i	  
!USES GERG 2008 MODEL	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k, MINTH(21), MAXTH(21)
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP, CHEM(21)
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM3/ MINTH, MAXTH
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL            
      REAL*8 AR1(21), AR0(21), AR(21,21)
      REAL*8 AZ(21),  COMPR,  DALPHA_TR, ALPHAR_SUM
      REAL*8 DROX(21), DTRX(21),TMA, TMB, TMC        
       
      TMA = 0D0
      TMB = 0D0
      TMC = 0D0    
      DALPHA_TR=0D0
                             
      CALL MOLX(MOL, XIT, XI)   
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
      CALL ZOTPX(PRESS, TEMP, COMPR, XI)
      CALL ARROTPX(PRESS, TEMP, ALPHAR_SUM, XI)                       
      CALL D8OTPX(PRESS, TEMP, DALPHA_TR, XI) 
             
        DO 5 i=1,21
         DROX(i) = 0D0
         DTRX(i) = 0D0
         AR1(i) = 0D0      
         AR0(i)=0D0
         AZ(i) =0D0
         CHEM(i) = 0D0     
5       CONTINUE 
        
        DO 6 i=1,21
         DO 6 j=1,21
         AR(i,j) = 0D0               
6       CONTINUE 
    
!ALPHAR_0i for component i      
       DO 10 i=1,21
        DO 10 k = 1, K1POL(i)             
        AR0(i) =  AR0(i) + N_0(i,k) * DL ** D_0(i,k) 
     &   * TAU ** t_0(i,k) 
10     CONTINUE      
       
       DO 20 i=1,21
        DO 20 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        AR0(i) =AR0(i)+  N_0(i,k) * (DL**D_0(i,k))
     &   * (TAU ** t_0(i,k)) * DEXP(-DL**C_0(i,k))        
20     CONTINUE              


      DO 25 i=1,19
        IF (XI(i).NE.0D0) THEN
      AZ(i) = AZ(i) + ( DLOG(  DL / (INVRHOR*RHOC(i)) ) ) +(RSTAR/R)* (
     &   N0_0(i,1) + N0_0(i,2) * (TEMPC(i)/ TEMP) +
     &   N0_0(i,3) * DLOG (TEMPC(i)/ TEMP)   )   
        ENDIF         
 25   CONTINUE
 
       DO 30 i=1,19
         DO 30 k=MINTH(i)+1, MAXTH(i),2       
        IF (XI(i).NE.0D0) THEN
      AZ(i) = AZ(i)  + (RSTAR/R)*  (
     &   -N0_0(i,k)*DLOG(DABS(DCOSH(THETA_0(i,k)*TEMPC(i)/TEMP))) ) 
        ENDIF       
 30   CONTINUE
 
       DO 35 i=1,19
        DO 35 k=MINTH(i), MAXTH(i),2       
        IF (XI(i).NE.0D0) THEN
      AZ(i) = AZ(i)  + (RSTAR/R)* ( N0_0(i,k)* 
     &   DLOG(DABS(DSINH(THETA_0(i,k)*TEMPC(i)/TEMP))) ) 
        ENDIF                  
 35   CONTINUE
 
       DO 38 i=20,21
        IF (XI(i).NE.0D0) THEN
      AZ(i) = AZ(i)  +  (
     &   DLOG(  DL / (INVRHOR*RHOC(i)) ) )  + (RSTAR/R)* (
     &   N0_0(i,1) + N0_0(i,2) * TEMPC(i)/ TEMP +
     &   N0_0(i,3) * DLOG (TEMPC(i)/ TEMP) )  
        ENDIF
 38   CONTINUE    
         
!ALPHAR_ij for component i       
      DO 40  i=1,21
        DO 40  j=1,21   
         DO 40 k=1,K2POL(i,j)
        IF (j.NE.i) THEN   
         AR(i,j)= AR(i,j) + N(i,j,k) 
     &   * (DL**D(i,j,k)) * (TAU**T(i,j,k))
        ENDIF   
40     CONTINUE
    
       DO 50  i=1,21
        DO 50  j=1,21
          DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        IF (j.NE.i) THEN
         AR(i,j) = AR(i,j) + N(i,j,k) 
     &   *(DL**D(i,j,k))* (TAU**T(i,j,k))
     &   * DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0
     &  -BT(i,j,k)*(DL-GM(i,j,k)))
        ENDIF
50     CONTINUE
       
!Alpha_r derivative with respect to x for component i 

        DO 60 i=1,21
         DO 60 k=1,21
          IF (k.NE.i) THEN
           AR1(i) = AR1(i) + XI(k)*F(i,k) * AR(i,k) 
!        PRINT *, AX(i)
          ENDIF
60      CONTINUE 

        DO 70 i=1,21
           AR1(i) = AR0(i) + AR1(i)
!        PRINT *, AX(i)
70      CONTINUE   
     
!d (1/rhor)/ d xi for component i
        DO 90 i=1,21
         DO 90 k=1, i-1
          IF (BTV(k,i)**2D0 * XI(k)+XI(i).NE.0D0) THEN
          DROX(i)= DROX(i) + 2D0*BTV(k,i) * GMV(k,i) 
     &   * (   (1D0/8D0)*  (   (  1D0/RHOC(k)**(1D0/3D0) + 
     &   1D0/RHOC(i)**(1D0/3D0)   )**3D0   )) * ( XI(k)*( XI(k)+XI(i) )
     &   / ( BTV(k,i)**2D0 * XI(k)+XI(i) ) +
     &   XI(k)*XI(i) * (  1D0 - ( XI(k)+XI(i) ) / 
     &   ( BTV(k,i)**2D0 * XI(k)+XI(i) )  ) * 
     &   (1D0/ ( BTV(k,i)**2D0 * XI(k)+XI(i) ) )    )
         ENDIF                 
90      CONTINUE
        
        DO 100 i=1,21
         DO 100 k= i+1, 21 
          IF (BTV(i,k)**2D0 * XI(i)+XI(k).NE.0D0) THEN
          DROX(i)= DROX(i) + 2D0*BTV(i,k) * GMV(i,k) * 
     &   (   (1D0/8D0)*  (   (  1D0/RHOC(i)**(1D0/3D0) + 
     &   1D0/RHOC(k)**(1D0/3D0) )**3D0   ) ) * ( XI(k)*( XI(i)+XI(k) )
     &   / ( BTV(i,k)**2D0 * XI(i)+XI(k) ) +
     &   XI(i)*XI(k) * (  1D0 - ( XI(i)+XI(k) ) / 
     &   ( BTV(i,k)**2D0 * XI(i)+XI(k) )  ) * 
     &   (1D0/ ( BTV(i,k)**2D0 * XI(i)+XI(k) ) )    )
         ENDIF      
100     CONTINUE
     
        DO 110 i=1,21
         DROX(i)= DROX(i) + 2D0*XI(i)*(1D0/RHOC(i))
110     CONTINUE         

!d (Tr)/ d xi for component i
        DO 120 i=1,21
         DO 120 k=1, i-1
         IF (BTT(k,i)**2D0 * XI(k)+XI(i).NE.0D0) THEN
          DTRX(i)= DTRX(i) +  2D0*BTT(k,i) * GMT(k,i) 
     &   * (TEMPC(k)*TEMPC(i))**0.5D0  * (  XI(k)*( XI(k)+XI(i) )
     &   / ( BTT(k,i)**2D0 * XI(k)+XI(i) ) +
     &   XI(k)*XI(i) * (  1D0 - ( XI(k)+XI(i) ) / 
     &   ( BTT(k,i)**2D0 * XI(k)+XI(i) )  ) * 
     &   (1D0/ ( BTT(k,i)**2D0 * XI(k)+XI(i) ) )   ) 
         ENDIF                   
120      CONTINUE

        DO 130 i=1,21
         DO 130 k= i+1, 21
          IF (BTT(i,k)**2D0 * XI(i)+XI(k).NE.0D0) THEN         
          DTRX(i)= DTRX(i) +  2D0*BTT(i,k) * GMT(i,k) 
     &    * (TEMPC(i)*TEMPC(k))**0.5D0 *  (   XI(k)*( XI(i)+XI(k) )
     &   / ( BTT(i,k)**2D0 * XI(i)+XI(k) ) +
     &   XI(i)*XI(k) * (  1D0 - ( XI(i)+XI(k) ) / 
     &   ( BTT(i,k)**2D0 * XI(i)+XI(k) )  ) * 
     &   (1D0/ ( BTT(i,k)**2D0 * XI(i)+XI(k) ) )    )     
          ENDIF     
130     CONTINUE         

        DO 140 i=1,21
         DTRX(i)= DTRX(i) +  2D0*XI(i)*TEMPC(i)         
140     CONTINUE  

        DO 150 k=1,21
        TMA = TMA + XI(k) * (1D0/INVRHOR) * DROX(k)
        TMB = TMB + XI(k) * DTRX(k)
        TMC = TMC + XI(k) * AR1(k) 

150     CONTINUE 
 
        DO 180 i=1,21
         IF (XI(i).NE.0D0) THEN

        CHEM(i) = (
     &  ALPHAR_SUM + (COMPR-1D0) * (1D0 + (1D0/INVRHOR) 
     &  * DROX(i) - TMA) + DALPHA_TR * TAU * (1/TEMPR)*
     &  ( DTRX(i)-TMB) + AR1(i) - TMC   - DLOG(COMPR) +
     &  1 + DLOG(XI(i))   + AZ(i) ) * R*TEMP/1000d0  

!Chemical Potential: CHEM(i)
        ENDIF       
180     CONTINUE        
             
      END
      
      
!*******************************************************************************	  
	  SUBROUTINE FOTPX(PRESS, TEMP, FUG, CFUG, XI)
!*******************************************************************************
!COMPUTES FUGACITY AND FUGACITY COEFFICIENT OF COMPONENT i
!USES GERG 2008 MODEL	   
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j,k
      REAL*8 PRESS, DL, TAU, XI(21)
      REAL*8 TEMPR, INVRHOR, TEMP, FUG(21), CFUG(21)
      REAL*8 N0_0(21,7),THETA_0(21,7), RHOC(21),TEMPC(21),MCOMP(21)
	  REAL*8 F(21,21), N_0(21,24), C_0(21,24), D_0(21,24), T_0(21,24)
	  REAL*8 BTV(21,21), GMV(21,21), R, RSTAR
	  REAL*8 BTT(21,21), GMT(21,21), PRESSREF, TEMPREF
	  REAL*8 D(21,21,12), T(21,21,12), N(21,21,12), ETA(21,21,12)
	  REAL*8 EPS(21,21,12), BT(21,21,12), GM(21,21,12)	   
	  INTEGER  K1POL(21), K1EXP(21), K2POL(21,21), K2EXP(21,21)
	  COMMON /PARAM/ N0_0, THETA_0, F
	  COMMON /PARAM2/ MCOMP
	  COMMON /PARAM/ N_0, C_0, D_0, T_0
	  COMMON /PARAM1/ BTV, GMV, BTT, GMT
	  COMMON /PARAM5/ RHOC, TEMPC
	  COMMON /PARAM/ D, T, N, ETA 
	  COMMON /PARAM/ EPS, BT, GM
	  COMMON /CONST/ K1POL, K1EXP, K2POL, K2EXP
	  COMMON /CONST1/ PRESSREF, TEMPREF, R, RSTAR
	  COMMON /CONST2/ DL            
      REAL*8 AR0(21), AR(21)
      REAL*8 COMPR,  DALPHA_TR, ALPHAR_SUM
      REAL*8 DROX(21), DTRX(21),TMA, TMB, TMC        

      TMA = 0D0
      TMB = 0D0
      TMC = 0D0    
      DALPHA_TR=0D0
                            
      CALL MOLX(MOL, XIT, XI)      
      CALL REDFUNC(INVRHOR, TEMP, TEMPR, TAU, XI) 
      CALL ZOTPX(PRESS, TEMP, COMPR, XI) 
      CALL ARROTPX(PRESS, TEMP, ALPHAR_SUM, XI)     
      CALL D8OTPX(PRESS, TEMP, DALPHA_TR, XI) 
              
        DO 5 i=1,21
         DROX(i) = 0D0
         DTRX(i) = 0D0    
         AR0(i)=0D0
         FUG(i) = 0D0
         CFUG(i)=0D0        
5       CONTINUE 
   
        DO 6 i=1,21
         AR(i) = 0D0               
6       CONTINUE 
     
!ALPHAR_0i for component i      
       DO 10 i=1,21
        DO 10 k = 1, K1POL(i)             
        AR0(i) =  AR0(i) + N_0(i,k) * DL ** D_0(i,k) 
     &   * TAU ** t_0(i,k) 
10     CONTINUE      

       DO 20 i=1,21
        DO 20 k = K1POL(i)+1, K1POL(i)+K1EXP(i)
        AR0(i) =AR0(i)+  N_0(i,k) * (DL**D_0(i,k))
     &   * (TAU ** t_0(i,k)) * DEXP(-DL**C_0(i,k))        
20     CONTINUE              
    
!ALPHAR_ij for component i       
      DO 40  i=1,21
        DO 40  j=1,21   
         DO 40 k=1,K2POL(i,j)
        IF (j.NE.i) THEN   
         AR(i)= AR(i) + XI(j)*F(i,j)*N(i,j,k) 
     &   * (DL**D(i,j,k)) * (TAU**T(i,j,k))
        ENDIF    
40     CONTINUE
    
       DO 50  i=1,21
        DO 50  j=1,21
          DO 50 k=K2POL(i,j)+1, K2POL(i,j) + K2EXP(i,j)
        IF (j.NE.i) THEN
         AR(i) = AR(i) + XI(j)*F(i,j)*N(i,j,k) 
     &   *(DL**D(i,j,k))* (TAU**T(i,j,k))
     &   * DEXP(-ETA(i,j,k)*(DL-EPS(i,j,k))**2D0
     &  -BT(i,j,k)*(DL-GM(i,j,k)))
        ENDIF
50     CONTINUE
            
!Alpha_r derivative with respect to x for component i
        DO 70 i=1,21
           AR(i) = AR0(i) + AR(i)
70      CONTINUE   
     
!d (1/rhor)/ d xi for component i
        DO 90 i=1,21
         DO 90 k=1, i-1
          IF (BTV(k,i)**2D0 * XI(k)+XI(i).NE.0D0) THEN
          DROX(i)= DROX(i) + 2D0*BTV(k,i) * GMV(k,i) 
     &   * (   (1D0/8D0)*  (   (  1D0/RHOC(k)**(1D0/3D0) + 
     &   1D0/RHOC(i)**(1D0/3D0)   )**3D0   )) * ( XI(k)*( XI(k)+XI(i) )
     &   / ( BTV(k,i)**2D0 * XI(k)+XI(i) ) +
     &   XI(k)*XI(i) * (  1D0 - ( XI(k)+XI(i) ) / 
     &   ( BTV(k,i)**2D0 * XI(k)+XI(i) )  ) * 
     &   (1D0/ ( BTV(k,i)**2D0 * XI(k)+XI(i) ) )    )
         ENDIF                 
90      CONTINUE
        
        DO 100 i=1,21
         DO 100 k= i+1, 21 
          IF (BTV(i,k)**2D0 * XI(i)+XI(k).NE.0D0) THEN
          DROX(i)= DROX(i) + 2D0*BTV(i,k) * GMV(i,k) * 
     &   (   (1D0/8D0)*  (   (  1D0/RHOC(i)**(1D0/3D0) + 
     &   1D0/RHOC(k)**(1D0/3D0) )**3D0   ) ) * ( XI(k)*( XI(i)+XI(k) )
     &   / ( BTV(i,k)**2D0 * XI(i)+XI(k) ) +
     &   XI(i)*XI(k) * (  1D0 - ( XI(i)+XI(k) ) / 
     &   ( BTV(i,k)**2D0 * XI(i)+XI(k) )  ) * 
     &   (1D0/ ( BTV(i,k)**2D0 * XI(i)+XI(k) ) )    )
         ENDIF    
100     CONTINUE
     
        DO 110 i=1,21
         DROX(i)= DROX(i) + 2D0*XI(i)*(1D0/RHOC(i))
110     CONTINUE         
             
!d (Tr)/ d xi for component i
        DO 120 i=1,21
         DO 120 k=1, i-1
         IF (BTT(k,i)**2D0 * XI(k)+XI(i).NE.0D0) THEN
          DTRX(i)= DTRX(i) +  2D0*BTT(k,i) * GMT(k,i) 
     &   * (TEMPC(k)*TEMPC(i))**0.5D0  * (  XI(k)*( XI(k)+XI(i) )
     &   / ( BTT(k,i)**2D0 * XI(k)+XI(i) ) +
     &   XI(k)*XI(i) * (  1D0 - ( XI(k)+XI(i) ) / 
     &   ( BTT(k,i)**2D0 * XI(k)+XI(i) )  ) * 
     &   (1D0/ ( BTT(k,i)**2D0 * XI(k)+XI(i) ) )   ) 
         ENDIF                   
120      CONTINUE

        DO 130 i=1,21
         DO 130 k= i+1, 21
          IF (BTT(i,k)**2D0 * XI(i)+XI(k).NE.0D0) THEN         
          DTRX(i)= DTRX(i) +  2D0*BTT(i,k) * GMT(i,k) 
     &    * (TEMPC(i)*TEMPC(k))**0.5D0 *  (   XI(k)*( XI(i)+XI(k) )
     &   / ( BTT(i,k)**2D0 * XI(i)+XI(k) ) +
     &   XI(i)*XI(k) * (  1D0 - ( XI(i)+XI(k) ) / 
     &   ( BTT(i,k)**2D0 * XI(i)+XI(k) )  ) * 
     &   (1D0/ ( BTT(i,k)**2D0 * XI(i)+XI(k) ) )    )     
          ENDIF     
130     CONTINUE         

        DO 140 i=1,21
         DTRX(i)= DTRX(i) +  2D0*XI(i)*TEMPC(i)         
140     CONTINUE  
           
        DO 150 k=1,21
         TMA = TMA + XI(k) * (1D0/INVRHOR) * DROX(k)
         TMB = TMB + XI(k) * DTRX(k)
         TMC = TMC + XI(k) * AR(k) 
150     CONTINUE 
!THIS CALCULATES FUGACITY f_i           
        DO 180 i=1,21
         IF (XI(i).NE.0D0) THEN
       FUG(i) = XI(i)* (DL / INVRHOR) * R* TEMP * DEXP (
     &  ALPHAR_SUM + (COMPR-1D0) * (1D0 + (1D0/INVRHOR) 
     &  * DROX(i) - TMA) + DALPHA_TR * TAU * (1/TEMPR)*
     &  ( DTRX(i)-TMB) + AR(i) - TMC  ) / 1000D0
           
!THIS CALCULATES Fi_i=DEXP(LOG(Fi_i)) WHICH ALSO EQUALS f_i / (x_i * P)
        CFUG(i) = DEXP(
     &  ALPHAR_SUM + (COMPR-1D0) * (1D0 + (1D0/INVRHOR) 
     &  * DROX(i) - TMA) + DALPHA_TR * TAU * (1/TEMPR)*
     &  ( DTRX(i)-TMB) + AR(i) - TMC   - DLOG(COMPR)     )
     
        ENDIF       
180     CONTINUE  
      END
!*******************************************************************************	

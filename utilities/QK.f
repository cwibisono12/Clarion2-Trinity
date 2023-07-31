      SUBROUTINE QK(E,QD2,QD4,R, D, T )
C
C
C		THIS ROUTINE CALCULATES THE ANGULAR DISTRIBUTIONS
C		ATTENUATION COEFFICIENTS AS DEFINED BY ROSE FOR
C		CYLINDRICAL GE(LI) DETECTORS
C
		EMEV=E/1000
		ELN=LOG(EMEV)
		EL1=ELN
		EL2=ELN**2
		EL3=EL1*EL2
		EL4=EL2**2
		EL5=EL4*EL1
		TLN=-1.1907-.5372*EL1-.043*EL2+.0218*EL3+
     1           .0765*EL4+.0095*EL5
		TAU=EXP(TLN)
C
		Z1=R/(D+T)
		Z2=R/D
		ALPHA=ATAN(Z1)
		GAMMA=ATAN(Z2)
C
		BL=0
		BU=ALPHA
		DELX1=(BU-BL)/1000
		SUM1=0
		SUM2=0
		SUM3=0
C
		DO 100 I=0,1000
C
		   IF(I.NE.0) GOTO 10
		   A=1.0
		   GO TO 50
10		   IF(I.NE.1000) GOTO 20
		   A=1.0
		   GO TO 50
20		   J=MOD(I,2)
		   IF(J.EQ.0) GOTO 30
		   A=4.0
		   GOTO 50
30                 A=2.0
50		   BETA=BL+I*DELX1
C
		   COSB=COS(BETA)
		   SINB=SIN(BETA)
		   SECB=1.0/COSB
		   C2=COSB**2
		   C4=COSB**4
		   FAC1=-1*TAU*T*SECB
		   EX1=EXP(FAC1)
C
		  TERM1=.5*(3*C2-1)*(1-EX1)*SINB*A*DELX1
		  TERM2=.125*A*(35*C4-30*C2+3)*(1-EX1)*SINB*DELX1
		TERM3=A*(1-EX1)*SINB*DELX1
C
		  SUM1=SUM1+TERM1
		  SUM2=SUM2+TERM2
		  SUM3=SUM3+TERM3
C
100		CONTINUE
C
		ANS1=SUM1/3
		ANS2=SUM2/3
		ANS3=SUM3/3
C
C
		BL=ALPHA
		BU=GAMMA
		DELX2=(BU-BL)/1000
		SUM4=0
		SUM5=0
		SUM6=0
C
		DO 105 I=0,1000
C
		   IF(I.NE.0) GO TO 60
		   A=1.0
		   GOTO 90
60		IF(I.NE.1000) GOTO 70
		   A=1.0
		   GOTO 90
70		   J2=MOD(I,2)
		   IF (J2.EQ.0) GOTO 80
		   A=4.0
		   GOTO 90
80		   A=2.0
90		   BETA=BL+I*DELX2
C
		  COSB=COS(BETA)
		  SINB=SIN(BETA)
		  SECB=1.0/COSB
		  CSCB=1.0/SINB
		  FAC2=-1*TAU*(R*CSCB-D*SECB)
		  EX2=EXP(FAC2)
		  C2=COSB**2
		  C4=COSB**4
C
		TERM4=.5*A*(3*C2-1)*(1-EX2)*SINB*DELX2
		TERM5=.125*A*(35*C4-30*C2+3)*(1-EX2)*SINB*DELX2
		TERM6=A*(1-EX2)*SINB*DELX2
C
		   SUM4=SUM4+TERM4
		   SUM5=SUM5+TERM5
		   SUM6=SUM6+TERM6
C
105		CONTINUE
C
	        ANS4=SUM4/3
		ANS5=SUM5/3
		ANS6=SUM6/3
C
C		CALCULATE Q2 AND Q4
C
		QD2=(ANS1+ANS4)/(ANS3+ANS6)
		QD4=(ANS2+ANS5)/(ANS3+ANS6)
C
C	      WRITE(*,*) 'QD2:', QD2
C            WRITE(*,*) 'QD4:', QD4
C
		RETURN
C
		END
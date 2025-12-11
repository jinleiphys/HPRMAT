c------------------------------------------------------------------------------
c angular_momentum.f - Wigner 3j and 6j coefficient calculations
c Extracted from Pierre Descouvemont's R-matrix code
c------------------------------------------------------------------------------

c*CLEBS - Clebsch-Gordan coefficients via 3j symbols
      SUBROUTINE CLEBS (L1,L2,L3,M1,M2,M3,Q)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FT(0:1000)
      SAVE
      DATALMEM,(FT(I),I=0,10)/9,2*1.0d0,2.0d0,6.0d0,24.0d0,
     1 120.0d0,720.0d0,5040.0d0,40320.0d0,362880.0d0,3628800.0d0/
      IS=(L1+L2+M1-M2)/2
      K3=-M3
      Q1=L3+1
   12 Q=0
      I1=L1
      I2=L2
      I3=L3
      K1=M1
      K2=M2
      IF(K1+K2+K3.NE.0)RETURN
      L=I1+I2+I3
      IF(MOD(L,2).EQ.1)GOTO8
      L=L/2
      IF(L.LE.LMEM)GOTO6
      DO10I=LMEM,L
   10 FT(I+1)=(I+1)*FT(I)
      LMEM=L
    6 J1=ABS(K1)
      J2=ABS(K2)
      J3=ABS(K3)
      IF(I1.LT.ABS(I2-I3).OR.I1.GT.I2+I3)RETURN
      IF(J1+J2.EQ.0)GOTO11
      J1=I1-J1
      J2=I2-J2
      J3=I3-J3
      IF(J1)8,2,1
    2 IF(I1.NE.0)GOTO1
      IF(J2.LT.0)GOTO8
   13 IF(J3.LT.0)GOTO8
    4 Q=SQRT(Q1/(I2+1))
      IS=IS+(I2-K2)/2
      IF(MOD(IS,2).EQ.1)Q=-Q
      RETURN
    1 IF(J2.GT.J1)GOTO3
      IF(J2.LT.0)GOTO8
      IS=IS+L
      J1=J2
      K1=K2
      K2=M1
      I1=I2
      I2=L1
      IF(I1.EQ.0)GOTO13
    3 IF(J3.GT.J1)GOTO5
      IF(J3.LT.0)GOTO8
      IS=IS+L
      J1=K3
      K3=K1
      K1=J1
      I3=I1
      I1=L3
      J1=J3
      IF(I1.EQ.0)GOTO4
    5 IF(K1.GE.0)GOTO9
      K1=-K1
      K2=-K2
      K3=-K3
      IS=IS+L
    9 CONTINUE
      Q1=Q1*FT(L-I3)/FT(L-I1)/FT(L-I2)/FT(L+1)
      I2=(I2+K2)/2
      I3=(I3+K3)/2
      K2=I2-K2
      K3=I3-K3
      J1=J1/2
      I1=J1+K1
      J2=I3-K2
      J3=MAX(J2,0)
      IS=IS+I1+K2
      X=0
      DO7I=J3,J1
    7 X=-X+FT(I1+I)*FT(I2+I3-I)/FT(J1-I)/FT(I3-I)/FT(I-J2)/FT(I)
      Q=X*       SQRT(Q1*FT(J1)*FT(K2)*FT(K3)*FT(I3)/FT(I1)/FT(I2))
      IF(MOD(IS,2).EQ.1)Q=-Q
      RETURN
    8 Q=0
      PRINT1010,L1,L2,L3,M1,M2,M3
 1010 FORMAT('ERROR 3J',2(3X,3I3))
      RETURN
   11 IF(MOD(L,2).EQ.1)RETURN
      I1=L-I1
      I2=L-L2
      I3=L-L3
      Q=SQRT(FT(I1)*FT(I2)*FT(I3)/FT(L+1)*Q1)
      I1=I1/2
      I2=I2/2
      I3=I3/2
      L =L/2
      Q=Q*FT(L )/FT(I1)/FT(I2)/FT(I3)
      IF(MOD(L +IS,2).EQ.1)Q=-Q
      RETURN
      ENTRY TROISJI (L1,L2,L3,M1,M2,M3,Q)
      ENTRY TROISJ (L1,L2,L3,M1,M2,M3,Q)
      K3=M3
      Q1=1
      IS=0
      GOTO12
      END

c*SIXJ - Wigner 6j coefficient
      SUBROUTINE SIXJ(J1,J2,J3,L1,L2,L3,Q)
c Computes the 6j coefficient
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION M(7),M1(4),M2(4),M3(4),FT(0:1000)
      SAVE
      DATA LMEM,(FT(I),I=0,10)/9,2*1.0d0,2.0d0,6.0d0,24.0d0,
     1 120.0d0,720.0d0,5040.0d0,40320.0d0,362880.0d0,3628800.0d0/
      ENTRYSIXJI(J1,J2,J3,L1,L2,L3,Q)
      I1=J1
      I2=J2
      I3=J3
      K1=L1
      K2=L2
      K3=L3
      IS=0
   24 Q=0
      M(1)=I1+I2+I3
      M(2)=I1+K2+K3
      M(3)=K1+I2+K3
      M(4)=K1+K2+I3
      DO 17 I=1,4
      IF(MOD(M(I),2).EQ.1) GO TO 8
   17 CONTINUE
      L=MAX(I1+I2+K1+K2,I1+I3+K1+K3,I2+I3+K2+K3)
      L=L/2
      IF(L.LE.LMEM) GO TO 6
      DO 10 I=LMEM,L
   10 FT(I+1)=FT(I)*(I+1)
      LMEM=L
    6 IF(I1.LT.ABS(I2-I3).OR.I1.GT.I2+I3) RETURN
      IF(I1.LT.ABS(K2-K3).OR.I1.GT.K2+K3) RETURN
      IF(K1.LT.ABS(I2-K3).OR.K1.GT.I2+K3) RETURN
      IF(K1.LT.ABS(K2-I3).OR.K1.GT.K2+I3) RETURN
      IF(I1) 8,2,1
    2 IF(I2.LT.0) GO TO 8
    9 IF(I3.LT.0) GO TO 8
   14 IF(K1.LT.0) GO TO 8
   19 IF(K2.LT.0) GO TO 8
   23 IF(K3.LT.0) GO TO 8
   27 Q=SQRT(1.0d0/(I2+1)/(K2+1))
      IS=(I2+K2+K1)/2+IS
      IF(MOD(IS,2).EQ.1) Q=-Q
      RETURN
    1 IF(I1.GT.1) GO TO 3
      IF(I2.LT.0) RETURN
   12 IF(I3.LT.0) RETURN
   16 IF(K1.LT.0) RETURN
   21 IF(K2.LT.0) RETURN
   25 IF(K3.LT.0) RETURN
   28 IF(I2.LT.I3) GO TO 4
      IC=I2
      I2=I3
      I3=IC
      IC=K2
      K2=K3
      K3=IC
    4 IF(K2.GT.K3) GO TO 5
      I11=I1+K1+I2-K2
      I11=I11/2
      I12=I11-I2+K2
      Q=SQRT(I11*I12*1.0d0/I3/(I3+1)/K3/(K3+1))
      IS =I11+K2+IS
      IF(MOD(IS,2).EQ.1) Q=-Q
      RETURN
    5 I11=K3-K1+I2
      I11=I11/2+1
      I12=I11+K1+1
      Q=SQRT(I11*I12*1.0d0/I3/(I3+1)/K2/(K2+1))
      IS =I12-1+IS
      IF(MOD(IS ,2).EQ.1) Q=-Q
      RETURN
    3 IF(I2.GE.I1) GO TO 7
      IF(I2.LT.0) GO TO 8
      IC=I2
      I2=I1
      I1=IC
      IC=K1
      K1=K2
      K2=IC
      IF(I1.EQ.0) GO TO 9
      IF(I1.EQ.1) GO TO 12
    7 IF(I3.GE.I1) GO TO 13
      IF(I3.LT.0) GO TO 8
      IC=I3
      I3=I1
      I1=IC
      IC=K3
      K3=K1
      K1=IC
      IF(I1.EQ.0) GO TO 14
      IF(I1.EQ.1) GO TO 16
   13 IF(K1.GE.I1) GO TO 18
      IF(K1.LT.0) GO TO 8
      IC=K1
      K1=I1
      I1=IC
      IC=K2
      K2=I2
      I2=IC
      IF(I1.EQ.0) GO TO 19
      IF(I1.EQ.1) GO TO 21
   18 IF(K2.GE.I1) GO TO 22
      IF(K2.LT.0) GO TO 8
      IC=K2
      K2=I1
      I1=IC
      IC=K1
      K1=I2
      I2=IC
      IF (I1.EQ.0) GO TO 23
      IF(I1.EQ.1) GO TO 25
   22 IF(K3.GE.I1) GO TO 26
      IF(K3.LT.0) GO TO 8
      IC=K3
      K3=I1
      I1=IC
      IC=K1
      K1=I3
      I3=IC
      IF(I1.EQ.0) GO TO 27
      IF(I1.EQ.1) GO TO 28
   26 M1(4)=I3
      M1(1)=I3
      M1(3)=K3
      M1(2)=K3
      M2(2)=I1
      M2(1)=I1
      M2(4)=K1
      M2(3)=K1
      M3(3)=I2
      M3(1)=I2
      M3(4)=K2
      M3(2)=K2
      M(1)=I1+I2+I3
      M(2)=I1+K2+K3
      M(3)=K1+I2+K3
      M(4)=K1+K2+I3
      Q1=1
      DO 11 I=1,4
      M(I)=M(I)/2
   11 Q1=FT(M(I)-M1(I))*FT(M(I)-M2(I))*FT(M(I)-M3(I))*Q1/FT(M(I)+1)
      Q1=SQRT(Q1)
      M1(1)=I1+K1
      M1(2)=I2+K2
      M1(3)=I3+K3
      IC=M1(1)+M1(2)
      M(5)=IC/2
      IC=M1(2)+M1(3)
      M(6)=IC/2
      IC=M1(1)+M1(3)
      M(7)=IC/2
      MAXZ=MIN(M(5),M(6),M(7))
      MINZ=MAX(M(1),M(2),M(3),M(4))
      X=0
      DO 15 I=MINZ,MAXZ
      Q2=1
      DO 20 J=1,7
      IJ=I-M(J)
      IF(J.GT.4) IJ=-IJ
   20 Q2=Q2*FT(IJ)
      Q2=FT(I+1)/Q2
   15 X=-X+Q2
      Q=X*Q1
      IS=MAXZ+IS
      IF(MOD(IS,2).EQ.1) Q=-Q
      RETURN
    8 PRINT 1010,J1,J2,J3,L1,L2,L3
 1010 FORMAT('ERROR 6J',2(3X,3I3))
      RETURN
      ENTRY RACAH(J1,J2,J3,L1,L2,L3,Q)
      IS=(J1+J2+J3+L1)/2
      I1=J1
      I2=J2
      I3=L2
      K1=L1
      K2=J3
      K3=L3
      GO TO 24
      END

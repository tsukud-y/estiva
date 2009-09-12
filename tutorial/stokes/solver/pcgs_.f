*
      SUBROUTINE PCGS(D,A,IA,N,N1,NL,B,EPS,ITR,S,X,DD,P,Q,R,R0,E,H,W,
     *                M,IER)
***********************************************************************
*  CONJUGATE GRADIENT SQURED METHOD WITH INCOMPLETE LU DECOMPOSITION. *
*                                                                     *
*  PARAMETERS: SAME AS ROUTINE PCG EXCEPT A AND IA.                   *
*   ON ENTRY:                                                         *
*     A      NON-ZERO ELEMENTS OF THE LOWER TRIANGULAR PART OF THE    *
*            MATRIX A INTO 1-ST TO NL-TH POSITION OF THE ARRAY A.     *
*            NON-ZERO ELEMENTS OF THE UPPER TRIANGULAR PART OF THE    *
*            MATRIX A INTO NL+1-ST TO 2*NL-TH POSITION OF THE ARRAY A.*
*     IA     COLUMN INDEX OF CORRESPONDING ELEMENT IN THE ARRAY A.    *
*   OTHERS:  WORKING PARAMETERS.                                      *
*                                                                     *
*  COPYRIGHT:     TSUTOMU OGUNI       FEB. 1 1993      VER. 2         *
***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(N1,2*NL), B(N), X(0:N), W(0:N), P(0:N), D(N), M(2*N)
     *  , Q(0:N), R(0:N), DD(0:N), R0(N), E(N), H(N), IA(N1,2*NL)
C
      IER = 0
      IF ((N1 .LT. N) .OR. (S .LT. 0.0)) THEN
       WRITE(*,*) '(SUBR. PCGS) INVALID ARGUMENT. ',N,N1,NL,S
       IER = 2
       RETURN
      ENDIF
C
      TH = 1.0D0
      IF (S .GT. 0.0 .AND. S .LT. 1.0) THEN
       TH = S
       S = 1.0D0
      ENDIF
      DO 5 I=1,2*N
    5  M(I) = 0
      DO 6 I=1,N
       DD(I) = 0.0D0
       DO 7 J=1,NL
        IF ( IA(I,J) .NE. 0) M(I) = M(I) + 1
    7  CONTINUE
       DO 8 J=NL+1,2*NL
        IF (IA(I,J) .NE. 0) M(I+N) = M(I+N) + 1
    8  CONTINUE
    6 CONTINUE   
       DD(0) = 0.0D0
       X(0) = 0.0D0
       P(0) = 0.0D0
       Q(0) = 0.0D0
       R(0) = 0.0D0
       W(0) = 0.0D0
C  INCOMPLETE CHOLESKY DECOMPOSITION
      IF (S .NE. 0.0) THEN   
       DD(1) = 1.0D0 / (S * D(1))
       DO 15 I=2,N
        SS = S * D(I)
        DO 16 K=1,M(I)
         SW = 0.0D0
         NN = IA(I,K) 
         IF (NN .NE. 0) THEN  
          DO 17 J=NL+1,NL+M(NN+N)
           IF (IA(NN,J) .NE. I) THEN
            SW = SW + A(NN,J) * TH
           ELSE
            SW = SW + A(NN,J)
           ENDIF
   17     CONTINUE
         ENDIF
   16    SS = SS - A(I,K) * SW * DD(NN)       
   15   DD(I) = 1.0D0 / SS
      ELSE
       DD(1) = 1.0D0 / D(1)
       DO 20 I=2,N
        SS = D(I)
        DO 21 K=1,M(I)
         NN = IA(I,K)
         DO 22 J=NL+1,NL+M(NN+N)
          IF(IA(NN,J) .EQ. I) SS = SS - A(I,K) * A(NN,J) * DD(NN)
   22    CONTINUE      
   21   CONTINUE
   20   DD(I) = 1.0D0 / SS
      ENDIF
      DO 30 I=1,N
       Q(I) = D(I) * X(I)
       DO 32 J=1,M(I)
   32   Q(I) = Q(I) + A(I,J) * X(IA(I,J))
       DO 34 J=NL+1,NL+M(I+N)
   34   Q(I) = Q(I) + A(I,J) * X(IA(I,J))
   30 CONTINUE 
      DO 40 I=1,N
   40  R(I) = B(I) - Q(I)
      DO 50 I=1,N
       DO 52 J=1,M(I)
   52   R(I) = R(I) - A(I,J) * R(IA(I,J))
   50  R(I) = R(I) * DD(I) 
      DO 60 I=N,1,-1
       SW = 0.0D0
       DO 62 J=NL+1,NL+M(I+N)
   62   SW = SW + A(I,J) * R(IA(I,J))
   60  R(I) = R(I) - DD(I) * SW
      C1 = 0.0D0
      DO 70 I=1,N
       R0(I) = R(I)
       P(I) = R(I)
       E(I) = R(I)
   70  C1 = C1 + R(I) * R(I)
C  ITERATION PHASE
      DO 200 K=1,ITR
       DO 80 I=1,N
        Q(I) = D(I) * P(I)
        DO 85 J=1,M(I)
   85    Q(I) = Q(I) + A(I,J) * P(IA(I,J))
        DO 87 J=NL+1,NL+M(N+I)
   87    Q(I) = Q(I) + A(I,J) * P(IA(I,J))
   80  CONTINUE
       DO 90 I=1,N
        DO 95 J=1,M(I)
   95    Q(I) = Q(I) - A(I,J) * Q(IA(I,J))
   90   Q(I) = DD(I) * Q(I)
       DO 100 I=N,1,-1
        SW = 0.0D0
        DO 105 J=NL+1,NL+M(I+N)
  105    SW = SW + A(I,J) * Q(IA(I,J))
  100   Q(I) = Q(I) - DD(I) * SW
       C2 = 0.0D0
       DO 110 I=1,N
  110   C2 = C2 + Q(I) * R0(I)
       IF (C2 .EQ. 0.0) THEN
        IER = 3
        ITR = K
        GO TO 300
       ENDIF
       ALPHA = C1 / C2
       C3 = 0.0D0
       X1 = 0.0D0
       X2 = 0.0D0
       DO 120 I=1,N
        H(I) = E(I) - ALPHA * Q(I)
  120  CONTINUE
       DO 130 I=1,N
  130   W(I) = E(I) + H(I)
       DO 140 I=1,N
        Q(I) = D(I) * W(I)
        DO 142 J=1,M(I)
  142    Q(I) = Q(I) + A(I,J) * W(IA(I,J))
        DO 144 J=NL+1,NL+M(I+N)
  144    Q(I) = Q(I) + A(I,J) * W(IA(I,J))
  140  CONTINUE
       DO 150 I=1,N
        DO 155 J=1,M(I)
  155    Q(I) = Q(I) - A(I,J) * Q(IA(I,J))
  150   Q(I) = DD(I) * Q(I)
       DO 160 I=N,1,-1
        SW = 0.0D0
        DO 165 J=NL+1,NL+M(I+N)
  165    SW = SW + A(I,J) * Q(IA(I,J))
  160   Q(I) = Q(I) - DD (I) * SW   
       DO 170 I=1,N
        Y = X(I)
        R(I) = R(I) - ALPHA * Q(I)
        X(I) = X(I) + ALPHA * W(I)
        C3 = C3 + R(I) * R0(I)
        X1 = X1 + Y * Y
  170   X2 = X2 + (X(I) - Y)**2
       IF (X1 .NE. 0.0) THEN
        RES = DSQRT(X2 / X1) 
        IF (RES .LE. EPS) THEN
         ITR = K
         IER = 0
         GO TO 300
        ENDIF
       ENDIF 
       IF (C1 .EQ. 0.0) THEN
        IER = 4
        ITR = K
        GO TO 300
       ENDIF
       BETA = C3 / C1
       C1 = C3
       DO 180 I=1,N
        E(I) = R(I) + BETA * H(I)
        P(I) = E(I) + BETA * (H(I) + BETA * P(I))
  180  CONTINUE
C
  200 CONTINUE
      IER = 1
      WRITE(*,*) '(SUBR. MLUCGS) NO CONVERGENCE. '
  300 CONTINUE
      EPS = RES
      IF (TH .NE. 1.0D0) S = TH
      RETURN
C  END OF PCGS
      END

*
      SUBROUTINE PCG(D,A,IA,N,N1,NL,B,EPS,ITR,S,X,DD,P,Q,R,M,IER)
**********************************************************************
*  ICCG AND MICCG METHOD FOR FINITE ELEMENT METHOD.                  *
*                                                                    *
*   PARAMETERS:                                                      *
*    ON ENTRY:                                                       *
*      D      1-DIM. ARRAY CONTAINING DIAGONAL ELEMENTS OF MATRIX A. *
*      A      2-DIM. ARRAY CONTAINING NON-AERO ELEMENTS OF LOWER     *
*             PART OF THE MATRIX A EXCEPT DIAGONALS.                 *
*      IA     2-DIM. ARRAY CONTAINING THE COLUMN INDEX OF ELEMENTS   *
*             IN TH ARRAY A.                                         *
*      N      THE ORDER OF THE MATRIX A.                             *
*      N1     THE LEADING DIMENSION OF THE ARRAY A.                  *
*      NL     MORE THAN MAXIMUM NUMBER OF NON-ZERO ELEMENTS IN EACH  *
*             ROW OF THE ARRAY A.                                    *
*      B      1-DIM. ARRAY CONTAINING RIGHT HAND SIDE VECTOR.        *
*      EPS    THE TOLERANCE FOR CONVERGENCE.                         *
*      ITR    MAXIMUM NUMBER OF ITERATIONS.                          *
*      S      THE SIGMA VALUE WHICH SPECIFIES A METHOD TO BE USED.   *
*   ON RETURN:                                                       *
*      X      1-DIM. ARRAY CONTAINING THE SOLUTION VECTOR.           *
*      EPS    RELATIVE ERROR AT RETURN.                              *
*      ITR    NUMBER OF ITERATIONS AT RETURN.                        *
*      IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.               *
*   OTHERS:   WORKING PARAMETERS:                                    *
*                                                                    *
*  COPYRIGHT:   TSUTOMU OGUNI    FEB. 1 1993   VER. 2                *
**********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION D(N), A(N1,2*NL), IA(N1,2*NL), B(N), X(0:N)
     *  , DD(0:N), R(N), P(0:N), Q(0:N), M(2*N)
C
      IER = 0
      IF (N1 .LT. N .OR. S .LT. 0.0) THEN
       WRITE(*,*) '(SUBR. PCG) INVALID ARGUMENT. ', N1, N, NL, S
       IER = 2
       RETURN
      ENDIF
C  INITIALIZATION
      TH = 1.0D0
      IF (S .GT. 0.0 .AND. S .LT. 1.0) THEN
       TH = S
       S = 1.0D0
      ENDIF
      DO 5 I=1,2*N
    5  M(I) = 0
      DO 10 I=1,N
       DD(I) = 0.0D0
       DO 10 J=1,NL
        K = IA(I,J)
        IF (K .NE. 0) THEN
         M(N+K) = M(N+K) + 1
         IA(K,M(N+K)+NL) = I
         A(K,M(N+K)+NL) = A(I,J)
         M(I) = M(I) + 1
        ENDIF 
   10 CONTINUE
      X(0) = 0.0D0
      DD(0) = 0.0D0
      P(0) = 0.0D0
      Q(0) = 0.0D0
      IF (S .NE. 0.0) THEN
       DO 12 I=1,N
        W = D(I) * S
        DO 14 K=1,M(I) 
         NN = IA(I,K)
         SS = A(I,K)
         IF (NN .NE. 0) THEN
          DO 16 J=NL+1,NL+M(NN+N)
           IF (IA(NN,J) .NE. I) THEN
            SS = SS + A(NN,J) * TH
           ENDIF
   16     CONTINUE       
         ENDIF
   14    W = W - A(I,K) * SS * DD(NN)
   12   DD(I) = 1.0D0 / W
      ELSE
       DO 30 I=1,N
        SS = D(I)
        DO 20 K=1,M(I)
   20    SS = SS - A(I,K)**2 * DD(IA(I,K))
   30   DD(I) = 1.0D0 / SS
      ENDIF
      DO 50 I=1,N
       Q(I) = D(I) * X(I)
       DO 40 J=1,M(I)
   40   Q(I) = Q(I) + A(I,J) * X(IA(I,J))
       DO 45 J=NL+1,NL+M(I+N)
   45   Q(I) = Q(I) + A(I,J) * X(IA(I,J))
   50 CONTINUE
      DO 60 I=1,N
   60  R(I) = B(I) - R(I)
      DO 80 I=1,N
       P(I) = R(I)
       DO 70 K=1,M(I)
   70   P(I) = P(I) - A(I,K) * P(IA(I,K))
   80  P(I) = DD(I) * P(I)
      DO 100 I=N,1,-1
       SS = 0.0D0
       DO 90 K=NL+1,NL+M(I+N)
   90   SS = SS +A(I,K) * P(IA(I,K))
  100  P(I) = P(I) - DD(I) * SS
      C1 = 0.0D0
      DO 110 I=1,N
  110  C1 = C1 + R(I) * P(I)
C  ITERATION PHASE
      DO 250 L=1,ITR
       DO 140 I=1,N
        Q(I) = D(I) * P(I)
        DO 130 J=1,M(I)
  130    Q(I) = Q(I) + A(I,J) * P(IA(I,J))
        DO 135 J=NL+1,NL+M(I+N)
  135       Q(I) = Q(I) + A(I,J) * P(IA(I,J))
  140  CONTINUE
      C2 = 0.0D0
      DO 150 I=1,N
  150  C2 = C2 + P(I) * Q(I)
      IF (C2 .EQ. 0.0) THEN
       IER = 3
       ITR = L
       GO TO 300
      ENDIF
      ALPHA = C1 / C2
      X1 = 0.0D0
      X2 = 0.0D0
      DO 160 I=1,N
       Y = X(I)
       X(I) = X(I) + ALPHA * P(I)
       R(I) = R(I) - ALPHA * Q(I)
       X1 = X1 + Y * Y
  160  X2 = X2 + (X(I) - Y)**2
      DO 190 I=1,N
       Q(I) = R(I)
       DO 180 K=1,M(I)
  180   Q(I) = Q(I) - A(I,K) * Q(IA(I,K))
  190  Q(I) = DD(I) * Q(I)
      DO 220 I=N,1,-1
       DO 220 K=NL+1,NL+M(I+N)
  220   Q(I) = Q(I) - DD(I) * A(I,K) * Q(IA(I,K))
      IF (X1 .NE. 0.0) THEN
       RES = DSQRT(X2 / X1)
       IF (RES .LT. EPS) THEN
        ITR = L
        EPS = RES
        IER = 0
        IF (TH .NE. 1.0) S = TH
        RETURN
       ENDIF
      ENDIF
      C3 = 0.0D0
      DO 230 I=1,N
  230  C3 = C3 + R(I) * Q(I)
      IF (C1 .EQ. 0.0) THEN
       ITR = L
       IER = 4
       GO TO 300
      ENDIF
      BETA = C3 / C1
      C1 = C3
      DO 240 I=1,N
       P(I) = Q(I) + BETA * P(I)
  240 CONTINUE
C
  250 CONTINUE
      IER = 1
      WRITE(*,*) '(SUBR. PCG) NO CONVERGENCE.'
  300 CONTINUE
      EPS = RES
      IF (TH .NE. 1.0) S = TH
      RETURN
C  END OF PCG
      END

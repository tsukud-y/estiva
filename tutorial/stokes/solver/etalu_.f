*
      SUBROUTINE ETALU(A,IA,L,L2,B,NWK,M,EPS,NUM,X,E,IE,WE,NE,IP,IS,
     *                 IER)
***********************************************************************
*  GAUSS METHOD FOR A NON-SYMMETRIC SPARSE MATRIX.                    *
*                                                                     *
*  PARAMETERS:                                                        *
*   ON ENTRY:                                                         *
*     A      THE ARRAY WHICH CONTAINS NON-ZERO ELEMENTS IN COLUMN-WISE*
*     IA     THE ARRAY WHICH HAS CORRESPONDING ROW INDEX WITH ARRAY A.*
*     L      THE LEADING DIMENSION OF THE ARRAY A.                    *
*     L2     THE LEADING DIMENSION OF THE ARRAY E.                    *
*     B      THE RIGHT HAND SIDE VECTOR.                              *
*     NWK    ACCUMULATED SUM OF NON-ZERO ELEMENTS IN EACH COLUMN OF A.*
*     M      THE ORDER OF THE MATRIX A.                               *
*     EPS    THE TOLERANCE FOR PIVOTAL ELEMENTS.                      *
*   ON RETURN:                                                        *
*     X      THE SOLUTION VECTOR.                                     *
*     E      THE ARRAY WHICH CONTAINS NON ZERO ELEMENTS IN ETA-VECTORS*
*            IN COLUMN-WISE.                                          *
*     IE     THE ARRAY WHICH HAS CORRESPONDING ROW INDEX WITH ARRAY E.*
*     NETA   ACCUMULATED SUM OF NON-ZERO ELEMENTS IN EACH COLUMN OF E.*
*            EACH COLUMN HAS TWO ENTRIES FOR U-ETA AND L-ETA.         *
*     IP     THE ARRAY WHICH HAS PIVOTAL ROWS.                        *
*     IS     THE NUMBER OF DEGENERATED COLUMNS.                       *
*     IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.                 *
*   OTHERS:  WORKING PARAMETERS.                                      *
*                                                                     *
*  COPYRIFHT:    TSUTOMU OGUNI      SEP. 1 1992      VER. 2           *
***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(L),IA(L),X(M),NWK(0:M),B(M),E(L2),IE(L2),
     *           WE(M),NE(-1:2*M),IP(M)
C
      IF (L .LT. M) THEN
       WRITE(*,*) '(SUBR. ETALU) INVALID ARGUMENT. ', L, M
       IER = 2
       RETURN
      ENDIF
      EPS1 = EPS * 1.0D-2
      IER = 0
      IS = 0
      NE(-1) = 0
      NE(0) = 0
      NUM = 0
      DO 10 I=1,M
   10  IP(I) = 0
C
      DO 90 K=1,M
       DO 20 I=1,M
   20   WE(I) = 0.0D0
       DO 30 I=NWK(K-1)+1,NWK(K)
   30   WE(IA(I)) = A(I)
       IF (K .NE. 1) THEN
        JK = K - 1
        CALL UPDTL(E,IE,NE,WE,IP,L2,JK,M)
C  COMPUTATION OF ETA
       ENDIF
C  SELECTION OF PIVOT
       IROW = 0
       VAL = 0.0D0
       DO 60 I=K,M
        WORK = DABS(WE(I))
        IF( WORK .GT. VAL) THEN
         VAL = WORK
         IROW = I
        ENDIF
   60  CONTINUE
       IF (VAL .LE. EPS) THEN
        WRITE(*,*) '(SUBR. ETALU) STOP WITH ZERO PIVOT', K
        IER = 1
        RETURN
       ENDIF
       IP(K) = IROW
       IF (IROW .NE. K) THEN
        WORK = WE(K)
        WE(K) = WE(IROW)
        WE(IROW) = WORK
       ENDIF
C      WRITE(*,*) (WE(I),I=1,M)
C  GENERATION OF ETA
       PIV = - 1.0D0 / WE(K)
       WE(K) = - 1.0D0
       DO 70 I=K,M
        WE(I) = WE(I) * PIV
        IF (DABS(WE(I)) .GE. EPS1) THEN
         NUM = NUM + 1
         E(NUM) = WE(I)
         IE(NUM) = I
        ENDIF
   70  CONTINUE
       NE(2*K-1) = NUM
       IF (K .EQ. 1) THEN
        NE(2) = NE(1)
       ELSE
        DO 80 I=1,K-1
         IF (DABS(WE(I)) .GE. EPS1) THEN
          NUM = NUM + 1
          E(NUM) = - WE(I)
          IE(NUM) = I
         ENDIF
   80   CONTINUE
        NE(2*K) = NUM
       ENDIF
C
   90 CONTINUE
C  COMPUTATION OF X
      DO 100 K=1,M
  100  X(K) = B(K)
      CALL UPDTL(E,IE,NE,X,IP,L2,M,M)
      CALL UPDTU(E,IE,NE,X,L2,M)
C
      RETURN
      END
C
      SUBROUTINE UPDTL(E,IE,NE,WE,IP,L2,K,M)
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION E(L2),IE(L2),NE(-1:2*M),WE(M),IP(M)
C
      DO 20 J=1,K
       IROW = IP(J)
       IF (IROW .NE. J) THEN
        WORK = WE(J)
        WE(J) = WE(IROW)
        WE(IROW) = WORK
       ENDIF
       WORK = WE(J)
       WE(J) = 0.0D0
       DO 10 I=NE(2*J-2)+1,NE(2*J-1)
   10   WE(IE(I)) = WE(IE(I)) + E(I) * WORK
   20 CONTINUE
      RETURN
      END
C
      SUBROUTINE UPDTU(E,IE,NE,WE,L2,M)
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION E(L2),IE(L2),NE(-1:2*M),WE(M)
C
      DO 40 J=M,2,-1
       IF (NE(2*J-1) .NE. NE(2*J)) THEN
        WORK = WE(J)
        DO 50 I=NE(2*J-1)+1,NE(2*J)
   50    WE(IE(I)) = WE(IE(I)) + E(I) * WORK
       ENDIF
   40 CONTINUE
      RETURN
C  END OF ETALU
      END

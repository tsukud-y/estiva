*
      SUBROUTINE STWART(IA,L,MA,M,R,C,IR,IC,JROW,JCOL,IP,JP,
     *                  KERNS,MEND,IW,LG,IER)
***********************************************************************
*  STWART METHOD OF BLOCKING FOR NON-SYMMETRIC SPARSE MATRIX.         *
*                                                                     *
*  PARAMETERS:                                                        *
*   ON ENTRY:                                                         *
*     IA     THE ARRAY WHICH CONTAINS ROW INDEX OF NON-ZERO ELEMENTS  *
*            OF THE MATRIX.                                           *
*     L      THE LEADING DIMENSION OF THE ARRAY A.                    *
*     MA     ACCUMULATED SUM OF NUMBERS OF NON-ZERO ELEMENTS IN       *
*            EACH COLUMN OF THE MATRIX A.                             *
*     M      THE ORDER OF THE MATRIX A.                               *
*   ON RETURN:                                                        *
*     JROW   THE INFORMATION ABOUT CHANGE OF ROWS.                    *
*     JCOL   THE INFORMATION ABOUT CHANGE OF COLUMNS.                 *
*     IW     BLOCK INDEX OF EACH COLUMN.                              *
*     LG     NUMBER OF BLOCKS OF THE KERNEL.                          *
*     IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.                 *
*   OTHERS:  WORKING PARAMETERS.                                      *
*                                                                     *
*  COPYRIGHT:     TSUTOMU OGUNI     SEP. 1 1991        VER. 1         *
***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       INTEGER*4 R, C
       DIMENSION IA(L),MA(0:M),R(M),C(M),IR(M),IP(M),JP(M),
     *   IC(M),JROW(M),JCOL(M),IW(M)
C
      IER = 0
      DO 1 I=1,M
       R(I) = 0
       C(I) = 0
       IW(I) = 0
       JROW(I) = 0
       JCOL(I) = 0
       IR(I) = 0
    1  IC(I) = 0
C
      CALL SINGLE(IA,L,IR,IC,R,C,MA,M,MEND,KERNS,JCOL,JROW)
      CALL FORFUL(IA,L,MA,M,IP,JP,IR,IC,KERNS,IER)
C
      KERN = KERNS
      MM = M
      LG = 1
    5 CONTINUE
C
      DO 10 I=1,M
       R(I) = 0
   10  C(I) = 0
      DO 20 J=1,M
       IF (JCOL(J) .EQ. 0) THEN
        DO 30 N=MA(J-1)+1,MA(J)
         IF (JROW(IA(N)) .EQ. 0) THEN
          R(IA(N)) = R(IA(N)) + 1
          C(J) = C(J) + 1
         ENDIF
   30   CONTINUE
       ENDIF
   20 CONTINUE
C  SEARCH OF MINIMUM ROW COUNT
      MIN = 100000
      IROW = 0
      DO 40 I=1,M
       IF (JROW(I) .EQ. 0) THEN
        IF (R(I) .LT. MIN) THEN
         MIN = R(I)
         IROW = I
        ENDIF
       ENDIF
   40 CONTINUE
      IF (IROW .EQ. 0) THEN
       WRITE(*,*) '(SUBR. STWART) STOP AT. ', KERN
       IER = 1
       RETURN
      ENDIF
      DO 50 J=1,M
       IF (JCOL(J) .EQ. 0) THEN
        DO 55 N=MA(J-1)+1,MA(J)
         IF (IA(N) .EQ. IROW) THEN
          IW(J) = LG
          KERN = KERN + 1
         ENDIF
   55   CONTINUE
       ENDIF
   50 CONTINUE
C
   57 CONTINUE
      IS = 0
      DO 60 J=1,M
       IF (IW(J) .EQ. LG) THEN
        IROW = JP(J)
        DO 70 K=1,M
         IF (JCOL(K) .EQ. 0) THEN
          IF (IW(K) .EQ. 0) THEN
           DO 80 N=MA(K-1)+1,MA(K)
            IF (IA(N) .EQ. IROW) THEN
             IW(K) = LG
             KERN = KERN + 1
             IS = IS + 1
            ENDIF
   80      CONTINUE
          ENDIF
         ENDIF
   70   CONTINUE
       ENDIF
   60 CONTINUE
C
      IF (IS .NE. 0) GO TO 57
      DO 90 J=1,M
       IF (IW(J) .EQ. LG) THEN
        JCOL(J) = MM
        JROW(JP(J)) = MM
        MM = MM - 1
       ENDIF
   90 CONTINUE
      IF (LG .LT. M) THEN
       IF (KERN .LE. M) THEN
        LG = LG + 1
        GO TO 5
       ENDIF
      ENDIF
      IER = 0
      RETURN
C  END OF STWART
      END

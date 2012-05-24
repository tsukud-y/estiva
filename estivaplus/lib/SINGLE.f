*
      SUBROUTINE SINGLE(IA,L,IR,IC,R,C,MA,M,MEND,KERNS,JCOL,JROW)
***********************************************************************
*  DETECTION OF SINGLETONS OF A SPARSE MATRIX.                        *
*                                                                     *
*  PARAMETERS:                                                        *
*   ON ENTRY:                                                         *
*     SAME AS MARKWZ ROUTINE.                                         *
*   ON RETURN:                                                        *
*     IR     ROW INDEX OF ROW SINGLETONS IN THE FIRST MEND ENTRIES    *
*            IN THE ARRAY IR. ROW INDEX OF COLUMN SINGLETONS IN       *
*            (KERNS-MEND-1) ENTRIES IN THE ARRAY IR.                  *
*     IC     COLUMN INDEX OF ROW SINGLETONS IN THE FIRST MEND ENTRIES *
*            IN THE ARRAY IC. COLUMN INDEX OF COLUMN SINGLETON IN     *
*            (KERNS-MEND-1) ENTRIES IN THE ARRAY IC.                  *
*   OTHERS:  THE WORKING PARAMETERS.                                  *
*                                                                     *
*  COPYRIGHT:      TSUTOMU OGUNI   SEP. 1 1992       VER. 2           *
***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       INTEGER*4 R, C
       DIMENSION IA(L),IC(M),IR(M),R(M),C(M),MA(0:M),JCOL(M),JROW(M)
C  SINGLETON
C  ROW SINGLETON
      DO 100 I=1,M
       JROW(I) = 0
       JCOL(I) = 0
       IC(I) = 0
       IR(I) = 0
       C(I) = MA(I) - MA(I-1)
       DO 110 N=MA(I-1)+1,MA(I)
  110   R(IA(N)) = R(IA(N)) + 1
  100 CONTINUE
      KERNS = 1
      DO 150 K=1,M
       IROW = 0
       ICOL = 0 
       DO 160 I=1,M
        IF (R(I) .EQ. 1) THEN
         IROW = I
         GO TO 170
        ENDIF
  160   CONTINUE
  170  CONTINUE    
       IF (IROW .EQ. 0) GO TO 180
       IR(KERNS) = IROW
       JROW(IROW) = KERNS
C
       DO 200 J=1,M
        IF (JCOL(J) .EQ. 0) THEN
         DO 190 N=MA(J-1)+1,MA(J)
          IF (IA(N) .EQ. IROW) THEN
           IC(KERNS) = J
           ICOL = J
           GO TO 200
          ENDIF
  190    CONTINUE
        ENDIF
  200  CONTINUE 
       JCOL(ICOL) = KERNS
       KERNS = KERNS + 1
       DO 210 N=MA(ICOL-1)+1,MA(ICOL)
        R(IA(N)) = R(IA(N)) - 1
  210  CONTINUE
  150 CONTINUE
C
  180 CONTINUE 
      MEND = KERNS - 1
C  COLUMN SINGLETON
      DO 120 K=1,M
       IROW = 0
       ICOL = 0
       DO 130 J=1,M
        IF (JCOL(J) .EQ. 0) THEN
         IF (C(J) .EQ. 1) THEN
          ICOL = J
          GO TO 135
         ENDIF
        ENDIF
  130  CONTINUE
  135 CONTINUE
       IF (ICOL .EQ. 0) RETURN
        JCOL(ICOL) = KERNS
        IC(KERNS) = ICOL
C
        DO 195 N=MA(ICOL-1)+1,MA(ICOL)
         IF (JROW(IA(N)) .EQ. 0) IROW = IA(N)
  195  CONTINUE
       IR(KERNS) = IROW
       JROW(IROW) = KERNS     
       KERNS = KERNS + 1
       DO 175 J=1,M
        IF (JCOL(J) .EQ. 0) THEN
         DO 185 N=MA(J-1)+1,MA(J)
          IF (IA(N) .EQ. IROW) C(J) = C(J) - 1
  185    CONTINUE
        ENDIF
  175  CONTINUE
  120 CONTINUE
C
      RETURN
C  END OF SINGLE
      END

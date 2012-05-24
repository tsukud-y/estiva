*
      SUBROUTINE FORFUL(IA,L,MA,M,IP,JP,IR,IC,KERNS,IER)
***********************************************************************
*  FORD-FULKERSON METHOD FOR ORDERING.                                *
*                                                                     *
*  PARAMETERS:                                                        *
*   ON ENTRY:  SAME AS STWART ROUTINE.                                *
*   ON RETURN:                                                        *
*     IP     COLUMN INDEX IN EACH I-TH ROW. ON THE CASE OF SINGLETON, *
*            MINUS SIGN.                                              *
*     JP     ROW INDEX IN EACH J-TH COLUMN.                           *
*     IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.                 *
*   OTHERS:  WORKING PARAMETERS.                                      *
*                                                                     *
*  COPYRIGHT:    TSUTOMU OGUNI      SEP. 1 1991         VER. 1        *
***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION IA(L), MA(0:M), IP(M), JP(M), IR(M), IC(M)
C
      IER = 0
      DO 1 K=1,M
       JP(K) = 0
    1  IP(K) = 0
      DO 3 K=1,KERNS-1
       IP(IR(K)) = - IC(K)
    3  JP(IC(K)) = IR(K)
      I = 1
      J = 0
    2 CONTINUE
      IF (JP(I) .EQ. 0) THEN
       DO 10 K=MA(I-1)+1,MA(I)
        IF (IP(IA(K)) .EQ. 0) THEN
         NK = IA(K)
         IP(NK) = I
         JP(I) = NK
         GO TO 20
        ENDIF
   10  CONTINUE
       GO TO 50
      ENDIF
   20 IF (I .GE. M) RETURN
      I = I + 1
      GO TO 2
   50 CONTINUE
      NN = 0
   52 NN = NN + 1
      IF (NN .GT. (MA(I) - MA(I-1))) THEN
       WRITE(*,*) '(SUBR. FORFUL) ERROR STOP. ', I, J
       IER = 1
       RETURN
      ENDIF 
      NH = IA(MA(I-1)+NN)
      J = IP(NH)
      IF (J .LE. 0) GO TO 52
      IP(NH) = I
      JP(I) = NH
   55 CONTINUE
      DO 60 K=MA(J-1)+1,MA(J)
       IF (IP(IA(K)) .EQ. 0) THEN
        NL = IA(K)
        IP(NL) = J
        JP(J) = NL
        GO TO 20
       ENDIF
   60 CONTINUE
      NN = 0
   62 NN = NN + 1
      IF (NN .GT. (MA(J) - MA(J-1))) THEN
       WRITE(*,*) '(SUBR. FORFUL) ERROR STOP. ', I, J
       IER = 1
       RETURN
      ENDIF
      NP = IA(MA(J-1)+NN)
      ND = IP(NP)
      IF (ND .LE. 0) GO TO 62
      IP(NP) = J
      JP(J) = NP
      J = ND
      GO TO 55
C  END OF FORFUL
      END

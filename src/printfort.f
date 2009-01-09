      SUBROUTINE PRINTFORT
C     ********************
C     -------------------------------------------------------------
     & ( CHAINE , TAILLE )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  : 
C  ---------
C      IMPRESSION D'UNE CHAINE DE CARACTERES (ISSUE D'UNE FONCTION C)
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "couplings_f.h"
C
C***********************************************************************
C
      CHARACTER     CHAINE(*)
      INTEGER       TAILLE
      CHARACTER     CHLOC*16384
      CHARACTER     NOM*64
      INTEGER       II

      TAILLE = MIN(TAILLE, 16384 - 1)
C
      DO II = 1, TAILLE
         CHLOC(II:II) = CHAINE(II)
      ENDDO
C
      WRITE(IFILE, 1000) CHLOC(1:TAILLE)
C
      RETURN

 1000 FORMAT(A, $)
C
      END
c@z

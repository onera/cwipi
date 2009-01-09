C
C***********************************************************************
C
C COUPLINGS_INIT_F 
C
C***********************************************************************
C
      
      SUBROUTINE COUPLINGS_INIT_F 
     &(GLOBALCOMM, OUTPUTUNIT, APPLINAME, APPLICOMM)
      
      IMPLICIT NONE

      include "couplings_f.h"

      INTEGER GLOBALCOMM
      INTEGER OUTPUTUNIT
      CHARACTER*(*) APPLINAME
      INTEGER APPLICOMM
      INTEGER L1
      INTEGER I

      L1 = LEN(APPLINAME)
      IFILE =  OUTPUTUNIT

      CALL COUPLINGS_INIT_CF (GLOBALCOMM, OUTPUTUNIT, 
     &                       APPLINAME, L1, 
     &                       APPLICOMM)

      END

C
C***********************************************************************
C
C COUPLINGS_ADD_LOCAL_INT_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_ADD_LOCAL_INT_CONTROL_PARAMETER_F
     &(NAME, INITIALVALUE)

      IMPLICIT NONE

      CHARACTER*(*) NAME
      INTEGER L
      INTEGER INITIALVALUE

      L = LEN(NAME)

      CALL COUPLINGS_ADD_LOCAL_INT_CONTROL_PARAMETER_CF
     &(NAME, L, INITIALVALUE) 

      END

C
C***********************************************************************
C
C COUPLINGS_ADD_LOCAL_DOUBLE_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_ADD_LOCAL_DOUBLE_CONTROL_PARAMETER_F
     &(NAME, INITIALVALUE)

      IMPLICIT NONE

      CHARACTER*(*) NAME
      INTEGER L
      DOUBLE PRECISION INITIALVALUE

      L = LEN(NAME)

      CALL COUPLINGS_ADD_LOCAL_DOUBLE_CONTROL_PARAMETER_CF
     &(NAME, L, INITIALVALUE) 

      END

C
C***********************************************************************
C
C COUPLINGS_SET_LOCAL_INT_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_SET_LOCAL_INT_CONTROL_PARAMETER_F
     &(NAME, INITIALVALUE)

      IMPLICIT NONE

      CHARACTER*(*) NAME
      INTEGER L
      INTEGER INITIALVALUE

      L = LEN(NAME)

      CALL COUPLINGS_SET_LOCAL_INT_CONTROL_PARAMETER_CF
     &(NAME, L, INITIALVALUE) 

      END

C
C***********************************************************************
C
C COUPLINGS_SET_LOCAL_DOUBLE_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_SET_LOCAL_DOUBLE_CONTROL_PARAMETER_F
     &(NAME, INITIALVALUE)

      IMPLICIT NONE

      CHARACTER*(*) NAME
      INTEGER L
      DOUBLE PRECISION INITIALVALUE

      L = LEN(NAME)

      CALL COUPLINGS_SET_LOCAL_DOUBLE_CONTROL_PARAMETER_CF
     &(NAME, L, INITIALVALUE) 

      END
C
C***********************************************************************
C
C COUPLINGS_GET_LOCAL_INT_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_GET_LOCAL_INT_CONTROL_PARAMETER_F
     &(NAME, VALUE)

      IMPLICIT NONE

      CHARACTER*(*) NAME
      INTEGER L
      INTEGER VALUE

      L = LEN(NAME)

      CALL COUPLINGS_GET_LOCAL_INT_CONTROL_PARAMETER_CF
     &(NAME, L, VALUE) 

      END

C
C***********************************************************************
C
C COUPLINGS_GET_LOCAL_DOUBLE_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_GET_LOCAL_DOUBLE_CONTROL_PARAMETER_F
     &(NAME, VALUE)

      IMPLICIT NONE

      CHARACTER*(*) NAME
      INTEGER L
      DOUBLE PRECISION VALUE

      L = LEN(NAME)

      CALL COUPLINGS_GET_LOCAL_DOUBLE_CONTROL_PARAMETER_CF
     &(NAME, L, VALUE) 

      END

C
C***********************************************************************
C
C COUPLINGS_DELETE_LOCAL_INT_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_DELETE_LOCAL_INT_CONTROL_PARAMETER_F
     &(NAME)

      IMPLICIT NONE

      CHARACTER*(*) NAME
      INTEGER L

      L = LEN(NAME)

      CALL COUPLINGS_DELETE_LOCAL_INT_CONTROL_PARAMETER_CF
     &(NAME, L) 

      END

C
C***********************************************************************
C
C COUPLINGS_DELETE_LOCAL_DOUBLE_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_DELETE_LOCAL_DOUBLE_CONTROL_PARAMETER_F
     &(NAME)

      IMPLICIT NONE

      CHARACTER*(*) NAME
      INTEGER L

      L = LEN(NAME)

      CALL COUPLINGS_DELETE_LOCAL_DOUBLE_CONTROL_PARAMETER_CF
     &(NAME, L) 

      END

C
C***********************************************************************
C
C COUPLINGS_GET_DISTANT_INT_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_GET_DISTANT_INT_CONTROL_PARAMETER_F
     &(APPLINAME, PARAMNAME, VALUE)

      IMPLICIT NONE

      CHARACTER*(*) APPLINAME
      CHARACTER*(*) PARAMNAME
      INTEGER L1, L2

      INTEGER VALUE

      L1 = LEN(APPLINAME)
      L2 = LEN(PARAMNAME)

      CALL COUPLINGS_GET_DISTANT_INT_CONTROL_PARAMETER_CF
     &(APPLINAME, L1, PARAMNAME, L2, VALUE) 

      END

C
C***********************************************************************
C
C COUPLINGS_GET_DISTANT_DOUBLE_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_GET_DISTANT_DOUBLE_CONTROL_PARAMETER_F
     &(APPLINAME, PARAMNAME, VALUE)

      IMPLICIT NONE

      CHARACTER*(*) APPLINAME
      CHARACTER*(*) PARAMNAME
      INTEGER L1, L2

      DOUBLE PRECISION VALUE

      L1 = LEN(APPLINAME)
      L2 = LEN(PARAMNAME)

      CALL COUPLINGS_GET_DISTANT_DOUBLE_CONTROL_PARAMETER_CF
     &(APPLINAME, L1, PARAMNAME, L2, VALUE) 

      END

C
C***********************************************************************
C
C COUPLINGS_SYNCHRONISE_CONTROL_PARAMETER_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_SYNCHRONISE_CONTROL_PARAMETER_F
     &(APPLINAME)

      IMPLICIT NONE

      CHARACTER*(*) APPLINAME
      
      INTEGER L

      L = LEN(APPLINAME)

      CALL COUPLINGS_SYNCHRONISE_CONTROL_PARAMETER_CF
     &(APPLINAME, L)

      END

 

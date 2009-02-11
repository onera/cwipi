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
 
C
C***********************************************************************
C
C COUPLINGS_CREATE_COUPLING_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_CREATE_COUPLING_F
     &(COUPLINGNAME,
     & CPLAPPLI,
     & DIM,
     & TOLERANCE,
     & MESHT,
     & SOLVERT,
     & OUTPUTFREQ,
     & OUTPUTFMT,
     & OUTPUTFMTOPT)

      IMPLICIT NONE

      CHARACTER*(*) COUPLINGNAME, CPLAPPLI
      INTEGER DIM, MESHT, SOLVERT
      DOUBLE PRECISION TOLERANCE, OUTPUTFREQ
      CHARACTER*(*) OUTPUTFMT, OUTPUTFMTOPT

      INTEGER LCOUPLINGNAME, LCPLAPPLI
      INTEGER LOUTPUTFMT, LOUTPUTFMTOPT

      LCOUPLINGNAME = LEN(COUPLINGNAME)
      LCPLAPPLI     = LEN(CPLAPPLI)
      LOUTPUTFMT    = LEN(OUTPUTFMT)
      LOUTPUTFMTOPT = LEN(OUTPUTFMTOPT)

      CALL COUPLINGS_CREATE_COUPLING_CF(COUPLINGNAME,
     &                                  LCOUPLINGNAME,
     &                                  CPLAPPLI,
     &                                  LCPLAPPLI,
     &                                  DIM,
     &                                  TOLERANCE,
     &                                  MESHT,
     &                                  SOLVERT,
     &                                  OUTPUTFREQ,
     &                                  OUTPUTFMT,
     &                                  LOUTPUTFMT,
     &                                  OUTPUTFMTOPT,
     &                                  LOUTPUTFMTOPT)

      END
C
C***********************************************************************
C
C COUPLINGS_SET_POINTS_TO_LOCATE_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_SET_POINTS_TO_LOCATE_F
     &(COUPLINGNAME,
     & NPTS,
     & COORDS)

      IMPLICIT NONE

      CHARACTER*(*) COUPLINGNAME
      INTEGER LCOUPLINGNAME
      
      INTEGER NPTS
      DOUBLE PRECISION COORDS(3*NPTS)
      
      LCOUPLINGNAME    = LEN(COUPLINGNAME)
      
      CALL COUPLINGS_SET_POINTS_TO_LOCATE_CF(COUPLINGNAME, 
     &                                       LCOUPLINGNAME, 
     &                                       NPTS,
     &                                       COORDS)

      END
C
C***********************************************************************
C
C COUPLINGS_DEFINE_MESH_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_DEFINE_MESH_F
     &(COUPLINGNAME,
     & NVERTEX,
     & NELTS,
     & COORDS,
     & CONNECINDEX,
     & CONNEC)

      IMPLICIT NONE

      CHARACTER*(*) COUPLINGNAME
      INTEGER LCOUPLINGNAME

      INTEGER NELTS, NVERTEX, CONNECINDEX(NELTS+1), CONNEC(*)
      DOUBLE PRECISION COORDS(3*NVERTEX)

      LCOUPLINGNAME    = LEN(COUPLINGNAME)

      CALL COUPLINGS_DEFINE_MESH_CF(COUPLINGNAME,
     &                              LCOUPLINGNAME,
     &                              NVERTEX,
     &                              NELTS,
     &                              COORDS,
     &                              CONNECINDEX,
     &                              CONNEC) 

      END

C
C***********************************************************************
C
C COUPLINGS_ADD_POLYHEDRA_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_ADD_POLYHEDRA_F
     &(COUPLINGNAME,
     & NELTS,
     & FACEIDX,
     & CELLTOFACE,
     & FACECONNECIDX,
     & FACECONNEC)

      IMPLICIT NONE

      CHARACTER*(*) COUPLINGNAME
      INTEGER LCOUPLINGNAME, NELTS
      INTEGER FACEIDX(NELTS), CELLTOFACE(*)
      INTEGER FACECONNECIDX(*), FACECONNEC(*)

      LCOUPLINGNAME = LEN(COUPLINGNAME)

      CALL COUPLINGS_ADD_POLYHEDRA_CF (COUPLINGNAME,
     &                                 LCOUPLINGNAME,
     &                                 NELTS,
     &                                 FACEIDX,
     &                                 CELLTOFACE,
     &                                 FACECONNECIDX,
     &                                 FACECONNEC)

      END

C
C***********************************************************************
C
C COUPLINGS_EXCHANGE_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_EXCHANGE_F
     &(COUPLINGNAME,
     & EXCHANGENAME,
     & EXCHANGEDIM, 
     & NSTEP, 
     & TIMEVALUE,
     & SENDINGFIELDNAME,
     & SENDINGFIELD,
     & RECEIVINGFIELDNAME,
     & RECEIVINGFIELD,
     & STATUS)

      IMPLICIT NONE

      CHARACTER*(*) COUPLINGNAME, EXCHANGENAME, SENDINGFIELDNAME
      CHARACTER*(*) RECEIVINGFIELDNAME
      INTEGER LCOUPLINGNAME, LEXCHANGENAME, LSENDINGFIELDNAME
      INTEGER LRECEIVINGFIELDNAME
      INTEGER EXCHANGEDIM, NSTEP, STATUS
      DOUBLE PRECISION TIMEVALUE
      DOUBLE PRECISION SENDINGFIELD(*), RECEIVINGFIELD(*)

      LCOUPLINGNAME       = LEN(COUPLINGNAME)
      LEXCHANGENAME       = LEN(EXCHANGENAME)
      LSENDINGFIELDNAME   = LEN(SENDINGFIELDNAME)
      LRECEIVINGFIELDNAME = LEN(RECEIVINGFIELDNAME)
      
      CALL COUPLINGS_EXCHANGE_CF(COUPLINGNAME,
     &                           LCOUPLINGNAME,
     &                           EXCHANGENAME,
     &                           LEXCHANGENAME,
     &                           EXCHANGEDIM, 
     &                           NSTEP, 
     &                           TIMEVALUE,
     &                           SENDINGFIELDNAME,
     &                           LSENDINGFIELDNAME,
     &                           SENDINGFIELD,
     &                           RECEIVINGFIELDNAME,
     &                           LRECEIVINGFIELDNAME,
     &                           RECEIVINGFIELD,
     &                           STATUS)

      END

C
C***********************************************************************
C
C COUPLINGS_EXCHANGE_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_EXCHANGE_F
     &(COUPLINGNAME,
     & EXCHANGENAME,
     & EXCHANGEDIM, 
     & NSTEP, 
     & TIMEVALUE,
     & SENDINGFIELDNAME,
     & SENDINGFIELD,
     & RECEIVINGFIELDNAME,
     & RECEIVINGFIELD,
     & PTINTERPOLATIONFCT,
     & STATUS)

      IMPLICIT NONE

      CHARACTER*(*) COUPLINGNAME, EXCHANGENAME, SENDINGFIELDNAME
      CHARACTER*(*) RECEIVINGFIELDNAME
      INTEGER LCOUPLINGNAME, LEXCHANGENAME, LSENDINGFIELDNAME
      INTEGER LRECEIVINGFIELDNAME
      INTEGER EXCHANGEDIM, NSTEP, STATUS
      DOUBLE PRECISION TIMEVALUE
      DOUBLE PRECISION SENDINGFIELD(*), RECEIVINGFIELD(*)

      LCOUPLINGNAME       = LEN(COUPLINGNAME)
      LEXCHANGENAME       = LEN(EXCHANGENAME)
      LSENDINGFIELDNAME   = LEN(SENDINGFIELDNAME)
      LRECEIVINGFIELDNAME = LEN(RECEIVINGFIELDNAME)
      
      CALL COUPLINGS_EXCHANGE_CF(COUPLINGNAME,
     &                           LCOUPLINGNAME,
     &                           EXCHANGENAME,
     &                           LEXCHANGENAME,
     &                           EXCHANGEDIM, 
     &                           NSTEP, 
     &                           TIMEVALUE,
     &                           SENDINGFIELDNAME,
     &                           LSENDINGFIELDNAME,
     &                           SENDINGFIELD,
     &                           RECEIVINGFIELDNAME,
     &                           LRECEIVINGFIELDNAME,
     &                           RECEIVINGFIELD,
     &                           STATUS)

      END

C
C***********************************************************************
C
C COUPLINGS_DELETE_COUPLING_F
C
C***********************************************************************
C

      SUBROUTINE COUPLINGS_DELETE_COUPLING_F(COUPLINGNAME)

      IMPLICIT NONE

      CHARACTER*(*) COUPLINGNAME 
      INTEGER LCOUPLINGNAME

      LCOUPLINGNAME       = LEN(COUPLINGNAME)

      CALL COUPLINGS_DELETE_COUPLING_CF(COUPLINGNAME, LCOUPLINGNAME)

      END

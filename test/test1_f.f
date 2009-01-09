      PROGRAM TESTF
      
      IMPLICIT NONE

      include "mpif.h"
      include "couplings_f.h"

      INTEGER LOCALCOM
      INTEGER IRANK
      CHARACTER *4 PROC
      INTEGER CODE
      INTEGER IIUNIT
      INTEGER IVALUE
      DOUBLE PRECISION DVALUE

      CALL MPI_INIT(CODE)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, CODE)
C
C CREATION DU FICHIER DE SORTIE LISTING
C (En parallele, un fichier listing par processus)
      WRITE(PROC,'(i4.4)') IRANK
      IIUNIT = 9
      OPEN (UNIT=IIUNIT, FILE='listing_code_fortran_'//PROC,
     &      FORM='FORMATTED', STATUS='UNKNOWN')
C
C INITIALISATION DE L'INTERFACE DE COUPLAGE
      CALL COUPLINGS_INIT_F (MPI_COMM_WORLD, 
     &                       IIUNIT, 
     &                       "CodeFortran", 
     &                       LOCALCOM)

      WRITE(IIUNIT,*)
      WRITE(IIUNIT,*) "Dump apres initialisation"
      WRITE(IIUNIT,*) "-------------------------"
      WRITE(IIUNIT,*)
      CALL COUPLINGS_DUMP_APPLICATION_PROPERTIES_F
C
C AJOUT DE PARAMETRES DE CONTROLE
      CALL COUPLINGS_ADD_LOCAL_INT_CONTROL_PARAMETER_F(
     &     "niterF", 10)

      CALL COUPLINGS_ADD_LOCAL_DOUBLE_CONTROL_PARAMETER_F(
     &     "physicalTimeF", 1.123D0)

      WRITE(IIUNIT,*)
      WRITE(IIUNIT,*) "Dump apres ajout de parametres"
      WRITE(IIUNIT,*) "------------------------------"
      WRITE(IIUNIT,*)
      CALL COUPLINGS_DUMP_APPLICATION_PROPERTIES_F
C
C MODIFICATION DES PARAMETRES DE CONTROLE
      CALL COUPLINGS_GET_LOCAL_INT_CONTROL_PARAMETER_F(
     &     "niterF", IVALUE)

      IVALUE = IVALUE + 1

      CALL COUPLINGS_SET_LOCAL_INT_CONTROL_PARAMETER_F(
     &     "niterF", IVALUE)

      CALL COUPLINGS_GET_LOCAL_DOUBLE_CONTROL_PARAMETER_F(
     &     "physicalTimeF", DVALUE)

      DVALUE = DVALUE + 0.1D0

      CALL COUPLINGS_SET_LOCAL_DOUBLE_CONTROL_PARAMETER_F(
     &     "physicalTimeF", DVALUE)

      WRITE(IIUNIT,*)
      WRITE(IIUNIT,*) "Dump apres modification de parametres"
      WRITE(IIUNIT,*) "-------------------------------------"
      WRITE(IIUNIT,*)
      CALL COUPLINGS_DUMP_APPLICATION_PROPERTIES_F
C
C ECHANGE DES PARAMETRES DE CONTROLE
      CALL COUPLINGS_SYNCHRONISE_CONTROL_PARAMETER_F("CodeC")
      WRITE(IIUNIT,*)
      WRITE(IIUNIT,*) "Dump apres synchronisation"
      WRITE(IIUNIT,*) "--------------------------"
      CALL COUPLINGS_DUMP_APPLICATION_PROPERTIES_F
C
C SUPPRESSION DES PARAMETRES DE CONTROLE
      CALL COUPLINGS_DELETE_LOCAL_INT_CONTROL_PARAMETER_F("niterF")
      CALL COUPLINGS_DELETE_LOCAL_DOUBLE_CONTROL_PARAMETER_F(
     &     "physicalTimeF")

      WRITE(IIUNIT,*)
      WRITE(IIUNIT,*) "Dump apres suppression des parametres"
      WRITE(IIUNIT,*) "-------------------------------------"
      WRITE(IIUNIT,*)
      CALL COUPLINGS_DUMP_APPLICATION_PROPERTIES_F

      CALL MPI_BARRIER(MPI_COMM_WORLD, CODE)
      CALL MPI_FINALIZE(CODE)
      END

!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2022  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

module pdm_io

  use pdm

  implicit none

  !
  ! Parametres
  ! ----------

  !
  ! Types de sufixe

  integer (kind = pdm_l_num_s), parameter :: pdm_io_suff_auto        = 0 ! Suffixe detemine automatiquement
  integer (kind = pdm_l_num_s), parameter :: pdm_io_suff_man         = 1 ! Suffixe fourni par l'utilisateur

  !
  ! Endianness

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_BIGENDIAN        = 0 ! Contenu en bigendian
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_LITTLEENDIAN     = 1 ! Contenu en little endian
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_NATIVE           = 3 ! Contenu natif machine

  !
  ! Types de données

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_INT    = 0 ! Type de donnee int
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_LONG   = 1 ! Type de donnee long
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_DOUBLE = 2 ! Type de donnee double
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_FLOAT  = 3 ! Type de donnee float
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_CHAR   = 4 ! Type de donnee char

  !
  ! Types d'entrees/sorties paralleles

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_ACCES_MPIIO_EO   = 0 ! Acces parallele avec MPIIO (explicit offset)
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_ACCES_MPIIO_IP   = 1 ! Acces parallele avec MPIIO (individual pointer)
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_ACCES_MPI_SIMPLE = 2 ! Acces parallele sans MPI
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_ACCES_SEQ        = 3 ! Acces 1 fichier par processus

  !
  ! Mode d'acces lecture, ecriture, lecture/ecriture

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_MODE_LECTURE  = 0    ! Acces en lecture
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_MODE_ECRITURE = 1    ! Acces en ecriture
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_MODE_AJOUT    = 2    ! Acces en lecture/ecriture

  !
  ! Indique si le nombre de sous-variables d'une variable est constant ou variable
  ! selon le numero de la donnee

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_N_COMPOSANTE_CONSTANT = 0
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_N_COMPOSANTE_VARIABLE = 1

  !
  ! Indique si les donnees lues ou écrites sont rangees par blocs contigus en memoire
  ! ou respectent une indirection (donnees entrelacees)

  integer(kind = pdm_l_num_s), parameter :: PDM_IO_DONNEES_BLOC       = 0
  integer(kind = pdm_l_num_s), parameter :: PDM_IO_DONNEES_ENTRELACEE = 1

  !
  ! Indique si le fichier contient une entete IOCEDRE

  integer(kind = pdm_l_num_s), parameter :: PDM_IO_ENTETE_ON  = 0
  integer(kind = pdm_l_num_s), parameter :: PDM_IO_ENTETE_OFF = 1

  !
  ! Indique le format du fichier

  integer(kind = pdm_l_num_s), parameter :: PDM_IO_FMT_TXT = 0
  integer(kind = pdm_l_num_s), parameter :: PDM_IO_FMT_BIN = 1

  !
  ! Active ou non le backup d'un fichier preexistant

  integer(kind = pdm_l_num_s), parameter :: PDM_IO_BACKUP_ON  = 0
  integer(kind = pdm_l_num_s), parameter :: PDM_IO_BACKUP_OFF = 1

  interface

  !>
  !! \brief Set the file position indicator
  !!
  !! \param [in] fichier         Pointer to \ref PDM_io_fichier_t object
  !! \param [in] offset          Adress
  !! \param [in] seek            Origin type
  !!
  !!

  subroutine PDM_io_seek (fichier, &
                          offset,  &
                          seek)    &
  bind (c, name='PDM_io_seek')
    use iso_c_binding
    implicit none

    type(c_ptr),     value :: fichier
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: offset
#else
    integer(c_int),  value :: offset
#endif
    integer(c_int),  value :: seek

  end subroutine PDM_io_seek


  !>
  !! \brief Lecture globale : Le processus maitre accede seul au fichier et redistribue
  !! l'information a l'ensemble des processus du communicateur
  !!
  !! \param [in]  fichier         Pointer to \ref PDM_io_fichier_t object
  !! \param [in]  taille_donnee   Taille unitaire de la donnee
  !! \param [in]  n_donnees       Nombre de donnees a lire
  !! \param [out] donnees         Donnees lues
  !!
  !!

  subroutine PDM_io_lecture_globale (fichier,       &
                                     taille_donnee, &
                                     n_donnees,     &
                                     donnees)       &
  bind (c, name='PDM_io_lecture_globale')
    use iso_c_binding
    implicit none

    type(c_ptr), value          :: fichier
    integer(c_int),  intent(in) :: taille_donnee
#ifdef PDM_LONG_G_NUM
    integer(c_long), intent(in) :: n_donnees
#else
    integer(c_int),  intent(in) :: n_donnees
#endif
    type(c_ptr)                 :: donnees

  end subroutine PDM_io_lecture_globale


  !>
  !! \brief Ecriture globale : Le processus maitre accede seul au fichier
  !!
  !! \param [in]  fichier         Pointer to \ref PDM_io_fichier_t object
  !! \param [in]  taille_donnee   Taille unitaire de la donnee
  !! \param [in]  n_donnees       Nombre de donnees a ecrire
  !! \param [in]  donnees         Donnees ecrites
  !!
  !!

  subroutine PDM_io_ecriture_globale (fichier,       &
                                      taille_donnee, &
                                      n_donnees,     &
                                      donnees)       &
  bind (c, name='PDM_io_ecriture_globale')
    use iso_c_binding
    implicit none

    type(c_ptr), value          :: fichier
    integer(c_int),  intent(in) :: taille_donnee
#ifdef PDM_LONG_G_NUM
    integer(c_long), intent(in) :: n_donnees
#else
    integer(c_int),  intent(in) :: n_donnees
#endif
    type(c_ptr)                 :: donnees

  end subroutine PDM_io_ecriture_globale


  !>
  !! \brief Fermeture du fichier sans destruction de la structure PDM_io associee a
  !! l'unite
  !!
  !! \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
  !!
  !!

  subroutine PDM_io_close (fichier) &
  bind (c, name='PDM_io_close')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: fichier

  end subroutine PDM_io_close


  !>
  !! \brief Destruction de la structure PDM_io associee a l'unite
  !!
  !! \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
  !!
  !!

  subroutine PDM_io_detruit (fichier) &
  bind (c, name='PDM_io_detruit')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: fichier

  end subroutine PDM_io_detruit


  !>
  !! \brief Affiche les informations sur le fichier
  !!
  !! \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
  !!
  !!

  subroutine PDM_io_dump (fichier) &
  bind (c, name='PDM_io_dump')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: fichier

  end subroutine PDM_io_dump


  !>
  !! \brief Retourne le communicateur du fichier
  !!
  !! \param [in]  fichier     Pointer to \ref PDM_io_fichier_t object
  !! \param [out] f_comm      Communicateur MPI
  !!
  !!

  subroutine PDM_io_get_comm (fichier, &
                              f_comm)  &
  bind (c, name='PDM_io_get_comm')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: fichier
    integer(c_int)     :: f_comm

  end subroutine PDM_io_get_comm


  !>
  !! \brief Active le swap endian
  !!
  !! \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
  !!
  !!

  subroutine PDM_io_swap_endian_on (fichier) &
  bind (c, name='PDM_io_swap_endian_on')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: fichier

  end subroutine PDM_io_swap_endian_on


  !>
  !! \brief Désactive le swap endian
  !!
  !! \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
  !!
  !!

  subroutine PDM_io_swap_endian_off (fichier) &
  bind (c, name='PDM_io_swap_endian_off')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: fichier

  end subroutine PDM_io_swap_endian_off


  !>
  !! \brief Initialise une phase d'écriture parallèle de tableaux de données associées
  !! aux numéros de variables PDM
  !! Chaque tableau a ses propres caractéristiques :
  !!         - taille de données
  !!         - nombre de donnée
  !!         - indirection (numérotation absolue)
  !!
  !! \param [in] unite              Unite du fichier
  !! \param [in] t_rangement        Type de rangement
  !! \param [in] num_var_cedre_max  Numéro max de variable PDM
  !! \param [in] n_partition_local  Nombre de partitions locales
  !!
  !!

  subroutine PDM_io_tab_ecr_debut (unite,             &
                                   t_rangement,       &
                                   num_var_cedre_max, &
                                   n_partition_local) &
  bind (c, name='PDM_io_tab_ecr_debut')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: unite
    integer(c_int), value :: t_rangement
    integer(c_int), value :: num_var_cedre_max
    integer(c_int), value :: n_partition_local

  end subroutine PDM_io_tab_ecr_debut


  !>
  !! \brief Definition d'une variable en ecriture
  !!
  !! \param [in] num_var_cedre          Numéro de variable PDM
  !! \param [in] num_indirection_cedre  Numéro d'indirection PDM
  !! \param [in] t_n_composantes        Type de tailles composantes (PDM_
  !! \param [in] n_composantes          Nombre de composantes pour chaque
  !! \param [in] taille_donnee          Taille unitaire de la donnnee
  !!

  subroutine PDM_io_tab_ecr_def_var (num_var_cedre,         &
                                     num_indirection_cedre, &
                                     t_n_composantes,       &
                                     n_composantes,         &
                                     taille_donnee)         &
  bind (c, name='PDM_io_tab_ecr_def_var')
    use iso_c_binding
    implicit none

    integer(c_int), value :: num_var_cedre
    integer(c_int), value :: num_indirection_cedre
    integer(c_int), value :: t_n_composantes
    integer(c_int), value :: n_composantes
    integer(c_int), value :: taille_donnee

  end subroutine PDM_io_tab_ecr_def_var


  !>
  !! \brief Finalise une phase d'écriture parallèle de tableaux de données associées
  !! aux numéros de variables PDM. Cette fonction déclenche réellement
  !! les écritures
  !!

  subroutine PDM_io_tab_ecr_fin () &
  bind (c, name='PDM_io_tab_ecr_fin')
    use iso_c_binding
    implicit none

  end subroutine PDM_io_tab_ecr_fin


  !>
  !! \brief Initialise une phase de lecture parallèle de tableaux de données associées
  !! aux numéros de variables PDM
  !! Chaque tableau a ses propres caractéristiques :
  !!         - taille de données
  !!         - nombre de donnée
  !!         - indirection (numérotation absolue)
  !!
  !! \param [in] unite              Unite du fichier
  !! \param [in] t_rangement        Type de rangement
  !! \param [in] num_var_cedre_max  Numéro max de variable PDM
  !! \param [in] n_partition_local  Nombre de partitions locales
  !!

  subroutine PDM_io_tab_lec_debut (unite,             &
                                   t_rangement,       &
                                   num_var_cedre_max, &
                                   n_partition_local) &
  bind (c, name='PDM_io_tab_lec_debut')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: unite
    integer(c_int), value :: t_rangement
    integer(c_int), value :: num_var_cedre_max
    integer(c_int), value :: n_partition_local

  end subroutine PDM_io_tab_lec_debut


  !>
  !! \brief Definition d'une variable en ecriture
  !!
  !! \param [in] num_var_cedre          Numéro de variable PDM
  !! \param [in] num_indirection_cedre  Numéro d'indirection PDM
  !! \param [in] t_n_composantes        Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
  !! \param [in] n_composantes          Nombre de composantes pour chaque donnee
  !! \param [in] taille_donnee          Taille unitaire de la donnnee
  !!

  subroutine PDM_io_tab_lec_def_var (num_var_cedre,         &
                                     num_indirection_cedre, &
                                     t_n_composantes,       &
                                     n_composantes,         &
                                     taille_donnee)         &
  bind (c, name='PDM_io_tab_lec_def_var')
    use iso_c_binding
    implicit none

    integer(c_int), value :: num_var_cedre
    integer(c_int), value :: num_indirection_cedre
    integer(c_int), value :: t_n_composantes
    integer(c_int), value :: n_composantes
    integer(c_int), value :: taille_donnee

  end subroutine PDM_io_tab_lec_def_var


  !>
  !! \brief Finalise une phase de lecture parallèle de tableaux de données associées
  !! aux numéros de variables PDM. Cette fonction déclenche réellement
  !! les écritures
  !!

  subroutine PDM_io_tab_lec_fin () &
  bind (c, name='PDM_io_tab_lec_fin')
    use iso_c_binding
    implicit none

  end subroutine PDM_io_tab_lec_fin

  end interface

contains


!>
!! \brief Ouverture d'un fichier pour acces parallele
!!
!! \param [in]  nom             Nom du fichier
!! \param [in]  fmt             Fichier text ou binaire
!! \param [in]  suff_t          Type de suffixe (manuel ou automatique)
!! \param [in]  suff_u          Suffixe (si suffixe manuel)
!! \param [in]  s_backup        Active le backup d'un fichier preexistant en mode ecriture
!! \param [in]  accesio         Type (parallele avec mpiio, parallele sans mpiio, sequentiel)
!! \param [in]  mode            Mode d'acces (lecture, ecriture, lecture/ecriture)
!! \param [in]  pdm_mpi_comm    Communicateur lie au fichier
!! \param [out] unite           Unite du fichier
!! \param [out] ierr            Indique si le fichier est de type PDM_io ou non (uniquement pour une ouverture en lecture)
!!
!!

subroutine PDM_io_open (nom,                &
                        fmt,                &
                        suff_t,             &
                        suff_u,             &
                        s_backup,           &
                        acces,              &
                        mode,               &
                        endian,             &
                        comm,               &
                        prop_noeuds_actifs, &
                        unite,              &
                        ierr)
  use iso_c_binding
  implicit none

  character (len=*) :: nom
  integer           :: fmt
  integer           :: suff_t
  character (len=*) :: suff_u
  integer           :: s_backup
  integer           :: acces
  integer           :: mode
  integer           :: endian
  integer           :: comm
  double precision  :: prop_noeuds_actifs
  type(c_ptr)       :: unite
  integer           :: ierr

  integer(c_int)    :: c_fmt
  integer(c_int)    :: c_suff_t
  integer(c_int)    :: c_s_backup
  integer(c_int)    :: c_acces
  integer(c_int)    :: c_mode
  integer(c_int)    :: c_endian
  integer(c_int)    :: c_comm
  real(c_double)    :: c_prop_noeuds_actifs
  integer(c_int)    :: c_ierr

  interface
    subroutine PDM_io_open_c (nom,                &
                              fmt,                &
                              suff_t,             &
                              suff_u,             &
                              s_backup,           &
                              acces,              &
                              mode,               &
                              endian,             &
                              comm,               &
                              prop_noeuds_actifs, &
                              unite,              &
                              ierr)               &
    bind (c, name='PDM_io_open')
      use iso_c_binding
      implicit none

      character(c_char)     :: nom(*)
      integer(c_int), value :: fmt
      integer(c_int), value :: suff_t
      character(c_char)     :: suff_u(*)
      integer(c_int), value :: s_backup
      integer(c_int), value :: acces
      integer(c_int), value :: mode
      integer(c_int), value :: endian
      integer(c_int), value :: comm
      real(c_double), value :: prop_noeuds_actifs
      type(c_ptr)           :: unite
      integer(c_int)        :: ierr

    end subroutine PDM_io_open_c
  end interface

  c_fmt      = fmt
  c_suff_t   = suff_t
  c_s_backup = s_backup
  c_acces    = acces
  c_mode     = mode
  c_endian   = endian

  c_comm = PDM_MPI_Comm_f2c(comm)

  c_prop_noeuds_actifs = prop_noeuds_actifs

  call PDM_io_open_c (nom//C_NULL_CHAR,      &
                      c_fmt,                 &
                      c_suff_t,              &
                      suff_u//C_NULL_CHAR,   &
                      c_s_backup,            &
                      c_acces,               &
                      c_mode,                &
                      c_endian,              &
                      c_comm,                &
                      c_prop_noeuds_actifs,  &
                      unite,                 &
                      c_ierr)

  ierr = c_ierr

end subroutine PDM_io_open


!>
!! \brief Return the current file position
!!
!! \param [in]  fichier         Pointer to \ref PDM_io_fichier_t object
!! \param [out] offset          Current position in file
!!
!!

subroutine PDM_io_tell (fichier, &
                        offset)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: fichier
  integer(pdm_g_num_s), intent(out) :: offset

  interface
    function PDM_io_tell_c (fichier) &
    result (offset)                  &
    bind (c, name='PDM_io_tell')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: fichier
#ifdef PDM_LONG_G_NUM
      integer(c_long)    :: offset
#else
      integer(c_int)     :: offset
#endif

    end function PDM_io_tell_c
  end interface

  offset = PDM_io_tell_c (fichier)

end subroutine PDM_io_tell


!>
!! \brief Lecture parallele de blocs de donnees suivie d'une redistribution des
!! des donnees suivant l'indirection
!!
!! \param [in]  fichier          Pointer to \ref PDM_io_fichier_t object
!! \param [in]  t_n_composantes  Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
!! \param [in]  n_composantes    Nombre de composantes pour chaque donnee
!! \param [in]  taille_donnee    Taille unitaire de la donnee
!! \param [in]  n_donnees        Nombre de donnees a lire
!! \param [in]  indirection      Indirection de redistribition des donnees
!! \param [out] donnees          Donnees lues
!!
!!

subroutine PDM_io_lec_par_entrelacee (fichier,         &
                                      t_n_composantes, &
                                      n_composantes,   &
                                      taille_donnee,   &
                                      n_donnees,       &
                                      indirection,     &
                                      donnees)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  integer(pdm_l_num_s)          :: t_n_composantes
  integer(pdm_l_num_s), pointer :: n_composantes(:)
  integer(pdm_l_num_s)          :: taille_donnee
  integer(pdm_l_num_s)          :: n_donnees
  integer(pdm_g_num_s), pointer :: indirection(:)
  type(c_ptr)                   :: donnees

  integer(c_int)                :: c_t_n_composantes
  integer(c_int)                :: c_taille_donnee
  integer(c_int)                :: c_n_donnees
  type(c_ptr)                   :: c_n_composantes = C_NULL_PTR
  type(c_ptr)                   :: c_indirection   = C_NULL_PTR

  interface
    subroutine PDM_io_lec_par_entrelacee_c (fichier,         &
                                            t_n_composantes, &
                                            n_composantes,   &
                                            taille_donnee,   &
                                            n_donnees,       &
                                            indirection,     &
                                            donnees)         &
    bind (c, name='PDM_io_lec_par_entrelacee')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: fichier
      integer(c_int)     :: t_n_composantes
      type(c_ptr), value :: n_composantes
      integer(c_int)     :: taille_donnee
      integer(c_int)     :: n_donnees
      type(c_ptr), value :: indirection
      type(c_ptr)        :: donnees

    end subroutine PDM_io_lec_par_entrelacee_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_taille_donnee   = taille_donnee
  c_n_donnees       = n_donnees

  c_n_composantes = c_loc(n_composantes)
  c_indirection   = c_loc(indirection  )

  call PDM_io_lec_par_entrelacee_c (fichier,           &
                                    c_t_n_composantes, &
                                    c_n_composantes,   &
                                    c_taille_donnee,   &
                                    c_n_donnees,       &
                                    c_indirection,     &
                                    donnees)

end subroutine PDM_io_lec_par_entrelacee


!>
!! \brief Lecture parallele de blocs de donnees
!! Les blocs doivent etre ranges par ordre croissant suivant la numerotation
!! des processus
!!
!! \param [in]  fichier          Pointer to \ref PDM_io_fichier_t object
!! \param [in]  t_n_composantes  Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
!! \param [in]  n_composantes    Nombre de composantes pour chaque donnee
!! \param [in]  taille_donnee    Taille unitaire de la donnee
!! \param [in]  n_donnees        Nombre de donnees a lire
!! \param [in]  debut_bloc       Adresse relative du debut de bloc
!! \param [out] donnees          Donnees lues
!!
!!

subroutine PDM_io_lec_par_bloc (fichier,         &
                                t_n_composantes, &
                                n_composantes,   &
                                taille_donnee,   &
                                n_donnees,       &
                                debut_bloc,      &
                                donnees)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  integer(pdm_l_num_s)          :: t_n_composantes
  integer(pdm_l_num_s), pointer :: n_composantes(:)
  integer(pdm_l_num_s)          :: taille_donnee
  integer(pdm_l_num_s)          :: n_donnees
  integer(pdm_g_num_s)          :: debut_bloc
  type(c_ptr)                   :: donnees

  integer(c_int)                :: c_t_n_composantes
  integer(c_int)                :: c_taille_donnee
  integer(c_int)                :: c_n_donnees
  type(c_ptr)                   :: c_n_composantes = C_NULL_PTR
#ifdef PDM_LONG_G_NUM
  integer(c_long)               :: c_debut_bloc
#else
  integer(c_int)                :: c_debut_bloc
#endif

  interface
    subroutine PDM_io_lec_par_bloc_c (fichier,         &
                                      t_n_composantes, &
                                      n_composantes,   &
                                      taille_donnee,   &
                                      n_donnees,       &
                                      debut_bloc,      &
                                      donnees)         &
    bind (c, name='PDM_io_lec_par_bloc')
      use iso_c_binding
      implicit none

      type(c_ptr),     value :: fichier
      integer(c_int)         :: t_n_composantes
      type(c_ptr),     value :: n_composantes
      integer(c_int)         :: taille_donnee
      integer(c_int)         :: n_donnees
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: debut_bloc
#else
      integer(c_int),  value :: debut_bloc
#endif
      type(c_ptr)            :: donnees

    end subroutine PDM_io_lec_par_bloc_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_taille_donnee   = taille_donnee
  c_n_donnees       = n_donnees
  c_debut_bloc      = debut_bloc

  c_n_composantes = c_loc(n_composantes)

  call PDM_io_lec_par_bloc_c (fichier,           &
                              c_t_n_composantes, &
                              c_n_composantes,   &
                              c_taille_donnee,   &
                              c_n_donnees,       &
                              c_debut_bloc,      &
                              donnees)

end subroutine PDM_io_lec_par_bloc


!>
!! \brief Tri des donnees suivant l'indirection puis ecriture parallele des blocs de
!! donnees
!!
!! \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
!! \param [in] t_n_composantes   Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
!! \param [in] n_composantes     Nombre de composantes pour chaque donnee
!! \param [in] taille_donnee     Taille unitaire de la donnee
!! \param [in] n_donnees         Nombre de donnees a ecrire
!! \param [in] indirection       Indirection de redistribition des donnees
!! \param [in] donnees           Donnees a ecrire
!!
!!

subroutine PDM_io_ecr_par_entrelacee (fichier,         &
                                      t_n_composantes, &
                                      n_composantes,   &
                                      taille_donnee,   &
                                      n_donnees,       &
                                      indirection,     &
                                      donnees)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  integer(pdm_l_num_s)          :: t_n_composantes
  integer(pdm_l_num_s), pointer :: n_composantes(:)
  integer(pdm_l_num_s)          :: taille_donnee
  integer(pdm_l_num_s)          :: n_donnees
  integer(pdm_g_num_s), pointer :: indirection(:)
  type(c_ptr)                   :: donnees

  integer(c_int)                :: c_t_n_composantes
  integer(c_int)                :: c_taille_donnee
  integer(c_int)                :: c_n_donnees
  type(c_ptr)                   :: c_n_composantes = C_NULL_PTR
  type(c_ptr)                   :: c_indirection   = C_NULL_PTR

  interface
    subroutine PDM_io_ecr_par_entrelacee_c (fichier,         &
                                            t_n_composantes, &
                                            n_composantes,   &
                                            taille_donnee,   &
                                            n_donnees,       &
                                            indirection,     &
                                            donnees)         &
    bind (c, name='PDM_io_ecr_par_entrelacee')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: fichier
      integer(c_int)     :: t_n_composantes
      type(c_ptr), value :: n_composantes
      integer(c_int)     :: taille_donnee
      integer(c_int)     :: n_donnees
      type(c_ptr), value :: indirection
      type(c_ptr)        :: donnees

    end subroutine PDM_io_ecr_par_entrelacee_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_taille_donnee   = taille_donnee
  c_n_donnees       = n_donnees

  c_n_composantes = c_loc(n_composantes)
  c_indirection   = c_loc(indirection  )

  call PDM_io_ecr_par_entrelacee_c (fichier,           &
                                    c_t_n_composantes, &
                                    c_n_composantes,   &
                                    c_taille_donnee,   &
                                    c_n_donnees,       &
                                    c_indirection,     &
                                    donnees)

end subroutine PDM_io_ecr_par_entrelacee


!>
!! \brief Ecriture parallele de blocs de donnees
!! Les blocs doivent etre rangés par ordre croissant suivant la numérotation
!! des processus
!!
!! \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
!! \param [in] t_n_composantes   Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
!! \param [in] n_composantes     Nombre de composantes pour chaque donnee
!! \param [in] taille_donnee     Taille unitaire de la donnee
!! \param [in] debut_bloc        Adresse relative du debut de bloc
!! \param [in] n_donnees         Nombre de donnees a lire
!! \param [in] donnees           Donnees a ecrire
!!
!!

subroutine PDM_io_ecr_par_bloc (fichier,         &
                                t_n_composantes, &
                                n_composantes,   &
                                taille_donnee,   &
                                n_donnees,       &
                                debut_bloc,      &
                                donnees)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  integer(pdm_l_num_s)          :: t_n_composantes
  integer(pdm_l_num_s), pointer :: n_composantes(:)
  integer(pdm_l_num_s)          :: taille_donnee
  integer(pdm_l_num_s)          :: n_donnees
  integer(pdm_g_num_s)          :: debut_bloc
  type(c_ptr)                   :: donnees

  integer(c_int)                :: c_t_n_composantes
  integer(c_int)                :: c_taille_donnee
  integer(c_int)                :: c_n_donnees
  type(c_ptr)                   :: c_n_composantes = C_NULL_PTR
#ifdef PDM_LONG_G_NUM
  integer(c_long)               :: c_debut_bloc
#else
  integer(c_int)                :: c_debut_bloc
#endif

  interface
    subroutine PDM_io_ecr_par_bloc_c (fichier,         &
                                      t_n_composantes, &
                                      n_composantes,   &
                                      taille_donnee,   &
                                      n_donnees,       &
                                      debut_bloc,      &
                                      donnees)         &
    bind (c, name='PDM_io_ecr_par_bloc')
      use iso_c_binding
      implicit none

      type(c_ptr),     value :: fichier
      integer(c_int)         :: t_n_composantes
      type(c_ptr),     value :: n_composantes
      integer(c_int)         :: taille_donnee
      integer(c_int)         :: n_donnees
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: debut_bloc
#else
      integer(c_int),  value :: debut_bloc
#endif
      type(c_ptr)            :: donnees

    end subroutine PDM_io_ecr_par_bloc_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_taille_donnee   = taille_donnee
  c_n_donnees       = n_donnees
  c_debut_bloc      = debut_bloc

  c_n_composantes = c_loc(n_composantes)

  call PDM_io_ecr_par_bloc_c (fichier,           &
                              c_t_n_composantes, &
                              c_n_composantes,   &
                              c_taille_donnee,   &
                              c_n_donnees,       &
                              c_debut_bloc,      &
                              donnees)

end subroutine PDM_io_ecr_par_bloc


!>
!! \brief Retourne le temps cumule d'acces aux fichiers
!!
!! \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
!! \param [out] t_cpu             Temps CPU
!! \param [out] t_elapsed         Temps elapsed
!!
!!

subroutine PDM_io_get_timer_fichier (fichier,   &
                                     t_cpu,     &
                                     t_elapsed)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  double precision, intent(out) :: t_cpu
  double precision, intent(out) :: t_elapsed

  real(c_double)                :: c_t_cpu
  real(c_double)                :: c_t_elapsed

  interface

    subroutine PDM_io_get_timer_fichier_c (fichier,   &
                                           t_cpu,     &
                                           t_elapsed) &
    bind (c, name='PDM_io_get_timer_fichier')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: fichier
      real(c_double)     :: t_cpu
      real(c_double)     :: t_elapsed

    end subroutine PDM_io_get_timer_fichier_c

  end interface

    call PDM_io_get_timer_fichier_c (fichier,     &
                                     c_t_cpu,     &
                                     c_t_elapsed)

    t_cpu     = c_t_cpu
    t_elapsed = c_t_elapsed

end subroutine PDM_io_get_timer_fichier



!>
!! \brief Retourne le temps cumule pour le swap des donnees
!!
!! \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
!! \param [out] t_cpu             Temps CPU
!! \param [out] t_elapsed         Temps elapsed
!!
!!

subroutine PDM_io_get_timer_swap_endian (fichier,   &
                                         t_cpu,     &
                                         t_elapsed)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  double precision, intent(out) :: t_cpu
  double precision, intent(out) :: t_elapsed

  real(c_double)                :: c_t_cpu
  real(c_double)                :: c_t_elapsed

  interface

    subroutine PDM_io_get_timer_swap_endian_c (fichier,   &
                                               t_cpu,     &
                                               t_elapsed) &
    bind (c, name='PDM_io_get_timer_swap_endian')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: fichier
      real(c_double)     :: t_cpu
      real(c_double)     :: t_elapsed

    end subroutine PDM_io_get_timer_swap_endian_c

  end interface

    call PDM_io_get_timer_swap_endian_c (fichier,     &
                                         c_t_cpu,     &
                                         c_t_elapsed)

    t_cpu     = c_t_cpu
    t_elapsed = c_t_elapsed

end subroutine PDM_io_get_timer_swap_endian



!>
!! \brief Retourne le temps cumule pour la distribution des donnees
!!
!! \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
!! \param [out] t_cpu             Temps CPU
!! \param [out] t_elapsed         Temps elapsed
!!
!!

subroutine PDM_io_get_timer_distrib (fichier,   &
                                     t_cpu,     &
                                     t_elapsed)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  double precision, intent(out) :: t_cpu
  double precision, intent(out) :: t_elapsed

  real(c_double)                :: c_t_cpu
  real(c_double)                :: c_t_elapsed

  interface

    subroutine PDM_io_get_timer_distrib_c (fichier,   &
                                           t_cpu,     &
                                           t_elapsed) &
    bind (c, name='PDM_io_get_timer_distrib')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: fichier
      real(c_double)     :: t_cpu
      real(c_double)     :: t_elapsed

    end subroutine PDM_io_get_timer_distrib_c

  end interface

    call PDM_io_get_timer_distrib_c (fichier,     &
                                     c_t_cpu,     &
                                     c_t_elapsed)

    t_cpu     = c_t_cpu
    t_elapsed = c_t_elapsed

end subroutine PDM_io_get_timer_distrib



!>
!! \brief Retourne le temps cumule total
!!
!! \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
!! \param [out] t_cpu             Temps CPU
!! \param [out] t_elapsed         Temps elapsed
!!
!!

subroutine PDM_io_get_timer_total (fichier,   &
                                   t_cpu,     &
                                   t_elapsed)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  double precision, intent(out) :: t_cpu
  double precision, intent(out) :: t_elapsed

  real(c_double)                :: c_t_cpu
  real(c_double)                :: c_t_elapsed

  interface

    subroutine PDM_io_get_timer_total_c (fichier,   &
                                         t_cpu,     &
                                         t_elapsed) &
    bind (c, name='PDM_io_get_timer_total')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: fichier
      real(c_double)     :: t_cpu
      real(c_double)     :: t_elapsed

    end subroutine PDM_io_get_timer_total_c

  end interface

    call PDM_io_get_timer_total_c (fichier,     &
                                   c_t_cpu,     &
                                   c_t_elapsed)

    t_cpu     = c_t_cpu
    t_elapsed = c_t_elapsed

end subroutine PDM_io_get_timer_total


!>
!! \brief Swap endian pour conversion little endian <-> big endian
!!
!! \param [in]  taille_donnee   Taille unitaire de la donnee
!! \param [in]  n_donnees       Nombre de donnees
!! \param [in]  donnees         Donnees
!! \param [out] resultats       Resultat
!!
!!

subroutine PDM_io_swap_endian (taille_donnee, &
                               n_donnees,     &
                               donnees,       &
                               resultats)
  use iso_c_binding
  implicit none

  integer, intent(in) :: taille_donnee
  integer, intent(in) :: n_donnees
  type(c_ptr), value  :: donnees
  type(c_ptr)         :: resultats

  integer(c_size_t)   :: c_taille_donnee
  integer(c_size_t)   :: c_n_donnees

  interface
    subroutine PDM_io_swap_endian_c (taille_donnee, &
                                     n_donnees,     &
                                     donnees,       &
                                     resultats)     &
    bind (c, name='PDM_io_swap_endian')
      use iso_c_binding
      implicit none

      integer(c_size_t), value :: taille_donnee
      integer(c_size_t), value :: n_donnees
      type(c_ptr),       value :: donnees
      type(c_ptr)              :: resultats

    end subroutine PDM_io_swap_endian_c
  end interface

  c_taille_donnee = taille_donnee
  c_n_donnees     = n_donnees

  call PDM_io_swap_endian_c (c_taille_donnee, &
                             c_n_donnees,     &
                             donnees,         &
                             resultats)

end subroutine PDM_io_swap_endian


!>
!! \brief Définit le format de la donnée indviduelle pour la sortie text
!!
!! \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
!! \param [in]  n_char_fmt        Nombre de caractères du format
!! \param [in]  data_type         Type de donnees
!! \param [in]  fmt               Format
!!
!!

subroutine PDM_io_fmt_donnee_set (fichier,    &
                                  n_char_fmt, &
                                  data_type,  &
                                  fmt)
  use iso_c_binding
  implicit none

  type(c_ptr), value :: fichier
  integer            :: n_char_fmt
  integer            :: data_type
  character(len=*)   :: fmt

  integer(c_int)     :: c_n_char_fmt
  integer(c_int)     :: c_data_type

  interface
    subroutine PDM_io_fmt_donnee_set_c (fichier,    &
                                        n_char_fmt, &
                                        data_type,  &
                                        fmt)        &
    bind (c, name='PDM_io_fmt_donnee_set')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: fichier
      integer(c_int), value :: n_char_fmt
      integer(c_int), value :: data_type
      character(c_char)     :: fmt(*)

    end subroutine PDM_io_fmt_donnee_set_c
  end interface

  c_n_char_fmt = n_char_fmt
  c_data_type  = data_type

  call PDM_io_fmt_donnee_set_c (fichier,          &
                                c_n_char_fmt,     &
                                c_data_type,      &
                                fmt//C_NULL_CHAR)

end subroutine PDM_io_fmt_donnee_set



!>
!! \brief Create a directory
!!
!! \param [in]  path   Path to new directory
!! \param [out] code   0 if successful, -1 else
!!
!!

subroutine PDM_io_mkdir (path, &
                         code)
  use iso_c_binding
  implicit none

  character(len=*), intent(in)  :: path
  integer,          intent(out) :: code

  interface
    function PDM_io_mkdir_c (path) &
    result (code)                  &
    bind (c, name='PDM_io_mkdir')
      use iso_c_binding
      implicit none

      character(c_char) :: path(*)
      integer(c_int)    :: code

    end function PDM_io_mkdir_c
  end interface

  code = PDM_io_mkdir_c (path//C_NULL_CHAR)

end subroutine PDM_io_mkdir



!>
!! \brief Calcul de la taille totale d'un champ de donnees
!!
!! \param [in]  fichier          Pointer to \ref PDM_io_fichier_t object
!! \param [in]  t_n_composantes  Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
!! \param [in]  n_composantes    Nombre de composantes pour chaque donnee
!! \param [in]  n_donnees        Nombre de donnees
!! \param [in]  indirection      Indirection de redistribition des donnees
!! \param [out] taille           Taille totale d'un champ de donnees
!!
!!

subroutine PDM_io_n_donnees_get (fichier,         &
                                 t_n_composantes, &
                                 n_composantes,   &
                                 n_donnees,       &
                                 indirection,     &
                                 taille)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: fichier
  integer                       :: t_n_composantes
  integer(pdm_l_num_s), pointer :: n_composantes(:)
  integer(pdm_l_num_s)          :: n_donnees
  integer(pdm_g_num_s), pointer :: indirection(:)
  integer(pdm_g_num_s)          :: taille

  integer(c_int)                :: c_t_n_composantes
  type(c_ptr)                   :: c_n_composantes = C_NULL_PTR
  integer(c_int)                :: c_n_donnees
  type(c_ptr)                   :: c_indirection   = C_NULL_PTR
#ifdef PDM_LONG_G_NUM
  integer(c_long)               :: c_taille
#else
  integer(c_int)                :: c_taille
#endif

  interface
    function PDM_io_n_donnees_get_c (fichier,         &
                                     t_n_composantes, &
                                     n_composantes,   &
                                     n_donnees,       &
                                     indirection)     &
    result (taille)                                   &
    bind (c, name='PDM_io_n_donnees_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: fichier
      integer(c_int), value :: t_n_composantes
      type(c_ptr),    value :: n_composantes
      integer(c_int), value :: n_donnees
      type(c_ptr),    value :: indirection
#ifdef PDM_LONG_G_NUM
      integer(c_long)       :: taille
#else
      integer(c_int)        :: taille
#endif

    end function PDM_io_n_donnees_get_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_n_donnees       = n_donnees

  c_n_composantes = c_loc(n_composantes)
  c_indirection   = c_loc(indirection)

  c_taille = PDM_io_n_donnees_get_c (fichier,           &
                                     c_t_n_composantes, &
                                     c_n_composantes,   &
                                     c_n_donnees,       &
                                     c_indirection)

  taille = c_taille

end subroutine PDM_io_n_donnees_get


!>
!! \brief Ajoute une partie des donnees dans un tableau associés à une variable
!! PDM
!!
!! \param [in] num_var_cedre          Numéro de variable PDM
!! \param [in] i_part                 indice de partition
!! \param [in] n_composantes          Nombre de composantes pour chaque donnee
!! \param [in] n_donnees              Nombre de donnees a lire
!! \param [in] indirection            Indirection de redistribition des donnees
!! \param [in] donnees                Donnees a écrire
!!
!!

subroutine PDM_io_tab_ecr_ajout_donnees (num_var_cedre, &
                                         i_part,        &
                                         n_composantes, &
                                         n_donnees,     &
                                         indirection,   &
                                         donnees)
  use iso_c_binding
  implicit none

  integer(pdm_l_num_s)          :: num_var_cedre
  integer(pdm_l_num_s)          :: i_part
  integer(pdm_l_num_s), pointer :: n_composantes(:)
  integer(pdm_l_num_s)          :: n_donnees
  integer(pdm_g_num_s), pointer :: indirection(:)
  type(c_ptr), value            :: donnees

  type(c_ptr)                   :: c_n_composantes = C_NULL_PTR
  type(c_ptr)                   :: c_indirection   = C_NULL_PTR

  interface
    subroutine PDM_io_tab_ecr_ajout_donnees_c (num_var_cedre, &
                                               i_part,        &
                                               n_composantes, &
                                               n_donnees,     &
                                               indirection,   &
                                               donnees)       &
    bind (c, name='PDM_io_tab_ecr_ajout_donnees')
      use iso_c_binding
      implicit none

      integer(c_int), value :: num_var_cedre
      integer(c_int), value :: i_part
      type(c_ptr),    value :: n_composantes
      integer(c_int), value :: n_donnees
      type(c_ptr),    value :: indirection
      type(c_ptr),    value :: donnees

    end subroutine PDM_io_tab_ecr_ajout_donnees_c
  end interface

  c_n_composantes = c_loc(n_composantes)
  c_indirection   = c_loc(indirection  )

  call PDM_io_tab_ecr_ajout_donnees_c (num_var_cedre,   &
                                       i_part,          &
                                       c_n_composantes, &
                                       n_donnees,       &
                                       c_indirection,   &
                                       donnees)

end subroutine PDM_io_tab_ecr_ajout_donnees


!>
!! \brief Ajoute une partie des donnees dans un tableau associés à une variable PDM
!!
!! \param [in] num_var_cedre          Numéro de variable PDM
!! \param [in] i_part                 indice de partition
!! \param [in] n_composantes          Nombre de composantes pour chaque donnee
!! \param [in] n_donnees              Nombre de donnees a lire
!! \param [in] indirection            Indirection de redistribition des donnees
!! \param [in] donnees                Donnees a écrire
!!
!!

subroutine PDM_io_tab_lec_ajout_donnees (num_var_cedre, &
                                         i_part,        &
                                         n_composantes, &
                                         n_donnees,     &
                                         indirection,   &
                                         donnees)
  use iso_c_binding
  implicit none

  integer(pdm_l_num_s)          :: num_var_cedre
  integer(pdm_l_num_s)          :: i_part
  integer(pdm_l_num_s), pointer :: n_composantes(:)
  integer(pdm_l_num_s)          :: n_donnees
  integer(pdm_g_num_s), pointer :: indirection(:)
  type(c_ptr), value            :: donnees

  type(c_ptr)                   :: c_n_composantes = C_NULL_PTR
  type(c_ptr)                   :: c_indirection   = C_NULL_PTR

  interface
    subroutine PDM_io_tab_lec_ajout_donnees_c (num_var_cedre, &
                                               i_part,        &
                                               n_composantes, &
                                               n_donnees,     &
                                               indirection,   &
                                               donnees)       &
    bind (c, name='PDM_io_tab_lec_ajout_donnees')
      use iso_c_binding
      implicit none

      integer(c_int), value :: num_var_cedre
      integer(c_int), value :: i_part
      type(c_ptr),    value :: n_composantes
      integer(c_int), value :: n_donnees
      type(c_ptr),    value :: indirection
      type(c_ptr),    value :: donnees

    end subroutine PDM_io_tab_lec_ajout_donnees_c
  end interface

  c_n_composantes = c_loc(n_composantes)
  c_indirection   = c_loc(indirection  )

  call PDM_io_tab_lec_ajout_donnees_c (num_var_cedre,   &
                                       i_part,          &
                                       c_n_composantes, &
                                       n_donnees,       &
                                       c_indirection,   &
                                       donnees)

end subroutine PDM_io_tab_lec_ajout_donnees

end module pdm_io

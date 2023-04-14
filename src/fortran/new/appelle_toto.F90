subroutine appelle_toto(fun, i) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  interface
    subroutine fun(i) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: i
    ! integer, intent(in) :: i
    end subroutine fun
  end interface
  
  ! type(c_funptr) :: c_fun
  ! procedure(fun), pointer :: f_fun
  integer, value, intent(in) :: i
  
  print *, "retour Fortran"
  
  ! call c_f_procpointer(c_fun, f_fun)

  ! print *, "c_f_procpointer OK", loc(f_fun)
  
  ! call f_fun(i)
  print *, "i = ", i
  call fun(i)

end subroutine appelle_toto

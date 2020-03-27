module mVariable

   use mBase

   ! Declarations
   implicit none

  type variable
     real(nk),dimension(:),allocatable:: ax, ay, b, cx, cy, d , var
     integer  :: CL_var_bas,   CL_var_haut,   CL_var_droite,   CL_var_gauche
     real(nk) :: var_bas   ,   var_haut   ,   var_droite   ,   var_gauche
     real(nk) :: q_var_bas ,   q_var_haut ,   q_var_droite ,   q_var_gauche
     real(nk) :: eps_zero = 1.0e-13
  end type variable


contains



end module mVariable

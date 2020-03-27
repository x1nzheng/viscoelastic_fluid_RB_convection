
#include "definitions.h"



program Navier_Stokes

  use mBase
  use mPression
  use mVitesse
  use mSolver
  use mTimestep
  use mChampsprof
  use mOutils
  use mPoste_traitement
  use mHyperbolic_part
  use mViscoelastic_models

  implicit none

  integer :: kkk

!  integer :: it_write
  integer :: point_choisi_espace , pce_x, pce_y, pce_k
  real(nk):: t1, t2
!!!-----------------
  eps_zero = 1.0e-30
!!!-----------------

  call Lecteur_donnees
!!!--------------------------adimensionnement des CLs-------------------------------
  call CL_Physique
  call adimensionnement_CL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Rayleigh = Rayleigh*8._nk
!  Weissenberg = Weissenberg/4._nk*sqrt(Rayleigh)
!  print*, Rayleigh, Weissenberg
!  pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!-----------------------------------
#if (TYPE_ECOULEMENT==3 && VISCOELASTIC_MODELS > 0)
  Ma2 = sqrt(Rayleigh) * Weissenberg / Prandtl / (1.- Beta)
  Reynolds= sqrt(Rayleigh) / Prandtl
#endif 
!!!------------------------------------

#if ( VISCOELASTIC_MODELS == 0 )
  var%W(T )%nomb_ad = 1._nk  /sqrt(Rayleigh)    ;  var%W(T )%rapp_ad = 0.
  var%W(u1)%nomb_ad = Prandtl/sqrt(Rayleigh)    ;  var%W(u1)%rapp_ad = 0.
  var%W(u2)%nomb_ad = Prandtl/sqrt(Rayleigh)    ;  var%W(u2)%rapp_ad = Prandtl
#if (DIMENSION_GEO == 3)     
  var%W(u3)%nomb_ad = Prandtl/sqrt(Rayleigh)    ;  var%W(u3)%rapp_ad = 0.
#endif

#elif ( VISCOELASTIC_MODELS > 0 )
      var%W(T )%nomb_ad = 1._nk       /sqrt(Rayleigh)    ;  var%W(T )%rapp_ad = 0.
      var%W(u1)%nomb_ad = Beta*Prandtl/sqrt(Rayleigh)    ;  var%W(u1)%rapp_ad = 0.
      var%W(u2)%nomb_ad = Beta*Prandtl/sqrt(Rayleigh)    ;  var%W(u2)%rapp_ad = Prandtl
#if (DIMENSION_GEO == 3)     
      var%W(u3)%nomb_ad = Beta*Prandtl/sqrt(Rayleigh)    ;  var%W(u3)%rapp_ad = 0.
#endif
#endif 
  var%W(P )%nomb_ad = 1.;  var%W(P )%rapp_ad = 0. 
#if (DIMENSION_GEO == 2 && VISCOELASTIC_MODELS > 0)
  do i=t11, t22       
     var%W(i)%nomb_ad = 1._nk / Weissenberg
  end do
#elif (DIMENSION_GEO == 3 && VISCOELASTIC_MODELS > 0)
  do i=t11, t33        
     var%W(i)%nomb_ad = 1._nk / Weissenberg
  end do
#endif


!!!----------------------------------------------------------------------------------
!!!----------------------------------------------------------------------------------
!!!----------------------------------------------------------------------------------
  allocate(time_size(nbr_it_temps) , grid_size(nbr_it_espace))
  call allocate_Erreur_var( Er )

  do k=1, 3
     call allocate_Erreur_var( Er_norme%norme(k) )
  end do

  call openfile 
 do kkk = 1, nbr_it_espace

#if (DIMENSION_GEO == 2)     
     do i=1,p-1
        var%W(i)%NNN    =  N1      * N2
        var%W(i)%NI     =  N1
        var%W(i)%NJ     =  N2
     end do
     var%W(p)%NNN    = (N1-1)*(N2-1)
     var%W(p)%NI     =  N1-1
     var%W(p)%NJ     =  N2-1
#elif (DIMENSION_GEO == 3)     
     do i=1,p-1
        var%W(i)%NNN    =  N1       *N2       *N3
        var%W(i)%NI     =  N1
        var%W(i)%NJ     =  N2
        var%W(i)%NK     =  N3
     end do
     var%W(p)%NNN  = (N1-1)*(N2-1)*(N3-1) 
     var%W(p)%NI     =  N1-1
     var%W(p)%NJ     =  N2-1    
     var%W(p)%NK     =  N3-1 
#endif

     do i=1,p-1
        var%W(i)%indice =  i
     end do

     do i=1,p
        call allocate_var(var%W(i))
     end do
     call ALLOUER
     call allocate_2
#if (DIMENSION_GEO == 3 && ENRG  == 1 )
     allocate( T_moy(N1) , Vitesse_u1_moy(N1) )
     T_moy = 0. ; Vitesse_u1_moy = 0. 
#endif

     call grille_et_pas_maillage
     call reperage_point_choisi(point_choisi_espace, pce_x, pce_y, pce_k)
     call discrete_coef
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call resolution_semi_implicite ( var )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if (DIMENSION_GEO == 2)   
     if (nbr_it_espace /=1) call calcul_erreur_espace_ponctuelle(Er, var, kkk, point_choisi_espace, cord_x1_0, cord_x2_0, 0, 0) 
#elif (DIMENSION_GEO == 3)   
     if (nbr_it_espace /=1) call calcul_erreur_espace_ponctuelle(Er, var, kkk, point_choisi_espace, &
          cord_x1_0, cord_x2_0,  1 ,  0, cord_x3_0)     
#endif 

     if ((nbr_it_temps /=1 .or. sonde ==1) .and. nbr_it_espace == 1 ) call closefile
     if ( nbr_it_temps /=1) call erreur_en_champs(Er_norme, 1, 0, 1, 0)

     do i=1,p
        call deallocate_var(var%W(i))
     end do
     call DEALLOUER
     call deallocate_2
     call dediscrete_coef


#if (DIMENSION_GEO == 3 && ENRG  == 1 )
     deallocate( T_moy , Vitesse_u1_moy )
#endif
  end do



  if ((nbr_it_espace /= 1 .or. sonde ==1) .and. nbr_it_temps == 1 ) call closefile
  if (nbr_it_espace  /=1) call erreur_en_champs(Er_norme, 0, 1, 1, 0)


  close(11)
  deallocate(time_size, grid_size)
  call deallocate_Erreur_var(Er)
  do k=1, 3
     call deallocate_Erreur_var(Er_norme%norme(k))
  end do
!!!----------------------------------------------------------------------  
!!!----------------------------------------------------------------------
!!!----------------------------------------------------------------------

contains































  subroutine resolution_semi_implicite ( var )
    implicit none

    type(vecteur_variable):: var
    integer :: i, j, k
    integer :: kk
    real(nk):: eps

!!!-------------------------------------------------------------------------------
    do kk = 1, nbr_it_temps  
       pas_write  = pas_write0

       do i= 1, p
          var%W(i)%var = 0.  ; var%W(i)%var0 = 0.  ; var%W(i)%var1 = 0.  ; var%W(i)%var10 = 0.          
       end do
#if (DIMENSION_GEO == 2 && VISCOELASTIC_MODELS > 0 )
       do j=3, 5
          var%W(j)%b = 1.0_nk 
       end do
#elif (DIMENSION_GEO == 3 && VISCOELASTIC_MODELS > 0 )
       do j=4, 9
          var%W(j)%b = 1.0_nk 
       end do
#endif

#if ( CONTINUE_CAL == 1 )
#if ( DIMENSION_GEO == 2 )
     call Read_initial_data_UT( var%W(u1), var%W(u2), var%W(T), var%W(t11), var%W(t12), var%W(t22))
     call Read_initial_data_P( var%W(P) )
     call champs_uvpT  (var%W(u1), var%W(u2), var%W(p), var%W(T), &
                        var%W(t11), var%W(t12), var%W(t22), 5, 1, it)  
#if ( CONTINUE_CAL == 1 )
#if ( EXP_TREATMENT == 1 )
     Tau_trans = 1.
     var%W(t11)%var0 = var%W(t11)%var0*exp(-dt/weissenberg)
     var%W(t12)%var0 = var%W(t12)%var0*exp(-dt/weissenberg)
     var%W(t22)%var0 = var%W(t22)%var0*exp(-dt/weissenberg)
#endif
     call calcaul_termes_NL_partie_hyperbolique (var%W(u1)%var0, var%W(u2)%var0, var%W(t11)%var0, &
                                                 var%W(t12)%var0, var%W(t22)%var0)
#endif
#elif (DIMENSION_GEO == 3)
     call Read_initial_data_UT( var%W(u1), var%W(u2), var%W(T), &
                 var%W(t11), var%W(t12), var%W(t22), var%W(u3), var%W(t13), var%W(t23), var%W(t33))
     call Read_initial_data_P( var%W(P) )
     call champs_uvpT  (var%W(u1), var%W(u2), var%W(p), var%W(T),var%W(t11), var%W(t12), &
                        var%W(t22), 5, 1, it, var%W(u3), var%W(t13), var%W(t23), var%W(t33)) 
                        
     call calcaul_termes_NL_partie_hyperbolique (var%W(u1)%var0, var%W(u2)%var0, var%W(t11)%var0, &
                                                 var%W(t12)%var0, var%W(t22)%var0, var%W(u3)%var0, &
                                                 var%W(t13)%var0, var%W(t23)%var0, var%W(t33)%var0)
#endif
#endif

#if ( ENRG == 1 )
#if ( CONTINUE_CAL == 0 )
       call initialisation_champ_T(var%W(T ))
#endif  
       call UT_points_ficts(var%W(u1))
       call UT_points_ficts(var%W(u2))
       call  P_points_ficts(var%W(p ))
       call UT_points_ficts(var%W(T )) 
       call Tau_points_ficts
#if (DIMENSION_GEO ==3)
       call UT_points_ficts(var%W(u3))
#endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! DICHOTOMIE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          it  = 1 ; eps = 1000.
          do while  ( eps > eps_zero .and. it < Total_step)  
#if ( DIMENSION_GEO == 2 )
       call timestep (var%W(u1), var%W(u2), var%W(p), var%W(T), var%W(t11), var%W(t12), var%W(t22), dt, it)
       if ( mod(it, 1) == 0) then
        call test_div_U_2D(var%W(u1)%var, var%W(u2)%var)
        call convergence_stationnaire(var%W(u1), var%W(u2),var%W(t11),var%W(t12),var%W(t22),&
             var%W(P), var%W(T), eps, kk, kkk)
       endif
#elif ( DIMENSION_GEO == 3) 
       call timestep (var%W(u1), var%W(u2), var%W(p), var%W(T), var%W(t11), var%W(t12), var%W(t22), dt, it, &
            var%W(u3), var%W(t13), var%W(t23), var%W(t33))
       if ( mod(it, 1) == 0) then
        call test_div_U_3D(var%W(u1)%var, var%W(u2)%var, var%W(u3)%var)
        call convergence_stationnaire(var%W(u1), var%W(u2),var%W(t11),var%W(t12),var%W(t22),&
             var%W(P), var%W(T), eps, kk, kkk, var%W(u3), var%W(t13), var%W(t23), var%W(t33))
       endif
#endif       
!!!-----------------
          if ((nbr_it_temps /=1 .or. nbr_it_espace /= 1) .and. it == pas_choisi_temps ) then
          if (nbr_it_temps /=1) call calcul_erreur_temps_ponctuelle(Er, var, kk, pce_x, pce_y, pce_k) 
             exit
          end if


          if (it == pas_write) then
#if (DIMENSION_GEO == 2)                                             
          call champs_uvpT  (var%W(u1), var%W(u2), var%W(p), var%W(T), &
                             var%W(t11), var%W(t12), var%W(t22), 5, 1, it)  
#elif (DIMENSION_GEO ==3)
          call champs_uvpT  (var%W(u1), var%W(u2), var%W(p), var%W(T),var%W(t11), var%W(t12), &
                        var%W(t22), 5, 1, it, var%W(u3), var%W(t13), var%W(t23), var%W(t33))   
#endif          
             pas_write = pas_write + pas_write0 
          end if
!!!!!!!!!!!!!call calcul_nombre_de_Nuesselt_with_constant_temperature   
          call Calculs_Average_Nusselt ( var%W(T)%var )            
!!!!!!!!!!!!!call calcul_nombre_de_Nuesselt_with_wall_heat_flux
!!!--------------------------------------------
          it= it + 1                         
       end do
!!!--------------------------------------------
!       if (nbr_it_temps /=1) then  
!          dt               = dt/2._nk
!          pas_choisi_temps = pas_choisi_temps*2
!       end if
    end do
#if (DIMENSION_GEO == 2)                                             
     call champs_uvpT  (var%W(u1), var%W(u2), var%W(p), var%W(T), &
                        var%W(t11), var%W(t12), var%W(t22), 5, 1, it)   
#elif (DIMENSION_GEO == 3)    
     call champs_uvpT  (var%W(u1), var%W(u2), var%W(p), var%W(T),var%W(t11), var%W(t12), &
                        var%W(t22), 5, 1, it, var%W(u3), var%W(t13), var%W(t23), var%W(t33))  
#endif   
    
  end subroutine Resolution_semi_implicite



end program Navier_Stokes


#include "definitions.h"

module mOutils
  use mBase
  implicit none

contains



  subroutine convergence_stationnaire (u1, u2, t11, t12, t22, P, T, eps, kk1, kkk1, u3, t13, t23, t33)
    implicit none
    integer :: i, j, k
    type(variable),intent(inout)         :: u1, u2, t11, t12, t22,  P, T
    type(variable),intent(inout),optional:: u3, t13, t23, t33
    real(nk),intent(out)                 :: eps
    integer ,intent(in )                 :: kk1, kkk1
    real(nk):: somme_u1, somme_u2, somme_P, somme_T, eps_u1, eps_u2, eps_P, eps_T,&
         somme_t11, somme_t12, somme_t22, eps_t11, eps_t12, eps_t22
 


#if(  DIMENSION_GEO == 2)
    real(nk),dimension(1:N1,1:N2):: t11_transvaser, t22_transvaser
    real(nk):: min_t11=0., min_t22=0., max_t11=0., max_t22=0.
    integer::  x1_x2_min_t11_x,  x1_x2_min_t22_x, x1_x2_max_t11_x, x1_x2_max_t22_x
    integer::  x1_x2_min_t11_y,  x1_x2_min_t22_y, x1_x2_max_t11_y, x1_x2_max_t22_y

#elif ( DIMENSION_GEO == 3)
    real(nk),dimension(1:N1,1:N2,1:N3):: t11_transvaser, t22_transvaser, t33_transvaser
    real(nk):: min_t11=0., min_t22=0., min_t33=0., max_t11=0., max_t22=0., max_t33=0.
    integer::  x1_x2_x3_min_t11_x,  x1_x2_x3_min_t22_x, x1_x2_x3_min_t33_x
    integer::  x1_x2_x3_min_t11_y,  x1_x2_x3_min_t22_y, x1_x2_x3_min_t33_y
    integer::  x1_x2_x3_min_t11_z,  x1_x2_x3_min_t22_z, x1_x2_x3_min_t33_z

#endif    

#if (DIMENSION_GEO == 3)
    real(nk):: somme_u3, eps_u3, somme_t13, somme_t23, somme_t33, eps_t13, eps_t23, eps_t33
#endif   


    real(nk):: max_valeur_propre_x1=0., max_valeur_propre_x2=0. 


    !! transvasement et recherche des positions de valeurs minimales


#if(  DIMENSION_GEO == 2)

    min_t11 = t11%var(1,1) ; max_t11 = t11%var(1,1)
    min_t22 = t22%var(1,1) ; max_t22 = t22%var(1,1)
    x1_x2_min_t11_x = 1
    x1_x2_min_t11_y = 1
    x1_x2_min_t22_x = 1    
    x1_x2_min_t22_y = 1 

    do j= 1, N2
       do i=1, N1 
          t11_transvaser(i,j) =  t11%var(i,j)
          t22_transvaser(i,j) =  t22%var(i,j)          
          if (min_t11 >  t11%var(i,j) ) then
             min_t11         = t11%var(i,j)
             x1_x2_min_t11_x = i
             x1_x2_min_t11_y = j
          end if
          if (max_t11 <  t11%var(i,j) ) then
             max_t11         = t11%var(i,j)
             x1_x2_max_t11_x = i
             x1_x2_max_t11_y = j
          end if

          if (min_t22 >  t22%var(i,j) ) then
             min_t22         = t22%var(i,j)
             x1_x2_min_t22_x = i   
             x1_x2_min_t22_y = j
          end if
          if (max_t22 <  t22%var(i,j) ) then
             max_t22         = t22%var(i,j)
             x1_x2_max_t22_x = i
             x1_x2_max_t22_y = j
         end if

       end do
    end do

#elif ( DIMENSION_GEO == 3)

    i = 1 
    j = 1
    k = 1
    min_t11 = t11%var(i,j,k) 
    min_t22 = t22%var(i,j,k) 
    min_t33 = t33%var(i,j,k) 
    x1_x2_x3_min_t11_x = 1
    x1_x2_x3_min_t11_y = 1
    x1_x2_x3_min_t11_z = 1
    x1_x2_x3_min_t22_x = 1
    x1_x2_x3_min_t22_y = 1
    x1_x2_x3_min_t22_z = 1
    x1_x2_x3_min_t33_x = 1
    x1_x2_x3_min_t33_y = 1
    x1_x2_x3_min_t33_z = 1

    do k= 1, N3
       do j= 1, N2
          do i=1, N1   

             t11_transvaser(i,j,k) =  t11%var(i,j,k)
             t22_transvaser(i,j,k) =  t22%var(i,j,k)
             t33_transvaser(i,j,k) =  t33%var(i,j,k)

             if (min_t11 >  t11%var(i,j,k) ) then
                min_t11          = t11%var(i,j,k)
                x1_x2_x3_min_t11_x = i
                x1_x2_x3_min_t11_y = j
                x1_x2_x3_min_t11_z = k
             end if

             if (min_t22 >  t22%var(i,j,k) ) then
                min_t22       = t22%var(i,j,k)
                x1_x2_x3_min_t22_x = i
                x1_x2_x3_min_t22_y = j
                x1_x2_x3_min_t22_z = k
             end if

             if (min_t33 >  t33%var(i,j,k) ) then
                min_t33       = t33%var(i,j,k)
                x1_x2_x3_min_t33_x = i
                x1_x2_x3_min_t33_y = j
                x1_x2_x3_min_t33_z = k
             end if

          end do
       end do
    end do

#endif

    somme_u1   = 0.0
    somme_u2   = 0.0
    somme_t11  = 0.0
    somme_t12  = 0.0
    somme_t22  = 0.0
    somme_P    = 0.0
    somme_T    = 0.0
#if (DIMENSION_GEO == 3)
    somme_u3   = 0.0
    somme_t13  = 0.0
    somme_t23  = 0.0
    somme_t33  = 0.0
#endif



#if (DIMENSION_GEO == 2)
    do j= 1, N2
       do i=1, N1     
          somme_u1 = somme_u1 + ( abs( u1%var(i,j) - u1%var10(i,j) ) ) 
          somme_u2 = somme_u2 + ( abs( u2%var(i,j) - u2%var10(i,j) ) ) 
          somme_T  = somme_T  + ( abs(  T%var(i,j) -  T%var10(i,j) ) ) 
#if(VISCOELASTIC_MODELS > 0)
          somme_t11 = somme_t11 + ( abs( t11%var(i,j) - t11%var10(i,j) ) ) 
          somme_t12 = somme_t12 + ( abs( t12%var(i,j) - t12%var10(i,j) ) ) 
          somme_t22 = somme_t22 + ( abs( t22%var(i,j) - t22%var10(i,j) ) ) 
#endif      
       end do
    end do
    do j= 1, N2-1
       do i=1, N1-1 
          somme_P  = somme_P  + (P%var(i,j) -  P%var10(i,j)) * (P%var(i,j) -  P%var10(i,j))
       end do
    end do
#elif (DIMENSION_GEO == 3)
    do k=1, N3
       do j= 1, N2
          do i=1, N1          
             somme_u1 = somme_u1 + ( abs( u1%var(i,j,k) - u1%var10(i,j,k) ) ) 
             somme_u2 = somme_u2 + ( abs( u2%var(i,j,k) - u2%var10(i,j,k) ) ) 
             somme_T  = somme_T  + ( abs(  T%var(i,j,k) -  T%var10(i,j,k) ) ) 
             somme_u3 = somme_u3 + ( abs( u3%var(i,j,k) - u3%var10(i,j,k) ) ) 
#if(VISCOELASTIC_MODELS > 0)
             somme_t11 = somme_t11 + ( abs( t11%var(i,j,k) - t11%var10(i,j,k) ) )
             somme_t12 = somme_t12 + ( abs( t12%var(i,j,k) - t12%var10(i,j,k) ) )
             somme_t13 = somme_t13 + ( abs( t13%var(i,j,k) - t13%var10(i,j,k) ) )
             somme_t22 = somme_t22 + ( abs( t22%var(i,j,k) - t22%var10(i,j,k) ) )
             somme_t23 = somme_t23 + ( abs( t23%var(i,j,k) - t23%var10(i,j,k) ) )
             somme_t33 = somme_t33 + ( abs( t33%var(i,j,k) - t33%var10(i,j,k) ) )
#endif      
          end do
       end do
    end do

    eps_u3  = somme_u3  /u1%NNN    
    eps_t13 = somme_t13 /u1%NNN    
    eps_t23 = somme_t23 /u1%NNN    
    eps_t33 = somme_t33 /u1%NNN    

    do k=1, N3 - 1
       do j= 1, N2 - 1 
          do i=1, N1 - 1
             somme_P  = somme_P  +  (P%var(i,j,k) -  P%var10(i,j,k)) * (P%var(i,j,k) -  P%var10(i,j,k))  
          end do
       end do
    end do
#endif

    eps_u1  = somme_u1  / u1%NNN
    eps_u2  = somme_u2  / u1%NNN
    eps_t11 = somme_t11 / u1%NNN
    eps_t12 = somme_t12 / u1%NNN
    eps_t22 = somme_t22 / u1%NNN
    eps_P   = SQRT( somme_P )  !/  P%NNN
    eps_T   = somme_T   /  T%NNN




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  PRINT INFORMATION ON TERMINAL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",' u1',nbr_it_temps,kk1,nbr_it_espace,kkk1,it,u1%niveau_conver,eps_u1
    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",' u2',nbr_it_temps,kk1,nbr_it_espace,kkk1,it,u2%niveau_conver,eps_u2
    u1%var10 = u1%var
    u2%var10 = u2%var

#if (DIMENSION_GEO == 3)
    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",' u3',nbr_it_temps,kk1,nbr_it_espace,kkk1,it,u3%niveau_conver,eps_u3
    u3%var10 = u3%var
#endif


#if(VISCOELASTIC_MODELS > 0)
    t11%var10 = t11%var
    t12%var10 = t12%var
    t22%var10 = t22%var

    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",'t11',nbr_it_temps,kk1,nbr_it_espace,kkk1,&
         it,t11%niveau_conver,eps_t11
    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",'t12',nbr_it_temps,kk1,nbr_it_espace,kkk1,&
         it,t11%niveau_conver,eps_t12            
    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",'t22',nbr_it_temps,kk1,nbr_it_espace,kkk1,&
         it,t11%niveau_conver,eps_t22

#if (DIMENSION_GEO == 3)
    t13%var10 = t13%var
    t23%var10 = t23%var
    t33%var10 = t33%var

    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",'t13',nbr_it_temps,kk1,nbr_it_espace,kkk1,&
         it,t13%niveau_conver,eps_t13
    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",'t23',nbr_it_temps,kk1,nbr_it_espace,kkk1,&
         it,t23%niveau_conver,eps_t23
    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",'t33',nbr_it_temps,kk1,nbr_it_espace,kkk1,&
         it,t33%niveau_conver,eps_t33
#endif

#endif

    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",' P ',nbr_it_temps,kk1,nbr_it_espace,kkk1,it,P%niveau_conver,eps_P
    P%var10  = P%var

#if ( ENRG  == 1)    
    if (print_erreur_temps ==1 ) print"(a,4(i4),i7,2(E24.10))",' T ',nbr_it_temps,kk1,nbr_it_espace,kkk1,it,T%niveau_conver,eps_T
    T%var10  = T%var
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (print_erreur_temps ==1 )  print*, ' '
    if (print_erreur_temps ==1 ) print"(a,E17.10,5x, a,E17.10)",'div_U_coins_coins=', div_U_sur_coins_avec_coins,&
         'div_U_centres_Uinter=', div_U_sur_centres_avec_Uinter , 'div_U_centres_coins=', div_U_sur_centres_avec_coins 
#if ( EXP_TREAMENT == 0)
#if(VISCOELASTIC_MODELS > 0  &&  DIMENSION_GEO == 2)
    if (print_erreur_temps ==1 ) print"(3(a,E17.10,5x))", '1/Ma2', 1./Ma2, &
        'min tau_11   ', minval(t11_transvaser), 'max tau_11   ', maxval(t11_transvaser)
    if (print_erreur_temps ==1 ) print"(3(a,E17.10,5x))", '1/Ma2', 1./Ma2, &
        'min tau_22   ', minval(t22_transvaser), 'max tau_22   ', maxval(t22_transvaser)
          
#elif (VISCOELASTIC_MODELS > 0  && DIMENSION_GEO == 3)
    if (print_erreur_temps ==1 ) print"(4(a,E12.4,5x))", '1/Ma2', 1./Ma2, &
        'min tau_11   ', minval(t11_transvaser), 'max tau_11   ', maxval(t11_transvaser)
    if (print_erreur_temps ==1 ) print"(4(a,E12.4,5x))", '1/Ma2', 1./Ma2, &
        'min tau_22   ', minval(t22_transvaser), 'max tau_22   ', maxval(t22_transvaser)
    if (print_erreur_temps ==1 ) print"(4(a,E12.4,5x))", '1/Ma2', 1./Ma2, &
        'min tau_33   ', minval(t33_transvaser), 'max tau_33   ', maxval(t33_transvaser)
#endif
#elif ( EXP_TREAMENT == 1 )
#if(VISCOELASTIC_MODELS > 0  &&  DIMENSION_GEO == 2)
    if (print_erreur_temps ==1 ) print"(3(a,E17.10,5x))", '1/Ma2', 1./Ma2, &
        'min tau_11   ', minval(t11_transvaser)/Tau_trans, 'max tau_11   ', maxval(t11_transvaser)/Tau_trans
    if (print_erreur_temps ==1 ) print"(3(a,E17.10,5x))", '1/Ma2', 1./Ma2, &
        'min tau_22   ', minval(t22_transvaser)/Tau_trans, 'max tau_22   ', maxval(t22_transvaser)/Tau_trans
          
#elif (VISCOELASTIC_MODELS > 0  && DIMENSION_GEO == 3)
    if (print_erreur_temps ==1 ) print"(4(a,E12.4,5x))", '1/Ma2', 1./Ma2, &
        'min tau_11   ', minval(t11_transvaser), 'max tau_11   ', maxval(t11_transvaser)
    if (print_erreur_temps ==1 ) print"(4(a,E12.4,5x))", '1/Ma2', 1./Ma2, &
        'min tau_22   ', minval(t22_transvaser)/Tau_trans, 'max tau_22   ', maxval(t22_transvaser)/Tau_trans
    if (print_erreur_temps ==1 ) print"(4(a,E12.4,5x))", '1/Ma2', 1./Ma2, &
        'min tau_33   ', minval(t33_transvaser)/Tau_trans, 'max tau_33   ', maxval(t33_transvaser)/Tau_trans
#endif
#endif

    if (print_erreur_temps ==1 )  print*




#if (DIMENSION_GEO == 2 && ENRG  == 1 && VISCOELASTIC_MODELS == 0 )    
    eps = max(abs(eps_u1), abs(eps_u2), abs(eps_T))    
#elif (DIMENSION_GEO == 2 && ENRG  == 1 && VISCOELASTIC_MODELS > 0)    
    eps = max(abs(eps_u1), abs(eps_u2), abs(eps_T), abs(eps_t11), abs(eps_t12), abs(eps_t22))    
#elif (DIMENSION_GEO == 2 && ENRG  == 0  && VISCOELASTIC_MODELS == 0 )    
    eps = max(abs(eps_u1), abs(eps_u2) , abs(eps_P))    
#elif (DIMENSION_GEO == 2 && ENRG  == 0  && VISCOELASTIC_MODELS > 0 )    
    eps = max(abs(eps_u1), abs(eps_u2), abs(eps_t11), abs(eps_t12), abs(eps_t22), abs(eps_P))    
#endif


#if (DIMENSION_GEO == 3 && ENRG  == 1 && VISCOELASTIC_MODELS == 0)    
    eps = max(abs(eps_u1), abs(eps_u2), abs(eps_u3), abs(eps_T)) !, abs(eps_P))    
#elif (DIMENSION_GEO == 3 && ENRG  == 1  && VISCOELASTIC_MODELS > 0)    
    eps = max(abs(eps_u1), abs(eps_u2), abs(eps_u3),  abs(eps_T), abs(eps_t11), abs(eps_t12), &
         abs(eps_t13), abs(eps_t22), abs(eps_t23), abs(eps_t33)) !, abs(eps_P))
#elif (DIMENSION_GEO == 3 && ENRG  == 0  && VISCOELASTIC_MODELS == 0)    
    eps = max(abs(eps_u1), abs(eps_u2), abs(eps_u3)) !, abs(eps_P))    
#elif (DIMENSION_GEO == 3 && ENRG  == 0  && VISCOELASTIC_MODELS > 0)    
    eps = max(abs(eps_u1), abs(eps_u2), abs(eps_u3), abs(eps_t11), abs(eps_t12), abs(eps_t13),&
         abs(eps_t22), abs(eps_t23), abs(eps_t33))!, abs(eps_P))    
#endif    
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
!    call Calculs_Average_Nusselt ( T )

#if ( DIMENSION_GEO == 2 )

   call kinetic_energy_budget( u1, u2, t11, t12, t22, T, P )

      write(30,*) it, it*dt, u1 %var(N1/2+1,N2/2+1), u2 %var(N1/2+1,N2/2+1),&
                             t11%var(N1/2+1,N2/2+1), t12%var(N1/2+1,N2/2+1),&
                             t22%var(N1/2+1,N2/2+1), T  %var(N1/2+1,N2/2+1)
            
      write(31,*) it, it*dt, u1 %var(N1/8*3+1,N2/2+1), u2 %var(N1/8*3+1,N2/2+1),&
                             t11%var(N1/8*3+1,N2/2+1), t12%var(N1/8*3+1,N2/2+1),&
                             t22%var(N1/8*3+1,N2/2+1), T  %var(N1/8*3+1,N2/2+1)

      write(32,*) it, it*dt, u1 %var(N1/8+1,N2/2+1), u2 %var(N1/8+1,N2/2+1),&
                             t11%var(N1/8+1,N2/2+1), t12%var(N1/8+1,N2/2+1),&
                             t22%var(N1/8+1,N2/2+1), T  %var(N1/8+1,N2/2+1)

      write(33,*) it, it*dt, u1 %var(N1/4+1,N2/4+1), u2 %var(N1/4+1,N2/4+1),&
                             t11%var(N1/4+1,N2/4+1), t12%var(N1/4+1,N2/4+1),&
                             t22%var(N1/4+1,N2/4+1), T  %var(N1/4+1,N2/4+1)
#elif ( DIMENSION_GEO == 3 )
      write(30,*) it, it*dt, u1 %var(N1/2+1,N2/2+1,N3/2+1), u2 %var(N1/2+1,N2/2+1,N3/2+1),&
                             u3 %var(N1/2+1,N2/2+1,N3/2+1), t11%var(N1/2+1,N2/2+1,N3/2+1),&
                             t12%var(N1/2+1,N2/2+1,N3/2+1), t13%var(N1/2+1,N2/2+1,N3/2+1),&
                             t22%var(N1/2+1,N2/2+1,N3/2+1), t23%var(N1/2+1,N2/2+1,N3/2+1),&
                             t33%var(N1/2+1,N2/2+1,N3/2+1), T  %var(N1/2+1,N2/2+1,N3/2+1)
            
      write(31,*) it, it*dt, u1 %var(N1/4+1,N2/4+1,N3/4+1), u2 %var(N1/4+1,N2/4+1,N3/4+1),&
                             u3 %var(N1/4+1,N2/4+1,N3/4+1), t11%var(N1/4+1,N2/4+1,N3/4+1),&
                             t12%var(N1/4+1,N2/4+1,N3/4+1), t13%var(N1/4+1,N2/4+1,N3/4+1),&
                             t22%var(N1/4+1,N2/4+1,N3/4+1), t23%var(N1/4+1,N2/4+1,N3/4+1),&
                             t33%var(N1/4+1,N2/4+1,N3/4+1), T  %var(N1/4+1,N2/4+1,N3/4+1)

      write(32,*) it, it*dt, u1 %var(N1/4+1,N2/2+1,N3/4+1), u2 %var(N1/4+1,N2/2+1,N3/4+1),&
                             u3 %var(N1/4+1,N2/2+1,N3/4+1), t11%var(N1/4+1,N2/2+1,N3/4+1),&
                             t12%var(N1/4+1,N2/2+1,N3/4+1), t13%var(N1/4+1,N2/2+1,N3/4+1),&
                             t22%var(N1/4+1,N2/2+1,N3/4+1), t23%var(N1/4+1,N2/2+1,N3/4+1),&
                             t33%var(N1/4+1,N2/2+1,N3/4+1), T  %var(N1/4+1,N2/2+1,N3/4+1)

      write(33,*) it, it*dt, u1 %var(N1/4+1,N2/4+1,N3/2+1), u2 %var(N1/4+1,N2/4+1,N3/2+1),&
                             u3 %var(N1/4+1,N2/4+1,N3/2+1), t11%var(N1/4+1,N2/4+1,N3/2+1),&
                             t12%var(N1/4+1,N2/4+1,N3/2+1), t13%var(N1/4+1,N2/4+1,N3/2+1),&
                             t22%var(N1/4+1,N2/4+1,N3/2+1), t23%var(N1/4+1,N2/4+1,N3/2+1),&
                             t33%var(N1/4+1,N2/4+1,N3/2+1), T  %var(N1/4+1,N2/4+1,N3/2+1)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

#if( DIMENSION_GEO == 2 && VISCOELASTIC_MODELS > 0)

    if (pas_temps_adptatif == 0) then 
       if (it > 1)  write(23,*) it, it*dt, x1_x2_min_t11_x, x1_x2_min_t11_y, minval(t11_transvaser),&       
                                           x1_x2_max_t11_x, x1_x2_max_t11_y, maxval(t11_transvaser),&       
                                           x1_x2_min_t22_x, x1_x2_min_t22_y, minval(t22_transvaser),&       
                                           x1_x2_max_t22_x, x1_x2_max_t22_y, maxval(t22_transvaser)
            
    else if (pas_temps_adptatif == 1) then        
       temps_cumule = temps_cumule + dt 
       if (it > 1)  write(23,*) it, it*dt, x1_x2_min_t11_x, x1_x2_min_t11_y, minval(t11_transvaser),&       
                                           x1_x2_max_t11_x, x1_x2_max_t11_y, maxval(t11_transvaser),&        
                                           x1_x2_min_t22_x, x1_x2_min_t22_y, minval(t22_transvaser),&        
                                           x1_x2_max_t22_x, x1_x2_max_t22_y, maxval(t22_transvaser)
    end if


#endif


  end subroutine convergence_stationnaire















  subroutine niveau_de_convergence(var)
    implicit none
    integer:: i, j, k
    type(variable),intent(inout):: var
    real(nk) :: somme, RH

#if (DIMENSION_GEO == 2)
    if (var%indice == p) then
       somme = 0.0
       do j=1, var%NJ
          do i=1, var%NI                             

             RH = abs(                           &
                  var%ax1(i,j)*var%var0(i-1,j)  +&
                  var%cx1(i,j)*var%var0(i+1,j)  +&
                  var%ax2(i,j)*var%var0(i,j-1)  +&
                  var%cx2(i,j)*var%var0(i,j+1)  +&
                  var%b  (i,j)*var%var0(i,j)    -&
                  var%d  (i,j)                  )
             somme = somme + RH*RH
          end do
       end do
    else if (var%indice /= p) then  
       somme = 0.0

       do j=1, var%NJ
          do i=1, var%NI
          
             RH=  abs(                          &
                  var%ax1(i,j)*var%var(i-1,j)   +&
                  var%cx1(i,j)*var%var(i+1,j)   +&
                  var%ax2(i,j)*var%var(i,j-1)   +&
                  var%cx2(i,j)*var%var(i,j+1)   +&
                  var%b  (i,j)*var%var(i,j)     -&
                  var%d  (i,j)                   )                  
             somme = somme + RH*RH
          end do
       end do
    end if
    var%niveau_conver =  sqrt(somme) !/real(var%N1*var%N2, nk)


#elif (DIMENSION_GEO == 3)

    if (var%indice == p) then
       somme = 0.0
       do k=1, var%NK
          do j=1, var%NJ
             do i=1, var%NI        
                RH = abs(                                 &
                     var%ax1(i,j,k)*var%var0(i-1,j,k)    +&
                     var%cx1(i,j,k)*var%var0(i+1,j,k)    +&
                     var%ax2(i,j,k)*var%var0(i,j-1,k)    +&
                     var%cx2(i,j,k)*var%var0(i,j+1,k)    +&
                     var%ax3(i,j,k)*var%var0(i,j,k-1)    +&
                     var%cx3(i,j,k)*var%var0(i,j,k+1)    +&
                     var%b  (i,j,k)*var%var0(i,j,k)      -&
                     var%d  (i,j,k)                      )
                somme = somme + RH*RH
             end do
          end do
       end do


    else if (var%indice /= p) then 

       somme = 0.0
       do k=1,var%NK
          do j=1, var%NJ
             do i=1, var%NI
                RH = abs(                                 &
                     var%ax1(i,j,k)*var%var(i-1,j,k)     +&
                     var%cx1(i,j,k)*var%var(i+1,j,k)     +&
                     var%ax2(i,j,k)*var%var(i,j-1,k)     +&
                     var%cx2(i,j,k)*var%var(i,j+1,k)     +&
                     var%ax3(i,j,k)*var%var(i,j,k-1)     +&
                     var%cx3(i,j,k)*var%var(i,j,k+1)     +&
                     var%b  (i,j,k)*var%var(i,j,k)       -&
                     var%d  (i,j,k)                      )
                somme = somme + RH*RH 
             end do
          end do
       end do
    end if
    var%niveau_conver = sqrt(somme)  !!  /real(var%N1*var%N2, nk)
#endif

  end subroutine niveau_de_convergence









#if (DIMENSION_GEO == 2)
  subroutine test_div_U_2D   (u1, u2 )
    implicit none   
    integer:: i, j
    real(nk),dimension(0:N1+1,0:N2+1),intent(in)      :: u1, u2  
    real(nk),dimension(1:N1-1,1:N2-1)                :: div_U_centres   
    real(nk),dimension(1:N1,1:N2)                    :: div_U_coins    
    real(nk) :: somme




!!!!!!--->  div_U_sur_coins_avec_coins
    if (it== pas_write ) open (82,file="div_U_coins")         

    somme =0.0
    div_U_coins = 0.0     
    do j= 1, N2
       do i=1, N1                   
          div_U_coins(i,j) = CENTRE_2(1._nk, u1, i, j, 1) + CENTRE_2(1._nk, u2, i, j, 2)        
          if (j== 1 .or. j == N2 .or. i == 1 .or. i == N1)  div_U_coins(i,j) = 0
          if (it== pas_write ) then 
             write(82,*) div_U_coins(i,j)                    
          end if
          somme = somme + div_U_coins(i,j)**2.                              
       end do
       if (it== pas_write )  write(82,*)"  "
    end do
    if (it== pas_write ) close(82) 
    div_U_sur_coins_avec_coins = sqrt(somme)


!!!!!!--->  div_U_sur_centres_avec_Uinter
    if (it== pas_write )  open (82,file="div_U_centres")  

    somme =0.0    
    div_U_centres = 0.0 
    do j= 1, N2-1
       do i=1, N1-1

          div_U_centres(i,j) = ( var_intf_u1(i+1,j) - var_intf_u1(i,j) )/dx1_R(i,j) + &
                               ( var_intf_u2(i,j+1) - var_intf_u2(i,j) )/dx2_R(i,j)

          if (it== pas_write ) then 
             write(82,*) div_U_centres(i,j)                    
          end if

          somme = somme + div_U_centres(i,j)**2
       end do
       if (it== pas_write ) write(82,*)"  "
    end do
    if (it== pas_write ) close(82)      
    div_U_sur_centres_avec_Uinter = sqrt(somme)


   end subroutine test_div_U_2D
#endif







#if (DIMENSION_GEO == 3)
  subroutine test_div_U_3D   (u1, u2, u3)
   implicit none  
    integer:: i, j, k
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(in)   :: u1, u2, u3  
    real(nk),dimension(1:N1-1,1:N2-1,1:N3-1)              :: div_U_centres   
    real(nk),dimension(1:N1,1:N2,1:N3)                    :: div_U_coins    
    real(nk) :: somme
    ! div centres avec coins 
    somme =0.0
    div_U_centres = 0.0  
    do k= 1,N3-1
       do j= 1,N2-1
          do i=1, N1-1
             div_U_centres(i,j,k) =                                                               &
                 (CENTRE_2(1._nk, u1, i  ,j  ,1  ,k  ) + CENTRE_2(1._nk, u1, i  ,j+1,1  ,k  )     &
                 +CENTRE_2(1._nk, u1, i  ,j  ,1  ,k+1) + CENTRE_2(1._nk, u1, i  ,j+1,1  ,k+1))/4. &
               + (CENTRE_2(1._nk, u1, i  ,j  ,2  ,k  ) + CENTRE_2(1._nk, u1, i+1,j  ,2  ,k  )     &
                 +CENTRE_2(1._nk, u1, i  ,j  ,2  ,k+1) + CENTRE_2(1._nk, u1, i+1,j  ,2  ,k+1))/4. &
               + (CENTRE_2(1._nk, u1, i  ,j  ,3  ,k  ) + CENTRE_2(1._nk, u1, i+1,j  ,3  ,k  )     &
                 +CENTRE_2(1._nk, u1, i  ,j+1,3  ,k  ) + CENTRE_2(1._nk, u1, i+1,j+1,3  ,k  ))/4.

             somme = somme +  div_U_centres(i,j,k)*div_U_centres(i,j,k)
          end do
       end do
    end do
    div_U_sur_centres_avec_coins = sqrt(somme) 



    !div coins avec coins 
    somme =0.0
    div_U_coins = 0.0  
    do k= 1,N3
       do j= 1,N2
          do i=1, N1

             div_U_coins(i,j,k) = CENTRE_2( 1._nk, u1, i, j, 1, k) + CENTRE_2( 1._nk, u1, i, j, 2, k) &             
                                + CENTRE_2( 1._nk, u1, i, j, 3, k)

             if ( i== 1 .or. i == N1   .or. j== 1 .or. j == N2   .or. k== 1 .or. k == N3   ) div_U_coins(i,j,k) = 0.
             if ( i== 2 .or. i == N1-1 .or. j== 2 .or. j == N2-1 .or. k== 2 .or. k == N3-1 ) div_U_coins(i,j,k) = 0.

             somme = somme +  div_U_coins(i,j,k)*div_U_coins(i,j,k)
          end do
       end do
    end do
    div_U_sur_coins_avec_coins = sqrt(somme) 



    ! div centres avec interfaces
    somme =0.0
    div_U_centres = 0.0  
    do k= 1, N3-1
       do j= 1, N2-1
          do i=1, N1-1

             div_U_centres(i,j,k) =                                    & 
                  ( var_intf_u1(i+1,j,k) - var_intf_u1(i,j,k)  )/dx1_R(i,j,k) + &
                  ( var_intf_u2(i,j+1,k) - var_intf_u2(i,j,k)  )/dx2_R(i,j,k) + &
                  ( var_intf_u3(i,j,k+1) - var_intf_u3(i,j,k)  )/dx3_R(i,j,k)    

             somme = somme +  div_U_centres(i,j,k)*div_U_centres(i,j,k) 
          end do
       end do
    end do
    div_U_sur_centres_avec_Uinter = sqrt(somme) 


  end subroutine test_div_U_3D
#endif











  subroutine reperage_point_choisi (point_choisi_espace, pce_x, pce_y, pce_k)
    implicit none
    integer :: i, j, k
    integer,intent(in   ):: point_choisi_espace
    integer,intent(inout):: pce_x, pce_y, pce_k

#if (DIMENSION_GEO == 2)     
    do j=1,N2
       do i=1,N1
          ll = i + (j-1)*N1
          if ( point_choisi_espace == ll ) then
             pce_x = i
             pce_Y = j 
          end if
       end do
    end do
#elif (DIMENSION_GEO == 3)          
    do K=1, N3
       do j=1,N2
          do i=1,N1
            ll = i + (j-1)*N1 + (k-1)*N1*N2
            if ( point_choisi_espace == ll ) then
               pce_x = i
               pce_Y = j 
               pce_k = k
             end if
          end do
       end do
    end do
#endif
  end subroutine reperage_point_choisi

















  subroutine calcul_erreur_espace_ponctuelle(Er, var, kkk, point_choisi_espace, cord_x1_0, cord_x2_0, &
       approx_O2, approx_O4, cord_x3_0)
    implicit none

    integer :: i, j, l, k, s
    integer :: pce_x, pce_y, pce_k
    type(Erreur_var)      :: Er
    type(vecteur_variable):: var
    integer,intent(in   ) :: kkk, point_choisi_espace
    integer,intent(inout) :: cord_x1_0, cord_x2_0
    integer,intent(inout),optional ::   cord_x3_0

    integer,intent(in   ) :: approx_O2, approx_O4
    real(nk),dimension(p,4)::approx
    integer:: cord=0, cord_x1x2, cord_x1x2x3, nbr_surface, nbr_ligne, nbr_points
    approx =0.



!!!!----------------------------------------------------------------------------------------------------
    if (kkk==1) then 
       nbr_points   = N1
       nbr_ligne    = 1
       nbr_surface  = 1
       i            = 1
       do while (nbr_points < point_choisi_espace ) 
          nbr_points  = nbr_points + N1 
          nbr_ligne   = nbr_ligne  +  1 

       end do
    end if
!!!!----------------------------------------------------------------------------------------------------
    if (kkk==1) then
#if (DIMENSION_GEO == 2)
       cord_x1_0 =  point_choisi_espace - (nbr_ligne-1)*N1
#elif (DIMENSION_GEO == 3)
       cord_x1_0 =  point_choisi_espace - (nbr_ligne-1)*N1 - (nbr_surface -1)*N1*N2 
       cord_x3_0 =  nbr_surface       
#endif 
       cord_x2_0 =  nbr_ligne
    elseif (kkk==2)  then
       cord_x1_0 = cord_x1_0*2 - 1   
       cord_x2_0 = cord_x2_0*2 - 1  
#if (DIMENSION_GEO == 3)
       cord_x3_0 = cord_x3_0*2 - 1
#endif           
    else
       cord_x1_0 = cord_x1_0*2  -1
       cord_x2_0 = cord_x2_0*2  -1
#if (DIMENSION_GEO == 3)
       cord_x3_0 = cord_x3_0*2  -1
#endif  
    end if


#if (DIMENSION_GEO == 2)
    cord_x1x2 = (cord_x2_0 -1)*N1 + cord_x1_0
    call reperage_point_choisi(cord_x1x2, pce_x, pce_y, pce_k)
    if (kkk == 1 ) print*, point_choisi_espace, cord_x1x2, N1
    if (kkk /= 1 ) print*, point_choisi_espace, cord_x1x2, N1

#elif (DIMENSION_GEO == 3)
    cord_x1x2x3  = (cord_x3_0 -1)*N1*N2 + (cord_x2_0 -1)*N1 + cord_x1_0
    call reperage_point_choisi(cord_x1x2x3, pce_x, pce_y, pce_k)
#endif   

!!!!----------------------------------------------------------------------------------------------------

#if (DIMENSION_GEO == 2)
    do i=1, p-1
       if (kkk==1) then      
          Er%var(i)%val_var_espace(kkk) = var%W(i)%var(pce_x, pce_y)
       else
          if (approx_O2 == 1) then 
             Er%var(i)%val_var_espace(kkk) = (                            &
                  var%W(i)%var(pce_x, pce_y) +&
                  var%W(i)%var(pce_x+1, pce_y) +&
                  var%W(i)%var(pce_x, pce_y+1) +&
                  var%W(i)%var(pce_x+1, pce_y+1) )/4._nk 
          elseif (approx_O4 == 1) then
             cord = i+(j-1)*N1 -1 - N1
             do j=1, 4         
                approx(i,j) = (                  &
                     -var%W(i)%var(pce_x-1, pce_y-1) +&
                     9._nk*var%W(i)%var(pce_x, pce_y-1) +&
                     9._nk*var%W(i)%var(pce_x+1, pce_y-1) -&
                     var%W(i)%var(pce_x+2, pce_y-1)  )/16._nk
                cord = cord + N1+2
             end do
             Er%var(i)%val_var_espace(kkk) = (-&
                  approx(i,1)           +&
                  9._nk*approx(i,2)           +&
                  9._nk*approx(i,3)           -&
                  approx(i,4)            )/16._nk             
          else
             Er%var(i)%val_var_espace(kkk) = var%W(i)%var(pce_x, pce_y) 
          end if
       endif
    end do

    print*, (Er%var(i)%val_var_espace(kkk),i=1,1), 'ttttttttttt'
!    grid_size(kkk) = dx2
    call calc_Erreur_espace_ponctuelle (Er,kkk)






    do J=1,N2, pas
       do i=1,N1, pas  !!!! SANS PAS POUR LE POISEUILLE AVEC N1 FIXE
          write(70,*) i,j, (var%W(s)%var(i,j), s=1,p-1)
       end do
    end do
    write(70,*)
    pas = pas * 2


#elif (DIMENSION_GEO == 3)



#endif     

  end subroutine calcul_erreur_espace_ponctuelle





  subroutine calc_Erreur_espace_ponctuelle (Er, kkk)
    implicit none
    integer :: i, j, k, l, ll
    type(Erreur_var):: Er
    integer ,intent(in) :: kkk


    print"(a, i4,3x,a,(E13.6) )",' kkk    dx            var   ordre      N1=',N1,'dt=', dt 
    if (kkk > 2) then
       do i=1,p
          Er%var(i)%ord_espace(kkk)= (Er%var(i)%val_var_espace(kkk-2)-Er%var(i)%val_var_espace(kkk-1) ) /&
               ( Er%var(i)%val_var_espace(kkk-1)-Er%var(i)%val_var_espace(kkk) ) /2.
          print"(i4,(E16.6), i4, 7(f20.15))", kkk, i ,Er%var(i)%ord_espace(kkk) , &
               ( Er%var(i)%val_var_espace(kkk-1)-Er%var(i)%val_var_espace(kkk) )
       end do
       print*,'   '
    end if

  end subroutine calc_Erreur_espace_ponctuelle








  subroutine calcul_erreur_temps_ponctuelle(Er, var, kk, pce_x, pce_y, pce_k )
    implicit none
    integer :: i, j, l, k, s
    integer :: pce_x, pce_y, pce_k
    type(Erreur_var):: Er
    type(vecteur_variable):: var  
    integer ,intent(in) ::  kk


#if (DIMENSION_GEO == 2)
    do i=1,p
       Er%var(i)%val_var_temps(kk) = var%W(i)%var(pce_x, pce_y)
       print*,i, kk, Er%var(i)%val_var_temps(kk)

    end do
#elif (DIMENSION_GEO ==3)
    do i=1,p
       Er%var(i)%val_var_temps(kk) = var%W(i)%var(pce_x, pce_y, pce_k)
       print*,i, kk, Er%var(i)%val_var_temps(kk)

    end do
#endif
    !Er%var(12)%val_var_temps(kk) = T_var (int(M-10)*(N+2*NG) + int(real(N)/3) )
    time_size(kk) = dt


    print*,' kk    dt            var   ordre' 
    if (kk > 2) then 
       do i=1,p
          Er%var(i)%ord_temps(kk)= (Er%var(i)%val_var_temps(kk-2)-Er%var(i)%val_var_temps(kk-1) ) /&
               ( Er%var(i)%val_var_temps(kk-1)-Er%var(i)%val_var_temps(kk) ) /2.
          print"(i4,(E16.6), i4, 7(f12.6))", kk, dt, i ,Er%var(i)%ord_temps(kk) 
       end do
       print*, ' '
    end if


!!!--------------------------------------------------------------------
#if (DIMENSION_GEO == 2)
    do J=1,N2
       do i=1,N1
          write(69,*) i,j, (var%W(s)%var(i,j), s=1,p-1)
       end do
    end do
    write(69,*)
    do J=1,N2-1
       do i=1,N1-1
          write(68,*) i,j,var%W(p)%var(i,j)
       end do
    end do
    write(68,*)
#elif (DIMENSION_GEO == 3)
    do k=1, N3
       do J=1,N2
          do i=1,N1
             write(69,*) i,j,k, (var%W(s)%var(i,j,k), s=1,p-1) 
          end do
       end do
    end do
    write(69,*)
    do k=1, N3-1
       do J=1,N2-1
          do i=1,N1-1
             write(68,*) i,j,k, var%W(p)%var(i,j,k)
          end do
       end do
    end do
    write(68,*)
#endif     
  end subroutine calcul_erreur_temps_ponctuelle











  subroutine openfile
    implicit none        

    open (30,file="monitor_P1")
    open (31,file="monitor_P2")
    open (32,file="monitor_P3")
    open (33,file="monitor_P4")

    open (34,file="Total_Nu_time_series")
#if(VISCOELASTIC_MODELS > 0)
    open (21,file="viscoelastic_parts")
    open (22,file="kinetic_energy_budget")
    open (23,file="min_t11_t22")
    
#endif
  end subroutine openfile


  subroutine closefile
    implicit none
    close(30)
    close(31)
    close(32)
    close(33)

    close(34)
#if(VISCOELASTIC_MODELS > 0)
    close(21)
    close(22)
    close(23)
#endif

  end subroutine closefile













  subroutine erreur_en_champs(Er_norme, erreur_en_temps, erreur_en_espace, RICH, SOL_ANA)
    implicit none

    integer :: g, l,o, s, ii , jj, kk, iii, ll , ppp
    type(Erreur_var_norme):: Er_norme 
    integer,intent(in):: erreur_en_temps, erreur_en_espace, RICH, SOL_ANA
    real(nk),dimension(var%W(1)%NNN,2,p-1):: variable   
    real(nk),dimension(var%W(P)%NNN,2    ):: variable_P 
    integer :: N1, N2, N1_P, N2_P, pas_choisi_temps, indice,  nN1, nN2, nN3, nN123
    real(nk):: dt, dx1, dx2  , indice1    , indice2  , indice3
    character(len=15) :: schema, variable_utau, variable_u123 
    integer :: N12, N22, N32
#if (DIMENSION_GEO == 3)
    integer :: N3, N3_P
#endif 


    !real(nk),dimension( var%W(1)%N2*var%W(1)%N3 , 2, 4):: u1tau   
    real(nk),dimension( 81*81 , 2, 4):: u1tau   
    real(nk),dimension( 81*81 ,    4):: diff

    real(nk),dimension(:,:,:),allocatable:: u123
    real(nk),dimension(:,:,:),allocatable:: diffu123
    real(nk),dimension(:,:  ),allocatable:: u123_2
    real(nk),dimension(:,:,:),allocatable:: Erreur_u123, ordre_u123




      variable =0._nk ; variable_P=0._nk 


    if (erreur_en_espace==1.and.RICH==1) then
       print*, '  '
       print*, '  '
       print*, 'Erreur spatialle en norme avec RICH'

       do ppp=1, 4

          if ( ppp == 1 ) then 
             print*,  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             print*,  '!!!!!!!!!!!!!!!!!! VARIBLES !!!!!!!!!!!!!!!!!!!!!!'
             print*,  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             open(11,file="champs_UT_en_fonction_dx",status="old")
          elseif ( ppp == 2 ) then 
             print*,  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             print*,  '!!!!!!!!!!!!!!!!!! Partie Hyperbolique !!!!!!!!!!!!!!!!!!!!!!'
             print*,  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             open(11,file="champs_Partie_Hyperbolique_en_fonction_dx",status="old")
          elseif ( ppp == 3 )  then 
             print*,  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             print*,  '!!!!!!!!!!!!!!!!!! Partie Parabolique  !!!!!!!!!!!!!!!!!!!!!!'
             print*,  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             open(11,file="champs_Partie_Parabolique_en_fonction_dx",status="old")
          elseif ( ppp == 4 )  then 
             print*,  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             print*,  '!!!!!!!!!!!!!!!!!! Gradient Pression !!!!!!!!!!!!!!!!!!!!!!'
             print*,  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             open(11,file="champs_Gradient_Pression_en_fonction_dx",status="old")
          end if


          !open(unit=10,file="champ_P_en_fonction_dt",status="old")
#if(DIMENSION_GEO == 2)
          read(11,*) N1, N2, pas_choisi_temps, dx1, dx2
          !read(10,*) N1_P, N2_P, pas_choisi_temps, dt
          print*, ' kkk   norme   u1              u2              t11             t12             t22',&
               '             T               P' 
#elif(DIMENSION_GEO == 3)
          read(11,*) N1, N2, N3, pas_choisi_temps, dt
          !read(10,*) N1_P, N2_P, N3_P, pas_choisi_temps, dt
          print*, ' kkk   norme   u1              u2              u3              t11             t12',&
               '                  t12             t22             t23             t33'              ,&
               '                   T               P'
#endif 

          do j=1,nbr_it_espace

#if(DIMENSION_GEO == 2)
             do i=1, N1*N2
                if ( ppp == 1 .or. ppp == 4 ) then
                   read(11,*) indice, ( variable(i,1,k), k=1,p-1)
                else if ( ppp == 2 ) then 
                   read(11,*) indice, ( variable(i,1,k), k=1,5)
                else  
                   read(11,*) indice, ( variable(i,1,k), k=1,2)
                endif

             end do
             read(11,*)             

#elif(DIMENSION_GEO == 3)   
             do i=1, N1*N2*N3
                if ( ppp == 1 .or. ppp == 4 ) then
                   read(11,*) indice, ( variable(i,1,k), k=1,p-1)                   
                else if ( ppp == 2 ) then  
                   read(11,*) indice, ( variable(i,1,k), k=1,9)
                else
                   read(11,*) indice, ( variable(i,1,k), k=1,3)
                end if
             end do
             read(11,*)
#endif



             if(j==1) then
                !print"(i4)",j
             else    

                do i=1, p-1
                   Er_norme%norme(1)%var(i)%val_var_espace(j)=     sum(abs(variable(:,1,i)-variable(:,2,i))   )               
                   Er_norme%norme(2)%var(i)%val_var_espace(j)=sqrt(sum(abs(variable(:,1,i)-variable(:,2,i))**2)) !&
                   ! /sum(variable(:,2,i)**2))  
                   Er_norme%norme(3)%var(i)%val_var_espace(j)=  maxval(abs(variable(:,1,i)-variable(:,2,i))   )  
                end do



#if(DIMENSION_GEO == 2)
                open (82,file="Error_espace_1")
      
                       write(82,*) Er_norme%norme(1)%var(1)%val_var_espace(1), Er_norme%norme(1)%var(2)%val_var_espace(1), &
                                   Er_norme%norme(1)%var(3)%val_var_espace(1), Er_norme%norme(1)%var(4)%val_var_espace(1), &
                                   Er_norme%norme(1)%var(5)%val_var_espace(1), Er_norme%norme(1)%var(6)%val_var_espace(1), &
                                   Er_norme%norme(1)%var(7)%val_var_espace(1)

                       write(82,*)"  "
!                          write(83,*)"  "         

                close(82)
#elif(DIMENSION_GEO == 3)   

#endif

#if(DIMENSION_GEO == 2)
                do k=1,3
                   print"(i4, 3x, a, i1, 2x, 7(E16.6))",j,'L',k, (Er_norme%norme(k)%var(i)%val_var_espace(j),i=1,p) 
                end do
#elif(DIMENSION_GEO == 3)   
                do k=1,3
                   print"(i4, 3x, a, i1, 2x, 11(E16.6))",j,'L',k, (Er_norme%norme(k)%var(i)%val_var_espace(j),i=1,p) 
                end do
#endif

                print*,  ' '         
             end if
             variable  (:,2,:) = variable  (:,1,:) 
             dx1 = dx1/2.0_nk
          end do

          do j=2, nbr_it_espace-1          
             do i=1, p
                do k=1,3           
                   Er_norme%norme(K)%var(i)%ord_espace(j) = &
                        log(Er_norme%norme(K)%var(i)%val_var_espace(j)/Er_norme%norme(K)%var(i)%val_var_espace(j+1))/log(2._nk)
                end do
             end do
          end do

          print*, '  '
          print*,'ordre de l erreur en espace en norme'

#if(DIMENSION_GEO == 2)
          print*, ' kk   norme   u1          u2          t11         t12         t22         T           P'
          do j=2, nbr_it_espace-1    
             do k=1,3
                print"(i4, 3x, a, i1, 2x, 7(f12.6))",j,'L',k, (Er_norme%norme(k)%var(i)%ord_espace(j),i=1,p) 
             end do
             print*, '  '
          end do
#elif(DIMENSION_GEO == 3)
          print*,' kk   norme   u1          u2          u3          t11         t12         t12         t22',&
               '         t23         t33         T           P'
          do j=2, nbr_it_espace-1          
             do k=1,3
                print"(i4, 3x, a, i1, 2x, 11(f12.6))",j,'L',k, (Er_norme%norme(k)%var(i)%ord_espace(j),i=1,p) 
             end do
             print*, '  ' 
          end do
#endif



!!!! pour rapport latex       
          do i=1, p                    
#if(DIMENSION_GEO == 2)
             if (i==1) print*, '$u_1$'
             if (i==2) print*, '$u_2$'
             if (i==3) print*, '$t_{11}$'
             if (i==4) print*, '$t_{12}$'
             if (i==5) print*, '$t_{22}$'
             if (i==6) print*, '$T$'
             if (i==7) print*, '$P$'
#elif(DIMENSION_GEO == 3)
             if (i==1) print*, '$u_1$'
             if (i==2) print*, '$u_2$'
             if (i==3) print*, '$u_3$'
             if (i==4) print*, '$t_{11}$'
             if (i==5) print*, '$t_{12}$'
             if (i==6) print*, '$t_{13}$'
             if (i==7) print*, '$t_{22}$'
             if (i==8) print*, '$t_{23}$'
             if (i==9) print*, '$t_{33}$'
             if (i==10) print*, '$T$'
             if (i==11) print*, '$P$'
#endif           
             do j=1,nbr_it_espace
                if (j==1 ) then 

                   print "(3(4a, a, 4a, a), a)", &
                        ('  ','&' ,'  ','&','     -      ', &
                        '  ','&' ,'  ','&','     -      ', k=1,3 ), '  \\'

                else if ( j==2 ) then 

                   print "(3(4a, E12.4, 4a, a), a)", &
                        ('  ','&' ,'  ','&', Er_norme%norme(k)%var(i)%val_var_espace(j), &
                        '  ','&' ,'  ','&','     -      ', k=1,3 ), '  \\' 

                else
                   print "(3(4a, E12.4, 4a, f12.4), a)", &
                        ('  ','&' ,'  ','&', Er_norme%norme(k)%var(i)%val_var_espace(j), &
                        '  ','&' ,'  ','&', Er_norme%norme(k)%var(i)%ord_espace(j-1), k=1,3 ), '  \\' 
                end if
             end do
          end do

       end do

    end if






    if (erreur_en_temps==1.and.RICH==1) then
       print*, '  '
       print*, '  '
       print*, 'Erreur temporelle en norme avec RICH'

       open(11,file="champs_UT_en_fonction_dt",status="old")
       open(unit=10,file="champ_P_en_fonction_dt",status="old")
#if(DIMENSION_GEO == 2)
       read(11,*) N1, N2, pas_choisi_temps, dt
       read(10,*) N1_P, N2_P, pas_choisi_temps, dt
       print*, ' kk   norme   u1              u2              t11             t12             t22',&
            '             T               P' 
#elif(DIMENSION_GEO == 3)
       read(11,*) N1, N2, N3, pas_choisi_temps, dt
       read(10,*) N1_P, N2_P, N3_P, pas_choisi_temps, dt
       print*, ' kk   norme   u1              u2              u3              t11             t12',&
            '             t12             t22             t23             t33'              ,&
            '             T               P'
#endif 

       do j=1,nbr_it_temps

          do i=1, var%W(1)%NNN
             read(11,*) indice, ( variable(i,1,k), k=1,p-1)
          end do
          read(11,*)

          do i=1, var%W(P)%NNN
             read(10,*) indice, variable_P(i,1)
          end do
          read(10,*)

          if(j==1) then
             !print"(i4)",j
          else           
             do i=1, p-1
                Er_norme%norme(1)%var(i)%val_var_temps(j)=     sum(abs(variable(:,1,i)-variable(:,2,i))   )               
                Er_norme%norme(2)%var(i)%val_var_temps(j)=sqrt(sum(abs(variable(:,1,i)-variable(:,2,i))**2))
                !/sum(variable(:,2,i)**2))  
                Er_norme%norme(3)%var(i)%val_var_temps(j)=  maxval(abs(variable(:,1,i)-variable(:,2,i))   )  
             end do
             Er_norme%norme(1)%var(P)%val_var_temps(j)  =     sum(abs(variable_P(:,1)-variable_P(:,2))   )               
             Er_norme%norme(2)%var(P)%val_var_temps(j)  =sqrt(sum(abs(variable_P(:,1)-variable_P(:,2))**2))
             !/sum(variable_P(:,2)**2))  
             Er_norme%norme(3)%var(P)%val_var_temps(j)  =  maxval(abs(variable_P(:,1)-variable_P(:,2))   )

#if(DIMENSION_GEO == 2)
             do k=1,3
                print"(i4, 3x, a, i1, 2x, 7(E16.6))",j,'L',k, (Er_norme%norme(k)%var(i)%val_var_temps(j),i=1,p) 
             end do
#elif(DIMENSION_GEO == 3)   
             do k=1,3
                print"(i4, 3x, a, i1, 2x, 11(E16.6))",j,'L',k, (Er_norme%norme(k)%var(i)%val_var_temps(j),i=1,p) 
             end do
#endif
             print*,  ' '         
          end if
          variable  (:,2,:) = variable  (:,1,:) 
          variable_P(:,2  ) = variable_P(:,1  )
          dt = dt/2.0_nk
       end do

       do j=2, nbr_it_temps-1          
          do i=1, p
             do k=1,3
                Er_norme%norme(K)%var(i)%ord_temps(j) = &
                     log(Er_norme%norme(K)%var(i)%val_var_temps(j)/Er_norme%norme(K)%var(i)%val_var_temps(j+1))/log(2._nk)
             end do
          end do
       end do

       print*, '  '
       print*,'ordre de l erreur temporelle en norme'

#if(DIMENSION_GEO == 2)
       print*, ' kk   norme   u1          u2          t11         t12         t22         T           P'
       do j=2, nbr_it_temps-1    
          do k=1,3
             print"(i4, 3x, a, i1, 2x, 7(f12.6))",j,'L',k, (Er_norme%norme(k)%var(i)%ord_temps(j),i=1,p) 
          end do
          print*, '  '
       end do
#elif(DIMENSION_GEO == 3)
       print*,' kk   norme   u1          u2          u3          t11         t12         t12         t22',&
            '         t23         t33         T           P'
       do j=2, nbr_it_temps-1          
          do k=1,3
             print"(i4, 3x, a, i1, 2x, 11(f12.6))",j,'L',k, (Er_norme%norme(k)%var(i)%ord_temps(j),i=1,p) 
          end do
          print*, '  ' 
       end do
#endif


!!!! pour rapport latex       
       do i=1, p                    
#if(DIMENSION_GEO == 2)
          if (i==1) print*, '$u_1$'
          if (i==2) print*, '$u_2$'
          if (i==3) print*, '$t_{11}$'
          if (i==4) print*, '$t_{12}$'
          if (i==5) print*, '$t_{22}$'
          if (i==6) print*, '$T$'
          if (i==7) print*, '$P$'
#elif(DIMENSION_GEO == 3)
          if (i==1) print*, '$u_1$'
          if (i==2) print*, '$u_2$'
          if (i==3) print*, '$u_3$'
          if (i==4) print*, '$t_{11}$'
          if (i==5) print*, '$t_{12}$'
          if (i==6) print*, '$t_{13}$'
          if (i==7) print*, '$t_{22}$'
          if (i==8) print*, '$t_{23}$'
          if (i==9) print*, '$t_{33}$'
          if (i==10) print*, '$T$'
          if (i==11) print*, '$P$'
#endif           
          do j=1,nbr_it_temps
             if (j==1 ) then 

                print "(3(4a, a, 4a, a), a)", &
                     ('  ','&' ,'  ','&','     -      ', &
                     '  ','&' ,'  ','&','     -      ', k=1,3 ), '  \\'

             else if ( j==2 ) then 

                print "(3(4a, E12.4, 4a, a), a)", &
                     ('  ','&' ,'  ','&', Er_norme%norme(k)%var(i)%val_var_temps(j), &
                     '  ','&' ,'  ','&','     -      ', k=1,3 ), '  \\' 

             else
                print "(3(4a, E12.4, 4a, f12.4), a)", &
                     ('  ','&' ,'  ','&', Er_norme%norme(k)%var(i)%val_var_temps(j), &
                     '  ','&' ,'  ','&', Er_norme%norme(k)%var(i)%ord_temps(j-1), k=1,3 ), '  \\' 
             end if
          end do
       end do

    end if






    if (erreur_en_temps==1 .and. SOL_ANA==1 .and. T_G_T == 1 ) then
       print*, '  '
       print*, '  '
       print*, 'Erreur temporelle en norme avec la solution analytique'

       open(unit=11,file="champs_UT_en_fonction_dt",status="old")
       open(unit=10,file="champ_P_en_fonction_dt"  ,status="old")
#if(DIMENSION_GEO == 2)
       read(11,*) N1, N2, pas_choisi_temps, dt
       read(10,*) N1_P, N2_P, pas_choisi_temps, dt
       print*, ' kk   norme   u1              u2              t11             t12             t22',&
            '             T               P' 
#elif(DIMENSION_GEO == 3)
       read(11,*) N1, N2, N3, pas_choisi_temps, dt
       read(10,*) N1_P, N2_P, N3_P, pas_choisi_temps, dt
       print*, ' kk   norme   u1              u2              u3              t11             t12',&
            '             t12             t22             t23             t33'              ,&
            '             T               P'
#endif 

       do j=1,nbr_it_temps

          do i=1, var%W(1)%NNN
             read(11,*) indice, (variable(i,1,k), k=1,p-1)
          end do
          read(11,*)

          do i=1, var%W(P)%NNN
             read(10,*) indice, variable_P(i,1)
          end do
          read(10,*)

          if(j==1) then
             !print"(i4)",j
          else           
             do i=1, p-1
                Er_norme%norme(1)%var(i)%val_var_temps(j)=     sum(abs(variable(:,1,i)-variable(:,2,i))   )               
                Er_norme%norme(2)%var(i)%val_var_temps(j)=sqrt(sum(abs(variable(:,1,i)-variable(:,2,i))**2)/sum(variable(:,2,i)**2))  
                Er_norme%norme(3)%var(i)%val_var_temps(j)=  maxval(abs(variable(:,1,i)-variable(:,2,i))   )  
             end do
             Er_norme%norme(1)%var(P)%val_var_temps(j)  =     sum(abs(variable_P(:,1)-variable_P(:,2))   )               
             Er_norme%norme(2)%var(P)%val_var_temps(j)  =sqrt(sum(abs(variable_P(:,1)-variable_P(:,2))**2)/sum(variable_P(:,2)**2))  
             Er_norme%norme(3)%var(P)%val_var_temps(j)  =  maxval(abs(variable_P(:,1)-variable_P(:,2))   )

#if(DIMENSION_GEO == 2)
             do k=1,3
                print"(i4, 3x, a, i1, 2x, 7(E16.6))",j,'L',k, (Er_norme%norme(k)%var(i)%val_var_temps(j),i=1,p) 
             end do
#elif(DIMENSION_GEO == 3)   
             do k=1,3
                print"(i4, 3x, a, i1, 2x, 11(E16.6))",j,'L',k, (Er_norme%norme(k)%var(i)%val_var_temps(j),i=1,p) 
             end do
#endif
             print*,  '  '           
          end if
          variable  (:,2,:) = variable  (:,1,:) 
          variable_P(:,2  ) = variable_P(:,1  )
          dt = dt/2.0_nk
       end do

       do j=2, nbr_it_temps-1          
          do i=1, p
             do k=1,3
                Er_norme%norme(K)%var(i)%ord_temps(j) = &
                     Er_norme%norme(K)%var(i)%val_var_temps(j)/Er_norme%norme(K)%var(i)%val_var_temps(j+1)/2._nk
             end do
          end do
       end do

       print*,'  '
       print*,'ordre de l erreur temporelle en norme'

#if(DIMENSION_GEO == 2)
       print*, ' kk   norme   u1          u2          t11         t12         t22         T           P'
       do j=2, nbr_it_temps-1    
          do k=1,3
             print"(i4, 3x, a, i1, 2x, 7(f12.6))",j,'L',k, (Er_norme%norme(k)%var(i)%ord_temps(j),i=1,p) 
          end do
          print*, '   '
       end do
#elif(DIMENSION_GEO == 3)
       print*,' kk   norme   u1          u2          u3          t11         t12         t12         t22',&
            '         t23         t33         T           P'
       do j=2, nbr_it_temps-1          
          do k=1,3
             print"(i4, 3x, a, i1, 2x, 11(f12.6))",j,'L',k, (Er_norme%norme(k)%var(i)%ord_temps(j),i=1,p) 
          end do
          print*,  '  '
       end do
#endif

    end if







#if (DIMENSION_GEO == 3)
    if (erreur_en_temps==0 .and. erreur_en_espace==0 .and. RICH == 0 .and. SOL_ANA==1)  then                              

!!$       N12=41
!!$       N22=41
!!$       N32=41

       N12=81
       N22=81
       N32=81

       open(11,file="solution_Poiseuille_2D_4x80x80_2",status="old")                            
       do j=1, N22
          do k=1, N32
             l = (j-1)*N32 + k             
             read(11,*) indice2, indice3, u1tau(l,1,1), u1tau(l,1,2), u1tau(l,1,3), u1tau(l,1,4)    
             if (j==1 .or. j==N22 .or. k==1 .or. k==N32) then  
                u1tau(l,1,1) = 0.               
             else if  ( k==1 .or. k==N32) then  
                u1tau(l,1,3) = 0.
             else if  ( j==1 .or. j==N22) then  
                u1tau(l,1,4) = 0.                
             end if

          end do
          read(11,*)
       end do
       close(11)

       do i=1, 4

          if( i==1 )  variable_utau  = "u_1"
          if( i==2 )  variable_utau  = "tau_11"
          if( i==3 )  variable_utau  = "tau_12"
          if( i==4 )  variable_utau  = "tau_13"
          print*, variable_utau

          do ii=1, 8

             if( ii==1 )  then 

                schema  = "   UPWIND1"
                open(11,file="champs_UT_x2x3_up1_f90_2",status="old")                            

             elseif( ii==2 ) then  

                schema  = "   UPWIND2"
                open(11,file="champs_UT_x2x3_up2_f90_2",status="old")                            

             elseif( ii==3 )  then 

                schema  = "   WENO3  "
                open(11,file="champs_UT_x2x3_w3_f90_2",status="old")                             

             elseif( ii==4 ) then 

                schema  = "   HOUC3  "
                open(11,file="champs_UT_x2x3_h3_f90_2",status="old")                            

             else if( ii==5 ) then 

                schema  = "   WENO5  "
                open(11,file="champs_UT_x2x3_w5_f90_2",status="old")                             

             elseif( ii== 6 ) then 

                schema  = "   HOUC5  "
                open(11,file="champs_UT_x2x3_h5_f90_2",status="old")                            

             elseif ( ii==7 ) then  

                schema  = "   WENO7  "
                open(11,file="champs_UT_x2x3_w7_f90_2",status="old")                            

             elseif ( ii==8 ) then 

                schema  = "   HOUC7  "
                open(11,file="champs_UT_x2x3_h7_f90_2",status="old")                            

             endif

             do j=1, N22
                do k=1, N32
                   l = (j-1)*N32 + k             
                   read(11,*) indice2, indice3, u1tau(l,2,1), indice1, indice1, u1tau(l,2,2), u1tau(l,2,3), u1tau(l,2,4)                                       
                   if (j==1 .or. j==N22 .or. k==1 .or. k==N32) then  
                      u1tau(l,2,1) = 0.               
                   else if  ( k==1 .or. k==N32) then  
                      u1tau(l,2,3) = 0.
                   else if  ( j==1 .or. j==N22) then  
                      u1tau(l,2,4) = 0.                
                   end if
                end do
                read(11,*)
             end do
             close(11)

!!$             do l=1, N22*N32
!!$                diff(l,i)= abs( u1tau(l,1,i)  - u1tau(l,2,i)  )
!!$             end do
             do j=1, N22
                do k=1, N32
                   l = (j-1)*N32 + k             
                   diff(l,i)= abs( u1tau(l,1,i)  - u1tau(l,2,i)  )
                end do
             end do




!!$             print"(a,3(a,E12.4),a)", schema , ' & & ',  sum(diff(:,i)),  &
!!$                  ' & & ', sqrt(sum(diff(:,i)*diff(:,i))),  ' & & ', maxval(diff(:,i)),' \\'  


!!$
!!$             print"(a,3(a,E12.4),a)", schema,  ' & ',  sum(diff(:,i)),  &
!!$                  ' &', sqrt(sum(diff(:,i)*diff(:,i))),  ' & ', maxval(diff(:,i)),' &'  

             print"(3(a,E12.4),a)", ' &',  sum(diff(:,i)),  &
                  ' &', sqrt(sum(diff(:,i)*diff(:,i))),  ' &', maxval(diff(:,i)),' \\'    

          end do
          print*, 
       end do
       stop        
    endif
#endif










#if (DIMENSION_GEO == 3)
    if (erreur_en_temps==0 .and. erreur_en_espace==0 .and. RICH == 0 .and. SOL_ANA==0)  then                              





       nN1 = 11 ; nN2 = 11 ; nN3 = 11 ;  nN123= nN1*nN2*nN3
       allocate( u123( nN123, 2, 3 ), diffu123(nN123,3, 4 ) )
       allocate( Erreur_u123(3,4,3),ordre_u123(3, 3, 3) ) 

       pas = 2
       do ii = 1, 4






          if( ii==1 )  then                                 

             open(11,file="champs_3D_11x11x11",status="old")                            
             open(12,file="champs_3D_21x21x21",status="old")   

             do k= 1, nN3
                do j= 1, nN2
                   do i=1, nN1   
                      l = (k-1)*nN1*nN2 + (j-1)*nN1 + i                       
                      read(11,*)  indice1, indice2, indice3, u123(l,1,1), u123(l,1,2), u123(l,1,3)                   
                   enddo
                   read(11,*) 
                enddo
                read(11,*) 
             enddo
             close(11)                                    

          elseif( ii==2 )  then                 
             open(12,file="champs_3D_41x41x41",status="old")                                
          elseif( ii==3 )  then                 
             open(12,file="champs_3D_81x81x81",status="old")                                
          elseif( ii==4 )  then                 
             open(12,file="champs_3D_161x161x161",status="old")                                             
          endif



          nN1 = nN1*2-1
          nN2 = nN2*2-1
          nN3 = nN3*2-1

          allocate( u123_2(nN1*nN2*nN3,3) )                              




          do k= 1, nN3
             do j= 1, nN2
                do i=1, nN1   
                   l = (k-1)*nN1*nN2 + (j-1)*nN1 + i
                   read(12,*)  indice1, indice2, indice3, u123_2(l,1), u123_2(l,2), u123_2(l,3)                       
                enddo
                read(12,*)
             enddo
             read(12,*)
          enddo
          close(12)



          kk =0    
          do k= 1,nN3,pas
             kk = kk + 1
             jj=0
             do j= 1,nN2,pas
                jj = jj + 1 
                iii=0
                do i=1,nN1,pas
                   iii= iii + 1

                   l = (k-1)*nN1*nN2 + (j-1)*nN1 + i                   
                   ll= (kk-1)*11*11 + (jj-1)*11 + iii

                   u123(ll,2,1) = u123_2(l,1)
                   u123(ll,2,2) = u123_2(l,2)
                   u123(ll,2,3) = u123_2(l,3)

                enddo
             enddo
          enddo
          pas = pas * 2

          do i=1, 3
             if( i==1 )  variable_u123  = "u_1"
             if( i==2 )  variable_u123  = "u_2"
             if( i==3 )  variable_u123  = "u_3"                

             do l=1, nN123
                diffu123(l,i,ii)= abs( u123(l,1,i)  -  u123(l,2,i)  )                    
             end do

             Erreur_u123(i,ii,1) =   sum(diffu123(:,i,ii))  
             Erreur_u123(i,ii,2) =   sqrt(sum(diffu123(:,i,ii)*diffu123(:,i,ii)))
             Erreur_u123(i,ii,3) =   maxval(diffu123(:,i,ii))

!!$             print"(a,3(a,E12.4),a)", variable_u123 , ' & & ',  sum(diffu123(:,i,ii)),  &
!!$                  ' & & ', sqrt(sum(diffu123(:,i,ii)*diffu123(:,i,ii))),  ' & & ', maxval(diffu123(:,i,ii)),' \\' 
          end do

          u123(:,1,1) =   u123(:,2,1) 
          u123(:,1,2) =   u123(:,2,2) 
          u123(:,1,3) =   u123(:,2,3)             

          deallocate(u123_2 )              
       end do

       ! calcul ordre
       do i=1, 3 ! les trois variables              
          do ii=1, 3  ! nombre de maillages-1                   
             do k=1,3 ! pour les trois normes                                 
                ordre_u123(i,ii,k)= log(Erreur_u123(i,ii,1)/Erreur_u123(i,ii+1,1))/log(2._nk)                 
             end do
          end do
       end do




!!!! pour rapport latex       
       do i=1, 3                    
          if (i==1) print*, '$u_1$'
          if (i==2) print*, '$u_2$'
          if (i==3) print*, '$u_3$'


          do ii=1,5
             if (ii==1 ) then 

                print "(3(4a, a, 4a, a), a)", &
                     ('  ','&' ,'  ','&','     -      ', &
                     '  ','&' ,'  ','&','     -      ', k=1,3 ), '  \\'

             else if ( ii==2 ) then 

                print "(3(4a, E12.4, 4a, a), a)", &
                     ('  ','&' ,'  ','&',Erreur_u123(i,ii-1,k) , &
                     '  ','&' ,'  ','&','     -      ', k=1,3 ), '  \\' 

             else
                print "(3(4a, E12.4, 4a, f12.4), a)", &
                     ('  ','&' ,'  ','&', Erreur_u123(i,ii-1,k), &
                     '  ','&' ,'  ','&', ordre_u123(i,ii-2,k), k=1,3 ), '  \\' 
             end if
          end do
       end do

       deallocate(u123, diffu123 ) 
       deallocate( Erreur_u123, ordre_u123 ) 
       stop 
    endif

#endif       

  end subroutine erreur_en_champs


















  subroutine CL_Physique 
    implicit none
    character(len=10):: paroi


#if    (TYPE_ECOULEMENT == 3) 
    open(47,file="CL_C_D_C",status="old")
    read(47,*)
#endif 

    read(47,*) paroi                                                 , &
         var%W(u1)%CL_var_bas, var%W(u1)%var_bas, var%W(u1)%q_var_bas, &
         var%W(u2)%CL_var_bas, var%W(u2)%var_bas, var%W(u2)%q_var_bas, &
         var%W(T )%CL_var_bas, var%W(T )%var_bas, var%W(T )%q_var_bas, &
         var%W(P )%CL_var_bas, var%W(P )%var_bas, var%W(P )%q_var_bas, &
#if (DIMENSION_GEO == 3)  
    var%W(u3)%CL_var_bas, var%W(u3)%var_bas, var%W(u3)%q_var_bas 
#elif (DIMENSION_GEO == 2)
    paroi
#endif 



    read(47,*) paroi                                                          , &
         var%W(u1)%CL_var_droite, var%W(u1)%var_droite, var%W(u1)%q_var_droite, &
         var%W(u2)%CL_var_droite, var%W(u2)%var_droite, var%W(u2)%q_var_droite, &
         var%W(T )%CL_var_droite, var%W(T )%var_droite, var%W(T )%q_var_droite, &
         var%W(P )%CL_var_droite, var%W(P )%var_droite, var%W(P )%q_var_droite, &
#if (DIMENSION_GEO == 3)  
    var%W(u3)%CL_var_droite, var%W(u3)%var_droite, var%W(u3)%q_var_droite 
#elif (DIMENSION_GEO == 2)
    paroi
#endif 

    read(47,*) paroi                                                    , &
         var%W(u1)%CL_var_haut, var%W(u1)%var_haut, var%W(u1)%q_var_haut, &
         var%W(u2)%CL_var_haut, var%W(u2)%var_haut, var%W(u2)%q_var_haut, &
         var%W(T )%CL_var_haut, var%W(T )%var_haut, var%W(T )%q_var_haut, &
         var%W(P )%CL_var_haut, var%W(P )%var_haut, var%W(P )%q_var_haut, &
#if (DIMENSION_GEO == 3)  
    var%W(u3)%CL_var_haut, var%W(u3)%var_haut, var%W(u3)%q_var_haut 
#elif (DIMENSION_GEO == 2)
    paroi
#endif 

    read(47,*) paroi                                                          , &
         var%W(u1)%CL_var_gauche, var%W(u1)%var_gauche, var%W(u1)%q_var_gauche, &
         var%W(u2)%CL_var_gauche, var%W(u2)%var_gauche, var%W(u2)%q_var_gauche, &
         var%W(T )%CL_var_gauche, var%W(T )%var_gauche, var%W(T )%q_var_gauche, &
         var%W(P )%CL_var_gauche, var%W(P )%var_gauche, var%W(P )%q_var_gauche, &
#if (DIMENSION_GEO == 3)  
    var%W(u3)%CL_var_gauche, var%W(u3)%var_gauche, var%W(u3)%q_var_gauche 
#elif (DIMENSION_GEO == 2)
    paroi
#endif 

#if (DIMENSION_GEO == 3)  
    read(47,*) paroi                                                      , &
         var%W(u1)%CL_var_avant, var%W(u1)%var_avant, var%W(u1)%q_var_avant, &
         var%W(u2)%CL_var_avant, var%W(u2)%var_avant, var%W(u2)%q_var_avant, &
         var%W(T )%CL_var_avant, var%W(T )%var_avant, var%W(T )%q_var_avant, &
         var%W(P )%CL_var_avant, var%W(P )%var_avant, var%W(P )%q_var_avant, &
         var%W(u3)%CL_var_avant, var%W(u3)%var_avant, var%W(u3)%q_var_avant 

    read(47,*) paroi                                                            , &
         var%W(u1)%CL_var_arriere, var%W(u1)%var_arriere, var%W(u1)%q_var_arriere, &
         var%W(u2)%CL_var_arriere, var%W(u2)%var_arriere, var%W(u2)%q_var_arriere, &
         var%W(T )%CL_var_arriere, var%W(T )%var_arriere, var%W(T )%q_var_arriere, &
         var%W(P )%CL_var_arriere, var%W(P )%var_arriere, var%W(P )%q_var_arriere, &
         var%W(u3)%CL_var_arriere, var%W(u3)%var_arriere, var%W(u3)%q_var_arriere 
#endif 
    close(47)    



  end subroutine CL_Physique









  subroutine adimensionnement_CL
    implicit none
    real :: DeltT, aveT

    if (C_D_C == 1 ) then 

    DeltT              =   var%W(T)%var_bas - var%W(T)%var_haut
    aveT               =  (var%W(T)%var_bas + var%W(T)%var_haut)/2._nk
    
    var%W(T)%var_bas   = (var%W(T)%var_bas  - aveT)/DeltT
    var%W(T)%q_var_bas =  var%W(T)%q_var_bas

    var%W(T)%var_haut  = (var%W(T)%var_haut - aveT)/DeltT
    var%W(T)%q_var_haut=  var%W(T)%q_var_haut    

    var%W(T)%var_droite   = (var%W(T)%var_droite  - aveT)/DeltT
    var%W(T)%q_var_droite =  var%W(T)%q_var_droite  

    var%W(T)%var_gauche   = (var%W(T)%var_gauche  - aveT)/DeltT
    var%W(T)%q_var_gauche =  var%W(T)%q_var_gauche    
#if (DIMENSION_GEO == 3)     
    var%W(T)%var_avant    = (var%W(T)%var_avant  - aveT)/DeltT
    var%W(T)%q_var_avant  =  var%W(T)%q_var_avant

    var%W(T)%var_arriere  = (var%W(T)%var_arriere - aveT)/DeltT
    var%W(T)%q_var_arriere=  var%W(T)%q_var_arriere
#endif
    end if

  end subroutine adimensionnement_CL




REAL FUNCTION UPWIND_1 (coef, var_in, i, j, indice, k) 
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j
    INTEGER, INTENT(IN) :: indice
    INTEGER, OPTIONAL, INTENT(IN) :: k
    real(nk), INTENT(IN)    :: coef
#if (DIMENSION_GEO == 2)
    real(nk),dimension( 0:N1+1,0:N2+1 ),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
    real(nk),dimension( 0:N1+1, 0:N2+1, 0:N3+1 ),intent(in):: var_in
#endif
   if(indice == 1) then

#if (DIMENSION_GEO == 2)
    UPWIND_1 = coef*( var_in(i  ,j      )-var_in(i-1,j      ) ) /dx1_L(i,j)
#elif (DIMENSION_GEO == 3)
    UPWIND_1 = coef*( var_in(i  ,j  ,k  )-var_in(i-1,j  ,k  ) ) /dx1_L(i,j,k)
#endif
   else if (indice == 2) then 

#if (DIMENSION_GEO == 2)
    UPWIND_1 = coef*( var_in(i  ,j      )-var_in(i  ,j-1    ) ) /dx2_L(i,j)
#elif (DIMENSION_GEO == 3)
    UPWIND_1 = coef*( var_in(i  ,j  ,k  )-var_in(i  ,j-1,k  ) ) /dx2_L(i,j,k)
#endif

#if(DIMENSION_GEO == 3)
   else if (indice == 3) then 
    UPWIND_1 = coef*( var_in(i  ,j  ,k  )-var_in(i  ,j  ,k-1) ) /dx3_L(i,j,k)
#endif
   end if
  END FUNCTION UPWIND_1

  REAL FUNCTION BACKWIND_1 (coef, var_in, i, j, indice, k) 
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j
    INTEGER, INTENT(IN) :: indice
    INTEGER, OPTIONAL, INTENT(IN) :: k
    real(nk), INTENT(IN)    :: coef
#if (DIMENSION_GEO == 2)
    real(nk),dimension( 0:N1+1,0:N2+1 ),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
    real(nk),dimension( 0:N1+1, 0:N2+1, 0:N3+1 ),intent(in):: var_in
#endif

   if(indice == 1) then
#if (DIMENSION_GEO == 2)
      BACKWIND_1 = coef*(var_in(i+1,j  )-var_in(i,j  ))/dx1_R(i,j  )
#elif (DIMENSION_GEO == 3)
      BACKWIND_1 = coef*(var_in(i+1,j,k)-var_in(i,j,k))/dx1_R(i,j,k)
#endif
   else if (indice == 2) then 

#if (DIMENSION_GEO == 2)
      BACKWIND_1 = coef*(var_in(i,j+1  )-var_in(i,j  ))/dx2_R(i,j  )
#elif (DIMENSION_GEO == 3)
      BACKWIND_1 = coef*(var_in(i,j+1,k)-var_in(i,j,k))/dx2_R(i,j,k)
#endif

#if(DIMENSION_GEO == 3)
   else if (indice == 3) then 
      BACKWIND_1 = coef*( var_in(i,j,k+1)-var_in(i,j,k))/dx3_R(i,j,k)
#endif
   end if
  END FUNCTION BACKWIND_1

  REAL FUNCTION UPWIND_2 (coef, var_in, i, j, indice, k) 
!!!-------------l-2---x2----l-1---x1---l
!!!              1           2         3
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j
    INTEGER, INTENT(IN) :: indice
    INTEGER, OPTIONAL, INTENT(IN) :: k
    real(nk), INTENT(IN)    :: coef
#if (DIMENSION_GEO == 2)
    real(nk),dimension( 0:N1+1,0:N2+1 ),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
    real(nk),dimension( 0:N1+1, 0:N2+1, 0:N3+1 ),intent(in):: var_in
#endif

   if(indice == 1) then
#if (DIMENSION_GEO == 2)
    UPWIND_2 = coef*( UPW2_1_A(i,j  )*var_in(i-2,j  )-UPW2_1_B(i,j  )*var_in(i-1,j  )+UPW2_1_C(i,j  )*var_in(i,j  ) )
#elif (DIMENSION_GEO == 3)
    UPWIND_2 = coef*( UPW2_1_A(i,j,k)*var_in(i-2,j,k)-UPW2_1_B(i,j,k)*var_in(i-1,j,k)+UPW2_1_C(i,j,k)*var_in(i,j,k) )
#endif
   else if (indice == 2) then 
#if (DIMENSION_GEO == 2)
    UPWIND_2 = coef*( UPW2_2_A(i,j  )*var_in(i,j-2  )-UPW2_2_B(i,j  )*var_in(i,j-1  )+UPW2_2_C(i,j  )*var_in(i,j  ) )
#elif (DIMENSION_GEO == 3)
    UPWIND_2 = coef*( UPW2_2_A(i,j,k)*var_in(i,j-2,k)-UPW2_2_B(i,j,k)*var_in(i,j-1,k)+UPW2_2_C(i,j,k)*var_in(i,j,k) )
#endif

#if(DIMENSION_GEO == 3)
   else if (indice == 3) then 
    UPWIND_2 = coef*( UPW2_3_A(i,j,k)*var_in(i,j,k-2)-UPW2_3_B(i,j,k)*var_in(i,j,k-1)+UPW2_3_C(i,j,k)*var_in(i,j,k) )
#endif
   end if
  END FUNCTION UPWIND_2

  REAL FUNCTION BACKWIND_2 (coef, var_in, i, j, indice, k) 
!!!-------------l-----x1----l+1---x2---l+2
!!!             1            2          3
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j
    INTEGER, INTENT(IN) :: indice
    INTEGER, OPTIONAL, INTENT(IN) :: k
    real(nk), INTENT(IN)    :: coef
#if (DIMENSION_GEO == 2)
    real(nk),dimension( 0:N1+1,0:N2+1 ),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
    real(nk),dimension( 0:N1+1, 0:N2+1, 0:N3+1 ),intent(in):: var_in
#endif
   if(indice == 1) then
    

#if (DIMENSION_GEO == 2)
    BACKWIND_2 = coef*(- BACKW2_1_A(i,j  )*var_in(i+2,j  )+BACKW2_1_B(i,j  )*var_in(i+1,j  )-BACKW2_1_C(i,j  )*var_in(i,j  ) )
#elif (DIMENSION_GEO == 3)
    BACKWIND_2 = coef*(- BACKW2_1_A(i,j,k)*var_in(i+2,j,k)+BACKW2_1_B(i,j,k)*var_in(i+1,j,k)-BACKW2_1_C(i,j,k)*var_in(i,j,k) )
#endif
   else if (indice == 2) then 

#if (DIMENSION_GEO == 2)
    BACKWIND_2 = coef*(- BACKW2_2_A(i,j  )*var_in(i,j+2  )+BACKW2_2_B(i,j  )*var_in(i,j+1  )-BACKW2_2_C(i,j  )*var_in(i,j  ) )
#elif (DIMENSION_GEO == 3)
    BACKWIND_2 = coef*(- BACKW2_2_A(i,j,k)*var_in(i,j+2,k)+BACKW2_2_B(i,j,k)*var_in(i,j+1,k)-BACKW2_2_C(i,j,k)*var_in(i,j,k) )
#endif

#if(DIMENSION_GEO == 3)
   else if (indice == 3) then 
    BACKWIND_2 = coef*(- BACKW2_3_A(i,j,k)*var_in(i,j,k+2)+BACKW2_3_B(i,j,k)*var_in(i,j,k+1)-BACKW2_3_C(i,j,k)*var_in(i,j,k) )
#endif
   end if
  END FUNCTION BACKWIND_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!-------------l-2----------l-1---------l---------l+1---------l+2
 !!!R_W_A1        1            2          3          4           5 
  REAL(nk) FUNCTION WENO_3_L(coef, var_in, i, j, indice, k)
    IMPLICIT NONE
   INTEGER, INTENT(IN) :: i, j
   INTEGER, INTENT(IN) :: indice
   INTEGER, OPTIONAL, INTENT(IN) :: k
   real(nk), INTENT(IN)    :: coef
   real(nk)                :: WENO_omega_0, WENO_omega_1
   real(nk)                :: q1, q2, q3

#if (DIMENSION_GEO == 2)
   real(nk),dimension(0:N1+1,0:N2+1),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
   real(nk),dimension(0:N1+1, 0:N2+1, 0:N3+1),intent(in):: var_in
#endif
    WENO_omega_0 = 1./3.
    WENO_omega_1 = 2./3.

   if (indice == 1) then
#if (DIMENSION_GEO == 2)
    q1 = (var_in(i-1,j) - var_in(i-2,j)) / dx1_L(i-1,j)
    q2 = (var_in(i  ,j) - var_in(i-1,j)) / dx1_L(i  ,j)
    q3 = (var_in(i+1,j) - var_in(i  ,j)) / dx1_R(i  ,j)
      WENO_3_L = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#elif (DIMENSION_GEO == 3)
    q1 = (var_in(i-1,j,k) - var_in(i-2,j,k)) / dx1_L(i-1,j,k)
    q2 = (var_in(i  ,j,k) - var_in(i-1,j,k)) / dx1_L(i  ,j,k)
    q3 = (var_in(i+1,j,k) - var_in(i  ,j,k)) / dx1_R(i  ,j,k)
      WENO_3_L = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#endif
   else if (indice == 2) then
#if (DIMENSION_GEO == 2)
    q1 = (var_in(i,j-1) - var_in(i,j-2)) / dx2_L(i,j-1)
    q2 = (var_in(i,j  ) - var_in(i,j-1)) / dx2_L(i,j  )
    q3 = (var_in(i,j+1) - var_in(i,j  )) / dx2_R(i,j  )
      WENO_3_L = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#elif (DIMENSION_GEO == 3)
    q1 = (var_in(i,j-1,k) - var_in(i,j-2,k)) / dx2_L(i,j-1,k)
    q2 = (var_in(i,j  ,k) - var_in(i,j-1,k)) / dx2_L(i,j  ,k)
    q3 = (var_in(i,j+1,k) - var_in(i,j  ,k)) / dx2_R(i,j  ,k)
      WENO_3_L = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#endif
#if (DIMENSION_GEO == 3)
   else if (indice == 3) then 
    q1 = (var_in(i,j,k-1) - var_in(i,j,k-2)) / dx3_L(i,j,k-1)
    q2 = (var_in(i,j,k  ) - var_in(i,j,k-1)) / dx3_L(i,j,k  )
    q3 = (var_in(i,j,k+1) - var_in(i,j,k  )) / dx3_R(i,j,k  )
      WENO_3_L = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#endif
   end if
!   return
  END FUNCTION WENO_3_L

  REAL(nk) FUNCTION WENO_3_R(coef, var_in, i, j, indice, k)
    IMPLICIT NONE
   INTEGER, INTENT(IN) :: i, j
   INTEGER, INTENT(IN) :: indice
   INTEGER, OPTIONAL, INTENT(IN) :: k
   real(nk), INTENT(IN)    :: coef
   real(nk)                :: WENO_omega_0, WENO_omega_1
   real(nk)                :: q1, q2, q3
#if (DIMENSION_GEO == 2)
   real(nk),dimension(0:N1+1,0:N2+1),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
   real(nk),dimension(0:N1+1, 0:N2+1, 0:N3+1),intent(in):: var_in
#endif
    WENO_omega_0 = 1./3.
    WENO_omega_1 = 2./3.

   if (indice == 1) then
#if (DIMENSION_GEO == 2)
    q1 = (var_in(i+2,j) - var_in(i+1,j)) / dx1_R(i+1,j)
    q2 = (var_in(i+1,j) - var_in(i  ,j)) / dx1_R(i  ,j)
    q3 = (var_in(i  ,j) - var_in(i-1,j)) / dx1_L(i  ,j)
      WENO_3_R = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#elif (DIMENSION_GEO == 3)
    q1 = (var_in(i+2,j,k) - var_in(i+1,j,k)) / dx1_R(i+1,j,k)
    q2 = (var_in(i+1,j,k) - var_in(i  ,j,k)) / dx1_R(i  ,j,k)
    q3 = (var_in(i  ,j,k) - var_in(i-1,j,k)) / dx1_L(i  ,j,k)
      WENO_3_R = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#endif
   else if (indice == 2) then
#if (DIMENSION_GEO == 2)
    q1 = (var_in(i,j+2) - var_in(i,j+1)) / dx2_R(i,j+1)
    q2 = (var_in(i,j+1) - var_in(i,j  )) / dx2_R(i,j  )
    q3 = (var_in(i,j  ) - var_in(i,j-1)) / dx2_L(i,j  )
      WENO_3_R = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#elif (DIMENSION_GEO == 3)
    q1 = (var_in(i,j+2,k) - var_in(i,j+1,k)) / dx2_R(i,j+1,k)
    q2 = (var_in(i,j+1,k) - var_in(i,j  ,k)) / dx2_R(i,j  ,k)
    q3 = (var_in(i,j  ,k) - var_in(i,j-1,k)) / dx2_L(i,j  ,k)
      WENO_3_R = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#endif
#if (DIMENSION_GEO == 3)
   else if (indice == 3) then 
    q1 = (var_in(i,j,k+2) - var_in(i,j,k+1)) / dx3_R(i,j,k+1)
    q2 = (var_in(i,j,k+1) - var_in(i,j,k  )) / dx3_R(i,j,k  )
    q3 = (var_in(i,j,k  ) - var_in(i,j,k-1)) / dx3_L(i,j,k  )
      WENO_3_R = coef*( WENO_omega_0*( -0.5_nk*q1 + 1.5_nk*q2 ) &
                      + WENO_omega_1*(  0.5_nk*q2 + 0.5_nk*q3 ) )
#endif
   end if
!   return
  END FUNCTION WENO_3_R




!  --------------l-3----------l-2-----------l-1-----------l-----------l+1-----------l+2-----------l+3-------
!  R_W_A2         7            6             5            1            2             3             4   
  REAL(nk) FUNCTION WENO_5_L(coef, var_in, i, j, indice, k)
    IMPLICIT NONE
   INTEGER, INTENT(IN) :: i, j
   INTEGER, INTENT(IN) :: indice
   INTEGER, OPTIONAL, INTENT(IN) :: k
   real(nk), INTENT(IN)    :: coef
   real(nk)                :: WENO_omega_0, WENO_omega_1, WENO_omega_2
   real(nk)                :: q1, q2, q3, q4, q5

#if (DIMENSION_GEO == 2)
   real(nk),dimension(0:N1+1,0:N2+1),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
   real(nk),dimension(0:N1+1, 0:N2+1, 0:N3+1),intent(in):: var_in
#endif
    WENO_omega_0 = 1./10.
    WENO_omega_1 = 6./10.
    WENO_omega_2 = 3./10.

   if (indice == 1) then
#if (DIMENSION_GEO == 2)
    q1 = (var_in(i-2,j) - var_in(i-3,j)) / dx1_L(i-2,j)
    q2 = (var_in(i-1,j) - var_in(i-2,j)) / dx1_L(i-1,j)
    q3 = (var_in(i  ,j) - var_in(i-1,j)) / dx1_L(i  ,j)
    q4 = (var_in(i+1,j) - var_in(i  ,j)) / dx1_R(i  ,j)
    q5 = (var_in(i+2,j) - var_in(i+1,j)) / dx1_R(i+1,j)
      WENO_5_L = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#elif (DIMENSION_GEO == 3)
    q1 = (var_in(i-2,j,k) - var_in(i-3,j,k)) / dx1_L(i-2,j,k)
    q2 = (var_in(i-1,j,k) - var_in(i-2,j,k)) / dx1_L(i-1,j,k)
    q3 = (var_in(i  ,j,k) - var_in(i-1,j,k)) / dx1_L(i  ,j,k)
    q4 = (var_in(i+1,j,k) - var_in(i  ,j,k)) / dx1_R(i  ,j,k)
    q5 = (var_in(i+2,j,k) - var_in(i+1,j,k)) / dx1_R(i+1,j,k)
      WENO_5_L = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#endif
   else if (indice == 2) then
#if (DIMENSION_GEO == 2)
    q1 = (var_in(i,j-2) - var_in(i,j-3)) / dx2_L(i,j-2)
    q2 = (var_in(i,j-1) - var_in(i,j-2)) / dx2_L(i,j-1)
    q3 = (var_in(i,j  ) - var_in(i,j-1)) / dx2_L(i,j  )
    q4 = (var_in(i,j+1) - var_in(i,j  )) / dx2_R(i,j  )
    q5 = (var_in(i,j+2) - var_in(i,j+1)) / dx2_R(i,j+1)
      WENO_5_L = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#elif (DIMENSION_GEO == 3)
    q1 = (var_in(i,j-2,k) - var_in(i,j-3,k)) / dx2_L(i,j-2,k)
    q2 = (var_in(i,j-1,k) - var_in(i,j-2,k)) / dx2_L(i,j-1,k)
    q3 = (var_in(i,j  ,k) - var_in(i,j-1,k)) / dx2_L(i,j  ,k)
    q4 = (var_in(i,j+1,k) - var_in(i,j  ,k)) / dx2_R(i,j  ,k)
    q5 = (var_in(i,j+2,k) - var_in(i,j+1,k)) / dx2_R(i,j+1,k)
      WENO_5_L = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#endif
#if (DIMENSION_GEO == 3)
   else if (indice == 3) then 
    q1 = (var_in(i,j,k-2) - var_in(i,j,k-3)) / dx3_L(i,j,k-2)
    q2 = (var_in(i,j,k-1) - var_in(i,j,k-2)) / dx3_L(i,j,k-1)
    q3 = (var_in(i,j,k  ) - var_in(i,j,k-1)) / dx3_L(i,j,k  )
    q4 = (var_in(i,j,k+1) - var_in(i,j,k  )) / dx3_R(i,j,k  )
    q5 = (var_in(i,j,k+2) - var_in(i,j,k+1)) / dx3_R(i,j,k+1)
      WENO_5_L = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#endif
   end if
!   return
  END FUNCTION WENO_5_L

  REAL(nk) FUNCTION WENO_5_R(coef, var_in, i, j, indice, k)
    IMPLICIT NONE
   INTEGER, INTENT(IN) :: i, j
   INTEGER, INTENT(IN) :: indice
   INTEGER, OPTIONAL, INTENT(IN) :: k
   real(nk), INTENT(IN)    :: coef
   real(nk)                :: WENO_omega_0, WENO_omega_1, WENO_omega_2
   real(nk)                :: q1, q2, q3, q4, q5
#if (DIMENSION_GEO == 2)
   real(nk),dimension(0:N1+1,0:N2+1),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
   real(nk),dimension(0:N1+1, 0:N2+1, 0:N3+1),intent(in):: var_in
#endif
    WENO_omega_0 = 1./10.
    WENO_omega_1 = 6./10.
    WENO_omega_2 = 3./10.

   if (indice == 1) then
#if (DIMENSION_GEO == 2)
    q1 = (var_in(i+3,j) - var_in(i+2,j)) / dx1_R(i+2,j)
    q2 = (var_in(i+2,j) - var_in(i+1,j)) / dx1_R(i+1,j)
    q3 = (var_in(i+1,j) - var_in(i  ,j)) / dx1_R(i  ,j)
    q4 = (var_in(i  ,j) - var_in(i-1,j)) / dx1_L(i  ,j)
    q5 = (var_in(i-1,j) - var_in(i-2,j)) / dx1_L(i-1,j)
      WENO_5_R = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#elif (DIMENSION_GEO == 3)
    q1 = (var_in(i+3,j,k) - var_in(i+2,j,k)) / dx1_R(i+2,j,k)
    q2 = (var_in(i+2,j,k) - var_in(i+1,j,k)) / dx1_R(i+1,j,k)
    q3 = (var_in(i+1,j,k) - var_in(i  ,j,k)) / dx1_R(i  ,j,k)
    q4 = (var_in(i  ,j,k) - var_in(i-1,j,k)) / dx1_L(i  ,j,k)
    q5 = (var_in(i-1,j,k) - var_in(i-2,j,k)) / dx1_L(i-1,j,k)
      WENO_5_R = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#endif
   else if (indice == 2) then
#if (DIMENSION_GEO == 2)
    q1 = (var_in(i,j+3) - var_in(i,j+2)) / dx2_R(i,j+2)
    q2 = (var_in(i,j+2) - var_in(i,j+1)) / dx2_R(i,j+1)
    q3 = (var_in(i,j+1) - var_in(i,j  )) / dx2_R(i,j  )
    q4 = (var_in(i,j  ) - var_in(i,j-1)) / dx2_L(i,j  )
    q5 = (var_in(i,j-1) - var_in(i,j-2)) / dx2_L(i,j-1)
      WENO_5_R = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#elif (DIMENSION_GEO == 3)
    q1 = (var_in(i,j+3,k) - var_in(i,j+2,k)) / dx2_R(i,j+2,k)
    q2 = (var_in(i,j+2,k) - var_in(i,j+1,k)) / dx2_R(i,j+1,k)
    q3 = (var_in(i,j+1,k) - var_in(i,j  ,k)) / dx2_R(i,j  ,k)
    q4 = (var_in(i,j  ,k) - var_in(i,j-1,k)) / dx2_L(i,j  ,k)
    q5 = (var_in(i,j-1,k) - var_in(i,j-2,k)) / dx2_L(i,j-1,k)
      WENO_5_R = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#endif
#if (DIMENSION_GEO == 3)
   else if (indice == 3) then 
    q1 = (var_in(i,j,k+3) - var_in(i,j,k+2)) / dx3_R(i,j,k+2)
    q2 = (var_in(i,j,k+2) - var_in(i,j,k+1)) / dx3_R(i,j,k+1)
    q3 = (var_in(i,j,k+1) - var_in(i,j,k  )) / dx3_R(i,j,k  )
    q4 = (var_in(i,j,k  ) - var_in(i,j,k-1)) / dx3_L(i,j,k  )
    q5 = (var_in(i,j,k-1) - var_in(i,j,k-2)) / dx3_L(i,j,k-1)
      WENO_5_R = coef*( WENO_omega_0*( 2._nk*q1 - 7._nk*q2 + 11._nk*q3 )/6._nk &
                      + WENO_omega_1*(-      q2 + 5._nk*q3 +  2._nk*q4 )/6._nk &
                      + WENO_omega_2*( 2._nk*q3 + 5._nk*q4 -        q5 )/6._nk )
#endif
   end if
  END FUNCTION WENO_5_R


  REAL(nk) FUNCTION CENTRE_2(coef, var_in, i, j, indice, k)
!!!----------l-1---x1---l----x2---l+1-----
!!!           1         2          3          
    IMPLICIT NONE
   INTEGER, INTENT(IN) :: i, j
   INTEGER, INTENT(IN) :: indice
   INTEGER, OPTIONAL, INTENT(IN) :: k
   real(nk), INTENT(IN)    :: coef
#if (DIMENSION_GEO == 2)
   real(nk),dimension(0:N1+1,0:N2+1),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
   real(nk),dimension(0:N1+1, 0:N2+1, 0:N3+1),intent(in):: var_in
#endif

   if (indice == 1) then
#if (DIMENSION_GEO == 2)
    centre_2 = coef*( CEN2_1_A(i,j  )*var_in(i+1,j  ) - CEN2_1_B(i,j  )*var_in(i,j  ) - CEN2_1_C(i,j  )*var_in(i-1,j  ) )
#elif (DIMENSION_GEO == 3)
    centre_2 = coef*( CEN2_1_A(i,j,k)*var_in(i+1,j,k) - CEN2_1_B(i,j,k)*var_in(i,j,k) - CEN2_1_C(i,j,k)*var_in(i-1,j,k) )
#endif
   else if (indice == 2) then
      
#if (DIMENSION_GEO == 2)
    centre_2 = coef*( CEN2_2_A(i,j  )*var_in(i,j+1  ) - CEN2_2_B(i,j  )*var_in(i,j  ) - CEN2_2_C(i,j  )*var_in(i,j-1  ) )
#elif (DIMENSION_GEO == 3)
    centre_2 = coef*( CEN2_2_A(i,j,k)*var_in(i,j+1,k) - CEN2_2_B(i,j,k)*var_in(i,j,k) - CEN2_2_C(i,j,k)*var_in(i,j-1,k) )
#endif
#if (DIMENSION_GEO == 3)
   else if (indice == 3) then 
    centre_2 = coef*( CEN2_3_A(i,j,k)*var_in(i,j,k+1) - CEN2_3_B(i,j,k)*var_in(i,j,k) - CEN2_3_C(i,j,k)*var_in(i,j,k-1) )
#endif
   end if
!   return
END FUNCTION CENTRE_2

  REAL FUNCTION CENTRE_4 (coef, var_in, i, j, indice, k) 
!!!-------------l-2----x1- --l-1----x2---l----x3---l+1----x4---l+2
!!!R_W_A1        1            2          3          4           5  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j
    INTEGER, INTENT(IN) :: indice
    INTEGER, OPTIONAL, INTENT(IN) :: k
    real(nk), INTENT(IN)    :: coef
#if (DIMENSION_GEO == 2)
    real(nk),dimension( 0:N1+1,0:N2+1 ),intent(in):: var_in
#elif (DIMENSION_GEO == 3)
    real(nk),dimension( 0:N1+1, 0:N2+1, 0:N3+1 ),intent(in):: var_in
#endif
   if(indice == 1) then

#if (DIMENSION_GEO == 2)
    centre_4 = coef*(-CEN4_1_A(i,j  )*var_in(i-2,j  )+CEN4_1_B(i,j  )*var_in(i-1,j  )-CEN4_1_C(i,j  )*var_in(i,j  )&
                     -CEN4_1_D(i,j  )*var_in(i+1,j  )+CEN4_1_E(i,j  )*var_in(i+2,j  ))
#elif (DIMENSION_GEO == 3)
    centre_4 = coef*(-CEN4_1_A(i,j,k)*var_in(i-2,j,k)+CEN4_1_B(i,j,k)*var_in(i-1,j,k)-CEN4_1_C(i,j,k)*var_in(i,j,k)&
                     -CEN4_1_D(i,j,k)*var_in(i+1,j,k)+CEN4_1_E(i,j,k)*var_in(i+2,j,k))
#endif
   else if (indice == 2) then 

#if (DIMENSION_GEO == 2)
    centre_4 = coef*(-CEN4_2_A(i,j  )*var_in(i,j-2  )+CEN4_2_B(i,j  )*var_in(i,j-1  )-CEN4_2_C(i,j  )*var_in(i,j  )&
                     -CEN4_2_D(i,j  )*var_in(i,j+1  )+CEN4_2_E(i,j  )*var_in(i,j+2  ))
#elif (DIMENSION_GEO == 3)
    centre_4 = coef*(-CEN4_2_A(i,j,k)*var_in(i,j-2,k)+CEN4_2_B(i,j,k)*var_in(i,j-1,k)-CEN4_2_C(i,j,k)*var_in(i,j,k)&
                     -CEN4_2_D(i,j,k)*var_in(i,j+1,k)+CEN4_2_E(i,j,k)*var_in(i,j+2,k))
#endif

#if(DIMENSION_GEO == 3)
   else if (indice == 3) then 
    centre_4 = coef*(-CEN4_3_A(i,j,k)*var_in(i,j,k-2)+CEN4_3_B(i,j,k)*var_in(i,j,k-1)-CEN4_3_C(i,j,k)*var_in(i,j,k)&
                     -CEN4_3_D(i,j,k)*var_in(i,j,k+1)+CEN4_3_E(i,j,k)*var_in(i,j,k+2) )
#endif
   end if
  END FUNCTION CENTRE_4

  

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine Read_initial_data_UT( u1, u2, T, t11, t12, t22, u3, t13, t23, t33)
   implicit none
   integer :: IOstatus
   integer :: i, j, k
   type(variable)              :: u1, u2, t11, t12, t22, T
   type(variable), optional    :: u3, t13, t23, t33
   real(nk):: Flanki

   u1 %var = 0. ; u2 %var = 0. ; T  %var = 0. ; t11%var = 0. ; t12%var = 0.  ; t22%var = 0.
   u1 %var0= 0. ; u2 %var0= 0. ; T  %var0= 0. ; t11%var0= 0. ; t12%var0= 0.  ; t22%var0= 0
#if (DIMENSION_GEO ==3)
   u3 %var = 0. ; t13%var = 0. ; t23%var = 0. ; t33%var = 0
   u3 %var0= 0. ; t13%var0= 0. ; t23%var0= 0. ; t33%var = 0.
#endif
#if (DIMENSION_GEO ==2)
#if (VISCOELASTIC_MODELS == 0)
   open (82,file="champs_UT", Form = "UNFORMATTED", access="SEQUENTIAL" , status="old",action='read')      
      do j= 1,N2
         do i= 1,N1      
            read(82,IOSTAT=IOstatus) Flanki, Flanki, &
            u1%var (i,j), u2%var (i,j), Flanki, Flanki, Flanki, T%var (i,j), &
            u1%var0(i,j), u2%var0(i,j), Flanki, Flanki, Flanki, T%var0(i,j)
         end do   
            read(82,IOSTAT=IOstatus)
      end do
   close(82)
#else 
   open (82,file="champs_UT", Form = "UNFORMATTED", access="SEQUENTIAL" , status="old",action='read')      
      do j= 1,N2
         do i= 1,N1         
            read(82,IOSTAT=IOstatus) Flanki, Flanki, &
            u1%var (i,j), u2%var (i,j), t11%var (i,j), t12%var (i,j), t22%var (i,j), T%var (i,j), &
            u1%var0(i,j), u2%var0(i,j), t11%var0(i,j), t12%var0(i,j), t22%var0(i,j), T%var0(i,j)
         end do   
            read(82,IOSTAT=IOstatus)
      end do
   close(82)
#endif
#elif (DIMENSION_GEO ==3)
   open (82,file="champs_UT", Form = "UNFORMATTED", access="SEQUENTIAL" , status="old",action='read')      
   do k=1,N3
      do j= 1,N2
         do i= 1,N1         
            read(82,IOSTAT=IOstatus) Flanki, Flanki, &
            u1 %var (i,j,k), u2 %var (i,j,k), u3 %var (i,j,k), t11%var (i,j,k), t12%var (i,j,k), &
            t13%var (i,j,k), t22%var (i,j,k), t23%var (i,j,k), t33%var (i,j,k), T  %var (i,j,k), &
            u1 %var0(i,j,k), u2 %var0(i,j,k), u3 %var0(i,j,k), t11%var0(i,j,k), t12%var0(i,j,k), &
            t13%var0(i,j,k), t22%var0(i,j,k), t23%var0(i,j,k), t33%var0(i,j,k), T  %var0(i,j,k)
         end do   
            read(82,IOSTAT=IOstatus)
      end do
         read(82,IOSTAT=IOstatus)
   end do
   close(82)
#endif
  end subroutine Read_initial_data_UT



  subroutine Read_initial_data_P( P)
   implicit none
   integer :: IOstatus
   integer :: i, j, k
   type(variable):: P
   real(nk):: Flanki

   P%var = 0. ; P%var0 = 0. 
#if (DIMENSION_GEO ==2)
   open (81,file="champ_P", Form = "UNFORMATTED", access="SEQUENTIAL" , status="old",action='read')      
     do j= 1,N2-1
         do i= 1,N1-1     
            read(81,IOSTAT=IOstatus) Flanki, Flanki, P%var(i,j), P%var0(i,j)
         end do   
            read(81,IOSTAT=IOstatus)
     end do
   close(81)
#elif (DIMENSION_GEO ==3)
   open (81,file="champ_P", Form = "UNFORMATTED", access="SEQUENTIAL" , status="old",action='read') 
   do k=1,N3-1     
     do j= 1,N2-1
         do i= 1,N1-1     
            read(81,IOSTAT=IOstatus) Flanki, Flanki, P%var(i,j,k), P%var0(i,j,k)
         end do   
            read(81,IOSTAT=IOstatus)
     end do
        read(81,*,IOSTAT=IOstatus)
   end do
   close(81)
#endif
   
  end subroutine Read_initial_data_P


  subroutine Calculs_Average_Nusselt(T)
   implicit none
   integer  :: i, j
#if ( DIMENSION_GEO == 2)
   real(nk), dimension(0:N1+1, 0:N2+1), intent(in):: T
   real(nk)  :: Total_Nu_paroi1, Total_Nu_paroi2
   real(nk), dimension(1:N1) :: Nu_paroi1, Nu_paroi2
#elif (DIMENSION_GEO ==3)
   real(nk), dimension(0:N1+1, 0:N2+1, 0:N3+1), intent(in):: T
   real(nk),dimension(:,:),allocatable  :: Nu_paroi1, Nu_paroi2
#endif
   real(nk)::  xi, yi, zi
   real(nk)::  Peak_Nu, Slope_Nu
   integer ::  size

#if ( DIMENSION_GEO == 2 )
   Nu_paroi1(:) = 0.0 ; Nu_paroi2(:) = 0.0 ; xi = 0.;Total_Nu_paroi1=0.;Total_Nu_paroi2=0.
     do i = 1, N1
        Nu_paroi1(i)  = - ( -3._nk*T(i,1 ) + 4._nk * T(i,1+1 ) - T(i,1+2 ) )/2._nk/dx2_R(i,1 )
        Nu_paroi2(i)  =   (  3._nk*T(i,N2) - 4._nk * T(i,N2-1) + T(i,N2-2) )/2._nk/dx2_L(i,N2) 
        xi = xi + dx1_R(i,1)
     end do
   
     do i = 1, N1-1
        Total_Nu_paroi1 =Total_Nu_paroi1 + dx1_R(i,1 )*(Nu_paroi1(i)+Nu_paroi1(i+1))*0.5_nk
        Total_Nu_paroi2 =Total_Nu_paroi2 + dx1_R(i,N2)*(Nu_paroi2(i)+Nu_paroi2(i+1))*0.5_nk
     enddo; Total_Nu_paroi1=Total_Nu_paroi1/2.; Total_Nu_paroi2=Total_Nu_paroi2/2.
     Time_period_integral_Nu = Time_period_integral_Nu + Total_Nu_paroi1*dt
       write(34,*)  it*dt, Total_Nu_paroi1, Total_Nu_paroi2

      Peak_Nu = (Total_Nu_paroi1 - Total_Nu_paroi1_0)*(Total_Nu_paroi1_0 - Total_Nu_paroi1_00)
      Slope_Nu = Total_Nu_paroi1_0 - Total_Nu_paroi1_00

      if (Peak_Nu < 0. .and. Slope_Nu > 0. ) then
         if ( var%W(2)%var(N1/8, N2/2) > 0. ) then
!!!!!!!! Nusselt number period and time periodic integral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    !!!   Here is Two period integrall for Nu, and period of Nu is half of that for velocity
            Period_t = dt*it - Period_t
            open (25,file="socillating_pattern_period_amplitude")
            write(25,*) Period_t, Time_period_integral_Nu/Period_t/L2
            close(25)
            Period_t = dt*it ; Time_period_integral_Nu = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call Calculs_Local_Nusselt (T, it, 11)
         else if ( var%W(2)%var(N1/8, N2/2) < 0. ) then
            call Calculs_Local_Nusselt (T, it, 12)
         end if
      else if (Peak_Nu < 0. .and. Slope_Nu < 0.) then
         if ( var%W(2)%var(N1/8, N2/2) > 0. ) then
            call Calculs_Local_Nusselt (T, it, 21)
         else if ( var%W(2)%var(N1/8, N2/2) < 0. ) then
            call Calculs_Local_Nusselt (T, it, 22)   
         end if
      end if
      Total_Nu_paroi1_00 = Total_Nu_paroi1_0 ; Total_Nu_paroi1_0 = Total_Nu_paroi1
#elif ( DIMENSION == 3 )
  Nu_paroi1(:,:) = 0.0 ; Nu_paroi2(:,:) = 0.0 ; xi = 0.; zi = 0.;Total_Nu_paroi1=0.;Total_Nu_paroi2=0.
     do i = 1, N1
      do k = 1, N3
        Nu_paroi1(i,k)  = - ( -3._nk*T(i,1 ,k) + 4._nk * T(i,1+1 ,k) - T(i,1+2 ,k) )/2._nk/dx2_R(i,1 ,k)
        Nu_paroi2(i,k)  =   (  3._nk*T(i,N2,k) - 4._nk * T(i,N2-1,k) + T(i,N2-2,k) )/2._nk/dx2_L(i,N2,k) 
        zi = zi + dx3_L(i,1,k)
      end do
       xi = xi + dx1_R(i,1,k)
     end do
   
     do i = 1, N1-1
      do k = 1, N3-1
        Total_Nu_paroi1 =Total_Nu_paroi1 + dx1_R(i,1 ,k)*dx2_R(i,1 ,k)*&
                        (Nu_paroi1(i,k)+Nu_paroi1(i+1,k)+Nu_paroi1(i,k+1)+Nu_paroi1(i+1,k+1))*0.25_nk
        Total_Nu_paroi2 =Total_Nu_paroi2 + dx1_R(i,N2,k)*dx2_R(i,N2,k)*&
                        (Nu_paroi2(i,k)+Nu_paroi2(i+1,k)+Nu_paroi2(i,k+1)+Nu_paroi2(i+1,k+1))*0.25_nk
     enddo;enddo; Total_Nu_paroi1=Total_Nu_paroi1/L1/L3; Total_Nu_paroi2=Total_Nu_paroi2/L1/L3
     Time_period_integral_Nu = Time_period_integral_Nu + Total_Nu_paroi1*dt
       write(34,*)  it*dt, Total_Nu_paroi1, Total_Nu_paroi2

      Peak_Nu = (Total_Nu_paroi1 - Total_Nu_paroi1_0)*(Total_Nu_paroi1_0 - Total_Nu_paroi1_00)
      Slope_Nu = Total_Nu_paroi1_0 - Total_Nu_paroi1_00

      if (Peak_Nu < 0. .and. Slope_Nu > 0. ) then
         if ( var%W(2)%var(N1/8, N2/2) > 0. ) then
!!!!!!!! Nusselt number period and time periodic integral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    !!!   Here is Two period integrall for Nu, and period of Nu is half of that for velocity
            Period_t = dt*it - Period_t
            open (25,file="socillating_pattern_period_amplitude")
            write(25,*) Period_t, Time_period_integral_Nu
            close(25)
            Period_t = dt*it ; Time_period_integral_Nu = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call Calculs_Local_Nusselt (T, it, 11)
         else if ( var%W(2)%var(N1/8, N2/2) < 0. ) then
            call Calculs_Local_Nusselt (T, it, 12)
         end if
      else if (Peak_Nu < 0. .and. Slope_Nu < 0.) then
         if ( var%W(2)%var(N1/8, N2/2) > 0. ) then
            call Calculs_Local_Nusselt (T, it, 21)
         else if ( var%W(2)%var(N1/8, N2/2) < 0. ) then
            call Calculs_Local_Nusselt (T, it, 22)   
         end if
      end if
      Total_Nu_paroi1_00 = Total_Nu_paroi1_0 ; Total_Nu_paroi1_0 = Total_Nu_paroi1
#endif

   end subroutine Calculs_Average_Nusselt

   subroutine Calculs_Local_Nusselt(T, it, indice)
#if ( DIMENSION_GEO == 2)
   real(nk), dimension(0:N1+1, 0:N2+1), intent(in):: T
#elif (DIMENSION_GEO ==3)
   real(nk), dimension(0:N1+1, 0:N2+1, 0:N3+1), intent(in):: T
#endif
    integer, intent(in)       :: it
    integer, optional, intent(in):: indice
    character(9) :: cTemp
    real(nk)  :: Total_Nu_paroi1, Total_Nu_paroi2
    real(nk), dimension(1:N1) :: Nu_paroi1, Nu_paroi2
    real(nk)::                             xi, yi, zi
    
    write( cTemp,'(i9)' ) it
#if ( DIMENSION_GEO == 2 )
    if (present(indice)) then
      if      ( indice == 11) then
         open (81,file="Top1_Peak_Nu_local_paroi_bas_haut") 
      else if ( indice == 12) then
         open (81,file="Top2_Peak_Nu_local_paroi_bas_haut") 
      else if ( indice == 21) then
         open (81,file="Bottom1_Peak_Nu_local_paroi_bas_haut") 
      else if ( indice == 22) then
         open (81,file="Bottom2_Peak_Nu_local_paroi_bas_haut") 
      endif
    else 
      open (81,file=trim(adjustl(cTemp))//"_Nu_local_paroi_bas_haut") 
    endif
     Nu_paroi1(:) = 0.0 ; Nu_paroi2(:) = 0.0 ; xi = - dx1_R(1,1) ;Total_Nu_paroi1=0.;Total_Nu_paroi2=0.
     do i = 1, N1
        Nu_paroi1(i)  = - ( -3._nk*T(i,1 ) + 4._nk * T(i,1 +1) - T(i,1 +2) )/2._nk/dx2_R(i,1 +1)
        Nu_paroi2(i)  =   (  3._nk*T(i,N2) - 4._nk * T(i,N2-1) + T(i,N2-2) )/2._nk/dx2_L(i,N2-1) 
        xi = xi + dx1_R(i,1) 
        write(81,*)  xi, Nu_paroi1(i), Nu_paroi2(i)
     end do
    close(81)
#elif ( DIMENSION_GEO == 3 )
    if (present(indice)) then
      if      ( indice == 11) then
         open (81,file="Top1_Peak_Nu_local_paroi_bas_haut") 
      else if ( indice == 12) then
         open (81,file="Top2_Peak_Nu_local_paroi_bas_haut") 
      else if ( indice == 21) then
         open (81,file="Bottom1_Peak_Nu_local_paroi_bas_haut") 
      else if ( indice == 22) then
         open (81,file="Bottom2_Peak_Nu_local_paroi_bas_haut") 
      endif
    else 
      open (81,file=trim(adjustl(cTemp))//"_Nu_local_paroi_bas_haut") 
    endif
     Nu_paroi1(:) = 0.0 ; Nu_paroi2(:) = 0.0 ; xi = - dx1_R(1,1,1) ;Total_Nu_paroi1=0.;Total_Nu_paroi2=0.
     do i = 1, N1
      do k = 1, N3
        Nu_paroi1(i)  = - ( -3._nk*T(i,1 ,k) + 4._nk * T(i,1 +1,k) - T(i,1 +2,k) )/2._nk/dx2_R(i,1 +1, k)
        Nu_paroi2(i)  =   (  3._nk*T(i,N2,k) - 4._nk * T(i,N2-1,k) + T(i,N2-2,k) )/2._nk/dx2_L(i,N2-1, k) 
        zi = zi + dx3_R(i,1,k)
        write(81,*)  xi, Nu_paroi1(i), Nu_paroi2(i)
      end do
       xi = xi + dx1_R(i,1,k) 
     end do
    close(81)
#endif
   end subroutine Calculs_Local_Nusselt

#if (DIMENSION_GEO ==2)
  subroutine kinetic_energy_budget( u1, u2, t11, t12, t22, T, P)
   implicit none
   integer :: i, j, k
   type(variable):: u1, u2, t11, t12, t22, T, P
   real(nk), dimension(0:N1+1,0:N2+1) :: D11, D12, D21, D22
   real(nk)                       :: Global_E, Global_D, Global_F, Global_G, Global_V, Diff_E
   real(nk)                       :: Total_normal_tau, Total_tau, Diff_Tau

   D11(:,:)=0.; D22(:,:)=0.
   Global_D=0.; Global_F=0.; Global_G=0.; Global_V=0.

   
#if (DIMENSION_GEO ==2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! global kinetic energy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 0, N1+1
      do j = 0, N2+1
        E_global(i,j) = 0.5_nk*(u1%var(i,j)*u1%var(i,j) + u2%var(i,j)*u2%var(i,j))
      enddo
    enddo
    do i = 1, N1-1
      do j = 1, N2-1
        Global_E = Global_E + dx1_R(i,j)*dx2_R(i,j) * &
                  (E_global(i,j) + E_global(i+1,j) + E_global(i,j+1) + E_global(i+1,j+1))*0.25_nk
      enddo
    enddo
    Diff_E = (Global_E - Last_Global_E) / dt
    Last_Global_E = Global_E
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! kinetic diffusion !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, N1
      do j = 1, N2 
         D11(i,j) = .25_nk*(P%var(i,j) + P%var(i-1,j) + P%var(i,j-1) + P%var(i-1,j-1))*u1 %var(i,j)    &
                   - u1%nomb_ad * (E_global(i+1,j) - E_global(i-1,j))/dx1_L(i,j)/2.


         D22(i,j) = .25_nk*(P%var(i,j) + P%var(i-1,j) + P%var(i,j-1) + P%var(i-1,j-1))*u2 %var(i,j)    &
                   - u1%nomb_ad * (E_global(i+1,j) - E_global(i,j-1))/dx2_L(i,j)/2.
      enddo
    enddo
    do i = 2, N1-1
      do j = 2, N2-1
          D_diffusion(i,j) = - (D11(i+1,j) - D11(i-1,j))/dx1_L(i,j)/2. - (D22(i,j+1) - D22(i,j-1))/dx2_L(i,j)/2.
      enddo
    enddo
!!!!!!!!!!!!! i = 1
    do j =2, N2-1
        i = 1
      D_diffusion(1,j)  = - (D11(i+1,j) - D11(i  ,j))/dx1_R(i,j) - (D22(i,j) - D22(i,j-1))/dx2_R(i,j)  
!!!!!!!!!!!!! i = N1
        i = N1
      D_diffusion(N1,j) = - (D11(i  ,j) - D11(i-1,j))/dx1_L(i,j) - (D22(i,j) - D22(i,j-1))/dx2_L(i,j)  
    enddo
!!!!!!!!!!!!! j = 1
    do i =2, N1-1
        j = 1
      D_diffusion(i,1)  = - (D11(i,j) - D11(i-1,j))/dx1_L(i,j) - (D22(i,j+1) - D22(i ,j  ))/dx2_R(i,j)  
!!!!!!!!!!!!! j = N2
        j = N2
      D_diffusion(i,N2) = - (D11(i,j) - D11(i-1,j))/dx1_L(i,j) - (D22(i,j  ) - D22(i,j-1))/dx2_L(i,j)  
    enddo
          
!!!!!!!!!!!!! corners
      D_diffusion( 1, 1) = - (D11(2 ,1 ) - D11(1   ,1 ))/dx1_R(1 ,1 ) - (D22(1 ,2 ) - D22(1 ,1   ))/dx2_R(1 ,1 )  
      D_diffusion( 1,N2) = - (D11(2 ,N2) - D11(1   ,N2))/dx1_R(1 ,N2) - (D22(1 ,N2) - D22(1 ,N2-1))/dx2_L(1 ,N2)  
      D_diffusion(N1, 1) = - (D11(N1,1 ) - D11(N1-1,1 ))/dx1_L(N1,1 ) - (D22(N1, 2) - D22(N1, 1  ))/dx2_R(N1,1 )  
      D_diffusion(N1,N2) = - (D11(N1,N2) - D11(N1-1,N2))/dx1_L(N1,N2) - (D22(N1,N2) - D22(N1,N2-1))/dx2_L(N1,N2)  

    do i = 1, N1-1
      do j = 1, N2-1
         Global_D = Global_D + dx1_R(i,j)*dx2_R(i,j) * &
         (D_diffusion(i,j) + D_diffusion(i+1,j) + D_diffusion(i,j+1) + D_diffusion(i+1,j+1))*0.25_nk
      enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! thermal energy input !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, N1
      do j = 1, N2 
         F_thermal(i,j) = u2%rapp_ad*T%var(i,j)*u2%var(i,j)
      enddo
    enddo
    do i = 1, N1-1
      do j = 1, N2-1
         Global_F = Global_F + dx1_R(i,j)*dx2_R(i,j) * &
                   (F_thermal(i,j) + F_thermal(i+1,j) + F_thermal(i,j+1) + F_thermal(i+1,j+1))*0.25_nk
      enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! viscous dissipation of kinetic energy !!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 2, N1-1
      do j = 2, N2-1
         V_dissipation(i,j) = - u1%nomb_ad * ( CENTRE_2(1._nk, u1%var, i, j, 1)*CENTRE_2(1._nk, u1%var, i, j, 1) &
                                             + CENTRE_2(1._nk, u1%var, i, j, 2)*CENTRE_2(1._nk, u1%var, i, j, 2) &
                                             + CENTRE_2(1._nk, u2%var, i, j, 1)*CENTRE_2(1._nk, u2%var, i, j, 1) &
                                             + CENTRE_2(1._nk, u2%var, i, j, 2)*CENTRE_2(1._nk, u2%var, i, j, 2) )
      enddo
    enddo

    do i = 2, N1-1
      j = 1
         V_dissipation(i,j) = - u1%nomb_ad * ( CENTRE_2(1._nk, u1%var, i, j, 1)*CENTRE_2(1._nk, u1%var, i, j, 1) &
                                             + BACKWIND_2(1._nk, u1%var, i, j, 2)*BACKWIND_2(1._nk, u1%var, i, j, 2) &
                                             + CENTRE_2(1._nk, u2%var, i, j, 1)*CENTRE_2(1._nk, u2%var, i, j, 1) &
                                             + BACKWIND_2(1._nk, u2%var, i, j, 2)*BACKWIND_2(1._nk, u2%var, i, j, 2) )
      j = N2
         V_dissipation(i,j) = - u1%nomb_ad * ( CENTRE_2(1._nk, u1%var, i, j, 1)*CENTRE_2(1._nk, u1%var, i, j, 1) &
                                             + UPWIND_2(1._nk, u1%var, i, j, 2)*UPWIND_2(1._nk, u1%var, i, j, 2) &
                                             + CENTRE_2(1._nk, u2%var, i, j, 1)*CENTRE_2(1._nk, u2%var, i, j, 1) &
                                             + UPWIND_2(1._nk, u2%var, i, j, 2)*UPWIND_2(1._nk, u2%var, i, j, 2) )
    enddo

    do j = 2, N2-1
      i = 1
         V_dissipation(i,j) = - u1%nomb_ad * ( BACKWIND_2(1._nk, u1%var, i, j, 1)*BACKWIND_2(1._nk, u1%var, i, j, 1) &
                                             + CENTRE_2(1._nk, u1%var, i, j, 2)*CENTRE_2(1._nk, u1%var, i, j, 2) &
                                             + BACKWIND_2(1._nk, u2%var, i, j, 1)*BACKWIND_2(1._nk, u2%var, i, j, 1) &
                                             + CENTRE_2(1._nk, u2%var, i, j, 2)*CENTRE_2(1._nk, u2%var, i, j, 2) )
      i = N1
         V_dissipation(i,j) = - u1%nomb_ad * ( UPWIND_2(1._nk, u1%var, i, j, 1)*UPWIND_2(1._nk, u1%var, i, j, 1) &
                                             + CENTRE_2(1._nk, u1%var, i, j, 2)*CENTRE_2(1._nk, u1%var, i, j, 2) &
                                             + UPWIND_2(1._nk, u2%var, i, j, 1)*UPWIND_2(1._nk, u2%var, i, j, 1) &
                                             + CENTRE_2(1._nk, u2%var, i, j, 2)*CENTRE_2(1._nk, u2%var, i, j, 2) )
    enddo

        V_dissipation(1, 1) = - u1%nomb_ad * ( BACKWIND_2(1._nk, u1%var, 1, 1, 1)*BACKWIND_2(1._nk, u1%var, 1, 1, 1) &
                                             + BACKWIND_2(1._nk, u1%var, 1, 1, 2)*BACKWIND_2(1._nk, u1%var, 1, 1, 2) &
                                             + BACKWIND_2(1._nk, u2%var, 1, 1, 1)*BACKWIND_2(1._nk, u2%var, 1, 1, 1) &
                                             + BACKWIND_2(1._nk, u2%var, 1, 1, 2)*BACKWIND_2(1._nk, u2%var, 1, 1, 2) )
        V_dissipation(1,N2) = - u1%nomb_ad * ( BACKWIND_2(1._nk, u1%var, 1,N2, 1)*BACKWIND_2(1._nk, u1%var, 1,N2, 1) &
                                             +   UPWIND_2(1._nk, u1%var, 1,N2, 2)*  UPWIND_2(1._nk, u1%var, 1,N2, 2) &
                                             + BACKWIND_2(1._nk, u2%var, 1,N2, 1)*BACKWIND_2(1._nk, u2%var, 1,N2, 1) &
                                             +   UPWIND_2(1._nk, u2%var, 1,N2, 2)*  UPWIND_2(1._nk, u2%var, 1,N2, 2) )
        V_dissipation(N1,1) = - u1%nomb_ad * (   UPWIND_2(1._nk, u1%var,N1, 1, 1)*  UPWIND_2(1._nk, u1%var,N1, 1, 1) &
                                             + BACKWIND_2(1._nk, u1%var,N1, 1, 2)*BACKWIND_2(1._nk, u1%var,N1, 1, 2) &
                                             +   UPWIND_2(1._nk, u2%var,N1, 1, 1)*  UPWIND_2(1._nk, u2%var,N1, 1, 1) &
                                             + BACKWIND_2(1._nk, u2%var,N1, 1, 2)*BACKWIND_2(1._nk, u2%var,N1, 1, 2) )
       V_dissipation(N1,N2) = - u1%nomb_ad * (   UPWIND_2(1._nk, u1%var,N1,N2, 1)*  UPWIND_2(1._nk, u1%var,N1,N2, 1) &
                                             +   UPWIND_2(1._nk, u1%var,N1,N2, 2)*  UPWIND_2(1._nk, u1%var,N1,N2, 2) &
                                             +   UPWIND_2(1._nk, u2%var,N1,N2, 1)*  UPWIND_2(1._nk, u2%var,N1,N2, 1) &
                                             +   UPWIND_2(1._nk, u2%var,N1,N2, 2)*  UPWIND_2(1._nk, u2%var,N1,N2, 2) )

    do i = 1, N1-1
      do j = 1, N2-1 
         Global_V = Global_V + dx1_R(i,j)*dx2_R(i,j) * &
         (V_dissipation(i,j) + V_dissipation(i+1,j) + V_dissipation(i,j+1) + V_dissipation(i+1,j+1))*0.25_nk
      enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! energy exchange of flow structure and ploymer !!!!!!!!!!!!!!!!!!!
    do i = 1, N1
      do j = 1, N2 
         G_exchange(i,j) =   ( CENTRE_2(u1%var(i,j), t11%var, i, j, 1) +  CENTRE_2(u1%var(i,j), t12%var, i, j, 2) &
                            +  CENTRE_2(u2%var(i,j), t12%var, i, j, 1) +  CENTRE_2(u2%var(i,j), t22%var, i, j, 2) )
      enddo
    enddo

    do i = 2, N1-1
      j = 1
         G_exchange(i,j) =   ( CENTRE_2(u1%var(i,j), t11%var, i, j, 1) +  BACKWIND_2(u1%var(i,j), t12%var, i, j, 2) &
                            +  CENTRE_2(u2%var(i,j), t12%var, i, j, 1) +  BACKWIND_2(u2%var(i,j), t22%var, i, j, 2) )
      j = N2
         G_exchange(i,j) =   ( CENTRE_2(u1%var(i,j), t11%var, i, j, 1) +    UPWIND_2(u1%var(i,j), t12%var, i, j, 2) &
                            +  CENTRE_2(u2%var(i,j), t12%var, i, j, 1) +    UPWIND_2(u2%var(i,j), t22%var, i, j, 2) )
    enddo
    do j = 2, N2-1
      i = 1
         G_exchange(i,j) =   ( BACKWIND_2(u1%var(i,j), t11%var, i, j, 1) +  CENTRE_2(u1%var(i,j), t12%var, i, j, 2) &
                            +  BACKWIND_2(u2%var(i,j), t12%var, i, j, 1) +  CENTRE_2(u2%var(i,j), t22%var, i, j, 2) )
      i = N1
         G_exchange(i,j) =   (   UPWIND_2(u1%var(i,j), t11%var, i, j, 1) +  CENTRE_2(u1%var(i,j), t12%var, i, j, 2) &
                            +    UPWIND_2(u2%var(i,j), t12%var, i, j, 1) +  CENTRE_2(u2%var(i,j), t22%var, i, j, 2) )
    enddo
         G_exchange(1,1) =   ( BACKWIND_2(u1%var(i,j), t11%var, 1, 1, 1) +  BACKWIND_2(u1%var(i,j), t12%var, 1, 1, 2) &
                            +  BACKWIND_2(u2%var(i,j), t12%var, 1, 1, 1) +  BACKWIND_2(u2%var(i,j), t22%var, 1, 1, 2) )
        G_exchange(1,N2) =   ( BACKWIND_2(u1%var(i,j), t11%var, 1,N2, 1) +    UPWIND_2(u1%var(i,j), t12%var, 1,N2, 2) &
                            +  BACKWIND_2(u2%var(i,j), t12%var, 1,N2, 1) +    UPWIND_2(u2%var(i,j), t22%var, 1,N2, 2) )
        G_exchange(N1,1) =   (   UPWIND_2(u1%var(i,j), t11%var,N1, 1, 1) +  BACKWIND_2(u1%var(i,j), t12%var,N1, 1, 2) &
                            +    UPWIND_2(u2%var(i,j), t12%var,N1, 1, 1) +  BACKWIND_2(u2%var(i,j), t22%var,N1, 1, 2) )
       G_exchange(N1,N2) =   (   UPWIND_2(u1%var(i,j), t11%var,N1,N2, 1) +    UPWIND_2(u1%var(i,j), t12%var,N1,N2, 2) &
                            +    UPWIND_2(u2%var(i,j), t12%var,N1,N2, 1) +    UPWIND_2(u2%var(i,j), t22%var,N1,N2, 2) )
    do i = 1, N1-1
      do j = 1, N2-1 
         Global_G = Global_G + dx1_R(i,j)*dx2_R(i,j) * &
         (G_exchange(i,j) + G_exchange(i+1,j) + G_exchange(i,j+1) + G_exchange(i+1,j+1))*0.25_nk
      enddo
    enddo
 
    write(22,*) it, it*dt, Global_E, Diff_E, Global_D, Global_F, Global_V, Global_G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VISCOELASTIC CONSTITUTIVE EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VISCOELASTIC CONSTITUTIVE EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VISCOELASTIC CONSTITUTIVE EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VISCOELASTIC CONSTITUTIVE EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, N1-1
      do j = 1, N2-1 
         Total_normal_tau = Total_normal_tau + dx1_R(i,j)*dx2_R(i,j) * (&
             (t11%var(i,j) + t11%var(i+1,j) + t11%var(i,j+1) + t11%var(i+1,j+1))*0.25_nk &
           + (t22%var(i,j) + t22%var(i+1,j) + t22%var(i,j+1) + t22%var(i+1,j+1))*0.25_nk )
      enddo
    enddo   
    Diff_Tau = ( Total_normal_tau  - Last_Total_Tau ) / dt 
    Last_Total_Tau = Total_normal_tau
 



    do i = N1/2+1, N1-1
      do j = 1, N2-1 
         Total_tau = Total_tau + dx1_R(i,j)*dx2_R(i,j) * (&
             (t11%var(i,j) + t11%var(i+1,j) + t11%var(i,j+1) + t11%var(i+1,j+1))*0.25_nk &
           + (t22%var(i,j) + t22%var(i+1,j) + t22%var(i,j+1) + t22%var(i+1,j+1))*0.25_nk &
           + (t12%var(i,j) + t12%var(i+1,j) + t12%var(i,j+1) + t12%var(i+1,j+1))*0.25_nk )
      enddo
    enddo  
    write(21,*) it, it*dt, Total_normal_tau, Diff_Tau, Total_tau*2

#elif (DIMENSION_GEO ==3)

#endif
   
  end subroutine kinetic_energy_budget

#endif


  SUBROUTINE PRINT_2D(ISTART,JSTART,IFIN,JFIN,NI,NJ,PHI,HEAD)
  INTEGER, INTENT(IN) :: ISTART,JSTART,NI,NJ,IFIN,JFIN
  INTEGER:: ISKIP,JSKIP,ISTA,IEND,I,JJ,J
  REAL(nk), ALLOCATABLE, DIMENSION(:) :: STORE
  REAL(nk), DIMENSION(1:NI,1:NJ), INTENT(IN) :: PHI
  CHARACTER(LEN=20) :: HEAD
!
  ALLOCATE(STORE(ISTART:IFIN))
!
  ISKIP=1; IF(IFIN-ISTART > 80) ISKIP=2
  JSKIP=1; IF(JFIN-JSTART > 80) JSKIP=2
  WRITE(6,110)HEAD
  ISTA=ISTART-12*ISKIP
! 100 CONTINUE
  DO  ! BOUCLE ETERNELLE
  ISTA=ISTA+12*ISKIP
  IEND=ISTA+12*ISKIP-1
  IEND=MIN0(IFIN,IEND)
  WRITE(6,112)
  DO JJ=JSTART,JFIN,JSKIP
     J=JSTART+JFIN-JJ
     DO I=ISTA,IEND
        STORE(I)=PHI(I,J)
        IF(ABS(STORE(I)).LT.1.E-20) STORE(I)=0.
     ENDDO
     WRITE(6,113)  J,STORE(ISTA:IEND:ISKIP)
  ENDDO
  WRITE(6,114) (I,I=ISTA,IEND,ISKIP)
!C------------------------------------------------
!     IF(IEND.LT.IFIN)GO TO 100
  IF(IEND.GE.IFIN) EXIT
  ENDDO ! BOUCLE ETERNELLE
110 FORMAT(" ",20("*-"),7X,A20,7X,20("-*"))
112 FORMAT("  J")
113 FORMAT(" ",I3,12(1PE10.2))
114 FORMAT(" I=   ",I6,11I10)
  DEALLOCATE(STORE)
  END SUBROUTINE PRINT_2D









end module mOutils
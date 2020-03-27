#include "definitions.h"
  

module mViscoelastic_models 
  
  use mBase
  use mVitesse

  implicit none
  
contains
  
  

  
  
  subroutine calcul_contraintes_viscoelatique(tau, d, const1, const2) 
    implicit none
    integer :: i, j, k
    real(nk), intent(in):: const1, const2
#if (DIMENSION_GEO == 2)    
    real(nk),dimension(0:N1+1,0:N2+1),intent(inout):: d
    real(nk),dimension(0:N1+1,0:N2+1),intent(out):: tau
#elif (DIMENSION_GEO == 3)    
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(out):: tau  
#endif 

#if (DIMENSION_GEO == 2) 
    do j= 1, N2
       do i=1, N1    
#if ( EXP_TREAMENT == 0 )      
          tau(i,j) = d(i,j)/(const2 + const1)
#elif ( EXP_TREAMENT == 1 )
          tau(i,j) = d(i,j)/const2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! transform back
!          tau(i,j) = tau(i,j)/Tau_trans
#endif
       end do
    end do    

#elif (DIMENSION_GEO == 3)    
    do k = 1, N3 
       do j= 1, N2          
          do i=1, N1       
             tau(i,j,k) = d(i,j,k)/(const2 + const1)
          end do
       end do
    end do        
#endif 
   
   
                          
  end subroutine calcul_contraintes_viscoelatique

  
  
 


  subroutine coefs_NL_PPT_tr_tau_tau(d, const, t11, t12, t22, indice, t13, t23, t33)
    implicit none
    integer :: i, j, k
    real(nk),intent(in):: const 
    integer,intent(in) :: indice 
    real(nk) ::  const1
#if (DIMENSION_GEO == 2)
    real(nk),dimension(0:N1+1,0:N2+1),intent(inout):: d    
    real(nk), dimension(0:N1+1,0:N2+1), intent(in):: t11, t12, t22
    real(nk), dimension(0:N1+1,0:N2+1), optional :: t13, t23, t33
#elif (DIMENSION_GEO == 3)   
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in):: t11, t12, t22, t13, t23, t33
#endif
         
#if ( EXP_TREAMENT == 0 )
    const1 = const * epsilon_PPT * sqrt(Rayleigh)/(1.-Beta)/Prandtl   
#elif ( EXP_TREAMENT == 1 )
    const1 = const * epsilon_PPT * sqrt(Rayleigh)/(1.-Beta)/Prandtl / Tau_trans  
#endif
#if (DIMENSION_GEO == 2)

    do j= 1, N2
       do i=1, N1          
          d(i,j)= d(i,j) - const1*( t11(i,j) + t22(i,j) ) * var%W(indice)%var(i,j)
       end do
    end do    
   
#elif (DIMENSION_GEO == 3)   
    do k = 1, N3 
       do j= 1, N2          
          do i=1, N1             
             d(i,j,k)= d(i,j,k) - const1*( t11(i,j,k) + t22(i,j,k) + t33(i,j,k) ) * var%W(indice)%var(i,j,k)
          end do
       end do
    end do
#endif

  end subroutine coefs_NL_PPT_tr_tau_tau

  subroutine coefs_NL_PPT_Xi(d, const, u1, u2, t11, t12, t22, indice, u3, t13, t23, t33)
   implicit none
   integer            :: i, j, k
   real(nk),intent(in):: const 
   integer,intent(in) :: indice
   real(nk)           :: const1
#if (DIMENSION_GEO == 2)
   real(nk), dimension(0:N1+1,0:N2+1),intent(inout):: d    
   real(nk), dimension(0:N1+1,0:N2+1), intent(in)  :: u1, u2, t11, t12, t22
   real(nk), dimension(0:N1+1,0:N2+1), optional    :: u3, t13, t23, t33
#elif (DIMENSION_GEO == 3)   
   real(nk), dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
   real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in)  :: u1, u2, u3, t11, t12, t22, t13, t23, t33
#endif
   const1 = const*xi_PPT

#if (DIMENSION_GEO == 2)

   do j= 2, N2-1
      do i=2, N1-1        
       if ( indice == 3 ) then
            d(i,j)= d(i,j) - const1*( CENTRE_2( 2*t11(i,j), u1, i, j, 1) &
                                    + CENTRE_2(   t12(i,j), u1, i, j, 2) &
                                    + CENTRE_2(   t12(i,j), u2, i, j, 1)) 
       endif
       if ( indice == 4 ) then
            d(i,j)= d(i,j) - const1*( CENTRE_2(   t12(i,j), u1, i, j, 1) &
                                    + CENTRE_2(   t12(i,j), u2, i, j, 2) &
                                    + CENTRE_2(.5*t22(i,j), u1, i, j, 2) &
                                    + CENTRE_2(.5*t22(i,j), u2, i, j, 1) &
                                    + CENTRE_2(.5*t11(i,j), u1, i, j, 2) &
                                    + CENTRE_2(.5*t11(i,j), u2, i, j, 1) )
       endif
       if ( indice == 5 ) then
            d(i,j)= d(i,j) - const1*( CENTRE_2( 2*t22(i,j), u2, i, j, 2) &
                                    + CENTRE_2(   t12(i,j), u1, i, j, 2) &
                                    + CENTRE_2(   t12(i,j), u2, i, j, 1) ) 
       endif
      end do
   end do    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  vertical wall
   do j=2, N2-1
       i=1
       if ( indice == 3 ) then
            d(i,j)= d(i,j) - const1*( BACKWIND_2( 2*t11(i,j), u1, i, j, 1) &
                                    + CENTRE_2  (   t12(i,j), u1, i, j, 2) &
                                    + BACKWIND_2(   t12(i,j), u2, i, j, 1) ) 
       endif
       if ( indice == 4 ) then
            d(i,j)= d(i,j) - const1*( BACKWIND_2(   t12(i,j), u1, i, j, 1) &
                                    + CENTRE_2  (   t12(i,j), u2, i, j, 2) &
                                    + CENTRE_2  (.5*t22(i,j), u1, i, j, 2) &
                                    + BACKWIND_2(.5*t22(i,j), u2, i, j, 1) &
                                    + CENTRE_2  (.5*t11(i,j), u1, i, j, 2) &
                                    + BACKWIND_2(.5*t11(i,j), u2, i, j, 1) )
       endif
       if ( indice == 5 ) then
            d(i,j)= d(i,j) - const1*( CENTRE_2  ( 2*t22(i,j), u2, i, j, 2) &
                                    + CENTRE_2  (   t12(i,j), u1, i, j, 2) &
                                    + BACKWIND_2(   t12(i,j), u2, i, j, 1) ) 
       endif

       i=N1
       if ( indice == 3 ) then
            d(i,j)= d(i,j) - const1*( UPWIND_2( 2*t11(i,j), u1, i, j, 1) &
                                    + CENTRE_2(   t12(i,j), u1, i, j, 2) &
                                    + UPWIND_2(   t12(i,j), u2, i, j, 1) ) 
       endif
       if ( indice == 4 ) then
            d(i,j)= d(i,j) - const1*( UPWIND_2(   t12(i,j), u1, i, j, 1) &
                                    + CENTRE_2(   t12(i,j), u2, i, j, 2) &
                                    + CENTRE_2(.5*t22(i,j), u1, i, j, 2) &
                                    + UPWIND_2(.5*t22(i,j), u2, i, j, 1) &
                                    + CENTRE_2(.5*t11(i,j), u1, i, j, 2) &
                                    + UPWIND_2(.5*t11(i,j), u2, i, j, 1) )
       endif
       if ( indice == 5 ) then
            d(i,j)= d(i,j) - const1*( CENTRE_2( 2*t22(i,j), u2, i, j, 2) &
                                    + CENTRE_2(   t12(i,j), u1, i, j, 2) &
                                    + UPWIND_2(   t12(i,j), u2, i, j, 1) ) 
       endif
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  horizontal wall
   do i=2, N1-1
      j=1
      if ( indice == 3 ) then
            d(i,j)= d(i,j) - const1*( CENTRE_2  ( 2*t11(i,j), u1, i, j, 1) &
                                    + BACKWIND_2(   t12(i,j), u1, i, j, 2) &
                                    + CENTRE_2  (   t12(i,j), u2, i, j, 1) ) 
       endif
       if ( indice == 4 ) then
            d(i,j)= d(i,j) - const1*( CENTRE_2  (   t12(i,j), u1, i, j, 1) &
                                    + BACKWIND_2(   t12(i,j), u2, i, j, 2) &
                                    + BACKWIND_2(.5*t22(i,j), u1, i, j, 2) &
                                    + CENTRE_2  (.5*t22(i,j), u2, i, j, 1) &
                                    + BACKWIND_2(.5*t11(i,j), u1, i, j, 2) &
                                    + CENTRE_2  (.5*t11(i,j), u2, i, j, 1) )
       endif
       if ( indice == 5 ) then
            d(i,j)= d(i,j) - const1*( BACKWIND_2( 2*t22(i,j), u2, i, j, 2) &
                                    + BACKWIND_2(   t12(i,j), u1, i, j, 2) &
                                    + CENTRE_2  (   t12(i,j), u2, i, j, 1) ) 
       endif

       j=N2
       if ( indice == 3 ) then
             d(i,j)= d(i,j) - const1*( CENTRE_2( 2*t11(i,j), u1, i, j, 1) &
                                     + UPWIND_2(   t12(i,j), u1, i, j, 2) &
                                     + CENTRE_2(   t12(i,j), u2, i, j, 1) ) 
       endif
       if ( indice == 4 ) then
            d(i,j)= d(i,j) - const1*( CENTRE_2(   t12(i,j), u1, i, j, 1) &
                                    + UPWIND_2(   t12(i,j), u2, i, j, 2) &
                                    + UPWIND_2(.5*t22(i,j), u1, i, j, 2) &
                                    + CENTRE_2(.5*t22(i,j), u2, i, j, 1) &
                                    + UPWIND_2(.5*t11(i,j), u1, i, j, 2) &
                                    + CENTRE_2(.5*t11(i,j), u2, i, j, 1) )
       endif
       if ( indice == 5 ) then
            d(i,j)= d(i,j) - const1*( UPWIND_2( 2*t22(i,j), u2, i, j, 2) &
                                    + UPWIND_2(   t12(i,j), u1, i, j, 2) &
                                    + CENTRE_2(   t12(i,j), u2, i, j, 1) ) 
       endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CORNERS
   if ( indice == 3 ) then
             d(1,1)= d(1,1) - const1*( BACKWIND_2( 2*t11(1,1), u1, 1, 1, 1) &
                                     + BACKWIND_2(   t12(1,1), u1, 1, 1, 2) &
                                     + BACKWIND_2(   t12(1,1), u2, 1, 1, 1) ) 

             d(1,N2)= d(1,N2) - const1*( BACKWIND_2( 2*t11(1,N2), u1, 1, N2, 1) &
                                       + UPWIND_2  (   t12(1,N2), u1, 1, N2, 2) &
                                       + BACKWIND_2(   t12(1,N2), u2, 1, N2, 1) ) 

             d(N1,1)= d(N1,1) - const1*( UPWIND_2  ( 2*t11(N1,1), u1, N1, 1, 1) &
                                       + BACKWIND_2(   t12(N1,1), u1, N1, 1, 2) &
                                       + UPWIND_2  (   t12(N1,1), u2, N1, 1, 1) ) 

             d(N1,N2)= d(N1,N2) - const1*( UPWIND_2( 2*t11(N1,N2), u1, N1, N2, 1) &
                                         + UPWIND_2(   t12(N1,N2), u1, N1, N2, 2) &
                                         + UPWIND_2(   t12(N1,N2), u2, N1, N2, 1) ) 
   endif
   if ( indice == 4 ) then
            d(1,1)= d(1,1) - const1*( BACKWIND_2(   t12(1,1), u1, 1, 1, 1) &
                                    + BACKWIND_2(   t12(1,1), u2, 1, 1, 2) &
                                    + BACKWIND_2(.5*t22(1,1), u1, 1, 1, 2) &
                                    + BACKWIND_2(.5*t22(1,1), u2, 1, 1, 1) &
                                    + BACKWIND_2(.5*t11(1,1), u1, 1, 1, 2) &
                                    + BACKWIND_2(.5*t11(1,1), u2, 1, 1, 1) )

            d(1,N2)= d(1,N2) - const1*( BACKWIND_2(   t12(1,N2), u1, 1, N2, 1) &
                                      + UPWIND_2  (   t12(1,N2), u2, 1, N2, 2) &
                                      + UPWIND_2  (.5*t22(1,N2), u1, 1, N2, 2) &
                                      + BACKWIND_2(.5*t22(1,N2), u2, 1, N2, 1) &
                                      + UPWIND_2  (.5*t11(1,N2), u1, 1, N2, 2) &
                                      + BACKWIND_2(.5*t11(1,N2), u2, 1, N2, 1) )

            d(N1,1)= d(N1,1) - const1*( UPWIND_2  (   t12(N1,1), u1, N1, 1, 1) &
                                     +  BACKWIND_2(   t12(N1,1), u2, N1, 1, 2) &
                                     +  BACKWIND_2(.5*t22(N1,1), u1, N1, 1, 2) &
                                     +  UPWIND_2  (.5*t22(N1,1), u2, N1, 1, 1) &
                                     +  BACKWIND_2(.5*t11(N1,1), u1, N1, 1, 2) &
                                     +  UPWIND_2  (.5*t11(N1,1), u2, N1, 1, 1) )

            d(N1,N2)= d(N1,N2) - const1*( UPWIND_2(   t12(N1,N2), u1, N1, N2, 1) &
                                       +  UPWIND_2(   t12(N1,N2), u2, N1, N2, 2) &
                                       +  UPWIND_2(.5*t22(N1,N2), u1, N1, N2, 2) &
                                       +  UPWIND_2(.5*t22(N1,N2), u2, N1, N2, 1) &
                                       +  UPWIND_2(.5*t11(N1,N2), u1, N1, N2, 2) &
                                       +  UPWIND_2(.5*t11(N1,N2), u2, N1, N2, 1) )
   endif
   if ( indice == 5 ) then
            d(1,1)= d(1,1) - const1*( BACKWIND_2( 2*t22(1,1), u2, 1, 1, 2) &
                                    + BACKWIND_2(   t12(1,1), u1, 1, 1, 2) &
                                    + BACKWIND_2(   t12(1,1), u2, 1, 1, 1) ) 

            d(1,N2)= d(1,N2) - const1*( UPWIND_2  ( 2*t22(1,N2), u2, 1, N2, 2) &
                                      + UPWIND_2  (   t12(1,N2), u1, 1, N2, 2) &
                                      + BACKWIND_2(   t12(1,N2), u2, 1, N2, 1) ) 

            d(N1,1)= d(N1,1) - const1*(  BACKWIND_2( 2*t22(N1,1), u2, N1, 1, 2) &
                                      +  BACKWIND_2(   t12(N1,1), u1, N1, 1, 2) &
                                      +  UPWIND_2  (   t12(N1,1), u2, N1, 1, 1) ) 

            d(N1,N2)= d(N1,N2) - const1*(  UPWIND_2( 2*t22(N1,N2), u2, N1, N2, 2) &
                                        +  UPWIND_2(   t12(N1,N2), u1, N1, N2, 2) &
                                        +  UPWIND_2(   t12(N1,N2), u2, N1, N2, 1) ) 
   endif


#elif (DIMENSION_GEO == 3)   
   do k = 2, N3-1
      do j= 2, N2-1         
         do i=2, N1-1      
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
            endif 
         end do
      end do
   end do

   do k = 2, N3-1
      do j= 2, N2-1         
         i= 1
         if ( indice == 4 ) then
            d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2( 2*t11(i,j,k), u1, i, j, 1, k) &
                                    + BACKWIND_2(   t12(i,j,k), u2, i, j, 1, k) &
                                    + BACKWIND_2(   t13(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
         endif
         if ( indice == 5 ) then
            d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t12(i,j,k), u1, i, j, 1, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                    + BACKWIND_2(.5*t22(i,j,k), u2, i, j, 1, k) &
                                    + BACKWIND_2(.5*t23(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                    + BACKWIND_2(.5*t11(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
         endif
         if ( indice == 6 ) then
            d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t13(i,j,k), u1, i, j, 1, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                    + BACKWIND_2(.5*t23(i,j,k), u2, i, j, 1, k) &
                                    + BACKWIND_2(.5*t33(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                    + BACKWIND_2(.5*t11(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 7 ) then
            d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t12(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 8 ) then
            d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(.5*t13(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                    + BACKWIND_2(.5*t12(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 9 ) then
            d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t13(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
         endif 
         i= N1
         if ( indice == 4 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                    + UPWIND_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
         endif
         if ( indice == 5 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                    + UPWIND_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                    + UPWIND_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
         endif
         if ( indice == 6 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                    + UPWIND_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                    + UPWIND_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 7 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 8 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                    + UPWIND_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 9 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
         endif 
      end do
   end do

   do k = 2, N3-1
      do i= 2, N1-1         
         j= 1
         if ( indice == 4 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                    + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
         endif
         if ( indice == 5 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                    + BACKWIND_2(.5*t22(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(.5*t11(i,j,k), u1, i, j, 2, k) &
                                    + BACKWIND_2(   t12(i,j,k), u2, i, j, 2, k) &
                                    + BACKWIND_2(.5*t13(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
         endif
         if ( indice == 6 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                    + BACKWIND_2(.5*t23(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(.5*t12(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 7 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + BACKWIND_2( 2*t22(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                    + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                    + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 8 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                    + BACKWIND_2(   t23(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                    + BACKWIND_2(.5*t13(i,j,k), u1, i, j, 2, k) &
                                    + BACKWIND_2(.5*t33(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(.5*t22(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 9 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
         endif 
         j= N2
         if ( indice == 4 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                    + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
         endif
         if ( indice == 5 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                    + UPWIND_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                    + UPWIND_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                    + UPWIND_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
         endif
         if ( indice == 6 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                    + UPWIND_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 7 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                    + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                    + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 8 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                    + UPWIND_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                    + UPWIND_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 9 ) then
            d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
         endif 
      end do
   end do

   do i=2, N1-1 
      do j= 2, N2-1         
         k = 1
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t13(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                       + BACKWIND_2( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + BACKWIND_2(   t13(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 3, k) )
            endif 
            k = N3
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                       + UPWIND_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + UPWIND_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 3, k) )
            endif 
      end do
   end do

   do i=2, N1-1 
      j = 1       
         k = 1
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t12(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t13(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) &
                                       + BACKWIND_2( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + BACKWIND_2(   t13(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 3, k) )
            endif 
      j = 1 
         k = N3
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t12(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) &
                                       + UPWIND_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + UPWIND_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 3, k) )
            endif 
      j = N2 
         k = 1
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                       + UPWIND_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t13(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                       + BACKWIND_2( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + BACKWIND_2(   t13(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 3, k) )
            endif 

      j = N2 
         k = N3
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                       + UPWIND_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( CENTRE_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                       + UPWIND_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + UPWIND_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 3, k) )
            endif 
   end do

   do k=2, N3-1 
      i = 1       
         j = 1
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + BACKWIND_2(   t12(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(   t13(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t12(i,j,k), u1, i, j, 1, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t12(i,j,k), u2, i, j, 2, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t13(i,j,k), u1, i, j, 1, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t12(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(   t23(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + BACKWIND_2(.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t13(i,j,k), u3, i, j, 1, k) &
                                       + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
            endif 
      i = 1 
         j = N2
            if ( indice == 4 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2( 2*t11(i,j,k), u1, i, j, 1, k) &
                                       + BACKWIND_2(   t12(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(   t13(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
            endif
            if ( indice == 5 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t12(i,j,k), u1, i, j, 1, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t22(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u2, i, j, 1, k) &
                                       + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
            endif
            if ( indice == 6 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t13(i,j,k), u1, i, j, 1, k) &
                                       + UPWIND_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                       + BACKWIND_2(.5*t23(i,j,k), u2, i, j, 1, k) &
                                       + BACKWIND_2(.5*t33(i,j,k), u3, i, j, 1, k) &
                                       + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                       + BACKWIND_2(.5*t11(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 7 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t12(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 8 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(.5*t13(i,j,k), u2, i, j, 1, k) &
                                       + UPWIND_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                       + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                       + UPWIND_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                       + UPWIND_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                       + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                       + BACKWIND_2(.5*t12(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
            endif
            if ( indice == 9 ) then
               d(i,j,k)= d(i,j,k) - const1*( BACKWIND_2(   t13(i,j,k), u3, i, j, 1, k) &
                                       + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                       + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                       + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                       + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
            endif 
      i = N1 
         j = 1
         if ( indice == 4 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                    + UPWIND_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
         endif
         if ( indice == 5 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                    + BACKWIND_2(.5*t22(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                    + UPWIND_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(.5*t11(i,j,k), u1, i, j, 2, k) &
                                    + BACKWIND_2(   t12(i,j,k), u2, i, j, 2, k) &
                                    + BACKWIND_2(.5*t13(i,j,k), u3, i, j, 2, k) &
                                    + UPWIND_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
         endif
         if ( indice == 6 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                    + BACKWIND_2(.5*t23(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                    + UPWIND_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                    + UPWIND_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(.5*t12(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 7 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + BACKWIND_2( 2*t22(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                    + BACKWIND_2(   t12(i,j,k), u1, i, j, 2, k) &
                                    + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 8 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                    + BACKWIND_2(   t23(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                    + BACKWIND_2(.5*t13(i,j,k), u1, i, j, 2, k) &
                                    + BACKWIND_2(.5*t33(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                    + UPWIND_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(.5*t22(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 9 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + BACKWIND_2(   t23(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
         endif 

      i = N1 
         j = N2
         if ( indice == 4 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  ( 2*t11(i,j,k), u1, i, j, 1, k) &
                                    + UPWIND_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) )
         endif
         if ( indice == 5 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t12(i,j,k), u1, i, j, 1, k) &
                                    + UPWIND_2  (.5*t22(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t23(i,j,k), u1, i, j, 3, k) &
                                    + UPWIND_2  (.5*t22(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (.5*t23(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (.5*t11(i,j,k), u1, i, j, 2, k) &
                                    + UPWIND_2  (   t12(i,j,k), u2, i, j, 2, k) &
                                    + UPWIND_2  (.5*t13(i,j,k), u3, i, j, 2, k) &
                                    + UPWIND_2  (.5*t11(i,j,k), u2, i, j, 1, k) &
                                    + CENTRE_2  (.5*t13(i,j,k), u2, i, j, 3, k) )
         endif
         if ( indice == 6 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t13(i,j,k), u1, i, j, 1, k) &
                                    + UPWIND_2  (.5*t23(i,j,k), u1, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u1, i, j, 3, k) &
                                    + UPWIND_2  (.5*t23(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (.5*t33(i,j,k), u3, i, j, 1, k) &
                                    + CENTRE_2  (.5*t11(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u3, i, j, 3, k) &
                                    + UPWIND_2  (.5*t11(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (.5*t12(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 7 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t12(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  ( 2*t22(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) &
                                    + UPWIND_2  (   t12(i,j,k), u1, i, j, 2, k) &
                                    + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 8 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (.5*t13(i,j,k), u2, i, j, 1, k) &
                                    + UPWIND_2  (   t23(i,j,k), u2, i, j, 2, k) &
                                    + CENTRE_2  (.5*t33(i,j,k), u2, i, j, 3, k) &
                                    + UPWIND_2  (.5*t13(i,j,k), u1, i, j, 2, k) &
                                    + UPWIND_2  (.5*t33(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  (.5*t12(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (.5*t22(i,j,k), u2, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u3, i, j, 3, k) &
                                    + UPWIND_2  (.5*t12(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (.5*t22(i,j,k), u3, i, j, 2, k) )
         endif
         if ( indice == 9 ) then
            d(i,j,k)= d(i,j,k) - const1*( UPWIND_2  (   t13(i,j,k), u3, i, j, 1, k) &
                                    + UPWIND_2  (   t23(i,j,k), u3, i, j, 2, k) &
                                    + CENTRE_2  ( 2*t33(i,j,k), u3, i, j, 3, k) &
                                    + CENTRE_2  (   t13(i,j,k), u1, i, j, 3, k) &
                                    + CENTRE_2  (   t23(i,j,k), u2, i, j, 3, k) )
         endif 
   end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! corners
! (1, 1, 1)
   if ( indice == 4 ) then
      d(1,1,1)= d(1,1,1) - const1*( BACKWIND_2( 2*t11(1,1,1), u1, 1, 1, 1, 1) &
                              + BACKWIND_2(   t12(1,1,1), u2, 1, 1, 1, 1) &
                              + BACKWIND_2(   t13(1,1,1), u3, 1, 1, 1, 1) &
                              + BACKWIND_2(   t12(1,1,1), u1, 1, 1, 2, 1) &
                              + BACKWIND_2(   t13(1,1,1), u1, 1, 1, 3, 1) )
   endif
   if ( indice == 5 ) then
      d(1,1,1)= d(1,1,1) - const1*( BACKWIND_2(   t12(1,1,1), u1, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t22(1,1,1), u1, 1, 1, 2, 1) &
                              + BACKWIND_2(.5*t23(1,1,1), u1, 1, 1, 3, 1) &
                              + BACKWIND_2(.5*t22(1,1,1), u2, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t23(1,1,1), u3, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t11(1,1,1), u1, 1, 1, 2, 1) &
                              + BACKWIND_2(   t12(1,1,1), u2, 1, 1, 2, 1) &
                              + BACKWIND_2(.5*t13(1,1,1), u3, 1, 1, 2, 1) &
                              + BACKWIND_2(.5*t11(1,1,1), u2, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t13(1,1,1), u2, 1, 1, 3, 1) )
   endif
   if ( indice == 6 ) then
      d(1,1,1)= d(1,1,1) - const1*( BACKWIND_2(   t13(1,1,1), u1, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t23(1,1,1), u1, 1, 1, 2, 1) &
                              + BACKWIND_2(.5*t33(1,1,1), u1, 1, 1, 3, 1) &
                              + BACKWIND_2(.5*t23(1,1,1), u2, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t33(1,1,1), u3, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t11(1,1,1), u1, 1, 1, 3, 1) &
                              + BACKWIND_2(.5*t12(1,1,1), u2, 1, 1, 3, 1) &
                              + BACKWIND_2(   t13(1,1,1), u3, 1, 1, 3, 1) &
                              + BACKWIND_2(.5*t11(1,1,1), u3, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t12(1,1,1), u3, 1, 1, 2, 1) )
   endif
   if ( indice == 7 ) then
      d(1,1,1)= d(1,1,1) - const1*( BACKWIND_2(   t12(1,1,1), u2, 1, 1, 1, 1) &
                              + BACKWIND_2( 2*t22(1,1,1), u2, 1, 1, 2, 1) &
                              + BACKWIND_2(   t23(1,1,1), u2, 1, 1, 3, 1) &
                              + BACKWIND_2(   t12(1,1,1), u1, 1, 1, 2, 1) &
                              + BACKWIND_2(   t23(1,1,1), u3, 1, 1, 2, 1) )
   endif
   if ( indice == 8 ) then
      d(1,1,1)= d(1,1,1) - const1*( BACKWIND_2(.5*t13(1,1,1), u2, 1, 1, 1, 1) &
                              + BACKWIND_2(   t23(1,1,1), u2, 1, 1, 2, 1) &
                              + BACKWIND_2(.5*t33(1,1,1), u2, 1, 1, 3, 1) &
                              + BACKWIND_2(.5*t13(1,1,1), u1, 1, 1, 2, 1) &
                              + BACKWIND_2(.5*t33(1,1,1), u3, 1, 1, 2, 1) &
                              + BACKWIND_2(.5*t12(1,1,1), u1, 1, 1, 3, 1) &
                              + BACKWIND_2(.5*t22(1,1,1), u2, 1, 1, 3, 1) &
                              + BACKWIND_2(   t23(1,1,1), u3, 1, 1, 3, 1) &
                              + BACKWIND_2(.5*t12(1,1,1), u3, 1, 1, 1, 1) &
                              + BACKWIND_2(.5*t22(1,1,1), u3, 1, 1, 2, 1) )
   endif
   if ( indice == 9 ) then
      d(1,1,1)= d(1,1,1) - const1*( BACKWIND_2(   t13(1,1,1), u3, 1, 1, 1, 1) &
                              + BACKWIND_2(   t23(1,1,1), u3, 1, 1, 2, 1) &
                              + BACKWIND_2( 2*t33(1,1,1), u3, 1, 1, 3, 1) &
                              + BACKWIND_2(   t13(1,1,1), u1, 1, 1, 3, 1) &
                              + BACKWIND_2(   t23(1,1,1), u2, 1, 1, 3, 1) )
   endif

! (1, N2, 1)
   if ( indice == 4 ) then
      d(1,N2,1)= d(1,N2,1) - const1*( BACKWIND_2( 2*t11(1,N2,1), u1, 1, N2, 1, 1) &
                              + BACKWIND_2(   t12(1,N2,1), u2, 1, N2, 1, 1) &
                              + BACKWIND_2(   t13(1,N2,1), u3, 1, N2, 1, 1) &
                              + UPWIND_2  (   t12(1,N2,1), u1, 1, N2, 2, 1) &
                              + BACKWIND_2(   t13(1,N2,1), u1, 1, N2, 3, 1) )
   endif
   if ( indice == 5 ) then
      d(1,N2,1)= d(1,N2,1) - const1*( BACKWIND_2(   t12(1,N2,1), u1, 1, N2, 1, 1) &
                              + UPWIND_2  (.5*t22(1,N2,1), u1, 1, N2, 2, 1) &
                              + BACKWIND_2(.5*t23(1,N2,1), u1, 1, N2, 3, 1) &
                              + BACKWIND_2(.5*t22(1,N2,1), u2, 1, N2, 1, 1) &
                              + BACKWIND_2(.5*t23(1,N2,1), u3, 1, N2, 1, 1) &
                              + UPWIND_2  (.5*t11(1,N2,1), u1, 1, N2, 2, 1) &
                              + UPWIND_2  (   t12(1,N2,1), u2, 1, N2, 2, 1) &
                              + UPWIND_2  (.5*t13(1,N2,1), u3, 1, N2, 2, 1) &
                              + BACKWIND_2(.5*t11(1,N2,1), u2, 1, N2, 1, 1) &
                              + BACKWIND_2(.5*t13(1,N2,1), u2, 1, N2, 3, 1) )
   endif
   if ( indice == 6 ) then
      d(1,N2,1)= d(1,N2,1) - const1*( BACKWIND_2(   t13(1,N2,1), u1, 1, N2, 1, 1) &
                              + UPWIND_2  (.5*t23(1,N2,1), u1, 1, N2, 2, 1) &
                              + BACKWIND_2(.5*t33(1,N2,1), u1, 1, N2, 3, 1) &
                              + BACKWIND_2(.5*t23(1,N2,1), u2, 1, N2, 1, 1) &
                              + BACKWIND_2(.5*t33(1,N2,1), u3, 1, N2, 1, 1) &
                              + BACKWIND_2(.5*t11(1,N2,1), u1, 1, N2, 3, 1) &
                              + BACKWIND_2(.5*t12(1,N2,1), u2, 1, N2, 3, 1) &
                              + BACKWIND_2(   t13(1,N2,1), u3, 1, N2, 3, 1) &
                              + BACKWIND_2(.5*t11(1,N2,1), u3, 1, N2, 1, 1) &
                              + UPWIND_2  (.5*t12(1,N2,1), u3, 1, N2, 2, 1) )
   endif
   if ( indice == 7 ) then
      d(1,N2,1)= d(1,N2,1) - const1*( BACKWIND_2(   t12(1,N2,1), u2, 1, N2, 1, 1) &
                              + UPWIND_2  ( 2*t22(1,N2,1), u2, 1, N2, 2, 1) &
                              + BACKWIND_2(   t23(1,N2,1), u2, 1, N2, 3, 1) &
                              + UPWIND_2  (   t12(1,N2,1), u1, 1, N2, 2, 1) &
                              + UPWIND_2  (   t23(1,N2,1), u3, 1, N2, 2, 1) )
   endif
   if ( indice == 8 ) then
      d(1,N2,1)= d(1,N2,1) - const1*( BACKWIND_2(.5*t13(1,N2,1), u2, 1, N2, 1, 1) &
                              + UPWIND_2  (   t23(1,N2,1), u2, 1, N2, 2, 1) &
                              + BACKWIND_2(.5*t33(1,N2,1), u2, 1, N2, 3, 1) &
                              + UPWIND_2  (.5*t13(1,N2,1), u1, 1, N2, 2, 1) &
                              + UPWIND_2  (.5*t33(1,N2,1), u3, 1, N2, 2, 1) &
                              + BACKWIND_2(.5*t12(1,N2,1), u1, 1, N2, 3, 1) &
                              + BACKWIND_2(.5*t22(1,N2,1), u2, 1, N2, 3, 1) &
                              + BACKWIND_2(   t23(1,N2,1), u3, 1, N2, 3, 1) &
                              + BACKWIND_2(.5*t12(1,N2,1), u3, 1, N2, 1, 1) &
                              + UPWIND_2  (.5*t22(1,N2,1), u3, 1, N2, 2, 1) )
   endif
   if ( indice == 9 ) then
      d(1,N2,1)= d(1,N2,1) - const1*( BACKWIND_2(   t13(1,N2,1), u3, 1, N2, 1, 1) &
                              + UPWIND_2  (   t23(1,N2,1), u3, 1, N2, 2, 1) &
                              + BACKWIND_2( 2*t33(1,N2,1), u3, 1, N2, 3, 1) &
                              + BACKWIND_2(   t13(1,N2,1), u1, 1, N2, 3, 1) &
                              + BACKWIND_2(   t23(1,N2,1), u2, 1, N2, 3, 1) )
   endif

! (1, 1, N3)
   if ( indice == 4 ) then
      d(1,1,N3)= d(1,1,N3) - const1*( BACKWIND_2( 2*t11(1,1,N3), u1, 1, 1, 1, N3) &
                              + BACKWIND_2(   t12(1,1,N3), u2, 1, 1, 1, N3) &
                              + BACKWIND_2(   t13(1,1,N3), u3, 1, 1, 1, N3) &
                              + BACKWIND_2(   t12(1,1,N3), u1, 1, 1, 2, N3) &
                              + UPWIND_2  (   t13(1,1,N3), u1, 1, 1, 3, N3) )
   endif
   if ( indice == 5 ) then
      d(1,1,N3)= d(1,1,N3) - const1*( BACKWIND_2(   t12(1,1,N3), u1, 1, 1, 1, N3) &
                              + BACKWIND_2(.5*t22(1,1,N3), u1, 1, 1, 2, N3) &
                              + UPWIND_2  (.5*t23(1,1,N3), u1, 1, 1, 3, N3) &
                              + BACKWIND_2(.5*t22(1,1,N3), u2, 1, 1, 1, N3) &
                              + BACKWIND_2(.5*t23(1,1,N3), u3, 1, 1, 1, N3) &
                              + BACKWIND_2(.5*t11(1,1,N3), u1, 1, 1, 2, N3) &
                              + BACKWIND_2(   t12(1,1,N3), u2, 1, 1, 2, N3) &
                              + BACKWIND_2(.5*t13(1,1,N3), u3, 1, 1, 2, N3) &
                              + BACKWIND_2(.5*t11(1,1,N3), u2, 1, 1, 1, N3) &
                              + UPWIND_2  (.5*t13(1,1,N3), u2, 1, 1, 3, N3) )
   endif
   if ( indice == 6 ) then
      d(1,1,N3)= d(1,1,N3) - const1*( BACKWIND_2(   t13(1,1,N3), u1, 1, 1, 1, N3) &
                              + BACKWIND_2(.5*t23(1,1,N3), u1, 1, 1, 2, N3) &
                              + UPWIND_2  (.5*t33(1,1,N3), u1, 1, 1, 3, N3) &
                              + BACKWIND_2(.5*t23(1,1,N3), u2, 1, 1, 1, N3) &
                              + BACKWIND_2(.5*t33(1,1,N3), u3, 1, 1, 1, N3) &
                              + UPWIND_2  (.5*t11(1,1,N3), u1, 1, 1, 3, N3) &
                              + UPWIND_2  (.5*t12(1,1,N3), u2, 1, 1, 3, N3) &
                              + UPWIND_2  (   t13(1,1,N3), u3, 1, 1, 3, N3) &
                              + BACKWIND_2(.5*t11(1,1,N3), u3, 1, 1, 1, N3) &
                              + BACKWIND_2(.5*t12(1,1,N3), u3, 1, 1, 2, N3) )
   endif
   if ( indice == 7 ) then
      d(1,1,N3)= d(1,1,N3) - const1*( BACKWIND_2(   t12(1,1,N3), u2, 1, 1, 1, N3) &
                              + BACKWIND_2( 2*t22(1,1,N3), u2, 1, 1, 2, N3) &
                              + UPWIND_2  (   t23(1,1,N3), u2, 1, 1, 3, N3) &
                              + BACKWIND_2(   t12(1,1,N3), u1, 1, 1, 2, N3) &
                              + BACKWIND_2(   t23(1,1,N3), u3, 1, 1, 2, N3) )
   endif
   if ( indice == 8 ) then
      d(1,1,N3)= d(1,1,N3) - const1*( BACKWIND_2(.5*t13(1,1,N3), u2, 1, 1, 1, N3) &
                              + BACKWIND_2(   t23(1,1,N3), u2, 1, 1, 2, N3) &
                              + UPWIND_2  (.5*t33(1,1,N3), u2, 1, 1, 3, N3) &
                              + BACKWIND_2(.5*t13(1,1,N3), u1, 1, 1, 2, N3) &
                              + BACKWIND_2(.5*t33(1,1,N3), u3, 1, 1, 2, N3) &
                              + UPWIND_2  (.5*t12(1,1,N3), u1, 1, 1, 3, N3) &
                              + UPWIND_2  (.5*t22(1,1,N3), u2, 1, 1, 3, N3) &
                              + UPWIND_2  (   t23(1,1,N3), u3, 1, 1, 3, N3) &
                              + BACKWIND_2(.5*t12(1,1,N3), u3, 1, 1, 1, N3) &
                              + BACKWIND_2(.5*t22(1,1,N3), u3, 1, 1, 2, N3) )
   endif
   if ( indice == 9 ) then
      d(1,1,N3)= d(1,1,N3) - const1*( BACKWIND_2(   t13(1,1,N3), u3, 1, 1, 1, N3) &
                              + BACKWIND_2(   t23(1,1,N3), u3, 1, 1, 2, N3) &
                              + UPWIND_2  ( 2*t33(1,1,N3), u3, 1, 1, 3, N3) &
                              + UPWIND_2  (   t13(1,1,N3), u1, 1, 1, 3, N3) &
                              + UPWIND_2  (   t23(1,1,N3), u2, 1, 1, 3, N3) )
   endif

! (1, N2, N3)
   if ( indice == 4 ) then
      d(1,N2,N3)= d(1,N2,N3) - const1*( BACKWIND_2( 2*t11(1,N2,N3), u1, 1, N2, 1, N3) &
                              + BACKWIND_2(   t12(1,N2,N3), u2, 1, N2, 1, N3) &
                              + BACKWIND_2(   t13(1,N2,N3), u3, 1, N2, 1, N3) &
                              + UPWIND_2  (   t12(1,N2,N3), u1, 1, N2, 2, N3) &
                              + UPWIND_2  (   t13(1,N2,N3), u1, 1, N2, 3, N3) )
   endif
   if ( indice == 5 ) then
      d(1,N2,N3)= d(1,N2,N3) - const1*( BACKWIND_2(   t12(1,N2,N3), u1, 1, N2, 1, N3) &
                              + UPWIND_2  (.5*t22(1,N2,N3), u1, 1, N2, 2, N3) &
                              + UPWIND_2  (.5*t23(1,N2,N3), u1, 1, N2, 3, N3) &
                              + BACKWIND_2(.5*t22(1,N2,N3), u2, 1, N2, 1, N3) &
                              + BACKWIND_2(.5*t23(1,N2,N3), u3, 1, N2, 1, N3) &
                              + UPWIND_2  (.5*t11(1,N2,N3), u1, 1, N2, 2, N3) &
                              + UPWIND_2  (   t12(1,N2,N3), u2, 1, N2, 2, N3) &
                              + UPWIND_2  (.5*t13(1,N2,N3), u3, 1, N2, 2, N3) &
                              + BACKWIND_2(.5*t11(1,N2,N3), u2, 1, N2, 1, N3) &
                              + UPWIND_2  (.5*t13(1,N2,N3), u2, 1, N2, 3, N3) )
   endif
   if ( indice == 6 ) then
      d(1,N2,N3)= d(1,N2,N3) - const1*( BACKWIND_2(   t13(1,N2,N3), u1, 1, N2, 1, N3) &
                              + UPWIND_2  (.5*t23(1,N2,N3), u1, 1, N2, 2, N3) &
                              + UPWIND_2  (.5*t33(1,N2,N3), u1, 1, N2, 3, N3) &
                              + BACKWIND_2(.5*t23(1,N2,N3), u2, 1, N2, 1, N3) &
                              + BACKWIND_2(.5*t33(1,N2,N3), u3, 1, N2, 1, N3) &
                              + UPWIND_2  (.5*t11(1,N2,N3), u1, 1, N2, 3, N3) &
                              + UPWIND_2  (.5*t12(1,N2,N3), u2, 1, N2, 3, N3) &
                              + UPWIND_2  (   t13(1,N2,N3), u3, 1, N2, 3, N3) &
                              + BACKWIND_2(.5*t11(1,N2,N3), u3, 1, N2, 1, N3) &
                              + UPWIND_2  (.5*t12(1,N2,N3), u3, 1, N2, 2, N3) )
   endif
   if ( indice == 7 ) then
      d(1,N2,N3)= d(1,N2,N3) - const1*( BACKWIND_2(   t12(1,N2,N3), u2, 1, N2, 1, N3) &
                              + UPWIND_2  ( 2*t22(1,N2,N3), u2, 1, N2, 2, N3) &
                              + UPWIND_2  (   t23(1,N2,N3), u2, 1, N2, 3, N3) &
                              + UPWIND_2  (   t12(1,N2,N3), u1, 1, N2, 2, N3) &
                              + UPWIND_2  (   t23(1,N2,N3), u3, 1, N2, 2, N3) )
   endif
   if ( indice == 8 ) then
      d(1,N2,N3)= d(1,N2,N3) - const1*( BACKWIND_2(.5*t13(1,N2,N3), u2, 1, N2, 1, N3) &
                              + UPWIND_2  (   t23(1,N2,N3), u2, 1, N2, 2, N3) &
                              + UPWIND_2  (.5*t33(1,N2,N3), u2, 1, N2, 3, N3) &
                              + UPWIND_2  (.5*t13(1,N2,N3), u1, 1, N2, 2, N3) &
                              + UPWIND_2  (.5*t33(1,N2,N3), u3, 1, N2, 2, N3) &
                              + UPWIND_2  (.5*t12(1,N2,N3), u1, 1, N2, 3, N3) &
                              + UPWIND_2  (.5*t22(1,N2,N3), u2, 1, N2, 3, N3) &
                              + UPWIND_2  (   t23(1,N2,N3), u3, 1, N2, 3, N3) &
                              + BACKWIND_2(.5*t12(1,N2,N3), u3, 1, N2, 1, N3) &
                              + UPWIND_2  (.5*t22(1,N2,N3), u3, 1, N2, 2, N3) )
   endif
   if ( indice == 9 ) then
      d(1,N2,N3)= d(1,N2,N3) - const1*( BACKWIND_2(   t13(1,N2,N3), u3, 1, N2, 1, N3) &
                              + UPWIND_2  (   t23(1,N2,N3), u3, 1, N2, 2, N3) &
                              + UPWIND_2  ( 2*t33(1,N2,N3), u3, 1, N2, 3, N3) &
                              + UPWIND_2  (   t13(1,N2,N3), u1, 1, N2, 3, N3) &
                              + UPWIND_2  (   t23(1,N2,N3), u2, 1, N2, 3, N3) )
   endif

! (N1, 1, 1)
   if ( indice == 4 ) then
      d(N1,1,1)= d(N1,1,1) - const1*( UPWIND_2  ( 2*t11(N1,1,1), u1, N1, 1, 1, 1) &
                              + UPWIND_2  (   t12(N1,1,1), u2, N1, 1, 1, 1) &
                              + UPWIND_2  (   t13(N1,1,1), u3, N1, 1, 1, 1) &
                              + BACKWIND_2(   t12(N1,1,1), u1, N1, 1, 2, 1) &
                              + BACKWIND_2(   t13(N1,1,1), u1, N1, 1, 3, 1) )
   endif
   if ( indice == 5 ) then
      d(N1,1,1)= d(N1,1,1) - const1*( UPWIND_2  (   t12(N1,1,1), u1, N1, 1, 1, 1) &
                              + BACKWIND_2(.5*t22(N1,1,1), u1, N1, 1, 2, 1) &
                              + BACKWIND_2(.5*t23(N1,1,1), u1, N1, 1, 3, 1) &
                              + UPWIND_2  (.5*t22(N1,1,1), u2, N1, 1, 1, 1) &
                              + UPWIND_2  (.5*t23(N1,1,1), u3, N1, 1, 1, 1) &
                              + BACKWIND_2(.5*t11(N1,1,1), u1, N1, 1, 2, 1) &
                              + BACKWIND_2(   t12(N1,1,1), u2, N1, 1, 2, 1) &
                              + BACKWIND_2(.5*t13(N1,1,1), u3, N1, 1, 2, 1) &
                              + UPWIND_2  (.5*t11(N1,1,1), u2, N1, 1, 1, 1) &
                              + BACKWIND_2(.5*t13(N1,1,1), u2, N1, 1, 3, 1) )
   endif
   if ( indice == 6 ) then
      d(N1,1,1)= d(N1,1,1) - const1*( UPWIND_2  (   t13(N1,1,1), u1, N1, 1, 1, 1) &
                              + BACKWIND_2(.5*t23(N1,1,1), u1, N1, 1, 2, 1) &
                              + BACKWIND_2(.5*t33(N1,1,1), u1, N1, 1, 3, 1) &
                              + UPWIND_2  (.5*t23(N1,1,1), u2, N1, 1, 1, 1) &
                              + UPWIND_2  (.5*t33(N1,1,1), u3, N1, 1, 1, 1) &
                              + BACKWIND_2(.5*t11(N1,1,1), u1, N1, 1, 3, 1) &
                              + BACKWIND_2(.5*t12(N1,1,1), u2, N1, 1, 3, 1) &
                              + BACKWIND_2(   t13(N1,1,1), u3, N1, 1, 3, 1) &
                              + UPWIND_2  (.5*t11(N1,1,1), u3, N1, 1, 1, 1) &
                              + BACKWIND_2(.5*t12(N1,1,1), u3, N1, 1, 2, 1) )
   endif
   if ( indice == 7 ) then
      d(N1,1,1)= d(N1,1,1) - const1*( UPWIND_2  (   t12(N1,1,1), u2, N1, 1, 1, 1) &
                              + BACKWIND_2( 2*t22(N1,1,1), u2, N1, 1, 2, 1) &
                              + BACKWIND_2(   t23(N1,1,1), u2, N1, 1, 3, 1) &
                              + BACKWIND_2(   t12(N1,1,1), u1, N1, 1, 2, 1) &
                              + BACKWIND_2(   t23(N1,1,1), u3, N1, 1, 2, 1) )
   endif
   if ( indice == 8 ) then
      d(N1,1,1)= d(N1,1,1) - const1*( UPWIND_2  (.5*t13(N1,1,1), u2, N1, 1, 1, 1) &
                              + BACKWIND_2(   t23(N1,1,1), u2, N1, 1, 2, 1) &
                              + BACKWIND_2(.5*t33(N1,1,1), u2, N1, 1, 3, 1) &
                              + BACKWIND_2(.5*t13(N1,1,1), u1, N1, 1, 2, 1) &
                              + BACKWIND_2(.5*t33(N1,1,1), u3, N1, 1, 2, 1) &
                              + BACKWIND_2(.5*t12(N1,1,1), u1, N1, 1, 3, 1) &
                              + BACKWIND_2(.5*t22(N1,1,1), u2, N1, 1, 3, 1) &
                              + BACKWIND_2(   t23(N1,1,1), u3, N1, 1, 3, 1) &
                              + UPWIND_2  (.5*t12(N1,1,1), u3, N1, 1, 1, 1) &
                              + BACKWIND_2(.5*t22(N1,1,1), u3, N1, 1, 2, 1) )
   endif
   if ( indice == 9 ) then
      d(N1,1,1)= d(N1,1,1) - const1*( UPWIND_2  (   t13(N1,1,1), u3, N1, 1, 1, 1) &
                              + BACKWIND_2(   t23(N1,1,1), u3, N1, 1, 2, 1) &
                              + BACKWIND_2( 2*t33(N1,1,1), u3, N1, 1, 3, 1) &
                              + BACKWIND_2(   t13(N1,1,1), u1, N1, 1, 3, 1) &
                              + BACKWIND_2(   t23(N1,1,1), u2, N1, 1, 3, 1) )
   endif

! (N1, 1, N3)
   if ( indice == 4 ) then
      d(N1,1,N3)= d(N1,1,N3) - const1*( UPWIND_2  ( 2*t11(N1,1,N3), u1, N1, 1, 1, N3) &
                              + UPWIND_2  (   t12(N1,1,N3), u2, N1, 1, 1, N3) &
                              + UPWIND_2  (   t13(N1,1,N3), u3, N1, 1, 1, N3) &
                              + BACKWIND_2(   t12(N1,1,N3), u1, N1, 1, 2, N3) &
                              + UPWIND_2  (   t13(N1,1,N3), u1, N1, 1, 3, N3) )
   endif
   if ( indice == 5 ) then
      d(N1,1,N3)= d(N1,1,N3) - const1*( UPWIND_2  (   t12(N1,1,N3), u1, N1, 1, 1, N3) &
                              + BACKWIND_2(.5*t22(N1,1,N3), u1, N1, 1, 2, N3) &
                              + UPWIND_2  (.5*t23(N1,1,N3), u1, N1, 1, 3, N3) &
                              + UPWIND_2  (.5*t22(N1,1,N3), u2, N1, 1, 1, N3) &
                              + UPWIND_2  (.5*t23(N1,1,N3), u3, N1, 1, 1, N3) &
                              + BACKWIND_2(.5*t11(N1,1,N3), u1, N1, 1, 2, N3) &
                              + BACKWIND_2(   t12(N1,1,N3), u2, N1, 1, 2, N3) &
                              + BACKWIND_2(.5*t13(N1,1,N3), u3, N1, 1, 2, N3) &
                              + UPWIND_2  (.5*t11(N1,1,N3), u2, N1, 1, 1, N3) &
                              + UPWIND_2  (.5*t13(N1,1,N3), u2, N1, 1, 3, N3) )
   endif
   if ( indice == 6 ) then
      d(N1,1,N3)= d(N1,1,N3) - const1*( UPWIND_2  (   t13(N1,1,N3), u1, N1, 1, 1, N3) &
                              + BACKWIND_2(.5*t23(N1,1,N3), u1, N1, 1, 2, N3) &
                              + UPWIND_2  (.5*t33(N1,1,N3), u1, N1, 1, 3, N3) &
                              + UPWIND_2  (.5*t23(N1,1,N3), u2, N1, 1, 1, N3) &
                              + UPWIND_2  (.5*t33(N1,1,N3), u3, N1, 1, 1, N3) &
                              + UPWIND_2  (.5*t11(N1,1,N3), u1, N1, 1, 3, N3) &
                              + UPWIND_2  (.5*t12(N1,1,N3), u2, N1, 1, 3, N3) &
                              + UPWIND_2  (   t13(N1,1,N3), u3, N1, 1, 3, N3) &
                              + UPWIND_2  (.5*t11(N1,1,N3), u3, N1, 1, 1, N3) &
                              + BACKWIND_2(.5*t12(N1,1,N3), u3, N1, 1, 2, N3) )
   endif
   if ( indice == 7 ) then
      d(N1,1,N3)= d(N1,1,N3) - const1*( UPWIND_2  (   t12(N1,1,N3), u2, N1, 1, 1, N3) &
                              + BACKWIND_2( 2*t22(N1,1,N3), u2, N1, 1, 2, N3) &
                              + UPWIND_2  (   t23(N1,1,N3), u2, N1, 1, 3, N3) &
                              + BACKWIND_2(   t12(N1,1,N3), u1, N1, 1, 2, N3) &
                              + BACKWIND_2(   t23(N1,1,N3), u3, N1, 1, 2, N3) )
   endif
   if ( indice == 8 ) then
      d(N1,1,N3)= d(N1,1,N3) - const1*( UPWIND_2  (.5*t13(N1,1,N3), u2, N1, 1, 1, N3) &
                              + BACKWIND_2(   t23(N1,1,N3), u2, N1, 1, 2, N3) &
                              +UPWIND_2  (.5*t33(N1,1,N3), u2, N1, 1, 3, N3) &
                              + BACKWIND_2(.5*t13(N1,1,N3), u1, N1, 1, 2, N3) &
                              + BACKWIND_2(.5*t33(N1,1,N3), u3, N1, 1, 2, N3) &
                              + UPWIND_2  (.5*t12(N1,1,N3), u1, N1, 1, 3, N3) &
                              + UPWIND_2  (.5*t22(N1,1,N3), u2, N1, 1, 3, N3) &
                              + UPWIND_2  (   t23(N1,1,N3), u3, N1, 1, 3, N3) &
                              + UPWIND_2  (.5*t12(N1,1,N3), u3, N1, 1, 1, N3) &
                              + BACKWIND_2(.5*t22(N1,1,N3), u3, N1, 1, 2, N3) )
   endif
   if ( indice == 9 ) then
      d(N1,1,N3)= d(N1,1,N3) - const1*( UPWIND_2  (   t13(N1,1,N3), u3, N1, 1, 1, N3) &
                              + BACKWIND_2(   t23(N1,1,N3), u3, N1, 1, 2, N3) &
                              + UPWIND_2  ( 2*t33(N1,1,N3), u3, N1, 1, 3, N3) &
                              + UPWIND_2  (   t13(N1,1,N3), u1, N1, 1, 3, N3) &
                              + UPWIND_2  (   t23(N1,1,N3), u2, N1, 1, 3, N3) )
   endif

! (N1, N2, 1)
   if ( indice == 4 ) then
      d(N1,N2,1)= d(N1,N2,1) - const1*( UPWIND_2  ( 2*t11(N1,N2,1), u1, N1, N2, 1, 1) &
                              + UPWIND_2  (   t12(N1,N2,1), u2, N1, N2, 1, 1) &
                              + UPWIND_2  (   t13(N1,N2,1), u3, N1, N2, 1, 1) &
                              + UPWIND_2  (   t12(N1,N2,1), u1, N1, N2, 2, 1) &
                              + BACKWIND_2(   t13(N1,N2,1), u1, N1, N2, 3, 1) )
   endif
   if ( indice == 5 ) then
      d(N1,N2,1)= d(N1,N2,1) - const1*( UPWIND_2  (   t12(N1,N2,1), u1, N1, N2, 1, 1) &
                              + UPWIND_2  (.5*t22(N1,N2,1), u1, N1, N2, 2, 1) &
                              + BACKWIND_2(.5*t23(N1,N2,1), u1, N1, N2, 3, 1) &
                              + UPWIND_2  (.5*t22(N1,N2,1), u2, N1, N2, 1, 1) &
                              + UPWIND_2  (.5*t23(N1,N2,1), u3, N1, N2, 1, 1) &
                              + UPWIND_2  (.5*t11(N1,N2,1), u1, N1, N2, 2, 1) &
                              + UPWIND_2  (   t12(N1,N2,1), u2, N1, N2, 2, 1) &
                              + UPWIND_2  (.5*t13(N1,N2,1), u3, N1, N2, 2, 1) &
                              + UPWIND_2  (.5*t11(N1,N2,1), u2, N1, N2, 1, 1) &
                              + BACKWIND_2(.5*t13(N1,N2,1), u2, N1, N2, 3, 1) )
   endif
   if ( indice == 6 ) then
      d(N1,N2,1)= d(N1,N2,1) - const1*( UPWIND_2  (   t13(N1,N2,1), u1, N1, N2, 1, 1) &
                              + UPWIND_2  (.5*t23(N1,N2,1), u1, N1, N2, 2, 1) &
                              + BACKWIND_2(.5*t33(N1,N2,1), u1, N1, N2, 3, 1) &
                              + UPWIND_2  (.5*t23(N1,N2,1), u2, N1, N2, 1, 1) &
                              + UPWIND_2  (.5*t33(N1,N2,1), u3, N1, N2, 1, 1) &
                              + BACKWIND_2(.5*t11(N1,N2,1), u1, N1, N2, 3, 1) &
                              + BACKWIND_2(.5*t12(N1,N2,1), u2, N1, N2, 3, 1) &
                              + BACKWIND_2(   t13(N1,N2,1), u3, N1, N2, 3, 1) &
                              + UPWIND_2  (.5*t11(N1,N2,1), u3, N1, N2, 1, 1) &
                              + UPWIND_2  (.5*t12(N1,N2,1), u3, N1, N2, 2, 1) )
   endif
   if ( indice == 7 ) then
      d(N1,N2,1)= d(N1,N2,1) - const1*( UPWIND_2  (   t12(N1,N2,1), u2, N1, N2, 1, 1) &
                              + UPWIND_2  ( 2*t22(N1,N2,1), u2, N1, N2, 2, 1) &
                              + BACKWIND_2(   t23(N1,N2,1), u2, N1, N2, 3, 1) &
                              + UPWIND_2  (   t12(N1,N2,1), u1, N1, N2, 2, 1) &
                              + UPWIND_2  (   t23(N1,N2,1), u3, N1, N2, 2, 1) )
   endif
   if ( indice == 8 ) then
      d(N1,N2,1)= d(N1,N2,1) - const1*( UPWIND_2  (.5*t13(N1,N2,1), u2, N1, N2, 1, 1) &
                              + UPWIND_2  (   t23(N1,N2,1), u2, N1, N2, 2, 1) &
                              + BACKWIND_2(.5*t33(N1,N2,1), u2, N1, N2, 3, 1) &
                              + UPWIND_2  (.5*t13(N1,N2,1), u1, N1, N2, 2, 1) &
                              + UPWIND_2  (.5*t33(N1,N2,1), u3, N1, N2, 2, 1) &
                              + BACKWIND_2(.5*t12(N1,N2,1), u1, N1, N2, 3, 1) &
                              + BACKWIND_2(.5*t22(N1,N2,1), u2, N1, N2, 3, 1) &
                              + BACKWIND_2(   t23(N1,N2,1), u3, N1, N2, 3, 1) &
                              + UPWIND_2  (.5*t12(N1,N2,1), u3, N1, N2, 1, 1) &
                              + UPWIND_2  (.5*t22(N1,N2,1), u3, N1, N2, 2, 1) )
   endif
   if ( indice == 9 ) then
      d(N1,N2,1)= d(N1,N2,1) - const1*( UPWIND_2  (   t13(N1,N2,1), u3, N1, N2, 1, 1) &
                              + UPWIND_2  (   t23(N1,N2,1), u3, N1, N2, 2, 1) &
                              + BACKWIND_2( 2*t33(N1,N2,1), u3, N1, N2, 3, 1) &
                              + BACKWIND_2(   t13(N1,N2,1), u1, N1, N2, 3, 1) &
                              + BACKWIND_2(   t23(N1,N2,1), u2, N1, N2, 3, 1) )
   endif

! (N1, N2, N3)
   if ( indice == 4 ) then
      d(N1,N2,N3)= d(N1,N2,N3) - const1*( UPWIND_2  ( 2*t11(N1,N2,N3), u1, N1, N2, 1, N3) &
                              + UPWIND_2  (   t12(N1,N2,N3), u2, N1, N2, 1, N3) &
                              + UPWIND_2  (   t13(N1,N2,N3), u3, N1, N2, 1, N3) &
                              + UPWIND_2  (   t12(N1,N2,N3), u1, N1, N2, 2, N3) &
                              + UPWIND_2  (   t13(N1,N2,N3), u1, N1, N2, 3, N3) )
   endif
   if ( indice == 5 ) then
      d(N1,N2,N3)= d(N1,N2,N3) - const1*( UPWIND_2  (   t12(N1,N2,N3), u1, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t22(N1,N2,N3), u1, N1, N2, 2, N3) &
                              + UPWIND_2  (.5*t23(N1,N2,N3), u1, N1, N2, 3, N3) &
                              + UPWIND_2  (.5*t22(N1,N2,N3), u2, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t23(N1,N2,N3), u3, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t11(N1,N2,N3), u1, N1, N2, 2, N3) &
                              + UPWIND_2  (   t12(N1,N2,N3), u2, N1, N2, 2, N3) &
                              + UPWIND_2  (.5*t13(N1,N2,N3), u3, N1, N2, 2, N3) &
                              + UPWIND_2  (.5*t11(N1,N2,N3), u2, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t13(N1,N2,N3), u2, N1, N2, 3, N3) )
   endif
   if ( indice == 6 ) then
      d(N1,N2,N3)= d(N1,N2,N3) - const1*( UPWIND_2  (   t13(N1,N2,N3), u1, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t23(N1,N2,N3), u1, N1, N2, 2, N3) &
                              + UPWIND_2  (.5*t33(N1,N2,N3), u1, N1, N2, 3, N3) &
                              + UPWIND_2  (.5*t23(N1,N2,N3), u2, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t33(N1,N2,N3), u3, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t11(N1,N2,N3), u1, N1, N2, 3, N3) &
                              + UPWIND_2  (.5*t12(N1,N2,N3), u2, N1, N2, 3, N3) &
                              + UPWIND_2  (   t13(N1,N2,N3), u3, N1, N2, 3, N3) &
                              + UPWIND_2  (.5*t11(N1,N2,N3), u3, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t12(N1,N2,N3), u3, N1, N2, 2, N3) )
   endif
   if ( indice == 7 ) then
      d(N1,N2,N3)= d(N1,N2,N3) - const1*( UPWIND_2  (   t12(N1,N2,N3), u2, N1, N2, 1, N3) &
                              + UPWIND_2  ( 2*t22(N1,N2,N3), u2, N1, N2, 2, N3) &
                              + UPWIND_2  (   t23(N1,N2,N3), u2, N1, N2, 3, N3) &
                              + UPWIND_2  (   t12(N1,N2,N3), u1, N1, N2, 2, N3) &
                              + UPWIND_2  (   t23(N1,N2,N3), u3, N1, N2, 2, N3) )
   endif
   if ( indice == 8 ) then
      d(N1,N2,N3)= d(N1,N2,N3) - const1*( UPWIND_2  (.5*t13(N1,N2,N3), u2, N1, N2, 1, N3) &
                              + UPWIND_2  (   t23(N1,N2,N3), u2, N1, N2, 2, N3) &
                              + UPWIND_2  (.5*t33(N1,N2,N3), u2, N1, N2, 3, N3) &
                              + UPWIND_2  (.5*t13(N1,N2,N3), u1, N1, N2, 2, N3) &
                              + UPWIND_2  (.5*t33(N1,N2,N3), u3, N1, N2, 2, N3) &
                              + UPWIND_2  (.5*t12(N1,N2,N3), u1, N1, N2, 3, N3) &
                              + UPWIND_2  (.5*t22(N1,N2,N3), u2, N1, N2, 3, N3) &
                              + UPWIND_2  (   t23(N1,N2,N3), u3, N1, N2, 3, N3) &
                              + UPWIND_2  (.5*t12(N1,N2,N3), u3, N1, N2, 1, N3) &
                              + UPWIND_2  (.5*t22(N1,N2,N3), u3, N1, N2, 2, N3) )
   endif
   if ( indice == 9 ) then
      d(N1,N2,N3)= d(N1,N2,N3) - const1*( UPWIND_2  (   t13(N1,N2,N3), u3, N1, N2, 1, N3) &
                              + UPWIND_2  (   t23(N1,N2,N3), u3, N1, N2, 2, N3) &
                              + UPWIND_2  ( 2*t33(N1,N2,N3), u3, N1, N2, 3, N3) &
                              + UPWIND_2  (   t13(N1,N2,N3), u1, N1, N2, 3, N3) &
                              + UPWIND_2  (   t23(N1,N2,N3), u2, N1, N2, 3, N3) )
   endif
#endif

 end subroutine coefs_NL_PPT_Xi
  
  
  
  
  



  subroutine coefs_NL_Giesekus_tau_tau(d, const, t11, t12, t22, indice, t13, t23, t33)
    implicit none
    integer :: i, j, k
    real(nk),intent(in):: const 
    integer,intent(in) :: indice
    real(nk) ::  const1     
#if (DIMENSION_GEO == 2)
    real(nk),dimension(0:N1+1,0:N2+1),intent(inout):: d
    real(nk), dimension(0:N1+1,0:N2+1), intent(in):: t11, t12, t22
    real(nk), dimension(0:N1+1,0:N2+1), optional :: t13, t23, t33
#elif (DIMENSION_GEO == 3)   
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in):: t11, t12, t22, t13, t23, t33
#endif
    const1 = const * alpha_G * sqrt(Rayleigh)/(1-Beta)/Prandtl    

#if (DIMENSION_GEO == 2)
    do j= 1, N2
       do i=1, N1          

         if ( indice==3 )  d(i,j)= d(i,j) - const1*( t11(i,j)**2. + t12(i,j)**2. )
             
         if ( indice==4 )  d(i,j)= d(i,j) - const1*( t11(i,j)*t12(i,j) + t12(i,j)*t22(i,j) )                  
                          
         if ( indice==5 )  d(i,j)= d(i,j) - const1*( t12(i,j)**2. + t22(i,j)**2. )

       end do
    end do    
#elif (DIMENSION_GEO == 3)   
    do k = 1, N3 
       do j= 1, N2          
          do i=1, N1             
                if ( indice==4 )  d(i,j,k)= d(i,j,k) - const1*( t11(i,j,k)**2. + t12(i,j,k)**2. + t13(i,j,k)**2. )
                
                if ( indice==5 )  d(i,j,k)= d(i,j,k) - const1*( t11(i,j,k)*t12(i,j,k) &
                                                              + t12(i,j,k)*t22(i,j,k) + t13(i,j,k)*t23(i,j,k) )
                
                if ( indice==6 )  d(i,j,k)= d(i,j,k) - const1*( t11(i,j,k)*t13(i,j,k) &
                                                              + t12(i,j,k)*t23(i,j,k) + t13(i,j,k)*t33(i,j,k) )
             
                if ( indice==7 )  d(i,j,k)= d(i,j,k) - const1*( t12(i,j,k)**2. + t22(i,j,k)**2. + t23(i,j,k)**2. )
                
                if ( indice==8 )  d(i,j,k)= d(i,j,k) - const1*( t12(i,j,k)*t13(i,j,k) &
                                                              + t22(i,j,k)*t23(i,j,k) + t23(i,j,k)*t33(i,j,k) )
                
                if ( indice==9 )  d(i,j,k)= d(i,j,k) - const1*( t13(i,j,k)**2. + t23(i,j,k)**2. + t33(i,j,k)**2. )
  
          end do
       end do
    end do
#endif

  end subroutine coefs_NL_Giesekus_tau_tau









end module mViscoelastic_models

#include "definitions.h"

module mPression

   use mBase
   use mVitesse

   ! Declarations
   implicit none




contains

  subroutine coefs_U_gradx1_Pres (d ,  P, const)
    implicit none
    integer :: i, j, k
    real(nk),intent(in):: const
#if   (DIMENSION_GEO == 2)
    real(nk),dimension(0:N1+1,0:N2+1),intent(inout)       :: d
    real(nk),dimension(0:N1,0:N2),intent(in   )           :: P
#elif (DIMENSION_GEO == 3)
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
    real(nk),dimension(0:N1,0:N2,0:N3),intent(in   )      :: P
#endif

#if (DIMENSION_GEO == 2)
    do j= 1, N2
       do i=1, N1    
          
#if (INCREMENTALE == 1)   
          d(i,j) =   d(i,j) - const * (  &
          2.*dx2_L(i,j)*( P(i,j  ) - P(i-1,j  ) )/(dx2_L(i,j)+dx2_R(i,j))/(dx1_L(i,j  )+dx1_R(i,j  )) &
        + 2.*dx2_R(i,j)*( P(i,j-1) - P(i-1,j-1) )/(dx2_L(i,j)+dx2_R(i,j))/(dx1_L(i,j-1)+dx1_R(i,j-1)) )
#endif                        
       end do   
    end do  
    
#elif (DIMENSION_GEO == 3)    
   
    do k= 1, N3
       do j= 1, N2
          do i= 1, N1   
      d(i,j,k) =  d(i,j,k) - const * (   &
       ( P(i  ,j  ,k  ) + P(i  ,j-1,k  ) + P(i  ,j  ,k-1) + P(i  ,j-1,k-1) )/4._nk   &
     - ( P(i-1,j  ,k  ) + P(i-1,j-1,k  ) + P(i-1,j  ,k-1) + P(i-1,j-1,k-1) )/4._nk ) &
     *2._nk/(dx1_L(i,j,k) + dx1_R(i,j,k))
!      dx2_L(i,j,k)*dx3_L(i,j,k)/(dx2_L(i,j,k)+dx2_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
!      2.*( P(i,j  ,k  ) - P(i-1,j  ,k  ) )/(dx1_L(i,j  ,k  )+dx1_R(i,j  ,k  )) &
!    + dx2_R(i,j,k)*dx3_L(i,j,k)/(dx2_L(i,j,k)+dx2_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
!    2.*( P(i,j-1,k  ) - P(i-1,j-1,k  ) )/(dx1_L(i,j-1,k  )+dx1_R(i,j-1,k  )) &
!    + dx2_L(i,j,k)*dx3_R(i,j,k)/(dx2_L(i,j,k)+dx2_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
!    2.*( P(i,j  ,k-1) - P(i-1,j  ,k-1) )/(dx1_L(i,j  ,k-1)+dx1_R(i,j  ,k-1)) &
!    + dx2_R(i,j,k)*dx3_R(i,j,k)/(dx2_L(i,j,k)+dx2_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
!    2.*( P(i,j-1,k-1) - P(i-1,j-1,k-1) )/(dx1_L(i,j-1,k-1)+dx1_R(i,j-1,k-1)) )
          end do
       end do
    end do
    
#endif 
  end subroutine coefs_U_gradx1_Pres


  subroutine coefs_U_gradx2_Pres (d ,  P, const)
   implicit none
   integer :: i, j, k
   real(nk),intent(in):: const
#if   (DIMENSION_GEO == 2)
   real(nk),dimension(0:N1+1,0:N2+1),intent(inout)       :: d
   real(nk),dimension(0:N1,0:N2),intent(in   )           :: P
#elif (DIMENSION_GEO == 3)
   real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
   real(nk),dimension(0:N1,0:N2,0:N3),intent(in   )      :: P
#endif

#if (DIMENSION_GEO == 2)

#if (INCREMENTALE == 1)   
    do j= 1, N2
       do i=1, N1 
          d(i,j) =   d(i,j) - const * (  &
          dx1_L(i,j)/(dx1_L(i,j)+dx1_R(i,j))*2.*( P(i  ,j  ) - P(i  ,j-1) )/(dx2_L(i  ,j  )+dx2_R(i  ,j  )) &
        + dx1_R(i,j)/(dx1_L(i,j)+dx1_R(i,j))*2.*( P(i-1,j  ) - P(i-1,j-1) )/(dx2_L(i-1,j  )+dx2_R(i-1,j  )))
       
       end do
    end do
#endif

#elif (DIMENSION_GEO == 3)
  
    do k= 1, N3
       do j= 1, N2
          do i=1, N1   

      d(i,j,k) =  d(i,j,k) - const * (   &
      ( P(i  ,j  ,k  ) + P(i-1,j  ,k  ) + P(i  ,j  ,k-1) + P(i-1,j  ,k-1) )/4._nk   &
    - ( P(i  ,j-1,k  ) + P(i-1,j-1,k  ) + P(i  ,j-1,k-1) + P(i-1,j-1,k-1) )/4._nk ) &
    *2._nk/(dx2_L(i,j,k) + dx2_R(i,j,k))
!      dx1_L(i,j,k)*dx3_L(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
!      2.*( P(i  ,j,k  ) - P(i  ,j-1,k  ) )/(dx2_L(i  ,j,k  )+dx2_R(i  ,j,k  )) &
!    + dx1_R(i,j,k)*dx3_L(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
!    2.*( P(i-1,j,k  ) - P(i-1,j-1,k  ) )/(dx2_L(i-1,j,k  )+dx2_R(i-1,j,k  )) &
!    + dx1_L(i,j,k)*dx3_R(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
!    2.*( P(i  ,j,k-1) - P(i  ,j-1,k-1) )/(dx2_L(i  ,j,k-1)+dx2_R(i  ,j,k-1)) &
!    + dx1_R(i,j,k)*dx3_R(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
!    2.*( P(i-1,j,k-1) - P(i-1,j-1,k-1) )/(dx2_L(i-1,j,k-1)+dx2_R(i-1,j,k-1)) )
          end do
       end do
    end do

#endif
  end subroutine coefs_U_gradx2_Pres

#if (DIMENSION_GEO == 3)
  subroutine coefs_U_gradx3_Pres (d, P, const)
  integer :: i, j, k
  real(nk),intent(in):: const
  real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
  real(nk),dimension(0:N1,0:N2,0:N3),intent(in   ):: P


   do k= 1, N3
      do j= 1, N2
         do i=1, N1   

     d(i,j,k) =  d(i,j,k) - const * (   &
     ( P(i  ,j  ,k  ) + P(i-1,j  ,k  ) + P(i  ,j-1,k  ) + P(i-1,j-1,k  ) )/4._nk   &
   - ( P(i  ,j  ,k-1) + P(i-1,j  ,k-1) + P(i  ,j-1,k-1) + P(i-1,j-1,k-1) )/4._nk ) &
   *2._nk/(dx3_L(i,j,k) + dx3_R(i,j,k))
!     dx1_L(i,j,k)*dx2_L(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx2_L(i,j,k)+dx2_R(i,j,k))*&
!     2.*( P(i  ,j  ,k) - P(i  ,j  ,k-1) )/(dx3_L(i  ,j  ,k)+dx3_R(i  ,j  ,k)) &
!   + dx1_R(i,j,k)*dx2_L(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx2_L(i,j,k)+dx2_R(i,j,k))*&
!   2.*( P(i-1,j  ,k) - P(i-1,j  ,k-1) )/(dx3_L(i-1,j  ,k)+dx3_R(i-1,j  ,k)) &
!   + dx1_L(i,j,k)*dx2_R(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx2_L(i,j,k)+dx2_R(i,j,k))*&
!   2.*( P(i  ,j-1,k) - P(i  ,j-1,k-1) )/(dx3_L(i  ,j-1,k)+dx3_R(i  ,j-1,k)) &
!   + dx1_R(i,j,k)*dx2_R(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx2_L(i,j,k)+dx2_R(i,j,k))*&
!   2.*( P(i-1,j-1,k) - P(i-1,j-1,k-1) )/(dx3_L(i-1,j-1,k)+dx3_R(i-1,j-1,k)) )
         end do
      end do
   end do

 end subroutine coefs_U_gradx3_Pres
#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!( u1%var, u2%var, P%var, P%var0, 3./2., u1%nomb_ad)
subroutine correction_Pression( u1, u2, P, Psi, const1, const2, u3 )
   implicit none
   real(nk),intent(in   ):: const1, const2
   integer :: i, j, k
#if (DIMENSION_GEO == 2)
   real(nk),dimension(0:N1+1,0:N2+1),intent(in)   :: u1, u2
   real(nk),dimension(0:N1  ,0:N2  ),intent(inout):: P, Psi 
   real(nk),intent(in),optional:: u3
#elif (DIMENSION_GEO == 3)
   real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(in)  :: u1, u2, u3
   real(nk),dimension(0:N1  ,0:N2  ,0:N3  ),intent(inout):: P, Psi 
#endif

#if (DIMENSION_GEO == 2)
   do j= 1, N2-1
      do i=1, N1-1
         P(i,j) = P(i,j)  + Psi(i,j) * const1   
      end do
   end do
  
#elif (DIMENSION_GEO == 3)

   do k= 1, N3-1 
      do j= 1, N2-1 
         do i=1, N1-1            
            P(i,j,k) = P(i,j,k)  + Psi(i,j,k) * const1 
         end do
      end do
   end do

#endif 
 end subroutine correction_Pression






  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!(u1%var, u2%var,P%var0, dt) 
  subroutine correction_vitesse(u1, u2, P, const, u3)
    implicit none
    real(nk),intent(in)   :: const
    integer:: i, j, k
#if (DIMENSION_GEO == 2)
    real(nk),dimension(0:N1+1,0:N2+1),intent(inout):: u1, u2
    real(nk),dimension(0:N1  ,0:N2  ),intent(in)   :: P
    real(nk),optional,intent(in)                   :: u3
#elif (DIMENSION_GEO == 3) 
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: u1, u2, u3
    real(nk),dimension(0:N1  ,0:N2  ,0:N3  ),intent(in)   :: P
#endif
    
#if (DIMENSION_GEO == 2)
    do j=1,N2-1
       do i=1, N1
          var_intf_u1(i,j) = var_intf_u1(i,j) - const*2.*( P(i,j) - P(i-1,j) )/(dx1_L(i,j)+dx1_R(i,j))
       end do
    end do
 
    do j=1,N2
       do i=1, N1-1
          var_intf_u2(i,j) = var_intf_u2(i,j) - const*2.*( P(i,j) - P(i,j-1) )/(dx2_L(i,j)+dx2_R(i,j))
       end do 
    end do


    do j= 1, N2
       do i=1, N1
          
          u1(i,j) = u1(i,j)  - const *                                   &
          ( dx2_L(i,j)/(dx2_L(i,j)+dx2_R(i,j))*2.*( P(i,j  ) - P(i-1,j  ) )/(dx1_L(i,j  )+dx1_R(i,j  )) &
          + dx2_R(i,j)/(dx2_L(i,j)+dx2_R(i,j))*2.*( P(i,j-1) - P(i-1,j-1) )/(dx1_L(i,j-1)+dx1_R(i,j-1)))  
         if (i==1  .or. i==N1 .or. j==1  .or. j==N2 ) u1(i,j) = 0.0

          u2(i,j) = u2(i,j) - const *                                   &
          ( dx1_L(i,j)/(dx1_L(i,j)+dx1_R(i,j))*2.*( P(i  ,j) - P(i  ,j-1) )/(dx2_L(i  ,j)+dx2_R(i  ,j)) &
          + dx1_R(i,j)/(dx1_L(i,j)+dx1_R(i,j))*2.*( P(i-1,j) - P(i-1,j-1) )/(dx2_L(i-1,j)+dx2_R(i-1,j)))
         if (i==1  .or. i==N1 .or. j==1  .or. j==N2 ) u2(i,j) = 0.0

                            
       end do
    end do

#elif (DIMENSION_GEO == 3) 
    do k= 1, N3-1
       do j= 1, N2-1
          do i=1, N1 
             var_intf_u1(i,j,k) = var_intf_u1(i,j,k) - const*2.*( P(i,j,k) - P(i-1,j,k) )/(dx1_L(i,j,k)+dx1_R(i,j,k))
          end do
       end do
    end do
    
    
  do k= 1, N3-1
     do j= 1, N2
        do i=1, N1-1 
           var_intf_u2(i,j,k) = var_intf_u2(i,j,k) - const*2.*( P(i,j,k) - P(i,j-1,k) )/(dx2_L(i,j,k)+dx2_R(i,j,k))
        end do     
     end do      
  end do



  do k= 1, N3
     do j= 1, N2-1
        do i=1, N1-1 
           var_intf_u3(i,j,k) = var_intf_u3(i,j,k) - const*2.*( P(i,j,k) - P(i,j,k-1) )/(dx3_L(i,j,k)+dx3_R(i,j,k))
        end do
     end do
  end do


    do k= 1, N3
       do j= 1, N2
          do i=1, N1   

         u1(i,j,k) = u1(i,j,k) - const * (   &
         dx2_L(i,j,k)*dx3_L(i,j,k)/(dx2_L(i,j,k)+dx2_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
         2.*( P(i,j  ,k  ) - P(i-1,j  ,k  ) )/(dx1_L(i,j  ,k  )+dx1_R(i,j  ,k  )) &
       + dx2_R(i,j,k)*dx3_L(i,j,k)/(dx2_L(i,j,k)+dx2_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
       2.*( P(i,j-1,k  ) - P(i-1,j-1,k  ) )/(dx1_L(i,j-1,k  )+dx1_R(i,j-1,k  )) &
       + dx2_L(i,j,k)*dx3_R(i,j,k)/(dx2_L(i,j,k)+dx2_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
       2.*( P(i,j  ,k-1) - P(i-1,j  ,k-1) )/(dx1_L(i,j  ,k-1)+dx1_R(i,j  ,k-1)) &
       + dx2_R(i,j,k)*dx3_R(i,j,k)/(dx2_L(i,j,k)+dx2_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
       2.*( P(i,j-1,k-1) - P(i-1,j-1,k-1) )/(dx1_L(i,j-1,k-1)+dx1_R(i,j-1,k-1)) )

          if (i==1  .or. i==N1 .or. j==1  .or. j==N2 .or. k==1 .or. k==N3) u1(i,j,k) = 0.0
 
            
         u2(i,j,k) = u2(i,j,k) - const * (   &
         dx1_L(i,j,k)*dx3_L(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
         2.*( P(i  ,j,k  ) - P(i  ,j-1,k  ) )/(dx2_L(i  ,j,k  )+dx2_R(i  ,j,k  )) &
       + dx1_R(i,j,k)*dx3_L(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
       2.*( P(i-1,j,k  ) - P(i-1,j-1,k  ) )/(dx2_L(i-1,j,k  )+dx2_R(i-1,j,k  )) &
       + dx1_L(i,j,k)*dx3_R(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
       2.*( P(i  ,j,k-1) - P(i  ,j-1,k-1) )/(dx2_L(i  ,j,k-1)+dx2_R(i  ,j,k-1)) &
       + dx1_R(i,j,k)*dx3_R(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx3_L(i,j,k)+dx3_R(i,j,k))*&
       2.*( P(i-1,j,k-1) - P(i-1,j-1,k-1) )/(dx2_L(i-1,j,k-1)+dx2_R(i-1,j,k-1)) )

          if (i==1  .or. i==N1 .or. j==1  .or. j==N2 .or. k==1 .or. k==N3) u2(i,j,k) = 0.0
             


         u3(i,j,k) = u3(i,j,k) - const * (   &
         dx1_L(i,j,k)*dx2_L(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx2_L(i,j,k)+dx2_R(i,j,k))*&
         2.*( P(i  ,j  ,k) - P(i  ,j  ,k-1) )/(dx3_L(i  ,j  ,k)+dx3_R(i  ,j  ,k)) &
       + dx1_R(i,j,k)*dx2_L(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx2_L(i,j,k)+dx2_R(i,j,k))*&
       2.*( P(i-1,j  ,k) - P(i-1,j  ,k-1) )/(dx3_L(i-1,j  ,k)+dx3_R(i-1,j  ,k)) &
       + dx1_L(i,j,k)*dx2_R(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx2_L(i,j,k)+dx2_R(i,j,k))*&
       2.*( P(i  ,j-1,k) - P(i  ,j-1,k-1) )/(dx3_L(i  ,j-1,k)+dx3_R(i  ,j-1,k)) &
       + dx1_R(i,j,k)*dx2_R(i,j,k)/(dx1_L(i,j,k)+dx1_R(i,j,k))/(dx2_L(i,j,k)+dx2_R(i,j,k))*&
       2.*( P(i-1,j-1,k) - P(i-1,j-1,k-1) )/(dx3_L(i-1,j-1,k)+dx3_R(i-1,j-1,k)) )

          if (i==1  .or. i==N1 .or. j==1  .or. j==N2 .or. k==1 .or. k==N3) u3(i,j,k) = 0.0
       
          end do
       end do
    end do
     
#endif 
  end subroutine correction_vitesse
  





  
  
  
  subroutine soutirage_champ_moyen(P)
    implicit none
#if (DIMENSION_GEO == 2)  
    real(nk),dimension(0:N1,0:N2),intent(inout):: P    
    real(nk) :: somme
    integer :: i,j 
    somme = 0.0

    do j = 1, N2-1
       do i =1, N1-1
          somme = somme +  P(i,j)          
       end do
    end do
    P = P - somme/real((N1-1)*(N2-1),nk)


#elif (DIMENSION_GEO == 3)
    real(nk),dimension(0:N1,0:N2,0:N3),intent(inout):: P    
    real(nk) :: somme
    integer :: i,j,k 
    somme = 0.0

    do k = 1, N3-1
       do j = 1, N2-1
          do i =1, N1-1
             somme = somme +  P(i,j,k)
          end do
       end do 
    end do
    P = P - somme/real((N1-1)*(N2-1)*(N3-1),nk)
#endif 
   end subroutine soutirage_champ_moyen
   
  









  

  
  

 










  subroutine calculs_CL_P ( var)
   implicit none
   type(variable),intent(inout):: var
   integer:: i, j, k
#if (DIMENSION_GEO == 2) 
   
   if ( var%CL_var_gauche == 2 ) then 
         AP(1,:)%CH(3) = AP(1,:)%CH(3) + AP(1,:)%CH(1); AP(1,:)%CH(1) = 0.
   end if 

   if ( var%CL_var_droite == 2 ) then  
         AP(N1-1,:)%CH(3) = AP(N1-1,:)%CH(3) + AP(N1-1,:)%CH(2); AP(N1-1,:)%CH(2) = 0.
   end if

   if ( var%CL_var_bas == 2 ) then
         AP(:,1)%CH(6) = AP(:,1)%CH(6) + AP(:,1)%CH(4); AP(:,1)%CH(4) = 0. 
   end if

   if ( var%CL_var_haut == 2 ) then 
         AP(:,N2-1)%CH(6) = AP(:,N2-1)%CH(6) + AP(:,N2-1)%CH(5); AP(:,N2-1)%CH(5) = 0.
   end if

#elif (DIMENSION_GEO == 3)    

    if ( var%CL_var_droite == 2 ) then 
       AP(:,:,1)%CH(9) = AP(:,:,1)%CH(9) + AP(:,:,1)%CH(7)
       AP(:,:,1)%CH(7) = 0.
    end if  
    
    if ( var%CL_var_gauche == 2 ) then      
       AP(:,:,N3-1)%CH(9) = AP(:,:,N3-1)%CH(9) + AP(:,:,N3-1)%CH(8)    
       AP(:,:,N3-1)%CH(8) = 0.     
    end if
    
    
    if ( var%CL_var_bas == 2 ) then 
       AP(:,1,:)%CH(6) = AP(:,1,:)%CH(6) + AP(:,1,:)%CH(4)      
       AP(:,1,:)%CH(4) = 0.
    end if

     
    if ( var%CL_var_haut == 2 )  then 
       AP(:,N2-1,:)%CH(6) = AP(:,N2-1,:)%CH(6) + AP(:,N2-1,:)%CH(5)     
       AP(:,N2-1,:)%CH(5) = 0.
    end if
    
    
    if ( var%CL_var_avant == 2 )  then 
       AP(1,:,:)%CH(3) = AP(1,:,:)%CH(3) + AP(1,:,:)%CH(1)      
       AP(1,:,:)%CH(1) = 0.          
    end if
    
    
    if ( var%CL_var_arriere == 2 )  then
       AP(N1-1,:,:)%CH(3) = AP(N1-1,:,:)%CH(3) + AP(N1-1,:,:)%CH(2)     
       AP(N1-1,:,:)%CH(2) = 0.
    end if
    
#endif    
  end subroutine calculs_CL_P













  
  subroutine coefs_P ( var )
   implicit none
   type(variable),intent(inout):: var
   integer:: i, j, k
   
#if (DIMENSION_GEO == 2)
   do j=1,N2-1
      do i=1, N1-1
         AP(i,j)%CH(1) =  2._nk / (.5_nk*(dx1_L(i  ,j)+dx1_R(i  ,j))*.5_nk*(dx1_L(i,j)+dx1_R(i,j)+dx1_L(i+1,j)+dx1_R(i+1,j)))
         AP(i,j)%CH(2) =  2._nk / (.5_nk*(dx1_L(i+1,j)+dx1_R(i+1,j))*.5_nk*(dx1_L(i,j)+dx1_R(i,j)+dx1_L(i+1,j)+dx1_R(i+1,j)))
         AP(i,j)%CH(3) = -2._nk / (.5_nk*(dx1_L(i  ,j)+dx1_R(i  ,j))*.5_nk*(dx1_L(i+1,j)+dx1_R(i+1,j)))

         AP(i,j)%CH(4) = -2._nk / (.5_nk*(dx2_L(i,j  )+dx2_R(i,j  ))*.5_nk*(dx2_L(i,j)+dx2_R(i,j)+dx2_L(i,j+1)+dx2_R(i,j+1)))
         AP(i,j)%CH(5) = -2._nk / (.5_nk*(dx2_L(i,j+1)+dx2_R(i,j+1))*.5_nk*(dx2_L(i,j)+dx2_R(i,j)+dx2_L(i,j+1)+dx2_R(i,j+1)))
         AP(i,j)%CH(6) = -2._nk / (.5_nk*(dx2_L(i,j  )+dx2_R(i,j  ))*.5_nk*(dx2_L(i,j+1)+dx2_R(i,j+1)))
      var%d(i,j) =  0._nk 
      end do
   end do
   
#elif (DIMENSION_GEO == 3)
   do k=1,N3-1
      do j=1,N2-1
         do i=1, N1-1
            AP(i,j,k)%CH(1) =  2._nk / (.5_nk*(dx1_L(i  ,j,k)+dx1_R(i  ,j,k))*&
                                        .5_nk*(dx1_L(i  ,j,k)+dx1_R(i,j,k)+dx1_L(i+1,j,k)+dx1_R(i+1,j,k)))
            AP(i,j,k)%CH(2) =  2._nk / (.5_nk*(dx1_L(i+1,j,k)+dx1_R(i+1,j,k))*&
                                        .5_nk*(dx1_L(i  ,j,k)+dx1_R(i,j,k)+dx1_L(i+1,j,k)+dx1_R(i+1,j,k)))
            AP(i,j,k)%CH(3) = -2._nk / (.5_nk*(dx1_L(i  ,j,k)+dx1_R(i  ,j,k))*.5_nk*(dx1_L(i+1,j,k)+dx1_R(i+1,j,k)))

            AP(i,j,k)%CH(4) = -2._nk / (.5_nk*(dx2_L(i,j  ,k)+dx2_R(i,j  ,k))*&
                                        .5_nk*(dx2_L(i,j  ,k)+dx2_R(i,j  ,k)+dx2_L(i,j+1,k)+dx2_R(i,j+1,k)))
            AP(i,j,k)%CH(5) = -2._nk / (.5_nk*(dx2_L(i,j+1,k)+dx2_R(i,j+1,k))*&
                                        .5_nk*(dx2_L(i,j  ,k)+dx2_R(i,j  ,k)+dx2_L(i,j+1,k)+dx2_R(i,j+1,k)))
            AP(i,j,k)%CH(6) = -2._nk / (.5_nk*(dx2_L(i,j  ,k)+dx2_R(i,j  ,k))*.5_nk*(dx2_L(i,j+1,k)+dx2_R(i,j+1,k)))

            AP(i,j,k)%CH(7) =  2._nk / (.5_nk*(dx3_L(i,j,k  )+dx3_R(i,j,k  ))*&
                                        .5_nk*(dx3_L(i,j,k  )+dx3_R(i,j,k  )+dx3_L(i,j,k+1)+dx3_R(i,j,k+1)))
            AP(i,j,k)%CH(8) =  2._nk / (.5_nk*(dx3_L(i,j,k+1)+dx3_R(i,j,k+1))*&
                                        .5_nk*(dx3_L(i,j,k  )+dx3_R(i,j,k  )+dx3_L(i,j,k+1)+dx3_R(i,j,k+1)))
            AP(i,j,k)%CH(9) = -2._nk / (.5_nk*(dx3_L(i,j,k  )+dx3_R(i,j,k  ))*.5_nk*(dx3_L(i,j,k+1)+dx3_R(i,j,k+1)))
         var%d(i,j,k) =  0._nk 
         end do
      end do 
   end do
  
#endif
 end subroutine coefs_P











  
  subroutine source_P (d, u1, u2, const, u3)
    implicit none
    real(nk),intent(in) :: const
    integer:: i, j, k
#if (DIMENSION_GEO == 2) 
    real(nk),dimension(0:N1+1,0:N2+1),intent(in) :: u1, u2
    real(nk),dimension(1:N1-1,1:N2-1),intent(out):: d 
    real(nk),optional :: u3
#elif (DIMENSION_GEO == 3)
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(in) :: u1, u2, u3
    real(nk),dimension(1:N1-1,1:N2-1,1:N3-1),intent(out):: d        
#endif
     

#if (DIMENSION_GEO == 2 )
    do j=1,N2-1
       do i=1,N1   
          var_intf_u1(i,j) = ( u1(i,j) + u1(i,j+1) )/2._nk
       end do
    end do

    do j=1,N2
       do i=1,N1-1
          var_intf_u2(i,j) = ( u2(i,j) + u2(i+1,j) )/2._nk
       end do
    end do

    do j= 1, N2-1
       do i=1, N1-1
         d(i,j) = 1._nk/const*                                         (  &
               (var_intf_u1(i+1,j) - var_intf_u1(i,j) ) / dx1_R(i  ,j  )  &
             + (var_intf_u2(i,j+1) - var_intf_u2(i,j) ) / dx2_R(i  ,j  )  )    
     end do 
    end do

#elif (DIMENSION_GEO == 3)

  do k= 1, N3-1
     do j= 1, N2-1
        do i=1, N1 
           var_intf_u1(i,j,k) =                    (&
                u1(i  ,j  ,k  ) + u1(i  ,j+1,k  ) + &
                u1(i  ,j  ,k+1) + u1(i  ,j+1,k+1) )/4._nk        
        end do
     end do
  end do


  do k= 1, N3-1
     do j= 1, N2
        do i=1, N1-1 
           var_intf_u2(i,j,k) =                    (&
                u2(i  ,j  ,k  ) + u2(i+1,j  ,k  ) + &
                u2(i  ,j  ,k+1) + u2(i+1,j  ,k+1) )/4._nk            
        end do     
     end do      
  end do

  do k= 1, N3
     do j= 1, N2-1
        do i=1, N1-1 
           var_intf_u3(i,j,k) =                    (&
                u3(i  ,j  ,k  ) + u3(i+1,j  ,k  ) + &
                u3(i  ,j+1,k  ) + u3(i+1,j+1,k  ) )/4._nk                         
        end do
     end do
  end do
  
  do k= 1, N3-1
     do j= 1, N2-1
        do i=1, N1-1
           d(i,j,k) = 1._nk/const *  ( & 
           (var_intf_u1(i+1,j  ,k  ) - var_intf_u1(i,j,k) ) / dx1_R(i  ,j  ,k  )  &
         + (var_intf_u2(i  ,j+1,k  ) - var_intf_u2(i,j,k) ) / dx2_R(i  ,j  ,k  )  &
         + (var_intf_u3(i  ,j  ,k+1) - var_intf_u3(i,j,k) ) / dx3_R(i  ,j  ,k  )  )                      
        end do
     end do
  end do

#endif 
  end subroutine source_P














 













  

  


  
  subroutine P_points_ficts (var, indice)
    implicit none
    
    type(variable),intent(inout):: var
    integer,optional,intent(in) :: indice
    integer:: i, j, k
    

#if (DIMENSION_GEO == 2)

       if (present(indice)) then !!! \Phi : (P^n+1 - P^n)            
          if ( var%CL_var_bas == 2 ) var%var0(1:N1-1,0) =   var%var0(1:N1-1,1) 
       else                      !!! Pression 
          if ( var%CL_var_bas == 2 ) var%var (1:N1-1,0) =   var%var (1:N1-1,1) - dx2_R(1,1) * var%q_var_bas          
       end if
          
       if (present(indice)) then
          if ( var%CL_var_haut == 2 ) var%var0(1:N1-1,N2) =   var%var0(1:N1-1,N2-1)
       else
          if ( var%CL_var_haut == 2 ) var%var (1:N1-1,N2) =   var%var (1:N1-1,N2-1) + dx2_L(1,N2) * var%q_var_haut          
       end if

       if (present(indice)) then
          if ( var%CL_var_gauche == 2 ) var%var0(0,1:N2-1) =   var%var0(1,1:N2-1) 
       else
          if ( var%CL_var_gauche == 2 ) var%var (0,1:N2-1) =   var%var (1,1:N2-1) - dx1_R(1,1) * var%q_var_gauche           
       end if
  
       if (present(indice)) then
          if ( var%CL_var_droite == 2 ) var%var0(N1,1:N2-1) =   var%var0(N1-1,1:N2-1) 
       else 
          if ( var%CL_var_droite == 2 ) var%var (N1,1:N2-1) =   var%var (N1-1,1:N2-1) + dx1_L(N1,1) * var%q_var_droite
       end if  
       

 
!!!----points sur les coins
!!!  bas-gauche
       if (present(indice)) then 
          var%var0(0,0) =   var%var0(1,1)
       else
          var%var(0,0)  =   var%var(1,1)
       end if
!!! bas-droit   
       if (present(indice)) then 
          var%var0(N1,0) =   var%var0(N1-1,1)
       else
          var%var(N1,0)  =   var%var(N1-1,1)
       end if         
!!! haut-droit 
       if (present(indice)) then 
          var%var0(N1,N2) =   var%var0(N1-1,N2-1)
       else
          var%var(N1,N2)  =   var%var(N1-1,N2-1)
       end if
!!!  haut-gauche
       if (present(indice)) then 
          var%var0(0,N2) =   var%var0(1,N2-1)
       else
          var%var(0,N2)  =   var%var(1,N2-1)
       end if

 
#elif (DIMENSION_GEO == 3)    
    
!!!  bas
    if (present(indice)) then
       if ( var%CL_var_bas == 2 ) var%var0(1:N1-1,0,1:N3-1) = var%var0(1:N1-1,1,1:N3-1)
    else
       if ( var%CL_var_bas == 2 ) var%var (1:N1-1,0,1:N3-1) = var%var (1:N1-1,1,1:N3-1) - dx2_R(1:N1-1,1,1:N3-1)*var%q_var_bas
    end if
!!!  haut      
    if (present(indice)) then
       if ( var%CL_var_haut == 2 ) var%var0(1:N1-1,N2,1:N3-1) = var%var0(1:N1-1,N2-1,1:N3-1)
    else
       if ( var%CL_var_haut == 2 ) var%var (1:N1-1,N2,1:N3-1) = var%var (1:N1-1,N2-1,1:N3-1) &
                                                              + dx2_L(1:N1-1,N2-1,1:N3-1)*var%q_var_haut  
    end if          


  !!!  droite           
    if (present(indice)) then !!! \Phi : (P^n+1 - P^n)
       if ( var%CL_var_droite == 2 ) var%var0(1:N1-1,1:N2-1,0) = var%var0(1:N1-1,1:N2-1,1)
    else                      !!! Pression 
       if ( var%CL_var_droite == 2 ) var%var (1:N1-1,1:N2-1,0) = var%var(1:N1-1,1:N2-1,1) - dx3_R(1:N1-1,1:N2-1,1)*var%q_var_droite
    end if
  !!!  gauche
    if (present(indice)) then
       if ( var%CL_var_gauche == 2 ) var%var0(1:N1-1,1:N2-1,N3) = var%var0(1:N1-1,1:N2-1,N3-1)
    else
       if ( var%CL_var_gauche == 2 ) var%var (1:N1-1,1:N2-1,N3) = var%var(1:N1-1,1:N2-1,N3-1) &
                                                                + dx3_L(1:N1-1,1:N2-1,N3-1)*var%q_var_gauche
    end if

!!!  avant       
    if (present(indice)) then
       if ( var%CL_var_avant == 2 ) var%var0(0,1:N2-1,1:N3-1) = var%var0(1,1:N2-1,1:N3-1)
    else 
       if ( var%CL_var_avant == 2 ) var%var (0,1:N2-1,1:N3-1) = var%var (1,1:N2-1,1:N3-1) - dx1_R(1,1:N2-1,1:N3-1)*var%q_var_avant  
    end if
!!!  arriere
    if (present(indice)) then  
       if ( var%CL_var_arriere == 2 ) var%var0(N1,1:N2-1,1:N3-1) = var%var0(N1-1,1:N2-1,1:N3-1)
    else
       if ( var%CL_var_arriere == 2 ) var%var (N1,1:N2-1,1:N3-1) = var%var (N1-1,1:N2-1,1:N3-1) &
                                                                 + dx1_L(N1-1,1:N2-1,1:N3-1)*var%q_var_arriere
    end if
     

    if (present(indice)) then 
      var%var0(1:N1-1, 0, 0) =   var%var0(1:N1-1,   1,   1)
      var%var0(1:N1-1,N2, 0) =   var%var0(1:N1-1,N2-1,   1)
      var%var0(1:N1-1, 0,N3) =   var%var0(1:N1-1,   1,N3-1)
      var%var0(1:N1-1,N2,N3) =   var%var0(1:N1-1,N2-1,N3-1)

      var%var0( 0, 1:N2-1, 0) =   var%var0(   1, 1:N2-1,   1)
      var%var0(N1, 1:N2-1, 0) =   var%var0(N1-1, 1:N2-1,   1)
      var%var0( 0, 1:N2-1,N3) =   var%var0(   1, 1:N2-1,N3-1)
      var%var0(N1, 1:N2-1,N3) =   var%var0(N1-1, 1:N2-1,N3-1)

      var%var0( 0, 0,1:N3-1) =   var%var0(   1,   1,1:N3-1)
      var%var0(N1, 0,1:N3-1) =   var%var0(N1-1,   1,1:N3-1)
      var%var0( 0,N2,1:N3-1) =   var%var0(   1,N2-1,1:N3-1)
      var%var0(N1,N2,1:N3-1) =   var%var0(N1-1,N2-1,1:N3-1)
    else
      var%var (1:N1-1, 0, 0) =   var%var (1:N1-1,   1,   1)
      var%var (1:N1-1, 0,N3) =   var%var (1:N1-1,   1,N3-1)
      var%var (1:N1-1,N2, 0) =   var%var (1:N1-1,N2-1,   1)
      var%var (1:N1-1,N2,N3) =   var%var (1:N1-1,N2-1,N3-1)

      var%var ( 0, 1:N2-1, 0) =   var%var (   1, 1:N2-1,   1)
      var%var (N1, 1:N2-1, 0) =   var%var (N1-1, 1:N2-1,   1)
      var%var ( 0, 1:N2-1,N3) =   var%var (   1, 1:N2-1,N3-1)
      var%var (N1, 1:N2-1,N3) =   var%var (N1-1, 1:N2-1,N3-1)

      var%var ( 0, 0,1:N3-1) =   var%var (   1,   1,1:N3-1)
      var%var (N1, 0,1:N3-1) =   var%var (N1-1,   1,1:N3-1)
      var%var ( 0,N2,1:N3-1) =   var%var (   1,N2-1,1:N3-1)
      var%var (N1,N2,1:N3-1) =   var%var (N1-1,N2-1,1:N3-1)
    end if



!!!----points sur les coins
    if (present(indice)) then 
      var%var0(0, 0 ,0 ) =   var%var0(   1,   1,   1)
      var%var0(N1,0 ,0 ) =   var%var0(N1-1,   1,   1)
      var%var0(0 ,N2,0 ) =   var%var0(   1,N2-1,   1)
      var%var0(N1,N2,0 ) =   var%var0(N1-1,N2-1,   1)

      var%var0(0, 0 ,N3) =   var%var0(   1,   1,N3-1)
      var%var0(N1,0 ,N3) =   var%var0(N1-1,   1,N3-1)
      var%var0(0 ,N2,N3) =   var%var0(   1,N2-1,N3-1)
      var%var0(N1,N2,N3) =   var%var0(N1-1,N2-1,N3-1)
    else
      var%var (0, 0 ,0 ) =   var%var (   1,   1,   1)
      var%var (N1,0 ,0 ) =   var%var (N1-1,   1,   1)
      var%var (0 ,N2,0 ) =   var%var (   1,N2-1,   1)
      var%var (N1,N2,0 ) =   var%var (N1-1,N2-1,   1)

      var%var (0, 0 ,N3) =   var%var (   1,   1,N3-1)
      var%var (N1,0 ,N3) =   var%var (N1-1,   1,N3-1)
      var%var (0 ,N2,N3) =   var%var (   1,N2-1,N3-1)
      var%var (N1,N2,N3) =   var%var (N1-1,N2-1,N3-1)
    end if
     
#endif
    
  end subroutine P_points_ficts








  SUBROUTINE EIGP(NI,NVP,EIGX,VPX,VPX1,AP_1,AP_2,AP_3,indice)
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)                               :: NI, NVP
  REAL(nk), DIMENSION(0:NI+1), INTENT(INOUT)        :: EIGX
  REAL(nk), DIMENSION(0:NI+1,0:NI+1), INTENT(INOUT) :: VPX, VPX1
  REAL(nk), DIMENSION(0:NI+1,0:NI+1)                :: WOKE
#if (DIMENSION_GEO == 2)
  REAL(nk), DIMENSION(0:N1,0:N2),     INTENT(IN)    :: AP_1, AP_2, AP_3
#elif (DIMENSION_GEO ==3)
  REAL(nk), DIMENSION(0:N1,0:N2,0:N3),INTENT(IN)    :: AP_1, AP_2, AP_3
#endif             
  REAL(nk), ALLOCATABLE, DIMENSION(:,:)             :: WOKE3, WOKE4, WOKE5
  REAL(nk), ALLOCATABLE, DIMENSION(:)               :: WOKE6
  INTEGER                                           :: I, J, NB, IFAIL, ILAENV
  INTEGER, OPTIONAL                                 :: indice

  !     INTERFACE
  !        SUBROUTINE MXM(A,NAR,B,NBR,C,NCC,ALPHA,BETA)
  !        INTEGER(KIND=4) :: NAR, NBR, NCC ! STOPPING ICI
  !        REAL(KIND=8), OPTIONAL, INTENT(IN) :: ALPHA, BETA
  !        REAL(KIND=8), DIMENSION(NAR,NBR), INTENT(IN) :: A
  !        REAL(KIND=8), DIMENSION(NBR,NCC), INTENT(IN) :: B
  !        REAL(KIND=8), DIMENSION(NAR,NCC), INTENT(INOUT) :: C
  !        END SUBROUTINE MXM
  !     END INTERFACE
  EIGX=0; VPX=0; VPX1=0; WOKE=0.
!  WOKE3(:,:)=0; WOKE4(:,:)=0; WOKE5(:,:)=0; WOKE6(:)=0
#if(DIMENSION_GEO ==2)
  DO i = 1 , NI-1
  WOKE(i  ,i+1) = AP_2(i  ,2)
  WOKE(i+1,i+1) = AP_3(i+1,2)
  WOKE(i+1,i  ) = AP_1(i+1,2)
  ENDDO
  WOKE(  1,  1) = AP_3(  1,2)
#elif(DIMENSION_GEO ==3)
  if ( indice == 1 ) then
   DO i = 1 , NI-1
      WOKE(i  ,i+1) = AP_2(i  ,2,2)
      WOKE(i+1,i+1) = AP_3(i+1,2,2)
      WOKE(i+1,i  ) = AP_1(i+1,2,2)
   ENDDO
      WOKE(1  ,1  ) = AP_3(1  ,2,2)
  else if ( indice == 3 ) then
   DO i = 1 , NI-1
      WOKE(i  ,i+1) = AP_2(2,2,i)
      WOKE(i+1,i+1) = AP_3(2,2,i+1)
      WOKE(i+1,i  ) = AP_1(2,2,i+1)
   ENDDO
      WOKE(1  ,1  ) = AP_3(2,2,  1)
  endif
#endif


  PRINT *,'NVP =',NVP
  !        CALL PRINT_2D(1,1,NI+2,NI+2,NI+2,NI+2,WOKE(0,0) ,HEADP)
  NB=ILAENV(1,'DGETRI',' ',NI+1,-1,-1,-1)
  IF(NB < 4) NB=4; NB = (NI+1)*(NB+1)
  ALLOCATE( WOKE6(1:NB) ); WOKE6 = 0.
  !
  !====67==1=========2=========3=========4=========5=========6=========7====
  !
  ALLOCATE(WOKE3(0:NI+1,0:NI+1), WOKE4(0:NI+1,0:NI+1), WOKE5(0:NI+1,0:NI+1))
  WOKE3 = 0.; WOKE4 = 0.; WOKE5 = 0.;
  IFAIL=0
  CALL DGEEV('N','V',NVP,WOKE(1:NVP,1:NVP),NVP,EIGX(1),WOKE3(1,0), &
             WOKE4(1:NVP,1:NVP),NVP,WOKE5(1:NVP,1:NVP),NVP,WOKE6,NB,IFAIL)
  DO J = 1 , NVP
     IF (WOKE3(J,0) /= 0.) THEN
       PRINT *,J,WOKE3(J,0)
       PRINT *,'VALEURS PROPRES COMPLEXES DANS S-PROG EIGGLP:4'
       STOP; ENDIF
  ENDDO
  !
  IF(IFAIL /= 0)THEN
     PRINT *,' CODE DE RETOUR DE DGEEV  ',IFAIL
     STOP;  ENDIF
  !
  VPX = WOKE5
  !
  CALL DLACPY('A',NVP,NVP,WOKE5(1:NVP,1:NVP),NVP,WOKE4(1:NVP,1:NVP),NVP)
  CALL DGETRF(NVP,NVP,WOKE4(1:NVP,1:NVP),NVP,WOKE3(1:NVP,1:NVP),IFAIL)
  !
  IF(IFAIL /= 0)THEN
     PRINT *,' CODE DE RETOUR DE DGETRF ',IFAIL
     STOP; ENDIF
  !
  CALL DGETRI(NVP,WOKE4(1:NVP,1:NVP),NVP,WOKE3(1:NVP,1:NVP),WOKE6,NB,IFAIL)
  VPX1 = WOKE4
  !
  IF(IFAIL /= 0)THEN
     PRINT *,' CODE DE RETOUR DE DGETRI ',IFAIL
     STOP; ENDIF
  !
  !
  4567    FORMAT(6(1PE12.4))
  !
   CALL MXM(VPX,NI+2,VPX1,NI+2,WOKE3,NI+2)
  !
   PRINT *,'EIGGLP : DIAGONALE DE A*INV(A) '
   PRINT 4567,(WOKE3(I,I),I = 0 , NI+1)
   PRINT 4567,(EIGX(I),I = 0 , NI+1)
   DEALLOCATE(WOKE3, WOKE4, WOKE5, WOKE6)

   IF ( present(indice) .and. indice == 3 ) THEN
   VPX  = TRANSPOSE(VPX)
   VPX1 = TRANSPOSE(VPX1)
   ENDIF
  !
  END SUBROUTINE EIGP
  





  SUBROUTINE RESOLP(NI,NJ,P,SP,AS,AP,AN,IDU,IFU,JDU,JFU,VPX,VPRX,VPRX1)
!(NI,NJ,Phi,DIV,AP(:,:)%CH(4), AP(:,:)%CH(6),AP(:,:)%CH(5),IDV,IFV,JDU,JFU,VPX(0),VPRX(0,0),VPRX1(0,0))
  IMPLICIT NONE
  INTEGER:: I,J, IND
  INTEGER, INTENT(IN) :: NI, NJ, IDU,IFU,JDU,JFU
  REAL(nk), DIMENSION(0:NI+1,0:NJ+1), INTENT(IN)    :: AS, AP, AN
  REAL(nk), DIMENSION(0:NI+1,0:NI+1), INTENT(IN)    :: VPRX,VPRX1
  REAL(nk), DIMENSION(0:NI+1       ), INTENT(IN)    :: VPX
  REAL(nk), DIMENSION(0:NI+1,0:NJ+1), INTENT(INOUT) :: P, SP
  REAL(nk), DIMENSION(:,:), ALLOCATABLE             :: WORK1, WORK2, WORK3
  REAL(nk), DIMENSION(:  ), ALLOCATABLE             :: WORK

  P(:,:) = 0.

  CALL MXM(VPRX1,NI+2,SP,NI+2,P,NJ+2)                 !!! VPRX1 * DIV = P!

  ALLOCATE( WORK(0:NI+1), WORK1(0:NI+1,0:NJ+1), WORK2(0:NI+1,0:NJ+1) )
  ALLOCATE( WORK3(0:NI+1,0:NJ+1) )

  WORK = 0.; WORK1 = 0.; WORK2 = 0.; WORK3 = AP

   DO I = 1 , NI
       WORK3(I,:) = WORK3(I,:) + VPX(I);
       IF(ABS(VPX(I)) < 1.E-7) then 
         WORK3(IND,JDU) = 1.E30     !IND = I
       end if
   ENDDO            !!! VPX + B = WORK3
  !

  DO I = IDU, IFU
  WORK1(I,JDU) = AN(I,JDU)/WORK3(I,JDU)
  WORK2(I,JDU) = P (I,JDU)/WORK3(I,JDU)   !!!   P / ( VPX + B ) = WORK2
  ENDDO
  !
  DO J = JDU+1 , JFU; DO I = IDU , IFU
     WORK(I) = 1./( WORK3(I,J)-AS(I,J)*WORK1(I,J-1) )
     WORK1(I,J) = AN(I,J)*WORK(I)
     WORK2(I,J) = ( P(I,J)+AS(I,J)*WORK2(I,J-1) )*WORK(I)
  ENDDO; ENDDO
  SP(IDU:IFU,JFU) = WORK2(IDU:IFU,JFU)
  !
  DO J = JFU-1 , JDU, -1; DO I = IDU , IFU
     SP(I,J) = WORK1(I,J)*SP(I,J+1)+WORK2(I,J)
  ENDDO; ENDDO

  CALL MXM(VPRX ,NI+2,SP,NI+2,P,NJ+2)
! !
  DEALLOCATE(WORK, WORK1, WORK2, WORK3)
  !
  END SUBROUTINE RESOLP


#if (DIMENSION_GEO ==3)
  SUBROUTINE RESOLP_3D(NI,NJ,NK,P,SP,AS,AP,AN,IDU,IFU,JDU,JFU,KDU,KFU,VPX,VPRX,VPRX1,VPZ,VPRZ,VPRZ1)
  !(NI,NJ,Phi,DIV,AP(:,:)%CH(4), AP(:,:)%CH(6),AP(:,:)%CH(5),IDV,IFV,JDU,JFU,VPX(0),VPRX(0,0),VPRX1(0,0))
  IMPLICIT NONE
  INTEGER:: I,J,K,IND,KND
  INTEGER, INTENT(IN) :: NI, NJ, NK, IDU,IFU,JDU,JFU,KDU,KFU
  REAL(8), DIMENSION(0:NI+1,0:NJ+1,0:NK+1), INTENT(IN) :: AS, AP, AN
  REAL(8), DIMENSION(0:NI+1,0:NI+1), INTENT(IN) :: VPRX,VPRX1
  REAL(8), DIMENSION(0:NI+1       ), INTENT(IN) :: VPX
  REAL(8), DIMENSION(0:NK+1,0:NK+1), INTENT(IN) :: VPRZ,VPRZ1
  REAL(8), DIMENSION(0:NK+1       ), INTENT(IN) :: VPZ
  REAL(8), DIMENSION(0:NI+1,0:NJ+1,0:NK+1), INTENT(INOUT) :: P, SP
  REAL(8), DIMENSION(0:NI+1,0:NJ+1,0:NK+1)                :: PP
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: WORK1, WORK2, WORK3
  REAL(8), DIMENSION(:,:  ), ALLOCATABLE :: WORK
  
  P = 0.; PP = 0.
    
    do K = KDU , KFU
  CALL MXM(VPRX1,NI+2,SP(:,:,K) ,NI+2,PP(:,:,K) ,NJ+2)                 !!! VPRX1 * DIV = PP
    end do

    do I = IDU , IFU
  CALL MXM(PP(I,:,:)   ,NJ+2,VPRZ1,NK+2,P(I,:,:) ,NK+2)                 !!! VPRX1 * DIV * VPXZ = P
    end do

  ALLOCATE( WORK(0:NI+1,0:NK+1), WORK1(0:NI+1,0:NJ+1,0:NK+1), WORK2(0:NI+1,0:NJ+1,0:NK+1) )
  ALLOCATE( WORK3(0:NI+1,0:NJ+1,0:NK+1) )

  WORK = 0.; WORK1 = 0.; WORK2 = 0.; WORK3 = AP

  DO K = KDU , KFU
    DO I = IDU , IFU
      WORK3(I,:,K) = WORK3(I,:,K) + VPX(I) + VPZ(K)
      IF(ABS(VPX(I)+VPZ(K)) < 1.E-7) THEN
         WORK3(I,JDU,K) = 1.E25; IND = I; KND = K
         print*, IND, KND
      ENDIF
    ENDDO
  ENDDO           !!! VPX + VPZ + B = WORK3
     
  DO I = IDU, IFU; DO K = KDU, KFU
  WORK1(I,JDU,K) = AN(I,JDU,K)/WORK3(I,JDU,K)
  WORK2(I,JDU,K) = P (I,JDU,K)/WORK3(I,JDU,K)   !!!   P / ( VPX + VPZ + B ) = WORK2
  ENDDO; ENDDO
  !
  DO J = JDU+1 , JFU; DO K = KDU , KFU; DO I = IDU , IFU
     WORK(I,K) = 1./( WORK3(I,J,K)-AS(I,J,K)*WORK1(I,J-1,K) )
     WORK1(I,J,K) = AN(I,J,K)*WORK(I,K)
     WORK2(I,J,K) = ( P (I,J,K)+AS(I,J,K)*WORK2(I,J-1,K) )*WORK(I,K)
  ENDDO; ENDDO; ENDDO
  !
  SP(IDU:IFU,JFU,KDU:KFU) = WORK2(IDU:IFU,JFU,KDU:KFU)

  DO J = JFU-1 , JDU, -1; DO K = KDU , KFU; DO I = IDU , IFU
     SP(I,J,K) = WORK2(I,J,K)+WORK1(I,J,K)*SP(I,J+1,K)
  ENDDO; ENDDO; ENDDO

  PP(:,:,:) = 0. ; P(:,:,:) = 0.
  DO K = KDU , KFU
  CALL MXM(VPRX ,NI+2,SP(:,:,K),NI+2,PP(:,:,K),NJ+2)
  ENDDO
   
  DO I = IDU , IFU
  CALL MXM(PP(I,:,:),NJ+2,VPRZ,NK+2,P(I,:,:),NK+2)
  ENDDO
  
  DEALLOCATE(WORK, WORK1, WORK2, WORK3)
  !
  END SUBROUTINE RESOLP_3D

#endif




!!!!!!!!!!!(VPRX ,NI+2,SP(:,:,K),NI+2,PP(:,:,K),NJ+2)
  SUBROUTINE MXM(A,NAR,B,NBR,C,NCC,ALPHA,BETA)
  IMPLICIT NONE
  INTEGER :: NAR, NBR, NCC ! STOPPING ICI
  REAL(nk), OPTIONAL, INTENT(IN) :: ALPHA, BETA
  REAL(nk) :: ONE, ZERO
  REAL(nk), DIMENSION(NAR,NBR), INTENT(IN) :: A
  REAL(nk), DIMENSION(NBR,NCC), INTENT(IN) :: B
  REAL(nk), DIMENSION(NAR,NCC), INTENT(INOUT) :: C
!
!=====MULTIPLICATION MATRICIELLE
!
  IF(PRESENT(ALPHA) .AND. PRESENT(BETA)) THEN
   CALL DGEMM('N','N',NAR,NCC,NBR,ALPHA,A,NAR,B,NBR,BETA,C,NAR)
  ELSE
   ONE = 1.0; ZERO = 0.
   CALL DGEMM('N','N',NAR,NCC,NBR,ONE,A,NAR,B,NBR,ZERO,C,NAR)
  ENDIF
!
  END SUBROUTINE MXM



end module mPression

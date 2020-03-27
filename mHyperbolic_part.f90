#include "definitions.h"
  
  
module mHyperbolic_part
  
  use mBase
  use mChampsprof
  use mOutils!, only: centre_2

  implicit none
 
contains






  
  subroutine calcaul_termes_NL_partie_hyperbolique(u1, u2, t11, t12, t22, A1_NL, A2_NL, dt, it, u3, t13, t23, t33, A3_NL)
    implicit none   
    integer  :: i, j, k 
    real(nk), intent(in) :: dt
    integer, intent(in)  :: it
#if ( DIMENSION_GEO == 2 )
    real(nk), dimension(0:N1+1,0:N2+1), intent(in) :: u1, u2, t11, t12, t22
    real(nk), dimension(1:N1,1:N2,1:5), intent(out):: A1_NL, A2_NL
    real(nk), dimension(0:N1+1,0:N2+1), optional   :: u3, t13, t23, t33, A3_NL
#elif ( DIMENSION_GEO == 3 )
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in):: u1, u2, t11, t12, t22, u3, t13, t23, t33
    real(nk), dimension(1:N1,1:N2,1:N3,1:9), intent(out) :: A1_NL, A2_NL, A3_NL
#endif

#if(EXP_TREATMENT == 1)
    Tau_trans = exp(mod(it, Re_zero)*dt / Weissenberg)
    if (  mod(it, Re_zero) == 0) then
        Tau_trans = exp( Re_zero*dt / Weissenberg )
    endif
#endif

#if ( QUASI_LINEAR == 1 )
#if ( DIMENSION_GEO == 2 )
      call valeurs_propres_A1_Lambda_L_R (u1, t11, t12, A(1)%Lambda, A(1)%L, A(1)%R)
      call valeurs_propres_A2_Lambda_L_R (u2, t12, t22, A(2)%Lambda, A(2)%L, A(2)%R)   

      call schema_WENO_A1 (u1, u2, t11, t12, t22, A(1)%Lambda, A(1)%L, A(1)%R, A1_NL)
      call schema_WENO_A2 (u1, u2, t11, t12, t22, A(2)%Lambda, A(2)%L, A(2)%R, A2_NL)
#elif (DIMENSION_GEO == 3)
      call valeurs_propres_A1_Lambda_L_R (u1, t11, t12, A(1)%Lambda, A(1)%L, A(1)%R, t13, t23)
      call valeurs_propres_A2_Lambda_L_R (u2, t12, t22, A(2)%Lambda, A(2)%L, A(2)%R, t13, t23)   
      call valeurs_propres_A3_Lambda_L_R (u3, t13, t23, A(3)%Lambda, A(3)%L, A(3)%R, t12, t33)   

      call schema_WENO_A1 (u1, u2, t11, t12, t22, A(1)%Lambda, A(1)%L, A(1)%R, A1_NL, u3, t13, t23, t33)  
      call schema_WENO_A2 (u1, u2, t11, t12, t22, A(2)%Lambda, A(2)%L, A(2)%R, A2_NL, u3, t13, t23, t33)  
      call schema_WENO_A3 (u1, u2, t11, t12, t22, A(3)%Lambda, A(3)%L, A(3)%R, A3_NL, u3, t13, t23, t33)  
#endif

#elif ( QUASI_LINEAR == 0 )

#if ( DIMENSION_GEO == 2 )
    do j=1, N2                 
       do i=1, N1
          call calcaul_termes_Convections_1(i, j, u1, u2, t11, t12, t22, A1_NL)
          call calcaul_termes_Convections_2(i, j, u1, u2, t11, t12, t22, A2_NL)
       end do
    end do
#elif ( DIMENSION_GEO == 3 )
   do k=1, N3
      do j= 1, N2                 
         do i= 1, N1
            call calcaul_termes_Convections_1(i, j, u1, u2 , t11, t12, t22, A1_NL, k, u3, t13, t23, t33)
            call calcaul_termes_Convections_2(i, j, u1, u2 , t11, t12, t22, A2_NL, k, u3, t13, t23, t33)
            call calcaul_termes_Convections_3(i, j, u1, u2 , t11, t12, t22, A3_NL, k, u3, t13, t23, t33)
         end do
      end do
   end do   
#endif

#endif 
  end subroutine calcaul_termes_NL_partie_hyperbolique




               




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     2D     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine valeurs_propres_A1_Lambda_L_R(u1, t11, t12, Lambda, L, R, t13, t23)
    implicit none   
    integer :: i, j, k
    real(nk):: const 
#if ( DIMENSION_GEO == 2 )
    real(nk), dimension(0:N1+1,0:N2+1),           intent(in)   :: u1, t11, t12
    real(nk), dimension(1:N1  ,1:N2  , 1:5),      intent(out)  :: Lambda
    real(nk), dimension(1:N1  ,1:N2  , 1:5, 1:5), intent(out)  :: L, R
    real(nk), optional                                         :: t13, t23
#elif  ( DIMENSION_GEO == 3 )
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1),    intent(in)   :: u1, t11, t12, t13, t23
    real(nk), dimension(1:N1,1:N2,1:N3,1:9),      intent(out)  :: Lambda
    real(nk), dimension(1:N1,1:N2,1:N3,1:9,1:9),  intent(out)  :: L, R
#endif
       
#if (DIMENSION_GEO == 2)
   do j= 1, N2                 
       do i=1, N1  

#if ( EXP_TREATMENT == 0 )
    if (t11(i,j) < - 1._nk/Ma2) then 
      Cent_or_HOUC_A(i,j) = 1
    else 
      Cent_or_HOUC_A(i,j) = 0 
      const = sqrt(t11(i,j) + 1._nk/Ma2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lambda
      Lambda(i,j,1) =  u1(i,j) - sqrt2*const
      Lambda(i,j,2) =  u1(i,j) + sqrt2*const
      Lambda(i,j,3) =  u1(i,j) -         const
      Lambda(i,j,4) =  u1(i,j) +         const
      Lambda(i,j,5) =  u1(i,j) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! L Matrice
      L(i,j,1,1) = 1.  
      L(i,j,1,2) = 1.  

      L(i,j,2,3) = 1.  
      L(i,j,2,4) = 1.  

      L(i,j,5,5) = 1.  

      L(i,j,3,1) = sqrt2*const   
      L(i,j,3,2) =-sqrt2*const   

      L(i,j,4,3) =       const   
      L(i,j,4,4) =-      const   

      L(i,j,5,3) = 2.* t12(i,j)/const   
      L(i,j,5,4) =-2.* t12(i,j)/const  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! R Matrice
      R(i,j,1,1) = .5_nk
      R(i,j,2,1) = .5_nk

      R(i,j,3,2) = .5_nk
      R(i,j,4,2) = .5_nk
 
      R(i,j,1,3) = .5/sqrt2/const   
      R(i,j,2,3) =-.5/sqrt2/const   

      R(i,j,3,4) = .5_nk   /const   
      R(i,j,4,4) =-.5_nk   /const 
      R(i,j,5,4) =- 2.* t12(i,j)/(t11(i,j) + 1._nk/Ma2)   

#elif ( EXP_TREATMENT == 1 )
    if (t11(i,j)/Tau_trans < - 1._nk/Ma2) then 
      Cent_or_HOUC_A(i,j) = 1
    else 
      Cent_or_HOUC_A(i,j) = 0 
      const = sqrt(t11(i,j)*Tau_trans + (Tau_trans**2)/Ma2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lambda
      Lambda(i,j,1) =  u1(i,j) - sqrt2*const / Tau_trans
      Lambda(i,j,2) =  u1(i,j) + sqrt2*const / Tau_trans
      Lambda(i,j,3) =  u1(i,j) -         const / Tau_trans
      Lambda(i,j,4) =  u1(i,j) +         const / Tau_trans
      Lambda(i,j,5) =  u1(i,j) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! L Matrice
      L(i,j,1,1) = 1.  
      L(i,j,1,2) = 1.  

      L(i,j,2,3) = 1.  
      L(i,j,2,4) = 1.  

      L(i,j,5,5) = 1.  

      L(i,j,3,1) = sqrt2*const   
      L(i,j,3,2) =-sqrt2*const   

      L(i,j,4,3) =       const   
      L(i,j,4,4) =-      const   

      L(i,j,5,3) = 2.*t12(i,j)*Tau_trans/const 
      L(i,j,5,4) =-2.*t12(i,j)*Tau_trans/const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! R Matrice
      R(i,j,1,1) = .5_nk
      R(i,j,2,1) = .5_nk

      R(i,j,3,2) = .5_nk
      R(i,j,4,2) = .5_nk
 
      R(i,j,1,3) = .5/sqrt2/const   
      R(i,j,2,3) =-.5/sqrt2/const   

      R(i,j,3,4) = .5_nk   /const   
      R(i,j,4,4) =-.5_nk   /const 
      R(i,j,5,4) =- 2.* t12(i,j)*Tau_trans/(t11(i,j)*Tau_trans + (Tau_trans**2)/Ma2) 

      R(i,j,5,5) =  1._nk
#endif
    end if

       end do
   end do
#elif (DIMENSION_GEO == 3)
   do k= 1, N3
      do j= 1, N2                 
         do i= 1, N1    
    
    if (t11(i,j,k) < - 1._nk/Ma2) then 
      Cent_or_HOUC_A(i,j,k) = 1
    else 
      Cent_or_HOUC_A(i,j,k) = 0
      const = sqrt(t11(i,j,k) + 1._nk/Ma2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lambda
      Lambda(i,j,k,1) =  u1(i,j,k) - sqrt2*const
      Lambda(i,j,k,2) =  u1(i,j,k) + sqrt2*const
      Lambda(i,j,k,3) =  u1(i,j,k) -       const
      Lambda(i,j,k,4) =  u1(i,j,k) -       const
      Lambda(i,j,k,5) =  u1(i,j,k) +       const
      Lambda(i,j,k,6) =  u1(i,j,k) +       const    
      Lambda(i,j,k,7) =  u1(i,j,k) 
      Lambda(i,j,k,8) =  u1(i,j,k) 
      Lambda(i,j,k,9) =  u1(i,j,k) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! L Matrice
      L(i,j,k,1,1) = 1.
      L(i,j,k,1,2) = 1.

      L(i,j,k,2,3) = 1.
      L(i,j,k,2,5) = 1.

      L(i,j,k,3,4) = 1.
      L(i,j,k,3,6) = 1.

      L(i,j,k,4,2) =-sqrt2*const
      L(i,j,k,4,1) = sqrt2*const

      L(i,j,k,5,3) = const
      L(i,j,k,5,5) =-const

      L(i,j,k,6,4) = const
      L(i,j,k,6,6) =-const

      L(i,j,k,7,3) = 2.*t12(i,j,k)/const
      L(i,j,k,7,5) =-2.*t12(i,j,k)/const
      L(i,j,k,7,7) = 1.

      L(i,j,k,8,1) =-t23(i,j,k)/sqrt2/const
      L(i,j,k,8,2) = t23(i,j,k)/sqrt2/const
      L(i,j,k,8,3) =    t13(i,j,k)/const
      L(i,j,k,8,4) =    t12(i,j,k)/const
      L(i,j,k,8,5) =-   t13(i,j,k)/const  
      L(i,j,k,8,6) =-   t12(i,j,k)/const  
      L(i,j,k,8,8) = 1.

      L(i,j,k,9,4) = 2.*t13(i,j,k)/const    
      L(i,j,k,9,6) =-2.*t13(i,j,k)/const
      L(i,j,k,9,9) = 1.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! R Matrice
      R(i,j,k,1,1)=.5_nk
      R(i,j,k,1,4)= 1./2.**(3./2.)/const

      R(i,j,k,2,1)=.5_nk
      R(i,j,k,2,4)=-1./2.**(3./2.)/const

      R(i,j,k,3,2)=.5_nk
      R(i,j,k,3,5)= 1/2./const

      R(i,j,k,4,3)=.5_nk
      R(i,j,k,4,6)= 1./2./const
    
      R(i,j,k,5,2)=.5_nk
      R(i,j,k,5,5)=-1/2./const
     
      R(i,j,k,6,3)=.5_nk
      R(i,j,k,6,6)=-1/2./const

      R(i,j,k,7,5)=-2.*t12(i,j,k)/const**2.
      R(i,j,k,7,7)= 1._nk

      R(i,j,k,8,4)= t23(i,j,k)/2./const**2.
      R(i,j,k,8,5)=-t13(i,j,k)   /const**2.
      R(i,j,k,8,6)=-t12(i,j,k)   /const**2.
      R(i,j,k,8,8)= 1._nk

      R(i,j,k,9,6)=-2.*t13(i,j,k)/const**2.
      R(i,j,k,9,9)= 1._nk
    end if
        end do 
      end do
   end do
#endif   


  end subroutine valeurs_propres_A1_Lambda_L_R

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine valeurs_propres_A2_Lambda_L_R(u2, t12, t22, Lambda, L, R, t13, t23)
    implicit none   
    integer :: i, j, k
#if ( DIMENSION_GEO == 2 )
    real(nk), dimension(0:N1+1,0:N2+1),           intent(in)   :: u2, t12, t22
    real(nk), dimension(1:N1  ,1:N2  , 1:5),      intent(out)  :: Lambda
    real(nk), dimension(1:N1  ,1:N2  , 1:5, 1:5), intent(out)  :: L, R
    real(nk), dimension(0:N1+1,0:N2+1), optional               :: t13, t23
#else 
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1),   intent(in)   :: u2, t12, t22, t13, t23
    real(nk), dimension(1:N1,1:N2,1:N3,1:9),     intent(out)  :: Lambda
    real(nk), dimension(1:N1,1:N2,1:N3,1:9,1:9), intent(out)  :: L, R
#endif
    real(nk):: const    

#if (DIMENSION_GEO == 2)
    
   do j= 1, N2                 
       do i=1, N1 

#if ( EXP_TREATMENT == 0 )
    if (t22(i,j) < - 1._nk/Ma2) then 
      Cent_or_HOUC_B(i,j) = 1
    else
      Cent_or_HOUC_B(i,j) = 0
      const = sqrt(t22(i,j) + 1._nk/Ma2) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lambda
      Lambda(i,j,1) =  u2(i,j) - sqrt2*const
      Lambda(i,j,2) =  u2(i,j) + sqrt2*const
      Lambda(i,j,3) =  u2(i,j) -       const
      Lambda(i,j,4) =  u2(i,j) +       const
      Lambda(i,j,5) =  u2(i,j) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! L Matrice
      L(i,j,2,1) = 1.  
      L(i,j,2,2) = 1. 

      L(i,j,1,3) = 1.  
      L(i,j,1,4) = 1.  

      L(i,j,3,5) = 1.  
      L(i,j,5,1) = sqrt2*const   
      L(i,j,5,2) =-sqrt2*const   

      L(i,j,4,3) =       const   
      L(i,j,4,4) =-      const   

      L(i,j,3,3) = 2.*t12(i,j)/const   
      L(i,j,3,4) =-2.*t12(i,j)/const   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! R Matrice
      R(i,j,1,2) = .5_nk
      R(i,j,2,2) = .5_nk  

      R(i,j,3,1) = .5_nk  
      R(i,j,4,1) = .5_nk  

      R(i,j,1,5) = .5/sqrt2/const   
      R(i,j,2,5) =-.5/sqrt2/const   

      R(i,j,3,4) = .5_nk   /const   
      R(i,j,4,4) =-.5_nk   /const   

      R(i,j,5,3) =  1._nk  
      R(i,j,5,4) =- 2.*t12(i,j)/(t22(i,j) + 1._nk/Ma2)    

#elif ( EXP_TREATMENT == 1 )
    if (t22(i,j)/Tau_trans < - 1._nk/Ma2) then 
      Cent_or_HOUC_B(i,j) = 1
    else
      Cent_or_HOUC_B(i,j) = 0
      const = sqrt(t22(i,j)*Tau_trans + (Tau_trans**2)/Ma2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lambda
      Lambda(i,j,1) =  u2(i,j) - sqrt2*const / Tau_trans
      Lambda(i,j,2) =  u2(i,j) + sqrt2*const / Tau_trans
      Lambda(i,j,3) =  u2(i,j) -       const / Tau_trans
      Lambda(i,j,4) =  u2(i,j) +       const / Tau_trans
      Lambda(i,j,5) =  u2(i,j) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! L Matrice
      L(i,j,2,1) = 1.  
      L(i,j,2,2) = 1. 

      L(i,j,1,3) = 1.  
      L(i,j,1,4) = 1.  

      L(i,j,3,5) = 1.  

      L(i,j,5,1) = sqrt2*const   
      L(i,j,5,2) =-sqrt2*const   

      L(i,j,4,3) =       const   
      L(i,j,4,4) =-      const   

      L(i,j,3,3) = 2.*t12(i,j)*Tau_trans/const
      L(i,j,3,4) =-2.*t12(i,j)*Tau_trans/const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! R Matrice
      R(i,j,1,2) = .5_nk
      R(i,j,2,2) = .5_nk  

      R(i,j,3,1) = .5_nk  
      R(i,j,4,1) = .5_nk  

      R(i,j,1,5) = .5/sqrt2/const   
      R(i,j,2,5) =-.5/sqrt2/const   

      R(i,j,3,4) = .5_nk   /const   
      R(i,j,4,4) =-.5_nk   /const   

      R(i,j,5,3) =  1._nk  
      R(i,j,5,4) =- 2.*t12(i,j)*Tau_trans/(t22(i,j)*Tau_trans + (Tau_trans**2)/Ma2)
#endif
    end if

       end do
   end do

#elif (DIMENSION_GEO == 3)
   do k= 1, N3
      do j= 1, N2                 
         do i= 1, N1    
      
    if (t22(i,j,k) < - 1._nk/Ma2) then 
      Cent_or_HOUC_B(i,j,k) = 1
    else
      Cent_or_HOUC_B(i,j,k) = 0
      const = sqrt(t22(i,j,k) + 1._nk/Ma2)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lambda
      Lambda(i,j,k,1) =  u2(i,j,k) - sqrt2*const
      Lambda(i,j,k,2) =  u2(i,j,k) + sqrt2*const
      Lambda(i,j,k,3) =  u2(i,j,k) -       const
      Lambda(i,j,k,4) =  u2(i,j,k) -       const
      Lambda(i,j,k,5) =  u2(i,j,k) +       const
      Lambda(i,j,k,6) =  u2(i,j,k) +       const
      Lambda(i,j,k,7) =  u2(i,j,k) 
      Lambda(i,j,k,8) =  u2(i,j,k)
      Lambda(i,j,k,9) =  u2(i,j,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! L Matrice
      L(i,j,k,1,3) = 1.
      L(i,j,k,1,5) = 1.

      L(i,j,k,2,1) = 1.
      L(i,j,k,2,2) = 1.
    
      L(i,j,k,3,4) = 1.
      L(i,j,k,3,6) = 1.

      L(i,j,k,4,3) = 2.*t12(i,j,k)/const
      L(i,j,k,4,5) =-2.*t12(i,j,k)/const
      L(i,j,k,4,7) = 1.

      L(i,j,k,5,3) = const
      L(i,j,k,5,5) =-const

      L(i,j,k,6,1) =-t13(i,j,k)/sqrt2/const
      L(i,j,k,6,2) = t13(i,j,k)/sqrt2/const
      L(i,j,k,6,3) = t23(i,j,k)/const
      L(i,j,k,6,4) = t12(i,j,k)/const  
      L(i,j,k,6,5) =-t23(i,j,k)/const
      L(i,j,k,6,6) =-t12(i,j,k)/const   
      L(i,j,k,6,8) = 1.

      L(i,j,k,7,1) = sqrt2*const
      L(i,j,k,7,2) =-sqrt2*const

      L(i,j,k,8,4) = const
      L(i,j,k,8,6) =-const

      L(i,j,k,9,4) = 2.* t23(i,j,k)/const 
      L(i,j,k,9,6) =-2.* t23(i,j,k)/const
      L(i,j,k,9,9) = 1.   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! R Matrice
      R(i,j,k,1,2)=.5_nk
      R(i,j,k,1,7)= 1._nk/2._nk**(3./2.)/const

      R(i,j,k,2,2)=.5_nk
      R(i,j,k,2,7)=-1./2.**(3./2.)/const

      R(i,j,k,3,1)= 0.5_nk
      R(i,j,k,3,5)= 1./2./const

      R(i,j,k,4,3)= 0.5_nk
      R(i,j,k,4,8)= 1./2./const

      R(i,j,k,5,1)= 0.5_nk
      R(i,j,k,5,5)=-1./2./const

      R(i,j,k,6,3)= 0.5_nk
      R(i,j,k,6,8)=-1._nk/2._nk/const

      R(i,j,k,7,4)= 1._nk
      R(i,j,k,7,5)=-2.*t12(i,j,k)/const**2

      R(i,j,k,8,5)=-t23(i,j,k)/const**2
      R(i,j,k,8,6)= 1._nk
      R(i,j,k,8,7)= t13(i,j,k)/2./const**2
      R(i,j,k,8,8)=-t12(i,j,k)/const**2

      R(i,j,k,9,8)=-2.*t23(i,j,k)/const**2
      R(i,j,k,9,9)= 1._nk
    end if
        end do 
      end do
   end do
#endif  
  end subroutine valeurs_propres_A2_Lambda_L_R

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if (DIMENSION_GEO == 3)
  subroutine valeurs_propres_A3_Lambda_L_R(u3, t13, t23, Lambda, L, R, t12, t33)
    implicit none   
    integer :: i, j, k
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1),   intent(in)   :: u3, t12, t13, t23, t33
    real(nk), dimension(1:N1,1:N2,1:N3,1:9),     intent(out)  :: Lambda
    real(nk), dimension(1:N1,1:N2,1:N3,1:9,1:9), intent(out)  :: L, R
    real(nk):: const    

   do k= 1, N3
      do j= 1, N2                 
         do i= 1, N1    
    if (t33(i,j,k) < - 1._nk/Ma2) then 
      Cent_or_HOUC_C(i,j,k) = 1
    else
      Cent_or_HOUC_C(i,j,k) = 0
      const = sqrt(t33(i,j,k) + 1._nk/Ma2)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lambda
      Lambda(i,j,k,1) =  u3(i,j,k) - sqrt2*const
      Lambda(i,j,k,2) =  u3(i,j,k) + sqrt2*const
      Lambda(i,j,k,3) =  u3(i,j,k) -       const
      Lambda(i,j,k,4) =  u3(i,j,k) -       const
      Lambda(i,j,k,5) =  u3(i,j,k) +       const
      Lambda(i,j,k,6) =  u3(i,j,k) +       const
      Lambda(i,j,k,7) =  u3(i,j,k) 
      Lambda(i,j,k,8) =  u3(i,j,k)
      Lambda(i,j,k,9) =  u3(i,j,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! L Matrice
      L(i,j,k,1,3) = 1.
      L(i,j,k,1,5) = 1.

      L(i,j,k,2,4) = 1.
      L(i,j,k,2,6) = 1.

      L(i,j,k,3,1) = 1.  
      L(i,j,k,3,2) = 1.  

      L(i,j,k,4,3) = 2.*t13(i,j,k)/const
      L(i,j,k,4,5) =-2.*t13(i,j,k)/const
      L(i,j,k,4,7) = 1.

      L(i,j,k,5,1) =-t12(i,j,k)/sqrt2/const  
      L(i,j,k,5,2) = t12(i,j,k)/sqrt2/const 
      L(i,j,k,5,3) =    t23(i,j,k)/const   
      L(i,j,k,5,4) =    t13(i,j,k)/const  
      L(i,j,k,5,5) =-   t23(i,j,k)/const
      L(i,j,k,5,6) =-   t13(i,j,k)/const
      L(i,j,k,5,8) = 1.

      L(i,j,k,6,3) =    const      
      L(i,j,k,6,5) =-   const   

      L(i,j,k,7,4) = 2.*t23(i,j,k)/const   
      L(i,j,k,7,6) =-2.*t23(i,j,k)/const
      L(i,j,k,7,9) = 1.

      L(i,j,k,8,4) = const    
      L(i,j,k,8,6) =-const

      L(i,j,k,9,1) = sqrt2*const    
      L(i,j,k,9,2) =-sqrt2*const  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! R Matrice
      R(i,j,k,1,3)=.5_nk
      R(i,j,k,1,9)= 1._nk/2._nk**(3./2.)/const

      R(i,j,k,2,3)=.5_nk
      R(i,j,k,2,9)=-1._nk/2._nk**(3./2.)/const

      R(i,j,k,3,1)=.5_nk
      R(i,j,k,3,6)= 1._nk/2._nk/const

      R(i,j,k,4,2)=.5_nk
      R(i,j,k,4,8)= 1._nk/2._nk/const

      R(i,j,k,5,1)=.5_nk
      R(i,j,k,5,6)=-1._nk/2._nk/const

      R(i,j,k,6,2)=.5_nk
      R(i,j,k,6,8)=-1._nk/2._nk/const

      R(i,j,k,7,4)= 1._nk
      R(i,j,k,7,6)=-2._nk*t13(i,j,k)/const**2

      R(i,j,k,8,5)= 1._nk
      R(i,j,k,8,6)=-t23(i,j,k)/const**2
      R(i,j,k,8,8)=-t13(i,j,k)/const**2
      R(i,j,k,8,9)= t12(i,j,k)/2./const**2

      R(i,j,k,9,7)= 1._nk
      R(i,j,k,9,8)=-2._nk*t23(i,j,k)/const**2
      
    end if
        end do 
      end do
   end do
  end subroutine valeurs_propres_A3_Lambda_L_R
#endif  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






  subroutine schema_WENO_A1 (u1, u2, t11, t12, t22, Lambda, L, R, NL, u3, t13, t23, t33)
    implicit none    
    integer  :: i,ii, jj, kk

#if (DIMENSION_GEO == 2)
    real(nk), dimension(0:N1+1,0:N2+1), intent(in), target:: u1, u2, t11, t12, t22
    real(nk), dimension(1:N1,1:N2,1:5,1:5), intent(in)    :: L, R
    real(nk), dimension(1:N1,1:N2,1:5    ), intent(in)    :: Lambda
    real(nk), dimension(1:N1,1:N2,1:5    ), intent(out)   :: NL
    real(nk), dimension(1:N1,1:N2,1:5    )                :: Lambda_dRWdx
    real(nk), optional                                    :: u3, t13, t23, t33
    real(nk),dimension(5)                                 :: D_W, L_W, R_W 
    type pointorW 
    real(nk), dimension(:,:), pointer:: var
    end type pointorW 
    type(pointorW), dimension(5) :: WW

    WW(1)%var => u1 ; WW(2)%var => u2 ; WW(3)%var => t11; WW(4)%var => t12; WW(5)%var => t22

#elif (DIMENSION_GEO ==3)
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1),intent(in),target:: u1, u2, t11, t12, t22, u3, t13, t23, t33
    real(nk), dimension(1:N1,1:N2,1:N3,1:9,1:9), intent(in)    :: L, R
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    ), intent(in)    :: Lambda
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    ), intent(out)   :: NL
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    )                :: Lambda_dRWdx
    real(nk),dimension(9)  :: D_W, L_W, R_W 
    type pointorW 
    real(nk), dimension(:,:,:), pointer:: var
    end type pointorW 
    type(pointorW), dimension(9) :: WW

    WW(1)%var => u1 ; WW(2)%var => u2 ; WW(3)%var => u3 ; WW(4)%var => t11; WW(5)%var => t12 
    WW(6)%var => t13; WW(7)%var => t22; WW(8)%var => t23; WW(9)%var => t33

#endif


#if (DIMENSION_GEO ==2)
   do jj = 1, N2 
    do ii = 1, N1

    if (Cent_or_HOUC_A(ii,jj)==0) then
 !!!-------------l-2----------l-1---------l---------l+1---------l+2
 !!!R_W_A1        1            2          3          4           5  
    L_W(:)= 0. 
    R_W(:)= 0.
#if (NL_TERMS == 10 )
    do i = 1,5 
        if ( ii == 1) then
            L_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 1 )                          !! backwind_o2
            R_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 1 )                          !! backwind_o2            
        else if ( ii == 2 ) then
            L_W(i) = UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! upwind_o1
            R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o3  
        else if ( ii > 2 .and. ii < N1 - 1 ) then
            L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o3 
            R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o3
        else if ( ii == N1 - 1 ) then 
            L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o3
            R_W(i) = BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 1 )                          !! backwind_o1
        else if ( ii == N1) then
            L_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! upwind_o2
            R_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! upwind_o2
        end if                     
    end do 
#elif (NL_TERMS == 11 )
    do i = 1,5 
        if ( ii == 1) then
            L_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 1 )                          !! backwind_o2
            R_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 1 )                          !! backwind_o2            
        else if ( ii == 2 ) then
            L_W(i) = UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! upwind_o1
            R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o3  
        else if ( ii == 3 ) then
            L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o3  
            R_W(i) = WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o5 
        else if ( ii > 3 .and. ii < N1 - 2 ) then
            L_W(i) = WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o5
            R_W(i) = WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o5
        else if ( ii == N1 - 2 ) then 
            L_W(i) = WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o5
            R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o3  
        else if ( ii == N1 - 1 ) then 
            L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! HOUC_o3
            R_W(i) = BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 1 )                          !! backwind_o1
        else if ( ii == N1) then
            L_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! upwind_o2
            R_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 1 )                          !! upwind_o2
        end if                     
    end do 

#endif
   D_W(:)= 0.
   if ( Lambda(ii,jj,1) > 0) then
      D_W(1)= L_W(1)
      D_W(3)= L_W(3)                  
   elseif ( Lambda(ii,jj,1) < 0) then 
      D_W(1)= R_W(1)
      D_W(3)= R_W(3)
   else
      D_W(1)= 0.
      D_W(3)= 0.
   end if
   if ( ii == N1 ) then
      Lambda_dRWdx(ii,jj,1) = Lambda(ii,jj,1)*( R(ii,jj,2,1)*D_W(1) + R(ii,jj,2,3)*D_W(3) )  
   else 
      Lambda_dRWdx(ii,jj,1) = Lambda(ii,jj,1)*( R(ii,jj,1,1)*D_W(1) + R(ii,jj,1,3)*D_W(3) ) 
   end if

   if ( Lambda(ii,jj,2) > 0) then
      D_W(1)= L_W(1)
      D_W(3)= L_W(3)                            
   elseif ( Lambda(ii,jj,2) < 0) then 
      D_W(1)= R_W(1)
      D_W(3)= R_W(3)
   else
      D_W(1)= 0.
      D_W(3)= 0.
   end if
   if ( ii == 1 ) then
      Lambda_dRWdx(ii,jj,2) = Lambda(ii,jj,2)*( R(ii,jj,1,1)*D_W(1) + R(ii,jj,1,3)*D_W(3) ) 
   else 
      Lambda_dRWdx(ii,jj,2) = Lambda(ii,jj,2)*( R(ii,jj,2,1)*D_W(1) + R(ii,jj,2,3)*D_W(3) ) 
   end if
 
   if ( Lambda(ii,jj,3) > 0) then
      D_W(2)= L_W(2)
      D_W(4)= L_W(4)                               
   elseif ( Lambda(ii,jj,3) < 0) then 
      D_W(2)= R_W(2)
      D_W(4)= R_W(4)     
   else
      D_W(2)= 0.
      D_W(4)= 0.
   end if
   if ( ii == N1 ) then
      Lambda_dRWdx(ii,jj,3) = Lambda(ii,jj,3)*( R(ii,jj,4,2)*D_W(2) + R(ii,jj,4,4)*D_W(4) )
   else 
      Lambda_dRWdx(ii,jj,3) = Lambda(ii,jj,3)*( R(ii,jj,3,2)*D_W(2) + R(ii,jj,3,4)*D_W(4) )
   end if
      
 
   if ( Lambda(ii,jj,4) > 0) then
      D_W(2)= L_W(2)
      D_W(4)= L_W(4)                       
   elseif ( Lambda(ii,jj,4) < 0) then 
      D_W(2)= R_W(2)
      D_W(4)= R_W(4)   
   else
      D_W(2)= 0.
      D_W(4)= 0.
   end if
   if ( ii == 1 ) then
      Lambda_dRWdx(ii,jj,4) = Lambda(ii,jj,4)*( R(ii,jj,3,2)*D_W(2) + R(ii,jj,3,4)*D_W(4) )
   else 
      Lambda_dRWdx(ii,jj,4) = Lambda(ii,jj,4)*( R(ii,jj,4,2)*D_W(2) + R(ii,jj,4,4)*D_W(4) )
   end if
      
 
   if ( Lambda(ii,jj,5) > 0) then
      D_W(4)= L_W(4)
      D_W(5)= L_W(5)                               
   elseif ( Lambda(ii,jj,5) < 0) then 
      D_W(4)= R_W(4)
      D_W(5)= R_W(5)     
   else
      D_W(4)= 0.
      D_W(5)= 0.
   end if
      Lambda_dRWdx(ii,jj,5) = Lambda(ii,jj,5)*( R(ii,jj,5,4)*D_W(4) + R(ii,jj,5,5)*D_W(5) )


        NL(ii,jj,1) = L(ii,jj,1,1)*Lambda_dRWdx(ii,jj,1) + L(ii,jj,1,2)*Lambda_dRWdx(ii,jj,2)
        NL(ii,jj,2) = L(ii,jj,2,3)*Lambda_dRWdx(ii,jj,3) + L(ii,jj,2,4)*Lambda_dRWdx(ii,jj,4)
        NL(ii,jj,3) = L(ii,jj,3,1)*Lambda_dRWdx(ii,jj,1) + L(ii,jj,3,2)*Lambda_dRWdx(ii,jj,2)
        NL(ii,jj,4) = L(ii,jj,4,3)*Lambda_dRWdx(ii,jj,3) + L(ii,jj,4,4)*Lambda_dRWdx(ii,jj,4)
        NL(ii,jj,5) = L(ii,jj,5,3)*Lambda_dRWdx(ii,jj,3) + L(ii,jj,5,4)*Lambda_dRWdx(ii,jj,4) &
                    + L(ii,jj,5,5)*Lambda_dRWdx(ii,jj,5) 
              else
               call calcaul_termes_Convections_1(ii, jj, u1, u2, t11, t12, t22, NL)
              end if
      end do
   end do
#elif (DIMENSION_GEO == 3)
   do kk=1, N3
      do jj= 1, N2                 
         do ii= 1, N1

    if (Cent_or_HOUC_A(ii,jj,kk)==0) then
 !!!-------------l-2----------l-1---------l---------l+1---------l+2
 !!!R_W_A1        1            2          3          4           5  
    L_W(:)= 0. 
    R_W(:)= 0. 
#if ( NL_TERMS == 10 )
    do i = 1,9
      if ( ii == 1) then
          L_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! backwind_o2
          R_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! backwind_o2
      else if ( ii == 2 ) then
          L_W(i) = UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! upwind_o1
          R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o3  
      else if ( ii > 2 .and. ii < N1 - 1 ) then
          L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o3 
          R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o3
      else if ( ii == N1 - 1 ) then
          L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o3
          R_W(i) = BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! backwind_o1
      else if ( ii == N1) then
          L_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! upwind_o2
          R_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! upwind_o2
      end if                     
    end do 
#elif  ( NL_TERMS == 11 )
    do i = 1,9
      if ( ii == 1) then
          L_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! backwind_o2
          R_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! backwind_o2                      
      else if ( ii == 2 ) then
          L_W(i) = UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! upwind_o1
          R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o3  
      else if ( ii == 3 ) then
          L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o3
          R_W(i) = WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o5
      else if ( ii > 3 .and. ii < N1 - 2 ) then
          L_W(i) = WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o5
          R_W(i) = WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o5
      else if ( ii == N1 - 2 ) then
          L_W(i) = WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o5
          R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o3
      else if ( ii == N1 - 1 ) then
          L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! HOUC_o3
          R_W(i) = BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! backwind_o1
      else if ( ii == N1) then
          L_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! upwind_o2
          R_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 1, kk )                            !! upwind_o2
      end if                     
    end do 
#endif
 
    D_W(:)= 0.
  if ( Lambda(ii,jj,kk,1) > 0) then
      D_W(1)= L_W(1)
      D_W(4)= L_W(4)                            
  elseif ( Lambda(ii,jj,kk,1) < 0) then 
      D_W(1)= R_W(1)
      D_W(4)= R_W(4)  
  else
      D_W(1)= 0.
      D_W(4)= 0.   
  end if
     Lambda_dRWdx(ii,jj,kk,1) = Lambda(ii,jj,kk,1)*( R(ii,jj,kk,1,1)*D_W(1) + R(ii,jj,kk,1,4)*D_W(4) )  
 
  if ( Lambda(ii,jj,kk,2) > 0) then
      D_W(1)= L_W(1)
      D_W(4)= L_W(4)                           
  elseif ( Lambda(ii,jj,kk,2) < 0) then 
      D_W(1)= R_W(1)
      D_W(4)= R_W(4)
  else
      D_W(1)= 0.
      D_W(4)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,2) = Lambda(ii,jj,kk,1)*( R(ii,jj,kk,2,1)*D_W(1) + R(ii,jj,kk,2,4)*D_W(4) )     
 
  if ( Lambda(ii,jj,kk,3) > 0) then
      D_W(2)= L_W(2)
      D_W(5)= L_W(5)                             
  elseif ( Lambda(ii,jj,kk,3) < 0) then 
      D_W(2)= R_W(2)
      D_W(5)= R_W(5)
  else
      D_W(2)= 0.
      D_W(5)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,3) = Lambda(ii,jj,kk,3)*( R(ii,jj,kk,3,2)*D_W(2) + R(ii,jj,kk,3,5)*D_W(5) )   
 
  if ( Lambda(ii,jj,kk,4) > 0) then
      D_W(3)= L_W(3)
      D_W(6)= L_W(6)                              
  elseif ( Lambda(ii,jj,kk,4) < 0) then 
      D_W(3)= R_W(3)
      D_W(6)= R_W(6)
  else
      D_W(3)= 0.
      D_W(6)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,4) = Lambda(ii,jj,kk,4)*( R(ii,jj,kk,4,3)*D_W(3) + R(ii,jj,kk,4,6)*D_W(6) )    
 
  if ( Lambda(ii,jj,kk,5) > 0) then
      D_W(2)= L_W(2)
      D_W(5)= L_W(5)                              
  elseif ( Lambda(ii,jj,kk,5) < 0) then 
      D_W(2)= R_W(2)
      D_W(5)= R_W(5)
  else
      D_W(2)= 0.
      D_W(5)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,5) = Lambda(ii,jj,kk,5)*( R(ii,jj,kk,5,2)*D_W(2) + R(ii,jj,kk,5,5)*D_W(5) ) 
 
  if ( Lambda(ii,jj,kk,6) > 0) then
      D_W(3)= L_W(3)
      D_W(6)= L_W(6)                             
  elseif ( Lambda(ii,jj,kk,6) < 0) then 
      D_W(3)= R_W(3)
      D_W(6)= R_W(6)
  else
      D_W(3)= 0.
      D_W(6)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,6) = Lambda(ii,jj,kk,6)*( R(ii,jj,kk,6,3)*D_W(3) + R(ii,jj,kk,6,6)*D_W(6) )      
 
  if ( Lambda(ii,jj,kk,7) > 0) then
      D_W(5)= L_W(5)
      D_W(7)= L_W(7)                             
  elseif ( Lambda(ii,jj,kk,7) < 0) then 
      D_W(5)= R_W(5)
      D_W(7)= R_W(7)
  else
      D_W(5)= 0.
      D_W(7)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,7) = Lambda(ii,jj,kk,7)*( R(ii,jj,kk,7,5)*D_W(7) + R(ii,jj,kk,7,7)*D_W(7) )   
 
  if ( Lambda(ii,jj,kk,8) > 0) then
      D_W(4)= L_W(4)
      D_W(5)= L_W(5)
      D_W(6)= L_W(6)    
      D_W(8)= L_W(8)            
          
  elseif ( Lambda(ii,jj,kk,8) < 0) then 
      D_W(4)= R_W(4)
      D_W(5)= R_W(5)
      D_W(6)= R_W(6) 
      D_W(8)= R_W(8)
  else
      D_W(4)= 0.
      D_W(5)= 0.
      D_W(6)= 0. 
      D_W(8)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,8) = Lambda(ii,jj,kk,8)*( R(ii,jj,kk,8,4)*D_W(4) + R(ii,jj,kk,8,5)*D_W(5) &
                                                   + R(ii,jj,kk,8,6)*D_W(6) + R(ii,jj,kk,8,8)*D_W(8) )  
 
  if ( Lambda(ii,jj,kk,9) > 0) then
      D_W(6)= L_W(6)
      D_W(9)= L_W(9)                             
  elseif ( Lambda(ii,jj,kk,9) < 0) then 
      D_W(6)= R_W(6)
      D_W(9)= R_W(9)
  else
      D_W(6)= 0.
      D_W(9)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,9) = Lambda(ii,jj,kk,9)*( R(ii,jj,kk,9,6)*D_W(6) + R(ii,jj,kk,9,9)*D_W(9) )  

        NL(ii,jj,kk,1) =L(ii,jj,kk,1,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,1,2)*Lambda_dRWdx(ii,jj,kk,2)
        NL(ii,jj,kk,2) =L(ii,jj,kk,2,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,2,5)*Lambda_dRWdx(ii,jj,kk,5)
        NL(ii,jj,kk,3) =L(ii,jj,kk,3,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,3,6)*Lambda_dRWdx(ii,jj,kk,6)
        NL(ii,jj,kk,4) =L(ii,jj,kk,4,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,4,2)*Lambda_dRWdx(ii,jj,kk,2)
        NL(ii,jj,kk,5) =L(ii,jj,kk,5,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,5,5)*Lambda_dRWdx(ii,jj,kk,5)
        NL(ii,jj,kk,6) =L(ii,jj,kk,6,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,6,6)*Lambda_dRWdx(ii,jj,kk,6)
        NL(ii,jj,kk,7) =L(ii,jj,kk,7,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,7,5)*Lambda_dRWdx(ii,jj,kk,5) &
                       +L(ii,jj,kk,7,7)*Lambda_dRWdx(ii,jj,kk,7)
        NL(ii,jj,kk,8) =L(ii,jj,kk,8,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,8,2)*Lambda_dRWdx(ii,jj,kk,2) &
                       +L(ii,jj,kk,8,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,8,4)*Lambda_dRWdx(ii,jj,kk,4) &
                       +L(ii,jj,kk,8,5)*Lambda_dRWdx(ii,jj,kk,5) +L(ii,jj,kk,8,6)*Lambda_dRWdx(ii,jj,kk,6) &
                       +L(ii,jj,kk,8,8)*Lambda_dRWdx(ii,jj,kk,8)
        NL(ii,jj,kk,9) =L(ii,jj,kk,9,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,9,6)*Lambda_dRWdx(ii,jj,kk,6) &
                       +L(ii,jj,kk,9,9)*Lambda_dRWdx(ii,jj,kk,9)
    else                   
       call calcaul_termes_Convections_1(ii, jj, u1, u2, t11, t12, t22, NL, kk, u3, t13, t23, t33)
    endif

        end do
     end do
   end do
#endif 
    
  end subroutine schema_WENO_A1







  subroutine schema_WENO_A2 (u1, u2, t11, t12, t22, Lambda, L, R, NL, u3, t13, t23, t33)
    implicit none   
    integer  :: i  
    integer  :: ii, jj, kk

#if (DIMENSION_GEO == 2)
    real(nk), dimension(0:N1+1,0:N2+1), intent(in), target:: u1, u2, t11, t12, t22
    real(nk), dimension(1:N1,1:N2,1:5,1:5), intent(in)    :: L, R
    real(nk), dimension(1:N1,1:N2,1:5    ), intent(in)    :: Lambda
    real(nk), dimension(1:N1,1:N2,1:5    ), intent(out)   :: NL
    real(nk), dimension(1:N1,1:N2,1:5    )                :: Lambda_dRWdx
    real(nk), optional                                    :: u3, t13, t23, t33
    real(nk),dimension(5)                                 :: D_W, L_W, R_W 
    type pointorW 
    real(nk), dimension(:,:), pointer:: var
    end type pointorW 
    type(pointorW), dimension(5) :: WW

#elif (DIMENSION_GEO ==3)
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in), target :: u1, u2, t11, t12, t22, u3, t13, t23, t33
    real(nk), dimension(1:N1,1:N2,1:N3,1:9,1:9), intent(in)       :: L, R
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    ), intent(in)       :: Lambda
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    ), intent(out)      :: NL
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    )                   :: Lambda_dRWdx
    real(nk),dimension(9)                                         :: D_W, L_W, R_W 
    type pointorW 
    real(nk), dimension(:,:,:), pointer:: var
    end type pointorW 
    type(pointorW), dimension(9) :: WW

#endif


#if (DIMENSION_GEO == 2)
    WW(1)%var => u1 ; WW(2)%var => u2 ; WW(3)%var => t11; WW(4)%var => t12; WW(5)%var => t22
#elif (DIMENSION_GEO ==3)
    WW(1)%var => u1 ; WW(2)%var => u2 ; WW(3)%var => u3 ; WW(4)%var => t11; WW(5)%var => t12 
    WW(6)%var => t13; WW(7)%var => t22; WW(8)%var => t23; WW(9)%var => t33
#endif


#if (DIMENSION_GEO ==2)
   do jj = 1, N2 
    do ii = 1, N1

    if (Cent_or_HOUC_B(ii,jj)==0) then
 !!!!-----------l-2(2*NG+N1)---------l-(2*NG+N1)---------l---------l+(2*NG+N1)--------l+2(2*NG+N1)
 !!!!R_W_A2          1                    2              3             4                   5      
    L_W(:)= 0. 
    R_W(:)= 0. 
#if ( NL_TERMS == 10 )
    do i = 1,5 
        if ( jj == 1) then 
            L_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 2 )                            !! backwind_o2
            R_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 2 )                            !! backwind_o2
        else if ( jj == 2 ) then
            L_W(i) = UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! upwind_o1
            R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o3 
        else if ( jj > 2 .and. jj < N2 - 1 ) then
            L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o3 
            R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o3 
        else if ( jj == N2 - 1 ) then
            L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o3 
            R_W(i) = BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 2 )                            !! backwind_o1
        else if ( jj == N2 ) then
            L_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! upwind_o2
            R_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! upwind_o2
        end if                     
    end do 
#elif ( NL_TERMS == 11 )
    do i = 1,5 
        if ( jj == 1) then 
            L_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 2 )                            !! backwind_o2
            R_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 2 )                            !! backwind_o2
        else if ( jj == 2 ) then
            L_W(i) = UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! upwind_o1
            R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o3 
        else if ( jj == 3 ) then
            L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o3 
            R_W(i) = WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o5
        else if ( jj > 3 .and. jj < N2 - 2 ) then
            L_W(i) = WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o5 
            R_W(i) = WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o5 
        else if ( jj == N2 - 2 ) then
            L_W(i) = WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o5 
            R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o3 
        else if ( jj == N2 - 1 ) then
            L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! HOUC_o3 
            R_W(i) = BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 2 )                            !! backwind_o1
        else if ( jj == N2 ) then
            L_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! upwind_o2
            R_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 2 )                            !! upwind_o2
        end if                     
    end do 
#endif
 
    D_W(:)= 0.
    if ( Lambda(ii,jj,1) > 0) then
       D_W(2)= L_W(2)
       D_W(5)= L_W(5)                             
    elseif ( Lambda(ii,jj,1) < 0) then 
       D_W(2)= R_W(2)
       D_W(5)= R_W(5)      
    else
       D_W(2)= 0.
       D_W(5)= 0.
    end if
    if ( jj == N2 ) then
      Lambda_dRWdx(ii,jj,1) = Lambda(ii,jj,1)*( R(ii,jj,2,2)*D_W(2) + R(ii,jj,2,5)*D_W(5) )
    else 
      Lambda_dRWdx(ii,jj,1) = Lambda(ii,jj,1)*( R(ii,jj,1,2)*D_W(2) + R(ii,jj,1,5)*D_W(5) )  
    end if
  
    if ( Lambda(ii,jj,2) > 0) then
       D_W(2)= L_W(2)
       D_W(5)= L_W(5)                 
           
    elseif ( Lambda(ii,jj,2) < 0) then 
       D_W(2)= R_W(2)
       D_W(5)= R_W(5) 
    else
       D_W(2)= 0.
       D_W(5)= 0.      
    end if
    if ( jj == 1 ) then
      Lambda_dRWdx(ii,jj,2) = Lambda(ii,jj,2)*( R(ii,jj,1,2)*D_W(2) + R(ii,jj,1,5)*D_W(5) ) 
    else 
      Lambda_dRWdx(ii,jj,2) = Lambda(ii,jj,2)*( R(ii,jj,2,2)*D_W(2) + R(ii,jj,2,5)*D_W(5) )  
    end if
  
    if ( Lambda(ii,jj,3) > 0) then
       D_W(1)= L_W(1)
       D_W(4)= L_W(4)                              
    elseif ( Lambda(ii,jj,3) < 0) then 
       D_W(1)= R_W(1)
       D_W(4)= R_W(4)
    else
       D_W(1)= 0.
       D_W(4)= 0.
    end if
    if ( jj == N2 ) then
      Lambda_dRWdx(ii,jj,3) = Lambda(ii,jj,3)*( R(ii,jj,4,1)*D_W(1) + R(ii,jj,4,4)*D_W(4) )  
    else 
      Lambda_dRWdx(ii,jj,3) = Lambda(ii,jj,3)*( R(ii,jj,3,1)*D_W(1) + R(ii,jj,3,4)*D_W(4) )   
    end if
  
    if ( Lambda(ii,jj,4) > 0) then
       D_W(1)= L_W(1)
       D_W(4)= L_W(4)                              
    elseif ( Lambda(ii,jj,4) < 0) then 
       D_W(1)= R_W(1)
       D_W(4)= R_W(4)
    else
       D_W(1)= 0.
       D_W(4)= 0.
    end if
    if ( jj == 1 ) then
      Lambda_dRWdx(ii,jj,4) = Lambda(ii,jj,4)*( R(ii,jj,3,1)*D_W(1) + R(ii,jj,3,4)*D_W(4) )
    else 
      Lambda_dRWdx(ii,jj,4) = Lambda(ii,jj,4)*( R(ii,jj,4,1)*D_W(1) + R(ii,jj,4,4)*D_W(4) )   
    end if
  
    if ( Lambda(ii,jj,5) > 0) then
       D_W(3)= L_W(3)
       D_W(4)= L_W(4)                             
    elseif ( Lambda(ii,jj,5) < 0) then 
       D_W(3)= R_W(3)
       D_W(4)= R_W(4)    
    else
       D_W(3)= 0.
       D_W(4)= 0.
    end if
       Lambda_dRWdx(ii,jj,5) = Lambda(ii,jj,5)*( R(ii,jj,5,3)*D_W(3) + R(ii,jj,5,4)*D_W(4) ) 


        NL(ii,jj,1) = L(ii,jj,1,3)*Lambda_dRWdx(ii,jj,3) + L(ii,jj,1,4)*Lambda_dRWdx(ii,jj,4)
        NL(ii,jj,2) = L(ii,jj,2,1)*Lambda_dRWdx(ii,jj,1) + L(ii,jj,2,2)*Lambda_dRWdx(ii,jj,2)
        NL(ii,jj,3) = L(ii,jj,3,3)*Lambda_dRWdx(ii,jj,3) + L(ii,jj,3,4)*Lambda_dRWdx(ii,jj,4) &
                    + L(ii,jj,3,5)*Lambda_dRWdx(ii,jj,5)
        NL(ii,jj,4) = L(ii,jj,4,3)*Lambda_dRWdx(ii,jj,3) + L(ii,jj,4,4)*Lambda_dRWdx(ii,jj,4)
        NL(ii,jj,5) = L(ii,jj,5,1)*Lambda_dRWdx(ii,jj,1) + L(ii,jj,5,2)*Lambda_dRWdx(ii,jj,2)
    else
        call calcaul_termes_Convections_2(ii, jj, u1, u2, t11, t12, t22, NL)
    end if

      end do
   end do
#elif (DIMENSION_GEO == 3)
   do kk=1, N3
      do jj= 1, N2                 
         do ii= 1, N1


    if (Cent_or_HOUC_B(ii,jj,kk)==0) then
 !!!!-----------l-2(2*NG+N1)---------l-(2*NG+N1)---------l---------l+(2*NG+N1)--------l+2(2*NG+N1)
 !!!!R_W_A2          1                    2              3             4                   5         
      L_W(:)= 0. 
      R_W(:)= 0. 
#if ( NL_TERMS == 10 )
      do i = 1,9
         if ( jj == 1) then 
             L_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! backwind_o2
             R_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! backwind_o2 
         else if ( jj == 2 ) then
             L_W(i) = UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! upwind_o1
             R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o3 
         else if ( jj > 2 .and. jj < N2 - 1 ) then
             L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o3 
             R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o3 
         else if ( jj == N2 - 1 ) then
             L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o3 
             R_W(i) = BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! backwind_o1
         else if ( jj == N2 ) then
             L_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! upwind_o2
             R_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! upwind_o2
         end if                     
      end do 
#elif ( NL_TERMS == 11 )
      do i = 1,9
         if ( jj == 1) then 
             L_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! backwind_o2
             R_W(i) = BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! backwind_o2 
         else if ( jj == 2 ) then
             L_W(i) = UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! upwind_o1
             R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o3 
         else if ( jj == 3 ) then
             L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o3
             R_W(i) = WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o5
         else if ( jj > 3 .and. jj < N2 - 2 ) then
             L_W(i) = WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o5
             R_W(i) = WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o5
         else if ( jj == N2 - 2 ) then
             L_W(i) = WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o5
             R_W(i) = WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o3
         else if ( jj == N2 - 1 ) then
             L_W(i) = WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! HOUC_o3 
             R_W(i) = BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! backwind_o1
         else if ( jj == N2 ) then
             L_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! upwind_o2
             R_W(i) = UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 2, kk )                         !! upwind_o2
         end if                     
      end do 
#endif
 
    D_W(:)= 0.
    if ( Lambda(ii,jj,kk,1) > 0) then
      D_W(2)= L_W(2)
      D_W(7)= L_W(7)                              
  elseif ( Lambda(ii,jj,kk,1) < 0) then 
      D_W(2)= R_W(2)
      D_W(7)= R_W(7)
  else
      D_W(2)= 0.
      D_W(7)= 0. 
  end if
     Lambda_dRWdx(ii,jj,kk,1) = Lambda(ii,jj,kk,1)*( R(ii,jj,kk,1,2)*D_W(2) + R(ii,jj,kk,1,7)*D_W(7) )   
 
  if ( Lambda(ii,jj,kk,2) > 0) then
      D_W(2)= L_W(2)
      D_W(7)= L_W(7)                              
  elseif ( Lambda(ii,jj,kk,2) < 0) then 
      D_W(2)= R_W(2)
      D_W(7)= R_W(7) 
  else
      D_W(2)= 0.
      D_W(7)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,2) = Lambda(ii,jj,kk,2)*( R(ii,jj,kk,2,2)*D_W(2) + R(ii,jj,kk,2,7)*D_W(7) ) 
  
  if ( Lambda(ii,jj,kk,3) > 0) then
      D_W(1)= L_W(1)
      D_W(5)= L_W(5)                              
  elseif ( Lambda(ii,jj,kk,3) < 0) then 
      D_W(1)= R_W(1)
      D_W(5)= R_W(5)  
  else
      D_W(1)= 0.
      D_W(5)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,3) = Lambda(ii,jj,kk,3)*( R(ii,jj,kk,3,1)*D_W(1) + R(ii,jj,kk,3,5)*D_W(5) )   
  
  if ( Lambda(ii,jj,kk,4) > 0) then
      D_W(3)= L_W(3)
      D_W(8)= L_W(8)                          
  elseif ( Lambda(ii,jj,kk,4) < 0) then 
      D_W(3)= R_W(3)
      D_W(8)= R_W(8) 
  else
      D_W(3)= 0.
      D_W(8)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,4) = Lambda(ii,jj,kk,4)* ( R(ii,jj,kk,4,3)*D_W(3) + R(ii,jj,kk,4,8)*D_W(8) )     
  
  if ( Lambda(ii,jj,kk,5) > 0) then
      D_W(1)= L_W(1)
      D_W(5)= L_W(5)                             
  elseif ( Lambda(ii,jj,kk,5) < 0) then 
      D_W(1)= R_W(1)
      D_W(5)= R_W(5)  
  else
      D_W(1)= 0.
      D_W(5)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,5) = Lambda(ii,jj,kk,5)*( R(ii,jj,kk,5,1)*D_W(1) + R(ii,jj,kk,5,5) *D_W(5) )  
  
  if ( Lambda(ii,jj,kk,6) > 0) then
      D_W(3)= L_W(3)
      D_W(8)= L_W(8)                             
  elseif ( Lambda(ii,jj,kk,6) < 0) then 
      D_W(3)= R_W(3)
      D_W(8)= R_W(8)       
  else
      D_W(3)= 0.
      D_W(8)= 0.
  end if
     Lambda_dRWdx(ii,jj,kk,6) = Lambda(ii,jj,kk,6)*( R(ii,jj,kk,6,3)*D_W(3) + R(ii,jj,kk,6,8)*D_W(8) )   
 
  if ( Lambda(ii,jj,kk,7) > 0) then
      D_W(4)= L_W(4)
      D_W(5)= L_W(5)                 
  elseif ( Lambda(ii,jj,kk,7) < 0) then 
      D_W(4)= R_W(4)
      D_W(5)= R_W(5)         
  else
      D_W(4)= 0.
      D_W(5)= 0. 
  end if
     Lambda_dRWdx(ii,jj,kk,7) = Lambda(ii,jj,kk,7)*( R(ii,jj,kk,7,4)*D_W(4) + R(ii,jj,kk,7,5)*D_W(5) )   
 
  if ( Lambda(ii,jj,kk,8) > 0) then
      D_W(5)= L_W(5)
      D_W(6)= L_W(6)
      D_W(7)= L_W(7)
      D_W(8)= L_W(8)                            
  elseif ( Lambda(ii,jj,kk,8) < 0) then 
      D_W(5)= R_W(5)
      D_W(6)= R_W(6)
      D_W(7)= R_W(7)
      D_W(8)= R_W(8)      
  else
      D_W(5)= 0.
      D_W(6)= 0.
      D_W(7)= 0.
      D_W(8)= 0.    
  end if
     Lambda_dRWdx(ii,jj,kk,8) = Lambda(ii,jj,kk,8)*( R(ii,jj,kk,8,5)*D_W(5) + R(ii,jj,kk,8,6)*D_W(6) &
                 + R(ii,jj,kk,8,7)*D_W(7) + R(ii,jj,kk,8,8)*D_W(8) )     
  
  if ( Lambda(ii,jj,kk,9) > 0) then
      D_W(8)= L_W(8)
      D_W(9)= L_W(9)                         
  elseif ( Lambda(ii,jj,kk,9) < 0) then 
      D_W(8)= R_W(8)
      D_W(9)= R_W(9)     
  else
      D_W(8)= 0.
      D_W(9)= 0. 
  end if
     Lambda_dRWdx(ii,jj,kk,9) = Lambda(ii,jj,kk,9)*( R(ii,jj,kk,9,8)*D_W(8) + R(ii,jj,kk,9,9)*D_W(9) )  
 

        NL(ii,jj,kk,1) =L(ii,jj,kk,1,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,1,5)*Lambda_dRWdx(ii,jj,kk,5)
        NL(ii,jj,kk,2) =L(ii,jj,kk,2,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,2,2)*Lambda_dRWdx(ii,jj,kk,2)
        NL(ii,jj,kk,3) =L(ii,jj,kk,3,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,3,6)*Lambda_dRWdx(ii,jj,kk,6)
        NL(ii,jj,kk,4) =L(ii,jj,kk,4,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,4,5)*Lambda_dRWdx(ii,jj,kk,5) &
                       +L(ii,jj,kk,4,7)*Lambda_dRWdx(ii,jj,kk,7)
        NL(ii,jj,kk,5) =L(ii,jj,kk,5,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,5,5)*Lambda_dRWdx(ii,jj,kk,5)
        NL(ii,jj,kk,6) =L(ii,jj,kk,6,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,6,2)*Lambda_dRWdx(ii,jj,kk,2) &
                       +L(ii,jj,kk,6,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,6,4)*Lambda_dRWdx(ii,jj,kk,4) &
                       +L(ii,jj,kk,6,5)*Lambda_dRWdx(ii,jj,kk,5) +L(ii,jj,kk,6,6)*Lambda_dRWdx(ii,jj,kk,6) &
                       +L(ii,jj,kk,6,8)*Lambda_dRWdx(ii,jj,kk,8)
        NL(ii,jj,kk,7) =L(ii,jj,kk,7,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,7,2)*Lambda_dRWdx(ii,jj,kk,2)
        NL(ii,jj,kk,8) =L(ii,jj,kk,8,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,8,6)*Lambda_dRWdx(ii,jj,kk,6)
        NL(ii,jj,kk,9) =L(ii,jj,kk,9,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,9,6)*Lambda_dRWdx(ii,jj,kk,6) &
                       +L(ii,jj,kk,9,9)*Lambda_dRWdx(ii,jj,kk,9)
    else                    
     call calcaul_termes_Convections_2(ii, jj, u1, u2, t11, t12, t22, NL, kk, u3, t13, t23, t33)
    endif

        end do
     end do
   end do
#endif 
    
  end subroutine schema_WENO_A2





#if (DIMENSION_GEO ==3)
  subroutine schema_WENO_A3 (u1, u2, t11, t12, t22, Lambda, L, R, NL, u3, t13, t23, t33)
    implicit none   
    integer  :: i  
    integer  :: ii, jj, kk
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1),intent(in),target:: u1, u2, t11, t12, t22, u3, t13, t23, t33
    real(nk), dimension(1:N1,1:N2,1:N3,1:9,1:9), intent(in)    :: L, R
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    ), intent(in)    :: Lambda
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    ), intent(out)   :: NL
    real(nk), dimension(1:N1,1:N2,1:N3,1:9    )                :: Lambda_dRWdx
    real(nk),dimension(9)  :: D_W, L_W, R_W 
    type pointorW 
    real(nk), dimension(:,:,:), pointer:: var
    end type pointorW 
    type(pointorW), dimension(9) :: WW

    WW(1)%var => u1 ; WW(2)%var => u2 ; WW(3)%var => u3 ; WW(4)%var => t11; WW(5)%var => t12 
    WW(6)%var => t13; WW(7)%var => t22; WW(8)%var => t23; WW(9)%var => t33


   do kk=1, N3
      do jj= 1, N2                 
         do ii= 1, N1

    if (Cent_or_HOUC_C(ii,jj,kk)==0) then
 !!!!-------------l-2*(2*NG+N1)*(2*NG+N2)----------l-(2*NG+N1)*(2*NG+N2)---------l---------l+(2*NG+N1)*(2*NG+N2)-----------l+2*(2*NG+N1)*(2*NG+N2)
 !!!!R_W_A3                  1                              2                    3                    4                                 5
     L_W(:)= 0. 
     R_W(:)= 0. 
#if ( NL_TERMS == 10 )
    do i = 1,9
      if ( kk == 1) then 
          L_W(i) =  BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! backwind_o2
          R_W(i) =  BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! backwind_o2   
      else if ( kk == 2 ) then
          L_W(i) =  UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! centré_o2
          R_W(i) =  WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o3 
      else if ( kk > 2 .and. kk < N3 - 1 ) then
          L_W(i) =  WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o3 
          R_W(i) =  WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o3 
      else if ( kk == N3 - 1 ) then
          L_W(i) =  WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o3 
          R_W(i) =  BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! centré_o2
      else if ( kk == N3 ) then
          L_W(i) =  UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! upwind_o2
          R_W(i) =  UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! upwind_o2
      end if                     
    end do 
#elif ( NL_TERMS == 11 )
    do i = 1,9
      if ( kk == 1) then 
          L_W(i) =  BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! backwind_o2
          R_W(i) =  BACKWIND_2( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! backwind_o2   
      else if ( kk == 2 ) then
          L_W(i) =  UPWIND_1  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! upwind_o1
          R_W(i) =  WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o3 
      else if ( kk == 3 ) then
          L_W(i) =  WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o3
          R_W(i) =  WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o5
      else if ( kk > 3 .and. kk < N3 - 2 ) then
          L_W(i) =  WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o5
          R_W(i) =  WENO_5_R  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o5
      else if ( kk == N3 - 2 ) then
          L_W(i) =  WENO_5_L  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o5
          R_W(i) =  WENO_3_R  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o3 
      else if ( kk == N3 - 1 ) then
          L_W(i) =  WENO_3_L  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! HOUC_o3 
          R_W(i) =  BACKWIND_1( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! backwind_o1
      else if ( kk == N3 ) then
          L_W(i) =  UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! upwind_o2
          R_W(i) =  UPWIND_2  ( 1._nk, WW(i)%var, ii, jj, 3, kk )                            !! upwind_o2
      end if                     
    end do 
#endif
 
   D_W(:)= 0.
   if (Lambda(ii,jj,kk,1) > 0) then
     D_W(3)= L_W(3)
     D_W(9)= L_W(9)                            
 elseif (Lambda(ii,jj,kk,1) < 0) then 
     D_W(3)= R_W(3)
     D_W(9)= R_W(9)       
 else
     D_W(3)= 0.
     D_W(9)= 0.   
 end if
    Lambda_dRWdx(ii,jj,kk,1) =Lambda(ii,jj,kk,1)*(R(ii,jj,kk,1,3)*D_W(3) +R(ii,jj,kk,1,9)*D_W(9) )  
     
 if (Lambda(ii,jj,kk,2) > 0) then
     D_W(3)= L_W(3)
     D_W(9)= L_W(9)                             
 elseif (Lambda(ii,jj,kk,2) < 0) then 
     D_W(3)= R_W(3)
     D_W(9)= R_W(9)       
 else
     D_W(3)= 0.
     D_W(9)= 0.  
 end if
    Lambda_dRWdx(ii,jj,kk,2) =Lambda(ii,jj,kk,2)*(R(ii,jj,kk,2,3)*D_W(3) +R(ii,jj,kk,2,9)*D_W(9) )  
     
 if (Lambda(ii,jj,kk,3) > 0) then
     D_W(1)= L_W(1)
     D_W(6)= L_W(6)                            
 elseif (Lambda(ii,jj,kk,3) < 0) then 
     D_W(1)= R_W(1)
     D_W(6)= R_W(6)         
 else
     D_W(1)= 0.
     D_W(6)= 0.
 end if
    Lambda_dRWdx(ii,jj,kk,3) =Lambda(ii,jj,kk,3)*(R(ii,jj,kk,3,1)*D_W(1) +R(ii,jj,kk,3,6)*D_W(6) ) 
     
 if (Lambda(ii,jj,kk,4) > 0) then
     D_W(2)= L_W(2)
     D_W(8)= L_W(8)                         
 elseif (Lambda(ii,jj,kk,4) < 0) then 
     D_W(2)= R_W(2)
     D_W(8)= R_W(8)        
 else
     D_W(2)= 0.
     D_W(8)= 0. 
 end if
    Lambda_dRWdx(ii,jj,kk,4) =Lambda(ii,jj,kk,4)*(R(ii,jj,kk,4,2)*D_W(2) +R(ii,jj,kk,4,8)*D_W(8) )  
     
 if (Lambda(ii,jj,kk,5) > 0) then
     D_W(1)= L_W(1)
     D_W(6)= L_W(6)                             
 elseif (Lambda(ii,jj,kk,5) < 0) then 
     D_W(1)= R_W(1)
     D_W(6)= R_W(6)        
 else
     D_W(1)= 0.
     D_W(6)= 0.      
 end if
    Lambda_dRWdx(ii,jj,kk,5) =Lambda(ii,jj,kk,5)*(R(ii,jj,kk,5,1)*D_W(1) +R(ii,jj,kk,5,6)*D_W(6) )  
     
 if (Lambda(ii,jj,kk,6) > 0) then
     D_W(2)= L_W(2)
     D_W(8)= L_W(8)                             
 elseif (Lambda(ii,jj,kk,6) < 0) then 
     D_W(2)= R_W(2)
     D_W(8)= R_W(8)      
 else
     D_W(2)= 0.
     D_W(8)= 0. 
 end if
    Lambda_dRWdx(ii,jj,kk,6) =Lambda(ii,jj,kk,6)*(R(ii,jj,kk,6,2)*D_W(2) +R(ii,jj,kk,6,8)*D_W(8) )  
     
 if (Lambda(ii,jj,kk,7) > 0) then
     D_W(4)= L_W(4)
     D_W(6)= L_W(6)                           
 elseif (Lambda(ii,jj,kk,7) < 0) then 
     D_W(4)= R_W(4)
     D_W(6)= R_W(6)  
 else
     D_W(4)= 0.
     D_W(6)= 0.
 end if
    Lambda_dRWdx(ii,jj,kk,7) =Lambda(ii,jj,kk,7)*(R(ii,jj,kk,7,4)*D_W(4) +R(ii,jj,kk,7,6)*D_W(6) )       
     
 if (Lambda(ii,jj,kk,8) > 0) then
     D_W(5)= L_W(5)
     D_W(6)= L_W(6)
     D_W(8)= L_W(8)
     D_W(9)= L_W(9)                            
 elseif (Lambda(ii,jj,kk,8) < 0) then 
     D_W(5)= R_W(5)
     D_W(6)= R_W(6)
     D_W(8)= R_W(8)
     D_W(9)= R_W(9)         
 else
     D_W(5)= 0.
     D_W(6)= 0.
     D_W(8)= 0.
     D_W(9)= 0.   
 end if
    Lambda_dRWdx(ii,jj,kk,8) =Lambda(ii,jj,kk,8)*(R(ii,jj,kk,8,5)*D_W(5) +R(ii,jj,kk,8,6)*D_W(6) &
                                                 +R(ii,jj,kk,8,8)*D_W(8) +R(ii,jj,kk,8,9)*D_W(9) )  
     
 if (Lambda(ii,jj,kk,9) > 0) then
     D_W(7)= L_W(7)
     D_W(8)= L_W(8)                             
 elseif (Lambda(ii,jj,kk,9) < 0) then 
     D_W(7)= R_W(7)
     D_W(8)= R_W(8)     
 else
     D_W(7)= 0.
     D_W(8)= 0.    
 end if
   Lambda_dRWdx(ii,jj,kk,9) =Lambda(ii,jj,kk,9)*(R(ii,jj,kk,9,7)*D_W(7) +R(ii,jj,kk,9,8)*D_W(8) )    


 NL(ii,jj,kk,1) =L(ii,jj,kk,1,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,1,5)*Lambda_dRWdx(ii,jj,kk,5)
 NL(ii,jj,kk,2) =L(ii,jj,kk,2,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,2,6)*Lambda_dRWdx(ii,jj,kk,6)
 NL(ii,jj,kk,3) =L(ii,jj,kk,3,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,3,2)*Lambda_dRWdx(ii,jj,kk,2)
 NL(ii,jj,kk,4) =L(ii,jj,kk,4,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,4,5)*Lambda_dRWdx(ii,jj,kk,5) &
                +L(ii,jj,kk,4,7)*Lambda_dRWdx(ii,jj,kk,7)
 NL(ii,jj,kk,5) =L(ii,jj,kk,5,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,5,2)*Lambda_dRWdx(ii,jj,kk,2) &
                +L(ii,jj,kk,5,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,5,4)*Lambda_dRWdx(ii,jj,kk,4) &
                +L(ii,jj,kk,5,5)*Lambda_dRWdx(ii,jj,kk,5) +L(ii,jj,kk,5,6)*Lambda_dRWdx(ii,jj,kk,6) &
                +L(ii,jj,kk,5,8)*Lambda_dRWdx(ii,jj,kk,8)
 NL(ii,jj,kk,6) =L(ii,jj,kk,6,3)*Lambda_dRWdx(ii,jj,kk,3) +L(ii,jj,kk,6,5)*Lambda_dRWdx(ii,jj,kk,5)
 NL(ii,jj,kk,7) =L(ii,jj,kk,7,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,7,6)*Lambda_dRWdx(ii,jj,kk,6) &
                +L(ii,jj,kk,7,9)*Lambda_dRWdx(ii,jj,kk,9)
 NL(ii,jj,kk,8) =L(ii,jj,kk,8,4)*Lambda_dRWdx(ii,jj,kk,4) +L(ii,jj,kk,8,6)*Lambda_dRWdx(ii,jj,kk,6)
 NL(ii,jj,kk,9) =L(ii,jj,kk,9,1)*Lambda_dRWdx(ii,jj,kk,1) +L(ii,jj,kk,9,2)*Lambda_dRWdx(ii,jj,kk,2)
           else                     
              call calcaul_termes_Convections_3(ii, jj, u1, u2, t11, t12, t22, NL, kk, u3, t13, t23, t33)
           endif
        end do
     end do
  end do   
  
  end subroutine schema_WENO_A3
#endif















 
  subroutine coefs_NL_partie_hyperbolique(d, const, indice, A1_NL, A2_NL, A3_NL)
    implicit none
    integer :: i, j, k
    real(nk),intent(in):: const 
    integer,intent(in) :: indice
#if (DIMENSION_GEO == 2)
    real(nk), dimension(0:N1+1,0:N2+1),intent(inout)  :: d
    real(nk), dimension(1:N1,1:N2,1:5), intent(in)    :: A1_NL, A2_NL
    real(nk),  optional ::  A3_NL
#elif (DIMENSION_GEO == 3)   
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
    real(nk), dimension(1:N1,1:N2,1:N3,1:9),  intent(in)   :: A1_NL, A2_NL, A3_NL
#endif

#if (DIMENSION_GEO == 2)
    do j= 1, N2
       do i=1, N1          
           d(i,j)= d(i,j) -  const*(A1_NL(i,j,indice)  + A2_NL(i,j,indice) )
       end do
    end do    
#elif (DIMENSION_GEO == 3)   
    do k = 1, N3 
       do j= 1, N2          
          do i=1, N1             
              d(i,j,k)= d(i,j,k) - const*(A1_NL(i,j,k,indice)  &
                                        + A2_NL(i,j,k,indice)  + A3_NL(i,j,k,indice) )
          end do
       end do
    end do
#endif 
    
  end subroutine coefs_NL_partie_hyperbolique









  subroutine  Tau_points_ficts
    implicit none
    integer :: i, j, k

#if (DIMENSION_GEO == 2)       
!!!--------------------------------------!! Paroi_bas !!--------------------------------------------
     !Neumann
     var%W(t11 )%var(1:N1,0) =  var%W(t11 )%var(1:N1,2)
     var%W(t12 )%var(1:N1,0) =  var%W(t12 )%var(1:N1,2)
     var%W(t22 )%var(1:N1,0) =  var%W(t22 )%var(1:N1,2)

!!!--------------------------------------!! Paroi_haut !!--------------------------------------------
     !Neumann
     var%W(t11 )%var(1:N1,N2+1) =  var%W(t11 )%var(1:N1,N2-1)
     var%W(t12 )%var(1:N1,N2+1) =  var%W(t12 )%var(1:N1,N2-1)
     var%W(t22 )%var(1:N1,N2+1) =  var%W(t22 )%var(1:N1,N2-1)
   
!!!--------------------------------------!! Paroi_gauche !!--------------------------------------------
     !Neumann
     var%W(t11)%var(0,1:N2) =  var%W(t11)%var(2,1:N2)
     var%W(t12)%var(0,1:N2) =  var%W(t12)%var(2,1:N2)
     var%W(t22)%var(0,1:N2) =  var%W(t22)%var(2,1:N2)       

!!!--------------------------------------!! Paroi_droite !!--------------------------------------------
     !Neumann
     var%W(t11)%var(N1+1,1:N2) =  var%W(t11)%var(N1-1,1:N2)
     var%W(t12)%var(N1+1,1:N2) =  var%W(t12)%var(N1-1,1:N2)
     var%W(t22)%var(N1+1,1:N2) =  var%W(t22)%var(N1-1,1:N2) 

     var%W(t11)%var(0,0) = var%W(t11)%var(2,2)
     var%W(t12)%var(0,0) = var%W(t12)%var(2,2)
     var%W(t22)%var(0,0) = var%W(t22)%var(2,2)

     var%W(t11)%var(0,N2+1) = var%W(t11)%var(2,N2-1)
     var%W(t12)%var(0,N2+1) = var%W(t12)%var(2,N2-1)
     var%W(t22)%var(0,N2+1) = var%W(t22)%var(2,N2-1)

     var%W(t11)%var(N1+1,N2+1) = var%W(t11)%var(N1-1,N2-1)
     var%W(t12)%var(N1+1,N2+1) = var%W(t12)%var(N1-1,N2-1)
     var%W(t22)%var(N1+1,N2+1) = var%W(t22)%var(N1-1,N2-1)

     var%W(t11)%var(N1+1,0) = var%W(t11)%var(N1-1,2)
     var%W(t12)%var(N1+1,0) = var%W(t12)%var(N1-1,2)
     var%W(t22)%var(N1+1,0) = var%W(t22)%var(N1-1,2)

#elif (DIMENSION_GEO == 3)
         
!!!--------------------------------------!! Paroi_bas !!--------------------------------------------
     !Neumann
     var%W(t11 )%var(1:N1,0,1:N3) =  var%W(t11 )%var(1:N1,2,1:N3) 
     var%W(t12 )%var(1:N1,0,1:N3) =  var%W(t12 )%var(1:N1,2,1:N3)
     var%W(t13 )%var(1:N1,0,1:N3) =  var%W(t13 )%var(1:N1,2,1:N3) 
     var%W(t22 )%var(1:N1,0,1:N3) =  var%W(t22 )%var(1:N1,2,1:N3)
     var%W(t23 )%var(1:N1,0,1:N3) =  var%W(t23 )%var(1:N1,2,1:N3)
     var%W(t33 )%var(1:N1,0,1:N3) =  var%W(t33 )%var(1:N1,2,1:N3)

!!!--------------------------------------!! Paroi_haut !!--------------------------------------------
     !Neumann
     var%W(t11 )%var(1:N1,N2+1,1:N3) =  var%W(t11 )%var(1:N1,N2-1,1:N3)
     var%W(t12 )%var(1:N1,N2+1,1:N3) =  var%W(t12 )%var(1:N1,N2-1,1:N3)
     var%W(t13 )%var(1:N1,N2+1,1:N3) =  var%W(t13 )%var(1:N1,N2-1,1:N3)
     var%W(t22 )%var(1:N1,N2+1,1:N3) =  var%W(t22 )%var(1:N1,N2-1,1:N3)
     var%W(t23 )%var(1:N1,N2+1,1:N3) =  var%W(t23 )%var(1:N1,N2-1,1:N3)
     var%W(t33 )%var(1:N1,N2+1,1:N3) =  var%W(t33 )%var(1:N1,N2-1,1:N3)
    
!!!--------------------------------------!! Paroi_gauche !!--------------------------------------------
     !Neumann
     var%W(t11)%var(0,1:N2,1:N3) =  var%W(t11)%var(2,1:N2,1:N3)
     var%W(t12)%var(0,1:N2,1:N3) =  var%W(t12)%var(2,1:N2,1:N3)
     var%W(t13)%var(0,1:N2,1:N3) =  var%W(t13)%var(2,1:N2,1:N3)
     var%W(t22)%var(0,1:N2,1:N3) =  var%W(t22)%var(2,1:N2,1:N3)
     var%W(t23)%var(0,1:N2,1:N3) =  var%W(t23)%var(2,1:N2,1:N3)
     var%W(t33)%var(0,1:N2,1:N3) =  var%W(t33)%var(2,1:N2,1:N3)      

!!!--------------------------------------!! Paroi_droite !!--------------------------------------------
     !Neumann
     var%W(t11)%var(N1+1,1:N2,1:N3) =  var%W(t11)%var(N1-1,1:N2,1:N3)
     var%W(t12)%var(N1+1,1:N2,1:N3) =  var%W(t12)%var(N1-1,1:N2,1:N3)
     var%W(t13)%var(N1+1,1:N2,1:N3) =  var%W(t13)%var(N1-1,1:N2,1:N3)
     var%W(t22)%var(N1+1,1:N2,1:N3) =  var%W(t22)%var(N1-1,1:N2,1:N3)
     var%W(t23)%var(N1+1,1:N2,1:N3) =  var%W(t23)%var(N1-1,1:N2,1:N3)
     var%W(t33)%var(N1+1,1:N2,1:N3) =  var%W(t33)%var(N1-1,1:N2,1:N3) 
   
 !!!--------------------------------------!! Paroi_avant !!--------------------------------------------
     !Neumann
     var%W(t11)%var(1:N1,1:N2,0) =  var%W(t11)%var(1:N1,1:N2,2)
     var%W(t12)%var(1:N1,1:N2,0) =  var%W(t12)%var(1:N1,1:N2,2)
     var%W(t13)%var(1:N1,1:N2,0) =  var%W(t13)%var(1:N1,1:N2,2)   
     var%W(t22)%var(1:N1,1:N2,0) =  var%W(t22)%var(1:N1,1:N2,2)
     var%W(t23)%var(1:N1,1:N2,0) =  var%W(t23)%var(1:N1,1:N2,2)
     var%W(t33)%var(1:N1,1:N2,0) =  var%W(t33)%var(1:N1,1:N2,2)     
 
 !!!--------------------------------------!! Paroi_arriere !!------------------------------------------
     !Neumann
     var%W(t11)%var(1:N1,1:N2,N3+1) =  var%W(t11)%var(1:N1,1:N2,N3-1)
     var%W(t12)%var(1:N1,1:N2,N3+1) =  var%W(t12)%var(1:N1,1:N2,N3-1)
     var%W(t13)%var(1:N1,1:N2,N3+1) =  var%W(t13)%var(1:N1,1:N2,N3-1)
     var%W(t22)%var(1:N1,1:N2,N3+1) =  var%W(t22)%var(1:N1,1:N2,N3-1)
     var%W(t23)%var(1:N1,1:N2,N3+1) =  var%W(t23)%var(1:N1,1:N2,N3-1)
     var%W(t33)%var(1:N1,1:N2,N3+1) =  var%W(t33)%var(1:N1,1:N2,N3-1) 

   var%W(t11 )%var(1:N1,0,0) =  var%W(t11 )%var(1:N1,2,2)
   var%W(t12 )%var(1:N1,0,0) =  var%W(t12 )%var(1:N1,2,2) 
   var%W(t13 )%var(1:N1,0,0) =  var%W(t13 )%var(1:N1,2,2) 
   var%W(t22 )%var(1:N1,0,0) =  var%W(t22 )%var(1:N1,2,2) 
   var%W(t23 )%var(1:N1,0,0) =  var%W(t23 )%var(1:N1,2,2) 
   var%W(t33 )%var(1:N1,0,0) =  var%W(t33 )%var(1:N1,2,2) 

   var%W(t11 )%var(1:N1,N2+1,0) =  var%W(t11 )%var(1:N1,N2-1,2)
   var%W(t12 )%var(1:N1,N2+1,0) =  var%W(t12 )%var(1:N1,N2-1,2)
   var%W(t13 )%var(1:N1,N2+1,0) =  var%W(t13 )%var(1:N1,N2-1,2)
   var%W(t22 )%var(1:N1,N2+1,0) =  var%W(t22 )%var(1:N1,N2-1,2)
   var%W(t23 )%var(1:N1,N2+1,0) =  var%W(t23 )%var(1:N1,N2-1,2)
   var%W(t33 )%var(1:N1,N2+1,0) =  var%W(t33 )%var(1:N1,N2-1,2)

   var%W(t11 )%var(1:N1,0,N3+1) =  var%W(t11 )%var(1:N1,2,N3-1)
   var%W(t12 )%var(1:N1,0,N3+1) =  var%W(t12 )%var(1:N1,2,N3-1) 
   var%W(t13 )%var(1:N1,0,N3+1) =  var%W(t13 )%var(1:N1,2,N3-1)
   var%W(t22 )%var(1:N1,0,N3+1) =  var%W(t22 )%var(1:N1,2,N3-1)
   var%W(t23 )%var(1:N1,0,N3+1) =  var%W(t23 )%var(1:N1,2,N3-1) 
   var%W(t33 )%var(1:N1,0,N3+1) =  var%W(t33 )%var(1:N1,2,N3-1)

   var%W(t11 )%var(1:N1,N2+1,N3+1) =  var%W(t11 )%var(1:N1,N2-1,N3-1)
   var%W(t12 )%var(1:N1,N2+1,N3+1) =  var%W(t12 )%var(1:N1,N2-1,N3-1)
   var%W(t13 )%var(1:N1,N2+1,N3+1) =  var%W(t13 )%var(1:N1,N2-1,N3-1)
   var%W(t22 )%var(1:N1,N2+1,N3+1) =  var%W(t22 )%var(1:N1,N2-1,N3-1)
   var%W(t23 )%var(1:N1,N2+1,N3+1) =  var%W(t23 )%var(1:N1,N2-1,N3-1)
   var%W(t33 )%var(1:N1,N2+1,N3+1) =  var%W(t33 )%var(1:N1,N2-1,N3-1)


   var%W(t11 )%var(0,1:N2,0) =  var%W(t11 )%var(2,1:N2,2)
   var%W(t12 )%var(0,1:N2,0) =  var%W(t12 )%var(2,1:N2,2)
   var%W(t13 )%var(0,1:N2,0) =  var%W(t13 )%var(2,1:N2,2)
   var%W(t22 )%var(0,1:N2,0) =  var%W(t22 )%var(2,1:N2,2)
   var%W(t23 )%var(0,1:N2,0) =  var%W(t23 )%var(2,1:N2,2)
   var%W(t33 )%var(0,1:N2,0) =  var%W(t33 )%var(2,1:N2,2)

   var%W(t11 )%var(N1+1,1:N2,0) =  var%W(t11 )%var(N1-1,1:N2,2)
   var%W(t12 )%var(N1+1,1:N2,0) =  var%W(t12 )%var(N1-1,1:N2,2)
   var%W(t13 )%var(N1+1,1:N2,0) =  var%W(t13 )%var(N1-1,1:N2,2)
   var%W(t22 )%var(N1+1,1:N2,0) =  var%W(t22 )%var(N1-1,1:N2,2)
   var%W(t23 )%var(N1+1,1:N2,0) =  var%W(t23 )%var(N1-1,1:N2,2)
   var%W(t33 )%var(N1+1,1:N2,0) =  var%W(t33 )%var(N1-1,1:N2,2)

   var%W(t11 )%var(0,1:N2,N3+1) =  var%W(t11 )%var(2,1:N2,N3-1)
   var%W(t12 )%var(0,1:N2,N3+1) =  var%W(t12 )%var(2,1:N2,N3-1)
   var%W(t13 )%var(0,1:N2,N3+1) =  var%W(t13 )%var(2,1:N2,N3-1)
   var%W(t22 )%var(0,1:N2,N3+1) =  var%W(t22 )%var(2,1:N2,N3-1)
   var%W(t23 )%var(0,1:N2,N3+1) =  var%W(t23 )%var(2,1:N2,N3-1)
   var%W(t33 )%var(0,1:N2,N3+1) =  var%W(t33 )%var(2,1:N2,N3-1)

   var%W(t11 )%var(N1+1,1:N2,N3+1) =  var%W(t11 )%var(N1-1,1:N2,N3-1)
   var%W(t12 )%var(N1+1,1:N2,N3+1) =  var%W(t12 )%var(N1-1,1:N2,N3-1)
   var%W(t13 )%var(N1+1,1:N2,N3+1) =  var%W(t13 )%var(N1-1,1:N2,N3-1)
   var%W(t22 )%var(N1+1,1:N2,N3+1) =  var%W(t22 )%var(N1-1,1:N2,N3-1)
   var%W(t23 )%var(N1+1,1:N2,N3+1) =  var%W(t23 )%var(N1-1,1:N2,N3-1)
   var%W(t33 )%var(N1+1,1:N2,N3+1) =  var%W(t33 )%var(N1-1,1:N2,N3-1)


   var%W(t11 )%var(0,0,1:N3) =  var%W(t11 )%var(2,2,1:N3)
   var%W(t12 )%var(0,0,1:N3) =  var%W(t12 )%var(2,2,1:N3)
   var%W(t13 )%var(0,0,1:N3) =  var%W(t13 )%var(2,2,1:N3)
   var%W(t22 )%var(0,0,1:N3) =  var%W(t22 )%var(2,2,1:N3)
   var%W(t23 )%var(0,0,1:N3) =  var%W(t23 )%var(2,2,1:N3)
   var%W(t33 )%var(0,0,1:N3) =  var%W(t33 )%var(2,2,1:N3)
   
   var%W(t11 )%var(0,N2+1,1:N3) =  var%W(t11 )%var(2,N2-1,1:N3)
   var%W(t12 )%var(0,N2+1,1:N3) =  var%W(t12 )%var(2,N2-1,1:N3)
   var%W(t13 )%var(0,N2+1,1:N3) =  var%W(t13 )%var(2,N2-1,1:N3)
   var%W(t22 )%var(0,N2+1,1:N3) =  var%W(t22 )%var(2,N2-1,1:N3)
   var%W(t23 )%var(0,N2+1,1:N3) =  var%W(t23 )%var(2,N2-1,1:N3)
   var%W(t33 )%var(0,N2+1,1:N3) =  var%W(t33 )%var(2,N2-1,1:N3)

   var%W(t11 )%var(N1+1,0,1:N3) =  var%W(t11 )%var(N1-1,2,1:N3)
   var%W(t12 )%var(N1+1,0,1:N3) =  var%W(t12 )%var(N1-1,2,1:N3)
   var%W(t13 )%var(N1+1,0,1:N3) =  var%W(t13 )%var(N1-1,2,1:N3)
   var%W(t22 )%var(N1+1,0,1:N3) =  var%W(t22 )%var(N1-1,2,1:N3)
   var%W(t23 )%var(N1+1,0,1:N3) =  var%W(t23 )%var(N1-1,2,1:N3)
   var%W(t33 )%var(N1+1,0,1:N3) =  var%W(t33 )%var(N1-1,2,1:N3)

   var%W(t11 )%var(N1+1,N2+1,1:N3) =  var%W(t11 )%var(N1-1,N2-1,1:N3)
   var%W(t12 )%var(N1+1,N2+1,1:N3) =  var%W(t12 )%var(N1-1,N2-1,1:N3)
   var%W(t13 )%var(N1+1,N2+1,1:N3) =  var%W(t13 )%var(N1-1,N2-1,1:N3)
   var%W(t22 )%var(N1+1,N2+1,1:N3) =  var%W(t22 )%var(N1-1,N2-1,1:N3)
   var%W(t23 )%var(N1+1,N2+1,1:N3) =  var%W(t23 )%var(N1-1,N2-1,1:N3)
   var%W(t33 )%var(N1+1,N2+1,1:N3) =  var%W(t33 )%var(N1-1,N2-1,1:N3)


   var%W(t11 )%var(N1+1,N2+1,N3+1) =  var%W(t11 )%var(N1-1,N2-1,N3-1)
   var%W(t12 )%var(N1+1,N2+1,N3+1) =  var%W(t12 )%var(N1-1,N2-1,N3-1)
   var%W(t13 )%var(N1+1,N2+1,N3+1) =  var%W(t13 )%var(N1-1,N2-1,N3-1)
   var%W(t22 )%var(N1+1,N2+1,N3+1) =  var%W(t22 )%var(N1-1,N2-1,N3-1)
   var%W(t23 )%var(N1+1,N2+1,N3+1) =  var%W(t23 )%var(N1-1,N2-1,N3-1)
   var%W(t33 )%var(N1+1,N2+1,N3+1) =  var%W(t33 )%var(N1-1,N2-1,N3-1)

   var%W(t11 )%var(0,N2+1,N3+1) =  var%W(t11 )%var(2,N2-1,N3-1)
   var%W(t12 )%var(0,N2+1,N3+1) =  var%W(t12 )%var(2,N2-1,N3-1)
   var%W(t13 )%var(0,N2+1,N3+1) =  var%W(t13 )%var(2,N2-1,N3-1)
   var%W(t22 )%var(0,N2+1,N3+1) =  var%W(t22 )%var(2,N2-1,N3-1)
   var%W(t23 )%var(0,N2+1,N3+1) =  var%W(t23 )%var(2,N2-1,N3-1)
   var%W(t33 )%var(0,N2+1,N3+1) =  var%W(t33 )%var(2,N2-1,N3-1)

   var%W(t11 )%var(N1+1,0,N3+1) =  var%W(t11 )%var(N1-1,2,N3-1)
   var%W(t12 )%var(N1+1,0,N3+1) =  var%W(t12 )%var(N1-1,2,N3-1)
   var%W(t13 )%var(N1+1,0,N3+1) =  var%W(t13 )%var(N1-1,2,N3-1)
   var%W(t22 )%var(N1+1,0,N3+1) =  var%W(t22 )%var(N1-1,2,N3-1)
   var%W(t23 )%var(N1+1,0,N3+1) =  var%W(t23 )%var(N1-1,2,N3-1)
   var%W(t33 )%var(N1+1,0,N3+1) =  var%W(t33 )%var(N1-1,2,N3-1)

   var%W(t11 )%var(N1+1,N2+1,0) =  var%W(t11 )%var(N1-1,N2-1,2)
   var%W(t12 )%var(N1+1,N2+1,0) =  var%W(t12 )%var(N1-1,N2-1,2)
   var%W(t13 )%var(N1+1,N2+1,0) =  var%W(t13 )%var(N1-1,N2-1,2)
   var%W(t22 )%var(N1+1,N2+1,0) =  var%W(t22 )%var(N1-1,N2-1,2)
   var%W(t23 )%var(N1+1,N2+1,0) =  var%W(t23 )%var(N1-1,N2-1,2)
   var%W(t33 )%var(N1+1,N2+1,0) =  var%W(t33 )%var(N1-1,N2-1,2)

   var%W(t11 )%var(0,0,N3+1) =  var%W(t11 )%var(2,2,N3-1)
   var%W(t12 )%var(0,0,N3+1) =  var%W(t12 )%var(2,2,N3-1)
   var%W(t13 )%var(0,0,N3+1) =  var%W(t13 )%var(2,2,N3-1)
   var%W(t22 )%var(0,0,N3+1) =  var%W(t22 )%var(2,2,N3-1)
   var%W(t23 )%var(0,0,N3+1) =  var%W(t23 )%var(2,2,N3-1)
   var%W(t33 )%var(0,0,N3+1) =  var%W(t33 )%var(2,2,N3-1)

   var%W(t11 )%var(0,N2+1,0) =  var%W(t11 )%var(2,N2-1,2)
   var%W(t12 )%var(0,N2+1,0) =  var%W(t12 )%var(2,N2-1,2)
   var%W(t13 )%var(0,N2+1,0) =  var%W(t13 )%var(2,N2-1,2)
   var%W(t22 )%var(0,N2+1,0) =  var%W(t22 )%var(2,N2-1,2)
   var%W(t23 )%var(0,N2+1,0) =  var%W(t23 )%var(2,N2-1,2)
   var%W(t33 )%var(0,N2+1,0) =  var%W(t33 )%var(2,N2-1,2)

   var%W(t11 )%var(N1+1,0,0) =  var%W(t11 )%var(N1-1,2,2)
   var%W(t12 )%var(N1+1,0,0) =  var%W(t12 )%var(N1-1,2,2)
   var%W(t13 )%var(N1+1,0,0) =  var%W(t13 )%var(N1-1,2,2)
   var%W(t22 )%var(N1+1,0,0) =  var%W(t22 )%var(N1-1,2,2)
   var%W(t23 )%var(N1+1,0,0) =  var%W(t23 )%var(N1-1,2,2)
   var%W(t33 )%var(N1+1,0,0) =  var%W(t33 )%var(N1-1,2,2)

   var%W(t11 )%var(0,0,0) =  var%W(t11 )%var(2,2,2)
   var%W(t12 )%var(0,0,0) =  var%W(t12 )%var(2,2,2)
   var%W(t13 )%var(0,0,0) =  var%W(t13 )%var(2,2,2)
   var%W(t22 )%var(0,0,0) =  var%W(t22 )%var(2,2,2)
   var%W(t23 )%var(0,0,0) =  var%W(t23 )%var(2,2,2)
   var%W(t33 )%var(0,0,0) =  var%W(t33 )%var(2,2,2)
#endif

  end subroutine Tau_points_ficts


        



  subroutine calcaul_termes_Convections_1(i, j, u1, u2, t11, t12, t22, NL, k, u3, t13, t23, t33)
   implicit none   
   integer, intent(in)         :: i, j
   integer, intent(in),optional:: k
#if ( DIMENSION_GEO == 2 )
    real(nk), dimension(0:N1+1,0:N2+1), intent(in)   :: u1, u2, t11, t12, t22
    real(nk), dimension(1:N1,1:N2,1:5), intent(out)  :: NL
    real(nk), dimension(0:N1+1,0:N2+1), optional     :: u3, t13, t23, t33
#else 
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in)   :: u1, u2, t11, t12, t22, u3, t13, t23, t33
    real(nk), dimension(1:N1,1:N2,1:N3,1:9),   intent(out)  :: NL
#endif
   

#if ( EXP_TREATMENT == 0)
#if (DIMENSION_GEO == 2)

    if (i ==1 ) then 
    NL(i,j,1) =   BACKWIND_2( u1 (i,j)          ,u1 ,i,j,1) -   BACKWIND_2( 1._nk             ,t11,i,j,1) 

    NL(i,j,2) =   BACKWIND_2( u1 (i,j)          ,u2 ,i,j,1) -   BACKWIND_2( 1._nk             ,t12,i,j,1)

    NL(i,j,3) =   BACKWIND_2( u1 (i,j)          ,t11,i,j,1) - 2*BACKWIND_2( 1._nk/Ma2+t11(i,j),u1 ,i,j,1)     

    NL(i,j,4) = - BACKWIND_2( 1._nk/Ma2+t11(i,j),u2 ,i,j,1) +   BACKWIND_2( u1 (i,j)          ,t12,i,j,1)

    NL(i,j,5) =   BACKWIND_2( u1 (i,j)          ,t22,i,j,1) - 2*BACKWIND_2( t12(i,j)          ,u2 ,i,j,1)
    
    else if (i == 2 ) then
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) -   CENTRE_2( 1._nk               ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) -   CENTRE_2( 1._nk               ,t12,i,j,1)
      if (u1 (i,j) > 0._nk) then
       NL(i,j,3) =   UPWIND_1( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   UPWIND_1( u1 (i,j)            ,t12,i,j,1)

       NL(i,j,5) =   UPWIND_1( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)
      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_R( u1 (i,j)            ,t12,i,j,1) 

       NL(i,j,5) =   WENO_3_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)
      end if

    else if ( i == 3 ) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk               ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk               ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_L( u1 (i,j)            ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_R( u1 (i,j)            ,t12,i,j,1)                

       NL(i,j,5) =   WENO_3_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)
      end if 
#elif (CENTRAL_D == 2)
    NL(i,j,1) =   CENTRE_4( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_4( 1._nk               ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_4( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_4( 1._nk               ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_L( u1 (i,j)            ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)            ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_5_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_5_R( u1 (i,j)            ,t12,i,j,1)                

       NL(i,j,5) =   WENO_5_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)            ,u2 ,i,j,1)
      end if 

#endif

    else if ( i > 3 .and. i < N1 - 2 ) then

#if ( CENTRAL_D == 1 )
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk               ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk               ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_L( u1 (i,j)            ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_R( u1 (i,j)            ,t12,i,j,1)                

       NL(i,j,5) =   WENO_3_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)
      end if 
#elif( CENTRAL_D == 2 )
    NL(i,j,1) =   CENTRE_4( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_4( 1._nk               ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_4( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_4( 1._nk               ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_5_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_5_L( u1 (i,j)            ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_5_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)            ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_5_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_5_R( u1 (i,j)            ,t12,i,j,1)                

       NL(i,j,5) =   WENO_5_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)            ,u2 ,i,j,1)
      end if 
#endif

    else if ( i == N1-2 ) then

#if ( CENTRAL_D == 1 )
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk               ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk               ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_L( u1 (i,j)            ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_R( u1 (i,j)            ,t12,i,j,1)                

       NL(i,j,5) =   WENO_3_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)
      end if 

#elif( CENTRAL_D == 2 )
    NL(i,j,1) =   CENTRE_4( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_4( 1._nk               ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_4( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_4( 1._nk               ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_5_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_5_L( u1 (i,j)            ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_5_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)            ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_R( u1 (i,j)            ,t12,i,j,1)                

       NL(i,j,5) =   WENO_3_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)            ,u2 ,i,j,1)
      end if 
#endif

    else if ( i == N1-1 ) then
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk               ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk               ,t12,i,j,1)
      if (u1 (i,j) > 0.) then
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)
 
       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   WENO_3_L( u1 (i,j)            ,t12,i,j,1)
 
       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)            ,u2 ,i,j,1)
      else 
       NL(i,j,3) =   BACKWIND_1( u1 (i,j)          ,t11,i,j,1) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)
 
       NL(i,j,4) = - CENTRE_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   BACKWIND_1( u1 (i,j)          ,t12,i,j,1)
 
       NL(i,j,5) =   BACKWIND_1( u1 (i,j)          ,t22,i,j,1) - 2*CENTRE_2(   t12(i,j)          ,u2 ,i,j,1)
      end if

    else if (i == N1) then

    NL(i,j,1) =   UPWIND_2( u1 (i,j)            ,u1 ,i,j,1) -   UPWIND_2( 1._nk               ,t11,i,j,1) 

    NL(i,j,2) =   UPWIND_2( u1 (i,j)            ,u2 ,i,j,1) -   UPWIND_2( 1._nk               ,t12,i,j,1)

    NL(i,j,3) =   UPWIND_2( u1 (i,j)            ,t11,i,j,1) - 2*UPWIND_2( 1._nk/Ma2 + t11(i,j),u1 ,i,j,1)     

    NL(i,j,4) = - UPWIND_2( 1._nk/Ma2 + t11(i,j),u2 ,i,j,1) +   UPWIND_2( u1 (i,j)            ,t12,i,j,1)

    NL(i,j,5) =   UPWIND_2( u1 (i,j)            ,t22,i,j,1) - 2*UPWIND_2( t12(i,j)            ,u2 ,i,j,1)
    end if 

#elif (DIMENSION_GEO == 3 )

    if (i ==1 ) then 
    NL(i,j,k,1) =   BACKWIND_2( u1 (i,j,k)            ,u1 ,i,j,1,k) -   BACKWIND_2( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   BACKWIND_2( u1 (i,j,k)            ,u2 ,i,j,1,k) -   BACKWIND_2( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   BACKWIND_2( u1 (i,j,k)            ,u3 ,i,j,1,k) -   BACKWIND_2( 1._nk                 ,t13,i,j,1,k)

    NL(i,j,k,4) =-2*BACKWIND_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   BACKWIND_2( u1 (i,j,k)            ,t11,i,j,1,k)

    NL(i,j,k,5) = - BACKWIND_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   BACKWIND_2( u1 (i,j,k)            ,t12,i,j,1,k)

    NL(i,j,k,6) = - BACKWIND_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   BACKWIND_2( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                  + BACKWIND_2( u1 (i,j,k)            ,t13,i,j,1,k)
                  
    NL(i,j,k,7) =-2*BACKWIND_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   BACKWIND_2( u1 (i,j,k)            ,t22,i,j,1,k)

    NL(i,j,k,8) =   BACKWIND_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   BACKWIND_2( t13(i,j,k)            ,u2 ,i,j,1,k) &
                  - BACKWIND_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   BACKWIND_2( u1 (i,j,k)            ,t23,i,j,1,k) 

    NL(i,j,k,9) =-2*BACKWIND_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   BACKWIND_2( u1 (i,j,k)            ,t33,i,j,1,k)
    
    else if (i == 2 ) then
    NL(i,j,k,1) =   CENTRE_2( u1 (i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   CENTRE_2( u1 (i,j,k)            ,u2 ,i,j,1,k) -   CENTRE_2( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   CENTRE_2( u1 (i,j,k)            ,u3 ,i,j,1,k) -   CENTRE_2( 1._nk                 ,t13,i,j,1,k)
      if (u1 (i,j,k) > 0._nk) then
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   UPWIND_1  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   UPWIND_1  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + UPWIND_1  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   UPWIND_1  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   UPWIND_1  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   UPWIND_1  ( u1 (i,j,k)            ,t33,i,j,1,k)
      else 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_R  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t33,i,j,1,k)
      end if

    else if ( i == 3 ) then

#if ( CENTRAL_D == 1 )
    NL(i,j,k,1) =   CENTRE_2( u1 (i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   CENTRE_2( u1 (i,j,k)            ,u2 ,i,j,1,k) -   CENTRE_2( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   CENTRE_2( u1 (i,j,k)            ,u3 ,i,j,1,k) -   CENTRE_2( 1._nk                 ,t13,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_L  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t33,i,j,1,k)

      else 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_R( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t33,i,j,1,k)
      end if 
#elif ( CENTRAL_D == 2 )
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   CENTRE_4( u1 (i,j,k)            ,u3 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t13,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_4( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_L( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_4( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_4( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_4( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t33,i,j,1,k)

      else 
        NL(i,j,k,4) =-2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_4( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_5_R  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_4( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_4( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_4( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t33,i,j,1,k)
      end if 

#endif

    else if ( i > 3 .and. i < N1 - 2 ) then

#if ( CENTRAL_D == 1 )
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   CENTRE_4( u1 (i,j,k)            ,u3 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t13,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_L  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t33,i,j,1,k)

      else 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_R  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t33,i,j,1,k)
      end if 
#elif( CENTRAL_D == 2 )
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   CENTRE_4( u1 (i,j,k)            ,u3 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t13,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_4( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_5_L( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_4( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_4( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_4( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t33,i,j,1,k)
      else 
        NL(i,j,k,4) =-2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_4( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_5_R  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_4( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_4( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_4( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_5_R  ( u1 (i,j,k)            ,t33,i,j,1,k)
      end if 
#endif

    else if ( i == N1-2 ) then

#if ( CENTRAL_D == 1 )
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   CENTRE_4( u1 (i,j,k)            ,u3 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t13,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_L  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t33,i,j,1,k)
      else 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_R  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t33,i,j,1,k)
      end if 

#elif( CENTRAL_D == 2 )
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_4( 1._nk                 ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_4( 1._nk                 ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_4( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_5_L( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_4( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_4( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_4( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_5_L  ( u1 (i,j,k)            ,t33,i,j,1,k)

      else 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_R  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_R  ( u1 (i,j,k)            ,t33,i,j,1,k)
      end if 
#endif

    else if ( i == N1-1 ) then
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   CENTRE_4( u1 (i,j,k)            ,u3 ,i,j,1,k) -   CENTRE_4( 1._nk                 ,t13,i,j,1,k)
      if (u1 (i,j,k) > 0.) then
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + WENO_3_L  ( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   WENO_3_L  ( u1 (i,j,k)            ,t33,i,j,1,k)
      else 
        NL(i,j,k,4) =-2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   BACKWIND_1( u1 (i,j,k)            ,t11,i,j,1,k)

        NL(i,j,k,5) = - CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   BACKWIND_1( u1 (i,j,k)            ,t12,i,j,1,k)
    
        NL(i,j,k,6) = - CENTRE_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                      + BACKWIND_1( u1 (i,j,k)            ,t13,i,j,1,k)
                      
        NL(i,j,k,7) =-2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   BACKWIND_1( u1 (i,j,k)            ,t22,i,j,1,k)
    
        NL(i,j,k,8) =   CENTRE_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2  ( t13(i,j,k)            ,u2 ,i,j,1,k) &
                      - CENTRE_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   BACKWIND_1( u1 (i,j,k)            ,t23,i,j,1,k) 
    
        NL(i,j,k,9) =-2*CENTRE_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   BACKWIND_1( u1 (i,j,k)            ,t33,i,j,1,k)
      end if

    else if (i == N1) then

    NL(i,j,k,1) =   UPWIND_2( u1 (i,j,k)            ,u1 ,i,j,1,k) -   UPWIND_2( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   UPWIND_2( u1 (i,j,k)            ,u2 ,i,j,1,k) -   UPWIND_2( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   UPWIND_2( u1 (i,j,k)            ,u3 ,i,j,1,k) -   UPWIND_2( 1._nk                 ,t13,i,j,1,k)

    NL(i,j,k,4) =-2*UPWIND_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k) +   UPWIND_2( u1 (i,j,k)            ,t11,i,j,1,k)

    NL(i,j,k,5) = - UPWIND_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) +   UPWIND_2( u1 (i,j,k)            ,t12,i,j,1,k)

    NL(i,j,k,6) = - UPWIND_2( t13(i,j,k)            ,u1 ,i,j,1,k) -   UPWIND_2( 1._nk/Ma2 + t11(i,j,k),u3 ,i,j,1,k) &
                  + UPWIND_2( u1 (i,j,k)            ,t13,i,j,1,k)
                  
    NL(i,j,k,7) =-2*UPWIND_2( t12(i,j,k)            ,u2 ,i,j,1,k) +   UPWIND_2( u1 (i,j,k)            ,t22,i,j,1,k)

    NL(i,j,k,8) =   UPWIND_2( t23(i,j,k)            ,u1 ,i,j,1,k) -   UPWIND_2( t13(i,j,k)            ,u2 ,i,j,1,k) &
                  - UPWIND_2( t12(i,j,k)            ,u3 ,i,j,1,k) +   UPWIND_2( u1 (i,j,k)            ,t23,i,j,1,k) 

    NL(i,j,k,9) =-2*UPWIND_2( t13(i,j,k)            ,u3 ,i,j,1,k) +   UPWIND_2( u1 (i,j,k)            ,t33,i,j,1,k)
    end if 
#endif   

#elif ( EXP_TREATMENT == 1 )

#if (DIMENSION_GEO == 2)

    if (i ==1 ) then 
    NL(i,j,1) =   BACKWIND_2( u1 (i,j)            ,u1 ,i,j,1) -   BACKWIND_2( 1._nk/Tau_trans      ,t11,i,j,1) 

    NL(i,j,2) =   BACKWIND_2( u1 (i,j)            ,u2 ,i,j,1) -   BACKWIND_2( 1._nk/Tau_trans      ,t12,i,j,1)

    NL(i,j,3) =   BACKWIND_2( u1 (i,j)            ,t11,i,j,1) - 2*BACKWIND_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)     

    NL(i,j,4) = - BACKWIND_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  BACKWIND_2( u1 (i,j)             ,t12,i,j,1)

    NL(i,j,5) =   BACKWIND_2( u1 (i,j)            ,t22,i,j,1) - 2*BACKWIND_2( t12(i,j)             ,u2 ,i,j,1)
    
    else if (i == 2 ) then
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans     ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans     ,t12,i,j,1)
      if (u1 (i,j) > 0._nk) then
       NL(i,j,3) =   UPWIND_1( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  UPWIND_1( u1 (i,j)             ,t12,i,j,1)

       NL(i,j,5) =   UPWIND_1( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)
      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_3_R( u1 (i,j)             ,t12,i,j,1) 

       NL(i,j,5) =   WENO_3_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)
      end if

    else if ( i == 3 ) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans          ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans          ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_3_L( u1 (i,j)             ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_3_R( u1 (i,j)             ,t12,i,j,1)                

       NL(i,j,5) =   WENO_3_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)
      end if 
#elif (CENTRAL_D == 2)
    NL(i,j,1) =   CENTRE_4( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_4( 1._nk/Tau_trans         ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_4( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_4( 1._nk/Tau_trans         ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_3_L( u1 (i,j)             ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)             ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_5_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_5_R( u1 (i,j)             ,t12,i,j,1)                

       NL(i,j,5) =   WENO_5_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)             ,u2 ,i,j,1)
      end if 

#endif

    else if ( i > 3 .and. i < N1 - 2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans         ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans         ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_3_L( u1 (i,j)             ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_3_R( u1 (i,j)             ,t12,i,j,1)                

       NL(i,j,5) =   WENO_3_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)
      end if 
#elif( CENTRAL_D == 2)
    NL(i,j,1) =   CENTRE_4( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_4( 1._nk/Tau_trans         ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_4( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_4( 1._nk/Tau_trans         ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_5_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_5_L( u1 (i,j)             ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_5_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)             ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_5_R( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_4(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_5_R( u1 (i,j)             ,t12,i,j,1)                

       NL(i,j,5) =   WENO_5_R( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)             ,u2 ,i,j,1)
      end if 
#endif

    else if ( i == N1-2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans         ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans         ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_3_L( u1 (i,j)             ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +   WENO_3_L( u1 (i,j)             ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_3_L( u1 (i,j)             ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)             ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +   WENO_3_R( u1 (i,j)             ,t12,i,j,1)                

       NL(i,j,5) =   WENO_3_R( u1 (i,j)             ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)
      end if 

#elif( CENTRAL_D == 2 )
    NL(i,j,1) =   CENTRE_4( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_4( 1._nk/Tau_trans         ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_4( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_4( 1._nk/Tau_trans         ,t12,i,j,1)
      if (u1 (i,j) > 0.) then 
       NL(i,j,3) =   WENO_5_L( u1 (i,j)             ,t11,i,j,1) - 2*CENTRE_4(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +   WENO_5_L( u1 (i,j)             ,t12,i,j,1)                 

       NL(i,j,5) =   WENO_5_L( u1 (i,j)             ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)             ,u2 ,i,j,1)

      else 
       NL(i,j,3) =   WENO_3_R( u1 (i,j)             ,t11,i,j,1) - 2*CENTRE_4(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)

       NL(i,j,4) = - CENTRE_4(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +   WENO_3_R( u1 (i,j)             ,t12,i,j,1)                

       NL(i,j,5) =   WENO_3_R( u1 (i,j)             ,t22,i,j,1) - 2*CENTRE_4( t12(i,j)             ,u2 ,i,j,1)
      end if 
#endif

    else if ( i == N1-1 ) then
    NL(i,j,1) =   CENTRE_2( u1 (i,j)            ,u1 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans         ,t11,i,j,1)

    NL(i,j,2) =   CENTRE_2( u1 (i,j)            ,u2 ,i,j,1) - CENTRE_2( 1._nk/Tau_trans         ,t12,i,j,1)
      if (u1 (i,j) > 0.) then
       NL(i,j,3) =   WENO_3_L( u1 (i,j)            ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)
 
       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  WENO_3_L( u1 (i,j)             ,t12,i,j,1)
 
       NL(i,j,5) =   WENO_3_L( u1 (i,j)            ,t22,i,j,1) - 2*CENTRE_2( t12(i,j)             ,u2 ,i,j,1)
      else 
       NL(i,j,3) =   BACKWIND_1( u1 (i,j)          ,t11,i,j,1) - 2*CENTRE_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)
 
       NL(i,j,4) = - CENTRE_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  BACKWIND_1( u1 (i,j)           ,t12,i,j,1)
 
       NL(i,j,5) =   BACKWIND_1( u1 (i,j)          ,t22,i,j,1) - 2*CENTRE_2(   t12(i,j)           ,u2 ,i,j,1)
      end if

    else if (i == N1) then

    NL(i,j,1) =   UPWIND_2( u1 (i,j)            ,u1 ,i,j,1) -   UPWIND_2( 1._nk/Tau_trans      ,t11,i,j,1) 

    NL(i,j,2) =   UPWIND_2( u1 (i,j)            ,u2 ,i,j,1) -   UPWIND_2( 1._nk/Tau_trans      ,t12,i,j,1)

    NL(i,j,3) =   UPWIND_2( u1 (i,j)            ,t11,i,j,1) - 2*UPWIND_2(Tau_trans/Ma2+t11(i,j),u1 ,i,j,1)     

    NL(i,j,4) = - UPWIND_2(Tau_trans/Ma2+t11(i,j),u2 ,i,j,1) +  UPWIND_2( u1 (i,j)             ,t12,i,j,1)

    NL(i,j,5) =   UPWIND_2( u1 (i,j)            ,t22,i,j,1) - 2*UPWIND_2( t12(i,j)             ,u2 ,i,j,1)
    end if 

#elif (DIMENSION_GEO == 3)
   

    if (i ==1 ) then 
    NL(i,j,k,1) =   BACKWIND_2( u1 (i,j,k)            ,u1 ,i,j,1,k) - BACKWIND_2( 1._nk                 ,t11,i,j,1,k) 

    NL(i,j,k,2) =   BACKWIND_2( u1 (i,j,k)            ,u2 ,i,j,1,k) - BACKWIND_2( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   BACKWIND_2( u1 (i,j,k)            ,t11,i,j,1,k) - 2*BACKWIND_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)     

    NL(i,j,k,4) = - BACKWIND_2( t12(i,j,k)            ,u1 ,i,j,1,k) -   BACKWIND_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                  + BACKWIND_2( u1 (i,j,k)            ,t12,i,j,1,k)

    NL(i,j,k,5) =   BACKWIND_2( u1 (i,j,k)            ,t22,i,j,1,k) - 2*BACKWIND_2( t12(i,j,k)            ,u2 ,i,j,1,k)
    
    else if (i == 2 ) then
    NL(i,j,k,1) =   CENTRE_2( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_2( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0._nk) then
       NL(i,j,k,3) =   UPWIND_1( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + UPWIND_1( u1 (i,j,k)            ,t12,i,j,1,k)

       NL(i,j,k,5) =   UPWIND_1( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)
      else 
       NL(i,j,k,3) =   WENO_3_R( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k)  -  CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_R( u1 (i,j,k)            ,t12,i,j,1,k) 

       NL(i,j,k,5) =   WENO_3_R( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)
      end if

    else if ( i == 3 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_2( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
       NL(i,j,k,3) =   WENO_3_L( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_L( u1 (i,j,k)            ,t12,i,j,1,k)                 
       NL(i,j,k,5) =   WENO_3_L( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)

      else 
       NL(i,j,k,3) =   WENO_3_R( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k)  -  CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_R( u1 (i,j,k)            ,t12,i,j,1,k)                

       NL(i,j,k,5) =   WENO_3_R( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)
      end if 
#elif (CENTRAL_D == 2)
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_4( 1._nk                     ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_4( 1._nk                     ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
       NL(i,j,k,3) =   WENO_3_L( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_4( t12(i,j,k)            ,u1 ,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_L( u1 (i,j,k)            ,t12,i,j,1,k)                 

       NL(i,j,k,5) =   WENO_3_L( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k)

      else 
       NL(i,j,k,3) =   WENO_5_R( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_4( t12(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_5_R( u1 (i,j,k)            ,t12,i,j,1,k)                

       NL(i,j,k,5) =   WENO_5_R( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k)
      end if 

#endif

    else if ( i > 3 .and. i < N1 - 2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_2( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
       NL(i,j,k,3) =   WENO_3_L( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_L( u1 (i,j,k)            ,t12,i,j,1,k)                 

       NL(i,j,k,5) =   WENO_3_L( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)

      else 
       NL(i,j,k,3) =   WENO_3_R( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_R( u1 (i,j,k)            ,t12,i,j,1,k)                

       NL(i,j,k,5) =   WENO_3_R( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)
      end if 
#elif( CENTRAL_D == 2)
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_4( 1._nk                     ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_4( 1._nk                     ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
       NL(i,j,k,3) =   WENO_5_L( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_4( t12(i,j,k)            ,u1 ,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) &
                     + WENO_5_L( u1 (i,j,k)            ,t12,i,j,1,k)                 

       NL(i,j,k,5) =   WENO_5_L( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k)

      else 
       NL(i,j,k,3) =   WENO_5_R( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_4( t12(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_5_R( u1 (i,j,k)            ,t12,i,j,1,k)                

       NL(i,j,k,5) =   WENO_5_R( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k)
      end if 
#endif

    else if ( i == N1-2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_2( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
       NL(i,j,k,3) =   WENO_3_L( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_L( u1 (i,j,k)            ,t12,i,j,1,k)                 

       NL(i,j,k,5) =   WENO_3_L( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)

      else 
       NL(i,j,k,3) =   WENO_3_R( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_R( u1 (i,j,k)            ,t12,i,j,1,k)                

       NL(i,j,k,5) =   WENO_3_R( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)
      end if 

#elif( CENTRAL_D == 2 )
    NL(i,j,k,1) =   CENTRE_4( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_4( 1._nk                     ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_4( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_4( 1._nk                     ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0.) then 
       NL(i,j,k,3) =   WENO_5_L( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_4( t12(i,j,k)            ,u1 ,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_5_L( u1 (i,j,k)            ,t12,i,j,1,k)                 

       NL(i,j,k,5) =   WENO_5_L( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k)

      else 
       NL(i,j,k,3) =   WENO_3_R( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)

       NL(i,j,k,4) = - CENTRE_4( t12(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_4( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                     + WENO_3_R( u1 (i,j,k)            ,t12,i,j,1,k)                

       NL(i,j,k,5) =   WENO_3_R( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_4( t12(i,j,k)            ,u2 ,i,j,1,k)
      end if 
#endif

    else if ( i == N1-1 ) then
    NL(i,j,k,1) =   CENTRE_2( u1 (i,j,k)            ,u1 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t11,i,j,1,k)

    NL(i,j,k,2) =   CENTRE_2( u1 (i,j,k)            ,u2 ,i,j,1,k) - CENTRE_2( 1._nk                     ,t12,i,j,1,k)
      if (u1 (i,j,k) > 0.) then
       NL(i,j,k,3) =   WENO_3_L( u1 (i,j,k)            ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)
 
       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) &
                     + WENO_3_L( u1 (i,j,k)            ,t12,i,j,1,k)
 
       NL(i,j,k,5) =   WENO_3_L( u1 (i,j,k)            ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)
      else 
       NL(i,j,k,3) =   BACKWIND_1( u1 (i,j,k)          ,t11,i,j,1,k) - 2*CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)
 
       NL(i,j,k,4) = - CENTRE_2( t12(i,j,k)            ,u1 ,i,j,1,k) -   CENTRE_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k) &
                     + BACKWIND_1( u1 (i,j,k)          ,t12,i,j,1,k)
 
       NL(i,j,k,5) =   BACKWIND_1( u1 (i,j,k)          ,t22,i,j,1,k) - 2*CENTRE_2( t12(i,j,k)            ,u2 ,i,j,1,k)
      end if

    else if (i == N1) then

    NL(i,j,k,1) =   UPWIND_2( u1 (i,j,k)            ,u1 ,i,j,1,k) -   UPWIND_2( 1._nk                  ,t11,i,j,1,k) 

    NL(i,j,k,2) =   UPWIND_2( u1 (i,j,k)            ,u2 ,i,j,1,k) -   UPWIND_2( 1._nk                 ,t12,i,j,1,k)

    NL(i,j,k,3) =   UPWIND_2( u1 (i,j,k)            ,t11,i,j,1,k) - 2*UPWIND_2( 1._nk/Ma2 + t11(i,j,k),u1 ,i,j,1,k)     

    NL(i,j,k,4) = - UPWIND_2( t12(i,j,k)            ,u1 ,i,j,1,k) -   UPWIND_2( 1._nk/Ma2 + t11(i,j,k),u2 ,i,j,1,k)  &
                  + UPWIND_2( u1 (i,j,k)            ,t12,i,j,1,k)

    NL(i,j,k,5) =   UPWIND_2( u1 (i,j,k)            ,t22,i,j,1,k) - 2*UPWIND_2( t12(i,j,k)            ,u2 ,i,j,1,k)
    end if 
#endif   
#endif
 end subroutine calcaul_termes_Convections_1

        




 subroutine calcaul_termes_Convections_2(i, j, u1, u2, t11, t12, t22, NL, k, u3, t13, t23, t33)
   implicit none   
   integer, intent(in)         :: i, j
   integer, intent(in),optional:: k
#if ( DIMENSION_GEO == 2 )
    real(nk), dimension(0:N1+1,0:N2+1), intent(in)   :: u1, u2, t11, t12, t22
    real(nk), dimension(1:N1,1:N2,1:5), intent(out)  :: NL
    real(nk), dimension(0:N1+1,0:N2+1), optional     :: u3, t13, t23, t33
#else 
    real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in)   :: u1, u2, t11, t12, t22, u3, t13, t23, t33
    real(nk), dimension(1:N1,1:N2,1:N3,1:9),   intent(out)  :: NL
#endif

#if ( EXP_TREATMENT == 0 )
#if (DIMENSION_GEO == 2)
    if (j == 1) then
    NL(i,j,1) =   BACKWIND_2( u2 (i,j)            ,u1 ,i,j,2) - BACKWIND_2( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   BACKWIND_2( u2 (i,j)            ,u2 ,i,j,2) - BACKWIND_2( 1._nk                   ,t22,i,j,2)

    NL(i,j,3) =   BACKWIND_2( u2 (i,j)            ,t11,i,j,2) - 2*BACKWIND_2( t12(i,j)            ,u1 ,i,j,2)
                      
    NL(i,j,4) =   BACKWIND_2( u2 (i,j)            ,t12,i,j,2) -   BACKWIND_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)
                          
    NL(i,j,5) =   BACKWIND_2( u2 (i,j)            ,t22,i,j,2) - 2*BACKWIND_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)
    else if (j == 2 ) then

    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk                   ,t22,i,j,2)
      if (u2 (i,j) > 0.) then 

       NL(i,j,3) =   UPWIND_1( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)
                      
       NL(i,j,4) =   UPWIND_1( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)   
                          
       NL(i,j,5) =   UPWIND_1( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)
      else 
       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)
                      
       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)  
                          
       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)
      end if

    else if ( j == 3 ) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk                   ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)  

       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)  

       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)                 
      end if
                      
#elif ( CENTRAL_D == 2 )
    NL(i,j,1) =   CENTRE_4( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_4( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_4( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_4( 1._nk                   ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)  

       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_5_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_5_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_5_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)                 
      end if

#endif

    else if ( j > 3 .and. j < N2-2) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk                   ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)                 
      end if
#elif (CENTRAL_D == 2)   
    NL(i,j,1) =   CENTRE_4( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_4( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_4( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_4( 1._nk                   ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_5_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_5_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)  

       NL(i,j,5) =   WENO_5_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_5_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_5_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_5_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)                 
      end if           

#endif
    
    else if ( j == N2-2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk                   ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -    CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)

       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)                 
      end if
#elif (CENTRAL == 2)              
    NL(i,j,1) =   CENTRE_4( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_4( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_4( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_4( 1._nk                   ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_5_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_5_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_5_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)            ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)

       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)                 
      end if
#endif

    else if ( j == N2-1 ) then
    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk                   ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk                   ,t22,i,j,2)
      if (u2 (i,j) > 0.) then 
 
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)
                       
       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2) 
                           
       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)               
      else 
       NL(i,j,3) =   BACKWIND_1( u2 (i,j)          ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)
                       
       NL(i,j,4) =   BACKWIND_1( u2 (i,j)          ,t12,i,j,2) -   CENTRE_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)
                           
       NL(i,j,5) =   BACKWIND_1( u2 (i,j)          ,t22,i,j,2) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)
      end if
    else if (j == N2) then
    NL(i,j,1) =   UPWIND_2( u2 (i,j)            ,u1 ,i,j,2) -   UPWIND_2( 1._nk               ,t12,i,j,2)
 
    NL(i,j,2) =   UPWIND_2( u2 (i,j)            ,u2 ,i,j,2) -   UPWIND_2( 1._nk               ,t22,i,j,2)

    NL(i,j,3) =   UPWIND_2( u2 (i,j)            ,t11,i,j,2) - 2*UPWIND_2( t12(i,j)            ,u1 ,i,j,2)
                      
    NL(i,j,4) =   UPWIND_2( u2 (i,j)            ,t12,i,j,2) -   UPWIND_2( 1._nk/Ma2 + t22(i,j),u1 ,i,j,2)
                          
    NL(i,j,5) =   UPWIND_2( u2 (i,j)            ,t22,i,j,2) - 2*UPWIND_2( 1._nk/Ma2 + t22(i,j),u2 ,i,j,2)

    end if

#elif (DIMENSION_GEO == 3)


    if (j == 1) then
    NL(i,j,k,1) =   BACKWIND_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   BACKWIND_2( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   BACKWIND_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   BACKWIND_2( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   BACKWIND_2( u2 (i,j,k)            ,u3 ,i,j,2,k) -   BACKWIND_2( 1._nk                 ,t23,i,j,2,k)

    NL(i,j,k,4) =-2*BACKWIND_2( t12(i,j,k)            ,u1 ,i,j,2,k) +   BACKWIND_2( u2 (i,j,k)            ,t11,i,j,2,k)
                    
    NL(i,j,k,5) = - BACKWIND_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   BACKWIND_2( u2 (i,j,k)            ,t12,i,j,2,k)

    NL(i,j,k,6) = - BACKWIND_2( t23(i,j,k)            ,u1 ,i,j,2,k) -   BACKWIND_2( t12(i,j,k)            ,u3 ,i,j,2,k) &
                  + BACKWIND_2( u2 (i,j,k)            ,t13,i,j,2,k)

    NL(i,j,k,7) =-2*BACKWIND_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   BACKWIND_2( u2 (i,j,k)            ,t22,i,j,2,k)

    NL(i,j,k,8) = - BACKWIND_2( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   BACKWIND_2( u2 (i,j,k)            ,t23,i,j,2,k)

    NL(i,j,k,9) =-2*BACKWIND_2( t23(i,j,k)            ,u3 ,i,j,2,k) +   BACKWIND_2( u2 (i,j,k)            ,t33,i,j,2,k)

    else if (j == 2 ) then

    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   CENTRE_2( u2 (i,j,k)            ,u3 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t23,i,j,2,k)
      if (u2 (i,j,k) > 0.) then 

        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   UPWIND_1  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   UPWIND_1  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + UPWIND_1  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   UPWIND_1  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   UPWIND_1  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   UPWIND_1  ( u2 (i,j,k)            ,t33,i,j,2,k)
    
      else 
        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_R  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t33,i,j,2,k)
      end if

    else if ( j == 3 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   CENTRE_2( u2 (i,j,k)            ,u3 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t23,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_L  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t33,i,j,2,k)
      else 

        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_R  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t33,i,j,2,k)
      end if
                      
#elif ( CENTRAL_D == 2 )
    NL(i,j,k,1) =   CENTRE_4( u2 (i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_4( u2 (i,j,k)            ,u2 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   CENTRE_4( u2 (i,j,k)            ,u3 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t23,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
        NL(i,j,k,4) =-2*CENTRE_4  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_4  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_L  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_4  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t33,i,j,2,k)
      else 

        NL(i,j,k,4) =-2*CENTRE_4  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_4  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_5_R  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_4  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t33,i,j,2,k)
      end if

#endif

    else if ( j > 3 .and. j < N2-2) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   CENTRE_2( u2 (i,j,k)            ,u3 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t23,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_L  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t33,i,j,2,k)
      else 

        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_R  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t33,i,j,2,k)
      end if
#elif (CENTRAL_D == 2)   
    NL(i,j,k,1) =   CENTRE_4( u2 (i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_4( u2 (i,j,k)            ,u2 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   CENTRE_4( u2 (i,j,k)            ,u3 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t23,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
        NL(i,j,k,4) =-2*CENTRE_4  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_4  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_5_L  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_4  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t33,i,j,2,k)
      else 

        NL(i,j,k,4) =-2*CENTRE_4  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_4  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_5_R  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_4  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_5_R  ( u2 (i,j,k)            ,t33,i,j,2,k)
      end if           

#endif
    
    else if ( j == N2-2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   CENTRE_2( u2 (i,j,k)            ,u3 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t23,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_L  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t33,i,j,2,k)
      else 

        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_R  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t33,i,j,2,k)  
      end if
#elif (CENTRAL == 2)              
    NL(i,j,k,1) =   CENTRE_4( u2 (i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_4( u2 (i,j,k)            ,u2 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   CENTRE_4( u2 (i,j,k)            ,u3 ,i,j,2,k) -   CENTRE_4( 1._nk                 ,t23,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
        NL(i,j,k,4) =-2*CENTRE_4  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_4  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_5_L  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_4  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_5_L  ( u2 (i,j,k)            ,t33,i,j,2,k)
      else 

        NL(i,j,k,4) =-2*CENTRE_4  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_4  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_R  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_4  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_4  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_R  ( u2 (i,j,k)            ,t33,i,j,2,k)
      end if
#endif

    else if ( j == N2-1 ) then
    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   CENTRE_2( u2 (i,j,k)            ,u3 ,i,j,2,k) -   CENTRE_2( 1._nk                 ,t23,i,j,2,k)
      if (u2 (i,j,k) > 0.) then 
 
        NL(i,j,k,4) =-2*CENTRE_2  ( t12(i,j,k)            ,u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - CENTRE_2  ( t23(i,j,k)            ,u1 ,i,j,2,k) -   CENTRE_2  ( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + WENO_3_L  ( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - CENTRE_2  ( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*CENTRE_2  ( t23(i,j,k)            ,u3 ,i,j,2,k) +   WENO_3_L  ( u2 (i,j,k)            ,t33,i,j,2,k)
      else 
        NL(i,j,k,4) =-2*BACKWIND_1( t12(i,j,k)            ,u1 ,i,j,2,k) +   BACKWIND_1( u2 (i,j,k)            ,t11,i,j,2,k)
                    
        NL(i,j,k,5) = - BACKWIND_1( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   BACKWIND_1( u2 (i,j,k)            ,t12,i,j,2,k)
    
        NL(i,j,k,6) = - BACKWIND_1( t23(i,j,k)            ,u1 ,i,j,2,k) -   BACKWIND_1( t12(i,j,k)            ,u3 ,i,j,2,k) &
                      + BACKWIND_1( u2 (i,j,k)            ,t13,i,j,2,k)
    
        NL(i,j,k,7) =-2*BACKWIND_1( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   BACKWIND_1( u2 (i,j,k)            ,t22,i,j,2,k)
    
        NL(i,j,k,8) = - BACKWIND_1( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   BACKWIND_1( u2 (i,j,k)            ,t23,i,j,2,k)
    
        NL(i,j,k,9) =-2*BACKWIND_1( t23(i,j,k)            ,u3 ,i,j,2,k) +   BACKWIND_1( u2 (i,j,k)            ,t33,i,j,2,k)
      end if
    else if (j == N2) then
    NL(i,j,k,1) =   UPWIND_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   UPWIND_2( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   UPWIND_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   UPWIND_2( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   UPWIND_2( u2 (i,j,k)            ,u3 ,i,j,2,k) -   UPWIND_2( 1._nk                 ,t23,i,j,2,k)

    NL(i,j,k,4) =-2*UPWIND_2( t12(i,j,k)            ,u1 ,i,j,2,k) +   UPWIND_2( u2 (i,j,k)            ,t11,i,j,2,k)
                    
    NL(i,j,k,5) = - UPWIND_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) +   UPWIND_2( u2 (i,j,k)            ,t12,i,j,2,k)

    NL(i,j,k,6) = - UPWIND_2( t23(i,j,k)            ,u1 ,i,j,2,k) -   UPWIND_2( t12(i,j,k)            ,u3 ,i,j,2,k) &
                  + UPWIND_2( u2 (i,j,k)            ,t13,i,j,2,k)

    NL(i,j,k,7) =-2*UPWIND_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k) +   UPWIND_2( u2 (i,j,k)            ,t22,i,j,2,k)

    NL(i,j,k,8) = - UPWIND_2( 1._nk/Ma2 + t22(i,j,k),u3 ,i,j,2,k) +   UPWIND_2( u2 (i,j,k)            ,t23,i,j,2,k)

    NL(i,j,k,9) =-2*UPWIND_2( t23(i,j,k)            ,u3 ,i,j,2,k) +   UPWIND_2( u2 (i,j,k)            ,t33,i,j,2,k)


    end if
   
             
#endif   

#elif( EXP_TREATMENT == 1 )

#if (DIMENSION_GEO == 2)
    if (j == 1) then
    NL(i,j,1) =   BACKWIND_2( u2 (i,j)            ,u1 ,i,j,2) - BACKWIND_2( 1._nk/Tau_trans        ,t12,i,j,2)

    NL(i,j,2) =   BACKWIND_2( u2 (i,j)            ,u2 ,i,j,2) - BACKWIND_2( 1._nk/Tau_trans        ,t22,i,j,2)

    NL(i,j,3) =   BACKWIND_2( u2 (i,j)            ,t11,i,j,2) - 2*BACKWIND_2( t12(i,j)             ,u1 ,i,j,2)
                      
    NL(i,j,4) =   BACKWIND_2( u2 (i,j)            ,t12,i,j,2) -   BACKWIND_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)
                          
    NL(i,j,5) =   BACKWIND_2( u2 (i,j)            ,t22,i,j,2) - 2*BACKWIND_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)
    else if (j == 2 ) then

    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t22,i,j,2)
      if (u2 (i,j) > 0.) then 

       NL(i,j,3) =   UPWIND_1( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)
                      
       NL(i,j,4) =   UPWIND_1( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2) 
                          
       NL(i,j,5) =   UPWIND_1( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)
      else 
       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)
                      
       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)
                          
       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)
      end if

    else if ( j == 3 ) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)  

       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)             ,u1 ,i,j,2)
 
       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)                 
      end if
                      
#elif ( CENTRAL_D == 2 )
    NL(i,j,1) =   CENTRE_4( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_4( 1._nk/Tau_trans         ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_4( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_4( 1._nk/Tau_trans         ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_5_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_5_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)

       NL(i,j,5) =   WENO_5_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)                 
      end if

#endif

    else if ( j > 3 .and. j < N2-2) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)

       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)

       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)                 
      end if
#elif (CENTRAL_D == 2)   
    NL(i,j,1) =   CENTRE_4( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_4( 1._nk/Tau_trans         ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_4( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_4( 1._nk/Tau_trans         ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_5_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_5_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)

       NL(i,j,5) =   WENO_5_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_5_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_5_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)

       NL(i,j,5) =   WENO_5_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)                 
      end if           

#endif
    
    else if ( j == N2-2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)

       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)                 
      end if
#elif (CENTRAL == 2)              
    NL(i,j,1) =   CENTRE_4( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_4( 1._nk/Tau_trans         ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_4( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_4( 1._nk/Tau_trans         ,t22,i,j,2)
      if (u2 (i,j) > 0. ) then
       NL(i,j,3) =   WENO_5_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_5_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_5_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)            
      else 

       NL(i,j,3) =   WENO_3_R( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_4( t12(i,j)             ,u1 ,i,j,2)

       NL(i,j,4) =   WENO_3_R( u2 (i,j)            ,t12,i,j,2) -   CENTRE_4(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2) 

       NL(i,j,5) =   WENO_3_R( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_4(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)                 
      end if
#endif

    else if ( j == N2-1 ) then
    NL(i,j,1) =   CENTRE_2( u2 (i,j)            ,u1 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t12,i,j,2)

    NL(i,j,2) =   CENTRE_2( u2 (i,j)            ,u2 ,i,j,2) - CENTRE_2( 1._nk/Tau_trans         ,t22,i,j,2)
      if (u2 (i,j) > 0.) then 
 
       NL(i,j,3) =   WENO_3_L( u2 (i,j)            ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)             ,u1 ,i,j,2)
                       
       NL(i,j,4) =   WENO_3_L( u2 (i,j)            ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)
                           
       NL(i,j,5) =   WENO_3_L( u2 (i,j)            ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)               
      else 
       NL(i,j,3) =   BACKWIND_1( u2 (i,j)          ,t11,i,j,2) - 2*CENTRE_2( t12(i,j)            ,u1 ,i,j,2)
                       
       NL(i,j,4) =   BACKWIND_1( u2 (i,j)          ,t12,i,j,2) -   CENTRE_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)  
                           
       NL(i,j,5) =   BACKWIND_1( u2 (i,j)          ,t22,i,j,2) - 2*CENTRE_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)
      end if
    else if (j == N2) then
    NL(i,j,1) =   UPWIND_2( u2 (i,j)            ,u1 ,i,j,2) -   UPWIND_2( 1._nk/Tau_trans      ,t12,i,j,2)
 
    NL(i,j,2) =   UPWIND_2( u2 (i,j)            ,u2 ,i,j,2) -   UPWIND_2( 1._nk/Tau_trans      ,t22,i,j,2)

    NL(i,j,3) =   UPWIND_2( u2 (i,j)            ,t11,i,j,2) - 2*UPWIND_2( t12(i,j)             ,u1 ,i,j,2)
                      
    NL(i,j,4) =   UPWIND_2( u2 (i,j)            ,t12,i,j,2) -   UPWIND_2(Tau_trans/Ma2+t22(i,j),u1 ,i,j,2)
                          
    NL(i,j,5) =   UPWIND_2( u2 (i,j)            ,t22,i,j,2) - 2*UPWIND_2(Tau_trans/Ma2+t22(i,j),u2 ,i,j,2)

    end if

#elif (DIMENSION_GEO == 3)



    if (j == 1) then
    NL(i,j,k,1) =   BACKWIND_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   BACKWIND_2( 1._nk                 ,t12,i,j,2,k)

    NL(i,j,k,2) =   BACKWIND_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   BACKWIND_2( 1._nk                 ,t22,i,j,2,k)

    NL(i,j,k,3) =   BACKWIND_2( u2 (i,j,k)            ,t11,i,j,2,k) - 2*BACKWIND_2( t12(i,j,k)            ,u1 ,i,j,2,k)
                      
    NL(i,j,k,4) =   BACKWIND_2( u2 (i,j,k)            ,t12,i,j,2,k) -   BACKWIND_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k) &
                  - BACKWIND_2( t12(i,j,k)            ,u2 ,i,j,2,k)
                          
    NL(i,j,k,5) =   BACKWIND_2( u2 (i,j,k)            ,t22,i,j,2,k) - 2*BACKWIND_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)
    else if (j == 2 ) then

    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t22,i,j,2,k)
      if (u2 (i,j,k) > 0.) then 

       NL(i,j,k,3) =   UPWIND_1( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)
                      
       NL(i,j,k,4) =   UPWIND_1( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)
                          
       NL(i,j,k,5) =   UPWIND_1( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)
      else 
       NL(i,j,k,3) =   WENO_3_R( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)
                      
       NL(i,j,k,4) =   WENO_3_R( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)
                          
       NL(i,j,k,5) =   WENO_3_R( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)
      end if

    else if ( j == 3 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t22,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
       NL(i,j,k,3) =   WENO_3_L( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_3_L( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_3_L( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)            
      else 

       NL(i,j,k,3) =   WENO_3_R( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_3_R( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_3_R( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)                 
      end if
                      
#elif ( CENTRAL_D == 2 )
    NL(i,j,k,1) =   CENTRE_4( u2 (i,j,k)            ,u1 ,i,j,2,k) - CENTRE_4( 1._nk                     ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_4( u2 (i,j,k)            ,u2 ,i,j,2,k) - CENTRE_4( 1._nk                     ,t22,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
       NL(i,j,k,3) =   WENO_3_L( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_4( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_3_L( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_4( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_3_L( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)            
      else 

       NL(i,j,k,3) =   WENO_5_R( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_4( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_5_R( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_4( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_5_R( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)                 
      end if

#endif

    else if ( j > 3 .and. j < N2-2) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t22,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
       NL(i,j,k,3) =   WENO_3_L( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_3_L( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_3_L( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)            
      else 

       NL(i,j,k,3) =   WENO_3_R( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_3_R( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_3_R( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)                 
      end if
#elif (CENTRAL_D == 2)   
    NL(i,j,k,1) =   CENTRE_4( u2 (i,j,k)            ,u1 ,i,j,2,k) - CENTRE_4( 1._nk                     ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_4( u2 (i,j,k)            ,u2 ,i,j,2,k) - CENTRE_4( 1._nk                     ,t22,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
       NL(i,j,k,3) =   WENO_5_L( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_4( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_5_L( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_4( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_5_L( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)            
      else 

       NL(i,j,k,3) =   WENO_5_R( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_4( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_5_R( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_4( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_5_R( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)                 
      end if           

#endif
    
    else if ( j == N2-2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t22,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
       NL(i,j,k,3) =   WENO_3_L( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_3_L( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_3_L( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)            
      else 

       NL(i,j,k,3) =   WENO_3_R( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_3_R( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_3_R( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)                 
      end if
#elif (CENTRAL == 2)              
    NL(i,j,k,1) =   CENTRE_4( u2 (i,j,k)            ,u1 ,i,j,2,k) - CENTRE_4( 1._nk                     ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_4( u2 (i,j,k)            ,u2 ,i,j,2,k) - CENTRE_4( 1._nk                     ,t22,i,j,2,k)
      if (u2 (i,j,k) > 0. ) then
       NL(i,j,k,3) =   WENO_5_L( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_4( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_5_L( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_4( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_5_L( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)            
      else 

       NL(i,j,k,3) =   WENO_3_R( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_4( t12(i,j,k)            ,u1 ,i,j,2,k)

       NL(i,j,k,4) =   WENO_3_R( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_4( t12(i,j,k)            ,u2 ,i,j,2,k)

       NL(i,j,k,5) =   WENO_3_R( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_4( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)                 
      end if
#endif

    else if ( j == N2-1 ) then
    NL(i,j,k,1) =   CENTRE_2( u2 (i,j,k)            ,u1 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t12,i,j,2,k)

    NL(i,j,k,2) =   CENTRE_2( u2 (i,j,k)            ,u2 ,i,j,2,k) - CENTRE_2( 1._nk                     ,t22,i,j,2,k)
      if (u2 (i,j,k) > 0.) then 
 
       NL(i,j,k,3) =   WENO_3_L( u2 (i,j,k)            ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)
                       
       NL(i,j,k,4) =   WENO_3_L( u2 (i,j,k)            ,t12,i,j,2,k) -   CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)
                           
       NL(i,j,k,5) =   WENO_3_L( u2 (i,j,k)            ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)               
      else 
       NL(i,j,k,3) =   BACKWIND_1( u2 (i,j,k)          ,t11,i,j,2,k) - 2*CENTRE_2( t12(i,j,k)            ,u1 ,i,j,2,k)
                       
       NL(i,j,k,4) =   BACKWIND_1( u2 (i,j,k)          ,t12,i,j,2,k) - CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                     - CENTRE_2( t12(i,j,k)            ,u2 ,i,j,2,k)
                           
       NL(i,j,k,5) =   BACKWIND_1( u2 (i,j,k)          ,t22,i,j,2,k) - 2*CENTRE_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)
      end if
    else if (j == N2) then
    NL(i,j,k,1) =   UPWIND_2( u2 (i,j,k)            ,u1 ,i,j,2,k) -   UPWIND_2( 1._nk                ,t12,i,j,2,k)
 
    NL(i,j,k,2) =   UPWIND_2( u2 (i,j,k)            ,u2 ,i,j,2,k) -   UPWIND_2( 1._nk                ,t22,i,j,2,k)

    NL(i,j,k,3) =   UPWIND_2( u2 (i,j,k)            ,t11,i,j,2,k) - 2*UPWIND_2( t12(i,j,k)            ,u1 ,i,j,2,k)
                      
    NL(i,j,k,4) =   UPWIND_2( u2 (i,j,k)            ,t12,i,j,2,k) -   UPWIND_2( 1._nk/Ma2 + t22(i,j,k),u1 ,i,j,2,k)   &
                  - UPWIND_2( t12(i,j,k)            ,u2 ,i,j,2,k)
                          
    NL(i,j,k,5) =   UPWIND_2( u2 (i,j,k)            ,t22,i,j,2,k) - 2*UPWIND_2( 1._nk/Ma2 + t22(i,j,k),u2 ,i,j,2,k)

    end if
   
             
#endif   
#endif
 end subroutine calcaul_termes_Convections_2


#if(DIMENSION_GEO ==3)    
 subroutine calcaul_termes_Convections_3(i, j, u1, u2, t11, t12, t22, NL, k, u3, t13, t23, t33)
   implicit none   
   integer, intent(in)         :: i, j, k
   real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in)   :: u1, u2, t11, t12, t22, u3, t13, t23, t33
   real(nk), dimension(1:N1,1:N2,1:N3,1:9), intent(out)    :: NL

   
   if (k ==1 ) then 
    NL(i,j,k,1) =   BACKWIND_2( u3 (i,j,k)            ,u1 ,i,j,3,k) - BACKWIND_2( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   BACKWIND_2( u3 (i,j,k)            ,u2 ,i,j,3,k) - BACKWIND_2( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   BACKWIND_2( u3 (i,j,k)            ,u3 ,i,j,3,k) - BACKWIND_2( 1._nk                     ,t33,i,j,3,k)     

    NL(i,j,k,4) =-2*BACKWIND_2( t13(i,j,k)            ,u1 ,i,j,3,k) + BACKWIND_2( u3 (i,j,k)                ,t11,i,j,3,k)

    NL(i,j,k,5) = - BACKWIND_2( t23(i,j,k)            ,u1 ,i,j,3,k) - BACKWIND_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                  + BACKWIND_2( t12(i,j,k)            ,u3 ,i,j,3,k) + BACKWIND_2( u3 (i,j,k)                ,t12,i,j,3,k)

    NL(i,j,k,6) = - BACKWIND_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - BACKWIND_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                  + BACKWIND_2( u3 (i,j,k)            ,t13,i,j,3,k)

    NL(i,j,k,7) =-2*BACKWIND_2( t23(i,j,k)            ,u2 ,i,j,3,k) + BACKWIND_2( u3 (i,j,k)                ,t22,i,j,3,k)

    NL(i,j,k,8) = - BACKWIND_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + BACKWIND_2( u3 (i,j,k)                ,t23,i,j,3,k)

    NL(i,j,k,9) =-2*BACKWIND_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + BACKWIND_2( u3 (i,j,k)                ,t33,i,j,3,k)  
    
    else if (k == 2 ) then
    NL(i,j,k,1) =   CENTRE_2( u3 (i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   CENTRE_2( u3 (i,j,k)            ,u2 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   CENTRE_2( u3 (i,j,k)            ,u3 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t33,i,j,3,k)     
      if (u1 (i,j,k) > 0._nk) then
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + UPWIND_1( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + UPWIND_1( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + UPWIND_1( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + UPWIND_1( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + UPWIND_1( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + UPWIND_1( u3 (i,j,k)                ,t33,i,j,3,k)  
      else 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_R( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t33,i,j,3,k)  
      end if

    else if ( k == 3 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u3 (i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   CENTRE_2( u3 (i,j,k)            ,u2 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   CENTRE_2( u3 (i,j,k)            ,u3 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t33,i,j,3,k)  
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_L( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t33,i,j,3,k)  

      else 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_R( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t33,i,j,3,k)  
      end if 
#elif (CENTRAL_D == 2)
    NL(i,j,k,1) =   CENTRE_4( u3 (i,j,k)            ,u1 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   CENTRE_4( u3 (i,j,k)            ,u2 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   CENTRE_4( u3 (i,j,k)            ,u3 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t33,i,j,3,k)  
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_4( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_4( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_4( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_L( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_4( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t33,i,j,3,k)  

      else 
        NL(i,j,k,4) =-2*CENTRE_4( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_4( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_4( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_5_R( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_4( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t33,i,j,3,k)  
      end if 

#endif

    else if ( k > 3 .and. k < N3 - 2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u3 (i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   CENTRE_2( u3 (i,j,k)            ,u2 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   CENTRE_2( u3 (i,j,k)            ,u3 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t33,i,j,3,k)  
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_L( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t33,i,j,3,k)  

      else 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_R( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t33,i,j,3,k)  
      end if 
#elif( CENTRAL_D == 2)
    NL(i,j,k,1) =   CENTRE_4( u3 (i,j,k)            ,u1 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   CENTRE_4( u3 (i,j,k)            ,u2 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   CENTRE_4( u3 (i,j,k)            ,u3 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t33,i,j,3,k)  
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_4( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_4( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_4( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_5_L( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_4( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t33,i,j,3,k)  

      else 
        NL(i,j,k,4) =-2*CENTRE_4( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_4( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_4( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_5_R( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_4( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_5_R( u3 (i,j,k)                ,t33,i,j,3,k)  
      end if 
#endif

    else if ( k == N3-2 ) then

#if (CENTRAL_D == 1)
    NL(i,j,k,1) =   CENTRE_2( u3 (i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   CENTRE_2( u3 (i,j,k)            ,u2 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   CENTRE_2( u3 (i,j,k)            ,u3 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t33,i,j,3,k)  
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_L( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t33,i,j,3,k)  
      else 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_R( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t33,i,j,3,k)  
      end if 

#elif( CENTRAL_D == 2 )
    NL(i,j,k,1) =   CENTRE_4( u3 (i,j,k)            ,u1 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   CENTRE_4( u3 (i,j,k)            ,u2 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   CENTRE_4( u3 (i,j,k)            ,u3 ,i,j,3,k) - CENTRE_4( 1._nk                     ,t33,i,j,3,k)  
      if (u1 (i,j,k) > 0.) then 
        NL(i,j,k,4) =-2*CENTRE_4( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_4( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_4( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_4( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_5_L( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_4( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_4( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_5_L( u3 (i,j,k)                ,t33,i,j,3,k) 

      else 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_R( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_R( u3 (i,j,k)                ,t33,i,j,3,k)  
      end if 
#endif

    else if ( k == N3-1 ) then
    NL(i,j,k,1) =   CENTRE_2( u1 (i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t11,i,j,3,k)

    NL(i,j,k,2) =   CENTRE_2( u1 (i,j,k)            ,u2 ,i,j,3,k) - CENTRE_2( 1._nk                     ,t12,i,j,3,k)
      if (u1 (i,j,k) > 0.) then
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                      + WENO_3_L( u3 (i,j,k)            ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + WENO_3_L( u3 (i,j,k)                 ,t33,i,j,3,k)  
      else 
        NL(i,j,k,4) =-2*CENTRE_2( t13(i,j,k)            ,u1 ,i,j,3,k) + BACKWIND_1( u3 (i,j,k)               ,t11,i,j,3,k)

        NL(i,j,k,5) = - CENTRE_2( t23(i,j,k)            ,u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                 ,u2 ,i,j,3,k) &
                      + CENTRE_2( t12(i,j,k)            ,u3 ,i,j,3,k) + BACKWIND_1( u3 (i,j,k)               ,t12,i,j,3,k)
    
        NL(i,j,k,6) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - CENTRE_2( t13(i,j,k)                 ,u3 ,i,j,3,k) &
                      + BACKWIND_1( u3 (i,j,k)          ,t13,i,j,3,k)
    
        NL(i,j,k,7) =-2*CENTRE_2( t23(i,j,k)            ,u2 ,i,j,3,k) + BACKWIND_1( u3 (i,j,k)               ,t22,i,j,3,k)
    
        NL(i,j,k,8) = - CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + BACKWIND_1( u3 (i,j,k)               ,t23,i,j,3,k)
    
        NL(i,j,k,9) =-2*CENTRE_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + BACKWIND_1( u3 (i,j,k)                 ,t33,i,j,3,k)  
      end if

    else if (k == N3) then

    NL(i,j,k,1) =   UPWIND_2( u3 (i,j,k)            ,u1 ,i,j,3,k) - UPWIND_2( 1._nk                     ,t13,i,j,3,k) 

    NL(i,j,k,2) =   UPWIND_2( u3 (i,j,k)            ,u2 ,i,j,3,k) - UPWIND_2( 1._nk                     ,t23,i,j,3,k)

    NL(i,j,k,3) =   UPWIND_2( u3 (i,j,k)            ,u3 ,i,j,3,k) - UPWIND_2( 1._nk                     ,t33,i,j,3,k)     

    NL(i,j,k,4) =-2*UPWIND_2( t13(i,j,k)            ,u1 ,i,j,3,k) + UPWIND_2( u3 (i,j,k)                ,t11,i,j,3,k)

    NL(i,j,k,5) = - UPWIND_2( t23(i,j,k)            ,u1 ,i,j,3,k) - UPWIND_2( t13(i,j,k)                ,u2 ,i,j,3,k) &
                  + UPWIND_2( t12(i,j,k)            ,u3 ,i,j,3,k) + UPWIND_2( u3 (i,j,k)                ,t12,i,j,3,k)

    NL(i,j,k,6) = - UPWIND_2( 1._nk/Ma2 + t33(i,j,k),u1 ,i,j,3,k) - UPWIND_2( t13(i,j,k)                ,u3 ,i,j,3,k) &
                  + UPWIND_2( u3 (i,j,k)            ,t13,i,j,3,k)

    NL(i,j,k,7) =-2*UPWIND_2( t23(i,j,k)            ,u2 ,i,j,3,k) + UPWIND_2( u3 (i,j,k)                ,t22,i,j,3,k)

    NL(i,j,k,8) = - UPWIND_2( 1._nk/Ma2 + t33(i,j,k),u2 ,i,j,3,k) + UPWIND_2( u3 (i,j,k)                ,t23,i,j,3,k)

    NL(i,j,k,9) =-2*UPWIND_2( 1._nk/Ma2 + t33(i,j,k),u3 ,i,j,3,k) + UPWIND_2( u3 (i,j,k)                ,t33,i,j,3,k)  
    end if 
             
   
 end subroutine calcaul_termes_Convections_3
#endif






      end module mHyperbolic_part
    
      
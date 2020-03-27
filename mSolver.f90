#include "definitions.h"

module mSolver

   use mBase
   use mPression
   use mVitesse

   ! Declarations
   implicit none

contains

  

#if (DIMENSION_GEO == 2)
  subroutine solver_ADI (vari)
    implicit none
    integer:: i, j, k
    type(variable),intent(inout)       :: vari
    real(nk),dimension(0:N1+1,0:N2+1)  :: R
    real(nk),dimension(0:N1+1,0:N2+1)  :: vari_cal1, vari_cal2

    R        = 0.
    vari_cal1= 0.; vari_cal2= 0.
    call coefs_UT_CL_ADI (vari) 
   
!    CALL PRINT_2D(1,1,N1,N2,N1,N2,vari%ax1 ,"     vari%ax1       ")
!    call PRINT_2D(1,1,N1,N2,N1,N2,vari%cx1 ,"     vari%cx1       ")
!    call PRINT_2D(1,1,N1,N2,N1,N2,vari%ax2 ,"     vari%ax2       ")
!    call PRINT_2D(1,1,N1,N2,N1,N2,vari%cx2 ,"     vari%cx2       ")
!    call PRINT_2D(1,1,N1,N2,N1,N2,vari%b ,  "     vari%b         ")
!   pause

       do j= 1, N2
          do i=1, N1
             R(i,j) =               vari%d(i,j)  - ( + & 
                  vari%ax1(i,j)*vari%var(i-1,j)      + &
                  vari%cx1(i,j)*vari%var(i+1,j)      + &
                  vari%ax2(i,j)*vari%var(i,j-1)      + &
                  vari%cx2(i,j)*vari%var(i,j+1)      + &
                  vari%b  (i,j)*vari%var(i,j)      )    
          end do
       end do
       
       do j=1,N2
          call thomas( vari%ax1(1:N1,j), vari%b1(1:N1,j), vari%cx1(1:N1,j), R(1:N1,j),         N1, vari_cal1(1:N1,j) )
       end do

       do i=1,N1
          call thomas( vari%ax2(i,1:N2), vari%b2(i,1:N2), vari%cx2(i,1:N2), vari_cal1(i,1:N2), N2, vari_cal2(i,1:N2) )
       end do
       vari%var(1:N1,1:N2) = vari%var(1:N1,1:N2) + vari_cal2(1:N1,1:N2)

    
  end subroutine solver_ADI












  subroutine  coefs_UT_CL_ADI (var) 
   implicit none
   type(variable),intent(inout):: var
   integer :: i, j, k
!-------------- bas------------
   do i= 2, N1-1 
       j = 1
      if ( var%CL_var_bas ==1 ) call CL_Dirichlet(var%var_bas , var%d(i,j) , var%b(i,j)  ,&
                                                  var%ax1(i,j), var%b1(i,j), var%cx1(i,j),&
                                                  var%ax2(i,j), var%b2(i,j), var%cx2(i,j),&
                                                  var= var%var(i,j))
      if ( var%CL_var_bas ==2 ) call CL_Neumann (var%d(i,j), var%cx2(i,j), var%ax2(i,j), dx2_R(i,j), var%q_var_bas)        

!-------------- haut------------
       j = N2 
       if ( var%CL_var_haut ==1 ) call CL_Dirichlet(var%var_haut, var%d(i,j) , var%b(i,j)  ,&
                                                    var%ax1(i,j), var%b1(i,j), var%cx1(i,j),&
                                                    var%ax2(i,j), var%b2(i,j), var%cx2(i,j),&
                                                    var= var%var(i,j))
       if ( var%CL_var_haut ==2 ) call CL_Neumann  (var%d(i,j), var%ax2(i,j), var%cx2(i,j), dx2_L(i,j), var%q_var_haut)
   end do

!-------------- gauche------------
   do j = 2, N2-1
       i = 1
       if ( var%CL_var_gauche ==1 ) call CL_Dirichlet(var%var_gauche,var%d(i,j), var%b(i,j)  ,&
                                                      var%ax1(i,j), var%b1(i,j), var%cx1(i,j),&
                                                      var%ax2(i,j), var%b2(i,j), var%cx2(i,j),&
                                                      var= var%var(i,j))
       if ( var%CL_var_gauche ==2 ) call CL_Neumann  (var%d(i,j), var%cx1(i,j), var%ax1(i,j), dx1_R(i,j), var%q_var_gauche)

!-------------- droite------------
       i = N1
      if ( var%CL_var_droite ==1 ) call CL_Dirichlet(var%var_droite,var%d(i,j), var%b(i,j)  ,&
                                                     var%ax1(i,j), var%b1(i,j), var%cx1(i,j),&
                                                     var%ax2(i,j), var%b2(i,j), var%cx2(i,j),&
                                                     var= var%var(i,j))
      if ( var%CL_var_droite ==2 ) call CL_Neumann  (var%d(i,j), var%ax1(i,j), var%cx1(i,j), dx1_L(i,j), var%q_var_droite)
   end do

!!!--- points sur les coins
!!!--- droite-bas ---
    i = N1; j = 1
   call CL_coin (                                &
        var%CL_var_bas   , var%CL_var_droite   , &
        var%var_bas      , var%var_droite      , &
        var%q_var_bas    , var%q_var_droite    , &
        dx2_R(i,j)       , dx1_L(i,j)          , &
        var%ax2(i,j)     , var%cx1(i,j)        , &
        var%cx2(i,j)     , var%ax1(i,j)        , &
        var%d(i,j)       , var%b(i,j)          , &
        var%ax1(i,j), var%b1(i,j), var%cx1(i,j), &
        var%ax2(i,j), var%b2(i,j), var%cx2(i,j), &
        varia= var%var(i,j) )

!!!--- droite-haut ---
    i = N1; j = N2
   call CL_coin (                                &
        var%CL_var_droite, var%CL_var_haut     , &
        var%var_droite   , var%var_haut        , &
        var%q_var_droite , var%q_var_haut      , &
        dx1_L(i,j)       , dx2_L(i,j)          , &
        var%cx1(i,j)     , var%cx2(i,j)        , &
        var%ax1(i,j)     , var%ax2(i,j)        , &
        var%d(i,j)       , var%b(i,j)          , &
        var%ax1(i,j), var%b1(i,j), var%cx1(i,j), &
        var%ax2(i,j), var%b2(i,j), var%cx2(i,j), &
        varia= var%var(i,j) )  
   
!!!--- gauche-haut ---
    i = 1; j = N2
   call CL_coin (                                &
        var%CL_var_haut  , var%CL_var_gauche   , &
        var%var_haut     , var%var_gauche      , &
        var%q_var_haut   , var%q_var_gauche    , &
        dx2_L(i,j)       , dx1_R(i,j)          , &
        var%cx2(i,j)     , var%ax1(i,j)        , &
        var%ax2(i,j)     , var%cx1(i,j)        , &
        var%d(i,j)       , var%b(i,j)          , &
        var%ax1(i,j), var%b1(i,j), var%cx1(i,j), &
        var%ax2(i,j), var%b2(i,j), var%cx2(i,j), &
        varia= var%var(i,j) )

   
!!!--- gauche-bas ---
    i = 1; j = 1
   call CL_coin (                                &
        var%CL_var_gauche, var%CL_var_bas      , &
        var%var_gauche   , var%var_bas         , &
        var%q_var_gauche , var%q_var_bas       , &
        dx1_R(i,j)       , dx2_R(i,j)          , &
        var%ax1(i,j)     , var%ax2(i,j)        , &
        var%cx1(i,j)     , var%cx2(i,j)        , &
        var%d(i,j)       , var%b(i,j)          , &
        var%ax1(i,j), var%b1(i,j), var%cx1(i,j), &
        var%ax2(i,j), var%b2(i,j), var%cx2(i,j), &
        varia= var%var(i,j) )


 end subroutine coefs_UT_CL_ADI

#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if (DIMENSION_GEO == 3)
 subroutine solver_ADI_3D  (vari)
   implicit none
   integer:: i, j, k
   type(variable),intent(inout)              :: vari
   real(nk),dimension(0:N1+1,0:N2+1,0:N3+1)  :: R
   real(nk),dimension(0:N1+1,0:N2+1,0:N3+1)  :: vari_cal1
   real(nk),dimension(0:N1+1,0:N2+1,0:N3+1)  :: vari_cal2   
   real(nk),dimension(0:N1+1,0:N2+1,0:N3+1)  :: vari_cal3   
  
   R        = 0.
   vari_cal1= 0.; vari_cal2= 0.; vari_cal3= 0.

   call coefs_UT_CL_ADI_3D(vari)
   
   do k= 1, N3
      do j= 1, N2
         do i=1, N1               
             R(i,j,k) =              vari%d(i,j,k)  - (   + & 
                 vari%ax1(i,j,k)*vari%var(i-1,j,k)        + &
                 vari%cx1(i,j,k)*vari%var(i+1,j,k)        + &
                 vari%ax2(i,j,k)*vari%var(i,j-1,k)        + &
                 vari%cx2(i,j,k)*vari%var(i,j+1,k)        + &
                 vari%ax3(i,j,k)*vari%var(i,j,k-1)        + &
                 vari%cx3(i,j,k)*vari%var(i,j,k+1)        + &
                 vari%b  (i,j,k)*vari%var(i,j,k)      )
        end do
      end do
   end do

   do k=1,N3
      do j=1,N2
         call thomas( vari%ax1(1:N1,j,k), vari%b1(1:N1,j,k), vari%cx1(1:N1,j,k),         R(1:N1,j,k), N1, vari_cal1(1:N1,j,k))
      end do
   end do

   do k=1,N3
      do i=1,N1
         call thomas( vari%ax2(i,1:N2,k), vari%b2(i,1:N2,k), vari%cx2(i,1:N2,k), vari_cal1(i,1:N2,k), N2, vari_cal2(i,1:N2,k) )
      end do
   end do

   do j=1,N2
      do i=1,N1
         call thomas( vari%ax3(i,j,1:N3), vari%b3(i,j,1:N3), vari%cx3(i,j,1:N3), vari_cal2(i,j,1:N3), N3, vari_cal3(i,j,1:N3) )
      end do
   end do
   vari%var(1:N1,1:N2,1:N3) = vari%var(1:N1,1:N2,1:N3) + vari_cal3(1:N1,1:N2,1:N3)
 end subroutine solver_ADI_3D





  subroutine  coefs_UT_CL_ADI_3D (var)
    implicit none
    type(variable),intent(inout):: var
    integer :: i, j, k
    
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   droite
    do j = 2, N2-1
      do i= 2, N1-1
          k = 1          
         if ( var%CL_var_droite ==1 )                                          &
              call CL_Dirichlet(var%var_droite, var%d(i,j,k),  var%b(i,j,k),   &
                                var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
                                var= var%var(i,j,k)  )   
         if ( var%CL_var_droite ==2 ) call CL_Neumann (var%d(i,j,k), var%cx3(i,j,k), var%ax3(i,j,k), dx3_R(i,j,k), var%q_var_droite)   
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   gauche
          k = N3                
         if ( var%CL_var_gauche ==1 )                                       & 
              call CL_Dirichlet(var%var_gauche, var%d(i,j,k), var%b(i,j,k), &
                                var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
                                var= var%var(i,j,k)  )    
         if ( var%CL_var_gauche ==2 ) call CL_Neumann (var%d(i,j,k), var%ax3(i,j,k), var%cx3(i,j,k), dx3_L(i,j,k), var%q_var_gauche) 
      end do 
    end do   
     

 !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   bas
    do k = 2, N3-1
      do i = 2, N1-1
          j = 1            
         if ( var%CL_var_bas ==1 )                                             & 
              call CL_Dirichlet(var%var_bas   , var%d(i,j,k) , var%b(i,j,k)  , &
                                var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
                                var= var%var(i,j,k)  )     
         if ( var%CL_var_bas ==2 ) call CL_Neumann (var%d(i,j,k), var%cx2(i,j,k), var%ax2(i,j,k), dx2_R(i,j,k), var%q_var_bas)     
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   haut
          j = N2          
         if ( var%CL_var_haut ==1 )                                            & 
              call CL_Dirichlet(var%var_haut  , var%d(i,j,k) , var%b(i,j,k)  , &
                                var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
                                var= var%var(i,j,k)  )      
         if ( var%CL_var_haut ==2 ) call CL_Neumann (var%d(i,j,k), var%ax2(i,j,k), var%cx2(i,j,k), dx2_L(i,j,k), var%q_var_haut)
      end do 
    end do
 
 
 
     
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   avant
    do k = 2, N3-1
      do j = 2, N2-1
          i = 1
            if ( var%CL_var_avant ==1 )                                         & 
               call CL_Dirichlet(var%var_avant , var%d(i,j,k) , var%b(i,j,k)  , &
                                 var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                 var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                 var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
                                 var= var%var(i,j,k)  )   
            if ( var%CL_var_avant ==2 ) call CL_Neumann (var%d(i,j,k), var%cx1(i,j,k), var%ax1(i,j,k), dx1_R(i,j,k),&
                                                         var%q_var_avant)              
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   arriere
          i = N1
          if ( var%CL_var_arriere ==1 )                                         & 
               call CL_Dirichlet(var%var_arriere,var%d(i,j,k) , var%b(i,j,k)  , &
                                 var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                 var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                 var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
                                 var= var%var(i,j,k)  )   
         if ( var%CL_var_arriere ==2 ) call CL_Neumann (var%d(i,j,k), var%ax1(i,j,k), var%cx1(i,j,k), dx1_L(i,j,k),&
                                                        var%q_var_arriere)          
      end do 
    end do
 
 
  
 
 
 !!!---- points sur les segments croisés (intersection entre deux faces)
     !droite-bas
     do i = 1, N1
         j = 1
          k = 1
        call CL_coin (                                      &
             var%CL_var_droite   , var%CL_var_bas         , &
             var%var_droite      , var%var_bas            , &
             var%q_var_droite    , var%q_var_bas          , &
             dx3_R(i,j,k)        , dx2_R(i,j,k)           , &
             var%ax3(i,j,k)      , var%ax2(i,j,k)         , &
             var%cx3(i,j,k)      , var%cx2(i,j,k)         , &
             var%d(i,j,k)        , var%b(i,j,k)           , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
   
      
     !droite-haut
     do i = 1, N1
         j = N2
          k = 1
        call CL_coin (                                      &
             var%CL_var_droite   , var%CL_var_haut        , &
             var%var_droite      , var%var_haut           , &
             var%q_var_droite    , var%q_var_haut         , &
             dx3_R(i,j,k)        , dx2_L(i,j,k)           , &
             var%ax3(i,j,k)      , var%cx2(i,j,k)         , &
             var%cx3(i,j,k)      , var%ax2(i,j,k)         , &
             var%d(i,j,k)        , var%b(i,j,k)           , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
     !gauche-bas
     do i = 1, N1
         j = 1
          k = N3
        call CL_coin (                                      &
             var%CL_var_gauche   , var%CL_var_bas         , &
             var%var_gauche      , var%var_bas            , &
             var%q_var_gauche    , var%q_var_bas          , &
             dx3_L(i,j,k)        , dx2_R(i,j,k)           , &
             var%cx3(i,j,k)      , var%ax2(i,j,k)         , &
             var%ax3(i,j,k)      , var%cx2(i,j,k)         , &
             var%d(i,j,k)        , var%b(i,j,k)           , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
 
 
     !gauche-haut
     do i = 1, N1
          j = N2
           k = N3          
         call CL_coin (                                     &
             var%CL_var_gauche   , var%CL_var_haut        , &
             var%var_gauche      , var%var_haut           , &
             var%q_var_gauche    , var%q_var_haut         , &
             dx3_L(i,j,k)        , dx2_L(i,j,k)           , &
             var%cx3(i,j,k)      , var%cx2(i,j,k)         , &
             var%ax3(i,j,k)      , var%ax2(i,j,k)         , &
             var%d(i,j,k)        , var%b(i,j,k)           , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
 
     !droite-arriere
     do j = 1, N2
         i = N1
          k = 1 
        call CL_coin (                                      &
             var%CL_var_droite   , var%CL_var_arriere     , &
             var%var_droite      , var%var_arriere        , &
             var%q_var_droite    , var%q_var_arriere      , &
             dx3_R(i,j,k)        , dx1_L(i,j,k)           , &
             var%ax3(i,j,k)      , var%cx1(i,j,k)         , &
             var%cx3(i,j,k)      , var%ax1(i,j,k)         , &
             var%d(i,j,k)        , var%b(i,j,k)           , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                     )
     end do 
      
 
     !droite-avant
     do j = 1, N2
         i = 1
          k = 1        
        call CL_coin (                                      &
             var%CL_var_droite   , var%CL_var_avant       , &
             var%var_droite      , var%var_avant          , &
             var%q_var_droite    , var%q_var_avant        , &
             dx3_R(i,j,k)        , dx1_R(i,j,k)           , &
             var%ax3(i,j,k)      , var%ax1(i,j,k)         , &
             var%cx3(i,j,k)      , var%cx1(i,j,k)         , &
             var%d(i,j,k)        , var%b(i,j,k)           , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do   
     
     !gauche-arriere
     do j = 1, N2
         i = N1
          k = N3
        call CL_coin (                                      &
             var%CL_var_gauche   , var%CL_var_arriere     , &
             var%var_gauche      , var%var_arriere        , &
             var%q_var_gauche    , var%q_var_arriere      , &
             dx3_L(i,j,k)        , dx1_L(i,j,k)           , &
             var%cx3(i,j,k)      , var%cx1(i,j,k)         , &
             var%ax3(i,j,k)      , var%ax1(i,j,k)         , &
             var%d(i,j,k)        , var%b(i,j,k)           , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                     )
     end do
     
     !gauche-avant
     do j = 1, N2
         i = 1
          k = N3
        call CL_coin (                                      &
             var%CL_var_gauche   , var%CL_var_avant       , &
             var%var_gauche      , var%var_avant          , &
             var%q_var_gauche    , var%q_var_avant        , &
             dx3_L(i,j,k)        , dx1_R(i,j,k)           , &
             var%cx3(i,j,k)      , var%ax1(i,j,k)         , &
             var%ax3(i,j,k)      , var%cx1(i,j,k)         , &
             var%d(i,j,k)         , var%b(i,j,k)          , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
     
     
 
     !avant-haut
     do k = 1, N3
         i = 1
          j = N2     
        call CL_coin (                                      &
             var%CL_var_avant    , var%CL_var_haut        , &
             var%var_avant       , var%var_haut           , &
             var%q_var_avant     , var%q_var_haut         , &
             dx1_R(i,j,k)        , dx2_L(i,j,k)           , &
             var%ax1(i,j,k)      , var%cx2(i,j,k)         , &
             var%cx1(i,j,k)      , var%ax2(i,j,k)         , &
             var%d(i,j,k)         , var%b(i,j,k)          , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
     
     !avant-bas
     do k = 1, N3
         i = 1
          j = 1
        call CL_coin (                                      &
             var%CL_var_avant    , var%CL_var_bas         , &
             var%var_avant       , var%var_bas            , &
             var%q_var_avant     , var%q_var_bas          , &
             dx1_R(i,j,k)        , dx2_R(i,j,k)           , &
             var%ax1(i,j,k)      , var%ax2(i,j,k)         , &
             var%cx1(i,j,k)      , var%cx2(i,j,k)         , &
             var%d(i,j,k)         , var%b(i,j,k)          , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
     
     !arriere-haut
     do k = 1, N3
         i = N1
          j = N2
        call CL_coin (                                      &
             var%CL_var_arriere  , var%CL_var_haut        , &
             var%var_arriere     , var%var_haut           , &
             var%q_var_arriere   , var%q_var_haut         , &
             dx1_L(i,j,k)        , dx2_L(i,j,k)           , &
             var%cx1(i,j,k)      , var%cx2(i,j,k)         , &
             var%ax1(i,j,k)      , var%ax2(i,j,k)         , &
             var%d(i,j,k)         , var%b(i,j,k)          , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
     
     !arriere-bas
     do k = 1, N3
         i = N1
          j = 1 
        call CL_coin (                                      &
             var%CL_var_arriere  , var%CL_var_bas         , &
             var%var_arriere     , var%var_bas            , &
             var%q_var_arriere   , var%q_var_bas          , &
             dx1_L(i,j,k)        , dx2_R(i,j,k)           , &
             var%cx1(i,j,k)      , var%ax2(i,j,k)         , &
             var%ax1(i,j,k)      , var%cx2(i,j,k)         , &
             var%d  (i,j,k)      , var%b  (i,j,k)         , &
             var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
             var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
             var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k), &
             varia= var%var(i,j,k)                        )
     end do
 
 
  endsubroutine coefs_UT_CL_ADI_3D

#endif

  
end module mSolver


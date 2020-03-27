#include "definitions.h"

module mVitesse

  use mBase
  use mOutils
   
  ! Declarations
  implicit none






contains
!!!!                                      (T , 2. / 3. * dt* T%nomb_ad, 1.)
  subroutine coefs_UT_diffusion_implicite ( var, const, const_2 )
    implicit none
    type(variable),intent(inout):: var
    real(nk),intent(in)         :: const , const_2
    integer                     :: i, j, k
#if (DIMENSION_GEO == 2)     

    do j= 1, N2
       do i=1, N1 
          var%ax1(i,j) =         - const * 2.0_nk / ( dx1_L(i,j)*( dx1_L(i,j)+dx1_R(i,j) ) )
          var%b1 (i,j) = const_2 + const * 2.0_nk / ( dx1_L(i,j)*  dx1_R(i,j)              )
          var%cx1(i,j) =         - const * 2.0_nk / ( dx1_R(i,j)*( dx1_L(i,j)+dx1_R(i,j) ) )
          var%ax2(i,j) =         - const * 2.0_nk / ( dx2_L(i,j)*( dx2_L(i,j)+dx2_R(i,j) ) )
          var%b2 (i,j) = const_2 + const * 2.0_nk / ( dx2_L(i,j)*  dx2_R(i,j)              )
          var%cx2(i,j) =         - const * 2.0_nk / ( dx2_R(i,j)*( dx2_L(i,j)+dx2_R(i,j) ) )
          var%b  (i,j) = const_2 + const *(2.0_nk / ( dx1_L(i,j)*  dx1_R(i,j)              ) &
                                         + 2.0_nk / ( dx2_L(i,j)*  dx2_R(i,j)              ))                 
       end do
    end do

#elif (DIMENSION_GEO == 3)

       do k = 1, N3 
          do j= 1, N2          
             do i=1, N1             
                var%ax1(i,j,k) =         - const * 2.0_nk / ( dx1_L(i,j,k)*( dx1_L(i,j,k)+dx1_R(i,j,k) ) )
                var%b1 (i,j,k) = const_2 + const * 2.0_nk / ( dx1_L(i,j,k)*  dx1_R(i,j,k)              )
                var%cx1(i,j,k) =         - const * 2.0_nk / ( dx1_R(i,j,k)*( dx1_L(i,j,k)+dx1_R(i,j,k) ) )
                var%ax2(i,j,k) =         - const * 2.0_nk / ( dx2_L(i,j,k)*( dx2_L(i,j,k)+dx2_R(i,j,k) ) )
                var%b2 (i,j,k) = const_2 + const * 2.0_nk / ( dx2_L(i,j,k)*  dx2_R(i,j,k)              )
                var%cx2(i,j,k) =         - const * 2.0_nk / ( dx2_R(i,j,k)*( dx2_L(i,j,k)+dx2_R(i,j,k) ) )
                var%ax3(i,j,k) =         - const * 2.0_nk / ( dx3_L(i,j,k)*( dx3_L(i,j,k)+dx3_R(i,j,k) ) )
                var%b3 (i,j,k) = const_2 + const * 2.0_nk / ( dx3_L(i,j,k)*  dx3_R(i,j,k)              )
                var%cx3(i,j,k) =         - const * 2.0_nk / ( dx3_R(i,j,k)*( dx3_L(i,j,k)+dx3_R(i,j,k) ) )
                var%b  (i,j,k) = const_2 + const *(2.0_nk / ( dx1_L(i,j,k)*  dx1_R(i,j,k)              ) &
                                                 + 2.0_nk / ( dx2_L(i,j,k)*  dx2_R(i,j,k)              ) &  
                                                 + 2.0_nk / ( dx3_L(i,j,k)*  dx3_R(i,j,k)              )) 
             end do
          end do
       end do
  
#endif 
  end subroutine coefs_UT_diffusion_implicite






  subroutine coefs_UT_convection(d, var_in, u1, u2, const, u3)
    implicit none
    integer :: i, j, k
    real(nk),intent(in):: const 
#if (DIMENSION_GEO == 2)
    real(nk),dimension(0:N1+1,0:N2+1),intent(inout):: d
    real(nk),dimension(0:N1+1,0:N2+1),intent(in)   :: var_in, u1, u2
    real(nk), optional, intent(in):: u3
#elif (DIMENSION_GEO == 3)
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(inout):: d
    real(nk),dimension(0:N1+1,0:N2+1,0:N3+1),intent(in)   :: var_in, u1, u2, u3
#endif

#if (DIMENSION_GEO == 2)
    do j= 1, N2
       do i=1, N1
          d(i,j) =  d(i,j) - const * ( &
        + CENTRE_2(u1(i,j), var_in, i, j, 1) &
        + CENTRE_2(u2(i,j), var_in, i, j, 2)  )
       end do
    end do
#elif (DIMENSION_GEO == 3)
       do k= 1, N3
          do j= 1, N2
             do i=1, N1
                d(i,j,k) =  d(i,j,k) - const * ( &
                  CENTRE_2(u1(i,j,k), var_in, i, j, 1, k) &
                + CENTRE_2(u2(i,j,k), var_in, i, j, 2, k) &
                + CENTRE_2(u3(i,j,k), var_in, i, j, 3, k) )                                                         
             end do
          end do
       end do
#endif    
  end subroutine coefs_UT_convection













  subroutine  coefs_UT_CL (var, Phi)
    implicit none
    type(variable),intent(inout)                :: var
    real(nk),dimension(0:N1+1,0:N2+1),intent(in):: Phi
    integer                                     :: i, j, k
 

#if (DIMENSION_GEO == 2)
!!!--------bas--------
    do i= 2, N1-1
       j=1      
       if ( var%CL_var_bas ==1 ) call CL_Dirichlet(var%var_bas , var%d(i,j) , var%b(i,j)  ,&
                                                   var%ax1(i,j), var%b1(i,j), var%cx1(i,j),&
                                                   var%ax2(i,j), var%b2(i,j), var%cx2(i,j))
       if ( var%CL_var_bas ==2 ) call CL_Neumann  (var%d(i,j), var%cx2(i,j), var%ax2(i,j), dx2_R(i,j), var%q_var_bas)             
!!!--------haut--------
       j=N2
      if ( var%CL_var_haut ==1 ) call CL_Dirichlet(var%var_haut, var%d(i,j) , var%b(i,j)  ,&
                                                   var%ax1(i,j), var%b1(i,j), var%cx1(i,j),&
                                                   var%ax2(i,j), var%b2(i,j), var%cx2(i,j))
      if ( var%CL_var_haut ==2 ) call CL_Neumann  (var%d(i,j), var%ax2(i,j), var%cx2(i,j), dx2_L(i,j), var%q_var_haut) 
    end do
!!!--------gauche--------
    do j= 2, N2-1
       i= 1
       if ( var%CL_var_gauche ==1 ) call CL_Dirichlet(var%var_gauche,var%d(i,j), var%b(i,j)  ,&
                                                      var%ax1(i,j), var%b1(i,j), var%cx1(i,j),&
                                                      var%ax2(i,j), var%b2(i,j), var%cx2(i,j))
       if ( var%CL_var_gauche ==2 ) call CL_Neumann  (var%d(i,j), var%cx1(i,j), var%ax1(i,j), dx1_R(i,j), var%q_var_gauche)      
!!!--------droite--------
      i= N1
       if ( var%CL_var_droite ==1 ) call CL_Dirichlet(var%var_droite,var%d(i,j), var%b(i,j)  ,&
                                                      var%ax1(i,j), var%b1(i,j), var%cx1(i,j),&
                                                      var%ax2(i,j), var%b2(i,j), var%cx2(i,j))
       if ( var%CL_var_droite ==2 ) call CL_Neumann  (var%d(i,j), var%ax1(i,j), var%cx1(i,j), dx1_L(i,j), var%q_var_droite)                          
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!--- points sur les coins ---
!!!--- droite-bas ---
      i= N1; j= 1
    call CL_coin (                                       &
         var%CL_var_bas   , var%CL_var_droite          , &
         var%var_bas      , var%var_droite             , &
         var%q_var_bas    , var%q_var_droite           , &
         dx2_R(i,j)       , dx1_L(i,j)                 , &
         var%ax2(i,j)     , var%cx1(i,j)               , &
         var%cx2(i,j)     , var%ax1(i,j)               , &
         var%d(i,j)       , var%b(i,j)                 , &
         var%ax1(i,j)     , var%b1(i,j) , var%cx1(i,j) , &
         var%ax2(i,j)     , var%b2(i,j) , var%cx2(i,j) )
    
!!!--- droite-haut ---
      i= N1; j= N2
    call CL_coin (                                       &
         var%CL_var_droite, var%CL_var_haut            , &
         var%var_droite   , var%var_haut               , &
         var%q_var_droite , var%q_var_haut             , &
         dx1_L(i,j)       , dx2_L(i,j)                 , &
         var%cx1(i,j)     , var%cx2(i,j)               , &
         var%ax1(i,j)     , var%ax2(i,j)               , &
         var%d(i,j)       , var%b(i,j)                 , &
         var%ax1(i,j)     , var%b1(i,j) , var%cx1(i,j) , &
         var%ax2(i,j)     , var%b2(i,j) , var%cx2(i,j) ) 
!!!--- gauche-haut ---
      i= 1; j= N2
    call CL_coin (                                       &
         var%CL_var_haut  , var%CL_var_gauche          , &
         var%var_haut     , var%var_gauche             , &
         var%q_var_haut   , var%q_var_gauche           , &
         dx2_L(i,j)       , dx1_R(i,j)                 , &
         var%cx2(i,j)     , var%ax1(i,j)               , &
         var%ax2(i,j)     , var%cx1(i,j)               , & 
         var%d(i,j)       , var%b(i,j)                 , &
         var%ax1(i,j)     , var%b1(i,j) , var%cx1(i,j) , &
         var%ax2(i,j)     , var%b2(i,j) , var%cx2(i,j) )
!!!--- droite-bas ---
      i= 1; j= 1
    call CL_coin (                                       &
         var%CL_var_gauche, var%CL_var_bas             , &
         var%var_gauche   , var%var_bas                , &
         var%q_var_gauche , var%q_var_bas              , &
         dx1_R(i,j)       , dx2_R(i,j)                 , &
         var%ax1(i,j)     , var%ax2(i,j)               , &
         var%cx1(i,j)     , var%cx2(i,j)               , &
         var%d(i,j)       , var%b(i,j)                 , &
         var%ax1(i,j)     , var%b1(i,j) , var%cx1(i,j) , &
         var%ax2(i,j)     , var%b2(i,j) , var%cx2(i,j) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
#elif (DIMENSION_GEO == 3)
 
         
    do j = 2,N2-1
       do i=2,N1-1
           k=1
          if ( var%CL_var_droite ==1 )                                          & 
               call CL_Dirichlet(var%var_droite, var%d(i,j,k) , var%b(i,j,k)  , &
                                 var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                 var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                 var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  )  
          if ( var%CL_var_droite ==2 ) call CL_Neumann (var%d(i,j,k), var%cx3(i,j,k), &
                                                        var%ax3(i,j,k), dx3_R(i,j,k), var%q_var_droite)   
                    
           k=N3
          if ( var%CL_var_gauche ==1 )                                          & 
               call CL_Dirichlet(var%var_gauche, var%d(i,j,k) , var%b(i,j,k)  , &
                                 var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                 var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                 var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  )           
          if ( var%CL_var_gauche ==2 ) call CL_Neumann (var%d(i,j,k), var%ax3(i,j,k), &
                                                        var%cx3(i,j,k), dx3_L(i,j,k), var%q_var_gauche) 
                    
       end do 
    end do   
    
    
    do k=2,N3-1
       do i=2,N1-1
           j=1
          if ( var%CL_var_bas ==1 )                                             & 
               call CL_Dirichlet(var%var_bas, var%d(i,j,k) , var%b(i,j,k)     , &
                                 var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                 var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                 var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  )               
          if ( var%CL_var_bas ==2 ) call CL_Neumann (var%d(i,j,k), var%cx2(i,j,k), &
                                                     var%ax2(i,j,k), dx2_R(i,j,k), var%q_var_bas)     
           
           j=N2     
          if ( var%CL_var_haut ==1 )                                            & 
               call CL_Dirichlet(var%var_haut, var%d(i,j,k) , var%b(i,j,k)    , &
                                 var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                 var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                 var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  )           
          if ( var%CL_var_haut ==2 ) call CL_Neumann (var%d(i,j,k), var%ax2(i,j,k), &
                                                      var%cx2(i,j,k), dx2_L(i,j,k), var%q_var_haut)
       end do 
    end do



    
    
    do k=2,N3-1
       do j=2,N2-1
           i=1 
          if (var%indice /= 1) then                          
             if ( var%CL_var_avant ==1 )                                           & 
                  call CL_Dirichlet(var%var_avant , var%d(i,j,k) , var%b(i,j,k)  , &
                                    var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                    var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                    var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  )        
             if ( var%CL_var_avant ==2 ) call CL_Neumann (var%d(i,j,k), var%cx1(i,j,k), &
                                                          var%ax1(i,j,k), dx1_R(i,j,k), var%q_var_avant) 
             
          else if (var%indice == 1) then  
             
             if ( var%CL_var_avant ==1 )                                          & 
                  call CL_Dirichlet(var%var_avant, var%d(i,j,k) , var%b(i,j,k)  , &
                                   var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                   var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                   var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  )        
             if ( var%CL_var_avant ==2 ) call CL_Neumann (var%d(i,j,k), var%cx1(i,j,k), &
                                                          var%ax1(i,j,k), dx1_R(i,j,k), var%q_var_avant)              

          endif

           i=N1      
          if ( var%CL_var_arriere ==1 )                                         & 
               call CL_Dirichlet(var%var_arriere,var%d(i,j,k) , var%b(i,j,k) , &
                                 var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
                                 var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
                                 var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  )            
          if ( var%CL_var_arriere ==2 ) call CL_Neumann (var%d(i,j,k), var%ax1(i,j,k), &
                                                         var%cx1(i,j,k), dx1_L(i,j,k), var%q_var_arriere)          
          
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
         var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  )  

  
    !droite-haut
      j = N2
       k = 1
     call CL_coin (                                      &
          var%CL_var_droite   , var%CL_var_haut        , &
          var%var_droite      , var%var_haut           , &
          var%q_var_droite    , var%q_var_haut         , &
          dx3_R(i,j,k)        , dx2_L(i,j,k)            , &
          var%ax3(i,j,k)      , var%cx2(i,j,k)         , &
          var%cx3(i,j,k)      , var%ax2(i,j,k)         , &
          var%d(i,j,k)        , var%b(i,j,k)           , &
          var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
          var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
          var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 

    !gauche-bas

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
          var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 
     
    !gauche-haut

      j = N2
       k = N3    
     call CL_coin (                                      &
          var%CL_var_gauche   , var%CL_var_haut        , &
          var%var_gauche      , var%var_haut           , &
          var%q_var_gauche    , var%q_var_haut         , &
          dx3_L(i,j,k)        , dx2_L(i,j,k)           , &
          var%cx3(i,j,k)      , var%cx2(i,j,k)         , &
          var%ax3(i,j,k)      , var%ax2(i,j,k)         , &
          var%d(i,j,k)        , var%b(i,j,k)           , &
          var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
          var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
          var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 
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
            var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 

    !droite-avant
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
            var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 

    !gauche-arriere
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
            var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 

    !gauche-avant
        i = 1
         k = N3
       call CL_coin (                                      &
            var%CL_var_gauche   , var%CL_var_avant       , &
            var%var_gauche      , var%var_avant          , &
            var%q_var_gauche    , var%q_var_avant        , &
            dx3_L(i,j,k)        , dx1_R(i,j,k)           , &
            var%cx3(i,j,k)      , var%ax1(i,j,k)         , &
            var%ax3(i,j,k)      , var%cx1(i,j,k)         , &
            var%d(i,j,k)        , var%b(i,j,k)           , &
            var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
            var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
            var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 
    end do
    
    !avant-haut
    do k = 1, N3
        i = 1
         j = N2
       call CL_coin (                                      &
            var%CL_var_avant   , var%CL_var_haut         , &
            var%var_avant      , var%var_haut            , &
            var%q_var_avant    , var%q_var_haut          , &
            dx1_R(i,j,k)       , dx2_L(i,j,k)            , &
            var%ax1(i,j,k)     , var%cx2(i,j,k)          , &
            var%cx1(i,j,k)     , var%ax2(i,j,k)          , &
            var%d(i,j,k)       , var%b(i,j,k)            , &
            var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
            var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
            var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 
    
    !avant-bas
        i = 1
         j = 1
       call CL_coin (                                      &
            var%CL_var_avant   , var%CL_var_bas          , &
            var%var_avant      , var%var_bas             , &
            var%q_var_avant    , var%q_var_bas           , &
            dx1_R(i,j,k)       , dx2_R(i,j,k)            , &
            var%ax1(i,j,k)     , var%ax2(i,j,k)          , &
            var%cx1(i,j,k)     , var%cx2(i,j,k)          , &
            var%d(i,j,k)       , var%b(i,j,k)            , &
            var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
            var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
            var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 
    
    !arriere-haut
        i = N1
         j = N2
       call CL_coin (                                      &
            var%CL_var_arriere   , var%CL_var_haut       , &
            var%var_arriere      , var%var_haut          , &
            var%q_var_arriere    , var%q_var_haut        , &
            dx1_L(i,j,k)         , dx2_L(i,j,k)          , &
            var%cx1(i,j,k)       , var%cx2(i,j,k)        , &
            var%ax1(i,j,k)       , var%ax2(i,j,k)        , &
            var%d(i,j,k)         , var%b(i,j,k)          , &
            var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
            var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
            var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 
    
    !arriere-bas
        i = N1
         j = 1
       call CL_coin (                                      &
            var%CL_var_arriere   , var%CL_var_bas        , &
            var%var_arriere      , var%var_bas           , &
            var%q_var_arriere    , var%q_var_bas         , &
            dx1_L(i,j,k)         , dx2_R(i,j,k)          , &
            var%cx1(i,j,k)       , var%ax2(i,j,k)        , &
            var%ax1(i,j,k)       , var%cx2(i,j,k)        , &
            var%d(i,j,k)         , var%b(i,j,k)          , &
            var%ax1(i,j,k), var%b1(i,j,k), var%cx1(i,j,k), &
            var%ax2(i,j,k), var%b2(i,j,k), var%cx2(i,j,k), &
            var%ax3(i,j,k), var%b3(i,j,k), var%cx3(i,j,k)  ) 
    end do
    
#endif
     

  end subroutine coefs_UT_CL

  

!!!!!---------------------------------- 
  subroutine CL_Dirichlet (val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3, b3, cx3, var) 
    implicit none
    
    real(nk),intent(in   ):: val
    real(nk),intent(inout):: d, b, b1, b2, ax1, cx1, ax2, cx2
    real(nk),intent(inout),optional::  b3, ax3, cx3, var
    d   =  val  
    b   =  0._nk
    b1  =  1._nk
    b2  =  1._nk
    ax1 =  0._nk
    cx1 =  0._nk
    ax2 =  0._nk
    cx2 =  0._nk
#if(DIMENSION_GEO == 3)
    b3  =  1._nk
    ax3 =  0._nk
    cx3 =  0._nk
#endif
    if (present(var)) then         
       var = val
       d   = 0._nk
    end if
  end subroutine CL_Dirichlet
!!!!!-------------------------------------
  !(var%d(i,j,k), var%cx1(i,j,k), var%ax1(i,j,k), dx1, var%q_var_gauche) 
  subroutine CL_Neumann (d, ac, coef, dx, q) 
    implicit none
    real(nk),intent(inout):: d, ac, coef
    real(nk),intent(in)::dx, q
         
    d    =  d     +  coef * 2._nk * dx * q
    ac   =  ac    +  coef
    coef =  0.0_nk     
    
  end subroutine CL_Neumann






!!!!!---------------------- 





  
  subroutine CL_coin (CL_1 , CL_2 , &
                      val_1, val_2, &
                      q_1  , q_2  , &
                      dx_1 , dx_2 , &
                      coef1, coef2, &
                      ac1  , ac2  , &
                      d    , b    , &
                      ax1  , b1 , cx1  , &
                      ax2  , b2 , cx2  , &
                      ax3  , b3 , cx3  , &
                      CL_3 , val_3, &
                      q_3  , dx_3 , &
                      coef3, ac3  , &
                      varia       )
    implicit none
    
    integer ,intent(in)   ::  CL_1,  CL_2 
    real(nk),intent(in)   ::  val_1, val_2, q_1, q_2, dx_1, dx_2
    real(nk),intent(inout)::  coef1, coef2
    real(nk),intent(inout)::  ac1  , ac2
    real(nk),intent(inout)::  d, b, ax1, b1, cx1, ax2, b2, cx2
 
    real(nk),intent(inout),optional:: ax3, b3, cx3
    
    integer ,intent(in)   ,optional:: CL_3
    real(nk),intent(in)   ,optional:: val_3, q_3, dx_3
    real(nk),intent(inout),optional:: coef3, ac3, varia
    
    real(nk) :: val=0., val1=0., val2=0., val3=0.



#if (DIMENSION_GEO == 2) 
      
    if ( CL_1  == 2 .and.  CL_2 == 2) then
       call CL_Neumann ( d, ac1, coef1, dx_1, q_1)
       call CL_Neumann ( d, ac2, coef2, dx_2, q_2)
    else if ( CL_1  == 1 .and.  CL_2 == 1) then  
       if (val_1 /= val_2 ) then
          val = max(val_1, val_2)
       else 
          val =  val_1
       end if
       if (present(varia)) then 
          call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, var=varia)       
       else
          call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2)       
       end if
    else if ( CL_1 /=  CL_2 ) then  
       if ( CL_1 == 1 ) val = val_1 
       if ( CL_2 == 1 ) val = val_2
       if (present(varia)) then 
          call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, var=varia)       
       else
          call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2)       
       end if
    end if

#elif (DIMENSION_GEO == 3) 

    if (present(CL_3)) then 
       if ( CL_1  == 2 .and.  CL_2 == 2 .and. CL_3 == 2 ) then
          call CL_Neumann ( d, ac1, coef1, dx_1, q_1)
          call CL_Neumann ( d, ac2, coef2, dx_2, q_2)
          call CL_Neumann ( d, ac3, coef3, dx_3, q_3)       
       else if ( CL_1  == 1 .and. CL_2 == 1 .and. CL_3 == 1) then  
          if (val_1 /= val_2 .or. val_1 /= val_3 .or. val_2 /= val_3) then
             val = max(val_1, val_2, val_3)
          else 
             val =  val_1
          end if
          if (present(varia)) then
             call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3,  b3, cx3, var=varia)
          else
             call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3,  b3, cx3)
          end if
       else if ( CL_1 /=  CL_2 .or. CL_1 /=  CL_3 .or. CL_2 /=  CL_3 ) then  
          if ( CL_1 == 1 ) val1 = val_1 
          if ( CL_2 == 1 ) val2 = val_2       
          if ( CL_3 == 1 ) val3 = val_3
          val = max(val1, val2, val3)
          if (present(varia)) then
             call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3,  b3, cx3, var=varia)
          else
             call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3,  b3, cx3)
          end if
       end if
       
    else
       
       if ( CL_1  == 2 .and.  CL_2 == 2) then
          call CL_Neumann ( d, ac1, coef1, dx_1, q_1)
          call CL_Neumann ( d, ac2, coef2, dx_2, q_2)
       else if ( CL_1  == 1 .and.  CL_2 == 1) then  
          val = max(val_1, val_2)
          if (present(varia)) then 
             call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3,  b3, cx3, var=varia)
          else
             call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3,  b3, cx3)
          end if
       else if ( CL_1 /=  CL_2 ) then  
          if ( CL_1 == 1 ) val = val_1 
          if ( CL_2 == 1 ) val = val_2
          if (present(varia)) then 
             call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3,  b3, cx3, var=varia)
          else
             call CL_Dirichlet(val, d, b, ax1, b1, cx1, ax2, b2, cx2, ax3,  b3, cx3)
          end if
       end if
    
    end if

#endif
    
  end subroutine CL_coin
!!!!!----------------------   
  

  

  

  
  subroutine UT_points_ficts ( var )
    implicit none
    type(variable),intent(inout):: var
    integer:: i, j, k

#if (DIMENSION_GEO == 2) 
   do i= 1, N1
!!!--------------------------------------!! Paroi_bas !!--------------------------------------------
      j=1
      if ( var%CL_var_bas == 1 ) var%var(i,j-1) = - var%var(i,j+1) +  2. * var%var_bas
      if ( var%CL_var_bas == 2 ) var%var(i,j-1) =   var%var(i,j+1) -  2. * dx2_R(i,j) * var%q_var_bas

!!!--------------------------------------!! Paroi_haut !!--------------------------------------------
      j=N2
      if ( var%CL_var_haut == 1 ) var%var(i,j+1) = - var%var(i,j-1) +   2. * var%var_haut  
      if ( var%CL_var_haut == 2 ) var%var(i,j+1) =   var%var(i,j-1) +   2. * dx2_L(i,j) * var%q_var_haut  
   end do
      
   do j= 1, N2
!!!--------------------------------------!! Paroi_gauche !!--------------------------------------------
      i=1
      if ( var%CL_var_gauche == 1 ) var%var(i-1,j) = - var%var(i+1,j) +  2. * var%var_gauche
      if ( var%CL_var_gauche == 2 ) var%var(i-1,j) =   var%var(i+1,j) -  2. * dx1_R(i,j) * var%q_var_gauche   

!!!--------------------------------------!! Paroi_droite !!--------------------------------------------
      i=N1
      if ( var%CL_var_droite == 1 ) var%var(i+1,j) = - var%var(i-1,j) +  2. * var%var_droite   
      if ( var%CL_var_droite == 2 ) var%var(i+1,j) =   var%var(i-1,j) +  2. * dx1_L(i,j) * var%q_var_droite
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#elif (DIMENSION_GEO == 3)

!!!--------------------------------------!! Paroi_droite !!--------------------------------------------        
    if ( var%CL_var_droite == 1 ) var%var(:,:,0) = - var%var(:,:,2) +  2. * var%var_droite   
    if ( var%CL_var_droite == 2 ) var%var(:,:,0) =   var%var(:,:,2) +  2. * dx3_R(:,:,1) * var%q_var_droite
!!!--------------------------------------!! Paroi_gauche !!--------------------------------------------       
    if ( var%CL_var_gauche == 1 ) var%var(:,:,N3+1) = - var%var(:,:,N3-1) +  2. * var%var_gauche
    if ( var%CL_var_gauche == 2 ) var%var(:,:,N3+1) =   var%var(:,:,N3-1) -  2. * dx3_L(:,:,N3) * var%q_var_gauche
   
!!!------------------------------------------!! Paroi_bas !!---------------------------------------------
    if ( var%CL_var_bas == 1 ) var%var(:,0,:) = - var%var(:,2,:) + 2. * var%var_bas       
    if ( var%CL_var_bas == 2 ) var%var(:,0,:) =   var%var(:,2,:) - 2. * dx2_R(:,1,:) * var%q_var_bas
!!!------------------------------------------!! Paroi_haut !!--------------------------------------------
    if ( var%CL_var_haut == 1 ) var%var(:,N2+1,:) = - var%var(:,N2-1,:) +  2. * var%var_haut  
    if ( var%CL_var_haut == 2 ) var%var(:,N2+1,:) =   var%var(:,N2-1,:) +  2. * dx2_L(:,N2,:) * var%q_var_haut  
   
!!!------------------------------------------!! Paroi_avant !!--------------------------------------------               
    if ( var%CL_var_avant == 1 ) var%var(0,:,:) = - var%var(2,:,:) +  2. * var%var_avant
    if ( var%CL_var_avant == 2 ) var%var(0,:,:) =   var%var(2,:,:) -  2. * dx1_R(1,:,:) * var%q_var_avant
!!!------------------------------------------!! Paroi_arriere !!-----------------------------------------
    if ( var%CL_var_arriere == 1 ) var%var(N1+1,:,:) = - var%var(N1-1,:,:) +  2. * var%var_arriere
    if ( var%CL_var_arriere == 2 ) var%var(N1+1,:,:) =   var%var(N1-1,:,:) +  2. * dx1_L(N1,:,:) * var%q_var_arriere 

    
#endif     

  end subroutine UT_points_ficts




    




  subroutine initialisation_champ_T(T)
   integer:: j
    type(variable), intent(inout):: T
    real(nk):: DeltT

    DeltT = T%var_bas - T%var_haut

#if (  DIMENSION_GEO   == 2  )
#if ( UNIFORM_GRID == 0)
    do j = 0, N2-1
      T%var(:,j) = (T%var_bas - yj/L2*DeltT)
    enddo
    T%var (:,1 ) = T%var_bas
    T%var (:,N2) = T%var_haut
    T%var0(:,:) = T%var(:,:)
#elif ( UNIFORM_GRID == 1)
    do j = 2, N2-1
      T%var(:,j) = T%var_bas - DeltT/(N2-1)*(j-1)
    enddo; T%var0(:,:) = T%var(:,:)
    T%var (N1/2  , N2/2) = T%var (N1/2  , N2/2) + 5e-3
    T%var (N1/8  , N2/2) = T%var (N1/8  , N2/2) - 2e-3
    T%var (N1/8*7, N2/2) = T%var (N1/8*7, N2/2) - 2e-3

#endif

#elif( DIMENSION_GEO  == 3  )
    do j = 2, N2-1
      T%var(:,j,:) = T%var_bas - DeltT/(N2-1)*(j-1)
    enddo; T%var0(:,:,:) = T%var(:,:,:)
    T%var (N1/2,N2/2,N3/2) = T%var0(N1/2,N2/2,N3/2) + 5e-2
    T%var (N1/2,N2/2,N3/8) = T%var0(N1/2,N2/2,N3/8) - 2e-2

#endif
  end subroutine initialisation_champ_T
end module mVitesse
















 

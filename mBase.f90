#include "definitions.h"

module mBase

  
  implicit none

!!$                                                                    2D                         3D  
!!$                      .. . . . . . . . . . .+                  ----------------------------------------
!!$                    . .                  . .                    W(1) = u1                  W(1)  = u1
!!$                 .    .                .   .                    W(2) = u2                  W(2)  = u2
!!$               .      .              .     .                                               W(3)  = u3  
!!$             ._____________________.       .                    W(3) = tau11               W(4)  = tau11
!!$             |        .            |       .  x2                W(4) = tau12               W(5)  = tau12
!!$             |        .            |       .                                               W(6)  = tau13
!!$             |        .            |       .                    W(5) = tau22               W(7)  = tau22
!!$             |        .            |       .                                               W(8)  = tau23
!!$             |        . . . . . . .|. . . ..                                               W(9)  = tau33
!!$             |       .             |      .                     W(6) = T                   W(10) = T
!!$             |     .               |    .                       W(7) = P                   W(11) = P 
!!$             |   .                 |   .    x1                             
!!$             | .                   | .  
!!$             |_____________________|.+        
!!$             +        
!!$                      x3 

  ! Declarations 
#if ( DIMENSION_GEO == 3 ) 
  integer :: u1=1, u2=2, u3=3, t11=4, t12=5, t13=6, t22=7, t23=8, t33=9, T=10, P=11
#elif( DIMENSION_GEO == 2 ) 
  integer :: u1=1, u2=2, t11=3, t12=4, t22=5, T=6, P=7
#endif
  
  integer, public, parameter :: sk=selected_real_kind(  6, 37   )
  integer, public, parameter :: nk=selected_real_kind( 15, 307  )
  integer, public, parameter :: pk=selected_real_kind( 33, 4931 )

  real(nk),parameter::pii=acos(-1.0_nk), sqrt2= sqrt(2.)
  integer :: i, j, l, k, s, ll, ii, jj 
  integer :: it, nbr_it_espace, nbr_it_temps, print_erreur_temps, pas_choisi_temps, pas_write
  integer :: pas=1

  ! grid size
  integer :: NNN
  integer :: N1
  integer :: N2
  real(nk):: dt, L1, L2
  integer :: cord_x1_0=0, cord_x2_0=0
  
#if ( DIMENSION_GEO == 3)
  integer :: N3
  real(nk):: L3
  integer :: cord_x3_0 = 0
#endif

#if ( DIMENSION_GEO == 2 )
  real(nk), dimension(:,:), allocatable    :: dx1_L, dx1_1_L, dx1_2_L
  real(nk), dimension(:,:), allocatable    :: dx1_R, dx1_1_R, dx1_2_R
  real(nk), dimension(:,:), allocatable    :: dx2_L, dx2_1_L, dx2_2_L
  real(nk), dimension(:,:), allocatable    :: dx2_R, dx2_1_R, dx2_2_R
  real(nk), dimension(:,:), allocatable    :: E_global, D_diffusion, F_thermal, V_dissipation, G_exchange
  real(nk)                                 :: Last_Global_E, Last_Total_Tau

#elif ( DIMENSION_GEO == 3 )
  real(nk), dimension(:,:,:), allocatable    :: dx1_L, dx1_1_L, dx1_2_L
  real(nk), dimension(:,:,:), allocatable    :: dx1_R, dx1_1_R, dx1_2_R
  real(nk), dimension(:,:,:), allocatable    :: dx2_L, dx2_1_L, dx2_2_L
  real(nk), dimension(:,:,:), allocatable    :: dx2_R, dx2_1_R, dx2_2_R
  real(nk), dimension(:,:,:), allocatable    :: dx3_L, dx3_1_L, dx3_2_L
  real(nk), dimension(:,:,:), allocatable    :: dx3_R, dx3_1_R, dx3_2_R

  real(nk), dimension(:,:,:), allocatable    :: E_global, D_diffusion, F_thermal, V_dissipation, G_exchange
  real(nk)                                   :: Last_Global_E, Last_Total_Tau

#endif
 
  ! grande structure des variables
#if (DIMENSION_GEO == 2)
  type variable
     real(nk),dimension(:,:),allocatable:: ax1, ax2, b, b1, b2, cx1, cx2, d , d0 , var, var0, var1, var10
     integer  :: CL_var_bas,   CL_var_haut,   CL_var_droite,   CL_var_gauche
     real(nk) ::    var_bas,      var_haut,      var_droite,      var_gauche
     real(nk) ::  q_var_bas,    q_var_haut,    q_var_droite,    q_var_gauche
     real(nk) ::    nomb_ad,       rapp_ad
     ! que pour le solver 
     integer  ::        NNN,            NI,               NJ
     integer  ::         N1,            N2
     integer  ::     indice 
     real(nk) ::     niveau_conver
  end type variable
#endif
   
#if (DIMENSION_GEO == 3)
  type variable
     real(nk),dimension(:,:,:),allocatable:: ax1, ax2, ax3, b1, b2, b3, cx1, cx2, cx3, b
     real(nk),dimension(:,:,:),allocatable:: d , d0 , var, var0, var1, var10
     integer  :: CL_var_bas,   CL_var_haut,   CL_var_droite,   CL_var_gauche,       CL_var_avant,   CL_var_arriere   
     real(nk) ::    var_bas,      var_haut,      var_droite,      var_gauche,          var_avant,      var_arriere
     real(nk) ::  q_var_bas,    q_var_haut,    q_var_droite,    q_var_gauche,        q_var_avant,    q_var_arriere
     real(nk) ::    nomb_ad,       rapp_ad
     ! que pour le solver 
     integer  ::        NNN,            NI,               NJ,              NK
     integer  ::         N1,            N2,               N3
     integer  ::     indice 
     real(nk) ::     niveau_conver
  end type variable
#endif
#if (DIMENSION_GEO == 2)
  real(nk),dimension(:,:),allocatable:: var_intf_u1, var_intf_u2, var_intf_u1_10, var_intf_u2_10
#elif (DIMENSION_GEO == 3) 
  real(nk),dimension(:,:,:),allocatable:: var_intf_u1, var_intf_u2, var_intf_u3
  real(nk),dimension(:,:,:),allocatable:: var_intf_u1_10, var_intf_u2_10, var_intf_u3_10
#endif

  ! pour l'etude des érreurs
  real(nk),dimension(:),allocatable:: time_size , grid_size
  type etude_ordre
     real(nk),dimension(:),allocatable:: val_var_espace, val_var_temps, ord_espace, ord_temps 
  end type etude_ordre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(nk) :: Total_Nu_paroi1_0, Total_Nu_paroi1_00
  real(nk) :: Time_period_integral_Nu, Period_t1=0., Period_t=0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#if ( EXP_TREATMENT == 1 )
  real(nk) :: Tau_trans=1.
  integer  :: Re_zero=0.
!#endif



#if (DIMENSION_GEO ==2)
  type vecteur_variable
     type(variable), dimension(7)  :: W     
  end type vecteur_variable
  type Erreur_var
       type(etude_ordre),dimension(7 ) :: var
    end type Erreur_var
#endif
#if (DIMENSION_GEO ==3)
  type vecteur_variable
     type(variable), dimension(11) :: W
  end type vecteur_variable
  type Erreur_var
       type(etude_ordre),dimension(11) :: var 
    end type Erreur_var
#endif
  type Erreur_var_norme
     type(erreur_var),dimension(3) :: norme
  end type Erreur_var_norme
  
  ! Declarations 
  type(vecteur_variable):: var
  type(Erreur_var      ):: Er
  type(Erreur_var_norme):: Er_norme  
  
!  var(7)%UI = N1-1; var(7)%UJ = N2-1
  
  ! hyperbolic part : termes non-linéaires  
  ! sans température
  ! Lambda : matrice des valeurs propres 
  ! L : matrice des vecteurs propres à gauche (Left)
  ! R : matrice des vecteurs propres à droite (Right)
  type val_vec_pro  
#if (DIMENSION_GEO == 2)
     real(nk),dimension(:,:,:,:),allocatable:: L, R
     real(nk),dimension(:,:,:)  ,allocatable:: Lambda, NL, NL0, Lambda_dRWdx
#elif (DIMENSION_GEO == 3)
     real(nk),dimension(:,:,:,:,:),allocatable::  L, R
     real(nk),dimension(:,:,:,:)  ,allocatable:: Lambda, NL, NL0, Lambda_dRWdx
#endif 
  end type val_vec_pro

#if (DIMENSION_GEO == 2)
     type(val_vec_pro),dimension(2):: A
#elif (DIMENSION_GEO == 3)
     type(val_vec_pro),dimension(3):: A
#endif 


  
!#if ( VISCOELASTIC_MODELS > 0 )
!#if (DIMENSION_GEO == 2)
!  type(mat_coefs),dimension(:,:),allocatable::mat 
!#elif (DIMENSION_GEO == 3)
!  type(mat_coefs),dimension(:,:,:),allocatable::mat 
!#endif
!#endif 

  integer :: sonde
  real(nk):: Peclet, Grashof, Prandtl, Rayleigh
  real(nk):: Ma2, Reynolds, Weissenberg, Beta,  CFL, epsilon_PPT , alpha_G,  xi_PPT
  real(nk):: div_U_sur_coins_avec_coins=0., div_U_sur_centres_avec_Uinter=0., eps_zero, &
       div_U_sur_centres_avec_coins=0.
  

#if (DIMENSION_GEO == 2)
       integer,dimension(:,:)  ,allocatable:: Cent_or_HOUC_A, Cent_or_HOUC_B
#elif (DIMENSION_GEO == 3)
       integer,dimension(:,:,:),allocatable:: Cent_or_HOUC_A, Cent_or_HOUC_B, Cent_or_HOUC_C
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if (DIMENSION_GEO == 2)
       integer,dimension(:,:)  ,allocatable:: UPW2_1_A, UPW2_1_B, UPW2_1_C, BACKW2_1_A, BACKW2_1_b, BACKW2_1_C
       integer,dimension(:,:)  ,allocatable:: UPW2_2_A, UPW2_2_B, UPW2_2_C, BACKW2_2_A, BACKW2_2_b, BACKW2_2_C
       integer,dimension(:,:)  ,allocatable:: CEN2_1_A, CEN2_1_B, CEN2_1_C, CEN2_2_A, CEN2_2_B, CEN2_2_C
       integer,dimension(:,:)  ,allocatable:: CEN4_1_A, CEN4_1_B, CEN4_1_C, CEN4_1_D, CEN4_1_E
       integer,dimension(:,:)  ,allocatable:: CEN4_2_A, CEN4_2_B, CEN4_2_C, CEN4_2_D, CEN4_2_E
#elif (DIMENSION_GEO == 3)
       integer,dimension(:,:,:)  ,allocatable:: UPW2_1_A, UPW2_1_B, UPW2_1_C, BACKW2_1_A, BACKW2_1_b, BACKW2_1_C
       integer,dimension(:,:,:)  ,allocatable:: UPW2_2_A, UPW2_2_B, UPW2_2_C, BACKW2_2_A, BACKW2_2_b, BACKW2_2_C
       integer,dimension(:,:,:)  ,allocatable:: UPW2_3_A, UPW2_3_B, UPW2_3_C, BACKW2_3_A, BACKW2_3_b, BACKW2_3_C
       integer,dimension(:,:,:)  ,allocatable:: CEN2_1_A, CEN2_1_B, CEN2_1_C, CEN2_2_A, CEN2_2_B, CEN2_2_C
       integer,dimension(:,:,:)  ,allocatable:: CEN2_3_A, CEN2_3_B, CEN2_3_C
       integer,dimension(:,:,:)  ,allocatable:: CEN4_1_A, CEN4_1_B, CEN4_1_C, CEN4_1_D, CEN4_1_E
       integer,dimension(:,:,:)  ,allocatable:: CEN4_2_A, CEN4_2_B, CEN4_2_C, CEN4_2_D, CEN4_2_E
       integer,dimension(:,:,:)  ,allocatable:: CEN4_3_A, CEN4_3_B, CEN4_3_C, CEN4_3_D, CEN4_3_E
#endif


#if (TYPE_ECOULEMENT == 3) 
  integer ::  T_G_T = 0, C_E =0, C_D_C = 1, Poiseuille = 0
#endif 


 integer :: ittt
 integer :: pas_temps_adptatif
 integer :: pas_write0, Total_step

 real(nk):: temps_cumule=0.
 real(nk):: eps_iteratif_P
 real(nk):: dp_x
 real(nk)::work
 real(nk):: resultats_int2
 real(nk),dimension(:),allocatable:: T_moy, Vitesse_u1_moy


 
 

#if (DIMENSION_GEO ==2)
 TYPE SHAMPOING
 REAL(nk) :: CH(6)
 END TYPE SHAMPOING
#elif (DIMENSION_GEO ==3)
 TYPE SHAMPOING
 REAL(nk) :: CH(9)
 END TYPE SHAMPOING
#endif
#if (DIMENSION_GEO == 2)
 TYPE(SHAMPOING), DIMENSION(:,:), ALLOCATABLE :: AP
 TYPE(SHAMPOING), DIMENSION(:,:), ALLOCATABLE :: CH_P
 REAL(nk), DIMENSION(:,:), ALLOCATABLE :: VPRX, VPRX1, OPERPX
 REAL(nk), DIMENSION(:,:), ALLOCATABLE :: DIV
 REAL(nk), DIMENSION(:), ALLOCATABLE :: VPX
 INTEGER :: IDU, IFU, JDU, JFU
#endif

#if (DIMENSION_GEO == 3)
 TYPE(SHAMPOING), DIMENSION(:,:,:), ALLOCATABLE :: AP
 TYPE(SHAMPOING), DIMENSION(:,:,:), ALLOCATABLE :: CH_P
 REAL(nk), DIMENSION(:,:), ALLOCATABLE :: VPRX, VPRX1, OPERPX
 REAL(nk), DIMENSION(:,:), ALLOCATABLE :: VPRZ, VPRZ1, OPERPZ
 REAL(nk), DIMENSION(:,:,:), ALLOCATABLE :: DIV
 REAL(nk), DIMENSION(:), ALLOCATABLE :: VPX
 REAL(nk), DIMENSION(:), ALLOCATABLE :: VPZ
 INTEGER :: IDU, IFU, JDU, JFU, KDU, KFU
#endif

  !!!!   time points
  real(nk) :: time1,  time2,  time3,  time4,  time5,  time6,  time7,  time8,  time9,  time10

contains

  SUBROUTINE ALLOUER
    IMPLICIT NONE
    INTEGER :: I
#if (DIMENSION_GEO ==2)
    ALLOCATE ( AP(0:N1,0:N2) )
    DO i = 1 , 6
      AP(:,:)%CH(i) = 0.
    ENDDO; 

    ALLOCATE ( CH_P(0:N1,0:N2) )
    DO i = 1 , 6
      CH_P(:,:)%CH(i) = 0.
    ENDDO
! ===== VALEURS ET VECTEURS PROPRES
    ALLOCATE ( VPX  (0:N1) ); VPX = 0.
    ALLOCATE ( VPRX (0:N1,0:N1), VPRX1(0:N1,0:N1) ); VPRX1= 0.; VPRX = 0.
    ALLOCATE ( OPERPX(0:N1,0:N1) ); OPERPX = 0.; 
    ALLOCATE ( DIV(0:N1,0:N2) ); DIV = 0.; 


    allocate( dx1_L     (0:N1+1,0:N2+1) )
    allocate( dx1_1_L   (0:N1+1,0:N2+1) )
    allocate( dx1_2_L   (0:N1+1,0:N2+1) )
    allocate( dx1_R     (0:N1+1,0:N2+1) )
    allocate( dx1_1_R   (0:N1+1,0:N2+1) )
    allocate( dx1_2_R   (0:N1+1,0:N2+1) )
    allocate( dx2_L     (0:N1+1,0:N2+1) )
    allocate( dx2_1_L   (0:N1+1,0:N2+1) )
    allocate( dx2_2_L   (0:N1+1,0:N2+1) )
    allocate( dx2_R     (0:N1+1,0:N2+1) )
    allocate( dx2_1_R   (0:N1+1,0:N2+1) )
    allocate( dx2_2_R   (0:N1+1,0:N2+1) )
#endif

#if (DIMENSION_GEO == 3)
    ALLOCATE ( AP(0:N1,0:N2,0:N3) )
    DO i = 1 , 9
      AP(:,:,:)%CH(i) = 0.
    ENDDO; 

    ALLOCATE ( CH_P(0:N1,0:N2,0:N3) )
    DO i = 1 , 6
      CH_P(:,:,:)%CH(i) = 0.
    ENDDO
! ===== VALEURS ET VECTEURS PROPRES
    ALLOCATE ( VPX  (0:N1),      VPZ  (0:N3) ); VPX = 0. ; VPZ = 0.
    ALLOCATE ( VPRX (0:N1,0:N1), VPRX1(0:N1,0:N1) ); VPRX1= 0.; VPRX = 0.
    ALLOCATE ( VPRZ (0:N3,0:N3), VPRZ1(0:N3,0:N3) ); VPRZ1= 0.; VPRZ = 0.
    ALLOCATE ( OPERPX(0:N1,0:N1) ); OPERPX = 0.; 
    ALLOCATE ( OPERPZ(0:N3,0:N3) ); OPERPZ = 0.; 
    ALLOCATE ( DIV(0:N1,0:N2,0:N3) ); DIV = 0.; 


    allocate( dx1_L     (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx1_1_L   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx1_2_L   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx1_R     (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx1_1_R   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx1_2_R   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx2_L     (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx2_1_L   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx2_2_L   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx2_R     (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx2_1_R   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx2_2_R   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx3_L     (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx3_1_L   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx3_2_L   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx3_R     (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx3_1_R   (0:N1+1,0:N2+1,0:N3+1) )
    allocate( dx3_2_R   (0:N1+1,0:N2+1,0:N3+1) )

#endif

  END SUBROUTINE ALLOUER

  SUBROUTINE DEALLOUER
    IMPLICIT NONE
    INTEGER :: I

#if (DIMENSION_GEO == 2)
     DEALLOCATE ( AP )
     DEALLOCATE ( CH_P )
     DEALLOCATE ( VPX )
     DEALLOCATE ( VPRX, VPRX1 )
     DEALLOCATE ( OPERPX )
     DEALLOCATE ( DIV )

     
#elif (DIMENSION_GEO == 3)     
     DEALLOCATE ( AP )
     DEALLOCATE ( CH_P )
     DEALLOCATE ( VPX, VPZ )
     DEALLOCATE ( VPRX, VPRX1, VPRZ, VPRZ1 )
     DEALLOCATE ( OPERPX, OPERPZ )
     DEALLOCATE ( DIV )

#endif


#if (DIMENSION_GEO ==2)

    deallocate( dx1_L     )
    deallocate( dx1_1_L   )
    deallocate( dx1_2_L   )
    deallocate( dx1_R     )
    deallocate( dx1_1_R   )
    deallocate( dx1_2_R   )
    deallocate( dx2_L     )
    deallocate( dx2_1_L   )
    deallocate( dx2_2_L   )
    deallocate( dx2_R     )
    deallocate( dx2_1_R   )
    deallocate( dx2_2_R   )
#endif

#if (DIMENSION_GEO == 3)

    deallocate( dx1_L     )
    deallocate( dx1_1_L   )
    deallocate( dx1_2_L   )
    deallocate( dx1_R     )
    deallocate( dx1_1_R   )
    deallocate( dx1_2_R   )
    deallocate( dx2_L     )
    deallocate( dx2_1_L   )
    deallocate( dx2_2_L   )
    deallocate( dx2_R     )
    deallocate( dx2_1_R   )
    deallocate( dx2_2_R   )
    deallocate( dx3_L     )
    deallocate( dx3_1_L   )
    deallocate( dx3_2_L   )
    deallocate( dx3_R     )
    deallocate( dx3_1_R   )
    deallocate( dx3_2_R   )

#endif

  END SUBROUTINE DEALLOUER


  subroutine allocate_var( var )
    implicit none
    type(variable):: var


#if ( DIMENSION_GEO == 2 ) 
    allocate( var%ax1  (0:var%NI+1,0:var%NJ+1) )
    allocate( var%ax2  (0:var%NI+1,0:var%NJ+1) )
    allocate( var%cx1  (0:var%NI+1,0:var%NJ+1) )
    allocate( var%cx2  (0:var%NI+1,0:var%NJ+1) )
    allocate( var%b    (0:var%NI+1,0:var%NJ+1) )
    allocate( var%b1   (0:var%NI+1,0:var%NJ+1) )
    allocate( var%b2   (0:var%NI+1,0:var%NJ+1) )
    allocate( var%d    (0:var%NI+1,0:var%NJ+1) )
!!!--
    allocate( var%d0   (0:var%NI+1,0:var%NJ+1) )
!!!--
    allocate( var%var  (0:var%NI+1,0:var%NJ+1) )
    allocate( var%var0 (0:var%NI+1,0:var%NJ+1) )
    allocate( var%var1 (0:var%NI+1,0:var%NJ+1) )
    allocate( var%var10(0:var%NI+1,0:var%NJ+1) )
!!!--
#endif
#if ( DIMENSION_GEO == 3 )
    allocate( var%ax1  (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%ax2  (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%cx1  (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%cx2  (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%ax3  (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%cx3  (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%b    (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%b1   (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%b2   (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%b3   (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%d    (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
!!!--
    allocate( var%d0   (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
!!!--
    allocate( var%var  (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%var0 (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%var1 (0:var%NI+1,0:var%NJ+1,0:var%NK+1) )
    allocate( var%var10(0:var%NI+1,0:var%NJ+1,0:var%NK+1) )

!!!--
#endif
  end subroutine allocate_var
  
  
  
  subroutine deallocate_var( var )
    implicit none
    type(variable):: var
    deallocate( var%ax1 )
    deallocate( var%ax2 )
    deallocate( var%cx1 )
    deallocate( var%cx2 )
    deallocate( var%b   )
    deallocate( var%b1  )
    deallocate( var%b2  )
    deallocate( var%d   )
!!!--
    deallocate( var%d0  )
!!!--
    deallocate( var%var, var%var0, var%var1, var%var10 )

    
#if ( DIMENSION_GEO == 3 ) 
    deallocate( var%ax3 )
    deallocate( var%cx3 )
    deallocate( var%b3  )
#endif    
  end subroutine deallocate_var



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine grille_et_pas_maillage
    implicit none
    integer :: i, j, k
    real(nk):: Sita, a, b, coef_y=0.95_nk, yj=0._nk, yj_0=0._nk
    real(nk)::             coef_x=0.95_nk, xi=0._nk, xi_0=0._nk
#if ( UNIFORM_GRID == 0)
    do i = 0, N1-1
      xi_0 = xi
      a = i ; b = N1-1
      Sita = -1._nk + 2._nk * (a/b)
      xi = 1._nk / coef_x * tanh(.5_nk * Sita * log((1+coef_x)/(1-coef_x)))
      xi = xi/ 2._nk + 0.5_nk
      xi = xi * L1
      
    dx1_L(i+1,:)   = xi - xi_0
    dx1_1_L(i+1,:) = 1._nk / dx1_L(i+1,:)
    dx1_2_L(i+1,:) = 1._nk / dx1_L(i+1,:)**2  
    enddo; dx1_L(1   ,:) = dx1_L(2, :)
           dx1_L(N1+1,:) = dx1_L(N1,:)

    do i = 0, N1
      dx1_R  (i,:)   =       dx1_L(i+1,:)
      dx1_1_R(i,:) = 1._nk / dx1_R(i,:)
      dx1_2_R(i,:) = 1._nk / dx1_R(i,:)**2
    enddo; dx1_R(N1,:)    = dx1_L(N1,:)
           dx1_R(N1+1,:)  = dx1_R(N1,:)
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j = 0, N2-1
      yj_0 = yj
      a = j ; b = N2-1
      Sita = -1._nk + 2._nk * (a/b)
      yj = 1._nk / coef_y * tanh(.5_nk * Sita * log((1+coef_y)/(1-coef_y)))
      yj = yj / 2._nk + 0.5_nk
      yj = yj * L2
      
    dx2_L(:,j+1)   = yj - yj_0
    dx2_1_L(:,j+1) = 1._nk / dx2_L(:,j+1)
    dx2_2_L(:,j+1) = 1._nk / dx2_L(:,j+1)**2  
    enddo; dx2_L(:, 1  ) = dx2_L(:, 2)
           dx2_L(:,N2+1) = dx2_L(:,N2)
    
    do j = 0, N2
    dx2_R  (:,j)   =       dx2_L(:,j+1)
    dx2_1_R(:,j) = 1._nk / dx2_R(:,j)
    dx2_2_R(:,j) = 1._nk / dx2_R(:,j)**2
    enddo; dx2_R(:,N2  )    = dx2_L(:,N2)
           dx2_R(:,N2+1)    = dx2_L(:,N2)

#elif ( UNIFORM_GRID == 1 )
#if ( DIMENSION_GEO == 2 )

    dx1_L(:,:)   = L1    / real(N1-1,nk)
    dx1_R(:,:)   = L1    / real(N1-1,nk)
    dx1_1_L(:,:) = 1._nk / dx1_L(:,:)
    dx1_1_R(:,:) = 1._nk / dx1_R(:,:)
    dx1_2_L(:,:) = 1._nk / dx1_L(:,:)**2
    dx1_2_R(:,:) = 1._nk / dx1_R(:,:)**2

    dx2_L(:,:)   = L2    / real(N2-1,nk)
    dx2_R(:,:)   = L2    / real(N2-1,nk)
    dx2_1_L(:,:) = 1._nk / dx2_L(:,:)
    dx2_1_R(:,:) = 1._nk / dx2_R(:,:)
    dx2_2_L(:,:) = 1._nk / dx2_L(:,:)**2
    dx2_2_R(:,:) = 1._nk / dx2_R(:,:)**2
#elif ( DIMENSION_GEO == 3 )
    dx1_L(:,:,:)   = L1    / real(N1-1,nk)
    dx1_R(:,:,:)   = L1    / real(N1-1,nk)
    dx1_1_L(:,:,:) = 1._nk / dx1_L(:,:,:)
    dx1_1_R(:,:,:) = 1._nk / dx1_R(:,:,:)
    dx1_2_L(:,:,:) = 1._nk / dx1_L(:,:,:)**2
    dx1_2_R(:,:,:) = 1._nk / dx1_R(:,:,:)**2

    dx2_L(:,:,:)   = L2    / real(N2-1,nk)
    dx2_R(:,:,:)   = L2    / real(N2-1,nk)
    dx2_1_L(:,:,:) = 1._nk / dx2_L(:,:,:)
    dx2_1_R(:,:,:) = 1._nk / dx2_R(:,:,:)
    dx2_2_L(:,:,:) = 1._nk / dx2_L(:,:,:)**2
    dx2_2_R(:,:,:) = 1._nk / dx2_R(:,:,:)**2

    dx3_L(:,:,:)   = L3    / real(N3-1,nk)
    dx3_R(:,:,:)   = L3    / real(N3-1,nk)
    dx3_1_L(:,:,:) = 1._nk / dx3_L(:,:,:)
    dx3_1_R(:,:,:) = 1._nk / dx3_R(:,:,:)
    dx3_2_L(:,:,:) = 1._nk / dx3_L(:,:,:)**2
    dx3_2_R(:,:,:) = 1._nk / dx3_R(:,:,:)**2
#endif
#endif


  end subroutine grille_et_pas_maillage
     






  subroutine thomas2dx (a,b,c,d,nx,ny,x)
    implicit none
    integer:: i,ii,j,nx,ny
    !integer, parameter::np=selected_real_kind(p=18,r=307)
    real(nk),dimension(nx,ny),intent(in):: a,b,c,d
    real(nk),dimension(nx,ny),intent(out):: x
    real(nk),dimension(nx,ny)::e,s

    do ii =  1 , ny
       s(1,ii)=-c(1,ii)/b(1,ii)
       e(1,ii)= d(1,ii)/b(1,ii)

       do i=2,nx
          s(i,ii)=-(c(i,ii)/(a(i,ii)*s(i-1,ii)+b(i,ii)))
          e(i,ii)=(d(i,ii)-a(i,ii)*e(i-1,ii))/(a(i,ii)*s(i-1,ii)+b(i,ii))
       end do

       x(nx,ii)= e(nx,ii)
       do j= nx-1, 1, -1
          x(j,ii) = s(j,ii)*x(j+1,ii)+ e(j ,ii)
       end do
    end do
  end subroutine thomas2dx


  subroutine thomas2dy (a,b,c,d,nx,ny,x)
    implicit none
    integer:: i,ii,j,nx,ny
    real(nk),dimension(nx,ny),intent(in):: a,b,c,d
    real(nk),dimension(nx,ny),intent(out):: x
    real(nk),dimension(nx,ny)::e,s

    do ii =  1 , nx
       s(ii,1) = -c(ii,1)/b(ii,1)
       e(ii,1) = d(ii,1)/b(ii,1)

       do i=2,ny
          s(ii,i)=-(c(ii,i)/(a(ii,i)*s(ii,i-1)+b(ii,i)))
          e(ii,i)=(d(ii,i)-a(ii,i)*e(ii,i-1))/(a(ii,i)*s(ii,i-1)+b(ii,i))
       end do

       x(ii,ny)= e(ii,ny)
       do j= ny-1, 1, -1
          x(ii,j) = s(ii,j)*x(ii,j+1)+ e(ii,j)
       end do
    end do

  end subroutine thomas2dy




  subroutine thomas(a,b,c,d,N,x)
    implicit none

    integer :: i,N
    real(nk), dimension(N), intent(in)  :: a,b,c,d
    real(nk), dimension(N), intent(out) :: x
    real(nk), dimension(N) :: s,e

    s(1)=c(1)/b(1) 
    e(1)=d(1)/b(1) 

    do i=2,N
       s(i)= c(i)                /(( b(i) - a(i)*s(i-1)))!+tiny( b(i) - a(i)*s(i-1)))
       e(i)= (d(i) - a(i)*e(i-1))/(( b(i) - a(i)*s(i-1)))!+tiny( b(i) - a(i)*s(i-1)))
    end do

    x(N)=e(N)

    do i=N-1,1,-1
       x(i)= e(i) - s(i)*x(i+1)
    end do
  end subroutine thomas







  subroutine allocate_Erreur_var( Er )
    implicit none
    type(Erreur_var):: Er
    do i=1,p
       allocate( Er%var(i)%ord_espace(nbr_it_espace))
       Er%var(i)%ord_espace(nbr_it_espace) = 0._nk

       allocate( Er%var(i)%ord_temps(nbr_it_temps))
       Er%var(i)%ord_temps(nbr_it_temps) = 0._nk


       allocate( Er%var(i)%val_var_espace(nbr_it_espace))
       Er%var(i)%val_var_espace(nbr_it_espace) = 0._nk

       allocate( Er%var(i)%val_var_temps(nbr_it_temps))
       Er%var(i)%val_var_temps(nbr_it_temps) = 0._nk
    end do
  end subroutine allocate_Erreur_var







  subroutine deallocate_Erreur_var( Er )
    implicit none
    type(Erreur_var):: Er
    do i=1,p
       deallocate( Er%var(i)%ord_espace)
       deallocate( Er%var(i)%ord_temps)
       deallocate( Er%var(i)%val_var_espace)
       deallocate( Er%var(i)%val_var_temps)
    end do
  end subroutine deallocate_Erreur_var



  subroutine allocate_2
    implicit none
    integer:: i

#if (DIMENSION_GEO == 2)
    do i = 1, 2
      allocate( A(i)%L           (1:N1,1:N2,1:5,1:5) )
      allocate( A(i)%R           (1:N1,1:N2,1:5,1:5) )
      allocate( A(i)%Lambda      (1:N1,1:N2,1:5) )
      allocate( A(i)%NL          (1:N1,1:N2,1:5) )
      allocate( A(i)%NL0         (1:N1,1:N2,1:5) )
      allocate( A(i)%Lambda_dRWdx(1:N1,1:N2,1:5) )
    end do
#elif (DIMENSION_GEO == 3)
    do i = 1, 3
      allocate( A(i)%L           (1:N1,1:N2,1:N3,1:9,1:9) )
      allocate( A(i)%R           (1:N1,1:N2,1:N3,1:9,1:9) )
      allocate( A(i)%Lambda      (1:N1,1:N2,1:N3,1:9) )
      allocate( A(i)%NL          (1:N1,1:N2,1:N3,1:9) )
      allocate( A(i)%NL0         (1:N1,1:N2,1:N3,1:9) )
      allocate( A(i)%Lambda_dRWdx(1:N1,1:N2,1:N3,1:9) )
    end do
#endif
    
#if (DIMENSION_GEO == 2) 
     allocate( var_intf_u1   (1:N1,1:N2-1), var_intf_u2   (1:N1-1,1:N2) )
     allocate( var_intf_u1_10(1:N1,1:N2-1), var_intf_u2_10(1:N1-1,1:N2) )
     allocate( Cent_or_HOUC_A(1:N1,1:N2) )
     allocate( Cent_or_HOUC_B(1:N1,1:N2) )

     allocate(  E_global       (0:N1+1,0:N2+1) )
     allocate(  D_diffusion    (1:N1  ,1:N2  ) )
     allocate(  F_thermal      (1:N1  ,1:N2  ) )
     allocate(  V_dissipation  (1:N1  ,1:N2  ) )
     allocate(  G_exchange     (1:N1  ,1:N2  ) )
#elif (DIMENSION_GEO == 3)
     allocate(var_intf_u1   (1:N1,1:N2-1,1:N3-1), var_intf_u2   (1:N1-1,1:N2,1:N3-1), var_intf_u3   (1:N1-1,1:N2-1,1:N3))
     allocate(var_intf_u1_10(1:N1,1:N2-1,1:N3-1), var_intf_u2_10(1:N1-1,1:N2,1:N3-1), var_intf_u3_10(1:N1-1,1:N2-1,1:N3))

     allocate( Cent_or_HOUC_A (1:N1,1:N2,1:N3) )
     allocate( Cent_or_HOUC_B (1:N1,1:N2,1:N3) )
     allocate( Cent_or_HOUC_C (1:N1,1:N2,1:N3) )

     allocate(  E_global       (0:N1+1,0:N2+1,0:N3+1) )
     allocate(  D_diffusion    (1:N1  ,1:N2  ,1:N3  ) )
     allocate(  F_thermal      (1:N1  ,1:N2  ,1:N3    ) )
     allocate(  V_dissipation  (1:N1  ,1:N2  ,1:N3    ) )
     allocate(  G_exchange     (1:N1  ,1:N2  ,1:N3    ) )
#endif
  end subroutine allocate_2

  subroutine deallocate_2
    implicit none
    integer:: i
    
#if (DIMENSION_GEO == 2)
    do i =1, 2
    deallocate( A(i)%L            )
    deallocate( A(i)%R            )
    deallocate( A(i)%Lambda       )
    deallocate( A(i)%NL           )
    deallocate( A(i)%NL0          )
    deallocate( A(i)%Lambda_dRWdx )
    end do
#elif (DIMENSION_GEO == 3)
    do i = 1, 3
    deallocate( A(i)%L            )
    deallocate( A(i)%R            )
    deallocate( A(i)%Lambda       )
    deallocate( A(i)%NL           )
    deallocate( A(i)%NL0          )
    deallocate( A(i)%Lambda_dRWdx )
    end do
#endif

#if (DIMENSION_GEO == 2) 
    deallocate(var_intf_u1, var_intf_u2)
    deallocate(var_intf_u1_10, var_intf_u2_10)
    deallocate(Cent_or_HOUC_A, Cent_or_HOUC_B )

    deallocate(  E_global       )
    deallocate(  D_diffusion    )
    deallocate(  F_thermal      )
    deallocate(  V_dissipation  )
    deallocate(  G_exchange     )
#elif (DIMENSION_GEO == 3)
    deallocate(var_intf_u1, var_intf_u2, var_intf_u3)
    deallocate(var_intf_u1_10, var_intf_u2_10, var_intf_u3_10)
    deallocate(Cent_or_HOUC_A, Cent_or_HOUC_B, Cent_or_HOUC_C )
#endif

  end subroutine deallocate_2

  subroutine discrete_coef
    implicit none
    real(nk) :: F, G, H
    
#if (DIMENSION_GEO == 2) 
     allocate( UPW2_1_A (1:N1,1:N2) ); allocate( UPW2_2_A (1:N1,1:N2) )
     allocate( UPW2_1_B (1:N1,1:N2) ); allocate( UPW2_2_B (1:N1,1:N2) )
     allocate( UPW2_1_C (1:N1,1:N2) ); allocate( UPW2_2_C (1:N1,1:N2) )

     allocate( BACKW2_1_A (1:N1,1:N2) ); allocate( BACKW2_2_A (1:N1,1:N2) )
     allocate( BACKW2_1_B (1:N1,1:N2) ); allocate( BACKW2_2_B (1:N1,1:N2) )
     allocate( BACKW2_1_C (1:N1,1:N2) ); allocate( BACKW2_2_C (1:N1,1:N2) )

     allocate( CEN2_1_A (1:N1,1:N2) ); allocate( CEN2_2_A (1:N1,1:N2) )
     allocate( CEN2_1_B (1:N1,1:N2) ); allocate( CEN2_2_B (1:N1,1:N2) )
     allocate( CEN2_1_C (1:N1,1:N2) ); allocate( CEN2_2_C (1:N1,1:N2) )

     allocate( CEN4_1_A (1:N1,1:N2) ); allocate( CEN4_2_A (1:N1,1:N2) )
     allocate( CEN4_1_B (1:N1,1:N2) ); allocate( CEN4_2_B (1:N1,1:N2) )
     allocate( CEN4_1_C (1:N1,1:N2) ); allocate( CEN4_2_C (1:N1,1:N2) )
     allocate( CEN4_1_D (1:N1,1:N2) ); allocate( CEN4_2_D (1:N1,1:N2) )
     allocate( CEN4_1_E (1:N1,1:N2) ); allocate( CEN4_2_E (1:N1,1:N2) )

     do j = 1, N2 
      do i = 1, N1
!!!!!!!!!!!!   UPWIND 2
        UPW2_1_A(i,j) = dx1_L(i,j) / ( dx1_L(i-1,j)*(dx1_L(i,j) + dx1_L(i-1,j)) )
        UPW2_1_B(i,j) = dx1_L(i,j) / ( dx1_L(i-1,j)*(dx1_L(i,j) + dx1_L(i-1,j)) ) &
                      + dx1_L(i-1,j)/( dx1_L(i  ,j)*(dx1_L(i,j) + dx1_L(i-1,j)) ) &
                      + 2. / (dx1_L(i,j) + dx1_L(i-1,j)) 
        UPW2_1_C(i,j) = dx1_L(i-1,j)/( dx1_L(i  ,j)*(dx1_L(i,j) + dx1_L(i-1,j)) ) &
                      + 2. / (dx1_L(i,j) + dx1_L(i-1,j))

        UPW2_2_A(i,j) = dx2_L(i,j) / ( dx2_L(i,j-1)*(dx2_L(i,j) + dx2_L(i,j-1)) )
        UPW2_2_B(i,j) = dx2_L(i,j) / ( dx2_L(i,j-1)*(dx2_L(i,j) + dx2_L(i,j-1)) )  &
                      + dx2_L(i,j-1)/( dx2_L(i,j  )*(dx2_L(i,j) + dx2_L(i,j-1)) )  &
                      + 2. / (dx2_L(i,j) + dx2_L(i,j-1)) 
        UPW2_2_C(i,j) = dx2_L(i,j-1)/( dx2_L(i,j  )*(dx2_L(i,j) + dx2_L(i,j-1)) )  &
                      + 2. / (dx2_L(i,j) + dx2_L(i,j-1))

!!!!!!!!!!!!   BACKWIND 2
        BACKW2_1_A(i,j) = dx1_R(i,j) / ( dx1_R(i+1,j)*(dx1_R(i,j) + dx1_R(i+1,j)) )
        BACKW2_1_B(i,j) = dx1_R(i,j) / ( dx1_R(i+1,j)*(dx1_R(i,j) + dx1_R(i+1,j)) )  &
                        + dx1_R(i+1,j)/( dx1_R(i  ,j)*(dx1_R(i,j) + dx1_R(i+1,j)) )  &
                        + 2. / (dx1_R(i,j) + dx1_R(i+1,j)) 
        BACKW2_1_C(i,j) = dx1_R(i+1,j)/( dx1_R(i  ,j)*(dx1_R(i,j) + dx1_R(i+1,j)) )  &
                        + 2. / (dx1_R(i,j) + dx1_R(i+1,j))

        BACKW2_2_A(i,j) = dx2_R(i,j) / ( dx2_R(i,j+1)*(dx2_R(i,j) + dx2_R(i,j+1)) )
        BACKW2_2_B(i,j) = dx2_R(i,j) / ( dx2_R(i,j+1)*(dx2_R(i,j) + dx2_R(i,j+1)) )  &
                        + dx2_R(i,j+1)/( dx2_R(i,j  )*(dx2_R(i,j) + dx2_R(i,j+1)) )  &
                        + 2. / (dx2_R(i,j) + dx2_R(i,j+1)) 
        BACKW2_2_C(i,j) = dx2_R(i,j+1)/( dx2_R(i,j  )*(dx2_R(i,j) + dx2_R(i,j+1)) )  &
                        + 2. / (dx2_R(i,j) + dx2_R(i,j+1))

!!!!!!!!!!!!   CENTRE 2
        CEN2_1_A(i,j) = dx1_L(i,j) / ( dx1_R(i,j)*(dx1_L(i,j) + dx1_R(i,j)) )
        CEN2_1_B(i,j) =(dx1_L(i,j) - dx1_R(i,j)) / (dx1_L(i,j) * dx1_R(i,j))
        CEN2_1_C(i,j) = dx1_R(i,j) / ( dx1_L(i,j)*(dx1_L(i,j) + dx1_R(i,j)) )

        CEN2_2_A(i,j) = dx2_L(i,j) / ( dx2_R(i,j)*(dx2_L(i,j) + dx2_R(i,j)) )
        CEN2_2_B(i,j) =(dx2_L(i,j) - dx2_R(i,j)) / (dx2_L(i,j) * dx2_R(i,j))
        CEN2_2_C(i,j) = dx2_R(i,j) / ( dx2_L(i,j)*(dx2_L(i,j) + dx2_R(i,j)) )

!!!!!!!!!!!!   CENTRE 4
      F = ((dx1_L(i-1,j)+dx1_L(i,j))**2)*((dx1_R(i,j)+dx1_R(i+1,j))**2)*(dx1_L(i-1,j)+dx1_L(i,j)+dx1_R(i,j)+dx1_R(i+1,j))
      G = (dx1_L(i,j)**2)*(dx1_R(i,j)**2)*(dx1_L(i,j)+dx1_R(i,j))
        CEN4_1_A(i,j) = G*(dx1_R(i,j)+dx1_R(i+1,j))**2
        CEN4_1_B(i,j) = F*dx1_R(i,j)**2
        CEN4_1_C(i,j) = G*( ((dx1_L(i-1,j)+dx1_L(i,j))**2) - ((dx1_R(i,j)+dx1_R(i+1,j))**2) ) - F*( dx1_L(i,j)**2 - dx1_R(i,j)**2 )
        CEN4_1_D(i,j) = F*dx1_L(i,j)**2
        CEN4_1_E(i,j) = G*((dx1_L(i-1,j)+dx1_L(i,j))**2)

        H = G*(dx1_L(i-1,j)+dx1_L(i,j))*(dx1_R(i,j)+dx1_R(i+1,j))*(dx1_L(i-1,j)+dx1_L(i,j)+dx1_R(i,j)+dx1_R(i+1,j))   &
        - F*dx1_L(i,j)*dx1_R(i,j)*(dx1_L(i,j)+dx1_R(i,j))
        CEN4_1_A(i,j) = CEN4_1_A(i,j) / H ; CEN4_1_B(i,j) = CEN4_1_B(i,j) / H
        CEN4_1_C(i,j) = CEN4_1_C(i,j) / H ; CEN4_1_D(i,j) = CEN4_1_D(i,j) / H
        CEN4_1_E(i,j) = CEN4_1_E(i,j) / H 


      F = ((dx2_L(i,j-1)+dx2_L(i,j))**2)*((dx2_R(i,j)+dx2_R(i,j+1))**2)*(dx2_L(i,j-1)+dx2_L(i,j)+dx2_R(i,j)+dx2_R(i,j+1))
      G = (dx2_L(i,j)**2)*(dx2_R(i,j)**2)*(dx2_L(i,j)+dx2_R(i,j))
        CEN4_2_A(i,j) = G*(dx2_R(i,j)+dx2_R(i,j+1))**2
        CEN4_2_B(i,j) = F*dx2_R(i,j)**2
        CEN4_2_C(i,j) = G*( ((dx2_L(i,j-1)+dx2_L(i,j))**2) - ((dx2_R(i,j)+dx2_R(i,j+1))**2) ) - F*( dx2_L(i,j)**2 - dx2_R(i,j)**2 )
        CEN4_2_D(i,j) = F*dx2_L(i,j)**2
        CEN4_2_E(i,j) = G*((dx2_L(i,j-1)+dx2_L(i,j))**2)

        H = G*(dx2_L(i,j-1)+dx2_L(i,j))*(dx2_R(i,j)+dx2_R(i,j+1))*(dx2_L(i,j-1)+dx2_L(i,j)+dx2_R(i,j)+dx2_R(i,j+1))   &
        - F*dx2_L(i,j)*dx2_R(i,j)*(dx2_L(i,j)+dx2_R(i,j))
        CEN4_2_A(i,j) = CEN4_2_A(i,j) / H ; CEN4_2_B(i,j) = CEN4_2_B(i,j) / H
        CEN4_2_C(i,j) = CEN4_2_C(i,j) / H ; CEN4_2_D(i,j) = CEN4_2_D(i,j) / H
        CEN4_2_E(i,j) = CEN4_2_E(i,j) / H 
      enddo 
     enddo
#elif (DIMENSION_GEO == 3)
     allocate( UPW2_1_A (1:N1,1:N2,1:N3) ); allocate( UPW2_2_A (1:N1,1:N2,1:N3) ); allocate( UPW2_3_A (1:N1,1:N2,1:N3) )
     allocate( UPW2_1_B (1:N1,1:N2,1:N3) ); allocate( UPW2_2_B (1:N1,1:N2,1:N3) ); allocate( UPW2_3_B (1:N1,1:N2,1:N3) )
     allocate( UPW2_1_C (1:N1,1:N2,1:N3) ); allocate( UPW2_2_C (1:N1,1:N2,1:N3) ); allocate( UPW2_3_C (1:N1,1:N2,1:N3) )

     allocate( BACKW2_1_A (1:N1,1:N2,1:N3) ); allocate( BACKW2_2_A (1:N1,1:N2,1:N3) ); allocate( BACKW2_3_A (1:N1,1:N2,1:N3) )
     allocate( BACKW2_1_B (1:N1,1:N2,1:N3) ); allocate( BACKW2_2_B (1:N1,1:N2,1:N3) ); allocate( BACKW2_3_B (1:N1,1:N2,1:N3) )
     allocate( BACKW2_1_C (1:N1,1:N2,1:N3) ); allocate( BACKW2_2_C (1:N1,1:N2,1:N3) ); allocate( BACKW2_3_C (1:N1,1:N2,1:N3) )

     allocate( CEN2_1_A (1:N1,1:N2,1:N3) ); allocate( CEN2_2_A (1:N1,1:N2,1:N3) ); allocate( CEN2_3_A (1:N1,1:N2,1:N3) )
     allocate( CEN2_1_B (1:N1,1:N2,1:N3) ); allocate( CEN2_2_B (1:N1,1:N2,1:N3) ); allocate( CEN2_3_B (1:N1,1:N2,1:N3) )
     allocate( CEN2_1_C (1:N1,1:N2,1:N3) ); allocate( CEN2_2_C (1:N1,1:N2,1:N3) ); allocate( CEN2_3_C (1:N1,1:N2,1:N3) )

     allocate( CEN4_1_A (1:N1,1:N2,1:N3) ); allocate( CEN4_2_A (1:N1,1:N2,1:N3) ); allocate( CEN4_3_A (1:N1,1:N2,1:N3) )
     allocate( CEN4_1_B (1:N1,1:N2,1:N3) ); allocate( CEN4_2_B (1:N1,1:N2,1:N3) ); allocate( CEN4_3_B (1:N1,1:N2,1:N3) )
     allocate( CEN4_1_C (1:N1,1:N2,1:N3) ); allocate( CEN4_2_C (1:N1,1:N2,1:N3) ); allocate( CEN4_3_C (1:N1,1:N2,1:N3) )
     allocate( CEN4_1_D (1:N1,1:N2,1:N3) ); allocate( CEN4_2_D (1:N1,1:N2,1:N3) ); allocate( CEN4_3_D (1:N1,1:N2,1:N3) )
     allocate( CEN4_1_E (1:N1,1:N2,1:N3) ); allocate( CEN4_2_E (1:N1,1:N2,1:N3) ); allocate( CEN4_3_E (1:N1,1:N2,1:N3) )

    do k = 1, N3
      do j = 1, N2 
        do i = 1, N1
!!!!!!!!!!!!   UPWIND 2
      UPW2_1_A(i,j,k) = dx1_L(i,j,k) / ( dx1_L(i-1,j,k)*(dx1_L(i,j,k) + dx1_L(i-1,j,k)) )
      UPW2_1_B(i,j,k) = dx1_L(i,j,k) / ( dx1_L(i-1,j,k)*(dx1_L(i,j,k) + dx1_L(i-1,j,k)) ) &
                      + dx1_L(i-1,j,k)/( dx1_L(i  ,j,k)*(dx1_L(i,j,k) + dx1_L(i-1,j,k)) ) &
                      + 2. / (dx1_L(i,j,k) + dx1_L(i-1,j,k)) 
      UPW2_1_C(i,j,k) = dx1_L(i-1,j,k)/( dx1_L(i  ,j,k)*(dx1_L(i,j,k) + dx1_L(i-1,j,k)) ) &
                      + 2. / (dx1_L(i,j,k) + dx1_L(i-1,j,k))
  
      UPW2_2_A(i,j,k) = dx2_L(i,j,k) / ( dx2_L(i,j-1,k)*(dx2_L(i,j,k) + dx2_L(i,j-1,k)) )
      UPW2_2_B(i,j,k) = dx2_L(i,j,k) / ( dx2_L(i,j-1,k)*(dx2_L(i,j,k) + dx2_L(i,j-1,k)) )  &
                      + dx2_L(i,j-1,k)/( dx2_L(i,j  ,k)*(dx2_L(i,j,k) + dx2_L(i,j-1,k)) )  &
                      + 2. / (dx2_L(i,j,k) + dx2_L(i,j-1,k)) 
      UPW2_2_C(i,j,k) = dx2_L(i,j-1,k)/( dx2_L(i,j  ,k)*(dx2_L(i,j,k) + dx2_L(i,j-1,k)) )  &
                      + 2. / (dx2_L(i,j,k) + dx2_L(i,j-1,k))

      UPW2_3_A(i,j,k) = dx3_L(i,j,k) / ( dx3_L(i,j,k-1)*(dx3_L(i,j,k) + dx3_L(i,j,k-1)) )
      UPW2_3_B(i,j,k) = dx3_L(i,j,k) / ( dx3_L(i,j,k-1)*(dx3_L(i,j,k) + dx3_L(i,j,k-1)) ) &
                      + dx3_L(i,j,k-1) / ( dx3_L(i,j,k)*(dx3_L(i,j,k) + dx3_L(i,j,k-1)) ) &
                      + 2. / (dx3_L(i,j,k) + dx3_L(i,j,k-1)) 
      UPW2_3_C(i,j,k) = dx3_L(i,j,k-1) / ( dx3_L(i,j,k)*(dx3_L(i,j,k) + dx3_L(i,j,k-1)) ) &
                      + 2. / (dx3_L(i,j,k) + dx3_L(i,j,k-1))
  
!!!!!!!!!!!!   BACKWIND 2
      BACKW2_1_A(i,j,k) = dx1_R(i,j,k) / ( dx1_R(i+1,j,k)*(dx1_R(i,j,k) + dx1_R(i+1,j,k)) )
      BACKW2_1_B(i,j,k) = dx1_R(i,j,k) / ( dx1_R(i+1,j,k)*(dx1_R(i,j,k) + dx1_R(i+1,j,k)) )  &
                        + dx1_R(i+1,j,k)/( dx1_R(i  ,j,k)*(dx1_R(i,j,k) + dx1_R(i+1,j,k)) )  &
                        + 2. / (dx1_R(i,j,k) + dx1_R(i+1,j,k)) 
      BACKW2_1_C(i,j,k) = dx1_R(i+1,j,k)/( dx1_R(i  ,j,k)*(dx1_R(i,j,k) + dx1_R(i+1,j,k)) )  &
                        + 2. / (dx1_R(i,j,k) + dx1_R(i+1,j,k))
  
      BACKW2_2_A(i,j,k) = dx2_R(i,j,k) / ( dx2_R(i,j+1,k)*(dx2_R(i,j,k) + dx2_R(i,j+1,k)) )
      BACKW2_2_B(i,j,k) = dx2_R(i,j,k) / ( dx2_R(i,j+1,k)*(dx2_R(i,j,k) + dx2_R(i,j+1,k)) )  &
                        + dx2_R(i,j+1,k)/( dx2_R(i,j  ,k)*(dx2_R(i,j,k) + dx2_R(i,j+1,k)) )  &
                        + 2. / (dx2_R(i,j,k) + dx2_R(i,j+1,k)) 
      BACKW2_2_C(i,j,k) = dx2_R(i,j+1,k)/( dx2_R(i,j  ,k)*(dx2_R(i,j,k) + dx2_R(i,j+1,k)) )  &
                        + 2. / (dx2_R(i,j,k) + dx2_R(i,j+1,k))

      BACKW2_3_A(i,j,k) = dx3_R(i,j,k) / ( dx3_R(i,j,k+1)*(dx3_R(i,j,k) + dx3_R(i,j,k+1)) )
      BACKW2_3_B(i,j,k) = dx3_R(i,j,k) / ( dx3_R(i,j,k+1)*(dx3_R(i,j,k) + dx3_R(i,j,k+1)) ) &
                        + dx3_R(i,j,k+1)/( dx3_R(i,j,k  )*(dx3_R(i,j,k) + dx3_R(i,j,k+1)) ) &
                        + 2. / (dx3_R(i,j,k) + dx3_R(i,j,k+1)) 
      BACKW2_3_C(i,j,k) = dx3_R(i,j,k+1)/( dx3_R(i,j,k  )*(dx3_R(i,j,k) + dx3_R(i,j,k+1)) ) &
                        + 2. / (dx3_R(i,j,k) + dx3_R(i,j,k+1))

!!!!!!!!!!!!   CNETRE 2
      CEN2_1_A(i,j,k) = dx1_L(i,j,k) / ( dx1_R(i,j,k)*(dx1_L(i,j,k) + dx1_R(i,j,k)) )
      CEN2_1_B(i,j,k) =(dx1_L(i,j,k) - dx1_R(i,j,k)) / (dx1_L(i,j,k) * dx1_R(i,j,k))
      CEN2_1_C(i,j,k) = dx1_R(i,j,k) / ( dx1_L(i,j,k)*(dx1_L(i,j,k) + dx1_R(i,j,k)) )

      CEN2_2_A(i,j,k) = dx2_L(i,j,k) / ( dx2_R(i,j,k)*(dx2_L(i,j,k) + dx2_R(i,j,k)) )
      CEN2_2_B(i,j,k) =(dx2_L(i,j,k) - dx2_R(i,j,k)) / (dx2_L(i,j,k) * dx2_R(i,j,k))
      CEN2_2_C(i,j,k) = dx2_R(i,j,k) / ( dx2_L(i,j,k)*(dx2_L(i,j,k) + dx2_R(i,j,k)) )

      CEN2_3_A(i,j,k) = dx3_L(i,j,k) / ( dx3_R(i,j,k)*(dx3_L(i,j,k) + dx3_R(i,j,k)) )
      CEN2_3_B(i,j,k) =(dx3_L(i,j,k) - dx3_R(i,j,k)) / (dx3_L(i,j,k) * dx3_R(i,j,k))
      CEN2_3_C(i,j,k) = dx3_R(i,j,k) / ( dx3_L(i,j,k)*(dx3_L(i,j,k) + dx3_R(i,j,k)) )

!!!!!!!!!!!!   CENTRE 4
      F = ((dx1_L(i-1,j,k)+dx1_L(i,j,k))**2)*((dx1_R(i,j,k)+dx1_R(i+1,j,k))**2) &
                                            *(dx1_L(i-1,j,k)+dx1_L(i,j,k)+dx1_R(i,j,k)+dx1_R(i+1,j,k))
      G = (dx1_L(i,j,k)**2)*(dx1_R(i,j,k)**2)*(dx1_L(i,j,k)+dx1_R(i,j,k))
      CEN4_1_A(i,j,k) = G*(dx1_R(i,j,k)+dx1_R(i+1,j,k))**2
      CEN4_1_B(i,j,k) = F*dx1_R(i,j,k)**2
      CEN4_1_C(i,j,k) = G*( ((dx1_L(i-1,j,k)+dx1_L(i,j,k))**2) - ((dx1_R(i,j,k)+dx1_R(i+1,j,k))**2) ) &
                      - F*( dx1_L(i,j,k)**2 - dx1_R(i,j,k)**2 )
      CEN4_1_D(i,j,k) = F*dx1_L(i,j,k)**2
      CEN4_1_E(i,j,k) = G*((dx1_L(i-1,j,k)+dx1_L(i,j,k))**2)
    
      H = G*(dx1_L(i-1,j,k)+dx1_L(i,j,k))*(dx1_R(i,j,k)+dx1_R(i+1,j,k))*(dx1_L(i-1,j,k)+dx1_L(i,j,k)+dx1_R(i,j,k)+dx1_R(i+1,j,k)) &
        - F*dx1_L(i,j,k)*dx1_R(i,j,k)*(dx1_L(i,j,k)+dx1_R(i,j,k))
      CEN4_1_A(i,j,k) = CEN4_1_A(i,j,k) / H ; CEN4_1_B(i,j,k) = CEN4_1_B(i,j,k) / H
      CEN4_1_C(i,j,k) = CEN4_1_C(i,j,k) / H ; CEN4_1_D(i,j,k) = CEN4_1_D(i,j,k) / H
      CEN4_1_E(i,j,k) = CEN4_1_E(i,j,k) / H 
    
    
      F = ((dx2_L(i,j-1,k)+dx2_L(i,j,k))**2)*((dx2_R(i,j,k)+dx2_R(i,j+1,k))**2)&
                                            *(dx2_L(i,j-1,k)+dx2_L(i,j,k)+dx2_R(i,j,k)+dx2_R(i,j+1,k))
      G = (dx2_L(i,j,k)**2)*(dx2_R(i,j,k)**2)*(dx2_L(i,j,k)+dx2_R(i,j,k))
      CEN4_2_A(i,j,k) = G*(dx2_R(i,j,k)+dx2_R(i,j+1,k))**2
      CEN4_2_B(i,j,k) = F*dx2_R(i,j,k)**2
      CEN4_2_C(i,j,k) = G*( ((dx2_L(i,j-1,k)+dx2_L(i,j,k))**2) - ((dx2_R(i,j,k)+dx2_R(i,j+1,k))**2) ) &
                      - F*( dx2_L(i,j,k)**2 - dx2_R(i,j,k)**2 )
      CEN4_2_D(i,j,k) = F*dx2_L(i,j,k)**2
      CEN4_2_E(i,j,k) = G*((dx2_L(i,j-1,k)+dx2_L(i,j,k))**2)
    
      H = G*(dx2_L(i,j-1,k)+dx2_L(i,j,k))*(dx2_R(i,j,k)+dx2_R(i,j+1,k))*(dx2_L(i,j-1,k)+dx2_L(i,j,k)+dx2_R(i,j,k)+dx2_R(i,j+1,k)) &
        - F*dx2_L(i,j,k)*dx2_R(i,j,k)*(dx2_L(i,j,k)+dx2_R(i,j,k))
      CEN4_2_A(i,j,k) = CEN4_2_A(i,j,k) / H ; CEN4_2_B(i,j,k) = CEN4_2_B(i,j,k) / H
      CEN4_2_C(i,j,k) = CEN4_2_C(i,j,k) / H ; CEN4_2_D(i,j,k) = CEN4_2_D(i,j,k) / H
      CEN4_2_E(i,j,k) = CEN4_2_E(i,j,k) / H 

      F = ((dx3_L(i,j,k-1)+dx3_L(i,j,k))**2)*((dx3_R(i,j,k)+dx3_R(i,j,k+1))**2)&
                                            *(dx3_L(i,j,k-1)+dx3_L(i,j,k)+dx3_R(i,j,k)+dx3_R(i,j,k+1))
      G = (dx3_L(i,j,k)**2)*(dx3_R(i,j,k)**2)*(dx3_L(i,j,k)+dx3_R(i,j,k))
      CEN4_2_A(i,j,k) = G*(dx3_R(i,j,k)+dx3_R(i,j,k+1))**2
      CEN4_2_B(i,j,k) = F*dx3_R(i,j,k)**2
      CEN4_2_C(i,j,k) = G*( ((dx3_L(i,j,k-1)+dx3_L(i,j,k))**2) - ((dx3_R(i,j,k)+dx3_R(i,j,k+1))**2) ) &
                      - F*( dx3_L(i,j,k)**2 - dx3_R(i,j,k)**2 )
      CEN4_2_D(i,j,k) = F*dx3_L(i,j,k)**2
      CEN4_2_E(i,j,k) = G*((dx3_L(i,j-1,k)+dx3_L(i,j,k))**2)
    
      H = G*(dx3_L(i,j,k-1)+dx3_L(i,j,k))*(dx3_R(i,j,k)+dx3_R(i,j,k+1))*(dx3_L(i,j,k-1)+dx3_L(i,j,k)+dx3_R(i,j,k)+dx3_R(i,j,k+1)) &
        - F*dx3_L(i,j,k)*dx3_R(i,j,k)*(dx3_L(i,j,k)+dx3_R(i,j,k))
      CEN4_3_A(i,j,k) = CEN4_3_A(i,j,k) / H ; CEN4_2_B(i,j,k) = CEN4_3_B(i,j,k) / H
      CEN4_3_C(i,j,k) = CEN4_3_C(i,j,k) / H ; CEN4_2_D(i,j,k) = CEN4_3_D(i,j,k) / H
      CEN4_3_E(i,j,k) = CEN4_3_E(i,j,k) / H 

        enddo
      enddo
    enddo

#endif
  end subroutine discrete_coef

  subroutine dediscrete_coef
    implicit none
    
#if (DIMENSION_GEO == 2) 
     deallocate( UPW2_1_A ); deallocate( UPW2_2_A )
     deallocate( UPW2_1_B ); deallocate( UPW2_2_B )
     deallocate( UPW2_1_C ); deallocate( UPW2_2_C )

     deallocate( BACKW2_1_A ); deallocate( BACKW2_2_A )
     deallocate( BACKW2_1_B ); deallocate( BACKW2_2_B )
     deallocate( BACKW2_1_C ); deallocate( BACKW2_2_C )

     deallocate( CEN2_1_A ); deallocate( CEN2_2_A )
     deallocate( CEN2_1_B ); deallocate( CEN2_2_B )
     deallocate( CEN2_1_C ); deallocate( CEN2_2_C )

     deallocate( CEN4_1_A ); deallocate( CEN4_2_A )
     deallocate( CEN4_1_B ); deallocate( CEN4_2_B )
     deallocate( CEN4_1_C ); deallocate( CEN4_2_C )
     deallocate( CEN4_1_D ); deallocate( CEN4_2_D )
     deallocate( CEN4_1_E ); deallocate( CEN4_2_E )
#elif (DIMENSION_GEO == 3)
     deallocate( UPW2_1_A ); deallocate( UPW2_2_A ); deallocate( UPW2_3_A )
     deallocate( UPW2_1_B ); deallocate( UPW2_2_B ); deallocate( UPW2_3_B )
     deallocate( UPW2_1_C ); deallocate( UPW2_2_C ); deallocate( UPW2_3_C )

     deallocate( BACKW2_1_A ); deallocate( BACKW2_2_A ); deallocate( BACKW2_3_A )
     deallocate( BACKW2_1_B ); deallocate( BACKW2_2_B ); deallocate( BACKW2_3_B )
     deallocate( BACKW2_1_C ); deallocate( BACKW2_2_C ); deallocate( BACKW2_3_C )

     deallocate( CEN2_1_A ); deallocate( CEN2_2_A ); deallocate( CEN2_3_A )
     deallocate( CEN2_1_B ); deallocate( CEN2_2_B ); deallocate( CEN2_3_B )
     deallocate( CEN2_1_C ); deallocate( CEN2_2_C ); deallocate( CEN2_3_C )

     deallocate( CEN4_1_A ); deallocate( CEN4_2_A ); deallocate( CEN4_3_A )
     deallocate( CEN4_1_B ); deallocate( CEN4_2_B ); deallocate( CEN4_3_B )
     deallocate( CEN4_1_C ); deallocate( CEN4_2_C ); deallocate( CEN4_3_C )
     deallocate( CEN4_1_D ); deallocate( CEN4_2_D ); deallocate( CEN4_3_D )
     deallocate( CEN4_1_E ); deallocate( CEN4_2_E ); deallocate( CEN4_3_E )
#endif
  end subroutine dediscrete_coef

  subroutine  Lecteur_donnees
    implicit none
    character(len=100):: work
    open(11,file="donnees",status="old")
    !open(11,file="test_boucle",status="old")
    
    read(11,*) 
    read(11,*) 
    read(11,*) work, N1 ;  read(11,*) work, N2
#if (DIMENSION_GEO == 3) 
    read(11,*) work, N3
#else
    read(11,*) 
#endif 
    read(11,*) 
    read(11,*) work, L1 ;  read(11,*) work, L2
#if (DIMENSION_GEO == 3) 
    read(11,*) work, L3
#else
    read(11,*) 
#endif 
    read(11,*); read(11,*); read(11,*); read(11,*) 
    read(11,*) work, nbr_it_espace
    read(11,*) work, nbr_it_temps
    read(11,*) work, pas_choisi_temps
    read(11,*) work, print_erreur_temps
    read(11,*) work, sonde
    read(11,*) work, dt
    read(11,*); read(11,*)
    read(11,*) work, Re_zero
    read(11,*); read(11,*); read(11,*); read(11,*) 
    read(11,*) work, Prandtl
    read(11,*) work, Reynolds
    read(11,*) work, Grashof
    read(11,*) work, Rayleigh
    read(11,*) ;     Peclet   = Prandtl * Reynolds 
#if ( VISCOELASTIC_MODELS > 0)
    read(11,*) work, Beta
    read(11,*) work, Weissenberg
    read(11,*) work, epsilon_PPT
    read(11,*) work, xi_PPT
    read(11,*) work, alpha_G
#else
    read(11,*);  read(11,*);  read(11,*);  read(11,*);  read(11,*)
#endif 
    read(11,*); read(11,*); read(11,*); read(11,*) 
    read(11,*) work, CFL
    read(11,*) work, pas_temps_adptatif
    read(11,*) work, pas_write0
    read(11,*) work, Total_step
    close(11) 


    N1 = N1+1
    N2 = N2+1
#if (DIMENSION_GEO == 3) 
    N3 = N3+1
#endif 

  end subroutine Lecteur_donnees
  
  

end module mBase

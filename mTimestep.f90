#include "definitions.h"

module mTimestep

  use mBase
  use mPression
  use mVitesse
  use mSolver
  use mOutils
  use mChampsprof
  use mHyperbolic_part
  use mViscoelastic_models


  !Declarations
  implicit none

contains
  
  
  subroutine timestep (u1, u2, P, T, t11, t12, t22, dt, it, u3, t13, t23, t33)
    implicit none
    type(variable):: u1, u2, P, T, t11, t12, t22
    type(variable),optional :: u3, t13, t23, t33
    real(nk), intent(in) :: dt
    integer, intent(in)  :: it


#if (SCHEMA_TEMPS ==1)
        call timestep_Euler1_semi_implicite(u1, u2, P, T, t11, t12, t22, dt, it, u3, t13, t23, t33)
#elif (SCHEMA_TEMPS ==2)
    if (it == 1) then 
#if (CONTINUE_CAL == 1) 
        call timestep_BDF2_Goda_semi_implicite(u1, u2, P, T, t11, t12, t22, dt, it, u3, t13, t23, t33)
#else 
        call timestep_Euler1_semi_implicite(u1, u2, P, T, t11, t12, t22, dt, it, u3, t13, t23, t33)
#endif 
    else 
        call timestep_BDF2_Goda_semi_implicite(u1, u2, P, T, t11, t12, t22, dt, it, u3, t13, t23, t33)
    end if
#endif

  end subroutine timestep
  





  subroutine timestep_Euler1_semi_implicite (u1, u2, P, T, t11, t12, t22, dt, it, u3, t13, t23, t33)
    implicit none
    type(variable):: u1, u2, P, T, t11, t12, t22
    type(variable),optional :: u3, t13, t23, t33
    real(nk), intent(in) :: dt
    integer, intent(in)  :: it
    logical, save :: ittt = .true.
    
    !calcul des termes non linéaire d'un écoulement viscoélastique 
#if(VISCOELASTIC_MODELS > 0)
#if ( DIMENSION_GEO == 2 )
#if(EXP_TREATMENT == 1)
    call calcaul_termes_NL_partie_hyperbolique  (u1%var0, u2%var0, t11%var0, t12%var0, t22%var0, &
                                                 A(1)%NL, A(2)%NL, dt, it)
#endif
    call calcaul_termes_NL_partie_hyperbolique  (u1%var, u2%var, t11%var, t12%var, t22%var, &
                                                 A(1)%NL, A(2)%NL, dt, it)
#elif( VISCOELASTIC_MODELS == 3 )
    call calcaul_termes_NL_partie_hyperbolique  (u1%var, u2%var, t11%var, t12%var, t22%var, A(1)%NL, &
                                                 A(2)%NL, dt, it, u3%var, t13%var, t23%var, t33%var, A(3)%NL)
#endif
#endif

    !calcul des contraintes viscoélastiques
#if(DIMENSION_GEO == 2 && VISCOELASTIC_MODELS > 0)
    
    t11%d = t11%var
    t12%d = t12%var
    t22%d = t22%var
    
    call coefs_NL_partie_hyperbolique(t11%d, dt, 3, A(1)%NL, A(2)%NL)
    call coefs_NL_partie_hyperbolique(t12%d, dt, 4, A(1)%NL, A(2)%NL)
    call coefs_NL_partie_hyperbolique(t22%d, dt, 5, A(1)%NL, A(2)%NL)
    
#if( VISCOELASTIC_MODELS == 2 )
    call coefs_NL_Giesekus_tau_tau(t11%d, dt, t11%var, t12%var, t22%var, 3)
    call coefs_NL_Giesekus_tau_tau(t12%d, dt, t11%var, t12%var, t22%var, 4)
    call coefs_NL_Giesekus_tau_tau(t22%d, dt, t11%var, t12%var, t22%var, 5)  
#endif 

#if( VISCOELASTIC_MODELS == 3 )
    call coefs_NL_PPT_tr_tau_tau(t11%d, dt, t11%var, t12%var, t22%var, 3)
    call coefs_NL_PPT_tr_tau_tau(t12%d, dt, t11%var, t12%var, t22%var, 4)
    call coefs_NL_PPT_tr_tau_tau(t22%d, dt, t11%var, t12%var, t22%var, 5)

    call coefs_NL_PPT_Xi(t11%d, dt, u1%var, u2%var, t11%var, t12%var, t22%var, 3)
    call coefs_NL_PPT_Xi(t12%d, dt, u1%var, u2%var, t11%var, t12%var, t22%var, 4)
    call coefs_NL_PPT_Xi(t22%d, dt, u1%var, u2%var, t11%var, t12%var, t22%var, 5)
#endif 


    call calcul_contraintes_viscoelatique( t11%var, t11%d, dt*t11%nomb_ad, 1._nk)
    call calcul_contraintes_viscoelatique( t12%var, t12%d, dt*t12%nomb_ad, 1._nk)
    call calcul_contraintes_viscoelatique( t22%var, t22%d, dt*t22%nomb_ad, 1._nk)

    call niveau_de_convergence(t11)
    call niveau_de_convergence(t12)
    call niveau_de_convergence(t22)
    call Tau_points_ficts

#elif(DIMENSION_GEO == 3 && VISCOELASTIC_MODELS > 0)

    t11%d = t11%var 
    t12%d = t12%var
    t13%d = t13%var        
    t22%d = t22%var    
    t23%d = t23%var
    t33%d = t33%var
    
    call coefs_NL_partie_hyperbolique(t11%d, dt, 4, A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique(t12%d, dt, 5, A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique(t13%d, dt, 6, A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique(t22%d, dt, 7, A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique(t23%d, dt, 8, A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique(t33%d, dt, 9, A(1)%NL, A(2)%NL, A(3)%NL)
    
#if( VISCOELASTIC_MODELS == 2 )
    call coefs_NL_Giesekus_tau_tau(t11%d, dt, t11%var, t12%var, t22%var, 4, t13%var, t23%var, t33%var)
    call coefs_NL_Giesekus_tau_tau(t12%d, dt, t11%var, t12%var, t22%var, 5, t13%var, t23%var, t33%var)
    call coefs_NL_Giesekus_tau_tau(t13%d, dt, t11%var, t12%var, t22%var, 6, t13%var, t23%var, t33%var)
    call coefs_NL_Giesekus_tau_tau(t22%d, dt, t11%var, t12%var, t22%var, 7, t13%var, t23%var, t33%var)
    call coefs_NL_Giesekus_tau_tau(t23%d, dt, t11%var, t12%var, t22%var, 8, t13%var, t23%var, t33%var)
    call coefs_NL_Giesekus_tau_tau(t33%d, dt, t11%var, t12%var, t22%var, 9, t13%var, t23%var, t33%var)
#endif 
#if( VISCOELASTIC_MODELS == 3 )

    call coefs_NL_PPT_tr_tau_tau(t11%d, dt, t11%var, t12%var, t22%var, 4, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t12%d, dt, t11%var, t12%var, t22%var, 5, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t13%d, dt, t11%var, t12%var, t22%var, 6, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t22%d, dt, t11%var, t12%var, t22%var, 7, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t23%d, dt, t11%var, t12%var, t22%var, 8, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t33%d, dt, t11%var, t12%var, t22%var, 9, t13%var, t23%var, t33%var)

    call coefs_NL_PPT_Xi(t11%d, dt, t11%var, t12%var, t22%var, 4, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t12%d, dt, t11%var, t12%var, t22%var, 5, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t13%d, dt, t11%var, t12%var, t22%var, 6, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t22%d, dt, t11%var, t12%var, t22%var, 7, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t23%d, dt, t11%var, t12%var, t22%var, 8, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t33%d, dt, t11%var, t12%var, t22%var, 9, t13%var, t23%var, t33%var)
#endif 
    
    call calcul_contraintes_viscoelatique(t11%var, t11%d, dt*t11%nomb_ad, 1._nk)
    call calcul_contraintes_viscoelatique(t12%var, t12%d, dt*t12%nomb_ad, 1._nk)
    call calcul_contraintes_viscoelatique(t13%var, t13%d, dt*t13%nomb_ad, 1._nk)
    call calcul_contraintes_viscoelatique(t22%var, t22%d, dt*t22%nomb_ad, 1._nk)
    call calcul_contraintes_viscoelatique(t23%var, t23%d, dt*t23%nomb_ad, 1._nk)
    call calcul_contraintes_viscoelatique(t33%var, t33%d, dt*t33%nomb_ad, 1._nk)
        
    call niveau_de_convergence(t11)
    call niveau_de_convergence(t12)
    call niveau_de_convergence(t13)
    call niveau_de_convergence(t22)
    call niveau_de_convergence(t23)
    call niveau_de_convergence(t33)   
    call Tau_points_ficts

#endif

#if ( DIMENSION_GEO == 2 && EXP_TREATMENT == 1)
    if (  mod(it, Re_zero) == 0) then

        t11%var0(:,:) = t11%var0(:,:)*exp( -Re_zero*dt/Weissenberg )
        t12%var0(:,:) = t12%var0(:,:)*exp( -Re_zero*dt/Weissenberg )
        t22%var0(:,:) = t22%var0(:,:)*exp( -Re_zero*dt/Weissenberg )

        t11%var (:,:) = t11%var (:,:)*exp( -Re_zero*dt/Weissenberg )
        t12%var (:,:) = t12%var (:,:)*exp( -Re_zero*dt/Weissenberg )
        t22%var (:,:) = t22%var (:,:)*exp( -Re_zero*dt/Weissenberg )
    endif
#endif

    if ( ittt ) then
#if (ENRG ==1)
    call coefs_UT_diffusion_implicite(T , dt*T%nomb_ad, 1._nk)
#endif
    call coefs_UT_diffusion_implicite(u1, dt*u1%nomb_ad, 1._nk)
    call coefs_UT_diffusion_implicite(u2, dt*u2%nomb_ad, 1._nk)
#if (DIMENSION_GEO ==3)
    call coefs_UT_diffusion_implicite(u3, dt*u2%nomb_ad, 1._nk)
#endif
#if (DIMENSION_GEO ==2)
    call coefs_P(P)  
    CALL EIGP((N1-1),(N1-1),VPX,VPRX,VPRX1,AP(:,:)%CH(1),AP(:,:)%CH(2),AP(:,:)%CH(3)) 
#elif (DIMENSION_GEO ==3)
    call coefs_P(P)
    CALL EIGP((N1-1),(N1-1),VPX,VPRX,VPRX1,AP(:,:,:)%CH(1),AP(:,:,:)%CH(2),AP(:,:,:)%CH(3),1) 
    CALL EIGP((N3-1),(N3-1),VPZ,VPRZ,VPRZ1,AP(:,:,:)%CH(7),AP(:,:,:)%CH(8),AP(:,:,:)%CH(9),3) 
#endif
      ittt = .false.
    end if    
    
   !avancement de la température 
#if (ENRG  == 1)
    T%d = T%var
    
#if (DIMENSION_GEO == 2)
    call coefs_UT_convection( T%d,  T%var, u1%var, u2%var, dt)  
    call solver_ADI(T)!, dt*T%nomb_ad) 
#elif ( DIMENSION_GEO == 3)
    call coefs_UT_convection( T%d,  T%var, u1%var, u2%var, dt, u3%var)
    call solver_ADI_3D(T)
#endif
    call UT_points_ficts(T)
    call niveau_de_convergence(T)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !avancement de la vitesse    
    u1%d = u1%var
    u2%d = u2%var

#if (ARCHIMEDE == 1)
    u2%d = u2%d  +  u2%rapp_ad *T%var *dt      
#endif

#if (DIMENSION_GEO == 2)
#if(VISCOELASTIC_MODELS == 0)
    call coefs_UT_convection(u1%d, u1%var, u1%var, u2%var, dt )
    call coefs_UT_convection(u2%d, u2%var, u1%var, u2%var, dt )
#else
    call coefs_NL_partie_hyperbolique(u1%d, dt, 1, A(1)%NL, A(2)%NL)
    call coefs_NL_partie_hyperbolique(u2%d, dt, 2, A(1)%NL, A(2)%NL)
    
#endif
    
#elif (DIMENSION_GEO == 3)        
    u3%d = u3%var
#if(VISCOELASTIC_MODELS == 0)
    call coefs_UT_convection( u1%d, u1%var, u1%var, u2%var, dt, u3%var )
    call coefs_UT_convection( u2%d, u2%var, u1%var, u2%var, dt, u3%var )
    call coefs_UT_convection( u3%d, u3%var, u1%var, u2%var, dt, u3%var )
#else
    call coefs_NL_partie_hyperbolique( u1%d, dt, 1, A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique( u2%d, dt, 2, A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique( u3%d, dt, 3, A(1)%NL, A(2)%NL, A(3)%NL)
#endif
#endif

    call coefs_U_gradx1_Pres(u1%d, P%var, dt)
    call coefs_U_gradx2_Pres(u2%d, P%var, dt)

#if (DIMENSION_GEO == 2)
    call solver_ADI(u1)!, dt*u1%nomb_ad)
    call solver_ADI(u2)!, dt*u2%nomb_ad)
#elif (DIMENSION_GEO == 3)
    call coefs_U_gradx3_Pres(u3%d, P%var, dt)

    call solver_ADI_3D(u1)
    call solver_ADI_3D(u2)
    call solver_ADI_3D(u3)
#endif

    call UT_points_ficts(u1)
    call UT_points_ficts(u2)
    call niveau_de_convergence(u1)
    call niveau_de_convergence(u2)
    
#if (DIMENSION_GEO == 3)  
    call UT_points_ficts(u3)
    call niveau_de_convergence(u3)
#endif

#if (DIMENSION_GEO == 2) 
    call source_P(P%d(1:N1-1,1:N2-1), u1%var, u2%var, dt)
    call RESOLP((N1-1), (N2-1), P%var0(0:N1,0:N2), P%d(0:N1,0:N2),AP(0:N1,0:N2)%CH(4), &
                AP(0:N1,0:N2)%CH(6),AP(0:N1,0:N2)%CH(5),1,(N1-1),1,(N2-1),VPX,VPRX,VPRX1)
    
#elif (DIMENSION_GEO == 3)
    call source_P(P%d(1:N1-1,1:N2-1,1:N3-1), u1%var, u2%var, dt, u3%var)
    call RESOLP_3D((N1-1), (N2-1), (N3-1), P%var0(0:N1,0:N2,0:N3), P%d(0:N1,0:N2,0:N3),       &
                 AP(0:N1,0:N2,0:N3)%CH(4), AP(0:N1,0:N2,0:N3)%CH(6),AP(0:N1,0:N2,0:N3)%CH(5), &
                 1,(N1-1),1,(N2-1),1,(N3-1),VPX,VPRX,VPRX1,VPZ,VPRZ,VPRZ1  )
#endif
    print*, 'P%var0'
    print*, P%var0
    pause

    call P_points_ficts(P,0)
    call niveau_de_convergence(P)

#if ( DIMENSION_GEO == 2)
    print*, 't11%var'
    print*, t11%var 
    pause
    print*, 't12%var'
    print*, t12%var 
    pause
    print*, 't22%var'
    print*, t22%var 
    pause
    print*, 'T%var'
    print*, T%var 
    pause
    print*, 'u1%var'
    print*, u1%var 
    pause
    print*, 'u2%var'
    print*, u2%var 
    pause
    print*, 'P%var0'
    print*, P%var0 
    pause
#elif (DIMENSION_GEO == 3)
    print*, 't11%var(:,:,0)'
    print*, t11%var(:,:,0)
    print*, 't11%var(:,:,1)'
    print*, t11%var(:,:,1)
    print*, 't11%var(:,:,2)'
    print*, t11%var(:,:,2)
    print*, 't11%var(:,:,3)'
    print*, t11%var(:,:,3)
    print*, 't11%var(:,:,4)'
    print*, t11%var(:,:,4)
    print*, 't11%var(:,:,5)'
    print*, t11%var(:,:,5)
    print*, 't11%var(:,:,6)'
    print*, t11%var(:,:,6)
    print*, 't11%var(:,:,7)'
    print*, t11%var(:,:,7)
    print*, 't11%var(:,:,8)'
    print*, t11%var(:,:,8)
    print*, 't11%var(:,:,9)'
    print*, t11%var(:,:,9)
    print*, 't11%var(:,:,10)'
    print*, t11%var(:,:,10)
    print*, 't11%var(:,:,11)'
    print*, t11%var(:,:,11)
    pause
    print*, 't12%var(:,:,0)'
    print*, t12%var(:,:,0)
    print*, 't12%var(:,:,1)'
    print*, t12%var(:,:,1)
    print*, 't12%var(:,:,2)'
    print*, t12%var(:,:,2)
    print*, 't12%var(:,:,3)'
    print*, t12%var(:,:,3)
    print*, 't12%var(:,:,4)'
    print*, t12%var(:,:,4)
    print*, 't12%var(:,:,5)'
    print*, t12%var(:,:,5)
    print*, 't12%var(:,:,6)'
    print*, t12%var(:,:,6)
    print*, 't12%var(:,:,7)'
    print*, t12%var(:,:,7)
    print*, 't12%var(:,:,8)'
    print*, t12%var(:,:,8)
    print*, 't12%var(:,:,9)'
    print*, t12%var(:,:,9)
    print*, 't12%var(:,:,10)'
    print*, t12%var(:,:,10)
    print*, 't12%var(:,:,11)'
    print*, t12%var(:,:,11)
    pause
    print*, 't22%var(:,:,0)'
    print*, t22%var(:,:,0)
    print*, 't22%var(:,:,1)'
    print*, t22%var(:,:,1)
    print*, 't22%var(:,:,2)'
    print*, t22%var(:,:,2)
    print*, 't22%var(:,:,3)'
    print*, t22%var(:,:,3)
    print*, 't22%var(:,:,4)'
    print*, t22%var(:,:,4)
    print*, 't22%var(:,:,5)'
    print*, t22%var(:,:,5)
    print*, 't22%var(:,:,6)'
    print*, t22%var(:,:,6)
    print*, 't22%var(:,:,7)'
    print*, t22%var(:,:,7)
    print*, 't22%var(:,:,8)'
    print*, t22%var(:,:,8)
    print*, 't22%var(:,:,9)'
    print*, t22%var(:,:,9)
    print*, 't22%var(:,:,10)'
    print*, t22%var(:,:,10)
    print*, 't22%var(:,:,11)'
    print*, t22%var(:,:,11)
    pause
    print*, 'T%var(:,:,0)'
    print*, T%var(:,:,0)
    print*, 'T%var(:,:,1)'
    print*, T%var(:,:,1)
    print*, 'T%var(:,:,2)'
    print*, T%var(:,:,2)
    print*, 'T%var(:,:,3)'
    print*, T%var(:,:,3)
    print*, 'T%var(:,:,4)'
    print*, T%var(:,:,4)
    print*, 'T%var(:,:,5)'
    print*, T%var(:,:,5)
    print*, 'T%var(:,:,6)'
    print*, T%var(:,:,6)
    print*, 'T%var(:,:,7)'
    print*, T%var(:,:,7)
    print*, 'T%var(:,:,8)'
    print*, T%var(:,:,8)
    print*, 'T%var(:,:,9)'
    print*, T%var(:,:,9)
    print*, 'T%var(:,:,10)'
    print*, T%var(:,:,10)
    print*, 'T%var(:,:,11)'
    print*, T%var(:,:,11)
    pause
    print*, 'u1%var(:,:,0)'
    print*, u1%var(:,:,0) 
    print*, 'u1%var(:,:,1)'
    print*, u1%var(:,:,1) 
    print*, 'u1%var(:,:,2)'
    print*, u1%var(:,:,2) 
    print*, 'u1%var(:,:,3)'
    print*, u1%var(:,:,3) 
    print*, 'u1%var(:,:,4)'
    print*, u1%var(:,:,4) 
    print*, 'u1%var(:,:,5)'
    print*, u1%var(:,:,5) 
    print*, 'u1%var(:,:,6)'
    print*, u1%var(:,:,6) 
    print*, 'u1%var(:,:,7)'
    print*, u1%var(:,:,7) 
    print*, 'u1%var(:,:,8)'
    print*, u1%var(:,:,8) 
    print*, 'u1%var(:,:,9)'
    print*, u1%var(:,:,9)  
    print*, 'u1%var(:,:,10)'
    print*, u1%var(:,:,10) 
    print*, 'u1%var(:,:,11)'
    print*, u1%var(:,:,11)  
    pause
    print*, 'u2%var(:,:,0)'
    print*, u2%var(:,:,0) 
    print*, 'u2%var(:,:,1)'
    print*, u2%var(:,:,1) 
    print*, 'u2%var(:,:,2)'
    print*, u2%var(:,:,2) 
    print*, 'u2%var(:,:,3)'
    print*, u2%var(:,:,3) 
    print*, 'u2%var(:,:,4)'
    print*, u2%var(:,:,4) 
    print*, 'u2%var(:,:,5)'
    print*, u2%var(:,:,5) 
    print*, 'u2%var(:,:,6)'
    print*, u2%var(:,:,6) 
    print*, 'u2%var(:,:,7)'
    print*, u2%var(:,:,7) 
    print*, 'u2%var(:,:,8)'
    print*, u2%var(:,:,8) 
    print*, 'u2%var(:,:,9)'
    print*, u2%var(:,:,9)
    print*, 'u2%var(:,:,10)'
    print*, u2%var(:,:,10) 
    print*, 'u2%var(:,:,11)'
    print*, u2%var(:,:,11)
    pause
    print*, 'P%var0(:,:,0)'
    print*, P%var0(:,:,0)
    print*, 'P%var0(:,:,1)'
    print*, P%var0(:,:,1)
    print*, 'P%var0(:,:,2)'
    print*, P%var0(:,:,2)
    print*, 'P%var0(:,:,3)'
    print*, P%var0(:,:,3)
    print*, 'P%var0(:,:,4)'
    print*, P%var0(:,:,4)
    print*, 'P%var0(:,:,5)'
    print*, P%var0(:,:,5)
    print*, 'P%var0(:,:,6)'
    print*, P%var0(:,:,6)
    print*, 'P%var0(:,:,7)'
    print*, P%var0(:,:,7)
    print*, 'P%var0(:,:,8)'
    print*, P%var0(:,:,8)
    print*, 'P%var0(:,:,9)'
    print*, P%var0(:,:,9)
    print*, 'P%var0(:,:,10)'
    print*, P%var0(:,:,10)
    pause
#endif
    !etape de correction
#if (DIMENSION_GEO == 2)
    call correction_Pression(u1%var, u2%var, P%var, P%var0, 1._nk, u1%nomb_ad)
    call correction_vitesse (u1%var, u2%var,        P%var0, dt)       
#elif (DIMENSION_GEO == 3)

    call correction_Pression(u1%var, u2%var, P%var, P%var0, 1._nk, u1%nomb_ad, u3%var) 
    call correction_vitesse (u1%var, u2%var,        P%var0, dt               , u3%var)
    call UT_points_ficts(u3)
#endif
    
    call  P_points_ficts(P )
    call UT_points_ficts(u1)
    call UT_points_ficts(u2)      

     end subroutine timestep_Euler1_semi_implicite
  
  



























  subroutine timestep_BDF2_Goda_semi_implicite (u1, u2, P, T, t11, t12, t22, dt, it, u3, t13, t23, t33)
    implicit none
    type(variable):: u1, u2, P, T, t11, t12, t22
    type(variable),optional :: u3, t13, t23, t33
    real(nk), intent(in) :: dt
    integer, intent(in)  :: it
    logical, save        :: ittt = .true.

    !calcul des termes non linéaire d'un écoulement viscoélastique
#if(VISCOELASTIC_MODELS > 0)
#if(DIMENSION_GEO == 2)
#if(EXP_TREATMENT == 1)
    call calcaul_termes_NL_partie_hyperbolique  (u1%var0, u2%var0, t11%var0, t12%var0, t22%var0, &
                                                 A(1)%NL, A(2)%NL, dt, it)
#endif
    do i = 1, 5  
       A(1)%NL0(:,:,i) = A(1)%NL(:,:,i)
       A(2)%NL0(:,:,i) = A(2)%NL(:,:,i)
    end do

    call calcaul_termes_NL_partie_hyperbolique  (u1%var, u2%var, t11%var, t12%var, t22%var, &
                                                 A(1)%NL, A(2)%NL, dt, it)
#elif(DIMENSION_GEO == 3) 
    do i = 1, 9
       A(1)%NL0(:,:,:,i) = A(1)%NL(:,:,:,i)
       A(2)%NL0(:,:,:,i) = A(2)%NL(:,:,:,i)
       A(3)%NL0(:,:,:,i) = A(3)%NL(:,:,:,i)
    end do

    call calcaul_termes_NL_partie_hyperbolique  (u1%var, u2%var, t11%var, t12%var, t22%var, A(1)%NL, &
                                                 A(2)%NL, dt, it, u3%var, t13%var, t23%var, t33%var, A(3)%NL)
#endif
#endif

    !calcul des contraintes viscoélastiques
#if(DIMENSION_GEO == 2 && VISCOELASTIC_MODELS > 0)
    
    t11%d = 4.*t11%var - t11%var0 
    t12%d = 4.*t12%var - t12%var0
    t22%d = 4.*t22%var - t22%var0

    call coefs_NL_partie_hyperbolique(t11%d, 4.*dt, 3, A(1)%NL, A(2)%NL)
    call coefs_NL_partie_hyperbolique(t12%d, 4.*dt, 4, A(1)%NL, A(2)%NL)
    call coefs_NL_partie_hyperbolique(t22%d, 4.*dt, 5, A(1)%NL, A(2)%NL)

    call coefs_NL_partie_hyperbolique(t11%d, -2.*dt, 3, A(1)%NL0, A(2)%NL0)
    call coefs_NL_partie_hyperbolique(t12%d, -2.*dt, 4, A(1)%NL0, A(2)%NL0)
    call coefs_NL_partie_hyperbolique(t22%d, -2.*dt, 5, A(1)%NL0, A(2)%NL0)
    
#if( VISCOELASTIC_MODELS == 2 )
    call coefs_NL_Giesekus_tau_tau(t11%d, 4.*dt, t11%var, t12%var, t22%var, 3)
    call coefs_NL_Giesekus_tau_tau(t12%d, 4.*dt, t11%var, t12%var, t22%var, 4)
    call coefs_NL_Giesekus_tau_tau(t22%d, 4.*dt, t11%var, t12%var, t22%var, 5)  

    call coefs_NL_Giesekus_tau_tau(t11%d, -2.*dt, t11%var0, t12%var0, t22%var0, 3)
    call coefs_NL_Giesekus_tau_tau(t12%d, -2.*dt, t11%var0, t12%var0, t22%var0, 4)
    call coefs_NL_Giesekus_tau_tau(t22%d, -2.*dt, t11%var0, t12%var0, t22%var0, 5)  
#endif 
#if( VISCOELASTIC_MODELS == 3 )
    call coefs_NL_PPT_tr_tau_tau(t11%d, 4.*dt, t11%var, t12%var, t22%var, 3)
    call coefs_NL_PPT_tr_tau_tau(t12%d, 4.*dt, t11%var, t12%var, t22%var, 4)
    call coefs_NL_PPT_tr_tau_tau(t22%d, 4.*dt, t11%var, t12%var, t22%var, 5)

    call coefs_NL_PPT_tr_tau_tau(t11%d, -2.*dt, t11%var0, t12%var0, t22%var0, 3)
    call coefs_NL_PPT_tr_tau_tau(t12%d, -2.*dt, t11%var0, t12%var0, t22%var0, 4)
    call coefs_NL_PPT_tr_tau_tau(t22%d, -2.*dt, t11%var0, t12%var0, t22%var0, 5)

    call coefs_NL_PPT_Xi(t11%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 3)
    call coefs_NL_PPT_Xi(t12%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 4)
    call coefs_NL_PPT_Xi(t22%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 5)

    call coefs_NL_PPT_Xi(t11%d, -2.*dt, u1%var0, u2%var0, t11%var0, t12%var0, t22%var0, 3)
    call coefs_NL_PPT_Xi(t12%d, -2.*dt, u1%var0, u2%var0, t11%var0, t12%var0, t22%var0, 4)
    call coefs_NL_PPT_Xi(t22%d, -2.*dt, u1%var0, u2%var0, t11%var0, t12%var0, t22%var0, 5)
#endif

    t11%var0 = t11%var  
    t12%var0 = t12%var
    t22%var0 = t22%var
    call calcul_contraintes_viscoelatique(t11%var, t11%d, 2.*dt*t11%nomb_ad, 3._nk)
    call calcul_contraintes_viscoelatique(t12%var, t12%d, 2.*dt*t12%nomb_ad, 3._nk)
    call calcul_contraintes_viscoelatique(t22%var, t22%d, 2.*dt*t22%nomb_ad, 3._nk)

!!!!!
    call niveau_de_convergence(t11)
    call niveau_de_convergence(t12)
    call niveau_de_convergence(t22)

#elif(DIMENSION_GEO == 3 && VISCOELASTIC_MODELS > 0)

    
    t11%d = 4.*t11%var - t11%var0; t12%d = 4.*t12%var - t12%var0
    t13%d = 4.*t13%var - t13%var0; t22%d = 4.*t22%var - t22%var0
    t23%d = 4.*t23%var - t23%var0; t33%d = 4.*t33%var - t33%var0
    
    call coefs_NL_partie_hyperbolique(t11%d, 4.*dt, 4, A(1)%NL , A(2)%NL , A(3)%NL )
    call coefs_NL_partie_hyperbolique(t12%d, 4.*dt, 5, A(1)%NL , A(2)%NL , A(3)%NL )
    call coefs_NL_partie_hyperbolique(t13%d, 4.*dt, 6, A(1)%NL , A(2)%NL , A(3)%NL )    
    call coefs_NL_partie_hyperbolique(t22%d, 4.*dt, 7, A(1)%NL , A(2)%NL , A(3)%NL )
    call coefs_NL_partie_hyperbolique(t23%d, 4.*dt, 8, A(1)%NL , A(2)%NL , A(3)%NL )
    call coefs_NL_partie_hyperbolique(t33%d, 4.*dt, 9, A(1)%NL , A(2)%NL , A(3)%NL )
    
    call coefs_NL_partie_hyperbolique(t11%d,-2.*dt, 4, A(1)%NL0, A(2)%NL0, A(3)%NL0)
    call coefs_NL_partie_hyperbolique(t12%d,-2.*dt, 5, A(1)%NL0, A(2)%NL0, A(3)%NL0)
    call coefs_NL_partie_hyperbolique(t13%d,-2.*dt, 6, A(1)%NL0, A(2)%NL0, A(3)%NL0)
    call coefs_NL_partie_hyperbolique(t22%d,-2.*dt, 7, A(1)%NL0, A(2)%NL0, A(3)%NL0)
    call coefs_NL_partie_hyperbolique(t23%d,-2.*dt, 8, A(1)%NL0, A(2)%NL0, A(3)%NL0)
    call coefs_NL_partie_hyperbolique(t33%d,-2.*dt, 9, A(1)%NL0, A(2)%NL0, A(3)%NL0)

#if( VISCOELASTIC_MODELS == 2 )
    call coefs_NL_Giesekus_tau_tau(t11%d, 4.*dt, t11%var , t12%var , t22%var , 4, t13%var , t23%var , t33%var )
    call coefs_NL_Giesekus_tau_tau(t12%d, 4.*dt, t11%var , t12%var , t22%var , 5, t13%var , t23%var , t33%var )
    call coefs_NL_Giesekus_tau_tau(t13%d, 4.*dt, t11%var , t12%var , t22%var , 6, t13%var , t23%var , t33%var )
    call coefs_NL_Giesekus_tau_tau(t22%d, 4.*dt, t11%var , t12%var , t22%var , 7, t13%var , t23%var , t33%var )
    call coefs_NL_Giesekus_tau_tau(t23%d, 4.*dt, t11%var , t12%var , t22%var , 8, t13%var , t23%var , t33%var )
    call coefs_NL_Giesekus_tau_tau(t33%d, 4.*dt, t11%var , t12%var , t22%var , 9, t13%var , t23%var , t33%var )

    call coefs_NL_Giesekus_tau_tau(t11%d,-2.*dt, t11%var0, t12%var0, t22%var0, 4, t13%var0, t23%var0, t33%var0)
    call coefs_NL_Giesekus_tau_tau(t12%d,-2.*dt, t11%var0, t12%var0, t22%var0, 5, t13%var0, t23%var0, t33%var0)
    call coefs_NL_Giesekus_tau_tau(t13%d,-2.*dt, t11%var0, t12%var0, t22%var0, 6, t13%var0, t23%var0, t33%var0)
    call coefs_NL_Giesekus_tau_tau(t22%d,-2.*dt, t11%var0, t12%var0, t22%var0, 7, t13%var0, t23%var0, t33%var0)
    call coefs_NL_Giesekus_tau_tau(t23%d,-2.*dt, t11%var0, t12%var0, t22%var0, 8, t13%var0, t23%var0, t33%var0)
    call coefs_NL_Giesekus_tau_tau(t33%d,-2.*dt, t11%var0, t12%var0, t22%var0, 9, t13%var0, t23%var0, t33%var0)
#endif 

#if( VISCOELASTIC_MODELS == 3 )
    call coefs_NL_PPT_tr_tau_tau(t11%d, 4.*dt, t11%var, t12%var, t22%var, 4, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t12%d, 4.*dt, t11%var, t12%var, t22%var, 5, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t13%d, 4.*dt, t11%var, t12%var, t22%var, 6, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t22%d, 4.*dt, t11%var, t12%var, t22%var, 7, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t23%d, 4.*dt, t11%var, t12%var, t22%var, 8, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_tr_tau_tau(t33%d, 4.*dt, t11%var, t12%var, t22%var, 9, t13%var, t23%var, t33%var)

    call coefs_NL_PPT_tr_tau_tau(t11%d,-2.*dt, t11%var0, t12%var0, t22%var0, 4, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_tr_tau_tau(t12%d,-2.*dt, t11%var0, t12%var0, t22%var0, 5, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_tr_tau_tau(t13%d,-2.*dt, t11%var0, t12%var0, t22%var0, 6, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_tr_tau_tau(t22%d,-2.*dt, t11%var0, t12%var0, t22%var0, 7, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_tr_tau_tau(t23%d,-2.*dt, t11%var0, t12%var0, t22%var0, 8, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_tr_tau_tau(t33%d,-2.*dt, t11%var0, t12%var0, t22%var0, 9, t13%var0, t23%var0, t33%var0)

    call coefs_NL_PPT_Xi(t11%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 4, u3 %var, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t12%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 5, u3 %var, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t13%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 6, u3 %var, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t22%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 7, u3 %var, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t23%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 8, u3 %var, t13%var, t23%var, t33%var)
    call coefs_NL_PPT_Xi(t33%d, 4.*dt, u1%var, u2%var, t11%var, t12%var, t22%var, 9, u3 %var, t13%var, t23%var, t33%var)

    call coefs_NL_PPT_Xi(t11%d,-2.*dt, u1 %var0, u2 %var0, t11%var0, t12%var0, t22%var0, 4, u3 %var0, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_Xi(t12%d,-2.*dt, u1 %var0, u2 %var0, t11%var0, t12%var0, t22%var0, 5, u3 %var0, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_Xi(t13%d,-2.*dt, u1 %var0, u2 %var0, t11%var0, t12%var0, t22%var0, 6, u3 %var0, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_Xi(t22%d,-2.*dt, u1 %var0, u2 %var0, t11%var0, t12%var0, t22%var0, 7, u3 %var0, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_Xi(t23%d,-2.*dt, u1 %var0, u2 %var0, t11%var0, t12%var0, t22%var0, 8, u3 %var0, t13%var0, t23%var0, t33%var0)
    call coefs_NL_PPT_Xi(t33%d,-2.*dt, u1 %var0, u2 %var0, t11%var0, t12%var0, t22%var0, 9, u3 %var0, t13%var0, t23%var0, t33%var0)
#endif

    t11%var0 = t11%var; t12%var0 = t12%var; t13%var0 = t13%var
    t22%var0 = t22%var; t23%var0 = t23%var; t33%var0 = t33%var
    call calcul_contraintes_viscoelatique(t11%var, t11%d, 2.*dt*t11%nomb_ad, 3._nk)
    call calcul_contraintes_viscoelatique(t12%var, t12%d, 2.*dt*t12%nomb_ad, 3._nk)
    call calcul_contraintes_viscoelatique(t13%var, t13%d, 2.*dt*t13%nomb_ad, 3._nk)
    call calcul_contraintes_viscoelatique(t22%var, t22%d, 2.*dt*t22%nomb_ad, 3._nk)
    call calcul_contraintes_viscoelatique(t23%var, t23%d, 2.*dt*t23%nomb_ad, 3._nk)
    call calcul_contraintes_viscoelatique(t33%var, t33%d, 2.*dt*t33%nomb_ad, 3._nk)
        
    call niveau_de_convergence(t11)
    call niveau_de_convergence(t12)
    call niveau_de_convergence(t13)
    call niveau_de_convergence(t22)
    call niveau_de_convergence(t23)
    call niveau_de_convergence(t33)

    call Tau_points_ficts

#endif
    
#if ( DIMENSION == 2 )
#if ( EXP_TREATMENT == 1)
    if (  mod(it, Re_zero) == 0) then

        t11%var0(:,:) = t11%var0(:,:)*exp( -Re_zero*dt/Weissenberg )
        t12%var0(:,:) = t12%var0(:,:)*exp( -Re_zero*dt/Weissenberg )
        t22%var0(:,:) = t22%var0(:,:)*exp( -Re_zero*dt/Weissenberg )

        t11%var (:,:) = t11%var (:,:)*exp( -Re_zero*dt/Weissenberg )
        t12%var (:,:) = t12%var (:,:)*exp( -Re_zero*dt/Weissenberg )
        t22%var (:,:) = t22%var (:,:)*exp( -Re_zero*dt/Weissenberg )
    endif
#endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( ittt ) then
#if (ENRG ==1)
        call coefs_UT_diffusion_implicite( T, 2./3.*dt*T %nomb_ad, 1._nk)
#endif
        call coefs_UT_diffusion_implicite(u1, 2./3.*dt*u1%nomb_ad, 1._nk)
        call coefs_UT_diffusion_implicite(u2, 2./3.*dt*u2%nomb_ad, 1._nk)
#if (DIMENSION_GEO ==3)
        call coefs_UT_diffusion_implicite(u3, 2./3.*dt*u3%nomb_ad, 1._nk)
#endif
#if (DIMENSION_GEO ==2)
        call coefs_P(P)  
        call calculs_CL_P(P)
        CALL EIGP((N1-1),(N1-1),VPX,VPRX,VPRX1,AP(:,:)%CH(1),AP(:,:)%CH(2),AP(:,:)%CH(3)) 
#elif (DIMENSION_GEO ==3)
        call coefs_P(P)  
        call calculs_CL_P(P)
        CALL EIGP((N1-1),(N1-1),VPX,VPRX,VPRX1,AP(:,:,:)%CH(1),AP(:,:,:)%CH(2),AP(:,:,:)%CH(3),1) 
        CALL EIGP((N3-1),(N3-1),VPZ,VPRZ,VPRZ1,AP(:,:,:)%CH(7),AP(:,:,:)%CH(8),AP(:,:,:)%CH(9),3) 
#endif
          ittt = .false.
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !avancement de la température 
#if (ENRG  == 1)
    T%d = ( 4._nk*T%var - T%var0 )/3._nk
#if (DIMENSION_GEO == 2)    
    call coefs_UT_convection( T%d,  T%var , u1%var , u2%var , 4./3.*dt)
    call coefs_UT_convection( T%d,  T%var0, u1%var0, u2%var0,-2./3.*dt)
#elif ( DIMENSION_GEO == 3)
    call coefs_UT_convection( T%d,  T%var , u1%var , u2%var ,  4./3.*dt, u3%var )
    call coefs_UT_convection( T%d,  T%var0, u1%var0, u2%var0, -2./3.*dt, u3%var0)
#endif
    T%var0 = T%var
#if (DIMENSION_GEO == 2)
    call solver_ADI(T)
#elif (DIMENSION_GEO == 3)
    call solver_ADI_3D(T)
#endif
    call UT_points_ficts(T)
    call niveau_de_convergence(T)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !avancement de la vitesse    
    u1%d =( 4. *u1%var - u1%var0 )/3._nk
    u2%d =( 4. *u2%var - u2%var0 )/3._nk
#if (ARCHIMEDE == 1)
    u2%d = u2%d  +  ( 2.*dt*u2%rapp_ad *T%var )/3._nk   
#endif
#if (DIMENSION_GEO == 2)
#if(VISCOELASTIC_MODELS == 0)
    call coefs_UT_convection(u1%d, u1%var, u1%var, u2%var, 4./3.*dt )
    call coefs_UT_convection(u2%d, u2%var, u1%var, u2%var, 4./3.*dt )

    call coefs_UT_convection(u1%d, u1%var0, u1%var0, u2%var0, -2./3.*dt )
    call coefs_UT_convection(u2%d, u2%var0, u1%var0, u2%var0, -2./3.*dt )    
#else
    call coefs_NL_partie_hyperbolique(u1%d, 4./3.*dt, 1, A(1)%NL, A(2)%NL)
    call coefs_NL_partie_hyperbolique(u2%d, 4./3.*dt, 2, A(1)%NL, A(2)%NL)
    
    call coefs_NL_partie_hyperbolique(u1%d, -2./3.*dt, 1, A(1)%NL0, A(2)%NL0)
    call coefs_NL_partie_hyperbolique(u2%d, -2./3.*dt, 2, A(1)%NL0, A(2)%NL0)
#endif
#elif (DIMENSION_GEO == 3)        
    u3%d = ( 4. *u3%var - u3%var0 )/3._nk
#if(VISCOELASTIC_MODELS == 0)
    call coefs_UT_convection(u1%d, u1%var, u1%var, u2%var, 4./3.*dt, u3%var)
    call coefs_UT_convection(u2%d, u2%var, u1%var, u2%var, 4./3.*dt, u3%var)
    call coefs_UT_convection(u3%d, u3%var, u1%var, u2%var, 4./3.*dt, u3%var)
    
    call coefs_UT_convection(u1%d, u1%var0, u1%var0, u2%var0, -2./3.*dt, u3%var0)
    call coefs_UT_convection(u2%d, u2%var0, u1%var0, u2%var0, -2./3.*dt, u3%var0)
    call coefs_UT_convection(u3%d, u3%var0, u1%var0, u2%var0, -2./3.*dt, u3%var0)
#else
    call coefs_NL_partie_hyperbolique(u1%d, 4./3.*dt, 1 , A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique(u2%d, 4./3.*dt, 2 , A(1)%NL, A(2)%NL, A(3)%NL)
    call coefs_NL_partie_hyperbolique(u3%d, 4./3.*dt, 3 , A(1)%NL, A(2)%NL, A(3)%NL)
    
    call coefs_NL_partie_hyperbolique(u1%d, -2./3.*dt, 1 , A(1)%NL0, A(2)%NL0, A(3)%NL0)
    call coefs_NL_partie_hyperbolique(u2%d, -2./3.*dt, 2 , A(1)%NL0, A(2)%NL0, A(3)%NL0)
    call coefs_NL_partie_hyperbolique(u3%d, -2./3.*dt, 3 , A(1)%NL0, A(2)%NL0, A(3)%NL0)
#endif
#endif
    call coefs_U_gradx1_Pres(u1%d, P%var, 2./3.*dt)
    call coefs_U_gradx2_Pres(u2%d, P%var, 2./3.*dt) 

    u1%var0 = u1%var
    u2%var0 = u2%var
#if (DIMENSION_GEO == 2)
    call solver_ADI(u1)
    call solver_ADI(u2)
#elif (DIMENSION_GEO == 3)
    call coefs_U_gradx3_Pres(u3%d, P%var, 2./3.*dt)

    u3%var0 = u3%var
    call solver_ADI_3D(u1)
    call solver_ADI_3D(u2)
    call solver_ADI_3D(u3)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call UT_points_ficts(u1)
    call UT_points_ficts(u2)
    call niveau_de_convergence(u1)
    call niveau_de_convergence(u2)
#if (DIMENSION_GEO == 3)  
    call UT_points_ficts(u3)
    call niveau_de_convergence(u3)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if (DIMENSION_GEO == 2)
    call source_P(P%d(1:N1-1,1:N2-1), u1%var, u2%var, dt)
    call RESOLP((N1-1), (N2-1), P%var0(:,:), P%d(:,:),AP(:,:)%CH(4), &
      AP(:,:)%CH(6),AP(:,:)%CH(5),1,(N1-1),1,(N2-1),VPX,VPRX,VPRX1)
#elif (DIMENSION_GEO == 3)

    call source_P(P%d(1:N1-1,1:N2-1,1:N3-1), u1%var, u2%var, dt, u3%var)
    call RESOLP_3D((N1-1), (N2-1), (N3-1), P%var0(:,:,:), P%d(:,:,:),AP(:,:,:)%CH(4), &
     AP(:,:,:)%CH(6),AP(:,:,:)%CH(5),1,(N1-1),1,(N2-1),1,(N3-1),VPX,VPRX,VPRX1,VPZ,VPRZ,VPRZ1  )
#endif

    call  P_points_ficts(P,0)
    call niveau_de_convergence(P)

#if (DIMENSION_GEO == 2)     
    call correction_Pression(u1%var, u2%var, P%var, P%var0, 3._nk/2._nk, u1%nomb_ad)
    call correction_vitesse (u1%var, u2%var,        P%var0, dt)            
#elif (DIMENSION_GEO == 3)
    call correction_Pression(u1%var, u2%var, P%var, P%var0, 3._nk/2._nk, u1%nomb_ad, u3%var)
    call correction_vitesse (u1%var, u2%var,        P%var0, dt                     , u3%var)

    call UT_points_ficts(u3)
#endif
    call  P_points_ficts(P )
    call UT_points_ficts(u1)
    call UT_points_ficts(u2)  

  end subroutine timestep_BDF2_Goda_semi_implicite
  
end module mTimestep

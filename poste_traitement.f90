program poste_traitement 
  implicit none

  integer :: i,ii, j,jj, k,kk, l,ll, s, nbr_delta_t 
  integer :: N, M,NM, indice, RICH, RICH_O2, RICH_O4, SOL_ANA, erreur_en_temps, erreur_en_espace
  
  
  real(nk),dimension(:),allocatable:: u_1, v_1, P_1, u_approx, v_approx, P_approx
  real(nk),dimension(:),allocatable:: u_2, v_2, P_2
  real(nk),dimension(:),allocatable:: u  , v  , P , u_ana, v_ana, Pres_ana
  real(nk),dimension(:),allocatable:: erreur_temps_P, erreur_temps_u,erreur_temps_v 
  real(nk):: Re, delta_x, delta_y, Lx, Ly
  
  integer :: nbr_it, pas_temps_choisi, it 
  real(nk):: delta_t, CFL
  type norme
     real(nk),dimension(:),allocatable:: erreur_temps_P, erreur_temps_u,erreur_temps_v 
  end type norme
  type(norme):: L1,L2,Lmax

  
!!!!!-----------------
!!!!!-----------------
  erreur_en_temps = 0
  RICH            = 1
  SOL_ANA         = 1
!!!!!----------------
  erreur_en_espace = 1
  RICH_O2          = 0 ; RICH_O4  = 0   
  SOL_ANA          = 1
!!!!!----------------
!!!!!----------------

  
  if (erreur_en_temps==1.and.RICH==1) then
     print*, 'En utilisant l extrapolation de Richadson'
!!!!--------------------------------------------------------------
     open(11,file="champs_uvp_en_fonction_delta_t",status="old")
     
     read(11,*) N, M, nbr_it, pas_temps_choisi,  delta_t,  Re
     read(11,*) 
     print*,'N=',N,'M=',M,'nbr_it=', nbr_it,'pas_temps_choisi=',pas_temps_choisi,'Delta_t=', delta_t,'Re=', Re
     NM = N*M  
 
      
!!$     read(11,*) N,M;  print*,N,M;  NM = N*M
!!$     read(11,*) nbr_it
!!$     print*, nbr_it
!!$     read(11,*) delta_t
!!$     read(11,*)
!!$     read(11,*)
!!!!---------------------------------
     allocate(L1%erreur_temps_u(nbr_it-1),L1%erreur_temps_v(nbr_it-1),L1%erreur_temps_P(nbr_it-1))
     allocate(L2%erreur_temps_u(nbr_it-1),L2%erreur_temps_v(nbr_it-1),L2%erreur_temps_P(nbr_it-1))
     allocate(Lmax%erreur_temps_u(nbr_it-1),Lmax%erreur_temps_v(nbr_it-1),Lmax%erreur_temps_P(nbr_it-1))

     allocate(u_1(NM),v_1(NM), P_1(NM),&
          u_2(NM),v_2(NM), P_2(NM),&
           u  (NM),v  (NM), P  (NM))

     u_1 = 0.0; v_1 = 0.0; P_1 = 0.0
     u_2 = 0.0; v_2 = 0.0; P_2 = 0.0
     u   = 0.0; v   = 0.0; P   = 0.0
!!$!!!!------------------------------
     open(21,file="deference_en_norme")
     open(22,file="ordre_en_norme_L1")
     open(23,file="ordre_en_norme_Lmax")
     open(24,file="ordre_en_norme_L2")
     write(21,*) '  t                         ',&
          'L1u                       ','L1v                       ','L1P                       ',&
          'Lmaxu                     ','Lmaxv                     ','LmaxP                     ',&
          'L2u                       ','L2v                       ','L2P                       '
!!!!---------------------------------
     do j=1,nbr_it
        do i=1,NM
           read(11,*) indice, u(i), v(i), P(i)         
        end do
        read(11,*)
        read(11,*)

        u_1 = u_2 ; v_1 = v_2 ; P_1 = P_2     
        u_2 = u   ; v_2 = v   ; P_2 = P

        if(j==1) then
           print*,j
        else           
           L1%erreur_temps_u(j-1) = sum( abs(u_1-u_2) ) 
           L1%erreur_temps_v(j-1) = sum( abs(v_1-v_2) )  
           L1%erreur_temps_P(j-1) = sum( abs(P_1-P_2) )

           Lmax%erreur_temps_u(j-1) = maxval( abs(u_1-u_2) ) 
           Lmax%erreur_temps_v(j-1) = maxval( abs(v_1-v_2) ) 
           Lmax%erreur_temps_P(j-1) = maxval( abs(P_1-P_2) )

           L2%erreur_temps_u(j-1) = sqrt( sum( abs(u_1 - u_2)**2) /sum( u_2**2 ) ) 
           L2%erreur_temps_v(j-1) = sqrt( sum( abs(v_1 - v_2)**2) /sum( v_2**2 ) )
           L2%erreur_temps_P(j-1) = sqrt( sum( abs(P_1 - P_2)**2) /sum( P_2**2 ) )
           
           print*,j,'L1  ', L1%erreur_temps_u(j-1)  , L1%erreur_temps_v(j-1)  , L1%erreur_temps_P(j-1)  
           print*,j,'L2  ', L2%erreur_temps_u(j-1)  , L2%erreur_temps_v(j-1)  , L2%erreur_temps_P(j-1)  
           print*,j,'Lmax', Lmax%erreur_temps_u(j-1), Lmax%erreur_temps_v(j-1), Lmax%erreur_temps_P(j-1)
           print*, 
!!!!--------------------------
        end if
        delta_t = delta_t/2.0_nk
     end do
     
     do i=2, nbr_it-1
        print*, L1%erreur_temps_u(i-1)/L1%erreur_temps_u(i)/2.0, &
                L1%erreur_temps_v(i-1)/L1%erreur_temps_v(i)/2.0, &
                L1%erreur_temps_P(i-1)/L1%erreur_temps_P(i)/2.0 

        print*, L2%erreur_temps_u(i-1)/L2%erreur_temps_u(i)/2.0, &
                L2%erreur_temps_v(i-1)/L2%erreur_temps_v(i)/2.0, &
                L2%erreur_temps_P(i-1)/L2%erreur_temps_P(i)/2.0 

        print*, Lmax%erreur_temps_u(i-1)/Lmax%erreur_temps_u(i)/2.0, &
                Lmax%erreur_temps_v(i-1)/Lmax%erreur_temps_v(i)/2.0, &
                Lmax%erreur_temps_P(i-1)/Lmax%erreur_temps_P(i)/2.0 
        print*,
      end do
        
     
     deallocate(L1%erreur_temps_u,L1%erreur_temps_v,L1%erreur_temps_P)
     deallocate(L2%erreur_temps_u,L2%erreur_temps_v,L2%erreur_temps_P)
     deallocate(Lmax%erreur_temps_u,Lmax%erreur_temps_v,Lmax%erreur_temps_P)

     deallocate(u_1, v_1, P_1)
     deallocate(u_2, v_2, P_2)
     deallocate(u  , v  , P  )

     close(21)
     close(22)
     close(23)
     close(24)
     close(11)
  end if

!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------

  if (erreur_en_temps==1.and.SOL_ANA == 1) then 
     print*,'confrontation avec les solutions anlytiques'
     open(11,file="champs_uvp_en_fonction_delta_t",status="old")
     read(11,*) N, M, nbr_it, pas_temps_choisi,  delta_t,  Re
     read(11,*) 
     print*,'N=',N,'M=',M,'nbr_it=', nbr_it,'pas_temps_choisi=',pas_temps_choisi,'Delta_t=', delta_t,'Re=', Re
     NM = N*M  
     it = 1
      
     allocate(u(NM), v(NM), P(NM))
     allocate(u_ana(NM), v_ana(NM), Pres_ana(NM))
     allocate(L1%erreur_temps_u(nbr_it),  L1%erreur_temps_v(nbr_it),  L1%erreur_temps_P(nbr_it))
     allocate(L2%erreur_temps_u(nbr_it),  L2%erreur_temps_v(nbr_it),  L2%erreur_temps_P(nbr_it))
     allocate(Lmax%erreur_temps_u(nbr_it),Lmax%erreur_temps_v(nbr_it),Lmax%erreur_temps_P(nbr_it))
     u  =0.0; v  =0.0; P  =0.0
    
     do j=1,nbr_it
        do i=1,NM
           read(11,*)indice, u(i), v(i), P(i)                    
        end do
        read(11,*)
        read(11,*)


        print*, u_ana(1), u(1) 
        L1%erreur_temps_u(j) = sum(abs( u - u_ana    ))  / (N*M)
        L1%erreur_temps_v(j) = sum(abs( v - v_ana    ))  / (N*M)
        L1%erreur_temps_P(j) = sum(abs( P - Pres_ana ))  / (N*M)
        
        Lmax%erreur_temps_u(j) = maxval(abs( u - u_ana   ))  
        Lmax%erreur_temps_v(j) = maxval(abs( v - v_ana   )) 
        Lmax%erreur_temps_P(j) = maxval(abs( P - Pres_ana))  
        
        
        L2%erreur_temps_u(j) = sqrt( sum( ( u - u_ana   )**2 )/sum( u_ana**2    )) 
        L2%erreur_temps_v(j) = sqrt( sum( ( v - v_ana   )**2 )/sum( v_ana**2    ))
        L2%erreur_temps_P(j) = sqrt( sum( ( P - Pres_ana)**2 )/sum( Pres_ana**2 ))
        
        print*, it,'L1  '  , L1%erreur_temps_u(j)  , L1%erreur_temps_v(j)  , L1%erreur_temps_P(j)  
        print*, it,'L2  '  , L2%erreur_temps_u(j)  , L2%erreur_temps_v(j)  , L2%erreur_temps_P(j)  
        print*, it,'Lmax'  , Lmax%erreur_temps_u(j), Lmax%erreur_temps_v(j), Lmax%erreur_temps_P(j)
        print*,           
         
        it = it + 1
          
     end do
     
    call plot_erreur
     
     do i=1, nbr_it-1
        print*,  log(L1%erreur_temps_u(i)/L1%erreur_temps_u(i+1))/log(2.0_nk),& 
                 log(L1%erreur_temps_v(i)/L1%erreur_temps_v(i+1))/log(2.0_nk),&
                 log(L1%erreur_temps_P(i)/L1%erreur_temps_P(i+1))/log(2.0_nk)
        
        print*,  log(L2%erreur_temps_u(i)/L2%erreur_temps_u(i+1))/log(2.0_nk),&
                 log(L2%erreur_temps_v(i)/L2%erreur_temps_v(i+1))/log(2.0_nk),&
                 log(L2%erreur_temps_P(i)/L2%erreur_temps_P(i+1))/log(2.0_nk)
        
        print*,  log(Lmax%erreur_temps_u(i)/Lmax%erreur_temps_u(i+1))/log(2.0_nk),&
                 log(Lmax%erreur_temps_v(i)/Lmax%erreur_temps_v(i+1))/log(2.0_nk),&
                 log(Lmax%erreur_temps_P(i)/Lmax%erreur_temps_P(i+1))/log(2.0_nk)
             
        print*,
     end do

     
         
      deallocate(u  , v  , P  )
      deallocate(u_ana  , v_ana  , Pres_ana  )
      deallocate(L1%erreur_temps_u,L1%erreur_temps_v,L1%erreur_temps_P)
      deallocate(L2%erreur_temps_u,L2%erreur_temps_v,L2%erreur_temps_P)
      deallocate(Lmax%erreur_temps_u,Lmax%erreur_temps_v,Lmax%erreur_temps_P)
   end if
  

!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------
!!!!---------------------------------------

if (erreur_en_espace==1.and.SOL_ANA == 1) then 

     open(11,file="champs_uvp_en_fonction_delta_x",status="old")
     read(11,*) nbr_it, pas_temps_choisi, delta_t, CFL
     read(11,*) indice,N,  M, delta_x, delta_y,  Re
     read(11,*)
     read(11,*)
     print*,'nbr_it=',nbr_it,'N=',N,'M=',M,'nbr_it=', nbr_it,'Delta_x=', delta_x,'Re=', Re, 'CFL=', CFL,&
          'pas_temps_choisi=', pas_temps_choisi, 'delta_t=', delta_t
  
     it = 1
      
     allocate(L1%erreur_temps_u(nbr_it),  L1%erreur_temps_v(nbr_it),  L1%erreur_temps_P(nbr_it))
     allocate(L2%erreur_temps_u(nbr_it),  L2%erreur_temps_v(nbr_it),  L2%erreur_temps_P(nbr_it))
     allocate(Lmax%erreur_temps_u(nbr_it),Lmax%erreur_temps_v(nbr_it),Lmax%erreur_temps_P(nbr_it))

     do j=1,nbr_it
        NM = N*M
        allocate(u(NM), v(NM), P(NM))
        allocate(u_ana(NM), v_ana(NM), Pres_ana(NM))
        u  =0.0; v  =0.0; P  =0.0
        u_ana=0.0; v_ana=0.0; Pres_ana=0.0
        
        
        
        do i=1,N*M
           read(11,*)indice, u(i), v(i), P(i)                    
        end do
        print*,N !, u(2*n-1), u_ana(2*n-1), u(2*n-1)- u_ana(2*n-1)

        read(11,*)
        read(11,*)
        read(11,*)
        read(11,*)
        
        L1%erreur_temps_u(j) = sum(abs( u - u_ana    ))  / (N*M)
        L1%erreur_temps_v(j) = sum(abs( v - v_ana    ))  / (N*M)
        L1%erreur_temps_P(j) = sum(abs( P - Pres_ana ))  / (N*M)
        
        Lmax%erreur_temps_u(j) = maxval(abs( u - u_ana   ))  
        Lmax%erreur_temps_v(j) = maxval(abs( v - v_ana   )) 
        Lmax%erreur_temps_P(j) = maxval(abs( P - Pres_ana))  
        
        
        L2%erreur_temps_u(j) = sqrt( sum( ( u - u_ana   )**2 )/sum( u_ana**2    )) 
        L2%erreur_temps_v(j) = sqrt( sum( ( v - v_ana   )**2 )/sum( v_ana**2    ))
        L2%erreur_temps_P(j) = sqrt( sum( ( P - Pres_ana)**2 )/sum( Pres_ana**2 ))
        
        print*, it,'L1  '  , L1%erreur_temps_u(j)  , L1%erreur_temps_v(j)  , L1%erreur_temps_P(j)  
        print*, it,'L2  '  , L2%erreur_temps_u(j)  , L2%erreur_temps_v(j)  , L2%erreur_temps_P(j)  
        print*, it,'Lmax'  , Lmax%erreur_temps_u(j), Lmax%erreur_temps_v(j), Lmax%erreur_temps_P(j)
        print*,           
         
        it = it + 1
        N= N*2
        M= M*2
        deallocate(u  , v  , P  )
        deallocate(u_ana  , v_ana  , Pres_ana  )      
     end do

   
     
     do i=1, nbr_it-1
        print*,  log(L1%erreur_temps_u(i)/L1%erreur_temps_u(i+1))/log(2.0_nk),& 
             log(L1%erreur_temps_v(i)/L1%erreur_temps_v(i+1))/log(2.0_nk),&
             log(L1%erreur_temps_P(i)/L1%erreur_temps_P(i+1))/log(2.0_nk)
        
        print*,  log(L2%erreur_temps_u(i)/L2%erreur_temps_u(i+1))/log(2.0_nk),&
             log(L2%erreur_temps_v(i)/L2%erreur_temps_v(i+1))/log(2.0_nk),&
             log(L2%erreur_temps_P(i)/L2%erreur_temps_P(i+1))/log(2.0_nk)
        
        print*,  log(Lmax%erreur_temps_u(i)/Lmax%erreur_temps_u(i+1))/log(2.0_nk),&
             log(Lmax%erreur_temps_v(i)/Lmax%erreur_temps_v(i+1))/log(2.0_nk),&
             log(Lmax%erreur_temps_P(i)/Lmax%erreur_temps_P(i+1))/log(2.0_nk)
             
        print*,
     end do

     
         
    
      deallocate(L1%erreur_temps_u,L1%erreur_temps_v,L1%erreur_temps_P)
      deallocate(L2%erreur_temps_u,L2%erreur_temps_v,L2%erreur_temps_P)
      deallocate(Lmax%erreur_temps_u,Lmax%erreur_temps_v,Lmax%erreur_temps_P)
   end if


!!!!----------------------------------------------------------------------
!!!!----------------------------------------------------------------------
!!!!----------------------------------------------------------------------
!!!!---------------------------------------------------------------------- 
 contains
!!$
!!!!----------------------------------------------------------------------
!!!!----------------------------------------------------------------------
!!!!----------------------------------------------------------------------
!!!!----------------------------------------------------------------------  

  subroutine plot_erreur
    implicit none
    
    integer :: i, j, l, ll
    real(nk),dimension(:),allocatable  :: x_som, y_som
    allocate(x_som(N*M), y_som(N*M))
    
    Lx  = pii!/2.0_nk 
    Ly  = pii!/2.0_nk
    
    delta_x   = Lx/N 
    delta_y   = Ly/M
    open (82,file="champs_erreur")
    do j=1,M
       do i=1,N
          ll= (j-1   )*(N     ) +i 
          x_som(ll) = (i-1)*delta_x  + delta_x/2.0_nk  +  pii/2.0_nk
          y_som(ll) = (j-1)*delta_y  + delta_y/2.0_nk  +  pii/2.0_nk               

          write(82,*) x_som(ll),y_som(ll), abs(  abs( u_ana(ll) )- abs(u(ll))),&
              abs(  abs( v_ana(ll) )- abs(v(ll)) ), &
              abs( abs(Pres_ana(ll)) - abs(P(ll) ) )    
       end do
       write(82,*) 
    end do
          
    
    close(82)
    
    deallocate(x_som,y_som)
  end subroutine plot_erreur
!!!!----------------------------------------------------------------------
!!!!----------------------------------------------------------------------
!!!!----------------------------------------------------------------------
!!!!----------------------------------------------------------------------  
  

  
end program poste_traitement
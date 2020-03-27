program poste_traitement_fonction_viscometrique      
  implicit none

  integer, parameter :: nk=selected_real_kind( 15, 307  )
  integer :: i, j, l, k, s, ll, ii, jj 
  integer :: N1, N2, N3   
  type variable
     real(nk),dimension(:),allocatable:: var
  end type variable
  type(variable), dimension(20) :: W
  real(nk),dimension(:),allocatable::x3_som, x2_som
  character(len=100):: work
  integer :: NG
  real(nk):: work1,  Reynolds, Beta, dx3, dx2
  real(nk):: viscosity_total

!!$  W(1)  = u1
!!$  W(2)  = u2
!!$  W(3)  = u3  
!!$  W(4)  = tau11
!!$  W(5)  = tau12
!!$  W(6)  = tau13
!!$  W(7)  = tau22
!!$  W(8)  = tau23
!!$  W(9)  = tau33
!!$  W(10) = T
!!$  W(11) = N1
!!$  W(12) = N2
!!$  W(13) = Shear rate  
!!$  W(14) = planar elongational rate 
!!$  W(15) = planar elongational rate_du2/dx2 
!!$  W(16) = planar elongational rate_du3/dx3
!!$  W(15) = Shear viscosity 
!!$  W(16) = planar elongational viscosity


!!!!!!!!!!!!!!!!!!$ lecture N1,2,3 taille maillage   
  open(11,file="donnees",status="old")  
  read(11,*) 
  read(11,*) 
  read(11,*) work, N1 ;  read(11,*) work, N2 ;  read(11,*) work, N3
  do i=1, 18
     read(11,*) 
  end do
  read(11,*) work, Reynolds
  read(11,*) 
  read(11,*) 
  read(11,*) 
  read(11,*) work, Beta
  close(11) 

  N1 = N1+1
  N2 = N2+1
  N3 = N3+1  

  do i=1, 20
     allocate( W(i)%var( N3*N2 )) 
  end do
  allocate(x3_som(N2*N3), x2_som(N2*N3))    

!!!!!!!!!!!!!!!!$ lecture : sectional cross fields  
  open(82,file="champs_UT_x2x3_1_2")    
  do j= 1,N2
     do k= 1,N3                 
        l = (j-1)*N3 + k                
        read(82,*) x3_som(l), x2_som(l) , ( W(ii)%var(l), ii=1, 10), work1,  W(13)%var(l)        
     end do
  end do
  close(82)
 
  dx3 =  x3_som(2) 
  dx2 =  x2_som(N3+1)



  
!!!!!!!!!!!!!!!!!!!!!!!!!!!$ profile viscosité
  
  open (81,file="prof_visco_x2")


  write(81,*)   &
       '  x2                        ' , & 
       'u_1                       ' , & 
       'u_2                       ' , & 
       'u_3                       ' , &  
       't_11                      ' , & 
       't_12                      ' , & 
       't_13                      ' , & 
       't_22                      ' , & 
       't_23                      ' , & 
       't_33  	              ' , & 
       'visco                      ' , & 
       'gamma                     ' , & 
       'N1                        ' , & 
       'N2                        ' , & 
       'eps_point                 ' , & 
       'visco_elong               ' , & 
       'eps_point_x2              ' , & 
       'eps_point_x3              ' , & 
       'phi1_coefN1               ' , & 
       'phi2_coefN2               '

  k= (N3-1)/2 + 1
  do j= 1,N2         
     
     
     ll = (j-2)*N3 + k      !!!! pour coef phi au centre du canal           
          
     l = (j-1)*N3 + k                
     
     
     if ( j == 1 .or. j == N2 ) then                
        
        if ( j == 1 ) then                                     
           W(15)%var(l)  =   (   4.* W(2)%var(l+N3) - W(2)%var(l+2*N3)  )  /2./dx2 
        end if
        if ( j == N2 ) then                                       
           W(15)%var(l)  =   ( - 4.* W(2)%var(l-N3) + 4.* W(2)%var(l-2*N3)  ) /2./dx2  
        end if

     else                        
        W(15)%var(l)   =  ( W(2)%var(l+N3) - W(2)%var(l-N3) ) /2./dx2              
     end if

     W(16)%var(l)   =  ( W(3)%var(l+1) - W(3)%var(l-1) ) /2./dx3              

     W(14)%var(l)= sqrt( W(15)%var(l)**2. +  W(16)%var(l)**2. )
     
     
     if ( j ==  ((N2-1)/2+1) ) then                   
        write(81,*)  x2_som(l), & !  1 
             ( W(ii)%var(l),ii=1,9), & ! 2-10 
             1._nk , & ! eta   11
             W(13)%var(ll) , &   ! 12 
             W(4)%var(l) -  W(7)%var(l)  , &  ! N1   13   
             abs(W(7)%var(ll) -  W(9)%var(ll)) , &    ! N2   14
             W(14)%var(l)  ,&    ! 15             
             3._nk  ,&  !  eta_e    16 
             W(15)%var(l) , &  ! 17
             W(16)%var(l) , &      ! 18     
             ( W(4)%var(ll) -  W(7)%var(ll) )/  W(13)%var(ll)**2. , &  ! 19   phi_1 vs N1                  
             abs(W(7)%var(ll) -  W(9)%var(ll))/  W(13)%var(ll)**2.   !  20  phi_2 vs N2     
     else        
        
        viscosity_total =  beta +  abs( W(5)%var(l))/ W(13)%var(l)
        
        if ( viscosity_total > 1._nk )  viscosity_total = 1._nk 
         
        
        write(81,*)  x2_som(l), & !  1 
             ( W(ii)%var(l),ii=1,9), & ! 2-10 
             viscosity_total, & ! eta   11  
             W(13)%var(l) , &   ! 12 
             W(4)%var(l) -  W(7)%var(l)    , &  ! N1   13   
             abs(W(7)%var(l) -  W(9)%var(l))    , &    ! N2   14
             W(14)%var(l) ,&    ! 15             
             3.*beta + ( W(7)%var(l) -  W(9)%var(l) ) / W(14)%var(l),&   !  eta_e    16 
             W(15)%var(l) , &  ! 17
             W(16)%var(l) , &    ! 18
             ( W(4)%var(l) -  W(7)%var(l) )/  W(13)%var(l)**2. , &  ! 19   phi_1 vs N1                  
             abs(W(7)%var(l) -  W(9)%var(l))/  W(13)%var(l)**2.   !  20  phi_2 vs N2 
     end if
     
     
  end do
  close(81)
  
  
  
  
  
  
  print*,  'je suis la' , Beta, Reynolds
  stop  
  























   end program poste_traitement_fonction_viscometrique

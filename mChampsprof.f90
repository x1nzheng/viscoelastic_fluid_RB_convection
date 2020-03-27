
#include "definitions.h"

module mChampsprof
  use mBase
  use mPression
  
  !Declarations
  implicit none

contains


  

  
  
  subroutine profil_plans (u1, u2)
    implicit none    
    real(nk),dimension(0:N1+1,0:N2+1),intent(in)          :: u1, u2 
    real(nk)                                              :: xi, yj, zk
#if (DIMENSION_GEO == 2)
    integer :: i, j, k

    real(nk),dimension(N2):: u11, tau_11, tau_12
    real(nk),dimension(N1):: u22
    real(nk)::u1_max = 0._nk , u2_max = 0._nk, x1_som_max=0._nk, x2_som_max=0._nk 
                 u11  = 0._nk  ; u22 = 0._nk ; tau_11 = 0._nk ; tau_12 = 0._nk
     
    open (81,file="prof_med_u1") 

    i= N1/2
    u1_max = u1(i,1)
    do j= 1,N2     
        yj = 0.              
          u11(j) = u1(i,j)
       tau_11(j) = var%W(t11)%var(i,j)
       tau_12(j) = var%W(t12)%var(i,j)

#if(VISCOELASTIC_MODELS == 0)
       write(81,*) yj, u1(i,j)
#elif(VISCOELASTIC_MODELS > 0)
       write(81,*) yj, u1(i,j),var%W(T)%var(i,j), var%W(t11)%var(i,j), var%W(t12)%var(i,j), var%W(t22)%var(i,j)
#endif 
    
       if(u1_max < u1(i,j) ) then 
          u1_max = u1(i,j)
          x2_som_max = j
       end if
       yj= yj + dx2_L(i,j)
    end do
    close(81)
    if (print_erreur_temps ==1 )  print*, ' '
    if (print_erreur_temps ==1 )  print"(a,(f10.6), 4X, a, (f10.6),4X, a, 2(f10.6))",&
         'Max_u1=',maxval(u11),'Min_u1=',minval(u11),'x2_som_max=',x2_som_max,u1_max
     
    
    
    open (81,file="prof_med_u2") 
    j=N2/2 
    i=1
    u2_max = u2(i,j)
    do i= 1,N1         
       xi = 0.
       u22(i) = u2(i,j)
       write(81,*) i, u2(i,j),var%W(T)%var(i,j)
       if(u2_max < u2(i,j) ) then 
          u2_max = u2(i,j)
          x1_som_max = i
       end if
       xi= xi + dx1_L(i,j)
    end do
    close(81)
    if (print_erreur_temps ==1 )   print"(a,(f10.6), 4X, a, (f10.6),4X, a, (f10.6))",&
         'Max_u2=', maxval(u22), 'Min_u2=',  minval(u22), 'x1_som_max=',  x1_som_max
    if (print_erreur_temps ==1 )   print*, ' '
    
#elif (DIMENSION_GEO == 3)
    

#endif 
  end subroutine profil_plans
  

  

  












  subroutine champs_uvpT (u1, u2, P, T, t11, t12, t22, indice1, indice2, it, u3, t13, t23, t33 )
    implicit none  
    integer :: i, j, k, ii, k1, k2, j1, j2
    integer,intent(in):: indice1, indice2, it
    character(9) :: cTemp
    real(nk)  :: Total_Nu_paroi1, Total_Nu_paroi2
    real(nk), dimension(1:N1) :: Nu_paroi1, Nu_paroi2
    integer ::  size
    type(variable), intent(in):: u1, u2, T, t11, t12, t22
    type(variable), intent(in):: P
    type(variable), optional, intent(in):: u3, t13, t23, t33
    real(nk), dimension(:,:),allocatable:: vortic
    real(nk)::                             xi, yi, zi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

#if (DIMENSION_GEO == 2)    
    call VTK_TEM(u1%var, u2%var, T%var, t11%var, t12%var, t22%var, it)
!!! recording the velocity temperature and stress on the y direction central line
    write( cTemp,'(i9)') it

    open (1,file=trim(adjustl(cTemp))//"_champs_central_yline_uT")
        yi = 0.
    do j= 1,N2
        i= N1/2
         write(1,*)  yi, u1 %var(i+1,j), u2 %var(i+1,j), t11%var(i+1,j), t12%var(i+1,j), t22%var(i+1,j), T  %var(i+1,j) 
         yi= yi + dx2_R(i,j)
    end do
    close(1)

    open (1,file=trim(adjustl(cTemp))//"_champs_central_xline_uT")
        xi = 0.
    do i= 1,N1
        j= N2/2
         write(1,*) xi, u1 %var(i,j+1), u2 %var(i,j+1), t11%var(i,j+1), t12%var(i,j+1),&
                        t22%var(i,j+1), T  %var(i,j+1)
         xi = xi + dx1_R(i,j)    
    end do
    close(1)

    open (81,file=trim(adjustl(cTemp))//"_champ_P",  Form = 'UNFORMATTED' )
!    open (81,file=trim(adjustl(cTemp))//"_champ_P" )
        yi = 0.5_nk*dx2_R(1,1)
    do j= 1,N2-1
       xi = 0.5_nk*dx1_R(1,1)
       do i= 1,N1-1         
          write(81) xi, yi, P%var(i,j), P%var0(i,j) 
          xi = xi + 0.5_nk*(dx1_L(i+1,j)+dx1_R(i+1,j))
       end do 
       yi = yi + 0.5_nk*(dx2_L(i,j+1)+dx2_R(i,j+1))
       write(81)"  "
    end do
    close(81)


    open (82,file=trim(adjustl(cTemp))//"_champs_UT",  Form = 'UNFORMATTED' )
!    open (82,file=trim(adjustl(cTemp))//"_champs_UT" )
        yi = 0. 
       do j= 1,N2
         xi = 0.
          do i= 1,N1         
            write(82) xi, yi, &
            u1%var (i,j), u2%var (i,j), t11%var (i,j), t12%var (i,j), t22%var (i,j), T%var (i,j), &
            u1%var0(i,j), u2%var0(i,j), t11%var0(i,j), t12%var0(i,j), t22%var0(i,j), T%var0(i,j)
            xi = xi + dx1_L(i,j)
          end do
           write(82)"  "      
       yi = yi + dx2_L(i,j)
       end do
    close(82)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open (82,file="champs_central_yline_uT") 
        yi = 0.
    do j= 1,N2
        i= N1/2
         write(82,*) yi, u1 %var(i+1,j), u2 %var(i+1,j), t11%var(i+1,j), t12%var(i+1,j),&
                        t22%var(i+1,j), T  %var(i+1,j)
         write(82,*)"  "    
         yi= yi + dx2_R(i,j)
    end do
    close(82)

    open (82,file="champs_central_xline_uT")
         xi = 0.
    do i= 1,N1
        j= N2/2
         write(82,*) xi, u1 %var(i,j+1), u2 %var(i,j+1), t11%var(i,j+1), t12%var(i,j+1),&
                        t22%var(i,j+1), T  %var(i,j+1)
         write(82,*)"  "  
         xi = xi + dx1_R(i,j)    
    end do
    close(82)

    open (82,file="champs_UT",  Form = 'UNFORMATTED' )
!    open (82,file="champs_UT" )
       yi = 0.
        do j= 1,N2
         xi = 0.
         do i= 1,N1         
            write(82) xi, yi, &
            u1%var (i,j), u2%var (i,j), t11%var (i,j), t12%var (i,j), t22%var (i,j), T%var (i,j), &
            u1%var0(i,j), u2%var0(i,j), t11%var0(i,j), t12%var0(i,j), t22%var0(i,j), T%var0(i,j)
            xi = xi + dx1_R(i,j)
         end do
            write(82)"  "     
            yi = yi + dx2_R(i,j)
        end do
    close(82)


      open (81,file="champ_P",  Form = 'UNFORMATTED' )
!      open (82,file="champ_P" )
        yi = 0.5_nk*dx2_R(1,1)
    do j= 1,N2-1
       xi = 0.5_nk*dx1_R(1,1)
       do i= 1,N1-1         
          write(81) xi, yi, P%var(i,j), P%var0(i,j) 
          xi = xi + 0.5_nk*(dx1_L(i+1,j)+dx1_R(i+1,j))
       end do 
       yi = yi + 0.5_nk*(dx2_L(i,j+1)+dx2_R(i,j+1))
       write(81)"  "
    end do
    close(81)

   call fonction_courant(var%W(1)%var,var%W(2)%var)
   call Calculs_Local_Nusselt (T%var, it)
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#elif (DIMENSION_GEO == 3)

    call VTK_TEM(u1%var, u2%var, T%var, t11%var, t12%var, t22%var, it, u3%var, t13%var, t23%var, t33%var)
!!!!! vitesse, \tau_ij et T    
    pause 
    open (82,file="champs_UT_x1x2_1_2") 
    open (83,file="champs_UT_x1x2_1_4") 
    open (84,file="champs_UT_x1x2_3_4")        
    k = N3/2 + 1              
    k1= N3/4 + 1 
    k2= N3/4*3 + 1 
    
    do j= 1,N2
       do i= 1,N1         
          write(82,*) i, j,(var%W(ii)%var(i,j,k ), ii=1, 10) 
          write(83,*) i, j,(var%W(ii)%var(i,j,k1), ii=1, 10) 
          write(84,*) i, j,(var%W(ii)%var(i,j,k2), ii=1, 10)                               
       end do
       write(82,*)"  "          
       write(83,*)"  "          
       write(84,*)"  "          
    end do
    close(82)
    close(83)
    close(84)






    open (82,file="champs_UT_x1x3_1_2") 
    open (83,file="champs_UT_x1x3_1_4") 
    open (84,file="champs_UT_x1x3_3_4")        
    j = N2/2 + 1              
    j1= N2/4 + 1 
    j2= N2/4*3 + 1
    
    do k= 1,N3
       do i= 1,N1         
          write(82,*) i, k, (var%W(ii)%var(i,j,k ), ii=1, 10) 
          write(83,*) i, k, (var%W(ii)%var(i,j1,k), ii=1, 10) 
          write(84,*) i, k, (var%W(ii)%var(i,j2,k), ii=1, 10) 
          
       end do
       write(82,*)"  "          
       write(83,*)"  "          
       write(84,*)"  "          
    end do
    close(82)
    close(83)
    close(84)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! sur l'ensemble du domaine !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open (82,file="champs_3D")     
    do k= 1, N3
       do j= 1, N2
          do i=1, N1   

#if(VISCOELASTIC_MODELS == 0)
             write(82,*) i, j, k, u1%var(i,j,k), u2%var(i,j,k), var%W(3)%var(i,j,k), var%W(10)%var(i,j,k)  
#elif(VISCOELASTIC_MODELS > 0)

             write(82,*) i, j, k, u1%var(i,j,k), u2%var(i,j,k), u3%var(i,j,k), &
              t11%var(i,j,k), t12%var(i,j,k), t13%var(i,j,k), &
              t22%var(i,j,k), t23%var(i,j,k), t33%var(i,j,k) , &
              T%var(i,j,k), P%var(i,j,k) 
#endif 

          end do
          write(82,*)"  "          
       end do
       write(82,*)"  " 
    end do
    close(82)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open (82,file="champs_3D_P")     
    do k= 1, N3-1
       do j= 1, N2-1
          do i=1, N1-1   
             
             write(82,*) i, j, k, P%var(i,j,k), P%var0(i,j,k) 
             
          end do
          write(82,*)"  "          
       end do
       write(82,*)"  "
    end do
    close(82)

#endif
    
  end subroutine champs_uvpT

  




  SUBROUTINE VTK_TEM(u1, u2, T, t11, t12, t22, it, u3, t13, t23, t33)
     IMPLICIT NONE

#if (DIMENSION_GEO == 2)
   integer, intent(in)                             :: it
   real(nk), dimension(0:N1+1,0:N2+1), intent(in)  :: u1, u2, T, t11, t12, t22
   real(nk), dimension(0:N1+1,0:N2+1), optional    :: u3, t13, t23, t33
#elif (DIMENSION_GEO == 3)   
   integer, intent(in)                                    :: it
   real(nk), dimension(0:N1+1,0:N2+1,0:N3+1), intent(in)  :: u1, u2, T, t11, t12, t22, u3, t13, t23, t33
#endif

     character(9) :: cTemp
     character(len=1), parameter :: eol=achar(10)
      !  character(len=7), parameter :: REAL_PREC = " float "
     character(len=8), parameter :: REAL_PREC = " double "
     integer, parameter          :: REALSIZE = SELECTED_REAL_KIND(12)
!       character(*), parameter :: bin_convert=OUTPUT_CONVERT
    integer, parameter           :: nunit=500
    integer                      :: n,m,nbr_dset
    integer                      :: err,reclen,rec_unit
    character(len=500)           :: header_main
    character(len=200)           :: header_grid, header_data
    character(len=100)           :: header_x,header_y,header_z
    character(len=100)           :: tmp_buffer 
!      character(len=12)       :: c_time
!      character(len=20)       :: itoas,name
    character(len=20)            :: nxbuf, nybuf, nzbuf, nxyzbuf
!      real(kind=4) , dimension(:,:,:,:), allocatable       :: vector
!      character(len=200),dimension(:), allocatable :: header_dset
    real(kind=REALSIZE) , dimension(:,:,:,:), allocatable       :: vector
    real(kind=REALSIZE) , dimension(:,:,:,:), allocatable       :: energy
    character(len=200) :: header_u, header_v, header_w, header_t
    character(len=200) :: header_t11, header_t12, header_t22
#if ( DIMENSION_GEO == 3 )
    character(len=200) ::  header_t13 , header_t23, header_t33  
#endif 
    character(len=200) :: header_E_E, header_D_E, header_F_E, header_V_E, header_G_E
    real(kind=REALSIZE), dimension(:), allocatable    :: XL, YL, ZL
    integer                 :: NX, NY, NZ, i, j



#if ( DIMENSION_GEO == 2 )
    write( cTemp,'(i9)') it
!-------------------------------------------------------------
! prepare coordinates
!-------------------------------------------------------------

       NX = N1; ALLOCATE(XL(1:NX)); XL(1) = 0.
       NY = N2; ALLOCATE(YL(1:NY)); YL(1) = 0.
       NZ = 1 ; ALLOCATE(ZL(1:NZ)); ZL(1:NZ)=0. 

       do i =2, NX
       XL(i) = XL(i-1) + dx1_R(i-1,1) 
       enddo
       do j =2, NY
       YL(j) = YL(j-1) + dx2_R(1,j-1) 
       enddo

!-------------------------------------------------------------
! Allocation
!-------------------------------------------------------------

       allocate(vector(1:7,1:NX,1:NY,1:NZ))
       allocate(energy(1:5,1:NX,1:NY,1:NZ))

!-------------------------------------------------------------
! prepare header
!-------------------------------------------------------------

       write(tmp_buffer,'(i8," ",i8," ",i8)') NX,NY,NZ
       write(nxbuf,'(i8)') NX; write(nybuf,'(i8)') NY;
       write(nzbuf,'(i8)') NZ; write(nxyzbuf,'(i10)') NX*NY*NZ;
!      write(c_time,'(E12.6)') time

!--- Main header

    header_main="# vtk DataFile Version 3.0"//eol//&
                   "vtk data for time "//eol//&
                   "BINARY"//eol
!                  "vtk data for time "//trim(c_time)//eol//&
!--- grid header : rectilinear grid

    header_grid="DATASET RECTILINEAR_GRID"//eol//&
                "DIMENSIONS "//trim(tmp_buffer)//eol
!   header_x="X_COORDINATES "//itoas(nx)//REAL_PREC//eol
!   header_y="Y_COORDINATES "//itoas(ny)//REAL_PREC//eol
!   header_z="Z_COORDINATES "//itoas(nz)//REAL_PREC//eol
    header_x="X_COORDINATES "//trim(nxbuf)//REAL_PREC//eol
    header_y="Y_COORDINATES "//trim(nybuf)//REAL_PREC//eol
    header_z="Z_COORDINATES "//trim(nzbuf)//REAL_PREC//eol

    !---  data header

!   header_data="POINT_DATA "//trim(itoas(NX*NY*NZ))//eol
    header_data="POINT_DATA "//trim(nxyzbuf)//eol

!-------------------------------------------------------------
!     Create VTK binary file
!-------------------------------------------------------------

!--- header_main, header_grid and coordinates
       header_u ="SCALARS  U"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_v ="SCALARS  V"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_w ="SCALARS  W"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t ="SCALARS  T"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t11 ="SCALARS T11"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t12 ="SCALARS T12"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t22 ="SCALARS T22"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_E_E ="SCALARS E_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_D_E ="SCALARS D_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_F_E ="SCALARS F_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_V_E ="SCALARS V_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_G_E ="SCALARS G_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

               
    reclen= len_trim(header_main) + len_trim(header_grid) + len_trim(header_x) + &
            NX*REALSIZE + len(eol) + len_trim(header_y) + NY*REALSIZE + len(eol) + &
            len_trim(header_z) + NZ*REALSIZE + len(eol) + len_trim(header_data) + &
            len_trim(header_u) + len_trim(header_v) + len_trim(header_w) + &
            len_trim(header_t11) + len_trim(header_t12) + len_trim(header_t22)+ &
            len_trim(header_t) + NX*NY*NZ*7*REALSIZE + len(eol)*7 + &
            len_trim(header_E_E)*5 + NX*NY*NZ*5*REALSIZE + len(eol)*5

        vector(1,1:NX,1:NY,1) = u1 (1:NX,1:NY)
        vector(2,1:NX,1:NY,1) = u2 (1:NX,1:NY)
        vector(3,1:NX,1:NY,1) = 0.
        vector(4,1:NX,1:NY,1) = t11(1:NX,1:NY)
        vector(5,1:NX,1:NY,1) = t12(1:NX,1:NY)
        vector(6,1:NX,1:NY,1) = t22(1:NX,1:NY)
        vector(7,1:NX,1:NY,1) = T  (1:NX,1:NY)

        energy(1,1:NX,1:NY,1) = E_global     (1:NX,1:NY)
        energy(2,1:NX,1:NY,1) = D_diffusion  (1:NX,1:NY)
        energy(3,1:NX,1:NY,1) = F_thermal    (1:NX,1:NY)
        energy(4,1:NX,1:NY,1) = V_dissipation(1:NX,1:NY)
        energy(5,1:NX,1:NY,1) = G_exchange   (1:NX,1:NY)
!
    open(nunit, file=trim(adjustl(cTemp))//'_sortie.vtk', form='UNFORMATTED', iostat=err, status='replace', & 
                access='direct',  convert    = 'BIG_ENDIAN',  recl=reclen )

       write(nunit,rec=1) trim(header_main), trim(header_grid), trim(header_x), &
           XL(1:NX),eol,trim(header_y),YL(1:NY),eol,trim(header_z),ZL(1:NZ),eol, &
           trim(header_data), &
           trim(header_u)    ,vector(1,1:NX,1:NY,1:NZ),eol, &
           trim(header_v)    ,vector(2,1:NX,1:NY,1:NZ),eol, &
           trim(header_w)    ,vector(3,1:NX,1:NY,1:NZ),eol, &
           trim(header_t11)  ,vector(4,1:NX,1:NY,1:NZ),eol, &
           trim(header_t12)  ,vector(5,1:NX,1:NY,1:NZ),eol, &
           trim(header_t22)  ,vector(6,1:NX,1:NY,1:NZ),eol, &
           trim(header_t)    ,vector(7,1:NX,1:NY,1:NZ),eol, &
           trim(header_E_E)  ,energy(1,1:NX,1:NY,1:NZ),eol, &
           trim(header_D_E)  ,energy(2,1:NX,1:NY,1:NZ),eol, &
           trim(header_F_E)  ,energy(3,1:NX,1:NY,1:NZ),eol, &
           trim(header_V_E)  ,energy(4,1:NX,1:NY,1:NZ),eol, &
           trim(header_G_E)  ,energy(5,1:NX,1:NY,1:NZ),eol
           
    close(nunit)
!
#elif ( DIMENSION_GEO == 3 )
     
    write( cTemp,'(i9)') it
!-------------------------------------------------------------
! prepare coordinates
!-------------------------------------------------------------

       NX = N1; ALLOCATE(XL(1:NX)); XL(1) = 0.
       NY = N2; ALLOCATE(YL(1:NY)); YL(1) = 0.
       NZ = N3; ALLOCATE(ZL(1:NZ)); ZL(1) = 0. 

       do i =2, NX
       XL(i) = XL(i-1) + dx1_R(i-1,1,1) 
       enddo
       do j =2, NY
       YL(j) = YL(j-1) + dx2_R(1,j-1,1) 
       enddo
       do k =2, NZ
       ZL(k) = ZL(k-1) + dx3_R(1,1,k-1) 
       enddo

!-------------------------------------------------------------
! Allocation
!-------------------------------------------------------------

       allocate(vector(1:10,1:NX,1:NY,1:NZ))
       allocate(energy(1:5,1:NX,1:NY,1:NZ))

!-------------------------------------------------------------
! prepare header
!-------------------------------------------------------------

       write(tmp_buffer,'(i8," ",i8," ",i8)') NX,NY,NZ
       write(nxbuf,'(i8)') NX; write(nybuf,'(i8)') NY;
       write(nzbuf,'(i8)') NZ; write(nxyzbuf,'(i10)') NX*NY*NZ;
!      write(c_time,'(E12.6)') time

!--- Main header

    header_main="# vtk DataFile Version 3.0"//eol//&
                   "vtk data for time "//eol//&
                   "BINARY"//eol
!                  "vtk data for time "//trim(c_time)//eol//&
!--- grid header : rectilinear grid

    header_grid="DATASET RECTILINEAR_GRID"//eol//&
                "DIMENSIONS "//trim(tmp_buffer)//eol
!   header_x="X_COORDINATES "//itoas(nx)//REAL_PREC//eol
!   header_y="Y_COORDINATES "//itoas(ny)//REAL_PREC//eol
!   header_z="Z_COORDINATES "//itoas(nz)//REAL_PREC//eol
    header_x="X_COORDINATES "//trim(nxbuf)//REAL_PREC//eol
    header_y="Y_COORDINATES "//trim(nybuf)//REAL_PREC//eol
    header_z="Z_COORDINATES "//trim(nzbuf)//REAL_PREC//eol

    !---  data header

!   header_data="POINT_DATA "//trim(itoas(NX*NY*NZ))//eol
    header_data="POINT_DATA "//trim(nxyzbuf)//eol

!-------------------------------------------------------------
!     Create VTK binary file
!-------------------------------------------------------------

!--- header_main, header_grid and coordinates
       header_u ="SCALARS  U"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_v ="SCALARS  V"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_w ="SCALARS  W"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t ="SCALARS  T"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t11 ="SCALARS T11"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t12 ="SCALARS T12"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t13 ="SCALARS T13"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t22 ="SCALARS T22"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t23 ="SCALARS T23"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol
       header_t33 ="SCALARS T33"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_E_E ="SCALARS E_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_D_E ="SCALARS D_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_F_E ="SCALARS F_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_V_E ="SCALARS V_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

                   header_G_E ="SCALARS G_E"//REAL_PREC//eol//&
                   "LOOKUP_TABLE default"//eol

               
    reclen= len_trim(header_main) + len_trim(header_grid) + len_trim(header_x) + &
            NX*REALSIZE + len(eol) + len_trim(header_y) + NY*REALSIZE + len(eol) + &
            len_trim(header_z) + NZ*REALSIZE + len(eol) + len_trim(header_data) + &
            len_trim(header_u) + len_trim(header_v) + len_trim(header_w) + &
            len_trim(header_t11) + len_trim(header_t12) + len_trim(header_t13) + &
            len_trim(header_t22) + len_trim(header_t23) + len_trim(header_t33) + &
            len_trim(header_t) + NX*NY*NZ*10*REALSIZE + len(eol)*10 + &
            len_trim(header_E_E)*5 + NX*NY*NZ*5*REALSIZE + len(eol)*5

        vector(1,1:NX,1:NY,1:NZ) = u1 (1:NX,1:NY,1:NZ)
        vector(2,1:NX,1:NY,1:NZ) = u2 (1:NX,1:NY,1:NZ)
        vector(3,1:NX,1:NY,1:NZ) = u3 (1:NX,1:NY,1:NZ)
        vector(4,1:NX,1:NY,1:NZ) = t11(1:NX,1:NY,1:NZ)
        vector(5,1:NX,1:NY,1:NZ) = t12(1:NX,1:NY,1:NZ)
        vector(6,1:NX,1:NY,1:NZ) = t13(1:NX,1:NY,1:NZ)
        vector(7,1:NX,1:NY,1:NZ) = t22(1:NX,1:NY,1:NZ)
        vector(8,1:NX,1:NY,1:NZ) = t23(1:NX,1:NY,1:NZ)
        vector(9,1:NX,1:NY,1:NZ) = t33(1:NX,1:NY,1:NZ)
        vector(10,1:NX,1:NY,1:NZ) = T (1:NX,1:NY,1:NZ)

        energy(1,1:NX,1:NY,1:NZ) = E_global     (1:NX,1:NY,1:NZ)
        energy(2,1:NX,1:NY,1:NZ) = D_diffusion  (1:NX,1:NY,1:NZ)
        energy(3,1:NX,1:NY,1:NZ) = F_thermal    (1:NX,1:NY,1:NZ)
        energy(4,1:NX,1:NY,1:NZ) = V_dissipation(1:NX,1:NY,1:NZ)
        energy(5,1:NX,1:NY,1:NZ) = G_exchange   (1:NX,1:NY,1:NZ)
!
    open(nunit, file=trim(adjustl(cTemp))//'_sortie.vtk', form='UNFORMATTED', iostat=err, status='replace', & 
                access='direct',  convert    = 'BIG_ENDIAN',  recl=reclen )

       write(nunit,rec=1) trim(header_main), trim(header_grid), trim(header_x), &
           XL(1:NX),eol,trim(header_y),YL(1:NY),eol,trim(header_z),ZL(1:NZ),eol, &
           trim(header_data), &
           trim(header_u)    ,vector(1,1:NX,1:NY,1:NZ),eol, &
           trim(header_v)    ,vector(2,1:NX,1:NY,1:NZ),eol, &
           trim(header_w)    ,vector(3,1:NX,1:NY,1:NZ),eol, &
           trim(header_t11)  ,vector(4,1:NX,1:NY,1:NZ),eol, &
           trim(header_t12)  ,vector(5,1:NX,1:NY,1:NZ),eol, &
           trim(header_t13)  ,vector(6,1:NX,1:NY,1:NZ),eol, &
           trim(header_t22)  ,vector(7,1:NX,1:NY,1:NZ),eol, &
           trim(header_t23)  ,vector(8,1:NX,1:NY,1:NZ),eol, &
           trim(header_t33)  ,vector(9,1:NX,1:NY,1:NZ),eol, &
           trim(header_t)    ,vector(10,1:NX,1:NY,1:NZ),eol, &
           trim(header_E_E)  ,energy(1,1:NX,1:NY,1:NZ),eol, &
           trim(header_D_E)  ,energy(2,1:NX,1:NY,1:NZ),eol, &
           trim(header_F_E)  ,energy(3,1:NX,1:NY,1:NZ),eol, &
           trim(header_V_E)  ,energy(4,1:NX,1:NY,1:NZ),eol, &
           trim(header_G_E)  ,energy(5,1:NX,1:NY,1:NZ),eol
           
    close(nunit)
!
#endif
  END SUBROUTINE VTK_TEM






  subroutine fonction_courant(u1, u2)
   implicit none
   
   integer :: i, j
   real(nk),dimension(0:N1+1,0,N2+1),intent(in):: u1, u2
   type(variable)  :: FC
   real(nk) :: somme, RH
   character(9) :: cTemp
   real(nk)::                             xi, yi, zi
#if (DIMENSION_GEO == 2)
 REAL(nk), DIMENSION(:,:), ALLOCATABLE :: FCVPRX, FCVPRX1
 REAL(nk), DIMENSION(:), ALLOCATABLE :: FCVPX
#endif


#if (DIMENSION_GEO == 2)
   FC%NNN    = var%W(7)%NNN 
   FC%NI     = var%W(7)%NI
   FC%NJ     = var%W(7)%NJ
!   FC%indice = 20 !var%W(7)%indice
#elif (DIMENSION_GEO == 3)
   FC%NNN = var%W(11)%NNN
   FC%NI     = var%W(7)%NI
   FC%NJ     = var%W(7)%NJ
   FC%NK     = var%W(7)%NK
#endif 
#if (DIMENSION_GEO == 2)
    ALLOCATE( FCVPX(0:N1) ); FCVPX = 0.
    ALLOCATE ( FCVPRX (0:N1,0:N1), FCVPRX1(0:N1,0:N1) ); FCVPRX1= 0.; FCVPRX = 0.

   call allocate_var(FC)
   
   
!!$!!!---CL-physique
   FC%CL_var_bas    = 2  ;   FC%var_bas    = 0.     ;   FC%q_var_bas    = 0.
   FC%CL_var_haut   = 2  ;   FC%var_haut   = 0.     ;   FC%q_var_haut   = 0.
   FC%CL_var_droite = 2  ;   FC%var_droite = 0.     ;   FC%q_var_droite = 0.   
   FC%CL_var_gauche = 2  ;   FC%var_gauche = 0.     ;   FC%q_var_gauche = 0.


!!$!!!---COEFFICIENT 
   do j=1, FC%NJ
      do i=1, FC%NI 

         FC%ax1(i,j) =  2._nk / (.5_nk*(dx1_L(i  ,j)+dx1_R(i  ,j))*.5_nk*(dx1_L(i,j)+dx1_R(i,j)+dx1_L(i+1,j)+dx1_R(i+1,j)))
         FC%cx1(i,j) =  2._nk / (.5_nk*(dx1_L(i+1,j)+dx1_R(i+1,j))*.5_nk*(dx1_L(i,j)+dx1_R(i,j)+dx1_L(i+1,j)+dx1_R(i+1,j)))
         FC%b1 (i,j) = -2._nk / (.5_nk*(dx1_L(i  ,j)+dx1_R(i  ,j))*.5_nk*(dx1_L(i+1,j)+dx1_R(i+1,j)))

         FC%ax2(i,j) = -2._nk / (.5_nk*(dx2_L(i,j  )+dx2_R(i,j  ))*.5_nk*(dx2_L(i,j)+dx2_R(i,j)+dx2_L(i,j+1)+dx2_R(i,j+1)))
         FC%cx2(i,j) = -2._nk / (.5_nk*(dx2_L(i,j+1)+dx2_R(i,j+1))*.5_nk*(dx2_L(i,j)+dx2_R(i,j)+dx2_L(i,j+1)+dx2_R(i,j+1)))
         FC%b2 (i,j) = -2._nk / (.5_nk*(dx2_L(i,j  )+dx2_R(i,j  ))*.5_nk*(dx2_L(i,j+1)+dx2_R(i,j+1)))
         FC%d  (i,j) =  0._nk 

      end do
   end do

    
!!$!!---CL
   if ( FC%CL_var_bas == 2 ) then
         FC%b2(:,1) =  FC%b2(:,1) -  FC%ax2(:,1)       
         FC%ax2(:,1) =  0.0_nk 
   end if

   if ( FC%CL_var_haut == 2 ) then                 
         FC%b2(:,N2-1) =  FC%b2(:,N2-1) -  FC%cx2(:,N2-1)       
         FC%cx2(:,N2-1) =  0.0_nk
   end if


   if ( FC%CL_var_gauche == 2 ) then 
         FC%b1(1,:)   =  FC%b1(1,:)  -  FC%ax1(1,:)       
         FC%ax1(1,:) =  0.0_nk
   end if 
   


   if ( FC%CL_var_droite == 2 ) then  
         FC%b1(N1-1,:)   =  FC%b1(N1-1,:) +  FC%cx1(N1-1,:)       
         FC%cx1(N1-1,:) =  0.0_nk  
   end if

   CALL EIGP((N1-1),(N1-1),FCVPX,FCVPRX,FCVPRX1,FC%ax1(:,:),FC%cx1(:,:),FC%b1(:,:)) 


!!$!!!---VORTICITE

         do j=1,N2-1
            do i=1,N1   
               var_intf_u1(i,j) = ( var%W(1)%var(i,j) +var%W(1)%var(i,j+1) )/2._nk
            end do
         end do
     
         do j=1,N2
            do i=1,N1-1
               var_intf_u2(i,j) = ( var%W(2)%var(i,j) + var%W(2)%var(i+1,j) )/2._nk
            end do
         end do
     
         do j= 1, N2-1
            do i=1, N1-1
               FC%d(i,j) =                                                 - ( &
                    (var_intf_u1(i+1,j) - var_intf_u1(i,j) ) / dx1_R(i  ,j  )  &
                  + (var_intf_u2(i,j+1) - var_intf_u2(i,j) ) / dx2_R(i  ,j  )  )    
          end do 
         end do
   
!!$!!---solver

   call RESOLP(FC%NI, FC%NJ, FC%var(:,:), FC%d(:,:),FC%ax2(:,:), &
               FC%b2(:,:),FC%cx2(:,:),1,FC%NI,1,FC%NJ,FCVPX,FCVPRX,FCVPRX1)

   somme = 0.
   RH    = 0.

   do j=1, FC%NJ
      do i=1, FC%NI                                    
         
         RH = abs(                               &
              FC%ax1(i,j)*FC%var(i-1,j)             +&
              FC%cx1(i,j)*FC%var(i+1,j)             +&
              FC%ax2(i,j)*FC%var(i,j-1)             +&
              FC%cx2(i,j)*FC%var(i,j+1)             +&
              FC%b(i,j)  *FC%var(i,j)               -&
              FC%d(i,j)                           )
         somme = somme + RH*RH        
      end do 
   end do
   print*, SQRT(somme) 
    

!!$!!---export
!!$!!!---VORTICITE

   do j=1,N2-1
      do i=1,N1   
         var_intf_u1(i,j) = ( var%W(1)%var(i,j) +var%W(1)%var(i,j+1) )/2._nk
      end do
   end do

   do j=1,N2
      do i=1,N1-1
         var_intf_u2(i,j) = ( var%W(2)%var(i,j) + var%W(2)%var(i+1,j) )/2._nk
      end do
   end do

   do j= 1, N2-1
      do i=1, N1-1
         FC%d(i,j) =                                                 - ( &
              (var_intf_u1(i+1,j) - var_intf_u1(i,j) ) / dx1_R(i  ,j  )  &
            + (var_intf_u2(i,j+1) - var_intf_u2(i,j) ) / dx2_R(i  ,j  )  )    
    end do 
   end do
   



   write( cTemp,'(i9)') it
   open (81,file=trim(adjustl(cTemp))//"_vortic")       
   yi = 0.
   do j = 1, N2 
     xi = 0.
     do i = 1, N1
         write(81,*) xi, yi, FC%var(i,j), FC%d(i,j)
         xi = xi + dx1_R(i,j)
      enddo
      write(81,*)"  "
      yi = yi + dx2_R(i,j)
   enddo
   close(81)

   open (81,file="champs_vortic")
   yi = 0.
   do j = 1, N2 
     xi = 0.
     do i = 1, N1
         
         write(81,*) xi, yi, FC%var(i,j), FC%d(i,j)
         xi = xi + dx1_R(i,j)
      enddo
      write(81,*)"  "
      yi = yi + dx2_R(i,j)
   enddo
   close(81)

   call deallocate_var(FC) 

     DEALLOCATE ( FCVPX )
     DEALLOCATE ( FCVPRX, FCVPRX1 )

#elif (DIMENSION_GEO == 3)
#endif

 end subroutine fonction_courant

end module mChampsprof






















program FiniteVolume

!  This program solves a FiniteVolume scalar problem 
   implicit none

   interface
      function flux(u_)
      real, dimension(2), intent (in) :: u_
      real, dimension(2)              :: flux
      end function flux
   end interface
   
   interface
      function U_init(x_)
      real, intent (in) :: x_
      real, dimension(2)              :: U_init
      end function U_init
   end interface

   interface
      function godunov(u_,v_)
      real, dimension(2), intent (in) :: u_,v_
      real, dimension(2)              :: godunov
      end function godunov
   end interface

   interface
      function correct_Flux(u_,v_,dt,dx)
      real, dimension(2), intent (in) :: u_,v_
      real     :: dt,dx
      real, dimension(2)              :: correct_Flux
      end function correct_Flux
   end interface

! loop int   
   integer:: i,j,k
   character(len=1) :: c1

!  constants
   real, parameter      :: pi = acos(-1.)

!  domaine spatial
   integer              :: nx
   real                 :: dx, xd, xf
   real, dimension(:),     Pointer :: X
!  domaine temporelle
   integer, parameter   :: nt=500
   real                 :: T
   real                 :: t_=   0.0,  dt=0.0

!  Flux Numeriques
   real, dimension(:,:),   Pointer :: F,F_c

!  Solution scalaire
   real, dimension(:,:),   Pointer :: U
   real,dimension (:,:),   Pointer :: sol

!  diverses valeurs numériques nécessaire
   real     :: vitesse, cfl, u_hat
   integer  :: n =0

!  error variables
   real     :: err_L1=0, err_L2=0, err_Li=0


!  file parameter
   integer, parameter   :: numfile_sol=1, numfile_data=2, numfile_param=3
   integer              :: n_imp=0, n_imp_max
   real                 :: t_imp
   character(len=32)    :: nomfile_sol = 'file_sol.txt',   nomfile_data = 'file_data.txt', nomfile_param = 'param.txt'
   character(len=32)    :: save_format
! =======================================================================================
! =======================================================================================
! =======================================================================================

   open(unit=numfile_param, file=nomfile_param, form ='formatted', status ='old')

   read(numfile_param,  *) nx;  
   read(numfile_param,  *) xd;  
   read(numfile_param,  *) xf;      
   read(numfile_param,  *) T;     
   read(numfile_param,  *) cfl;   
   read(numfile_param,  *) n_imp_max;     

   dx = real((xf-xd))/nx

   allocate(X(0:nx-1));   X(0:nx-1) = (/  (xd+ (i+0.5)*dx, i = 0,nx-1)  /)
   allocate(F(2,0:nx), F_c(2,0:nx))
   allocate(U(2,0:nx-1))
   allocate(sol(3,nx))

   t_imp=T/real(n_imp_max)

   do i=0,nx-1;   U(:,i) = U_init(X(i));   end do

   print *, "init"
   
   write(save_format, '("(" i5 "(f10.6, f12.8, f12.8 /))")') nx 

   open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')   
   write(unit= numfile_data, fmt='("nt = "i5)') nt
   write(unit= numfile_data, fmt='("nx = "i5)') nx
   write(unit= numfile_data, fmt='("save_max =" i5)')  int(T/t_imp)
   
   open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')

!  boucle while sur le temps 
   do while (t_<real(T))
      n = n +1

      vitesse =0
      do i=0,nx-2
         
      ! u_hat = (sqrt(U(1,i  ))/(sqrt(U(1,i))+ sqrt(U(1,i+1))))*(U(2,i  )/U(1,i  )) &
      !       + (sqrt(U(1,i+1))/(sqrt(U(1,i))+ sqrt(U(1,i+1))))*(U(2,i+1)/U(1,i+1)) 

      !    if(  abs(u_hat) >vitesse )   vitesse = abs(u_hat)
         if(abs(U(2,i)/U(1,i)) >vitesse )   vitesse = abs(U(2,i)/U(1,i))
      end do

      if(vitesse >1e-20) then
         dt = min(cfl * dx /(2* vitesse), T-t_)
      else
         print *,'hey'
         dt = cfl*dx 
      end if
      

      do i=0,nx
         if(U(1,i)<1e-6) U(1,i) = 1e-6  
         if(i==0) then
            F(:,0)  = godunov(U(:,0),     U(:,0)) 
            ! F_c(:,0)= correct_Flux(U(:,0),U(:,0), dt,dx)

         else if(i==nx) then
            F(:,nx) = godunov(U(:,nx-1),  U(:,nx-1)) 
            ! F_c(:,nx)= correct_Flux(U(:,nx-1),U(:,nx-1), dt,dx)

         else if(i /= 0 .and. i/=nx) then 
            F(:,i)  = godunov(U(:,i-1),   U(:,i))
            ! F_c(:,i)= correct_Flux(U(:,i-1),U(:,i), dt,dx)

         end if
         
        
      end do

      U(:,:) = U(:,:)- ((dt/dx)* (F  (:,1:nx)-F  (:,0:nx-1))) 
                     ! - ((dt/dx) *(F_c(:,1:nx)-F_c(:,0:nx-1)))

      if(t_ >=  n_imp*t_imp)  then
         print *, "loop : ",n,", n_imp",n_imp,", time :",t_," ; ","dt : ",dt, ";"
         n_imp = n_imp +1
         sol(1,:)=X(0:nx-1)
         sol(2,:)=U(1,0:nx-1)
         sol(3,:)=U(2,0:nx-1)
         write(unit=numfile_sol,  fmt=save_format) sol
         write(unit=numfile_data, fmt='("time_save =" f10.6)')  t_
      end if
      
      t_ = t_ +dt


   end do
   
   close(unit=numfile_data)
   close(unit=numfile_param)
   close(unit=numfile_sol)

   print *, "program complete !"

end program FiniteVolume
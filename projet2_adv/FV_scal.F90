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
   function U_exa(x,t,wl,wr)
      implicit none
      real, intent(in)    :: x,t,wl,wr
      real, dimension(2)  :: U_exa
      end function U_exa
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
   real                 :: t_=0.0,  dt=0.0

!  Flux Numeriques
   real, dimension(:,:),   Pointer :: F,F_c

!  Solution scalaire
   real, dimension(:,:),   Pointer :: U, U_ex
   real,dimension (:,:),   Pointer :: sol

!  diverses valeurs numériques nécessaire
   real     :: vitesse, cfl, u_hat
   integer  :: n =0

!  error variables
   real     :: err_L1=0, err_L2=0, err_Li=0
   real     :: wl=1.0 ,wr=4.0


!  file parameter
   integer, parameter   :: numfile_sol=1, numfile_data=2, numfile_param=3, numfile_err =4, numfile_conv=5
   integer              :: n_imp=0, n_imp_max
   real                 :: t_imp

   character(len=32)    :: nomfile_sol = 'file_sol.txt',   nomfile_data = 'file_data.txt', nomfile_param = 'param.txt', nomfile_err = 'error_file.txt', &
                           nomfile_conv = 'convergence_err.txt'
   character(len=128)    :: save_format
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

   allocate(F(2,0:nx),     F_c(2,0:nx))
   allocate(U(2,0:nx-1),  U_ex(2,0:nx-1))

   allocate(sol(5,nx))

   t_imp=T/real(n_imp_max)

   do i=0,nx-1;   U(:,i) = U_init(X(i));   end do

   print *, "init"
   
   write(save_format, '("(" i5 "(f10.6, f12.8, f12.8, f12.8, f12.8 /))")') nx 

   open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')   
   write(unit= numfile_data, fmt='("nt = "i5)') nt
   write(unit= numfile_data, fmt='("nx = "i5)') nx
   write(unit= numfile_data, fmt='("save_max =" i5)')  int(T/t_imp)
   
   open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')
   open(unit=numfile_err,  file=nomfile_err, form ='formatted', status ='old')

!  boucle while sur le temps 
   do while (t_<real(T))
      n = n +1

      vitesse =0
      do i=0,nx-2
         
      u_hat = (sqrt(U(1,i  ))/(sqrt(U(1,i))+ sqrt(U(1,i+1))))*(U(2,i  )/U(1,i  )) &
            + (sqrt(U(1,i+1))/(sqrt(U(1,i))+ sqrt(U(1,i+1))))*(U(2,i+1)/U(1,i+1)) 

         if(  abs(u_hat) >vitesse )   vitesse = abs(u_hat)
         ! if(abs(U(2,i)/U(1,i)) >vitesse )   vitesse = abs(U(2,i)/U(1,i))
      end do

      if(vitesse >1e-20) then
         dt = min(cfl * dx /(2* vitesse), T-t_)
      else
         print *,'hey'
         dt = cfl*dx 
      end if
      ! print *,"TIME =",t_, 'dt=',dt, 'dx=',dx

      call correct_Flux(U,dt,dx,nx,F) 

      if(t_ >=  n_imp*t_imp)  then
         print *, "loop : ",n,", n_imp",n_imp,", time :",t_," ; ","dt : ",dt, ";"

         
         do i=0,nx-1
            U_ex(:,i) = U_exa(X(i),t_-1.,wl,wr) 
         end do

         n_imp = n_imp +1
         sol(1,:)=X(0:nx-1)
         sol(2,:)=U(1,0:nx-1);         sol(3,:)=U(2,0:nx-1)
         sol(4,:)=U_ex(1,0:nx-1);      sol(5,:)=U_ex(2,0:nx-1)
         write(unit=numfile_sol,  fmt=save_format) sol

         if(sum( abs(U-U_ex))*dx > err_L1) err_L1 = sum( abs(U-U_ex))*dx 
         if(sum( (U-U_ex)**2)*dx > err_L2) err_L2 = sum( (U-U_ex)**2)*dx 

         
         write(unit=numfile_err, fmt='(" --------------- at time : "f10.6" ----------------- ")') t_ 
         write(unit=numfile_err, fmt='("err_L1 :" f16.10 )') sum( abs(U-U_ex))*dx 
         write(unit=numfile_err, fmt='("err_L2 :" f16.10 )') sum( (U-U_ex)**2)*dx   
         write(unit=numfile_data, fmt='("time_save =" f10.6)')  t_
      end if
      
      t_ = t_ +dt


   end do
   
   close(unit=numfile_data)
   close(unit=numfile_err)
   close(unit=numfile_param)
   close(unit=numfile_sol)

   open(unit=numfile_conv,  file=nomfile_conv, form ='formatted', status ='old', position='append')
   write(unit=numfile_conv, fmt='("=====================")') 
   write(unit=numfile_conv, fmt='("for nx = "i5" we have error :")' ) nx
   write(unit=numfile_conv, fmt='("err_L1 :" f16.10 )') sum( abs(U-U_ex))*dx 
   write(unit=numfile_conv, fmt='("err_L2 :" f16.10 )') sum( (U-U_ex)**2)*dx  
   write(unit=numfile_conv, fmt='("=====================")') 

   close(unit=numfile_conv)

   print *, "program complete !"

end program FiniteVolume
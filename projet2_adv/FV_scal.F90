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
      subroutine  Q_init(X,W,arg_string,arg_real)
      real, dimension(:),   intent(in)    :: X
      real, dimension(:),   intent(in), optional :: arg_real    ! (nb valeurs, nb cell, val pb)
      character(len=32),    intent(in), optional :: arg_string  ! (nom pb)
      real, dimension(:,:), intent(inout)  :: W
      end subroutine Q_init
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
   real                 :: t_=-1.0,  dt=0.0

!  Flux Numeriques
   real, dimension(:,:),   Pointer :: F,F_c

!  Solution scalaire
   real, dimension(:,:),   Pointer :: Q, Q_ex
   real,dimension (:,:),   Pointer :: sol

!  diverses valeurs numériques nécessaire
   real     :: vitesse, cfl, u_hat
   integer  :: n =0

!  error variables
   real     :: err_L1=0, err_L2=0, err_Li=0


!  file parameter
   integer, parameter   :: numfile_sol=1, numfile_data=2, numfile_param=3, numfile_err =4, numfile_conv=5
   integer              :: n_imp_max, n_imp=1
   real, dimension(:),  Pointer   :: t_imp
   real, dimension(8)   :: t_imp_ref= (/-1.0,   -0.5, 0.0,  0.5,  1.0,  1.5,  3.5,  6.0 /)

   character(len=32)    :: nomfile_sol = 'file_sol.txt',   nomfile_data = 'file_data.txt', nomfile_param = 'param.txt', nomfile_err = 'error_file.txt', &
                           nomfile_conv = 'convergence_err.txt'
   character(len=128)    :: save_format

!  Problem def
   integer :: nb_val =2
   real  :: wl=1.0 ,wr=4.0
   real  :: ul=1.0 ,ur=-2.0
   real  :: rhoL=2.0,rhoR=2.0
   character(len=32)    :: pb_name = "gaz"
! =======================================================================================
! =======================================================================================
! =======================================================================================

   ! declaration et callibrations des variables 

   open(unit=numfile_param, file=nomfile_param, form ='formatted', status ='old')

   read(numfile_param,  *) nx;  
   read(numfile_param,  *) xd;  
   read(numfile_param,  *) xf;      
   read(numfile_param,  *) T;     
   read(numfile_param,  *) cfl;   
   read(numfile_param,  *) n_imp_max;   
   read(numfile_param,  *) wL;  
   read(numfile_param,  *) wR;    
   read(numfile_param,  *) uL;  
   read(numfile_param,  *) uR;   
   read(numfile_param,  *) rhoL;  
   read(numfile_param,  *) rhoR;    
   

   dx = real((xf-xd))/nx

   allocate(X(0:nx-1));   X(0:nx-1) = (/  (xd+ (i+0.5)*dx, i = 0,nx-1)  /)

   allocate(F(2,0:nx),     F_c(2,0:nx))
   allocate(Q(2,0:nx-1),  Q_ex(2,0:nx-1))

   allocate(sol(5,nx))

   call Q_init(X,Q, arg_string= pb_name, arg_real= (/ real(nb_val), real(nx), wl, wr, rhoL, uL, rhoR, uR /) )
   
   if(n_imp_max<=8) then

      do i=1,8
         if(t_imp_ref(i)>T) then
            n_imp_max=i-1
            exit
         end if
      end do

   

      print *, n_imp_max
      allocate(t_imp(n_imp_max))
      do i=1,n_imp_max; t_imp(i)= t_imp_ref(i); end do

   end if

   print *, n_imp_max

   ! preparations des sorties

   print *, "init"
   
   write(save_format, '("(" i5 "(f10.6, f12.8, f12.8, f12.8, f12.8 /))")') nx 

   open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')   
   write(unit= numfile_data, fmt='("nt = "i5)') nt
   write(unit= numfile_data, fmt='("nx = "i5)') nx
   write(unit= numfile_data, fmt='("save_max =" i5)')  n_imp_max
   
   open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')
   open(unit=numfile_err,  file=nomfile_err, form ='formatted', status ='old')

!  boucle while sur le temps 
   do while (t_<real(T))
      n = n +1

      vitesse =0
      do i=0,nx-2
         
      u_hat = (sqrt(Q(1,i  ))/(sqrt(Q(1,i))+ sqrt(Q(1,i+1))))*(Q(2,i  )/Q(1,i  )) &
            + (sqrt(Q(1,i+1))/(sqrt(Q(1,i))+ sqrt(Q(1,i+1))))*(Q(2,i+1)/Q(1,i+1)) 

         if(  abs(u_hat) >vitesse )   vitesse = abs(u_hat)
         ! if(abs(Q(2,i)/Q(1,i)) >vitesse )   vitesse = abs(Q(2,i)/Q(1,i))
      end do

      if(vitesse >1e-20) then
         dt = min(cfl * dx /(2* vitesse), abs(T-t_))
      else
         print *,'hey', dt
         dt = cfl*dx 
      end if
      
      t_ = t_ +dt

      call Update_LeVeque(Q,dt,dx,nx,F) 

      if(t_ >=  t_imp(n_imp) ) then
         print *, "loop : ",n,", n_imp",n_imp,", time :",t_," ; ","dt : ",dt, ";"

         ! do i=0,nx-1
         !    Q_ex(:,i) = U_exa(X(i),t_,wl,wr) 
         ! end do

         n_imp = n_imp +1


         sol(1,:)=X(0:nx-1)
         sol(2,:)=Q(1,0:nx-1);         sol(3,:)=Q(2,0:nx-1)
         sol(4,:)=Q_ex(1,0:nx-1);      sol(5,:)=Q_ex(2,0:nx-1)
         write(unit=numfile_sol,  fmt=save_format) sol

         if(sum( abs(Q-Q_ex))*dx > err_L1) err_L1 = sum( abs(Q-Q_ex))*dx 
         if(sum( (Q-Q_ex)**2)*dx > err_L2) err_L2 = sum( (Q-Q_ex)**2)*dx 

         write(unit=numfile_err, fmt='(" --------------- at time : "f10.6" ----------------- ")') t_ 
         write(unit=numfile_err, fmt='("err_L1 :" f16.10 )') sum( abs(Q-Q_ex))*dx 
         write(unit=numfile_err, fmt='("err_L2 :" f16.10 )') sum( (Q-Q_ex)**2)*dx   
         write(unit=numfile_data, fmt='("time_save =" f10.6)')  t_
      end if
      
   end do
   
   close(unit=numfile_data)
   close(unit=numfile_err)
   close(unit=numfile_param)
   close(unit=numfile_sol)

   open(unit=numfile_conv,  file=nomfile_conv, form ='formatted', status ='old', position='append')
   write(unit=numfile_conv, fmt='("=====================")') 
   write(unit=numfile_conv, fmt='("for nx = "i5" we have error :")' ) nx
   write(unit=numfile_conv, fmt='("err_L1 :" f16.10 )') sum( abs(Q-Q_ex))*dx 
   write(unit=numfile_conv, fmt='("err_L2 :" f16.10 )') sum( (Q-Q_ex)**2)*dx  
   write(unit=numfile_conv, fmt='("=====================")') 

   close(unit=numfile_conv)

   print *, "program complete !"

end program FiniteVolume
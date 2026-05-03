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
      subroutine  Q_init(X,Q,arg_string,arg_real)
      real, dimension(:),   intent(in)    :: X
      real, dimension(:),   intent(in), optional :: arg_real    ! (nb valeurs, nb cell, val pb)
      character(len=32),    intent(in), optional :: arg_string  ! (nom pb)
      real, dimension(:,:), intent(inout)  :: Q
      end subroutine Q_init

      
      subroutine Update_LeVeque(Q,dt,dx,nx) 
         implicit none
         integer,intent(in)  :: nx
         real, dimension(2,1:nx), intent(inout)  :: Q
         real, intent(in) ::  dt,dx
      end subroutine Update_LeVeque
      
      subroutine Q_exa(X,t,Q_ex,arg_string,arg_real)

         implicit none
         real, dimension(:),  intent(in)   :: X
         real,                intent(in)   :: t
         real, dimension(:,:),intent(in), Pointer :: Q_ex
         real, dimension(:),  intent(in), optional :: arg_real    ! (nb valeurs, nb cell, val pb)
         character(len=32),   intent(in), optional :: arg_string  ! (nom pb)
      end subroutine Q_exa

   
      subroutine Update(Q,X,dt,nx, arg_string)
         implicit none
         integer,                intent(in)  :: nx
         real, dimension(2,0:nx-1),intent(inout)  :: Q
         real, dimension(0:nx+1),intent(in)  :: X

         real,                   intent(in) ::  dt
         character(len=32),      intent(in),optional :: arg_string
      end subroutine
   end interface
   



! loop int   
   integer:: i,j,k
   character(len=1) :: blank

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

!  Solution scalaire
   real, dimension(:,:),   Pointer :: Q, Q_ex
   real,dimension (:,:),   Pointer :: sol
   integer :: is_allocated

!  diverses valeurs numériques nécessaire
   real     :: vitesse, cfl, u_hat
   integer  :: n =0

!  error variables
   real     :: err_L1=0, err_L2=0, err_Li=0

!  probleme variable
   character(len=32)    :: methode_update

!  file parameter
   integer, parameter   :: numfile_sol=1, numfile_data=2, numfile_param=3, numfile_err =4, numfile_conv=5
   integer              :: n_imp_max, n_imp=1
   ! real, dimension(:),  Pointer   :: t_imp
   real, dimension(8)   :: t_imp= (/-1.0,   -0.5, 0.0,  0.5,  1.0,  1.5,  3.5,  6.0 /)

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

      print *, pb_name

      open(unit=numfile_param, file=nomfile_param, form ='formatted', status ='old')

      read(numfile_param,  *) nx;  
      read(numfile_param,  *) xd;  
      read(numfile_param,  *) xf;      
      read(numfile_param,  *) T;     
      read(numfile_param,  *) cfl;   
      read(numfile_param,  *) n_imp_max;   
      read(numfile_param,  *) blank;
      read(numfile_param,  *) methode_update;      
      read(numfile_param,  *) wL;  
      read(numfile_param,  *) wR;    
      read(numfile_param,  *) uL;  
      read(numfile_param,  *) uR;   
      read(numfile_param,  *) rhoL;  
      read(numfile_param,  *) rhoR;    
   
   print *, methode_update

   dx = real((xf-xd))/nx

   ! allocation des pointeurs 

      allocate(X(0:nx),  stat=is_allocated);   X(1:nx) = (/  (xd+ i*dx, i = 1,nx-1)  /); X(0)=xd; X(nx)=xf
      print *, "is X alloc :",is_allocated
      allocate(Q(2,nx),  stat=is_allocated);    print *, "is Q alloc :",is_allocated
      allocate(Q_ex(2,nx),stat=is_allocated);   print *, "is Qex alloc :",is_allocated
      allocate(sol(5,nx),stat=is_allocated);    print *, "is sol alloc :",is_allocated

   ! initialisation des solutions numérique et exactes
      call Q_init(X,Q,     arg_string= pb_name, arg_real= (/ real(nb_val), real(nx), wl, wr, rhoL, uL, rhoR, uR /))
      call Q_exa(X,t_,Q_ex,arg_string= pb_name, arg_real= (/ real(nb_val), real(nx), wl, wr, rhoL, uL, rhoR, uR /))

   ! preparations des sorties
   
      ! format de sauvegarde ->file_sol.txt
      write(save_format, '("(" i5 "(f10.6, f12.8, f12.8, f12.8, f12.8 /))")') nx 

      open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')   
      write(unit= numfile_data, fmt='("nt = "i5)') nt
      write(unit= numfile_data, fmt='("nx = "i5)') nx
      write(unit= numfile_data, fmt='("save_max =" i5)')  n_imp_max
      
      open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')
      open(unit=numfile_err,  file=nomfile_err, form ='formatted', status ='old')

!  boucle while sur le temps 

   print *,"init"

   do while (t_<real(T))
      n = n +1

      if (.not. associated(Q)) print *, "the Sol is desassociated " 

      ! critère CFL 
         vitesse =0
         do i=0,nx-2
            u_hat = (sqrt(Q(1,i  ))/(sqrt(Q(1,i))+ sqrt(Q(1,i+1))))*(Q(2,i  )/Q(1,i  )) &
                  + (sqrt(Q(1,i+1))/(sqrt(Q(1,i))+ sqrt(Q(1,i+1))))*(Q(2,i+1)/Q(1,i+1)) 

            if(  abs(u_hat) >vitesse )   vitesse = abs(u_hat)
         end do

         if(vitesse >1e-20) then
            dt = min(cfl * dx /(2* vitesse), abs(T-t_))
         else
            dt = cfl*dx/2 
            print *,'pas de temps trop petit, forcage du pas de temps', dt
         end if

   ! sauvegarde 
      if(t_ >=  t_imp(n_imp) ) then
         print *, "loop : ",n,", n_imp",n_imp,", time :",t_," ; ","dt : ",dt, ";"

         ! calcul de la solution exacte
         call Q_exa(X,t_,Q_ex)

         n_imp = n_imp +1

         ! sauvegarde dans le fichier sol
            sol(1,:)=X(1:nx);
            sol(2,1:nx)=Q(1,1:nx);         
            sol(3,1:nx)=Q(2,1:nx);
            sol(4,:)=Q_ex(1,1:nx);      
            sol(5,:)=Q_ex(2,1:nx);
            write(unit=numfile_sol,  fmt=save_format) sol
         ! calcul des erreurs numérique
            if(sum( abs(Q-Q_ex))*dx > err_L1) then
               do i=1,nx ; if(abs(Q(1,i))<= max(rhoL, rhoR)) err_L1 = err_L1+ ( (abs(Q(1,i)-Q_ex(1,i))))*dx ;   end do
            end if

            if(maxval(abs(Q-Q_ex))>err_Li) then
               do i=1,nx ; if(abs(Q(1,i))<= max(rhoL, rhoR))  err_Li = max(err_Li, abs(Q(1,i)-Q_ex(1,i))); end do
            end if

            if(sum( (Q-Q_ex)**2)*dx > err_L2) then
               do i=1,nx ; if(abs(Q(1,i))<= max(rhoL, rhoR))  err_L2 =err_L2 + ( (abs((Q(1,i)-Q_ex(1,i))))**2)*dx ;    end do
            end if
         ! sauvegarde des erreurs numériques
            write(unit=numfile_err, fmt='(" --------------- at time : "f10.6" ----------------- ")') t_ 
            write(unit=numfile_err, fmt='("err_L1 :" f16.10 )') err_L1
            write(unit=numfile_err, fmt='("err_L2 :" f16.10 )') err_L2
            write(unit=numfile_err, fmt='("err_Li :" f16.10 )') err_Li
            write(unit=numfile_data, fmt='("time_save =" f10.6)')  t_
      end if
      
      ! calcul de la solution numérique au pas de temps suivant
      t_ = t_ +dt
      call Update(Q=Q, X=X,dt=dt,nx =nx, arg_string = methode_update) 
      
   end do
   
   close(unit=numfile_data)
   close(unit=numfile_err)
   close(unit=numfile_param)
   close(unit=numfile_sol)

   ! sauvegarde de l'erreur numérique dans le fichier de convergence
   open(unit=numfile_conv,  file=nomfile_conv, form ='formatted', status ='old', position='append')
   write(unit=numfile_conv, fmt='("=====================")') 
   write(unit=numfile_conv, fmt='("for nx = "i5" we have error :")' ) nx
   write(unit=numfile_conv, fmt='("err_L1 :" f16.10 )') err_L1
   write(unit=numfile_conv, fmt='("err_L2 :" f16.10 )') err_L2
   write(unit=numfile_conv, fmt='("err_Li :" f16.10 )') err_Li
   write(unit=numfile_conv, fmt='("=====================")') 
   close(unit=numfile_conv)

   print *, "program complete !"

end program FiniteVolume
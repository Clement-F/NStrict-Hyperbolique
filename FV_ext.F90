
function flux(Q) result(f)
   implicit none
   real, dimension(2), intent(in)  :: Q
   real, dimension(2)              :: f

   f(1) = Q(2) 
   f(2) = (Q(2)/Q(1))*Q(2)
   return
end function flux

function godunov(Ql,Qr)
   implicit none
   real, dimension(2), intent(in)  :: Ql,Qr
   real, dimension(2)              :: godunov

   interface
   function flux(Q) result(f)
      implicit none
      real, dimension(2), intent(in)  :: Q
      real, dimension(2)              :: f
   end function flux
   end interface

   real              :: u_hat,uL,uR,rhoL,rhoR

   rhoL = Ql(1); rhoR = Qr(1)
   uL = Ql(2)/Ql(1); uR = Qr(2)/Qr(1)

   if(uL<0 .and. 0<uR) then
      godunov = 0.
   else 
      u_hat = uL*(sqrt(rhoL)/(sqrt(rhoL) + sqrt(rhoR))) &
            + uR*(sqrt(rhoR)/(sqrt(rhoL) + sqrt(rhoR)))

      if(u_hat>0) then 
         godunov = flux(Ql)

      else if(u_hat<0) then
         godunov = flux(Qr)

      else if(u_hat==0) then
         godunov = (flux(Ql)+flux(Qr))*0.5

      else 
         print *, "ERROR", Ql, Qr
         godunov = 0.

      end if
   end if

   ! print *, u_,v_,godunov
   return
end function godunov

subroutine Update_LeVeque(Q,dt,dx,nx) 
   implicit none
   integer,intent(in)  :: nx
   real, dimension(2,1:nx), intent(inout)  :: Q

   real, intent(in) ::  dt,dx
   
   interface
      function godunov(u_,v_)
         implicit none
         real, dimension(2), intent(in)  :: u_,v_
         real, dimension(2)              :: godunov
      end function godunov
      
      function flux(Q) result(f)
         implicit none
         real, dimension(2), intent(in)  :: Q
         real, dimension(2)              :: f
      end function flux
   end interface

   real   :: s1,s2
   real, dimension(2,nx) :: F
   real, dimension(2,nx) :: z1,z2, z1c,z2c
   real, dimension(2,nx)   :: s
   real, dimension(2)   :: F_=0, F_pred=0
   real                 :: u_hat, theta1, theta2, phi1,phi2
   real                 :: rhoL,rhoR,uL,uR
   integer :: i


   do i=1,nx
      if(Q(1,i)<1e-20) Q(1,i) = 1e-20 

      if(i==1) then
         F(:,1)  = godunov(Q(:,1),     Q(:,1)) 

      else if(i==nx) then
         F(:,nx) = godunov(Q(:,nx),  Q(:,nx)) 

      else if(i /= 1 .and. i/=nx) then 
         F(:,i)  = godunov(Q(:,i),   Q(:,i+1))

      end if
   end do

   ! print *, "godunov calc"
   
   z1=0; z2=0; s=0

   do i=1,nx-1
      rhoL=Q(1,i);         rhoR  =Q(1,i+1)
      uL  =Q(2,i)/Q(1,i);  uR    =Q(2,i+1)/Q(1,i+1)


      if(uL<0 .and. 0<uR) then
         z1(:,i) = -flux(Q(:,i));  z2(:,i)= flux(Q(:,i+1))
         s(1,i) = uL;  s(2,i) =uR

      else  
         u_hat = uL*(sqrt(rhoL)/(sqrt(rhoL) + sqrt(rhoR))) &
               + uR*(sqrt(rhoR)/(sqrt(rhoL) + sqrt(rhoR)))

         if(u_hat<0) then
            ! z1(:,i) = flux(Q(:,i))-flux(Q(:,i-1)); z2(:,i) =0
            z1(:,i) = Q(:,i+1) - Q(:,i); z2(:,i) =0

            s(1,i) = u_hat; s(2,i) = u_hat

         else 
            ! z1(:,i) = 0; z2(:,i) = flux(Q(:,i))-flux(Q(:,i-1))
            z1(:,i) = 0; z2(:,i) = Q(:,i+1)-Q(:,i)

            s(1,i) = u_hat; s(2,i) = u_hat

         end if

      end if      
      if(Q(1,i)<0.) call exit(0)

   end do

   do i=2,nx-1
      
   theta1 =0; theta2=0; phi1=0; phi2=0

      if(z1(1,i)/=z1(1,i)) print *,"z1",i
      if(z2(1,i)/=z2(1,i)) print *,"z2",i
      if(s(1,i) /=s(1,i) ) print *,"s",i
      if(F_(1)  /=F_(1)  ) print *,"F",i
      if(Q(1,i) /=Q(1,i) ) print *,"Q",i
      if(phi1 /= phi1)  print *,'Phi1',phi1
      if(phi2 /= phi2)  print *,'Phi2',phi2

      if(s(1,i)>0 .or. s(2,i)>0) then 
         ! print *,'upwind'
         if(abs(z1(1,i))> 1e-10 .or. abs(z1(2,i)) > 1e-10) theta1 = (z1(1,i)*z1(1,i-1) + z1(2,i)*z1(2,i-1))/(z1(1,i)*z1(1,i) + z1(2,i)*z1(2,i))
         if(abs(z2(1,i))> 1e-10 .or. abs(z2(2,i)) > 1e-10) theta2 = (z2(1,i)*z2(1,i-1) + z2(2,i)*z2(2,i-1))/(z2(1,i)*z2(1,i) + z2(2,i)*z2(2,i))

      else 
         ! print *,'downwind'
         if(abs(z1(1,i))> 1e-10 .or. abs(z1(2,i)) > 1e-10) theta1 = (z1(1,i)*z1(1,i+1) + z1(2,i)*z1(2,i+1))/(z1(1,i)*z1(1,i) + z1(2,i)*z1(2,i))
         if(abs(z2(1,i))> 1e-10 .or. abs(z2(2,i)) > 1e-10) theta2 = (z2(1,i)*z2(1,i+1) + z2(2,i)*z2(2,i+1))/(z2(1,i)*z2(1,i) + z2(2,i)*z2(2,i))
      
      end if

      if(theta1 /=theta1 ) print *,i,"theta1",z1(:,i),z1(:,i-1)
      if(theta2 /=theta2 ) print *,i,"theta2",z2(:,i),z2(:,i-1)

      phi1 = max(0.,min(1.,(theta1)));  
      phi2 = max(0.,min(1.,(theta2)));  

      ! print *, theta1, phi1
      ! print *, theta2, phi2

      ! F_ = 0.5*(    sign(1.,s(1,i)) * (1-(dt/dx) *abs(s(2,i))))*(z1(:,i)*phi1)
      ! F_ = F_+ 0.5*(sign(1.,s(2,i)) * (1-(dt/dx) *abs(s(1,i))))*(z2(:,i)*phi2)

      F_ = 0.5*(    abs(s(1,i)) * (1-(dt/dx) *abs(s(2,i))))*(z1(:,i)*phi1)
      F_ = F_+ 0.5*(abs(s(2,i)) * (1-(dt/dx) *abs(s(1,i))))*(z2(:,i)*phi2)

      ! if(F_(1) /=0 .or. F_(2) /=0) print *,"F_ :",F_
      
      if(F_(1) /= F_(1)) call exit(0)

      if( (Q(1,i)-((dt/dx)* (F(1,i)-F(1,i-1))) - ((dt/dx)* (F_(1) -F_pred(1)) ))<0. )  then
         z1c(:,i-1)=0; z2c(:,i-1)=0;
         z1c(:,i)  =0; z2c(:,i)  =0;
         F_pred =0; 
         ! print *, "nul",i
      else 
         z1c(:,i)=z1(:,i)*phi1; z2c(:,i)=z2(:,i)*phi2;
         F_pred = F_
      end if
   end do

   ! print *, "limited flux calc"

   do i=2,nx-1
      ! F_ = 0.5*(    sign(real(1.),s(1,i)) * (1-(dt/dx) *abs(s(2,i))))*z1c(:,i)
      ! F_ = F_+ 0.5*(sign(real(1.),s(2,i)) * (1-(dt/dx) *abs(s(1,i))))*z2c(:,i)

      
      F_ = 0.5*(    abs(s(1,i)) * (1-(dt/dx) *abs(s(2,i))))*(z1c(:,i)*phi1)
      F_ = F_+ 0.5*(abs(s(2,i)) * (1-(dt/dx) *abs(s(1,i))))*(z2c(:,i)*phi2)

      F(:,i) = F(:,i) + F_

      Q(:,i)= Q(:,i)- ((dt/dx)* (F(:,i)-F(:,i-1)))
   end do

   ! print *,"updated state"



end subroutine Update_LeVeque

subroutine Q_init(X,Q,arg_string,arg_real)
   implicit none
   real, dimension(:),  intent(in)    :: X
   real, dimension(:),  intent(in), optional :: arg_real    ! (nb valeurs, nb cell, val pb)
   character(len=32),   intent(in), optional :: arg_string  ! (nom pb)
   real, dimension(:,:),intent(inout) :: Q

   integer :: i=0

   ! saving info
   integer,          save :: nx
   character(len=32),save :: nom_pb
   real             ,save :: wL,wR,rhoL,rhoR,uL,uR
   real :: x_
   
   if(present(arg_real)) nx = arg_real(2)
   if(present(arg_string)) nom_pb = arg_string

   
   Q=0
   print *, "Q_init"
   if(nom_pb=="gaz") then
      
      if(present(arg_real)) then; 
         wL  = arg_real(3);wR= arg_real(4)
         rhoL= arg_real(5);uL= arg_real(6)
         rhoR= arg_real(7);uR= arg_real(8)
      end if

      do i=1,nx
         x_ = (X(i)+X(i+1))/2
         if(x_>-wL-uL      .and. x_<-uL ) then
            Q(1,i) = rhoL; Q(2,i)= rhoL*uL

         else if (x_>-uR .and. x_<-uR+wR) then
            Q(1,i) = rhoR; Q(2,i)= rhoR*uR

         else
            Q(1,i)= 1e-20; Q(2,i)=1e-20
         end if
      end do
   end if

   ! print *,Q

   return 
end subroutine Q_init

subroutine Q_exa(X,t,Q_ex,arg_string,arg_real)

   implicit none
   real, dimension(:),  intent(in)   :: X
   real,                intent(in)   :: t
   real, dimension(:),  intent(in), optional :: arg_real    ! (nb valeurs, nb cell, val pb)
   character(len=32),   intent(in), optional :: arg_string  ! (nom pb)

   real, dimension(:,:),intent(out), Pointer :: Q_ex

   real  :: T1, T2,  R
   real  :: Xhi, u_hat
   integer ::i

   integer,          save :: nx
   character(len=32),save :: nom_pb
   real             ,save :: wL,wR,rhoL,rhoR,uL,uR
   real :: x_
   
   if(present(arg_real)) nx = arg_real(2)
   if(present(arg_string)) nom_pb = arg_string

   allocate(Q_ex(2,nx))

   if(nom_pb=="gaz") then
      if(present(arg_real)) then; 
         wL  = arg_real(3);wR= arg_real(4)
         rhoL= arg_real(5);uL= arg_real(6)
         rhoR= arg_real(7);uR= arg_real(8)
      end if

      u_hat = uL*(sqrt(rhoL)/(sqrt(rhoL) + sqrt(rhoR))) &
            + uR*(sqrt(rhoR)/(sqrt(rhoL) + sqrt(rhoR)))
      R = rhoL/rhoR

      T1 = wl/(uL-u_hat) 
      T2 = (wr**2 +2*wr*wl*R + (wl**2) *R )/(2*wl*R*(uL - uR))

      Q_ex=0

      do i=1,nx
      
         x_ = (X(i)+X(i+1))/2
         if(t<0.0) then
            if((x_+1.) -(t+1) <0. .and. (x_+1.+wl) -(t+1)>0.  ) then
               Q_ex(1,i) =rhoL; Q_ex(2,i) =rhoL*ul
            else if((x_-1.) +(t+1) >0. .and. (x_-1.-wr) +(t+1)<0.  ) then
               Q_ex(1,i) =rhoR; Q_ex(2,i) =rhoR*ur
            end if
         else if(t<T1) then
            Xhi = u_hat*t
            if((x_-xhi)<0. .and. (x_+1.+wl) -(t+1)>0.) then 
               Q_ex(1,i) =rhoL; Q_ex(2,i) =rhoL*ul
            else if((x_-xhi)>0. .and. (x_-1.-wr) +(t+1)<0.) then
               Q_ex(1,i) =rhoR; Q_ex(2,i) =rhoR*ur
            end if

         else if(t<T2) then
            Xhi = uR*t -wl*R + sqrt(2*wl*R*(uL-uR)*t +wl**2 *(R**2-R))
            if((x_-xhi)<0. .and. (x_+1.+wl) -(t+1)>0.) then 
               Q_ex(:,i) = 0
            else if((x_-xhi)>0. .and. (x_-1.-wr) +(t+1)<0.) then
               Q_ex(1,i) =rhoR; Q_ex(2,i) =rhoR*ur
            end if
         end if
      end do
   end if
   return

end subroutine Q_exa


function Lax_Friedrichs(qg,qd)
   implicit none
   real, dimension(2), intent(in)  ::qg,qd
   real, dimension(2)              :: Lax_Friedrichs

   interface
   function flux(Q) result(f)
      implicit none
      real, dimension(2), intent(in)  :: Q
      real, dimension(2)              :: f
   end function flux
   end interface
   real              :: beta

   
   beta = max( abs(qg(2)/qg(1)), abs(qd(2)/qd(1)))

   Lax_Friedrichs = 0.5* (flux(qg)+flux(qd) - beta*(qd -qg))
   return 

end function Lax_Friedrichs


function minmod(x,y,z)

   implicit none
   real, intent(in)    :: x,y,z
   real                :: minmod

   minmod = min(0.,max(x,y,z)) + max(0.,min(x,y,z))
   return

end function

subroutine Update(Q,X,dt,nx, arg_string)
   implicit none
   integer,                intent(in)  :: nx
   real, dimension(2,1:nx),intent(inout)  :: Q  ! Q_i
   real, dimension(0:nx+1),intent(in)  :: X     

   real,                   intent(in) ::  dt
   character(len=32),      intent(in),optional :: arg_string

   character(len=32) :: methode_update
   integer :: i 
   real, dimension(2,0:nx) :: F           ! F_{i+1/2}
   real, dimension(:,:), Pointer :: Q_int ! Q_{t n+1/2,i}
   real, dimension(:,:), Pointer :: delta
   real                  :: dx
   real  :: minmod
   real  :: alpha 

    interface
      function godunov(qg,qd)
         real,dimension(2), intent (in) :: qg,qd
         real,dimension(2)              :: godunov
      end function godunov

      function Lax_Friedrichs(qg,qd)
         implicit none
         real,dimension(2), intent(in)  :: qg,qd
         real,dimension(2)              :: Lax_Friedrichs
      end function Lax_Friedrichs

   end interface




   dx = X(2)-X(1)

   if(present(arg_string)) methode_update = arg_string

   if(methode_update == "godunov") then
   

      do i=0,nx
         if(Q(1,i)<1e-20) Q(1,i) = 1e-20 
         if(i==0) then;          F(:,0)  = godunov(Q(:,1),     Q(:,1)) 
         else if(i==nx) then;    F(:,nx) = godunov(Q(:,nx),  Q(:,nx)) 
         else if(i /= 1 .and. i/=nx) then 
            F(:,i)  = godunov(Q(:,i),   Q(:,i+1))
         end if
      end do

      do i=1,nx
         Q(:,i)= Q(:,i)- ((dt/dx)* (F(:,i)-F(:,i-1))) 
      end do
      return 

   else if(methode_update == "LeVeque") then

      call Update_LeVeque(Q,dt,dx,nx)
      return 

   else if(methode_update == "limitation") then

      allocate(Q_int(2,1:nx))
      allocate(delta(2,1:nx))
      alpha = 0.5

      delta(1,1)  = (Q(1,2)-Q(1,1))/(2*dx);      delta(2,1)  = (Q(2,2) -Q(2,1))/(2*dx);  
      delta(1,nx) = (Q(1,nx)-Q(1,nx-1))/(2*dx);  delta(2,nx) = (Q(2,nx)-Q(2,nx-1))/(2*dx);
      do i=2,nx-1
         delta(1,i) = (Q(1,i+1)-Q(1,i-1))/(2*dx);
         delta(2,i) = (Q(2,i+1)-Q(2,i-1))/(2*dx)
      end do

      delta(1,1)    = minmod(delta(1,1),   2*alpha* (Q(1,1)-Q(1,1))/dx,      2*alpha* (Q(1,2)-Q(1,1))/dx)
      delta(2,1)    = minmod(delta(2,1),   2*alpha* (Q(2,1)-Q(2,1))/dx,      2*alpha* (Q(2,1)-Q(2,1))/dx)

      delta(1,nx)   = minmod(delta(1,nx),  2*alpha* (Q(1,nx)-Q(1,nx-1))/dx,  2*alpha* (Q(1,nx)-Q(1,nx))/dx)
      delta(2,nx)   = minmod(delta(2,nx),  2*alpha* (Q(2,nx)-Q(2,nx-1))/dx,  2*alpha* (Q(2,nx)-Q(2,nx))/dx)

      do i=2,nx-1
         delta(1,i) = minmod(delta(1,i),2*alpha* (Q(1,i)-Q(1,i-1))/dx,2*alpha* (Q(1,i+1)-Q(1,i))/dx)
         delta(2,i) = minmod(delta(2,i),2*alpha* (Q(2,i)-Q(2,i-1))/dx,2*alpha* (Q(2,i+1)-Q(2,i))/dx)
      end do

      do i=0,nx
        if(i==0) then;        F(:,0)  = Lax_Friedrichs(Q(:,1) ,Q(:,1))
        else if (i==nx) then; F(:,nx) = Lax_Friedrichs(Q(:,nx),Q(:,nx))

        else;                 F(:,i)  = Lax_Friedrichs(Q(:,i-1)+0.5*dx*delta(:,i-1)    ,Q(:,i)-0.5*dx*delta(:,i))
        end if
      end do

      do i=1,nx
         Q_int(:,i)= Q(:,i)- ((dt/(2*dx))* (F(:,i)-F(:,i-1)))
      end do 

      delta(1,1)  = (Q_int(1,2)-Q_int(1,1))/(2*dx);      delta(2,1)  = (Q_int(2,2) -Q_int(2,1))/(2*dx);  
      delta(1,nx) = (Q_int(1,nx)-Q_int(1,nx-1))/(2*dx);  delta(2,nx) = (Q_int(2,nx)-Q_int(2,nx-1))/(2*dx);
      do i=2,nx-1
         delta(1,i) = (Q_int(1,i+1)-Q_int(1,i-1))/(2*dx);
         delta(2,i) = (Q_int(2,i+1)-Q_int(2,i-1))/(2*dx)
      end do

      delta(1,1)    = minmod(delta(1,1),   2*alpha* (Q_int(1,1)-Q_int(1,1))/dx,      2*alpha* (Q_int(1,2)-Q_int(1,1))/dx)
      delta(2,1)    = minmod(delta(2,1),   2*alpha* (Q_int(2,1)-Q_int(2,1))/dx,      2*alpha* (Q_int(2,1)-Q_int(2,1))/dx)

      delta(1,nx)   = minmod(delta(1,nx),  2*alpha* (Q_int(1,nx)-Q_int(1,nx-1))/dx,  2*alpha* (Q_int(1,nx)-Q_int(1,nx))/dx)
      delta(2,nx)   = minmod(delta(2,nx),  2*alpha* (Q_int(2,nx)-Q_int(2,nx-1))/dx,  2*alpha* (Q_int(2,nx)-Q_int(2,nx))/dx)

      do i=2,nx-1
         delta(1,i) = minmod(delta(1,i),2*alpha* (Q_int(1,i)-Q_int(1,i-1))/dx,2*alpha* (Q_int(1,i+1)-Q_int(1,i))/dx)
         delta(2,i) = minmod(delta(2,i),2*alpha* (Q_int(2,i)-Q_int(2,i-1))/dx,2*alpha* (Q_int(2,i+1)-Q_int(2,i))/dx)
      end do

      do i=0,nx
        if(i==0) then;        F(:,0)  = Lax_Friedrichs(Q_int(:,1) ,Q_int(:,1))
        else if (i==nx) then; F(:,nx) = Lax_Friedrichs(Q_int(:,nx),Q_int(:,nx))

        else;                 F(:,i)  = Lax_Friedrichs(Q_int(:,i-1)+0.5*dx*delta(:,i-1)    ,Q_int(:,i)-0.5*dx*delta(:,i))
        end if
      end do

      do i=2,nx
         Q(:,i)= Q(:,i)- ((dt/dx)* (F(:,i)-F(:,i-1)))
      end do 

   
   
   end if

end subroutine
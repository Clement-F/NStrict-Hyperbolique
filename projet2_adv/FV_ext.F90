
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
      print *, "delta choc"

      else 
         print *, "ERROR", Ql, Qr

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
   end interface

   interface
   function flux(Q) result(f)
      implicit none
      real, dimension(2), intent(in)  :: Q
      real, dimension(2)              :: f
   end function flux
   end interface

   real   :: s1,s2
   real*8, dimension(2,nx) :: F
   real, dimension(2,nx) :: z1,z2,s, z1c,z2c
   real*8, dimension(2)   :: F_=0, F_pred=0
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
            z1(:,i) = flux(Q(:,i))-flux(Q(:,i-1)); z2(:,i) =0
            ! z1(:,i) = Q(:,i) - Q(:,i-1); z2(:,i) =0
            s(1,i) = u_hat; s(2,i) = u_hat

         else 
            z1(:,i) = 0; z2(:,i) = flux(Q(:,i))-flux(Q(:,i-1))
            ! z1(:,i) = 0; z2(:,i) = Q(:,i)-Q(:,i-1)
            s(1,i) = u_hat; s(2,i) = u_hat

         end if

      end if      
      if(Q(1,i)<0.) call exit(0)

      ! print *,i,rhoL,rhoR, ul,ur

   end do

   ! print *, "correct flux calc"

   do i=2,nx-1
      
   theta1 =0; theta2=0; phi1=0; phi2=0

      if(z1(1,i)/=z1(1,i)) print *,"z1",i
      if(z2(1,i)/=z2(1,i)) print *,"z2",i
      if(s(1,i) /=s(1,i) ) print *,"s",i
      if(F_(1)  /=F_(1)  ) print *,"F",i
      if(Q(1,i) /=Q(1,i) ) print *,"Q",i
      if(phi1 /= phi1)  print *,'Phi1',phi1
      if(phi2 /= phi2)  print *,'Phi2',phi2

      if(abs(z1(1,i))> 1e-10 .or. abs(z1(2,i)) > 1e-10) theta1 = (z1(1,i)*z1(1,i-1) + z1(2,i)*z1(2,i-1))/(z1(1,i)*z1(1,i) + z1(2,i)*z1(2,i))
      if(abs(z2(1,i))> 1e-10 .or. abs(z2(2,i)) > 1e-10) theta2 = (z2(1,i)*z2(1,i-1) + z2(2,i)*z2(2,i-1))/(z2(1,i)*z2(1,i) + z2(2,i)*z2(2,i))

      if(theta1 /=theta1 ) print *,i,"theta1",z1(:,i),z1(:,i-1)
      if(theta2 /=theta2 ) print *,i,"theta2",z2(:,i),z2(:,i-1)

      phi1 = max(0.,min(1.,theta1));  
      phi2 = max(0.,min(1.,theta2));  

      F_ = 0.5*(    sign(1.,s(1,i)) * (1-(dt/dx) *abs(s(2,i))))*(z1(:,i)*phi1)
      F_ = F_+ 0.5*(sign(1.,s(2,i)) * (1-(dt/dx) *abs(s(1,i))))*(z2(:,i)*phi2)

      if(F_(1) /=0 .or. F_(2) /=0) print *,"F_ :",F_
      
      if(F_(1) /= F_(1)) call exit(0)

      if( (Q(1,i)-((dt/dx)* (F(1,i)-F(1,i-1))) - ((dt/dx)* (F_(1) -F_pred(1)) ))<0. )  then
         z1c(:,i-1)=0; z2c(:,i-1)=0;
         z1c(:,i)  =0; z2c(:,i)  =0;
         F_pred =0; print *, "nul",i
      else 
         z1c(:,i)=z1(:,i)*phi1; z2c(:,i)=z2(:,i)*phi2;
         F_pred = F_
      end if
   end do

   ! print *, "limited flux calc"

   do i=2,nx-1
      F_ = 0.5*(    sign(1.,s(1,i)) * (1-(dt/dx) *abs(s(2,i))))*z1c(:,i)
      F_ = F_+ 0.5*(sign(1.,s(2,i)) * (1-(dt/dx) *abs(s(1,i))))*z2c(:,i)
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
         if(X(i)>-wL-uL      .and. X(i)<-uL ) then
            Q(1,i) = rhoL; Q(2,i)= rhoL*uL

         else if (X(i)>-uR .and. X(i)<-uR+wR) then
            Q(1,i) = rhoR; Q(2,i)= rhoR*uR

         else
            Q(1,i)= 1e-20; Q(2,i)=1e-20
         end if
      end do
   end if

   ! print *,Q

   return 
end subroutine Q_init

! subroutine Q_exa(X,t,arg_string,arg_real)

!    implicit none
!    real, dimension(:),  intent(in)   :: X
!    real,                intent(in)   :: t
!    real, dimension(:),  intent(in), optional :: arg_real    ! (nb valeurs, nb cell, val pb)
!    character(len=32),   intent(in), optional :: arg_string  ! (nom pb)

!    real, dimension(:,:), Pointer :: Q

!    real  :: T1, T2,  R
!    real  :: Xhi, u_hat
!    integer ::i

!    integer,          save :: nx
!    character(len=32),save :: nom_pb
!    real             ,save :: wL,wR,rhoL,rhoR,uL,uR

   

!    print *,"init"
   
!    if(present(arg_real)) nx = arg_real(2)
!    if(present(arg_string)) nom_pb = arg_string

!    allocate(Q(2,nx))

!    if(nom_pb=="gaz") then
!       print *,"pb gaz"
!       if(present(arg_real)) then; 
!          wL  = arg_real(3);wR= arg_real(4)
!          rhoL= arg_real(5);uL= arg_real(6)
!          rhoR= arg_real(7);uR= arg_real(8)
!       end if


!       u_hat = uL*(sqrt(rhoL)/(sqrt(rhoL) + sqrt(rhoR))) &
!             + uR*(sqrt(rhoR)/(sqrt(rhoL) + sqrt(rhoR)))
!       R = rhoL/rhoR

!       T1 = wl/(uL-u_hat) 
!       T2 = (wr**2 +2*wr*wl*R + (wl**2) *R )/(2*wl*R*(uL - uR))

!       Q=0

!       do i=1,nx
!       ! print *, "assign ",X(i)
!          if(t<0.0) then
!             if((X(i)+1.) -(t+1) <0. .and. (X(i)+1.+wl) -(t+1)>0.  ) Q =ul
!             if((X(i)-1.) +(t+1) >0. .and. (X(i)-1.-wr) +(t+1)<0.  ) Q =ur
            
!          else if(t<T1) then
!             Xhi = u_hat*t
!             if((X(i)-xhi)<0. .and. (X(i)+1.+wl) -(t+1)>0.) then 
!                Q = ul
!             else if((X(i)-xhi)>0. .and. (X(i)-1.-wr) +(t+1)<0.) then
!                Q = ur
!             end if

!          else if(t<T2) then
!             Xhi = uR*t -wl*R + sqrt(2*wl*R*(uL-uR)*t +wl**2 *(R**2-R))
!             if((X(i)-xhi)<0. .and. (X(i)+1.+wl) -(t+1)>0.) then 
!                Q = 0
!             else if((X(i)-xhi)>0. .and. (X(i)-1.-wr) +(t+1)<0.) then
!                Q = ur
!             end if
!          end if
!       end do
!    end if
!    return

! end subroutine Q_exa
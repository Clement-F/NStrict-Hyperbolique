
function flux(u) result(f)
   implicit none
   real, dimension(2), intent(in)  :: u
   real, dimension(2)              :: f

   f(1) = u(2) 
   f(2) = (u(2)/u(1))*u(2)
   return
end function flux

function godunov(u_,v_)
   implicit none
   real, dimension(2), intent(in)  :: u_,v_
   real, dimension(2)              :: godunov

   interface
   function flux(u) result(f)
      implicit none
      real, dimension(2), intent(in)  :: u
      real, dimension(2)              :: f
   end function flux
   end interface

   real              :: u_hat

   if((u_(2)/u_(1))<0 .and. 0<(v_(2)/v_(1))) then
      godunov = 0.
      ! print *, "flux nul"
   else 
      u_hat = (u_(2)/u_(1)) * (sqrt(u_(1))/(sqrt(u_(1)) + sqrt(v_(1)))) &
            + (v_(2)/v_(1)) * (sqrt(v_(1))/(sqrt(u_(1)) + sqrt(v_(1))))

      if(u_hat>0) then 
         godunov = flux(u_)

      else if(u_hat<0) then
         godunov = flux(v_)

      else if(u_hat==0) then
         godunov = (flux(u_)+flux(v_))*0.5
      print *, "delta choc"

      else 
         ! print *, "ERROR"

      end if
   end if

   ! print *, u_,v_,godunov
   return
end function godunov

subroutine correct_Flux(U,dt,dx,nx,F) 
   implicit none
   integer,intent(in)  :: nx
   real, dimension(2,0:nx-1), intent(inout)  :: U
   real, dimension(2,0:nx),   intent(out) :: F

   real, intent(in) ::  dt,dx
   
   interface
   function godunov(u_,v_)
      implicit none
      real, dimension(2), intent(in)  :: u_,v_
      real, dimension(2)              :: godunov
   end function godunov
   end interface

   interface
   function flux(u) result(f)
      implicit none
      real, dimension(2), intent(in)  :: u
      real, dimension(2)              :: f
   end function flux
   end interface

   real   :: s1,s2
   real, dimension(2)   :: z1,z2, z1_pred, z2_pred
   real, dimension(2)   :: F_, F_pred
   real                 :: u_hat, theta1, theta2, phi
   integer :: i

   
   do i=0,nx
      if(U(1,i)<1e-6) U(1,i) = 1e-6  
      if(i==0) then
         F(:,0)  = godunov(U(:,0),     U(:,0)) 

      else if(i==nx) then
         F(:,nx) = godunov(U(:,nx-1),  U(:,nx-1)) 

      else if(i /= 0 .and. i/=nx) then 
         F(:,i)  = godunov(U(:,i-1),   U(:,i))

      end if
   end do


   ! U(:,:) = U(:,:)- ((dt/dx)* (F  (:,1:nx)-F  (:,0:nx-1))) 
   do i=0,nx-1
      s1=0;s2=0;
      if((U(2,i-1)/U(1,i-1))<0 .and. 0<(U(2,i)/U(1,i))) then
         z1 = -flux(U(:,i-1));   z2= flux(U(:,i))
         s1 = U(2,i-1)/U(1,i-1); s2 = U(2,i)/U(1,i)

         ! print *,'z1 =',z1,'z2 =',z2
      else  
         u_hat = (U(2,i-1) /U(1,i-1))  * (sqrt(U(1,i-1))/(sqrt(U(1,i-1)) + sqrt(U(1,i)))) &
               + (U(2,i)   /U(1,i))    * (sqrt(U(1,i))  /(sqrt(U(1,i-1)) + sqrt(U(1,i))))

         if(u_hat<0) then
            z1 = flux(U(:,i))-flux(U(:,i-1)); z2 =0
            ! print *,'z1 =',z1
            s1 = u_hat; s2 = u_hat

         else 
            z1 = 0; z2 = flux(U(:,i))-flux(U(:,i-1))
            ! print *,'z2 =', z2
            s1 = u_hat; s2 = u_hat

         end if

      end if
      theta1 =0; theta2=0; phi=0
      if(z1(1) > 1e-6) theta1 = (z1(1)*z1_pred(1) + z1(2)*z1_pred(2))/(z1(1)*z1(1) + z1(2)*z1(2))
      if(z2(1) > 1e-6) theta2 = (z2(1)*z2_pred(1) + z2(2)*z2_pred(2))/(z2(1)*z2(1) + z2(2)*z2(2))

      ! print*, 'thetas = ',theta1, theta2
      ! print*, (z1(1)*z1_pred(1) + z1(2)*z1_pred(2))
      ! print*, (z1(1)*z1(1) + z1(2)*z1(2))
      phi = max(0.,min(1.,theta1));  z1 = phi*z1    
      ! print *,phi
      phi = max(0.,min(1.,theta2));  z2 = phi*z2    
      ! print *,phi

      F_ = 0.5*(    sign(1.,s1) * (1-(dt/dx) *abs(s1)))*z1
      F_ = F_+ 0.5*(sign(1.,s2) * (1-(dt/dx) *abs(s2)))*z2
   !    ! print *, sign(1.,s1), sign(1.,s2)
      ! print *, 'U=  ',U(:,i-1), U(:,i)
   !    print *, 'F,s=',F_, s1, s2
   !    print *, 'z=  ',z1, z2 
   !    print *, 'phi=',phi, theta1, theta2
      
   !    if(F_(1) /= F_(1)) call exit(0)

      if( (U(1,i)-((dt/dx)* (F(1,i+1)-F(1,i))) - ((dt/dx)* (F_(1)-F_pred(1))))<0. )  then

         U(:,i)= U(:,i)- ((dt/dx)* (F(:,i)-F(:,i-1)))
   !       print *, U(1,i)
   !       print *, "Ordre 1"
         F_pred =0; z1_pred =0; z2_pred =0

      else 
         U(:,i) = U(:,i)-((dt/dx)* (F(:,i+1)-F(:,i))) - ((dt/dx)* (F_-F_pred))
   !       print *, U(1,i)
   !       print *, "Ordre 2"
         z1_pred = z1; z2_pred = z2
         F_pred = F_
      end if

      
   !   
         ! U(:,i) = U(:,i)-((dt/dx)* (F(:,i+1)-F(:,i)))

   !    print *, 'U=  ',U(:,i-1), U(:,i)
   !    print *,'----------------'

   !    if(U(1,i)<0.) call exit(0)

   end do




end subroutine correct_Flux

function U_init(x) result(U)
   implicit none
   real, intent(in)    :: x
   real, dimension(2)  :: U
   
   u=0
   if(x>-2. .and. x<-1.) then
      u(1) = 2; u(2)=2

   else if (x>1. .and. x<5.) then
      u(1) = 1; u(2)=-1

   else
      u(1)= 1e-6; u(2)=1e-6
   end if

   return 
end function U_init

function U_exa(x,t,wl,wr) result(u)

   implicit none
   real, intent(in)    :: x,t,wl,wr
   real, dimension(2)  :: u
   real  :: T1, T2,  R
   real, dimension(2)  :: ul,ur
   real  :: Xhi, u_hat

   
   interface
      function U_init(x) 
         implicit none
         real, intent(in)  :: x
         real, dimension(2):: U_init
      end function U_init
   end interface


   ul = U_init(-wl-1 +1e-2);   ur = U_init( wr+1-1e-6)

   ! print *,wl, ul, wr, ur

   u_hat = (ul(2)/ul(1)) * (sqrt(ul(1))/(sqrt(ul(1)) + sqrt(ur(1)))) &
         + (ur(2)/ur(1)) * (sqrt(ur(1))/(sqrt(ul(1)) + sqrt(ur(1))))
   R = ul(1)/ur(1)

   T1 = wl/(   (ul(2)/ul(1))-u_hat) 
   T2 = (wr**2 +2*wr*wl*R + (wl**2) *R )/(2*wl*R*(ul(2)/ul(1) - ur(2)/ur(1)))

   u=0


   if(t<0.0) then
      if((x+1.) -(t+1) <0. .and. (x+1.+wl) -(t+1)>0.  ) u =ul
      if((x-1.) +(t+1) >0. .and. (x-1.-wr) +(t+1)<0.  ) u =ur
      
   else if(t<T1) then
      Xhi = u_hat*t
      if((x-xhi)<0. .and. (x+1.+wl) -(t+1)>0.) then 
         u = ul
      else if((x-xhi)>0. .and. (x-1.-wr) +(t+1)<0.) then
         u = ur
      end if

   else if(t<T2) then
      Xhi = ur(2)/ur(1)*t -wl*R + sqrt(2*wl*R*(ul(2)/ul(1)-ur(2)/ur(1))*t +wl**2 *(R**2-R))
      if((x-xhi)<0. .and. (x+1.+wl) -(t+1)>0.) then 
         u = 0
      else if((x-xhi)>0. .and. (x-1.-wr) +(t+1)<0.) then
         u = ur
      end if
   end if

   return

end function U_exa
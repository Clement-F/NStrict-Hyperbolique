
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

function correct_Flux(u_,v_,dt,dx) result(f)
   implicit none
   real, dimension(2), intent(in)  :: u_,v_
   real ::  dt,dx
   real, dimension(2)              :: f
   
   interface
   function flux(u) result(f)
      implicit none
      real, dimension(2), intent(in)  :: u
      real, dimension(2)              :: f
   end function flux
   end interface

   real   :: s1,s2
   real, dimension(2)   :: z1,z2
   real                 :: u_hat

   if((u_(2)/u_(1))<0 .and. 0<(v_(2)/v_(1))) then
      z1 = -flux(u_); z2= flux(v_)
      s1 = u_(2)/u_(1); s2 = v_(2)/v_(1)
   else  
      u_hat = (u_(2)/u_(1)) * (sqrt(u_(1))/(sqrt(u_(1)) + sqrt(v_(1)))) &
            + (v_(2)/v_(1)) * (sqrt(v_(1))/(sqrt(u_(1)) + sqrt(v_(1))))

      if(u_hat<0) then
         z1 = flux(v_)-flux(u_); z2 =0
         s1 = u_hat; s2 = u_hat

      else 
         z1 = 0; z2 = flux(v_)-flux(u_)
         s1 = u_hat; s2 = u_hat

      end if

   end if

   z1 = min(z1,u_,v_)

   f = 0.5*(    sign(1.,s1) * (1-(dt/dx) *abs(s1)))*z1
   f = f + 0.5*(sign(1.,s2) * (1-(dt/dx) *abs(s2)))*z2

   ! print *, u_,v_,f,s1,s2,z1,z2, u_hat

   return

end function correct_Flux

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
      ! Xhi = ur(2)/ur(1)*t -wl*R + sqrt(2*wl*R*(ul(2)/ul(1)-ur(2)/ur(1))*t +wl**2 *(R**2-R))
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
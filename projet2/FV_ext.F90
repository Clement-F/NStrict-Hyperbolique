
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
      print *, "flux nul"
   else 
      u_hat = (sqrt(u_(1))/(sqrt(u_(1)) + sqrt(v_(1))))*(u_(2)/u_(1)) &
            + (sqrt(v_(1))/(sqrt(u_(1)) + sqrt(v_(1))))*(v_(2)/v_(1)) 

      if(u_hat>0) then 
         godunov = flux(u_)

      else if(u_hat<0) then
         godunov = flux(v_)

      else if(u_hat==0) then
         godunov = (flux(u_)+flux(v_))*0.5
      print *, "delta choc"

      else 
         print *, "ERROR"

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

   real   :: s1,s2
   real, dimension(2)   :: w1,w2
   real                 :: u_hat

   if(u_(2)/u_(1)<0 .and. 0<v_(2)/v_(1)) then
      w1 = -u_; w2= v_
      s1 = u_(2)/u_(1); s2 = v_(2)/v_(1)
   else  
      u_hat = (sqrt(u_(1))/(sqrt(u_(1)) + sqrt(v_(1))))*(u_(2)/u_(1)) &
            + (sqrt(v_(1))/(sqrt(u_(1)) + sqrt(v_(1))))*(v_(2)/v_(1)) 

      if(u_hat<0) then
         w1 = v_-u_; w2 =0
         s1 = u_hat; s2 = u_hat

      else 
         w1 = 0; w2 = v_-u_
         s1 = u_hat; s2 = u_hat

      end if

   end if

   f = 0.5*(abs(s1)* (1-(dt/dx) *abs(s1)))*w1
   f = f + 0.5*(abs(s2)* (1-(dt/dx) *abs(s2)))*w2

   print *, u_,v_,f,s1,s2,w1,w2, u_hat

   return

end function correct_Flux

function U_init(x) result(U)
   implicit none
   real, intent(in)    :: x
   real, dimension(2)  :: U
   
   u=0
   if(x>-2 .and. x<-1) then
      u(1) = 2; u(2)=2

   else if (x>1 .and. x<5) then
      u(1) = 1; u(2)=-1

   else
      u(1)= 1e-6; u(2)=1e-6
   end if

   return 
end function U_init


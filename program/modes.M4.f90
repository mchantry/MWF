!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
module modes
  !*************************************************************************
  use parameters
  use rp_emulator
  implicit none
  save
  
contains

  function udotgradu(v,vx,vz) result(ans)
    type(rpe_var) :: ans(i_KK)
    type(rpe_var),intent(in) :: v(i_KK),vx(i_KK),vz(i_KK)
    ans(1) =v(1)*vx(1) + (v(2)*vx(2))/rp2 + (v(3)*vx(3))/rp2 + (v(4)*vx(4))/rp2 + v(8)*vz(1) + (v(9)*vz(2))/rp2 + (v(10)*vz(3))/rp2 + (v(11)*vz(4))/rp2 + (d_beta*v(2)*v(5))/rp2 - d_beta*v(3)*v(6) + (3*d_beta*v(4)*v(7))/rp2
    ans(2) =v(1)*vx(2) + v(2)*vx(1) - (v(2)*vx(3))/rp2 - (v(3)*vx(2))/rp2 + (v(3)*vx(4))/rp2 + (v(4)*vx(3))/rp2 + v(8)*vz(2) + v(9)*vz(1) - (v(9)*vz(3))/rp2 - (v(10)*vz(2))/rp2 + (v(10)*vz(4))/rp2 + (v(11)*vz(3))/rp2 + (d_beta*v(2)*v(6))/rp2 - d_beta*v(3)*v(5) + d_beta*v(3)*v(7) - (3*d_beta*v(4)*v(6))/rp2
    ans(3) =v(1)*vx(3) - (v(2)*vx(2))/rp2 + v(3)*vx(1) + (v(2)*vx(4))/rp2 + (v(4)*vx(2))/rp2 + v(8)*vz(3) - (v(9)*vz(2))/rp2 + v(10)*vz(1) + (v(9)*vz(4))/rp2 + (v(11)*vz(2))/rp2 + (d_beta*v(2)*v(5))/rp2 + (d_beta*v(2)*v(7))/rp2 + (3*d_beta*v(4)*v(5))/rp2
    ans(4) =v(1)*vx(4) + (v(2)*vx(3))/rp2 + (v(3)*vx(2))/rp2 + v(4)*vx(1) + v(8)*vz(4) + (v(9)*vz(3))/rp2 + (v(10)*vz(2))/rp2 + v(11)*vz(1) + (d_beta*v(2)*v(6))/rp2 - d_beta*v(3)*v(5)
    ans(5) =v(1)*vx(5) + (v(2)*vx(6))/rp2 + (v(3)*vx(5))/rp2 + (v(3)*vx(7))/rp2 + (v(4)*vx(6))/rp2 + v(8)*vz(5) + (v(9)*vz(6))/rp2 + (v(10)*vz(5))/rp2 + (v(10)*vz(7))/rp2 + (v(11)*vz(6))/rp2 + (d_beta*v(5)*v(6))/rp2 - (d_beta*v(6)*v(7))/rp2
    ans(6) =v(1)*vx(6) + (v(2)*vx(5))/rp2 - (v(2)*vx(7))/rp2 + (v(4)*vx(5))/rp2 + v(8)*vz(6) + (v(9)*vz(5))/rp2 - (v(9)*vz(7))/rp2 + (v(11)*vz(5))/rp2 - (d_beta*v(5)*v(5))/rp2 - d_beta*v(5)*v(7)
    ans(7) =v(1)*vx(7) - (v(2)*vx(6))/rp2 + (v(3)*vx(5))/rp2 + v(8)*vz(7) - (v(9)*vz(6))/rp2 + (v(10)*vz(5))/rp2 + (3*d_beta*v(5)*v(6))/rp2
    ans(8) =v(1)*vx(8) + (v(2)*vx(9))/rp2 + (v(3)*vx(10))/rp2 + (v(4)*vx(11))/rp2 + v(8)*vz(8) + (v(9)*vz(9))/rp2 + (v(10)*vz(10))/rp2 + (v(11)*vz(11))/rp2 + (d_beta*v(5)*v(9))/rp2 - d_beta*v(6)*v(10) + (3*d_beta*v(7)*v(11))/rp2
    ans(9) =v(1)*vx(9) + v(2)*vx(8) - (v(2)*vx(10))/rp2 - (v(3)*vx(9))/rp2 + (v(3)*vx(11))/rp2 + (v(4)*vx(10))/rp2 + v(8)*vz(9) + v(9)*vz(8) - (v(9)*vz(10))/rp2 - (v(10)*vz(9))/rp2 + (v(10)*vz(11))/rp2 + (v(11)*vz(10))/rp2 - d_beta*v(5)*v(10) + (d_beta*v(6)*v(9))/rp2 - (3*d_beta*v(6)*v(11))/rp2 + d_beta*v(7)*v(10)
    ans(10) =v(1)*vx(10) - (v(2)*vx(9))/rp2 + v(3)*vx(8) + (v(2)*vx(11))/rp2 + (v(4)*vx(9))/rp2 + v(8)*vz(10) - (v(9)*vz(9))/rp2 + v(10)*vz(8) + (v(9)*vz(11))/rp2 + (v(11)*vz(9))/rp2 + (d_beta*v(5)*v(9))/rp2 + (3*d_beta*v(5)*v(11))/rp2 + (d_beta*v(7)*v(9))/rp2
    ans(11) =v(1)*vx(11) + (v(2)*vx(10))/rp2 + (v(3)*vx(9))/rp2 + v(4)*vx(8) + v(8)*vz(11) + (v(9)*vz(10))/rp2 + (v(10)*vz(9))/rp2 + v(11)*vz(8) - d_beta*v(5)*v(10) + (d_beta*v(6)*v(9))/rp2
  end function udotgradu

  function udotgradu2(u,ux,uz) result(ans)
    type(rpe_var) :: ans(i_KK)
    type(rpe_var),intent(in) :: u(i_KK),ux(i_KK),uz(i_KK)
    ans(1:i_K0)= nluw(u(1:i_K0),ux(1:i_K0)) + nluyv(u(1:i_K0),u(i_K0+1:2*i_K0-1)) + nluw(u(2*i_K0:i_KK),uz(1:i_K0))
    ans(i_K0+1:2*i_K0-1)= nluv(u(1:i_K0),ux(i_K0+1:2*i_K0-1)) + nlvvy(u(i_K0+1:2*i_K0-1))        + nluv(u(2*i_K0:i_KK),uz(i_K0+1:2*i_K0-1))
    ans(2*i_K0:i_KK)=nluw(u(1:i_K0),ux(2*i_K0:i_KK))+ nluyv(u(2*i_K0:i_KK),u(i_K0+1:2*i_K0-1))+ nluw(u(2*i_K0:i_KK),uz(2*i_K0:i_KK))
  end function udotgradu2
  
  function nluw(u,w) result(ans)
    type(rpe_var) :: ans(0:i_K1)
    type(rpe_var),intent(in) :: u(0:i_K1),w(0:i_K1)
    ans(0)=u(0)*w(0) + (u(1)*w(1))/rp2 + (u(2)*w(2))/rp2 + (u(3)*w(3))/rp2
    ans(1)=u(0)*w(1) + u(1)*w(0) - (u(1)*w(2))/rp2 - (u(2)*w(1))/rp2 + (u(2)*w(3))/rp2 + (u(3)*w(2))/rp2	
    ans(2)=u(0)*w(2) - (u(1)*w(1))/rp2 + u(2)*w(0) + (u(1)*w(3))/rp2 + (u(3)*w(1))/rp2
    ans(3)=u(0)*w(3) + (u(1)*w(2))/rp2 + (u(2)*w(1))/rp2 + u(3)*w(0)
    RETURN
  end function nluw

  function nluyv(u,v) result(ans)
    type(rpe_var) :: ans(0:i_K1)
    type(rpe_var),intent(in) :: u(0:i_K1),v(1:i_K1)
    ans(0)=(d_beta*u(1)*v(1))/rp2 - d_beta*u(2)*v(2) + (3*d_beta*u(3)*v(3))/rp2
    ans(1)=(d_beta*u(1)*v(2))/rp2 - d_beta*u(2)*v(1) + d_beta*u(2)*v(3) - (3*d_beta*u(3)*v(2))/rp2
    ans(2)=(d_beta*u(1)*v(1))/rp2 + (d_beta*u(1)*v(3))/rp2 + (3*d_beta*u(3)*v(1))/rp2
    ans(3)=(d_beta*u(1)*v(2))/rp2 - d_beta*u(2)*v(1)
    RETURN
  end function nluyv
  
  function nluv(u,v) result(ans)
    type(rpe_var) :: ans(1:i_K1)
    type(rpe_var),intent(in) :: u(0:i_K1),v(1:i_K1)
    ans(1)=u(0)*v(1) + (u(1)*v(2))/rp2 + (u(2)*v(1))/rp2 + (u(2)*v(3))/rp2 + (u(3)*v(2))/rp2
    ans(2)=u(0)*v(2) + (u(1)*v(1))/rp2 - (u(1)*v(3))/rp2 + (u(3)*v(1))/rp2
    ans(3)=u(0)*v(3) - (u(1)*v(2))/rp2 + (u(2)*v(1))/rp2
  end function nluv
  
  function nlvvy(v) result(ans)
    type(rpe_var) :: ans(1:i_K1)
    type(rpe_var),intent(in) :: v(1:i_K1)
    ans(1)=(d_beta*v(1)*v(2))/rp2 - (d_beta*v(2)*v(3))/rp2
    ans(2)=- (d_beta*v(1)*v(1))/rp2 - d_beta*v(1)*v(3)
    ans(3)=(3*d_beta*v(1)*v(2))/rp2
  end function nlvvy
  
  function nlvv(v) result(ans)
    type(rpe_var) :: ans(0:i_K1)
    type(rpe_var),intent(in) :: v(1:i_K1)
    ans(0)=v(1)*v(1)/rp2 + v(2)*v(2)/rp2 + v(3)*v(3)/rp2
    ans(1)=v(1)*v(2) - v(2)*v(3)
    ans(2)=v(1)*v(1)/rp2 + v(3)*v(1)
    ans(3)=v(1)*v(2)
  end function nlvv
  
  type(rpe_var) function velU(V,y) 
    type(rpe_var), intent(in) :: V(i_KK),y
    velU= V(1) &
         +sin(d_beta*y) * V(2) &
         +cos(2d0*d_beta*y) * V(3)&
         +sin(3d0*d_beta*y) * V(4)
    return
  end function velU
  
  type(rpe_var) function velW(V,y)
    type(rpe_var), intent(in) :: V(i_KK),y
    velW = V(8)&
         +sin(d_beta*y) * V(9)&
         +cos(2d0*d_beta*y) * V(10)&
         +sin(3d0*d_beta*y) * V(11)
    return
  end function velW

  type(rpe_var) function velV(V,y)
    type(rpe_var), intent(in) :: V(i_KK),y
    velV= cos(d_beta*y) * V(5)&
         +sin(2d0*d_beta*y) * V(6)&
         +cos(3d0*d_beta*y) * V(7)
    return
  end function velV

end module modes


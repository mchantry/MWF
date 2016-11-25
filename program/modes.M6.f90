!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
module modes
  !*************************************************************************
  use parameters
  use rp_emulator
  !   use timestep
  implicit none
  save
contains

  function udotgradu(u,ux,uz) result(ans)
    type(rpe_var) :: ans(i_KK)
    type(rpe_var),intent(in) :: u(i_KK),ux(i_KK),uz(i_KK)
    ans(1:i_K0)= nluw(u(1:i_K0),ux(1:i_K0)) + nluyv(u(1:i_K0),u(i_K0+1:2*i_K0-1)) + nluw(u(2*i_K0:i_KK),uz(1:i_K0))
    ans(i_K0+1:2*i_K0-1)= nluv(u(1:i_K0),ux(i_K0+1:2*i_K0-1)) + nlvvy(u(i_K0+1:2*i_K0-1))        + nluv(u(2*i_K0:i_KK),uz(i_K0+1:2*i_K0-1))
    ans(2*i_K0:i_KK)=nluw(u(1:i_K0),ux(2*i_K0:i_KK))+ nluyv(u(2*i_K0:i_KK),u(i_K0+1:2*i_K0-1))+ nluw(u(2*i_K0:i_KK),uz(2*i_K0:i_KK))
  end function udotgradu
  
  function nluw(u,w) result(ans)
    type(rpe_var) :: ans(0:i_K1)
    type(rpe_var),intent(in) :: u(0:i_K1),w(0:i_K1)
    ans(0) = u(0)*w(0) + (u(1)*w(1))/2 + (u(2)*w(2))/2 + (u(3)*w(3))/2 + (u(4)*w(4))/2 + (u(5)*w(5))/2
    ans(1) = u(0)*w(1) + u(1)*w(0) - (u(1)*w(2))/2 - (u(2)*w(1))/2 + (u(2)*w(3))/2 + (u(3)*w(2))/2 - (u(3)*w(4))/2 - (u(4)*w(3))/2 + (u(4)*w(5))/2 + (u(5)*w(4))/2
    ans(2) = u(0)*w(2) - (u(1)*w(1))/2 + u(2)*w(0) + (u(1)*w(3))/2 + (u(3)*w(1))/2 + (u(2)*w(4))/2 + (u(4)*w(2))/2 + (u(3)*w(5))/2 + (u(5)*w(3))/2
    ans(3) = u(0)*w(3) + (u(1)*w(2))/2 + (u(2)*w(1))/2 + u(3)*w(0) - (u(1)*w(4))/2 - (u(4)*w(1))/2 + (u(2)*w(5))/2 + (u(5)*w(2))/2
    ans(4) = u(0)*w(4) - (u(1)*w(3))/2 + (u(2)*w(2))/2 - (u(3)*w(1))/2 + u(4)*w(0) + (u(1)*w(5))/2 + (u(5)*w(1))/2
    ans(5) = u(0)*w(5) + (u(1)*w(4))/2 + (u(2)*w(3))/2 + (u(3)*w(2))/2 + (u(4)*w(1))/2 + u(5)*w(0)
  end function nluw
  
  function nlvv(v) result(ans)
    type(rpe_var) :: ans(0:i_K1)
    type(rpe_var),intent(in) :: v(1:i_K1)
    
    ans(0) =v(1)**2/2 + v(2)**2/2 + v(3)**2/2 + v(4)**2/2 + v(5)**2/2 
    ans(1) = v(1)*v(2) - v(2)*v(3) + v(3)*v(4) - v(4)*v(5) 
    ans(2) = v(1)**2/2 + v(3)*v(1) + v(2)*v(4) + v(3)*v(5) 
    ans(3) = v(1)*v(2) + v(1)*v(4) - v(2)*v(5) 
    ans(4) = - v(2)**2/2 + v(1)*v(3) + v(1)*v(5) 
    ans(5) = v(1)*v(4) + v(2)*v(3)
  end function nlvv
  
  
  function nluv(u,v) result(ans)
    type(rpe_var) :: ans(1:i_K1)
    type(rpe_var),intent(in) :: u(0:i_K1),v(1:i_K1)
    
    ans(1) = u(0)*v(1) + (u(1)*v(2))/2 + (u(2)*v(1))/2 + (u(2)*v(3))/2 + (u(3)*v(2))/2 + (u(3)*v(4))/2 + (u(4)*v(3))/2 + (u(4)*v(5))/2 + (u(5)*v(4))/2 
    ans(2) = u(0)*v(2) + (u(1)*v(1))/2 - (u(1)*v(3))/2 + (u(3)*v(1))/2 + (u(2)*v(4))/2 - (u(4)*v(2))/2 - (u(3)*v(5))/2 + (u(5)*v(3))/2 
    ans(3) = u(0)*v(3) - (u(1)*v(2))/2 + (u(2)*v(1))/2 + (u(1)*v(4))/2 + (u(4)*v(1))/2 + (u(2)*v(5))/2 + (u(5)*v(2))/2 
    ans(4) = u(0)*v(4) + (u(1)*v(3))/2 + (u(2)*v(2))/2 + (u(3)*v(1))/2 - (u(1)*v(5))/2 + (u(5)*v(1))/2 
    ans(5) = u(0)*v(5) - (u(1)*v(4))/2 + (u(2)*v(3))/2 - (u(3)*v(2))/2 + (u(4)*v(1))/2
  end function nluv
  
  function nluyv(u,v) result(ans)
    type(rpe_var) :: ans(0:i_K1)
    type(rpe_var),intent(in) :: u(0:i_K1),v(1:i_K1)
    
    ans(0) = (d_beta*u(1)*v(1))/2 - d_beta*u(2)*v(2) + (3*d_beta*u(3)*v(3))/2 - 2*d_beta*u(4)*v(4) + (5*d_beta*u(5)*v(5))/2 
    ans(1) = (d_beta*u(1)*v(2))/2 - d_beta*u(2)*v(1) + d_beta*u(2)*v(3) - (3*d_beta*u(3)*v(2))/2 + (3*d_beta*u(3)*v(4))/2 - 2*d_beta*u(4)*v(3) + 2*d_beta*u(4)*v(5) - (5*d_beta*u(5)*v(4))/2 
    ans(2) = (d_beta*u(1)*v(1))/2 + (d_beta*u(1)*v(3))/2 + (3*d_beta*u(3)*v(1))/2 - d_beta*u(2)*v(4) - 2*d_beta*u(4)*v(2) + (3*d_beta*u(3)*v(5))/2 + (5*d_beta*u(5)*v(3))/2 
    ans(3) = (d_beta*u(1)*v(2))/2 - d_beta*u(2)*v(1) + (d_beta*u(1)*v(4))/2 - 2*d_beta*u(4)*v(1) + d_beta*u(2)*v(5) - (5*d_beta*u(5)*v(2))/2 
    ans(4) = (d_beta*u(1)*v(3))/2 + d_beta*u(2)*v(2) + (3*d_beta*u(3)*v(1))/2 + (d_beta*u(1)*v(5))/2 + (5*d_beta*u(5)*v(1))/2 
    ans(5) = (d_beta*u(1)*v(4))/2 - d_beta*u(2)*v(3) + (3*d_beta*u(3)*v(2))/2 - 2*d_beta*u(4)*v(1)
  end function nluyv
  
  
  function nlvvy(v) result(ans)
    type(rpe_var) :: ans(1:i_K1)
    type(rpe_var),intent(in) :: v(1:i_K1)
    
    ans(1) = (d_beta*v(1)*v(2))/2 - (d_beta*v(2)*v(3))/2 + (d_beta*v(3)*v(4))/2 - (d_beta*v(4)*v(5))/2 
    ans(2) = - (d_beta*v(1)**2)/2 - d_beta*v(1)*v(3) - d_beta*v(2)*v(4) - d_beta*v(3)*v(5) 
    ans(3) = (3*d_beta*v(1)*v(2))/2 + (3*d_beta*v(1)*v(4))/2 - (3*d_beta*v(2)*v(5))/2 
    ans(4) = d_beta*v(2)**2 - 2*d_beta*v(1)*v(3) - 2*d_beta*v(1)*v(5) 
    ans(5) = (5*d_beta*v(1)*v(4))/2 + (5*d_beta*v(2)*v(3))/2
  end function nlvvy
  
  type(rpe_var) function velU(V,y) 
    type(rpe_var), intent(in) :: V(i_KK),y
    velU= V(1) &
         +sin(d_beta*y) * V(2) &
         +cos(2d0*d_beta*y) * V(3)&
         +sin(3d0*d_beta*y) * V(4)&
         +cos(4d0*d_beta*y) * V(5)&
         +sin(5d0*d_beta*y) * V(6)
    return
  end function velU
  
  type(rpe_var) function velW(V,y)
    type(rpe_var), intent(in) :: V(i_KK),y
    velW = V(12)&
         +sin(d_beta*y) * V(13)&
         +cos(2d0*d_beta*y) * V(14)&
         +sin(3d0*d_beta*y) * V(15)&
         +cos(4d0*d_beta*y) * V(16)&
         +sin(5d0*d_beta*y) * V(17)
    return
  end function velW

  type(rpe_var) function velV(V,y)
    type(rpe_var), intent(in) :: V(i_KK),y
    velV= cos(d_beta*y) * V(7)&
         +sin(2d0*d_beta*y) * V(8)&
         +cos(3d0*d_beta*y) * V(9)&
         +sin(4d0*d_beta*y) * V(10)&
         +cos(5d0*d_beta*y) * V(11)
    return
  end function velV

end module modes


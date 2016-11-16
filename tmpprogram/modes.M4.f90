!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
module modes
  !*************************************************************************
  use parameters
  implicit none
  save
  
contains
  
  function nluw(u,w) result(ans)
    double precision :: ans(0:i_K1)
    double precision,intent(in) :: u(0:i_K1),w(0:i_K1)
    ans(0)=u(0)*w(0) + (u(1)*w(1))/2 + (u(2)*w(2))/2 + (u(3)*w(3))/2
    ans(1)=u(0)*w(1) + u(1)*w(0) - (u(1)*w(2))/2 - (u(2)*w(1))/2 + (u(2)*w(3))/2 + (u(3)*w(2))/2	
    ans(2)=u(0)*w(2) - (u(1)*w(1))/2 + u(2)*w(0) + (u(1)*w(3))/2 + (u(3)*w(1))/2
    ans(3)=u(0)*w(3) + (u(1)*w(2))/2 + (u(2)*w(1))/2 + u(3)*w(0)
    RETURN
  end function nluw

  function nluyv(u,v) result(ans)
    double precision :: ans(0:i_K1)
    double precision,intent(in) :: u(0:i_K1),v(1:i_K1)
    ans(0)=(d_beta*u(1)*v(1))/2 - d_beta*u(2)*v(2) + (3*d_beta*u(3)*v(3))/2
    ans(1)=(d_beta*u(1)*v(2))/2 - d_beta*u(2)*v(1) + d_beta*u(2)*v(3) - (3*d_beta*u(3)*v(2))/2
    ans(2)=(d_beta*u(1)*v(1))/2 + (d_beta*u(1)*v(3))/2 + (3*d_beta*u(3)*v(1))/2
    ans(3)=(d_beta*u(1)*v(2))/2 - d_beta*u(2)*v(1)
    RETURN
  end function nluyv
  
  function nluv(u,v) result(ans)
    double precision :: ans(1:i_K1)
    double precision,intent(in) :: u(0:i_K1),v(1:i_K1)
    ans(1)=u(0)*v(1) + (u(1)*v(2))/2 + (u(2)*v(1))/2 + (u(2)*v(3))/2 + (u(3)*v(2))/2
    ans(2)=u(0)*v(2) + (u(1)*v(1))/2 - (u(1)*v(3))/2 + (u(3)*v(1))/2
    ans(3)=u(0)*v(3) - (u(1)*v(2))/2 + (u(2)*v(1))/2
  end function nluv
  
  function nlvvy(v) result(ans)
    double precision :: ans(1:i_K1)
    double precision,intent(in) :: v(1:i_K1)
    ans(1)=(d_beta*v(1)*v(2))/2 - (d_beta*v(2)*v(3))/2
    ans(2)=- (d_beta*v(1)**2)/2 - d_beta*v(1)*v(3)
    ans(3)=(3*d_beta*v(1)*v(2))/2
  end function nlvvy
  
  function nlvv(v) result(ans)
    double precision :: ans(0:i_K1)
    double precision,intent(in) :: v(1:i_K1)
    ans(0)=v(1)**2/2 + v(2)**2/2 + v(3)**2/2
    ans(1)=v(1)*v(2) - v(2)*v(3)
    ans(2)=v(1)**2/2 + v(3)*v(1)
    ans(3)=v(1)*v(2)
  end function nlvv
  
  double precision function velU(V,y) 
    double precision, intent(in) :: V(i_KK),y
    velU= V(1) &
         +dsin(d_beta*y) * V(2) &
         +dcos(2d0*d_beta*y) * V(3)&
         +dsin(3d0*d_beta*y) * V(4)
    return
  end function velU
  
  double precision function velW(V,y)
    double precision, intent(in) :: V(i_KK),y
    velW = V(8)&
         +dsin(d_beta*y) * V(9)&
         +dcos(2d0*d_beta*y) * V(10)&
         +dsin(3d0*d_beta*y) * V(11)
    return
  end function velW

  double precision function velV(V,y)
    double precision, intent(in) :: V(i_KK),y
    velV= dcos(d_beta*y)     * V(5)&
         +dsin(2d0*d_beta*y) * V(6)&
         +dcos(3d0*d_beta*y) * V(7)
    return
  end function velV

end module modes


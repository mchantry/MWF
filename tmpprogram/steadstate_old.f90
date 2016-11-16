!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
module steadystate
  !*************************************************************************
  use mpif
  use variables
  use transform
  use velocity
  implicit none
  save

  type (phys) :: up,upx,upz
contains

  subroutine ss_mpt2vec(m,v)
    type(mpt), intent(in) :: m
    double complex,intent(out) :: v(i_KMN)
    integer :: i
    _loop_kmn_vars
    i = 0
    do n = 0, i_NN1 
      do m = 0, i_M1
         do k = 1, i_K
          i = i + 1
          v(i)=dcmplx(m%Re(k,m,n),m%Im(k,m,n))
        end do
      end do
    end do
  end subroutine ss_mpt2vec

  subroutine ss_vec2mpt(v,m)
    type(mpt), intent(out) :: m
    double complex, intent(in) :: v(i_KMN)
    integer :: i
    _loop_kmn_vars
    i = 0
    do n = 0, i_NN1 
      do m = 0, i_M1
         do k = 1, i_K
          i = i + 1
          m%Re(k,m,n) = dble(v(i))
          m%Im(k,m,n) = dimag(v(i))
        end do
      end do
    end do
  end subroutine ss_vec2mpt

  subroutine ssw_RHS(u,v)
    double complex, intent(in) :: u(i_KMN)
    double complex, intent(out) :: v(i_KMN)
    type(mpt), intent(in) :: uu
    type(mpt), intent(out) :: vv
    call ss_vec2mpt(u,uu)
    call ss_RHS(uu,vv)
    call ss_mpt2vec(vv,v)
  end subroutine ssw_RHS

  subroutine ssw_LHS(u,v)
    double complex, intent(in) :: u(i_KMN)
    double complex, intent(out) :: v(i_KMN)
    type(mpt), intent(in) :: uu
    type(mpt), intent(out) :: vv
    call ss_vec2mpt(u,uu)
    call ss_LHS(uu,vv)
    call ss_mpt2vec(vv,v)
  end subroutine ssw_LHS

  subroutine ss_RHS(u,v)
    type(mpt), intent(in) :: u
    type(mpt), intent(out) :: v

    vel_c = u
    call vel_TS()
    v = vel_c - u
  end subroutine ss_RHS

  subroutine ss_l_nonlinear(uin,v)
    type(mpt), intent(in) :: uin
    type(mpt), intent(out) :: v
    type(spec) :: u,ux,uz
    type(phys) :: p,px,pz,ans
    _loop_mn_vars
    
    call var_mpt2spec(vel_c,u)
    call var_spec_grad(u,ux,uz)
    call tra_spec2phys(u,p)
    call tra_spec2phys(ux,px)
    call tra_spec2phys(uz,pz)
    
    do m = 0, var_M%pH1
       do n = 0,i_N1
          ans%Re(:,n,m)=linear_udotgradu(up%Re(:,n,m),upx%Re(:,n,m),upz%Re(:,n,m),p%Re(:,n,m),px%Re(:,n,m),pz%Re(:,n,m))
       end do
    end do
    
    call tra_phys2spec(ans,u)
    call var_spec2mpt(u,v)

  end subroutine ss_l_nonlinear

    subroutine ss_LHS()
    type(mpt), intent(in) :: u
    type(mpt), intent(out) :: v    
    _loop_kmn_vars
    
    call ss_l_nonlinear(u,v)
    _loop_kmn_begin
    v%Re(k,m,n) = v%Re(k,m,n) * d_dt
    v%Im(k,m,n) = v%Im(k,m,n) * d_dt
    v%Re(k,m,n) = v%Re(k,m,n) + rhs%Re(k,m,n)*u%Re(k,m,n)
    v%Im(k,m,n) = v%Im(k,m,n) + rhs%Re(k,m,n)*u%Im(k,m,n)
    v%Re(k,m,n) = v%Re(k,m,n) / lhs%Re(k,m,n)
    v%Im(k,m,n) = v%Im(k,m,n) / lhs%Re(k,m,n)
    

    v%Re(k,m,n) = v%Re(k,m,n) - u%Re(k,m,n)
    v%Im(k,m,n) = v%Im(k,m,n) - u%Im(k,m,n)  
    _loop_kmn_end

    if (var_N%pH0 == 0) then
       vel_c2%Re(1,0,0)=0d0
    end if

  end subroutine ss_LHS

  subroutine newton()

    integer :: mxmv = 1000
    integer, parameter :: l = 1
    double precision :: UU(i_2KMN),b(i_2KMN),x(i_2KMN), work(i_2KMN,2*l+5))
    integer :: lwd(i_2KMN*(2*l+5)), info
    double precision :: error

!    Take UU guess
  while (TRUE) then 
    call ssw_RHS(UU,b)
    call ss_size(b,error)
    if (error < big_tol) then
      call ss_save(UU)
      return
    end if
    x=0d0
    call bicgstab2 (.TRUE.,l, i_KMN, x, b, ssw_LHS, .FALSE., small_tol,
     &      'rel',mxmv, work, ldw, info)
    if (info /=0) then
        print*, info
        return
    else  
      UU = UU - x 
    end if 
  end do
    
  end subroutine newton

end module steadystate


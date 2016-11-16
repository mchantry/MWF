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
  type (mpt) :: um,umx,umz
contains
#ifndef _MPI
  subroutine ss_mpt2vec(m,v)
    type(mpt), intent(in) :: m
    double complex,intent(out) :: v(i_2KMN)
    integer :: i
    _loop_kmn_vars
    i = 0
    _loop_kmn_begin
          i = i + 1
          v(i)=m%Re(k,m,n)
          i = i + 1
          v(i)=m%Im(k,m,n)
    _loop_kmn_end
  end subroutine ss_mpt2vec

#else
  subroutine ss_mpt2vec(m,v)
    type(mpt), intent(in) :: m
    double complex,intent(out) :: v(i_2KMN)
    integer :: i
    type(mpt) :: c1
    _loop_kmn_vars
  if(mpi_rnk==0) then
    i=0
    _loop_kmn_begin
      i=i+1
      v(i)=m%Re(k,m,n)
      i = i + 1
      v(i)=m%Im(k,m,n)    
    _loop_kmn_end
    do r = 1, mpi_sze-1
      pN0 = var_N%pH0_(r)
      pN1 = var_N%pH1_(r)
      mpi_tg = r
      call mpi_recv( c1%Re(1,0,0), i_M*(pN1+1)*i_K, mpi_double_precision,  &
        r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
      call mpi_recv( c1%Im(1,0,0), i_M*(pN1+1)*i_K, mpi_double_precision,  &
        r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
!      i=(pN0-1)*i_M*i_K*2 !Unnecessary
      _loop_kmn_begin
        i = i + 1
        v(i)=m%Re(k,m,n)
        i = i + 1
        v(i)=m%Im(k,m,n)    
      _loop_kmn_end
    end do
      else
         mpi_tg = mpi_rnk
         call mpi_send( m%Re(1,0,0), i_M*(var_N%pH1+1)*i_K, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( m%Im(1,0,0), i_M*(var_N%pH1+1)*i_K, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
      end if
  end subroutine ss_vec2mpt
#endif      

  subroutine ss_vec2mpt(v,m)
    type(mpt), intent(out) :: m
    double precision, intent(in) :: v(i_2KMN)
    integer :: i
    _loop_kmn_vars
    i = 0
    do n = 0, i_NN1 
      do m = 0, i_M1
         do k = 1, i_K
          i = i + 1
          m%Re(k,m,n) = v(i)
          i = i + 1
          m%Im(k,m,n) = v(i)
        end do
      end do
    end do
  end subroutine ss_vec2mpt

  subroutine ss_update_mpt(v)
    double precision, intent(in) :: v(i_2KMN)
    call ss_vec2mpt(v,um)
    call vel_addlam(um)
    call var_spec_grad(um,umx,umz)
    call tra_spec2phys(um,up)
    call tra_spec2phys(umx,upx)
    call tra_spec2phys(umz,upz)
  end subroutine ss_update_mpt

  subroutine ssw_RHS(u,v)
    double complex, intent(in) :: u(i_2KMN)
    double complex, intent(out) :: v(i_2KMN)
    type(mpt), intent(in) :: uu
    type(mpt), intent(out) :: vv
    call ss_vec2mpt(u,uu)
    call ss_RHS(uu,vv)
    call ss_mpt2vec(vv,v)
  end subroutine ssw_RHS

  subroutine ssw_LHS(u,v)
    double complex, intent(in) :: u(i_2KMN)
    double complex, intent(out) :: v(i_2KMN)
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

  subroutine ss_LHS(u,v)
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

  subroutine ss_save(u)
    double precision, intent(in) :: u(i_2KMN)
    type(mpt) :: m
    call ss_vec2mpt(u,m)
    call io_save_state(m)
  end subroutine ss_save

  subroutine newton()
    integer :: mxmv = 1000
    integer, parameter :: l = 1
    double precision :: UU(i_2KMN),b(i_2KMN),x(i_2KMN), work(i_2KMN,2*l+5))
    integer :: lwd(i_2KMN*(2*l+5)), info
    double precision :: error
!    Take UU guess
    do while (TRUE) then 
      call ss_update_mpt(UU) !Updates the MPT and Phys versions of UU
       call ssw_RHS(UU,b)
       call ss_size(b,error)
       if (error < big_tol) then
          call ss_save(UU)
          return
       end if
       x=0d0
       call bicgstab2 (.TRUE.,l, i_2KMN, x, b, ssw_LHS, .FALSE., small_tol,
       &      'rel',mxmv, work, ldw, info)
#ifdef _MPI
      call mpi_bcast(info,1,mpi_integer, 0,mpi_comm_world,mpi_er) !CONFIRM THE WORKING OF THIS CALL
#endif
       if (info /=0) then
          print*, info
#ifdef _MPI
        call mpi_barrier(mpi_comm_world, mpi_er)
        call mpi_finalize(mpi_er)
#else
          return
#endif
        end if
       else  
          UU = UU - x 
       end if
    end do
  end subroutine newton

end module steadystate


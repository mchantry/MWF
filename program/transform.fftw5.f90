#include "../parallel.h"
!*************************************************************************
! Transformation using the FFT, de-aliased.
! Written for FFTW version 3 -- documentation at  www.fftw.org
! See FFTW reference, sections 2.3, 4.3.3, 4.4.1-2
!
! \sum_{-Ma}^{Ma} == \sum_0^{3M-1}
! logical size of DFT n==3M=Th
!
!*************************************************************************
 module transform
!*************************************************************************
   use parameters
   use variables
   use mpif
   use FFTW3
!   include 'fftw3.f03'

   implicit none
   save
   !i_NN = (i_N /2 + 1)
   !i_Na = (i_3N/2 + 1) = (3*i_NN)/2 + 1
   integer, parameter, private :: i_Na = i_3N/2 !+ 1!(3*i_NN)/2 != 3*(i_N/2 + 1)/2
   complex,     private :: state_inp(i_KK,0:i_3M-1)!,state_res(i_KK,0:i_3M-1)
   complex,     private :: state_mid(i_KK,0:i_Na)
   REAL(KIND=RKD),   private :: state_phy(i_KK,0:i_3N-1)
!   integer*8,          private :: plan_inp2mid, plan_mid2inp, plan_mid2phy, plan_phy2mid
   TYPE(c_ptr), private :: plan_inp2mid, plan_mid2inp, plan_mid2phy, plan_phy2mid
   complex,     private ::  T(i_KK,0:i_3M-1,0:i_Np-1)
   complex,     private :: Ts(i_KK,0:i_NN1,0:i_Mp-1)

 contains

!------------------------------------------------------------------------
! Setup plans for 2d transforms.  
!------------------------------------------------------------------------
   subroutine tra_precompute
      integer, parameter :: flag=32  !see fftw3.f
      integer :: sgn, n(1), howmany, inembed(1),ouembed(1),istride
      n = (/i_3M/)
      howmany = i_KK
      inembed = (/i_3M/) 
      istride=i_KK
      sgn=1
      plan_inp2mid = fftwf_plan_many_dft(1, n, howmany,&
            state_inp, inembed, istride, &
            1,  state_inp, inembed, istride, 1, sgn, flag) !To Real transform
      sgn=-1
      plan_mid2inp = fftwf_plan_many_dft(1, n, howmany,&
            state_inp, inembed, istride, &
            1,  state_inp, inembed, istride, 1, sgn, flag) !To Real transform

      n = (/i_3N/)
      howmany = i_KK
      inembed = (/i_Na+1/)!i_Na+1 
      ouembed = (/i_3N/) 
      istride=i_KK
      plan_mid2phy = fftwf_plan_many_dft_c2r(1, n, howmany,state_mid, inembed, istride, &
         1,  state_phy, ouembed, istride, 1, flag) !To Real transform
      plan_phy2mid = fftwf_plan_many_dft_r2c(1, n, howmany,state_phy, ouembed, istride, &
         1,  state_mid, inembed, istride, 1, flag) !To Real transform

   end subroutine tra_precompute

!------------------------------------------------------------------------
!  Convert spectral to real space
!------------------------------------------------------------------------
   subroutine tra_spec2phys(s, p)
      type (spec), intent(in)  :: s
      type (phys), intent(out) :: p
      integer :: n,m,mm,ms
      				! for each r_n ...      
      do n = 0,var_N%pH1
         state_inp=0d0
         do m = -(i_M-i_MM),i_MM1
            mm=modulo(m,i_3M)
            ms=modulo(m,i_M)
            state_inp(:,mm)=cmplx(s%Re(:,ms,n),s%Im(:,ms,n))
         end do
         call fftwf_execute_dft(plan_inp2mid,state_inp,state_inp)
         T(:,:,n)=state_inp   
      end do
      call tra_T2Ts()
      do m = 0,var_M%pH1
         state_mid(:,i_NN:)=0d0
         state_mid(:,0:i_NN1)=Ts(:,:,m)
         call fftwf_execute_dft_c2r(plan_mid2phy,state_mid,state_phy)
         p%Re(:,:,m) = state_phy
      end do

   end subroutine tra_spec2phys


   subroutine tra_spectest(s)
      type (spec), intent(inout)  :: s
      type (spec) :: st
      integer :: n,m,mm,ms
                  ! for each r_n ...      
      do n = 0,var_N%pH1
         state_inp=0d0
         do m = -(i_M-i_MM),i_MM1
            mm=modulo(m,i_3M)
            ms=modulo(m,i_M)
            if (n==1) print*,m,mm,ms
            state_inp(:,mm)=cmplx(s%Re(:,ms,n),s%Im(:,ms,n))
         end do
         if (n==1) print*,'BEFORE _-----------------'
         if (n==1) print*, real(state_inp(1,:))
         call fftwf_execute_dft(plan_inp2mid,state_inp,state_inp)
         state_inp=state_inp/(i_3M)
         if (n==1) print*,'AFTER --------------------'
         call fftwf_execute_dft(plan_mid2inp,state_inp,state_inp)
         if (n==1) print*, real(state_inp(1,:))
         do m = - (i_M - i_MM),i_MM1
            mm = modulo(m,i_3M)
            ms = modulo(m,i_M)
            s%Re(:,ms,n)= real(state_inp(:,mm))
            s%Im(:,ms,n)=imag(state_inp(:,mm))
         end do
      end do
   end subroutine tra_spectest


!------------------------------------------------------------------------
!  Convert real to spectral space
!------------------------------------------------------------------------
   subroutine tra_phys2spec(p, s)
      type (spec), intent(out)  :: s
      type (phys), intent(in) :: p
      integer :: n,m,mm,ms
      REAL(KIND=RKD) :: scale_
                  ! scale, FFTW 4.7.2
      scale_ = 1d0 / real(i_3M*i_3N)
      do m = 0,var_M%pH1
         state_phy = scale_ * p%Re(:,:,m)
         call fftwf_execute_dft_r2c(plan_phy2mid,state_phy,state_mid)
         Ts(:,:,m) = state_mid(:,0:i_NN1)
      end do
      call tra_Ts2T
      do n = 0,var_N%pH1
         state_inp = T(:,:,n)
         call fftwf_execute_dft(plan_mid2inp,state_inp,state_inp)
         do m = - (i_M - i_MM),i_MM1
            mm = modulo(m,i_3M)
            ms = modulo(m,i_M)
            s%Re(:,ms,n)= real(state_inp(:,mm))
            s%Im(:,ms,n)=imag(state_inp(:,mm))
         end do
      end do
   end subroutine tra_phys2spec

!------------------------------------------------------------------------
!  Transforms
!------------------------------------------------------------------------

#ifndef _MPI
   subroutine tra_Ts2T()
     _loop_mn_vars
     do n = 0, i_NN1
        do m = 0, i_3M-1   
           T(:,m,n) = Ts(:,n,m)
     _loop_mn_end

   end subroutine tra_Ts2T!tra_Ts2Tvar_tran2spec

#else
   subroutine tra_Ts2T()!var_tran2spec(c,s)
      REAL(KIND=RKD) :: bsend(i_KK,2*i_Np*i_Mp,0:_Np-1)
      REAL(KIND=RKD) :: brecv(i_KK,2*i_Np*i_Mp,0:_Np-1)
      integer :: stp, dst,src, n,m,l,j

      do stp = 0, mpi_sze-1
         src  = modulo(mpi_sze-stp+mpi_rnk, mpi_sze)         
         mpi_tg = stp 
         call mpi_irecv( brecv(1,1,stp), 2*i_KK*(var_M%pH1_(src)+1)*(var_N%pH1+1), &
            mpi_real, src, mpi_tg, mpi_comm_world,  &
            mpi_rq(stp), mpi_er)
      end do
      do stp = 0, mpi_sze-1
         dst  = modulo(stp+mpi_rnk, mpi_sze)
         l = 1
         do n = var_N%ph0_(dst), var_N%ph0_(dst)+var_N%ph1_(dst)
            do m = 0, var_M%pH1
               bsend(:,l,  stp)   =  real(Ts(:,n,m))
               bsend(:,l+1,stp) = imag(Ts(:,n,m))
               l = l + 2
            end do
         end do
         mpi_tg = stp
         call mpi_isend( bsend(1,1,stp), 2*i_KK*(var_N%pH1_(dst)+1)*(var_M%pH1+1),  &
            mpi_real, dst, mpi_tg, mpi_comm_world,  &
            mpi_rq(mpi_sze+stp), mpi_er)
      end do

      do stp = 0, mpi_sze-1
         src  = modulo(mpi_sze-stp+mpi_rnk, mpi_sze)      
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
         do n = 0, var_N%pH1
            do m = var_M%pH0_(src), var_M%pH0_(src)+var_M%pH1_(src)
               do j = 1,i_KK
                  T(j,m,n)= cmplx(brecv(j,l,stp),brecv(j,l+1,stp))
               end do 
               l = l + 2
            end do
         end do
      end do

      do stp = 0, mpi_sze-1
         call mpi_wait( mpi_rq(mpi_sze+stp), mpi_st, mpi_er)
      end do
   end subroutine tra_Ts2T
#endif

#ifndef _MPI
   subroutine tra_T2Ts()
     _loop_mn_vars
      do m = 0, i_3M-1    
        do n = 0, i_NN1
     Ts(:,n,m)=T(:,m,n) 
     _loop_mn_end
     
   end subroutine tra_T2Ts!
   
#else
   subroutine tra_T2Ts()
      REAL(KIND=RKD) :: bsend(i_KK,2*i_Np*i_Mp,0:_Np-1)
      REAL(KIND=RKD) :: brecv(i_KK,2*i_Np*i_Mp,0:_Np-1)
      integer :: stp, dst,src, n,m,l,j,lk
      do stp = 0, mpi_sze-1
         src  = modulo(mpi_sze-stp+mpi_rnk, mpi_sze)         
         mpi_tg = stp 
            call mpi_irecv( brecv(1,1,stp), 2*i_KK*(var_N%pH1_(src)+1)*(var_M%pH1+1), &
            mpi_real, src, mpi_tg, mpi_comm_world,  &
            mpi_rq(stp), mpi_er)
      end do
      do stp = 0, mpi_sze-1
         dst  = modulo(stp+mpi_rnk, mpi_sze)
         l = 1
         lk=l+i_KK1
         do m = var_M%ph0_(dst), var_M%ph0_(dst)+var_M%ph1_(dst)
            do n = 0, var_N%pH1
               bsend(:,l,stp)   =  real(T(:,m,n))
               bsend(:,l+1,stp) = imag(T(:,m,n))
               l=l+2
            end do
         end do
         mpi_tg = stp
         call mpi_isend( bsend(1,1,stp), 2*i_KK*(var_M%ph1_(dst)+1)*(var_N%pH1+1),  &
            mpi_real, dst, mpi_tg, mpi_comm_world,  &
            mpi_rq(mpi_sze+stp), mpi_er)
      end do

      do stp = 0, mpi_sze-1
         src  = modulo(mpi_sze-stp+mpi_rnk, mpi_sze)      
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
         do m = 0, var_M%pH1
            do n = var_N%pH0_(src), var_N%pH0_(src)+var_N%pH1_(src)
               do j=1,i_KK
                  Ts(j,n,m)= cmplx(brecv(j,l,stp),brecv(j,l+1,stp))
               end do 
               l = l + 2
            end do
         end do
      end do

      do stp = 0, mpi_sze-1
         call mpi_wait( mpi_rq(mpi_sze+stp), mpi_st, mpi_er)
      end do
   end subroutine tra_T2Ts
#endif
!*************************************************************************
 end module transform
!*************************************************************************

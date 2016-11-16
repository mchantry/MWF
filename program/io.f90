!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
#include "../parallel.h"
 module io
!**************************************************************************
   use velocity
   use mpif
   use netcdf
   use turb
   implicit none
   save

   character(200)        :: io_statefile
   integer               :: io_save1, io_save2
   integer,     private  :: io_KE, io_HI
   type (mpt),  private  :: c1,c2,c3
   character(200) :: oloc='RUNNING'
 contains
 
!--------------------------------------------------------------------------
!  initialiser fn
!--------------------------------------------------------------------------
   subroutine io_precompute()
      io_statefile = 'state.cdf.in'
      io_save1 = 0
      io_save2 = 0
      io_KE    = 20
      io_hi    = 0      
      if (s_HIS) &
           io_hi    = 30
   end subroutine io_precompute 
 

!--------------------------------------------------------------------------
!  Open files written to every ... steps runtime
!--------------------------------------------------------------------------
   ! subroutine io_openfiles()
   !    character(10), save :: s = 'unknown', a = 'sequential', p='append'
   !    if(mpi_rnk/=0) return
   !    if(io_KE/=0)  open(io_KE,status=s,access=a,position=p, file='vel_energy.dat')
   !    if(io_hi/=0)  open(io_hi,status=s,access=a,position=p, file='vel_history.dat')
   !    s = 'old'
   ! end subroutine io_openfiles

   subroutine io_openfiles()
      character(10), save :: s = 'replace', a = 'sequential', p='append'
      if(mpi_rnk/=0) return
      if(io_KE/=0)  open(io_KE,status=s,position=p, file='vel_energy.dat')
      if(io_hi/=0)  open(io_hi,status=s,position=p, file='vel_history.dat')
      s = 'old'
!      a = 'stream'
   end subroutine io_openfiles


!--------------------------------------------------------------------------
!  Close files written to during runtime
!--------------------------------------------------------------------------
   subroutine io_closefiles()
      if(mpi_rnk/=0) return
      if(io_KE/=0) close(io_KE)
      if(io_hi/=0) close(io_hi)
   end subroutine io_closefiles


   !-------------------------------------------------------------------------
   subroutine io_clk_time(t)
     real, intent(out) :: t
#ifdef _CPUTIME
     call cpu_time(t)
#else
     integer, save :: ct,ctrt,ctmx,ct_=0,wrap=0
     call system_clock(ct,ctrt,ctmx)
     if(ct<ct_) wrap = wrap + 1
     ct_ = ct
     t = (real(ct)+real(ctmx)*wrap)/real(ctrt)
#endif
   end subroutine io_clk_time
    
!-------------------------------------------------------------------------


!--------------------------------------------------------------------------
!  Write to files
!--------------------------------------------------------------------------
   subroutine io_write2files()
      type (phys) :: vp
      if(modulo(tim_step,i_save_rate1)==0) then
!         if (tim_step > 0) then
            call io_save_state()
            call io_save_Uspec()
!         end if
         io_save1 = io_save1+1
      endif

      if(modulo(tim_step,i_save_rate2)==0) then
         call vel_mpt2phys(vel_c,vp)
         if(io_KE/=0) call io_write_energy(vp)
         if(io_hi/=0) call io_write_history(vp)
         if(s_tur) then
            call turb_sample(vp)
            if (modulo(tim_step,i_WT*i_MT*i_save_rate2)==0) &
                 call turb_save_state()
         end if
         io_save2 = io_save2+1
      end if
      
      if(modulo(tim_step,i_save_rate2*50)==0) then
         call io_closefiles()
         call io_openfiles()
      end if

   end subroutine io_write2files


!--------------------------------------------------------------------------
!  Load state - start from previous solution
!--------------------------------------------------------------------------
   subroutine io_load_state()
      integer :: e, f, i, rd
      double precision :: d

      if(mpi_rnk==0)      print*, "Loading: ",io_statefile
      e=nf90_open(io_statefile, nf90_nowrite, f)      
      if(e/=nf90_noerr) then
         if(mpi_rnk==0) open(99, file='PRECOMPUTING')
         if(mpi_rnk==0) close(99, status='delete')
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'state file not found!'
      end if

      e=nf90_get_att(f,nf90_global,'t', d)
      if(d_time<0d0) tim_t = d
      if(mpi_rnk==0 .and. dabs(tim_t-d)>1d-8)  &
         print*,' t    :',d,' --> ', tim_t

      e=nf90_get_att(f,nf90_global,'Re', d)
      if(mpi_rnk==0 .and. dabs(d_Re-d)>1d-8)  &
         print*,' Re   :',d,' --> ', d_Re
      e=nf90_get_att(f,nf90_global,'alpha', d)
      if(mpi_rnk==0 .and. dabs(d_alpha-d)>1d-8)  &
         print*,' alpha:',d,' --> ', d_alpha
      e=nf90_get_att(f,nf90_global,'gamma', d)
      if(mpi_rnk==0 .and. dabs(d_gamma-d)>1d-8)  &
         print*,' gamma:',d,' --> ', d_gamma

      call io_load_mpt(f,'mpt',vel_c)

      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

   end subroutine io_load_state


!--------------------------------------------------------------------------
!  Load coll variable
!--------------------------------------------------------------------------
   subroutine io_load_mpt(f,nm, a)
      integer,          intent(in)  :: f
      character(*),     intent(in)  :: nm
      type (mpt),      intent(out) :: a
      integer :: N1,M1, mm_,mm,nn,e,i
      integer :: K__, M__, N__,K1,K2
      integer :: m,n
      logical :: error=.false.
          
      e=nf90_inq_varid(f,nm, i)
      if(e/=nf90_noerr)  print*, 'Field '//nm//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_get_att(f,i, 'K',  K__)
      e=nf90_get_att(f,i, 'MM',  M__)
      e=nf90_get_att(f,i, 'NN',  N__)
      if(K__ /=i_K)  error=.true.
      if(M__ /=i_MM) error=.true.
      if(N__ /=i_NN) error=.true.

      if(mpi_rnk==0) then
         if(K__ /=i_K)  print*, nm, ' K :', K__, ' --> ',i_K 
         if(M__ /=i_MM)  print*, nm, ' MM :', M__, ' --> ',i_MM 
         if(N__ /=i_NN)  print*, nm, ' NN :', N__,' --> ',i_NN 
      end if
      if (error) then
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'wrong state size'
      end if
      
      a%Re = 0d0
      a%Im = 0d0

      e=nf90_get_var(f,i, a%Re(:,:,:),start=(/1,1,var_N%pH0+1,1/))
      e=nf90_get_var(f,i, a%Im(:,:,:),start=(/1,1,var_N%pH0+1,2/))
      
    end subroutine io_load_mpt

   subroutine io_load_mpt_old(f,nm, a)
      integer,          intent(in)  :: f
      character(*),     intent(in)  :: nm
      type (mpt),      intent(out) :: a
      integer :: N1,M1, mm_,mm,nn,e,i
      integer :: K__, M__, N__,K1,K2
      integer :: m,n
          
      e=nf90_inq_varid(f,nm, i)
      if(e/=nf90_noerr)  print*, 'Field '//nm//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_get_att(f,i, 'K',  K__)
      e=nf90_get_att(f,i, 'MM',  M__)
      e=nf90_get_att(f,i, 'NN',  N__)
      if(mpi_rnk==0) then
         if(K__ /=i_K)  print*, nm, ' K :', K__, ' --> ',i_K 
         if(M__ /=i_MM)  print*, nm, ' MM :', M__, ' --> ',i_MM 
         if(N__ /=i_NN)  print*, nm, ' NN :', N__,' --> ',i_NN 
      end if

      a%Re = 0d0
      a%Im = 0d0

      N1 = min(N__,i_NN)-1
      M1 = min(M__,i_MM)-1
      K1 = min(i_K0,(K__+1)/2)
      K2 = (K__+1)/2+1
      do n = 0, N__-1 
         if(n<var_N%pH0 .or. n>var_N%pH0+var_N%pH1)  cycle
         nn = n - var_N%pH0
!         print*,n,mpi_rnk,nn
         do m = -M__,M__
            if (abs(m) > i_MM1) cycle
            mm_ = modulo(m,2*(M__-1))
            mm  = modulo(m,i_M)
            if (i_K == K__) then
               e=nf90_get_var(f,i, a%Re(1:i_K,mm,nn),start=(/1,mm_+1,n+1,1/))
               e=nf90_get_var(f,i, a%Im(1:i_K,mm,nn),start=(/1,mm_+1,n+1,2/))
            else 
               e=nf90_get_var(f,i, a%Re(1:K1,mm,nn)            ,start=(/1,mm_+1,n+1,1/))
               e=nf90_get_var(f,i, a%Re(i_K0+1:i_K0+K1-1,mm,nn),start=(/K2,mm_+1,n+1,1/))
               e=nf90_get_var(f,i, a%Im(1:K1,mm,nn)            ,start=(/1,mm_+1,n+1,2/))
               e=nf90_get_var(f,i, a%Im(i_K0+1:i_K0+K1-1,mm,nn),start=(/K2,mm_+1,n+1,2/))
            end if
         end do
      end do
    end subroutine io_load_mpt_old

!--------------------------------------------------------------------------
!  Save state
!--------------------------------------------------------------------------
   subroutine io_save_state()
      character(4) :: cnum
      integer :: e, f
      integer :: md,nd,kd, ReImd, dims(4)
      integer :: ss 

      write(cnum,'(I4.4)') io_save1

      if(mpi_rnk==0) then
         print*, ' saving state'//cnum//'  t=', tim_t
         e=nf90_create('state'//cnum//'.cdf.dat', nf90_clobber, f)

         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         e=nf90_put_att(f, nf90_global, 'gamma', d_gamma)

         e=nf90_def_dim(f, 'N', i_NN, nd)
         e=nf90_def_dim(f, 'M', i_M, md)
         e=nf90_def_dim(f, 'K', i_K, kd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         dims = (/kd,md,nd,ReImd/)
         call io_define_mpt(f, 'mpt', dims, ss)

         e=nf90_enddef(f)
      end if

      call io_save_mpt(f,ss, vel_c)

      if(mpi_rnk==0)  &
         e=nf90_close(f)

   end subroutine io_save_state

!--------------------------------------------------------------------------
!  Save coll variable
!--------------------------------------------------------------------------
   subroutine io_define_mpt(f,nm,dims, id)
      integer,      intent(in) :: f, dims(4)
      character(*), intent(in) :: nm
      integer, intent(out) :: id
      integer :: e
      e=nf90_def_var(f, nm, nf90_double, dims, id)
      e=nf90_put_att(f, id,  'K', i_K)      
      e=nf90_put_att(f, id,  'MM', i_MM)
      e=nf90_put_att(f, id,  'NN', i_NN)
   end subroutine io_define_mpt

   subroutine io_save_mpt(f,id,a)
      integer,     intent(in) :: f, id
      type (mpt), intent(in) :: a
      integer :: e
      
#ifndef _MPI
      e=nf90_put_var(f,id,a%Re(1:i_K,0:i_M1,0:i_NN1), start=(/1,1,1,1/))
      e=nf90_put_var(f,id,a%Im(1:i_K,0:i_M1,0:i_NN1), start=(/1,1,1,2/))

#else
      integer :: r, pN0,pN1

      if(mpi_rnk==0) then
         e=nf90_put_var(f,id,a%Re(1:i_K,0:i_M1,0:var_N%pH1), start=(/1,1,1,1/))
         e=nf90_put_var(f,id,a%Im(1:i_K,0:i_M1,0:var_N%pH1), start=(/1,1,1,2/))         
         do r = 1, mpi_sze-1
            pN0 = var_N%pH0_(r)
            pN1 = var_N%pH1_(r)
            mpi_tg = r
            call mpi_recv( c1%Re(1,0,0), i_M*(pN1+1)*i_K, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            call mpi_recv( c1%Im(1,0,0), i_M*(pN1+1)*i_K, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            e=nf90_put_var(f,id,c1%Re(1:i_K,0:i_M1,0:pN1), start=(/1,1,pN0+1,1/))
            e=nf90_put_var(f,id,c1%Im(1:i_K,0:i_M1,0:pN1), start=(/1,1,pN0+1,2/))
         end do
      else
         mpi_tg = mpi_rnk
         call mpi_send( a%Re(1,0,0), i_M*(var_N%pH1+1)*i_K, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%Im(1,0,0), i_M*(var_N%pH1+1)*i_K, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
      end if
#endif      
   end subroutine io_save_mpt

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  save spectrum !var_mpt2spec
!--------------------------------------------------------------------------
   subroutine io_save_Uspec()
     
     type(spec) :: tmp
     double precision :: n_(0:i_NN1), m_(0:i_MM1)
     double precision :: n__(0:i_NN1), m__(0:i_MM1)
     double precision :: d, dRe, dIm
     character(4) :: cnum
     integer :: i
     _loop_kmn_vars
10   format(i4,1e20.12)
     
     call var_mpt2spec(vel_c,tmp)
     n_ = 0d0
     m_ = 0d0
     _loop_mn_begin
     mm = m
     if (m > i_MM1) mm = i_M - m
     do k=1,i_KK
        dRe=tmp%Re(k,m,n)
        dIm=tmp%Im(k,m,n) 
        d = dsqrt(dRe*dRe+dIm*dIm)
        n_(nn)  = max(d, n_(nn))
        m_(mm)  = max(d, m_(mm))
     end do
     _loop_mn_end
     
#ifdef _MPI
     call mpi_barrier(mpi_comm_world, mpi_er)
     call mpi_allreduce(n_(0), n__(0), i_NN, mpi_double_precision,  &
          mpi_max, mpi_comm_world, mpi_er)
     n_ = n__
     call mpi_allreduce(m_(0), m__(0), i_MM, mpi_double_precision,  &
          mpi_max, mpi_comm_world, mpi_er)
     m_ = m__
#endif
     if(mpi_rnk/=0) return
     write(cnum,'(I4.4)') io_save1
     open(11, status='unknown', file='vel_spec'//cnum//'.dat')
     write(11,*) '# t = ', tim_t
     write(11,*) '# m'
     do i = 0, i_MM1
        write(11,10) i, m_(i)      
     end do
     write(11,*)
     write(11,*) '# n'
     do i = 0, i_NN1
        write(11,10) i, n_(i)      
     end do
     close(11)
     
   end subroutine io_save_Uspec


!--------------------------------------------------------------------------
!  save spectrum !var_mpt2spec
!--------------------------------------------------------------------------
   subroutine io_save_spectrum()
      double precision :: n_(0:i_NN1), m_(0:i_MM1)
      double precision :: n__(0:i_NN1), m__(0:i_MM1)
      double precision :: d, dRe, dIm
      character(4) :: cnum
      integer :: i
      _loop_kmn_vars
   10 format(i4,1e20.12)
      
      n_ = 0d0
      m_ = 0d0
      _loop_kmn_begin
      mm = m
      if (m > i_MM1) mm = i_M - m
!      print*,mpi_rnk,k,mm,nn
      dRe=vel_c%Re(k,m,n)
      dIm=vel_c%Im(k,m,n) 
      d = dsqrt(dRe*dRe+dIm*dIm)
      n_(nn)  = max(d, n_(nn))
      m_(mm)  = max(d, m_(mm))
      _loop_kmn_end
      
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
      call mpi_allreduce(n_(0), n__(0), i_NN, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      n_ = n__
      call mpi_allreduce(m_(0), m__(0), i_MM, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      m_ = m__
#endif
      if(mpi_rnk/=0) return
      write(cnum,'(I4.4)') io_save1
      open(11, status='unknown', file='vel_spec'//cnum//'.dat')
      write(11,*) '# t = ', tim_t
      write(11,*) '# m'
      do i = 0, i_MM1
         write(11,10) i, m_(i)      
      end do
      write(11,*)
      write(11,*) '# n'
      do i = 0, i_NN1
         write(11,10) i, n_(i)      
      end do
      close(11)

   end subroutine io_save_spectrum

 
!--------------------------------------------------------------------------
!  write to energy file
!--------------------------------------------------------------------------
   subroutine io_write_energy(sp)
      type (phys), intent(in) :: sp
      double precision :: E

      call vel_energy(sp,E)
      
      if(mpi_rnk/=0) return
      write(io_KE,'(2e20.12)')  tim_t, E
      
      if(E>d_minE .or. tim_t<20d0) return
      print*, 'io_write_energy: Relaminarised!'
      open(99,file=oloc)
      close(99, status='delete')

   end subroutine io_write_energy

!--------------------------------------------------------------------------
!  write to history file
!--------------------------------------------------------------------------
   subroutine io_write_history(sp)
      type (phys), intent(in) :: sp
      double precision :: H(4,i_H)
      integer :: i

      if(mpi_rnk/=0) return
      call vel_history(sp,0d0,H)
      do i=1,i_H
         write(io_HI,'(6e20.12)')  tim_t, 0d0,H(1,i),H(2,i),H(3,i),H(4,i)
      end do

   end subroutine io_write_history


   Subroutine io_writeVTK_xz(V,y,cnum,xs,zs)
  
     type(phys),intent(in) :: V
     double precision, intent(in) :: y
     integer, optional :: xs,zs
     character(4) :: cnum
     integer :: io,ix,iz,M_,N_
     character(10) :: s = 'unknown', a = 'sequential'
     
     if (.not. present(zs)) then
        xs=1
        zs=1
     end if

     M_=i_3M/xs
     N_=i_3N/zs

     if(mpi_rnk==0) then
     
        open(io,status=s,access=a, file='XZ'//cnum//'.vtk')
        
        write(io,'(A)') "# vtk DataFile Version 3.0"
        write(io,'(A)') "M x 2 x N"
        write(io,'(A)') "ASCII"
        write(io,'(A)') "DATASET RECTILINEAR_GRID"
        write(io,'(A,3i5)') "DIMENSIONS",M_+1,1,N_+1
        
        write(io,'(A,i5,A)') "X_COORDINATES ",M_+1, " FLOAT"
        do ix=0,i_3M,xs
           write(io,'(f9.5)') d_Lx*dble(ix)/(i_3M)
        end do
        
        write(io,'(A,i5,A)') "Y_COORDINATES ",1, " FLOAT"
        write(io,'(f9.5)') 0d0
        
        write(io,'(A,i5,A)') "Z_COORDINATES ",N_+1, " FLOAT"
        do ix=0,i_3N,zs
           write(io,'(f9.5)') d_Lz*dble(ix)/(i_3N)
        end do
        
        write(io,*) 
        
        write(io,'(A,i10)') "POINT_DATA",(N_+1)*(M_+1) 
        write(io,'(A)') "VECTORS velocity FLOAT"

        do iz=0,i_3N-1,zs
           do ix=0,i_3M-1,xs
              write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
           end do
           ix=0
           write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
        end do
        iz=0
        do ix=0,i_3M-1,xs
           write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
        end do
        ix=0
        write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
        
        close(io)
     END if
  
   end subroutine io_writeVTK_xz

   Subroutine io_writeVTK_txz(V,cnum,xs,zs,ct)
  
     type(phys),intent(in) :: V
     double precision, intent(in) :: ct
     integer, optional :: xs,zs
     character(4) :: cnum
     integer :: io,ix,iz,M_,N_,t
     character(10) :: s = 'unknown', a = 'sequential'
     double precision :: ke

     if (.not. present(zs)) then
        xs=1
        zs=1
     end if

     M_=i_3M/xs
     N_=i_3N/zs
     t=0

     if(mpi_rnk==0) then
     
        open(io,status=s,access=a, file='TXZ'//cnum//'.vtk')
        
        write(io,'(A)') "# vtk DataFile Version 3.0"
        write(io,'(A)') "M x 2 x N"
        write(io,'(A)') "ASCII"
        write(io,'(A)') "DATASET RECTILINEAR_GRID"
        write(io,'(A,3i5)') "DIMENSIONS",M_+1,1,N_+1
        
        write(io,'(A,i5,A)') "X_COORDINATES ",M_+1, " FLOAT"
        do ix=0,i_3M,xs
           write(io,'(f9.5)') d_Lx*dble(ix)/(i_3M)
        end do
        
        write(io,'(A,i5,A)') "Y_COORDINATES ",1, " FLOAT"
        write(io,'(f9.5)') 0d0
        
        write(io,'(A,i5,A)') "Z_COORDINATES ",N_+1, " FLOAT"
        do ix=0,i_3N,zs
           write(io,'(f9.5)') d_Lz*dble(ix)/(i_3N)
        end do
        
        write(io,*) 
        
        write(io,'(A,i10)') "POINT_DATA",(N_+1)*(M_+1) 
        write(io,'(A)') "SCALAR turb FLOAT"

        do iz=0,i_3N-1,zs
           do ix=0,i_3M-1,xs
              ke = epos(V%Re(:,iz,ix))
              if (ke > ct) t = t + 1
              write(io,'(3e13.5)') ke
           end do
           ix=0
           write(io,'(3e13.5)') epos(V%Re(:,iz,ix))
        end do
        iz=0
        do ix=0,i_3M-1,xs
           write(io,'(3e13.5)') epos(V%Re(:,iz,ix))
        end do
        ix=0
        write(io,'(3e13.5)') epos(V%Re(:,iz,ix))
        
        close(io)
     END if
     print*,"Turbulent fraction :", dble(t)/(N_+1*M_+1)
  
   end subroutine io_writeVTK_txz


   Subroutine io_writeVTK(V,cnum,xs,zs)
  
     type(phys),intent(in) :: V
     integer, optional :: xs,zs
     character(4) :: cnum
     integer :: io,ix,iy,iz,M_,N_
     double precision :: y
     character(10) :: s = 'unknown', a = 'sequential'
     
     if (.not. present(zs)) then
        xs=1
        zs=1
     end if
     
     M_=i_3M/xs
     N_=i_3N/zs
     
     if(mpi_rnk==0) then
        
        open(io,status=s,access=a, file='XYZ'//cnum//'.vtk')
        
        write(io,'(A)') "# vtk DataFile Version 3.0"
        write(io,'(A)') "Lx x 2 x Lz"
        write(io,'(A)') "ASCII"
        write(io,'(A)') "DATASET RECTILINEAR_GRID"
        write(io,'(A,3i5)') "DIMENSIONS",M_+1,i_Ny,N_+1
        
        write(io,'(A,i5,A)') "X_COORDINATES ",M_+1, " FLOAT"
        do ix=0,i_3M,xs
           write(io,'(f9.5)') d_Lx*dble(ix)/(i_3M)
        end do
        
        write(io,'(A,i5,A)') "Y_COORDINATES ",i_Ny, " FLOAT"
        do iy=0,i_Ny-1
           write(io,'(f9.5)') -1d0 + 2d0*dble(iy)/(i_Ny-1)
        end do
        
        
        write(io,'(A,i5,A)') "Z_COORDINATES ",N_+1, " FLOAT"
        do ix=0,i_3N,zs
           write(io,'(f9.5)') d_Lz*dble(ix)/(i_3N)
        end do
        
        write(io,*) 
        
        write(io,'(A,i10)') "POINT_DATA",(N_+1)*(M_+1)*i_Ny
        write(io,'(A)') "VECTORS velocity FLOAT"
        
        do iz=0,i_3N-1,zs
           do iy=0,i_Ny-1
              y=-1d0 + 2d0*dble(iy)/(i_Ny-1)
              do ix=0,i_3M-1,xs
                 write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
              end do
              ix=0
              write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
           end do
        end do
        iz=0
        do iy=0,i_Ny-1
           y=-1d0 + 2d0*dble(iy)/(i_Ny-1)
           do ix=0,i_3M-1,xs
              write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
           end do
           ix=0
           write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
        end do
        close(io)
     END if
  
   end subroutine io_writeVTK
   

!**************************************************************************
 end module io
!**************************************************************************


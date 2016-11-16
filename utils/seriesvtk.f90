#include "../parallel.h"
program state2vtk

  use io
  use mpif

  implicit none

  character(4) :: cnum
  integer :: i, i0, i1,nx,nz
  type(phys) :: pv
  
  print*, 'initialising...'
  call mpi_precompute()
  call par_precompute()
  call var_precompute()
  call tra_precompute()
  call vel_precompute()
  call  io_precompute()
  
  if (mpi_sze == 1) then
     print*,"Reading Initial and final state number"
     read(*,*) i0,i1
     print*, "Save every nx,nz points"
     read(*,*) nx,nz
  else
     open(99,status='old',file='series.in')
     read(99,*) i0,i1
     read(99,*) nx,nz
     close(99)
  end if
  do i=i0,i1
     write(cnum,'(I4)') i
     if(i<1000) cnum(1:1) = '0'
     if(i< 100) cnum(2:2) = '0'
     if(i<  10) cnum(3:3) = '0'
     print*, cnum
     io_statefile='state'//cnum//'.cdf.dat'
     call io_load_state()
     call vel_mpt2phys(vel_c,pv)
!     if (mpi_sze == 1) then
        call io_writevtk_xz(pv,0d0,cnum,nx,nz)
 !    else
 !       call iom_writevtk_xz(pv,0d0,cnum,nx,nz)
!     end if
  end do

end program state2vtk

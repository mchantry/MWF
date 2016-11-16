#include "../parallel.h"
program state2vtk

  use io

  implicit none

  character(4) :: cnum
  integer :: i, i0, i1,nx,nz
  type(phys) :: pv
  double precision :: ct
!  _loop_mn_vars

  print*, 'initialising...'
  call mpi_precompute()
  call par_precompute()
  call var_precompute()
  call tra_precompute()
  call vel_precompute()
  call  io_precompute()
  
  print*,"Reading Initial and final state number"
  read(*,*) i0,i1
  print*, "Save every nx,nz points"
  read(*,*) nx,nz
  print*, "Critical level for turbulent"
  read(*,*) ct
  do i=i0,i1
     write(cnum,'(I4)') i
     if(i<1000) cnum(1:1) = '0'
     if(i< 100) cnum(2:2) = '0'
     if(i<  10) cnum(3:3) = '0'
     print*, cnum
     io_statefile='state'//cnum//'.cdf.dat'
     call io_load_state()
     call vel_mpt2phys(vel_c,pv)
     call io_writevtk_txz(pv,cnum,nx,nz,ct)
  end do

end program state2vtk

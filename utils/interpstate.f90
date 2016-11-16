#include "../parallel.h"
program state2vtk

  use io

  implicit none

  type(mpt) :: vel_1
  double precision :: rscl
9 format(A)

  print*, 'initialising...'
  call mpi_precompute()
  call par_precompute()
  call var_precompute()
  call tra_precompute()
  call vel_precompute()
  call  io_precompute()
  
  if(mpi_sze/=1) stop 'set _Np 1'

  print*,"Reading state1 name"
  read(*,9) io_statefile
  call io_load_state()  
  call var_mpt_copy(vel_c,vel_1)
  print*,"Reading state2 name"
  read(*,9) io_statefile
  call io_load_state()
  print*, "Read rescale factor: u := (1-d) u1 + d u2 ;  0 -> state1,  1 -> state2"
  read(*,*) rscl

  vel_c%Re = (1d0-rscl)*vel_1%Re + rscl*vel_c%Re
  vel_c%Im = (1d0-rscl)*vel_1%Im + rscl*vel_c%Im

  call io_save_state()

end program state2vtk

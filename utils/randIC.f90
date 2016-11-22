!*************************************************************************
! main.f90 (executable)
!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
PROGRAM double
!*************************************************************************
  use variables
  use mpif
  use velocity
  use io
  implicit none

  double precision :: scl=1d-8

  print*, 'initialising...'
  call mpi_precompute()
  call par_precompute()
  call var_precompute()
  call tra_precompute()
  call  io_precompute()

  call var_randmpt(vel_c)
  call var_maskmpt(vel_c)
  call io_save_state()
!*************************************************************************
 END PROGRAM double
!*************************************************************************


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

  double precision :: scl=1d-12

  call mpi_precompute()
  if(mpi_rnk==0)  print*, 'initialising...'
  call par_precompute()
  call var_precompute()
!  call tra_precompute()
  call  io_precompute()
  if(mpi_rnk==0)  print*, 'initialised'

  call io_load_state()
  if(mpi_rnk==0)  print*, 'loaded state'
  call var_quadX(vel_c)
  if(mpi_rnk==0)  print*, 'shifted wavenumbers'
  call var_randadd(vel_c,scl)
  if(mpi_rnk==0)  print*, 'added noise'
  call io_save_state()
  if(mpi_rnk==0)  print*, 'saved state'
  call io_save_Uspec()
  if(mpi_rnk==0)  print*, 'saved spectrum'

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
      call mpi_finalize(mpi_er)
#endif

!*************************************************************************
 END PROGRAM double
!*************************************************************************


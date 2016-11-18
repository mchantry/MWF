#include "../parallel.h"
program fftwtest

  use variables
  use velocity

  implicit none

  type(spec) :: tsp,csp
  type(phys) :: phy,cphy
  double precision :: rscl
9 format(A)

  print*, 'initialising...'
  call mpi_precompute()
  call par_precompute()
  call var_precompute()
  call tra_precompute()
  
  if(mpi_sze/=1) stop 'set _Np 1'

  call var_randspec(tsp)
  call var_spec_copy(tsp,csp)
  call tra_spec2phys(tsp,phy)

  cphy%Re = phy%Re

  call tra_phys2spec(phy,tsp)

  call tra_spec2phys(tsp,phy)

  call var_spec_sub(tsp,csp)

  cphy%Re = cphy%Re - phy%Re 

  print*,maxval(csp%Re)

  print*,maxval(cphy%Re)
  print*,maxval(tsp%Re)


end program fftwtest

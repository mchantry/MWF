#include "../parallel.h"
 module parameters
!***************************************************************************
   implicit none
   save

   INTEGER, PARAMETER :: RKD = SELECTED_REAL_KIND(6,37)

   REAL(KIND=RKD)            ::  d_Re         = 200d0!67.8d0!875d0
   
   !NUMBER OF MODES TO USE IE HIGHEST WAVENUMBER + 1
   integer,          parameter :: i_MM	        = 32!4096!512!2048!4096
   integer,          parameter :: i_NN          = 32!1024!512!2048!512!2048!4096
   integer,          parameter :: i_K0          = 4
   
   REAL(KIND=RKD), parameter :: d_PI          = 3.1415926535897931d0
   REAL(KIND=RKD), parameter :: d_Lx          = 16d0!8192d0!1024d0!4096d0!960d0!16d0
   REAL(KIND=RKD), parameter :: d_Lz          = 16d0!2048d0!1024d0!4096d0!1024d0!4096d0!960d0!1024d0
   REAL(KIND=RKD), parameter :: d_alpha       = 2d0*d_PI/d_Lx!0.5d0
   REAL(KIND=RKD), parameter :: d_gamma       = 2d0*d_PI/d_Lz!0.5d0

   logical,          parameter :: s_reflect     = .false.!.TRUE.!.FALSE. 
   logical,          parameter :: s_uvreflect   = .FALSE.
   logical,          parameter :: s_wreflect    = .false.

   integer, parameter          :: i_maxtstep    = 1d8
   integer, parameter          :: i_save_rate1  = 5000!12500!1d8!50000!100000!5000
   integer, parameter          :: i_save_rate2  = 50!25!25!50 
   REAL(KIND=RKD), parameter :: d_maxt        = -1d0
   REAL(KIND=RKD), parameter :: d_cpuhours    = 9.7d0!1.7d0!19.6d0!0.98d0!1d99 !90d0
   REAL(KIND=RKD), parameter :: d_time        = -1d0 !0d0
   REAL(KIND=RKD), parameter :: d_thdeg       = 0d0!24d0
   REAL(KIND=RKD), parameter :: d_dt          = 0.01d0!0.02d0
   REAL(KIND=RKD), parameter :: d_minE        = 1d-6

   REAL(KIND=RKD), parameter :: d_HYPO        = 0d0!1d-6 !1d-01
   integer, parameter          :: i_PHYPO       = 2
   REAL(KIND=RKD), parameter :: d_drag        = 1d-2 !1d-01
   REAL(KIND=RKD), parameter :: d_vdrag       = 0d0
   logical,          parameter :: s_dragall     = .true.!.false.

   logical, parameter          :: s_HIS         = .FALSE. !HISTORY DATA
   integer, parameter          :: i_H           = 32
   logical, parameter          :: s_TUR         = .TRUE.! TURBULENT DATA
   integer, parameter          :: i_DM          = 6!60!!0!12!60!48!6!6!3!96 !Single point in M
   integer, parameter          :: i_DN          = 6!6!12!6!3!9  !equal to dz of 1.875
   integer, parameter          :: i_MT          = 10!Number of times to average over
   integer, parameter          :: i_WT          = 10!20  ! Gaps of MT to make between output
   !---------------------------------------------------------------------------
   !  Fixed parameters
   !---------------------------------------------------------------------------
   integer,	     parameter :: i_N           = 2*(i_NN-1)
   integer,          parameter :: i_K           = 2*i_K0-1
   integer,          parameter :: i_KK          = 3*i_K0-1
   integer,          parameter :: i_M           = 2*(i_MM-1)
   integer,          parameter :: i_Ny          = 15
   REAL(KIND=RKD), parameter :: d_beta        = d_PI/2d0
   REAL(KIND=RKD), parameter :: d_theta       = d_thdeg/180d0*d_PI

   integer,          parameter :: i_3M          = 3*i_MM
   integer,          parameter :: i_3N          = 3*i_NN
   integer,          parameter :: i_SM          = i_3M/i_DM
   integer,          parameter :: i_SN          = i_3N/i_DN
   
   integer,          parameter :: i_Mp          = (_Np + i_3M - 1)/_Np
   integer,          parameter :: i_Np          = (_Np + i_NN- 1)/_Np
   integer,          parameter :: i_SMp         = (_Np + i_SM - 1)/_Np

   integer,          parameter :: i_NN1  = i_NN-1
   integer,          parameter :: i_N1  = i_N-1
   integer,          parameter :: i_M1  = i_M-1
   integer,          parameter :: i_MM1  = i_MM-1
   integer,          parameter :: i_KK1  = i_KK-1
   integer,          parameter :: i_K1  = i_K0-1

   REAL(KIND=RKD) :: tim_t
   integer          :: tim_step
!---------------------------------------------------------------------------
!  check params are ok
!---------------------------------------------------------------------------

contains

   subroutine par_precompute()
      integer :: itmp
      if (d_time > 0d0) then
         tim_t=d_time
      else
         tim_t=0d0
      end if
      tim_step=0
      itmp=mod(i_M,_Np)
!      if(itmp /= 0) stop 'mpi_precompute: incorrect num procs, M'
      itmp=mod(i_NN,_Np)
!      if(itmp /= 0) stop 'mpi_precompute: incorrect num procs, N'

!      print*, "i_Mp : ",i_Mp
!      print*, "i_Np : ",i_Np


   end subroutine par_precompute
 

!***************************************************************************
 end module parameters
!***************************************************************************

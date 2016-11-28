#include "../parallel.h"
 module parameters
!***************************************************************************
   use rp_emulator
   implicit none
   save

!   INTEGER, PARAMETER :: RKD = KIND(1.0)!SELECTED_REAL_KIND(6,37)
   INTEGER, PARAMETER :: RKD = KIND(1.0d0)!SELECTED_REAL_KIND(13,300)
   !Default number of bits
   INTEGER, parameter :: i_NB = 15 

   REAL(KIND=RKD)            ::  td_Re         = 200d0
   !NUMBER OF MODES TO USE IE HIGHEST WAVENUMBER + 1
   integer,          parameter :: i_MM	        = 32
   integer,          parameter :: i_NN          = 32
   integer,          parameter :: i_K0          = 4
   
   REAL(KIND=RKD), parameter :: td_PI          = 3.1415926535897931d0
   REAL(KIND=RKD), parameter :: td_Lx          = 16d0
   REAL(KIND=RKD), parameter :: td_Lz          = 16d0
   REAL(KIND=RKD), parameter :: td_alpha       = 2d0*td_PI/td_Lx
   REAL(KIND=RKD), parameter :: td_gamma       = 2d0*td_PI/td_Lz

   logical,          parameter :: s_reflect     = .false.!.TRUE.!.FALSE. 
   logical,          parameter :: s_uvreflect   = .FALSE.
   logical,          parameter :: s_wreflect    = .false.

   integer, parameter          :: i_maxtstep    = 150000!1d8
   integer, parameter          :: i_save_rate1  = 5000!12500!1d8!50000!100000!5000
   integer, parameter          :: i_save_rate2  = 50!25!25!50 
   REAL(KIND=RKD), parameter :: td_maxt        = -1d0
   REAL(KIND=RKD), parameter :: d_cpuhours    = 9.7d0!1.7d0!19.6d0!0.98d0!1d99 !90d0
   REAL(KIND=RKD), parameter :: td_time        = -1d0 !0d0
   REAL(KIND=RKD), parameter :: td_thdeg       = 0d0!24d0
   REAL(KIND=RKD), parameter :: td_dt          = 0.01d0!0.02d0
   REAL(KIND=RKD), parameter :: td_minE        = 1d-6

   REAL(KIND=RKD), parameter :: td_HYPO        = 0d0!1d-6 !1d-01
   integer, parameter          :: i_PHYPO       = 2
   REAL(KIND=RKD), parameter :: td_drag        = 1d-2 !1d-01
   REAL(KIND=RKD), parameter :: td_vdrag       = 0d0
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
   REAL(KIND=RKD),   parameter :: td_beta        = td_PI/2d0
   REAL(KIND=RKD),   parameter :: td_theta       = td_thdeg/180d0*td_PI

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

   type(rpe_var) :: tim_t
   integer          :: tim_step

! RPE versions
   type(rpe_var)   :: rp0             
   type(rpe_var)   :: rp2             
   type(rpe_var)   :: d_Re            
   type(rpe_var)   :: d_PI            
   type(rpe_var)   :: d_Lx            
   type(rpe_var)   :: d_Lz            
   type(rpe_var)   :: d_alpha         
   type(rpe_var)   :: d_gamma         
   type(rpe_var)   :: d_maxt
   type(rpe_var)   :: d_time          
   type(rpe_var)   :: d_thdeg         
   type(rpe_var)   :: d_dt            
   type(rpe_var)   :: d_minE          
   type(rpe_var)   :: d_HYPO          
   type(rpe_var)   :: d_drag          
   type(rpe_var)   :: d_vdrag         
   type(rpe_var)   :: d_beta          
   type(rpe_var)   :: d_theta         
   
!---------------------------------------------------------------------------
!  check params are ok
!---------------------------------------------------------------------------

contains

   subroutine par_precompute()
      integer :: itmp
      if (d_time > 0d0) then
         tim_t=d_time!rpe_literal(d_time,i_nb)
      else
         tim_t=rpe_literal(0.0,i_nb)
      end if
      tim_step=0
      itmp=mod(i_M,_Np)
      itmp=mod(i_NN,_Np)

      RPE_ACTIVE = .TRUE.
      RPE_DEFAULT_SBITS = i_nb
      RPE_IEEE_HALF = .TRUE.

      rp0 = rpe_literal(0.0,i_nb)            
      rp2 = rpe_literal(2.0,i_nb)            
      d_Re = rpe_literal(td_Re,i_nb)
      d_PI = rpe_literal(td_PI,i_nb)            
      d_Lx = rpe_literal(td_Lx,i_nb)            
      d_Lz = rpe_literal(td_Lz,i_nb)            
      d_alpha = rpe_literal(td_alpha,i_nb)         
      d_gamma = rpe_literal(td_gamma,i_nb)         
      d_time = rpe_literal(td_time,i_nb)
      d_thdeg = rpe_literal(td_thdeg,i_nb)         
      d_dt = rpe_literal(td_dt,i_nb)            
      d_minE = rpe_literal(td_minE,i_nb)          
      d_HYPO = rpe_literal(td_HYPO,i_nb)          
      d_drag = rpe_literal(td_drag,i_nb)          
      d_vdrag = rpe_literal(td_vdrag,i_nb)         
      d_beta = rpe_literal(td_beta,i_nb)          
      d_theta = rpe_literal(td_theta,i_nb)         
      
   end subroutine par_precompute
 

!***************************************************************************
 end module parameters
!***************************************************************************

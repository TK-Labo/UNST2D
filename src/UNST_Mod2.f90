! Define Global Variables
! Coded by TK Labo

module d1riv_globals_mod

! ===========================================
!  1d-riv (development) - only single UNST2D
! ===========================================

! -- common --
!----------------
! constant param
!----------------
real(8), parameter :: rivgg = 9.8d0
real(8), parameter :: rivpi = 3.14d0
real(8), parameter :: rivth = 1.0d-3
real(8), parameter :: thd_1d = 1.0d-4
real(8), parameter :: tha_1d = 1.0d-6
real(8), parameter :: tha2_1d = 1.0d-4

!------------
! time param
!------------
integer spout  ! dispout in spin up
real(8) unsttime_r
integer d1_spin_upn
real(8) rivdt

!-------------
! num of data
!-------------
integer ndan  ! num of cross-section
integer nriv  ! num of river
integer, allocatable :: riv_ndan(:)

!---------------
! cross section
!---------------
real(8), allocatable :: h_1d(:), ho_1d(:)  ! water head
real(8), allocatable :: vv_1d(:)           ! velocity
real(8), allocatable :: q_1d(:), qo_1d(:)  ! volume
real(8), allocatable :: a_1d(:), r_1d(:)   ! area, radius
real(8), allocatable :: b_1d(:), rn_1d(:)  ! width, roughness

real(8), allocatable :: h_1dmax(:), vv_1dmax(:)  ! max h, vv

! -- fractional step -- TK.Labo
real(8) dz
real(8) :: first_depth_1d
! main array
integer, allocatable :: ntype_1d(:)  ! (-*) -> (*)
real(8), allocatable :: kp_1d(:), dx_1d(:)  ! kp, dx
real(8), allocatable :: depth_1d(:)  ! depth
real(8), allocatable :: rbed_1d(:), rzmax_1d(:)  ! river bed & head
! coords
integer :: composit_rn  ! composite roughness cofficient flag (1:activate)
integer, allocatable :: dan_record(:)  ! num of record point
real(8), allocatable :: dan_x(:,:), dan_z(:,:)  ! xy record coords
real(8), allocatable :: dan_rn(:,:)  ! roughness of coords
! table
integer :: max_tbl
integer, allocatable :: n_tbl(:)  ! table num
real(8), allocatable :: h_table(:,:), a_table(:,:), r_table(:,:)  ! h, a, r
real(8), allocatable :: b_table(:,:), rn_table(:,:)  ! b, rn
! boundary
integer nbound_1d  ! boundary num
integer, allocatable :: b_idx_1d(:)  ! bound id
real(8), allocatable :: b_dt_1d(:), b_data_1d(:,:)
real(8), allocatable :: up_h_1d(:)  ! upstream(dummy cross section h)
real(8), allocatable :: up_q_1d(:)  ! upstream(dummy cross section q)
integer, allocatable :: b_upme_1d(:)  ! mesh
integer, allocatable :: b_dome_1d(:)  ! mesh
! else inflow
integer :: nsubflow_1d
! else inflow array
integer, allocatable :: subflow_input(:,:)  ! out(1,:) -> in(2,:)
real(8), allocatable :: subq_all(:)  ! sum inflow
real(8), allocatable :: tributaryq_1d(:)  ! tributary
real(8), allocatable :: breakq_1d(:), breakq_l_1d(:), breakq_r_1d(:)  ! break
real(8), allocatable :: weirq_1d(:), weirq_l_1d(:), weirq_r_1d(:)  ! weir(overflow)
real(8), allocatable :: pumpq_1d(:), sluiceq_1d(:)  ! pump, sluice
real(8), allocatable :: weir_dx_1d(:)  ! weir width
! 1d <-> 2d
! weir(overflow)
integer, allocatable :: ncnct_1d2d(:,:)  ! connect mesh num(1:left, 2:right)
integer, allocatable :: cnct_1d2d_idx(:,:,:)  ! connect mesh ids(1:left, 2:right)
real(8), allocatable :: crown_1d(:,:)  ! crown elevation(1:left 2:right)
real(8), allocatable :: inland_1d(:,:)  ! inland elevation(average) (1:left, 2:right)
real(8), allocatable :: w_alpha_1d(:,:)  ! coef
real(8), allocatable :: w_angle_1d(:)  ! coef
! break(development)
integer, allocatable :: bktype_1d(:,:)  ! breakpoint flag

character(len=100) finit, fd1out, fd1mx
integer fd1out_unit, fd1mx_unit

end module d1riv_globals_mod
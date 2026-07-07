! Define Global Variables
! Coded by K.Kawaike and TK Labo
! Released on July 7th 2025

module unst_globals_mod

! ================
!  constant param
! ================
real(8), parameter :: pi = 3.14d0   ! pi
real(8), parameter :: gg = 9.8d0    ! gravitational acceleration
real(8), parameter :: th = 1.0d-3   ! Limit depth of movement
real(8), parameter :: fita = 0.5d0  !

! =============
!  time & step 
! =============
integer mstep
integer lkout, lpout
real(8) unstdt, unsttime, timmax
real(8) dt2, dtq, dtrain

! ======
!  grid
! ======
! -- num of data --
integer mesh, link, node
! -- mesh param & condition --
integer, allocatable :: ko(:), melink(:,:), menode(:,:)  ! components(link & node)
real(8), allocatable :: node_dx(:,:), node_dy(:,:)       ! components(link)
real(8), allocatable :: smesh(:)  ! area
real(8), allocatable :: xmesh(:), ymesh(:)  ! centroid coords
real(8), allocatable :: rtuv_x(:,:), rtuv_y(:,:)  ! interpolarate conf
integer, allocatable :: inf(:)  ! landuse type id
real(8), allocatable :: baseo(:), mn(:), lambda(:)
real(8), allocatable :: rnof(:)  ! rainfall correction factor(No select value v.1.0.5)
    ! - - -
real(8), allocatable :: unsth(:), ho(:), hmax(:)   ! depth
real(8), allocatable :: umm(:), uum(:), uummax(:)  ! flux & speed
real(8), allocatable :: vnm(:), vvm(:), vvmmax(:)  ! flux & speed
real(8), allocatable :: qr_sum(:)
! -- link param & condition --
integer, allocatable :: limesh(:,:), linode(:,:)  ! components(mesh & link)
real(8), allocatable :: scv(:)  !
real(8), allocatable :: rthl(:,:), ux(:), uy(:)  ! interpolarate conf
real(8), allocatable :: dl(:), rnx(:), rbeta(:), umbeta(:), vnbeta(:)
real(8), allocatable :: blink(:)  ! length  v.1.0.5
    ! - - -
real(8), allocatable :: hl(:)  ! depth
real(8), allocatable :: um(:), umo(:), uu(:)  ! flux & speed
real(8), allocatable :: vn(:), vno(:), vv(:)  ! flux & speed
integer, allocatable :: lhan(:), lhano(:)   ! interpolarate switch
! -- node param --
real(8), allocatable :: dnox(:), dnoy(:)

! =====
!  qin
! =====
! -- switch --
integer inls
! -- num of data --
integer iqnum, iqin
! -- qin param -- 
integer, allocatable :: inl(:)  ! link id or mesh id
integer, allocatable :: lkyokai_dir(:)
real(8), allocatable :: lkyokai_dx(:), lkyokai_dy(:)  ! link length
real(8), allocatable :: qin(:,:)

! ======
!  rain
! ======
real(8), allocatable :: rain(:,:,:)  ! rain data
integer, allocatable :: urain_i(:), urain_j(:)  ! match mesh_id-rain_ij

! =============================
!  unst用鉛直浸透モデル(infilt)
! =============================
real(8), allocatable, save :: unst_soildepth(:)
real(8), allocatable, save :: unst_gammaa(:)
real(8), allocatable, save :: unst_ksv(:)
real(8), allocatable, save :: unst_faif(:)
real(8), allocatable, save :: unst_infilt_limit(:)
real(8), allocatable, save :: unst_gampt_ff(:), unst_gampt_f(:)

! ======
!  else
! ======
real(8) unst_error_v, unst_dis_v ! numerical error & discharge q

! ======================================
!  flow resistance caused by vegetation
! ======================================
! -- switch --
integer plantFN, plantDa
! -- param & condition --
integer, allocatable :: plantF_array(:), plantN_array(:)       !植生帯（ヨシ帯）
real(8), allocatable :: plant_D_array(:), plant_a_array(:)     !植生帯（樹林帯）
real(8), allocatable :: plant_hv_array(:), vr_cd(:)            !植生帯（樹林帯）
real(8), allocatable :: dk_val(:), plant_lambda(:)             !植生帯（樹林帯）
real(8), allocatable :: vr_hv(:)                               !植生帯（樹林帯）

real(8), allocatable :: rri_x(:), rri_y(:)

! ==========
!  Paddydam
! ==========
! -- switch --
integer paddydam
! -- param & condition --
integer, allocatable :: paddyid(:), device(:), phid(:)
integer, allocatable :: pqout_idx(:), pdrain(:), min_pmeshid(:)
integer wtyp, paddy, ca, nhp
real(8), allocatable :: pqh(:), min_dist(:), totalqp(:)
real(8), allocatable :: dr_dist(:), dhp(:,:), dhj(:,:)
real(8), allocatable :: psmesh(:), paddy_q(:), etp(:)            !paddydat
real(8) lh, wh1, wh2, ww1, ww2, dld, unstdh, p_data, ph, dd, outa    !paddydat
integer, allocatable :: ttp(:), orifice_num(:), paddycount_num(:)    !paddydat
integer, allocatable :: paddyid2mesh(:,:), drain2phidx(:)
integer max_paddycount_num

! ==========
!  Drainage
! ==========
! -- switch --
integer drainarea
! -- param & condition --
integer, allocatable :: inf_dr(:), tripTime(:)                   !下水道圃場
integer dr_no
integer, allocatable :: drp(:), drc(:)                       !下水道・圃場モデル
real(8), allocatable :: drr(:), drr_dist(:)                     !下水道・圃場モデル
real(8), allocatable :: vol_dr(:), vol(:)                       !下水道・圃場モデル

! =================
!  Line Embankment
! =================
integer mmorid, morid
integer, allocatable :: infl(:)
real(8), allocatable :: zbbk(:)

! ===========
!  Discharge
! ===========
! -- switch --
integer dsmesh
! -- num of data --
integer dsmenum, ids2mesh, idsdepth
! -- param & condition --
real(8), allocatable :: ds_dt(:)
integer, allocatable :: ds_inf(:), ds_upme(:), ds2me(:)
real(8), allocatable :: sum_rtuv(:), sum_rtuv_x(:), sum_rtuv_y(:)
real(8), allocatable :: ds_wl(:,:)

! ==================
!  file path & unit
! ==================
! -- file path --
character(len=50) :: &
    fdhp, fpaddyh, fpaddyq, &
    fh, fhmx, fuum, fvvm, &
    fuumx, fvvmx, fstorage, fq, &
    fvminus, fkyokaiq
! -- file path (use from RRI) --
character(len=50) :: &
    fdsmesh, fplantF, fplantN, fplantD, fplanta, &
    fpaddy, fpqout, fpaddy_param, &
    finf_dr, fdrain, fmorid, fd1riv_cntl, &
    fplanthv, fplantcd
! -- unit --
integer :: &
    fdhp_unit, fpaddyh_unit, fpaddyq_unit, &
    fh_unit, fhmx_unit, fuum_unit, fvvm_unit, &
    fuumx_unit, fvvmx_unit, fstorage_unit, fq_unit

! ======================
!  1d-riv (development)
! ======================
! -- switch --
integer d1riv
! -- num of data --

! -- param & condition --
real(8), allocatable :: riv_eq(:)

end module unst_globals_mod
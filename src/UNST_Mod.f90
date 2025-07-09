module globals

!------------------UNST MODEL----------------
real(8) pi
real(8) gg, unstdt, dt2, unsttime, dtq, th, fita, dtrain, precond_time
real(8), allocatable :: unsth(:), ho(:), hl(:), hmax(:), uummax(:), vvmmax(:)
real(8), allocatable :: um(:), umo(:), umm(:), uu(:), vn(:), vno(:), vnm(:), vv(:)
real(8), allocatable :: baseo(:)
real(8), allocatable :: dnox(:), dnoy(:)
real(8), allocatable :: scv(:), rthl(:,:), ux(:), uy(:), smesh(:)
real(8), allocatable :: mn(:), rnof(:), lambda(:), rbeta(:), umbeta(:), vnbeta(:)
real(8), allocatable :: xmesh(:), ymesh(:), rtuv(:, :)
real(8), allocatable :: rtuv_x(:,:), rtuv_y(:,:)
real(8), allocatable :: uum(:), vvm(:)
integer sep_rtuv
real(8) ocpy, unstbeta
real(8), allocatable :: qin(:,:), qin1(:)   !RRI_UNST-2Dはqin1(:) >> qinu(:,:), qinv(:,:)
integer, allocatable :: limesh(:,:), linode(:,:)
integer, allocatable :: inf(:)    
integer, allocatable :: lhan(:), lhano(:)
integer, allocatable :: ko(:), menode(:,:), melink(:,:)
integer mesh, link, node, mstep
integer iqnum, iqin
integer, allocatable :: inl(:)
real(8), allocatable :: rnx(:), dl(:)
integer plantFN, plantDa, paddydam, drainarea, dsmesh, precond_hot

real(8) timmax
real(8), allocatable :: node_dx(:,:), node_dy(:,:), lkyokai_dx(:), lkyokai_dy(:), lkyokai_vect(:)
integer lkout, lpout
integer inls, nx_rain, ny_rain, t_rain !UNST-2D only
real(8), allocatable :: rri_x(:), rri_y(:)

! 樹林帯モデル
integer, allocatable :: plantF_array(:), plantN_array(:)      !樹林帯（ヨシなど）
real(8), allocatable :: plant_D_array(:), plant_a_array(:)                         !樹林帯(防備林)
real(8), allocatable :: dk_val(:), plant_lambda(:)             !樹林帯(防備林)

! 田んぼダムモデル
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

! 下水道・圃場モデル
integer, allocatable :: inf_dr(:), tripTime(:)                   !下水道圃場             
integer dr_no
integer, allocatable :: drp(:), drc(:)                       !下水道・圃場モデル
real(8), allocatable :: drr(:), drr_dist(:)                     !下水道・圃場モデル
real(8), allocatable :: vol_dr(:), vol(:)                       !下水道・圃場モデル


real(8), allocatable :: rain(:,:,:)    
integer, allocatable :: urain_i(:), urain_j(:)
integer tt_max_rain !UNST-2D only

!盛り土
integer mmorid, morid
integer, allocatable :: nmorili(:), infl(:), lnode1(:), lnode2(:), neib1(:), neib2(:)
real(8), allocatable :: zbbk(:)

!水収支
real(8), allocatable :: qr_sum(:)

! 1driv k.kawaike - UNST-2D original
real(8) cslmd
integer d1riv, msctn
integer knode, knode1, knode2, ksetu
real(8), allocatable :: hn(:), qn(:), hs(:), qs(:), us(:), as(:), bs(:)
real(8), allocatable :: htes(:), fs(:), cap(:), cbm(:), qys(:), sap(:), sbm(:)
integer, allocatable :: iptn(:), ihnum(:)
real(8), allocatable :: hb(:), rn(:), dx2(:)
real(8), allocatable :: ha(:,:), hr(:,:)
integer, allocatable :: mnd(:), nd(:,:)
real(8), allocatable :: qhyd(:,:), hhyd(:,:)
integer, allocatable :: kkousi(:), kkasen(:)
integer, parameter :: str_type = 0  ! 2507 廃止予定
integer count1, count2, count3, count4
real(8), allocatable :: unstc(:), co(:), cr(:), str(:), stro(:), strmx(:)
real(8) qout, dvr, vrain, sv0, v_minus_all, v_minus(100), phi(100)
real(8) svc, v_cminus, v_cplus, v_cextre !tracer
! weir(1d-2d) d.baba - UNST-2D original
real(8), allocatable :: f_qs(:)
real(8), allocatable :: dq_l(:), dq_r(:)
real(8), allocatable :: lcrown(:), rcrown(:)
real(8), allocatable :: wir_alpha(:), wir_angle(:)
! tabele d.baba - UNST-2D original
integer ndan_num
real(8) max_dz
integer, allocatable :: ndan(:)
real(8), allocatable :: dcoord_x(:,:), dcoord_z(:,:)
real(8), allocatable :: min_z(:), max_z(:)
real(8), allocatable :: ah(:,:), hss(:,:)

! 下流端からの排水
integer dsmenum, idsfilter2, idsdepth
real(8), allocatable :: dsdt(:)
integer, allocatable :: dsinf(:), dsupper(:), dsfilter2(:)
real(8), allocatable :: sum_rtuv(:), sum_rtuv_x(:), sum_rtuv_y(:)
real(8), allocatable :: dsdepth(:,:)

! 遺伝的アルゴリズム  2507 廃止予定
integer, allocatable :: genes(:), paddycluster(:)
integer ga, paddyclass

! unst用鉛直浸透モデル(infilt)
real(8), allocatable, save :: unst_soildepth(:)
real(8), allocatable, save :: unst_gammaa(:)
real(8), allocatable, save :: unst_ksv(:)
real(8), allocatable, save :: unst_faif(:)
real(8), allocatable, save :: unst_infilt_limit(:)
real(8), allocatable, save :: unst_gampt_ff(:), unst_gampt_f(:)

end module globals

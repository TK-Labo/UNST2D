!***********************************************************************
! UNST_Main.f90
! Coded by K.Kawaike, Added by TK Labo
! Released on July 7th 2025
!***********************************************************************

!-----------------------------------------------------------------------
! UNSTメインプログラム/UNST2D main program
!-----------------------------------------------------------------------
program UNST2D
    use unst_globals_mod
    use unst_read
    use unst_cal_sub
    use d1riv_globals_mod
    use unst_d1riv_read
    use unst_1d_main
    use unst_wrfile
    implicit none
    integer i, me, li, k

    write(*,*) 'UNST version1_0_5_1'  ! version 1.0.5.1 (release 26/07)

    ! データ読み込み など/data read etc
    call unst_rdat  ! メイン読み込み/main read
    call open_unst_output_files  ! 出力ファイルを開く/open output files
    ! オプション読み込み
    if(dsmesh==1) call dsmeshdat  ! 境界流出設定ファイル(dsmesh.dat)の読み込み/read the boundary outflow condition file (dsmesh.dat)
    if(d1riv==1) call d1rivdat(timmax, dt2, mesh, baseo, fd1riv_cntl)  ! 1次元河道の読み込み/read 1D-river
    if(plantFN==1) call plantFNdat  ! 植生ファイル(倒伏あり)の読み込み/read the vegetation file (considering deformation due to drag force)
    if(plantDa==1) call plantDadat  ! 植生ファイル(倒伏なし;樹木竹林)の読み込み/read the vegetation file (tree or bamboo)
    if(paddydam==1) call paddydat  ! 田んぼファイル(田んぼダム含む)の読み込み/read the paddy or paddydam file 
    if(drainarea==1) call draindat  ! 下水道圃場ファイルの読み込み/read the drainage area and the parameters file
    if(morid==1) call moriddat  ! 線盛土ファイルの読み込み/read the line embankment file
    
    ! 各種変数の初期化initiald condition
    call unst_initiald
    if(d1riv==1) call d1rivinitiald(dt2)  ! 1次元河道用/for 1D-river
    if(paddydam==1) call paddyinitiald  ! 田んぼ用/for paddy
    if(drainarea==1) call draininitiald  ! 下水道圃場用/for drainage

    unsttime = 0.0d0  ! calculation time
    mstep = 0  ! for output counter

    ! zero-time output
    call diskwrite  ! disk output
    call dispwrite  ! display output
    if(paddydam==1) call paddywrite  ! disk output  
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                  Loop Start
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    do while (unsttime + unstdt <= timmax)
        ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !  Equation of motion (cal flux)
        ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        call flux   ! cal link flux
        call lkyokai ! apply boundary condition

        ! flow flux(velocity) ahead of the flood inundation front set
        !$omp parallel do default(shared),private(me,li,k)
        do me = 1, mesh
            if(unsth(me) >= th) cycle  ! early return(dry mesh)
            do k = 1, ko(me)
                li = melink(k, me)
                if((um(li)*node_dy(k, me) - vn(li)*node_dx(k, me)) > 0.0d0) then
                    um(li) = 0.0d0
                    vn(li) = 0.0d0
                endif
            enddo
        enddo
        !$omp end parallel do

        if(d1riv==1) then
            call calc_2d_to_1d_inflow
            call calc_1d_to_2d_outflow
            call d1riv_main(unsttime)
        endif

        ! update time(1)
        unsttime = unsttime + unstdt
        mstep = mstep + 1

        ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !  Equation of continuity (cal waterdepth)
        ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        call suisin  ! cal mesh water depth

        ! cal velocity
        call velocity

        ! update max record
        !$omp parallel do default(shared),private(me)
        do me = 1, mesh
            hmax(me) = max(hmax(me), unsth(me))
            uummax(me) = max(uummax(me), abs(uum(me)))
            vvmmax(me) = max(vvmmax(me), abs(vvm(me)))
        enddo
        !$omp end parallel do

        ! update variable（new >>> old）
        call replace

        ! reset subflow
        if(d1riv==1) then
            riv_eq = 0.0d0
            !$omp parallel do default(shared),private(i)
            do i = 1, ndan
                h_1dmax(i) = max(h_1dmax(i), h_1d(i))
                vv_1dmax(i) = max((vv_1dmax(i)), abs(vv_1d(i)))
            enddo
            !$omp end parallel do
        endif

        ! update time(2)
        unsttime = unsttime + unstdt
        mstep = mstep + 1

        ! output result
        if(mod(mstep, lkout) == 0) call diskwrite   ! disk
        if(mod(mstep, lpout) == 0) call dispwrite   ! display
        if(paddydam==1 .and. mod(mstep, lkout) == 0) call paddywrite  ! paddydam

    enddo

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                  Loop End
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    call wrhmax
    write(*, 1999) unsttime, timmax

1999 format('      - normal end -  time=', f8.0, '  timmax=', f8.0)

!
    ! deallocate variables
    call close_unst_output_files
    
    deallocate(baseo, dnox, dnoy, smesh, scv, rthl, ux, uy, xmesh, ymesh, rtuv_x, rtuv_y)
    deallocate(limesh, linode, inf, ko, menode, melink, inl, qin, lkyokai_dx, lkyokai_dy)
    deallocate(unsth, ho, hl, hmax, uummax, vvmmax)
    deallocate(um, umo, umm, uu, vn, vno, vnm, vv)
    deallocate(mn, rnof, lambda, rbeta)
    deallocate(uum, vvm, lhan, lhano, qr_sum, rnx, dl)
    if(plantFN==1) deallocate(plantF_array, plantN_array)
    if(plantDa==1) deallocate(dk_val)
    if(paddydam==1) deallocate(paddyid, pqout_idx, pdrain, min_pmeshid, device)
    if(paddydam==1) deallocate(orifice_num, min_dist, psmesh, dr_dist, dhp, phid)
    if(paddydam==1) deallocate(paddy_q, pqh, drain2phidx)
    if(drainarea==1) deallocate(inf_dr, drp, drr, drr_dist, dhj)
    if(dsmesh==1) then
        if(any(ds_inf==2)) deallocate(ds_wl)
        if(dsmesh==1) deallocate(ds_dt, ds_inf, ds2me, ds_upme)
    endif
    if(d1riv==1) deallocate(riv_ndan, h_1d, ho_1d, vv_1d, q_1d, qo_1d, a_1d, r_1d, b_1d, rn_1d)
    if(d1riv==1) deallocate(h_1dmax, ntype_1d, kp_1d, dx_1d, depth_1d, rbed_1d, rzmax_1d)
    if(d1riv==1) deallocate(dan_record, n_tbl, h_table, a_table, r_table, b_table, rn_table)
    if(d1riv==1) deallocate(b_idx_1d, b_dt_1d, b_data_1d, up_h_1d, up_q_1d, b_upme_1d, b_dome_1d)
    if(d1riv==1) deallocate(subflow_input, subq_all, tributaryq_1d, pumpq_1d, sluiceq_1d)
    if(d1riv==1) deallocate(breakq_1d, breakq_l_1d, breakq_r_1d)
    if(d1riv==1) deallocate(weirq_1d, weirq_l_1d, weirq_r_1d, weir_dx_1d)
    if(d1riv==1) deallocate(w_alpha_1d, w_angle_1d, bktype_1d)
    if(d1riv==1) deallocate(ncnct_1d2d, cnct_1d2d_idx, crown_1d, inland_1d)
    
end program UNST2D

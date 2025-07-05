!***********************************************************************
! UNST_Main.f90
! Coded by K.Kawaike and TK Labo
! Released on July 7th 2025 
!***********************************************************************

program UNST
    use globals
    implicit none
    integer i, me, li, k, iqmax, ihmax
    real(8) alpg2, alpbet, cslgb, alpha

    write(*,*) 'UNST version0_0_0'

    ! data read
    call unst_rdat
    if(plantFN==1) call plantFNdat
    if(plantDa==1) call plantDadat
    if(paddydam==1) call paddydat
    if(drainarea==1) call draindat
    if(morid==1) call moriddat
    if(d1riv==1) then
        alpha = 1.0d0
        unstbeta = 1.0d0
        cslmd = 1.0d0
        alpg2 = alpha/(2.0d0*gg)
        alpbet = (alpha-unstbeta)/(2.0d0*unstbeta)
        cslgb = cslmd*gg/unstbeta
        call d1rivdat
    endif

    if(d1riv==1) call bndry(iqmax,ihmax)
    
    ! initiald condition
    call unst_initiald
    if(paddydam==1) call paddyinitiald
    if(drainarea==1) call draininitiald

    unsttime = 0.0d0
    mstep = 0

    ! 田んぼダムの適用
    if(str_type == 0) then
        phi(61) = 0.20d0
        phi(63) = 0.20d0
    elseif(str_type == 1) then
        phi(61) = 0.05d0
        phi(63) = 0.20d0
    endif

    if(d1riv==1) call elm(alpg2, alpbet, cslgb)

    call diskwrite
    call dispwrite
    if(paddydam==1) call paddywrite
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                  loop start↓
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! ----------------------
    !  equation of motion
    ! ----------------------
    if(d1riv==1) goto 100
  1 call flux
    if(inls /= 2) call lkyokai

    ! hangling tips (sentan no toriatsukai)
    !$omp parallel do default(shared),private(me,li,k)
    do me = 1, mesh
        if(unsth(me) >= th) goto 200
        do k = 1, ko(me)
            li = melink(me, k)
            if((um(li)*node_dy(me,k) - vn(li)*node_dx(me,k)) > 0.0d0) then
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        enddo
200 enddo
    !$omp end parallel do

    !河川網の計算
    if(d1riv==1) then
        call links
        call inodes
        call bnodes(iqmax, ihmax)
        call elm(alpg2, alpbet, cslgb)
    endif
100 continue

    !横流入流量リセット
    if(d1riv==1) then
    !$omp parallel do default(shared),private(i)
    do i = 1, msctn
        qys(i) = 0.0d0
    enddo
    !$omp end parallel do
    endif

    ! time step       
    unsttime = unsttime + unstdt
    mstep = mstep + 1

    ! -----------------------------
    !  equation of continuity 
    ! -----------------------------
    call suisin

    ! calculate velocity
    call velocity

    ! calculate max h and v 
    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        hmax(me) = max(hmax(me), unsth(me))
        uummax(me) = max(uummax(me), abs(uum(me)))
        vvmmax(me) = max(vvmmax(me), abs(vvm(me)))
    enddo
    !$omp end parallel do
    
    ! replace data ( new >>> old )
    call replace

    ! time step 
    unsttime = unsttime + unstdt
    mstep = mstep + 1

    ! output
    if(mod(mstep, lkout) == 0) call diskwrite
    if(mod(mstep, lpout) == 0) call dispwrite
    if(paddydam==1 .and. mod(mstep, lkout) == 0) call paddywrite

    ! time judging
    if(unsttime + unstdt <= timmax) goto 1
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                  loop end
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    call wrhmax
    write(*, 1999) unsttime, timmax
1999 format('      - normal end -  time=', f8.0, '  timmax=', f8.0)

    close(91)
    close(93)
    close(94)
    close(95)
    close(96)
    close(97)
!
    ! deallocate variables
    deallocate(baseo, dnox, dnoy, smesh, scv, rthl, ux, uy, xmesh, ymesh, rtuv_x, rtuv_y)
    deallocate(limesh, linode, inf, ko, menode, melink, inl)
    if(inls == 1) deallocate(qin,lkyokai_dx,lkyokai_dy)
    if(inls == 0) deallocate(qin1,lkyokai_dx,lkyokai_dy)
    deallocate(unsth, ho, hl, hmax, uummax, vvmmax)
    deallocate(um, umo, umm, uu, vn, vno, vnm, vv)
    deallocate(mn, rnof, lambda, rbeta, umbeta, vnbeta)
    deallocate(uum, vvm, lhan, lhano, qr_sum,rnx,dl)
    deallocate(unstc, co, str, stro, strmx)
    if(plantFN==1) deallocate(plantF_array, plantN_array)
    if(plantDa==1) deallocate(plant_D_array,plant_a_array,dk_val)
    if(paddydam==1) deallocate(paddyid, pqout_idx, pdrain, min_pmeshid, device)
    if(paddydam==1) deallocate(orifice_num, min_dist, psmesh, dr_dist, dhp, phid)
    if(paddydam==1) deallocate(paddy_q, pqh)
    if(drainarea==1) deallocate(inf_dr, drp, drr, drr_dist, dhj)
    if(d1riv==1) deallocate(iptn, hb, rn, dx2, ihnum, ha, hr, mnd, nd)
    if(d1riv==1) deallocate(kkousi, kkasen, qhyd, hhyd)
    if(d1riv==1) deallocate(hn, qn, hs, us, as, bs, htes, fs, cap, cbm, qys, sap, sbm)
    if(dsmesh==1) deallocate(dsdt, dsinf, dsupper, dsfilter2)
    if(ga==1) deallocate(genes)
    
end program UNST

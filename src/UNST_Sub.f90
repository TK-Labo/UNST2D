! Calculation Subroutines
! Coded by K.Kawaike and TK Labo
! Released on July 7th 2025

!================================
! Equation of motion, continuity
!================================
!--------------------------------------------------------------------------
! UNST2D 主な変数/UNST2D Main variable
!   um(link): linkのx軸方向のフラックス(t=n+2)/link flux  x-axis(t=n+2)
!   umo(link): linkのx軸方向のフラックス(t=n)/link flux  x-axis(t=n)
!   uu(link): linkのx軸方向の流速/link velocity  x-axis
!   vn(link): linkのy軸方向のフラックス(t=n+2)/link flux  y-axis(t=n+2)
!   vno(link): linkのy軸方向のフラックス(t=n)/link flux  y-axis(t=n)
!   vv(link): linkのy軸方向の流速/link velocity  x-axis
!   hl(link): linkの(補間)水深/link depth
!   umm(mesh): meshのx軸方向のフラックス/mesh flux  x-axis
!   uum(mesh): meshのx軸方向の流速/mesh velocity  x-axis
!   vnm(mesh): meshのy軸方向のフラックス/mesh flux  y-axis
!   vvm(mesh): meshのy軸方向の流速/mesh velocity  y-axis
!   unsth(mesh): meshの水深(t=n+1)/mesh depth(t=n+1)
!   ho(mesh): meshの水深(t=n-1)/mesh depth(t=n-1)
!--------------------------------------------------------------------------
module unst_cal_sub
    use unst_globals_mod
    use, intrinsic :: ieee_arithmetic
    implicit none

contains
!-------------------------------------
! Calculate flux (equation of motion)
!-------------------------------------
!---------------------------------------------------------------------------------
! UNST2D 運動方程式/UNST2D equation of motion
! 変数/Variable:
! > 一時変数/temporary variable
!   hhe: linkに隣接するmesh(2)の水位/water level of the mesh(2) adjacent to the link
!   hhw: linkに隣接するmesh(1)の水位/water level of the mesh(1) adjacent to the link
!   hhep: linkに隣接するmesh(1)の有効水深(水深-移動限界水深th)/
!         effective depth of mesh (2) adjacent to the link 
!          (depth - limit depth of movement th)
!   hhwp: linkに隣接するmesh(1)の有効水深(水深-移動限界水深th)/
!         effective depth of mesh (1) adjacent to the link 
!          (depth - limit depth of movement th)
!   hhan: linkに隣接する2つのmeshの水位(水深)差/
!         the difference in water level (depth) between two meshes 
!          adjacent to the link
!   sgn: positive and negative
!   hh1: linkに隣接するmeshの高い方の水位と高い方の標高値の差/
!        the difference between the highest water level and 
!         the highest elevation value of the meshes adjacent to the link
!   u11: 対流項/convective term
!   v11: 対流項/convective term
!   u13: 圧力項/gravitational term
!   v13: 圧力項/gravitational term
!   sqx: sqrt(uu**2 + vv**2)
!   ram: 摩擦項/friction term
!---------------------------------------------------------------------------------
subroutine flux
    implicit none
    real(8) hhe, hhw, hhep, hhwp, hhan, sgn, hh1 !Discontinuous(Drop formula, Overflow formula)
    real(8) u11, v11, u13, v13, sqx, ram !Continuous
    integer li, k, me1, me2
    real(8) b, h1, h2, uvmn, vol1, vol2
    real(8) vr_hl, vr_u, vr_v

    !$omp parallel do default(shared) &
    !$omp private(hhe, hhw, hhep, hhwp, hhan, sgn, hh1, u11, v11, u13, v13, sqx, ram, li, k, me1, me2) &
    !$omp private(b, h1, h2, uvmn, vol1, vol2, vr_hl, vr_u, vr_v)
  do li = 1, link
    if(limesh(li, 2) == 0) cycle
    if((inf(limesh(li, 2)) == 0) .or. (inf(limesh(li, 1)) == 0) .or. &
       (rbeta(li) == 0.0d0) .or. &
       (unsth(limesh(li, 1)) <= th .and. unsth(limesh(li, 2)) <= th)) then
        um(li) = 0.0d0
        vn(li) = 0.0d0
        lhan(li) = 1
        cycle
    endif
    ! water level, effective water depth set
    hhe = unsth(limesh(li, 2)) + baseo(limesh(li, 2))
    hhw = unsth(limesh(li, 1)) + baseo(limesh(li, 1))
    hhep = unsth(limesh(li, 2)) - th
    hhwp = unsth(limesh(li, 1)) - th

    ! line embankment added by CTI r.nishizawa
    if (morid == 1) then
    if(infl(li) == 1) then
        b = blink(li)
        h1 = max(hhe, hhw) - zbbk(li)
        h2 = min(hhe, hhw) - zbbk(li)

        hhan = hhe - hhw
        sgn = sign(1.0d0, hhan)

        if (h1 > th) then
            if((h2/h1) <= 2.0d0/3.0d0 .and. h1 > th) then
                um(li) = -sgn * 0.35d0 * h1 * sqrt(2.0d0*gg*h1)*ux(li)
                vn(li) = -sgn * 0.35d0 * h1 * sqrt(2.0d0*gg*h1)*uy(li)
            elseif((h2/h1) > 2.0d0/3.0d0 .and. h1 > th) then
                um(li) = -sgn * 0.91d0 * h2 * sqrt(2.0d0*gg*(h1-h2))*ux(li)
                vn(li) = -sgn * 0.91d0 * h2 * sqrt(2.0d0*gg*(h1-h2))*uy(li)
            endif

            uvmn = sqrt(um(li)**2 + vn(li)**2) * b * dt2
            vol1 = unsth(limesh(li,1)) * smesh(limesh(li,1)) * (1.0d0 - lambda(limesh(li,1)))
            vol2 = unsth(limesh(li,2)) * smesh(limesh(li,2)) * (1.0d0 - lambda(limesh(li,2)))

            if (hhw > hhe .and. uvmn > vol1) then
                um(li) = um(li) * vol1/uvmn
                vn(li) = vn(li) * vol1/uvmn
            elseif (hhe > hhw .and. uvmn > vol2) then
                um(li) = um(li) * vol2/uvmn
                vn(li) = vn(li) * vol2/uvmn
            endif
        else
            um(li) = 0.0d0
            vn(li) = 0.0d0
        endif
        lhan(li) = 1
        cycle
    endif
    endif

    ! when the water surface is not continuous
    ! 段落ちの式
    if(hhe < baseo(limesh(li, 1))) then
        if ((inf(limesh(li,1))==71 .and. inf(limesh(li,2))<71) .or. &
            (inf(limesh(li,2))==71 .and. inf(limesh(li,1))<71)) then
            !田んぼダム適用 by k.yamamura(inf=71)
            if(unsth(limesh(li,1)) > lh + th) then
                um(li) = 0.544d0*(unsth(limesh(li, 1)) - lh)*sqrt(gg*(unsth(limesh(li, 1)) - lh))*ux(li)
                vn(li) = 0.544d0*(unsth(limesh(li, 1)) - lh)*sqrt(gg*(unsth(limesh(li, 1)) - lh))*uy(li)
            else
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        else
            ! 田んぼダム適用場所以外
            if(unsth(limesh(li, 1)) > th) then
                if(baseo(limesh(li,2))==-9999.0d0 .and. baseo(limesh(li,1))/=-9999.0d0) then
                    !自由流出
                    um(li) = unsth(limesh(li,1))*sqrt(2.0d0*gg*unsth(limesh(li,1)))*ux(li)
                    vn(li) = unsth(limesh(li,1))*sqrt(2.0d0*gg*unsth(limesh(li,1)))*uy(li)
                else
                    um(li) = 0.544d0*unsth(limesh(li, 1))*sqrt(gg*unsth(limesh(li, 1)))*ux(li)
                    vn(li) = 0.544d0*unsth(limesh(li, 1))*sqrt(gg*unsth(limesh(li, 1)))*uy(li)
                endif
            else
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        endif
        lhan(li) = 1
        cycle
    elseif (hhw < baseo(limesh(li, 2))) then
        if ((inf(limesh(li,1))==71 .and. inf(limesh(li,2))<71) .or. &
            (inf(limesh(li,2))==71 .and. inf(limesh(li,1))<71)) then
            ! 田んぼダム適用 by k.yamamura(inf=71)
            if(unsth(limesh(li,2)) > lh + th) then
                um(li) = -0.544d0*(unsth(limesh(li, 2)) - lh)*sqrt(gg*(unsth(limesh(li, 2)) - lh))*ux(li)
                vn(li) = -0.544d0*(unsth(limesh(li, 2)) - lh)*sqrt(gg*(unsth(limesh(li, 2)) - lh))*uy(li)
            else
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        else
            !　田んぼダム適用以外
            if(unsth(limesh(li,2)) > th) then
                if(baseo(limesh(li,1))==-9999.0d0 .and. baseo(limesh(li,2))/=-9999.0d0) then
                    !自由流出
                    um(li) = -unsth(limesh(li,2))*sqrt(2.0d0*gg*unsth(limesh(li,2)))*ux(li)
                    vn(li) = -unsth(limesh(li,2))*sqrt(2.0d0*gg*unsth(limesh(li,2)))*uy(li)
                else
                    um(li) = -0.544d0*unsth(limesh(li, 2))*sqrt(gg*unsth(limesh(li, 2)))*ux(li)
                    vn(li) = -0.544d0*unsth(limesh(li, 2))*sqrt(gg*unsth(limesh(li, 2)))*uy(li)
                endif
            else
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        endif
        lhan(li) = 1
        cycle

    ! 完全越流の式
    elseif(hhep*hhwp < 0.0d0) then
        if(unsth(limesh(li, 2)) > 0.0d0 .or. unsth(limesh(li, 1)) > 0.0d0) then
            hhan = hhep - hhwp
            sgn = sign(1.0d0, hhan)
            hh1 =   max(unsth(limesh(li, 2))+baseo(limesh(li, 2)), unsth(limesh(li, 1)) + baseo(limesh(li, 1))) &
                    - max(baseo(limesh(li, 2)), baseo(limesh(li, 1)))
            if ((inf(limesh(li,1))==71 .and. inf(limesh(li,2))<71) .or. &
                (inf(limesh(li,2))==71 .and. inf(limesh(li,1))<71)) then
                ! 田んぼダム適用
                if(hh1 > lh + th) then
                    um(li) = - sgn*0.35d0*(hh1-lh)*sqrt(2.0d0*gg*(hh1-lh))*ux(li)
                    vn(li) = - sgn*0.35d0*(hh1-lh)*sqrt(2.0d0*gg*(hh1-lh))*uy(li)
                else
                    um(li) = 0.0d0
                    vn(li) = 0.0d0
                endif
            else
                ! 田んぼダム適用以外
                um(li) = - sgn*0.35d0*hh1*sqrt(2.0d0*gg*hh1)*ux(li)
                vn(li) = - sgn*0.35d0*hh1*sqrt(2.0d0*gg*hh1)*uy(li)
            endif
        else
            um(li) = 0.0d0
            vn(li) = 0.0d0
        endif
        lhan(li) = 1
        cycle
    endif

    ! when the water surface is continuous
    ! convective terms
    u11 = 0.0d0
    v11 = 0.0d0
    do k = 1, ko(limesh(li, 1))
        me1 = 0
        me2 = 0
        if((lhano(melink(limesh(li, 1), k)) == 1) .or. (melink(limesh(li, 1), k) == li) .or. &
           (uu(melink(limesh(li, 1), k))==0.0d0 .and. vv(melink(limesh(li, 1), k))==0.0d0)     ) cycle
        if(uu(melink(limesh(li, 1), k))*node_dy(limesh(li,1),k) > 0.0d0) then
            me1 = limesh(li, 1)
        else
            if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 1)) me1 = limesh(melink(limesh(li, 1), k), 2)
            if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 2)) me1 = limesh(melink(limesh(li, 1), k), 1)
            ! water level
            if(me1 == 0 .and. dsmesh == 1) then
                if (ds_inf(limesh(li,1))== 2 .or. ds_inf(limesh(li,1))== 3) me1 = limesh(li, 1)
            endif
        endif

        if(vv(melink(limesh(li, 1), k))*node_dx(limesh(li,1),k) <0.0d0) then
            me2 = limesh(li, 1)
        else
            if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 1)) me2 = limesh(melink(limesh(li, 1), k), 2)
            if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 2)) me2 = limesh(melink(limesh(li, 1), k), 1)
            ! water level
            if(me2 == 0 .and. dsmesh == 1) then
                if (ds_inf(limesh(li,1))== 2 .or. ds_inf(limesh(li,1))== 3) me2 = limesh(li, 1)
            endif
        endif
        if (me1==0 .and. me2==0) cycle
        u11 = u11 + uu(melink(limesh(li, 1), k))*umm(me1)*node_dy(limesh(li,1),k) &
            - vv(melink(limesh(li, 1), k))*umm(me2)*node_dx(limesh(li,1),k)
        v11 = v11 + uu(melink(limesh(li, 1), k))*vnm(me1)*node_dy(limesh(li,1),k) &
            - vv(melink(limesh(li, 1), k))*vnm(me2)*node_dx(limesh(li,1),k)
    enddo

    do k = 1, ko(limesh(li, 2))
        me1 = 0
        me2 = 0
        if((lhano(melink(limesh(li, 2), k)) == 1) .or. (melink(limesh(li, 2), k) == li) .or. &
           (uu(melink(limesh(li, 2), k))==0.0d0 .and. vv(melink(limesh(li, 2), k))==0.0d0)     ) cycle
        if(uu(melink(limesh(li, 2), k))*node_dy(limesh(li,2),k) > 0.0d0) then
            me1 = limesh(li, 2)
        else
            if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 1)) me1 = limesh(melink(limesh(li, 2), k), 2)
            if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 2)) me1 = limesh(melink(limesh(li, 2), k), 1)
            ! water level
            if(me1 == 0 .and. dsmesh == 1) then
                if (ds_inf(limesh(li,2)) == 2 .or. ds_inf(limesh(li,2)) == 3) me1 = limesh(li, 2)
            endif
        endif

        if(vv(melink(limesh(li, 2), k))*node_dx(limesh(li,2),k) < 0.0d0) then
            me2 = limesh(li, 2)
        else
            if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 1)) me2 = limesh(melink(limesh(li, 2), k), 2)
            if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 2)) me2 = limesh(melink(limesh(li, 2), k), 1)
            ! water level
            if(me2 == 0 .and. dsmesh == 1) then
                if (ds_inf(limesh(li,2)) == 2 .or. ds_inf(limesh(li,2)) == 3) me2 = limesh(li, 2)
            endif
        endif
        if (me1==0 .and. me2==0) cycle
        u11 = u11 + uu(melink(limesh(li, 2), k))*umm(me1)*node_dy(limesh(li,2),k) &
            - vv(melink(limesh(li, 2), k))*umm(me2)*node_dx(limesh(li,2),k)
        v11 = v11 + uu(melink(limesh(li, 2), k))*vnm(me1)*node_dy(limesh(li,2),k) &
            - vv(melink(limesh(li, 2), k))*vnm(me2)*node_dx(limesh(li,2),k)
    enddo

    u11 = u11*dt2/scv(li)
    v11 = v11*dt2/scv(li)

    ! gravitational term
    if(dl(li) /= 0.0d0) then
        u13 = gg*hl(li)*dt2/dl(li)*ux(li)*(unsth(limesh(li, 2)) + baseo(limesh(li, 2)) &
            - unsth(limesh(li, 1)) - baseo(limesh(li, 1)))
        v13 = gg*hl(li)*dt2/dl(li)*uy(li)*(unsth(limesh(li, 2)) + baseo(limesh(li, 2)) &
            - unsth(limesh(li, 1)) - baseo(limesh(li, 1)))
    else
        u13 = 0.0d0
        v13 = 0.0d0
    endif

    ! shear term
    sqx = sqrt(uu(li)**2 + vv(li)**2)
    ram = gg*rnx(li)**2*sqx/hl(li)**1.333333
    if(hl(li) <= th) ram = 0.0d0

    ! defended forest term
    ! fixed  v.1.0.5
    if(plantDa==1) then
    !        plant_force = gg*(sqx**2.0d0)/(dk_val(li)**2.0d0)
    !        ram = ram + plant_force
        vr_hl = min(hl(li), vr_hv(li))
        vr_u = dk_val(li) * vr_hl * uu(li) * sqx * dt2
        vr_v = dk_val(li) * vr_hl * vv(li) * sqx * dt2
    else
        vr_u = 0.0d0
        vr_v = 0.0d0
    endif

    ! um, vn calculation
    um(li) = ((1.0d0 - dt2*ram*(1.0d0 - fita))*umo(li) - u11 - u13 + vr_u)/ &
            (1.0d0 + dt2*ram*fita)
    vn(li) = ((1.0d0 - dt2*ram*(1.0d0 - fita))*vno(li) - v11 - v13 + vr_v)/ &
            (1.0d0 + dt2*ram*fita)
    lhan(li) = 0
    enddo
    !$omp end parallel do

end subroutine flux

!------------------------
! equation of continuity
!------------------------
!-------------------------------------------------------------------------------------
! UNST2D 運動方程式/UNST2D equation of motion
! 変数/Variable:
! > 一時変数/temporary variable
!   f: linkの流出量(m^3)/link out-q
!   sumf: meshの流出量(m^3)/mesh out-q
!   tmp_wl: 水位(dsmesh,dsinf=2)/water level(dsmesh,dsinf=2)
!   unst_tmp_error: 水深計算における数値誤差/numerical errors in water depth calculations
!   unst_tmp_ds: 排水処理における排水量/wastewater volume in wastewater treatment
!   rr: meshの降水量(mm)/mesh rainfall(mm)
!   q(mesh): meshの流量/mesh q
!-------------------------------------------------------------------------------------
subroutine suisin
    implicit none
    integer me, k, ii, it, ilt, nn, ct, nt
    real(8) f, sumf, rr    !Equation of Continuity
    real(8), allocatable :: q(:)
    real(8) tmp_wl
    real(8) unst_tmp_error, unst_tmp_ds

    allocate(q(mesh))

    if (paddydam==1) paddy_q = 0.0d0
    unst_tmp_error = 0.0d0
    unst_tmp_ds = 0.0d0

    it = int(unsttime/dtrain) + 1

    if(drainarea==1) vol = 0.0d0
    if(drainarea==1) vol_dr = 0.0d0
    ct = int(unsttime / (dt2+0.000001d0))
    nt = int((unsttime + (dt2+0.000001d0)) / (dt2+0.000001d0))

    !$omp parallel do default(shared) &
    !$omp private(f,sumf,k,me,ilt,nn,tmp_wl,rr) &
    !$omp reduction(+:unst_tmp_error, unst_tmp_ds)
    do me = 1, mesh
        if(inf(me) == 0) cycle
        if(dsmesh==1) then
            ! water level(dstype:2) added by d.baba
            if (ds_inf(me)==2) then
                ii = int( unsttime / ds_dt(me)) + 1
                tmp_wl = ds_wl(me, ii) + (unsttime - ds_dt(me)*dble(ii - 1))/ds_dt(me)*(ds_wl(me, ii + 1) - ds_wl(me, ii))
                unsth(me) = max(tmp_wl - baseo(me), 0.000d0)
                cycle
            ! uniform flow(dstype:3) added by d.baba
            elseif (ds_inf(me)==3) then
                do k = 1, ko(me)
                    if (limesh(melink(me, k),1) /= ds_upme(me) .and. limesh(melink(me, k),2) /= ds_upme(me)) then
                        um(melink(me, k)) = um(melink(ds_upme(me), k))
                        vn(melink(me, k)) = vn(melink(ds_upme(me), k))
                    endif
                enddo
                unsth(me) = unsth(ds_upme(me))
                cycle
            endif
        endif
        ! added by k.yamamura
        rr = rain(it, urain_i(me), urain_j(me))

        ! inflow by flux
        sumf = 0.0d0
        do k = 1, ko(me)
            umbeta(melink(me, k)) = um(melink(me, k))*rbeta(melink(me, k))
            vnbeta(melink(me, k)) = vn(melink(me, k))*rbeta(melink(me, k))
            f = umbeta(melink(me, k))*node_dy(me,k) - vnbeta(melink(me, k))*node_dx(me,k)
            sumf = sumf + f
        enddo

        ! calculate depth
        q(me) = dt2*(-sumf + rr*rnof(me)*smesh(me) + riv_eq(me))  ! fixed

        ! Update Depth
        unsth(me) = ho(me) + q(me)/( (1.0d0-lambda(me)) * smesh(me) )  ! fixed
        if(paddydam==1 .and. inf(me) == 71) then
            if(unsth(me) >= etp(paddyid(me))) then
                unsth(me) = unsth(me) - etp(paddyid(me))
                q(me) = q(me) - (etp(paddyid(me)) * smesh(me))
            endif
        endif
        if(baseo(me)==-9999.0d0) then
            unst_tmp_ds = unst_tmp_ds + unsth(me) * smesh(me)  ! v.1.0.5
            unsth(me) = 0.0d0
        endif

        ! 浸透モデル
        if( unst_soildepth(me) > 0.d0 .and. unst_ksv(me) > 0.d0 ) call unst_infilt(me)

        if(unsth(me)<0.0d0) unst_tmp_error = unst_tmp_error - unsth(me) * smesh(me)
        unsth(me) = max(unsth(me), 0.0d0)

        !　calculate sewerage and fields depth
        ! added by k.yamamura
        if(drainarea==1) then
        if(inf_dr(me)>0) then
                if(unsth(me)*smesh(me)<=vol_dr(me)) then
                    vol(me) = 0.0d0
                else
                    vol(me) = vol_dr(me)
                    unsth(me) = (unsth(me)*smesh(me) - vol(me)) / smesh(me)
                endif
                !Time-delayed discharge to the outlet point
                !$omp critical
                ilt = tripTime(me)
                if(ilt>0) then
                    do nn = 1, dr_no
                        dhj(ilt,nn) = dhj(ilt,nn) + vol(me)
                    enddo
                endif
                !$omp end critical
        endif
        endif

        ! calculate paddy field flow
        if(paddydam==1) then
            if(paddyid(me)>0) then
                !$omp atomic
                paddy_q(paddyid(me)) = paddy_q(paddyid(me)) + (unsth(me) * smesh(me))
            endif
        endif

        ! check water balance
        if(paddydam==1 .and. inf(me)==71) then
            ! 田んぼダムを考慮する場合、
            ! 田んぼダムのメッシュ上の qr_sum の計算はここでは行いません。
            continue
        else
            ! 下水道・圃場の場合の処理です。
            if(drainarea==1) then
                if(inf_dr(me) > 0) qr_sum(me) = qr_sum(me) + (q(me) - vol(me))
            else
                ! 下水道・圃場でない場合は、自由流出する下流端メッシュでなければ
                ! q(me) を加算します。
                if(baseo(me)/=-9999.0d0) qr_sum(me) = qr_sum(me) + q(me)
            endif
        endif

    enddo
    !$omp end parallel do

    unst_error_v = unst_error_v + unst_tmp_error  ! v.1.0.5
    unst_dis_v = unst_dis_v + unst_tmp_ds  ! v.1.0.5

    ! calculate paddy field dam outflow
    if (paddydam==1) call paddyflow
    if (ct /= nt .and. paddydam==1) call timedelay_paddyflow

    ! calculate sewerage and fields outflow
    if (ct /= nt .and. drainarea==1)  call drainflow

    deallocate(q)

end subroutine suisin

!----------------------------
! Calculate velocity
!----------------------------
subroutine velocity
    implicit none
    integer li, me, k, i

    ! hl calculation
    !$omp parallel do default(shared),private(li)
    do li = 1, link
        if(limesh(li, 2) == 0) then
        hl(li) = unsth(limesh(li, 1))
        else
        hl(li) = unsth(limesh(li, 1))*rthl(li, 1) + unsth(limesh(li, 2))*rthl(li, 2)
        endif
    enddo
    !$omp end parallel do

    ! umm, vnm calculation
    !$omp parallel do default(shared),private(me, k)
    do me = 1, mesh
        umm(me) = 0.0d0
        vnm(me) = 0.0d0
    do k = 1, ko(me)
        umm(me) = umm(me) + um(melink(me, k))*rtuv_x(me, k)
        vnm(me) = vnm(me) + vn(melink(me, k))*rtuv_y(me, k)
    enddo
    enddo
    !$omp end parallel do

    if(dsmesh==1) then
    !$omp parallel do default(shared),private(i, k)
    do i = 1, ids2mesh
        do k = 1, ko(ds2me(i))
            if(limesh(melink(ds2me(i), k),2)==0) then
                um(melink(ds2me(i), k)) = umm(ds2me(i))
                vn(melink(ds2me(i), k)) = vnm(ds2me(i))
                if(unsth(ds2me(i)) < th) then
                if((umm(ds2me(i))*node_dy(ds2me(i),k) - &
                     vnm(ds2me(i))*node_dx(ds2me(i),k)) > 0.0d0) then
                        um(melink(ds2me(i), k)) = 0.0d0
                        vn(melink(ds2me(i), k)) = 0.0d0
                endif
                endif
            endif
        enddo
    enddo
    !$omp end parallel do
    endif

    ! uu, vv calculation
    !$omp parallel do default(shared),private(li)
    do li = 1, link
        if(hl(li) > th) then
        uu(li) = um(li)/hl(li)
        vv(li) = vn(li)/hl(li)
        else
        uu(li) = 0.0d0
        vv(li) = 0.0d0
        endif
    enddo
    !$omp end parallel do

    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        if(unsth(me) < th) then
        uum(me) = 0.0d0
        vvm(me) = 0.0d0
        else
        uum(me) = umm(me)/unsth(me)
        vvm(me) = vnm(me)/unsth(me)
        endif
    enddo
    !$omp end parallel do

end subroutine velocity

!-------------------------
! Replace (Update) data
!-------------------------
subroutine replace
    implicit none
    integer li, me

    !$omp parallel do default(shared),private(li)
    do li = 1, link
        umo(li) = um(li)
        vno(li) = vn(li)
        lhano(li) = lhan(li)
    enddo
    !$omp end parallel do

    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        ho(me) = unsth(me)
    enddo
    !$omp end parallel do

end subroutine replace

!--------------------------------------
! Culculate inflow on the boundares
!--------------------------------------
subroutine lkyokai
    implicit none
    real(8) s_qin ! UNST-2D original
    integer j, ii

    ii = int(unsttime/dtq) + 1
    if(inls == 1) then
        !$omp parallel do default(shared),private(j, s_qin)
        do j = 1, iqnum
            s_qin = qin(j, ii) + (unsttime - dtq*dble(ii - 1))/dtq*(qin(j, ii + 1) - qin(j, ii))
            ! only unst-2D (Not used RRI to UNST)
            select case(lkyokai_dir(j))
            case(1)  ! from left side
                um(inl(j)) = s_qin/lkyokai_dy(j)
            case(2)  ! from down side
                vn(inl(j)) = s_qin/lkyokai_dx(j)
            case(3)  ! from right side
                um(inl(j)) = -s_qin/lkyokai_dy(j)
            case(4)  ! from upper side
                vn(inl(j)) = -s_qin/lkyokai_dx(j)
            case default
                ! no qin
            end select
        enddo
        !$omp end parallel do
    elseif(inls == 0) then
        s_qin = qin(1,ii) + (unsttime - dtq*dble(ii - 1))/dtq*(qin(1,ii + 1) - qin(1,ii))
        ! only unst-2D (Not used RRI to UNST)
        select case(lkyokai_dir(1))
        case(1)  ! from left side
            um(inl(1)) = s_qin/lkyokai_dy(1)
        case(2)  ! from down side
            vn(inl(1)) = s_qin/lkyokai_dx(1)
        case(3)  ! from right side
            um(inl(1)) = -s_qin/lkyokai_dy(1)
        case(4)  ! from upper side
            vn(inl(1)) = -s_qin/lkyokai_dx(1)
        case default
            ! no qin
        end select
    endif

end subroutine lkyokai

!----------------------------------
! Calculate infiltration process
!----------------------------------
subroutine unst_infilt(me)
    implicit none
    
    real(8) unst_gampt_ff_temp
    integer me
    
    unst_gampt_f(me) = 0.d0
    unst_gampt_ff_temp = unst_gampt_ff(me)
    if( unst_gampt_ff_temp .le. 0.01d0 ) unst_gampt_ff_temp = 0.01d0

    ! unst_gampt_f(me) : infiltration capacity [m/s]
    ! unst_gampt_ff : accumulated infiltration depth [m]
    unst_gampt_f(me) = unst_ksv(me) * (1.d0 + unst_faif(me) * unst_gammaa(me) / unst_gampt_ff_temp)

    ! unst_gampt_f(me) : infiltration capacity -> infiltration rate [m/s]
    if( unst_gampt_f(me) .ge. unsth(me) / unstdt ) unst_gampt_f(me) = unsth(me) / unstdt

    ! unst_gampt_ff should not exceeds a certain level
    if( unst_infilt_limit(me) .ge. 0.d0 .and. unst_gampt_ff(me) .ge. unst_infilt_limit(me) ) unst_gampt_f(me) = 0.d0

    ! update unst_gampt_ff [m]
    unst_gampt_ff(me) = unst_gampt_ff(me) + unst_gampt_f(me) * unstdt

    ! hs : hs - infiltration rate * dt [m]
    unsth(me) = unsth(me) - unst_gampt_f(me) * unstdt
    if( unsth(me) .le. 0.d0 ) unsth(me) = 0.d0
    
end subroutine unst_infilt

end module unst_cal_sub

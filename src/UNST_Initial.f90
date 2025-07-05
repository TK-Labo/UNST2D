! initialize conditions
! coded by k.kawaike and TK Labo
! Released on July 7th 2025

subroutine unst_initiald
    use globals
    implicit none
    integer me, li, k, k2, infn, j
    real(8) plant_lambda_link, mesh_dx, mesh_dy
    integer i, idummy !unst-2D only

    ! 初期水位、流量 read
    do i = 1, msctn
        read(15, *) idummy, hs(i), qs(i), qys(i)
        ! f_qs(i) = qs(i) !河道流量の下限値(連続水面の維持)
    enddo
    close(15)

    ! 河川
    !$omp parallel do default(shared),private(i)
    do i = 1, msctn
        hn(i) = hs(i)
        qn(i) = qs(i)
        bs(i) = 0.0d0
        sap(i) = 0.0d0
        sbm(i) = 0.0d0
    enddo
    !$omp end parallel do

    ! coefficient of roughness by inf
    !$omp parallel do default(shared),private(me, k, k2)
    do me = 1, mesh
        if(str_type == 0) then
            if(inf(me) == 26) then
                unsth(me) = 1.0d0
                unstc(me) = cr(me)
            else
                unsth(me) = 0.0d0
                unstc(me) = 0.0d0
            endif
        elseif(str_type == 1) then
            if(inf(me) == 26) then
                unsth(me) = 0.0d0
                unstc(me) = 0.0d0
            else
                unsth(me) = 0.0d0
                unstc(me) = 0.0d0
            endif
        endif
        
        ho(me) = unsth(me)
        hmax(me) = unsth(me)
        co(me) = unstc(me)!tracer
        qr_sum(me) = 0.0d0

        str(me) = 0.0d0    !貯留量
        stro(me) = str(me)
        rnof(me) = 1.00d0

        ! if(inf(me) ==  1) rnof(me) = 1.00d0   ! その他
        ! if(inf(me) == 24) rnof(me) = 1.00d0   ! その他河川
        ! if(inf(me) == 63) rnof(me) = 1.00d0   ! 田んぼダムを考慮しない水田
        ! if(inf(me) == 71) rnof(me) = 1.00d0   ! 田んぼダム(added by k.yamamura)
        ! if(inf(me) == 25) rnof(me) = 1.00d0   ! ポンプ
        ! if(inf(me) == 26) rnof(me) = 1.00d0   ! ため池
        ! if(inf(me) == 31) rnof(me) = 1.00d0   ! 事務所
        ! if(inf(me) == 32) rnof(me) = 1.00d0   ! 住宅
        ! if(inf(me) ==  4) rnof(me) = 1.00d0   ! 学校
        ! if(inf(me) ==  5) rnof(me) = 1.00d0   ! 公園
        ! if(inf(me) ==  6) rnof(me) = 1.00d0   ! 農地
        ! if(inf(me) ==  7) rnof(me) = 1.00d0   ! 山地
        ! if(inf(me) ==  8) rnof(me) = 1.00d0   ! 道路

        ! 貯留可能量
        if(str_type == 0) then
            strmx(me) = 0.0d0
        elseif(str_type == 1) then
            strmx(me) = 0.0d0
            if(inf(me) ==  4) strmx(me) = smesh(me)*0.4d0*0.3d0  ! 学校
            if(inf(me) ==  5) strmx(me) = smesh(me)*0.6d0*0.2d0   ! 公園
            ! if(inf(me) ==  1) strmx(me) = 0.0d0   ! その他(住宅地扱い)
            ! if(inf(me) == 21) strmx(me) = 0.0d0   ! 堀川
            ! if(inf(me) == 22) strmx(me) = 0.0d0   ! 比津川
            ! if(inf(me) == 23) strmx(me) = 0.0d0   ! 中川
            ! if(inf(me) == 24) strmx(me) = 0.0d0   ! その他河川
            ! if(inf(me) == 25) strmx(me) = 0.0d0   ! ポンプ
            ! if(inf(me) == 26) strmx(me) = 0.0d0   ! ため池
            ! if(inf(me) == 31) strmx(me) = 0.0d0   ! 事業所・商店
            ! if(inf(me) == 32) strmx(me) = 0.0d0   ! 住宅地
            ! if(inf(me) ==  6) strmx(me) = 0.0d0   ! 農地
            ! if(inf(me) == 61) strmx(me) = 0.0d0   ! 水田
            ! if(inf(me) == 62) strmx(me) = 0.0d0   ! 水路との接続がない水田
            ! if(inf(me) == 63) strmx(me) = 0.0d0   ! 田んぼダムを考慮しない水田
            ! if(inf(me) ==  7) strmx(me) = 0.0d0   ! 山地
            ! if(inf(me) ==  8) strmx(me) = 0.0d0   ! 道路
        endif

        ! added by k.yamamura
        if(plantDa==1) then
            if (plant_a_array(me) /= 0.0d0 .and. plant_D_array(me) /= 0.0d0) then
                plant_lambda(me) = plant_a_array(me)*plant_D_array(me)
            else
                plant_lambda(me) = 0.0d0
            endif
        endif

        do k = 1, ko(me)
            k2 = mod(k, ko(me)) + 1
            node_dx(me, k) = dnox(menode(me, k2)) - dnox(menode(me, k))
            node_dy(me, k) = dnoy(menode(me, k2)) - dnoy(menode(me, k))
        enddo

    enddo
    !$omp end parallel do

    !$omp parallel do default(shared),private(li, plant_lambda_link)
    do li = 1, link
        um(li)  = 0.0d0
        umo(li) = 0.0d0
        vn(li)  = 0.0d0
        vno(li) = 0.0d0
!     ---------------------------------- 補間水深
        if(limesh(li, 2) /= 0) then
        hl(li) = unsth(limesh(li, 1))*rthl(li, 1) + unsth(limesh(li, 2))*rthl(li, 2)
        else
        hl(li) = unsth(limesh(li, 1))
        endif
        lhan(li) = 0
        lhano(li) = lhan(li)

        dl(li) = 0.0d0
        rnx(li) = 0.0d0
        if(limesh(li,2) /= 0) then
            mesh_dx = xmesh(limesh(li, 2)) - xmesh(limesh(li, 1))
            mesh_dy = ymesh(limesh(li, 2)) - ymesh(limesh(li, 1))
            dl(li) = sqrt(mesh_dx**2 + mesh_dy**2)
            rnx(li) = 0.5d0*(mn(limesh(li, 1)) + mn(limesh(li, 2)))
        endif

        if(plantDa==1) then
            if(limesh(li,2)/=0) then
                plant_lambda_link = 0.5d0*(plant_lambda(limesh(li, 1))+plant_lambda(limesh(li, 2)))
            else
                plant_lambda_link = plant_lambda(limesh(li, 1))
            endif
            dk_val(li) = 0.0d0
            ! 抗力係数CD = 1.2d0
            if(plant_lambda_link>0.0d0) dk_val(li) = 1.0d0/sqrt(plant_lambda_link*1.2d0/(2.0d0*gg))
        endif
    enddo
    !$omp end parallel do

    do infn = 1, 100
        v_minus(infn) = 0.0d0    !属性ごとの水量調整量(負の水深をゼロにするときに増やした水量)
    enddo
    v_minus_all = 0.0d0    !負の水深をゼロにするときに増やした総水量
    v_cminus = 0.0d0    !負の濃度をゼロにするときに増やした総水量
    v_cplus = 0.0d0     !1より大きい濃度を1にするときに減らした総水量
    v_cextre = 0.0d0    !水深→0による極大の濃度をゼロにするときに増やした総水量 

    vrain = 0.0d0         !積算降雨量(外力の合計)
    qout = 0.0d0          !排水量(排水機場)
    dvr = 0.0d0           !河川からの流入量
    svc = 0.0d0           !流入トレーサー量
    sv0 = 0.0d0

    !$omp parallel do default(shared),private(me),reduction(+:sv0) 
    do me = 1, mesh
        if(inf(me) == 0) goto 101
        sv0 = sv0 + unsth(me)*smesh(me)*(1.0d0 - lambda(me))    !初期の総水量
        svc = svc + unstc(me)*unsth(me)*smesh(me)*(1.0d0 - lambda(me))
101 enddo
    !$omp end parallel do

    ! 上流端のdx,dy計算
    lkyokai_dx = 0.0d0
    lkyokai_dy = 0.0d0
    lkyokai_vect = 0

    ! UNST-2D original
    if(inls == 1) then
        !$omp parallel do default(shared),private(j)
        do j = 1, iqnum
            lkyokai_dx(j) = dnox(linode(inl(j), 1)) - dnox(linode(inl(j), 2))
            lkyokai_dy(j) = dnoy(linode(inl(j), 1)) - dnoy(linode(inl(j), 2))
            if (abs(lkyokai_dx(j)) < abs(lkyokai_dy(j))) lkyokai_vect(j) = 1
        enddo
        !$omp end parallel do
    elseif(inls == 0) then
        lkyokai_dx(1) = dnox(linode(inl(1), 1)) - dnox(linode(inl(1), 2))
        lkyokai_dy(1) = dnoy(linode(inl(1), 1)) - dnoy(linode(inl(1), 2))
        if (abs(lkyokai_dx(1)) < abs(lkyokai_dy(1))) lkyokai_vect(1) = 1
    endif

end subroutine unst_initiald

! 開発中(断面データから一次元用水理テーブルを自動算出) by d.baba
subroutine d1riv_table
    use globals

    integer i, j, k
    integer ihnummax
    real(8) dp_step                     ! 水位間隔[m]
    logical point1_judge, point2_judge  ! 水面との位置関係を判定
    real(8) section_dx, section_dz      ! 端点距離のx,z成分
    real(8) intersection_x              ! 水面との交点(x座標)
    real(8) tmp_dan_a, tmp_dan_s        ! 河積と潤辺
    real(8), allocatable :: tmp_x(:),tmp_z(:)  ! 断面測点座標

    dp_step = 0.1d0
    ihnummax = int(max_dz / dp_step) + 1

    allocate( ha(msctn, ihnummax), hr(msctn, ihnummax))
    allocate( hss(msctn, ihnummax), ah(msctn, ihnummax))

    write(*,*) '一次元河道 - 水理テーブル計算中...'

    ha = 0.0d0  ! 河積
    hr = 0.0d0  ! 径深
    hss = 0.0d0  ! 潤辺
    ah = 0.0d0  ! 水位

    allocate (tmp_x(ndan_num), tmp_z(ndan_num))
    !$omp parallel do default(shared), private(i)
    ! 断面毎に算出
    do i = 1, msctn
        tmp_x = 0.0d0
        tmp_z = 0.0d0
        ! 対象断面の座標を取得
        do j = 1, ndan(i)
            tmp_x(j) = dcoord_x(i,j)
            tmp_z(j) = dcoord_z(i,j)
        end do
        ! データ数を取得
        ihnum(i) = int((max_z(i)-min_z(j))/dp_step)+1

        do j = 1, ihnum(i)
            ! 水位データを取得
            ah(i,j) = min_z(j) + dble(j-1) * dp_step
            ! 断面特性を計算
            do k = 1, ndan(j)-1
                ! 測点が水面下であるかを判定
                point1_judge = tmp_z(k) <= ah(i,j)
                point2_judge = tmp_z(k+1) <= ah(i,j)
                ! 点間の距離を算出
                section_dx = tmp_x(k+1) - tmp_x(k)
                section_dz = tmp_z(k+1) - tmp_z(k)
                ! 状況に応じて算出
                if (point1_judge .and. point2_judge) then
                    ! 両測点が水面下の場合
                    ! 河積[m2]
                    tmp_dan_a = 0.5d0 * section_dx * &
                                (2.0d0 * ah(i,j) - tmp_z(k) - tmp_z(k+1))
                    ha(i,j) = ha(i,j) + tmp_dan_a
                    ! 潤辺[m]
                    tmp_dan_s = sqrt(section_dx**2 + section_dz**2)
                    hss(i,j) = hss(i,j) + tmp_dan_s
                elseif (point1_judge .or. point2_judge) then
                    ! 片方の測点のみが水面下の場合
                    ! 傾きがある場合のみ算出（しきい値以下の場合は水面上とする）
                    if (abs(section_dz) < 0.000001d0) cycle
                    if (point1_judge) then
                    ! 左岸側が水面下の場合
                        ! 水面との交点
                        intersection_x = tmp_x(k) + &
                                         (ah(i,j) - tmp_z(k)) * section_dx / section_dz
                        ! 河積
                        tmp_dan_a = 0.5d0 * (intersection_x - tmp_x(k)) * &
                                    (ah(i,j) - tmp_z(k))
                        ! 潤辺
                        tmp_dan_s = sqrt((intersection_x - tmp_x(k))**2 + &
                                         (ah(i,j) - tmp_z(k))**2)
                    else
                    ! 右岸側が水面下の場合
                        ! 水面との交点
                        intersection_x = tmp_x(k+1) + &
                                         (ah(i,j) - tmp_z(k+1)) * section_dx / section_dz
                        ! 河積
                        tmp_dan_a = 0.5d0 * (tmp_x(k+1) - intersection_x) * &
                                    (ah(i,j) - tmp_z(k+1))
                        ! 潤辺
                        tmp_dan_s = sqrt((tmp_x(k+1) - intersection_x)**2 + &
                                         (ah(i,j) - tmp_z(k+1))**2)
                    endif
                    ! 河積
                    ha(i,j) = ha(i,j) + tmp_dan_a
                    ! 潤辺
                    hss(i,j) = hss(i,j) + tmp_dan_s
                endif
            enddo
            ! 径深
            if (hss(i,j) > 0.0d0) hr(i,j) = ha(i,j) / hss(i,j)
        enddo
    enddo
    !$omp end parallel do
    deallocate(tmp_x, tmp_z, hss, ah)
    
end subroutine d1riv_table

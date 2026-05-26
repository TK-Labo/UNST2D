! initialize conditions
! coded by k.kawaike and TK Labo
! Released on July 7th 2025

subroutine unst_initiald
    use unst_globals_mod
    use unst_1d_main
    implicit none
    integer me, li, k, k2, j
    real(8) plant_lambda_link, mesh_dx, mesh_dy
    real(8) max_baseo

    !-----------------------
    ! Initialize time param
    !-----------------------
    unsttime = 0.0d0
    mstep = 0

    !---------------------
    ! Initialize variable
    !---------------------
    ! -- 01 mesh data --
    umm = 0.0d0
    vnm = 0.0d0
    riv_eq = 0.0d0  ! v.1.0.5
    unst_error_v = 0.0d0  ! v1.0.5
    unst_dis_v = 0.0d0  ! v1.0.5

    ! shape param
    !$omp parallel do default(shared),private(me, k, k2)
    do me = 1, mesh
        do k = 1, ko(me)
            k2 = mod(k, ko(me)) + 1
            node_dx(k, me) = dnox(menode(k2, me)) - dnox(menode(k, me))
            node_dy(k, me) = dnoy(menode(k2, me)) - dnoy(menode(k, me))
        enddo
    enddo
    !$omp end parallel do

    ! infomation
    !$omp parallel do default(shared),private(me, k, k2)
    do me = 1, mesh
        unsth(me) = 0.0d0
        if (dsmesh==1) then
            if (ds_inf(me)==2) unsth(me) = ds_wl(me,1)
        endif
        ho(me) = unsth(me)
        hmax(me) = unsth(me)
        qr_sum(me) = unsth(me) * smesh(me)  ! fix  v.1.0.5
        rnof(me) = 1.00d0
    enddo
    !$omp end parallel do

    ! option
    if(plantDa==1) then
        !$omp parallel do default(shared),private(me)
        do me = 1, mesh
            if (plant_a_array(me) /= 0.0d0 .and. plant_D_array(me) /= 0.0d0) then
                plant_lambda(me) = plant_a_array(me)*plant_D_array(me)  ! a_s 遮断面積
            else
                plant_lambda(me) = 0.0d0
            endif
        enddo
        !$omp end parallel do
        deallocate(plant_a_array, plant_D_array)
    endif
    ! ---

    ! -- 02 link data --
    um  = 0.0d0
    umo = 0.0d0
    vn  = 0.0d0
    vno = 0.0d0
    lhan = 0
    lhano = 0

    ! shape param
    !$omp parallel do default(shared),private(li, mesh_dx, mesh_dy)
    do li = 1, link
        dl(li) = 0.0d0
        if(limesh(li,2) /= 0) then
            mesh_dx = xmesh(limesh(li, 2)) - xmesh(limesh(li, 1))
            mesh_dy = ymesh(limesh(li, 2)) - ymesh(limesh(li, 1))
            dl(li) = sqrt(mesh_dx**2 + mesh_dy**2)
            blink(li) = sqrt((dnox(linode(li, 1)) - dnox(linode(li, 2))) **2 &
                             + (dnoy(linode(li, 1)) - dnoy(linode(li, 2)) ** 2))
        endif
    enddo
    !$omp end parallel do

    ! infomation
    !$omp parallel do default(shared),private(li)
    do li = 1, link
        if(limesh(li, 2) /= 0) then
            hl(li) = unsth(limesh(li, 1))*rthl(li, 1) + unsth(limesh(li, 2))*rthl(li, 2)
        else
            hl(li) = unsth(limesh(li, 1))
        endif
        rnx(li) = 0.0d0
        if(limesh(li,2) /= 0) rnx(li) = 0.5d0*(mn(limesh(li, 1)) + mn(limesh(li, 2)))
    enddo
    !$omp end parallel do

    ! option : plantDa
    if(plantDa==1) then
        !$omp parallel do default(shared),private(li, plant_lambda_link)
        do li = 1, link
            if(limesh(li, 2)/=0) cycle
            if(plant_lambda(limesh(li, 1))>0.0d0 .and. plant_lambda(limesh(li, 2))>0.0d0) then
                plant_lambda_link = 0.5d0*(plant_lambda(limesh(li, 1))+plant_lambda(limesh(li, 2)))
                vr_hv(li) = 0.5d0*(plant_hv_array(limesh(li, 1))+plant_hv_array(limesh(li, 2)))
            elseif(plant_lambda(limesh(li, 1))>0.0d0) then
                plant_lambda_link = plant_lambda(limesh(li, 1))
                vr_hv(li) = plant_hv_array((limesh(li, 1)))
            elseif(plant_lambda(limesh(li, 2))>0.0d0) then
                plant_lambda_link = plant_lambda(limesh(li, 2))
                vr_hv(li) = plant_hv_array(limesh(li, 2))
            endif

            dk_val(li) = 0.0d0
            !if(plant_lambda_link>0.0d0) dk_val(li) = 1.0d0/sqrt(plant_lambda_link*plant_cd/(2.0d0*gg))
            ! fixed  v.1.0.5  樹高有無の境界で樹高および占有が半減しないように調整
            dk_val(li) = 0.5d0 * vr_cd(li) * plant_lambda_link
        enddo
        !$omp end parallel do
        deallocate(plant_lambda, plant_hv_array, vr_cd)
    endif

    ! option : plantFN v.1.0.5.1
    ! if(plantFN==1) then
    !     do li = 1, link
    !         if(limesh(li,2)/=0) cycle
    !         ! N
    !         if(plantN_array(limesh(li, 1))>0.0d0 .or. plantN_array(limesh(li, 2))>0.0d0) then
    !             vr_N(li) = 0.5d0 * (plantN_array(limesh(li, 1)) + plantN_array(limesh(li, 2)))
    !         endif
    !         ! Al
    !         if(plantAl_array(limesh(li, 1))>0.0d0 .or. plantAl_array(limesh(li, 2))>0.0d0) then
    !             vr_Al(li) = 0.5d0 * (plantAl_array(limesh(li, 1)) + plantAl_array(limesh(li, 2)))
    !         endif
    !         ! l
    !         if(plantl_array(limesh(li, 1))>0.0d0 .or. plantl_array(limesh(li, 2))>0.0d0) then
    !             vr_l(li) = 0.5d0 * (plantl_array(limesh(li, 1)) + plantl_array(limesh(li, 2)))
    !         endif
    !         ! d
    !         if(plantdd_array(limesh(li, 1))>0.0d0 .or. plantdd_array(limesh(li, 2))>0.0d0) then
    !             vr_dd(li) = 0.5d0 * (plantdd_array(limesh(li, 1)) + plantdd_array(limesh(li, 2)))
    !         endif
    !         vr_ld(li) = 0.5d0*vr_l(li)*vr_dd(li)
    !         ! C_d
    !         if(stemCd_array(limesh(li, 1))>0.0d0 .or. stemCd_array(limesh(li, 2))>0.0d0) then
    !             stem_cd(li) = 0.5d0 * (stemCd_array(limesh(li, 1)) + stemCd_array(limesh(li, 2)))
    !         endif
    !         ! C_dl
    !         if(leavesCdl_array(limesh(li, 1))>0.0d0 .or. leavesCdl_array(limesh(li, 2))>0.0d0) then
    !             leaves_cdl(li) = 0.5d0 * (leavesCdl_array(limesh(li, 1)) + leavesCdl_array(limesh(li, 2)))
    !         endif
    !         ! C_sl
    !         if(leavesCsl_array(limesh(li, 1))>0.0d0 .or. leavesCsl_array(limesh(li, 2))>0.0d0) then
    !             leaves_csl(li) = 0.5d0 * (leavesCsl_array(limesh(li, 1)) + leavesCsl_array(limesh(li, 2)))
    !         endif
    !     enddo
    ! endif

    ! option : morido
    if(mmorid==1) then
        !$omp parallel do default(shared),private(li,max_baseo)
        do li = 1, link
            if(infl(li)==1) then
                max_baseo = max(baseo(limesh(li, 1)), baseo(limesh(li, 2)))
                if(zbbk(li) <= max_baseo) then
                    infl(li) = 0
                    write(*,*) ' UNST2D - WARNING : invalid Embankment height linkid  >>> linkid', li
                endif
            endif
        enddo
        !$omp end parallel do
    endif
    ! ---

    ! -- 03 qin data --
    lkyokai_dir = 0
    lkyokai_dx = 0.0d0
    lkyokai_dy = 0.0d0
    !$omp parallel do default(shared),private(j, k)
    do j = 1, iqnum
            lkyokai_dx(j) = dnox(linode(inl(j), 1)) - dnox(linode(inl(j), 2))
            lkyokai_dy(j) = dnoy(linode(inl(j), 1)) - dnoy(linode(inl(j), 2))
            ! v.1.0.5
            if(abs(lkyokai_dx(j)) > abs(lkyokai_dy(j))) then
                if(lkyokai_dx(j) < 0) lkyokai_dir(j) = 2  ! from down side
                if(lkyokai_dx(j) > 0) lkyokai_dir(j) = 4  ! from upper side
            else
                if(lkyokai_dy(j) < 0) lkyokai_dir(j) = 3  ! from right side
                if(lkyokai_dy(j) > 0) lkyokai_dir(j) = 1  ! from left side
            endif
            lkyokai_dx(j) = abs(lkyokai_dx(j))
            lkyokai_dy(j) = abs(lkyokai_dy(j))
            if(lkyokai_dir(j)==0) write(*,*) ' UNST2D - WARNING : No length qin  >>> linkid', inl(j)
    enddo
    !$omp end parallel do

    if(d1riv==1) call d1riv_spinup

end subroutine unst_initiald

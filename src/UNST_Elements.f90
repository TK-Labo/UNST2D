! River Basin Disaster Resilience and Sustainability by All
! coded by k.yamamura
! added by d.baba
!

! paddy field dam
subroutine paddyinitiald
    use globals
    implicit none
    integer i, pa, k, ii
    real(8), allocatable :: hpa(:), gradp(:), speedp(:) 
    real(8), allocatable :: ttpA(:), ttpB(:)

    allocate(hpa(mesh), gradp(mesh), speedp(mesh))
    allocate(ttpA(mesh), ttpB(mesh))
    allocate(paddycount_num(paddy))

    ! initialize
    totalqp = 0.0d0
    pqh = 0.0d0
    
    ! area / paddy field & paddycount_num
    allocate(psmesh(paddy))
    max_paddycount_num = 0
    paddycount_num = 0
    psmesh = 0.0d0
    do i = 1, mesh
        if (paddyid(i) > 0) then
        psmesh(paddyid(i)) = psmesh(paddyid(i)) + smesh(i)
        paddycount_num(paddyid(i)) = paddycount_num(paddyid(i)) + 1
        max_paddycount_num = max(max_paddycount_num, paddycount_num(paddyid(i)))
        endif
    enddo

    ! paddyid2mesh
    allocate(paddyid2mesh(paddy, max_paddycount_num))
    paddyid2mesh = 0
    do i = 1, mesh
        if (paddyid(i) > 0) then
            do k = 1, paddycount_num(paddyid(i))
                if (paddyid2mesh(paddyid(i), k) == 0) then
                    paddyid2mesh(paddyid(i), k) = i
                    exit
                endif
            enddo
        endif
    enddo

    !depth of submergence
    etp = (etp / 86400.0d0 / 1000) * dt2
    
    ! calculate distance from paddy field to waterway
    hpa = 0.0d0
    !$omp parallel do default(shared),private(pa)
    do pa = 1, paddy
        if (min_dist(pa) > 0.001d0) then
            ! elevation difference btw paddy field and nearest waterway
            hpa(pa) = abs(baseo(pqout_idx(pa)) - baseo(min_pmeshid(pa)))
            ! gradient
            gradp(pa) = hpa(pa) / min_dist(pa)
            ! calculated velocity（クラーヘン式の洪水伝播速度）
            if (gradp(pa) >= 0.01d0) speedp(pa) = 3.5d0                             ! 1/100 <= gradient
            if (0.01d0 > gradp(pa) .and. gradp(pa) >= 0.005d0) speedp(pa) = 3.0d0   ! 1/100 ~ 1/200 gradient 
            if (gradp(pa) < 0.005d0) speedp(pa) = 2.1d0
            ttpA(pa) = min_dist(pa) / speedp(pa)
        else
            ttpA(pa) = 0.0d0
        endif
        !　elevation difference btw nearest waterway and draining point
        hpa(pa) = abs(baseo(pqout_idx(pa)) - baseo(pdrain(pa)))
        ! gradient
        gradp(pa) = hpa(pa) / dr_dist(pa)
        ! velocity
        if (gradp(pa) >= 0.01d0) speedp(pa) = 3.5d0                             ! 1/100 <= gradient
        if (0.01d0 > gradp(pa) .and. gradp(pa) >= 0.005d0) speedp(pa) = 3.0d0   ! 1/100 ~ 1/200 gradient 
        if (gradp(pa) < 0.005d0) speedp(pa) = 2.1d0
        ttpB(pa) = dr_dist(pa) / speedp(pa)
        ttp(pa) = int(ttpA(pa) + ttpB(pa))
    enddo
    !$omp end parallel do

    ! drain_meshid for paddynumber
    allocate(drain2phidx(paddy))
    drain2phidx = 0
    do pa = 1, paddy
        do ii = 1, nhp
            if (phid(ii) == pdrain(pa)) then
                drain2phidx(pa) = ii
                exit  ! 見つかったら内側のループを抜ける
            endif
        enddo
        if(drain2phidx(pa) == 0) then
            print *, "error: please enter pdrain", pdrain(pa)
            read(*,*)
        endif
    enddo

    outa = (((p_data / 2.0d0) ** 2.0d0) * pi)

    deallocate(hpa, gradp, speedp, ttpA, ttpB)

end subroutine paddyinitiald


! paddy field dam outflow calculation
subroutine timedelay_paddyflow
    use globals
    implicit none
    integer j, ii, k, idx
    real(8) temp_qin

    ! add time delay
    !$omp parallel do default(shared) private(k, idx, temp_qin)
    do k = 1, nhp
        idx = phid(k)
        temp_qin = dhp(1,k) / smesh(idx)
        unsth(idx) = max(unsth(idx) + temp_qin, 0.0d0)
        if(unsth(idx)/=0.0d0) qr_sum(idx) = qr_sum(idx) + dhp(1,k)
    enddo
    !$omp end parallel do
    do j = 1, nhp
        do ii = 1, 71999
            dhp(ii, j) = dhp(ii+1, j)
        enddo
    enddo

end subroutine timedelay_paddyflow      


subroutine paddyflow
    use globals
    implicit none
    integer pa, ii, ilt, k, j
    real(8) qp1, qp2, qp3, qp, total_depth, hhp

    !$omp parallel do default(shared) &
    !$omp private(pa, ii, ilt, k, j) &
    !$omp private(qp1, qp2, qp3, qp, total_depth, hhp)
    do pa = 1, paddy
        hhp = paddy_q(pa) / psmesh(pa)
        pqh(pa) = hhp
        if (wh1 < hhp .and. hhp < lh) then
            ! device=0(not paddy field dam)
            if (device(pa) == 0) then
                qp1 = &
                    0.6d0 * (2.0d0 / 3.0d0) * sqrt(2.0d0 * gg) * ww1 * ((hhp - wh1) ** (3.0d0 / 2.0d0))
            ! device=1(機能一体型)
            elseif (device(pa) == 1) then
                if (hhp <= (wh1 + wh2)) then
                    if (wtyp == 1) then
                        qp1 = &
                            0.6d0 * (2.0d0 / 3.0d0) * sqrt(2.0d0 * gg) * ww2 * ((hhp - wh1) ** (3.0d0 / 2.0d0))
                    elseif (wtyp == 2) then
                        qp1 = 0.6d0 * (8.0d0 / 15.0d0) * tan((ca * pi) / 180.0d0) * sqrt(2.0d0 * gg) &
                                * ((hhp - wh1) ** (5.0d0 / 2.0d0))
                    endif
                elseif (hhp > (wh1 + wh2)) then
                    if (wtyp == 1) then
                        qp1 = 0.6d0 * (2.0d0 / 3.0d0) * sqrt(2.0d0 * gg) * ww2 * (wh2 ** (3.0d0 / 2.0d0)) + &
                            0.6d0 * (2.0d0 / 3.0d0) * sqrt(2.0d0 * gg) * ww1 * ((hhp - (wh1 + wh2)) ** (3.0d0 / 2.0d0))
                    elseif (wtyp == 2) then
                        qp1 = 0.6d0 * (8.0d0 / 15.0d0) * tan((ca * pi) / 180.0d0) * sqrt(2.0d0 * gg) * (wh2 ** (5.0d0 / 2.0d0)) + &
                            0.6d0 * (2.0d0 / 3.0d0) * sqrt(2.0d0 * gg) * ww1 * ((hhp - (wh1 + wh2)) ** (3.0d0 / 2.0d0))
                    end if
                endif
            ! device=2(機能分離型)
            elseif (device(pa) == 2) then
                if (hhp <= (lh - dld)) then
                    qp1 = &
                        0.6d0 * (2.0d0 / 3.0d0) * sqrt(2.0d0 * gg) * ww1 * ((hhp - wh1) ** (3.0d0 / 2.0d0))
                    qp2 = 0.6d0 * (((dd / 2.0d0) ** 2.0d0) * pi) * sqrt(2.0d0 * gg * (hhp + unstdh))                      
                    qp1 = min(qp1, qp2)
                elseif (hhp >= (lh - dld)) then
                    qp1 = 0.6d0 * (((dd / 2.0d0) ** 2.0d0) * pi) * sqrt(2.0d0 * gg * (hhp + unstdh)) + &
                    0.6d0 * (2.0d0 / 3.0d0) * sqrt(2.0d0 * gg) * ww1 * ((hhp - (lh - dld)) ** (3.0d0 / 2.0d0))
                end if
            else
                print *, "error: please enter device"
                read(*,*)
            endif
            ! Orifice outflow
            qp3 = 0.6d0 * outa * sqrt(2.0d0 * gg * (hhp + ph))
            qp = min(qp1, qp3)
            totalqp(pa) = (qp * dt2) * orifice_num(pa)
            total_depth = (paddy_q(pa) - totalqp(pa)) / psmesh(pa)
            pqh(pa) =  max(total_depth, 0.0d0)
            if(pqh(pa) == 0.0d0) totalqp(pa) = paddy_q(pa)
        
            ilt = ttp(pa) !time delay（s）
            if(ilt > 0) then
                ii = drain2phidx(pa)
                !$omp atomic
                dhp(ilt, ii) = dhp(ilt, ii) + totalqp(pa)
            endif
        end if

        ! update paddy field water depth
        do k = 1, paddycount_num(pa)
            j = paddyid2mesh(pa, k)
            unsth(j) = pqh(pa)
            qr_sum(j) = qr_sum(j) + ((unsth(j) - ho(j)) * smesh(j))
        enddo
    enddo
    !$omp end parallel do

end subroutine paddyflow


! sewerage and fields are 
subroutine draininitiald
    use globals
    implicit none

    integer me
    real(8), allocatable :: heightA(:), heightB(:), distanceA(:), distanceB(:)
    real(8), allocatable :: gradA(:), gradB(:), speedA(:), speedB(:)
    real(8), allocatable :: tripA(:), tripB(:)
    !

    allocate(heightA(mesh), heightB(mesh), distanceA(mesh), distanceB(mesh))
    allocate(gradA(mesh), gradB(mesh), speedA(mesh), speedB(mesh))
    allocate(tripA(mesh), tripB(mesh))

    ! distance from area to drp
    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        vol_dr(me) = 0.0d0
        if (inf_dr(me) > 0) then

            vol_dr(me) = drr(inf_dr(me))*1.0d-3*smesh(me) / 3600.0d0 * dt2

            distanceA(me) = 0.0d0
            distanceB(me) = 0.0d0

            distanceA(me) = sqrt((xmesh(me) - xmesh(drc(inf_dr(me))))**2 &
                                + (ymesh(me) - ymesh(drc(inf_dr(me))))**2)
            ! current mesh to drc
            distanceB(me) = drr_dist(inf_dr(me))
            ! drc to drp					
            ! height
            heightA(me)= 0.0d0
            heightB(me)= 0.0d0
            heightA(me) = abs(baseo(me) - baseo(drc(inf_dr(me))))
            heightB(me) = abs(baseo(drc(inf_dr(me))) - baseo(drp(inf_dr(me))))

            ! gradient
            if(heightA(me) == 0.0d0) then
                gradA(me) = 0.0d0
            else
                gradA(me) = heightA(me) / distanceA(me)
            endif
            if(heightB(me) == 0.0d0) then
                gradB(me) = 0.0d0
            else 
                gradB(me) = heightB(me) / distanceB(me)
            endif
            ! velocity
            if (gradA(me) >= 0.01d0) speedA(me) = 3.5d0                             ! 1/100 <= gradient
            if (0.01d0 > gradA(me) .and. gradA(me) >= 0.005d0) speedA(me) = 3.0d0   ! 1/100 ~ 1/200 gradient
            if (gradA(me) < 0.005d0) speedA(me) = 2.1d0                             

            if (gradB(me) >= 0.01d0) speedB(me) = 3.5d0
            if (0.01d0 > gradB(me) .and. gradB(me) >= 0.005d0) speedB(me) = 3.0d0
            if (gradB(me) < 0.005d0) speedB(me) = 2.1d0

            ! time delay
            tripA(me) = distanceA(me) / speedA(me)
            tripB(me) = distanceB(me) / speedB(me)

            ! time delay = tripA + tripB
            tripTime(me) = int(tripA(me) + tripB(me))
        else
            tripTime(me) = 0
        endif
    enddo
    !$omp end parallel do

    deallocate(drc)
    deallocate(heightA, heightB, distanceA, distanceB)
    deallocate(gradA, gradB, speedA, speedB)
    deallocate(tripA, tripB)

end subroutine draininitiald


! drain outflow calculation
subroutine drainflow
    use globals
    implicit none
    integer(8) nn, idx, j, ii
    real(8) temp_qin

    ! add time delay
    !$omp parallel do default(shared) private(nn, idx, temp_qin)
    do nn = 1, dr_no
        idx = drp(nn)
        temp_qin = dhj(1,nn) / smesh(idx)
        unsth(idx) = max(unsth(idx) + temp_qin, 0.0d0)
        if(unsth(idx)/=0.0d0) qr_sum(idx) = qr_sum(idx) + dhj(1,nn)
    enddo
    !$omp end parallel do
    ! time delay set
    do j = 1, dr_no
        do ii = 1, 71999
            dhj(ii, j) = dhj(ii+1, j)
        enddo
    enddo   
!      
end subroutine drainflow 


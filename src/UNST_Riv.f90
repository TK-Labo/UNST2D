! 1D River Model
! coded by TK Labo (D.Baba)
! Released on July 7th 2025

! 上流端･下流端境界条件（ハイドログラフ） read
subroutine bndry(iqmax, ihmax)
    use globals
    integer i, j, iup, iqmax, idw, ihmax
    character*50 tdummy

    ! qhyd ：上流端流入流量
    read(14, *)
    read(14, *)
    read(14, *) iup, iqmax

    ! allocate
    allocate( qhyd(iup, iqmax) )

    if(iup==0) goto 201
    do i = 1, iup
        read(14, *)
        do j = 1, iqmax
            read(14, *) tdummy, tdummy, qhyd(i, j)
        enddo
    enddo

    ! hhyd ：下流端潮位 (o.p.+2.2m = t.p.+0.9m 一定)
201 read(14, *)
    read(14, *)
    read(14, *) idw, ihmax

    ! allocate
    allocate( hhyd(idw, ihmax) )

    if(idw==0) goto 202
    do i = 1, idw
        read(14, *)
    do j = 1, ihmax
        read(14, *) tdummy, tdummy, hhyd(i, j)
    enddo
    enddo
202 close(14)

end subroutine bndry


! ハイドログラフの内挿
subroutine hydro
    use globals
    real(8) dtp, qqq, hhh
    integer iqhy, ihhy, iqmax, ihmax, it2
    
    ! 流量ハイドログラフの内挿
    entry qhydr(qqq, iqhy, iqmax)
        dtp = time/3600.0d0
        it2 = int(dtp) + 1
        if(it2 >= iqmax) it2 = iqmax - 1
        dtp = dtp - dble(it2 - 1)
        qqq = (qhyd(iqhy, it2 + 1) - qhyd(iqhy, it2))*dtp + qhyd(iqhy, it2)
    return

    ! 水位ハイドログラフの内挿
    entry hhydr(hhh, ihhy, ihmax)
        dtp = time/3600.0d0
        it2 = int(dtp) + 1
        if(it2 >= ihmax) it2 = ihmax - 1
        dtp = dtp - dble(it2 - 1)
        hhh = (hhyd(ihhy, it2 + 1) - hhyd(ihhy, it2))*dtp + hhyd(ihhy, it2)
    return

end subroutine hydro


subroutine soltn
    use globals
    real(8) cpm, a, v, r, b, an
    real(8) qcs, bvc, dhh
    integer i, iqhy, ihhy, iqmax, ihmax, inodek, inode, id

    ! ノードをもたない断面の水位・流量
    entry links
        !$omp parallel do default(shared),private(cpm, a, v, i)
        do i = 1, msctn
            if(iptn(i) /= 0) goto 10
            cpm = cap(i) - cbm(i)
            a = as(i) + (cap(i)*sbm(i) - cbm(i)*sap(i))/cpm
            v = us(i) - cap(i)*cbm(i)/as(i)*(sap(i) - sbm(i))/cpm
            !$omp critical
            call invaz(i, a, hn(i))
            !$omp end critical
            qn(i) = a*v
     10 enddo
        !$omp end parallel do
    return

    ! 外部ノードに与える境界条件
    entry bnodes(iqmax, ihmax)
        iqhy = 0
        ihhy = 0
        do inodek = knode1 + 1, knode
            i = nd(inodek, 1)

            if(inodek == knode1 + 1) goto 32   ! 上流端流量（嵩見橋）
            if(inodek == knode1 + 2 .and. iptn(i) == -1) goto 31   ! 上流端水位（宍道湖）
            if(inodek == knode1 + 3 .and. iptn(i) ==  1) goto 36   ! 下流端水位（中海）
            !if(inodek == knode1 + 4 .and. iptn(i) ==  1) goto 37   ! 下流端流量（手貝水門上）
            !if(inodek == knode1 + 5 .and. iptn(i) == -1) goto 32   ! 上流端流量（手貝水門下）
            goto 30

            ! 上流端：水位
        31 ihhy = ihhy + 1
            call hhydr(hn(i), ihhy, ihmax)
            call abrz(i, hn(i), a, r, b)
            qn(i) = a*(us(i) + cap(i)/as(i)*(a - as(i) - sbm(i)))
            goto 30

            ! 上流端：流量
        32 iqhy = iqhy + 1
            call qhydr(qn(i),iqhy,iqmax)
            an = as(i) + (qn(i) - qs(i) + cap(i)*sbm(i))/(us(i) + cap(i))
            call invaz(i, an, hn(i))
            goto 30

            ! 上流端：水位-流量曲線
        33 continue
            !call uphq( time, inode, hs(i), qq, dqdh )
            !hn(i)=hs(i)+(qq-qs(i)+sbm(i)*cs(i))/((us(i)+cs(i))*bs(i)-dqdh)
            !call uphq( time, inode, hn(i), qn(i), dqdh)
            goto 30

            !下流端：水位
        36 ihhy = ihhy + 1
            call hhydr(hn(i), ihhy, ihmax)
            call abrz(i, hn(i), a, r, b)
            qn(i) = a*(us(i) + cbm(i)/as(i)*(a - as(i) - sap(i)))
            goto 30

            !下流端：流量
        37 iqhy = iqhy + 1
            call qhydr(qn(i), iqhy, iqmax)
            an = as(i) + (qn(i) - qs(i) + cbm(i)*sap(i))/(us(i) + cbm(i))
            call invaz(i, an, hn(i))
            goto 30

            ! 下流端：水位-流量曲線
        38 continue
            !call dnhq( time, inode, hs(i), qq, dqdh )
            !hn(i)=hs(i)+(qq-qs(i)-sap(i)*cs(i))/((us(i)-cs(i))*bs(i)-dqdh)
            !call dnhq( time, inode, hn(i), qn(i), dqdh)
            goto 30
     30 enddo
    return

    !内部ノード（分・合流点）の水位・流量
    entry inodes
    if(knode1 == 0) goto 40
    do inode = 1, knode1
        qcs = 0.0d0
        bvc = 0.0d0
        do id = 1, mnd(inode)
            i = nd(inode, id)
            ! 下流端ノード 
            if(iptn(i) == 1) then
                qcs = qcs + (qs(i) - cbm(i)*sap(i))
                bvc = bvc + (bs(i)*(us(i) + cbm(i)))
            !上流端ノード 
            elseif(iptn(i) == -1) then
                qcs = qcs - (qs(i) - cap(i)*sbm(i))
                bvc = bvc - (bs(i)*(us(i) + cap(i)))
            endif
        enddo

        dhh = -qcs/bvc

        do id = 1, mnd(inode)
            i = nd(inode, id)
            hn(i) = hs(i) + dhh
            call abrz(i, hn(i), a, r, b)
            ! 下流端ノード 
            if(iptn(i) == 1) then
                qn(i) = (qs(i) - cbm(i)*sap(i)) + (us(i) + cbm(i))*(a-as(i))
            ! 上流端ノード 
            elseif(iptn(i) == -1) then
                qn(i) = (qs(i) - cap(i)*sbm(i)) + (us(i) + cap(i))*(a - as(i))
            endif
        enddo
    enddo

 40 return

end subroutine soltn


subroutine elm(alpg2, alpbet, cslgb)
    use globals
    real(8) alpg2, alpbet, cslgb, a, r, b, cs
    integer i, j

    !$omp parallel do default(shared),private(a, r, b, cs, i)
    do i = 1, msctn
        !$omp critical
        call abrz(i, hn(i), a, r, b)
        !$omp end critical
        if(hn(i) - hb(i) <= 0.0d0) hn(i) = hb(i) + 1.0d-3
        hs(i) = hn(i)
        as(i) = a
        qs(i) = qn(i)
        us(i) = qn(i)/a
        bs(i) = b
        htes(i) = us(i)**2*alpg2 + hn(i)
        fs(i) = rn(i)**2*us(i)*abs(us(i))/r**1.3333333
        cs = sqrt((alpbet*us(i))**2 + cslgb*a/b)  
        cap(i) = alpbet*us(i) + cs
        cbm(i) = alpbet*us(i) - cs
    enddo
    !$omp end parallel do

    !$omp parallel do default(shared),private(i)
    do i = 1, msctn
        ! upstream side
        if(iptn(i) == -1) goto 16
        sap(i) = -dt2*((qs(i) - qs(i - 1))/dx2(i - 1) - qys(i - 1) &
            + bs(i)*cap(i)/cslmd*((htes(i) - htes(i - 1))/dx2(i - 1) + (fs(i - 1) + fs(i))*0.5d0))
        ! downstream side
    16  if(iptn(i) ==  1) goto 15
        sbm(i) = -dt2*((qs(i + 1) - qs(i))/dx2(i) - qys(i) &
            + bs(i)*cbm(i)/cslmd*((htes(i + 1) - htes(i))/dx2(i) + (fs(i) + fs(i+1))*0.5d0))
 15 enddo
    !$omp end parallel do

end subroutine elm


subroutine hahr
    use globals
    real(8) rh, hh, hd, a, r, b, dif
    integer i, idm, ii

    ! 断面iの水位が与えられたとき
    ! 断面積a, 径深r, 水路幅b を求める
    entry abrz(i, hh, a, r, b)
        hd = hh - hb(i)
        if(hd < th) goto 411
        idm = int(hd/0.1d0) + 1
        ! 断面積、径深
        if(idm <= ihnum(i)-1) then
            dif = hd - dble(idm - 1)*0.1d0
            a = ha(i, idm) + (ha(i, idm + 1) - ha(i, idm))/0.1d0*dif
            r = hr(i, idm) + (hr(i, idm + 1) - hr(i, idm))/0.1d0*dif
        else
            dif = hd - dble(ihnum(i))*0.1d0
            a = ha(i, ihnum(i)) + (ha(i, ihnum(i)) - ha(i, ihnum(i) - 1))/0.1d0*dif
            r = hr(i, ihnum(i)) + (hr(i, ihnum(i)) - hr(i, ihnum(i) - 1))/0.1d0*dif
        endif
        ! 川幅
        if(abs(a - as(i)) < 1.0d-5 .or. abs(hn(i) - hs(i)) < 1.0d-7) then
            b = bs(i)
            if(b == 0.0d0) b = a/hd
        else
            b = (a - as(i))/(hn(i) - hs(i))
        endif
    return
    !
    411 write(*,412) i, hh, time
    412 format(1h /1h ,'too low stage  i, hh=', i5, f10.4, 'time=', f10.2)
        a = 0.0d0
        b = 0.0d0
        r = 0.0d0
    return

    ! 断面iの面積a が与えられたとき
    ! 水位hh を求める
    entry invaz(i, a, hh)
        do ii = 1, ihnum(i)
            if(ha(i, ii) > a) then
                rh = dble(ii - 2)*0.1d0 + (a - ha(i, ii - 1))/(ha(i, ii) - ha(i, ii - 1))*0.1d0
                goto 410
            endif
        enddo
410 hh = rh + hb(i)

end subroutine hahr


! honma weir equation added by d.baba
subroutine honma_weir(h1,h2,rivb,q0)
    use globals
    real(8), intent(in) :: h1
    real(8), intent(in) :: h2
    real(8), intent(in) :: rivb
    real(8), intent(out) :: q0
    
    if (h1 == h2) then
        q0 = 0.0d0
    elseif(h2*3.0d0 < h1*2.0d0) then
        ! kanzen
        q0 = 0.35d0 * rivb * sqrt(2.0d0*gg*h1)
    else
        ! moguri
        q0 = 0.91d0 * rivb * h2 * sqrt(2*(h1-h2))
    endif
    
end subroutine honma_weir

! weir calculation added by d.baba
subroutine weir_equation
    use globals
    real(8) tmp_h1, tmp_h2, h1, h2, wir_inout, q0, weir_ratio

    !$omp parallel do default(shared)&
    !$omp private(tmp_h1, tmp_h2, h1, h2, wir_inout, q0, weir_ratio)
    do i = 1, ksetu
        ! left levee
        if (hs(kkasen(i)) < lcrown((kkasen(i))) .and. &
            (baseo(kkousi(i)) + unsth(kkousi(i)) < lcrown((kkasen(i))))) goto 420
        tmp_h1 = max(baseo(kkousi(i)) + unsth(kkousi(i)) - lcrown((kkasen(i))), 0.0d0)
        tmp_h2 = max(hs(kkasen(i)) - lcrown((kkasen(i))), 0.0d0)
        if (tmp_h1>tmp_h2) then
            h1 = tmp_h1
            h2 = tmp_h2
            wir_inout = -1.0d0
        else
            h1 = tmp_h2
            h2 = tmp_h1
            wir_inout = 1.0d0
        endif
        call honma_weir(h1, h2, dx2(kkasen(i)), q0)
        if (q0 > 0.0d0) then
            dq_l(kkasen(i)) = wir_inout * wir_alpha(kkasen(i)) * q0 * wir_angle(kkasen(i))
        else
            dq_l(kkasen(i)) = 0.0d0
        endif

        ! right levee
420     if (hs(kkasen(i)) < rcrown((kkasen(i))) .and. &
            (baseo(kkousi(i)) + unsth(kkousi(i)) < rcrown((kkasen(i))))) goto 421
        tmp_h1 = max(baseo(kkousi(i)) + unsth(kkousi(i)) - rcrown((kkasen(i))), 0.0d0)
        tmp_h2 = max(hs(kkasen(i)) - rcrown((kkasen(i))), 0.0d0)
        if (tmp_h1>tmp_h2) then
            h1 = tmp_h1
            h2 = tmp_h2
            wir_inout = -1.0d0
        else
            h1 = tmp_h2
            h2 = tmp_h1
            wir_inout = 1.0d0
        endif
        call honma_weir(h1, h2, dx2(kkasen(i)), q0)
        if (q0 > 0.0d0) then
            dq_r(kkasen(i)) = wir_inout * wir_alpha(kkasen(i)) * q0 * wir_angle(kkasen(i))
        else
            dq_r(kkasen(i)) = 0.0d0
        endif
        
421     qys(kkasen(i)) = (dq_l(kkasen(i))+dq_r(kkasen(i))) * dt2
        ! hosei
        if (qn(kkasen(i)) - f_qs(kkasen(i)) < qys(kkasen(i))) then
            weir_ratio = (qn(kkasen(i)) - f_qs(kkasen(i))) / qys(kkasen(i))
            dq_l(kkasen(i)) = dq_l(kkasen(i)) * weir_ratio
            dq_r(kkasen(i)) = dq_r(kkasen(i)) * weir_ratio
            qys(kkasen(i)) = dq_l(kkasen(i))+dq_r(kkasen(i))
        endif
    enddo
    !$omp end parallel do

end subroutine weir_equation


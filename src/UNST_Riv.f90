! 1D River Model
! Coded by TK Labo (D.Baba)
! Released on July 12th 2025

!================================
! Equation of motion, continuity
!================================
!--------------------------------------------------------------------------
! 1D Model 主な変数/1D Model Main variable
!   leap-flogに合わせてn -> n+2 で表記
!   h_1d(ndan): 断面水位(t=n+2)/cross-sectional depth of a river(t=n+2)
!   ho_1d(ndan): 断面水位(t=n)/cross-sectional depth of a river(t=n)
!   a_1d(ndan): 断面河積/cross-sectional area of a river
!   r_1d(ndan): 断面径深/cross-sectional hydraulic radius of a river
!   rn_1d(ndan): 断面粗度/cross-sectional roughness of a river
!   vv_1d(ndan): 断面流速/cross-sectional velocity of a river
!   q_1d(ndan): 断面流量/cross-sectional flow rate of a river(t=n+2)
!   qo_1d(ndan): 断面流量/cross-sectional flow rate of a river(t=n)
!   subq_all(ndan): 断面横流入量/cross-sectional sub flow rate of a river
!--------------------------------------------------------------------------

module d1riv_cal_sub
    use d1riv_globals_mod

contains
!=====================================
! method of fractional step
!=====================================
!-------------------------------------
! equation of motion
!-------------------------------------
subroutine fractional_step_robust
    implicit none
    integer :: n, i
    real(8) :: term_adv, term_pres, term_fric_coef
    real(8) :: tmp_q, current_h, current_a, current_r, current_rn
    real(8) :: dh_dx

    !$omp parallel do default(shared) &
    !$omp private(n, i, term_adv, term_pres, term_fric_coef, tmp_q) &
    !$omp private(current_h, current_a, current_r, current_rn, dh_dx)
    do n = 2, ndan
        if((ntype_1d(n)==1) .or. (ntype_1d(n)==-2) .or. &
            (ntype_1d(n)==3) .or. (ntype_1d(n)==100)) cycle 
        ! --- dry ---
        if (a_1d(n) < tha_1d .or. h_1d(n) - rbed_1d(n) <= thd_1d) then
            q_1d(n) = 0.0d0
            vv_1d(n) = 0.0d0
            cycle
        endif

        ! set local
        current_a  = a_1d(n)
        current_h  = h_1d(n)
        current_r  = r_1d(n)
        current_rn = rn_1d(n)

        ! --- Step 1: 移流項/Advection ---
        term_adv = 0.0d0
        if(vv_1d(n) > 0.0d0) then
            term_adv = (vv_1d(n) * (q_1d(n) - q_1d(n-1)) + &
                        q_1d(n) * (vv_1d(n) - vv_1d(n-1))    ) / dx_1d(n)
        elseif(vv_1d(n) < 0.0d0) then
            if(ntype_1d(n)/=-1) then
                term_adv = (vv_1d(n) * (q_1d(n+1) - q_1d(n)) + &
                            q_1d(n) * (vv_1d(n+1) - vv_1d(n))   ) / dx_1d(n+1)
            else
                ! upstream
                do i = 1, nbound_1d
                    if(n == b_idx_1d(i)) then
                        term_adv = (vv_1d(n) * (up_q_1d(i) - q_1d(n))) / dx_1d(n)
                        exit
                    endif
                enddo
            endif
        endif

        ! q*
        tmp_q = q_1d(n) - term_adv * rivdt

        ! --- Step 2: 圧力項/Pressure & 摩擦項/Friction ---
        term_pres = 0.0d0
        term_fric_coef = 0.0d0
        if((current_a >= tha_1d) .and. (current_r > 1.0d-6)) then
            
            ! 1. Pressure Explicit
            dh_dx = (h_1d(n) - h_1d(n-1)) / dx_1d(n)
            term_pres = rivgg * current_a * dh_dx

            ! 2. Friction Semi-Implicit           
            ! K
            if (abs(vv_1d(n)) > 1.0d-10) then
                 term_fric_coef = rivgg * (current_rn**2) * abs(vv_1d(n)) / (current_r**(4.0d0/3.0d0))
            else
                 term_fric_coef = 0.0d0
            endif

            q_1d(n) = (tmp_q - term_pres * rivdt) / (1.0d0 + term_fric_coef * rivdt)

        else
            q_1d(n) = tmp_q
        endif

    end do
    !$omp end parallel do

end subroutine fractional_step_robust

!----------
! sub-flow
!----------
!--------------------------------------------------------------------------
! フラクショナルステップ法 横流入/method of fractional step sub flow
!   tributaryq_1d(ndan): 分流量/cross-sectional tributary q
!   breakq_1d(ndan): 破堤流量/cross-sectional break q
!   weirq_1d(ndan): 越流量/cross-sectional weir q(left+right)
!   weirq_l_1d(ndan): 左岸越流量/cross-sectional left weir q
!   weirq_r_1d(ndan): 右岸越流量/cross-sectional right weir q
!   pumpq_1d(ndan): ポンプ流量/cross-sectional pump q
!   sluiceq_1d(ndan): 樋門流量/cross-sectional sluice q
!--------------------------------------------------------------------------
subroutine sub_flow_1d
    implicit none

    integer :: n, i
    integer :: ups_d1, cnt_d1, sub_d1, main_d1, dns_d1  ! cross section id
    real(8) :: term_adv, term_pres, term_fric_coef
    real(8) :: tmp_q
    real(8) :: coef

    ! sum sub qin
    !$omp parallel do default(shared), private(n)
    do n = 1, ndan
        subq_all(n) = subq_all(n) + &
                      (tributaryq_1d(n) + breakq_1d(n) + &
                       pumpq_1d(n) + sluiceq_1d(n)) / rivdt
    enddo
    !$omp end parallel do

    ! Reset sub qin ---
    tributaryq_1d = 0.0d0
    breakq_1d = 0.0d0
    weirq_1d = 0.0d0
    weirq_l_1d = 0.0d0
    weirq_r_1d = 0.0d0
    pumpq_1d = 0.0d0
    sluiceq_1d = 0.0d0
    ! -----------------

    ! initialize ---
    ups_d1 = 0
    cnt_d1 = 0
    dns_d1 = 0
    main_d1 = 0
    sub_d1 = 0
    tmp_q = 0.0d0
    ! --------------

    if(nsubflow_1d > 0) then
        do i = 1, nsubflow_1d
            ! initialize ---
            term_adv = 0.0d0
            term_pres = 0.0d0
            term_fric_coef = 0.0d0
            ! --------------

            ! set cross section id ---
            if(ntype_1d(subflow_input(1,i))==-2 .or. ntype_1d(subflow_input(1,i))==-4) then
                ! tributary
                if(ntype_1d(subflow_input(1,i))==-4 .and. &
                   h_1d(subflow_input(1,i)) < h_1d(subflow_input(2,i))) then
                    q_1d(subflow_input(1,i)) = 0.0d0
                    cycle
                endif
                ups_d1 = subflow_input(1,i) + 1  ! sub upper (upstream)
                cnt_d1 = subflow_input(1,i)  ! sub (center)
                sub_d1 = subflow_input(1,i)  ! sub
                main_d1 = subflow_input(2,i)  ! main
                dns_d1 = subflow_input(2,i)  ! main (downstream)
                coef = -1.0d0
            elseif(ntype_1d(subflow_input(1,i))==-3) then
                ! distributary
                ups_d1 = subflow_input(1,i) + 1  ! main upper (upstream)
                cnt_d1 = subflow_input(1,i)  ! main (center)
                main_d1 = subflow_input(1,i)  ! main
                sub_d1 = subflow_input(2,i)  ! sub
                dns_d1 = subflow_input(2,i) - 1  ! sub (downstream)
                coef = 1.0d0
            else
                cycle
            endif
            ! ------------------------

            ! Step 1 ---
            if(a_1d(cnt_d1) >= tha2_1d) then
                if(vv_1d(sub_d1) > 0.0d0) then
                    term_adv = (vv_1d(sub_d1) * (q_1d(sub_d1) - q_1d(dns_d1)) + &
                            q_1d(sub_d1) * (vv_1d(sub_d1) - vv_1d(dns_d1))) / dx_1d(sub_d1)
                elseif(vv_1d(sub_d1) < 0.0d0) then
                    term_adv = (vv_1d(sub_d1) * (q_1d(ups_d1) - q_1d(sub_d1)) + &
                            q_1d(sub_d1) * (vv_1d(ups_d1) - vv_1d(sub_d1))) / dx_1d(ups_d1)
                endif
            endif
        
            tmp_q = q_1d(sub_d1) - term_adv * rivdt
            ! ----------
        
            ! Step 2 ---
            if((a_1d(cnt_d1) >= tha2_1d) .and. (r_1d(cnt_d1) > 1.0d-6)) then
                term_pres = rivgg * a_1d(cnt_d1) * (h_1d(cnt_d1) - h_1d(dns_d1)) / dx_1d(sub_d1)
                if(abs(vv_1d(sub_d1)) > 1.0d-10) then
                    term_fric_coef = rivgg * (rn_1d(cnt_d1)**2) * &
                                        abs(vv_1d(sub_d1)) / (r_1d(cnt_d1)**(4.d0/3.d0))
                endif
                q_1d(sub_d1) = (tmp_q - term_pres * rivdt) / (1.0d0 + term_fric_coef * rivdt)
            else
                q_1d(sub_d1) = tmp_q
            endif
            subq_all(main_d1) = subq_all(main_d1) + q_1d(sub_d1) * coef
            tributaryq_1d(main_d1) = tributaryq_1d(main_d1) + q_1d(sub_d1) * coef
            ! ----------
        enddo
    endif

end subroutine sub_flow_1d

! ------------------------
!  equation of continuity 
! ------------------------
subroutine continuous_1d
    implicit none

    integer n, i
    real(8) :: tmp_a

    !$omp parallel do default(shared), private(n, i, tmp_a)
    do n = 1, ndan
        tmp_a = 0.0d0
        if(ntype_1d(n)/=-1) then
            ! else upstream
            tmp_a = a_1d(n) - (q_1d(n+1) - q_1d(n) - subq_all(n)) * rivdt/((dx_1d(n)+dx_1d(n+1))*0.5d0)
        elseif(ntype_1d(n)==100) then
            ! downstream(connect 1d-2d)
            if(a_1d(n) > tha_1d) then
                vv_1d(n) = q_1d(n) / a_1d(n)
            else
                vv_1d(n) = 0.0d0
            endif
            cycle
        elseif(ntype_1d(n)==1) then
            ! downstream(nomal)
            if(a_1d(n) > tha_1d) then
                vv_1d(n) = q_1d(n) / a_1d(n)
            else
                vv_1d(n) = 0.0d0
            endif
            cycle
        else
            ! upstream (boundary-q)
            do i = 1, nbound_1d
                if(n == b_idx_1d(i)) then
                    tmp_a = a_1d(n) - (up_q_1d(i) - q_1d(n) - subq_all(n)) * rivdt/dx_1d(n)
                    exit
                endif
            enddo
        endif

        if(tmp_a >= tha_1d) then
            a_1d(n) = tmp_a
        else
            a_1d(n) = tha_1d
        endif

        call interp_a_1d(n)        

        if(a_1d(n) > tha_1d) then
            vv_1d(n) = q_1d(n) / a_1d(n)
        else
            vv_1d(n) = 0.0d0
        endif

    enddo
    !$omp end parallel do

    subq_all = 0.0d0  ! reset

end subroutine continuous_1d

!-------------------------------------
! interpolalate
!-------------------------------------
! --------------
!  a_1d -> else
! --------------
subroutine interp_a_1d(n)
    implicit none

    integer, intent(in) :: n
    integer :: i
    real(8) :: coef
    
    if(a_1d(n) >= a_table(n_tbl(n),n)) then
        coef = (a_1d(n) - a_table(n_tbl(n)-1,n)) / &
                (a_table(n_tbl(n),n) - a_table(n_tbl(n)-1,n))
        h_1d(n) = h_table(n_tbl(n)-1,n) + &
                  (h_table(n_tbl(n),n) - h_table(n_tbl(n)-1,n)) * coef
        r_1d(n) = r_table(n_tbl(n)-1,n) + &
                  (r_table(n_tbl(n),n) - r_table(n_tbl(n)-1,n)) * coef
        rn_1d(n) = h_table(n_tbl(n)-1,n) + &
                  (rn_table(n_tbl(n),n) - rn_table(n_tbl(n)-1,n)) * coef
        b_1d(n) = b_table(n_tbl(n)-1,n) + &
                  (b_table(n_tbl(n),n) - b_table(n_tbl(n)-1,n)) * coef
    elseif(a_1d(n) <= a_table(1,n)) then
        write(*,*)  ' UNST2D - WARNING : 1D-area csid ', unsttime_r, kp_1d(n), a_1d(n)
        h_1d(n) = rbed_1d(n) + thd_1d
        r_1d(n) = 0.0d0
        rn_1d(n) = 0.0d0
    else
        do i = 2, n_tbl(n)
            if(a_1d(n) < a_table(i,n)) then
                coef = (a_1d(n) - a_table(i-1,n)) / &
                        (a_table(i,n) - a_table(i-1,n))
                h_1d(n) = h_table(i-1,n) + &
                            (h_table(i,n) - h_table(i-1,n)) * coef
                r_1d(n) = r_table(i-1,n) + &
                            (r_table(i,n) - r_table(i-1,n)) * coef
                rn_1d(n) = rn_table(i-1,n) + &
                            (rn_table(i,n) - rn_table(i-1,n)) * coef
                b_1d(n) = b_table(i-1,n) + &
                            (b_table(i,n) - b_table(i-1,n)) * coef
                exit
            endif
        enddo
    endif

end subroutine interp_a_1d

! --------------
!  h_1d -> else
! --------------
subroutine interp_h_1d(n)
    implicit none

    integer, intent(in) :: n
    integer :: i
    real(8) :: coef
    
    if(h_1d(n) >= h_table(n_tbl(n),n)) then
        coef = (h_1d(n) - h_table(n_tbl(n)-1,n)) / &
                (h_table(n_tbl(n),n) - h_table(n_tbl(n)-1,n))
        a_1d(n) = a_table(n_tbl(n)-1,n) + &
                  (a_table(n_tbl(n),n) - a_table(n_tbl(n)-1,n)) * coef
        r_1d(n) = r_table(n_tbl(n)-1,n) + &
                  (r_table(n_tbl(n),n) - r_table(n_tbl(n)-1,n)) * coef
        rn_1d(n) = h_table(n_tbl(n)-1,n) + &
                  (rn_table(n_tbl(n),n) - rn_table(n_tbl(n)-1,n)) * coef
        b_1d(n) = b_table(n_tbl(n)-1,n) + &
                  (b_table(n_tbl(n),n) - b_table(n_tbl(n)-1,n)) * coef
    elseif(h_1d(n) < h_table(1,n)+thd_1d) then
        write(*,*) ' UNST2D - WARNING : 1D-depth csid ', unsttime_r, kp_1d(n), h_1d(n)
        a_1d(n) = 0.0d0
        r_1d(n) = 0.0d0
        rn_1d(n) = 0.0d0
    else
        do i = 2, n_tbl(n)
            if(h_1d(n) < h_table(i,n)) then
                coef = (h_1d(n) - h_table(i-1,n)) / &
                        (h_table(i,n) - h_table(i-1,n))
                a_1d(n) = a_table(i-1,n) + &
                            (a_table(i,n) - a_table(i-1,n)) * coef
                r_1d(n) = r_table(i-1,n) + &
                            (r_table(i,n) - r_table(i-1,n)) * coef
                rn_1d(n) = rn_table(i-1,n) + &
                            (rn_table(i,n) - rn_table(i-1,n)) * coef
                b_1d(n) = b_table(i-1,n) + &
                            (b_table(i,n) - b_table(i-1,n)) * coef
                exit
            endif
        enddo
    endif
        
end subroutine interp_h_1d

!----------------------------------------------
! Culculate inflow or outflow on the boundares
!----------------------------------------------
! downstream
subroutine h_bound_1d
    implicit none

    integer :: n, i
    real(8) :: tmp_h

    do n = 1, nbound_1d
        if(ntype_1d(b_idx_1d(n))/=1 .or. b_dt_1d(n)<=0.0d0 .or. b_dome_1d(n)==0) cycle
        i = int(unsttime_r / b_dt_1d(n)) + 1
        tmp_h = b_data_1d(i,n) + (unsttime_r - b_dt_1d(n) * dble(i - 1)) / &
                b_dt_1d(n) * (b_data_1d(i+1,n) - b_data_1d(i,n))
        h_1d(b_idx_1d(n)) = max(tmp_h, rbed_1d(b_idx_1d(n))+thd_1d)
        call interp_h_1d(b_idx_1d(n))
    enddo
    
end subroutine h_bound_1d

subroutine q_bound_1d
    implicit none

    integer :: n, i
    real(8) :: tmp_q

    do n = 1, nbound_1d
        if(ntype_1d(b_idx_1d(n))/=-1 .or. b_dt_1d(n)<=0.0d0 .or. ntype_1d(n)==100) cycle
        i = int(unsttime_r / b_dt_1d(n)) + 1
        tmp_q = -1.d0 * (b_data_1d(i,n) + (unsttime_r - b_dt_1d(n) * dble(i - 1)) / &
                b_dt_1d(n) * (b_data_1d(i+1,n) - b_data_1d(i,n)))
        up_q_1d(n) = min(tmp_q, -0.1d0)
    enddo
    
end subroutine q_bound_1d

end module d1riv_cal_sub

module unst_cnct_1d2d
    use unst_globals_mod
    use d1riv_globals_mod

contains
!=====================================
! connect 1d - 2d
!=====================================
! ---------------------
!  honma weir equation added by d.baba, supported by CJS
! ---------------------
subroutine honma_weir(h1,h2,rivb,q0)
    implicit none
    real(8), intent(in) :: h1
    real(8), intent(in) :: h2
    real(8), intent(in) :: rivb
    real(8), intent(out) :: q0
    real(8) arg

    arg = 0.0d0
    q0 = 0.0d0
    if (h1 == h2) then
        q0 = 0.0d0
    elseif(h2*3.0d0 < h1*2.0d0) then
        ! kanzen
        q0 = 0.35d0 * rivb * h1 * sqrt(2.0d0*rivgg*h1)
    else
        arg = h1 - h2
        ! moguri
        if(arg>0.0d0) then
            q0 = 0.91d0 * rivb * h2 * sqrt(2.0d0*rivgg*(h1-h2))
        else
            q0 = 0.0d0
        endif
    endif

end subroutine honma_weir

! ------------------
!  weir calculation
! ------------------
subroutine weir_equation
    implicit none
    real(8) inlandh, tmp_h1, tmp_h2, h1, h2, wir_inout, q0
    real(8) :: total_out_vol, available_vol, limit_ratio
    integer :: n, i

    ! tmp_h1: 2d excess water level
    ! tmp_h2: 1d excess water level

    weirq_1d = 0.0d0
    weirq_l_1d = 0.0d0
    weirq_r_1d = 0.0d0
    !$omp parallel do default(shared) &
    !$omp private(i, n, inlandh, tmp_h1, tmp_h2, h1, h2) &
    !$omp private(wir_inout, q0, total_out_vol, available_vol, limit_ratio)
    do n = 1, ndan
        ! left levee
        if(ncnct_1d2d(1,n)/=0 .and. w_alpha_1d(1,n) > 0.0d0 .and. bktype_1d(1,n)/=1) then
            inlandh = 0.0d0
            q0 = 0.0d0
            do i = 1, ncnct_1d2d(1,n)
                inlandh = inlandh + unsth(cnct_1d2d_idx(1,i,n)) + baseo(cnct_1d2d_idx(1,i,n))
            enddo
            inlandh = inlandh / ncnct_1d2d(1,n)
            if (crown_1d(1,n) < h_1d(n) .or. crown_1d(1,n) < inlandh) then
                tmp_h1 = max(inlandh - crown_1d(1,n), 0.0d0)
                tmp_h2 = max(h_1d(n) - crown_1d(1,n), 0.0d0)
                if (tmp_h1>tmp_h2) then
                    ! 2d -> 1d (weirq_l_1d>0)
                    h1 = tmp_h1
                    h2 = tmp_h2
                    wir_inout = 1.0d0
                else
                    ! 2d <- 1d (weirq_l_1d<0)
                    h1 = tmp_h2
                    h2 = tmp_h1
                    wir_inout = -1.0d0
                endif
                if ((inlandh - inland_1d(1,n)) > th .or. (h_1d(n) - rbed_1d(n)) > thd_1d ) then
                    call honma_weir(h1, h2, weir_dx_1d(n), q0)
                    if (q0 > 0.0d0) weirq_l_1d(n) = wir_inout * w_alpha_1d(1,n) * q0 * w_angle_1d(n)  ! * rivdt
                else
                    weirq_l_1d(n) = 0.0d0
                endif
            else
                weirq_l_1d(n) = 0.0d0
            endif
        endif

        ! right levee
        if(ncnct_1d2d(2,n)/=0 .and. w_alpha_1d(2,n) > 0.0d0 .and. bktype_1d(2,n)/=1) then
            inlandh = 0.0d0
            q0 = 0.0d0
            do i = 1, ncnct_1d2d(2,n)
                inlandh = inlandh + unsth(cnct_1d2d_idx(2,i,n)) + baseo(cnct_1d2d_idx(2,i,n))
            enddo
            inlandh = inlandh / ncnct_1d2d(2,n)
            if (crown_1d(2,n) < h_1d(n) .or. crown_1d(2,n) < inlandh) then
                tmp_h1 = max(inlandh - crown_1d(2,n), 0.0d0)
                tmp_h2 = max(h_1d(n) - crown_1d(2,n), 0.0d0)
                if (tmp_h1>tmp_h2) then
                    ! 2d -> 1d (weirq_r_1d>0)
                    h1 = tmp_h1
                    h2 = tmp_h2
                    wir_inout = 1.0d0
                else
                    ! 2d <- 1d (weirq_r_1d<0)
                    h1 = tmp_h2
                    h2 = tmp_h1
                    wir_inout = -1.0d0
                endif
                if ((inlandh - inland_1d(2,n)) > th .or. (h_1d(n) - rbed_1d(n)) > thd_1d ) then
                    call honma_weir(h1, h2, weir_dx_1d(n), q0)
                    if (q0 > 0.0d0) weirq_r_1d(n) = wir_inout * w_alpha_1d(2,n) * q0 * w_angle_1d(n)  ! * rivdt
                else
                    weirq_r_1d(n) = 0.0d0
                endif
            else
                weirq_r_1d(n) = 0.0d0
            endif
        endif

        weirq_1d(n) = weirq_l_1d(n) + weirq_r_1d(n)
        
        ! Volume check
        if (weirq_1d(n) < 0.0d0) then
            total_out_vol = abs(weirq_1d(n)) * rivdt
            available_vol = max(a_1d(n), 0.0d0) * dx_1d(n) * 0.8d0
            if (total_out_vol > available_vol) then
                if (total_out_vol > 1.0d-10) then ! ゼロ割防止
                    limit_ratio = available_vol / total_out_vol
                    weirq_l_1d(n) = weirq_l_1d(n) * limit_ratio
                    weirq_r_1d(n) = weirq_r_1d(n) * limit_ratio
                    weirq_1d(n)   = weirq_1d(n)   * limit_ratio
                else
                    weirq_l_1d(n) = 0.0d0
                    weirq_r_1d(n) = 0.0d0
                    weirq_1d(n)   = 0.0d0
                endif
            endif
        endif
        if (ncnct_1d2d(1,n)/=0 .and. cnct_1d2d_idx(1,1,n)/=0) then
            riv_eq(cnct_1d2d_idx(1,1,n)) = riv_eq(cnct_1d2d_idx(1,1,n)) - weirq_l_1d(n)
        endif
        if ((ncnct_1d2d(2,n)==0 .or. cnct_1d2d_idx(2,1,n)==0)) then
            riv_eq(cnct_1d2d_idx(2,1,n)) = riv_eq(cnct_1d2d_idx(2,1,n)) - weirq_r_1d(n)
        endif
    enddo
    !$omp end parallel do
    
end subroutine weir_equation

! ----------------------------------------------
!  inflow 1D -> 2D
! ----------------------------------------------
subroutine calc_2d_to_1d_inflow
    implicit none
    integer :: n, ii, me_2d
    real(8) :: h_1d_v, h_2d_v, h_diff, dist_x, flux_q
    real(8) :: manning_coef, area_eff, radius_eff
    
    do n = 1, nbound_1d
        ii = b_idx_1d(n)
        if(ntype_1d(ii)/=-1) cycle  ! 上流端のみ
        
        me_2d = b_upme_1d(n)
        if(me_2d == 0) cycle
        
        h_1d_v = h_1d(ii)
        h_2d_v = unsth(me_2d) + baseo(me_2d)
        
        flux_q = 0.0d0
        
        ! 2D(unsth+baseo) > 1D(h_1d)
        if (h_2d_v > h_1d_v + th) then  ! safety(th)
            h_diff = h_2d_v - h_1d_v
            dist_x = dx_1d(ii) * 0.5d0
            
            if(a_1d(ii) > tha_1d) then
                area_eff = a_1d(ii)
                radius_eff = r_1d(ii)
                manning_coef = rn_1d(ii)
            else
                ! dry
                area_eff = a_table(2, ii)
                radius_eff = r_table(2, ii)
                manning_coef = rn_table(2, ii)
            endif
            
            ! Manning
            if(manning_coef > 0.0d0 .and. radius_eff > 0.0d0) then
                flux_q = area_eff * (1.0d0/manning_coef) * &
                         (radius_eff**(2.0d0/3.0d0)) * sqrt(h_diff/dist_x)
                
                ! 発散する場合は有効, 原則コメントアウト
                ! limit(default 70%)
                !if(unsth(me_2d) > th) then
                !    flux_q = min(flux_q, 0.8d0*unsth(me_2d)*smesh(me_2d)/rivdt)
                !endif
            endif
        endif
        
        ! set
        up_q_1d(n) = -flux_q
        
        riv_eq(me_2d) = riv_eq(me_2d) - flux_q
        
    enddo
    
end subroutine calc_2d_to_1d_inflow

! --------------------
!  discharge 1d -> 2d
! --------------------
subroutine calc_1d_to_2d_outflow
    implicit none
    integer :: n, ii, me_2d
    real(8) :: h_1d_v, h_2d_v, h_diff, dist_x
    real(8) :: area_eff, radius_eff, flux_q
    
    do n = 1, nbound_1d
        ii = b_idx_1d(n)
        if(ntype_1d(ii)/=100) cycle
        
        me_2d = b_dome_1d(n)
        if(me_2d == 0) cycle
        
        h_1d_v = h_1d(ii)
        h_2d_v = unsth(me_2d) + baseo(me_2d)
        
        flux_q = 0.0d0
        
        ! 1D(h_1d) > 2D(unsth+baseo)
        if (h_1d_v > h_2d_v) then
            h_diff = h_1d_v - h_2d_v
            dist_x = dx_1d(ii)
            if(dist_x<=0.0d0) dist_x = dx_1d(ii+1)
            
            if(h_1d_v - rbed_1d(ii) > thd_1d .and. a_1d(ii) > tha_1d) then
                ! Manning
                flux_q = -a_1d(ii) * (1.0d0/rn_1d(ii)) * &
                         (r_1d(ii)**(2.0d0/3.0d0)) * sqrt(h_diff/dist_x)
                
                ! 発散する場合は有効, 原則コメントアウト
                ! discharge limit(default 70%)
                !flux_q = min(flux_q, 0.7d0*a_1d(ii)*dx_1d(ii)/rivdt)
            endif
            
        ! 2D(h_1d) > 1D(unsth+baseo)
        else if (h_2d_v > h_1d_v) then
            h_diff = h_2d_v - h_1d_v
            dist_x = dx_1d(ii)
            if(dist_x<=0.0d0) dist_x = dx_1d(ii+1)
            
            if(unsth(me_2d) > th) then
                area_eff = min(a_1d(ii), b_1d(ii)*unsth(me_2d))
                if(area_eff < tha_1d) area_eff = tha_1d
                
                radius_eff = area_eff / (b_1d(ii) + 2.0d0*unsth(me_2d))
                
                flux_q = area_eff * (1.0d0/rn_1d(ii)) * &
                         (radius_eff**(2.0d0/3.0d0)) * sqrt(h_diff/dist_x)
                
                ! 発散する場合は有効, 原則コメントアウト
                ! discharge limit(default 70%)
                !flux_q = max(flux_q, -0.7d0*unsth(me_2d)*smesh(me_2d)/rivdt)
            endif
        endif
        
        ! set boundary data
        b_data_1d(1, n) = h_2d_v
        b_data_1d(2, n) = h_2d_v
         
        riv_eq(me_2d) = riv_eq(me_2d) - flux_q
        
    enddo
    
end subroutine calc_1d_to_2d_outflow

end module unst_cnct_1d2d

module unst_1d_main
    use d1riv_globals_mod
    use d1riv_cal_sub
    use unst_cnct_1d2d
    
contains    

!=====================================
! Main
!=====================================
subroutine d1riv_main(unsttime)
    implicit none
    real(8), intent(in) :: unsttime

    unsttime_r = unsttime

    call h_bound_1d
    call fractional_step_robust
    call sub_flow_1d
    call q_bound_1d
    call continuous_1d
    call weir_equation

end subroutine d1riv_main

subroutine d1riv_spinup
    implicit none
    integer i

    do i = 1, ndan
        h_1d(i) = ho_1d(i)
        q_1d(i) = qo_1d(i)
        call interp_h_1d(i)
    enddo
    if (d1_spin_upn>0) then
        unsttime_r = 0.0d0
        write(*,*) 'Sipn up 1d river -- ', 0, '/', d1_spin_upn
        do i = 1, d1_spin_upn
            call h_bound_1d
            call fractional_step_robust
            call sub_flow_1d
            call q_bound_1d
            call continuous_1d

            if(mod(i,spout)==0) then
                write(*,*) 'Sipn up 1d river -- ', i, '/', d1_spin_upn
            endif
        enddo
        write(*,*) '  Sipn up 1d river fin -- ', d1_spin_upn, '/', d1_spin_upn
    endif

end subroutine d1riv_spinup

end module unst_1d_main
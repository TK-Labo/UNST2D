! Read 1D River  ! release 25/12 v.1.0.5
!--------------------------------------------------------------------------------
! UNST2D一次元河道(水路)読み込みサブルーチン/UNST2D 1d river read subroutine
! 変数/Variable:
! > グローバル変数/global variable
!   ndan: 断面総数/num of cross-section
!   ntype_1d(ndan): 断面タイプ/cross-sectional type
!                   =0:通常/nomal,
!                   =1:下流端/downstream, =-1:上流端/upstream
!                   =2:合流本川/tributary-main, =-2:合流支川/tributary-sub
!                   =3:分流支川/distributary-sub, =-3:分流本川/distributary-main
!                   =4:合流(逆流防止樋門あり)本川/tributary(stop back water)-main
!                   =-4:合流(逆流防止樋門あり)支川/tributary(stop back water)-sub
!                   =100:下流端(2次元との接続)/downstream(connect-2D)
! > ローカル変数/local variable
!   
!--------------------------------------------------------------------------------

module unst_d1riv_read
    use d1riv_globals_mod

contains

subroutine d1rivdat(timmax, dtq, mesh, baseo, frivcntl)
    implicit none
    real(8), intent(in) :: dtq, timmax
    integer, intent(in) :: mesh
    real(8), intent(in) :: baseo(mesh)
    character(len=50), intent(in) :: frivcntl
    integer idummy, j, i, ios, n, nb_data, k, max_ncnct_1d2d
    integer oudan_format
    integer read_count, read_count2, tmp_ndan
    integer xznum
    real(8) rdummy, rpre_d, rpre_h, rpre_m, rpre_s
    real(8) sum_inland_bs_l, sum_inland_bs_r
    character(len=100) line
    integer, allocatable :: ud_idx_1d(:,:)  ! up-down bound cross section id
    integer, allocatable :: subflow_type(:)  ! type
    integer frivcntl_unit, &
            fjudan_unit, foudan_unit, fbound_unit, f1d2d_unit
    character(len=100) fjudan, foudan, fbound, f1d2d, foudan_table

    !=======================
    ! Step1 : Read cntl.dat
    !=======================
    ! input param
    open(newunit=frivcntl_unit, file=frivcntl, action = 'read')
    read(frivcntl_unit, '(7x, 4f6.1)') rpre_d, rpre_h, rpre_m, rpre_s
    read(frivcntl_unit, '(7x, f11.2)') rivdt  ! common
    read(frivcntl_unit, '(7x, i11)') nriv           ! common
    read(frivcntl_unit, *)
    !---

    ! input file path
    read(frivcntl_unit, '(A)') fjudan  ! common()
    read(frivcntl_unit, '(A)') fbound  ! common(boundary condition)
    read(frivcntl_unit, '(A)') finit   ! common(initialize)
    read(frivcntl_unit, '(A)') f1d2d   ! common(connect 1d -2d)
    read(frivcntl_unit, *)
    read(frivcntl_unit, *)
    !! oudan
    read(frivcntl_unit, '(7x, i11)') oudan_format
    read(frivcntl_unit, '(7x, f11.2)') dz  ! common(table dz)
    read(frivcntl_unit, '(A)') foudan  ! common()
    read(frivcntl_unit, '(A)') foudan_table
    read(frivcntl_unit, *)
    read(frivcntl_unit, *)

    read(frivcntl_unit, '(7x, i11)') nsubflow_1d  ! common()
    if(nsubflow_1d > 0) then
        allocate(subflow_type(nsubflow_1d))
        allocate(subflow_input(2,nsubflow_1d))
        do n = 1, nsubflow_1d
            read(frivcntl_unit, *, iostat=ios) subflow_type(n), &
                        subflow_input(1,n), subflow_input(2,n)
            if(ios/=0) read(frivcntl_unit, *, iostat=ios) subflow_type(n), &
                                                            subflow_input(1,n)
        enddo
    endif
    read(frivcntl_unit, *)
    read(frivcntl_unit, *)

    read(frivcntl_unit, '(A)') fd1out   ! output
    read(frivcntl_unit, '(A)') fd1mx   ! output

    write(*,*)
    write(*,*) '== 1D river Parameters =='
    write(*,'(a19,f10.2,7a,f10.2,7a,f10.2,7a,f10.2,7a)') ' - 1d spinup time:', rpre_d, ' day /', &
                                    rpre_h, ' hour /', rpre_m, ' min /' , rpre_s, ' sec '
    write(*,'(a18,i5)') ' - num of river :', nriv
    write(*,'(a18,i5)') ' - num of joint :', nsubflow_1d
    write(*,*)
    write(*,*) '* 1d River Data *'
    write(*,'(a16, a)') ' - judan      :', fjudan
    write(*,'(a16, a)') ' - bound      :', fbound
    write(*,'(a16, a)') ' - init       :', finit
    write(*,'(a16, a)') ' - cntn_1d2d  :', f1d2d
    write(*,'(a16, a)') ' - oudan      :', foudan
    if(oudan_format==1) write(*,'(a18,f10.2, a3)') '   - dz pitch   :', dz, '(m)'
    if(oudan_format==1) write(*,'(a18, a)') '   - oudan      :', foudan_table
    write(*,*)
    write(*,*) '* UNST Output *'
    write(*,'(a16, a)') ' - d1out      :', fd1out
    write(*,'(a16, a)') ' - d1outmax   :', fd1mx
    write(*,*)

    allocate(ud_idx_1d(2,nriv))
    allocate(riv_ndan(nriv))

    close(frivcntl_unit)

    !============================
    ! Step2 : Set constant param
    !============================
    d1_spin_ups = rpre_d*86400 +  rpre_h*3600 + rpre_m*60 + rpre_s

    read_count = 0
    read_count2 = 1
    ndan = 0

    !==========================
    ! Step3 : Read input files
    !==========================
    !------------------
    ! open input files
    !------------------
    open(newunit= fjudan_unit, file = fjudan, action = 'read')
    open(newunit= fbound_unit, file = fbound, action = 'read')
    open(newunit= f1d2d_unit, file = f1d2d, action = 'read')
    if(oudan_format==0) then
        open(newunit= foudan_unit, file = foudan, action = 'read')
    else
        open(newunit= foudan_unit, file = foudan_table, action = 'read')
    endif

    !---------------
    ! 01 judan data
    !---------------
    ! -- read --
    tmp_ndan = 0
    do i = 1, nriv
        ud_idx_1d(1,i) = tmp_ndan + 1    ! downstream
        read(fjudan_unit, *, iostat=ios) riv_ndan(i)
        tmp_ndan = tmp_ndan + riv_ndan(i)
        ud_idx_1d(2,i) = tmp_ndan  ! upstream
    enddo
    close(fjudan_unit)
    ndan = tmp_ndan
    write(*,*) 'num of cross section - ', ndan

    ! -- allocate --
    allocate(ntype_1d(ndan))
    allocate(h_1d(ndan), vv_1d(ndan), q_1d(ndan))
    allocate(h_1dmax(ndan), vv_1dmax(ndan))
    allocate(kp_1d(ndan), dx_1d(ndan))
    allocate(a_1d(ndan), r_1d(ndan), b_1d(ndan), rn_1d(ndan))
    allocate(depth_1d(ndan))
    allocate(rbed_1d(ndan), rzmax_1d(ndan))
    allocate(dan_record(ndan))
    ! - - -
    allocate(n_tbl(ndan))
    ! - - -
    allocate(subq_all(ndan))
    allocate(tributaryq_1d(ndan))
    allocate(breakq_1d(ndan), breakq_l_1d(ndan), breakq_r_1d(ndan))
    allocate(weirq_1d(ndan), weirq_l_1d(ndan), weirq_r_1d(ndan))
    allocate(pumpq_1d(ndan), sluiceq_1d(ndan))
    ! - - -
    allocate(weir_dx_1d(ndan))
    allocate(crown_1d(2,ndan), inland_1d(2,ndan))
    allocate(w_alpha_1d(2,ndan), w_angle_1d(ndan))
    allocate(ncnct_1d2d(2,ndan))
    ! - - -
    allocate(bktype_1d(2,ndan))  ! development
    bktype_1d = 0
    ! ------------

    ! set upstream-downstream id
    ntype_1d = 0
    ntype_1d(ud_idx_1d(1,:)) = 1  ! default bound-h
    ntype_1d(ud_idx_1d(2,:)) = -1  ! default bound-q

    ! update subflow id
    if(nsubflow_1d > 0) then
        do n = 1, nsubflow_1d
            select case(subflow_type(n))
            case (100)  ! downstream -> 2d
                ntype_1d(subflow_input(1,n)) = 100
            case (1)  ! tributary (sub(-2) -> main(2))
                if(ntype_1d(subflow_input(2,n)) /= -1) then
                    ntype_1d(subflow_input(1,n)) = -2  ! sub river
                    ntype_1d(subflow_input(2,n)) = 2  ! main river
                else
                    write(*,*) ' UNST2D - WARNING : subflow id ', n, ' - subflow type is invaild value, not subflow applyed.'
                endif
            case (2)  ! distributary  (main(-3) -> sub(3))
                if(ntype_1d(subflow_input(1,n)) /= 1) then
                    ntype_1d(subflow_input(1,n)) = -3  ! main river
                    ntype_1d(subflow_input(2,n)) = 3  ! sub river
                else
                    write(*,*) ' UNST2D - WARNING : subflow id ', n, ' - subflow type is invaild value, not subflow applyed.'
                endif
            case (3)  ! tributary sluice(stop back water sub(-4) -> main(4))
                write(*,*) ' UNST2D - WARNING : subflow id ', n, ' - subflow type is development, Comming soon...'
                if(ntype_1d(subflow_input(2,n)) /= -1) then
                    ntype_1d(subflow_input(1,n)) = -4  ! sub river
                    ntype_1d(subflow_input(2,n)) = 4  ! main river
                else
                    write(*,*) ' UNST2D - WARNING : subflow id ', n, ' - subflow type is invaild value, not subflow applyed.'
                endif
            case (4)  ! pump
                write(*,*) ' UNST2D - WARNING : subflow id ', n, ' - subflow type is development.'
            case (-100)  ! RRI-upstream (only RRI-UNST2D)  add v.1.0.6
                ntype_1d(subflow_input(2,n)) = -100
            case default
                write(*,*) ' UNST2D - WARNING : subflow id ', n, ' - subflow type is invaild value, not subflow applyed.'
            end select
        enddo
        deallocate(subflow_type)
        deallocate(ud_idx_1d)
    endif

    !---------------
    ! 02 oudan data
    !---------------
    if(oudan_format==0) then
        read_count = 0
        xznum = 0
        do n = 1, ndan
            read(foudan_unit, '(A)', iostat=ios) line
            if(ios/=0) exit
            read_count = 0
            line = trim(adjustl(line))  ! line clean
            call read_oudan_header(line, kp_1d(n), dan_record(n))
            do i = 1, dan_record(n)
                read(foudan_unit,*)
            enddo
            xznum = max(xznum, dan_record(n))
        enddo
        rewind(foudan_unit)

        ! -- allocate --
        allocate(dan_x(xznum,ndan), dan_z(xznum,ndan), dan_rn(xznum,ndan))

        dan_rn = 0.030d0
        rbed_1d = 9999.d0
        rzmax_1d = -9999.0d0
        n_tbl = 0
        max_tbl = 0
        do n = 1, ndan
            read(foudan_unit,*)  ! skip header
            do i = 1, dan_record(n)
                read(foudan_unit, '(A)') line
                line = trim(adjustl(line))
                call read_oudan_data(line, dan_x(i,n), dan_z(i,n), dan_rn(i,n))
                rbed_1d(n) = min(rbed_1d(n), dan_z(i,n))
                rzmax_1d(n) = max(rzmax_1d(n), dan_z(i,n))
            enddo
            crown_1d(1,n) = dan_z(1,n)  ! left levee
            crown_1d(2,n) = dan_z(dan_record(n),n) ! right levee
            n_tbl(n) = int((rzmax_1d(n)-rbed_1d(n))/dz)+1
            max_tbl = max(max_tbl, n_tbl(n))
        enddo

        call d1riv_table(foudan_table)
    
    else
        read(foudan_unit, *) max_tbl
        allocate(h_table(max_tbl, ndan), a_table(max_tbl, ndan))
        allocate(r_table(max_tbl, ndan), b_table(max_tbl, ndan))
        allocate(rn_table(max_tbl, ndan))
        do n = 1, ndan
            read(foudan_unit, *) kp_1d(n), dx_1d(n), n_tbl(n), crown_1d(1,n), crown_1d(2,n), w_angle_1d(n)
            do i = 1, n_tbl(n)
                read(foudan_unit, *) h_table(i,n), a_table(i,n), r_table(i,n), b_table(i,n), rn_table(i,n)
            enddo
            rbed_1d(n) = h_table(1,n)
        enddo
    endif
    close(foudan_unit)

    !------------------
    ! 03 boundary data
    !------------------
    nbound_1d = 0
    max_nb = 0
    do n = 1, ndan
        select case(ntype_1d(n))
        case (-1, 1, 100)
            nbound_1d = nbound_1d + 1
        case default
            cycle
        end select
    enddo

    allocate(b_dt_1d(nbound_1d))

    b_dt_1d = timmax
    do n = 1, nbound_1d
        read(fbound_unit, '(A)', iostat=ios) line
        if(ios/=0) exit
        line = trim(adjustl(line))  ! line clean
        read(line, *) idummy, nb_data
        read(fbound_unit,*)
        max_nb = max(max_nb, nb_data+1)
    enddo
    rewind(fbound_unit)

    if(any(ntype_1d==100)) max_nb = max(int(timmax / dtq)+1, max_nb)

    allocate(b_idx_1d(nbound_1d))
    allocate(b_data_1d(max_nb, nbound_1d))
    allocate(up_q_1d(nbound_1d), up_h_1d(nbound_1d))
    allocate(b_upme_1d(nbound_1d), b_dome_1d(nbound_1d))

    b_idx_1d = 0
    b_data_1d = 0.0d0
    b_upme_1d = 0
    b_dome_1d = 0
    do n = 1, nbound_1d
        read(fbound_unit,*) b_idx_1d(n), nb_data
        select case(ntype_1d(b_idx_1d(n)))
        case(-1, 1, 100)
            if(nb_data>0) then
                b_dt_1d(n) = timmax / dble(nb_data)  ! set time
                read(fbound_unit,*) (b_data_1d(i,n), i=1, nb_data)
            else
                nb_data = 1
                b_dt_1d(n) = timmax
                if(ntype_1d(b_idx_1d(n)) == -1) then
                    read(fbound_unit, *) b_upme_1d(n)
                    b_data_1d(1,n) = 0.10d0
                elseif(ntype_1d(b_idx_1d(n))==1 .or. ntype_1d(b_idx_1d(n))==100) then
                    read(fbound_unit, *) b_dome_1d(n)
                    b_data_1d(1,n) = rbed_1d(b_idx_1d(n)) + thd_1d                    
                endif
            endif
            b_data_1d(nb_data+1,n) = b_data_1d(nb_data,n)
        case default  ! development
            write(*,*) ' UNST2D - WARNING : BOUNDARY DATA - id', n, ' - data num is invaild value'
        end select
    enddo
    close(fbound_unit)

    !----------------------
    ! connect 1d - 2d data
    !----------------------
    ncnct_1d2d = 0
    max_ncnct_1d2d = 0
    read(f1d2d_unit, *)
    do n = 1, ndan
        read(f1d2d_unit, '(A)') line
        line = trim(adjustl(line))
        read(line, *) idummy, w_alpha_1d(1,n), w_alpha_1d(2,n), &
                      ncnct_1d2d(1,n), ncnct_1d2d(2,n)
        max_ncnct_1d2d = max(max_ncnct_1d2d,ncnct_1d2d(1,n),ncnct_1d2d(2,n))
    enddo
    rewind(f1d2d_unit)

    allocate(cnct_1d2d_idx(2,max(max_ncnct_1d2d,1),ndan))
    cnct_1d2d_idx = 0
    
    if (max_ncnct_1d2d/=0) then
        read(f1d2d_unit, *)
        do n = 1, ndan
            if (ncnct_1d2d(1,n)==0 .and. ncnct_1d2d(2,n)==0) then
                read(f1d2d_unit, *)
            else
                read(f1d2d_unit, *) idummy, rdummy, rdummy, idummy, idummy, &
                                    (cnct_1d2d_idx(1,i,n), i=1, ncnct_1d2d(1,n)), &
                                    (cnct_1d2d_idx(2,j,n), j=1, ncnct_1d2d(1,n))
                sum_inland_bs_r = 0.0d0
                sum_inland_bs_l = 0.0d0
                do k = 1, ncnct_1d2d(1,n)
                    if(cnct_1d2d_idx(1,1,n)==0) then
                        w_alpha_1d(1,n) = 0.0d0
                        exit
                    endif
                    sum_inland_bs_l = sum_inland_bs_l + baseo(cnct_1d2d_idx(1,k,n))
                enddo
                do k = 1, ncnct_1d2d(2,n)
                    if(cnct_1d2d_idx(2,1,n)==0) then
                        w_alpha_1d(2,n) = 0.0d0
                        exit
                    endif
                    sum_inland_bs_r = sum_inland_bs_r + baseo(cnct_1d2d_idx(2,k,n))
                enddo
                if(cnct_1d2d_idx(1,1,n)/=0) inland_1d(1,n) = sum_inland_bs_l / ncnct_1d2d(1,n)
                if(cnct_1d2d_idx(2,1,n)/=0) inland_1d(2,n) = sum_inland_bs_r / ncnct_1d2d(2,n)

                if(crown_1d(1,n) < inland_1d(1,n) .and. w_alpha_1d(1,n) /= 0.0d0) then
                    write(*,*) ' UNST2D -  WARNING : corss section id', n, ' - inland elevation higher than river crown(left side).'
                    crown_1d(1,n) = inland_1d(1,n)
                endif
                if(crown_1d(2,n) < inland_1d(2,n) .and. w_alpha_1d(2,n) /= 0.0d0) then
                    write(*,*) ' UNST2D - WARNING : corss section id', n, ' - inland elevation higher than river crown(right side).'
                    crown_1d(2,n) = inland_1d(2,n)
                endif
            endif
        enddo
    else
        write(*,*) "No Connect 1d-2d"
    endif

    close(f1d2d_unit)

end subroutine d1rivdat

subroutine d1riv_table(foudan_table)
    implicit none

    character(len=100), intent(in) :: foudan_table

    integer n, i, j
    integer foudan_table_unit
    real(8) tmp_dx, tmp_dz, intercect_x  ! length
    real(8) tmp_a, tmp_s, tmp_rn_s  ! addition
    logical l_point_wd, r_point_wd  ! check

    real(8), allocatable :: s_table(:,:)  ! s

    real(8) tmp_angle

    ! allocate ---
    allocate(h_table(max_tbl, ndan), a_table(max_tbl, ndan))
    allocate(r_table(max_tbl, ndan), b_table(max_tbl, ndan))
    allocate(rn_table(max_tbl, ndan))
    allocate(s_table(max_tbl, ndan))
    ! ------------

    write(*,*) '    * Calculating 1D hydraulic tables ...'

    ! initialize
    h_table = 0.0d0
    a_table = 0.0d0
    r_table = 0.0d0
    b_table = 0.0d0
    rn_table = 0.030d0  ! deault val
    s_table = 0.0d0

    do n = 1, ndan
        do i = 1, n_tbl(n)
            h_table(i,n) = rbed_1d(n) + dble(i-1) * dz
            if(i==1) cycle
            tmp_a = 0.0d0
            tmp_s = 0.0d0
            tmp_rn_s = 0.0d0
            ! coord j - j+1 -------------------------------
            do j = 1, dan_record(n)-1
                ! comparison of water level and coords
                l_point_wd = dan_z(j,n) <= h_table(i,n)
                r_point_wd = dan_z(j+1,n) <= h_table(i,n)
                ! set distance
                tmp_dx = dan_x(j+1,n) - dan_x(j,n)
                tmp_dz = dan_z(j+1,n) - dan_z(j,n)
            
                if(l_point_wd .and. r_point_wd) then
                    ! below water level (left and right)
                    ! area of river[m2]
                    tmp_a = 0.5d0 * tmp_dx * &
                            (2.0d0 * h_table(i,n) - &
                             dan_z(j,n) - dan_z(j+1,n))
                    ! wetted perimeter[m]
                    tmp_s = sqrt(tmp_dx**2 + tmp_dz**2)
                    ! width
                    b_table(i,n) = b_table(i,n) + tmp_dx
                    ! roughness
                    if(composit_rn==1) then
                        tmp_rn_s = tmp_rn_s + &
                                   (dan_rn(j,n)**(3.0d0/2.0d0)) * tmp_s
                    endif
                
                elseif(l_point_wd .or. r_point_wd) then
                    if(abs(tmp_dz) < 0.00001d0) cycle
                    if(l_point_wd) then
                        ! below water level (only left)
                        intercect_x = dan_x(j,n) + &
                                      (h_table(i,n) - dan_z(j,n)) * &
                                      (tmp_dx / tmp_dz)
                        ! area of river[m2]
                        tmp_a = 0.5d0 * (intercect_x - dan_x(j,n)) * &
                                (h_table(i,n) - dan_z(j,n))
                        ! wetted perimeter[m]
                        tmp_s = sqrt((intercect_x - dan_x(j,n))**2 + &
                                     (h_table(i,n) - dan_z(j,n))**2)
                        ! width
                        b_table(i,n) = b_table(i,n) + (intercect_x - dan_x(j,n))
                        ! roughness
                        if(composit_rn==1) then
                            tmp_rn_s = tmp_rn_s + &
                                       (dan_rn(j,n)**(3.0d0/2.0d0)) * tmp_s
                        endif

                    else
                        ! below water level (only right)
                        intercect_x = dan_x(j+1,n) + &
                                      (h_table(i,n) - dan_z(j+1,n)) * &
                                      (tmp_dx / tmp_dz)
                        ! area of river[m2]
                        tmp_a = 0.5d0 * (dan_x(j+1,n) - intercect_x) * &
                                (h_table(i,n) - dan_z(j+1,n))
                        ! wetted perimeter[m]
                        tmp_s = sqrt((dan_x(j+1,n) - intercect_x)**2 + &
                                     (h_table(i,n) - dan_z(j+1,n))**2)
                        ! width
                        b_table(i,n) = b_table(i,n) + (dan_x(j+1,n) - intercect_x)
                        ! roughness
                        if(composit_rn==1) then
                            tmp_rn_s = tmp_rn_s + &
                                       (dan_rn(j,n)**(3.0d0/2.0d0)) * tmp_s
                        endif

                    endif
                else
                    ! above water level
                    cycle
                endif
            
                a_table(i,n) = a_table(i,n) + tmp_a
                s_table(i,n) = s_table(i,n) + tmp_s
                if(composit_rn==1) then
                    if(s_table(i,n) > 0.0d0) then
                        if(composit_rn==1) then
                            rn_table(i,n) = (tmp_rn_s / s_table(i,n))**(2.0d0/3.0d0)
                        elseif(composit_rn==0) then
                            rn_table(i,n) = dan_rn(1,n)
                        endif
                    else
                        write(*,*) ' UNST2D - WARNING : id ', n, ' - roughness is invaild value, replace 0.03d0.' 
                    endif
                endif
            enddo
            ! ---------------------------------------------
            ! hydraulic radius
            if(s_table(i,n) > 0.0d0) r_table(i,n) = a_table(i,n) / s_table(i,n)
        enddo
    enddo
    
    deallocate(s_table)
    deallocate(dan_x, dan_z, dan_rn)

    ! set dx and angle
    dx_1d = 0.0d0
    weir_dx_1d = 0.0d0
    w_angle_1d = 0.0d0
    do n = 1, ndan
        tmp_angle = 99999.0d0  ! river bed angle (1/tmp_angle)
        select case(ntype_1d(n))
        case (1,100)  ! downstream
            dx_1d(n) = kp_1d(n)
            weir_dx_1d(n) = (kp_1d(n+1) - kp_1d(n))*0.5d0  ! reference next
            ! same next
            if (rbed_1d(n) < rbed_1d(n+1)) tmp_angle = (kp_1d(n+1) - kp_1d(n)) / (rbed_1d(n+1) - rbed_1d(n))
            if((1.0d0 / tmp_angle) > (1.0d0/12000.0d0)) then
                w_angle_1d(n) = dcos(155.0d0 - 38.0d0 * log10(1.0d0/tmp_angle) * rivpi / 180.0d0)
            else
                w_angle_1d(n) = 1.0d0
            endif
        case (3)  ! distributary
            dx_1d(n) = kp_1d(n) - kp_1d(n-1)
            weir_dx_1d(n) = dx_1d(n)  ! same previous(test)
            if (rbed_1d(n) > rbed_1d(n-1)) tmp_angle = dx_1d(n) / (rbed_1d(n) - rbed_1d(n-1))
            if((1.0d0 / tmp_angle) > (1.0d0/12000.0d0)) then
                w_angle_1d(n) = dcos(155.0d0 - 38.0d0 * log10(1.0d0/tmp_angle) * rivpi / 180.0d0)
            else
                w_angle_1d(n) = 1.0d0
            endif
        case (-1)  ! upstream
            dx_1d(n) = kp_1d(n) - kp_1d(n-1)
            weir_dx_1d(n) = dx_1d(n)*0.5d0
            if (rbed_1d(n) > rbed_1d(n-1)) tmp_angle = dx_1d(n) / (rbed_1d(n) - rbed_1d(n-1))
            if((1.0d0 / tmp_angle) > (1.0d0/12000.0d0)) then
                w_angle_1d(n) = dcos(155.0d0 - 38.0d0 * log10(1.0d0/tmp_angle) * rivpi / 180.0d0)
            else
                w_angle_1d(n) = 1.0d0
            endif
        case (-2)  ! tributary
            dx_1d(n) = kp_1d(n)
            weir_dx_1d(n) = ((kp_1d(n+1) - kp_1d(n)) + dx_1d(n))*0.5d0
            do i = 1, nbound_1d
                if(n == b_idx_1d(i)) then
                    if (rbed_1d(n-1) > rbed_1d(subflow_input(2,i))) then
                        tmp_angle = dx_1d(n)/ (rbed_1d(n-1) - rbed_1d(subflow_input(2,i)))
                    endif
                    exit
                endif
            enddo
            if((1.0d0 / tmp_angle) > (1.0d0/12000.0d0)) then
                w_angle_1d(n) = dcos(155.0d0 - 38.0d0 * log10(1.0d0/tmp_angle) * rivpi / 180.0d0)
            else
                w_angle_1d(n) = 1.0d0
            endif
        case default  ! else 
            dx_1d(n) = kp_1d(n) - kp_1d(n-1)
            weir_dx_1d(n) = ((kp_1d(n+1) - kp_1d(n)) + dx_1d(n))*0.5d0
            if (rbed_1d(n) > rbed_1d(n-1)) tmp_angle = dx_1d(n) / (rbed_1d(n) - rbed_1d(n-1))
            if((1.0d0 / tmp_angle) > (1.0d0/12000.0d0)) then
                w_angle_1d(n) = dcos(155.0d0 - 38.0d0 * log10(1.0d0/tmp_angle) * rivpi / 180.0d0)
            else
                w_angle_1d(n) = 1.0d0
            endif
        end select
    enddo

    ! check
    open(newunit=foudan_table_unit, file = foudan_table, action = 'write')
    write(foudan_table_unit,'(i10)') max_tbl 
    write(foudan_table_unit,*) '          WL(m)  AREA(m2) RADIUS(m)  WIDTH(m) ROUGHNESS '
    do n = 1, ndan
        write(foudan_table_unit,'(2f10.3, i5, 3f10.2)') kp_1d(n), dx_1d(n), n_tbl(n), crown_1d(1,n), crown_1d(2,n), w_angle_1d(n)
        do i = 1, n_tbl(n)
            write(foudan_table_unit,'(5x, 6f10.2)') h_table(i,n), a_table(i,n), r_table(i,n), b_table(i,n), rn_table(i,n)
        enddo
    enddo
    close(foudan_table_unit)

    write(*,*) '    * Finish Calculating 1D hydraulic tables.'
    
end subroutine d1riv_table

subroutine read_oudan_header(line, kp, r)
    implicit none
    character(len=100), intent(in) :: line
    integer, intent(out) :: r
    real(8), intent(out) :: kp

    integer :: i
    integer :: j, tmp_idx
    character(len=100), allocatable:: tmp_value(:)

    allocate(tmp_value(16))

    kp = 0.0d0  ! kp
    r = 0

    i = 0
    j = 1
    do while (i < 16)
        i = i + 1
        tmp_idx = index(line(j:), ',')
        if (tmp_idx == 0) then  ! 最後のデータ
            tmp_value(i) = line(j:)
            exit
        else
            tmp_value(i) = line(j:j+tmp_idx-2)
            j = j + tmp_idx
        endif
    enddo

    ! out: kp
    if (len_trim(adjustl(tmp_value(1))) /= 0) read(tmp_value(1), *) kp
    kp = kp * 1000.0d0  ! km -> m
    ! out: num of point
    if (len_trim(adjustl(tmp_value(7))) /= 0) read(tmp_value(7), *) r
    
    deallocate(tmp_value)
    
end subroutine read_oudan_header

subroutine read_oudan_data(line, xc, zc, r)
    implicit none
    character(len=100), intent(in) :: line
    real(8), intent(out) :: xc, zc, r
    integer ios, idummy
    character(len=100) :: tmp_value

    xc = 0.0d0  ! x coord
    zc = 0.0d0  ! y coord
    r = 0.0d0  ! roughness

    tmp_value = trim(adjustl(line))
    read(tmp_value,*,iostat=ios) idummy, xc, zc, r
    if(ios/=0) then
        read(tmp_value,*,iostat=ios) idummy, xc, zc
        r = 0.030d0
    endif
    
end subroutine read_oudan_data

subroutine d1rivinitiald
    implicit none
    integer finit_unit
    integer i, ii, n, idummy, init_type
    real(8) init_depth

    h_1d = 0.0d0
    vv_1d = 0.0d0
    q_1d = 0.0d0
    r_1d = 0.0d0
    a_1d = 0.0d0
    b_1d = 0.0d0
    h_1dmax = 0.0d0
    vv_1dmax = 0.0d0
    subq_all = 0.0d0
    tributaryq_1d = 0.0d0
    breakq_1d = 0.0d0
    breakq_r_1d = 0.0d0
    breakq_l_1d = 0.0d0
    weirq_1d = 0.0d0
    weirq_l_1d = 0.0d0
    weirq_r_1d = 0.0d0
    pumpq_1d = 0.0d0
    sluiceq_1d = 0.0d0
    up_q_1d = 0.0d0

    open(newunit=finit_unit, file = finit, action = 'read')
    ii = 1
    read(finit_unit, *) init_type
    if(init_type==0) then
        ! only fractional
        do n = 1, nriv
            read(finit_unit, *) idummy, init_depth
            do i = ii, ii+riv_ndan(n)-1
                h_1d(i) = rbed_1d(i) + init_depth
            enddo
            ii = ii + riv_ndan(n)
        enddo
    elseif(init_type==1) then
        do n = 1, ndan
            read(finit_unit, *) idummy, h_1d(n), q_1d(n)
        enddo
    endif
    close(finit_unit)

end subroutine d1rivinitiald

end module unst_d1riv_read

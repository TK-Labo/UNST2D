! Write data
! Coded by K.Kawaike and TK Labo
! Release on July 7th 2025

!=======================================
! Set Output Format / Open & Close File
!=======================================
module unst_write_procedures
    use unst_globals_mod
    use d1riv_globals_mod

contains
    !-----------------------
    ! Write 2D nomal output
    !-----------------------
    subroutine write_array_data(unit_num, data, fmt_data)
        implicit none
        integer, intent(in) :: unit_num
        real(8), intent(in) :: data(:)
        character(*), intent(in) :: fmt_data
        integer :: me

        write(unit_num, '(a,f8.0,a)') ' time=', unsttime, '(s)'
        write(unit_num, fmt_data) (data(me), me = 1, mesh)  ! meshもグローバル変数
    end subroutine write_array_data

    !-----------------------
    ! Write paddy output(1)
    !-----------------------
    subroutine write_paddyarray_data(unit_num, data, fmt_data)
        implicit none
        integer, intent(in) :: unit_num
        real(8), intent(in) :: data(:)
        character(*), intent(in) :: fmt_data
        integer :: i

        write(unit_num, '(a,f8.0,a)') ' time=', unsttime, '(s)'
        write(unit_num, fmt_data) (data(i), i = 1, paddy)  ! meshもグローバル変数
    end subroutine write_paddyarray_data

    !-----------------------
    ! Write paddy output(2)
    !-----------------------
    subroutine write_paddy_data(unit_num)
        implicit none
        integer, intent(in) :: unit_num
        integer :: me, i

        write(unit_num, 1040) unsttime
        do me = 1, nhp
            do i = 1, 72000
                if(dhp(i,me) > 0.0d0) write(unit_num, 1041) i, dhp(i,me), me
            enddo
        enddo
    1040 format(' time=', f8.0, '(s)')
    1041 format(i8, f10.5, i8)
    end subroutine write_paddy_data

    !-----------------------
    ! Write 1D River output
    !-----------------------
    subroutine write_multi_data(unit_num,dnum,d1,d2,d3,d4,d5,d6,d7)
        implicit none
        integer, intent(in) :: unit_num, dnum
        real(8), intent(in) :: d1(dnum), d2(dnum), d3(dnum), d4(dnum), d5(dnum), &
                                d6(dnum), d7(dnum)
        integer :: i
        
        write(unit_num, 2111) unsttime
        do i = 1, dnum
            !write(unit_num, 2112) i, hs(i), qs(i), qys(i)
            ! d1: kp
            ! d2: river bed
            ! d3: water level
            ! d4: q
            ! d5: velocity
            ! d6: subq
            ! d7: area
            write(unit_num, 2112) i, d1(i), d3(i), d3(i)-d2(i), -d4(i), -d5(i), d6(i), d7(i)
        enddo
    2111 format(f12.1, '  time(sec)')
    2112 format(i5, f10.3,10f10.5)
    !2112 format(i4, 10f14.8)
    end subroutine write_multi_data

    !------------------
    ! Open Output file
    !------------------
    subroutine open_unst_output_files
        implicit none

        ! output file set
        open(newunit = fh_unit, file = fh, action = 'write')
        open(newunit = fhmx_unit, file = fhmx, action = 'write')
        open(newunit = fuum_unit, file = fuum, action = 'write')
        open(newunit = fvvm_unit, file = fvvm, action = 'write')
        open(newunit = fuumx_unit, file = fuumx, action = 'write')
        open(newunit = fvvmx_unit, file = fvvmx, action = 'write')
        open(newunit = fstorage_unit, file = fstorage, action = 'write')
        open(newunit = fq_unit, file = fq, action = 'write')

        if(paddydam==1) then
            open(newunit = fdhp_unit, file = fdhp, action = 'write')
            open(newunit = fpaddyh_unit, file = fpaddyh, action = 'write')
            open(newunit = fpaddyq_unit, file = fpaddyq, action = 'write')
        endif

        if(d1riv==1) then
            open(newunit = fd1out_unit, file = fd1out, action = 'write')
            open(newunit = fd1mx_unit, file = fd1mx, action = 'write')
            write(fd1out_unit, '(A)') &
                '  rid     kp(m)     wl(m)  depth(m)   q(m3/s) velo(m/s)  subq(m3)'
            write(fd1mx_unit, '(A)') &
                '  rid     kp(m)     wl(m)  depth(m)   q(m3/s) velo(m/s)  subq(m3)'
        endif

    end subroutine open_unst_output_files

    !-------------------
    ! Close Output file
    !-------------------
    subroutine close_unst_output_files
        implicit none

        close(fh_unit)
        close(fhmx_unit)
        close(fuum_unit)
        close(fvvm_unit)
        close(fuumx_unit)
        close(fvvmx_unit)
        close(fstorage_unit)
        close(fq_unit)

        if(paddydam==1) then
            close(fdhp_unit)
            close(fpaddyh_unit)
            close(fpaddyq_unit)
        endif

        if(d1riv==1) then
            close(fd1out_unit)
            close(fd1mx_unit)
        endif

    end subroutine close_unst_output_files

end module unst_write_procedures
    
!=====================
! Output Main Process
!=====================
module unst_wrfile
    use unst_globals_mod
    use unst_write_procedures
    implicit none

contains
    !--------------------
    ! Display Indication
    !--------------------
    subroutine dispwrite
        implicit none
        integer me
        real(8) sv, sa, saj
    
        ! inundation water level (Former subroutine sumqa (UNST_Sub.f90))  v.1.0.5
        sv = 0.0d0
        sa = 0.0d0
        saj = 0.0d0
        !$omp parallel do default(shared),private(me),reduction(+:sv, sa, saj)
        do me = 1, mesh
            if(inf(me) == 0) cycle
            sv = sv + unsth(me)*smesh(me)*(1.0d0 - lambda(me))
            if(unsth(me) > 1.0d-3) sa  = sa  + smesh(me)
            if(unsth(me) > 1.0d-3) saj = saj + smesh(me)*(1.0d0 - lambda(me))  ! flooded area
        enddo
        !$omp end parallel do

        write(*, 1000) unsttime, sv, unst_error_v, unst_dis_v
1000 format('UNST------   unsttime=', f8.0, '(s)', 10f18.4)
! 1005 format('   unsttime=', f8.0, '(s)', 2f18.4, 4i10)
! 1099 format('   unsttime=', f8.0, '(s)', 10f21.7)
    end subroutine dispwrite

    !----------------------
    ! Writing Data to File
    !----------------------
    subroutine diskwrite
        implicit none
        character(len=20) :: fmt1
        character(len=20) :: fmt2

        fmt1 = '(10f15.3)'
        fmt2 = '(10f15.7)'
    
        !$omp parallel
        !$omp sections
        !$omp section
        call write_array_data(fh_unit, unsth, fmt1)

        !$omp section
        call write_array_data(fuum_unit, uum, fmt1)

        !$omp section
        call write_array_data(fvvm_unit, vvm, fmt1)

        !$omp section
        call write_array_data(fstorage_unit, ((unsth*smesh) - qr_sum), fmt2)

        !$omp section
        call write_array_data(fq_unit, qr_sum, fmt1)

        !$omp end sections
        !$omp end parallel
        if(d1riv==1) call write_multi_data(fd1out_unit, ndan, kp_1d, rbed_1d, &
                                            h_1d, q_1d, vv_1d, subq_all, a_1d)
    end subroutine diskwrite

    !----------------------------
    ! Writing Paddy Data to File
    !----------------------------
    subroutine paddywrite
        implicit none
        character(len=20) :: fmt1, fmt2

        fmt1 = '(10f15.9)'
        fmt2 = '(10f15.7)'

        !$omp parallel
        !$omp sections
        !$omp section
        call write_paddy_data(fdhp_unit)

        !$omp section
        call write_paddyarray_data(fpaddyh_unit, pqh, fmt1)

        !$omp section
        call write_paddyarray_data(fpaddyq_unit, totalqp,  fmt2)
        !$omp end sections
        !$omp end parallel
    end subroutine paddywrite

    !-------------------------------------
    ! Maximum depth of inundation written
    ! Maximum water level written
    !-------------------------------------
    subroutine wrhmax
        implicit none
        character(len=20) :: fmt1

        fmt1 = '(10f15.3)'

        !$omp parallel
        !$omp sections
        !$omp section
        call write_array_data(fhmx_unit, hmax, fmt1)

        !$omp section
        call write_array_data(fuumx_unit, uummax, fmt1)

        !$omp section
        call write_array_data(fvvmx_unit, vvmmax, fmt1)
        !$omp end sections
        !$omp end parallel

        ! v.1.0.5
        if(d1riv==1) call write_multi_data(fd1mx_unit, ndan, kp_1d, rbed_1d, &
                                            h_1dmax, q_1d, vv_1dmax, subq_all, a_1d)
    end subroutine wrhmax

end module unst_wrfile

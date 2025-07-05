! write data
! Coded by K.Kawaike and TK Labo
! Release on July 7th 2025

module write_procedures
    use globals
    contains
    ! 配列全体書き込み用サブルーチン
    subroutine write_array_data(unit_num, data, fmt_data)
        implicit none
        integer, intent(in) :: unit_num
        real(8), intent(in) :: data(:)
        character(*), intent(in) :: fmt_data
        integer :: me
        
        write(unit_num, '(a,f8.0,a)') ' time=', unsttime, '(s)' 
        write(unit_num, fmt_data) (data(me), me = 1, mesh)  ! meshもグローバル変数
    end subroutine
    
    subroutine write_paddyarray_data(unit_num, data, fmt_data)
        implicit none
        integer, intent(in) :: unit_num
        real(8), intent(in) :: data(:)
        character(*), intent(in) :: fmt_data
        integer :: i
        
        write(unit_num, '(a,f8.0,a)') ' time=', unsttime, '(s)' 
        write(unit_num, fmt_data) (data(i), i = 1, paddy)  ! meshもグローバル変数
    end subroutine

    ! 複合データ書き込み用サブルーチン UNST-2D original
    subroutine write_multi_data(unit_num)
        implicit none
        integer, intent(in) :: unit_num
        integer :: i
        
        write(unit_num, 2111) unsttime
        do i = 1, msctn
            write(unit_num, 2112) i, hs(i), qs(i), qys(i)
        enddo
    2111 format(f12.1, '  time(sec)')
    2112 format(i4, 3f14.8)
    end subroutine

    ! paddy書き込み用サブルーチン
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
    end subroutine
end module write_procedures
    
subroutine wrfile
    use globals
    use write_procedures
    real(8) sv, sa, saj

    ! display indication
    entry dispwrite
    call sumqa(sv, sa, saj)
    write(*, 1000) unsttime, sv
    if(str_type == 1) then
    write(79, 1000) unsttime, sv_a, dvr, qout, sv_a + sumstr45, sv0 + vrain + dvr - qout, sv_a+sumstr45-sv0-vrain-dvr+qout, v_minus_all
    !write(71,1000) unsttime, sv, sa, saj
    ! write(71,1000) unsttime, sv, sv3, sv6, sv26, sumstr45
    ! write(73,1000) unsttime, sv3, v_minus(31)+v_minus(32)
    ! write(76,1005) unsttime, sv6, v_minus(61), count1, count2, count3, count4
    ! write(72,1000) unsttime, sv26, v_minus(26)
    ! write(74,1000) unsttime, sumstr45, sumstr4, sumstr5
    write(78,1000) unsttime, sv_a, v_minus_all, v_minus(21), v_minus(22), v_minus(23), v_minus(24), v_minus(25)
    write(77,1099) unsttime, tracer, svc, tracer-svc, v_cminus, v_cplus, v_cextre
    endif
1000 format('UNST------   unsttime=', f8.0, '(s)', 7f18.4)
! 1005 format('   unsttime=', f8.0, '(s)', 2f18.4, 4i10)
1099 format('   unsttime=', f8.0, '(s)', 10f21.7)
    return

    ! writing data to a file
    entry diskwrite
    !$omp parallel
    !$omp sections
    !$omp section
    call write_array_data(91, unsth, '(10f15.3)') 

    !$omp section
    if(d1riv==1) call write_multi_data(81)

    !$omp section
    call write_array_data(94, uum, '(10f15.3)')

    !$omp section
    call write_array_data(95, vvm, '(10f15.3)')

    !$omp section
    call write_array_data(101, ((unsth*smesh) - qr_sum), '(10e15.7)')

    !$omp section
    call write_array_data(102, qr_sum, '(10f15.3)')

    !$omp section
    if(str_type == 1) call write_array_data(103, unstc, '(10f7.4)')
    
    !$omp end sections
    !$omp end parallel
    return

    entry wrhmax
    !$omp parallel
    !$omp sections
    !$omp section
    call write_array_data(93, hmax,  '(10f15.3)')

    !$omp section
    call write_array_data(96, uummax,  '(10f15.3)')

    !$omp section
    call write_array_data(97, vvmmax,  '(10f15.3)')
    !$omp end sections
    !$omp end parallel
    return

    entry paddywrite
    !$omp parallel
    !$omp sections
    !$omp section
    call write_paddy_data(98)

    !$omp section
    call write_paddyarray_data(99, pqh,  '(10f15.9)')

    !$omp section
    call write_paddyarray_data(100, paddy_q,  '(10f10.5)')
    !$omp end sections
    !$omp end parallel
    return

    end subroutine wrfile

! data read
! coded by k.kawaike and TK Labo
! Release on July 7th 2025

module unst_read
    use unst_globals_mod
    implicit none
    
contains

subroutine unst_rdat
    implicit none

    integer no, me, li, k, i, j, tt, dummy, read_count, ios

    ! file paths
    character(len=50) :: &
        fnode, flink, fmesh, finf, fbs, fqin, fmesh2ij, &
        frn, fcntl, frain

    ! unit numbers
    integer fnode_unit, flink_unit, fmesh_unit, finf_unit, &
        fbs_unit, fqin_unit, fmesh2ij_unit, &
        frn_unit, fcntl_unit, frain_unit

    character(len=100) :: line

    integer, allocatable :: infmn(:)
    real(8), allocatable :: temp_mn(:)
    real(8), allocatable :: temp_soildepth(:), temp_gammaa(:), temp_ksv(:), temp_faif(:)
    integer n, numrn
    integer sep_rtuv

    real(8) dkout, dpout
    real(8) tday, thour, tmin, tsec
    integer ino_qin, t ! UNST-2D only
    integer nx_rain, ny_rain, tt_max_rain
    real(8) xllcorner_rain, yllcorner_rain, cellsize_rain_x, cellsize_rain_y

    real(8) rdummy ! temp variable UNST-2D only

    ! config file
    fcntl = 'cntl.dat'

    !------------------------------
    ! Step1 : Read UNST2D_cntl.dat
    !------------------------------
    ! parameter set
    open(newunit=fcntl_unit, file=fcntl, action = 'read')
    read(fcntl_unit, 1004) tday, thour, tmin, tsec
    read(fcntl_unit, 1002) unstdt
    read(fcntl_unit, 1002) dkout
    read(fcntl_unit, 1002) dpout
    read(fcntl_unit, 1002) dtrain
    read(fcntl_unit, 1002) dtq
    read(fcntl_unit, 1001) inls
    read(fcntl_unit, *)
    write(*,*)
    write(*,*) '== UNST Parameters =='
    if(inls==1) write(*,*) ' - Qin Mode : Multi'
    if(inls==0) write(*,*) ' - Qin Mode : One Point'
    write(*,*) '* Time & Step Parameters *'
    write(*,'(a18,4f6.1)') ' - Time         :', tday, thour, tmin, tsec
    write(*,'(a18,f10.2)') ' - unst dt      :', unstdt
    write(*,'(a18,f10.2, a3)') ' - disk out     :', dkout, '(s)'
    write(*,'(a18,f10.2, a3)') ' - disk out     :', dkout, '(s)'
    write(*,'(a18,f10.2, a3)') ' - dt rain      :', dtrain,'(s)'
    write(*,'(a18,f10.2, a3)') ' - dt qin       :', dtq,'(s)'
    write(*,*)

    ! required filename set
    read(fcntl_unit, 1009) fnode
    read(fcntl_unit, 1009) flink
    read(fcntl_unit, 1009) fmesh
    read(fcntl_unit, 1009) finf
    read(fcntl_unit, 1009) fbs
    read(fcntl_unit, 1009) fqin
    read(fcntl_unit, 1009) frain
    read(fcntl_unit, 1009) fmesh2ij
    read(fcntl_unit, 1009) frn
    read(fcntl_unit, *)
    write(*,*) '* Requred Data *'
    write(*,'(a16, a)') ' - node       :', fnode
    write(*,'(a16, a)') ' - link       :', flink
    write(*,'(a16, a)') ' - mesh       :', fmesh
    write(*,'(a16, a)') ' - inf        :', finf
    write(*,'(a16, a)') ' - bs         :', fbs
    write(*,'(a16, a)') ' - qin        :', fqin
    write(*,'(a16, a)') ' - rain       :', frain
    write(*,'(a16, a)') ' - mesh2ij    :', fmesh2ij
    write(*,'(a16, a)') ' - infilt(rn) :', frn
    write(*,*)

    read(fcntl_unit, *) xllcorner_rain
    read(fcntl_unit, *) yllcorner_rain
    read(fcntl_unit, *) cellsize_rain_x, cellsize_rain_y
    read(fcntl_unit, *)
    read(fcntl_unit, *)
    write(*,'(a18)') '* Rain Propaties *'
    write(*,'(a18,f10.2)') '   > xllcorner  :', xllcorner_rain
    write(*,'(a18,f10.2)') '   > yllcorner  :', yllcorner_rain
    write(*,'(a18,f10.2)') '   > cellsize-x :', cellsize_rain_x
    write(*,'(a18,f10.2)') '   > cellsize-y :', cellsize_rain_y
    write(*,*)

    ! option switch & filename set
    write(*,*) '* Option *'
    ! option : discharge
    read(fcntl_unit, 1003) dsmesh
    read(fcntl_unit, 1009) fdsmesh
    read(fcntl_unit, *)
    read(fcntl_unit, *)
    if(dsmesh==1) then
        write(*,*) '+ Discharge(dsmesh) Active'
        write(*,'(a16, a)') ' - dsmesh     :', fdsmesh
    endif

    ! option : 1D river
    read(fcntl_unit, 1003) d1riv
    read(fcntl_unit, 1009) fd1riv_cntl
    read(fcntl_unit, *)
    read(fcntl_unit, *)
    if(d1riv==1) then
        write(*,*) '+ 1D River(d1riv) Active'
        write(*,'(a16, a)') ' - d1riv_cntl :', fd1riv_cntl
    endif

    ! option : vegetation resistance ; considering lodging, ex)reed
    read(fcntl_unit, 1003) plantFN
    read(fcntl_unit, 1009) fplantF
    read(fcntl_unit, 1009) fplantN
    read(fcntl_unit, *)
    read(fcntl_unit, *)
    if(plantFN==1) then
        write(*,*) '+ Vegetation Resistance(plantFN) Active'
        write(*,'(a16, a)') ' - plantF     :', fplantF
        write(*,'(a16, a)') ' - plantN     :', fplantN
    endif

    ! option : vegetation resistance ; flood defence forest ex)tree, bamboo
    read(fcntl_unit, 1003) plantDa
    read(fcntl_unit, 1009) fplantD
    read(fcntl_unit, 1009) fplanta
    read(fcntl_unit, 1009) fplanthv
    read(fcntl_unit, 1009) fplantcd
    read(fcntl_unit, *)
    read(fcntl_unit, *)
    if(plantDa==1) then
        write(*,*) '+ Vegetation Resistance(plantDa) Active'
        write(*,'(a16, a)') ' - plantD     :', fplantD
        write(*,'(a16, a)') ' - planta     :', fplanta
        write(*,'(a16, a)') ' - planthv    :', fplanthv
        write(*,'(a16, a)') ' - plantcd    :', fplantcd
    endif

    ! option : paddydam
    read(fcntl_unit, 1003) paddydam
    read(fcntl_unit, 1009) fpaddy
    read(fcntl_unit, 1009) fpqout
    read(fcntl_unit, 1009) fpaddy_param
    read(fcntl_unit, *)
    read(fcntl_unit, 1009) fdhp
    read(fcntl_unit, 1009) fpaddyh
    read(fcntl_unit, 1009) fpaddyq
    read(fcntl_unit, *)
    read(fcntl_unit, *)
    if(paddydam==1) then
        write(*,*) '+ Paddydam Active'
        write(*,'(a16, a)') ' - paddy      :', fpaddy
        write(*,'(a16, a)') ' - pqout      :', fpqout
        write(*,'(a16, a)') ' - paddy_param:', fpaddy_param
        write(*,'(a24, a)') ' - (Output) dhp       :', fdhp
        write(*,'(a24, a)') ' - (Output) paddyh    :', fpaddyh
        write(*,'(a24, a)') ' - (Output) paddyq    :', fpaddyq
    endif

    ! option : drainage area
    read(fcntl_unit, 1003) drainarea
    read(fcntl_unit, 1009) finf_dr
    read(fcntl_unit, 1009) fdrain
    read(fcntl_unit, *)
    read(fcntl_unit, *)
    if(drainarea==1) then
        write(*,*) '+ Drainage Area(drainarea) Active'
        write(*,'(a16, a)') ' - inf_dr     :', finf_dr
        write(*,'(a16, a)') ' - drain      :', fdrain
    endif

    ! option : line embankment
    read(fcntl_unit, 1003) morid
    read(fcntl_unit, 1009) fmorid
    read(fcntl_unit, *)
    read(fcntl_unit, *)
    if(morid==1) then
        write(*,*) '+ Line Embankment(morid) Active'
        write(*,'(a16, a)') ' - morid      :', fmorid
    endif

    if(dsmesh==0 .and. plantFN==0 .and. plantDa==0 .and. paddydam==0 .and. &
       drainarea==0 .and. morid==0) then
        write(*,*) ' - No Options'
    endif
    write(*,*)

    ! output filename set
    read(fcntl_unit, 1009) fh
    read(fcntl_unit, 1009) fhmx
    read(fcntl_unit, 1009) fuum
    read(fcntl_unit, 1009) fvvm
    read(fcntl_unit, 1009) fuumx
    read(fcntl_unit, 1009) fvvmx
    read(fcntl_unit, 1009) fstorage
    read(fcntl_unit, 1009) fq

1001 format(7x, i11)
1002 format(7x, f11.2)
1003 format(17x, i1)
1004 format(7x, 4f6.1)
1009 format(a50)
    close(fcntl_unit)

    !----------------------------
    ! Step2 : Set constant param
    !----------------------------
    timmax = 3600.0d0*24.0d0*tday + 3600.d0*thour + 60.d0*tmin + tsec
    dt2 = 2.0d0*unstdt

    lpout = int(dpout/unstdt)  ! display write dt
    lkout = int(dkout/unstdt)  ! disk write dt

    !--------------------------
    ! Step3 : Read input files
    !--------------------------
    ! open input files
    open(newunit = fnode_unit, file = fnode, action = 'read')
    open(newunit = flink_unit, file = flink, action = 'read')
    open(newunit = fmesh_unit, file = fmesh, action = 'read')
    open(newunit = finf_unit, file = finf, action = 'read')
    open(newunit = fbs_unit, file = fbs, action = 'read')
    open(newunit = fqin_unit, file = fqin, action = 'read')
    open(newunit = fmesh2ij_unit, file = fmesh2ij, action = 'read')
    open(newunit = frn_unit, file = frn, action = 'read')
    open(newunit = frain_unit, file = frain, action = 'read')  ! only single unst

    write(*,*) '* grid data *'
    ! -- 01 grid : node (vertex coordinates) --
    read(fnode_unit, *) node
    allocate(dnox(node), dnoy(node))
    do no = 1, node
        read(fnode_unit, 1212) dnox(no), dnoy(no)
    enddo
1212 format(8x, 2f11.2)  ! format change v.1.0.5
    close(fnode_unit)
    write(*,*) ' - node : ', node, 'points'
    ! ---

    ! -- 01 grid : link (edge) --
    read(flink_unit, *) link
    allocate(limesh(link, 2), linode(link, 2))
    allocate(scv(link), rthl(link, 2))
    allocate(ux(link), uy(link))
    do li = 1, link
        read(flink_unit, 1222) limesh(li, 1), limesh(li, 2), linode(li, 1), linode(li, 2)
        read(flink_unit, 1223) scv(li), rthl(li, 1), rthl(li, 2)
        read(flink_unit, 1224) ux(li), uy(li)
    enddo
1222 format(8x, 2i8, 2x, 2i8)
1223 format(8x, f10.2, 2f10.5)
1224 format(18x, 2f10.5)
    close(flink_unit)
    write(*,*) ' - link : ', link, 'links'
    ! ---

    ! -- 01 grid : mesh --
    read(fmesh_unit, 1009, iostat = ios) line
    line = trim(adjustl(line))
    read(line, *, iostat = ios) mesh, sep_rtuv
    if(ios/=0) then
        read(line, *, iostat = ios) mesh
        sep_rtuv = 0
        if(ios/=0) stop 'UNST2D - ERROR-1001 : mesh could not be loaded.'
    endif
    write(*,*) ' - mesh : ', mesh, 'grids'
    write(*,*) '   > interpolation type:', sep_rtuv

    allocate(ko(mesh), menode(mesh, 6), melink(mesh, 6), smesh(mesh), xmesh(mesh), ymesh(mesh))
    allocate(rtuv_x(mesh, 6), rtuv_y(mesh, 6))
    do me = 1, mesh
        read(fmesh_unit, 1232) ko(me), (menode(me, k), k = 1, ko(me))
        read(fmesh_unit, 1233) (melink(me, k), k = 1, ko(me))
        read(fmesh_unit, 1234) smesh(me), xmesh(me), ymesh(me)
        if (sep_rtuv==0) then
            ! weight distance only
            read(fmesh_unit, 1235) (rtuv_x(me, k), k = 1, ko(me))
            rtuv_y(me, 1:ko(me)) = rtuv_x(me, 1:ko(me))
        elseif (sep_rtuv==1) then
            ! weight distance and angle added by d.baba
            read(fmesh_unit, 1235) (rtuv_x(me, k), k = 1, ko(me))
            read(fmesh_unit, 1235) (rtuv_y(me, k), k = 1, ko(me))
        endif
    enddo
1232 format(8x, i5, 5x, 20i8)
1233 format(13x, 20i8)
1234 format(6x, 3f15.2)
1235 format(13x, 20f10.5)
    close(fmesh_unit)
    ! ---

    ! -- 02 attribute(mesh) : inf (landuse infromation) --  ! update v.1.0.5
    allocate(inf(mesh), lambda(mesh))
    do me = 1, mesh
        read(finf_unit, 1009, iostat = ios) line
        if(ios/=0) exit
        read(line, *, iostat = ios) inf(me), lambda(me)
        if(ios/=0) then
            read(line, *, iostat = ios) inf(me)
            lambda(me) = 0.0d0
        endif
    enddo
    close(finf_unit)
    ! ---

    ! -- 02 attribute(mesh) : baseo (elevation) --
    allocate(baseo(mesh))
    do me = 1, mesh
        read(fbs_unit, *) dummy, baseo(me)
    enddo
    close(fbs_unit)
    ! ---

    ! -- 02 attribute(mesh) : rain --  ! update v.1.0.5
    read(frain_unit, *) t, nx_rain, ny_rain
    do i = 1, ny_rain
        read(frain_unit, *, iostat = ios) (rdummy, j=1, nx_rain)
    enddo
    rewind(frain_unit)
    tt_max_rain = int(timmax/dtrain)

    allocate(rain(0:tt_max_rain, 0:ny_rain, 0:nx_rain))
    do tt = 0, tt_max_rain
        read(frain_unit, *) t, nx_rain, ny_rain
        if(tt == tt_max_rain) write(*,*) 'rain_time = ', t
        do i = 1, ny_rain
            read(frain_unit, *) (rain(tt, i, j), j = 1, nx_rain)
        enddo
    enddo
    close(frain_unit)
    ! convert from (mm/h) to (m/s)
    rain = rain / 3600.d0 / 1000.d0
    
    allocate(urain_j(mesh), urain_i(mesh))
    allocate(rri_x(mesh), rri_y(mesh))
    do me = 1, mesh
     read(fmesh2ij_unit, *) rri_x(me), rri_y(me)
    end do
    close(fmesh2ij_unit)

    do me = 1, mesh
        urain_j(me) = int( (rri_x(me) - xllcorner_rain) / cellsize_rain_x) + 1
        urain_i(me) = ny_rain - int( (rri_y(me) - yllcorner_rain) / cellsize_rain_y)
        if(urain_i(me) <= 0 .or. urain_j(me) <= 0 &
           .or. urain_i(me) > ny_rain .or. urain_j(me) > nx_rain) then
            urain_i(me) = 0
            urain_j(me) = 0
        endif
    enddo
    if(any(urain_i==0) .or. any(urain_j==0)) then
        write(*,*) ' UNST2D - WARNIG : No overlap between UNST2D and rainfall area found.'
    endif
    ! ---

    ! -- 02 attribute : roughness --
    allocate( mn(mesh) )
    allocate( unst_soildepth(mesh), unst_gammaa(mesh), unst_ksv(mesh), unst_faif(mesh), unst_infilt_limit(mesh) )
    allocate( unst_gampt_ff(mesh), unst_gampt_f(mesh) )
    read(frn_unit, *) numrn
    read(frn_unit, *)
    
    allocate( infmn(numrn), temp_mn(numrn) )
    allocate( temp_soildepth(numrn), temp_gammaa(numrn), temp_ksv(numrn), temp_faif(numrn) )
    
    temp_soildepth = 0.0d0
    temp_gammaa = 0.0d0
    temp_ksv = 0.0d0
    temp_faif = 0.0d0
    
    do n = 1, numrn
        read(frn_unit, '(A)', iostat = ios) line
        if (ios /= 0) exit
        read_count = 0
        do i = 1, len_trim(line)
            if (line(i:i) == ' ' .and. line(i+1:i+1) /= ' ') then
                read_count = read_count + 1
            endif
        enddo
        read_count = read_count + 1
        if (read_count >= 2) then
            read(line, *) infmn(n), temp_mn(n)
        endif
        if (read_count >= 7) then
            read(line, *) infmn(n), temp_mn(n), temp_soildepth(n), temp_gammaa(n), temp_ksv(n), temp_faif(n)
        endif
    enddo
    close(frn_unit)

    !$omp parallel do default(shared),private(me, n)
    do me = 1, mesh
        unst_soildepth(me) = 0.0d0
        unst_gammaa(me) = 0.0d0
        unst_ksv(me) = 0.0d0
        unst_faif(me) = 0.0d0
        unst_infilt_limit(me) = 0.0d0
        unst_gampt_ff(me) = 0.0d0
        unst_gampt_f(me) = 0.0d0
        
        do n = 1, numrn
            if(inf(me) == infmn(n)) then
                mn(me) = temp_mn(n)
                unst_soildepth(me) = temp_soildepth(n)
                unst_gammaa(me) = temp_gammaa(n)
                unst_ksv(me) = temp_ksv(n)
                unst_faif(me) = temp_faif(n)
                if( unst_soildepth(me) .gt. 0.d0 .and. unst_ksv(me) .gt. 0.d0 ) unst_infilt_limit(me) = unst_soildepth(me) * unst_gammaa(me)
            endif
        enddo
    enddo
    !$omp end parallel do
    
    deallocate( infmn, temp_mn, temp_soildepth, temp_gammaa, temp_ksv, temp_faif )
    ! ---

    ! -- 03 attribute(link) : rbeta(passage rate) --  ! update v.1.0.5
    allocate(rbeta(link))
    rbeta = 1.0d0
    ! ---

    ! -- 03 attribute(link) : qin --  ! update v.1.0.5
    if(inls==1) then
        read(fqin_unit, *) iqnum, iqin
        allocate(inl(iqnum), qin(iqnum, iqin+1))
        allocate(lkyokai_dx(iqnum), lkyokai_dy(iqnum), lkyokai_dir(iqnum))

        qin = 0.0d0

        do i = 1, iqnum
            read(fqin_unit, *) inl(i)
        enddo
        read(fqin_unit, *)
        do i =1, iqnum
            read(fqin_unit, *) (qin(i,ino_qin),ino_qin=1,iqin) 
            if(i < iqnum) read(fqin_unit, *)
        enddo
    elseif(inls == 0) then
        read(fqin_unit, *) iqin
        allocate(inl(1), qin(1,iqin+1))
        allocate(lkyokai_dx(1), lkyokai_dy(1), lkyokai_dir(1))

        qin = 0.0d0

        read(fqin_unit,*) inl(1)
        read(fqin_unit,*)
        do i = 1, iqin
            read(fqin_unit, *) qin(1,i)
        enddo
    endif
    close(fqin_unit)
    ! ---

    ! allocate variables
    allocate(unsth(mesh), ho(mesh), hl(link), hmax(mesh), uummax(mesh), vvmmax(mesh))
    allocate(umm(0:mesh), vnm(0:mesh))
    allocate(uum(mesh), vvm(mesh))
    allocate(node_dx(mesh,6), node_dy(mesh,6))
    allocate(rnof(mesh), qr_sum(mesh), riv_eq(mesh))
    allocate(um(link), umo(link), uu(link), vn(link), vno(link), vv(link))
    allocate(umbeta(link), vnbeta(link))
    allocate(lhan(link), lhano(link))
    allocate(rnx(link), dl(link))
    allocate(blink(link))

end subroutine unst_rdat

!==================================================
! Read Vegetation Resistance ; considering lodging  ! development 25/12
!==================================================
subroutine plantFNdat
    implicit none
    integer me, fplantF_unit, fplantN_unit

    allocate(plantF_array(mesh), plantN_array(mesh))

    ! plant resistance / 1  (per one reed)
    open(newunit = fplantF_unit, file = fplantF, action = 'read')
    do me = 1, mesh
          read(fplantF_unit, *) plantF_array(me)
    enddo
    close(fplantF_unit)

    ! number of plants
    open(newunit = fplantN_unit, file = fplantN, action = 'read')
    do me = 1, mesh
        read(fplantN_unit, *) plantN_array(me)
    enddo
    close(fplantN_unit)

end subroutine plantFNdat

!===================================================
! Read Vegetation Resistance ; flood defence forest  ! fixed & development 25/12
!===================================================
!--------------------------------------------------------------------------
! UNST2D下水道・圃場読み込みサブルーチン/UNST2D drainage process read subroutine
! 変数/Variable:
! > グローバル変数/global variable
!   plant_a_array(mesh): 樹林密度(本/m3)/
!   plant_D_array(mesh): 胸高直径(m)/
!   plant_hv_array(mesh): mesh上の植生の有効樹高(m)/
!   vr_cd(link): 抗力係数
!   plant_lambda(mesh): 単位体積に占める植生の遮断面積/
!   vr_hv(link): link上の植生の有効樹高(m)/
!   dk_val(link): 植生抵抗の係数
!--------------------------------------------------------------------------
subroutine plantDadat
    implicit none
    integer me, li
    integer fplantD_unit, fplanta_unit, fplanthv_unit, fplantcd_unit

    allocate(plant_D_array(mesh), plant_a_array(mesh), plant_lambda(mesh))
    allocate(plant_hv_array(mesh), vr_cd(link))
    allocate(vr_hv(link), dk_val(link))

    ! breast height diarneter　(m)
    open(newunit = fplantD_unit, file = fplantD, action = 'read')
    do me = 1, mesh
        read(fplantD_unit, *) plant_D_array(me)
    enddo
    close(fplantD_unit)

    ! density (tree / m^2)
    open(newunit = fplanta_unit, file = fplanta, action = 'read')
    do me = 1, mesh
        read(fplanta_unit, *) plant_a_array(me)
    enddo
    close(fplanta_unit)

    ! height (m)
    open(newunit = fplanthv_unit, file = fplanthv, action = 'read')
    do me = 1, mesh
        read(fplanthv_unit, *) plant_hv_array(me)
    enddo
    close(fplanthv_unit)

    ! C_D
    open(newunit = fplantcd_unit, file = fplantcd, action = 'read')
    do li = 1, link
        read(fplanta_unit, *) vr_cd(me)
    enddo
    close(fplantcd_unit)
!
end subroutine plantDadat

!===============
! Read Paddydam
!===============
subroutine paddydat
    implicit none
    integer me, pa, ii
    integer fpaddy_unit, fpqout_unit, fpaddy_param_unit

    ! paddyid(unique id) read
    allocate(paddyid(mesh))
    open(newunit = fpaddy_unit, file = fpaddy, action = 'read')
    do me = 1, mesh
        read(fpaddy_unit,*) paddyid(me)
    enddo
    close(fpaddy_unit)

    ! distance from waterway (m)
    open(newunit = fpqout_unit, file = fpqout, action = 'read')
    read(fpqout_unit, *) paddy
    read(fpqout_unit, *)
    allocate(pqout_idx(paddy), pdrain(paddy), min_pmeshid(paddy), device(paddy))
    allocate(orifice_num(paddy), min_dist(paddy), etp(paddy))
    allocate(dr_dist(paddy), pqh(paddy), paddy_q(paddy), totalqp(paddy))
    if(paddy>0) then
    do pa = 1, paddy
          read(fpqout_unit, *) pqout_idx(pa), min_pmeshid(pa), min_dist(pa), orifice_num(pa), &
                         pdrain(pa), dr_dist(pa), device(pa), etp(pa)
    enddo
    read(fpqout_unit, *)
    read(fpqout_unit, *)
    read(fpqout_unit, *) nhp
    allocate(dhp(72000, nhp), phid(nhp))
    do ii = 1, nhp
          read(fpqout_unit, *) phid(ii)
    enddo
    endif
    close(fpqout_unit)

    ! paddy field parameters
    open(newunit = fpaddy_param_unit, file = fpaddy_param, action = 'read')
    read(fpaddy_param_unit, 1361) lh           ! 畦畔高
    read(fpaddy_param_unit, 1363) wh1          ! 分離型のセキ板の高さ
    read(fpaddy_param_unit, 1361) wh2          ! 一体型のセキ板の高さ
    read(fpaddy_param_unit, 1361) ww1          ! 落水口幅
    read(fpaddy_param_unit, 1361) ww2          ! 器具の切欠幅
    read(fpaddy_param_unit, 1362) ca           ! 器具の中心角
    read(fpaddy_param_unit, 1362) wtyp         ! セキのタイプ（1:四角セキ, 2:三角セキ）
    read(fpaddy_param_unit, 1361) dld          ! 畦畔天端と器具上端の高さの差
    read(fpaddy_param_unit, 1361) dd           ! 器具の穴の直径
    read(fpaddy_param_unit, 1361) unstdh       ! 田面から器具の穴中心までの高さ
    read(fpaddy_param_unit, 1361) p_data       ! 排水溝パイプ直径
    read(fpaddy_param_unit, 1361) ph           ! 田面からパイプ中心までの高さ
    close(fpaddy_param_unit)
1361 format(7x, f11.2)
1362 format(7x, i11)
1363 format(7x, f11.3)
!
    !allocate variables
    allocate(ttp(mesh))

end subroutine paddydat

!==============================
! Read Sewerage and Fields Are 
!==============================
!--------------------------------------------------------------------------
! UNST2D下水道・圃場読み込みサブルーチン/UNST2D drainage process read subroutine
! 変数/Variable:
! > グローバル変数/global variable
!   inf_dr: 施設分類固有id/facility classification unique id
!   dr_no: 施設の数/num of drainage facility
!--------------------------------------------------------------------------
subroutine draindat
    implicit none
    integer me, nn, drn
    integer finf_dr_unit, fdrain_unit

    ! inf_dr(sewerage and fields id) data read
    ! (0:not 1~:sewerage and fields are)
    open(newunit = finf_dr_unit, file = finf_dr, action = 'read')
    allocate(inf_dr(mesh))
    do me = 1, mesh
    read(finf_dr_unit, *) inf_dr(me)
    enddo
    close(finf_dr_unit)

    ! drainage system planning data read
    ! (keikaku haisui noryoku)
    open(newunit = fdrain_unit, file = fdrain, action = 'read')
    read(fdrain_unit, *) dr_no
    if(dr_no>0) then
    read(fdrain_unit, *)
    allocate(drp(dr_no), drc(dr_no))
    allocate(drr(dr_no), drr_dist(dr_no))
    allocate(vol_dr(mesh), vol(mesh))
    allocate(dhj(72000, dr_no))
    do nn = 1, dr_no
          read(fdrain_unit, *)  drn, drr(nn), drp(nn), drc(nn), drr_dist(nn)
          if (nn == dr_no) exit
          if (drr(nn) < 0.0) then
                print *, 'Enter the planned treatment rainfall.'
          endif
    enddo
    endif
    close(fdrain_unit)
    ! allocate variables
    allocate(tripTime(mesh))

end subroutine draindat

!======================
! Read Line Embankment
!======================
!--------------------------------------------------------------------------------
! UNST2D線盛土・通過率処理読み込みサブルーチン
!       /UNST2D line embankment process read subroutine
! 変数/Variable:
! > グローバル変数/global variable
!   mmorid: 盛り土対象linkの数/num of target link
!   nmorili: 対象link/target link
!   infl(link): 線盛土link判定フラグ/Line embankment link judgment flag
!   zbbk: 盛土(地盤)高/Line bank ground level
!   rbeta(link): 通過率(0.0~1.0)/passage rate
!--------------------------------------------------------------------------------
subroutine moriddat
    implicit none
    integer mmo
    integer fmorid_unit, ios
    character(len=50) line
    integer, allocatable :: nmorili(:)  ! change  v.1.0.5

    open(newunit = fmorid_unit, file = fmorid, action = 'read')
    read(fmorid_unit, *) mmorid
    read(fmorid_unit, *)

    allocate( nmorili(mmorid), infl(link) )
    allocate( zbbk(link) )

    infl = 0  ! initialize

    do mmo = 1, mmorid
        read(fmorid_unit, '(A)', iostat = ios) line
        if(ios/=0) exit
            read(line, *, iostat = ios) nmorili(mmo), zbbk(nmorili(mmo)), rbeta(nmorili(mmo))
            if(ios/=0) then
                read(line, *, iostat = ios) nmorili(mmo), zbbk(nmorili(mmo))
                rbeta(nmorili(mmo)) = 0.0d0
                if(ios/=0) stop 'UNST2D - ERROR-1004 : morido could not be loaded.'
            endif
        if(zbbk(nmorili(mmo)) > 0.0d0) infl(nmorili(mmo)) = 1
    enddo
    close(fmorid_unit)

    deallocate(nmorili)

end subroutine moriddat

!================
! Read Discharge  ! v.1.0.5 separate
!================
!----------------------------------------------------------------------------
! UNST2D排水処理読み込みサブルーチン/UNST2D discharge process read subroutine
! 入力/Input:
!   lasth: RRIの計算時間/calculation time(/hour)
! 変数/Variable:
! > 一時変数/temporary variable
!   idswlmax: 水位データの最大数/max num of water level data
!   dsme(dsmenum): 排水対象meshid/target meshid
!   dstype(dsmenum): 処理タイプ/process type
!                    (=1:自由流出/free flow, 
!                     =2:水位付与/water level specification,
!                     =3:等流/uniform flow  *only continuous rectangular mesh)
!   dswl_num(dsmenum): 水位データ数/num of water level data
!   sum_rtuv_x(mesh): ウェイト再計算用総和変数/
!                     Summation variable for recalculating weights
!   sum_rtuv_y(mesh): ウェイト再計算用総和変数/
!                     Summation variable for recalculating weights
! > グローバル変数/global variable
!   dsmenum: 排水処理対象mesh数/num of target mesh
!   ids2mesh: 水位付与meshの数/num of water level specification mesh
!   dsinf(mesh): = dstype(i) (Expanded to mesh number array)
!   ds2me(dsmenum): dsinf=2 meshid array(the array is allocated extra)
!   dsdt(mesh): 水位データタイムステップ(水位付与meshのみ値を持つ)/
!               water level data time step(*only dsinf=2 mesh has value)
!   dswl(mesh,idswlmax): 水位データ(水位付与meshのみ値を持つ)/
!                        water level data(*only dsinf=2 mesh has value)
!   dsupme(mesh): 上流meshid(等流meshのみ値を持つ)/
!                 upstream meshid(*only dsinf=3 mesh has value)
!----------------------------------------------------------------------------
subroutine dsmeshdat
    implicit none
    integer i, j, k, ios, read_count
    integer fdsmesh_unit
    integer idswlmax
    character(len=100) :: line
    integer, allocatable :: dsme(:), dstype(:), dswl_num(:)

    ! downstream data read
    open(newunit = fdsmesh_unit, file = fdsmesh, action = 'read')
    read(fdsmesh_unit, *) dsmenum
    allocate( dsme(dsmenum), dstype(dsmenum), dswl_num(dsmenum))  ! temporary
    allocate( ds_dt(mesh), ds_inf(mesh), ds2me(dsmenum), ds_upme(mesh))  ! global
    ! apply initial v
    dstype = 1  ! common
    ds_inf = 0   ! common
    idswlmax = 0  ! dsinf=2
    ids2mesh = 0  ! dsinf=2
    ds_dt = timmax  ! dsinf=2
    dswl_num = 0  ! dsinf=2
    ds2me = 0  ! dsinf=2
    ds_upme = 0  ! dsinf=3
    read(fdsmesh_unit, *)
    do i = 1, dsmenum
        read(fdsmesh_unit, '(A)', iostat = ios) line
        if (ios /= 0) exit

        ! count data
        line = trim(adjustl(line))
        read_count = 0
        do j = 1, len_trim(line)
            if (line(j:j) == ' ' .and. line(j+1:j+1) /= ' ') then
                read_count = read_count + 1
            endif
        enddo

        if (read_count == 0) then
            read(line, *) dsme(i)
        elseif (read_count == 1) then
            read(line, *) dsme(i), dstype(i)
        elseif (read_count == 2) then
            read(line, *) dsme(i), dstype(i), dswl_num(i)
        else
            stop 'UNST2D - ERROR-1003 : dsmesh could not be loaded.'
        endif

        ds_inf(dsme(i)) = dstype(i)
        idswlmax = max(dswl_num(i), idswlmax, 0)
        ! count dstype2 & cal wl time step
        if(dstype(i)==2) then
            ids2mesh = ids2mesh + 1
            ds2me(ids2mesh) = dsme(i)
            if (dswl_num(i)<=0) then
                write(*,*) ' UNST2D - WARNIG : meshid' , dsme(i) , 'is illegal value, No discharge.'
                dstype(i) = 0
                ds_inf(dsme(i)) = 0
                ids2mesh = ids2mesh - 1
            else
                ds_dt(dsme(i)) = timmax/dswl_num(i)
            endif
        endif
    enddo

    ! dstype2: water level(m)
    if (any(dstype == 2)) then
        read(fdsmesh_unit, *)
        allocate(ds_wl(mesh,idswlmax+1))

        ds_wl = 0.0d0  ! apply initial v

        ! downstream mesh rtuv set
        allocate(sum_rtuv_x(mesh), sum_rtuv_y(mesh))
    endif

    do i = 1, dsmenum
        ! dstype1: free outflow mesh set
        if (dstype(i)==1) then
            baseo(dsme(i)) = -9999.0d0

        ! dstype2: water level(m)
        elseif (dstype(i)==2) then
            ! read water level(m)
            read(fdsmesh_unit, *) (ds_wl(dsme(i),j),j=1,dswl_num(i))
            ds_wl(dsme(i),dswl_num(i)+1) = ds_wl(dsme(i),dswl_num(i))
            ! cal downstream rtuv set
            sum_rtuv_x(dsme(i)) = 0.0d0
            sum_rtuv_y(dsme(i)) = 0.0d0
            do k = 1, ko(dsme(i))
                if (limesh(melink(dsme(i),k),2)/=0) then
                    sum_rtuv_x(dsme(i)) = sum_rtuv_x(dsme(i)) + rtuv_x(dsme(i),k)
                    sum_rtuv_y(dsme(i)) = sum_rtuv_y(dsme(i)) + rtuv_y(dsme(i),k)
                endif
            enddo
            do k = 1, ko(dsme(i))
                if (limesh(melink(dsme(i),k),2)/=0) then
                    rtuv_x(dsme(i),k) = rtuv_x(dsme(i), k)/sum_rtuv_x(dsme(i))
                    rtuv_y(dsme(i),k) = rtuv_y(dsme(i), k)/sum_rtuv_y(dsme(i))
                else
                    rtuv_x(dsme(i),k) = 0.0d0
                    rtuv_y(dsme(i),k) = 0.0d0
                endif
            enddo
            
        ! dstype3: uniform flow(beta)
        !! only square mesh -- copy upper mesh env. !!
        elseif (dstype(i)==3) then
            if (ko(dsme(i))==4) then
                do k = 1, ko(dsme(i))
                    ! set upper meshid
                    if (limesh(melink(dsme(i), k),2) == 0) then
                        if (limesh(melink(dsme(i), mod(k+1,4)+1),1) == dsme(i)) then
                            ds_upme(dsme(i)) = limesh(melink(dsme(i), mod(k+1,4)+1),2)
                        else
                            ds_upme(dsme(i)) = limesh(melink(dsme(i), mod(k+1,4)+1),1)
                        endif
                        !write(*,*) dsme(i) , '-' , dsupper(dsme(i))
                    endif
                enddo
            else
                write(*,*) ' UNST2D - WARINIG : meshid', dsme(i), 'is not square, No discharge.'
                ds_inf(dsme(i)) = 0
            endif
        endif
    enddo
    close(fdsmesh_unit)
    if (any(dstype == 2)) deallocate(sum_rtuv_x, sum_rtuv_y)
    deallocate(dsme, dstype, dswl_num)

end subroutine dsmeshdat
    
end module unst_read
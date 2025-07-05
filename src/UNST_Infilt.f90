! infiltration
subroutine unst_infilt(me)
    use globals
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
    
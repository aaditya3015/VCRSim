module fluidprop
    use syswide
    implicit none
    private
    public fluid_state
    interface
    end interface
    type fluid_state
        real(dp) :: p, h, x, t, d, s
        real(dp) :: hl, hv, rhol, rhov
        real(dp) :: cpl, cpv, kl, kv, mul, muv
    end type fluid_state
contains
    subroutine propph(p, h, state)
        real(dp), intent(in) :: p, h
        type(fluid_state), intent(out) :: state
        state%p = p; state%h = h
    end subroutine propph
end module fluidprop

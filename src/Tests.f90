module tests
    use syswide
    use moistaircalc
    use solvers
    implicit none
    private
    public test_moistaircalc
contains
    subroutine test_moistaircalc
        real(dp) :: sol
        print *, "Psat water @ 273.15 K calc, actual", psat_water(273.15_dp), 0.61115_dp
        print *, "Psat water @ 323.15 K calc, actual", psat_water(323.15_dp), 12.3513_dp
        print *, "Psat water @ 373.15 K calc, actual", psat_water(373.15_dp), 101.4180_dp
        print *, "dew point, calc vs actual", dewpoint(101.325_dp, 300.0_dp, 1, 1.0_dp), 300.0_dp
        sol = illinois1(samplefunc, 0.0_dp, 1.0_dp, 10e-5_dp)
        print *, sol, samplefunc(sol)
    end subroutine test_moistaircalc

    real(dp) function samplefunc(x) result(y)
        real(dp), intent(in) :: x
        !y = x*x*x - 27.0_dp
        y = x**4 - 2.0_dp*x**2 + 0.25_dp
    end function samplefunc
end module tests

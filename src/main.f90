
program hello
    use syswide
    use tests
    use solvers, only: illinois1, bisection, fzero, err, niter, iprint
    use auto_diff
    use moistaircalc, only: psat_water, psat_water_ad, dewpoint
    use fluidprop
    use Nelder_Mead, only: minim
    use Powell_Optimize, only: uobyqa
    implicit none
    !call test_moistaircalc
    real(dp) :: a, b, tol, y
    real(dp) :: step(2), p(2), var(2)
    integer(int) :: ifault
    iprint = 0
    a = 200.0_dp
    b = 380.0_dp
    tol = 1.0e-3_dp
    step = (/0.5_dp, 0.5_dp/)
    p = (/1.0_dp, 1.0_dp/)

    print *, "illinois1 algoritm - ", illinois1(func, a, b, tol), err, niter
    print *, "bisection algoritm - ", bisection(func, a, b, tol), err, niter
    print *, "dew point", dewpoint(100.325_dp, 313.0_dp, 1, 1.0_dp), err, niter
    !print *, "fzero method - ", fzero(func, a, b, tol), err, niter

    !call test_moistaircalc
    !sol1 = newton_ad(ad_func, 283.15_dp, 1.0e-3_dp)

    !minim(p, step, nop, func, maxfn, iprint, stopcr, nloop, iquad, &
    !simp, var, functn, ifault)

    !call minim(p, step, 2, y, 100, 1, 1.0e-3_dp, 6, 1, &
    !           1.0e-3_dp, var, functn, ifault)
    !print *, p
    !p = (/5.0_dp, 1.0_dp/)
!uobyqa(n, x, rhobeg, rhoend, iprint, maxfun)
    !call uobyqa(2, p, 0.5_dp, 1e-3_dp, 1, 100, calfun)
    !print *, p

contains
    real(dp) function func(x) result(y)
        real(dp), intent(in) :: x
        y = psat_water(x) - 20.0_dp
        !y = (x**(2.0_dp)) - 100.0_dp
    end function func
    function ad_func(x) result(y)
        type(auto_var) :: x, y
        y = psat_water_ad(x) - 8.0_dp
    end function ad_func
    SUBROUTINE functn(p, func)
        !use syswide
        IMPLICIT NONE
        !INTEGER, PARAMETER     :: dp = kind(0.0_d0)!SELECTED_REAL_KIND(12, 60)
        REAL(dp), INTENT(IN)  :: p(:)
        REAL(dp), INTENT(OUT) :: func
        func = (p(1) + 2.0_dp*p(2) - 7.0_dp)**2 + (2.0_dp*p(1) + p(2) - 5.0_dp)**2
    END SUBROUTINE functn
    SUBROUTINE calfun(n, x, f)
        use syswide
        IMPLICIT NONE
        !INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
        INTEGER, INTENT(IN)    :: n
        REAL(dp), INTENT(IN)  :: x(:)
        REAL(dp), INTENT(OUT) :: f
        f = (x(1) + 2.0_dp*x(2) - 7.0_dp)**2 + (2.0_dp*x(1) + x(2) - 5.0_dp)**2
    END SUBROUTINE calfun
end program

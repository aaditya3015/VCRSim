module solvers
    use syswide
    use auto_diff
    implicit none
    real(dp) :: err
    integer(int) :: niter, iprint
    private
    public brent, bisection, newton_ad, illinois, regfalsi, &
        illinois1, fzero, err, niter, iprint

    interface
        real(dp) function func(x)
            use syswide
            implicit none
            real(dp), intent(in) :: x
        end function
        function ad_func(x) result(y)
            use auto_diff
            type(auto_var) :: x, y
        end function ad_func
    end interface
contains

    real(dp) function bisection(f, aa, bb, abstol) result(sol)
        procedure(func) :: f
        real(dp), intent(in) :: aa, bb, abstol
        real(dp) :: a, fa, b, fb, c, fc
        integer(int), parameter :: maxiter = 30
        integer(int) :: i
        a = aa; b = bb  ! making local copies
        fa = f(a); fb = f(b)
        if (fa*fb .gt. 0.0_dp) then
            print *, "roots are not bracketed"
            sol = 0.0_dp
            return
        end if
        err = 1.0_dp
        i = 0
        solloop: do while (err .gt. abstol .and. i .lt. maxiter)
            err = abs(b - a)
            c = (a + b)/2.0_dp
            fc = f(c)
            if (fc*fb .gt. 0.0_dp) then
                b = c
                fb = fc
            else
                a = c
                fa = fc
            end if
            i = i + 1
            !print *, c, i
        end do solloop
        niter = i
        sol = c
    end function bisection

    real(dp) function illinois(f, aa, bb, abstol) result(sol)
        procedure(func) :: f
        real(dp), intent(in) :: aa, bb, abstol
        real(dp) :: a, fa, b, fb, c, fc
        logical :: pside
        integer(int), parameter :: maxiter = 30
        integer(int) :: i
        a = aa; b = bb  ! making local copies
        fa = f(a); fb = f(b)
        if (fa*fb .gt. 0.0_dp) then
            print *, "roots are not bracketed"
            sol = 0.0_dp
            return
        end if
        err = 1.0_dp
        pside = .true.
        i = 0
        solloop: do while (err .gt. abstol .and. i .lt. maxiter)
            err = abs(b - a)
            c = (fa*b - fb*a)/(fa - fb)
            fc = f(c)
            if (fb*fc .gt. 0.0_dp) then
                b = c; fb = fc
                if (.not. (pside)) then
                    fa = fa/2.0_dp
                    pside = .false.
                end if
            else !if (fa*fc .gt. 0.0_dp) then
                a = c; fa = fc
                if (pside) then
                    fb = fb/2.0_dp
                    pside = .true.
                end if
                !else
                !c = (a + b)/2.0_dp
            end if
            i = i + 1
            !print *, c, i, err
        end do solloop
        niter = i
        sol = c
    end function illinois

    real(dp) function illinois1(f, aa, bb, abstol) result(sol)
        procedure(func) :: f
        real(dp), intent(in) :: aa, bb, abstol
        real(dp) :: a, fa, b, fb, c, fc, reltol, slp
        integer(int), parameter :: maxiter = 30
        integer(int) :: i
        a = aa; b = bb  ! making local copies
        fa = f(a); fb = f(b)
        fa = fa; fb = fb
        if (fa*fb .gt. 0.0_dp) then
            print *, "roots are not bracketed"
            sol = 0.0_dp
            return
        end if
        reltol = 1.0e3_dp
        i = 0
        err = 1.0_dp
        solloop: do while (err .gt. abstol .and. i .lt. maxiter)
            err = abs(b - a)
            slp = abs((fb - fa)/(b - a))
            c = (fa*b - fb*a)/(fa - fb)
            !if (abs(c - a)/abs(b - a) .lt. reltol .or. abs(c - b)/abs(b - a) .lt. reltol) then
            if (slp .gt. reltol) then
                c = (a + b)/2.0_dp
                fc = f(c)
                if (fa*fc .lt. 0.0_dp) then
                    b = a; fb = fa
                end if
                a = c; fa = fc
            else
                fc = f(c)
                if (fa*fc .lt. 0.0_dp) then
                    b = a; fb = fa
                else
                    fb = 0.5_dp*fb
                end if
                a = c; fa = fc
            end if
            i = i + 1
            if (iprint .eq. 1) then
                print *, c, i, err, slp
            end if
        end do solloop
        niter = i
        sol = c
    end function illinois1

    real(dp) function regfalsi(f, aa, bb, abstol) result(sol)
        procedure(func) :: f
        real(dp), intent(in) :: aa, bb, abstol
        real(dp) :: a, fa, b, fb, c, fc, cprv, fcprv
        integer(int), parameter :: maxiter = 30
        integer(int) :: i
        real(dp) :: mm = 0.5_dp, nn = 0.5_dp
        a = aa; b = bb  ! making local copies
        fa = f(a); fb = f(b)
        if (fa*fb .gt. 0.0_dp) then
            print *, "roots are not bracketed"
            sol = 0.0_dp
            return
        end if
        err = 1.0_dp
        !c = (fb*a - fa*b)/(fb - fa)
        c = (a + b)/2.0_dp
        fc = f(c)
        i = 0
        solloop: do while (err .gt. abstol .and. i .lt. maxiter)
            err = abs(b - a)
            cprv = (fb*a - fa*b)/(fb - fa)
            fcprv = f(cprv)
            if (fc*fcprv .lt. 0.0_dp) then
                c = cprv
                fc = fcprv
            elseif (fcprv*fc .gt. 0.0_dp .and. fcprv*fb .gt. 0.0_dp) then
                c = (a*fb - nn*b*fa)/(fb - nn*fa)
                fc = f(c)
            elseif (fcprv*fc .gt. 0.0_dp .and. fcprv*fb .lt. 0.0_dp) then
                c = (mm*a*fb - b*fa)/(mm*fb - fa)
                fc = f(c)
            end if
            if (fc*fb .gt. 0.0_dp) then
                b = c
                fb = fc
            else
                a = c
                fa = fc
            end if
            i = i + 1
            !print *, c, i
        end do solloop
        niter = i
        sol = c
    end function regfalsi

    real(dp) function brent(f, aa, bb, abstol) result(sol)
        procedure(func) :: f
        real(dp), intent(in) :: aa, bb, abstol
        real(dp) :: a, fa, b, fb, c, fc, d, s, fs, tol
        logical :: mflag
        integer(int), parameter :: maxiter = 30
        integer(int) :: i
        !real(dp) :: abstol = 1e-3_dp
        tol = epsilon(1.0_dp)
        a = aa; b = bb  ! making local copies
        fa = f(a); fb = f(b)
        if (fa*fb .gt. 0.0_dp) then
            print *, "roots are not bracketed"
            sol = 0.0_dp
            return
        end if
        if (abs(fa) .lt. abs(fb)) then
            call swap(a, b)
            call swap(fb, fa)
        end if
        c = a
        fc = fa
        s = c
        fs = fc
        mflag = .true.
        i = 0
        err = 1.0_dp
        solloop: do while (abs(fs) .lt. 1e-12_dp .or. (err .gt. abstol .and. i .lt. maxiter))
            err = abs(b - a)
            !only inverse quadratic approximation
            if (abs(fa - fc) .gt. tol .and. abs(fb - fc) .gt. tol) then
                s = a*fb*fc/((fa - fb)*(fa - fc)) + b*fa*fc/((fb - fa) - (fb - fc)) + &
                    c*fa*fb/((fc - fa) - (fc - fb))
                !print *, "used quadratic"
            else
                s = (a*fb - b*fa)/(fb - fa)
            end if
            if ((s - b)*(s - 3.0_dp*a + b)/4.0_dp .gt. 0.0_dp .or. &
                (mflag .and. abs(s - b) .ge. abs(b - c)/2.0_dp) .or. &
                (.not. (mflag) .and. abs(s - b) .ge. abs(c - d)/2.0_dp) .or. &
                (mflag .and. abs(b - c) .lt. abstol) .or. &
                (.not. (mflag) .and. abs(c - d) .lt. abstol)) then
                s = (a + b)/2.0_dp
                mflag = .true.
            else
                mflag = .false.
            end if
            fs = f(s)
            d = c
            c = b
            if (fa*fs .lt. 0) then
                b = s
                fb = fs
            else
                a = s
                fa = fs
            end if
            if (abs(fa) .le. abs(fb)) then
                call swap(a, b)
                call swap(fa, fb)
            end if
            i = i + 1
            !print *, s, i
        end do solloop
        niter = i
        sol = s
    end function brent

    real(dp) function newton_ad(f, xguess, abstol) result(sol)
        procedure(ad_func) :: f
        real(dp), intent(in) :: xguess, abstol
        real(dp):: xold, xnew
        integer(int) :: i, imax = 20
        type(auto_var) :: x, y
        !call x%set(xguess)
        !y = f(x)
        err = 1.0_dp
        i = 1
        xnew = xguess
        do while ((abs(err) .gt. abstol) .or. (i .ge. imax))
            xold = xnew
            call x%set(xold)
            y = f(x)
            xnew = xold - y%get_value()/y%get_derivative()
            err = xnew - xold
            print *, xnew
        end do
        niter = i
        sol = xnew
    end function newton_ad

    real(dp) function fzero(f, aa, bb, abstol) result(sol)
        procedure(func) :: f
        real(dp), intent(in) :: aa, bb, abstol
        real(dp) :: tsol
        tsol = illinois1(f, aa, bb, abstol)
        !if (f(tsol) .gt. abstol .or. tsol .lt. abstol) then
        !    tsol = illinois1(f, aa, bb, abstol)
        !end if
        sol = tsol
    end function fzero

    subroutine swap(a, b)
        real(dp), intent(inout) :: a, b
        real(dp) :: c
        c = a
        a = b
        b = c
    end subroutine swap
end module solvers

module moistaircalc
    use syswide
    use solvers, only: illinois1, fzero
    implicit none
    private
    public psat_water, density_moistair, humidity_ratio, dewpoint, enthalpy_moistair, &
        moistairprop, psat_water_ad
    !interface moistair
    !procedure density_moistair
    !end interface moistair
contains
    function psat_water(T) result(pws)
        ! this function calculates the saturation vapour pressure of the water at a given temperature
        ! See ashrae handbook - fundamentals - chapter1 for reference
        real(dp) :: T, pws
        logical :: tneg, tpos
        real(dp), dimension(13) :: CF
        CF(1) = -5.6745359e3_dp
        CF(2) = 6.3925247_dp
        CF(3) = -9.6778430e-3_dp
        CF(4) = 6.2215701e-7_dp
        CF(5) = 2.0747825e-9_dp
        CF(6) = -9.4840240e-13_dp
        CF(7) = 4.1635019_dp
        CF(8) = -5.8002206e3_dp
        CF(9) = 1.3914993_dp
        CF(10) = -4.8640239e-2_dp
        CF(11) = 4.1764768e-5_dp
        CF(12) = -1.4452093e-8_dp
        CF(13) = 6.5459673_dp

        tneg = (T .gt. 173.15_dp .and. T .lt. 273.15_dp)
        tpos = (T .ge. 273.15_dp .and. T .lt. 473.15_dp)

        if (tneg) then
            pws = dexp(CF(1)/T + CF(2) + CF(3)*T + CF(4)*T**(2) + CF(5)*T**(3) + CF(6)*T**(4) + CF(7)*dlog(T))
        elseif (tpos) then
            pws = dexp(CF(8)/T + CF(9) + CF(10)*T + CF(11)*T**(2) + CF(12)*T**(3) + CF(13)*dlog(T))
        else
            pws = 0.0_dp*T
        end if
        pws = pws/1000.0_dp; 
    end function psat_water

    function psat_water_ad(T) result(pws)
        ! this function calculates the saturation vapour pressure of the water at a given temperature
        ! See ashrae handbook - fundamentals - chapter1 for reference
        use auto_diff
        type(auto_var) :: T, pws
        logical :: tneg, tpos
        real(dp), dimension(13) :: CF
        CF(1) = -5.6745359e3_dp
        CF(2) = 6.3925247_dp
        CF(3) = -9.6778430e-3_dp
        CF(4) = 6.2215701e-7_dp
        CF(5) = 2.0747825e-9_dp
        CF(6) = -9.4840240e-13_dp
        CF(7) = 4.1635019_dp
        CF(8) = -5.8002206e3_dp
        CF(9) = 1.3914993_dp
        CF(10) = -4.8640239e-2_dp
        CF(11) = 4.1764768e-5_dp
        CF(12) = -1.4452093e-8_dp
        CF(13) = 6.5459673_dp

        tneg = (T%get_value() .gt. 173.15_dp .and. T%get_value() .lt. 273.15_dp)
        tpos = (T%get_value() .ge. 273.15_dp .and. T%get_value() .lt. 473.15_dp)

        if (tneg) then
            pws = exp(CF(1)/T + CF(2) + CF(3)*T + CF(4)*T**(2) + CF(5)*T**(3) + CF(6)*T**(4) + CF(7)*log(T))
        elseif (tpos) then
            pws = exp(CF(8)/T + CF(9) + CF(10)*T + CF(11)*T**(2) + CF(12)*T**(3) + CF(13)*log(T))
        else
            pws = 0.0_dp*T
        end if
        pws = pws/1000.0_dp; 
    end function psat_water_ad

    real(dp) function density_moistair(p, t, iden, val) result(rho)
        !this function evaluates the dewpoint temperature of the most air
        ! p -> Total Pressure in kPa
        ! t-> dry bulb temperature in K
        ! identifier (iden) -> 1 for relative humidity
        !                   -> 2 for humidity humidity_ratio
        !                   -> 3 for wet bulbt temperature
        !                   -> 4 for dew point temparature
        ! val-> value of the identifier
        real(dp), intent(in) :: p
        real(dp), intent(in) :: t
        integer(int), intent(in) :: iden
        real(dp), intent(in) :: val
        real(dp) :: w, v
        w = humidity_ratio(p, t, iden, val)
        v = 0.287042_dp*(t)*(1.0_dp + 1.607858_dp*w)/p; 
        rho = 1.0_dp/v; 
    end function density_moistair

    real(dp) function humidity_ratio(p, t, iden, val) result(w)
        !this functions calculates the humidity ratio of the moistair
        ! p -> Total Pressure in kPa
        ! t-> dry bulb temperature in K
        ! identifier (iden) -> 1 for relative humidity
        !                   -> 2 for humidity_ratio
        !                   -> 3 for wet bulbt temperature
        !                   -> 4 for dew point temparature
        ! val -> value of the idening kg_vap/kg_dryair
        real(dp), intent(in) :: p, t, val
        integer(int), intent(in) :: iden
        real(dp) :: t0, rh, pws, pws_wbt, wbt, ws_wbt, dpt, w_const
        t0 = 273.15_dp
        w_const = 0.621945_dp
        select case (iden)
        case (1)
            rh = val
            pws = psat_water(t)
            w = w_const*rh*(pws)/(p - rh*pws)
        case (2)
            w = val
        case (3)
            wbt = val
            pws_wbt = psat_water(wbt)
            ws_wbt = w_const*(pws_wbt)/(p - pws_wbt)
            if (t >= t0) then
                w = ((2501_dp - 2.326_dp*(wbt - t0))*ws_wbt - &
                     1.006_dp*((t - t0) - (wbt - t0)))/ &
                    (2501.0_dp + 1.86_dp*(t - t0) - 4.186_dp*(wbt - t0))
            else
                w = ((2830_dp - 0.24_dp*(wbt - t0))*ws_wbt - &
                     1.006_dp*((t - t0) - (wbt - t0)))/ &
                    (2830.0_dp + 1.86_dp*(t - t0) - 2.1_dp*(wbt - t0))
            end if
        case (4)
            dpt = val
            pws = psat_water(dpt)
            w = w_const*(pws)/(p - pws)
        end select
    end function humidity_ratio

    real(dp) function dewpoint(p, t, iden, val) result(t_dpt)
        !this function evaluates the dewpoint temperature of the most air
        ! p -> Total Pressure in kPa
        ! t-> dry bulb temperature in K
        ! identifier (iden) -> 1 for relative humidity
        !                   -> 2 for humidity humidity_ratio
        !                   -> 3 for wet bulb temperature
        !                   -> 4 for dew point temparature
        ! val -> value of the idening kg_vap/kg_dryair
        implicit none
        real(dp), intent(in) :: p, t, val
        integer(int), intent(in) :: iden
        real(dp) :: w, p_w
        select case (iden)
        case (2)
            w = val
        case (1, 3)
            w = humidity_ratio(p, t, iden, val)
        case (4)
            t_dpt = val
        end select
        p_w = p/(0.621945_dp/w + 1.0_dp); 
        if (iden /= 4) then
            !t_dpt = 0.0_dp!fzero(@(temp) (psat_water(temp) - p_w),[t-100,t]);
            t_dpt = illinois1(dummy, t - 100, t, 1e-3_dp)
            !t_dpt = newton(ad_func, )
        end if
    contains
        real(dp) function dummy(temp) result(err)
            real(dp), intent(in) :: temp
            err = psat_water(temp) - p_w
        end function dummy
    end function dewpoint

    real(dp) function enthalpy_moistair(p, t, iden, val) result(enth)
        !this function evaluates the enthalpy of the most air
        ! p -> Total Pressure in kPa
        ! t-> dry bulb temperature in K
        ! identifier (iden) -> 1 for relative humidity
        !                   -> 2 for humidity humidity_ratio
        !                   -> 3 for wet bulb temperature
        !                   -> 4 for dew point temparature
        ! val -> value of the idening kg_vap/kg_dryair
        real(dp), intent(in) :: p, t, val
        integer(int), intent(in) :: iden
        real(dp), parameter :: t0 = 273.15_dp
        enth = 1.006_dp*(t - t0) + humidity_ratio(p, t, iden, val)*(2501.0_dp + 1.86_dp*(t - t0))
    end function enthalpy_moistair

    subroutine moistairprop(p, t, iden, val, w, wbt, dpt, enth)
        real(dp), intent(in) :: p, t, val
        real(dp), intent(out) :: w, wbt, dpt, enth
        integer(int), intent(in) :: iden
        real(dp) :: p_w
        w = humidity_ratio(p, t, iden, val)
        p_w = p/(0.621945_dp/w + 1.0_dp); 
        if (iden /= 4) then
            !t_dpt = 0.0_dp!fzero(@(temp) (psat_water(temp) - p_w),[t-100,t]);
            dpt = fzero(dummy, t - 100, t, 1e-3_dp)
        else
            dpt = val
        end if
        wbt = 0.0_dp
        enth = 0.0_dp
    contains
        real(dp) function dummy(temp) result(err)
            real(dp), intent(in) :: temp
            err = psat_water(temp) - p_w
        end function dummy
    end subroutine moistairprop
end module moistaircalc

module TPcorrelations
    use syswide
    implicit none
    private
    real(dp), parameter:: grav = 9.8
    public friedel, lockhart_martinelli, msh, gronnerud
contains

    real(dp) function friedel(x, G, D, dl, dv, mul, muv, sigma) result(dpdz)
        real(dp), intent(in) :: x, G, D, dl, dv, mul, muv, sigma
        real(dp) :: fl, fv, Rel, Rev
        Rel = G*D/mul; Rev = G*D/muv
        fl = 0.079_dp/(Rel**0.25_dp)
        fv = 0.079_dp/(Rev**0.25_dp)
        dpdz = 1.0_dp
    end function friedel

    real(dp) function lockhart_martinelli(x, G, D, dl, dv, mul, muv) result(dpdz)
        real(dp), intent(in) :: x, G, D, dl, dv, mul, muv
        real(dp) :: Rel, Rev, fl, fv, xtt, dpdzl, dpdzv, phi_ltt, phi_vtt
        real(dp) :: C = 20.0_dp
        Rel = G*D/mul; Rev = G*D/muv
        fl = 0.079_dp/(Rel**0.25_dp)
        fv = 0.079_dp/(Rev**0.25_dp)
        xtt = (1.0_dp/x - 1.0_dp)**(0.9_dp)*(dv/dl)**(0.5_dp)*(mul/muv)**(0.1_dp)
        if (Rel .gt. 4000.0_dp) then
            dpdzl = 2.0_dp*fl*G*G/(dl*D)
            phi_ltt = 1.0_dp + C/(xtt) + 1.0_dp/(xtt*xtt)
            dpdz = phi_ltt*phi_ltt*(dpdzl)
        else
            dpdzv = 2.0_dp*fv*G*G/(dv*D)
            phi_vtt = 1.0_dp + C*xtt + xtt*xtt
            dpdz = phi_vtt*phi_vtt*(dpdzv)
        end if
    end function lockhart_martinelli

    real(dp) function msh(x, G, D, dl, dv, mul, muv) result(dpdz)
        real(dp), intent(in) :: x, G, D, dl, dv, mul, muv
        real(dp) :: Rel, Rev, fl, fv, dpdzl, dpdzv, delta
        Rel = G*D/mul; Rev = G*D/muv
        fl = 0.079_dp/(Rel**0.25_dp)
        fv = 0.079_dp/(Rev**0.25_dp)
        dpdzl = 2.0_dp*fl*G*G/(dl*D)
        dpdzv = 2.0_dp*fv*G*G/(dv*D)
        delta = dpdzl + 2.0_dp*(dpdzv - dpdzl)*x
        dpdz = delta*(1.0_dp - x)**(1.0_dp/3.0_dp) + dpdzv*x**3.0_dp
    end function msh

    real(dp) function gronnerud(x, G, D, dl, dv, mul, muv) result(dpdz)
        real(dp), intent(in) :: x, G, D, dl, dv, mul, muv
        real(dp) :: Rel, fl, dpdzl, frl, ffr, dpdzfr, phi_l
        Rel = G*D/mul
        fl = 0.079_dp/(Rel**0.25_dp)
        dpdzl = 2.0_dp*fl*G*G/(dl*D)
        frl = G*G/(grav*D*dl*dl)
        ffr = frl**0.3_dp + 0.0055_dp*(dlog(1/frl))**2.0_dp
        dpdzfr = ffr*(x + 4.0_dp*(x**1.8_dp - x**10.0_dp*ffr**0.5_dp))
        phi_l = 1.0_dp + dpdzfr*(dl/dv*(muv/mul)**0.25_dp - 1.0_dp)
        dpdz = phi_l*dpdzl
    end function gronnerud

end module TPcorrelations

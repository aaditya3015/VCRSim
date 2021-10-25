module SPcorrelations
    use syswide
    implicit none
    private
    public dittus_boelter, sieder_tate, churchill
contains

    real(dp) pure function churchill(Re, e) result(ff)
        real(dp), intent(in) :: Re
        real(dp), intent(in), optional :: e
        real(dp) :: ee, A, B
        if (present(e)) then
            ee = e
        else
            ee = 1.52e-6_dp
            !print *, "assumed value of absolute roughness = 1.52 micron"
        end if
        A = (2.457_dp*log(1.0_dp/((7.0_dp/Re)**0.9_dp + (0.27_dp*ee/Re))))**16
        B = (37530.0_dp/Re)**16
        ff = 8.0_dp*((8.0_dp/Re)**12.0_dp + 1.0_dp/(A + B)**1.5_dp)**(1.0_dp/12.0_dp)
    end function churchill

    real(dp) pure function dittus_boelter(Re, Pr, iheat) result(Nu)
        real(dp), intent(in) :: Re, Pr
        integer(int), intent(in), optional :: iheat
        if (present(iheat)) then
            if (iheat .eq. 0) then
                Nu = 0.023_dp*Re**0.8_dp*Pr**0.4_dp
            else
                Nu = 0.023_dp*Re**0.8_dp*Pr**0.3_dp
            end if
        else
            Nu = 0.023_dp*Re**0.8_dp*Pr**0.4_dp
        end if
    end function dittus_boelter

    real(dp) function sieder_tate(Re, Pr) result(Nu)
        real(dp), intent(in) :: Re, Pr
        real(dp) :: fs
        fs = 1.0_dp/(1.58_dp*dlog(Re) - 3.28_dp)
        Nu = fs/2.0_dp*(Re - 1000.0_dp)*Pr/ &
             (1.0_dp + 12.7_dp*(Pr**(2.0_dp/3.0_dp) - 1.0_dp)*dsqrt(fs/2.0_dp))
    end function sieder_tate
end module SPcorrelations

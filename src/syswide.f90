module syswide
!module to apply properties system wide
    use, intrinsic :: ISO_C_BINDING
    implicit none
    !integer, parameter :: autodiffenabled = 1
    integer, parameter:: dp = kind(0.0d0) !(C_DOUBLE) !double precision
    integer, parameter:: int = kind(1) !integer
end module syswide

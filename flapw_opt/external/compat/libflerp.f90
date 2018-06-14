! libflerp (library FLEur wraP) allows more convenient access to some functions of Fleur.

module constants
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
end module

subroutine testfunc(a)
    integer a
    a = 3*a
end subroutine
subroutine array_test(a, n_a)
    use constants
    implicit none
    integer :: n_a
    real(dp) :: a(n_a,2)
    a = 2*a
end subroutine
subroutine array_out(a, n_a)
    use constants
    implicit none
    integer :: n_a, i
    real(dp) :: a(n_a, 2)
    do i = 1,n_a
        a(i, 1:2) = (/2*i, 2*i +1 /)
    end do
end subroutine

subroutine wrap_bessel(lmax, x, jl, djl)
    use m_sphbes
    use m_dsphbs
    implicit none
    integer, intent(in) :: lmax
    real, intent(in) :: x
    real, intent(out), dimension(0:lmax) :: jl, djl
    ! return both bessel j_l and its derivative at position x
    call sphbes(lmax, x, jl)
    call dsphbs(lmax, x, jl, djl)
end subroutine
    
subroutine wrap_ylm(lmax, r, Ylm)
    use m_ylm
    implicit none
    integer, intent(in) :: lmax
    real, intent(in) :: r(3)
    complex, intent(out) :: Ylm((lmax+1)**2)
    call ylm4(lmax, r, Ylm)
end subroutine

subroutine wrap_radfun( &
    l, E, dx, r_MT, &
    v_radmesh, n_radpts, &
    sol, sol_dE, &
    u_R, u_dr_R, u_dE_R, u_dr_dE_R, u_dE_norm)
    use constants
    use m_radfun
    implicit none
    integer, intent(in) :: l, n_radpts
    real(dp), intent(in) :: E, dx, r_MT
    real(dp), intent(in), dimension(n_radpts) :: v_radmesh
    real(dp), intent(out), dimension(n_radpts, 2) :: sol, sol_dE
    real(dp), intent(out) :: u_R, u_dr_R, u_dE_R, u_dr_dE_R, u_dE_norm
    ! dummy arguments
    real(dp) :: wronk
    integer :: node_u, node_d
    ! locals
    real(dp) :: r0

    r0 = r_MT * exp(dx * (1 - n_radpts)) ! first mesh point
    !sol = 1337.
    !sol_dE = 3117.
    !u_R = 10.
    !u_dr_R = 11.

    call radfun(l, &  ! radial quantum number
                E, &  ! energy
                v_radmesh, & ! potential on radial mesh
                n_radpts, &  ! number of radial grid points
                r0, & ! first point of radial mesh
                dx, & ! exp(dx) is the spacing between points
                n_radpts, & 
                ! jmtd == max(jri(:)), "relict" from fixed-size days
                sol, sol_dE, & ! solution & energy derivative on mesh
                ! function values at r_MT and two types of derivatives
                ! .._dr is the radial derivative
                ! .._dE is the energy derivative
                u_R, & 
                u_dr_R, &
                u_dE_R, &
                u_dr_dE_R, &  
                u_dE_norm, & ! <u_dE | u_dE>, i.e. norm 
                node_u, node_d, wronk) ! dummy arguments, unused
end subroutine

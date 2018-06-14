! note: compile with 8-byte reals!
module rad_help
    use binmat
    implicit none
contains
    subroutine readpot(fname, pot)
        character(len=*), intent(in) :: fname
        real, intent(out), allocatable :: pot(:)
        call read_mat(pot, fname)        
    end subroutine
end module
program main
    use m_radfun
    use binmat
    use rad_help
    implicit none
    integer :: node_d, node_u
    integer :: l
    integer :: n_radpts 
    real :: E, r0, dx, r_MT
    real, allocatable :: sol(:, :), sol_dE(:, :), v_radmesh(:)
    real :: u_at_R, &
            u_dr_at_R, &
            u_dE_at_R, &
            u_dr_dE_at_R, &
            u_dE_norm
    ! potential for l=0, dumped from Si_bulk
    character(len=*), parameter :: pot_fname = "vr_jri-521_l-0.bin"
    real :: wronk ! wronskian, unused

    l = 0  ! l-value for Ylm
    E = -0.145868090423439 ! -1. ! sensible value? (energy)
    ! from Si bulk example
    dx = 2.2e-2 ! from inp file, logarithmic increment
    r_MT = 2.16 ! a.u., bohr radii
    ! potential given on the radial mesh
    !allocate(v_radmesh(n_radpts))
    !v_radmesh = 0. ! TODO: needs to contain the l=0 component
    call readpot(pot_fname, v_radmesh)
    !n_radpts = 521 ! number of points in the radial mesh
    n_radpts = size(v_radmesh)
    ! on a radial mesh of the potential
    ! r0 = r_MT / exp(dx) ^ (n_radpts - 1)
    r0 = r_MT * exp(dx * (1 - n_radpts)) ! first mesh point
    ! solution and d(solution)/dE, energy derivative
    ! both have two components, the "big" and "small" part
    allocate(sol(n_radpts, 2), sol_de(n_radpts, 2))
    print*, "Calling radfun with ", &
            "l = ", l, &
            "E = ", E, &
            "n_radpts = ", n_radpts, &
            "r0 = ", r0, &
            "dx = ", dx
    
    call radfun(l, &
                E, &
                v_radmesh, &
                n_radpts, &
                r0, &
                dx, &
                n_radpts, & 
                ! jmtd == max(jri(:)), "relict" from fixed-size days
                sol, sol_dE, &
                u_at_R, &
                u_dr_at_R, &
                u_dE_at_R, &
                u_dr_dE_at_R, &
                u_dE_norm, &
                node_u, node_d, wronk)
    ! dump the solution
    print*, "Dumping solution (derivative) to sol.bin (sol_dE.bin)"
    call write_mat(sol, "sol.bin")
    call write_mat(sol_dE, "sol_dE.bin")
                


end program
subroutine mainwurst
    use m_radfun
    use binmat
    use rad_help
    implicit none
    integer :: node_d, node_u
    integer :: l
    integer :: n_radpts 
    real :: E, r0, dx, r_MT
    real, allocatable :: sol(:, :), sol_dE(:, :), v_radmesh(:)
    real :: u_at_R, &
            u_dr_at_R, &
            u_dE_at_R, &
            u_dr_dE_at_R, &
            u_dE_norm
    ! potential for l=0, dumped from Si_bulk
    character(len=*), parameter :: pot_fname = "vr_jri-521_l-0.bin"
    real :: wronk ! wronskian, unused

    l = 0  ! l-value for Ylm
    E = -0.145868090423439 ! -1. ! sensible value? (energy)
    ! from Si bulk example
    dx = 2.2e-2 ! from inp file, logarithmic increment
    r_MT = 2.16 ! a.u., bohr radii
    ! potential given on the radial mesh
    !allocate(v_radmesh(n_radpts))
    !v_radmesh = 0. ! TODO: needs to contain the l=0 component
    call readpot(pot_fname, v_radmesh)
    !n_radpts = 521 ! number of points in the radial mesh
    n_radpts = size(v_radmesh)
    ! on a radial mesh of the potential
    ! r0 = r_MT / exp(dx) ^ (n_radpts - 1)
    r0 = r_MT * exp(dx * (1 - n_radpts)) ! first mesh point
    ! solution and d(solution)/dE, energy derivative
    ! both have two components, the "big" and "small" part
    allocate(sol(n_radpts, 2), sol_de(n_radpts, 2))
    print*, "Calling radfun with ", &
            "l = ", l, &
            "E = ", E, &
            "n_radpts = ", n_radpts, &
            "r0 = ", r0, &
            "dx = ", dx
    
    call radfun(l, &
                E, &
                v_radmesh, &
                n_radpts, &
                r0, &
                dx, &
                n_radpts, & 
                ! jmtd == max(jri(:)), "relict" from fixed-size days
                sol, sol_dE, &
                u_at_R, &
                u_dr_at_R, &
                u_dE_at_R, &
                u_dr_dE_at_R, &
                u_dE_norm, &
                node_u, node_d, wronk)
    ! dump the solution
    print*, "Dumping solution (derivative) to sol.bin (sol_dE.bin)"
    call write_mat(sol, "sol.bin")
    call write_mat(sol_dE, "sol_dE.bin")
                


end subroutine

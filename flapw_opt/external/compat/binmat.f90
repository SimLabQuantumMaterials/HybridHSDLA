! Write a real or complex matrix into a binary file
!
! File format: small header with 3 entries:
! 1 byte: 'C' or 'R' for real or complex
! 4 bytes: integer, # of rows
! 4 bytes: integer, # of columns
! afterwards, row*col entries in column-major format, per entry:
! if real:
!   8 bytes: double, entry value
! if complex:
!   16 bytes: two doubles, first real, then imaginary part
! To accomplish this with Fortran IO, it needs `access=stream`, which is 
! available since F2013
module binmat
implicit none
private
public :: write_mat, read_mat
interface write_mat
    module procedure write_mat_R, write_mat_C, write_vec_C, write_vec_R
end interface
interface read_mat
    module procedure read_mat_R, read_mat_C, read_vec_R
end interface
integer, parameter :: dp = selected_real_kind(12)
contains !!!!!!!!!!!!!!!!!!!!!!!
subroutine write_vec_R(V, fname)
    real(kind=dp), intent(in) :: V(:)
    character(len=*), intent(in) :: fname
    integer :: length
    integer :: i
    length = size(V)
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'R', length, 1
    do i = 1, length
        write(10) V(i)
    end do
    close(10)
end subroutine
subroutine write_vec_C(V, fname)
    complex(kind=dp), intent(in) :: V(:)
    character(len=*), intent(in) :: fname
    integer :: length
    integer :: i
    length = size(V)
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'C', length, 1
    do i = 1, length
        write(10) V(i)
    end do
    close(10)
end subroutine
subroutine write_mat_R(M, fname)
    real(kind=dp), intent(in) :: M(:, :)
    character(len=*), intent(in) :: fname
    integer :: r, c ! size of M
    integer :: i, j
    r = size(M, 1)
    c = size(M, 2)
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'R', r, c
    do j = 1, c
        do i = 1, r
            write(10) M(i, j)
        end do
    end do
    close(10)
end subroutine
subroutine write_mat_C(M, fname)
    complex(kind=dp), intent(in) :: M(:, :)
    character(len=*), intent(in) :: fname
    integer :: r, c ! size of M
    integer :: i, j
    r = size(M, 1)
    c = size(M, 2)
    open(unit=10, file=fname, form='unformatted', access='stream')
    write(10) 'C', r, c
    do j = 1, c
        do i = 1, r
            write(10) M(i, j)
        end do
    end do
    close(10)
end subroutine
subroutine read_mat_R(M, fname)
    real(kind=dp), intent(out), allocatable :: M(:,:)
    character(len=*), intent(in) :: fname
    character(1) :: dtype
    integer :: r, c, i, j
    open(unit=10, file=fname, form='unformatted', access='stream')
    read(10) dtype, r, c
    !print*, 'Read: ', dtype, r, c
    if (dtype /= 'R') then
        print*, 'not a real matrix, type = ', dtype
        close(10)
        stop
    endif
    allocate(M(r, c))
    do j = 1, c
        do i = 1, r
            read(10) M(i, j)
        end do
    end do
    close(10)
end subroutine
subroutine read_mat_C(M, fname)
    complex(kind=dp), intent(out), allocatable :: M(:,:)
    character(len=*), intent(in) :: fname
    character(1) :: dtype
    integer :: r, c, i, j
    open(unit=10, file=fname, form='unformatted', access='stream')
    read(10) dtype, r, c
    !print*, 'Read: ', dtype, r, c
    if (dtype /= 'C') then
        print*, 'not a complex matrix, type = ', dtype
        close(10)
        stop
    endif
    allocate(M(r, c))
    do j = 1, c
        do i = 1, r
            read(10) M(i, j)
        end do
    end do
    close(10)
end subroutine
subroutine read_vec_R(V, fname)
    real(kind=dp), intent(out), allocatable :: V(:)
    character(len=*), intent(in) :: fname
    character(1) :: dtype
    integer :: length, cols
    integer :: i
    length = size(V)
    open(unit=10, file=fname, form='unformatted', access='stream')
    read(10) dtype, length, cols
    if (dtype /= 'R') then
        print*, 'not a real vector, type = ', dtype
        close(10)
        stop
    endif
    if (cols /= 1) then
        print*, 'not a vector, cols = ', cols
        close(10)
        stop
    endif
    allocate(V(length))
    do i = 1, length
        read(10) V(i)
    end do
    close(10)
end subroutine
end module

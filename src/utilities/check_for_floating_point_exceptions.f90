module check_for_floating_point_exceptions
    use ieee_arithmetic
    implicit none
    public :: check_fp_exceptions
    
    interface check_fp_exceptions
        module procedure check_fp_exceptions_array_dim2
        module procedure check_fp_exceptions_array_dim1
        module procedure check_fp_exceptions_scalar
    end interface check_fp_exceptions

contains

    ! Subroutine to check for exceptions in an array of dim2
    subroutine check_fp_exceptions_array_dim2(input, var_name)
    
        implicit none
        real(8), intent(in) :: input(:,:)     ! Assumed-shape 2D array
        character(len=*), intent(in) :: var_name
        integer :: i, j
        logical :: has_nan, has_inf, has_zero

        has_nan = .false.
        has_inf = .false.
        has_zero = .false.

        do i = 1, size(input, 1)
            do j = 1, size(input, 2)
                if (ieee_is_nan(input(i, j))) then
                    has_nan = .true.
                    write(404,*) 'NaN encountered in variable ', trim(var_name), ' at index (', i, ',', j, ')'
                    stop
                elseif (.not. ieee_is_finite(input(i, j))) then
                    has_inf = .true.
                    write(404,*) 'Infinite value encountered in variable ', trim(var_name), ' at index (', i, ',', j, ')'
                    stop
                elseif (abs(input(i, j)) == 0.0d0) then
                    has_zero = .true.
                    write(404,*) 'Absolute zero encountered in variable ', trim(var_name), ' at index (', i, ',', j, ')'
                end if
            end do
        end do

    end subroutine check_fp_exceptions_array_dim2


    ! Subroutine to check for exceptions in an array of dim1
    subroutine check_fp_exceptions_array_dim1(input, var_name)
        real(8), intent(in) :: input(:)     ! Assumed-shape array
        character(len=*), intent(in) :: var_name
        integer :: i
        logical :: has_nan, has_inf, has_zero

        has_nan = .false.
        has_inf = .false.
        has_zero = .false.

        do i = 1, size(input)
            if (ieee_is_nan(input(i))) then
                has_nan = .true.
                write(404,*) 'NaN encountered in variable ', trim(var_name), ' at index ', i
                stop
            elseif (.not. ieee_is_finite(input(i))) then
                has_inf = .true.
                write(404,*) 'Infinite value encountered in variable ', trim(var_name), ' at index ', i
                stop
            elseif (abs(input(i)) == 0.0d0) then
                has_zero = .true.
                write(404,*) 'Absolute zero encountered in variable ', trim(var_name), ' at index ', i
            end if
        end do

    end subroutine check_fp_exceptions_array_dim1

    ! Subroutine to check for exceptions in a scalar
    subroutine check_fp_exceptions_scalar(input, var_name)
        real(8), intent(in) :: input
        character(len=*), intent(in) :: var_name
        logical :: is_nan, is_inf, is_zero

        is_nan = ieee_is_nan(input)
        is_inf = .not. ieee_is_finite(input)
        is_zero = (abs(input) == 0.0d0)

        if (is_nan) then
            write(404,*) 'NaN encountered in variable ', var_name
            stop
        end if

        if (is_inf) then
            write(404,*) 'Infinite value encountered in variable ', var_name
            stop
        end if

        if (is_zero) then
            write(404,*) 'Absolute zero encountered in variable ', var_name
        end if

    end subroutine check_fp_exceptions_scalar


end module check_for_floating_point_exceptions

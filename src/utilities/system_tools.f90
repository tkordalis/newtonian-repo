module system_tools

    contains

    subroutine execute_command(cmd)
        character(len=100), intent(in) :: cmd
        integer :: status
        integer, external :: system

        ! Execute command
        status = system(cmd)

        if (status /= 0) then
        print *, "Error executing command: ", cmd
        end if

    end subroutine execute_command

    subroutine check_dir(dir)
        character(len=*), intent(in) :: dir
        character(len=100) :: command

        ! Check if directory exists
        command = "if [ ! -d " // trim(dir) // " ]; then mkdir " // trim(dir) // "; fi"

        call execute_command(command)

        ! command ="touch " // trim(dir) // "./results.dat" 

        ! call execute_command(command)

    end subroutine check_dir

end module system_tools
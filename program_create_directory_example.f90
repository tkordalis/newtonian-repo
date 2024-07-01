program create_directory_example
    implicit none
    
    character(len=100) :: command
    character(len=100) :: directory_name
    
    ! Directory name to check/create
    directory_name = "test_directory"
    
    ! Check if directory exists
    command = "if [ ! -d " // trim(directory_name) // " ]; then mkdir " // trim(directory_name) // "; fi"
    
    ! Execute the command
    call execute_command(command)
    
contains

    subroutine execute_command(cmd)
        character(len=100), intent(in) :: cmd
        integer :: status
        
        ! Execute command
        status = system(cmd)
        
        if (status /= 0) then
            print *, "Error executing command: ", cmd
        end if
        
    end subroutine execute_command

end program create_directory_example

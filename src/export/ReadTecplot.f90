

Module ReadTecplot
    Use Formats

    Private 
    Public  :: TecplotFileReader


    Type TecplotFileReader
        Character(len=:), Allocatable :: file
        Character(len=:), Allocatable :: tmp_folder

        contains 
            procedure :: getTitle 
            procedure :: getZonename
            procedure :: getVariable
            procedure :: getElements
            procedure :: getNumberOfNodes
            procedure :: getNumberOfElements
            procedure :: close
    End Type TecplotFileReader

    interface TecplotFileReader
        module procedure constructor
    end interface



    contains

    Function constructor(file) Result(this)
        Implicit None 
        Character(len=*), Intent(In) :: file 
        Type(TecplotFileReader)      :: this

        this%file = file
        call runPythonScript(this%file, this%tmp_folder)
    End Function constructor

    Function getNumberOfNodes(this) Result(output)
        Implicit None 
        Class(TecplotFileReader), Intent(In) :: this
        Integer                              :: output

        Character(len=100)                   :: line    
        Integer                              :: funit

        open(newunit = funit, file = pathJoin(this%tmp_folder,"numberOfNodes"), action='read')
        Read (funit,'(a)') line
        close(funit)

        output = toInt(line)
    End Function getNumberOfNodes


    Function getNumberOfElements(this) Result(output)
        Implicit None 
        Class(TecplotFileReader), Intent(In) :: this
        Integer                              :: output

        Character(len=100)                   :: line    
        Integer                              :: funit

        open(newunit = funit, file = pathJoin(this%tmp_folder,"numberOfElements"), action='read')
        Read (funit,'(a)') line
        close(funit)

        output = toInt(line)
    End Function getNumberOfElements

    Function getTitle(this) Result(output)
        Implicit None 
        Class(TecplotFileReader), Intent(In) :: this
        Character(len=:), Allocatable        :: output

        output = readTitle(this%tmp_folder)
    End Function getTitle

    Function getZonename(this) Result(output)
        Implicit None 
        Class(TecplotFileReader), Intent(In) :: this
        Character(len=:), Allocatable        :: output

        output = readZoneName(this%tmp_folder)
    End Function getZonename


    Function getVariable(this,var) Result(output)
        Implicit None 
        Class(TecplotFileReader), Intent(In) :: this 
        Character(len=*)        , Intent(In) :: var
        Real(8), Dimension(:), Allocatable   :: output

        output = readVariable(this%tmp_folder,var)
    End Function getVariable


    Function getElements(this) Result(output)
        Implicit None 
        Class(TecplotFileReader), Intent(In) :: this
        Integer, Dimension(:,:), Allocatable :: output 

        output = readElements(this%tmp_folder)
    End Function getElements


    Subroutine close(this)
        Implicit None
        Class(TecplotFileReader), Intent(In) :: this

        call removeTemporaryDirectory(this%tmp_folder)  
    End Subroutine close



!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Helpful Routines
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Subroutine runPythonScript(file, dest)
        Implicit None 
        Character(len=*),              Intent(In)   :: file
        Character(len=:), Allocatable, Intent(InOut):: dest

        Character(len=:), Allocatable :: filename
        Integer                       :: funit

        ! the output is a hidden folder
        dest = replace(".{}", "{}", getFileName(file))

        call removeTemporaryDirectory(dest)
        call dumpTecplotToFolders(file,dest)
    End Subroutine runPythonScript


    Function readElements(dest) Result(output)
        Implicit None 
        Character(len=*), Intent(In)        :: dest
        Integer, Dimension(:,:), Allocatable:: output

        Integer                             :: funit
        Integer                             :: i,j
        Integer                             :: nelem, nnodes
    
        open(funit, file = pathJoin(dest,"elements"), action="read")
        read(funit,*) nelem, nnodes

        allocate( output( nelem, nnodes ) )

        if (nnodes /= 3 ) then
            print*, "[Error] readElements : the method is not implemented for nnodes /= 3. Please change."
            stop 
        end if

        do i = 1, nelem
            read(funit,*) output(i,1), output(i,2), output(i,3)
        end do
        close(funit)

    End Function readElements


    Function readVariable(dest, var) Result(output)
        Implicit None 
        Character(len=*), Intent(In)       :: dest
        Character(len=*), Intent(In)       :: var
        Real(8), Dimension(:), Allocatable :: output

        Character(len=:), Allocatable      :: filename
        logical                            :: fileExists
        Integer                            :: funit

        Integer                            :: i
        Integer                            :: n
        !<><><><><><><><><><><><><><><><><><><><><><><><>
        ! Check if file exists 
        !<><><><><><><><><><><><><><><><><><><><><><><><>

        filename = pathJoin(dest,var)
        inquire(file = filename, exist = fileExists)    

        if (.not. fileExists) Then 
            print*, replace("[Error] ReadVariable : The variable * does not exist in the tecplot file.","*", var)
            stop
        End if

        !<><><><><><><><><><><><><><><><><><><><><><><><>
        ! Read file 
        !<><><><><><><><><><><><><><><><><><><><><><><><>

        open(newunit = funit, file = filename, action = "read")
        read (funit,*) n
        allocate( output(n) )
        do i = 1, n
        read (funit,*) output(i)
        end do

        close(funit)
    End Function readVariable

    Function readTitle(dest) Result(output)
        Implicit None 
        Character(len=*), Intent(In) :: dest
        Character(len=:), Allocatable:: output

        Integer, Parameter           :: buffer = 10000
        Character(len=buffer)        :: line
        Integer                      :: funit

        open(newunit = funit, file=pathJoin(dest,'title'), action="read")
        read(funit,'(a)') line
        close(funit)
        output = trim(adjustl(line))
    End Function readTitle

    Function readZoneName(dest) Result(output)
        Implicit None 
        Character(len=*), Intent(In) :: dest
        Character(len=:), Allocatable:: output

        Integer, Parameter           :: buffer = 10000
        Character(len=buffer)        :: line
        Integer                      :: funit

        open(newunit = funit, file = pathJoin(dest,'zonename'), action='read')
        read (funit,'(a)') line
        close(funit)
        output = trim(adjustl(line))
    End Function readZoneName




    Subroutine dumpTecplotToFolders(file, dest)
        Implicit None 
        Character(len=*), Intent(In) :: file 
        Character(len=*), Intent(In) :: dest

        Integer                      :: funit

        Open(newunit = funit, file = ".tmp.py")
        write(funit,"(a)") "import os" 
        write(funit,"(a)") "import binarytecplot as bt" 
        write(funit,"(a)") replace('filename = "*"',"*",file)
        write(funit,"(a)") 'tecline=bt.LoadTecplotFile(filename)'   

        write(funit,"(a)") replace('tecline.dumpToFolder("*")',"*",dest)
        call execute_command_line("python3.8 .tmp.py", wait = .true.)
        Close(funit, status='delete')
    End Subroutine dumpTecplotToFolders



    Function getFileName(filename) Result(filename_)
        ! the tmp folder is hidden and located as the running pwd.
        ! the pattern is .{filename}
        Implicit None 
        Character(len=*), Intent(In)  :: filename

        Character(len=:), Allocatable :: filename_
        Integer                       :: ipos

        if (index(filename,"/") == 0) then
            filename_ = filename
            return 
        end if

        filename_ = reverse(filename)
        filename_ = replace(filename_,"/", " ") 
        ipos      = index(filename_, " ")
        filename_ = reverse(filename_(1:ipos))
        filename_ = trim(adjustl(filename_))
    End Function getFileName


    Subroutine removeTemporaryDirectory(filename)
        Implicit None 
        Character(len=*), Intent(In) :: filename

        call execute_command_line("rm -rf "// filename, wait=.true. )
    End Subroutine removeTemporaryDirectory
        



End Module ReadTecplot





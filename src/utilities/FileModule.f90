!*********************************************************************
! Author: Dionisis Pettas
! Date  : 9/12/2020
! Reads the lines of a file
!*********************************************************************

Module FileModule
    Use Formats

    private 
    public :: fFile, ReadFile, OpenFile


    ! Max length of the line
    Integer, Parameter :: MAXLEN     = 1000
    Integer, Parameter :: BATCHLINES = 1000


    Type fFile
        Character(len=:), Allocatable                    :: filename
        Integer                                          :: nrows 
        Character(len=MAXLEN), Dimension(:), Allocatable :: line ! this will be more efficient with linked lists
        contains 
            procedure :: getFilename
            procedure :: getNumberOfRows
            procedure :: getLine
            procedure :: getLineThatContains
            procedure :: existInFile
            procedure :: push_back
            procedure :: findLineIdThatContains
            procedure :: replaceLinewith
            procedure :: dumpToFile
    End Type fFile


    Interface ReadFile
        module procedure constructor
    End Interface ReadFile

    Interface OpenFile
        module procedure openfile_
    end interface


    contains

    Function constructor(filename) Result(this)
        Implicit None 
        Character(len=*), Intent(In) :: filename
        Type(fFile)                  :: this

        Integer                      :: funit
        Character(len=MAXLEN)        :: line_
        Integer                      :: i

        this%filename = filename
        this%nrows    = countNumberOfRows(filename)

        allocate( this%line( this%nrows ) ) 
        
        open( newunit = funit, file = filename, action = "read")

        do i = 1, this%nrows
            read(funit,"(a)") line_     
            this%line(i) = trim(line_)
        end do  
        close(funit)
    End Function constructor

    Function openfile_(filename) Result(this)
        Implicit None 
        Type(fFile)                  :: this
        Character(len=*), Intent(In) :: filename

        this%filename = filename

        this%nrows = 0
        Allocate(this%line(BATCHLINES))
    End Function openfile_


    Subroutine push_back(this,str) 
        Implicit None 
        Class(fFile)                                     :: this 
        Character(len=*), Intent(In)                     :: str

        Character(len=MAXLEN), Dimension(:), Allocatable :: tmp

        if ( mod(this%nrows,BATCHLINES) == 0 .and. this%nrows /= 0) then 
            allocate  (tmp      , source = this%line  )
            deallocate(this%line)
            allocate  (this%line ( this%nrows + BATCHLINES ) )  
            this%line(1:this%nrows) = tmp(1:this%nrows)
            deallocate(tmp)
        end if


        this%nrows            = this%nrows + 1
        this%line(this%nrows) = str

    End Subroutine push_back



    Subroutine replaceLinewith(this, id, str) 
        Implicit None 
        Class(fFile)    , Intent(InOut) :: this
        Integer                         :: id
        Character(len=*), Intent(In)    :: str


        if ( id < 1 ) then 
            print*, "[Error] fFile: the line id should be larger than 1"
            stop 
        end if

        if ( id > this%nrows ) then 
            print*, "[Error] fFile: the line id is larger that the number of rows of the file"
            stop
        end if 
        
        this%line(id) = str
    End Subroutine replaceLinewith


    Subroutine dumpToFile(this, filename) 
        Implicit None 
        Class(fFile)    , Intent(In) :: this 
        Character(len=*), Intent(In) :: filename

        Integer                      :: funit
        Integer                      :: i

        open( newunit = funit, file = filename, action = "write")
        do i = 1, this%nrows
            write(funit,"(a)") trim(this%line(i))
        end do  
        close(funit)
    End Subroutine dumpToFile

    Function findLineIdThatContains(this, str) Result(id) 
        Implicit None 
        Class(fFile)    , Intent(In) :: this
        Character(len=*), Intent(In) :: str
        Integer                      :: id
        
        Character(len=MAXLEN)        :: line_1
        Character(len=MAXLEN)        :: line_2
        Integer                      :: i

        do i = 1, this%nrows
            line_1 = this%line(i)
            line_2 = replace(line_1, str,"")
            
            id = i 
            if ( line_1 /= line_2 ) return 
        end do
        id = -1 
    End Function findLineIdThatContains

    Function existInFile(this, str) Result(output)
        Implicit None 
        Class(fFile)    , Intent(In) :: this
        Character(len=*), Intent(In) :: str

        Logical                      :: output

        output = this%findLineIdThatContains(str) /= -1
    End Function existInFile


    Function getLineThatContains(this, str) Result(line)
        Implicit None 
        Class(fFile)    , Intent(In)  :: this
        Character(len=*), Intent(In)  :: str
        Character(len=:), Allocatable :: line
    
        Integer                       :: id

        id = this%findLineIdThatContains(str)
        if ( id == -1) then 
            print*, replace("[Error] fFile: there is no line in the file that contains *.","*",str)
            stop
        end if 

        line = this%getLine( id )
    End Function getLineThatContains
        

    Function getNumberOfRows(this) Result(output)
        Implicit None 
        Class(fFile) :: this
        Integer      :: output
        
        output = this%nrows

    End Function getNumberOfRows

    Function getFilename(this) Result(output)
        Implicit None 
        Class(fFile)    , Intent(In)  :: this
        Character(len=:), Allocatable :: output

        output = this%filename
    End Function getFilename



    Function getLine(this, id) Result(line)
        Implicit None 
        Class(fFile), Intent(In)      :: this
        Integer     , Intent(In)      :: id     
        Character(len=:), Allocatable :: line
        
        if ( id < 1 ) then 
            print*, "[Error] fFile: the line id should be larger than 1"
            stop 
        end if

        if ( id > this%nrows ) then 
            print*, "[Error] fFile: the line id is larger that the number of rows of the file"
            stop
        end if 

        line = trim(this%line(id))
    End Function getLine



    Function countNumberOfRows(filename) Result(rows)
        Implicit None 
        Character(len=*), Intent(In) :: filename
        Integer                      :: rows

        Integer                      :: funit
        Character(len=1)             :: tmp
        Integer                      :: io

        rows = 0
        open( newunit = funit, file = filename, action="Read")
        do 
            read(funit,"(a)", iostat = io)  tmp
            if ( io /= 0 ) exit
            rows = rows + 1
        end do
        close( funit )
    End Function countNumberOfRows

End Module FileModule


Module UNV
    Use Formats
    Use MeshGeneration, only: python

 Private 
 

 Public  :: unvFileReader


 Type unvFileReader
    Character(len=:)       , Allocatable :: filename
    Character(len=:)       , Allocatable :: foldername
    Real(8), Dimension(:,:), Allocatable :: nodes 
    Integer, Dimension(:,:), Allocatable :: elements
    Integer, Dimension(:,:), Allocatable :: surface_elements 
    Integer, Dimension(:)  , Allocatable :: faces

    contains 
        procedure :: readNodes 
        procedure :: readElements
        procedure :: readSurfaceElements
        procedure :: getNodes
        procedure :: getElements
        procedure :: getSurfaceElements
        procedure :: getBoundaryName
        procedure :: getBoundaryId
        procedure :: getNumberOfBoundaries
        procedure :: getNumberOfNodes
        procedure :: getNumberOfElements
        procedure :: getNumberOfSurfaceElements
        procedure :: getBoundary
        procedure :: info
        procedure :: close
        procedure :: FileNotFound
        procedure :: executePythonScript
        procedure :: dumpToTecplotFile
 End Type unvFileReader


 Interface unvFileReader 
    module procedure constructor
 End Interface unvFileReader


 contains


    Function constructor(filename) Result(this)
        Implicit None 
        Character(len=*), Intent(In) :: filename
        Type(unvFileReader)          :: this

        this%filename         = filename
        this%foldername       = concatenate(".",getFilenameFromPath(this%filename))

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Check if File exist
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        if ( .not. this%FileNotFound()) Then 
            print*, replace("[Error] unvFileReader : The file * does not exist.","*",this%filename)
            stop
        end if
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        call this%executePythonScript()
        call this%readNodes          ("nodes.dat")
    call this%readElements       ("elements.dat")
        Call this%readSurfaceElements("surface_elements.dat")   
    End Function constructor


    ! This method is called like :
    ! a = (object)%getNumberOfNodes()  
    ! (this) is ommited ALWAYS
    Function getNumberOfNodes(this) Result(output)
        Implicit None 
        Class(unvFileReader), Intent(In) :: this
        Integer                          :: output

        output = size(this%nodes,1)
    End Function getNumberOfNodes

    Function getNumberOfElements(this) Result(output)
        Implicit None 
        Class(unvFileReader), Intent(In) :: this
        Integer                          :: output

        output = size(this%elements,1)
    End Function getNumberOfElements


    Function getNumberOfSurfaceElements(this) Result(output)
        Implicit None 
        Class(unvFileReader), Intent(In) :: this
        Integer                          :: output

        output = size(this%surface_elements,1)
    End Function getNumberOfSurfaceElements

    Function getNumberOfBoundaries(this) Result (output)
        Implicit None 
        Class(unvFileReader), Intent(In) :: this
        Integer                          :: output
    
        Character(len=:), Allocatable    :: item
        Integer                          :: tmpfile
    
        item       = pathJoin( this%foldername, "bnd_*_elements.dat")

        Call Execute_Command_Line("ls -l "//item//"| wc -l > .boundaries.list")

        Open(newunit = tmpfile, file = ".boundaries.list")
        Read (tmpfile,*) output 
        Close(tmpfile, status = 'delete')   
        Deallocate(item)
    End Function getNumberOfBoundaries


    Function getBoundaryName(this, nBnd) Result(output)
        Implicit None 
        Class(unvFileReader), Intent(In) :: this
        Integer             , Intent(In) :: nBnd        
        Character(len=:), Allocatable    :: output

        Character(len=300)           :: line
        Character(len=:), Allocatable:: item
        Integer                      :: tmpfile
        Integer                      :: i

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Constrains
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        if (nBnd < 0 )&
        stop "[Error] getBoundaryName : The nBnd should be positive"

        if (nBnd > this%getNumberOfBoundaries())&
        stop "[Error] getBoundaryName : The number of boundaries exceeds the maximum value."
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        item = pathJoin(this%foldername,"bnd_*_elements.dat")
        call execute_command_line(replace("ls -v * > .boundaries.list","*", item), wait = .true.)   

        Open(newunit = tmpfile, file = ".boundaries.list")  
        ! skip lines        
        do i = 1, nBnd - 1 ; Read(tmpfile,*) ; end do       
        Read (tmpfile,'(a)') line   
        Close(tmpfile, status = 'delete' )

        line   = remove( trim(adjustl(line)), [this%foldername,"/"])
        output = remove( trim(adjustl(line)), ["bnd_", toStr(nBnd),"_", "_elements.dat"])   

    End Function getBoundaryName


    Function getBoundaryId(this, name) Result(output) 
        Implicit None 
        Class(unvFileReader), Intent(In) :: this
        Character(len=*)    , Intent(In) :: name
        Integer                          :: output

        Integer                          :: i 

        do i = 1, this%getNumberOfBoundaries()
            if (this%getBoundaryName(i) == name ) exit      
        end do
        output = i
    
    End Function getBoundaryId

    Subroutine getBoundary (this, name, elements, faces)
        Implicit None
        Class(unvFileReader)              , Intent(In) :: this
        Character(len=*)                  , Intent(In) :: name
        Integer, Dimension(:), Allocatable             :: elements 
        Integer, Dimension(:), Allocatable             :: faces

        Character(len=:)     , Allocatable :: bnd

        ! Read Boundary Elements
        bnd = pathJoin( this%foldername, replace("bnd_*_*_elements.dat","*", [ toStr(this%getBoundaryId(name)), name] ))
        call readBoundaryFile(bnd, elements)

        ! Read Boundary Faces
        bnd = pathJoin( this%foldername, replace("bnd_*_*_faces.dat","*", [toStr(this%getBoundaryId(name)), name] ) )
        call readBoundaryFile(bnd, faces)   
    End Subroutine getBoundary

  Subroutine executePythonScript (this)
        Implicit None 
        Class(unvFileReader), Intent(In) :: this
        Integer                          :: pyscript

        Open( newunit = pyscript, file = ".splitunv.py")

        Write(pyscript,'(a)') 'import unv'
        Write(pyscript,'(a)') replace('filename="*"'                           ,'*', this%filename   )
        Write(pyscript,'(a)') replace('f=unv.Reader("*")'                      ,'*', this%filename   )
        Write(pyscript,'(a)') replace('f.exportTofolder("*", zero_based=False)','*', this%foldername )
        ! Call execute_command_line( "python3.8 .splitunv.py > /dev/null", wait = .true.)
        Call execute_command_line( concatenate(python," ", ".splitunv.py"," ","> /dev/null"), wait = .true.)
        Close(pyscript, status = "delete")
    End Subroutine executePythonScript


    Subroutine readNodes(this, filename)
        Implicit None 
        Class(unvFileReader)             :: this
        Character(len=*)    , Intent(In) :: filename

        Integer                                             :: nodefile
        Integer                                             :: num
        Integer         , Parameter                         :: buffer = 300
        Character(len=buffer)                               :: line     
        Integer                                             :: i


        If ( Allocated(this%nodes) ) Deallocate(this%nodes)
 
        Open(newunit = nodefile, file = pathJoin(this%foldername, filename) ) 

        ! Read The Maximum Values of the file (num)
        Read(nodefile ,'(a)') line
        
        line = remove(line,"#")
        num  = toInt (line)

        Allocate( this%nodes(num,3) ) ! x, y, z

        Do i = 1, num ; Read(nodefile,*) this%nodes(i,1), & 
                                             this%nodes(i,2), &
                                                                         this%nodes(i,3)
        End Do 

        Close(nodefile)
    End Subroutine readNodes


    Subroutine readElements(this, filename)
        Implicit None 
        Class(unvFileReader)                                :: this
        Character(len=*), Intent(In)                        :: filename


        Integer                                             :: elementfile
        Integer                                             :: num
        Integer         , Parameter                         :: buffer = 300
        Character(len=buffer)                               :: line     
        Integer                                             :: i


        If ( Allocated( this%elements) ) Deallocate( this%elements)
 

        Open(newunit = elementfile, file = pathJoin(this%foldername,filename) ) 
        ! Read The Maximum Values of the file (num)
        Read(elementfile,'(a)') line
        num = toInt( remove(line,"#") )

        Allocate( this%elements(num,3) ) 

        Do i = 1, num
        Read(elementfile,*) this%elements(i,1), &
                        this%elements(i,2), &
                        this%elements(i,3)
        End Do 

        Close(elementfile)
    End Subroutine readElements

    Subroutine readSurfaceElements(this, filename)
        Implicit None 
        Class(unvFileReader), Intent(InOut)                 :: this 
        Character(len=*)    , Intent(In)                    :: filename

        Integer         , Parameter                         :: buffer = 300
        Character(len=buffer)                               :: line     

        Integer                                             :: i
        Integer                                             :: num
        Integer                                             :: elementfile


        If ( Allocated(this%surface_elements) ) Deallocate(this%surface_elements)
 
        Open(newunit = elementfile, file = pathJoin(this%foldername, filename)) 

        ! Read The Maximum Values of the file (num)
        Read(elementfile,'(a)') line
        num = toInt( remove(line,"#") )

        Allocate( this%surface_elements(num,2) ) 

        Do i = 1, num
        Read(elementfile,*) this%surface_elements(i,1), &
                        this%surface_elements(i,2)
        End Do 

        Close(elementfile)
    End Subroutine readSurfaceElements


    Function getElements(this) Result(output)
        Implicit None 
        Class(unvFileReader), Intent(In)     :: this 
        Integer, Dimension(:,:), Allocatable :: output

        allocate(output, source = this%elements)    
    End Function getElements

    Function getNodes(this) Result(output)
        Implicit None 
        Class(unvFileReader), Intent(In)     :: this 
        Real(8), Dimension(:,:), Allocatable :: output 

        allocate(output, source = this%nodes)
    End Function getNodes

    Function getSurfaceElements(this) Result(output)
        Implicit None 
        Class(unvFileReader), Intent(In)     :: this 
        Integer, Dimension(:,:), Allocatable :: output

        allocate(output, source = this%surface_elements)    
    End Function getSurfaceElements


    Subroutine dumpToTecplotFile(this, filename, coordinates_name)
        Implicit None 
        Class(unvFileReader), Intent(In)           :: this 
        Character(len=*)    , Intent(In), Optional :: filename
        Character(len=*)    , Intent(In), Optional :: coordinates_name

        Character(len=:)    , Allocatable          :: pythonfilename
        Integer                                    :: python_script
        Character(len=:)    , Allocatable          :: output
        Character(len=:)    , Allocatable          :: coordinates_name_

        if  ( present(filename) ) Then ; output = filename
        Else                           ; output = replace( this%filename, ".unv", ".plt") 
        End if

        if  ( present(coordinates_name) ) Then ; coordinates_name_ = coordinates_name
        Else                                   ; coordinates_name_ = "X Y"
        End if


        pythonfilename = concatenate(".python38_", replace(getFilenameFromPath(this%filename),".unv",".py"))

        open( newunit = python_script, file = pythonfilename )
        write(python_script,"(a)") "import unv"
        write(python_script,"(a)") replace("f = unv.Reader('*')","*", this%filename)
        write(python_script,"(a)") replace("f.toAsciiTecplot('*','*',renameCoords = '*')","*", [output,"from unv", coordinates_name_])

        ! call execute_command_line( replace("python3.8 *","*", pythonfilename) )
        call execute_command_line( concatenate(python," ", pythonfilename) )
        close(python_script,status = "delete")


        If (allocated(output)) deallocate(output)
    End Subroutine dumpToTecplotFile





    Subroutine close(this)
        Implicit None
        Class(unvFileReader), Intent(InOut) :: this
        Logical                             :: dir_exist

        If (Allocated(this%nodes)           ) Deallocate(this%nodes           )
    If (Allocated(this%elements)        ) Deallocate(this%elements        )
        If (Allocated(this%surface_elements)) Deallocate(this%surface_elements)
        If (Allocated(this%faces)           ) Deallocate(this%faces           )

        inquire ( directory = this%foldername , exist = dir_exist )

        If (dir_exist)        call Execute_Command_Line("rm -rf "//this%foldername)
        deallocate(this%filename)
        deallocate(this%foldername)
    End Subroutine close


    Subroutine info(this)
        ! use ConstitutiveEquationOfTheBubble
        Implicit None 
        Class(unvFileReader), Intent(In) :: this

        Character(len=:), Allocatable :: buffer
        Character(len=:), Allocatable :: files
        Character(len=:), Allocatable :: line 
        Integer                       :: ibnd
        character(25)                 :: dateNtime
        

        files = ""
        do ibnd = 1, this%getNumberOfBoundaries()
            files  = concatenate(files, this%getBoundaryName(ibnd))
            files  = concatenate(files,", ")
        end do 
        files = files(1:len(files)-2)
        call fdate(dateNtime)
        print*, dateNtime
        ! print*, "Constitutive Equation for the Bubbles: ", VorPV
        print*, "Unv Info ..."
        print*, concatenate("  Filename : "                          , this%filename  )
        print*, concatenate("  Exported Folder           : "         , this%foldername)
        print*, concatenate("  Number of Nodes           : "         , toStr(this%getNumberOfNodes()          ))
        print*, concatenate("  Number of Elements        : "         , toStr(this%getNumberOfElements()       ))
        print*, concatenate("  Number of surface elements: "         , toStr(this%getNumberOfSurfaceElements()))
        print*, concatenate("  Number of Boundaries      : "         , toStr(this%getNumberOfBoundaries()     ))
        print*, concatenate("  Boundary Names            : "         , files                                   )
        print*, replace    ("  Coordinate X : Min = *, Max = *  ","*", [toStr(minval(this%nodes(:,1))), &
                                                                        toStr(maxval(this%nodes(:,1)))]) 
        print*, replace    ("  Coordinate Y : Min = *, Max = *  ","*", [toStr(minval(this%nodes(:,2))), &
                                                                          toStr(maxval(this%nodes(:,2)))]) 

    End Subroutine info

    Function FileNotFound(this) Result(output)
        Implicit None 
        Class(unvFileReader)         :: this
    Logical                      :: output
    Inquire(file = this%filename, exist = output)
    End Function FileNotFound


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Helpful Functions
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    Subroutine readBoundaryFile(bndfile, array)
        Implicit None 
        Character(len=*), Intent(In)       :: bndfile 
        Integer, Dimension(:), Allocatable :: array

        Integer                            :: bndunit
        Character(len=100)                 :: header
        Integer                            :: n
        Integer                            :: i

        Open(newunit = bndunit, file = bndfile) 
        Read(bndunit,"(a)") header  
        header = remove(header,"#")     
        n      = toInt(trim(adjustl(header)))

        If ( ALlocated(array) ) Deallocate(array)
        Allocate(array (n) )

        do i = 1, n
            Read(bndunit,*) array(i)
        end do
        Close(bndunit)
    End Subroutine readBoundaryFile

End Module UNV




Module TecplotUtilities

    Type TecVar 
        Character(len=:)     , Allocatable :: name_ 
        Real(8), Dimension(:), Allocatable :: value_
    End Type TecVar

    interface TecVar
        module procedure TecVar_Default_Constructor 
        module procedure TecVar_Copy_Constructor
    end interface TecVar


    contains 

        Function TecVar_Default_Constructor (name_, value_) Result(output)
            Implicit None 
            Character(len=*)     , Intent(In) :: name_
            Real(8), Dimension(:), Intent(In) :: value_

            Type(TecVar)                      :: output

            output%name_  = name_ 
            if ( allocated(output%value_) ) Allocate(output%value_, mold = value_)
            output%value_ = value_
        End Function TecVar_Default_Constructor

        Function TecVar_Copy_Constructor (other) Result(this)
            Implicit None 
            Type(TecVar), Intent(In) :: other   
            Type(TecVar)             :: this
        
            this%name_ = other%name_ 
            this%value_= other%value_ 
        End Function TecVar_Copy_Constructor

End Module TecplotUtilities





module WriteTecplot
    Use Formats
    Use TecplotUtilities
    Private 
    Public :: Tecplot_File


    Type Tecplot_File   
        Character(len=:)            , Allocatable :: filename
        Character(len=:)            , Allocatable :: title
        Character(len=:)            , Allocatable :: zonename

        Integer                                   :: nnodes   = 0
        Integer                                   :: nelements= 0
        Character(len=:)            , Allocatable :: datapacking
        Character(len=:)            , Allocatable :: femtype
        Integer                                   :: nVars
        Type(TecVar), Dimension(:)  , Allocatable :: Var
        Integer     , Dimension(:,:), Allocatable :: elements
        contains    
            procedure :: addVariable 
            procedure :: addElements
            procedure :: VariableLine
            procedure :: setTitle
            procedure :: setZonename
            procedure :: header
            procedure :: writePointStructure
            procedure :: writeBlockStructure
            procedure :: dumpToFile 
            procedure :: dumpToFile_Binary
            procedure :: convertToBinary
            procedure :: getNumberOfVariables
            procedure :: clear
            final     :: deconstructor
    End Type Tecplot_File   

    Interface Tecplot_File
        module procedure Tecplot_File_Constructor
    End Interface Tecplot_File


    contains

        Function Tecplot_File_Constructor(datapacking, femtype) Result(this)
            Implicit None 
            Type(Tecplot_File)            :: this 
            Character(len=*), Intent(In)  :: datapacking
            Character(len=*), Intent(In)  :: femtype

            if (datapacking /= "POINT" .and. datapacking /="BLOCK" ) Then 
                Print*, "[Error] Tecplot_File."
                Print*, "The datapacking should be either POINT or BLOCK"
                stop
            end if


            Allocate( this%Var (30) )

            this%nVars     = 0
            this%filename  = ""     
            this%title     = ""
            this%nnodes    = 0
            this%nelements = 0
            this%datapacking = datapacking
            this%femtype     = femtype
        End Function Tecplot_File_Constructor

        Function getNumberOfVariables(this) Result(output) 
            Implicit None 
            Class(Tecplot_File), Intent(In) :: this
            Integer                         :: output

            output = this%nVars
        End Function getNumberOfVariables

        Subroutine addVariable(this, name_, value_, cellcentered)
            Implicit None 
            Class(Tecplot_File)                     :: this 
            Character(len=*)     , Intent(In)       :: name_ 
            Real(8), Dimension(:), Intent(In)       :: value_
            Logical, intent(in), optional           :: cellcentered
            Integer                                 :: i


            this%nVars        = this%nVars + 1
            if (present(cellcentered)) then
                this%nelements = size(value_)
            else
                this%nnodes       = size(value_)
            endif
            i                 = this%nVars
            if (this%nVars > 30 ) then 
                print*, "Tecplot_File: Maximum Value Reached."
                stop
            end if 

            this%Var(i)%name_ = name_ 
            this%Var(i)%value_= value_
            
        End Subroutine addVariable

        Subroutine clear(this)
            Implicit None 
            Class(Tecplot_File) :: this

            
            if (allocated(this%Var)      ) Deallocate(this%Var) 
            if (allocated(this%elements) ) Deallocate(this%elements)
        End Subroutine clear

        Subroutine deconstructor(this)
            Implicit None 
            Type(Tecplot_File) :: this 
            
            this%filename    = "" 
            this%nnodes      = 0
            this%nelements   = 0
            this%datapacking = ""
            this%femtype     = ""
            If ( allocated( this%Var)      ) Deallocate(this%Var)
            If ( allocated( this%elements) ) Deallocate(this%elements)
        End Subroutine deconstructor

            
        Subroutine addElements(this, elements) 
            Implicit None 
            Class(Tecplot_File)                 :: this 
            Integer, Dimension(:,:), Intent(In) :: elements

            this%nelements = size(elements,1)
            allocate( this%elements, mold = elements )
            this%elements  = elements   

        End Subroutine addElements


        Subroutine setTitle(this, title)
            Implicit None 
            Class(Tecplot_File), Intent(InOut) :: this
            Character(len=*)   , Intent(In)    :: title

            this%title = title
        End Subroutine setTitle


        Subroutine setZonename(this, zonename)
            Implicit None 
            Class(Tecplot_File), Intent(InOut) :: this
            Character(len=*)   , Intent(In)    :: zonename

            this%zonename = zonename
        End Subroutine setZonename

        Function VariableLine(this) Result(line)
            Implicit None 
            Class(Tecplot_File), Intent(In)  :: this 
            Character(len=:)   , Allocatable :: line 
            
            Integer                          :: i

            line = ""
            do i = 1, this%nVars
                line = line//'"'//this%Var(i)%name_//'",'
            end do 
            ! remove last comma
            line = line(1:len(line)-1)
            line = "VARIABLES = "//line
        End Function VariableLine


        Function header(this,zonename, datapacking, zonetype, cellcentered) Result(line)
            Implicit None
            Class(Tecplot_File)           :: this
            Character(len=*), Intent(In)  :: zonename
            Character(len=*), Intent(In)  :: datapacking
            Character(len=*), Intent(In)  :: zonetype
            Logical,Intent(In), optional  :: cellcentered
            Character(len=:), Allocatable :: line

            line = ""
            line = line//replace('ZONE T="*"'   ,"*", zonename)                         //","
            line = line//replace("DATAPACKING=*","*", datapacking)                      //","

            if ( present(cellcentered) .and. (cellcentered .eqv. .true.) ) then
                line = line//replace("VARLOCATION=([*,*]=NODAL, [*]=CELLCENTERED)", &
                                                   "*", [toStr(1),toStr(2),toStr(3)])   //","
            endif

            line = line//replace("N=*"          ,"*", toStr(this%nnodes) )              //","
            line = line//replace("E=*"          ,"*", toStr(this%nelements) )           //","
            line = line//replace("ZONETYPE=*"   ,"*", zonetype)
        End Function header


        Subroutine writeBlockStructure(this, unit) 
            Implicit None 
            Class(Tecplot_File), Intent(In) :: this 
            Integer            , Intent(In) :: unit

            Type(TecVar)                    :: variable
            Integer                         :: i,j 

            ! modyfying the block structure subroutine to write the two first variables 
            ! in a nodal form [ Z, R ] and the third in a cell-centered form [Stretch]
            ! writeBlockStructure should be used ONLY for the mesh quality files
            do i = 1, this%nVars - 1
                variable = this%var(i)
                write(unit,*) (variable%value_(j), j = 1, this%nnodes)
            end do

            variable = this%var(this%nVars)

            write(unit,*) (variable%value_(j), j = 1, this%nelements)
            
        End Subroutine writeBlockStructure

        Subroutine writePointStructure(this, unit) 
            Implicit None 
            Class(Tecplot_File), Intent(In) :: this 
            Integer            , Intent(In) :: unit

            Integer                         :: nvariables
            Real(8)                         :: val
            Character(len=:), Allocatable   :: line
            Integer                         :: node 
            Integer                         :: variable

            nvariables = this%nVars
            do node = 1, this%nnodes
                line = ""
                do variable = 1, nvariables
                    val  = this%Var(variable)%value_(node)
                    line = line//replace(" * ","*",toStr(val, fmt = "F20.12" )  )
                end do  
                write(unit,'(a)') line
            end do
    
        End Subroutine writePointStructure


        Subroutine dumpToFile(this, filename)
            Implicit None 
            Class(Tecplot_File)             :: this 
            Character(len=*)   , Intent(In) :: filename

            Integer                         :: tfile
            Type(TecVar)                    :: variable
            Integer                         :: v
            Integer                         :: i
            Integer                         :: j

            this%filename    = filename

            if ( containFolders(this%filename) ) call createFolder(this%filename)

            open(newunit = tfile, file = this%filename)
            write(tfile,"(a)") replace('TITLE = "*"',"*", this%title )

            write(tfile,"(a)") this%VariableLine()
    
            if      (this%datapacking == "BLOCK" )Then; write(tfile,"(a)") this%header(this%zonename, "BLOCK", "FETRIANGLE", cellcentered = .True.)
            else if (this%datapacking == "POINT" )Then; write(tfile,"(a)") this%header(this%zonename, "POINT", "FETRIANGLE")
            end  if 


            if      (this%datapacking == "BLOCK" )Then; call this%writeBlockStructure(tfile)
            else if (this%datapacking == "POINT" )Then; call this%writePointStructure(tfile)
            end  if 


            do j = 1, this%nelements
                write(tfile,"(*(i6,1x))") this%elements(j,:)
            end do 
            

            close(tfile)

        End Subroutine dumpToFile

        Subroutine convertToBinary(this)
            Implicit None 
            Class(Tecplot_File), Intent(In) :: this
            Character(len=:), Allocatable   :: cmd

            cmd = "preplot $1 $2 > /dev/null && mv $3 $4 "
            cmd = replace(cmd,"$1", this%filename     )
            cmd = replace(cmd,"$2", this%filename//"_")
            cmd = replace(cmd,"$3", this%filename//"_")
            cmd = replace(cmd,"$4", this%filename     )

            call execute_command_line(cmd, wait =.true.)
        End Subroutine convertToBinary

        Subroutine dumpToFile_Binary(this, filename, title)
            Implicit None 
            include 'tecio.f90'
            Class(Tecplot_File)             :: this 
            Character(len=*)   , Intent(In) :: title
            Character(len=*)   , Intent(In) :: filename

            Integer                         :: nnodes
            Integer                         :: nelements
            Integer                         :: nVars

            Integer(4)                      :: NPts      ! NPts   - NUMBER OF NODES
            Integer(4)                      :: NElm      ! NElm   - NUMBER OF ELEMENTS
            Integer(4)                      :: KMax      ! KMax   - NUMBER OF FACES // NOT USED
            Integer(4)                      :: StrandID  ! StrandID   - ASSOCIATED ZONES  // NOT USED
            Integer(4)                      :: ParentZn  ! Parent element / Not Used
            Integer(4)                      :: DIsDouble
            Real(8)                         :: SolTime
            Integer(4)                      :: Debug
            Integer(4)                      :: VIsDouble
            Integer(4)                      :: FileType

            Pointer(NullPtr, Null) 
            Integer(4)                      :: Null(*)
            Character(len=1)                :: NullChar
            Integer                         :: i, ivar, iii
            Character(len=:), Allocatable   :: varline 
            Real(4), Dimension(:), Allocatable :: tmp_var
            Integer :: element(3)
            Integer(4), Dimension(:,:), Allocatable :: elements

            nnodes    = this%nnodes
            nelements = this%nelements
            nVars     = this%nVars

            !*******************************************************
            ! Binary Structure
            !*******************************************************
            Debug     = 0 
            VIsDouble = 1
            FileType  = 0
            NullPtr   = 0
            
            NullChar  = Char(0)

            !*******************************************************
            !*******************************************************
            varline = ""
            do ivar = 1, this%getNumberOfVariables()
                varline = varline//this%Var(ivar)%name_//" "
            end do
            varline = trim(adjustl(varline) )
            !*******************************************************
            ! Open File
            !*******************************************************
            i = tecini112(   title //NullChar,&
                           varline //NullChar,&
                                         filename//NullChar,&
                                         "."     //NullChar,&
                                         FileType          ,&
                                         Debug             ,&
                                         VIsDouble)

            i = tecFil112(1)
            !*******************************************************
            ! Header
            !*******************************************************
            NPts    = nnodes 
            Nelm    = size(this%elements,1)
            KMax    = 0
            SolTime = 0.0
            StrandID= 0
            ParentZn= 0

            i = TecZne112(title//NullChar,&
                           2,             &  ! FETRIANGLE
                           NPts,          &
                           NElm,          &
                           KMax,          &
                           0,             &
                           0,             &
                           0,             &
                           SolTime,       &
                           StrandID,      &
                           ParentZn,      &
                           1,             &
                           0,             &
                           0,             &
                           0,             &
                           0,             &
                           0,             &
                           Null,          &
                           Null,          &
                           Null,          &
                           0)

            !*******************************************************
            ! Tecplot Vars
            !*******************************************************
            III       = NPts    
            DIsDouble = 1   
    
            allocate(tmp_var (NPts) )   
            do ivar = 1, this%getNumberOfVariables()
                tmp_var = real(this%Var(ivar)%value_,4)
                i = tecDat112(III, tmp_var, DIsDouble)
            end do
            Deallocate(tmp_var)



            !*******************************************************
            ! Elements
            !*******************************************************
            allocate( elements, mold = this%elements) 
            do i = 1, this%nelements
            print*, i

            elements(i,:) = this%elements(i,:)
            end do
            i = tecNod112( elements )
            !*******************************************************
            ! Close File
            !*******************************************************
            i = tecEnd112()
        End Subroutine dumpToFile_Binary

end module WriteTecplot





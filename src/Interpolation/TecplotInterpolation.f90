!*****************************************************************
! Author : Dionisis Pettas
! Date   : 10/12/2020
!*****************************************************************

                                                     
Module TecplotInterpolation
    Use Formats
    Use WriteTecplot
    Use ReadTecplot
    Use MeshGeneration


    private 
    public :: TecplotInterpolator

    Type TecplotInterpolator
        Real(8)           , Dimension(:)  , Allocatable :: x ! reference coordinates
        Real(8)           , Dimension(:)  , Allocatable :: y
        Integer           , Dimension(:,:), Allocatable :: elements
        Character(len=:)  ,                 Allocatable :: varnamelist
        Character(len=:)  ,                 Allocatable :: output_name
        Type(Tecplot_File)                              :: tec

        contains
            procedure :: referenceDomain
            procedure :: mapToReferenceDomain
            procedure :: interpolateInto
            procedure :: executeScript
            procedure :: setTecplotName
            procedure :: getInterpolatedValues
            procedure :: getElements
            procedure :: close
    End Type TecplotInterpolator

    Interface TecplotInterpolator
        module procedure constructor
    End Interface TecplotInterpolator

    Character(len=*), Parameter :: tecplot = "tecfocus"



    contains 

    Function constructor(x, y, elements) Result(this)
        Implicit None 
        Type(TecplotInterpolator) :: this
        Real(8), Dimension(:)    , Intent(In) :: x
        Real(8), Dimension(:)    , Intent(In) :: y
        Integer, Dimension(:,:)  , Intent(In) :: elements

        call this%referenceDomain(x, y, elements)
    End Function constructor

    Subroutine referenceDomain(this, x, y, elements)
        Implicit None
        Class(TecplotInterpolator), Intent(InOut) :: this
        Real(8), Dimension(:)     , Intent(In)    :: x
        Real(8), Dimension(:)     , Intent(In)    :: y
        Integer, Dimension(:,:)   , Intent(In)    :: elements
        
        if ( .not. allocated(this%x)        ) allocate( this%x       , source = x        )
        if ( .not. allocated(this%y)        ) allocate( this%y       , source = y        )
        if ( .not. allocated(this%elements) ) allocate( this%elements, source = elements )
    End Subroutine referenceDomain


    Subroutine mapToReferenceDomain(this, time, solution, variableName) 
        Implicit None 
        Class(TecplotInterpolator), Intent(InOut) :: this
        Real(8)                   , Intent(In)    :: time
        Real(8), Dimension(:,:)   , Intent(In)    :: solution

        Integer                                   :: id
        Integer                                   :: nunknowns
        Interface 
            Function variableName( vid ) Result(var)
                Implicit None 
                Integer         , Intent(In)  :: vid
                Character(len=:), Allocatable :: var
            End Function variableName
        End Interface
        


        nunknowns = size(solution,2)

        this%tec  = Tecplot_File( datapacking = "POINT", femtype = "FETRIANGLE")


        call this%tec%setZonename( toStr(time)   )

        ! Write Reference Domain
        call this%tec%addVariable ("x_tec"  , this%x  )
        call this%tec%addVariable ("y_tec"  , this%y  )

        this%varnamelist = ''
        this%varnamelist = concatenate(this%varnamelist, '"x_tec"',' ')
        this%varnamelist = concatenate(this%varnamelist, '"y_tec"',' ')
            
        do id = 1, nunknowns
        call this%tec%addVariable( variableName(id) , solution(:,id))
        this%varnamelist = concatenate(this%varnamelist, replace('"*"',"*", variableName(id)) ,' ')
        end do
        
        call this%tec%addElements( this%elements )

    !call this%tec%dumpToFile("new.plt" )
    !call this%tec%convertToBinary()
    !call this%tec%clear()
    
    End Subroutine mapToReferenceDomain

    Subroutine interpolateInto(this, mesh )
        Implicit None 
        Class(TecplotInterpolator)    :: this
        Type(SalomeMeshGeneration)    :: mesh

        Integer                       :: t_int
        Character(len=:), Allocatable :: new
        Character(len=:), Allocatable :: old
        Character(len=:), Allocatable :: out

        new = ".new.plt"
        old = ".old.plt"

        if ( this%output_name == "") then ; out = ".out.plt"
        else                              ; out = this%output_name
        end if 

        call mesh%convertToBinaryTecplot(new, renameCoords = "x_tec y_tec")
        call this%tec%dumpToFile        (old )
        call this%tec%convertToBinary()

        call this%executeScript(new, old, out)

        call execute_command_line( replace("rm -f *","*",new) )
        call execute_command_line( replace("rm -f *","*",old) )
        deallocate( new )
        deallocate( old )
        deallocate( out )
    End Subroutine interpolateInto

    Subroutine setTecplotName(this, filename) 
        Implicit None
        Class(TecplotInterpolator)   :: this 
        Character(len=*), Intent(In) :: filename

        this%output_name = filename
    End Subroutine setTecplotName




    Function getInterpolatedValues(this, getVariableName) Result(sol)
        Implicit None 
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        Interface
            Function getVariableName(id) Result(var)
                Implicit None 
                Integer         , Intent(In)  :: id
                Character(len=:), Allocatable :: var

            End Function getVariableName
        End Interface
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        Class(TecplotInterpolator), Intent(In) :: this
        Real(8), Dimension(:,:), Allocatable   :: sol

        Integer                                :: iukn
        Integer                                :: numberOfUnknowns
        Type(TecplotFileReader)                :: tf

        !*****************************************************************
        ! get the number of unknowns
        !*****************************************************************
        iukn = 1
        do 
            if ( getVariableName(iukn) == "") exit 
            iukn = iukn + 1
        end do
        numberOfUnknowns = iukn - 1

        !*****************************************************************
        ! Read Tecplot
        !*****************************************************************
        tf  = TecplotFileReader( this%output_name )

        allocate( sol( tf%getNumberOfNodes() , numberOfUnknowns) )

        do iukn = 1, numberOfUnknowns
        sol(:,iukn) = tf%getVariable( getVariableName(iukn) )
        end do
        call tf%close()
    End Function getInterpolatedValues


    Function getElements(this) Result(elements)
        Implicit None 
        Class(TecplotInterpolator)           :: this
        Type (TecplotFileReader)             :: tf
        Integer, Dimension(:,:), Allocatable :: elements

        !*****************************************************************
        ! Read Tecplot File 
        !*****************************************************************
        tf  = TecplotFileReader( this%output_name )

        allocate( elements, source = tf%getElements() ) 
        call tf%close()
    End Function getElements


    Subroutine close(this)
     Implicit None 
     Class(TecplotInterpolator) , Intent(InOut)  :: this

     call execute_command_line( replace("rm -f {}","{}", this%output_name) )

     if ( allocated(this%x)          ) deallocate(this%x)
     if ( allocated(this%y)          ) deallocate(this%y)
     if ( allocated(this%elements)   ) deallocate(this%elements)
     if ( allocated(this%varnamelist)) deallocate(this%varnamelist)
     if ( allocated(this%output_name)) deallocate(this%output_name)

    End Subroutine close




    Subroutine executeScript(this, new, old, out)
        Implicit None 
        Class(TecplotInterpolator), Intent(In)    :: this
        Character(len=*)          , Intent(In)    :: new
        Character(len=*)          , Intent(In)    :: old
        Character(len=*)          , Intent(In)    :: out

        Character(len=:)          , Allocatable   :: tecplotscript
        Character(len=:)          , Allocatable   :: new_
        Character(len=:)          , Allocatable   :: old_
        Character(len=:)          , Allocatable   :: out_
        Integer                                   :: t_int

        new_          = concatenate('"',new,'"')
        old_          = concatenate('"',old,'"')
        out_          = concatenate('"',out,'"')
        tecplotscript = ".inter.mcr"

        open ( newunit = t_int, file = tecplotscript)

        write(t_int,"(a)")"#!MC 1410"
        write(t_int,"(a)") replace("$!READDATASET  '* *'","*", [new_, old_])
        write(t_int,"(a)")"  READDATAOPTION = NEW"
        write(t_int,"(a)")"  RESETSTYLE = YES"
        write(t_int,"(a)")"  VARLOADMODE = BYNAME"
        write(t_int,"(a)")"  ASSIGNSTRANDIDS = YES"
        write(t_int,"(a)") replace("  VARLIST =  [*-*]","*",[toStr(3), toStr(this%tec%getNumberOfVariables())])
        write(t_int,"(a)")"$!FIELDLAYERS SHOWMESH = YES"
        write(t_int,"(a)")"$!REDRAWALL"
        write(t_int,"(a)")"$!KRIG"
        write(t_int,"(a)")"  SOURCEZONES =  [2]"
        write(t_int,"(a)")"  DESTINATIONZONE = 1"
        write(t_int,"(a)") replace("  VARLIST =  [*-*]","*",[toStr(3), toStr(this%tec%getNumberOfVariables())])
        write(t_int,"(a)")"  KRIGRANGE = 0.3"
        write(t_int,"(a)")"  KRIGZEROVALUE = 0"
        write(t_int,"(a)")"  KRIGDRIFT = LINEAR"
        write(t_int,"(a)")"  INTERPPTSELECTION = OCTANTNPOINTS"
        write(t_int,"(a)")"  INTERPNPOINTS = 8"
        write(t_int,"(a)") replace("$!WRITEDATASET  *", "*", out_)
        write(t_int,"(a)")"  INCLUDETEXT = NO"
        write(t_int,"(a)")"  INCLUDEGEOM = NO"
        write(t_int,"(a)")"  INCLUDEDATASHARELINKAGE = YES"
        write(t_int,"(a)")"  ZONELIST =  [1]"
        write(t_int,"(a)")"  BINARY = YES"
        write(t_int,"(a)")"  USEPOINTFORMAT = NO"
        write(t_int,"(a)")"  PRECISION = 9"
        write(t_int,"(a)")"  TECPLOTVERSIONTOWRITE = TECPLOTCURRENT"

        call execute_command_line( replace( concatenate(tecplot," -b -p * > /dev/null"), "*", tecplotscript) )
        close(t_int, status='delete')
    End Subroutine executeScript


End Module TecplotInterpolation




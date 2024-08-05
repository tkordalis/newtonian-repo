Module IO_module

    Type Solution
        Real(8)                              :: time
        Real(8), Dimension(:,:), Allocatable :: TL

        Integer                              :: Increment

        contains 
            procedure :: getSolutionVars
    End Type Solution

    contains


    Subroutine getSolutionVars(this, time, sol, increment, RemeshCounter)
        Implicit None 
        Class(Solution), Intent(In)          :: this
        Real(8)                              :: time
        Real(8), Dimension(:,:), Allocatable :: sol
        Real(8)                              :: Pressure_Bubble
        Integer                              :: increment
        Integer, optional                    :: RemeshCounter

        if ( allocated(sol) ) deallocate(sol) 

        allocate( sol , source = this%TL    )
        time              = this%time
        increment         = this%increment
        
        

    End Subroutine getSolutionVars




    Function ReadSolution (filename) Result(this)
        Use VariableMapping 
        Use ReadTecplot
        Use Formats
        Implicit None 
        Character(len=*)       , Intent(In)          :: filename 
        Type(Solution)                               :: this

        Type(TecplotFileReader)                      :: tf
        Character(len=:), Allocatable                :: var
        Character(len=100), Dimension(:),Allocatable :: strGlobal

    
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Open File (tf is (tecplotfile) is the obj that controls the readability
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         tf  = TecplotFileReader(filename)

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Allocate Solution Array
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        allocate( this%TL( tf%getNumberOfNodes(), getNumberOfUnknowns() ) )

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Read Field Variables
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        var = "Vr" ; this%TL(:, getVariableId(var) ) = tf%getVariable(var)
        var = "Vz" ; this%TL(:, getVariableId(var) ) = tf%getVariable(var)
        var = "P"  ; this%TL(:, getVariableId(var) ) = tf%getVariable(var)
        
        var = "Z"  ; this%TL(:, getVariableId(var) ) = tf%getVariable(var)
        var = "R"  ; this%TL(:, getVariableId(var) ) = tf%getVariable(var)
        
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! The Globals are constained in the title of the file 
        ! split the title to extract the variables 
        ! use replace 
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        print*, tf%getTitle()
        strGlobal              = split(tf%getTitle(), sep = ",")
        ! print*, strGlobal(1)
        ! print*, toInt   ( remove( strGlobal(1), ["Increment ", "="]) )
        ! print*, strGlobal(2)
        ! print*, remove( strGlobal(2), ["Remesh_counter ", "="])
        ! ! print*, toInt   ( remove( strGlobal(2), ["Remesh_Counter ", "="]) )
        ! print*, strGlobal(3)
        ! print*, toDouble( remove( strGlobal(3), ["Bubble_Pressure_1 ", "="]) )

        ! pause
        this%increment         = toInt   ( remove( strGlobal(1), ["Increment        ", "="]) )
        

      
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! The zone title (zonename) contains the dimensionless time of the simulation
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        this%time = toDouble(tf%getZonename())

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Close the file
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        call tf%close()

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Deallocate arrays
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        deallocate( strGlobal )
        if (allocated(var)) deallocate( var       )
    End Function ReadSolution

    Subroutine WriteTecplotFile(filename, title, zone, array, elements)
        Use VariableMapping
        Use Tecplot
        Implicit None
        Character(len=*)       , Intent(In) :: filename
        Character(len=*)       , Intent(In) :: title
        Character(len=*)       , Intent(In) :: zone
        Real(8), Dimension(:,:), Intent(In) :: array
        Integer, Dimension(:,:), Intent(In) :: elements

        Type(Tecplot_File)                  :: tecfile
        Integer                             :: ivar
        Integer                             :: unknowns


        unknowns = size(array,2)

        tecfile  = Tecplot_File( datapacking = "POINT", femtype = "FETRIANGLE")

        call tecfile%setTitle   (title)
        call tecfile%setZonename(zone )

        do ivar = 1, unknowns
        call tecfile%addVariable( getVariableName(ivar), array(:,ivar) )
        end do

        call tecfile%addElements( elements )

        call tecfile%dumpToFile( filename )
        call tecfile%convertToBinary()
        call tecfile%clear()
    End Subroutine WriteTecplotFile

    
! 888888888888888888888888888888888888888888888888888888888888888888888888888888888888


    Subroutine exportFiles(Solution, elements, time, Increment, datapack)
        Use Physical_module
        Use Formats
        Use Tecplot
        Use RemeshVariables
        Use FieldFunctions
        Use VariableMapping
        Use MESH_MODULE, only: Xm, Ym
        Use GLOBAL_ARRAYS_MODULE, only: TL
        USE ENUMERATION_MODULE,   only: NM_MESH
        Use BOUNDARY_ENUMERATION_MODULE, only: NBE, &
                                                        bnd1_elements, bnd1_faces, &
                                                        bnd2_elements, bnd2_faces, &
                                                        bnd3_elements, bnd3_faces, &
                                                        bnd4_elements, bnd4_faces, &
                                                        bnd5_elements, bnd5_faces, &
                                                        bnd6_elements, bnd6_faces, &
                                                        bnd7_elements, bnd7_faces, &
                                                        getBoundaryNodesOfTheElement
        Implicit None 
        Real(8), Dimension(:,:), Intent(In) :: Solution
        Integer, Dimension(:,:), Intent(In) :: elements
        Real(8),                 Intent(In) :: time
        Integer,                 Intent(In) :: Increment
        Character(len=*),        Intent(In) :: datapack


        Type(Tecplot_File)                  :: tecfile
        Character(len=:), Allocatable       :: title
        Character(len=:), Allocatable       :: filename
        
        Real(8), Dimension(:,:), Allocatable  :: stresses_nodes
        Integer                             :: i

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Write Tecplot File
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        select case(datapack)
            case("POINT")
                
                tecfile  = Tecplot_File( datapacking = datapack, femtype = "FETRIANGLE")
        
                title = "" 
                title = replace("Increment = *, time = *", "*", &
                       [ toStr(Increment), toStr(time) ] )
                
        
                call tecfile%setTitle   (title)
                call tecfile%setZonename(toStr(time))
                call tecfile%addVariable("Z"      , Solution(:, getVariableId("Z"  )))
                call tecfile%addVariable("R"      , Solution(:, getVariableId("R"  )))
                call tecfile%addVariable("X"      , Xm(:))
                call tecfile%addVariable("Y"      , Ym(:))
                call tecfile%addVariable("Vz"      , Solution(:, getVariableId("Vz"  )))
                call tecfile%addVariable("Vr"      , Solution(:, getVariableId("Vr"  )))
                call tecfile%addVariable("P"      , Solution(:, getVariableId("P"  )))
                

                call tecfile%addElements( elements )
        
                call tecfile%dumpToFile(replace("./tecplot/time_*.plt", "*", toStr(time) ) )
                call tecfile%convertToBinary()
                call tecfile%clear()

            case("BLOCK")

                tecfile  = Tecplot_File( datapacking = datapack, femtype = "FETRIANGLE")

                title = "Mesh Quality"

                call tecfile%setTitle   (title)
                call tecfile%setZonename(toStr(time))

                ! call tecfile%addVariable("Z"      , Solution(:, getVariableId("Z"                ))                       )
                ! call tecfile%addVariable("R"      , Solution(:, getVariableId("R"                ))                       )
                ! call tecfile%addVariable("Stretch", (180d0/3.141599265d0)*minimumAngleOfTriangle(Solution, elements), cellcentered = .True.)

                call tecfile%addElements( elements )

                call tecfile%dumpToFile(replace("./meshq/timeM_*.plt", "*", toStr(time) ) )
                call tecfile%convertToBinary()
                call tecfile%clear()

            case default
                write(*,*) 'Invalid datapack choice'
        end select
    End Subroutine exportFiles


End Module IO_module


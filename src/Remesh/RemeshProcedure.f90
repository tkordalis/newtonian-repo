module RemeshProcedure
    use RemeshVariables

    public  :: checkAndRemesh
    private :: Remesh


    integer :: BeforeRemesh_counter        = 0
    integer :: AfterRemesh_counter         = 0
    integer :: AfterRemeshFading1_counter  = 0
    integer :: AfterRemeshFading2_counter  = 0
    logical :: BeforeRemesh_logical        = .false.
    logical :: AfterRemesh_logical         = .false.
    logical :: AfterRemeshFading1_logical  = .false.
    logical :: AfterRemeshFading2_logical  = .false.

    contains


    Subroutine checkAndRemesh(Solution_, elements, x_mesh, y_mesh, Increment, ReallocateForRemesh )
        Use FieldFunctions,   only : minimumAngleOfTriangle
        Use TIME_INTEGRATION, only : DT, DT_constant
        Implicit none
        Real(8), dimension(:,:), intent(in) :: Solution_
        Integer, Dimension(:,:), Intent(In) :: elements
        Real(8), Dimension(:)  , Intent(In) :: x_mesh
        Real(8), Dimension(:)  , Intent(In) :: y_mesh
        Integer,                 Intent(in) :: Increment
        Logical,                 Intent(out):: ReallocateForRemesh
    
        Real(8), Dimension(:)  , Allocatable  :: Stretch, minimumAnglesAllTriangles
        Real(8)                               :: maxSkewness

        Real(8), parameter  :: theta_degrees_threshold     = 10.d0  ! degrees

        Real(8), parameter  :: remeshThreshold     = (3.1415926535d0/180.d0) * theta_degrees_threshold
        Integer, parameter  :: IterThreshold       = 2
        Integer, parameter  :: IterFading1Threshold = 1
        Integer, parameter  :: IterFading2Threshold = 1

        

        ReallocateForRemesh = .false.
        ! Stretch = relativeElementArea(Solution_, elements, x_mesh, y_mesh)
        minimumAnglesAllTriangles = minimumAngleOfTriangle(Solution_, elements)
        

        if ( BeforeRemesh_logical .eqv. .true. ) then 
            BeforeRemesh_counter = BeforeRemesh_counter + 1
        else
            BeforeRemesh_counter = 0
        endif

        if ( AfterRemesh_logical .eqv. .true. ) then 
            AfterRemesh_counter = AfterRemesh_counter + 1
        else
            AfterRemesh_counter = 0
        endif

        if ( AfterRemeshFading1_logical .eqv. .true. ) then 
            AfterRemeshFading1_counter = AfterRemeshFading1_counter + 1
        else
            AfterRemeshFading1_counter = 0
        endif

        if ( AfterRemeshFading2_logical .eqv. .true. ) then 
            AfterRemeshFading2_counter = AfterRemeshFading2_counter + 1
        else
            AfterRemeshFading2_counter = 0
        endif



        ! maxSkewness = abs(minval(Stretch))
        maxSkewness = abs(minval(minimumAnglesAllTriangles))

        if ( BeforeRemesh_logical .eqv. .True. ) then
            print*, 'Remesh will be performed shortly'
            ! print*, 'Max Skewness =', maxSkewness
            print*, 'Min Angle =', (180d0/3.14159265d0)*maxSkewness
            print*, 'Steps to remesh:', IterThreshold - BeforeRemesh_counter
            DT = Dt_constant
            ! DT = Dt_constant/real(IterThreshold)
        elseif ( AfterRemesh_logical .eqv. .True. ) then
            print*, 'Remesh has already been performed'
            print*, 'Steps after remesh:', AfterRemesh_counter
            DT = Dt_constant
            ! DT = Dt_constant/real(IterThreshold)
        elseif ( (AfterRemesh_logical .eqv. .false.) .and. (AfterRemeshFading1_logical .eqv. .true.) ) then
            print*, 'Restoring timestep to its initial value: Fade #1'
            print*, 'Step increasing timestep:', AfterRemeshFading1_counter
            DT = Dt_constant
            ! DT = Dt_constant/(0.5d0*real(IterThreshold))
        elseif ( (AfterRemeshFading1_logical .eqv. .false.) .and. (AfterRemeshFading2_logical .eqv. .true.) ) then
            print*, 'Restoring timestep to its initial value: Fade #2'
            print*, 'Step increasing timestep:', AfterRemeshFading2_counter
            DT = Dt_constant
        else 
            DT = Dt_constant
        endif


        ! If ( (maxSkewness .gt. remeshThreshold) .and. (BeforeRemesh_counter .eq. 0) ) then
        If ( (maxSkewness .lt. remeshThreshold) .and. (BeforeRemesh_counter .eq. 0) ) then
            BeforeRemesh_logical = .true.
            BeforeRemesh_counter = 1
        endif

        if ( BeforeRemesh_counter .ge. IterThreshold ) then
            BeforeRemesh_logical = .false.
            call Remesh
            ReallocateForRemesh  = .true.
            AfterRemesh_logical  = .true.
            AfterRemesh_counter  = 1
        endif


        if ( AfterRemesh_counter .ge. IterThreshold ) then
            AfterRemesh_logical = .False.
            AfterRemeshFading1_logical = .true.
        endif

        if ( AfterRemeshFading1_counter .ge. IterFading1Threshold ) then
            AfterRemeshFading1_logical = .false.
            AfterRemeshFading2_logical = .true.
        endif

        if ( AfterRemeshFading2_counter .ge. IterFading2Threshold ) then
            AfterRemeshFading2_logical = .false.
        endif
    
    end Subroutine checkAndRemesh
    


    Subroutine Remesh
            Use unv
            Use MeshGeneration
            Use Tecplot
            Use VariableMapping
            Use Formats
            Use ELEMENTS_MODULE
            Use TIME_INTEGRATION
            Use PHYSICAL_MODULE
            Use GLOBAL_ARRAYS_MODULE
            Use BOUNDARY_ENUMERATION_MODULE
            Use MESH_MODULE
            Use ENUMERATION_MODULE
            Use CONTINUATION_MODULE
            Use FLOW_ARRAYS_MODULE
            Use GLOBAL_ARRAYS_MODULE
            Use CSR_STORAGE
            Implicit None
            Character(len=:), Allocatable        :: new_mesh_name
            Type(unvFileReader)                  :: unvf
            Type(SalomeMeshGeneration)           :: new_mesh
            Type(TecplotInterpolator)            :: tecint
            Type(TecplotInterpolator)            :: tecint_o
            Type(TecplotInterpolator)            :: tecint_b
            Real(8), Dimension(:)  , Allocatable :: Z
            Real(8), Dimension(:)  , Allocatable :: R
            Real(8), Dimension(:,:), ALlocatable :: solution
            Real(8), Dimension(:,:), ALlocatable :: solution_o
            Real(8), Dimension(:,:), ALlocatable :: solution_b
            Integer, Dimension(:)  , Allocatable :: bnd3_e, bnd3_f
            Integer, Dimension(:)  , Allocatable :: globnodes
            integer :: ii 

            new_mesh_name = replace("remesh_{}.unv", "{}", toStr(time) )
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            ! Create New Mesh
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            Remesh_counter = Remesh_counter + 1

            if (StructMeshOnlyInTheFront) Remesh_counter_structuredInTheFront = Remesh_counter_structuredInTheFront + 1


            print*, "1. Create New Mesh"
            new_mesh = SalomeMeshGeneration("./mesh/remesh.py")

            
            ! call getBoundaryNodes(TL, bnd3_elements, bnd3_faces, "Z", Z)
            ! call getBoundaryNodes(TL, bnd3_elements, bnd3_faces, "R", R)
            ! call new_mesh%writeBoundaryNodes("Bubble1.dat", Z, R)
            call getVariableFromBoundaryNodes(TL, bnd3_elements, bnd3_faces, "Z", Z)
            call getVariableFromBoundaryNodes(TL, bnd3_elements, bnd3_faces, "R", R)
            call new_mesh%writeBoundaryNodes("Bubble1.dat", Z, R)
            deallocate(Z)
            deallocate(R)
            

            call getVariableFromBoundaryNodes(TL, bnd4_elements, bnd4_faces, "Z", Z)
            call getVariableFromBoundaryNodes(TL, bnd4_elements, bnd4_faces, "R", R)
            call new_mesh%writeBoundaryNodes("Ambient.dat", Z, R)
            deallocate(Z)
            deallocate(R)


            call new_mesh%setUnvFilename( new_mesh_name )
            call new_mesh%execute()

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            ! Intepolate
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            print*, "2. Interpolate Old Mesh to the New one"

            Z      = TL(:,getVariableId("Z"))
            R      = TL(:,getVariableId("R")) 

            tecint   = TecplotInterpolator( Z, R, NM_MESH )
            tecint_o = TecplotInterpolator( Z, R, NM_MESH )
            tecint_b = TecplotInterpolator( Z, R, NM_MESH )

            call tecint  %setTecplotName("solution.plt"  )
            call tecint_o%setTecplotName("solution_o.plt")
            call tecint_b%setTecplotName("solution_b.plt")

            call tecint  %mapToReferenceDomain(time , TL , getVariableName)
            call tecint_o%mapToReferenceDomain(timeo, TLo, getVariableName)
            call tecint_b%mapToReferenceDomain(timeo-dt, TLb, getVariableName)

            print*, "  i. interpolation of the solution"
            call tecint  %interpolateInto(new_mesh)
            print*, " ii. interpolation of the solution at the previous step"
            call tecint_o%interpolateInto(new_mesh)
            print*, "iii. interpolation of the solution two time steps before"
            call tecint_b%interpolateInto(new_mesh)

            print*, " iv. Store Field Variables"
            solution   = tecint  %getInterpolatedValues(getVariableName)
            solution_o = tecint_o%getInterpolatedValues(getVariableName)
            solution_b = tecint_b%getInterpolatedValues(getVariableName)

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            ! Close Files
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            print*, "3. Close Files"
            call new_mesh%close()
            call tecint  %close()
            call tecint_o%close()
            call tecint_b%close()

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            ! Release of Memory
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            print*, "4. Release Memory"

            PHASE_f     = -1           ! release internal memory
            call PARDISO&
               (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, DDUM_f, IDUM_f, IDUM_f,&
                IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, DDUM_f, DDUM_f, ERROR_f)

            call Allocate_Mesh_Arrays                 (.FALSE.)
            call Allocate_Enumeration_Indices         (.FALSE.)
            call Allocate_Continuation_Arrays         (.FALSE.)
            call Allocate_Flow_Arrays                 (.FALSE.)
            call Allocate_Global_Arrays               (.FALSE.)
            call Allocate_CSR_Arrays                  (.FALSE.)
            call Allocate_CSR_Connectivity            (.FALSE.)
            call Allocate_Boundary_Enumeration_Indices(.FALSE.)

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            ! Initialize the problem using the new mesh
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            print*, "5. Reinitialize the problem"
            unvf = unvFileReader( new_mesh_name ) 
            call unvf%info()

            call DISCRETIZATION_DATA( unvf )
            call ALLOCATE_MESH_ARRAYS(.TRUE.,NODTOL)
            call MESH ( unvf )
            
            call ALLOCATE_ENUMERATION_INDICES(.TRUE.,NEL_2d,NBF_1d,NBF_2d,NED_2d)
            call NNUM_ELEMENT( unvf )
            call SURROUNDING_NUMBERING
            call NNUM_f
            call ALLOCATE_BOUNDARY_ENUMERATION_INDICES(.TRUE.,NEL_2d)
            call DEFINE_BOUNDARY_NUMBERING( unvf )

            call unvf%close()

            call CSR_MAIN('NU')
            call CSR_MAIN('CO')
            call INITIALIZE_PARDISO
    
            call ALLOCATE_CONTINUATION_ARRAYS( .TRUE., NUNKNOWNS_f           )
            call ALLOCATE_FLOW_ARRAYS        ( .TRUE., NUNKNOWNS_f, NEX_f    )
            call ALLOCATE_GLOBAL_ARRAYS      ( .TRUE., NODTOL, NEL_2d, NEQ_f )

          !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
          ! Copy Solution to Working Arrays
          !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
          print*, "6. Copy to Working Arrays"
          
          ! move unv file to mesh
          call execute_command_line("mkdir -p mesh", wait = .true.)
          call execute_command_line( replace("mv * mesh","*",new_mesh_name), wait = .true. )

          TL = solution
          TLo= solution_o
          !TLb= solution_b
          TLb= TLo

          deallocate(Z)
          deallocate(R)
          deallocate(solution  )
          deallocate(solution_o)
          !deallocate(solution_b)
          deallocate(new_mesh_name)

    End Subroutine Remesh


end module RemeshProcedure

Module InflatedBubbleBoundary
    Use Boundary_Equations
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use DirichletBoundaries
    use MESH_MODULE, only: Xm, Ym
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    
    Private 

    Public :: InflatedBubble, NewInflatedBubble

    Type InflatedBubble
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes


        Integer                            :: gid

        Real(8)                            :: pressure
        Real(8)                            :: maxRcoordinate
        logical                            :: isMaxRcoordinateLargerThanThreshold

        contains
        procedure :: applyBoundaryConditions
        procedure :: pressureEquation
        procedure :: airFlowrateEquation
        procedure :: setProperties
        procedure :: setPressure
        procedure :: setMaxRcoordinateOfInterfaceAndLogicalOperator
        procedure :: getMaxRcoordinate
        procedure :: getMaxZcoordinate

        final     :: deconstructor
    End Type InflatedBubble

    Contains
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Constructor - Each time we remesh the object is deallocated and recreated
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function NewInflatedBubble( elements, faces) Result(This)
        Implicit None 
        Integer, Dimension(:), Intent(In) :: elements
        Integer, Dimension(:), Intent(In) :: faces

        Type(InflatedBubble)                      :: This

        
        if ( size(elements) /= size(faces) ) Then
            Print*, "[Error] NewInflatedBubble."
            Print*, "Incompartible size of elements and faces"
            Stop
        end if
            

        This%nelem = size(elements)
        ! Memory Allocation of the local variables
        If ( Allocated(This%elements) ) Allocate( This%elements ( This%nelem) )
        If ( Allocated(This%faces   ) ) Allocate( This%faces    ( This%nelem) )

        ! assign the elements
        This%elements = elements
        This%faces    = faces

        call getBoundaryNodesOfWholeBoundary( This%nelem, This%elements, This%faces, This%nodes )

    End Function NewInflatedBubble

    Subroutine setProperties(This, gid)
        Implicit None 
        Class(InflatedBubble)           :: This 
        Integer, Intent(In)     :: gid

        this%gid        = gid
    End Subroutine setProperties

    Subroutine setPressure(This, Pressure)
        Implicit None 
        Class(InflatedBubble)       :: This
        Real(8), Intent(In) :: Pressure

        This%Pressure = Pressure
    End Subroutine setPressure

    Function pressureEquation(this, pressure) Result(output)
        Use CSR_STORAGE,    Only: Ar_f 
        Implicit None
        Class(InflatedBubble) :: this
        Real(8)       :: output

        Real(8)       :: pressure
        
        output = this%pressure - pressure 

        Ar_f(this%gid,:) = 0.d0

        ! call loopOverElements(this%nelem, this%elements, this%faces, this%gid, SurfaceIntegration   ) ! first  constrain

    end Function pressureEquation


    Function airFlowrateEquation(this, time) Result(output)
    use physical_module, only: calculateFlowrate
        Implicit None
        Class(InflatedBubble) :: this
        Real(8)       :: output

        Real(8)       :: time
        Real(8)       :: flowrate
        Real(8)       :: int_u_dS
        
        flowrate = calculateFlowrate(time)
        
        int_u_dS = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, inflowVelocity )

        output = -flowrate + int_u_dS

        call loopOverElements(this%nelem, this%elements, this%faces, this%gid, inflowVelocity   )

    end Function airFlowrateEquation


    Subroutine applyBoundaryConditions(This, FlagNR)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL 
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Implicit None 
        Class(InflatedBubble)   , Intent(In)      :: This
        ! Class(InflatedBubble)   :: This
        Character(len=3), Intent(In)      :: FlagNR 

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_2
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_3
        Integer                              :: iel, element, face
        ! Integer                              :: i
        Integer                              :: ii

        ! clear rows in order to substitute kinematic bc - currently the row for the kinematic is the 
        ! R ( since i apply y-equidistibution in Z )
        call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfResidual)
        If (FlagNR == "NRP") &
        call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfJacobian)

        do iel = 1, this%nelem
            element =This%elements(iel)
            face    =This%faces   (iel) 


            call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

            ! do ii=1,size(TL_,1) 
            !     if (TL_(ii, getVariableId("R")) .gt. 1.15d0) this%isMaxRcoordinateLargerThanThreshold = .true.
            ! enddo

            call Kinematic                        (element, face, TL_, RES_1, .true.)
            
            call Theta_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_2, .true.)
            
            call Stresses                         (element, face, TL_, RES_3, .true., This%pressure )

            
    
            if (FlagNR == "NRP") Then
                call CalculateJacobianContributionsOf(Kinematic                        ,element, face, TL_, RES_1)

                call CalculateJacobianContributionsOf(Theta_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_2)
                
                call CalculateJacobianContributionsOf(Stresses                         ,element, face, TL_, RES_3, This%pressure)
                
                !Extra Unknowns
                call CalculateExtraJacobianContributionsOf(Stresses, element, face, TL_, RES_3, 1, This%pressure,   this%gid)
          end if
        end do 

        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyBoundaryConditions


    subroutine setMaxRcoordinateOfInterfaceAndLogicalOperator(This, StructuredMeshOnlyInTheFront, solution_vector) 
        Use GLOBAL_ARRAYS_MODULE,        Only: TL 
        Implicit none
        Class(InflatedBubble),intent(inout)  :: This
        integer                :: inode, node
        real(8), dimension(:,:), intent(in)    :: solution_vector
        logical, intent(inout)   :: StructuredMeshOnlyInTheFront


        StructuredMeshOnlyInTheFront = .false.

        do inode = 1,size(this%nodes)
            node = this%nodes(inode)
            ! print*, TL(node, getVariableId("R"))
            if (solution_vector(node, getVariableId("R")) .gt. 1.15d0) then 
                print*, 'global node and R coordinate =',node, solution_vector(node, getVariableId("R"))
                StructuredMeshOnlyInTheFront = .true.
                 print*, StructuredMeshOnlyInTheFront
                exit
            endif
        enddo
        
    end subroutine setMaxRcoordinateOfInterfaceAndLogicalOperator

    function getMaxRcoordinate(This)    result(output)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL 
        Implicit none
        Class(InflatedBubble),intent(inout)  :: This
        integer                :: inode, node
        real(8)                :: maxRcoordinate, output

        maxRcoordinate = maxval( TL(this%nodes,getVariableId("R")) )

        output = maxRcoordinate
        
    end function getMaxRcoordinate

    function getMaxZcoordinate(This)    result(output)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL 
        Implicit none
        Class(InflatedBubble),intent(inout)  :: This
        integer                :: inode, node
        real(8)                :: maxZcoordinate, output

        maxZcoordinate = maxval( TL(this%nodes,getVariableId("Z")) )

        output = maxZcoordinate
        
    end function getMaxZcoordinate


    function getPressureAtMaxZcoordinate(This)    result(output)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL 
        Implicit none
        Class(InflatedBubble),intent(inout)  :: This
        integer                :: inode, node
        real(8)                :: maxlocZcoordinate, output

        maxlocZcoordinate = maxloc( TL(this%nodes,getVariableId("Z")),1 )
        print*, TL(this%nodes(maxlocZcoordinate),getVariableId("Z"))

        output = maxlocZcoordinate
        
    end function getPressureAtMaxZcoordinate


    Subroutine deconstructor(This)
        Implicit None
        Type(InflatedBubble) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

end module InflatedBubbleBoundary
Module AmbientBoundary
    Use Boundary_Equations
    Use NumericalBoundaryJacobian 
    Use ExtraEquations
    Use DirichletBoundaries,         only: ApplyDirichletAtNode_
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    use MESH_MODULE, only: Xm, Ym
    use physical_module, only: initial_position, ambient_position


    private 

    public :: Ambient, NewAmbient

    Type Ambient
        ! Wall Properties
        
        Integer                            :: nelem 
        Integer, Dimension(:), Allocatable :: elements 
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes
        Real(8)                            :: minimumZcoord


        Real(8) :: datumPressure

        contains 
            procedure :: applyBoundaryConditions
            procedure :: setDatumPressure
            final     :: deconstructor
    End Type Ambient    


    contains

    Function NewAmbient( elements, faces) Result(This)
        Implicit None 
        Type(Ambient)                      :: This 
        Integer, Dimension(:), Allocatable :: elements 
        Integer, Dimension(:), Allocatable :: faces 
        real(8)                            :: a1 


        
        This%nelem = size(elements)
        ! Memory Allocation of the local variables
        If ( Allocated(This%elements) ) Allocate( This%elements ( This%nelem) )
        If ( Allocated(This%faces   ) ) Allocate( This%faces    ( This%nelem) )

        ! assign the elements
        This%elements = elements 
        This%faces    = faces

        call getBoundaryNodesOfWholeBoundary( This%nelem, This%elements, This%faces, This%nodes )

        ! this%minimumZcoord = Xm(   this%nodes(  minloc( Xm(this%nodes), dim=1 )  )   )
        ! this%maximumRcoord = Ym(   this%nodes(  maxloc( Ym(this%nodes), dim=1 )  )   )

        this%minimumZcoord = minval(  Xm( this%nodes )  )

        ambient_position = this%minimumZcoord 
        initial_position = this%minimumZcoord 

    End Function NewAmbient




    Subroutine applyBoundaryConditions(This, FlagNr)
        Use physical_module,             only: ambient_position, vm_ambient
        Use TIME_INTEGRATION,            only: time, dt
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ksi => Xm, Eta => Ym
        ! Use DirichletBoundaries, only : updateAllNodesOfTheBoundary, &
        !                             ClearRowsOfResidual, ClearRowsOfJacobian

        Implicit None 
        Class(Ambient)  , Intent(In)         :: This 
        Character(len=3), Intent(In)         :: FlagNr

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1, RES_2, RES_3
        Integer                              :: iel, node_counter, node
        Integer                              :: element 
        Integer                              :: face
        Integer                              :: i
        Integer                              :: ii
        Integer                              :: datumPressure_boundary_node
       

        
        do node_counter = 1, size(this%nodes)
            node = this%nodes(node_counter)

            call ApplyDirichletAtNode_(node, "Z", ambient_position, FlagNr )
            call ApplyDirichletAtNode_(node, "R", Eta(node), FlagNr )
            call ApplyDirichletAtNode_(node, "Vz", vm_ambient, FlagNr )
        enddo

        node = this%nodes(  maxloc( Ym(this%nodes), dim=1 )  )
        call ApplyDirichletAtNode_(node, "P", This%datumPressure, FlagNr )

        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyBoundaryConditions



    Subroutine setDatumPressure(This, datumPressure)
        Implicit None 
        Class(Ambient)       :: This
        Real(8), Intent(In) :: datumPressure

        This%datumPressure = datumPressure
    End Subroutine setDatumPressure



    Subroutine deconstructor(This) 
        Implicit None
        Type(Ambient) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module AmbientBoundary

Module AmbientHenryBoundary
    Use Boundary_EquationsDO
    Use NumericalBoundaryJacobian 
    Use ExtraEquations
    Use DirichletBoundaries,         only: ApplyDirichletAtNode_
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    use MESH_MODULE, only: Xm, Ym
    use physical_module, only: initial_position, position_o, position, ambient_position_o, ambient_position


    private 

    public :: AmbientHenry, NewAmbientHenry

    Type AmbientHenry
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
    End Type AmbientHenry    


    contains

    Function NewAmbientHenry( elements, faces) Result(This)
        Implicit None 
        Type(AmbientHenry)                      :: This 
        Integer, Dimension(:), Allocatable :: elements 
        Integer, Dimension(:), Allocatable :: faces 

        
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

        initial_position   = this%minimumZcoord  ;  position_o = initial_position 
        position           = position_o          ;  ambient_position_o = position 
        ambient_position   = ambient_position_o

    End Function NewAmbientHenry

    Subroutine applyBoundaryConditions(This, FlagNr, dVtankdt, Pressure_bc)
        Use physical_module,             only: ambient_position_o, vm_ambient, Rtank, pi, KoN
        Use TIME_INTEGRATION,            only: dt
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ym
        Use DirichletBoundaries, only : updateAllNodesOfTheBoundary, &
                                    ClearRowsOfResidual, ClearRowsOfJacobian

        Implicit None 
        Class(AmbientHenry)  , Intent(In)         :: This
        Character(len=3), Intent(In)         :: FlagNr
        Real(8),          Intent(In)         :: dVtankdt
        Real(8),          Intent(In)         :: Pressure_bc

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Integer                              :: node_counter, node

        ! call updateAllNodesOfTheBoundary('C',This%elements, This%faces, ClearRowsOfResidual)
        ! If (FlagNR == "NRP") &
        !     call updateAllNodesOfTheBoundary('C',This%elements, This%faces, ClearRowsOfJacobian)
        
        vm_ambient = dVtankdt / (pi*Rtank**2)
        ! vm_ambient = 0.d0
       
        do node_counter = 1, size(this%nodes)
            node = this%nodes(node_counter)
            call ApplyDirichletAtNode_(node, "Z", ambient_position_o + dt*vm_ambient, FlagNr )
            call ApplyDirichletAtNode_(node, "R", Ym(node), FlagNr )
            call ApplyDirichletAtNode_(node, "Vz", vm_ambient, FlagNr )
            
            call ApplyDirichletAtNode_(node, "P", Pressure_bc, FlagNr )
            ! call ApplyDirichletAtNode_(node, "C", KoN*Pressure_bc, FlagNr )
        enddo

        ! weak imposition of henry's law
        ! do iel = 1, this%nelem
        !     element =This%elements(iel)
        !     face    =This%faces   (iel)

        !     call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)
        !     call Henry                        (element, face, TL_, RES_1, .true., Pressure_bc)
        !     if (FlagNR == "NRP") &
        !             call CalculateJacobianContributionsOf(Henry  ,element, face, TL_, RES_1, Pressure_bc)
        !             ! call CalculateJacobianContributionsOf(Henry  ,element, face, TL_, RES_1, Pressure_bc, .true.)
        ! enddo

        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyBoundaryConditions


    Subroutine setDatumPressure(This, datumPressure)
        Implicit None 
        Class(AmbientHenry)       :: This
        Real(8), Intent(In) :: datumPressure

        This%datumPressure = datumPressure
    End Subroutine setDatumPressure



    Subroutine deconstructor(This) 
        Implicit None
        Type(AmbientHenry) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module AmbientHenryBoundary

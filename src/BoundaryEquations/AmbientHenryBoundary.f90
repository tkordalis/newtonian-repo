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


        Real(8) :: datumPressure, concentrationPressure

        contains 
            procedure :: applyBoundaryConditions
            procedure :: setDatumPressure
            procedure :: setConcentrationPressure
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
        this%minimumZcoord = minval(  Xm( this%nodes )  )

        initial_position   = this%minimumZcoord  ;  position_o = initial_position 
        position           = position_o          ;  ambient_position_o = position 
        ambient_position   = ambient_position_o

    End Function NewAmbientHenry

    Subroutine applyBoundaryConditions(This, FlagNr, Pressure_bc, naturalBCs)
        Use physical_module,             only: ambient_position_o, vm_ambient, Rtank, pi, KoN
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH

        Use TIME_INTEGRATION,            only: dt
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ym
        Use DirichletBoundaries, only : updateAllNodesOfTheBoundary, &
                                    ClearRowsOfResidual, ClearRowsOfJacobian

        Implicit None 
        Class(AmbientHenry)  , Intent(In)         :: This
        Character(len=3), Intent(In)         :: FlagNr
        Real(8),          Intent(In)         :: Pressure_bc
        logical,          Intent(In)         :: naturalBCs

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_stresses, RES_kinematic
        Integer                              :: node_counter, node
        Integer                              :: iel, element, face

        if (naturalBCs) then
            do iel = 1, this%nelem
                element =This%elements(iel)
                face    =This%faces   (iel)

                call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

                call Stresses                         (element, face, TL_, RES_stresses, .true., Pressure_bc )

                if (FlagNR == "NRP") &
                    call CalculateJacobianContributionsOf(Stresses      ,element, face, TL_, RES_stresses, Pressure_bc )
            enddo

        else
            call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfResidual)
            If (FlagNR == "NRP") then
                call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfJacobian)
            endif
            
            do iel = 1, this%nelem
                element =This%elements(iel)
                face    =This%faces   (iel)

                call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

                call Kinematic        (element, face, TL_, RES_kinematic, .true. )
                
                If (FlagNR == "NRP") &
                    call CalculateJacobianContributionsOf(Kinematic     ,element, face, TL_, RES_kinematic )
            enddo
           
            do node_counter = 1, size(this%nodes)
                node = this%nodes(node_counter)
                call ApplyDirichletAtNode_(node, "R", Ym(node), FlagNr )
                
                call ApplyDirichletAtNode_(node, "C", KoN*this%concentrationPressure, FlagNr )
            enddo
            
        endif

        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyBoundaryConditions


    Subroutine setDatumPressure(This, datumPressure)
        Implicit None 
        Class(AmbientHenry)       :: This
        Real(8), Intent(In) :: datumPressure

        This%datumPressure = datumPressure
    End Subroutine setDatumPressure

    Subroutine setConcentrationPressure(This, ConcentrationPressure)
        Implicit None 
        Class(AmbientHenry)       :: This
        Real(8), Intent(In) :: ConcentrationPressure

        This%ConcentrationPressure = ConcentrationPressure
    End Subroutine setConcentrationPressure




    Subroutine deconstructor(This) 
        Implicit None
        Type(AmbientHenry) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module AmbientHenryBoundary


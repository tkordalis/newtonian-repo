Module SymmetryDiffusionBoundary
    Use Boundary_EquationsDO 
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use DirichletBoundaries, only : ApplyDirichletAtNode_, updateAllNodesOfTheBoundary, &
                                    ClearRowsOfResidual, ClearRowsOfJacobian
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    use MESH_MODULE, only: Xm, Ym


    private 

    public :: SymmetryDiffusion, NewSymmetryDiffusion

    Type SymmetryDiffusion
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes

        Character(1)                       :: position

        contains 
            procedure :: applyBoundaryConditions
            procedure :: setPosition

            final     :: deconstructor
    End Type SymmetryDiffusion


    contains

    Function NewSymmetryDiffusion( elements, faces) Result(This)
        Implicit None 
        Type(SymmetryDiffusion)                     :: This 
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
    End Function NewSymmetryDiffusion


    Subroutine setPosition(This, pos)
        Implicit None 
        Class(SymmetryDiffusion)              :: This
        character(*), Intent(In)     :: pos

        this%position        = pos
    End Subroutine setPosition


    Subroutine applyBoundaryConditions(This, FlagNr, naturalBCs)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ym

        Implicit None 
        Class(SymmetryDiffusion) , Intent(In)         :: This 
        Character(len=3), Intent(In)         :: FlagNr
        logical,          Intent(In)         :: naturalBCs


        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_Xequid, RES_concentration
        Integer                              :: inode, node
        Integer                              :: iel, element, face

        
        if (naturalBCs) then
            do iel = 1, this%nelem
                element =This%elements(iel)
                face    =This%faces   (iel)

                call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

                call zeroConcentrationFlux  ( element, face, TL_, RES_concentration, .true. )

                if (FlagNR == "NRP") &
                    call CalculateJacobianContributionsOf( zeroConcentrationFlux  ,element, face, TL_, RES_concentration )
            enddo
        else

            call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfResidual)
            If (FlagNR == "NRP") &
                    call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfJacobian)
                    
            do iel = 1, this%nelem
                element =This%elements(iel)
                face    =This%faces   (iel)
                call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

                select case(this%position)
                case('X')
                    call X_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_Xequid, .true.)
                    if (FlagNR == "NRP") Then
                        call CalculateJacobianContributionsOf(X_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_Xequid)
                    end if
                case default
                        Print*, "[Error] : Equid. in symmetryDiffusion type wrong value of position."
                End Select
            end do
            

            do inode = 1, size(this%nodes)
                node = this%nodes(inode)
                call ApplyDirichletAtNode_(node, "Vr", 0.d0      , FlagNr )
                Call ApplyDirichletAtNode_(node, 'R'  , 0.d0     , FlagNr )
            end do
        endif

        If ( Allocated(TL_) ) Deallocate(TL_)

    End Subroutine applyBoundaryConditions


    Subroutine deconstructor(This)
        Implicit None
        Type(SymmetryDiffusion) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

end module SymmetryDiffusionBoundary

Module SymmetryBoundary
    Use Boundary_Equations 
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use DirichletBoundaries,         only: ApplyDirichletAtNode_
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    use MESH_MODULE, only: Xm, Ym


    private 

    public :: Symmetry, NewSymmetry

    Type Symmetry
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes

        Character(1)                       :: position

        contains 
            procedure :: applyBoundaryConditions
            procedure :: setPosition

            final     :: deconstructor
    End Type Symmetry


    contains

    Function NewSymmetry( elements, faces) Result(This)
        Implicit None 
        Type(Symmetry)                     :: This 
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
    End Function NewSymmetry


    Subroutine setPosition(This, pos)
        Implicit None 
        Class(Symmetry)              :: This
        character(*), Intent(In)     :: pos

        this%position        = pos
    End Subroutine setPosition


    Subroutine applyBoundaryConditions(This, FlagNr)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ksi => Xm, Eta => Ym

        Implicit None 
        Class(Symmetry) , Intent(In)         :: This 
        Character(len=3), Intent(In)         :: FlagNr

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Integer                              :: iel, inode, node
        Integer                              :: element 
        Integer                              :: face

        ! do iel = 1, this%nelem
        !     element =This%elements(iel)
        !     face    =This%faces   (iel)
        !     call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

        !     select case(this%position)
        !     case('X')
        !         call X_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_1, .true.)
        !         if (FlagNR == "NRP") Then
        !             call CalculateJacobianContributionsOf(X_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_1)
        !         end if
        !     case('Y')
        !         call Y_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_1, .true.)
        !         if (FlagNR == "NRP") Then
        !             call CalculateJacobianContributionsOf(Y_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_1)
        !         end if
        !     case default
        !             Print*, "[Error] : Equid. in symmetry type wrong value of position."
        !     End Select
        ! end do
        

        do inode = 1, size(this%nodes)
            node = this%nodes(inode)
            call ApplyDirichletAtNode_(node, "Vr", 0.d0      , FlagNr )
            ! Call ApplyDirichletAtNode_(node, 'R'  , 0.d0     , FlagNr )
            Call ApplyDirichletAtNode_(node, 'R'  , Eta(node), FlagNr )
        end do

        If ( Allocated(TL_) ) Deallocate(TL_)

    End Subroutine applyBoundaryConditions


    Subroutine deconstructor(This)
        Implicit None
        Type(Symmetry) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

end module SymmetryBoundary

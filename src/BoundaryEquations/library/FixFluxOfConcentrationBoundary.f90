Module FixFluxOfConcentrationBoundary
    Use Boundary_Equations
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use DirichletBoundaries,         only: ApplyDirichletAtNode_
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    use MESH_MODULE, only: Xm, Ym
    ! use physical_module, only: Rtank
    
    private 

    public :: FixFluxConcentrationValue, NewFixFluxConcentrationValue

    Type FixFluxConcentrationValue
        ! Wall Properties
        
        Integer                            :: nelem 
        Integer, Dimension(:), Allocatable :: elements 
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes
        Real(8)                            :: fluxValueAtTheBoundary
        contains 
            procedure :: setFluxValueAtTheBoundary
            procedure :: applyBoundaryConditions
            final     :: deconstructor
    End Type FixFluxConcentrationValue    


    contains

    Function NewFixFluxConcentrationValue( elements, faces) Result(This)
        Implicit None 
        Type(FixFluxConcentrationValue)                      :: This 
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

    End Function NewFixFluxConcentrationValue


    subroutine setFluxValueAtTheBoundary(this, value)
        implicit none
        Class(FixFluxConcentrationValue)  , Intent(InOut)         :: This 
        Real(8),                        Intent(In)         :: value


        this%fluxValueAtTheBoundary = value
    end subroutine setFluxValueAtTheBoundary




    Subroutine applyBoundaryConditions(This, FlagNr)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ksi => Xm, Eta => Ym

        Implicit None 
        Class(FixFluxConcentrationValue)  , Intent(In)         :: This 
        Character(len=3), Intent(In)         :: FlagNr

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Integer                              :: iel, inode, node
        Integer                              :: element 
        Integer                              :: face
        Integer                              :: i
        Integer                              :: ii

    
        do iel = 1, this%nelem
            element =This%elements(iel)
            face    =This%faces   (iel) 

            call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)
            call fixConcentrationFlux       (element, face, TL_, RES_1, .true., this%fluxValueAtTheBoundary )
            if (FlagNR == "NRP") Then
                call CalculateJacobianContributionsOf(fixConcentrationFlux      ,element, face, TL_, RES_1, this%fluxValueAtTheBoundary )
            end if
        end do 

        If ( Allocated(TL_) ) Deallocate(TL_)

    End Subroutine applyBoundaryConditions




    Subroutine deconstructor(This) 
        Implicit None
        Type(FixFluxConcentrationValue) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module FixFluxOfConcentrationBoundary

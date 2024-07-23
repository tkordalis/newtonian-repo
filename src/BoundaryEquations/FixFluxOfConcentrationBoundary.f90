Module FixConcentrationValueBoundary
    Use Boundary_Equations
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use DirichletBoundaries,         only: ApplyDirichletAtNode_
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    use MESH_MODULE, only: Xm, Ym
    ! use physical_module, only: Rtank
    
    private 

    public :: FixConcentrationValue, NewFixConcentrationValue

    Type FixConcentrationValue
        ! Wall Properties
        
        Integer                            :: nelem 
        Integer, Dimension(:), Allocatable :: elements 
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes
        Real(8)                            :: maximumRcoord
        contains 
            procedure :: applyBoundaryConditions
            final     :: deconstructor
    End Type FixConcentrationValue    


    contains

    Function NewFixConcentrationValue( elements, faces) Result(This)
        Implicit None 
        Type(FixConcentrationValue)                      :: This 
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

        this%maximumRcoord = maxval(Ym(this%nodes))

        ! Rtank = this%maximumRcoord

    End Function NewFixConcentrationValue



    Subroutine applyBoundaryConditions(This, value, FlagNr)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ksi => Xm, Eta => Ym

        Implicit None 
        Class(FixConcentrationValue)  , Intent(In)         :: This 
        Real(8),          Intent(In)         :: value
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
            call fixConcentrationFlux       (element, face, TL_, RES_1, .true.)
            if (FlagNR == "NRP") Then
                call CalculateJacobianContributionsOf(fixConcentrationFlux      ,element, face, TL_, RES_1)
            end if
        end do 

        If ( Allocated(TL_) ) Deallocate(TL_)

    End Subroutine applyBoundaryConditions




    Subroutine deconstructor(This) 
        Implicit None
        Type(FixConcentrationValue) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module FixConcentrationValueBoundary

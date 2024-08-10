Module FixWallConcentrationBoundary
    Use Boundary_EquationsDO
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use DirichletBoundaries,         only: ApplyDirichletAtNode_
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    use MESH_MODULE, only: Xm, Ym
    use physical_module, only: Rtank
    
    private 

    public :: FixWallConcentration, NewFixWallConcentration

    Type FixWallConcentration
        ! Wall Properties
        
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes
        Real(8)                            :: maximumRcoord
        contains 
            procedure :: applyBoundaryConditions
            final     :: deconstructor
    End Type FixWallConcentration    


    contains

    Function NewFixWallConcentration( elements, faces) Result(This)
        Implicit None 
        Type(FixWallConcentration)                      :: This 
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

        Rtank = this%maximumRcoord

    End Function NewFixWallConcentration



    Subroutine applyBoundaryConditions(This, FlagNr)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ksi => Xm, Eta => Ym

        Implicit None 
        Class(FixWallConcentration)  , Intent(In)         :: This 
        Character(len=3), Intent(In)         :: FlagNr

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Integer                              :: iel, inode, node
        Integer                              :: element 
        Integer                              :: face
        Integer                              :: i
        Integer                              :: ii

    
        do inode = 1, size(this%nodes)
            node = this%nodes(inode)
            call ApplyDirichletAtNode_(node, "Vz", 0.d0, FlagNr )
            call ApplyDirichletAtNode_(node, "Vr", 0.d0, FlagNr )
            call ApplyDirichletAtNode_(node, "Z", Ksi(node), FlagNr )
            call ApplyDirichletAtNode_(node, "R", Eta(node), FlagNr )

            call ApplyDirichletAtNode_(node, "C", 1.d0, FlagNr )
        end do 

    End Subroutine applyBoundaryConditions




    Subroutine deconstructor(This) 
        Implicit None
        Type(FixWallConcentration) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module FixWallConcentrationBoundary

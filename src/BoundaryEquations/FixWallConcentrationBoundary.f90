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



    Subroutine applyBoundaryConditions(This, FlagNr, naturalBCs)
        Use MESH_MODULE,                 only: Xm, Ym
        Implicit None 
        Class(FixWallConcentration)  , Intent(In)         :: This 
        Character(len=3), Intent(In)         :: FlagNr
        logical,          Intent(In)         :: naturalBCs


        Integer                              :: inode, node

        if (naturalBCs) then


        else
            do inode = 1, size(this%nodes)
                node = this%nodes(inode)
                call ApplyDirichletAtNode_(node, "Vz", 0.d0, FlagNr )
                call ApplyDirichletAtNode_(node, "Vr", 0.d0, FlagNr )
                call ApplyDirichletAtNode_(node, "Z", Xm(node), FlagNr )
                call ApplyDirichletAtNode_(node, "R", Ym(node), FlagNr )
            end do
        endif 

    End Subroutine applyBoundaryConditions




    Subroutine deconstructor(This) 
        Implicit None
        Type(FixWallConcentration) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module FixWallConcentrationBoundary

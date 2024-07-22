Module MovingWallBoundary
    Use Boundary_Equations
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use DirichletBoundaries,         only: ApplyDirichletAtNode_
    use boundary_enumeration_module, only: getBoundaryNodesOfWholeBoundary
    use MESH_MODULE, only: Xm, Ym
    ! use physical_module, only: Rtank
    
    private 

    public :: MovingWall, NewMovingWall

    Type MovingWall
        ! Wall Properties
        
        Integer                            :: nelem 
        Integer, Dimension(:), Allocatable :: elements 
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes
        Real(8)                            :: maximumRcoord
        contains 
            procedure :: applyBoundaryConditions
            procedure :: getLidVelocity
            final     :: deconstructor
    End Type MovingWall    


    contains

    Function NewMovingWall( elements, faces) Result(This)
        Implicit None 
        Type(MovingWall)                      :: This 
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

    End Function NewMovingWall



    Subroutine applyBoundaryConditions(This, FlagNr, time)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Use MESH_MODULE,                 only: Ksi => Xm, Eta => Ym

        Implicit None 
        Class(MovingWall)  , Intent(In)         :: This 
        Character(len=3),   Intent(In)         :: FlagNr
        Real(8),            Intent(In)         :: time

        Real(8)                              :: Ulid
        Integer                              :: iel, inode, node
        Integer                              :: element 
        Integer                              :: face
        Integer                              :: i
        Integer                              :: ii


        do inode = 1, size(this%nodes)
            node = this%nodes(inode)
            Ulid = this%getLidVelocity( Ksi(node), time )
            call ApplyDirichletAtNode_(node, "Vz", Ulid, FlagNr )
            call ApplyDirichletAtNode_(node, "Vr", 0.d0, FlagNr )
            call ApplyDirichletAtNode_(node, "Z", Ksi(node), FlagNr )
            call ApplyDirichletAtNode_(node, "R", Eta(node), FlagNr )
            
        end do 

    End Subroutine applyBoundaryConditions

    function getLidVelocity(This, Xnode,time) result(output)
        implicit none
        Class(MovingWall) :: This
        real(8), intent(in) :: Xnode
        real(8), intent(in) :: time
        real(8)             :: output

        output = ( 1.d0 - (2.d0*Xnode-1.d0)**14.d0 )
        ! print*, time, Xnode, output
        ! pause
    end function getLidVelocity


    Subroutine deconstructor(This) 
        Implicit None
        Type(MovingWall) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module MovingWallBoundary

Module InternalEquidistribution
    Use Boundary_Equations
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    
    Private 

    Public :: IntEquidistribution, NewIntEquidistribution

    Type IntEquidistribution
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces

        Character(1)                       :: position

        contains
        procedure :: applyEquidistribution
        procedure :: setPosition
        
        final     :: deconstructor
    End Type IntEquidistribution

    Contains
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Constructor - Each time we remesh the object is deallocated and recreated
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function NewIntEquidistribution( elements, faces) Result(This)
        Implicit None 
        Integer, Dimension(:), Intent(In) :: elements
        Integer, Dimension(:), Intent(In) :: faces

        Type(IntEquidistribution)                      :: This

        
        if ( size(elements) /= size(faces) ) Then
            Print*, "[Error] NewIntEquidistribution."
            Print*, "Incompartible size of elements and faces"
            Stop
        end if
            

        This%nelem = size(elements)
        ! Memory Allocation of the local variables
        If ( Allocated(This%elements) ) Allocate( This%elements ( This%nelem) )
        If ( Allocated(This%faces   ) ) Allocate( This%faces    ( This%nelem) )

        ! assign the elements
        This%elements = elements
        This%faces    = faces

    End Function NewIntEquidistribution

    Subroutine setPosition(This, pos)
        Implicit None 
        Class(IntEquidistribution)   :: This 
        character(*), Intent(In)     :: pos

        this%position        = pos
    End Subroutine setPosition


    Subroutine applyEquidistribution(This, FlagNR)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL 
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Implicit None 
        Class(IntEquidistribution)   , Intent(In)   :: This
        Character(len=3), Intent(In)           :: FlagNR 

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_2
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_3
        Integer                              :: iel, element, face

        do iel = 1, this%nelem
            element =This%elements(iel)
            face    =This%faces   (iel)
            call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

            select case(this%position)
            case('X')
                call X_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_2, .true.)
                if (FlagNR == "NRP") Then
                    call CalculateJacobianContributionsOf(X_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_2)
                end if
            case('Y')
                call Theta_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_2, .true.)
                if (FlagNR == "NRP") Then
                    call CalculateJacobianContributionsOf(Theta_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_2)
                end if
                
                ! call Y_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_2, .true.)
                ! if (FlagNR == "NRP") Then
                !     call CalculateJacobianContributionsOf(Y_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_2)
                ! end if
            ! case('O') ! Other than X or Y
            !     call Theta_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_2, .true.)
            !     if (FlagNR == "NRP") Then
            !         call CalculateJacobianContributionsOf(Theta_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_2)
            !     end if
            case default
                    Print*, "[Error] : Equid. type wrong value of position."
            End Select
        end do 
        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyEquidistribution

    Subroutine deconstructor(This)
        Implicit None
        Type(IntEquidistribution) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End module InternalEquidistribution


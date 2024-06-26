!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><2
! Module BubbleBoundary 
! Author : Dionisis Pettas
! Date   : 16.10.2020
! 
! Boundary Conditions Clustering for Bubble Definition with Cylindical 
! coordinates
!
! Future Idea Create an Abstract Class to inherit boundary conditions
! in the subclasses
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Module BubbleBoundary
    Use Boundary_Equations
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use DirichletBoundaries
    Private 

    Public :: Bubble, NewBubble


    Type Bubble 
        ! Bubble Properties
        Real(8)                            :: Rcenter = 0.d0
        Real(8)                            :: Zcenter

        ! FEM : Prop
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces
        ! Global Contrain
        Integer                            :: gid
        Real(8)                            :: pressure, InitialPressure
        Real(8)                            :: volume  , InitialVolume
        contains 
            ! constrains and boundary conditions
            procedure :: applyBoundaryConditions
            procedure :: volumeConservation
            procedure :: fixCentroid
            procedure :: PressureVolumeConservation
            ! setters
            procedure :: setProperties
            procedure :: setInitialPressure
            procedure :: setPressure
            procedure :: setInitialVolume
            procedure :: setVolume
            procedure :: setCentroid
            ! getters
            procedure :: getCentroid
            procedure :: getDragForce
            procedure :: getAspectRatio
            procedure :: getVolume

            final     :: deconstructor
    End Type Bubble



    Contains
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Constructor - Each time we remesh the object is deallocated and recreated
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function NewBubble( elements, faces) Result(This)
        Implicit None 
        Integer, Dimension(:), Intent(In) :: elements
        Integer, Dimension(:), Intent(In) :: faces

        Type(Bubble)                      :: This

        
        if ( size(elements) /= size(faces) ) Then
            Print*, "[Error] NewBubble."
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
    End Function NewBubble

    Subroutine setProperties(This, gid)
        Implicit None 
        Class(Bubble)           :: This 
        Integer, Intent(In)     :: gid

        this%gid        = gid
    End Subroutine setProperties
    Subroutine setInitialPressure(This, InitialPressure)
        Implicit None 
        Class(Bubble)       :: This
        Real(8), Intent(In) :: InitialPressure

        This%InitialPressure = InitialPressure
    End Subroutine setInitialPressure
    Subroutine setPressure(This, Pressure)
        Implicit None 
        Class(Bubble)       :: This
        Real(8), Intent(In) :: Pressure

        This%Pressure = Pressure
    End Subroutine setPressure
    Subroutine setInitialVolume(This, InitialVolume)
        Implicit None 
        Class(Bubble)       :: This
        Real(8), Intent(In) :: InitialVolume

        This%InitialVolume = InitialVolume
    End Subroutine setInitialVolume
    Subroutine setVolume(This)
        Implicit None 
        Class(Bubble)       :: This
        
        This%Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
    End Subroutine setVolume
    Subroutine setCentroid(this)
        Implicit None 
        Class(Bubble)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Centroid, Volume
        
        Volume   = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        Centroid = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Centroid = Centroid/Volume

        this%Zcenter = Centroid
    end Subroutine setCentroid
    
    
    
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                       volumeConservation                       
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function volumeConservation(this) Result(output)
        Implicit None
        Class(Bubble) :: this
        Real(8)       :: output

        Real(8)       :: Volume
        

        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )

        output = Volume - this%InitialVolume

        call loopOverElements(this%nelem, this%elements, this%faces, this%gid, SurfaceIntegration   ) ! first  constrain
    end Function volumeConservation


    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                   PressureVolumeConservation                   
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function PressureVolumeConservation(this) Result(output)
        Implicit None
        Class(Bubble) :: this
        Real(8)       :: output

        Real(8)       :: Volume
        

        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        
        output = this%pressure * Volume - this%InitialPressure * this%InitialVolume

        call loopOverElements(this%nelem, this%elements, this%faces, this%gid, SurfaceIntegration, this%pressure ) ! first  constrain
    end Function PressureVolumeConservation


    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                          fixCentroid                          
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function fixCentroid(this) Result(output)
        Use constrainJacobians, only: jacobianOfConstrainCentroid
        Implicit None 
        Class(Bubble)               :: this
        Real(8)                     :: output

        Real(8)                     :: Volume

        Real(8)                     :: Centroid

        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )

        
        Centroid = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Centroid = Centroid/Volume

        call jacobianOfConstrainCentroid(this%gid)

        output = Centroid
    End Function fixCentroid



    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                    applyBoundaryConditions                     
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Subroutine applyBoundaryConditions(This, FlagNR)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Implicit None 
        Class(Bubble)   , Intent(In)      :: This
        Character(len=3), Intent(In)      :: FlagNR 

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_2
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_3
        Integer                              :: iel 
        Integer                              :: element 
        Integer                              :: face
        Integer                              :: i
        Integer                              :: ii

        call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfResidual)
        If (FlagNR == "NRP") &
        call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfJacobian) 

        do iel = 1, this%nelem
            element =This%elements(iel)
            face    =This%faces   (iel) 

            call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)
            call Kinematic                        (element, face, TL_, RES_1, .true.)
            call Theta_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_2, .true.)
            call Stresses                         (element, face, TL_, RES_3, .true., This%pressure )

            
    
            if (FlagNR == "NRP") Then
                call CalculateJacobianContributionsOf(Kinematic                        ,element, face, TL_, RES_1)
                call CalculateJacobianContributionsOf(Theta_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_2)
                call CalculateJacobianContributionsOf(Stresses                         ,element, face, TL_, RES_3, This%pressure)
                !Extra Unknowns
                call CalculateExtraJacobianContributionsOf(Stresses, element, face, TL_, RES_3, 1, This%pressure,   this%gid)
          end if
        end do 

        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyBoundaryConditions



    Function getCentroid(this) Result(output)
        Implicit None 
        Class(Bubble)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Centroid, Volume

        Volume   = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        Centroid = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Centroid = Centroid/Volume

        output = Centroid
    end Function getCentroid
    Function getDragForce(this) Result(output)
        Implicit none
        Class(Bubble)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: DragForce

        DragForce = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, DragForceCalculation)

        output = DragForce
    end Function getDragForce
    Function getAspectRatio(this)    Result(AR)
        Use GLOBAL_ARRAYS_MODULE,      Only: TL
        Implicit none
        Class(Bubble), Intent(In)            :: this
        Real(8)                              :: output

        Real(8), Dimension(:)  , Allocatable :: Z_coord, R_coord
        Real(8)                              :: Height, Width, AR
        
        call getBoundaryNodes(TL, This%elements, This%faces, "Z", Z_coord)
        call getBoundaryNodes(TL, This%elements, This%faces, "R", R_coord)
        
        Height = maxval(Z_coord) - minval(Z_coord)
        Width  = maxval(R_coord) ! - minval(R_coord) -> commented out since it is always zero

        AR = Height / (2.d0*Width)

        output = AR
    End Function getAspectRatio
    Function getInitialVolume(This)  Result(output)
        Implicit None 
        Class(Bubble)       :: This
        Real(8)             :: output

        output = This%InitialVolume 
    End Function getInitialVolume
    Function getVolume(This)  Result(output)
        Implicit None 
        Class(Bubble)       :: This
        Real(8)             :: output

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
    End Function getVolume



    Subroutine deconstructor(This) 
        Implicit None
        Type(Bubble) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )

    End Subroutine deconstructor

End Module BubbleBoundary

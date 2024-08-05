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
Module StaticCSBubbleBoundaryNewtonian
    Use Boundary_Equations
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use constrainJacobians
    Use DirichletBoundaries
    Private 

    Public :: StaticCSBubbleNewtonian, NewStaticCSBubbleNewtonian

    Type StaticCSBubbleNewtonian 
        ! StaticCSBubbleNewtonian Properties
        Real(8)                            :: Zcenter_o

        ! FEM : Prop
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces
        ! Global Contrain
        Integer                            :: gid
        Real(8)                            :: pressure, pressure_o, InitialPressure
        Real(8)                            ::           volume_o  , InitialVolume
        Real(8)                            :: mol  , Initialmol, Initialmol_dim
        contains 
            ! constrains and boundary conditions
            procedure :: applyBoundaryConditions
            procedure :: volumeConservation
            ! procedure :: fixCentroid
            procedure :: PressureVolumeConservation
            procedure :: check_engine
            ! setters
            procedure :: setProperties
            procedure :: setInitialPressure
            procedure :: setPressure
            procedure :: setInitialVolume
            procedure :: setCentroid_o
            procedure :: setVolume_o
            procedure :: setInitialmol_dim

            ! getters
            procedure :: getCentroid
            procedure :: getVelocity
            procedure :: getDragForce
            procedure :: getAspectRatio
            procedure :: getVolume
            procedure :: getdVtankdt

            final     :: deconstructor
    End Type StaticCSBubbleNewtonian



    Contains
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Constructor - Each time we remesh the object is deallocated and recreated
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function NewStaticCSBubbleNewtonian( elements, faces) Result(This)
        Implicit None 
        Integer, Dimension(:), Intent(In) :: elements
        Integer, Dimension(:), Intent(In) :: faces

        Type(StaticCSBubbleNewtonian)                      :: This

        
        if ( size(elements) /= size(faces) ) Then
            Print*, "[Error] NewStaticCSBubbleNewtonian."
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
    End Function NewStaticCSBubbleNewtonian

    
    

    
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                       volumeConservation                       
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function volumeConservation(this) Result(output)
        Implicit None
        Class(StaticCSBubbleNewtonian) :: this
        Real(8)       :: output

        Real(8)       :: Volume
        

        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        ! print '(A11, 2x, f6.4)', "initial_V =" , this%InitialVolume
        ! print '(A11, 2x, f6.4)', "current_V = ", Volume
        output = Volume - this%InitialVolume
        ! print*, "equation = ", output


        call loopOverElements(this%nelem, this%elements, this%faces, this%gid, SurfaceIntegration   ) ! first  constrain
    end Function volumeConservation


    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                   PressureVolumeConservation                   
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function PressureVolumeConservation(this) Result(output)
        Implicit None
        Class(StaticCSBubbleNewtonian) :: this
        Real(8)       :: output

        Real(8)       :: Volume
        

        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )

        print '(A11, 2x, f11.4)', "initialPV =" , this%InitialPressure * this%InitialVolume
        print '(A11, 2x, f11.4)', "currentPV = ", this%pressure * Volume 
        print '(A11, 2x, f11.4)', "initial_P =" , this%InitialPressure
        print '(A11, 2x, f11.4)', "current_P = ", this%pressure
        print '(A11, 2x, f6.4)', "initial_V =" , this%InitialVolume
        print '(A11, 2x, f6.4)', "current_V = ", Volume
        
        output = this%pressure * Volume - this%InitialPressure * this%InitialVolume
        print*, "equation = ", output
        pause

        call loopOverElements(this%nelem, this%elements, this%faces, this%gid, SurfaceIntegration, this%pressure ) ! first  constrain
    end Function PressureVolumeConservation

    subroutine check_engine(this)
        use check_for_floating_point_exceptions
        implicit none
        Class(StaticCSBubbleNewtonian) :: this

        call check_fp_exceptions(this%pressure, "pressure")
        call check_fp_exceptions(this%InitialPressure, "InitialPressure")
        call check_fp_exceptions(this%InitialVolume, "InitialVolume")
        
        print*, "pressure =", this%pressure
        print*, "InitialPressure =", this%InitialPressure
        print*, "InitialVolume =", this%InitialVolume

    end subroutine check_engine
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                    applyBoundaryConditions                     
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Subroutine applyBoundaryConditions(This, FlagNR, kinematic_logical)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        Implicit None 
        Class(StaticCSBubbleNewtonian)   , Intent(In)      :: This
        Character(len=3), Intent(In)      :: FlagNR 
        logical, Intent(In), optional     :: kinematic_logical 

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_2
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_3
        Integer                              :: iel 
        Integer                              :: element 
        Integer                              :: face
        Integer                              :: i
        Integer                              :: ii



        if (present(kinematic_logical) .and. (kinematic_logical)) then
            call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfResidual)
            If (FlagNR == "NRP") &
                call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfJacobian)
        endif

        do iel = 1, this%nelem
            element =This%elements(iel)
            face    =This%faces   (iel)

            call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)
            if (present(kinematic_logical) .and. (kinematic_logical)) then
                call Kinematic                        (element, face, TL_, RES_1, .true.)
            else
                call Theta_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_2, .true.)
                call Stresses                         (element, face, TL_, RES_3, .true., This%pressure )
            endif

            
    
            if (FlagNR == "NRP") Then
                if (present(kinematic_logical) .and. (kinematic_logical)) then
                    call CalculateJacobianContributionsOf(Kinematic                        ,element, face, TL_, RES_1)
                else
                    call CalculateJacobianContributionsOf(Theta_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_2)
                    call CalculateJacobianContributionsOf(Stresses                         ,element, face, TL_, RES_3, This%pressure, .true. )
                    !Extra Unknowns
                    call CalculateExtraJacobianContributionsOf(Stresses, element, face, TL_, RES_3, 1, This%pressure,   this%gid)
                endif
          end if
        end do 

        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyBoundaryConditions


    Subroutine setProperties(This, gid)
        Implicit None 
        Class(StaticCSBubbleNewtonian)           :: This 
        Integer, Intent(In)     :: gid

        this%gid        = gid
    End Subroutine setProperties
    Subroutine setInitialPressure(This, InitialPressure)
        Implicit None 
        Class(StaticCSBubbleNewtonian)       :: This
        Real(8), Intent(In) :: InitialPressure

        This%InitialPressure = InitialPressure
    End Subroutine setInitialPressure
    Subroutine setInitialVolume(This)
        Implicit None 
        Class(StaticCSBubbleNewtonian)       :: This

        This%InitialVolume = this%getVolume()
    End Subroutine setInitialVolume


    Subroutine setPressure(This, Pressure)
        Implicit None 
        Class(StaticCSBubbleNewtonian)       :: This
        Real(8), Intent(In) :: Pressure

        This%Pressure = Pressure
    End Subroutine setPressure

    Subroutine setInitialmol_dim(This)
        Implicit None 
        Class(StaticCSBubbleNewtonian)       :: This

        This%Initialmol_dim = (Pchar*this%InitialPressure) * (length_char**3*this%InitialVolume) / ( 8.314 * (273.d0 + 20.d0) )
    End Subroutine setInitialmol_dim
    Subroutine setInitialmol(This)
        Implicit None 
        Class(StaticCSBubbleNewtonian)       :: This
        ! dimensionless equation is P*V=n

        This%Initialmol = (this%InitialPressure) * (this%InitialVolume)
    End Subroutine setInitialmol
    
    Subroutine setmol(This, mol)
        use physical_module, only: Pchar, length_char
        Implicit None 
        Class(StaticCSBubbleNewtonian)       :: This
        Real(8), Intent(In) :: mol

    End Subroutine setmol


    Subroutine setCentroid_o(this)
        Implicit None 
        Class(StaticCSBubbleNewtonian)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Centroid, Volume
        
        Volume   = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        Centroid = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Centroid = Centroid/Volume

        this%Zcenter_o = Centroid
    end Subroutine setCentroid_o

    Subroutine setVolume_o(this)
        Implicit None 
        Class(StaticCSBubbleNewtonian)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Centroid, Volume
        
        This%volume_o   = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
    end Subroutine setVolume_o


    ! ********************************************************************


    Function getCentroid(this) Result(output)
        Implicit None 
        Class(StaticCSBubbleNewtonian)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Centroid, Volume

        Volume   = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        Centroid = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Centroid = Centroid/Volume

        output = Centroid
    end Function getCentroid

    Function getVelocity(this) Result(output)
        use TIME_INTEGRATION, only:dt
        Implicit None 
        Class(StaticCSBubbleNewtonian)              :: this
        Real(8)                                     :: output

        output = ( this%getCentroid() - this%Zcenter_o ) / dt

    end Function getVelocity

    Function getDragForce(this) Result(output)
        Implicit none
        Class(StaticCSBubbleNewtonian)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: DragForce

        DragForce = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, DragForceCalculationNewtonian)

        output = DragForce
    end Function getDragForce
    Function getAspectRatio(this)    Result(AR)
        Use GLOBAL_ARRAYS_MODULE,      Only: TL
        Implicit none
        Class(StaticCSBubbleNewtonian), Intent(In)            :: this
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
        Class(StaticCSBubbleNewtonian)       :: This
        Real(8)             :: output

        output = This%InitialVolume 
    End Function getInitialVolume
    Function getVolume(This)  Result(output)
        Implicit None 
        Class(StaticCSBubbleNewtonian)       :: This
        Real(8)             :: output

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
    End Function getVolume

    Function getdVtankdt(This)  Result(output)
        use time_integration, only: dt
        Implicit None 
        Class(StaticCSBubbleNewtonian)       :: This
        Real(8)             :: volume
        Real(8)             :: output

        volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )

        output = ( volume - this%Volume_o ) / dt
    End Function getdVtankdt


    Subroutine deconstructor(This) 
        Implicit None
        Type(StaticCSBubbleNewtonian) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )
    End Subroutine deconstructor

End Module StaticCSBubbleBoundaryNewtonian


Module BubbleDiffusionStaticCSBoundary
    Use Boundary_EquationsDO
    Use NumericalBoundaryJacobian
    Use ExtraEquations
    Use constrainJacobians
    Use DirichletBoundaries
    Private 

    Public :: BubbleDiffusionStaticCS, NewBubbleDiffusionStaticCS

    Type BubbleDiffusionStaticCS 
        ! BubbleDiffusionStaticCS Properties
        Real(8)                            :: Zcenter_o

        ! FEM : Prop
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes

        ! Global Contrain
        Integer                            :: gidP, gidC
        Real(8)                            :: pressure, pressure_o, InitialPressure
        Real(8)                            ::           volume_o  , InitialVolume
        Real(8)                            :: mol,      mol_o     , Initialmol
        ! Real(8)                            :: mol  , Initialmol, Initialmol_dim
        contains 
            ! constrains and boundary conditions
            procedure :: applyBoundaryConditions
            procedure :: volumeConservation
            ! procedure :: fixCentroid
            ! procedure :: PressureVolumeConservation
            procedure :: PressureVolumeMolConservation
            procedure :: molBalance
            
            procedure :: check_engine
            ! setters
            procedure :: setProperties
            procedure :: setInitialPressure
            procedure :: setPressure
            procedure :: setInitialVolume
            procedure :: setCentroid_o
            procedure :: setVolume_o
            procedure :: setInitialmol
            procedure :: setmol_o
            procedure :: setmol

            ! getters
            procedure :: getPressure
            procedure :: getCentroid
            procedure :: getVelocity
            procedure :: getDragForce
            procedure :: getAspectRatio
            procedure :: getVolume
            procedure :: getdVtankdt
            procedure :: getmol


            final     :: deconstructor
    End Type BubbleDiffusionStaticCS



    Contains
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Constructor - Each time we remesh the object is deallocated and recreated
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function NewBubbleDiffusionStaticCS( elements, faces) Result(This)
        Implicit None 
        Integer, Dimension(:), Intent(In) :: elements
        Integer, Dimension(:), Intent(In) :: faces

        Type(BubbleDiffusionStaticCS)                      :: This

        
        if ( size(elements) /= size(faces) ) Then
            Print*, "[Error] NewBubbleDiffusionStaticCS."
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

        call getBoundaryNodesOfWholeBoundary( This%nelem, This%elements, This%faces, This%nodes )

    End Function NewBubbleDiffusionStaticCS

    
    

    
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                       volumeConservation                       
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function volumeConservation(this) Result(output)
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: Volume
        
        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        
        output = Volume - this%InitialVolume
        
        call loopOverElements(this%nelem, this%elements, this%faces, this%gidP, SurfaceIntegration   ) ! first  constrain
    end Function volumeConservation


    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                   PressureVolumeConservation                   
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function PressureVolumeConservation(this) Result(output)
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: Volume
        

        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        
        output = this%pressure * Volume - this%InitialPressure * this%InitialVolume
     
        call loopOverElements(this%nelem, this%elements, this%faces, this%gidP, SurfaceIntegration, this%pressure ) 
        
    end Function PressureVolumeConservation

    Function PressureVolumeMolConservation(this) Result(output)
        Use physical_module, only: IdN
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: Volume

        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        
        output = this%pressure * Volume - IdN*this%mol
     
        call loopOverElements(this%nelem, this%elements, this%faces, this%gidP, SurfaceIntegration, this%pressure ) 
        
    end Function PressureVolumeMolConservation


    Function molBalance(this) Result(output)
        Use physical_module, only: PeN
        Use time_integration, only: dt
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: totalMolFlux

        totalMolFlux = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_F )
        
        ! In Deen p. 27, eq. 2.2-2, there is the form of the macroscopic balance
        output = (this%mol - this%mol_o)/dt + totalMolFlux
     
        call loopOverElements(this%nelem, this%elements, this%faces, this%gidC, int_n_dot_F ) 
        
    end Function molBalance

    subroutine check_engine(this)
        use check_for_floating_point_exceptions
        implicit none
        Class(BubbleDiffusionStaticCS) :: this

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
        Class(BubbleDiffusionStaticCS)   , Intent(In)      :: This
        Character(len=3), Intent(In)      :: FlagNR 
        logical, Intent(In), optional     :: kinematic_logical 

        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_1
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_2
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_3
        Integer                              :: iel, element, face
        Integer                              :: node_counter, node


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
                ! call Kinematic                        (element, face, TL_, RES_1, .true.)
                call Kinematic_mass                   (element, face, TL_, RES_1, .true.)
            else
                call Theta_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_2, .true.)
                call Stresses                         (element, face, TL_, RES_3, .true., This%pressure )
                call weakHenry                        (element, face, TL_, RES_1, .true., This%pressure )

            endif

            
    
            if (FlagNR == "NRP") Then
                if (present(kinematic_logical) .and. (kinematic_logical)) then
                    call CalculateJacobianContributionsOf(Kinematic_mass                   ,element, face, TL_, RES_1)
                else
                    call CalculateJacobianContributionsOf(Theta_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_2)
                    call CalculateJacobianContributionsOf(Stresses                         ,element, face, TL_, RES_3, This%pressure )
                    !Extra Unknowns
                    call CalculateExtraJacobianContributionsOf(Stresses                    ,element, face, TL_, RES_3, 1, This%pressure,   this%gidP)
                    
                    call CalculateJacobianContributionsOf(weakHenry                        ,element, face, TL_, RES_3, This%pressure )
                    ! Extra Unknowns
                    call CalculateExtraJacobianContributionsOf(weakHenry                   ,element, face, TL_, RES_3, 1, This%pressure,   this%gidP)
                endif
          end if
        end do 

        ! do node_counter = 1, size(this%nodes)
        !     node = this%nodes(node_counter)
        !     call ApplyDirichletAtNode_(node, "C", KoN*This%pressure, FlagNr, this%gidP )
        ! enddo

        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyBoundaryConditions


    Subroutine setProperties(This, gidP, gidC)
        Implicit None 
        Class(BubbleDiffusionStaticCS)           :: This 
        Integer, Intent(In)             :: gidP
        Integer, Intent(In), optional   ::  gidC

        this%gidP        = gidP
        if (present(gidC)) this%gidC = gidC
    End Subroutine setProperties

    Subroutine setInitialPressure(This, InitialPressure)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: InitialPressure

        This%InitialPressure = InitialPressure
    End Subroutine setInitialPressure

    Subroutine setInitialVolume(This)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This

        This%InitialVolume = this%getVolume()
        This%volume_o = This%InitialVolume
    End Subroutine setInitialVolume

    Subroutine setInitialmol(This)
        use physical_module, only:IdN
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        ! dimensionless equation is P*V-IdN*n=0
        ! CAUTION -- the n has 4pi/3 embeded in it (due to the volume)

        This%Initialmol = (this%InitialPressure) * (this%InitialVolume) / IdN
        This%mol_o      = This%Initialmol
        This%mol        = This%mol_o
    End Subroutine setInitialmol


    Subroutine setPressure(This, Pressure)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: Pressure

        This%Pressure = Pressure
    End Subroutine setPressure

    Subroutine setmol(This, mol)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: mol

        This%mol = mol
    End Subroutine setmol


    Subroutine setCentroid_o(this)
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: Centroid, Volume
        
        Volume   = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        Centroid = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Centroid = Centroid/Volume

        this%Zcenter_o = Centroid
    end Subroutine setCentroid_o

    Subroutine setVolume_o(this)
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        
        This%volume_o   = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
    end Subroutine setVolume_o

    Subroutine setmol_o(This, mol_o)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: mol_o

        This%mol_o = mol_o
    End Subroutine setmol_o


    ! ********************************************************************


    function getPressure(this)  Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8)             :: output

        output = This%pressure 

    end function getPressure

    Function getmol(This)  Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8)             :: output

        output = This%mol 
    End Function getmol

    Function getCentroid(this) Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Centroid, Volume

        Volume   = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        Centroid = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Centroid = Centroid/Volume

        output = Centroid
    end Function getCentroid

    Function getVelocity(this) Result(output)
        use TIME_INTEGRATION, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: zcenter
        Real(8)                                     :: output
        zcenter= this%getCentroid()
        output = ( zcenter - this%Zcenter_o ) / dt
    end Function getVelocity

    Function getDragForce(this) Result(output)
        Implicit none
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: DragForce

        DragForce = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, DragForceCalculation)

        output = DragForce
    end Function getDragForce
    Function getAspectRatio(this)    Result(AR)
        Use GLOBAL_ARRAYS_MODULE,      Only: TL
        Implicit none
        Class(BubbleDiffusionStaticCS), Intent(In)            :: this
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
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8)             :: output

        output = This%InitialVolume 
    End Function getInitialVolume
    Function getVolume(This)  Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8)             :: output

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
    End Function getVolume

    Function getdVtankdt(This)  Result(output)
        use time_integration, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8)             :: volume
        Real(8)             :: output

        volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )

        output = ( volume - this%Volume_o ) / dt
    End Function getdVtankdt


    Subroutine deconstructor(This) 
        Implicit None
        Type(BubbleDiffusionStaticCS) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )
    End Subroutine deconstructor

End Module BubbleDiffusionStaticCSBoundary

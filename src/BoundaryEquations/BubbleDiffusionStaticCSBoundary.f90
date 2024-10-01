
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
        Real(8)                            :: Zcenter, Zcenter_o, InitialZcenter

        ! FEM : Prop
        Integer                            :: nelem
        Integer, Dimension(:), Allocatable :: elements
        Integer, Dimension(:), Allocatable :: faces
        Integer, Dimension(:), Allocatable :: nodes

        ! Global Contrain
        Integer                            :: gidP, gidC, gidV, gidU
        Real(8)                            :: pressure, pressure_o, InitialPressure
        Real(8)                            :: mol,      mol_o     , Initialmol
        Real(8)                            :: volume,   volume_o  , InitialVolume
        Real(8)                            :: velocity, velocity_o, Initialvelocity
        contains 
            ! constrains and boundary conditions
            procedure :: PressureVolumeMolConservation
            procedure :: molBalance
            procedure :: volumeEquation
            procedure :: velocityEquation

            procedure :: applyBoundaryConditions
            
            ! ======== setters ========
            procedure :: setProperties
            ! ---- Initial values ----
            procedure :: setInitialPressure
            procedure :: setInitialmol
            procedure :: setInitialVolume
            procedure :: setInitialVelocity
            procedure :: setInitialZcenter
            ! ---- Previous values ----
            procedure :: setPressure_o
            procedure :: setmol_o
            procedure :: setVolume_o
            procedure :: setVelocity_o
            procedure :: setZcenter_o
            ! ---- Current values ----
            procedure :: setPressure
            procedure :: setmol
            procedure :: setVolume
            procedure :: setVelocity

            ! ======== getters ========

            ! ---- Current values ----
            procedure :: getPressure
            procedure :: getmol
            procedure :: getVolume
            procedure :: getVelocity
            procedure :: getZcenter
            procedure :: calculateZcenter

            procedure :: getZcenter_o
            
            procedure :: getdVtankdt


            procedure :: check_engine
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
    !                       Extra Equations                     
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Function PressureVolumeMolConservation(this) Result(output)
        Use physical_module, only: IdN
        use CSR_STORAGE, only: Ar_f
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: Volume

        output = this%pressure * this%Volume - IdN*this%mol

        Ar_f(this%gidP,:) = 0.d0

        write(*,*) ' '
        write(*,*) '------------------ PressureVolumeMolConservation ------------------'
        write(*,'(5(a14,4x))') 'this%pressure' , 'this%Volume' , 'this%mol', 'output'
        write(*,'(5(e14.7,4x))') this%pressure , this%Volume , this%mol, output
        write(*,*) ' '
        pause
    end Function PressureVolumeMolConservation



    Function molBalance(this) Result(output)
        Use time_integration, only: dt
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: totalMolFlux

        totalMolFlux = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_F )
        
        ! In Deen p. 27, eq. 2.2-2, there is the form of the macroscopic balance
        output = (this%mol - this%mol_o)/dt + totalMolFlux
     
        call loopOverElements(this%nelem, this%elements, this%faces, this%gidC, int_n_dot_F ) 
        
        write(*,*) ' '
        write(*,*) '------------------ molBalance ------------------'
        write(*,'(5(a14,4x))') 'this%mol', 'this%mol_o' , 'totalMolFlux', 'output'
        write(*,'(5(e14.7,4x))') this%mol, this%mol_o , totalMolFlux, output
        write(*,*) ' '
        pause
    end Function molBalance



    Function volumeEquation(this) Result(output)
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: calculatedVolume

        calculatedVolume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )

        output = - this%volume + calculatedVolume

        call loopOverElements(this%nelem, this%elements, this%faces, this%gidV, SurfaceIntegration )

        write(*,*) ' '
        write(*,*) '------------------ volumeEquation ------------------'
        write(*,'(5(a14,4x))') 'this%volume', 'calculatedVolume', 'output'
        write(*,'(5(e14.7,4x))') this%volume, calculatedVolume, output
        write(*,*) ' '
        pause

    end Function volumeEquation



    Function velocityEquation(this) Result(output)
        use time_integration, only: dt
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: calculatedInt_Z_dV

        calculatedInt_Z_dV = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV )

        output = - this%velocity*this%volume*dt - this%volume*this%Zcenter_o + calculatedInt_Z_dV

        call loopOverElements(this%nelem, this%elements, this%faces, this%gidU, int_Z_dV ) 

        write(*,*) ' '
        write(*,*) '------------------ velocityEquation ------------------'
        write(*,'(5(a14,4x))') 'this%velocity', 'calculatedInt_Z_dV', 'output'
        write(*,'(5(e14.7,4x))') this%velocity, calculatedInt_Z_dV, output
        write(*,*) ' '
        pause

    end Function velocityEquation


   
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
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_kinematic
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_thetaEquid
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_stresses
        Integer                              :: iel, element, face
        Integer                              :: node_counter, node
        Real(8)                              :: Volume


        if (present(kinematic_logical) .and. (kinematic_logical)) then
            call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfResidual)
            call updateAllNodesOfTheBoundary('R',This%elements, This%faces, ClearRowsOfResidual)
            If (FlagNR == "NRP") then
                call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfJacobian)
                call updateAllNodesOfTheBoundary('R',This%elements, This%faces, ClearRowsOfJacobian)
            endif
        endif

        
        do iel = 1, this%nelem
            element =This%elements(iel)
            face    =This%faces   (iel)

            call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)
            if (present(kinematic_logical) .and. (kinematic_logical)) then
                call Kinematic_mass_gasInterface        (element, face, TL_, RES_kinematic, .true., This%mol, this%volume, This%velocity )
                call Theta_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_thetaEquid, .true.)
            else
                call Stresses                         (element, face, TL_, RES_stresses, .true., This%pressure )
            endif

    
            if (FlagNR == "NRP") Then
                if (present(kinematic_logical) .and. (kinematic_logical)) then
                    call CalculateJacobianContributionsOf(Kinematic_mass_gasInterface        ,element, face, TL_, RES_kinematic, This%mol, this%Volume, This%velocity )
                    
                    call CalculateJacobianContributionsOf(Theta_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_thetaEquid)
                    !Extra Unknown
                    call CalculateExtraJacobianContributionsOf(Kinematic_mass_gasInterface   ,element, face, TL_, RES_kinematic, 1, This%mol, this%Volume, This%velocity, this%gidC)
                    call CalculateExtraJacobianContributionsOf(Kinematic_mass_gasInterface   ,element, face, TL_, RES_kinematic, 2, This%mol, this%Volume, This%velocity, this%gidV)
                    call CalculateExtraJacobianContributionsOf(Kinematic_mass_gasInterface   ,element, face, TL_, RES_kinematic, 3, This%mol, this%Volume, This%velocity, this%gidU)
                    
                else
                    call CalculateJacobianContributionsOf(Stresses                         ,element, face, TL_, RES_stresses, This%pressure )
                    !Extra Unknown
                    call CalculateExtraJacobianContributionsOf(Stresses                    ,element, face, TL_, RES_stresses, 1, This%pressure,   this%gidP)
                endif
          end if
        end do 

        do node_counter = 1, size(this%nodes)
            node = this%nodes(node_counter)
            call ApplyDirichletAtNode_(node, "C", KoN*This%pressure, FlagNr, this%gidP )
            ! call ApplyDirichletAtNode_(node, "C", 1.d0, FlagNr, this%gidP )
            ! call ApplyDirichletAtNode_(node, "C", KoN*This%pressure_o, FlagNr )
        enddo

        If ( Allocated(TL_) ) Deallocate(TL_)
    End Subroutine  applyBoundaryConditions



    ! ====================== setters ====================== 
    Subroutine setProperties(This, gidP, gidC, gidV, gidU)
        Implicit None 
        Class(BubbleDiffusionStaticCS)           :: This 
        Integer, Intent(In)             :: gidP
        Integer, Intent(In), optional   ::  gidC, gidV, gidU

        this%gidP        = gidP
        if (present(gidC)) then 
            this%gidC = gidC
            this%gidV = gidV
            this%gidU = gidU
        endif
    End Subroutine setProperties

    
    ! ------------------ Initial values ------------------

    Subroutine setInitialPressure(This, InitialPressure)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: InitialPressure

        This%InitialPressure = InitialPressure
        This%Pressure_o      = InitialPressure
        This%Pressure        = InitialPressure
    End Subroutine setInitialPressure


    Subroutine setInitialmol(This)
        use physical_module, only:IdN, pi4o3
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        ! dimensionless equation is P*V-IdN*n=0

        This%Initialmol = (this%InitialPressure) * (this%InitialVolume) / IdN
        This%mol_o      = This%Initialmol
        This%mol        = This%mol_o
    End Subroutine setInitialmol

    Subroutine setInitialVolume(This)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This

        This%InitialVolume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        This%volume_o      = This%InitialVolume
        This%volume        = This%InitialVolume
    End Subroutine setInitialVolume

    Subroutine setInitialVelocity(This)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This

        This%InitialVelocity = 0.d0
        This%velocity_o = This%InitialVelocity
        This%velocity   = This%InitialVelocity
    End Subroutine setInitialVelocity

    Subroutine setInitialZcenter(This)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This

        This%InitialZcenter = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        This%InitialZcenter = This%InitialZcenter/This%InitialVolume
        This%Zcenter_o      = This%InitialZcenter
    End Subroutine setInitialZcenter


    ! ------------------ Previous values ------------------
     
    Subroutine setPressure_o(This, Pressure_o)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: Pressure_o

        This%Pressure_o = Pressure_o
    End Subroutine setPressure_o

    Subroutine setmol_o(This, mol_o)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: mol_o

        This%mol_o = mol_o
    End Subroutine setmol_o

    Subroutine setVolume_o(this, Volume_o)
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8), Intent(In) :: Volume_o
        
        This%volume_o   = Volume_o
    end Subroutine setVolume_o

    Subroutine setVelocity_o(This, Velocity_o)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: Velocity_o

        This%Velocity_o = Velocity_o
    End Subroutine setVelocity_o

    Subroutine setZcenter_o(This)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This

        This%Zcenter_o = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        This%Zcenter_o = This%Zcenter_o/This%volume_o
    End Subroutine setZcenter_o


    ! ------------------ Current values ------------------

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

    Subroutine setVolume(This, Volume)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: Volume

        This%Volume = Volume
    End Subroutine setVolume

    Subroutine setVelocity(This, Velocity)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8), Intent(In) :: Velocity

        This%Velocity = Velocity
    End Subroutine setVelocity

   

   


    ! ********************************************************************
    ! ====================== getters ====================== 


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

    Function getVolume(This)  Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8)             :: output

        output = this%volume
    End Function getVolume

    Function getVelocity(this) Result(output)
        use TIME_INTEGRATION, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: output
        output = this%velocity
    end Function getVelocity

    Function getZcenter(this) Result(output)
        use TIME_INTEGRATION, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: output
        output = this%Zcenter
    end Function getZcenter

    Function getZcenter_o(this) Result(output)
        use TIME_INTEGRATION, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: output
        output = this%Zcenter_o
    end Function getZcenter_o

    
    

    Function calculateZcenter(this) Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        Zcenter = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Zcenter = Zcenter/This%Volume

        output = Zcenter
    end Function calculateZcenter


    Function getdVtankdt(This)  Result(output)
        use time_integration, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8)             :: volume
        Real(8)             :: output


        output = ( this%volume - this%Volume_o ) / dt
    End Function getdVtankdt


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

    Subroutine deconstructor(This) 
        Implicit None
        Type(BubbleDiffusionStaticCS) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )
    End Subroutine deconstructor

End Module BubbleDiffusionStaticCSBoundary
  
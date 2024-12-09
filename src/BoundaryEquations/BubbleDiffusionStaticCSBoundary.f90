
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
            procedure :: calculateZcenter
           
            procedure :: calculatendotF
            procedure :: calculatendotUbubblemUmesh_z
            procedure :: calculatendotUbubblemUmesh_r
            procedure :: calculatendotUmUmesh_z
            procedure :: calculatendotUmUmesh_r
            procedure :: calculatendotgradC_z
            procedure :: calculatendotgradC_r
            procedure :: calculateSherwood
            procedure :: calculateReynolds
            procedure :: calculate_kL
            
            procedure :: printEachContributionOfKinematicBC

            procedure :: getZcenter_o
            
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
        ! write(*,*) '------------------ PressureVolumeMolConservation ------------------'
        ! write(*,'(5(a14,4x))') 'this%pressure' , 'this%Volume' , 'this%mol', 'output'
        ! write(*,'(5(e14.7,4x))') this%pressure , this%Volume , this%mol, output
        ! write(*,*) ' '
        ! pause
    end Function PressureVolumeMolConservation



    Function molBalance(this) Result(output)
        Use time_integration, only: dt
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: totalMolFlux

        totalMolFlux = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_F )
        
        ! In Deen p. 27, eq. 2.2-2, there is the form of the macroscopic balance
        ! output = (this%mol - this%mol_o)/dt + totalMolFlux
        ! Now I multiply with dt inside int_n_dot_F to avoid dividing with a small term
        output = (this%mol - this%mol_o) + totalMolFlux
     
        call loopOverElements(this%nelem, this%elements, this%faces, this%gidC, int_n_dot_F ) 
        
        ! write(*,*) ' '
        ! write(*,*) '------------------ molBalance ------------------'
        ! write(*,'(5(a14,4x))') 'this%mol', 'this%mol_o' , 'totalMolFlux', 'output'
        ! write(*,'(5(e14.7,4x))') this%mol, this%mol_o , totalMolFlux, output
        ! write(*,*) ' '
        ! pause
    end Function molBalance



    Function volumeEquation(this) Result(output)
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: calculatedVolume

        calculatedVolume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )

        output = - this%volume + calculatedVolume

        call loopOverElements(this%nelem, this%elements, this%faces, this%gidV, SurfaceIntegration )

        ! write(*,*) ' '
        ! write(*,*) '------------------ volumeEquation ------------------'
        ! write(*,'(5(a14,4x))') 'this%volume', 'calculatedVolume', 'output'
        ! write(*,'(5(e14.7,4x))') this%volume, calculatedVolume, output
        ! write(*,*) ' '
        ! pause
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

        ! write(*,*) ' '
        ! write(*,*) '------------------ velocityEquation ------------------'
        ! write(*,'(5(a14,4x))') 'this%velocity', 'calculatedInt_Z_dV', 'output'
        ! write(*,'(5(e14.7,4x))') this%velocity, calculatedInt_Z_dV, output
        ! write(*,*) ' '
        ! pause
    end Function velocityEquation


   
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !                    applyBoundaryConditions                     
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Subroutine applyBoundaryConditions(This, FlagNR, naturalBCs, kinematicBC)
        Use GLOBAL_ARRAYS_MODULE,        Only: TL
        Use ENUMERATION_MODULE,          Only: NM_MESH
        Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
        use physical_module, only: Pinitial, Pchar, deltaPN
        Implicit None 
        Class(BubbleDiffusionStaticCS)   , Intent(In)      :: This
        Character(len=3), Intent(In)         :: FlagNR 
        logical,          Intent(In)         :: naturalBCs
        logical,         Intent(In), optional:: kinematicBC


        Real(8), Dimension(:,:), Allocatable :: TL_
        Real(8), Dimension(NBF_2d,NEQ_f)     :: RES_kinematic, RES_thetaEquid, RES_stresses, RES_concentration
        Integer                              :: iel, element, face
        Integer                              :: inode, node
        Real(8)                              :: Volume


        if (naturalBCs) then
            do iel = 1, this%nelem
                element =This%elements(iel)
                face    =This%faces   (iel)

                call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

                call Stresses            ( element, face, TL_, RES_stresses, .true., This%pressure )

                if (FlagNR == "NRP") Then
                    call CalculateJacobianContributionsOf  ( Stresses  ,element, face, TL_, RES_stresses, This%pressure )
                    
                    !Extra Unknown
                    call CalculateExtraJacobianContributionsOf  ( Stresses ,element, face, TL_, RES_stresses, 1, This%pressure,   this%gidP )
                endif
            enddo
        else
            if (present(kinematicBC) .and. kinematicBC) then
                call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfResidual)
                If (FlagNR == "NRP") then
                    call updateAllNodesOfTheBoundary('Z',This%elements, This%faces, ClearRowsOfJacobian)
                endif

                do iel = 1, this%nelem
                    element =This%elements(iel)
                    face    =This%faces   (iel)

                    call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

                    call Kinematic_mass        (element, face, TL_, RES_kinematic, .true. )

                    If (FlagNR == "NRP") then
                        call CalculateJacobianContributionsOf(Kinematic_mass     ,element, face, TL_, RES_kinematic )
                    endif
                enddo
            else
                do iel = 1, this%nelem
                    element =This%elements(iel)
                    face    =This%faces   (iel)

                    call copyArrayToLocalValues(TL, nm_mesh(element,:), 1, TL_)

                    call Theta_EQUIDISTRIBUTION_RESIDUAL_f(element, face, TL_, RES_thetaEquid, .true.)
                    
                    If (FlagNR == "NRP") then
                        call CalculateJacobianContributionsOf(Theta_EQUIDISTRIBUTION_RESIDUAL_f,element, face, TL_, RES_thetaEquid)
                    endif
                enddo
                
                do inode = 1, size(this%nodes)
                    node = this%nodes(inode)
                    call ApplyDirichletAtNode_(node, "C", ( this%Pressure - (Pinitial/Pchar) ) / deltaPN, FlagNr, this%gidP )
                enddo
            endif
        endif
        
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

    Function calculateZcenter(this) Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        Zcenter = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_Z_dV)
        Zcenter = Zcenter/This%Volume

        output = Zcenter
    end Function calculateZcenter

    Function getZcenter_o(this) Result(output)
        use TIME_INTEGRATION, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: output
        output = this%Zcenter_o
    end Function getZcenter_o


    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !                       Auxiliary functions
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    Function calculatendotF(this) Result(output)
        use time_integration, only:dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_F)
        output = output/dt
    end Function calculatendotF

    Function calculatendotUbubblemUmesh_z(this) Result(output)
        use time_integration, only:dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_UbubblemUmesh_z)
    end Function calculatendotUbubblemUmesh_z

    Function calculatendotUbubblemUmesh_r(this) Result(output)
        use time_integration, only:dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_UbubblemUmesh_r)
    end Function calculatendotUbubblemUmesh_r

    Function calculatendotUmUmesh_z(this) Result(output)
        use time_integration, only:dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_UmUmesh_z)
    end Function calculatendotUmUmesh_z

    Function calculatendotUmUmesh_r(this) Result(output)
        use time_integration, only:dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_UmUmesh_r)
    end Function calculatendotUmUmesh_r
   
    Function calculatendotgradC_z(this) Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_gradC_z)

    end Function calculatendotgradC_z

    Function calculatendotgradC_r(this) Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        output = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_gradC_r)

    end Function calculatendotgradC_r

    Function calculateSherwood(this) Result(output)
        use physical_module, only: Pinitial, Pchar, deltaPN
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: output

        Real(8)                                     :: dummy

        dummy = PeN * ( - integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_gradC_r) - integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_n_dot_gradC_z) )/ integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_dS)

        output = dummy / ( (this%Pressure - (Pinitial/Pchar))/deltaPN )
    end Function calculateSherwood

    Function calculate_kL(this) Result(output)
        use physical_module, only: Pinitial, Pchar, deltaPN, Cchar, length_char, time_char, Pchar, KHenry
        use time_integration, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: output

        Real(8)                                     :: dummy

        output = (Cchar*length_char**3/time_char * (this%mol - this%mol_o)/dt) / (length_char**2 * integrateOverAllElementsOfTheBoundary (this%elements, this%faces, int_dS)) / (KHenry*Pchar*this%Pressure  - KHenry*Pinitial)

    end Function calculate_kL

    Function calculateReynolds(this) Result(output)
        Implicit None 
        Class(BubbleDiffusionStaticCS)              :: this
        Real(8)                                     :: output

        Real(8)                                     :: Zcenter

        output = rho * (velocity_char*this%velocity) * (length_char*(this%volume/pi4o3)**0.333333d0) / viscosity

    end Function calculateReynolds

    Function printEachContributionOfKinematicBC(this) Result(output)
        use time_integration, only:dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        output = integrateOverAllElementsOfTheBoundaryAndPrintEachContribution (this%elements, this%faces, int_n_dot_UbubblemUmesh_z, 300)
        output = integrateOverAllElementsOfTheBoundaryAndPrintEachContribution (this%elements, this%faces, int_n_dot_UbubblemUmesh_r, 301)
        output = integrateOverAllElementsOfTheBoundaryAndPrintEachContribution (this%elements, this%faces, int_n_dot_UmUmesh_z, 302)
        output = integrateOverAllElementsOfTheBoundaryAndPrintEachContribution (this%elements, this%faces, int_n_dot_UmUmesh_r, 303)
        output = integrateOverAllElementsOfTheBoundaryAndPrintEachContribution (this%elements, this%faces, int_n_dot_gradC_z, 304)
        output = integrateOverAllElementsOfTheBoundaryAndPrintEachContribution (this%elements, this%faces, int_n_dot_gradC_r, 305)
    end Function printEachContributionOfKinematicBC

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

    Subroutine deconstructor(This) 
        Implicit None
        Type(BubbleDiffusionStaticCS) :: This

        If (Allocated(This%elements) ) Deallocate( This%elements)
        If (Allocated(This%faces)    ) Deallocate( This%faces   )
    End Subroutine deconstructor

End Module BubbleDiffusionStaticCSBoundary
  
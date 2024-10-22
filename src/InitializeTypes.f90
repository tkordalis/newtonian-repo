Module BoundaryConditions
    Use BOUNDARY_ENUMERATION_MODULE
    Use BubbleDiffusionStaticCSBoundary
    Use FixWallConcentrationBoundary
    Use SymmetryDiffusionBoundary 
    Use AmbientHenryBoundary

    Type(FixWallConcentration)         :: wall
    Type(SymmetryDiffusion)            :: symmetryaxis
    Type(AmbientHenry)                 :: ambientinterf
    Type(BubbleDiffusionStaticCS)      :: bubble, bubble2

    contains
    Subroutine DefineTheBoundaries()
        use TIME_INTEGRATION, only : increment
        Implicit None

        wall            = NewFixWallConcentration    (bnd1_elements, bnd1_faces)

        symmetryaxis    = NewSymmetryDiffusion       (bnd2_elements, bnd2_faces)
        call symmetryaxis%setPosition('X')

        bubble          = NewBubbleDiffusionStaticCS (bnd3_elements, bnd3_faces)
        call bubble%setProperties( gidP=1, gidC=2, gidV=3, gidU=4 )

        bubble2         = NewBubbleDiffusionStaticCS (bnd4_elements, bnd4_faces)
        call bubble2%setProperties( gidP=5, gidC=6, gidV=7, gidU=8 )

        ambientinterf   = NewAmbientHenry            (bnd5_elements, bnd5_faces)
        call ambientinterf%setDatumPressure(Pambient_o_Pchar)
        call ambientinterf%setConcentrationPressure(Pinitial/Pchar)

        if (increment .gt. 1) then
           
            call bubble%setPressure( Pressure_bubble )
            call bubble%setmol     ( mol_bubble      )
            Call bubble%setVolume  ( Volume_bubble   )
            Call bubble%setVelocity( Velocity_bubble )

            call bubble%setPressure_o( Pressure_bubbleo )
            call bubble%setmol_o     ( mol_bubbleo      )
            Call bubble%setVolume_o  ( Volume_bubbleo   )
            Call bubble%setVelocity_o( Velocity_bubbleo )
            Call bubble%setZcenter_o()

            call bubble2%setPressure( Pressure_bubble2 )
            call bubble2%setmol     ( mol_bubble2      )
            Call bubble2%setVolume  ( Volume_bubble2   )
            Call bubble2%setVelocity( Velocity_bubble2 )

            call bubble2%setPressure_o( Pressure_bubble2o )
            call bubble2%setmol_o     ( mol_bubble2o      )
            Call bubble2%setVolume_o  ( Volume_bubble2o   )
            Call bubble2%setVelocity_o( Velocity_bubble2o )
            Call bubble2%setZcenter_o()
        endif
       
    End Subroutine DefineTheBoundaries
End Module BoundaryConditions



Module InitialConditions
    Use BoundaryConditions
    Use Physical_Module
    Use ELEMENTS_MODULE,      only : Nex_f
    Use GLOBAL_ARRAYS_MODULE, only : TL, TLo, TLb, TLp
    Use RemeshVariables


    contains
    subroutine setInitalConditions
        Implicit None
        ! integer :: i,j
        ! I should print in the title of the .plts the Remesh_counter to read and define the boundaries correctly regardless
        Remesh_counter = 0
        TLo(:,:) = 0.D0
        TLo(:,getVariableId("Z"))   = Xm
        TLo(:,getVariableId("R"))   = Ym
        TLo(:,getVariableId("P"))   = Pambient_o_Pchar + ratio_of_pressures*( initial_position - TLo(:,getVariableId("Z")) )
        TLo(:,getVariableId("C"))   = 1.d0

        ! TLo(bubble%nodes(:),getVariableId("C"))   = KoN*Pressure_Bubble
        TLb = TLo
        TL  = TLo
        TLp = TL

        call bubble%setInitialVolume()
        call bubble%setInitialZcenter()

        call bubble2%setInitialVolume()
        call bubble2%setInitialZcenter()

        Pressure_Bubbleo = Pambient_o_Pchar + ratio_of_pressures*( initial_position - bubble%InitialZcenter ) + 2.d0/BoN
        Pressure_Bubble  = Pressure_Bubbleo

        Pressure_Bubble2o = Pambient_o_Pchar + ratio_of_pressures*( initial_position - bubble2%InitialZcenter ) + 2.d0/BoN
        Pressure_Bubble2  = Pressure_Bubble2o

        call bubble%setInitialPressure(Pressure_Bubble)
        call bubble%setInitialmol()
        call bubble%setInitialVelocity()

        call bubble2%setInitialPressure(Pressure_Bubble2)
        call bubble2%setInitialmol()
        call bubble2%setInitialVelocity()



        Mol_Bubbleo = bubble%getmol()
        Mol_Bubble  = Mol_Bubbleo

        volume_bubbleo = bubble%getVolume()
        volume_bubble  = volume_bubbleo

        velocity_bubbleo = bubble%getVelocity()
        velocity_bubble  = velocity_bubbleo
        ! -----------------------------------
        Mol_Bubble2o = bubble2%getmol()
        Mol_Bubble2  = Mol_Bubble2o

        volume_bubble2o = bubble2%getVolume()
        volume_bubble2  = volume_bubble2o

        velocity_bubble2o = bubble2%getVelocity()
        velocity_bubble2  = velocity_bubble2o

    end subroutine setInitalConditions
end Module InitialConditions



module solveAllExtraConstraints
    Use BoundaryConditions
    Use ELEMENTS_MODULE, only: NEX_f
    Use RemeshVariables

    Implicit None 
    
    contains
    Subroutine applyBCs_solveExtraConstraints( FlagNR, BubblePressure, BubbleMol, BubbleVolume, BubbleVelocity, Bubble2Pressure, Bubble2Mol, Bubble2Volume, Bubble2Velocity )
        use CSR_STORAGE, only: Ah_f
        use FLOW_ARRAYS_MODULE, only: Be_f
        use Physical_module, only: PeN, IdN
        use TIME_INTEGRATION, only:time, dt
        use pressure_variation
        implicit none
        character(*),                    intent(in)  :: FlagNR
        Real(8),                         intent(in)  :: BubblePressure, Bubble2Pressure
        Real(8),                         intent(in)  :: BubbleMol, Bubble2Mol
        Real(8),                         intent(in)  :: BubbleVolume, Bubble2Volume
        Real(8),                         intent(in)  :: BubbleVelocity, Bubble2Velocity

            Ah_f = 0.d0
            Call bubble%setPressure( BubblePressure )  ;  Call bubble2%setPressure( Bubble2Pressure )
            Call bubble%setmol     ( BubbleMol      )  ;  Call bubble2%setmol     ( Bubble2Mol      )
            Call bubble%setVolume  ( BubbleVolume   )  ;  Call bubble2%setVolume  ( Bubble2Volume   )
            Call bubble%setVelocity( BubbleVelocity )  ;  Call bubble2%setVelocity( Bubble2Velocity )
            

            call wall%applyBoundaryConditions(FlagNR, naturalBCs = .true.)

            call symmetryaxis%applyBoundaryConditions(FlagNR, naturalBCs = .true.)
            
            Call bubble%applyBoundaryConditions(FlagNR, naturalBCs = .true.)
            Call bubble2%applyBoundaryConditions(FlagNR, naturalBCs = .true.)

            call ambientinterf%applyBoundaryConditions( FlagNR, PressureChamber(time), naturalBCs = .true. )
            

            Be_f(1) = bubble%PressureVolumeMolConservation()
            
            Ah_f(1,1) = bubble%getVolume()
            Ah_f(1,2) = -IdN
            Ah_f(1,3) = bubble%getPressure()
            Ah_f(1,4) = 0.d0


            Be_f(2) = bubble%molBalance()

            Ah_f(2,1) = 0.d0
            ! since i multiply the whole equation with dt, the dt goes inside the integral and nb is multiplied with 1
            Ah_f(2,2) = 1.d0
            Ah_f(2,3) = 0.d0
            Ah_f(2,4) = 0.d0


            Be_f(3) = bubble%volumeEquation()

            Ah_f(3,1) = 0.d0
            Ah_f(3,2) = 0.d0
            Ah_f(3,3) = -1.d0
            Ah_f(3,4) = 0.d0

            Be_f(4) = bubble%velocityEquation()

            Ah_f(4,1) = 0.d0
            Ah_f(4,2) = 0.d0
            Ah_f(4,3) = -dt*bubble%getVelocity() - bubble%getZcenter_o()
            Ah_f(4,4) = -dt*bubble%getVolume()


            Be_f(5) = bubble2%PressureVolumeMolConservation()
            
            Ah_f(5,5) = bubble2%getVolume()
            Ah_f(5,6) = -IdN
            Ah_f(5,7) = bubble2%getPressure()
            Ah_f(5,8) = 0.d0


            Be_f(6) = bubble2%molBalance()

            Ah_f(6,5) = 0.d0
            ! since i multiply the whole equation with dt, the dt goes inside the integral and nb is multiplied with 1
            Ah_f(6,6) = 1.d0
            Ah_f(6,7) = 0.d0
            Ah_f(6,8) = 0.d0


            Be_f(7) = bubble2%volumeEquation()

            Ah_f(7,5) = 0.d0
            Ah_f(7,6) = 0.d0
            Ah_f(7,7) = -1.d0
            Ah_f(7,8) = 0.d0

            Be_f(8) = bubble2%velocityEquation()

            Ah_f(8,5) = 0.d0
            Ah_f(8,6) = 0.d0
            Ah_f(8,7) = -dt*bubble2%getVelocity() - bubble2%getZcenter_o()
            Ah_f(8,8) = -dt*bubble2%getVolume()


            Call bubble%applyBoundaryConditions(FlagNR, naturalBCs = .false.)
            Call bubble2%applyBoundaryConditions(FlagNR, naturalBCs = .false.)
            call symmetryaxis%applyBoundaryConditions(FlagNR, naturalBCs = .false.)
            Call bubble%applyBoundaryConditions(FlagNR, naturalBCs = .false., kinematicBC = .true.)
            Call bubble2%applyBoundaryConditions(FlagNR, naturalBCs = .false., kinematicBC = .true.)
            call wall%applyBoundaryConditions(FlagNR, naturalBCs = .false.)

            call ambientinterf%applyBoundaryConditions( FlagNR, PressureChamber(time), naturalBCs = .false. )

    End Subroutine applyBCs_solveExtraConstraints
    


end module solveAllExtraConstraints



Module BubbleOutput
    Use BoundaryConditions
    Use ELEMENTS_MODULE, only : Nex_f
    Use system_tools, only: check_dir
    contains
    Subroutine openBubbleFiles
        Implicit None
        character(*), parameter :: fileplace  = "./1_results_dat/"
        character(18), dimension(14) :: title_results0
        character(18), dimension(12) :: title_results1
        integer :: i
        
        call check_dir(fileplace)
        Open(20,File=fileplace//'results_dimensionless_b1.dat')
        Open(21,File=fileplace//'results_dimensional_b1.dat')
        Open(30,File=fileplace//'results_dimensionless_b2.dat')
        Open(31,File=fileplace//'results_dimensional_b2.dat')
        title_results0 = [ 'time', 'pressure', 'mol','volume','velocity', 'displacement', 'int_ndotF', 'int_ndotgradC', 'int_ndotUbmUmesh', 'int_ndotUmUmesh', 'ChamberP', 'concentration', 'Reynolds', 'Sherwood' ]
        title_results1 = [ 'time', 'pressure', 'mol','volume','velocity', 'displacement', 'int_ndotF', 'ChamberP', 'Radius', 'concentration', 'Reynolds', 'Sherwood' ]
        do i=1,size(title_results0)
            write(20,'(A25,3x)', advance='no') title_results0(i)
            write(30,'(A25,3x)', advance='no') title_results0(i)
        enddo
        write(20,*)
        write(30,*)
        do i=1,size(title_results1)
            write(21,'(A25,3x)', advance='no') title_results1(i)
            write(31,'(A25,3x)', advance='no') title_results1(i)
        enddo
        write(21,*)
        write(31,*)
    end Subroutine openBubbleFiles

    Subroutine WriteBubbleFiles(TIME)
        use pressure_variation, only: PressureChamber
        Use BoundaryConditions
        Use Formats
        Implicit none
        Real(8), Intent(In) :: Time
        Real(8) :: dummy

        write(20,'(14(f26.16,3x))') TIME, bubble%getPressure(), bubble%getmol(), bubble%getVolume(), bubble%getVelocity(), &
                                    bubble%calculateZcenter(), bubble%calculatendotF(), bubble%calculatendotgradC_z()+bubble%calculatendotgradC_r(), &
                                    bubble%calculatendotUbubblemUmesh_z() + bubble%calculatendotUbubblemUmesh_r(), bubble%calculatendotUmUmesh_z() + bubble%calculatendotUmUmesh_r(), &
                                    PressureChamber(time), bubble%getmol()/bubble%getVolume(), bubble%calculateReynolds(), bubble%calculateSherwood()

        write(21,'(12(f26.16,3x))') TIME*time_char, Pchar*bubble%getPressure(), nchar*bubble%getmol(), length_char**3.d0*bubble%getVolume(), velocity_char*bubble%getVelocity(), &
                                    length_char*bubble%calculateZcenter(), nchar/time_char*bubble%calculatendotF(), Pchar*PressureChamber(time), length_char*(3.d0*bubble%getVolume()/4.d0/pi)**0.333333d0, &
                                    Cchar*bubble%getmol()/bubble%getVolume(), bubble%calculateReynolds(), bubble%calculateSherwood()


        write(30,'(14(f26.16,3x))') TIME, bubble2%getPressure(), bubble2%getmol(), bubble2%getVolume(), bubble2%getVelocity(), &
                                    bubble2%calculateZcenter(), bubble2%calculatendotF(), bubble2%calculatendotgradC_z()+bubble2%calculatendotgradC_r(), &
                                    bubble2%calculatendotUbubblemUmesh_z() + bubble2%calculatendotUbubblemUmesh_r(), bubble2%calculatendotUmUmesh_z() + bubble2%calculatendotUmUmesh_r(), &
                                    PressureChamber(time), bubble2%getmol()/bubble2%getVolume(), bubble2%calculateReynolds(), bubble2%calculateSherwood()

        write(31,'(12(f26.16,3x))') TIME*time_char, Pchar*bubble2%getPressure(), nchar*bubble2%getmol(), length_char**3.d0*bubble2%getVolume(), velocity_char*bubble2%getVelocity(), &
                                    length_char*bubble2%calculateZcenter(), nchar/time_char*bubble2%calculatendotF(), Pchar*PressureChamber(time), length_char*(3.d0*bubble2%getVolume()/4.d0/pi)**0.333333d0, &
                                    Cchar*bubble2%getmol()/bubble2%getVolume(), bubble2%calculateReynolds(), bubble2%calculateSherwood()

! dummy = bubble%printEachContributionOfKinematicBC()
    End Subroutine WriteBubbleFiles
end Module BubbleOutput
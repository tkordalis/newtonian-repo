Module BoundaryConditions
    Use BOUNDARY_ENUMERATION_MODULE
    Use BubbleDiffusionStaticCSBoundary
    Use FixWallConcentrationBoundary
    Use SymmetryDiffusionBoundary 
    Use AmbientHenryBoundary

    Type(FixWallConcentration)         :: wall
    Type(SymmetryDiffusion)            :: symmetryaxis
    Type(AmbientHenry)                 :: ambientinterf
    Type(BubbleDiffusionStaticCS)      :: bubble

    contains
    Subroutine DefineTheBoundaries()
        use TIME_INTEGRATION, only : increment
        Implicit None

        wall            = NewFixWallConcentration    (bnd1_elements, bnd1_faces)

        symmetryaxis    = NewSymmetryDiffusion       (bnd2_elements, bnd2_faces)
        call symmetryaxis%setPosition('X')

        bubble          = NewBubbleDiffusionStaticCS (bnd3_elements, bnd3_faces)
        call bubble%setProperties( gidP=1, gidC=2, gidV=3, gidU=4 )

        ambientinterf   = NewAmbientHenry            (bnd4_elements, bnd4_faces)
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

        Pressure_Bubbleo = Pambient_o_Pchar + ratio_of_pressures*( initial_position ) + 2.d0/BoN
        Pressure_Bubble  = Pressure_Bubbleo

        ! TLo(bubble%nodes(:),getVariableId("C"))   = KoN*Pressure_Bubble
        TLb = TLo
        TL  = TLo
        TLp = TL

        
        call bubble%setInitialPressure(Pressure_Bubble)
        call bubble%setInitialvolume()
        call bubble%setInitialmol()
        call bubble%setInitialVelocity()


        Mol_Bubbleo = bubble%getmol()
        Mol_Bubble  = Mol_Bubbleo

        volume_bubbleo = bubble%getVolume()
        volume_bubble  = volume_bubbleo

        velocity_bubbleo = bubble%getVelocity()
        velocity_bubble  = velocity_bubbleo

    end subroutine setInitalConditions
end Module InitialConditions



module solveAllExtraConstraints
    Use BoundaryConditions
    Use ELEMENTS_MODULE, only: NEX_f
    Use RemeshVariables

    Implicit None 
    
    contains
    Subroutine applyBCs_solveExtraConstraints( FlagNR, BubblePressure, BubbleMol, BubbleVolume, BubbleVelocity )
        use CSR_STORAGE, only: Ah_f
        use FLOW_ARRAYS_MODULE, only: Be_f
        use Physical_module, only: PeN, IdN
        use TIME_INTEGRATION, only:time, dt
        use pressure_variation
        implicit none
        character(*),                    intent(in)  :: FlagNR
        Real(8),                         intent(in)  :: BubblePressure
        Real(8),                         intent(in)  :: BubbleMol
        Real(8),                         intent(in)  :: BubbleVolume
        Real(8),                         intent(in)  :: BubbleVelocity

            Ah_f = 0.d0
            Call bubble%setPressure( BubblePressure )
            Call bubble%setmol     ( BubbleMol      )
            Call bubble%setVolume  ( BubbleVolume   )
            Call bubble%setVelocity( BubbleVelocity )
            

            call wall%applyBoundaryConditions(FlagNR, naturalBCs = .true.)

            call symmetryaxis%applyBoundaryConditions(FlagNR, naturalBCs = .true.)
            
            Call bubble%applyBoundaryConditions(FlagNR, naturalBCs = .true.)

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


            Call bubble%applyBoundaryConditions(FlagNR, naturalBCs = .false.)
            call symmetryaxis%applyBoundaryConditions(FlagNR, naturalBCs = .false.)
            Call bubble%applyBoundaryConditions(FlagNR, naturalBCs = .false., kinematicBC = .true.)
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
        character(18), dimension(14) :: title_results20
        character(18), dimension(12) :: title_results21
        integer :: i
        
        call check_dir(fileplace)
        Open(20,File=fileplace//'results_dimensionless.dat')
        Open(21,File=fileplace//'results_dimensional.dat')
        title_results20 = [ 'time', 'pressure', 'mol','volume','velocity', 'displacement', 'int_ndotF', 'int_ndotgradC', 'int_ndotUbmUmesh', 'int_ndotUmUmesh', 'ChamberP', 'concentration', 'Reynolds', 'Sherwood' ]
        title_results21 = [ 'time', 'pressure', 'mol','volume','velocity', 'displacement', 'int_ndotF', 'ChamberP', 'Radius', 'concentration', 'Reynolds', 'Sherwood' ]
        do i=1,size(title_results20)
            write(20,'(A25,3x)', advance='no') title_results20(i)
        enddo
        write(20,*)
        do i=1,size(title_results21)
            write(21,'(A25,3x)', advance='no') title_results21(i)
        enddo
        write(21,*)
    end Subroutine openBubbleFiles

    Subroutine WriteBubbleFiles(TIME)
        use pressure_variation, only: PressureChamber
        Use Physical_module, only: Pressure_Bubble
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

! dummy = bubble%printEachContributionOfKinematicBC()
    End Subroutine WriteBubbleFiles
end Module BubbleOutput
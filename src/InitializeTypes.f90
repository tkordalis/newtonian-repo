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
        Implicit None

        wall            = NewFixWallConcentration    (bnd1_elements, bnd1_faces)

        symmetryaxis    = NewSymmetryDiffusion       (bnd2_elements, bnd2_faces)
        call symmetryaxis%setPosition('X')

        bubble          = NewBubbleDiffusionStaticCS (bnd3_elements, bnd3_faces)
        call bubble%setProperties( gidP=1, gidC=1 )

        ambientinterf   = NewAmbientHenry            (bnd4_elements, bnd4_faces)
        call ambientinterf%setDatumPressure(Pambient_o_Pchar)
       
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
    TLb = TLo
    TL  = TLo
    TLp = TL

    Pressure_Bubbleo = Pambient_o_Pchar + ratio_of_pressures*( initial_position ) + 2.d0/BoN
    Pressure_Bubble  = Pressure_Bubbleo
    
    call bubble%setInitialPressure(Pressure_Bubble)
    call bubble%setInitialvolume()
    call bubble%setInitialmol()

    Mol_Bubbleo = bubble%getmol()
    Mol_Bubble   = Mol_Bubbleo

    

    print*, ' '
    ! print*, 'initialMol = ', bubble%Initialmol
    ! print*, 'initialMol = ', bubble%Initialmol/bubble%InitialVolume 
    print*, 'KoN * Pambient = ', PeN * KoN

    print*, ' '
    ! write(*,"(10X,A15,2X,F16.8)") "velocity_char ="    , velocity_char
    ! write(*,"(10X,A15,2X,F16.8)") "time_char ="    , time_char
    ! write(*,"(10X,A15,2X,E12.5)") "P ="    , bubble%InitialPressure* Pchar
    ! write(*,"(10X,A15,2X,E12.5)") "V ="    , bubble%InitialVolume * length_char**3.d0
    ! write(*,"(10X,A15,2X,E12.5)") "PV ="    , bubble%InitialPressure* Pchar * bubble%InitialVolume * length_char**3.d0
    ! write(*,"(10X,A15,2X,E12.5)") "RT ="    , Rgas*Tgas
    ! write(*,"(10X,A15,2X,E12.5)") "n ="    , bubble%InitialPressure* Pchar * bubble%InitialVolume * length_char**3.d0/(Rgas*Tgas)
    ! write(*,"(10X,A15,2X,E12.5)") "n_initMol ="    , bubble%Initialmol*cchar*length_char**3.d0
    ! write(*,"(10X,A15,2X,E12.5)") "n_initMol2 ="    , bubble%InitialPressure* bubble%InitialVolume/IdN*cchar*length_char**3.d0
    ! print*, 'initialMol = ', Mol_Bubbleo
    ! print*, 'cchar = ', cchar 
    ! print*, 'Pressure_Bubble = ', Pressure_Bubble*Pchar
    ! print*, 'NondimPressure_Bubble = ', Pressure_Bubble
    ! print*, 'Volume_Bubble = ', bubble%InitialVolume*length_char**3
    ! print*, 'NondimVolume_Bubble = ', bubble%InitialVolume
    ! print*, 'NondimPV = ', Pressure_Bubble*bubble%InitialVolume
    ! print*, 'Mol_Bubble = ', Pressure_Bubble*Pchar*bubble%InitialVolume*length_char**3/Rgas/Tgas
    ! print*, 'Mol_Bubble/4pi/3 = ', Pressure_Bubble*Pchar*bubble%InitialVolume*length_char**3/Rgas/Tgas/(4.d0*pi/3.d0)
    ! print*, '2--Mol_Bubble = ', Cchar*length_char**3
    ! print*, 'NondimMol_Bubble = ', Pressure_Bubble*bubble%InitialVolume/IdN
    ! print*, 'NondimMol_Bubble/4pi/3 = ', Pressure_Bubble*bubble%InitialVolume/IdN/4.1889d0
    ! pause
    
    ! do i=1, size(tlo,1)
    !     write(404,'(f16.8,3x)') (tlo(i,j),j=1,size(tlo,2))
    ! enddo
   
    ! call check_fp_exceptions()
    
    end subroutine setInitalConditions
end Module InitialConditions



module solveAllExtraConstraints
    Use BoundaryConditions
    Use ELEMENTS_MODULE, only: NEX_f
    Use RemeshVariables

    Implicit None 
    
    contains
    Subroutine applyBCs_solveExtraConstraints( FlagNR, BubblePressure, BubbleMol )
        use CSR_STORAGE, only: Ah_f
        use FLOW_ARRAYS_MODULE, only: Be_f
        use Physical_module, only: PeN, IdN
        use TIME_INTEGRATION, only:time, dt
        use pressure_variation
        implicit none
        character(*),                    intent(in)  :: FlagNR
        Real(8),                         intent(in)  :: BubblePressure
        Real(8),                         intent(in)  :: BubbleMol
        Real(8)                                       :: dVtankdt

            Ah_f = 0.d0
            Call bubble%setPressure( BubblePressure )
            Call bubble%setmol( BubbleMol )
            Call bubble%applyBoundaryConditions(FlagNR)
            Be_f(1) = bubble%PressureVolumeMolConservation()
            Ah_f(1,1) = bubble%getVolume()
            Ah_f(1,2) = -IdN

            Be_f(2) = bubble%molBalance()
            Ah_f(2,2) = 1.d0/dt

            ! Be_f(1) = bubble%volumeConservation()
            ! Ah_f(:,:) = 0.d0

            dVtankdt = bubble%getdVtankdt()
            


            call symmetryaxis%applyBoundaryConditions(FlagNR)
            Call bubble%applyBoundaryConditions(FlagNR, .true.)
            call wall%applyBoundaryConditions(FlagNR)

            call ambientinterf%applyBoundaryConditions( FlagNR, dVtankdt, PressureChamber(time) )

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
        character(18), dimension(7) :: title_results
        integer :: i
        
        call check_dir(fileplace)
        Open(20,File=fileplace//'results_dimensionless.dat')
        title_results = [ 'time', 'displacement', 'velocity', 'pressure', 'volume', 'ChamberP', 'mol' ]
        do i=1,size(title_results)
            write(20,'(A17,3x)', advance='no') title_results(i)
        enddo
        write(20,*) " " 
    end Subroutine openBubbleFiles

    Subroutine WriteBubbleFiles(TIME)
        use pressure_variation, only: PressureChamber
        Use Physical_module, only: Pressure_Bubble
        Use BoundaryConditions
        Use Formats
        Implicit none
        Real(8), Intent(In) :: Time

        write(20,'(7(f16.7,3x))') TIME, bubble%getCentroid(), bubble%getVelocity(), Pressure_Bubble, bubble%getVolume(), PressureChamber(time), bubble%getmol()
        
    End Subroutine WriteBubbleFiles
end Module BubbleOutput
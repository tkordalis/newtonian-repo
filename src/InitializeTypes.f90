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
        call bubble%setProperties( gidP=1 )

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
    ! TLo(:,getVariableId("C"))   = 1.d0
    TLb = TLo
    TL  = TLo
    TLp = TL

    Pressure_Bubbleo = Pambient_o_Pchar + ratio_of_pressures*( initial_position ) + 2.d0/BoN
    Pressure_Bubble  = Pressure_Bubbleo
    
    call bubble%setInitialPressure(Pressure_Bubble)
    call bubble%setInitialvolume()
    
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
    Subroutine applyBCs_solveExtraConstraints( FlagNR, Bubble1Pressure )
        use CSR_STORAGE, only: Ah_f
        use FLOW_ARRAYS_MODULE, only: Be_f
        use Physical_module, only: Pambient_o_Pchar, Rtank, pi
        use TIME_INTEGRATION, only:time
        use pressure_variation
        implicit none
        character(*),                    intent(in)  :: FlagNR
        Real(8),                         intent(in)  :: Bubble1Pressure
        Real(8)                                       :: dVtankdt


            Call bubble%setPressure( Bubble1Pressure )
            Call bubble%applyBoundaryConditions(FlagNR)
            ! Be_f(1) = bubble%PressureVolumeConservation()
            ! Ah_f(:,:) = bubble%getVolume()
            Be_f(1) = bubble%volumeConservation()
            Ah_f(:,:) = 0.d0

            dVtankdt = bubble%getdVtankdt()
            ! write(404,'(4(f16.9,2x))') time, dVtankdt, Rtank, dVtankdt / (pi*Rtank**2)
            ! call bubble%check_engine()


            call symmetryaxis%applyBoundaryConditions(FlagNR)
            Call bubble%applyBoundaryConditions(FlagNR, .true.)
            call wall%applyBoundaryConditions(FlagNR)

            call ambientinterf%applyBoundaryConditions( FlagNR, dVtankdt, PressureChamber(time) )
            ! call ambientinterf%applyBoundaryConditions( FlagNR, dVtankdt, Pambient_o_Pchar )

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
        character(18), dimension(6) :: title_results
        integer :: i
        
        call check_dir(fileplace)
        Open(20,File=fileplace//'results_dimensionless.dat')
        title_results = [ 'time', 'displacement', 'velocity', 'pressure', 'volume', 'ChamberP' ]
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

        write(20,'(6(f16.7,3x))') TIME, bubble%getCentroid(), bubble%getVelocity(), Pressure_Bubble, bubble%getVolume(), PressureChamber(time)
        
    End Subroutine WriteBubbleFiles
end Module BubbleOutput
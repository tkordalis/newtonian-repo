Module BoundaryConditions
  Use FixWallBoundary
  Use SymmetryBoundary 
  Use InflatedBubbleBoundary 
  Use AmbientBoundary
  Use InternalEquidistribution
  Use BOUNDARY_ENUMERATION_MODULE
  Use RemeshVariables

  Type(FixWall)                 :: Wall
  Type(Symmetry)                :: SymmetryAxis
  Type(InflatedBubble)          :: Bubble1
  Type(Ambient)                 :: topSurface
  Type(IntEquidistribution)     :: InternalEquidistributionZ
  Type(IntEquidistribution)     :: InternalEquidistributionR

  contains
    Subroutine DefineTheBoundaries()
      Implicit None

        Wall           = NewFixWall ( bnd1_elements, bnd1_faces )
        
        SymmetryAxis   = NewSymmetry( bnd2_elements, bnd2_faces )
        call SymmetryAxis%setPosition('X')
        
        Bubble1        = NewInflatedBubble ( bnd3_elements, bnd3_faces)
        call Bubble1%setProperties ( gid = 1 )
        
        topSurface     = NewAmbient ( bnd4_elements, bnd4_faces )
        call topSurface%setDatumPressure( Pambient_o_Pchar )
        

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
    ! I should print in the title of the .plts the Remesh_counter to read and define the boundaries correctly regardless
    Remesh_counter = 0
    TLo(:,:) = 0.D0
    TLo(:,getVariableId("Z"))   = Xm
    TLo(:,getVariableId("R"))   = Ym
    TLo(:,getVariableId("P"))   = Pambient_o_Pchar + ratio_of_pressures*( -initial_position + TLo(:,getVariableId("Z")) )
    TLb = TLo
    ! ----------------------------------------------------------------------
    ! ASSIGN NR INITIAL GUESS
    ! ----------------------------------------------------------------------
    TL  = TLo
    TLp = TL
    position         = initial_position
    ! Pressure_Bubbleo = 0.d0
    Pressure_Bubbleo = Pambient_o_Pchar
    Pressure_Bubble  = Pressure_Bubbleo
    
    call Bubble1%setPressure(Pressure_Bubble)
    
    end subroutine setInitalConditions
end Module InitialConditions



module solveAllExtraConstraints
    Use BoundaryConditions
    Use ELEMENTS_MODULE, only: NEX_f
    Use RemeshVariables

    Implicit None
    
    contains
    Subroutine applyBCs_and_solveExtraConstraints( FlagNR, JacGlobalEq_GlobUnkn, ExtraRes, Bubble1Pressure )
        use TIME_INTEGRATION, only: time
        implicit none
        character(*),                    intent(in)  :: FlagNR
        Real(8), dimension(NEX_f,NEX_f), intent(out) :: JacGlobalEq_GlobUnkn
        Real(8), dimension(NEX_f),       intent(out) :: ExtraRes
        Real(8),                         intent(in)  :: Bubble1Pressure



            if ( Remesh_counter .eq. 0 ) then
                call InternalEquidistributionZ%applyEquidistribution(FlagNR)
                call InternalEquidistributionR%applyEquidistribution(FlagNR)
            endif
            
            ! if (Remesh_counter_structuredInTheFront .eq. 0) then
            !     call InternalEquidistributionR%applyEquidistribution(FlagNR)
            ! endif


            Call Bubble1%setPressure( Bubble1Pressure )
            Call Bubble1%applyBoundaryConditions(FlagNR)
            ExtraRes(1) = Bubble1%airFlowrateEquation(time)
            JacGlobalEq_GlobUnkn(:,:) = 0.d0

            call SymmetryAxis%applyBoundaryConditions(FlagNR)
            
            call Wall%applyBoundaryConditions(FlagNR)
            
            call topSurface%applyBoundaryConditions(FlagNR)

    End Subroutine applyBCs_and_solveExtraConstraints

    


end module solveAllExtraConstraints



Module BubbleOutput
    Use BoundaryConditions
    Use ELEMENTS_MODULE, only : Nex_f
    Use system_tools, only: check_dir
    contains
    Subroutine openBubbleFiles
        Implicit None
        character(*), parameter :: fileplace  = "./1_results_dat/"
        
            call check_dir(fileplace)
            Open(20,File=fileplace//'Flowrate.dat')
            Open(21,File=fileplace//'Pressure_BubbleSB.dat')
            Open(22,File=fileplace//'MaxZ.dat')
            Open(23,File=fileplace//'MaxR.dat')
           
        
    end Subroutine openBubbleFiles            
    Subroutine calculateBubbleVariables(DT)
        Use Physical_module  
        Implicit none
        real(8), intent(in) :: DT
        Real(8)             :: displacement1_dt, displacement2_dt
        
       
    end Subroutine calculateBubbleVariables        
    Subroutine WriteBubbleFiles(TIME)
        Use Physical_module!, only: Pressure_Bubble, calculateFlowrate
        Use BoundaryConditions
        Use Formats
        Implicit none
        Real(8), Intent(In) :: Time
        Real(8)             :: flowrate           
       
        flowrate = calculateFlowrate(time)

        write(20,'(4(f20.15,3x))') TIME, flowrate,        time_char*time, volumetric_flowrate*flowrate
        write(21,'(4(f18.11,3x))') TIME, Pressure_Bubble, time_char*time, Pchar*Pressure_Bubble
        write(22,'(4(f18.11,3x))') TIME, Bubble1%getMaxZcoordinate(), time_char*time, length_char*Bubble1%getMaxZcoordinate()
        write(23,'(4(f18.11,3x))') TIME, Bubble1%getMaxRcoordinate(), time_char*time, length_char*Bubble1%getMaxRcoordinate()
       
    End Subroutine WriteBubbleFiles
end Module BubbleOutput
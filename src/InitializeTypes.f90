Module BoundaryConditions
  Use BOUNDARY_ENUMERATION_MODULE
  Use FixConcentrationValueBoundary
  Use FixFluxOfConcentrationBoundary

  Type(FixConcentrationValue)        :: Left
  Type(FixConcentrationValue)        :: Right
  ! Type(FixConcentrationValue)        :: Up
  ! Type(FixConcentrationValue)        :: Down
  ! Type(FixFluxConcentrationValue)    :: Left
  ! Type(FixFluxConcentrationValue)    :: Right
  Type(FixFluxConcentrationValue)    :: Up
  Type(FixFluxConcentrationValue)    :: Down
  
  contains
    Subroutine DefineTheBoundaries()
      Implicit None

        Down      = NewFixFluxConcentrationValue  (bnd1_elements, bnd1_faces)
        call Down%setFluxValueAtTheBoundary(0.d0)

        ! Down    = NewFixConcentrationValue  (bnd1_elements, bnd1_faces)
        ! call Down%setValueAtTheBoundary(1.d0)

        Up    =  NewFixFluxConcentrationValue (bnd2_elements, bnd2_faces)
        call Up%setFluxValueAtTheBoundary(0.d0)


        ! Up   = NewFixConcentrationValue  (bnd2_elements, bnd2_faces)
        ! call Up%setValueAtTheBoundary(1.d0)

        Left    = NewFixConcentrationValue  (bnd3_elements, bnd3_faces)
        call Left%setValueAtTheBoundary(1.d0)

        ! Left    = NewFixFluxConcentrationValue  (bnd3_elements, bnd3_faces)
        ! call Left%setFluxValueAtTheBoundary(-1.d0)

        Right   = NewFixConcentrationValue  (bnd4_elements, bnd4_faces)
        call Right%setValueAtTheBoundary(0.d0)

        ! Right   = NewFixFluxConcentrationValue  (bnd4_elements, bnd4_faces)
        ! call Right%setFluxValueAtTheBoundary(0.d0)

    
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
    ! TLo(:,getVariableId("Z"))   = Xm
    ! TLo(:,getVariableId("R"))   = Ym
    TLb = TLo
    ! ----------------------------------------------------------------------
    ! ASSIGN NR INITIAL GUESS
    ! ----------------------------------------------------------------------
    TL  = TLo
    TLp = TL
    
    end subroutine setInitalConditions
end Module InitialConditions



module solveAllExtraConstraints
    Use BoundaryConditions
    Use ELEMENTS_MODULE, only: NEX_f
    Use RemeshVariables

    Implicit None
    
    contains
    Subroutine applyBCs_and_solveExtraConstraints( FlagNR, JacGlobalEq_GlobUnkn, ExtraRes )
        use TIME_INTEGRATION, only: time
        implicit none
        character(*),                    intent(in)  :: FlagNR
        Real(8), dimension(NEX_f,NEX_f), intent(out), optional :: JacGlobalEq_GlobUnkn
        Real(8), dimension(NEX_f),       intent(out), optional :: ExtraRes

        call Up%applyBoundaryConditions(FlagNR)
        call Down%applyBoundaryConditions(FlagNR)
        call Left%applyBoundaryConditions(FlagNR)
        call Right%applyBoundaryConditions(FlagNR)
            
        
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
            
           
        
    end Subroutine openBubbleFiles
    Subroutine calculateVariables(DT)
        Use Physical_module  
        Implicit none
        real(8), intent(in) :: DT
        Real(8)             :: displacement1_dt, displacement2_dt
        
       
    end Subroutine calculateVariables
    Subroutine WriteBubbleFiles(TIME)
        Use Physical_module!, only: Pressure_Bubble, calculateFlowrate
        Use BoundaryConditions
        Use Formats
        Implicit none
        Real(8), Intent(In) :: Time
        Real(8)             :: flowrate           
       

        ! write(20,'(4(f20.15,3x))') TIME, flowrate,        time_char*time, volumetric_flowrate*flowrate
        
       
    End Subroutine WriteBubbleFiles
end Module BubbleOutput
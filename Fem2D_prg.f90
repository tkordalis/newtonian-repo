! ***********************************************************************
!   diffusion in a quadrangle: left wall has high c1, right wall low c2
!   upper and lower walls have a flux imposed
! ***********************************************************************
PROGRAM FEM2D
   Use unv
   Use FieldFunctions
   Use FileModule
   Use IO_module
   Use BoundaryConditions
   Use BOUNDARY_ENUMERATION_MODULE
   Use ELEMENTS_MODULE
   Use NRAPSHON_MODULE
   Use MESH_MODULE
   Use GAUSS_MODULE
   Use ENUMERATION_MODULE
   Use BOUNDARY_ENUMERATION_MODULE
   Use CSR_STORAGE
   Use FLOW_ARRAYS_MODULE
   Use GLOBAL_ARRAYS_MODULE
   Use PHYSICAL_MODULE
   Use CONTINUATION_MODULE
   Use TIME_INTEGRATION
   Use Formats
   Use Tecplot
   Use MeshGeneration
   ! Use RemeshProcedure
   Use BubbleOutput
   Use InitialConditions
   USE OMP_PARALLEL
   Use RemeshVariables
   use newton_bulk_call
   Implicit None
   character(25)                        :: dateNtime
   
   Type(unvFileReader)                  :: unvf
   Type(Solution)                       :: sol
   Type(Solution)                       :: sol_o
   Type(Solution)                       :: sol_b
   Logical                              :: ReallocateForRemesh = .false.

   Logical                              :: ReadSolutionFromFile = .false.
   ! Logical                            :: ReadSolutionFromFile = .true.

   

   call openBubbleFiles()
      
   ! Read Solution
   ! unvf  = unvFileReader("UnBounded.unv")
   unvf  = unvFileReader("Bounded.unv")

   if (ReadSolutionFromFile) Then
      sol_b = ReadSolution("./sol/time_3.2000.plt")
      sol_o = ReadSolution("./sol/time_3.4000.plt")
      sol   = ReadSolution("./sol/time_3.6000.plt")
   else
      Remesh_counter_structuredInTheFront = 0
      Remesh_counter = 0
   end if
   call unvf%info()

   call DISCRETIZATION_DATA( unvf )
   call ALLOCATE_MESH_ARRAYS(.TRUE.,NODTOL)
   call MESH ( unvf )

   CALL GAUSS_EVALUATION

   call ALLOCATE_ENUMERATION_INDICES(.TRUE.,NEL_2d,NBF_1d,NBF_2d,NED_2d)
   call NNUM_ELEMENT( unvf )
   call SURROUNDING_NUMBERING
   call NNUM_f
   call ALLOCATE_BOUNDARY_ENUMERATION_INDICES(.TRUE.,NEL_2d)
   call DEFINE_BOUNDARY_NUMBERING( unvf )

   call unvf%close()
  !  ! ----------------------------------------------------------------------
  !  ! CSR STORAGE
  !  ! ----------------------------------------------------------------------
   call CSR_MAIN('NU')
   call CSR_MAIN('CO')
   call INITIALIZE_PARDISO 
  !  ! ----------------------------------------------------------------------
  !  ! ALLOCATE MEMORY FOR GLOBAL SYSTEM OF EQUATIONS
  !  ! ----------------------------------------------------------------------
   call ALLOCATE_CONTINUATION_ARRAYS( .TRUE., NUNKNOWNS_f           )
   call ALLOCATE_FLOW_ARRAYS        ( .TRUE., NUNKNOWNS_f, NEX_f    )
   call ALLOCATE_GLOBAL_ARRAYS      ( .TRUE., NODTOL, NEL_2d, NEQ_f )



   Call DefineTheBoundaries()

   CALL DIMENSIONLESS_NUMBERS

   call set_num_threads
   call set_DT

  !  ! ----------------------------------------------------------------------
  !  ! ASSIGN INITIAL CONDITIONS
  !  ! ----------------------------------------------------------------------
   call setInitalConditions()

   INCREMENT      = 0
   INCREMENT_STEP = 0
   TIME           = 0.0D0
   DT             = Dt_constant
   Dto            = DT
   Dtb            = DT

   if (ReadSolutionFromFile) then
      ! call sol_b%getSolutionVars(timeb, TLb, increment, Remesh_counter, pressure_bubble  )
      ! call sol_o%getSolutionVars(timeo, TLo, increment, Remesh_counter, Pressure_Bubbleo )
      ! call sol  %getSolutionVars(time , TL , increment, Remesh_counter, pressure_bubble  )
   endif
   
   call WriteBubbleFiles(TIME)

   if (.not. ReadSolutionFromFile) call exportFiles( TL,  NM_MESH, time, Increment, "POINT" )

   vm_ambient = 0.d0

   LOOP_TIME_INTEGRATION: DO
      
      INCREMENT         = INCREMENT + 1
      INCREMENT_STEP    = INCREMENT_STEP + 1
      TIME              = TIME + DT
      

      Write(*,*)'======================================================================'
      Write(*,*)'======================================================================'
      Write(*,*) replace("TIME STEP  =    *         TIME  =    *           DT  =    *", &
                    "*", [toStr(INCREMENT),         toStr(TIME),            toStr(DT)])
      Write(*,*)'======================================================================'
      Write(*,*)''


      Write(*,*)''

      
    ! SOLVE THE NONLINEAR PDE SYSTEM
      Call NEWTON_RAPSHON_f(ReallocateForRemesh)

      
      call WriteBubbleFiles(TIME)
      

      if ( (mod(increment,1) .eq. 0) ) then
      call exportFiles( TL,  NM_MESH, time, Increment, "POINT" )
      ! call exportFiles( TL,  NM_MESH, time, Increment, "BLOCK" )  ! mesh quality plt 
      endif
      
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      ! Check Criteria for Remeshing  
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      ! call checkAndRemesh( TL, NM_MESH,  Xm, Ym, Increment, ReallocateForRemesh )

      call UPDATE_SOLUTION( INCREMENT )
      
      call fdate(dateNtime)
      print*, dateNtime
      
      

    ENDDO LOOP_TIME_INTEGRATION
    
    
   ! ----------------------------------------------------------------------
   !  RELEASE OF MEMORY
   ! ----------------------------------------------------------------------
   CALL ALLOCATE_MESH_ARRAYS                 (.FALSE.,NODTOL)
   CALL ALLOCATE_ENUMERATION_INDICES         (.FALSE.,NEL_2d,NBF_1d,NBF_2d,NED_2d)
   CALL ALLOCATE_CONTINUATION_ARRAYS         (.FALSE., NUNKNOWNS_f )
   CALL ALLOCATE_FLOW_ARRAYS                 (.FALSE.,NUNKNOWNS_f, NEX_f)
   CALL ALLOCATE_GLOBAL_ARRAYS               (.FALSE., NODTOL, NEL_2d, NEQ_f )
   CALL ALLOCATE_CSR_ARRAYS                  (.FALSE., 1, NUNKNOWNS_f, NZ_f, NEX_f )
   CALL ALLOCATE_CSR_CONNECTIVITY            (.FALSE.,NBF_2d,NEL_2d)
   CALL ALLOCATE_BOUNDARY_ENUMERATION_INDICES(.FALSE.,NEL_2d)

END PROGRAM FEM2D
! ---------------------------------------------------------------------2
! ---------------------------------------------------------------------2
!    FEM IS THE MAIN PROGRAM FOR THE SOLUTION OF A SYSTEM OF PDEs
!    IN A CUBIC DOMAIN X,Y,Z
!    IT USES TRILINEAR BASIS FUNCTIONS & HEXAHEDRAL ELEMENTS
! ----------------------------------------------------------------------
!    PROGRAM WRITTEN BY J.A.T. & A.N.B. (1985, MIT    )
!    REVISED BY Y.D.                    (2013, UPATRAS)
!    REVISED BY S.E.V.                  (2018, UPATRAS)
!    REVISED BY D.P.                    (2020, UPATRAS)
! ----------------------------------------------------------------------

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
   Use RemeshProcedure
   Use BubbleOutput
   Use InitialConditions
   USE OMP_PARALLEL
   Use RemeshVariables
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
   unvf  = unvFileReader("UnBounded.unv")


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


   CALL DIMENSIONLESS_NUMBERS
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


   call set_num_threads
   call set_DT

  !  ! ----------------------------------------------------------------------
  !  ! ASSIGN INITIAL CONDITIONS
  !  ! ----------------------------------------------------------------------
   Call DefineTheBoundaries()
   call setInitalConditions()

   INCREMENT      = 0
   INCREMENT_STEP = 0
   TIME           = 0.0D0
   DT             = Dt_constant
   Dto            = DT
   Dtb            = DT

   if (ReadSolutionFromFile) then
      call sol_b%getSolutionVars(timeb, TLb, increment, Remesh_counter, pressure_bubble  )
      call sol_o%getSolutionVars(timeo, TLo, increment, Remesh_counter, Pressure_Bubbleo )
      call sol  %getSolutionVars(time , TL , increment, Remesh_counter, pressure_bubble  )
   endif
   


   if (.not. ReadSolutionFromFile) call exportFiles( TL,  NM_MESH, time, Increment, "POINT" )


   LOOP_TIME_INTEGRATION: DO
      
      INCREMENT         = INCREMENT + 1
      INCREMENT_STEP    = INCREMENT_STEP + 1
      TIME              = TIME + DT
      
      vm_ambient        = -( calculateFlowrate(time) ) / ( (Rtank)**2 - 1.d0 )
      ambient_position  = ambient_position + vm_ambient*dt


      Write(*,*)'======================================================================'
      Write(*,*)'======================================================================'
      Write(*,*) replace("TIME STEP  =    *         TIME  =    *           DT  =    *", &
                    "*", [toStr(INCREMENT),         toStr(TIME),            toStr(DT)])
      Write(*,*)'======================================================================'
      Write(*,*)''
      print*, "vm_ambient=",vm_ambient
      print*, "ambient_position=",ambient_position
      print*, "Rtank=",Rtank
      Write(*,*)''

      
    ! SOLVE THE NONLINEAR PDE SYSTEM
      Call NEWTON_RAPSHON_f(ReallocateForRemesh)

      call calculateBubbleVariables(DT)
      
      call WriteBubbleFiles(TIME)
      

      if ( (mod(increment,10) .eq. 0) ) then
      call exportFiles( TL,  NM_MESH, time, Increment, "POINT" )
      call exportFiles( TL,  NM_MESH, time, Increment, "BLOCK" )  ! mesh quality plt 
      endif
      
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      ! Check Criteria for Remeshing  
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      call checkAndRemesh( TL, NM_MESH,  Xm, Ym, Increment, ReallocateForRemesh )

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








!---------------------------------------------------------------------
!---------------------------------------------------------------------
!                  SUBROUTINE   NEWTON_RAPSHON_f
!---------------------------------------------------------------------

SUBROUTINE NEWTON_RAPSHON_f(ReallocateForRemesh)
   Use BoundaryConditions
   Use IO_module
   Use Formats
   Use ELEMENTS_MODULE
   Use CONTINUATION_MODULE,  only: ArN
   Use NRAPSHON_MODULE
   Use FLOW_ARRAYS_MODULE
   Use GLOBAL_ARRAYS_MODULE, only: TL, TLo, TLp, DpL
   Use CSR_STORAGE
   Use OMP_PARALLEL
   
   Use ExtraEquations    
   Use DirichletBoundaries
   Use ArrayTools
   Use solveAllExtraConstraints
   ! use BoundaryConditions
   use enumeration_module, only: gntr
   use VariableMapping
   Use RemeshVariables
   Use RemeshProcedure, only: AfterRemesh_counter
   IMPLICIT NONE
   
   Logical,  Intent(in):: ReallocateForRemesh


   ! LOCAL VARIABLES
   INTEGER          :: I, J, IERROR, IEL, K, INFO, N_dense
   INTEGER          :: xNITER, BvNiter
   CHARACTER(LEN=1) :: CCL, FoM
   LOGICAL          :: EMERGENCY, LMSG, FG_WS, changeBvN
   logical          :: MNR_FAILED
   REAL(8)          :: xERROR_NR, xF
   real(8)          :: Res_Norm_First_Iteration
   
   INTEGER, DIMENSION(NEX_f)  :: IPVT
   Real(8) :: Volume
   Real(8) :: Centroid 
   Character(len=:)     , Allocatable :: filename
   Character(len=:)     , Allocatable :: title
   Character(len=:)     , Allocatable :: zone

   Integer, Dimension(:), Allocatable :: elements 
   Integer, Dimension(:), Allocatable :: faces 
   Integer                            :: nelements, ii, jj, kk
   Integer, Dimension(:), Allocatable :: all_nodes

   Integer, Dimension(:)  , Allocatable :: globnodes
   Real(8), Dimension(:)  , Allocatable :: Z
   Real(8), Dimension(:)  , Allocatable :: R

   integer              :: n_timer, n_timer_o, time, time_in_seconds, counter_print

   ! INITIALIZE NEWTON LOOP VARIABLES
   ERR_nr     = 0
   EMERGENCY  = .FALSE.
   MNR_FAILED = .FALSE.
   xF         = 1.D0
   Res_Norm_First_Iteration = 0.d0
   300 CONTINUE

   IF(EMERGENCY)THEN
     xNITER = 50*NITER
     xF     = 0.5D0*xF
     IF ( (xF.LT.1.d-2) .or. (ITER_f .ge. 130) ) THEN
       WRITE(*,*) 'VERY SMALL RELAXATION FACTOR, GENERAL STOP!'
       STOP
       
     ENDIF

     xERROR_NR          = 1.D+0*ERROR_NR
     TL                 = TLp
     Pressure_Bubble    = Pressure_Bubbleo
     WRITE(*,*) 'xF=', xF
   ELSE
     xNITER = NITER
     xF        = 1D0
     xERROR_NR = ERROR_NR
   ENDIF
   
   
   ITER_f = 0

   RES_NORM       = 1.D0
   COR_NORM_OLD_f = 1.D0
   COR_NORM_NEW_f = 1.D0

    
   if (ReallocateForRemesh) then
    call DefineTheBoundaries()
  endif

   n_timer_o = time()
   n_timer   = 0
   ITER_f = 0
   ! ----------------------------------------------------------------------
   !  NEWTON ITERATION LOOP
   !  UPDATE THE SOLUTION VECTOR TL AND CHECK CONVERGENCE
   ! ----------------------------------------------------------------------
   LOOP_NEWTON_RAPSHON_f: DO WHILE ((COR_NORM_NEW_f.GT.xERROR_NR) .or. (RES_NORM.GT.xERROR_NR))

     
     

     ITER_f = ITER_f + 1
     if (COR_NORM_NEW_f <= 1.d0) xF = 1.10d0 * xF
     if (xF             >= 1.d0) xF = 1.d0 

     CHECK_EMERGENCY: IF(EMERGENCY)THEN
                        FLAG_NR='NRP'
                      ENDIF CHECK_EMERGENCY
     IF (MNR_FAILED) FLAG_NR='NRP'
     if (AfterRemesh_counter .eq. 1) FLAG_NR='NRP' 

     
     ! INITIALIZE LINEAR SYSTEM'S MATRICES & VECTORS
     B_f  = 0.D0
     Be_f = 0.D0
     Bi_f = 0.D0
     Bw_f = 0.D0
     S_f  = 0.D0
     Se_f = 0.D0
     Sa_f = 0.D0
     Sb_f = 0.D0
     IF (FLAG_NR=='NRP') THEN
       A_f  = 0.D0
       Ah_f = 0.D0
       Ar_f = 0.D0
       Ac_f = 0.D0
       Ai_f = 0.D0
     ENDIF
  
      ! print*, bubble1%pressure
      ! pause
     !$OMP  PARALLEL DO NUM_THREADS(NTHREADS)&
     !$OMP& DEFAULT (SHARED)&
     !$OMP& PRIVATE (IEL)
      DO IEL = 1, NEL_2d
       CALL FLOW_EQUATIONS(IEL, FLAG_NR)
      ENDDO
     !$OMP END PARALLEL DO
     

     
      call applyBCs_and_solveExtraConstraints( FLAG_NR, Ah_f, Be_f, Pressure_Bubble )
     

     ! CALCULATE RESIDUAL NORM
     RES_NORM = DOT_PRODUCT(B_f,B_f)+DOT_PRODUCT(Be_f,Be_f)
     RES_NORM = DSQRT(RES_NORM)
     
     ! print *, DOT_PRODUCT(B_f,B_f)
     ! print *, DOT_PRODUCT(Be_f,Be_f)
     jj=1
     kk=0
     do ii = 1, size(B_f)
        if (jj .gt. NEQ_f) jj=1
        
        if (jj .eq. 1) then
          kk=kk+1
          ! print*,'global_node=', kk
        endif
        
        ! if ( getVariableName(jj) == 'P' ) then
        ! print*,'global_node=', kk
        ! print*, ' ' 
        ! print*,'(X,Y) =', Xm(kk), Ym(kk) 
        ! print*, ' ' 
        ! print*, 'variable','     ', 'residual' 
        ! print*, getVariableName(jj), '     ',B_f(ii) 
        ! print*, ' ' 
        ! print*, ' ' 
        ! print*, ' //////////////////////////////////////////////////////// ' 
        ! pause
        ! endif
        jj=jj+1
     enddo
     ! ----------------------------------------------------------------------
     !    LU DECOMPOSITION
     ! ----------------------------------------------------------------------
     IF(FLAG_NR=='NRP')THEN
       
       PHASE_f = 22  ! ONLY FACTORIZATION
       NRHS_f  = 1

       CALL PARDISO&
       (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
        IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, DDUM_f, DDUM_f, ERROR_f)

       IF(ERROR_f .NE. 0)THEN

         WRITE(*,*)'THE FOLLOWING ERROR_F WAS DETECTED IN LU DECOMPOSITION: ', ERROR_f
         WRITE(*,*)' '
         STOP
       ENDIF
       
     ENDIF
     
      
     !--------------------------------------------------------------
     !    CALCULATE Sa_f = A_f(-1).Ac_f
     !--------------------------------------------------------------
     PHASE_f     = 33  ! ONLY SOLUTION
     IPARM_f(8)  = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS
     NRHS_f      = Nex_f


     CALL PARDISO&
     (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
     IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, Ac_f, Sa_f, ERROR_f)
     

     !--------------------------------------------------------------
     !    CALCULATE Sb_f = A_f(-1).B_f
     !--------------------------------------------------------------
     PHASE_f     = 33  ! ONLY SOLUTION
     IPARM_f(8)  = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS
     NRHS_f      = 1

     CALL PARDISO&
     (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
     IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, B_f, Sb_f, ERROR_f)
     
     !--------------------------------------------------------------      
     !    SOLVE INTERMEDIATE SYSTEM Ai_f.Se_f = Bi_f
     !--------------------------------------------------------------
     DO I = 1, SIZE(Ai_f,1)
       DO J = 1, SIZE(Ai_f,2)
         Ai_f(I,J) = DOT_PRODUCT(Ar_f(I,:),Sa_f(:,J)) - Ah_f(I,J)
       ENDDO
     ENDDO
     

     DO I = 1, SIZE(Bi_f)
       Bi_f(I) = DOT_PRODUCT(Ar_f(I,:),Sb_f) - Be_f(I)
     ENDDO
     
     N_dense = Nex_f

     ! PERFORM LU DECOMPOSITION AND CALCULATE CONDITIONING NUMBER
     CALL DGETRF( N_dense, N_dense, Ai_f, N_dense, IPVT, INFO )
     IF (INFO.NE.0) PRINT*, '[Error] NEWTON_RAPSHON_f / Dense matrix is singular to working precision'
     IF (INFO.NE.0) STOP

     ! PERFORM BACK-SUBSTITUTION
     CALL DGETRS( 'N', N_dense, 1, Ai_f, N_dense, IPVT, Bi_f, N_dense, INFO )
     
     ! CORRECTION VECTOR
     Se_f = Bi_f
     
      
     !-------------------------------------------------------------
     !    CALCULATE WHOLE SOLUTION
     !-------------------------------------------------------------
     DO I = 1, SIZE(Bw_f)
       Bw_f(I) = B_f(I) - DOT_PRODUCT(Ac_f(I,:),Se_f)
     ENDDO

     PHASE_f     = 33  ! ONLY SOLUTION
     IPARM_f(8)  = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS

     ! CALL PARDISO&
     ! (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
     ! IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, Bw_f, S_f, ERROR_f)

     
     DO I = 1, SIZE(S_f)
       S_f(I) = Sb_f(I) - DOT_PRODUCT(Sa_f(I,:),Se_f)
     ENDDO
     
      n_timer         = time()
      time_in_seconds = n_timer - n_timer_o
      n_timer_o       = n_timer
     ! CHECK CONVERGENCE & CHOOSE KIND OF NR METHOD
     CALL CHECK_CONVERGENCE&
     (ITER_f, xNITER, FLAG_NR, CCL, FoM, COR_NORM_NEW_f, COR_NORM_OLD_f, S_f, Se_f,&
      NUNKNOWNS_f, NEX_f, RES_NORM, LMSG, MNR_FAILED, Res_Norm_First_Iteration, time_in_seconds)

     IF (LMSG) THEN

        EMERGENCY = .TRUE.
        
          GOTO 300

     ENDIF

     IF (CCL=='Y') CYCLE LOOP_NEWTON_RAPSHON_f
     ERR_nr = max(ERR_nr,COR_NORM_NEW_f)
     
     
     ! UPDATE SOLUTION
     K = 0
     DO I = 1, NODTOL
       DO J = 1, NEQ_f
         K = K + 1
         TL(I,J) = TL(I,J) - xF*S_f(K)
       ENDDO
     ENDDO
     
     
      Pressure_Bubble    = Pressure_Bubble     - xF * Se_f(1)
      call Bubble1%setPressure(Pressure_Bubble)

            
      filename = replace("Iteration_*.plt","*", toStr(ITER_f) ) 
      title = replace("Bubble_Pressure_1 = *", "*", toStr(Pressure_Bubble))

      zone     = toStr(ITER_f)

      call WriteTecplotFile(filename, title, zone, TL, NM_MESH)

   ENDDO LOOP_NEWTON_RAPSHON_f

   call Bubble1%setPressure(Pressure_Bubble)
   ! call Bubble1%setMaxRcoordinateOfInterfaceAndLogicalOperator(StructMeshOnlyInTheFront, TL)

   StructMeshOnlyInTheFront = .false.
   do k = 1, size(Bubble1%nodes)
        kk = Bubble1%nodes(k)
            if (TL(kk, getVariableId("R")) .gt. 1.15d0) then 
                print*, 'global node and R coordinate =',kk, TL(kk, getVariableId("R"))
                StructMeshOnlyInTheFront = .true.
                 print*, StructMeshOnlyInTheFront
                exit
            endif
        enddo

   ! If convergence is achieved, remove iteration_*.plt
   call execute_command_line("rm -f Iteration_*.plt")
    
   

END SUBROUTINE NEWTON_RAPSHON_f


Subroutine loopOverElements(nelements, elements, faces, gid, procedure_, globUnknown)
  Use CSR_STORAGE,               Only: Ar_f, Ah_f
  Use constrainJacobians,        only: jacobianOfConstrain

  Implicit None
  Integer,                       Intent(In) :: nelements
  Integer, Dimension(nelements), Intent(In) :: elements
  Integer, Dimension(nelements), Intent(In) :: faces
  Real(8),             intent(in), optional :: globUnknown
  Integer                                   :: gid
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  Interface
    Function procedure_ (nelem, ned) Result(out)
      Implicit None 
      Integer, Intent(In) :: nelem
      Integer, Intent(In) :: ned
      Real(8)             :: out
    End Function procedure_
  End Interface 
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  Integer                                   :: iel
  Integer                                   :: element
  Integer                                   :: face

     
  do iel = 1, nelements
    element = elements(iel)
    face    = faces   (iel)

    call jacobianOfConstrain(gid, element, face, procedure_)
  end do   
  
End Subroutine loopOverElements


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                 SUBROUTINE FLOW_EQUATIONS
!-----------------------------------------------------------------------


SUBROUTINE FLOW_EQUATIONS(IEL, FLAG_NR)
 
   Use PHYSICAL_MODULE
   Use ELEMENTS_MODULE,             Only: NBF_2d, NEQ_f
   Use ENUMERATION_MODULE,          Only: NM_MESH
   Use BOUNDARY_ENUMERATION_MODULE, Only: NBE
   Use GLOBAL_ARRAYS_MODULE,        Only: TL
   Use MESH_MODULE,                 Only: Xm, Ym, EPS_MESH
   Use Boundary_Equations 
   Use BulkEquations 
   Use NumericalBoundaryJacobian
   Implicit None
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   !  ARGUMENTS
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   Integer          :: IEL, ICH
   Character(len=3) :: FLAG_NR
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
   !  LOCAL VARIABLES
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   Integer :: INOD, II, IBND, NOD
   Logical :: FLAG
   Real(8) :: SGN, XX, YY
   Real(8), Dimension(NBF_2d,NEQ_f) :: TEMP_TL
   Real(8), Dimension(NBF_2d,NEQ_f) :: TEMP_RES
   Integer                          :: face


   ! COPY SOLUTION TO A LOCAL ARRAY
   TEMP_TL = 0.D0
   DO INOD = 1, NBF_2d
     II = NM_MESH(IEL,INOD)
     TEMP_TL(INOD,:) = TL(II,:)
   ENDDO


   CALL DOMI_RESIDUAL_f( IEL, TEMP_TL, TEMP_RES, .TRUE. )
       
   IF (FLAG_NR=='NRP') CALL DOMI_JACOBIAN_f( IEL, TEMP_TL, TEMP_RES )
   
 
   
END SUBROUTINE FLOW_EQUATIONS


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                 SUBROUTINE CHECK_CONVERGENCE
!-----------------------------------------------------------------------


SUBROUTINE CHECK_CONVERGENCE&
(ITER, MITER, FLAG_NR, CCL, FoM, RSUM_NEW, RSUM_OLD, B, Be, IDIM_B, IDIM_Be,&
       RES_NRM, LMSG, MNR_FAILED, Res_Norm_First_Iteration, time_in_seconds )

   USE NRAPSHON_MODULE,         only: ERROR_NR, NITER
   use TIME_INTEGRATION,        only: TIME
   IMPLICIT NONE

   ! ARGUMENTS
   LOGICAL,          INTENT(INOUT) :: LMSG, MNR_FAILED
   INTEGER,          INTENT(IN)    :: ITER, MITER
   CHARACTER(len=3), INTENT(INOUT) :: FLAG_NR
   CHARACTER(len=1), INTENT(OUT)   :: CCL
   CHARACTER(len=1), INTENT(IN)    :: FoM
   REAL(8),          INTENT(OUT)   :: RSUM_NEW
   REAL(8),          INTENT(INOUT) :: RSUM_OLD, RES_NRM, Res_Norm_First_Iteration

   INTEGER,                     INTENT(IN) :: IDIM_B, IDIM_Be
   REAL(8), DIMENSION(IDIM_B),  INTENT(IN) :: B
   REAL(8), DIMENSION(IDIM_Be), INTENT(IN) :: Be
   INTEGER   ,                  INTENT(IN) :: time_in_seconds


   ! LOCAL VARIABLES
   REAL(8)  :: ERROR_NEW
   INTEGER  :: I, IJ
   real(8), PARAMETER  :: Weak_Criterion_MNR   = 1.d-3
   real(8), PARAMETER  :: Strong_Criterion_MNR = 1.d-6
   real(8)  :: Criterion_MNR
   LMSG = .FALSE.

   CCL = 'N'
   

   RSUM_NEW  = DSQRT(DOT_PRODUCT(B,B)+DOT_PRODUCT(Be,Be))
   ERROR_NEW = MAX(MAXVAL(DABS(B)),MAXVAL(DABS(Be)))
   If (ITER==1)   Res_Norm_First_Iteration = RSUM_NEW
   Criterion_MNR = Weak_Criterion_MNR
   If (TIME.GT.10.89d0) Criterion_MNR = Strong_Criterion_MNR
   IF (RSUM_NEW.NE.RSUM_NEW) THEN
     WRITE(*,*) 'NaN ENCOUNTERED, GENERAL STOP'
     LMSG = .TRUE.
     RETURN
   ENDIF
   
   
   ! SELECT AMONG FULL AND MODIFIED NEWTON-RAPHSON
   KIND_OF_NEWTON_RAPHSON: SELECT CASE(FLAG_NR)
   
     CASE('NRP')
   ! ---------------------------------------------------------------------
     IF( RSUM_NEW .LT. 1.0D-1 )THEN

       IF( RSUM_NEW .GT. ERROR_NR ) THEN
         FLAG_NR = 'MNR'
       ELSE
         FLAG_NR = 'NRP'
       ENDIF

     ELSE

       FLAG_NR = 'NRP'

     ENDIF


     WRITE(*, 50)ITER, time_in_seconds
     WRITE(*, 51)RES_NRM, RSUM_NEW
    
      
     CASE('MNR')
   ! ---------------------------------------------------------------------
     IF( RSUM_NEW .GE. 1.0D-2 )THEN

       FLAG_NR = 'NRP'
       WRITE(*, 52)ITER, time_in_seconds
       WRITE(*, 51)RES_NRM, RSUM_NEW
       MNR_FAILED = .TRUE.
       CCL = 'Y'
       RETURN
      
     ELSE
       IF(RSUM_NEW/RSUM_OLD .GE.0.5D0)THEN
         MNR_FAILED = .TRUE.
         FLAG_NR = 'NRP'
       ELSE
         FLAG_NR = 'MNR'
       ENDIF

     ENDIF
         
     IF (( RSUM_NEW .LE. ERROR_NR ).AND.(Res_Norm_First_Iteration.GT.1.d-3)) THEN
       FLAG_NR = 'NRP'
     ENDIF

     IF( (NITER - 7) .LT. ITER ) MNR_FAILED = .TRUE.
     WRITE(*, 52)ITER, time_in_seconds
     WRITE(*, 51)RES_NRM, RSUM_NEW
     
     
     CASE DEFAULT
     WRITE(*,*)' INCORRECT CHOICE OF NEWTON RAPHSON FLAG '


     END SELECT KIND_OF_NEWTON_RAPHSON
     
      
    ! ONLY FULL NEWTON RAPHSON 
    ! FLAG_NR = 'NRP'


   IF( RSUM_NEW .GT. 5.0D+9 )THEN
     WRITE(*,*)'PROGRAM UNABLE TO CONVERGE'
     WRITE(*,*)'TOO LARGE NORMA !!!'
     WRITE(*,*)'OVER-RELAXATION NR ITERATIONS !!!'
     LMSG = .TRUE.
     RETURN
   ENDIF


   IF( ITER .GT. MITER )THEN
     WRITE(*,*)'PROGRAM UNABLE TO CONVERGE'
     WRITE(*,*)'TOO MANY ITERATIONS !!!'
     WRITE(*,*)'OVER-RELAXATION NR ITERATIONS !!!'
     LMSG = .TRUE.
     RETURN
   ENDIF
     

   RSUM_OLD = RSUM_NEW
     
     
   !  FORMAT STATEMENTS
   50 FORMAT(I3,'. ITERATION. FULL NEWTON RAPSHON - Elapsed Time = ', I4, ' seconds' )
   51 FORMAT(3x,'  RESIDUAL NORM = ', E12.4, 7x, 'CORRECTION NORM = ', E12.4)
   52 FORMAT(I3,'. ITERATION. MODIFIED NEWTON RAPSHON - Elapsed Time = ', I4, ' seconds')


END SUBROUTINE CHECK_CONVERGENCE
      






! ----------------------------------------------------------------------
! ----------------------------------------------------------------------






!----------------------------------------------------------------------
!----------------------------------------------------------------------
!               REAL(8) FUNCTION F_DX
!
!----------------------------------------------------------------------


REAL(8) FUNCTION F_DX ( X )

   IMPLICIT NONE
      
   ! ARGUMENTS
   REAL(8), INTENT(IN)  :: X
      
   ! LOCAL VARIABLES
   REAL(8) :: EP_RES = 1.0D-9
      
      
   F_DX = EP_RES*DMAX1(1.0D0,DABS(X))*SIGN(1.0D0,X)

      
END FUNCTION F_DX


  

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                   SUBROUTINE BASIS_2d
!----------------------------------------------------------------------
      
SUBROUTINE BASIS_2d&
( IG, X_loc, Y_loc, BFN_2d, DFDC_2d, DFDE_2d, XX, DXDC, DXDE, YY, DYDC,&
   DYDE, CJAC, AJAC, DFDX, DFDY, NGAUSS )

   USE ELEMENTS_MODULE, only: NBF_2d

   IMPLICIT NONE
   ! ARGUMENTS
   INTEGER, INTENT(IN)                            :: IG, NGAUSS
   REAL(8), INTENT(IN), DIMENSION(NBF_2d)         :: X_loc, Y_loc
   REAL(8), INTENT(IN), DIMENSION(NBF_2d, NGAUSS) :: BFN_2d, DFDC_2d, DFDE_2d
  
   REAL(8), INTENT(INOUT) :: XX, DXDC, DXDE
   REAL(8), INTENT(INOUT) :: YY, DYDC, DYDE
   REAL(8), INTENT(INOUT) :: CJAC, AJAC
  
   REAL(8), INTENT(INOUT), DIMENSION(NBF_2d) :: DFDX, DFDY

   ! LOCAL VARIABLES
   INTEGER :: I
   REAL(8) :: X,  Y
   REAL(8) :: Cx, Cy
   REAL(8) :: Ex, Ey
      
      
   ! CALCULATE PARTIAL DERIVATIVES OF TRANSFORMATION
   ! AT A POINT C,E IN NELEM

   XX = 0.D0 ; DXDC = 0.D0 ; DXDE = 0.D0

   DO I = 1, NBF_2d
     X    = X_loc(I)
     XX   = XX   +  BFN_2d(I,IG)*X
     DXDC = DXDC + DFDC_2d(I,IG)*X
     DXDE = DXDE + DFDE_2d(I,IG)*X
   ENDDO


   YY = 0.D0 ; DYDC = 0.D0 ; DYDE = 0.D0

   DO I = 1, NBF_2d
     Y    = Y_loc(I)
     YY   = YY   +  BFN_2d(I,IG)*Y
     DYDC = DYDC + DFDC_2d(I,IG)*Y
     DYDE = DYDE + DFDE_2d(I,IG)*Y
   ENDDO


   ! CALCULATE JACOBIAN OF THE TRANSFORMATION
   CJAC = DXDC*DYDE-DXDE*DYDC
   AJAC = DABS(CJAC)
  
   ! CALCULATE DERIVATIVES OF BASIS FUNCTIONS WRT X,Y COORDINATES
   DO I = 1, NBF_2d
     DFDX(I) = (DFDC_2d(I,IG)*DYDE-DFDE_2d(I,IG)*DYDC)/CJAC
     DFDY(I) = (DFDE_2d(I,IG)*DXDC-DFDC_2d(I,IG)*DXDE)/CJAC
   ENDDO
      
  
END SUBROUTINE BASIS_2d


 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!  SUBROUTINE UPDATE_SOLUTION
!--------------------------------------------------------------------------

SUBROUTINE UPDATE_SOLUTION( INCREMENT )

   USE PHYSICAL_MODULE,  only: Pressure_Bubble, Pressure_Bubbleo 


   USE GLOBAL_ARRAYS_MODULE
   USE TIME_INTEGRATION,     only: DT, TIME, DTo, DTb
   IMPLICIT NONE

   ! ARGUMENTS
   INTEGER, INTENT(IN)        :: INCREMENT

   ! LOCAL VARIABLES
   REAL(8) :: Lb, Lo, L

  if (INCREMENT.GT.2) then 
   CALL LAGRANGE_EXTRAPOLATION(TIME+dt, TIME-(dtb+DTo), TIME-dto, TIME, Lb, Lo, L)

   TLp = Lb*TLb + Lo*TLo + L*TL
   TLb = TLo
   TLo = TL
   TL  = TLp

  
  else 
   ! UPDATE VECTORS
   TLp = TL  + Dt*(TL-TLo)
   TLb = TLo
   TLo = TL
   TL  = TLp

 
  endif 

  Pressure_Bubbleo = Pressure_Bubble

END SUBROUTINE UPDATE_SOLUTION




SUBROUTINE LAGRANGE_EXTRAPOLATION(Xp, Xb, Xo, X, Lb, Lo, L)
  IMPLICIT NONE

  ! ARGUMENTS
  REAL(8), INTENT(IN)        :: Xp, Xb, Xo, X
  REAL(8), INTENT(OUT)       :: Lb, Lo, L

  
  Lb = (Xp-Xo)*(Xp-X) /((Xb-Xo)*(Xb-X))

  Lo = (Xp-Xb)*(Xp-X) /((Xo-Xb)*(Xo-X))

  L  = (Xp-Xb)*(Xp-Xo)/((X -Xb)*(X -Xo))
 

END SUBROUTINE LAGRANGE_EXTRAPOLATION


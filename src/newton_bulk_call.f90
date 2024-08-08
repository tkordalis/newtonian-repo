module newton_bulk_call
    use solution_check_update
    contains

    SUBROUTINE NEWTON_RAPSHON_f(ReallocateForRemesh)
        Use BoundaryConditions
        Use IO_module
        Use Formats
        Use ELEMENTS_MODULE
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
        ! Use RemeshProcedure, only: AfterRemesh_counter
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


        Character(len=:)     , Allocatable :: filename
        Character(len=:)     , Allocatable :: title
        Character(len=:)     , Allocatable :: zone

        Integer, Dimension(:), Allocatable :: elements 
        Integer, Dimension(:), Allocatable :: faces 
        Integer                            :: nelements, ii, jj, kk
        Integer, Dimension(:), Allocatable :: all_nodes





        integer              :: n_timer, n_timer_o, time, time_in_seconds, counter_print

        ! INITIALIZE NEWTON LOOP VARIABLES
        ERR_nr     = 0
        EMERGENCY  = .FALSE.
        MNR_FAILED = .FALSE.
        xF         = 1.D0
        Res_Norm_First_Iteration = 0.d0
        300 CONTINUE

        IF(EMERGENCY)THEN
            xNITER = 10*NITER
            xF     = 0.5D0*xF
            IF ( (xF.LT.1.d-2) .or. (ITER_f .ge. 50) ) THEN
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

        call bubble%setCentroid_o()
        call bubble%setVolume_o()

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
            ! if (AfterRemesh_counter .eq. 1) FLAG_NR='NRP' 


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
        !   !$OMP  PARALLEL DO NUM_THREADS(NTHREADS)&
        !   !$OMP& DEFAULT (SHARED)&
        !   !$OMP& PRIVATE (IEL)
            DO IEL = 1, NEL_2d
                ! CALL CONCENTRATION_EQUATION(IEL, FLAG_NR)
                CALL FLOW_EQUATIONS(IEL, FLAG_NR)
            ENDDO
        !      !$OMP END PARALLEL DO
                pause

              ! jj=1
              ! kk=0
              ! do ii = 1, size(B_f)
              !    if (jj .gt. NEQ_f) jj=1
              !    if (jj .eq. 1) then
              !      kk=kk+1
              !    endif
              !    !! if ( getVariableName(jj) == 'P' ) then
              !      ! if ( abs(B_f(ii)) .gt. 1.d0) then 
              !      print*,'global_node=', kk
              !      print*, ' ' 
              !      print '(A9,2x,f10.5,2x,f10.5,2x,f10.5)','(X,Y,R) =', Xm(kk), Ym(kk), sqrt(Xm(kk)**2 + Ym(kk)**2)
              !      print*, ' ' 
              !      print*, 'variable','     ', 'residual' 
              !      print*, getVariableName(jj), '     ',B_f(ii) 
              !      print*, ' ' 
              !      print*, ' ' 
              !      print*, ' ' 
              !      print*, ' //////////////////////////////////////////////////////// ' 
              !      pause
              !      ! endif
              !    !! endif
              !    jj=jj+1
              ! enddo

            call applyBCs_solveExtraConstraints( FLAG_NR, Pressure_bubble )

           

            ! CALCULATE RESIDUAL NORM
            RES_NORM = DOT_PRODUCT(B_f,B_f)+DOT_PRODUCT(Be_f,Be_f)
            RES_NORM = DSQRT(RES_NORM)

            ! print *, DOT_PRODUCT(B_f,B_f)
            ! print *, DOT_PRODUCT(Be_f,Be_f)
            ! pause


            
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

            ! Ac_f = 0.d0
            ! Ar_f = 0.d0
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

            ! PHASE_f     = 33  ! ONLY SOLUTION
            ! IPARM_f(8)  = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS


            ! ****************************************************
            ! ****************************************************
            ! In case a problem has no extra equations I use this
            ! ****************************************************
            !Solve for S_f = A_f(-1).Bw_f
            ! IF(NEX_f == 0)THEN
            !     CALL PARDISO(PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
            !     IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, B_f, S_f, ERROR_f)
            ! ELSE
            !     CALL PARDISO(PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
            !     IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, Bw_f, S_f, ERROR_f)
            ! END IF
            ! ****************************************************
            ! ****************************************************
            
            PHASE_f    = 33  ! ONLY SOLUTION
            IPARM_f(8) = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS

            CALL PARDISO&
            (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
            IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, Bw_f, S_f, ERROR_f)


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

            ! print *, DOT_PRODUCT(S_f,S_f)
            ! print *, DOT_PRODUCT(Se_f,Se_f)
            ! pause


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
         ! print*, ' '
         ! print*, '------------------------------------ '
         ! print*, ' '
         ! print '(a5, 3x,i3)', "node=", I
         ! print*, ' '
         ! print '(a7, 3x, f10.5, 3x, f10.5, 3x, f10.5)', '(X,Y,R) =', Xm(I), Ym(I), sqrt(Xm(I)**2+Ym(I)**2)
         ! print*, ' '
         ! print '(a9, 3x,a2, 3x, a11, e20.10)', "variable=", getVariableName(J), "correction=", s_f(k)
         ! print*, ' '
         ! pause 

                ENDDO
            ENDDO


            Pressure_bubble = Pressure_bubble - xF*Se_f(1)

            filename = replace("Iteration_*.plt","*", toStr(ITER_f) ) 

            title = replace("Solution_*", "*", toStr(iter_f))

            zone     = toStr(ITER_f)

            call WriteTecplotFile(filename, title, zone, TL, NM_MESH)
            ! pause

        ENDDO LOOP_NEWTON_RAPSHON_f

        call bubble%setPressure(Pressure_bubble)

        ! call Bubble1%setMaxRcoordinateOfInterfaceAndLogicalOperator(StructMeshOnlyInTheFront, TL)


        ! If convergence is achieved, remove iteration_*.plt
        call execute_command_line("rm -f Iteration_*.plt")

           

    END SUBROUTINE NEWTON_RAPSHON_f


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

       CALL DOMI_RESIDUAL_fluid( IEL, TEMP_TL, TEMP_RES, .TRUE. )
       ! IF (FLAG_NR=='NRP') CALL DOMI_JACOBIAN_f( IEL, TEMP_TL, TEMP_RES )
       IF (FLAG_NR=='NRP') call NumJacBulk(  DOMI_RESIDUAL_fluid  ,IEL, TEMP_TL, TEMP_RES)
        
       
    END SUBROUTINE FLOW_EQUATIONS


    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------


    SUBROUTINE CONCENTRATION_EQUATION(IEL, FLAG_NR)
     
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


       ! CALL DOMI_RESIDUAL_chemSpecies( IEL, TEMP_TL, TEMP_RES, .TRUE. )
           
       ! IF (FLAG_NR=='NRP') CALL DOMI_JACOBIAN_f( IEL, TEMP_TL, TEMP_RES )
       ! IF (FLAG_NR=='NRP') call NumJacBulk(  DOMI_RESIDUAL_chemSpecies  ,IEL, TEMP_TL, TEMP_RES)
     
       
    END SUBROUTINE CONCENTRATION_EQUATION




end module newton_bulk_call
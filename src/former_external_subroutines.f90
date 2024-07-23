module former_external_subroutines

    contains

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
            xNITER = 50*NITER
            xF     = 0.5D0*xF
            IF ( (xF.LT.1.d-2) .or. (ITER_f .ge. 130) ) THEN
                WRITE(*,*) 'VERY SMALL RELAXATION FACTOR, GENERAL STOP!'
                STOP

            ENDIF

            xERROR_NR          = 1.D+0*ERROR_NR
            TL                 = TLp

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
            !$OMP  PARALLEL DO NUM_THREADS(NTHREADS)&
            !$OMP& DEFAULT (SHARED)&
            !$OMP& PRIVATE (IEL)
            DO IEL = 1, NEL_2d
                CALL CONCENTRATION_EQUATION(IEL, FLAG_NR)
                ! CALL FLOW_EQUATIONS(IEL, FLAG_NR)
            ENDDO
            !$OMP END PARALLEL DO



            call applyBCs_and_solveExtraConstraints( FLAG_NR )


            ! CALCULATE RESIDUAL NORM
            RES_NORM = DOT_PRODUCT(B_f,B_f)+DOT_PRODUCT(Be_f,Be_f)
            RES_NORM = DSQRT(RES_NORM)

            ! print *, DOT_PRODUCT(B_f,B_f)
            ! print *, DOT_PRODUCT(Be_f,Be_f)
            jj=1
            kk=0
            ! do ii = 1, size(B_f)
            !    if (jj .gt. NEQ_f) jj=1

            !    if (jj .eq. 1) then
            !      kk=kk+1
            !      ! print*,'global_node=', kk
            !    endif

            !    ! if ( getVariableName(jj) == 'P' ) then
            !      print*,'global_node=', kk
            !      print*, ' ' 
            !      print*,'(X,Y) =', Xm(kk), Ym(kk) 
            !      print*, ' ' 
            !      print*, 'variable','     ', 'residual' 
            !      print*, getVariableName(jj), '     ',B_f(ii) 
            !      print*, ' ' 
            !      print*, ' ' 
            !      print*, ' ' 
            !      print*, ' //////////////////////////////////////////////////////// ' 
            !      pause
            !      ! if ( abs(B_f(ii)) .gt. 1.d+2) then 
            !      ! endif
            !    ! endif
            !    jj=jj+1
            ! enddo
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
            ! PHASE_f     = 33  ! ONLY SOLUTION
            ! IPARM_f(8)  = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS
            ! NRHS_f      = Nex_f


            ! CALL PARDISO&
            ! (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
            ! IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, Ac_f, Sa_f, ERROR_f)


            !--------------------------------------------------------------
            !    CALCULATE Sb_f = A_f(-1).B_f
            !--------------------------------------------------------------
            ! PHASE_f     = 33  ! ONLY SOLUTION
            ! IPARM_f(8)  = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS
            ! NRHS_f      = 1

            ! CALL PARDISO&
            ! (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
            ! IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, B_f, Sb_f, ERROR_f)

            !--------------------------------------------------------------      
            !    SOLVE INTERMEDIATE SYSTEM Ai_f.Se_f = Bi_f
            !--------------------------------------------------------------
            ! DO I = 1, SIZE(Ai_f,1)
            !   DO J = 1, SIZE(Ai_f,2)
            !     Ai_f(I,J) = DOT_PRODUCT(Ar_f(I,:),Sa_f(:,J)) - Ah_f(I,J)
            !   ENDDO
            ! ENDDO


            ! DO I = 1, SIZE(Bi_f)
            !   Bi_f(I) = DOT_PRODUCT(Ar_f(I,:),Sb_f) - Be_f(I)
            ! ENDDO

            N_dense = Nex_f

            ! PERFORM LU DECOMPOSITION AND CALCULATE CONDITIONING NUMBER
            ! CALL DGETRF( N_dense, N_dense, Ai_f, N_dense, IPVT, INFO )
            ! IF (INFO.NE.0) PRINT*, '[Error] NEWTON_RAPSHON_f / Dense matrix is singular to working precision'
            ! IF (INFO.NE.0) STOP

            ! PERFORM BACK-SUBSTITUTION
            ! CALL DGETRS( 'N', N_dense, 1, Ai_f, N_dense, IPVT, Bi_f, N_dense, INFO )

            ! CORRECTION VECTOR
            ! Se_f = Bi_f


            !-------------------------------------------------------------
            !    CALCULATE WHOLE SOLUTION
            !-------------------------------------------------------------
            ! DO I = 1, SIZE(Bw_f)
            !   Bw_f(I) = B_f(I) - DOT_PRODUCT(Ac_f(I,:),Se_f)
            ! ENDDO

            ! PHASE_f     = 33  ! ONLY SOLUTION
            ! IPARM_f(8)  = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS

            PHASE_f    = 33  ! ONLY SOLUTION
            IPARM_f(8) = 2   ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS

            !Solve for S_f = A_f(-1).Bw_f
            IF(NEX_f == 0)THEN
                CALL PARDISO(PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
                IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, B_f, S_f, ERROR_f)
            ELSE
                CALL PARDISO(PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
                IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, Bw_f, S_f, ERROR_f)
            END IF

            ! CALL PARDISO&
            ! (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
            ! IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, Bw_f, S_f, ERROR_f)


            ! DO I = 1, SIZE(S_f)
            !   S_f(I) = Sb_f(I) - DOT_PRODUCT(Sa_f(I,:),Se_f)
            ! ENDDO

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




            filename = replace("Iteration_*.plt","*", toStr(ITER_f) ) 

            title = replace("Solution_*", "*", toStr(iter_f))

            zone     = toStr(ITER_f)

            call WriteTecplotFile(filename, title, zone, TL, NM_MESH)


        ENDDO LOOP_NEWTON_RAPSHON_f


        ! call Bubble1%setMaxRcoordinateOfInterfaceAndLogicalOperator(StructMeshOnlyInTheFront, TL)


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


       CALL DOMI_RESIDUAL_f( IEL, TEMP_TL, TEMP_RES, .TRUE. )
           
       ! IF (FLAG_NR=='NRP') CALL DOMI_JACOBIAN_f( IEL, TEMP_TL, TEMP_RES )
       IF (FLAG_NR=='NRP') call NumJacBulk(  DOMI_RESIDUAL_f  ,IEL, TEMP_TL, TEMP_RES)
     
       
    END SUBROUTINE CONCENTRATION_EQUATION


    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !                 SUBROUTINE CHECK_CONVERGENCE
    !-----------------------------------------------------------------------


    SUBROUTINE CHECK_CONVERGENCE&
    (ITER, MITER, FLAG_NR, CCL, FoM, RSUM_NEW, RSUM_OLD, B, Be, IDIM_B, IDIM_Be,&
           RES_NRM, LMSG, MNR_FAILED, Res_Norm_First_Iteration, time_in_seconds )
        
        use check_for_floating_point_exceptions
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
        
        
        
        LMSG = .FALSE.

        CCL = 'N'
        

        
        RSUM_NEW  = DSQRT(DOT_PRODUCT(B,B)+DOT_PRODUCT(Be,Be))
        
        ERROR_NEW = MAX(MAXVAL(DABS(B)),MAXVAL(DABS(Be)))
        
        If (ITER==1)   Res_Norm_First_Iteration = RSUM_NEW
        
        
        ! call check_fp_exceptions(RSUM_NEW, "RSUM_NEW")
        ! call check_fp_exceptions(RSUM_OLD, "RSUM_OLD")
        ! call check_fp_exceptions(RES_NRM, "RES_NRM")
        ! call check_fp_exceptions(Res_Norm_First_Iteration, "Res_Norm_First_Iteration")
        ! call check_fp_exceptions(B, "B")
        ! call check_fp_exceptions(Be, "Be")
        ! call check_fp_exceptions(ERROR_NEW, "ERROR_NEW")
        
        

        IF (RSUM_NEW.NE.RSUM_NEW) THEN
            WRITE(*,*) 'NaN ENCOUNTERED, GENERAL STOP'
            LMSG = .TRUE.
            RETURN
        ENDIF
       
       
       ! SELECT AMONG FULL AND MODIFIED NEWTON-RAPHSON
        KIND_OF_NEWTON_RAPHSON: SELECT CASE(FLAG_NR)

            CASE('NRP')

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
        FLAG_NR = 'NRP'


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
          



      

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !                   SUBROUTINE BASIS_2d
    !----------------------------------------------------------------------
          



     

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !  SUBROUTINE UPDATE_SOLUTION
    !--------------------------------------------------------------------------

    SUBROUTINE UPDATE_SOLUTION( INCREMENT )




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


end module former_external_subroutines
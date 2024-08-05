module solution_check_update

    contains

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


end module solution_check_update
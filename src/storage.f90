module storage

    contains

    SUBROUTINE MATRIX_STORAGE_RESIDUAL( TEMP, NM, IDM, INEQ, B, IDIM_B )

        IMPLICIT NONE

       ! ARGUMENTS
       INTEGER, INTENT(IN)    :: IDM, INEQ, IDIM_B
       
       INTEGER, DIMENSION(IDM),      INTENT(IN)    :: NM
       REAL(8), DIMENSION(IDM,INEQ), INTENT(IN)    :: TEMP
       REAL(8), DIMENSION(IDIM_B),   INTENT(INOUT) :: B

       ! LOCAL VARIABLES
       INTEGER :: I, IEQ, IROW


       ! STORE THE RESIDUAL VECTOR IN THE GLOBAL VECTOR B
       DO I = 1, IDM

       ! ITERATE OVER IDIM NODES
         DO IEQ = 1, INEQ

         ! ITERATE OVER INEQ EQUATIONS
           IROW  = NM(I) + IEQ - 1
           !$OMP ATOMIC
           B(IROW) = B(IROW) + TEMP(I,IEQ)

         ENDDO
       ENDDO
       
    END SUBROUTINE MATRIX_STORAGE_RESIDUAL



    SUBROUTINE MATRIX_STORAGE_JACOBIAN&
    (TP, IDM, JDIM, JNEQ, INEQ, NM, IA, IIA, CSR, ICSR, A, IDIM_A)

       IMPLICIT NONE

       ! ARGUMENTS
       INTEGER, INTENT(IN) :: ICSR, IIA
       INTEGER, INTENT(IN) :: IDM, JDIM, INEQ, JNEQ
       INTEGER, INTENT(IN) :: IDIM_A

       INTEGER, DIMENSION(ICSR), INTENT(IN) :: CSR
       INTEGER, DIMENSION(IIA),  INTENT(IN) :: IA
       INTEGER, DIMENSION(IDM),  INTENT(IN) :: NM

       REAL(8), DIMENSION(IDM,JDIM,JNEQ,INEQ), INTENT(IN)    :: TP
       REAL(8), DIMENSION(IDIM_A),             INTENT(INOUT) :: A

       ! LOCAL VARIABLES
       INTEGER  :: I, J, L, IEQ, JEQ
       INTEGER  :: IROW, ICOL, IAD


       ! ITERATE OVER IDIM NODES
       DO I = 1, IDM
         IROW = NM(I)

       ! ITERATE OVER INEQ EQUATIONS
         DO IEQ = 1, INEQ
           IAD  = (IEQ-1)*(IA(IROW+1)-IA(IROW))

         ! ITERATE OVER JDIM BASIS FUNCTIONS
           DO J = 1, JDIM
             L = (I-1)*JDIM + J

           ! ITERATE OVER JNEQ EQUATIONS
             DO JEQ = 1, JNEQ
               ICOL = IAD + CSR(L) + JEQ - 1

             ! ITERATE OVER ALL EQUATIONS
               !$OMP ATOMIC
               A(ICOL) = A(ICOL) + TP(I,J,JEQ,IEQ)
               ! write(404,'(4(i3,2x),f20.10)') J, JEQ, I, IEQ, TP(I,J,JEQ,IEQ)

             ENDDO
           ENDDO

         ENDDO
       ENDDO
           
    END SUBROUTINE MATRIX_STORAGE_JACOBIAN

end module storage
 
Module VariableMapping

    contains 

      Function variableExist( Var ) Result(output)
        Implicit None 
        Character(len=*), Intent(In) :: Var
        Logical                      :: output

        output = getVariableId(Var) /= -1
      End Function variableExist

      Function getNumberOfUnknowns() Result(output)
        Implicit None 
        Integer :: output

        output = 0
        do 
          output = output + 1
          if (getVariableName(output) == "") then
            output = output - 1
            exit
          end if
        end do
      End Function getNumberOfUnknowns



      Function getVariableId( Var ) Result (Vid)
        Use Formats
        Implicit None 
        Character(len=*), Intent(In) :: Var
        Integer                      :: Vid 

        Select Case(Var)
        Case('C' ) ; Vid = 1
        ! Case('Vz' ) ; Vid = 2
        ! Case('P'  ) ; Vid = 3
        ! Case('Z'  ) ; Vid = 4
        ! Case('R'  ) ; Vid = 5
        
        Case Default; Vid = -1
        End Select 
      End Function getVariableId

      Function getVariableName( vid) Result(var)
        Implicit None 
        Integer, Intent(In)           :: vid
        Character(len=:), Allocatable :: var


        Select Case(vid)
        Case(1)      ; Var = 'C' 
        ! Case(2)      ; Var = 'Vz' 
        ! Case(3)      ; Var = 'P'  
        ! Case(4)      ; Var = 'Z'  
        ! Case(5)      ; Var = 'R'  
        
        Case Default ; Var = ''
        End Select 
      End Function getVariableName

End Module VariableMapping



MODULE PHYSICAL_MODULE
  
    Real(8), PARAMETER       :: pi        = 3.141592653589793D0
    Real(8), parameter       :: g_grav    =  9.81d0           ! m/s2:  gravitational acceleration
    

    Real(8), parameter       :: rho       = 1000.d0           ! kg/m3: density of water
    Real(8), parameter       :: length_char   =  1.d0         ! m
    


    Real(8), parameter       :: velocity_char       =  0.01d0  ! m/s
  
    Real(8), parameter       :: viscosity = 1.0d0


    !_______________________________________________________________________________

    Real(8), parameter       :: diffusivity = 1.0d-6  ! m2/s
    !_______________________________________________________________________________
  

    Real(8), parameter       :: surface_tension = 0.073d0
    
    ! Real(8), parameter       :: time_char           = length_char/velocity_char            ! s
    Real(8), parameter       :: time_char           = length_char**2/diffusivity            ! s

    !_______________________________________________________________________________

    Real(8), parameter       :: viscous_stress      = viscosity/time_char

    Real(8), parameter       :: inertial_stress     = rho * ( velocity_char )**2

    Real(8), parameter       :: gravity_stress      = rho*g_grav*length_char

    Real(8), parameter       :: capillary_stress    = surface_tension/length_char
    !_______________________________________________________________________________


    Real(8), parameter       :: Pchar             = viscous_stress
    
    
    Real(8), parameter       :: ratio_of_pressures = gravity_stress/Pchar

    Real(8)                  :: ReN  !  REYNOLDS NUMBER
    

    Real(8)                  :: eo1, eo2
    Real(8)                  :: es1, es2
    Real(8)                  :: e_bnd
  
    Contains


      Subroutine DIMENSIONLESS_NUMBERS

        Implicit None
      
        Integer :: I
        Real(8) :: DUM

      
        ! initial_position =  ho/length_char
        ! ambient_position =  ho/length_char

        ReN   =  inertial_stress/Pchar
        ! print*, "pi  =", pi
        ! print*, " "
        ! print*, "rho  =", rho
        ! print*, " "
        ! print*, "length_char  =", length_char
        ! print*, " "
        ! print*, "g_grav  =", g_grav
        ! print*, " "
        ! print*, "viscosity  =", viscosity
        ! print*, " "
        ! print*, "surface_tension  =", surface_tension
        ! print*, " "
        ! print*, "velocity_char  =", velocity_char
        ! print*, " "
        ! print*, "time_char  =", time_char
        ! print*, " "
        ! print*, "viscous_stress  =", viscous_stress
        ! print*, " "
        ! print*, "inertial_stress  =", inertial_stress
        ! print*, " "
        ! print*, "gravity_stress  =", gravity_stress
        ! print*, " "
        ! print*, "capillary_stress  =", capillary_stress
        ! print*, " "
        ! print*, "Pchar  =", Pchar
        ! print*, " "
        ! print*, "ratio_of_pressures  =", ratio_of_pressures
        ! print*, " "
        ! print*, "ReN  =", ReN
        ! print*, " "
        ! pause
        write(*,"(A6,2X,F16.8)") "Reff ="    , length_char
        write(*,"(A6,2X,F16.8)") "ReN  ="    , ReN
        
        eo1   = 0.0D0
        eo2   = 0.1D0
        e_bnd = - 1.0D+4
        

    END SUBROUTINE DIMENSIONLESS_NUMBERS

    function calculateFlowrate(time)  result(output)
      real(8), intent(in) :: time
      real(8)             :: output

      ! output = 1.d0 - exp(-1.d0*time)
      output = 1.d0

    end function calculateFlowrate


  END MODULE PHYSICAL_MODULE
  
! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------  

Module TIME_INTEGRATION 
  Implicit None
  
  Real(8) :: TIME
  Real(8) :: timeo
  Real(8) :: timeb
  Real(8) :: timep
  Real(8) :: INITIAL_TIME
  Real(8) :: FINAL_TIME
  
  Real(8) :: Dt,DTo, Dtb
  Real(8) :: Dt_constant
  Real(8) :: Dt_min
  Real(8) :: Dt_max
  integer :: INCREMENT_STEP, INCREMENT_TIME
  LOGICAL :: Adjust_Dt




  contains

  subroutine set_DT
      Dt_constant = 0.01d0
  end subroutine set_DT
  

  subroutine update_DT
        
  end subroutine update_DT

End Module TIME_INTEGRATION


! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

Module RemeshVariables
  Implicit none
  integer :: Remesh_counter
  Logical :: StructMeshOnlyInTheFront
  integer :: Remesh_counter_structuredInTheFront
end module RemeshVariables




Module ELEMENTS_MODULE
 Integer, Parameter:: NBF_1d = 2             ! NUMBER OF BASIS FUNCTIONS PER 1d ELEMENT
 Integer, Parameter:: NBF_2d = 3             ! NUMBER OF BASIS FUNCTIONS PER 2d ELEMENT
 Integer, Parameter:: NED_2d = 3             ! NUMBER OF EDGES PER 2d ELEMENT
    
 !   NUMBER OF EQUATIONS
 Integer, Parameter:: NEQ_f = 1             ! NUMBER OF PDEs SYSTEM TO SOLVE FOR FLOW
 Integer, Parameter:: NEX_f = 0 

 !   NUMBER OF ELEMENTS
 Integer :: NEL_1d                 ! TOTAL NUMBER OF 1D ELEMENTS
 Integer :: NEL_2d                 ! TOTAL NUMBER OF 2D ELEMENTS
 Integer :: NEL                    ! TOTAL NUMBER OF 1D & 2D
 
 !   NUMBER OF NODES
 Integer :: NODTOL                 ! TOTAL NUMBER OF NODES
 
 !   NUMBER OF UNKNOWNS
 Integer :: NUNKNOWNS_f            ! TOTAL NUMBER OF UNKNOWNS
    
    
 Contains


  Subroutine Discretization_Data ( unvfile ) 
    Use unv
    Implicit None
    Type(unvFileReader), Intent(In) :: unvfile
   
    NODTOL      = unvfile%getNumberOfNodes()
    NEL_1d      = unvfile%getNumberOfSurfaceElements()
    NEL_2d      = unvfile%getNumberOfElements()

    ! TOTAL NUMBER OF UNKNOWNS
    NUNKNOWNS_f = NODTOL*NEQ_f
  End Subroutine Discretization_Data

End Module Elements_Module
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------


 MODULE OMP_PARALLEL
   USE IFPORT, only: SETENVQQ    
         
   INTEGER :: NTHREADS
   
   contains

   subroutine set_num_threads

     logical :: success
     CHARACTER(LEN=100) :: STR_ITER_TMP, THREADS, FN

     NTHREADS = 12 ! number of threads that pardiso will use and the loop of the jacobian

     WRITE(STR_ITER_TMP,'(I4)') NTHREADS

     STR_ITER_TMP = TRIM(ADJUSTL(STR_ITER_TMP))
   
     FN = 'MKL_NUM_THREADS='//STR_ITER_TMP
     
     success = SETENVQQ(FN)

   end subroutine set_num_threads

         
 END MODULE OMP_PARALLEL
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------


  MODULE CONTINUATION_MODULE
    
    INTEGER            :: INCREMENT
      
    REAL(8)            :: ArN
    REAL(8)            :: ArNo, ArNb, ArNp

    REAL(8)            :: dArN = - 0.01D0
    
    REAL(8), PARAMETER :: INITIAL_ArN = 1.0D0
    REAL(8), PARAMETER :: FINAL_ArN   = 1.D+5
    
    REAL(8)            :: dS, S0SM1, SNRM

    REAL(8)            :: dDpL, dDpLo, dDpLodS
    REAL(8)            ::       dArNo, dArNodS
    
    REAL(8), ALLOCATABLE, DIMENSION(:) :: dTL, dTLo, dTLodS
    
    
    ! ----------------------------------------------------------------------
    contains
    ! ----------------------------------------------------------------------


    Subroutine Allocate_Continuation_Arrays(L, NUNKNOWNS_f_) 
      Implicit None 
      Logical, Intent(In)           :: L
      Integer, Intent(In), Optional :: NUNKNOWNS_f_
 
      If (L) Then
        allocate(    dTL(NUNKNOWNS_f_) )          ;  dTL    = 0.0D0
        allocate(   dTLo(NUNKNOWNS_f_) )          ;  dTLo   = 0.0D0
        allocate( dTLodS(NUNKNOWNS_f_) )          ;  dTLodS = 0.0D0
      Else 
        deallocate( dTL    )
        deallocate( dTLo   )
        deallocate( dTLodS )
      End If 
    End Subroutine Allocate_Continuation_Arrays
        
  END MODULE CONTINUATION_MODULE
  

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------


  MODULE FLOW_ARRAYS_MODULE

 
    REAL(8), ALLOCATABLE, DIMENSION(:)            :: B_f, Be_f, Bi_f, Bw_f
    REAL(8), ALLOCATABLE, DIMENSION(:)            :: S_f, Se_f, Sb_f
    REAL(8), ALLOCATABLE, DIMENSION(:,:)          :: Sa_f
    LOGICAL, ALLOCATABLE, DIMENSION(:)            :: Bc_f, Bce_f
    
    
    ! ----------------------------------------------------------------------
    contains
    ! ----------------------------------------------------------------------


    Subroutine Allocate_Flow_Arrays(L, NUNKNOWNS_f_, NEX_f_) 
      Implicit None 
      Logical, Intent(In)           :: L
      Integer, Intent(In), Optional :: NUNKNOWNS_f_,NEX_f_
 
      If (L) Then
        allocate( B_f(NUNKNOWNS_f_) )          ;  B_f   = 0.0D0
        allocate( Be_f(NEX_f_) )               ;  Be_f  = 0.0D0
        allocate( Bi_f(NEX_f_) )               ;  Bi_f  = 0.0D0
        allocate( Bw_f(NUNKNOWNS_f_) )         ;  Bw_f  = 0.0D0
        allocate( S_f(NUNKNOWNS_f_)  )         ;  S_f   = 0.0D0
        allocate( Se_f(NEX_f_)  )              ;  Se_f  = 0.0D0
        allocate( Sa_f(NUNKNOWNS_f_,NEX_f_)  ) ;  Sa_f  = 0.0D0
        allocate( Sb_f(NUNKNOWNS_f_)  )        ;  Sb_f  = 0.0D0
        allocate( Bc_f(NUNKNOWNS_f_)  )        ;  Bc_f  = .FALSE.
        allocate( Bce_f(NUNKNOWNS_f_) )        ;  Bce_f = .FALSE.
      Else 
        deallocate( B_f   ) 
        deallocate( Be_f  ) 
        deallocate( Bi_f  )
        deallocate( Bw_f  )
        deallocate( S_f   )
        deallocate( Se_f  )
        deallocate( Sa_f  )
        deallocate( Sb_f  )
        deallocate( Bc_f  )
        deallocate( Bce_f )
      ENDIF
 
    End Subroutine Allocate_Flow_Arrays
    
  END MODULE FLOW_ARRAYS_MODULE
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------


  MODULE GLOBAL_ARRAYS_MODULE

 
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: TL, TLo, TLb, TLp    
    REAL(8)                              :: DpL, DpLo, DpLb, DpLp
    
    
    ! ----------------------------------------------------------------------
    contains
    ! ----------------------------------------------------------------------


    Subroutine Allocate_Global_Arrays(L, NODTOL_,NEL_2d_,NEQ_f_) 
      Implicit None 
      Logical, Intent(In)           :: L
      Integer, Intent(In), Optional :: NODTOL_, NEQ_f_, NEL_2d_
 
      If (L) Then
        allocate( TL (NODTOL_,NEQ_f_) )              ;  TL  = 0.D0
        allocate( TLo(NODTOL_,NEQ_f_) )              ;  TLo = 0.D0
        allocate( TLb(NODTOL_,NEQ_f_) )              ;  TLb = 0.D0
        allocate( TLp(NODTOL_,NEQ_f_) )              ;  TLp = 0.D0

      Else 
        deallocate( TL  )
        deallocate( TLo )
        deallocate( TLb )
        deallocate( TLp )
      End if 
    End Subroutine Allocate_Global_Arrays
    
  END MODULE GLOBAL_ARRAYS_MODULE
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------


  MODULE NRAPSHON_MODULE

    INTEGER, PARAMETER :: NITER     = 1000
    REAL(8), PARAMETER :: ERROR_NR  = 5.d-9
    
    
    INTEGER            :: ITER_f
    REAL(8)            :: RES_NORM, COR_NORM_OLD_f, COR_NORM_NEW_f
    REAL(8)            :: ERR_nr
    CHARACTER(LEN=3)   :: FLAG_NR = 'NRP'

    contains
    REAL(8) FUNCTION F_DX ( X )
        IMPLICIT NONE

        REAL(8), INTENT(IN)  :: X
        REAL(8) :: EP_RES = 1.0D-9


        F_DX = EP_RES*DMAX1(1.0D0,DABS(X))*SIGN(1.0D0,X)

    END FUNCTION F_DX


  END MODULE NRAPSHON_MODULE
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------



MODULE ENUMERATION_MODULE
 
    ! NODES ARE ENUMERATED ALONG THE x- DIRECTION 
    INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: NM_MESH       ! FOR TRIANGULAR ELEMENTS
    INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: NM_f          ! FOR TRIANGULAR ELEMENTS
    
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: NME_MESH      ! FOR TRIANGULAR ELEMENTS
    
    ! SURROUNDING VARIABLES & ARRAYS
    INTEGER                            :: IDEsuN, IDNsuN, maxIDNsuN
    INTEGER, ALLOCATABLE, DIMENSION(:) :: EsuN1, NsuN1
    INTEGER, ALLOCATABLE, DIMENSION(:) :: EsuN2, NsuN2
    
    ! ----------------------------------------------------------------------
    contains
    ! ----------------------------------------------------------------------


    Subroutine Allocate_Enumeration_Indices(L, NEL_2d_, NBF_1d_, NBF_2d_, NED_2d_) 
      Implicit None 
      Logical, Intent(In)           :: L
      Integer, Intent(In), Optional :: NEL_2d_,NBF_1d_,NBF_2d_,NED_2d_
 
      If (L) Then
        Allocate(   NM_MESH(NEL_2d_, NBF_2d_)          )  ;  NM_MESH    = 0
        Allocate(      NM_f(NEL_2d_, NBF_2d_)          )  ;  NM_f       = 0
        Allocate(  NME_MESH(NEL_2d_, NED_2d_, NBF_1d_) )  ;  NME_MESH   = 0
      Else 
        Deallocate( NM_MESH   )
        Deallocate( NM_f      )
        Deallocate( NME_MESH  )
      End If 
    End Subroutine Allocate_Enumeration_Indices

 
    !----------------------------------------------------------------------- 
 

   Subroutine NNUM_ELEMENT( unvfile )
     use unv
     Implicit None
     Type(unvFileReader), Intent(In) :: unvfile
     
     NM_MESH = unvfile%getElements() 
   End Subroutine NNUM_ELEMENT

 
    !-----------------------------------------------------------------------

  
   INTEGER FUNCTION GNTR(GNOD)
   
     USE ELEMENTS_MODULE, only: NEQ_f

     IMPLICIT NONE
     
     ! ARGUMENTS
     INTEGER, INTENT(IN) :: GNOD

     ! GIVEN A GLOBAL NODE RETURNS THE ROW OF
     ! THE FIRST EQUATION IN THE RESIDUAL VECTOR
                  ! ------------
     GNTR = (GNOD-1)*NEQ_f + 1

   END FUNCTION GNTR

 
    !-----------------------------------------------------------------------

  
   SUBROUTINE NNUM_f
 
     USE ELEMENTS_MODULE, only: NEL_2d, NBF_2d

     IMPLICIT NONE
 
    ! LOCAL VARIABLES
     INTEGER :: IEL, INOD, II
     
     
     DO IEL = 1, NEL_2d
       DO INOD = 1, NBF_2d
         II = NM_MESH(IEL,INOD)
         
         NM_f(IEL,INOD) = GNTR(II)
         
       ENDDO
     ENDDO
     
     
   END SUBROUTINE NNUM_f

 
    !----------------------------------------------------------------------- 



   SUBROUTINE SURROUNDING_NUMBERING


     USE ELEMENTS_MODULE, only: NBF_2d, NODTOL, NEL_2d

       
     IMPLICIT NONE

     INTEGER :: IEL, INOD, JNOD, KNOD
     INTEGER :: I, J, K, L, II, JJ, KK, LL
     INTEGER :: IM
     INTEGER, DIMENSION(NBF_2d)         :: ISID, JSID 
     INTEGER, ALLOCATABLE, DIMENSION(:) :: LPN
     INTEGER, ALLOCATABLE, DIMENSION(:) :: ITEMP

    !-------------------------------------------------------------
    !-------------------------------------------------------------
    
    !-------------------------------------------------------------
    !    ELEMENTS SURROUNDING NODES
    !-------------------------------------------------------------
     if ( allocated(EsuN1) ) Deallocate(EsuN1)
     if ( allocated(EsuN2) ) Deallocate(EsuN2)
     if ( allocated(LPN)   ) Deallocate(LPN)
     if ( allocated(ITEMP) ) Deallocate(ITEMP)
     if ( allocated(NsuN1) ) Deallocate(NsuN1)
     if ( allocated(NsuN2) ) Deallocate(NsuN2)
     ! EsuN1 --> ARRAY HOLDING ELEMENTS THAT SURROUND EACH NODE
     ! EsuN2 --> INDEXING ARRAY


     ! ALLOCATE ARRAY EsuN2 & LPN
     ALLOCATE( EsuN2(NODTOL+1) )
     ALLOCATE( LPN(NODTOL) )
     ALLOCATE( ITEMP(NODTOL) )


     ! COUNT NUMBER OF ELEMENTS CONNECTED TO EACH NODE
     EsuN2 = 0

     DO IEL  = 1, NEL_2d
       DO INOD = 1, NBF_2d
         II = NM_MESH(IEL,INOD)
         EsuN2(II) = EsuN2(II) + 1
       ENDDO
     ENDDO

   
     ! RESHUFFLING
     DO II = 2, NODTOL+1
       EsuN2(II) = EsuN2(II) + EsuN2(II-1)
     ENDDO


     ! DEFINE IDEsuN
     IDEsuN = EsuN2(NODTOL+1)


     ! RESHUFFLING
     DO II = NODTOL+1, 2, -1
       EsuN2(II) = EsuN2(II-1)
     ENDDO
     EsuN2(1) = 0


     ! ALLOCATE ARRAY EsuN1
     ALLOCATE( EsuN1(IDEsuN) )

     EsuN1 = 0

     ! STORE THE ELEMENTS IN ARRAY EsuN1
     DO IEL  = 1, NEL_2d
       DO INOD = 1, NBF_2d
         II = NM_MESH(IEL,INOD)
         JJ = EsuN2(II)
         LOOP_IN: DO 
           JJ = JJ + 1
           IF(EsuN1(JJ)==0)THEN
             EsuN1(JJ) = IEL
             EXIT LOOP_IN
           ENDIF 
         ENDDO LOOP_IN
       ENDDO
     ENDDO


    !-------------------------------------------------------------
    !    NODES SURROUNDING NODES
    !-------------------------------------------------------------


     ! NsuN1 --> ARRAY HOLDING NODES THAT SURROUND EACH NODE
     ! NsuN2 --> INDEXING ARRAY


     ! DEFINE maxIDNsuN
     maxIDNsuN = (NBF_2d-1)*IDEsuN


     ! ALLOCATE ARRAY NsuN2
     ALLOCATE( NsuN2(NODTOL+1) )
     ALLOCATE( NsuN1(maxIDNsuN) )


     ! INITIALIZATION
     NsuN1 = 0
     NsuN2 = 0


     ! ASSEMBLE EsuN, NsuN
     JJ = 0

     LOOP_1: DO II = 1, NODTOL

       LPN = 0

       LOOP_2: DO KK = EsuN2(II)+1, EsuN2(II+1)
         IEL  = EsuN1(KK)

         LOOP_3: DO INOD = 1, NBF_2d
           JNOD = NM_MESH(IEL,INOD)
           KNOD = JNOD
           IF(LPN(JNOD)==0)THEN
             JJ = JJ + 1
             NsuN1(JJ) = KNOD 
             LPN(JNOD)   = INOD
           ENDIF
         ENDDO LOOP_3

       ENDDO LOOP_2

       NsuN2(II+1) = JJ
     ENDDO LOOP_1


     ! DEFINE IDNsuN
     IDNsuN = NsuN2(NODTOL+1)

     
     ! RESHUFFLING
     LOOPNODES: DO I = 1, NODTOL

       L = 0
       DO J = NsuN2(I)+1, NsuN2(I+1)
         L = L + 1
         ITEMP(L) = NsuN1(J)
       ENDDO

       LL = 0
       DO J = NsuN2(I)+1, NsuN2(I+1)
         LL = LL + 1

         JNOD = 0
         JNOD = MINVAL(ITEMP(1:L), MASK= ITEMP .GT. 0)

         DO KK = 1, L
           IF(JNOD == ITEMP(KK))THEN
             ITEMP(KK) = 0
           ENDIF
         ENDDO

         NsuN1(J) = JNOD
       ENDDO

     ENDDO LOOPNODES


    END SUBROUTINE SURROUNDING_NUMBERING

END MODULE ENUMERATION_MODULE
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------



 MODULE CSR_STORAGE
    

   INTEGER                              :: NZ_f, N_f, NRHS_f, MAXFCT_f, MNUM_f,&
                                           MTYPE_f, PHASE_f, MSGLVL_f, IDUM_f
   INTEGER                              :: ERROR_f = 0
   REAL(8)                              :: DDUM_f
   
   REAL(8), ALLOCATABLE, DIMENSION(:)   :: A_f         ! LINEAR SYSTEM'S MATRIX A
   REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Ah_f, Ar_f, Ac_f, Ai_f
   
   INTEGER, ALLOCATABLE, DIMENSION(:)   :: IA_f        ! INDEXING OF MATRIX A 
   INTEGER, ALLOCATABLE, DIMENSION(:)   :: CA_f        ! COLUMNS OF MATRIX A
   INTEGER, ALLOCATABLE, DIMENSION(:)   :: RA_f        ! ROWS OF MATRIX A
   INTEGER, ALLOCATABLE, DIMENSION(:)   :: DA_f        ! DIAGONAL OF MATRIX A
   INTEGER,              DIMENSION(64)  :: IPARM_f = 0
   INTEGER,              DIMENSION(64)  :: PT_f    = 0
   
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CSR_f
    
   ! ----------------------------------------------------------------------
   contains
   ! ----------------------------------------------------------------------


    Subroutine Allocate_CSR_Arrays(L, K,NUNKNOWNS_f_,NZ_f_,NEX_f_)
 
      Implicit None 
      Logical, Intent(In)           :: L
      Integer, Intent(In), Optional :: K, NUNKNOWNS_f_, NZ_f_,NEX_f_
 
      If (L) Then
        If (K.EQ.1) Then
          Allocate( IA_f(NUNKNOWNS_f_+1) )      ;  IA_f = 0
        Else
          Allocate(  A_f(NZ_f_) )               ;  A_f  = 0
          Allocate( Ah_f(NEX_f_,NEX_f_) )       ;  Ah_f = 0
          Allocate( Ar_f(NEX_f_,NUNKNOWNS_f_) ) ;  Ar_f = 0
          Allocate( Ac_f(NUNKNOWNS_f_,NEX_f_) ) ;  Ac_f = 0
          Allocate( Ai_f(NEX_f_,NEX_f_) )       ;  Ai_f = 0
          Allocate( CA_f(NZ_f_) )               ;  CA_f = 0
          Allocate( RA_f(NZ_f_) )               ;  RA_f = 0
          Allocate( DA_f(NUNKNOWNS_f_) )        ;  DA_f = 0
        End If
      Else 
        deallocate( IA_f )
        deallocate( CA_f )
        deallocate( A_f  )
        deallocate( Ah_f )
        deallocate( Ar_f )
        deallocate( Ac_f )
        deallocate( Ai_f )
        deallocate( RA_f )
        deallocate( DA_f )
      End If
 
    End Subroutine Allocate_CSR_Arrays

 
    !-----------------------------------------------------------------------


    Subroutine Allocate_CSR_Connectivity(L, NBF_2d_, NEL_2d_) 
      Implicit None 
      Logical, Intent(In)           :: L
      Integer, Intent(In), Optional :: NBF_2d_, NEL_2d_
 
      If (L) Then
        allocate( CSR_f(NEL_2d_,NBF_2d_*NBF_2d_) )  ;  CSR_f = 0
      Else 
        deallocate( CSR_f )
      End If 
    End Subroutine Allocate_CSR_Connectivity

 
    !-----------------------------------------------------------------------


   SUBROUTINE CSR_MAIN(ICHOICE)

     USE ELEMENTS_MODULE,    only: NEQ_f, NBF_2d, NODTOL, NEL_2d, NUNKNOWNS_f
     USE ENUMERATION_MODULE, only: NM_MESH, NM_f

     IMPLICIT NONE

     ! ARGUMENTS
     CHARACTER(LEN=2), INTENT(IN) :: ICHOICE
     
     ! LOCAL VARIABLES
     INTEGER :: IEL, I, J


     ! NUMBERING: NU ; CONNECTIVITY: CO
     SELECT CASE(ICHOICE)

       ! DEFINE NON-ZEROS, IA & JA
       CASE('NU')
       CALL CSR_NUMBERING

       ! CSR CONNECTIVITY FOR ELEMENTAL MATRIX STORAGE
       CASE('CO')
       CALL ALLOCATE_CSR_CONNECTIVITY(.TRUE.,NBF_2d,NEL_2d)

       ! DEFINE CSR FLOW CONNECTIVITY FOR ALL ELEMENTS
       LOOP_ALL_ELEMENTS_CSR_f: DO IEL = 1, NEL_2d
         
         CALL CSR_CONNECTIVITY&
         (IEL, NEL_2d, NM_f, NBF_2d, CSR_f, NBF_2d*NBF_2d, NZ_f, IA_f, &
          CA_f, RA_f, NUNKNOWNS_f+1)

         CALL CSR_DIAGONAL&
         (IEL, NEL_2d, NM_f, NBF_2d, NEQ_f, DA_f, NUNKNOWNS_f, RA_f,   &
          CA_f, NZ_f, IA_f, NUNKNOWNS_f+1)
          
       ENDDO LOOP_ALL_ELEMENTS_CSR_f

     END SELECT
     

   END SUBROUTINE CSR_MAIN
   
   
    !-----------------------------------------------------------------------


   SUBROUTINE CSR_NUMBERING

     USE ELEMENTS_MODULE,    only: NODTOL, NUNKNOWNS_f, NBF_2d, NEQ_f, NEX_f
     USE ENUMERATION_MODULE, only: IDNsuN, NsuN1, NsuN2, IDNsuN, NsuN1, NsuN2, GNTR

     IMPLICIT NONE

     INTEGER :: I, II, J, JJ, K, KK, L, MZ
     

     ! INITIALIZE NZ
     NZ_f = 0
     
     
     ! ALLOCATE IA ARRAY
     CALL ALLOCATE_CSR_ARRAYS(.TRUE., 1, NUNKNOWNS_f, NZ_f, NEX_f)


     ! ASSIGN VALUES TO IA ARRAY & NZ PARAMETER
     IA_f(1) = 1  ;  JJ = 0

     DO I = 1, NODTOL
       DO II = 1, NEQ_f

         JJ = JJ + 1
         MZ = 0

         DO J = NsuN2(I)+1, NsuN2(I+1)
           NZ_f = NZ_f + NEQ_f
           MZ   = MZ   + NEQ_f
         ENDDO

         IA_f(JJ+1) = IA_f(JJ) + MZ

       ENDDO
     ENDDO


     ! ALLOCATE JA, AA & RA
     CALL ALLOCATE_CSR_ARRAYS(.TRUE., 2, NUNKNOWNS_f, NZ_f, NEX_f)


     ! ASSIGN VALUES TO RA ARRAY
     DO I = 1, NUNKNOWNS_f
       DO J = IA_f(I), IA_f(I+1)-1
         RA_f(J) = I
       ENDDO
     ENDDO

     
     ! ASSIGN VALUES TO JA ARRAY
     KK = 0

     DO I = 1, NODTOL
       DO II = 1, NEQ_f

         DO J = NsuN2(I)+1, NsuN2(I+1)
           L = NsuN1(J)
           DO JJ = 1, NEQ_f
             KK = KK + 1
             CA_f(KK) = GNTR(L) + JJ - 1
           ENDDO
         ENDDO
       ENDDO
     ENDDO
     

   END SUBROUTINE CSR_NUMBERING
 
 
    !-----------------------------------------------------------------------


   SUBROUTINE CSR_CONNECTIVITY&
   (IEL_, NEL_2d_, NM_, NBF_2d_, CSR_, EJS_, NZ_, IA_, JA_, RA_, NUNP1_)
     
     IMPLICIT NONE

     ! ARGUMENTS
     INTEGER, INTENT(IN) :: IEL_, NEL_2d_, NBF_2d_, EJS_, NZ_, NUNP1_
     INTEGER, DIMENSION(NEL_2d_,NBF_2d_), INTENT(IN)  :: NM_
     INTEGER, DIMENSION(NZ_),            INTENT(IN)  :: RA_
     INTEGER, DIMENSION(NZ_),            INTENT(IN)  :: JA_
     INTEGER, DIMENSION(NUNP1_),         INTENT(IN)  :: IA_

     INTEGER, DIMENSION(NEL_2d_,EJS_),    INTENT(OUT) :: CSR_

     ! LOCAL VARIABLES
     INTEGER :: I, J, K, L, II, JJ
     

     L = 0

     DO I = 1, NBF_2d_
       II = NM_(IEL_,I)
       DO J = 1, NBF_2d_
         JJ = NM_(IEL_,J)

         L = L + 1

         DO K = IA_(II), IA_(II+1)-1
           IF(RA_(K)==II .AND. JA_(K)==JJ)THEN
             CSR_(IEL_,L) = K
           ENDIF
         ENDDO

       ENDDO
     ENDDO
     

   END SUBROUTINE CSR_CONNECTIVITY
   
   
    !-----------------------------------------------------------------------


   SUBROUTINE CSR_DIAGONAL&
   (IEL_, NEL_2d_, NM_f_, NBF_2d_, NEQ_f, DA_, NUNKNOWNS_, RA_, JA_, NZ_, IA_, NUNP1_)


     IMPLICIT NONE

    ! ARGUMENTS
     INTEGER, INTENT(IN) :: IEL_, NEL_2d_, NBF_2d_, NUNKNOWNS_, NZ_, NUNP1_
     INTEGER, INTENT(IN) :: NEQ_f

     INTEGER, DIMENSION(NZ_),            INTENT(IN)  :: RA_
     INTEGER, DIMENSION(NZ_),            INTENT(IN)  :: JA_
     INTEGER, DIMENSION(NUNP1_),         INTENT(IN)  :: IA_
     INTEGER, DIMENSION(NUNKNOWNS_),     INTENT(OUT) :: DA_

     INTEGER, DIMENSION(NEL_2d_,NBF_2d_), INTENT(IN)  :: NM_f_

     ! LOCAL VARIABLES
     INTEGER :: I, J, K, II, JJ, IEQ



     DO I = 1, NBF_2d_
       JJ = NM_f_(IEL_,I)

       DO IEQ = 0, NEQ_f-1
         II = JJ + IEQ

         DO K = IA_(II), IA_(II+1)-1
           IF(RA_(K)==II .AND. JA_(K)==II)THEN
             DA_(II) = K
           ENDIF
         ENDDO

       ENDDO

     ENDDO
     

   END SUBROUTINE CSR_DIAGONAL
 
 
    !-----------------------------------------------------------------------


   SUBROUTINE INITIALIZE_PARDISO
     
     USE ELEMENTS_MODULE, only: NUNKNOWNS_f
     USE OMP_PARALLEL,    only: NTHREADS
     
     IMPLICIT NONE
     
     
     ! INITIALIZE PARDISO VARIABLES

     ! FOR FLOW PROBLEM
     N_f =  NUNKNOWNS_f
     NRHS_f     = 1
     MAXFCT_f   = 1
     MNUM_f     = 1
     IPARM_f(3) = NTHREADS   ! SET OMP_NUM_THREADS PARAMETER
     IPARM_f(1) = 0
     


      ! SETUP PARDISO CONTROL PARAMETERS AND INITIALIZE THE SOLVERS
  
      ! FOR FLOW PROBLEM
      MTYPE_f = 1              ! REAL & STRUCTURALLY SYMMETRIC MATRIX
  
      
  
      ! REORDERING AND SYMBOLIC FACTORIZATION, THIS STEP ALSO ALLOCATES
      ! ALL MEMORY THAT IS NECESSARY FOR THE FACTORIZATION
  
      ! FOR FLOW PROBLEM
  
      PHASE_f  = 11            ! ONLY REORDERING AND SYMBOLIC FACTORIZATION
      MSGLVL_f = 0             ! WITHOUT STATISTICAL INFORMATION
  
  
      CALL PARDISO&
      (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
      IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, DDUM_f, DDUM_f, ERROR_f)
  
  
      ! PRINT*, 'REORDERING OF A_f COMPLETED ... '
  
  
      IF(ERROR_f .NE. 0)THEN
        PRINT*, 'THE FOLLOWING ERROR_f WAS DETECTED: ', ERROR_f
        PRINT*, 'DURING SYMBOLIC FACTORIZATION. SUB: PARDISO.'
        PRINT*, 'PROGRAM STOP.'
        PRINT*, ''
        STOP
      END IF
  
  
     ! PRINT*, 'NUMBER OF NONZEROS IN FACTORS   = ',IPARM_f(18)
     ! PRINT*, 'NUMBER OF FACTORIZATION MFLOPS  = ',IPARM_f(19)
      
      
   END SUBROUTINE INITIALIZE_PARDISO

 END MODULE CSR_STORAGE
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------



 MODULE MESH_MODULE
    

   ! GLOBAL MESH PARAMETERS
   REAL(8), ALLOCATABLE, DIMENSION(:) :: Xm, Ym
   
   INTEGER, PARAMETER :: NBD = 8  ! NUMBER OF DEFINED BOUNDARIES
   
   REAL(8), PARAMETER :: EPS_MESH = 1.D-4  ! MESH TOLERANCE
    
    
    ! ----------------------------------------------------------------------
    contains
    ! ----------------------------------------------------------------------


   Subroutine Allocate_Mesh_Arrays (L,NODTOL_) 
     Implicit None
     Logical, Intent(In)           :: L
     Integer, Intent(In), Optional :: NODTOL_
 
     If (L) Then
       Allocate( Xm(NODTOL_) )  ;  Xm = 0.D0
       Allocate( Ym(NODTOL_) )  ;  Ym = 0.D0
     Else 
       Deallocate( Xm )
       Deallocate( Ym )
     End If 
   End Subroutine Allocate_Mesh_Arrays

 
    !----------------------------------------------------------------------- 



   Subroutine MESH ( unvfile )  
     Use unv
     Implicit None
     Type(unvFileReader), Intent(In)      :: unvfile
     Real(8), Dimension(:,:), Allocatable :: nodes


     nodes = unvfile%getNodes()
     Xm    = nodes(:,1)
     Ym    = nodes(:,2)
     Deallocate(nodes)
   End Subroutine MESH

 
    !----------------------------------------------------------------------- 


 END MODULE MESH_MODULE
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------



MODULE BOUNDARY_ENUMERATION_MODULE
   Use VariableMapping 
   USE PHYSICAL_MODULE
   USE ELEMENTS_MODULE,    only: NUNKNOWNS_f, NEL_2d, NODTOL, NBF_2d,&
                                   NED_2d
   USE ENUMERATION_MODULE, only: NM_MESH, NME_MESH
   USE MESH_MODULE,        only: Xm, Ym, NBD, EPS_MESH
   Use RemeshVariables

  

   Integer, Dimension(:), Allocatable :: bnd1_elements, bnd1_faces
   Integer, Dimension(:), Allocatable :: bnd2_elements, bnd2_faces
   Integer, Dimension(:), Allocatable :: bnd3_elements, bnd3_faces
   Integer, Dimension(:), Allocatable :: bnd4_elements, bnd4_faces
   Integer, Dimension(:), Allocatable :: bnd5_elements, bnd5_faces
   Integer, Dimension(:), Allocatable :: bnd6_elements, bnd6_faces
   Integer, Dimension(:), Allocatable :: bnd7_elements, bnd7_faces
   Integer, Dimension(:), Allocatable :: bnd8_elements, bnd8_faces

   ! BOUNDARY INDICES
   TYPE INLO
     LOGICAL, DIMENSION(NBD) :: KOE  !  KOE(IBND) = FALSE, IF ELEMENT DOES NOT BELONG TO BOUNDARY # IBND
                                     !  KOE(IBND) = TRUE,  IF ELEMENT BELONGS TO BOUNDARY # IBND
                                         
     INTEGER, DIMENSION(NBD) :: BEE  !  BEE(IBND) = 0,   DEFAULT, THE ELEMENT DOES NOT BELONG TO A BOUNDARY
                                     !  BEE(IBND) = IED, 1<=IED<=NED_2d, EDGE # IED IS ATTACHED TO BOUNDARY # IBND
   END TYPE INLO
   
   TYPE(INLO), ALLOCATABLE, DIMENSION(:) :: NBE
    
    
    ! ----------------------------------------------------------------------
    contains
    ! ----------------------------------------------------------------------






    Subroutine Allocate_Boundary_Enumeration_Indices(L,NEL_2d_) 
      Implicit None 
      !     ARGUMENTS
      Logical, Intent(In)           :: L
      Integer, Intent(In), Optional :: NEL_2d_ 
      !     LOCAL VARIABLES
      Integer             :: IEL
 
      If (L) Then
        Allocate( NBE(NEL_2d_) )
        do IEL = 1, NEL_2d
          NBE(IEL)%KOE = .FALSE.
          NBE(IEL)%BEE = 0
        End Do
      Else 
        Deallocate( NBE )
        if (allocated(bnd1_elements)) deallocate(bnd1_elements)
        if (allocated(bnd1_faces   )) deallocate(bnd1_faces   )

        if (allocated(bnd2_elements)) deallocate(bnd2_elements)
        if (allocated(bnd2_faces   )) deallocate(bnd2_faces   )

        if (allocated(bnd3_elements)) deallocate(bnd3_elements)
        if (allocated(bnd3_faces   )) deallocate(bnd3_faces   )

        if (allocated(bnd4_elements)) deallocate(bnd4_elements)
        if (allocated(bnd4_faces   )) deallocate(bnd4_faces   )
        
        if (allocated(bnd5_elements)) deallocate(bnd5_elements)
        if (allocated(bnd5_faces   )) deallocate(bnd5_faces   )

        if (allocated(bnd6_elements)) deallocate(bnd6_elements)
        if (allocated(bnd6_faces   )) deallocate(bnd6_faces   )

        if (allocated(bnd7_elements)) deallocate(bnd7_elements)
        if (allocated(bnd7_faces   )) deallocate(bnd7_faces   )

        if (allocated(bnd8_elements)) deallocate(bnd8_elements)
        if (allocated(bnd8_faces   )) deallocate(bnd8_faces   )
      end if
 
    End Subroutine Allocate_Boundary_Enumeration_Indices

 
    !----------------------------------------------------------------------- 
   Subroutine update_NBE( ibnd, bnd_elements, faces)
      Implicit None
      Integer, Intent(In)                :: ibnd
      Integer, Dimension(:), Allocatable :: bnd_elements 
      Integer, Dimension(:), Allocatable :: faces

      Integer                            :: iel
      Integer                            :: element2d
      ! Store which of the elements have face in the boundaries

      do iel = 1, size(bnd_elements)
        element2d = bnd_elements(iel)

        nbe(element2d)%koe(ibnd) = .true. 
        nbe(element2d)%bee(ibnd) = faces(iel)
     end do
   End Subroutine Update_NBE

   Function getBoundaryNodesOfTheElement(bnd_element, face) Result(local_nodes)
      Implicit None 
      Integer, Intent(In)   :: bnd_element 
      Integer, Intent(In)   :: face 

      Integer, Dimension(2) :: local_nodes 

      if      (face == 1) Then ; local_nodes = NM_MESH(bnd_element,(/1,2/) )
      else if (face == 2) Then ; local_nodes = NM_MESH(bnd_element,(/2,3/) )
      else if (face == 3) Then ; local_nodes = NM_MESH(bnd_element,(/3,1/) )
      end  if 

   End Function getBoundaryNodesOfTheElement


  Subroutine WriteBoundaryNodesAt(unit, elements, faces)
    Implicit None 
    Integer, Intent(In)                :: unit 
    Integer, Dimension(:), Allocatable :: elements 
    Integer, Dimension(:), Allocatable :: faces 

    Integer                            :: iel 
    Integer                            :: nodes(2)


    do iel = 1, size(elements)
      nodes = getBoundaryNodesOfTheElement(elements(iel), faces(iel))
      
      Write(unit,*) Xm(nodes(1)), Ym(nodes(1))
      Write(unit,*) Xm(nodes(2)), Ym(nodes(2))
      Write(unit,*) ''
    end do 
  End Subroutine WriteBoundaryNodesAt

  Subroutine commitBoundary(unvfile, name, ibnd, elements, faces)
    Use formats
    Use unv
    Implicit None 
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Arguments
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Type(unvFileReader)               , Intent(In) :: unvfile
    Character(len=*)                  , Intent(In) :: name
    Integer                           , Intent(In) :: ibnd
    Integer, Dimension(:), Allocatable, Intent(In) :: elements
    Integer, Dimension(:), Allocatable, Intent(In) :: faces
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Local Variables 
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Print*,replace("Read Boundary * : *", "*", [toStr(ibnd), name])
    call unvfile%getBoundary(name, elements, faces) 
    call Update_NBE( unvfile%getBoundaryId(name), elements, faces)
  End Subroutine commitBoundary




  Subroutine DEFINE_BOUNDARY_NUMBERING ( unvfile )
    Use unv
    Implicit None
    Type(unvFileReader), Intent(In) :: unvfile  



    call commitBoundary(unvfile, 'BottomWall'        , 1, bnd1_elements, bnd1_faces)
    call commitBoundary(unvfile, 'TopWall'     , 2, bnd2_elements, bnd2_faces)
    call commitBoundary(unvfile, 'LeftWall'     , 3, bnd3_elements, bnd3_faces)
    call commitBoundary(unvfile, 'RightWall'     , 4, bnd4_elements, bnd4_faces)
    

     ! call WriteBoundaryNodesAt(101, bnd1_elements, bnd1_faces)
     ! call WriteBoundaryNodesAt(102, bnd2_elements, bnd2_faces)
     ! call WriteBoundaryNodesAt(103, bnd3_elements, bnd3_faces)
     ! call WriteBoundaryNodesAt(104, bnd4_elements, bnd4_faces)
     
  End Subroutine DEFINE_BOUNDARY_NUMBERING


    Subroutine getBoundaryNodes(Solution, bndElements, bndFaces, variable, output)
      Implicit None 
      Real(8), Dimension(:,:), Intent(In) :: Solution
      Integer, Dimension(:)  , Intent(In) :: bndElements
      Integer, Dimension(:)  , Intent(In) :: bndFaces
      Character(len=*)       , Intent(In) :: variable
      Real(8), Dimension(:)  , Allocatable, Intent(InOut):: output

      Integer                             :: ielem
      Integer                             :: bndTol
      Integer                             :: nbf_1d
      Integer                             :: element
      Integer                             :: face
      Integer, Dimension(:)  , Allocatable:: node
      Integer                             :: counter

      nbf_1d = 2
      bndTol = size(bndElements)

      
      if (allocated(output) ) deallocate(output)
      allocate( node   (         nbf_1d) )
      allocate( output (bndTol * nbf_1d) )

      counter = 0
      do ielem = 1, bndTol 
        element = bndElements(ielem)
        face    = bndFaces   (ielem)

        node    = getBoundaryNodesOfTheElement( element, face )
        
        counter = counter + 1
        output(counter) = Solution(node(1), getVariableId(variable) )
        counter = counter + 1
        output(counter) = Solution(node(2), getVariableId(variable) )

      end do 
      
      deallocate(node)
    End Subroutine getBoundaryNodes


  !----------------------------------------------------------------------- 

      Subroutine getVariableFromBoundaryNodes( Solution, bnd_elements, faces, variable, ValuesOnNodes )
        Use geometry, only: distance
        Implicit None
        Real(8), Dimension(:,:), Intent(In) :: Solution
        Integer, Dimension(:)  , Intent(In) :: bnd_elements
        Integer, Dimension(:)  , Intent(In) :: faces
        Character(len=*)       , Intent(In) :: variable
        Real(8), Dimension(:)  , Allocatable, Intent(Out):: ValuesOnNodes
        ! Integer, Dimension(:)  , Allocatable, Intent(Out):: GlobalNodes

        integer :: counter_i, counter_j, NodeMaxZlocation, NodeMinDist
        Integer, Dimension(:), Allocatable   :: nodesOnBoundary, nodesOnBoundary_temp
        Real(8), Dimension(:), Allocatable   :: Z_coordinates_temp,R_coordinates_temp, Dist
        Real(8), dimension(2) :: point_ref, point_next

        call getBoundaryNodesOfWholeBoundary( size(bnd_elements), bnd_elements, faces, nodesOnBoundary_temp )

        allocate( nodesOnBoundary   (size(nodesOnBoundary_temp)) )
        allocate( ValuesOnNodes     (size(nodesOnBoundary_temp)) )
        allocate( Z_coordinates_temp(size(nodesOnBoundary_temp)) )
        allocate( R_coordinates_temp(size(nodesOnBoundary_temp)) )
        allocate( Dist              (size(nodesOnBoundary_temp)) )

        do counter_i =1, size(Z_coordinates_temp)
          Z_coordinates_temp(counter_i) = Solution( nodesOnBoundary_temp(counter_i), getVariableId("Z") )
          R_coordinates_temp(counter_i) = Solution( nodesOnBoundary_temp(counter_i), getVariableId("R") )
        enddo

        
        nodesOnBoundary    = 0
        NodeMaxZlocation   = maxloc(Z_coordinates_temp, dim=1)
        nodesOnBoundary(1) = nodesOnBoundary_temp(NodeMaxZlocation)
        nodesOnBoundary_temp(NodeMaxZlocation) = -1

        outer_loop:do counter_i = 1, size(nodesOnBoundary_temp)


          point_ref = [ Solution( nodesOnBoundary(counter_i), getVariableId("Z") ) , &

                        Solution( nodesOnBoundary(counter_i), getVariableId("R") )  ]
          
          Dist = 40000.d0
          inner_loop:do counter_j = 1, size(nodesOnBoundary_temp)

            if ( ( nodesOnBoundary_temp(counter_j) .eq. nodesOnBoundary(counter_i) ) &
            .or. ( nodesOnBoundary_temp(counter_j) .lt. 0) ) cycle inner_loop

            point_next = [ Solution( nodesOnBoundary_temp(counter_j), getVariableId("Z") ), &
                           Solution( nodesOnBoundary_temp(counter_j), getVariableId("R") ) ]

            call distance(point_ref, point_next, Dist(counter_j))
          enddo inner_loop

          NodeMinDist = minloc(Dist, dim=1)

          nodesOnBoundary(counter_i+1) = nodesOnBoundary_temp(NodeMinDist)
         
          nodesOnBoundary_temp(NodeMinDist) = -1

        enddo outer_loop

        ValuesOnNodes = Solution(nodesOnBoundary(:), getVariableId(variable))
        ! GlobalNodes   = nodesOnBoundary
      End Subroutine getVariableFromBoundaryNodes



      Subroutine getBoundaryNodesOfWholeBoundary( nelements, bnd_elements, faces, all_nodes )
        
        Implicit None
        Integer,                            intent(in)  :: nelements
        Integer, Dimension(:),              intent(in)  :: bnd_elements, faces
        Integer, Dimension(:), Allocatable, intent(out) :: all_nodes

        Integer, Dimension(:), Allocatable :: unsorted_list
        Integer, Dimension(:), Allocatable :: sorted_list
        Integer, Dimension(:), Allocatable :: unsorted_list_zeroes

        Integer, Dimension(2) :: elementBoundaryNodes
        Integer :: counter, i, j, element, face, temp_var

        if (allocated(all_nodes)) deallocate(all_nodes)

        allocate( unsorted_list       ( 2*nelements    ) )
        allocate( unsorted_list_zeroes( 2*nelements    ) )


        unsorted_list  = 0   ;  unsorted_list_zeroes = 0

        counter = 1

        do i = 1, nelements
          element = bnd_elements(i)
          face    = faces(i)
          elementBoundaryNodes(:) = getBoundaryNodesOfTheElement(element, face)
          unsorted_list(counter:counter+1) = elementBoundaryNodes(:)
          counter = counter + 2
          ! Write(*,*) elementBoundaryNodes(1), elementBoundaryNodes(2)
          ! pause
        enddo

        unsorted_list_zeroes = unsorted_list
        
        loop_over_values: do i = 1, size(unsorted_list)-1
                    temp_var = unsorted_list(i)
                    if (unsorted_list_zeroes(i) .eq. 0) cycle loop_over_values
                    comparison_loop: do j = i+1, size(unsorted_list)   
                                if ( unsorted_list(j) .eq. temp_var ) then
                                  unsorted_list_zeroes(j) = 0
                                  exit comparison_loop
                                endif
                              enddo comparison_loop
                  enddo loop_over_values

        counter = 0

        do i = 1, size(unsorted_list_zeroes)
          
          if (unsorted_list_zeroes(i) .ne. 0) then
            counter = counter + 1
          endif
        enddo
        allocate( sorted_list(counter))
        allocate( all_nodes  (counter))
        sorted_list  = 0
        counter = 0
        do i = 1, size(unsorted_list_zeroes)
          
          if (unsorted_list_zeroes(i) .ne. 0) then
            counter = counter + 1
            sorted_list(counter) = unsorted_list_zeroes(i)
          endif
        enddo

        all_nodes(:) = sorted_list

        deallocate( unsorted_list        )
        deallocate( sorted_list          )
        deallocate( unsorted_list_zeroes )

      End Subroutine getBoundaryNodesOfWholeBoundary



END MODULE BOUNDARY_ENUMERATION_MODULE
 

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

  
Module DirichletBoundaries
  Use VariableMapping
  Use BOUNDARY_ENUMERATION_MODULE
  Use ENUMERATION_MODULE,          Only: GNTR
  Use FLOW_ARRAYS_MODULE,          Only: Bc_f, Bce_f, B_f, Be_f
  Use GLOBAL_ARRAYS_MODULE,        Only: TL
  Use ELEMENTS_MODULE,             Only: NUNKNOWNS_f
  Use CSR_STORAGE,                 Only: A_f, Ac_f, IA_f, CSR_f, NZ_f, Ar_f, DA_f
  
  Contains


      



    Subroutine updateAllNodesOfTheBoundary( FemValueName, &
                                            bnd_elements, &
                                            bnd_faces   , &
                                            procedure_)
      Implicit None 
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      ! Arguments
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      Character(len=*)     , Intent(In) :: FemValueName
      Integer, Dimension(:), Intent(In) :: bnd_elements
      Integer, Dimension(:), Intent(In) :: bnd_faces

      Interface 
        Subroutine procedure_ (node, FemValueName)
          Implicit None 
          Integer,          Intent(In) :: node
          Character(len=*), Intent(In) :: FemValueName
        End Subroutine procedure_
      End Interface
      
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      ! Local Variables 
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      Integer :: iel 
      Integer :: element 
      Integer :: face
      Integer :: node(2)

      do iel = 1, size(bnd_elements)
        element = bnd_elements(iel)
        face    = bnd_faces   (iel)

        node    = getBoundaryNodesOfTheElement( element, face )

        ! Since the elements are linear only 2 nodes lie in the surface
        call procedure_(node(1), FemValueName )
        call procedure_(node(2), FemValueName )
      end do
    End Subroutine updateAllNodesOfTheBoundary


    Function integrateOverAllElementsOfTheBoundary(bnd_elements, bnd_faces, procedure_) Result(output)
      Implicit None 
      Interface
        Function procedure_ (bnd_elements, bnd_faces) Result(out)
          Implicit None 
          Real(8)             :: out 
          Integer, Intent(In) :: bnd_elements
          Integer, Intent(In) :: bnd_faces 
        End Function procedure_
      End Interface

      Integer, Dimension(:), Intent(In) :: bnd_elements
      Integer, Dimension(:), Intent(In) :: bnd_faces 
      Real(8)                           :: output 
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
      ! Local Variables
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
      Integer :: iel 

      output = 0.d0 
      do iel = 1, size(bnd_elements) 
        output = output + procedure_(bnd_elements(iel), bnd_faces(iel))
      end do 
      
    End Function integrateOverAllElementsOfTheBoundary


    Subroutine ApplyDirichletAtNode_(node, FemValueName, value, FlagNr, gid )
      Implicit None 
      Integer,          Intent(In)           :: node
      Character(len=*), Intent(In)           :: FemValueName
      Real(8),          Intent(In)           :: value
      Character(len=3), Intent(In)           :: FlagNr
      Integer,          Intent(In), Optional :: gid 

      Integer                                :: irow 
      Integer                                :: v_id 
      Logical                                :: isAssociatedWithGlobalVar

      v_id = getVariableId(FemValueName)
      irow = gntr(node) - 1

      isAssociatedWithGlobalVar = .false.
      if (present(gid) .and. gid /= 0) isAssociatedWithGlobalVar = .true. 
        

      B_f (irow +v_id)   = TL(node,v_id) - value
      Bc_f(irow +v_id)   = .TRUE.

      if (FlagNr == "NRP") then 
        call updateJacobian(irow+v_id)
        if (isAssociatedWithGlobalVar)&
          call updateExtraJacobian(gid, irow+v_id, 1.d0 )
      end if  

    End Subroutine ApplyDirichletAtNode_


    Subroutine updateJacobian(irow)
      Implicit None 
      Integer, Intent(In) :: irow

      Integer             :: il 
      Integer             :: iu
      Integer             :: icol 

      il   = IA_f(irow)
      iu   = IA_f(irow+1)-1
      icol = DA_f(irow)
      A_f(il:iu) = 0.d0
      A_f(icol ) = 1.d0 
    End Subroutine updateJacobian
      
    Subroutine updateExtraJacobian(gid, irow, value)
      Implicit None 
      Integer, Intent(In) :: gid
      Integer, Intent(In) :: irow
      Real(8), Intent(In) :: value 


      Ac_f(irow,gid) = value
    End Subroutine updateExtraJacobian



    
    Subroutine ClearRowsOfResidual(node, FemValueName)
      Implicit None 
      Integer         , Intent(In) :: node
      Character(len=*), Intent(In) :: FemValueName

      Integer                      :: v_id
      Integer                      :: irow 

      Select Case(FemValueName)
        Case('Vr' ) ; V_id = 1
        Case('Vz' ) ; V_id = 2
        Case('P'  ) ; V_id = 3
        Case('Z'  ) ; V_id = 4
        Case('R'  ) ; V_id = 5
        Case('Srr') ; V_id = 6
        Case('Srz') ; V_id = 7
        Case('Szz') ; V_id = 8
        Case('Stt') ; V_id = 9
      End Select 

      irow               = gntr(node) - 1 + v_id
      
      B_f (irow)   = 0.0d0 
    End Subroutine ClearRowsOfResidual

    
    Subroutine ClearRowsOfJacobian(node, FemValueName)
      Implicit None
      Integer         , Intent(In) :: node
      Character(len=*), Intent(In) :: FemValueName

      Integer                      :: v_id
      Integer                      :: irow 
      Integer                      :: il  
      Integer                      :: iu

      Select Case(FemValueName)
        Case('Vr' ) ; V_id = 1
        Case('Vz' ) ; V_id = 2
        Case('P'  ) ; V_id = 3
        Case('Z'  ) ; V_id = 4
        Case('R'  ) ; V_id = 5
        Case('Srr') ; V_id = 6
        Case('Srz') ; V_id = 7
        Case('Szz') ; V_id = 8
        Case('Stt') ; V_id = 9
      End Select 

      irow  = gntr(node) - 1 + v_id

      il = ia_f(irow  )
      iu = ia_f(irow+1)-1
      
      A_f (il:iu) = 0.d0
      Ac_f(irow,:)= 0.d0

    End Subroutine ClearRowsOfJacobian

End Module DirichletBoundaries

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------



MODULE GAUSS_MODULE
  
    USE ELEMENTS_MODULE, only: NED_2d, NBF_2d  

    Type BasisFunctions 
      Integer                           :: nbf
      Real(8), Dimension(NBF_2d,NBF_2d) :: bfn 
      Real(8), Dimension(NBF_2d,NBF_2d) :: dfdx1
      Real(8), Dimension(NBF_2d,NBF_2d) :: dfdx2
    End Type BasisFunctions


    Type(BasisFunctions) :: BasisFunctionsAtNodes


    ! NUMBER OF GAUSS POINTS
    INTEGER, PARAMETER :: NGAUSS_1d = 3
    INTEGER, PARAMETER :: NGAUSS_2d = 4
    
    INTEGER, PARAMETER :: NGAUSS_LC = 6

    ! GAUSS POINTS AT UNITY ELEMENT
    REAL(8), DIMENSION(NGAUSS_1d)   :: GAPT_1d
    REAL(8), DIMENSION(NGAUSS_1d)   :: WO_1d
    
    REAL(8), DIMENSION(NGAUSS_2d,2) :: GAPT_2d
    REAL(8), DIMENSION(NGAUSS_2d)   :: WO_2d
    
    REAL(8), DIMENSION(NGAUSS_LC)   :: GAPT_LC
    REAL(8), DIMENSION(NGAUSS_LC)   :: WO_LC

   
    ! LINEAR BASIS FUNCTIONS AND DERIVATIVES AT UNITY ELEMENT IN 1D
    REAL(8), DIMENSION(NBF_2d,NGAUSS_1d,NED_2d) :: BFN_E, DFDC_E, DFDE_E

    ! LINEAR BASIS FUNCTIONS AND DERIVATIVES AT UNITY ELEMENT IN 2D
    REAL(8), DIMENSION(NBF_2d,NGAUSS_2d)        :: BFN_2d, DFDC_2d, DFDE_2d


    
    
    ! ----------------------------------------------------------------------
    contains
    ! ----------------------------------------------------------------------
    Subroutine getBasisFunctionsAtFace( face, bfn, dfdc, dfde )
      ! Returns the basis functions and its derivatives
      ! calculated at gauss points of the face
      Implicit None 
      Integer, Intent(In)                  :: face 
      Real(8), Dimension(:,:), Allocatable :: bfn 
      Real(8), Dimension(:,:), Allocatable :: dfdc 
      Real(8), Dimension(:,:), Allocatable :: dfde

      Integer                              :: nbf_2d_ 
      Integer                              :: npoints_

      If (Allocated(bfn ) ) Deallocate (bfn )
      If (Allocated(dfdc) ) Deallocate (dfdc)
      If (Allocated(dfde) ) Deallocate (dfde)

      nbf_2d_  = size(bfn_e,1)
      npoints_ = size(bfn_e,2)

      Allocate( bfn (nbf_2d_, npoints_) )
      Allocate( dfdc(nbf_2d_, npoints_) )
      Allocate( dfde(nbf_2d_, npoints_) )

      Select Case(face)

        Case(1); bfn(:,:) = bfn_e(:,:,1); dfdc = dfdc_e(:,:,1); dfde = dfde_e(:,:,1)

        Case(2); bfn(:,:) = bfn_e(:,:,3); dfdc = dfdc_e(:,:,3); dfde = dfde_e(:,:,3)

        Case(3); bfn(:,:) = bfn_e(:,:,2); dfdc = dfdc_e(:,:,2); dfde = dfde_e(:,:,2)

      Case Default
      Print*, 'Wrong value of face. Possible values 1,2,3'
      Print*, 'Current Value of face :',face
      stop
      End Select


    End Subroutine getBasisFunctionsAtFace


    Subroutine getNormalVectorAtFace(derivs, face, nr, nz, ds, normalize, forceOnObject)
      Implicit None 
      Real(8), Dimension(:), Intent(In):: derivs
      Integer, Intent(In)              :: face
      Real(8), Intent(InOut)           :: nr
      Real(8), Intent(InOut)           :: nz
      Real(8), Intent(InOut), Optional :: ds
      Logical, Intent(In)   , Optional :: normalize
      Logical, Intent(In)   , Optional :: forceOnObject


      Real(8)                          :: ds_
      Real(8)                          :: dzdc
      Real(8)                          :: dzde
      Real(8)                          :: drdc
      Real(8)                          :: drde

      dzdc = derivs(1) ; dzde = derivs(2)
      drdc = derivs(3) ; drde = derivs(4) 

      Select Case(face)
      Case(1); nr = - dzdc        ; nz =   drdc 
      Case(2); nr =  (dzdc-dzde)  ; nz = -(drdc-drde)
      Case(3); nr =   dzde        ; nz = - drde  
      Case Default
        Print*, "[Error] : getNormalVectorAtFace. Wrong Value of face. Possible values 1,2,3."
      End Select
      
      ds_ = dsqrt(nr*nr+nz*nz)
      if ( present(ds) ) ds = ds_

      if ( present(normalize) .and. normalize ) then
        nr = nr/ds
        nz = nz/ds  
      end if 
      
      if ( present(forceOnObject) .and. forceOnObject ) then
        nr = -nr
        nz = -nz
      endif
      
    End Subroutine getNormalVectorAtFace


    SUBROUTINE GAUSS_EVALUATION 

      CALL GAUSS_EVALUATION_1d
      CALL GAUSS_EVALUATION_2d
      
      CALL GAUSS_EVALUATION_LC

    END SUBROUTINE GAUSS_EVALUATION

 
    !----------------------------------------------------------------------- 


  SUBROUTINE GAUSS_EVALUATION_1d
         
     IMPLICIT NONE
     
   ! LOCAL VARIABLE
     INTEGER :: I, IG
     
     REAL(8), DIMENSION(NBF_2d) :: BFN, DFDC, DFDE



     GAPT_1d(1) = (1.d0 - DSQRT(3.D0/5.D0))*0.50
     GAPT_1d(2) = (1.d0)                   *0.50
     GAPT_1d(3) = (1.d0 + DSQRT(3.D0/5.D0))*0.50

     WO_1d(1)   = (5.D0/9.D0)*0.5d0
     WO_1d(2)   = (8.D0/9.D0)*0.5d0
     WO_1d(3)   = (5.D0/9.D0)*0.5d0

     IG = 0
     
     DO I = 1, NGAUSS_1d

       IG = IG  + 1

     ! EDGE 1 (C,E=0)
       CALL LINEAR_TRIANGLE( BFN, DFDC, DFDE, NBF_2d, GAPT_1d(I), 0.D0 )

       BFN_E (:,IG,1) = BFN
       DFDC_E(:,IG,1) = DFDC
       DFDE_E(:,IG,1) = DFDE

     ! EDGE 2 (C=0,E)
       CALL LINEAR_TRIANGLE( BFN, DFDC, DFDE, NBF_2d, 0.D0, GAPT_1d(I) )

       BFN_E(:,IG,2)  = BFN
       DFDC_E(:,IG,2) = DFDC
       DFDE_E(:,IG,2) = DFDE

     ! EDGE 3 (C,E=1-C)
       CALL LINEAR_TRIANGLE( BFN, DFDC, DFDE, NBF_2d, GAPT_1d(I), 1.D0-GAPT_1d(I) )

       BFN_E(:,IG,3)  = BFN
       DFDC_E(:,IG,3) = DFDC
       DFDE_E(:,IG,3) = DFDE
       
     ENDDO

  END SUBROUTINE GAUSS_EVALUATION_1d

 
  !----------------------------------------------------------------------- 



  SUBROUTINE GAUSS_EVALUATION_2d

   ! LOCAL VARIABLE
     INTEGER :: I, J, IG
     
     REAL(8), DIMENSION(NBF_2d) :: BFN, DFDC, DFDE


     GAPT_2d(1,1) = 0.3333333333333333D0
     GAPT_2d(1,2) = 0.3333333333333333D0
     
     GAPT_2d(2,1) = 0.2D0
     GAPT_2d(2,2) = 0.2D0
     
     GAPT_2d(3,1) = 0.2D0
     GAPT_2d(3,2) = 0.6D0
     
     GAPT_2d(4,1) = 0.6D0
     GAPT_2d(4,2) = 0.2D0

     WO_2d(1)     = (- 0.5625D0          )*0.5D0
     WO_2d(2)     = (0.5208333333333333D0)*0.5D0
     WO_2d(3)     = (0.5208333333333333D0)*0.5D0
     WO_2d(4)     = (0.5208333333333333D0)*0.5D0

     IG = 0

     DO I = 1, NGAUSS_2d

       IG = IG  + 1

     ! EVALUATE BASIS FUNCTIONS AND THEIR PARTIAL DERIVATIVES
     ! AT GAUSS LEGENDRE POINTS IN THE UNIT ELEMENT COORDINATES
       CALL LINEAR_TRIANGLE( BFN, DFDC, DFDE, NBF_2d, GAPT_2d(I,1), GAPT_2d(I,2) )

       BFN_2d (:,IG) = BFN
       DFDC_2d(:,IG) = DFDC              ! C-> ksi
       DFDE_2d(:,IG) = DFDE              ! E-> eta
     
     ENDDO
        

    ! Node 1: At x1 = 0.0 x2 = 0.0 
    call LINEAR_TRIANGLE(BFN, DFDC, DFDE, NBF_2d,  0.d0, 0.d0 )
    BasisFunctionsAtNodes%nbf        = NBF_2d
    BasisFunctionsAtNodes%bfn  (:,1) = bfn 
    BasisFunctionsAtNodes%dfdx1(:,1) = dfdc
    BasisFunctionsAtNodes%dfdx2(:,1) = dfde

    ! Node 2: At x1 = 1.0 x2 = 0.0 
    call LINEAR_TRIANGLE(BFN, DFDC, DFDE, NBF_2d,  1.d0, 0.d0 )
    BasisFunctionsAtNodes%nbf        = NBF_2d
    BasisFunctionsAtNodes%bfn  (:,2) = bfn 
    BasisFunctionsAtNodes%dfdx1(:,2) = dfdc
    BasisFunctionsAtNodes%dfdx2(:,2) = dfde

    ! Node 3: At x1 = 0.0 x2 = 1.0 
    call LINEAR_TRIANGLE(BFN, DFDC, DFDE, NBF_2d,  0.d0, 1.d0 )
    BasisFunctionsAtNodes%nbf        = NBF_2d
    BasisFunctionsAtNodes%bfn  (:,3) = bfn 
    BasisFunctionsAtNodes%dfdx1(:,3) = dfdc
    BasisFunctionsAtNodes%dfdx2(:,3) = dfde

  END SUBROUTINE GAUSS_EVALUATION_2d

 
  !----------------------------------------------------------------------- 


   SUBROUTINE GAUSS_EVALUATION_LC
         
     IMPLICIT NONE

     GAPT_LC(1) = 0.5D0*(1.D0-0.9324695142031520278123D0)
     GAPT_LC(2) = 0.5D0*(1.D0-0.6612093864662645136610D0)
     GAPT_LC(3) = 0.5D0*(1.D0-0.2386191860831969086305D0)
     GAPT_LC(4) = 0.5D0*(1.D0+0.2386191860831969086305D0)
     GAPT_LC(5) = 0.5D0*(1.D0+0.6612093864662645136610D0)
     GAPT_LC(6) = 0.5D0*(1.D0+0.9324695142031520278123D0)

     WO_LC(1)   = 0.5D0*0.1713244923791703450403
     WO_LC(2)   = 0.5D0*0.3607615730481386075698
     WO_LC(3)   = 0.5D0*0.4679139345726910473899
     WO_LC(4)   = 0.5D0*0.4679139345726910473899
     WO_LC(5)   = 0.5D0*0.3607615730481386075698
     WO_LC(6)   = 0.5D0*0.1713244923791703450403

   END SUBROUTINE GAUSS_EVALUATION_LC

 
  !----------------------------------------------------------------------- 


  SUBROUTINE LINEAR_TRIANGLE( BFN, DFDC, DFDE, NBF, C, E )

     IMPLICIT NONE

   ! ARGUMENTS
     INTEGER, INTENT(IN) :: NBF
     REAL(8), INTENT(IN) :: C, E
     
     REAL(8), INTENT(INOUT), DIMENSION(NBF) :: BFN, DFDC, DFDE

     BFN(1)   =  1.D0 - C - E
     BFN(2)   =  C
     BFN(3)   =  E

     DFDC(1) = -1.D0
     DFDC(2) =  1.D0
     DFDC(3) =  0.D0

     DFDE(1) = -1.D0
     DFDE(2) =  0.D0
     DFDE(3) =  1.D0

  END SUBROUTINE LINEAR_TRIANGLE

END MODULE GAUSS_MODULE
    

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module basis_calculations
  
  contains

  subroutine inverse_derivatives(dA, dA_inverse, Jac)
    implicit none
    real(8), dimension(:), intent(in)  :: dA
    real(8), dimension(:), intent(out) :: dA_inverse
    real(8),               intent(out) :: Jac
  
    Jac = dA(4) * dA(1) - dA(3) * dA(2)
    dA_inverse(1)  =   dA(4)/Jac
    dA_inverse(3)  = - dA(2)/Jac
    dA_inverse(2)  = - dA(3)/Jac
    dA_inverse(4)  =   dA(1)/Jac
  end subroutine inverse_derivatives


  subroutine basis_spatial_derivatives( gauss_point, Z_nodes, R_nodes, dBFNdZ, dBFNdR, Jac )
    Use ELEMENTS_MODULE,  only: NBF_2d, NEQ_f
    use GAUSS_MODULE,     only: NGAUSS_2d, BFN_2d, DFDC_2d, DFDE_2d
    implicit none
    integer,               intent(in)  :: gauss_point
    real(8), dimension(:), intent(in)  :: Z_nodes, R_nodes
    real(8), dimension(:), intent(out) :: dBFNdZ, dBFNdR
    Real(8)              , intent(out) :: Jac
    ! local variables
    integer                              :: II
    Real(8), Dimension(NBF_2d,NGAUSS_2d) :: BFN, dBFNdx1, dBFNdx2
    Real(8)                              :: R, dRdx1, dRdx2, dx1dR, dx2dR
    Real(8)                              :: Z, dZdx1, dZdx2, dx1dZ, dx2dZ
    Real(8), dimension(4)                :: outpt


    BFN     = BFN_2d
    dBFNdx1 = DFDC_2d
    dBFNdx2 = DFDE_2d

    Z     = 0.d0 ; dZdx1 = 0.d0 ; dZdx2 = 0.d0
    do II = 1, NBF_2d
      Z     = Z     + Z_nodes(II)*BFN    (II, gauss_point)
      dZdx1 = dZdx1 + Z_nodes(II)*dBFNdx1(II, gauss_point)
      dZdx2 = dZdx2 + Z_nodes(II)*dBFNdx2(II, gauss_point)
    enddo 

    R     = 0.d0 ; dRdx1 = 0.d0 ; dRdx2 = 0.d0
    do II = 1, NBF_2d
      R     = R     + R_nodes(II)*BFN    (II, gauss_point)
      dRdx1 = dRdx1 + R_nodes(II)*dBFNdx1(II, gauss_point)
      dRdx2 = dRdx2 + R_nodes(II)*dBFNdx2(II, gauss_point)
    enddo 


    call inverse_derivatives( [dZdx1, dZdx2, dRdx1, dRdx2], outpt, Jac )

    dx1dZ = outpt(1)
    dx2dZ = outpt(2)
    dx1dR = outpt(3)
    dx2dR = outpt(4)
    

    dBFNdZ  = 0.d0 ; dBFNdR  = 0.d0
    do II = 1, NBF_2d

      dBFNdZ(II) = dx1dZ * dBFNdx1(II, gauss_point)  +  dx2dZ * dBFNdx2(II, gauss_point)
      dBFNdR(II) = dx1dR * dBFNdx1(II, gauss_point)  +  dx2dR * dBFNdx2(II, gauss_point)

    enddo 
  end subroutine basis_spatial_derivatives


  
  subroutine basis_interpolation_chainrule( value_at_nodes_of_element, gauss_point, Z_nodes, R_nodes, VAR, dVARdZ, dVARdR)
    Use ELEMENTS_MODULE,  only: NBF_2d, NEQ_f
    use GAUSS_MODULE,     only: NGAUSS_2d, BFN_2d, DFDC_2d, DFDE_2d
    implicit none
    real(8), dimension(:), intent(in)  :: value_at_nodes_of_element, Z_nodes, R_nodes
    integer,               intent(in)  :: gauss_point
    real(8),               intent(out) :: VAR, dVARdZ, dVARdR !values of the variable on the gauss point

    ! local variables
    Real(8), Dimension(NBF_2d)           :: BFN, dBFNdZ, dBFNdR
    integer                              :: II
    Real(8)                              :: dVARdZ_dummy, dVARdR_dummy, Jac


    BFN     = BFN_2d (:, gauss_point)

    call basis_spatial_derivatives( gauss_point, Z_nodes, R_nodes, dBFNdZ, dBFNdR, Jac )


    VAR     = 0.d0 
    dVARdZ_dummy  = 0.d0 ; dVARdR_dummy  = 0.d0
    do II = 1, NBF_2d

      VAR          = VAR          + value_at_nodes_of_element(II)*BFN   (II)
      dVARdZ_dummy = dVARdZ_dummy + value_at_nodes_of_element(II)*dBFNdZ(II)
      dVARdR_dummy = dVARdR_dummy + value_at_nodes_of_element(II)*dBFNdR(II)

    enddo 

    dVARdZ  = dVARdZ_dummy
    dVARdR  = dVARdR_dummy

  end subroutine basis_interpolation_chainrule


  subroutine Hugn_calculation(velocity_z, velocity_r, gauss_point, Z_nodes, R_nodes, Hugn)
    Use ELEMENTS_MODULE,  only: NBF_2d, NEQ_f
    use GAUSS_MODULE,     only: NGAUSS_2d, BFN_2d, DFDC_2d, DFDE_2d
    Implicit none
    real(8),               intent(in)  :: velocity_z, velocity_r
    integer,               intent(in)  :: gauss_point
    real(8), dimension(:), intent(in)  :: Z_nodes, R_nodes
    real(8),               intent(out) :: Hugn

    ! local variables
    Real(8), Dimension(NBF_2d)           :: BFN, dBFNdZ, dBFNdR
    Real(8)                              :: dummy, Jac
    Real(8), dimension(4)                :: outpt
    integer                              :: II


    BFN     = BFN_2d (:, gauss_point)

    call basis_spatial_derivatives( gauss_point, Z_nodes, R_nodes, dBFNdZ, dBFNdR, Jac )


    Hugn  = 0.d0
    Dummy = 0.d0
     do II = 1, NBF_2d

       Dummy = Dummy + dabs( velocity_z* dBFNdZ(II) + velocity_r*dBFNdR(II) )

     enddo 

    Hugn  = DSQRT(velocity_z**2+velocity_r**2)/(Dummy + 1.d-8 )

  end subroutine Hugn_calculation

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

end module basis_calculations





module sr_representation      
  USE PHYSICAL_MODULE        
  implicit none        
  
  contains

  subroutine inverse_tensor( s, si )
    implicit none
    real(8), dimension(2,2), intent(in)   :: s
    real(8), dimension(2,2), intent(inout):: si
    real(8)                               :: dets
    dets = s(1,1)*s(2,2) - s(1,2)*s(1,2)
    si = 0.d0
    si(1,1) =   s(2,2)
    si(1,2) = - s(1,2)
    si(2,1) = - s(2,1)
    si(2,2) =   s(1,1)
    si = si/dets
  end subroutine inverse_tensor
  
  subroutine SR_CONFORMATION( s, gut, st, si)
    implicit none
    real(8), dimension(2,2), intent(in) :: s, gut
    real(8), dimension(2,2), intent(out):: st, si
    real(8)                             :: wmega
    real(8), dimension(2,2)             :: gu, w

    call inverse_tensor( s, si )
    
    w = 0.d0
    wmega = ( (s(1,2)*gut(1,1)+s(2,2)*gut(1,2)) - (s(1,1)*gut(2,1)+s(1,2)*gut(2,2)) )/(s(1,1)+s(2,2))
    w(1,2) =   wmega
    w(2,1) = - wmega

    !    define velocity gradient
    gu = transpose(gut)
    
    st = - matmul(s,gu) - matmul(w,s)
  end subroutine SR_CONFORMATION
  



  ! subroutine set_S_get_Stress( S, Stress)
  !   implicit none
  !   real(8), dimension(3,3), intent(in)  :: s
  !   real(8), dimension(3,3), intent(out) :: Stress

    
  !   real(8), dimension(3,3) :: conformation, unity   
    
  !   unity = 0.d0

  !   unity(1,1) = 1.d0  ;  unity(2,2) = 1.d0  ;  unity(3,3) = 1.d0 
    
  !   conformation = matmul(S,S)

  !   Stress = (conformation - unity) / WiN

  ! end subroutine set_S_get_Stress

end module sr_representation


     
MODULE LOG_REPRESENTATION
      
   USE PHYSICAL_MODULE
        
   IMPLICIT NONE
    
    
    ! ----------------------------------------------------------------------
    contains
    ! ----------------------------------------------------------------------


  SUBROUTINE LOG_CONFORMATION( S, GUT, EXP_S, EXP_mS, R, B, W )
        
     IMPLICIT NONE
           
   ! ARGUMENTS
     REAL(8), DIMENSION(2,2), INTENT(IN) :: S, GUT
     REAL(8), DIMENSION(2,2), INTENT(OUT):: EXP_S, EXP_mS, R, B, W
        
   ! LOCAL VARIABLES
     REAL(8)                             :: L1, L2, wmega
     REAL(8), DIMENSION(2)               :: EV1, EV2
     REAL(8), DIMENSION(2,2)             :: RT, M, MD     
        
   ! CALCULATE EIGENVALUES & EIGENVECTORS OF S TENSOR
     CALL EIGENSPECTRUM( S(1,1), S(1,2), S(2,2), L1, L2, EV1, EV2 )
        
   ! DEFINE R & R TRANSPOSE TENSORS
     R(1,1) = EV1(1)
     R(2,1) = EV1(2)
     R(1,2) = EV2(1)
     R(2,2) = EV2(2)
        
     RT = TRANSPOSE(R)
     
     IF (DABS(S(1,2)).LT.1.D-9) THEN
        
     ! DEFINE B TENSOR
       B = 0.5D0*(TRANSPOSE(GUT)+GUT)
          
     ! DEFINE W TENSOR
       W = 0.D0
       
     ELSE
        
     ! DEFINE M TENSOR & ITS DIAGONAL
       M = MATMUL(MATMUL(RT,GUT),R)
        
       MD      = 0.D0
       MD(1,1) = M(1,1)
       MD(2,2) = M(2,2)
        
     ! DEFINE B TENSOR
       B = MATMUL(MATMUL(R,MD),RT)
        
     ! DEFINE W TENSOR
       wmega = (DEXP(L2)*M(1,2)+DEXP(L1)*M(2,1))/(DEXP(L2)-DEXP(L1))
        
       W = 0.D0
       W(1,2) =  wmega
       W(2,1) = -wmega
        
       W = MATMUL(MATMUL(R,W),RT)
       
     ENDIF

   ! DEFINE EXPONENTIAL TENSOR-FUNCTION OF S
     EXP_S = 0.D0
     EXP_S(1,1) = DEXP(L1)
     EXP_S(2,2) = DEXP(L2)
        
     EXP_S = MATMUL(MATMUL(R,EXP_S),RT)
      
   ! DEFINE EXPONENTIAL TENSOR-FUNCTION OF S
     EXP_mS = 0.D0
     EXP_mS(1,1) = DEXP(-L1)
     EXP_mS(2,2) = DEXP(-L2)
        
     EXP_mS = MATMUL(MATMUL(R,EXP_mS),RT)
     
  END SUBROUTINE LOG_CONFORMATION

 
  !----------------------------------------------------------------------- 


  SUBROUTINE EIGENSPECTRUM( Axx, Ayx, Ayy, L1, L2, EV1, EV2 )
        
     IMPLICIT NONE
           
   ! ARGUMENTS
     REAL(8),               INTENT(IN)  :: Axx, Ayx, Ayy
     REAL(8),               INTENT(OUT) :: L1, L2
     REAL(8), DIMENSION(2), INTENT(OUT) :: EV1, EV2
        
   ! LOCAL VARIABLES
     REAL(8)                            :: TrA, Dis, Length
    
        
   ! DEFINE TRACE OF THE TENSOR
     TrA = Axx + Ayy
        
   ! DEFINE DISCRIMINANT OF THE CHARACTERISTIC EQUATION
     Dis = (Axx-Ayy)**2 + 4.D0*Ayx**2

   ! CALCULATE EIGENVALUES
     L1 = (TrA+DSQRT(Dis))/2.D0
        
     L2 = L1 - DSQRT(Dis)

     
   ! CALCULATE EIGENVECTORS
     IF ((Axx-Ayy).GT.0.D0) THEN
        
       EV1(1) = Axx - Ayy + DSQRT(Dis)
       EV1(2) = 2.D0*Ayx
       Length = DSQRT(DOT_PRODUCT(EV1,EV1))
       EV1    = EV1/Length
          
       EV2(1) = -EV1(2)
       EV2(2) =  EV1(1)
          
     ELSEIF ((Axx-Ayy).LT.0.D0) THEN
        
       EV1(1) =  2.D0*Ayx
       EV1(2) = -(Axx - Ayy) + DSQRT(Dis)
       Length = DSQRT(DOT_PRODUCT(EV1,EV1))
       EV1    = EV1/Length
          
       EV2(1) = -EV1(2)
       EV2(2) =  EV1(1)
          
     ELSEIF (Ayx.NE.0.D0) THEN
        
       EV1(1) =  SIGN( DSQRT(2.D0)/2.D0, Ayx )
       EV1(2) =  DSQRT(2.D0)/2.D0
          
       EV2(1) = -DSQRT(2.D0)/2.D0
       EV2(2) =  SIGN( DSQRT(2.D0)/2.D0, Ayx )          
          
     ELSE
        
       EV1(1) = 1.D0
       EV1(2) = 0.D0
          
       EV2(1) = 0.D0
       EV2(2) = 1.D0
          
     ENDIF

  END SUBROUTINE EIGENSPECTRUM


END MODULE LOG_REPRESENTATION
    




 

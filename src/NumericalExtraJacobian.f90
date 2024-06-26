Module NumericalBoundaryJacobian
  
  use storage,    only: MATRIX_STORAGE_JACOBIAN

  Interface CalculateJacobianContributionsOf
    Module Procedure NumericalJacobian_Simple
    Module Procedure NumericalJacobian_Stresses
  End Interface CalculateJacobianContributionsOf

  Interface CalculateExtraJacobianContributionsOf
    Module Procedure CalculateExtraJacobianContributionsOf_1_global
    Module Procedure CalculateExtraJacobianContributionsOf_2_global
  End Interface CalculateExtraJacobianContributionsOf
  
  Contains

  Subroutine NumericalJacobian_Simple(Equation, NELEM, NED, TEMP_TL, TEMP_RES )
      !********************************************************************
      ! NumJac : Computes the contributions of the Equation in the Jacobian
      !          matrix.
      !********************************************************************
      Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f, NUNKNOWNS_f
      Use ENUMERATION_MODULE,   Only: NM_f
      Use CSR_STORAGE,          Only: A_f, IA_f, CSR_f, NZ_f, Ac_f, Ar_f

      Implicit None
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      !  ARGUMENTS
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      Integer,                               Intent(In)          :: NELEM, NED
      Real(8), Dimension(NBF_2d,NEQ_f),      Intent(InOut)       :: TEMP_TL
      Real(8), Dimension(NBF_2d,NEQ_f),      Intent(In)          :: TEMP_RES
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      ! Interface Equation 
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      Interface 
        Subroutine Equation ( nelem, ned, temp_tl, temp_res, store )
          Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f
          Integer,                           Intent(In)           :: NELEM, NED
          Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)           :: TEMP_TL
          Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out)          :: TEMP_RES
          LOGICAL,                           Intent(In)           :: STORE
        End Subroutine Equation
      End Interface   

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      !  LOCAL VARIABLES
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      Integer                                                    :: II, JJ, IW, JW, I, J, INOD, IEQ, JNOD, JEQ
      Integer                                                    :: IROW, JCOL, ICOL, IAD, L
      Real(8)                                                    :: F_DX, EPS_JAC
      Integer, Dimension(NBF_2d)                                 :: NM
      Integer, Dimension(NBF_2d*NBF_2d)                          :: CSC
      Real(8), Dimension(NBF_2d,NEQ_f)                           :: dTEMP_RES
      Real(8), Dimension(NBF_2d,NBF_2d, NEQ_f, NEQ_f)            :: TEMP_JAC
      Real(8), Dimension(NBF_2d, NEQ_f)                          :: TEMP_DpL


      !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
      !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
      TEMP_JAC  = 0.D0
            
            
      !  ELEMENTAL JACOBIAN
      !  ITERATE OVER J NODES
      do JW = 1, NBF_2d   
        !    ITERATE OVER J EQUATIONS
        do JEQ = 1, NEQ_f         
          !      ADD A SMALL CHANGE TO VARIABLE  
          EPS_JAC         = F_DX( TEMP_TL(JW,JEQ) )    
          TEMP_TL(JW,JEQ) = TEMP_TL(JW,JEQ) + EPS_JAC
               
          !      RECOMPUTE RESIDUAL WITH A CHANGE IN VARIABLE
          CALL EQUATION( NELEM, NED, TEMP_TL, dTEMP_RES, .FALSE. )
               
          !      RETURN THE ORIGINAL VALUE TO THE UNKNOWN    
          TEMP_TL(JW,JEQ) = TEMP_TL(JW,JEQ) - EPS_JAC
               
          !      COMPUTE DERIVATIVE USING FINITE DIFFERENCES
          TEMP_JAC(:,JW,JEQ,:) = ( dTEMP_RES - TEMP_RES )/EPS_JAC

        end do
      end do


      !  STORE THE ELEMENT INTEGRATION MATRIX IN THE GLOBAL MATRIX A
      NM  = NM_f (NELEM,1:NBF_2d)
      CSC = CSR_f(NELEM,1:NBF_2d*NBF_2d)

      CALL MATRIX_STORAGE_JACOBIAN&
          (TEMP_JAC, NBF_2d, NBF_2d, NEQ_f, NEQ_f, NM, IA_f, NUNKNOWNS_f+1,&
                CSC, NBF_2d*NBF_2d, A_f, NZ_f)
  End Subroutine NumericalJacobian_Simple




  Subroutine NumericalJacobian_Stresses(Equation, nelem, ned, temp_tl, temp_res, globalValue)
    Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f, NUNKNOWNS_f
    Use ENUMERATION_MODULE,   Only: NM_f
    Use CSR_STORAGE,          Only: A_f, IA_f, CSR_f, NZ_f, Ac_f, Ar_f
    Implicit None
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !  ARGUMENTS
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Integer,                          Intent(In)            :: NELEM, NED
    Real(8), Dimension(NBF_2d,NEQ_f), Intent(InOut)         :: TEMP_TL
    Real(8), Dimension(NBF_2d,NEQ_f), Intent(In)            :: TEMP_RES
    Real(8)                         , Intent(In)            :: globalValue
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Interface Equation 
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Interface 
      Subroutine Equation ( nelem, ned, temp_tl, temp_res, store, globalValue )
        Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f
        Integer,                           Intent(In)           :: nelem
        Integer,                           Intent(In)           :: ned
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)           :: temp_tl
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out)          :: temp_res
        Logical,                           Intent(In)           :: store
        Real(8),                           Intent(In)           :: globalValue
      End Subroutine Equation
    End Interface   
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !  LOCAL VARIABLES
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Integer                                         :: II, JJ, IW, JW, I, J, INOD, IEQ, JNOD, JEQ
    Integer                                         :: IROW, JCOL, ICOL, IAD, L
    Real(8)                                         :: F_DX, EPS_JAC
    Integer, Dimension(NBF_2d)                      :: NM
    Integer, Dimension(NBF_2d*NBF_2d)               :: CSC
    Real(8), Dimension(NBF_2d,NEQ_f)                :: dTEMP_RES
    Real(8), Dimension(NBF_2d,NBF_2d, NEQ_f, NEQ_f) :: TEMP_JAC
    Real(8), Dimension(NBF_2d, NEQ_f)               :: TEMP_DpL
    Real(8)                                         :: tmpGlobalValue


    !  initialize working (temporary) areas for element integration
    !  before forming elemental jacobian and rhs vector
    temp_jac       = 0.d0
    tmpGlobalValue = globalValue 
            
    !**********************************************************
    !  elemental jacobian
    !  ** iterate over elements nodes
    !  ** iterate over equations 
    !  ** perturb the node and compute the residual of the 
    !     perturbed equation
    !  ** return the residual in its original value 
    !**********************************************************

    do jw = 1, nbf_2d
      do jeq = 1, neq_f
        eps_jac         = f_dx( temp_tl(jw,jeq) )
        temp_tl(jw,jeq) = temp_tl(jw,jeq) + eps_jac
                    
        call equation( nelem, ned, temp_tl, dtemp_res, .false., tmpGlobalValue )
                    
        temp_tl(jw,jeq) = temp_tl(jw,jeq) - eps_jac
                    
        temp_jac(:,jw,jeq,:) = ( dtemp_res - temp_res )/eps_jac
        
      enddo
    enddo


    !  store the element integration matrix in the global matrix a
    NM  = NM_f (NELEM,1:NBF_2d)
    CSC = CSR_f(NELEM,1:NBF_2d*NBF_2d)

    CALL MATRIX_STORAGE_JACOBIAN&
         (TEMP_JAC, NBF_2d, NBF_2d, NEQ_f, NEQ_f, NM, IA_f, NUNKNOWNS_f+1,&
          CSC, NBF_2d*NBF_2d, A_f, NZ_f)
  end subroutine NumericalJacobian_Stresses


  Subroutine CalculateExtraJacobianContributionsOf_1_global&
    ( Equation, nelem, ned, temp_tl, temp_res, gid, gVal1, store_id)

    Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f, NUNKNOWNS_f
    Use ENUMERATION_MODULE,   Only: NM_f
    Use CSR_STORAGE,          Only: A_f, IA_f, CSR_f, NZ_f, Ac_f, Ar_f
    Implicit None
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !  ARGUMENTS
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Integer,                          Intent(In)            :: NELEM, NED
    Real(8), Dimension(NBF_2d,NEQ_f), Intent(InOut)         :: TEMP_TL
    Real(8), Dimension(NBF_2d,NEQ_f), Intent(In)            :: TEMP_RES
    Integer,                          Intent(In)            :: gid
    Real(8)                         , Intent(In)            :: gVal1
    Integer                         , Intent(In)            :: store_id
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Interface Equation 
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Interface 
      Subroutine Equation ( nelem, ned, temp_tl, temp_res, store, gVal1 )
        Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f
        Integer,                           Intent(In)           :: nelem
        Integer,                           Intent(In)           :: ned
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)           :: temp_tl
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out)          :: temp_res
        Logical,                           Intent(In)           :: store
        Real(8),                           Intent(In)           :: gVal1
      End Subroutine Equation
    End Interface   
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !  LOCAL VARIABLES
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Integer                                         :: II, JJ, IW, JW, I, J, INOD, IEQ, JNOD, JEQ
    Integer                                         :: IROW, JCOL, ICOL, IAD, L
    Real(8)                                         :: F_DX 
    Real(8)                                         :: eps
    Integer, Dimension(NBF_2d)                      :: NM
    Integer, Dimension(NBF_2d*NBF_2d)               :: CSC
    Real(8), Dimension(NBF_2d,NEQ_f)                :: dTEMP_RES
    Real(8), Dimension(NBF_2d,NBF_2d, NEQ_f, NEQ_f) :: TEMP_JAC
    Real(8), Dimension(NBF_2d, NEQ_f)               :: TEMP_DpL
    Real(8)                                         :: gVal


    !  initialize working (temporary) areas for element integration
    !  before forming elemental jacobian and rhs vector
            

    !**********************************************************
    !  ELEMENTAL JACOBIAN FOR FLOW RATE IMPOSITION
    !  ADD A SMALL CHANGE TO VARIABLE  
    !**********************************************************
    if ( gid == 1 ) Then 

      gVal = gVal1 
      eps  = F_DX( gVal1 )

      gVal = gVal + eps 
      call equation( NELEM, NED, TEMP_TL, dTEMP_RES, .FALSE., gVal )
      gVal = gVal - eps 
    Else 
      Print*, "[Error] CalculateExtraJacobianContributionsOf (_1_global)."
      Print*, "        the gid should be 1"
    End if      
    
    !**********************************************************
    !  COMPUTE DERIVATIVE USING FINITE DIFFERENCES
    !**********************************************************
    temp_dpl = ( dtemp_res - temp_res )/eps

    !**********************************************************
    !  STORE THE ELEMENT INTEGRATION MATRIX IN THE GLOBAL MATRIX A
    !**********************************************************
    NM  = NM_f (NELEM,1:NBF_2d)
    CSC = CSR_f(NELEM,1:NBF_2d*NBF_2d)

    !**********************************************************
    !  STORE MATRIX Ac_f
    !**********************************************************
    Call updateJacobianExtraColumn ( store_id, nm, temp_dpl )
  end subroutine CalculateExtraJacobianContributionsOf_1_global


  Subroutine CalculateExtraJacobianContributionsOf_2_global&
    ( Equation, nelem, ned, temp_tl, temp_res, gid, gVal1, gVal2, store_id)

    Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f, NUNKNOWNS_f
    Use ENUMERATION_MODULE,   Only: NM_f
    Use CSR_STORAGE,          Only: A_f, IA_f, CSR_f, NZ_f, Ac_f, Ar_f
    Implicit None
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !  ARGUMENTS
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Integer,                          Intent(In)            :: NELEM, NED
    Real(8), Dimension(NBF_2d,NEQ_f), Intent(InOut)         :: TEMP_TL
    Real(8), Dimension(NBF_2d,NEQ_f), Intent(In)            :: TEMP_RES
    Integer,                          Intent(In)            :: gid
    Real(8)                         , Intent(In)            :: gVal1
    Real(8)                         , Intent(In)            :: gVal2
    Integer                         , Intent(In)            :: store_id
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Interface Equation 
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Interface 
      Subroutine Equation ( nelem, ned, temp_tl, temp_res, gVal1, gVal2, store )
        Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f
        Integer,                           Intent(In)           :: nelem
        Integer,                           Intent(In)           :: ned
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)           :: temp_tl
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out)          :: temp_res
        Real(8),                           Intent(In)           :: gVal1
        Real(8),                           Intent(In)           :: gVal2
        Logical,                           Intent(In)           :: store
      End Subroutine Equation
    End Interface   
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !  LOCAL VARIABLES
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Integer                                         :: II, JJ, IW, JW, I, J, INOD, IEQ, JNOD, JEQ
    Integer                                         :: IROW, JCOL, ICOL, IAD, L
    Real(8)                                         :: F_DX 
    Real(8)                                         :: eps
    Integer, Dimension(NBF_2d)                      :: NM
    Integer, Dimension(NBF_2d*NBF_2d)               :: CSC
    Real(8), Dimension(NBF_2d,NEQ_f)                :: dTEMP_RES
    Real(8), Dimension(NBF_2d,NBF_2d, NEQ_f, NEQ_f) :: TEMP_JAC
    Real(8), Dimension(NBF_2d, NEQ_f)               :: TEMP_DpL
    Real(8)                                         :: gVal


    !  initialize working (temporary) areas for element integration
    !  before forming elemental jacobian and rhs vector
            
    !**********************************************************
    !  ELEMENTAL JACOBIAN FOR FLOW RATE IMPOSITION
    !  ADD A SMALL CHANGE TO VARIABLE  
    !**********************************************************
    if ( gid == 1 ) Then 

      gVal = gVal1 
      eps  = F_DX( gVal1 )
      gVal = gVal + eps 

      call equation( NELEM, NED, TEMP_TL, dTEMP_RES, gVal, gVal2, .FALSE. )
      gVal = gVal - eps 
    else if (gid == 2) Then

      gVal = gVal2 
      eps  = F_DX( gVal2 )
      gVal = gVal + eps 

      call equation( NELEM, NED, TEMP_TL, dTEMP_RES, gVal1, gVal, .FALSE. )
      gVal = gVal - eps 
    else 
      Print*, "[Error] CalculateExtraJacobianContributionsOf_2_global"
      Print*, "Incorrect value of gid."
    End if      
    
    !**********************************************************
    !  COMPUTE DERIVATIVE USING FINITE DIFFERENCES
    !**********************************************************
    temp_dpl = ( dtemp_res - temp_res )/eps

    !**********************************************************
    !  STORE THE ELEMENT INTEGRATION MATRIX IN THE GLOBAL MATRIX A
    !**********************************************************
    NM  = NM_f (NELEM,1:NBF_2d)
    CSC = CSR_f(NELEM,1:NBF_2d*NBF_2d)

    !**********************************************************
    !  STORE MATRIX Ac_f
    !**********************************************************
    Call updateJacobianExtraColumn ( store_id, nm, temp_dpl )
  end subroutine CalculateExtraJacobianContributionsOf_2_global


  Subroutine NumJac( NELEM, NED, TEMP_TL, TEMP_RES, Equation )
      !********************************************************************
      ! NumJac : Computes the contributions of the Equation in the Jacobian
      !          matrix.
      !********************************************************************
      Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f, NUNKNOWNS_f
      Use ENUMERATION_MODULE,   Only: NM_f
      Use CSR_STORAGE,          Only: A_f, IA_f, CSR_f, NZ_f, Ac_f, Ar_f

      Implicit None
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      !  ARGUMENTS
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      Integer,                               Intent(In)          :: NELEM, NED
      Real(8), Dimension(NBF_2d,NEQ_f),      Intent(InOut)       :: TEMP_TL
      Real(8), Dimension(NBF_2d,NEQ_f),      Intent(In)          :: TEMP_RES
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      ! Interface Equation 
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      Interface 
        Subroutine Equation ( nelem, ned, temp_tl, temp_res, store )
          Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f
          Integer,                           Intent(In)           :: NELEM, NED
          Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)           :: TEMP_TL
          Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out)          :: TEMP_RES
          LOGICAL,                           Intent(In)           :: STORE
        End Subroutine Equation
      End Interface   

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      !  LOCAL VARIABLES
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      Integer                                                    :: II, JJ, IW, JW, I, J, INOD, IEQ, JNOD, JEQ
      Integer                                                    :: IROW, JCOL, ICOL, IAD, L
      Real(8)                                                    :: F_DX, EPS_JAC
      Integer, Dimension(NBF_2d)                                 :: NM
      Integer, Dimension(NBF_2d*NBF_2d)                          :: CSC
      Real(8), Dimension(NBF_2d,NEQ_f)                           :: dTEMP_RES
      Real(8), Dimension(NBF_2d,NBF_2d, NEQ_f, NEQ_f)            :: TEMP_JAC
      Real(8), Dimension(NBF_2d, NEQ_f)                          :: TEMP_DpL


      !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
      !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
      TEMP_JAC  = 0.D0
            
            
      !  ELEMENTAL JACOBIAN
      !  ITERATE OVER J NODES
      do JW = 1, NBF_2d   
        !    ITERATE OVER J EQUATIONS
        do JEQ = 1, NEQ_f         
          !      ADD A SMALL CHANGE TO VARIABLE  
          EPS_JAC         = F_DX( TEMP_TL(JW,JEQ) )    
          TEMP_TL(JW,JEQ) = TEMP_TL(JW,JEQ) + EPS_JAC
               
          !      RECOMPUTE RESIDUAL WITH A CHANGE IN VARIABLE
          CALL EQUATION( NELEM, NED, TEMP_TL, dTEMP_RES, .FALSE. )
               
          !      RETURN THE ORIGINAL VALUE TO THE UNKNOWN    
          TEMP_TL(JW,JEQ) = TEMP_TL(JW,JEQ) - EPS_JAC
               
          !      COMPUTE DERIVATIVE USING FINITE DIFFERENCES
          TEMP_JAC(:,JW,JEQ,:) = ( dTEMP_RES - TEMP_RES )/EPS_JAC

        end do
      end do


      !  STORE THE ELEMENT INTEGRATION MATRIX IN THE GLOBAL MATRIX A
      NM  = NM_f (NELEM,1:NBF_2d)
      CSC = CSR_f(NELEM,1:NBF_2d*NBF_2d)

      CALL MATRIX_STORAGE_JACOBIAN&
          (TEMP_JAC, NBF_2d, NBF_2d, NEQ_f, NEQ_f, NM, IA_f, NUNKNOWNS_f+1,&
                CSC, NBF_2d*NBF_2d, A_f, NZ_f)
  End Subroutine NumJac


  Subroutine updateJacobianExtraColumn(id, element, Contributions) 
    Use CSR_STORAGE,          Only: Ac_f
    Implicit None 
    Integer,                 Intent(In) :: id 
    Integer, Dimension(:)  , Intent(In) :: element
    Real(8), Dimension(:,:), Intent(In) :: Contributions

    Integer :: i 
    Integer :: j 
    Integer :: ieq
    Integer :: nbf_2d_
    Integer :: neq_f_

    nbf_2d_ = size(Contributions,1)
    neq_f_  = size(Contributions,2)

    do i = 1, nbf_2d_ 
      do ieq = 1, neq_f_ 
        j = element(i) + IEQ - 1
        Ac_f(j, id) = Ac_f(j, id) + Contributions(i,ieq)
      end do
    end do
  End Subroutine updateJacobianExtraColumn

End Module NumericalBoundaryJacobian


 
module constrainJacobians 
    use storage,    only: MATRIX_STORAGE_RESIDUAL
    contains

      Subroutine jacobianOfConstrain(gid, nelem, ned, Constrain )
        Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f, NUNKNOWNS_f
        Use ENUMERATION_MODULE,   Only: NM_f, NM_MESH
        Use GLOBAL_ARRAYS_MODULE, Only: TL
        Use CSR_STORAGE,          Only: A_f, IA_f, CSR_f, NZ_f, Ar_f
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        ! Interface
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Interface
            Function Constrain (nelem, ned) Result(out)
            Implicit None 
            Integer, Intent(In) :: nelem
            Integer, Intent(In) :: ned
            Real(8)             :: out
            End Function Constrain
        End Interface 
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer,       Intent(In) :: NELEM, NED
        Integer,       Intent(In) :: gid 
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        INTEGER :: II, JJ, IW, JW, I, J, INOD, IEQ, JNOD, JEQ
        INTEGER :: IROW, JCOL, ICOL, IAD, L
        REAL(8) :: F_DX, EPS_JAC
           
        INTEGER, DIMENSION(NBF_2d)        :: NM
        INTEGER, DIMENSION(NBF_2d*NBF_2d) :: CSC
           
        REAL(8) :: TEMP_RES
        REAL(8) :: dTEMP_RES
           
        REAL(8), DIMENSION(NBF_2d, NEQ_f) :: TEMP_JAC
        
        ! *****************************************
        ! the constrain is the volume of the bubble
        ! So temp_res is the unperturbed Vbubble 
        TEMP_RES = Constrain( nelem, ned )
        ! *****************************************
        
        
        ! INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        ! BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        TEMP_JAC  = 0.D0
           
        ! ELEMENTAL JACOBIAN
        ! ITERATE OVER J NODES
        DO JW = 1, NBF_2d
             
            ! ITERATE OVER J EQUATIONS
            DO JEQ = 1, NEQ_f
                
            ! ADD A SMALL CHANGE TO VARIABLE  
            EPS_JAC = F_DX( TL(NM_MESH(NELEM,JW),JEQ) )    
            TL(NM_MESH(NELEM,JW),JEQ) = TL(NM_MESH(NELEM,JW),JEQ) + EPS_JAC
                
            ! *************************************************************** 
            ! RECOMPUTE RESIDUAL WITH A CHANGE IN VARIABLE
            ! dtemp_res is the perturbed volume of the bubble with respect to 
            ! the bulk unknown
            dTEMP_RES = Constrain( nelem, ned)
            ! *************************************************************** 
        
        
            ! RETURN THE ORIGINAL VALUE TO THE UNKNOWN
            TL(NM_MESH(NELEM,JW),JEQ) = TL(NM_MESH(NELEM,JW),JEQ) - EPS_JAC
                
        
            ! *************************************************************** 
            ! COMPUTE DERIVATIVE USING FINITE DIFFERENCES
            ! temp_jac is THE DERIVATIVE OF THE VOLUME WITH RESPECT TO THE BULK UNKNOWN
            TEMP_JAC(JW,JEQ) = ( dTEMP_RES - TEMP_RES )/EPS_JAC
            ! *************************************************************** 
                
            ENDDO
        ENDDO
        
        
        ! STORE MATRIX Ar_f
        NM  = NM_f(NELEM,1:NBF_2d)
           
        DO IW = 1, NBF_2d
            DO IEQ = 1, NEQ_f
            JW = NM(IW) + IEQ - 1
            Ar_f(gid,JW) = Ar_f(gid,JW) + TEMP_JAC(IW,IEQ)
            ENDDO
        ENDDO
      End Subroutine jacobianOfConstrain
        
        
      Subroutine jacobianOfConstrainCentroid( gid )
        ! constrain is the integral Z dV of the interface
        Use DirichletBoundaries,         Only: integrateOverAllElementsOfTheBoundary
        Use BOUNDARY_ENUMERATION_MODULE, Only: bnd3_elements, bnd3_faces, bnd4_elements, bnd4_faces, getBoundaryNodesOfWholeBoundary
        Use ExtraEquations,              Only: SurfaceIntegration, int_Z_dV
        
        Use ELEMENTS_MODULE,      Only: NBF_2d, NEQ_f, NUNKNOWNS_f
        Use ENUMERATION_MODULE,   Only: NM_f, NM_MESH, GNTR
        Use GLOBAL_ARRAYS_MODULE, Only: TL
        Use CSR_STORAGE,          Only: A_f, IA_f, CSR_f, NZ_f, Ar_f
        Implicit None        
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer,       Intent(In) :: gid 
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        INTEGER :: II, JJ, IW, JW, I, J, INOD, IEQ, JNOD, JEQ
        INTEGER :: IROW, JCOL, ICOL, IAD, L
        REAL(8) :: F_DX, EPS_JAC
           
        INTEGER, DIMENSION(NBF_2d)        :: NM
        INTEGER, DIMENSION(NBF_2d*NBF_2d) :: CSC
           
        REAL(8) :: TEMP_RES
        REAL(8) :: dTEMP_RES
        REAL(8) :: volume_unperturbed  , volume_perturbed
        REAL(8) :: int_Z_dV_unperturbed, int_Z_dV_perturbed
           
        ! REAL(8), DIMENSION(NBF_2d, NEQ_f) :: TEMP_JAC
        REAL(8), DIMENSION(:, :), allocatable :: TEMP_JAC
        Integer, dimension(:   ), allocatable :: boundary_nodes

        REAL(8) :: derivative_int_dV, derivative_int_Z_dV
       

        if     (gid .eq. 1) then 
            call getBoundaryNodesOfWholeBoundary( size(bnd3_elements), bnd3_elements, bnd3_faces, boundary_nodes )
            volume_unperturbed   = integrateOverAllElementsOfTheBoundary( bnd3_elements, bnd3_faces, SurfaceIntegration )
            int_Z_dV_unperturbed = integrateOverAllElementsOfTheBoundary( bnd3_elements, bnd3_faces, int_Z_dV )
        elseif (gid .eq. 3) then 
            call getBoundaryNodesOfWholeBoundary( size(bnd4_elements), bnd4_elements, bnd4_faces, boundary_nodes )
            volume_unperturbed   = integrateOverAllElementsOfTheBoundary( bnd4_elements, bnd4_faces, SurfaceIntegration )
            int_Z_dV_unperturbed = integrateOverAllElementsOfTheBoundary( bnd4_elements, bnd4_faces, int_Z_dV )
        endif
        
        allocate( TEMP_JAC( size(boundary_nodes), NEQ_f ) )

        
        TEMP_JAC  = 0.D0
                
        DO JW = 1, size(boundary_nodes)
             
            ! ITERATE OVER J EQUATIONS
            DO JEQ = 1, NEQ_f
              
            derivative_int_dV   = 0.d0
            derivative_int_Z_dV = 0.d0

            ! ADD A SMALL CHANGE TO VARIABLE
            EPS_JAC = F_DX( TL(boundary_nodes(JW),JEQ) )
            TL(boundary_nodes(JW),JEQ) = TL(boundary_nodes(JW),JEQ) + EPS_JAC
                
            
            if (gid .eq. 1) then 
                volume_perturbed = integrateOverAllElementsOfTheBoundary( bnd3_elements, bnd3_faces, SurfaceIntegration )
                int_Z_dV_perturbed = integrateOverAllElementsOfTheBoundary( bnd3_elements, bnd3_faces, int_Z_dV )


            elseif (gid .eq. 3) then 
                volume_perturbed   = integrateOverAllElementsOfTheBoundary( bnd4_elements, bnd4_faces, SurfaceIntegration )
                int_Z_dV_perturbed = integrateOverAllElementsOfTheBoundary( bnd4_elements, bnd4_faces, int_Z_dV )

            endif
            ! *************************************************************** 
        
            ! RETURN THE ORIGINAL VALUE TO THE UNKNOWN
            TL(boundary_nodes(JW),JEQ) = TL(boundary_nodes(JW),JEQ) - EPS_JAC
                
        

            derivative_int_dV   = ( volume_perturbed   - volume_unperturbed)    / eps_jac

            derivative_int_Z_dV = ( int_Z_dV_perturbed - int_Z_dV_unperturbed ) / eps_jac

            ! TEMP_JAC(JW,JEQ) = derivative_int_Z_dV / volume_unperturbed - ( int_Z_dV_unperturbed / volume_unperturbed**2 ) * derivative_int_dV
            
            TEMP_JAC(JW,JEQ) = ( int_Z_dV_perturbed / volume_perturbed - int_Z_dV_unperturbed / volume_unperturbed ) / eps_jac
            ENDDO
        ENDDO
        
        ! STORE MATRIX Ar_f
        ! NM  = NM_f(NELEM,1:NBF_2d)
           
        DO IW = 1, size(boundary_nodes)
            DO IEQ = 1, NEQ_f
            ! STORE MATRIX Ar_f
            JW = GNTR(boundary_nodes(IW)) + IEQ - 1
                    
            ! i separated the jacobian contributions for the PV equations and
            ! centroid and this subroutine only calculates the centroid, so line 2 of Ar_f
            Ar_f(2,JW) = TEMP_JAC(IW,IEQ)
            ENDDO
        ENDDO
        
        
      END SUBROUTINE jacobianOfConstrainCentroid
        

end module constrainJacobians




Module Boundary_Equations
    Use ArrayTools, Only: copyArrayToLocalValues
    use storage,    only: MATRIX_STORAGE_RESIDUAL
    use basis_calculations, only: BASIS_2d


    Contains

    Subroutine Stresses( NELEM, NED, TEMP_TL, TEMP_RES, STORE, gVar )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use CONTINUATION_MODULE,     Only: INCREMENT
        Use ELEMENTS_MODULE,         Only: NBF_2d, NEL_2d, NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, BFN_E, &
                                            getBasisFunctionsAtFace, &
                                            getNormalVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer,                           Intent(In)  :: NELEM, NED
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
        Logical,                           Intent(In)  :: STORE
        Real(8),                           intent(in)  :: gVar 

    
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        ! FEM variables and their derivatives
        Real(8)                              :: R, dRdx1, dRdx2 
        Real(8)                              :: Z, dZdx1, dZdx2     
        ! Basis Functions and their derivatives
        Real(8), Dimension(:,:), Allocatable ::  bfn 
        Real(8), Dimension(:,:), Allocatable :: dbfndx1
        Real(8), Dimension(:,:), Allocatable :: dbfndx2 
        ! Jacobian of Transformation and the reverse derivatives
        Real(8)                              :: JacT
        Real(8)                              :: dx1dR
        Real(8)                              :: dx2dR
        Real(8)                              :: dx1dZ
        Real(8)                              :: dx2dZ
        ! Normal Vector Components
        Real(8)                              :: nr 
        Real(8)                              :: nz 
        Real(8)                              :: dS

        Integer, Dimension(NBF_2d)           :: NM 
        Real(8), Dimension(NEQ_f)            :: TERM_RES
        ! Basis Function 
        Real(8)                              :: BIFN, DBIR, DBIZ
        Integer :: KK, II, JJ, LL, IW, I
        Integer :: IROW, JCOL

        !*********************************************************************
        ! calculate the basis functions at the face of the triangle
        !*********************************************************************
        NM = NM_MESH(NELEM,:)
        call getBasisFunctionsAtFace(ned, bfn, dbfndx1, dbfndx2)

        !*********************************************************************
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !*********************************************************************
        TEMP_RES = 0.D0

        !*********************************************************************
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !*********************************************************************
        LOOP_GAUSS: DO KK = 1, NGAUSS_1d

            !*********************************************************************
            ! Calculate the variation of the the FEM variables
            ! in the parent element
            !*********************************************************************

            R = 0.d0; dRdx1 = 0.d0; dRdx2 = 0.d0
            Z = 0.d0; dZdx1 = 0.d0; dZdx2 = 0.d0
            do ii = 1, nbf_2d
                R     =  R    + TEMP_TL(ii, getVariableId("R")) *  bfn   (ii,kk)
                dRdx1 = dRdx1 + TEMP_TL(ii, getVariableId("R")) * dbfndx1(ii,kk)
                dRdx2 = dRdx2 + TEMP_TL(ii, getVariableId("R")) * dbfndx2(ii,kk)

                Z     =  Z    + TEMP_TL(ii, getVariableId("Z")) *  bfn   (ii,kk)
                dZdx1 = dZdx1 + TEMP_TL(ii, getVariableId("Z")) * dbfndx1(ii,kk)
                dZdx2 = dZdx2 + TEMP_TL(ii, getVariableId("Z")) * dbfndx2(ii,kk)    
            end do

            !*********************************************************************
            ! Calculate the Jacobian of Transformation
            !*********************************************************************
            JacT   = dRdx2 * dZdx1 - dRdx1 * dZdx2
            dx1dZ  =   dRdx2/JacT
            dx1dR  = - dZdx2/JacT
            dx2dZ  = - dRdx1/JacT
            dx2dR  =   dZdx1/JacT

            !*********************************************************************
            ! Calculate the normal vectors with respect to the face of the 
            ! triangle
            !*********************************************************************

            call getNormalVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, nr, nz, dS, normalize = .true.)
            


            ! if (kk==1 .and. store == .true.) then
            !     print*, 'Bubble pressure=', gvar; pause
            !     endif
            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------
            loop_residuals_f:DO IW = 1, NBF_2d
        
                    BIFN =  bfn   (iw,kk)
                    DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                    DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
        
                    ! -n * T = + P_bubble n + 1/Bo * (I-nn)\nabla \cdot \phi
        
                    TERM_RES     = 0.D0
                    ! TERM_RES(getVariableId("Vr"))  = nr*gVar*BIFN*R + (1.d0/BoN)*((1.D0-nr*nr)*DBIR + BIFN/R + (    -nr*nz)*DBIZ)*R
                    ! TERM_RES(getVariableId("Vz"))  = nz*gVar*BIFN*R + (1.d0/BoN)*((    -nr*nz)*DBIR          + (1.D0-nz*nz)*DBIZ)*R
        
                    !      FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_1d(KK) * dS
                              
            end do loop_residuals_f
        end do LOOP_GAUSS
      

        !---------------------------------------------------------------------
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !---------------------------------------------------------------------
        if ( STORE ) then 
            NM = NM_f(NELEM,1:NBF_2d)
            call MATRIX_STORAGE_RESIDUAL ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
        end if               

    end Subroutine Stresses

! ********************************************************************

    Subroutine Kinematic( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        Use CONTINUATION_MODULE,     Only: INCREMENT
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d, NEL_2d, NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, BFN_E, DFDC_E,&
                                       DFDE_E, &
                                                                         getBasisFunctionsAtFace, &
                                                                         getNormalVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb, DpL
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use MESH_MODULE,             Only: Xm, Ym
        Use LOG_REPRESENTATION
        Use TIME_INTEGRATION,        Only: Dt, TIME, DTo, DTb
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer,                           Intent(In)  :: NELEM, NED
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
        Logical,                           Intent(In)  :: STORE
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        ! FEM variables and their derivatives
        Real(8)                                        :: R ,  dRdx1,  dRdx2 
        Real(8)                                        :: Z ,  dZdx1,  dZdx2    
        Real(8)                                        :: Vr, dVrdx1, dVrdx2
        Real(8)                                        :: Vz, dVzdx1, dVzdx2

        ! Basis Functions and their derivatives
        Real(8), Dimension(:,:), Allocatable           ::  bfn 
        Real(8), Dimension(:,:), Allocatable           :: dbfndx1
        Real(8), Dimension(:,:), Allocatable           :: dbfndx2 
        ! Jacobian of Transformation and the reverse derivatives
        Real(8)                                        :: JacT
        Real(8)                                        :: dx1dR
        Real(8)                                        :: dx2dR
        Real(8)                                        :: dx1dZ
        Real(8)                                        :: dx2dZ
        ! Normal Vector Components
        Real(8)                                        :: nr 
        Real(8)                                        :: nz 
        ! Normal Vector Components
        Real(8)                                        :: tr
        Real(8)                                        :: tz
        ! Surface Arc Length
        Real(8)                                        :: dS
    

        Integer                                        :: KK, II,  IW
    
        Real(8)                                        :: Rb
        Real(8)                                        :: Zb
        Real(8)                                        :: BIFN, DBIR, DBIZ, SBFN
          
        Real(8)                                        :: Zo, Ro
        Real(8)                                        :: dRdt, dZdt
    
        Real(8)                                        :: Uelem, tsupg
                      
        Integer, Dimension(NBF_2d)                     :: NM 
        Real(8), Dimension(NEQ_f)                      :: TERM_RES
      

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        NM = NM_MESH(NELEM,:)
        call getBasisFunctionsAtFace(ned, bfn, dbfndx1, dbfndx2)
              
        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0

        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------
        LOOP_GAUSS: DO KK = 1, NGAUSS_1d

            !*********************************************************************
            ! Calculate the variation of the the FEM variables
            ! in the parent element
            !*********************************************************************
            Vr = 0.d0; dVrdx1 = 0.d0; dVrdx2 = 0.d0
            Vz = 0.d0; dVzdx1 = 0.d0; dVzdx2 = 0.d0
            R  = 0.d0;  dRdx1 = 0.d0;  dRdx2 = 0.d0
            Z  = 0.d0;  dZdx1 = 0.d0;  dZdx2 = 0.d0
            do ii = 1, nbf_2d
                R     =  R     + TEMP_TL(ii, getVariableId("R"))  *  bfn   (ii,kk)
                dRdx1 = dRdx1  + TEMP_TL(ii, getVariableId("R"))  * dbfndx1(ii,kk)
                dRdx2 = dRdx2  + TEMP_TL(ii, getVariableId("R"))  * dbfndx2(ii,kk)

                Z     =  Z     + TEMP_TL(ii, getVariableId("Z"))  *  bfn   (ii,kk)
                dZdx1 = dZdx1  + TEMP_TL(ii, getVariableId("Z"))  * dbfndx1(ii,kk)
                dZdx2 = dZdx2  + TEMP_TL(ii, getVariableId("Z"))  * dbfndx2(ii,kk)  

                Vr    =  Vr    + TEMP_TL(ii, getVariableId("Vr")) *  bfn   (ii,kk)
                dVrdx1= dVrdx1 + TEMP_TL(ii, getVariableId("Vr")) * dbfndx1(ii,kk)
                dVrdx2= dVrdx2 + TEMP_TL(ii, getVariableId("Vr")) * dbfndx2(ii,kk)  

                Vz    =  Vz    + TEMP_TL(ii, getVariableId("Vz")) *  bfn   (ii,kk)
                dVzdx1= dVzdx1 + TEMP_TL(ii, getVariableId("Vz")) * dbfndx1(ii,kk)
                dVzdx2= dVzdx2 + TEMP_TL(ii, getVariableId("Vz")) * dbfndx2(ii,kk)  
            end do

            !*********************************************************************
            ! Calculate the Jacobian of Transformation
            !*********************************************************************
            JacT   = dRdx2 * dZdx1 - dRdx1 * dZdx2
            dx1dZ  =   dRdx2/JacT
            dx1dR  = - dZdx2/JacT
            dx2dZ  = - dRdx1/JacT
            dx2dR  =   dZdx1/JacT

            !*********************************************************************
            ! Calculate the normal vectors with respect to the face of the 
            ! triangle
            !*********************************************************************

            call getNormalVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, nr, nz, dS, normalize = .true.)
            ! if (KK .eq. 1) then
            ! print*, nr, nz
            ! pause
            ! endif
            !*********************************************************************
            ! Calculate the Solution of the previous time steps
            ! to take into account the time variation
            !*********************************************************************
            Ro = 0.d0; Zo = 0.d0 
            Rb = 0.d0; Zb = 0.d0    

            do ii = 1, nbf_2d
                Ro = Ro + TLo(NM(ii), getVariableId("R")) * bfn(ii,kk)
                Zo = Zo + TLo(NM(ii), getVariableId("Z")) * bfn(ii,kk)

                Rb = Rb + TLb(NM(ii), getVariableId("R")) * bfn(ii,kk)
                Zb = Zb + TLb(NM(ii), getVariableId("Z")) * bfn(ii,kk)
            end do
            !*********************************************************************
            ! Calculate the time derivatives of the node
            !*********************************************************************

            dRdt = (R - Ro)/Dt
            dZdt = (Z - Zo)/Dt

            !*********************************************************************
            ! Calculate pspg values
            !*********************************************************************
            ! Uelem = dabs( Vr * tr + Vz * tz  )
            ! tsupg = dS/(Uelem + dS/dt)

            !*********************************************************************
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !*********************************************************************
            loop_residuals_f:do iw = 1, nbf_2d

                BIFN =  bfn   (iw,kk)
                DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
            
                SBFN         = BIFN !+ tsupg*(Term_R*tR+Term_Z*tZ)*DFDL(IW,KK)
                TERM_RES     = 0.D0
            
                TERM_RES(getVariableId("Z")) = SBFN * ( nR * (Vr-dRdt) + nZ * (Vz-dZdt) ) * R * dS
                ! TERM_RES(getVariableId("R")) = SBFN * ( nR * (Vr-dRdt) + nZ * (Vz-dZdt) ) * R * dS

                !      FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,:) = TEMP_RES(IW,:) + TERM_RES * WO_1d(KK)
            end do loop_residuals_f

        end do loop_gauss
      
        !*********************************************************************
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !*********************************************************************
        if ( STORE ) then 
            NM = NM_f(NELEM,1:NBF_2d)
            call MATRIX_STORAGE_RESIDUAL ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
        end if

    End Subroutine Kinematic

    Subroutine fixConcentrationFlux( NELEM, NED, TEMP_TL, TEMP_RES, STORE, value )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use CONTINUATION_MODULE,     Only: INCREMENT
        Use ELEMENTS_MODULE,         Only: NBF_2d, NEL_2d, NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, BFN_E, &
                                            getBasisFunctionsAtFace, &
                                            getNormalVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        use MESH_MODULE,            only: Xm, Ym
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer,                           Intent(In)  :: NELEM, NED
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
        Logical,                           Intent(In)  :: STORE
        Real(8),                           Intent(In)  :: value
        ! Real(8),                           intent(in)  :: gVar 

    
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
            !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        ! FEM variables and their derivatives
        Real(8)                              :: R, dRdx1, dRdx2 
        Real(8)                              :: Z, dZdx1, dZdx2     
        ! Basis Functions and their derivatives
        Real(8), Dimension(:,:), Allocatable ::  bfn 
        Real(8), Dimension(:,:), Allocatable :: dbfndx1
        Real(8), Dimension(:,:), Allocatable :: dbfndx2 
        ! Jacobian of Transformation and the reverse derivatives
        Real(8)                              :: JacT
        Real(8)                              :: dx1dR
        Real(8)                              :: dx2dR
        Real(8)                              :: dx1dZ
        Real(8)                              :: dx2dZ
        ! Normal Vector Components
        Real(8)                              :: nr 
        Real(8)                              :: nz 
        Real(8)                              :: dS

        Integer, Dimension(NBF_2d)           :: NM 
        Real(8), Dimension(NEQ_f)            :: TERM_RES
        ! Basis Function 
        Real(8)                              :: BIFN, DBIR, DBIZ
        Integer :: KK, II, JJ, LL, IW, I
        Integer :: IROW, JCOL

        !*********************************************************************
        ! calculate the basis functions at the face of the triangle
        !*********************************************************************
        NM = NM_MESH(NELEM,:)
        call getBasisFunctionsAtFace(ned, bfn, dbfndx1, dbfndx2)

        !*********************************************************************
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !*********************************************************************
        TEMP_RES = 0.D0

        !*********************************************************************
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !*********************************************************************
        LOOP_GAUSS: DO KK = 1, NGAUSS_1d

            !*********************************************************************
            ! Calculate the variation of the the FEM variables
            ! in the parent element
            !*********************************************************************

            R = 0.d0; dRdx1 = 0.d0; dRdx2 = 0.d0
            Z = 0.d0; dZdx1 = 0.d0; dZdx2 = 0.d0
            do ii = 1, nbf_2d
                Z     =  Z    + Xm( NM(ii) ) *  bfn   (ii,kk)
                dZdx1 = dZdx1 + Xm( NM(ii) ) * dbfndx1(ii,kk)
                dZdx2 = dZdx2 + Xm( NM(ii) ) * dbfndx2(ii,kk)    
                
                R     =  R    + Ym( NM(ii) ) *  bfn   (ii,kk)
                dRdx1 = dRdx1 + Ym( NM(ii) ) * dbfndx1(ii,kk)
                dRdx2 = dRdx2 + Ym( NM(ii) ) * dbfndx2(ii,kk)

            end do

            !*********************************************************************
            ! Calculate the Jacobian of Transformation
            !*********************************************************************
            JacT   = dRdx2 * dZdx1 - dRdx1 * dZdx2
            dx1dZ  =   dRdx2/JacT
            dx1dR  = - dZdx2/JacT
            dx2dZ  = - dRdx1/JacT
            dx2dR  =   dZdx1/JacT

            !*********************************************************************
            ! Calculate the normal vectors with respect to the face of the 
            ! triangle
            !*********************************************************************

            call getNormalVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, nr, nz, dS, normalize = .true.)
            


            ! if (kk==1 .and. store == .true.) then
            !     print*, 'Bubble pressure=', gvar; pause
            !     endif
            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------
            loop_residuals_f:DO IW = 1, NBF_2d
        
                    BIFN =  bfn   (iw,kk)
                    DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                    DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
        
                    ! -n * T = + P_bubble n + 1/Bo * (I-nn)\nabla \cdot \phi
        
                    TERM_RES     = 0.D0
                    TERM_RES(getVariableId("C"))  = value*BIFN
                    ! TERM_RES(getVariableId("Vr"))  = nr*gVar*BIFN*R + (1.d0/BoN)*((1.D0-nr*nr)*DBIR + BIFN/R + (    -nr*nz)*DBIZ)*R
                    ! TERM_RES(getVariableId("Vz"))  = nz*gVar*BIFN*R + (1.d0/BoN)*((    -nr*nz)*DBIR          + (1.D0-nz*nz)*DBIZ)*R
        
                    !      FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_1d(KK) * dS
                              
            end do loop_residuals_f
        end do LOOP_GAUSS
      

        !---------------------------------------------------------------------
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !---------------------------------------------------------------------
        if ( STORE ) then 
            NM = NM_f(NELEM,1:NBF_2d)
            call MATRIX_STORAGE_RESIDUAL ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
        end if               

    end Subroutine fixConcentrationFlux

! ********************************************************************

    Subroutine X_EQUIDISTRIBUTION_RESIDUAL_f( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        USE PHYSICAL_MODULE
        USE ELEMENTS_MODULE,     only: NBF_2d, NEL_2d, NEQ_f, NUNKNOWNS_f
        USE GAUSS_MODULE,        only: WO_1d, NGAUSS_1d, &
                                         getBasisFunctionsAtFace
        USE ENUMERATION_MODULE,  only: NM_MESH, NM_f
        USE FLOW_ARRAYS_MODULE,  only: B_f
        USE MESH_MODULE,         only: Ksi_ => Xm, Eta_ =>Ym
        Implicit None
                
        !**************************************************************************
        !  ARGUMENTS
        !**************************************************************************
        Integer,                           Intent(In)  :: NELEM, NED
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
        Logical,                           Intent(In)  :: STORE                       
        !**************************************************************************
        ! LOCAL VARIABLES
        !**************************************************************************
        Real(8)                              :: R  ,   dRdx1,   dRdx2
        Real(8)                              :: Z  ,   dZdx1,   dZdx2
        Real(8)                              :: Ksi, dKsidx1, dKsidx2
        Real(8)                              :: Eta, dEtadx1, dEtadx2
        ! Basis Functions and their derivatives
        Real(8), Dimension(:,:), Allocatable ::  bfn 
        Real(8), Dimension(:,:), Allocatable :: dbfndx1
        Real(8), Dimension(:,:), Allocatable :: dbfndx2 
        ! Jacobian of transformation
        Real(8)                              :: JacC
        Real(8)                              :: dx1dKsi
        Real(8)                              :: dx1dEta
        Real(8)                              :: dx2dKsi
        Real(8)                              :: dx2dEta
        ! Derivatives in Computational Space 
        Real(8)                              :: dRdEta
        Real(8)                              :: dRdKsi
        Real(8)                              :: dZdKsi
        Real(8)                              :: dZdEta
        ! Normal Vector Components
        Real(8)                              :: n_ksi
        Real(8)                              :: n_eta
        Real(8)                              :: dS
        ! Scale Factor 
        Real(8)                              :: QKsi
        ! Basis Functions on Computational Domain   
        Real(8)                              :: BIFN
        Real(8)                              :: DBIEta
        Real(8)                              :: DBIKsi
        ! Element indices
        Integer, Dimension(NBF_2d)           :: NM 
        ! Temporary Residual
        Real(8), Dimension(NEQ_f)            :: TERM_RES
        Integer                              :: KK, II, IW
    
        !*********************************************************************
        ! calculate the basis functions at the face of the triangle
        !*********************************************************************
        call getBasisFunctionsAtFace(ned, bfn, dbfndx1, dbfndx2)
                       
        !*********************************************************************
        ! initialize working (temporary) areas for element integration
        ! before forming elemental jacobian and rhs vector
        !*********************************************************************
        NM       = NM_MESH(NELEM,:)
        TEMP_RES = 0.D0
        !*********************************************************************
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !*********************************************************************

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d

            !*********************************************************************
            ! Calculate the variation of the the FEM variables
            ! in the parent element
            !*********************************************************************
            R   = 0.d0; dRdx1   = 0.d0; dRdx2   = 0.d0
            Z   = 0.d0; dZdx1   = 0.d0; dZdx2   = 0.d0
            Ksi = 0.d0; dKsidx1 = 0.d0; dKsidx2 = 0.d0 
            Eta = 0.d0; dEtadx1 = 0.d0; dEtadx2 = 0.d0

            do ii = 1, nbf_2d
                R       =  R      + TEMP_TL(ii,getVariableId("R")) *  bfn   (ii,kk)
                dRdx1   = dRdx1   + TEMP_TL(ii,getVariableId("R")) * dbfndx1(ii,kk)
                dRdx2   = dRdx2   + TEMP_TL(ii,getVariableId("R")) * dbfndx2(ii,kk)

                Z       =  Z      + TEMP_TL(ii,getVariableId("Z")) *  bfn   (ii,kk)
                dZdx1   = dZdx1   + TEMP_TL(ii,getVariableId("Z")) * dbfndx1(ii,kk)
                dZdx2   = dZdx2   + TEMP_TL(ii,getVariableId("Z")) * dbfndx2(ii,kk)
                
                Ksi     =  Ksi    + Ksi_(NM(ii)) * bfn    (ii,kk)
                dKsidx1 = dKsidx1 + Ksi_(NM(ii)) * dbfndx1(ii,kk)
                dKsidx2 = dKsidx2 + Ksi_(NM(ii)) * dbfndx2(ii,kk)

                Eta     =  Eta    + Eta_(NM(ii)) * bfn    (ii,kk)
                dEtadx1 = dEtadx1 + Eta_(NM(ii)) * dbfndx1(ii,kk)
                dEtadx2 = dEtadx2 + Eta_(NM(ii)) * dbfndx2(ii,kk)
            end do
            !*********************************************************************
            ! Calculate the Jacobian of Transformation
            ! The Jacobian of tranformation in referred on computational domain
            !*********************************************************************
            JacC   = dEtadx2 * dKsidx1 - dEtadx1 * dKsidx2
            dx1dKsi=   dEtadx2/JacC
            dx1dEta= - dKsidx2/JacC
            dx2dKsi= - dEtadx1/JacC
            dx2dEta=   dKsidx1/JacC

            !*********************************************************************
            ! Calculate the normal vectors with respect to the face of the 
            ! triangle
            !*********************************************************************

            Select Case(ned)
            Case(1); n_eta = -  dKsidx1          ; n_ksi =    dEtadx1
            Case(2); n_eta = -( dKsidx1-dKsidx2) ; n_ksi =  ( dEtadx1-dEtadx2)
            Case(3); n_eta = -(-dKsidx2)         ; n_ksi =  (-dEtadx2)  ! face 3 has reverse numbering
            Case Default
                Print*, "[Error] : Stresses. Wrong Value of face. Possible values 1,2,3."
            End Select

            dS = dsqrt(n_ksi**2 + n_eta**2)
            n_eta = n_eta/dS
            n_ksi = n_ksi/dS
                        
            !*********************************************************************
            ! Calculate Derivatives In the Computational Space
            !*********************************************************************
            dRdEta = dRdx1 * dx1dEta + dRdx2 * dx2dEta
            dRdKsi = dRdx1 * dx1dKsi + dRdx2 * dx2dKsi

            dZdEta = dZdx1 * dx1dEta + dZdx2 * dx2dEta
            dZdKsi = dZdx1 * dx1dKsi + dZdx2 * dx2dKsi

            !*********************************************************************
            ! Define Scale Factor
            !*********************************************************************
            QKsi = dsqrt(dZdKsi**2 + dRdKsi**2)

            !*********************************************************************
            ! Iterate over weighting Functions
            !*********************************************************************
            LOOP_RESIDUALS_f:DO IW = 1, NBF_2d

                BIFN    =   bfn  (iw,kk)
                DBIEta  = dbfndx1(iw,kk) * dx1dEta + dbfndx2(iw,kk) * dx2dEta
                DBIKsi  = dbfndx1(iw,kk) * dx1dKsi + dbfndx2(iw,kk) * dx2dKsi
                
                TERM_RES = 0.D0 
                TERM_RES(getVariableId("Z")) = e_bnd*DLOG(QKsi)*DBIKsi * dS

            ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f) * WO_1d(KK)
            ENDDO LOOP_RESIDUALS_f
        ENDDO LOOP_GAUSS
                
        !*********************************************************************
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !*********************************************************************
        if ( STORE ) then
            NM = NM_f(NELEM,1:NBF_2d)
            CALL MATRIX_STORAGE_RESIDUAL( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
        end if

        !*********************************************************************
        ! Deallocate The local arrays
        !*********************************************************************
        deallocate( bfn )
        deallocate( dbfndx1 )
        deallocate( dbfndx2 )

    End Subroutine X_EQUIDISTRIBUTION_RESIDUAL_f

! ********************************************************************

    Subroutine Y_EQUIDISTRIBUTION_RESIDUAL_f( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        USE PHYSICAL_MODULE
        USE ELEMENTS_MODULE,     only: NBF_2d, NEL_2d, NEQ_f, NUNKNOWNS_f
        USE GAUSS_MODULE,        only: WO_1d, NGAUSS_1d, &
                                         getBasisFunctionsAtFace
        USE ENUMERATION_MODULE,  only: NM_MESH, NM_f
        USE FLOW_ARRAYS_MODULE,  only: B_f
        USE MESH_MODULE,         only: Ksi_ => Xm, Eta_ =>Ym
        Implicit None
                
        !**************************************************************************
        !  ARGUMENTS
        !**************************************************************************
        Integer,                           Intent(In)  :: NELEM, NED
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
        Logical,                           Intent(In)  :: STORE                       
        !**************************************************************************
        ! LOCAL VARIABLES
        !**************************************************************************
        Real(8)                              :: R  ,   dRdx1,   dRdx2
        Real(8)                              :: Z  ,   dZdx1,   dZdx2
        Real(8)                              :: Ksi, dKsidx1, dKsidx2
        Real(8)                              :: Eta, dEtadx1, dEtadx2
        ! Basis Functions and their derivatives
        Real(8), Dimension(:,:), Allocatable ::  bfn 
        Real(8), Dimension(:,:), Allocatable :: dbfndx1
        Real(8), Dimension(:,:), Allocatable :: dbfndx2 
        ! Jacobian of transformation
        Real(8)                              :: JacC
        Real(8)                              :: dx1dKsi
        Real(8)                              :: dx1dEta
        Real(8)                              :: dx2dKsi
        Real(8)                              :: dx2dEta
        ! Derivatives in Computational Space 
        Real(8)                              :: dRdEta
        Real(8)                              :: dRdKsi
        Real(8)                              :: dZdKsi
        Real(8)                              :: dZdEta
        ! Normal Vector Components
        Real(8)                              :: n_ksi
        Real(8)                              :: n_eta
        Real(8)                              :: dS
        ! Scale Factor 
        Real(8)                              :: QEta    
        ! Basis Functions on Computational Domain   
        Real(8)                              :: BIFN
        Real(8)                              :: DBIEta
        Real(8)                              :: DBIKsi
        ! Element indices
        Integer, Dimension(NBF_2d)           :: NM 
        ! Temporary Residual
        Real(8), Dimension(NEQ_f)            :: TERM_RES
        Integer                              :: KK, II, IW
    
        !*********************************************************************
        ! calculate the basis functions at the face of the triangle
        !*********************************************************************
        call getBasisFunctionsAtFace(ned, bfn, dbfndx1, dbfndx2)
                       
        !*********************************************************************
        ! initialize working (temporary) areas for element integration
        ! before forming elemental jacobian and rhs vector
        !*********************************************************************
        NM       = NM_MESH(NELEM,:)
        TEMP_RES = 0.D0
        !*********************************************************************
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !*********************************************************************

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d

            !*********************************************************************
            ! Calculate the variation of the the FEM variables
            ! in the parent element
            !*********************************************************************
            R   = 0.d0; dRdx1   = 0.d0; dRdx2   = 0.d0
            Z   = 0.d0; dZdx1   = 0.d0; dZdx2   = 0.d0
            Ksi = 0.d0; dKsidx1 = 0.d0; dKsidx2 = 0.d0 
            Eta = 0.d0; dEtadx1 = 0.d0; dEtadx2 = 0.d0

            do ii = 1, nbf_2d
                R       =  R      + TEMP_TL(ii,getVariableId("R")) *  bfn   (ii,kk)
                dRdx1   = dRdx1   + TEMP_TL(ii,getVariableId("R")) * dbfndx1(ii,kk)
                dRdx2   = dRdx2   + TEMP_TL(ii,getVariableId("R")) * dbfndx2(ii,kk)

                Z       =  Z      + TEMP_TL(ii,getVariableId("Z")) *  bfn   (ii,kk)
                dZdx1   = dZdx1   + TEMP_TL(ii,getVariableId("Z")) * dbfndx1(ii,kk)
                dZdx2   = dZdx2   + TEMP_TL(ii,getVariableId("Z")) * dbfndx2(ii,kk)
                
                Ksi     =  Ksi    + Ksi_(NM(ii)) * bfn    (ii,kk)
                dKsidx1 = dKsidx1 + Ksi_(NM(ii)) * dbfndx1(ii,kk)
                dKsidx2 = dKsidx2 + Ksi_(NM(ii)) * dbfndx2(ii,kk)

                Eta     =  Eta    + Eta_(NM(ii)) * bfn    (ii,kk)
                dEtadx1 = dEtadx1 + Eta_(NM(ii)) * dbfndx1(ii,kk)
                dEtadx2 = dEtadx2 + Eta_(NM(ii)) * dbfndx2(ii,kk)
            end do
            !*********************************************************************
            ! Calculate the Jacobian of Transformation
            ! The Jacobian of tranformation in referred on computational domain
            !*********************************************************************
            JacC   = dEtadx2 * dKsidx1 - dEtadx1 * dKsidx2
            dx1dKsi=   dEtadx2/JacC
            dx1dEta= - dKsidx2/JacC
            dx2dKsi= - dEtadx1/JacC
            dx2dEta=   dKsidx1/JacC

            !*********************************************************************
            ! Calculate the normal vectors with respect to the face of the 
            ! triangle
            !*********************************************************************

            Select Case(ned)
            Case(1); n_eta = -  dKsidx1          ; n_ksi =    dEtadx1
            Case(2); n_eta = -( dKsidx1-dKsidx2) ; n_ksi =  ( dEtadx1-dEtadx2)
            Case(3); n_eta = -(-dKsidx2)         ; n_ksi =  (-dEtadx2)  ! face 3 has reverse numbering
            Case Default
                Print*, "[Error] : Stresses. Wrong Value of face. Possible values 1,2,3."
            End Select

            dS    = dsqrt(n_ksi**2 + n_eta**2)
            n_eta = n_eta/dS
            n_ksi = n_ksi/dS
                        
            !*********************************************************************
            ! Calculate Derivatives In the Computational Space
            !*********************************************************************
            dRdEta = dRdx1 * dx1dEta + dRdx2 * dx2dEta
            dRdKsi = dRdx1 * dx1dKsi + dRdx2 * dx2dKsi

            dZdEta = dZdx1 * dx1dEta + dZdx2 * dx2dEta
            dZdKsi = dZdx1 * dx1dKsi + dZdx2 * dx2dKsi

            !*********************************************************************
            ! Define Scale Factor
            !*********************************************************************
            QEta = dsqrt(dZdEta**2 + dRdEta**2)

            !*********************************************************************
            ! Iterate over weighting Functions
            !*********************************************************************
            LOOP_RESIDUALS_f:DO IW = 1, NBF_2d

                BIFN    =   bfn  (iw,kk)
                DBIEta  = dbfndx1(iw,kk) * dx1dEta + dbfndx2(iw,kk) * dx2dEta
                DBIKsi  = dbfndx1(iw,kk) * dx1dKsi + dbfndx2(iw,kk) * dx2dKsi
                
                TERM_RES = 0.D0 
                TERM_RES(getVariableId("R")) = e_bnd*DLOG(QEta)*DBIEta * dS

            ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f) * WO_1d(KK)
            ENDDO LOOP_RESIDUALS_f
        ENDDO LOOP_GAUSS
                
        !*********************************************************************
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !*********************************************************************
        if ( STORE ) then
            NM = NM_f(NELEM,1:NBF_2d)
            CALL MATRIX_STORAGE_RESIDUAL( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
        end if

        !*********************************************************************
        ! Deallocate The local arrays
        !*********************************************************************
        deallocate( bfn )
        deallocate( dbfndx1 )
        deallocate( dbfndx2 )

    End Subroutine Y_EQUIDISTRIBUTION_RESIDUAL_f

! ********************************************************************

    SUBROUTINE Theta_EQUIDISTRIBUTION_RESIDUAL_f( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        USE CONTINUATION_MODULE, only: INCREMENT
        USE PHYSICAL_MODULE
        USE ELEMENTS_MODULE,         only: NBF_2d, NEL_2d, NEQ_f, NUNKNOWNS_f
        USE ENUMERATION_MODULE,      only: NM_MESH, NM_f, NME_MESH 
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, BFN_E, DFDC_E,&
                                        DFDE_E, &
                                                                            getBasisFunctionsAtFace, &
                                                                            getNormalVectorAtFace
        USE GLOBAL_ARRAYS_MODULE,    only: TLo, TLb
        USE FLOW_ARRAYS_MODULE,      only: B_f
        USE MESH_MODULE,             only: Xm, Ym, EPS_MESH
        USE MESH_MODULE,             only: Ksi_ => Xm, Eta_ => Ym
        USE LOG_REPRESENTATION

        IMPLICIT NONE
    
        !  ARGUMENTS
        INTEGER,                           INTENT(IN)  :: NELEM, NED
        REAL(8), DIMENSION(NBF_2d, NEQ_f), INTENT(IN)  :: TEMP_TL
        REAL(8), DIMENSION(NBF_2d, NEQ_f), INTENT(OUT) :: TEMP_RES
        LOGICAL,                           INTENT(IN)  :: STORE
       
        !**************************************************************************
            !  LOCAL VARIABLES
        !**************************************************************************
            ! FEM variables
        Real(8) :: R  , dRdx1  , dRdx2
        Real(8) :: Z  , dZdx1  , dZdx2
        Real(8) :: Ksi, dKsidx1, dKsidx2
        Real(8) :: Eta, dEtadx1, dEtadx2
        ! Jacobian Of Transformation
        Real(8) :: JacC
        Real(8) :: dx1dEta
        Real(8) :: dx1dKsi
        Real(8) :: dx2dEta
        Real(8) :: dx2dKsi
        ! Basis Functions
        Real(8), Dimension(:,:), Allocatable ::  bfn
        Real(8), Dimension(:,:), Allocatable :: dbfndx1
        Real(8), Dimension(:,:), Allocatable :: dbfndx2
        Real(8), Dimension(NBF_2d,NGAUSS_1d)           :: dbfndS 

            ! Normal and Tangent Vectors
        Real(8) :: n_ksi, n_eta
        Real(8) :: t_ksi, t_eta
        Real(8) :: n_x1 , n_x2
        Real(8) :: t_x1 , t_x2
        Real(8) :: dS
        

        INTEGER :: KK, II, JJ, LL, IW, JW, I, J, INOD, IEQ, JNOD, JEQ, ISTEP, IMOD
        INTEGER :: IROW, JCOL, N1, N2

        REAL(8)  :: WET
        REAL(8)  :: X, Xo, Xb, dXdC, dXdE, dXdY0, dXdX0, dXdt, dX0dY, dX0dX, QX0, SX0
        REAL(8)  :: Y, Yo, Yb, dYdC, dYdE, dYdY0, dYdX0, dYdt, dY0dY, dY0dX, QY0, SY0
        REAL(8)  :: X0, dX0dC, dX0dE
        REAL(8)  :: Y0, dY0dC, dY0dE
        REAL(8)  :: CJAC, AJAC, JacT, dL, nx, ny, nl, CJAC0, AJAC0, dL0
        real(8)  :: DThetaDT
        REAL(8)  :: BIFN, DBIX,  DBIY
        REAL(8)  ::       DBIX0, DBIY0
        real(8)  :: w1, w2
    
        REAL(8)  :: U1, U2, U3
        REAL(8)  :: TrC
       
       
    
        INTEGER, DIMENSION(NBF_2d) :: NM 
        REAL(8), DIMENSION(NBF_2d) :: DFDX,   DFDY
        REAL(8), DIMENSION(NBF_2d) :: DFDX0,  DFDY0
        REAL(8), DIMENSION(NEQ_f)  :: TERM_RES
   
        Real(8), Dimension(:,:), Allocatable :: dfdc, dfde  
        Real(8), Dimension(:)  , Allocatable :: X_loc
        Real(8), Dimension(:)  , Allocatable :: Y_loc

        Real(8), Dimension(NBF_2d,NGAUSS_1d) ::  DFDL0
        REAL(8)                              :: dKSIdTHETA, dETAdTHETA, THETA
        real(8)                              :: Scale_Factor
        REAL(8)                              :: dTHETAdKSI, dTHETAdETA, dQdtheta


        w1 = 1.d0
        w2 = 2.d0 - w1
        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !--------------------------------------------------------------------- 
        call copyArrayToLocalValues(Xm, nm_mesh(nelem,:), X_loc )
        call copyArrayToLocalValues(Ym, nm_mesh(nelem,:), Y_loc )
        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        NM = NM_MESH(NELEM,:)
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)
        call getBasisFunctionsAtFace(ned, bfn, dbfndx1, dbfndx2)
        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0


        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------
        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
                R   = 0.d0; dRdx1   = 0.d0; dRdx2   = 0.d0
                Z   = 0.d0; dZdx1   = 0.d0; dZdx2   = 0.d0
                Ksi = 0.d0; dKsidx1 = 0.d0; dKsidx2 = 0.d0  
                Eta = 0.d0; dEtadx1 = 0.d0; dEtadx2 = 0.d0

                do ii = 1, nbf_2d
                    R      =  R       + TEMP_TL(ii,getVariableId("R")) *  bfn   (ii,kk)
                    dRdx1   = dRdx1   + TEMP_TL(ii,getVariableId("R")) * dbfndx1(ii,kk)
                    dRdx2   = dRdx2   + TEMP_TL(ii,getVariableId("R")) * dbfndx2(ii,kk)

                    Z      =  Z       + TEMP_TL(ii,getVariableId("Z")) *  bfn   (ii,kk)
                    dZdx1   = dZdx1   + TEMP_TL(ii,getVariableId("Z")) * dbfndx1(ii,kk)
                    dZdx2   = dZdx2   + TEMP_TL(ii,getVariableId("Z")) * dbfndx2(ii,kk)

                    Ksi    =  Ksi     + Ksi_(NM(ii))                   *  bfn   (ii,kk)
                    dKsidx1 = dKsidx1 + Ksi_(NM(ii))                   * dbfndx1(ii,kk)
                    dKsidx2 = dKsidx2 + Ksi_(NM(ii))                   * dbfndx2(ii,kk)

                    Eta    =  Eta     + Eta_(NM(ii))                   *  bfn   (ii,kk)
                    dEtadx1 = dEtadx1 + Eta_(NM(ii))                   * dbfndx1(ii,kk)
                    dEtadx2 = dEtadx2 + Eta_(NM(ii))                   * dbfndx2(ii,kk)
                end do
                !*********************************************************************
                ! Calculate the Jacobian of Transformation
                ! The Jacobian of tranformation in referred on computational domain
                !*********************************************************************
                JacC   = dEtadx2 * dKsidx1 - dEtadx1 * dKsidx2
                dx1dKsi=   dEtadx2/JacC
                dx1dEta= - dKsidx2/JacC
                dx2dKsi= - dEtadx1/JacC
                dx2dEta=   dKsidx1/JacC

            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
                CALL BASIS_2d&
                ( KK, TEMP_TL(:,getVariableId("Z")) , TEMP_TL(:,getVariableId("R")) , BFN, DFDC, DFDE, X, dXdC, dXdE, Y, dYdC, dYdE,&
                CJAC, AJAC, DFDX, DFDY, NGAUSS_1d )

                CALL BASIS_2d&
                ( KK, X_loc, Y_loc, BFN, DFDC, DFDE, X0, dX0dC, dX0dE, Y0, dY0dC, dY0dE,&
                CJAC0, AJAC0, DFDX0, DFDY0, NGAUSS_1d )

                !*********************************************************************
                ! Calculate the normal and tangent vectors of the parent element
                !*********************************************************************
                Select Case(ned)
                !---------------------------------------------------------------------
                Case(1); n_x1 =  0.d0            ; n_x2 = -1.d0
                        t_x1 = +1.d0            ; t_x2 =  0.d0
                !---------------------------------------------------------------------
                Case(2); n_x1 =  1.d0/dsqrt(2.d0); n_x2 =  1.d0/dsqrt(2.d0)
                        t_x1 = -1.d0/dsqrt(2.d0); t_x2 =  1.d0/dsqrt(2.d0)
                !---------------------------------------------------------------------
                Case(3); n_x1 =-(-1.d0)          ; n_x2 =   0.d0
                        t_x1 =   0.d0           ; t_x2 =-(-1.d0) ! The minus is for reverse numbering
                End Select
                !*********************************************************************
                ! Calculate the normal vectors with respect to the face of the 
                ! triangle
                !*********************************************************************
                Select Case(ned)
                Case(1); n_eta = -  dKsidx1          ; n_ksi =    dEtadx1
                Case(2); n_eta = -( dKsidx1-dKsidx2) ; n_ksi =  ( dEtadx1-dEtadx2)
                Case(3); n_eta = -(-dKsidx2)         ; n_ksi =  (-dEtadx2)  ! face 3 has reverse numbering
                Case Default
                    Print*, "[Error] : Stresses. Wrong Value of face. Possible values 1,2,3."
                End Select

                dS    = dsqrt(n_ksi**2 + n_eta**2)
                n_eta = n_eta/dS
                n_ksi = n_ksi/dS
     
                t_ksi = n_eta
                t_eta = n_ksi   
            ! DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n
            SELECT CASE(NED)
      
                CASE(1)
                    dL           =   DSQRT(w1*DXDC**2+w2*DYDC**2)
                    dL0          =   DSQRT(DX0DC**2+DY0DC**2)
                    DFDL0        =   DFDC/dL0
                    Scale_Factor = DSQRT(w1*DX0DC**2+w2*DY0DC**2)
                 
             
                CASE(3)
                    dL           =    DSQRT(w1*DXDE**2+w2*DYDE**2)
                    dL0          =    DSQRT(DX0DE**2+DY0DE**2)
                    DFDL0        =    DFDE/dL0
                    Scale_Factor = DSQRT(w1*DX0DE**2+w2*DY0DE**2) 
               
             
                CASE(2)
                    dL           =   DSQRT(w1*(DXDC-DXDE)**2+w2*(DYDC-DYDE)**2)
                    dL0          =   DSQRT((DX0DC-DX0DE)**2+(DY0DC-DY0DE)**2)
                    DFDL0        =   (DFDC-DFDE)/dL0
                    Scale_Factor = DSQRT(w1*(DX0DC-DX0DE)**2+w2*(DY0DC-DY0DE)**2) 
            END SELECT

            WET = WO_1d(KK)*dL0


            !---------------------------------------------------------------------
            !    CALCULATE DEPENDENT VARIABLE AND PARTIAL DERIVATIVES AT
            !    THE GAUSSIAN INTEGRATION POINTS
            !---------------------------------------------------------------------
            dXdX0  = 0.D0
            dYdX0  = 0.D0
            dXdY0  = 0.D0
            dYdY0  = 0.D0
          
            DO II = 1, NBF_2d
          
                JJ = NM_MESH(NELEM,II)
            
                dXdX0 = dXdX0 + TEMP_TL(II,getVariableId("Z")) * DFDX0(II)
            
                dYdX0 = dYdX0 + TEMP_TL(II,getVariableId("R")) * DFDX0(II)
            
                dXdY0 = dXdY0 + TEMP_TL(II,getVariableId("Z")) * DFDY0(II)
            
                dYdY0 = dYdY0 + TEMP_TL(II,getVariableId("R")) * DFDY0(II) 
            ENDDO

            !-----------------------------------------------------------------------
            !     DEFINE SCALE FACTOR
            !-----------------------------------------------------------------------
        
                dQdtheta = DLOG(dL/Scale_Factor)
        
        
            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------
            LOOP_RESIDUALS_f:DO IW = 1, NBF_2d
    
                BIFN  = BFN(IW,KK)
                DBIX0 = DFDX0(IW)
                DBIY0 = DFDY0(IW)
    
                TERM_RES = 0.D0
    
                TERM_RES(getVariableId("R")) = e_bnd*dQdtheta*DFDL0(IW,KK)
                ! TERM_RES(getVariableId("Z")) = e_bnd*dQdtheta*DFDL0(IW,KK)
    
                ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)*WET 
            
            ENDDO LOOP_RESIDUALS_f
      
        ENDDO LOOP_GAUSS
    

        !---------------------------------------------------------------------
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !---------------------------------------------------------------------
        IF ( STORE ) THEN
    
        NM = NM_f(NELEM,1:NBF_2d)
      
        CALL MATRIX_STORAGE_RESIDUAL&
        ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
      
        ENDIF
    
    END SUBROUTINE Theta_EQUIDISTRIBUTION_RESIDUAL_f

! ********************************************************************


! SUBROUTINE OUTFLOW_RESIDUAL_f( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
!       Use VariableMapping
!       USE CONTINUATION_MODULE,     only: INCREMENT
!       USE PHYSICAL_MODULE
!       USE ELEMENTS_MODULE,         only: NBF_2d, NEL_2d, NEQ_f, NUNKNOWNS_f
!       USE GAUSS_MODULE,            only: WO_1d, NGAUSS_1d, BFN_E, DFDC_E,&
!                                          DFDE_E, &
!                                                  getBasisFunctionsAtFace
!       USE ENUMERATION_MODULE,      only: NM_MESH, NM_f
!       USE GLOBAL_ARRAYS_MODULE,    only: TLo, TLb, DpL
!       USE FLOW_ARRAYS_MODULE,      only: B_f
!       USE MESH_MODULE,             only: Xm, Ym
!       USE sr_representation,       only: set_S_get_Stress
!       USE TIME_INTEGRATION,        only: Dt
!       IMPLICIT NONE
      
!     !  ARGUMENTS
!       INTEGER,                           INTENT(IN)  :: NELEM, NED
!       REAL(8), DIMENSION(NBF_2d, NEQ_f), INTENT(IN)  :: TEMP_TL
!       REAL(8), DIMENSION(NBF_2d, NEQ_f), INTENT(OUT) :: TEMP_RES
!       LOGICAL,                           INTENT(IN)  :: STORE
                  
!     !  LOCAL VARIABLES
!       INTEGER :: KK, II, JJ, LL, IW, JW, I, J, INOD, IEQ, JNOD, JEQ, ISTEP, IMOD
!       INTEGER :: IROW, JCOL

!       REAL(8)  :: WET
!       REAL(8)  :: Z, dZdC, dZdE
!       REAL(8)  :: R, dRdC, dRdE
!       REAL(8)  :: CJAC, AJAC, JacT, dL, nr, nz, nl
  
!       REAL(8)  :: BIFN, DBIR, DBIZ
!       REAL(8)  :: BJFN, DBJX, DBJY, DBJZ
      
!       REAL(8)  :: P
      
    
!       REAL(8)  :: Prr, Prz, Pzz
!       Real(8)  :: Trr, Trz, Tzz, Ttt
!       Real(8)  :: Srr, Srz, Szz, Stt

!       Real(8), Dimension(3,3)     :: S_tensor_o, S_tensor, Stress_tensor_o, Stress_tensor

      
!       INTEGER, DIMENSION(NBF_2d) :: NM 
!       REAL(8), DIMENSION(NBF_2d) :: X_loc, Y_loc
!       REAL(8), DIMENSION(NBF_2d) :: DFDR,  DFDZ
!       REAL(8), DIMENSION(NEQ_f)  :: TERM_RES
      
!       ! REAL(8), DIMENSION(NBF_2d,NGAUSS_1d) :: BFN, DFDC, DFDE
!       Real(8), Dimension(:,:), Allocatable ::  bfn
!       Real(8), Dimension(:,:), Allocatable :: dfdc, dfde 
!     !  VARIABLES for the EVP HB MODEL
      
      

!     !---------------------------------------------------------------------
!     !  COPY X VECTOR TO LOCAL VECTOR
!     !---------------------------------------------------------------------
!       DO II = 1, NBF_2d
!               JJ = NM_MESH(NELEM,II)
      
!               X_loc(II) = Xm(JJ)
!               Y_loc(II) = Ym(JJ)
!       ENDDO
      
!       call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)

      

!       TEMP_RES = 0.D0


   
!       LOOP_GAUSS: DO KK = 1, NGAUSS_1d


!           CALL BASIS_2d&
!           ( KK, TEMP_TL(:,getVariableId("Z")), TEMP_TL(:,getVariableId("R")), BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
!                   CJAC, AJAC, DFDZ, DFDR, NGAUSS_1d )

!               SELECT CASE(NED)
              
!                       CASE(1)
!                               dL = DSQRT(DRDC**2+DZDC**2)
!                               nz =   DRDC/dL
!                               nr = - DZDC/dL

              
!                       CASE(2)
!                               dL = DSQRT((DRDC-DRDE)**2+(DZDC-DZDE)**2)
!                               nr = - (DZDC-DZDE)/dL
!                               nz =   (DRDC-DRDE)/dL
              
!                       CASE(3)
                      
!                       ! I mutually changed case 2 and case 3 and it worked
!                       ! in the future i need to call subroutines to do the procedure automatically
                              
!                               dL = DSQRT(DRDE**2+DZDE**2)
!                               nr = + DZDE/dL
!                               nz = - DRDE/dL
!               END SELECT
!               ! dL = 1.d0
!               ! nr = 0.d0
!               ! nz = 1.d0
!                ! if (kk .eq. 1) then
!                ! print*, ned
!                ! print*, nr, nz
!                ! pause
!                ! endif
                
              
!               WET = WO_1d(KK)*dL
              
!         !---------------------------------------------------------------------
!         !    CALCULATE DEPENDENT VARIABLE AND PARTIAL DERIVATIVES AT
!         !    THE GAUSSIAN INTEGRATION POINTS
!         !---------------------------------------------------------------------

!               P   = 0.D0  ;    

!               Srr = 0.D0  ;
!               Srz = 0.D0  ;
!               Szz = 0.D0  ;
!               Stt = 0.D0  ;
              
!               DO II = 1, NBF_2d
              
!                       JJ = NM_MESH(NELEM,II)
                      
!                       P     = P     + TEMP_TL(II,getVariableId("P"))*BFN(II,KK)
                      
!                       Srr    = Srr  + TEMP_TL(II,getVariableId("Srr"))*BFN(II,KK)
                      
!                       Srz    = Srz  + TEMP_TL(II,getVariableId("Srz"))*BFN(II,KK)

!                       Szz    = Szz  + TEMP_TL(II,getVariableId("Szz"))*BFN(II,KK)
!               ENDDO

!             S_tensor = 0.d0
!             Stt      = 0.d0 ! Stt is not present in the current weak form so i dont care to retrieve its value
!                             ! i just give a nominal value 0 to calculate the matrix Stress
!             S_tensor(1,1) = Srr ; S_tensor(1,2) = Srz 
!             S_tensor(2,1) = Srz ; S_tensor(2,2) = Szz 
!             S_tensor(3,3) = Stt

!             call set_S_get_Stress( S_tensor  , Stress_tensor   )
!             Trr = Stress_tensor(1,1)
!             Trz = Stress_tensor(1,2)
!             Tzz = Stress_tensor(2,2)
!         !---------------------------------------------------------------------     
!         !    DEFINE EXTRA STRESS TENSOR
!         !---------------------------------------------------------------------
!               Prr = - P + Trr
!               Prz =       Trz
!               Pzz = - P + Tzz
    
!               LOOP_RESIDUALS_f:DO IW = 1, NBF_2d

!                       BIFN = BFN(IW,KK)
!                       DBIR = DFDR(IW)
!                       DBIZ = DFDZ(IW)

!                       TERM_RES     = 0.D0

!                       TERM_RES(getVariableId("Vr"))  = - ( nr*(Prr  ) + nz*(Prz  ))*BIFN*R

!                       TERM_RES(getVariableId("Vz"))  = - ( nr*(Prz )  + nz*(Pzz  ))*BIFN*R


!                       TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)*WET 
                      
!               ENDDO LOOP_RESIDUALS_f
              

!       ENDDO LOOP_GAUSS
      

   
!       IF ( STORE ) THEN
      
!               NM = NM_f(NELEM,1:NBF_2d)
              
!               CALL MATRIX_STORAGE_RESIDUAL&
!               ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
              
!       ENDIF

!     END SUBROUTINE OUTFLOW_RESIDUAL_f



end module Boundary_Equations

! ********************************************************************

    Subroutine Y_EQUIDISTRIBUTION_RESIDUAL_f( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        USE PHYSICAL_MODULE
        USE ELEMENTS_MODULE,     only: NBF_2d,  NEQ_f, NUNKNOWNS_f
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

            dS    = sqrt(n_ksi**2 + n_eta**2)
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
            QEta = sqrt(dZdEta**2 + dRdEta**2)

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

    Subroutine Kinematic_mass( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
                                                        getBasisFunctionsAtFace, &
                                                        getNormalVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use TIME_INTEGRATION,        Only: Dt
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
        Real(8)                                        :: C , dCdx1 , dCdx2, dCdZ, dCdR

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

        ! Surface Arc Length
        Real(8)                                        :: dS
    

        Integer                                        :: KK, II,  IW
    
        Real(8)                                        :: Rb
        Real(8)                                        :: Zb
        Real(8)                                        :: BIFN, DBIR, DBIZ, SBFN
          
        Real(8)                                        :: Zo, Ro
        Real(8)                                        :: dRdt, dZdt
    
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
            C  = 0.d0; dCdx1  = 0.d0;  dCdx2 = 0.d0
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

                C    =  C    + TEMP_TL(ii, getVariableId("C")) *  bfn   (ii,kk)
                dCdx1= dCdx1 + TEMP_TL(ii, getVariableId("C")) * dbfndx1(ii,kk)
                dCdx2= dCdx2 + TEMP_TL(ii, getVariableId("C")) * dbfndx2(ii,kk) 
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
            dCdZ = dCdx1 * dx1dZ + dCdx2 * dx2dZ
            dCdR = dCdx1 * dx1dR + dCdx2 * dx2dR
            ! Calculate pspg values
            !*********************************************************************
            ! Uelem = abs( Vr * tr + Vz * tz  )
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
            
                TERM_RES(getVariableId("Z")) = SBFN * ( ( nR * (Vr-dRdt) + nZ * (Vz-dZdt) )*C - (nR*dCdR + nZ*dCdZ)/PeN  )* R * dS

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

    End Subroutine Kinematic_mass

! ********************************************************************

    Subroutine Kinematic_mass_gasInterface( NELEM, NED, TEMP_TL, TEMP_RES, STORE, bmol, bvolume, bvelocity )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
                                                        getBasisFunctionsAtFace, &
                                                        getNormalVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use TIME_INTEGRATION,        Only: Dt 
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer,                           Intent(In)  :: NELEM, NED
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
        Logical,                           Intent(In)  :: STORE
        Real(8),                           intent(in)  :: bmol 
        Real(8),                           intent(in)  :: bvolume 
        Real(8),                           intent(in)  :: bvelocity 
        
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        ! FEM variables and their derivatives
        Real(8)                                        :: R ,  dRdx1,  dRdx2 
        Real(8)                                        :: Z ,  dZdx1,  dZdx2    
        Real(8)                                        :: Vr, dVrdx1, dVrdx2
        Real(8)                                        :: Vz, dVzdx1, dVzdx2
        Real(8)                                        :: C , dCdx1 , dCdx2, dCdZ, dCdR

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

        ! Surface Arc Length
        Real(8)                                        :: dS
    

        Integer                                        :: KK, II,  IW
    
        Real(8)                                        :: Rb
        Real(8)                                        :: Zb
        Real(8)                                        :: BIFN, DBIR, DBIZ, SBFN
          
        Real(8)                                        :: Zo, Ro
        Real(8)                                        :: dRdt, dZdt
    
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
            C  = 0.d0; dCdx1  = 0.d0;  dCdx2 = 0.d0
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

                C    =  C    + TEMP_TL(ii, getVariableId("C")) *  bfn   (ii,kk)
                dCdx1= dCdx1 + TEMP_TL(ii, getVariableId("C")) * dbfndx1(ii,kk)
                dCdx2= dCdx2 + TEMP_TL(ii, getVariableId("C")) * dbfndx2(ii,kk) 
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
            
            !*********************************************************************
            ! Calculate the Solution of the previous time steps
            ! to take into account the time variation
            !*********************************************************************
            Ro = 0.d0; Zo = 0.d0 

            do ii = 1, nbf_2d
                Ro = Ro + TLo(NM(ii), getVariableId("R")) * bfn(ii,kk)
                Zo = Zo + TLo(NM(ii), getVariableId("Z")) * bfn(ii,kk)
            end do
            !*********************************************************************
            ! Calculate the time derivatives of the node
            !*********************************************************************

            dRdt = (R - Ro)/Dt
            dZdt = (Z - Zo)/Dt

            !*********************************************************************
            dCdZ = dCdx1 * dx1dZ + dCdx2 * dx2dZ
            dCdR = dCdx1 * dx1dR + dCdx2 * dx2dR
            ! Calculate pspg values
            !*********************************************************************
            ! Uelem = abs( Vr * tr + Vz * tz  )
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
            
                TERM_RES(getVariableId("Z")) = SBFN * (  +( nZ * (bvelocity - dZdt ) + nR * (- dRdt ) )*bmol/bvolume &
                                                            - ( ( nR * (Vr-dRdt) + nZ * (Vz-dZdt) )*C  &
                                                                + (nR*dCdR + nZ*dCdZ)/PeN )  )* R * dS

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

    End Subroutine Kinematic_mass_gasInterface

! ********************************************************************

    Subroutine Henry( NELEM, NED, TEMP_TL, TEMP_RES, STORE, gVar )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
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
        Real(8)                              :: R, dRdx1, dRdx2, C
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
        Real(8)                              :: BIFN
        Integer :: KK, II, IW 

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

            R = 0.d0      ;     Z = 0.d0     ;     C = 0.d0
            dRdx1 = 0.d0  ; dRdx2 = 0.d0  ; dZdx1 = 0.d0  ; dZdx2 = 0.d0

            do ii = 1, nbf_2d
                R     =  R    + TEMP_TL(ii, getVariableId("R")) *  bfn   (ii,kk)
                dRdx1 = dRdx1 + TEMP_TL(ii, getVariableId("R")) * dbfndx1(ii,kk)
                dRdx2 = dRdx2 + TEMP_TL(ii, getVariableId("R")) * dbfndx2(ii,kk)

                Z     =  Z    + TEMP_TL(ii, getVariableId("Z")) *  bfn   (ii,kk)
                dZdx1 = dZdx1 + TEMP_TL(ii, getVariableId("Z")) * dbfndx1(ii,kk)
                dZdx2 = dZdx2 + TEMP_TL(ii, getVariableId("Z")) * dbfndx2(ii,kk) 
                
                C     =  C    + TEMP_TL(ii, getVariableId("C")) *  bfn   (ii,kk)
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
            

            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------

            loop_residuals_f:DO IW = 1, NBF_2d
        
                    BIFN =  bfn   (iw,kk)
        
                    TERM_RES     = 0.D0
                    TERM_RES(getVariableId("C"))  = (C - KoN*gVar)*BIFN*R
        
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

    end Subroutine Henry

! *****************************************************************

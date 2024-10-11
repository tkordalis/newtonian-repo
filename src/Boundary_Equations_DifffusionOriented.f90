Module Boundary_EquationsDO
    Use ArrayTools, Only: copyArrayToLocalValues
    use storage,    only: MATRIX_STORAGE_RESIDUAL
    use basis_calculations, only: BASIS_2d
    use check_for_floating_point_exceptions

    Contains

    Subroutine Stresses( NELEM, NED, TEMP_TL, TEMP_RES, STORE, gVar )
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
            !     print*, 'Bubble pressure=', gvar
            !     endif
            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------

            loop_residuals_f:DO IW = 1, NBF_2d
        
                    BIFN =  bfn   (iw,kk)
                    DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                    DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
        
                    ! -n * T = + P_bubble n + 1/Bo * (-nn)\nabla \cdot \phi
        
                    TERM_RES     = 0.D0
                    TERM_RES(getVariableId("Vr"))  = nr*gVar*BIFN*R + (1.d0/BoN)*((1.D0-nr*nr)*DBIR + BIFN/R + (    -nr*nz)*DBIZ)*R
                    TERM_RES(getVariableId("Vz"))  = nz*gVar*BIFN*R + (1.d0/BoN)*((    -nr*nz)*DBIR          + (1.D0-nz*nz)*DBIZ)*R
        
                    !      FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_1d(KK) * dS
                              
            end do loop_residuals_f
        end do LOOP_GAUSS
      
                ! print*, 'temp_res=', temp_res

        !---------------------------------------------------------------------
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !---------------------------------------------------------------------
        if ( STORE ) then 
            NM = NM_f(NELEM,1:NBF_2d)
            call MATRIX_STORAGE_RESIDUAL ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
        end if               

    end Subroutine Stresses

! ********************************************************************
    
    Subroutine Kinematic_mass_gasInterf_z( NELEM, NED, TEMP_TL, TEMP_RES, STORE, bmol, bvolume, bvelocity )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
                                                        getBasisFunctionsAtFace, &
                                                        getNormalVectorAtFace, &
                                                        getTangentVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use TIME_INTEGRATION,        Only: Dt , increment, time
        use basis_calculations,      only: BASIS_2d
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
        Real(8), Dimension(:,:), Allocatable :: dfdc, dfde  

        ! Jacobian of Transformation and the reverse derivatives
        Real(8)                                        :: JacT
        Real(8)                                        :: dx1dR
        Real(8)                                        :: dx2dR
        Real(8)                                        :: dx1dZ
        Real(8)                                        :: dx2dZ
        ! Normal Vector Components
        Real(8)                                        :: nr , tr 
        Real(8)                                        :: nz , tz 

        ! Surface Arc Length
        Real(8)                                        :: dS, dL
    

        Integer                                        :: KK, II,  IW
    
        Real(8)                                        :: Rb
        Real(8)                                        :: Zb
        Real(8)                                        :: BIFN, DBIR, DBIZ, SBFN
          
        Real(8)                                        :: Zo, Ro, bconcentration
        Real(8)                                        :: dRdt, dZdt, Uelem, tsupg
        Real(8)                                        :: X, dXdC, dXdE, Y, dYdC, dYdE, CJAC, AJAC
        REAL(8), DIMENSION(NBF_2d)                     :: DFDX,   DFDY
        Integer, Dimension(NBF_2d)                     :: NM 
        Real(8), Dimension(NEQ_f)                      :: TERM_RES
        Real(8), Dimension(NBF_2d,NGAUSS_1d)           ::  DFDL

      

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
        bconcentration = 0.d0!bmol/bvolume
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


            CALL BASIS_2d&
                ( KK, TEMP_TL(:,getVariableId("Z")) , TEMP_TL(:,getVariableId("R")) , BFN, DFDC, DFDE, X, dXdC, dXdE, Y, dYdC, dYdE, &
                CJAC, AJAC, DFDX, DFDY, NGAUSS_1d )

            !*********************************************************************
            ! Calculate the normal vectors with respect to the face of the 
            ! triangle
            !*********************************************************************

            call getNormalVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, nr, nz, dS, normalize = .true.)
            
            call getTangentVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, tr, tz, normalize = .true.)


             SELECT CASE(NED)
                CASE(1)
                    dL           =   sqrt(DXDC**2+DYDC**2)
                    DFDL         =   DFDC/dL
                CASE(3)
                    dL           =    sqrt(DXDE**2+DYDE**2)
                    DFDL         =    DFDE/dL
                CASE(2)
                    dL           =   sqrt((DXDC-DXDE)**2+(DYDC-DYDE)**2)
                    DFDL         =   (DFDC-DFDE)/dL
            END SELECT
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
            Uelem = abs( Vr * tr + Vz * tz  )
            ! Uelem = abs( (Vr-dRdt) * tr + (Vz-dZdt) * tz  )
            ! Uelem = abs( (-dRdt) * tr + (-dZdt) * tz  )
            ! Uelem = abs( bvelocity + 1.d-8 )
            tsupg = dS/(Uelem + dS/dt)*bmol/bvolume

            !*********************************************************************
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !*********************************************************************
            

            loop_residuals_f:do iw = 1, nbf_2d

                BIFN =  bfn   (iw,kk)
                DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
            
                SBFN         = BIFN + tsupg*( (Vr-dRdt)*tR + (Vz-dZdt)*tZ )*DFDL(IW,KK)
                ! SBFN         = BIFN + tsupg*( (- dRdt )*tR + (bvelocity - dZdt )*tZ )*DFDL(IW,KK)
                TERM_RES     = 0.D0
            
                TERM_RES(getVariableId("Z")) = SBFN*(- dZdt*(bconcentration-C) + (nz*bvelocity*bconcentration - (nZ*Vz + nR*Vr)*C + (nR*dCdR + nZ*dCdZ)/PeN )*nZ )* R * dS

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

    End Subroutine Kinematic_mass_gasInterf_z


    Subroutine Kinematic_mass_gasInterf_r( NELEM, NED, TEMP_TL, TEMP_RES, STORE, bmol, bvolume, bvelocity )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
                                                        getBasisFunctionsAtFace, &
                                                        getNormalVectorAtFace, &
                                                        getTangentVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use TIME_INTEGRATION,        Only: Dt , increment, time
        use basis_calculations,      only: BASIS_2d
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
        Real(8), Dimension(:,:), Allocatable :: dfdc, dfde  

        ! Jacobian of Transformation and the reverse derivatives
        Real(8)                                        :: JacT
        Real(8)                                        :: dx1dR
        Real(8)                                        :: dx2dR
        Real(8)                                        :: dx1dZ
        Real(8)                                        :: dx2dZ
        ! Normal Vector Components
        Real(8)                                        :: nr , tr 
        Real(8)                                        :: nz , tz 

        ! Surface Arc Length
        Real(8)                                        :: dS, dL
    

        Integer                                        :: KK, II,  IW
    
        Real(8)                                        :: Rb
        Real(8)                                        :: Zb
        Real(8)                                        :: BIFN, DBIR, DBIZ, SBFN
          
        Real(8)                                        :: Zo, Ro, bconcentration
        Real(8)                                        :: dRdt, dZdt, Uelem, tsupg
        Real(8)                                        :: X, dXdC, dXdE, Y, dYdC, dYdE, CJAC, AJAC
        REAL(8), DIMENSION(NBF_2d)                     :: DFDX,   DFDY
        Integer, Dimension(NBF_2d)                     :: NM 
        Real(8), Dimension(NEQ_f)                      :: TERM_RES
        Real(8), Dimension(NBF_2d,NGAUSS_1d)           ::  DFDL

      

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
        bconcentration = 0.d0!bmol/bvolume
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


            CALL BASIS_2d&
                ( KK, TEMP_TL(:,getVariableId("Z")) , TEMP_TL(:,getVariableId("R")) , BFN, DFDC, DFDE, X, dXdC, dXdE, Y, dYdC, dYdE, &
                CJAC, AJAC, DFDX, DFDY, NGAUSS_1d )

            !*********************************************************************
            ! Calculate the normal vectors with respect to the face of the 
            ! triangle
            !*********************************************************************

            call getNormalVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, nr, nz, dS, normalize = .true.)
            
            call getTangentVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, tr, tz, normalize = .true.)


             SELECT CASE(NED)
                CASE(1)
                    dL           =   sqrt(DXDC**2+DYDC**2)
                    DFDL         =   DFDC/dL
                CASE(3)
                    dL           =    sqrt(DXDE**2+DYDE**2)
                    DFDL         =    DFDE/dL
                CASE(2)
                    dL           =   sqrt((DXDC-DXDE)**2+(DYDC-DYDE)**2)
                    DFDL         =   (DFDC-DFDE)/dL
            END SELECT
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
            Uelem = abs( Vr * tr + Vz * tz  )
            ! Uelem = abs( (Vr-dRdt) * tr + (Vz-dZdt) * tz  )
            ! Uelem = abs( (-dRdt) * tr + (-dZdt) * tz  )
            ! Uelem = abs( bvelocity + 1.d-8 )
            tsupg = dS/(Uelem + dS/dt)*bmol/bvolume

            !*********************************************************************
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !*********************************************************************
            

            loop_residuals_f:do iw = 1, nbf_2d

                BIFN =  bfn   (iw,kk)
                DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
            
                SBFN         = BIFN + tsupg*( (Vr-dRdt)*tZ + (Vz-dZdt)*tR )*DFDL(IW,KK)
                ! SBFN         = BIFN + tsupg*( (Vr-dRdt)*tR + (Vz-dZdt)*tZ )*DFDL(IW,KK)
                ! SBFN         = BIFN + tsupg*( (- dRdt )*tR + (bvelocity - dZdt )*tZ )*DFDL(IW,KK)
                TERM_RES     = 0.D0
            
                TERM_RES(getVariableId("R")) = SBFN*(- dRdt*(bconcentration-C) + (nz*bvelocity*bconcentration - (nZ*Vz + nR*Vr)*C + (nR*dCdR + nZ*dCdZ)/PeN )*nR )* R * dS

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

    End Subroutine Kinematic_mass_gasInterf_r



    Subroutine Kinematic_mass_gasInterfaceSUPG( NELEM, NED, TEMP_TL, TEMP_RES, STORE, bmol, bvolume, bvelocity )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
                                                        getBasisFunctionsAtFace, &
                                                        getNormalVectorAtFace, &
                                                        getTangentVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use TIME_INTEGRATION,        Only: Dt , increment, time
        use basis_calculations,      only: BASIS_2d
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
        Real(8), Dimension(:,:), Allocatable :: dfdc, dfde  

        ! Jacobian of Transformation and the reverse derivatives
        Real(8)                                        :: JacT
        Real(8)                                        :: dx1dR
        Real(8)                                        :: dx2dR
        Real(8)                                        :: dx1dZ
        Real(8)                                        :: dx2dZ
        ! Normal Vector Components
        Real(8)                                        :: nr , tr 
        Real(8)                                        :: nz , tz 

        ! Surface Arc Length
        Real(8)                                        :: dS, dL
    

        Integer                                        :: KK, II,  IW
    
        Real(8)                                        :: Rb
        Real(8)                                        :: Zb
        Real(8)                                        :: BIFN, DBIR, DBIZ, SBFN
          
        Real(8)                                        :: Zo, Ro
        Real(8)                                        :: dRdt, dZdt, Uelem, tsupg
        Real(8)                                        :: X, dXdC, dXdE, Y, dYdC, dYdE, CJAC, AJAC
        REAL(8), DIMENSION(NBF_2d)                     :: DFDX,   DFDY
        Integer, Dimension(NBF_2d)                     :: NM 
        Real(8), Dimension(NEQ_f)                      :: TERM_RES
        Real(8), Dimension(NBF_2d,NGAUSS_1d)           ::  DFDL

      

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


            CALL BASIS_2d&
                ( KK, TEMP_TL(:,getVariableId("Z")) , TEMP_TL(:,getVariableId("R")) , BFN, DFDC, DFDE, X, dXdC, dXdE, Y, dYdC, dYdE, &
                CJAC, AJAC, DFDX, DFDY, NGAUSS_1d )

            !*********************************************************************
            ! Calculate the normal vectors with respect to the face of the 
            ! triangle
            !*********************************************************************

            call getNormalVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, nr, nz, dS, normalize = .true.)
            
            call getTangentVectorAtFace( [dZdx1, dZdx2, dRdx1, dRdx2] , &
                                         ned, tr, tz, normalize = .true.)


             SELECT CASE(NED)
                CASE(1)
                    dL           =   sqrt(DXDC**2+DYDC**2)
                    DFDL         =   DFDC/dL
                CASE(3)
                    dL           =    sqrt(DXDE**2+DYDE**2)
                    DFDL         =    DFDE/dL
                CASE(2)
                    dL           =   sqrt((DXDC-DXDE)**2+(DYDC-DYDE)**2)
                    DFDL         =   (DFDC-DFDE)/dL
            END SELECT
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
            Uelem = abs( Vr * tr + Vz * tz  )
            ! Uelem = abs( (Vr-dRdt) * tr + (Vz-dZdt) * tz  )
            ! Uelem = abs( (-dRdt) * tr + (-dZdt) * tz  )
            ! Uelem = abs( bvelocity + 1.d-8 )
            tsupg = dS/(Uelem + dS/dt)

            !*********************************************************************
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !*********************************************************************
            

            loop_residuals_f:do iw = 1, nbf_2d

                BIFN =  bfn   (iw,kk)
                DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
            
                SBFN         = BIFN !+ tsupg*( (Vr-dRdt)*tR + (Vz-dZdt)*tZ )*DFDL(IW,KK)
                ! SBFN         = BIFN + tsupg*( (- dRdt )*tR + (bvelocity - dZdt )*tZ )*DFDL(IW,KK)
                TERM_RES     = 0.D0
            
                ! TERM_RES(getVariableId("Z")) = ( SBFN *( -( nZ * (bvelocity - dZdt ) + nR * (- dRdt ) )*bmol/bvolume &
                !                                 + ( nZ * (Vz-dZdt) + nR * (Vr-dRdt) )*C ) &
                !                             - BIFN*(nR*dCdR + nZ*dCdZ)/PeN )* R * dS

                ! TERM_RES(getVariableId("Z")) = ( SBFN *( -( nZ * (bvelocity - dZdt ) + nR * (- dRdt ) )*max(0,increment-2)*bmol/bvolume &
                ! TERM_RES(getVariableId("Z")) = ( SBFN *( -( nZ * (bvelocity - dZdt ) + nR * (- dRdt ) )*(1.d0 -exp(-time))*bmol/bvolume &

                ! TERM_RES(getVariableId("Z")) = ( SBFN *( -( nZ * (bvelocity - dZdt ) + nR * (- dRdt ) )*bmol/bvolume &
                TERM_RES(getVariableId("Z")) = ( SBFN *( -( nZ * (bvelocity - dZdt ) + nR * (- dRdt ) )*bmol/bvolume &
                                                - ( nZ * (Vz-dZdt) + nR * (Vr-dRdt) )*C ) &
                                            + BIFN*(nR*dCdR + nZ*dCdZ)/PeN )* R * dS
                                                ! + ( nZ * (Vz-dZdt) + nR * (Vr-dRdt) )*(Cwater/Cchar+C) ) &
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

    End Subroutine Kinematic_mass_gasInterfaceSUPG

! ********************************************************************
    
    Subroutine Kinematic( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
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
            ! if ((KK .eq. 1) .and. (store)) then
            ! print*, Z, R
            ! print*, nz, nr 
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
            
                TERM_RES(getVariableId("Z")) = SBFN * ( nR * (Vr-dRdt) + nZ * (Vz-dZdt) ) * R * dS

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

! ********************************************************************

    Subroutine X_EQUIDISTRIBUTION_RESIDUAL_f( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
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

            dS = sqrt(n_ksi**2 + n_eta**2)
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
            QKsi = sqrt(dZdKsi**2 + dRdKsi**2)

            !*********************************************************************
            ! Iterate over weighting Functions
            !*********************************************************************
            LOOP_RESIDUALS_f:DO IW = 1, NBF_2d

                BIFN    =   bfn  (iw,kk)
                DBIEta  = dbfndx1(iw,kk) * dx1dEta + dbfndx2(iw,kk) * dx2dEta
                DBIKsi  = dbfndx1(iw,kk) * dx1dKsi + dbfndx2(iw,kk) * dx2dKsi
                
                TERM_RES = 0.D0 
                ! TERM_RES(getVariableId("Z")) = e_bnd*DLOG(QKsi)*DBIKsi * dS
                TERM_RES(getVariableId("Z")) = DLOG(QKsi)*DBIKsi * dS

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

    SUBROUTINE Theta_EQUIDISTRIBUTION_RESIDUAL_f( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        USE PHYSICAL_MODULE
        USE ELEMENTS_MODULE,         only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        USE ENUMERATION_MODULE,      only: NM_MESH, NM_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
                                                            getBasisFunctionsAtFace, &
                                                            getNormalVectorAtFace
        USE FLOW_ARRAYS_MODULE,      only: B_f
        USE MESH_MODULE,             only: Xm, Ym, EPS_MESH
        USE MESH_MODULE,             only: Ksi_ => Xm, Eta_ => Ym

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

            ! Normal and Tangent Vectors
        Real(8) :: n_ksi, n_eta
        Real(8) :: t_ksi, t_eta
        Real(8) :: n_x1 , n_x2
        Real(8) :: t_x1 , t_x2
        Real(8) :: dS
        

        INTEGER :: KK, II, JJ, IW

        REAL(8)  :: WET
        REAL(8)  :: X, dXdC, dXdE, dXdY0, dXdX0
        REAL(8)  :: Y, dYdC, dYdE, dYdY0, dYdX0
        REAL(8)  :: X0, dX0dC, dX0dE
        REAL(8)  :: Y0, dY0dC, dY0dE
        REAL(8)  :: CJAC, AJAC, dL, CJAC0, AJAC0, dL0
        REAL(8)  :: BIFN
        REAL(8)  ::       DBIX0, DBIY0
        real(8)  :: w1, w2
    
    
        INTEGER, DIMENSION(NBF_2d) :: NM 
        REAL(8), DIMENSION(NBF_2d) :: DFDX,   DFDY
        REAL(8), DIMENSION(NBF_2d) :: DFDX0,  DFDY0
        REAL(8), DIMENSION(NEQ_f)  :: TERM_RES
   
        Real(8), Dimension(:,:), Allocatable :: dfdc, dfde  
        Real(8), Dimension(:)  , Allocatable :: X_loc
        Real(8), Dimension(:)  , Allocatable :: Y_loc

        Real(8), Dimension(NBF_2d,NGAUSS_1d) ::  DFDL0
        real(8)                              :: Scale_Factor
        REAL(8)                              :: dQdtheta


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

            ! if (kk==1 .and. store) then 
            ! print*, 'theta_equid'
            ! print*, TEMP_TL(:,getVariableId("Z")) , TEMP_TL(:,getVariableId("R"))
            ! endif
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
                Case(2); n_x1 =  1.d0/sqrt(2.d0); n_x2 =  1.d0/sqrt(2.d0)
                        t_x1 = -1.d0/sqrt(2.d0); t_x2 =  1.d0/sqrt(2.d0)
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

                dS    = sqrt(n_ksi**2 + n_eta**2)
                n_eta = n_eta/dS
                n_ksi = n_ksi/dS
     
                t_ksi = n_eta
                t_eta = n_ksi   
            ! DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n
            SELECT CASE(NED)
      
                CASE(1)
                    dL           =   sqrt(w1*DXDC**2+w2*DYDC**2)
                    dL0          =   sqrt(DX0DC**2+DY0DC**2)
                    DFDL0        =   DFDC/dL0
                    Scale_Factor = sqrt(w1*DX0DC**2+w2*DY0DC**2)
                 
             
                CASE(3)
                    dL           =    sqrt(w1*DXDE**2+w2*DYDE**2)
                    dL0          =    sqrt(DX0DE**2+DY0DE**2)
                    DFDL0        =    DFDE/dL0
                    Scale_Factor = sqrt(w1*DX0DE**2+w2*DY0DE**2) 
               
             
                CASE(2)
                    dL           =   sqrt(w1*(DXDC-DXDE)**2+w2*(DYDC-DYDE)**2)
                    dL0          =   sqrt((DX0DC-DX0DE)**2+(DY0DC-DY0DE)**2)
                    DFDL0        =   (DFDC-DFDE)/dL0
                    Scale_Factor = sqrt(w1*(DX0DC-DX0DE)**2+w2*(DY0DC-DY0DE)**2) 
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
    
                ! TERM_RES(getVariableId("R")) = e_bnd*dQdtheta*DFDL0(IW,KK)
                TERM_RES(getVariableId("R")) = dQdtheta*DFDL0(IW,KK)
    
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

    Subroutine zeroConcentrationFlux( NELEM, NED, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
                                            getBasisFunctionsAtFace, &
                                            getNormalVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo
        Use TIME_INTEGRATION,        Only: Dt

        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer,                           Intent(In)  :: NELEM, NED
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
        Logical,                           Intent(In)  :: STORE

    
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
        Real(8)                              :: dS, Ro, Zo, dRdt, dZdt, Vr, Vz

        Integer, Dimension(NBF_2d)           :: NM 
        Real(8), Dimension(NEQ_f)            :: TERM_RES
        ! Basis Function 
        Real(8)                              :: BIFN, DBIR, DBIZ
        Real(8)                              :: C , dCdx1 , dCdx2, dCdZ, dCdR

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

            R  = 0.d0; dRdx1 = 0.d0; dRdx2 = 0.d0
            Z  = 0.d0; dZdx1 = 0.d0; dZdx2 = 0.d0
            Vr = 0.d0; Vz    = 0.d0
            C  = 0.d0; dCdx1 = 0.d0; dCdx2 = 0.d0
            do ii = 1, nbf_2d
                R     =  R    + TEMP_TL(ii, getVariableId("R")) *  bfn   (ii,kk)
                dRdx1 = dRdx1 + TEMP_TL(ii, getVariableId("R")) * dbfndx1(ii,kk)
                dRdx2 = dRdx2 + TEMP_TL(ii, getVariableId("R")) * dbfndx2(ii,kk)

                Z     =  Z    + TEMP_TL(ii, getVariableId("Z")) *  bfn   (ii,kk)
                dZdx1 = dZdx1 + TEMP_TL(ii, getVariableId("Z")) * dbfndx1(ii,kk)
                dZdx2 = dZdx2 + TEMP_TL(ii, getVariableId("Z")) * dbfndx2(ii,kk)    
                
                Vr     =  Vr    + TEMP_TL(ii, getVariableId("Vr")) *  bfn   (ii,kk)
                Vz     =  Vz    + TEMP_TL(ii, getVariableId("Vz")) *  bfn   (ii,kk)

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
            

            Ro = 0.d0; Zo = 0.d0 

            do ii = 1, nbf_2d
                Ro = Ro + TLo(NM(ii), getVariableId("R")) * bfn(ii,kk)
                Zo = Zo + TLo(NM(ii), getVariableId("Z")) * bfn(ii,kk)
            end do

            dRdt = (R - Ro)/Dt
            dZdt = (Z - Zo)/Dt

            dCdZ = dCdx1 * dx1dZ + dCdx2 * dx2dZ
            dCdR = dCdx1 * dx1dR + dCdx2 * dx2dR

            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------

            loop_residuals_f:DO IW = 1, NBF_2d
        
                    BIFN =  bfn   (iw,kk)
                    DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                    DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
        
        
                    TERM_RES     = 0.D0
                    TERM_RES(getVariableId("C"))  = ( BIFN * ( nR * (Vr-dRdt) + nZ * (Vz-dZdt) ) * C ) * R
        
                    !      FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_1d(KK)  * dS
                              
            end do loop_residuals_f
        end do LOOP_GAUSS
      
                ! print*, 'temp_res=', temp_res

        !---------------------------------------------------------------------
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !---------------------------------------------------------------------
        if ( STORE ) then 
            NM = NM_f(NELEM,1:NBF_2d)
            call MATRIX_STORAGE_RESIDUAL ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
        end if               

    end Subroutine zeroConcentrationFlux

! ********************************************************************

    Subroutine weakHenry( NELEM, NED, TEMP_TL, TEMP_RES, STORE, gVar )
        Use VariableMapping
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d,  NEQ_f, NUNKNOWNS_f
        Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d, &
                                            getBasisFunctionsAtFace, &
                                            getNormalVectorAtFace
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo
        Use TIME_INTEGRATION,        Only: Dt

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
        Real(8)                              :: dS, Ro, Zo, dRdt, dZdt, Vr, Vz

        Integer, Dimension(NBF_2d)           :: NM 
        Real(8), Dimension(NEQ_f)            :: TERM_RES
        ! Basis Function 
        Real(8)                              :: BIFN, DBIR, DBIZ
        Real(8)                              :: C , dCdx1 , dCdx2, dCdZ, dCdR

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

            R  = 0.d0; dRdx1 = 0.d0; dRdx2 = 0.d0
            Z  = 0.d0; dZdx1 = 0.d0; dZdx2 = 0.d0
            Vr = 0.d0; Vz    = 0.d0
            C  = 0.d0; dCdx1 = 0.d0; dCdx2 = 0.d0
            do ii = 1, nbf_2d
                R     =  R    + TEMP_TL(ii, getVariableId("R")) *  bfn   (ii,kk)
                dRdx1 = dRdx1 + TEMP_TL(ii, getVariableId("R")) * dbfndx1(ii,kk)
                dRdx2 = dRdx2 + TEMP_TL(ii, getVariableId("R")) * dbfndx2(ii,kk)

                Z     =  Z    + TEMP_TL(ii, getVariableId("Z")) *  bfn   (ii,kk)
                dZdx1 = dZdx1 + TEMP_TL(ii, getVariableId("Z")) * dbfndx1(ii,kk)
                dZdx2 = dZdx2 + TEMP_TL(ii, getVariableId("Z")) * dbfndx2(ii,kk)    
                
                Vr     =  Vr    + TEMP_TL(ii, getVariableId("Vr")) *  bfn   (ii,kk)
                Vz     =  Vz    + TEMP_TL(ii, getVariableId("Vz")) *  bfn   (ii,kk)

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
            

            Ro = 0.d0; Zo = 0.d0 

            do ii = 1, nbf_2d
                Ro = Ro + TLo(NM(ii), getVariableId("R")) * bfn(ii,kk)
                Zo = Zo + TLo(NM(ii), getVariableId("Z")) * bfn(ii,kk)
            end do

            dRdt = (R - Ro)/Dt
            dZdt = (Z - Zo)/Dt

            dCdZ = dCdx1 * dx1dZ + dCdx2 * dx2dZ
            dCdR = dCdx1 * dx1dR + dCdx2 * dx2dR

            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------

            loop_residuals_f:DO IW = 1, NBF_2d
        
                    BIFN =  bfn   (iw,kk)
                    DBIR = dbfndx1(iw,kk) * dx1dR + dbfndx2(iw,kk) * dx2dR
                    DBIZ = dbfndx1(iw,kk) * dx1dZ + dbfndx2(iw,kk) * dx2dZ
        
        
                    TERM_RES     = 0.D0
                    ! TERM_RES(getVariableId("C"))  = PeN * BIFN * ( nR * (Vr-dRdt) + nZ * (Vz-dZdt) ) * KoN * gVar * R/PeN
                    TERM_RES(getVariableId("C"))  = ( BIFN * ( nR * (Vr-dRdt) + nZ * (Vz-dZdt) ) * KoN * gVar - (nR*dCdR + nZ*dCdZ)/PeN ) * R
        
                    !      FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_1d(KK)  * dS
                              
            end do loop_residuals_f
        end do LOOP_GAUSS
      
                ! print*, 'temp_res=', temp_res

        !---------------------------------------------------------------------
        !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
        !---------------------------------------------------------------------
        if ( STORE ) then 
            NM = NM_f(NELEM,1:NBF_2d)
            call MATRIX_STORAGE_RESIDUAL ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
        end if               

    end Subroutine weakHenry

end module Boundary_EquationsDO

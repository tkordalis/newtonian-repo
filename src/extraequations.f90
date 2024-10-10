Module ExtraEquations
    Use ArrayTools,              Only: copyArrayToLocalValues
    Use VariableMapping,         Only: getVariableId
    Use GAUSS_MODULE,            Only: WO_1d, NGAUSS_1d       , &
                                        BFN_E, DFDC_E, DFDE_E  , &
                                        getBasisFunctionsAtFace, &
                                        getNormalVectorAtFace
    Use ELEMENTS_MODULE,         Only: NBF_2d, NEL_2d, NEQ_f, NUNKNOWNS_f
    Use ENUMERATION_MODULE,      Only: NM_MESH
    Use GLOBAL_ARRAYS_MODULE,    Only: TL, TLo, TLb
    Use FLOW_ARRAYS_MODULE,      Only: B_f
    Use MESH_MODULE,             Only: Xm, Ym
    Use PHYSICAL_MODULE,         Only: PeN
    use basis_calculations,      only: BASIS_2d




    Private
    Public :: SurfaceIntegration,   int_Z_dV,   int_n_dot_F, &
                int_n_dot_gradC_z,    int_n_dot_gradC_r,    int_n_dot_UmUmesh_r,    int_n_dot_UmUmesh_z, &
                                                         int_n_dot_UbubblemUmesh_z, int_n_dot_UbubblemUmesh_r

    Contains


    Function int_n_dot_F( NELEM, NED ) Result(TEMP_RES)
        use time_integration, only:dt, time
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  ARGUMENT
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer,      Intent(In)  :: NELEM, NED
        Real(8)                   :: TEMP_RES

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer                    :: KK, ii
        Real(8)                    :: WET
        Real(8)                    :: R, dRdC, dRdE
        Real(8)                    :: Z, dZdC, dZdE
        Real(8)                    :: CJAC, AJAC, dL, nR, nZ, Ro ,Zo ,Vz ,Vr ,C ,dCdZ ,dCdR, dRdt, dZdt

        REAL(8)                    :: TERM_1, term_2_, term_3_
        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc, Zo_loc, Ro_loc, Vz_loc, Vr_loc, C_loc

        Real(8), Parameter                   :: pi = 3.14159265359d0

        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !---------------------------------------------------------------------
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("Z")), nm_mesh(nelem,:), Zo_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("R")), nm_mesh(nelem,:), Ro_loc )

        call copyArrayToLocalValues( TL(:,getVariableId("Vz")), nm_mesh(nelem,:), Vz_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("Vr")), nm_mesh(nelem,:), Vr_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("C")), nm_mesh(nelem,:) , C_loc )

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        TERM_1   = 0.D0
        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
            CALL BASIS_2d&
            ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDZ, DFDR,  NGAUSS_1d )
            
            Ro = 0.d0 ; Zo = 0.d0 ; Vz = 0.d0 ; Vr = 0.d0 ; C = 0.d0 ; dCdZ = 0.d0 ; dCdR = 0.d0
            do ii = 1, nbf_2d
                    Ro    =  Ro    + Ro_loc(ii)  *  bfn   (ii,kk)

                    Zo    =  Zo    + Zo_loc(ii)  *  bfn   (ii,kk)

                    Vz    =  Vz    + Vz_loc(ii) *  bfn   (ii,kk)
                    Vr    =  Vr    + Vr_loc(ii) *  bfn   (ii,kk)

                    C    =  C    + C_loc(ii) *  bfn   (ii,kk)
                    dCdZ = dCdZ  + C_loc(ii) *  DFDZ(ii)
                    dCdR = dCdR  + C_loc(ii) *  DFDR(ii) 
            end do

            dZdt = (Z-Zo) / dt
            dRdt = (R-Ro) / dt

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n

            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,         &
            normalize = .true., forceOnObject = .true.) 
            WET = WO_1d(KK)*dL

            TERM_1  = TERM_1 + ( ( nR*(Vr-dRdt) + nZ*(Vz-dZdt) ) * C - (nR*dCdR + nZ*dCdZ)/PeN ) * (R*WET)
            

        ENDDO LOOP_GAUSS


        TEMP_RES = 2.d0*pi*TERM_1*dt

    End Function int_n_dot_F



    Function int_Z_dV( NELEM, NED ) Result(TEMP_RES)
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  ARGUMENT
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer,      Intent(In)  :: NELEM, NED
        Real(8)                   :: TEMP_RES

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer                    :: KK
        Real(8)                    :: WET
        Real(8)                    :: R, dRdC, dRdE
        Real(8)                    :: Z, dZdC, dZdE
        Real(8)                    :: CJAC, AJAC, dL, nR, nZ

        REAL(8)                    :: TERM_1
        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc

        Real(8), Parameter                   :: pi = 3.14159265359d0

        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !---------------------------------------------------------------------
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        TERM_1   = 0.D0

        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------
        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
            CALL BASIS_2d&
            ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDZ, DFDR, NGAUSS_1d )

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n

            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,         &
            normalize = .true., forceOnObject = .true.) 
            WET = WO_1d(KK)*dL


            TERM_1  = TERM_1 + (Z)*(nr*(R) + (nZ)*(Z) )*R*WET

        ENDDO LOOP_GAUSS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Reference: Karapetsas, Photeinos, Dimakopoulos, Tsamopoulos
        ! Dynamics and motion of a gas bubble in a viscoplastic medium under
        ! acoustic excitation, 2019
        ! int(z)dV = (1/4) int(z n.r)dS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        TEMP_RES = 2.d0*pi*TERM_1 / 4.d0

    End Function int_Z_dV



    Function SurfaceIntegration ( NELEM, NED) Result(TEMP_RES)
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        integer,  intent(In)  :: NELEM, NED
        Real(8)               :: TEMP_RES   
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        Real(8), Parameter         :: PI =  3.1415926535897932384626433D0
        Integer :: KK

        Real(8)  :: WET
        Real(8)  :: R, dRdC, dRdE
        Real(8)  :: Z, dZdC, dZdE
        Real(8)  :: CJAC, AJAC, dL, nR, nZ

        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: bfn, dfdc, dfde
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc


        ! Copy X,Y coordinates to local vectors
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )

        ! Copy Basis Functions at the Corresponding face
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)

        ! Initialize temporary residual
        TEMP_RES = 0.D0

        ! Integrate
        loop_gauss: do KK = 1, NGAUSS_1d

            ! if (kk==1 ) then 
            ! print*, 'SurfaceIntegration'
            ! print*, Z_loc(1) , R_loc(1)
            ! endif
            call BASIS_2d&
            ( KK , Z_loc, R_loc, &
            BFN, DFDC , DFDE , &
            Z  , dZdC , dZdE , &
            R  , dRdC , dRdE , &
            CJAC, AJAC, DFDZ, DFDR, NGAUSS_1d)

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n
            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,        &
            normalize = .true., forceOnObject = .true.)
            WET = WO_1d(KK)*dL

            TEMP_RES = TEMP_RES + (  nr * R  + nz * Z  )*R*WET

        end do loop_gauss

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! Reference: Karapetsas, Photeinos, Dimakopoulos, Tsamopoulos
        ! Dynamics and motion of a gas bubble in a viscoplastic medium under
        ! acoustic excitation, 2019
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        TEMP_RES = 2.d0* PI * TEMP_RES / 3.d0

    End Function SurfaceIntegration

! ********************************************************************
! ********************************************************************

    Function int_n_dot_UbubblemUmesh_z( NELEM, NED ) Result(TEMP_RES)
        use time_integration, only:dt, time
        use PHYSICAL_MODULE, only: mol_bubble, volume_bubble, velocity_bubble
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  ARGUMENT
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer,      Intent(In)  :: NELEM, NED
        Real(8)                   :: TEMP_RES

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer                    :: KK, ii
        Real(8)                    :: WET
        Real(8)                    :: R, dRdC, dRdE
        Real(8)                    :: Z, dZdC, dZdE
        Real(8)                    :: CJAC, AJAC, dL, nR, nZ, Ro ,Zo ,Vz ,Vr ,C ,dCdZ ,dCdR, dRdt, dZdt

        REAL(8)                    :: TERM_1, term_2_, term_3_
        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc, Zo_loc, Ro_loc, Vz_loc, Vr_loc, C_loc

        Real(8), Parameter                   :: pi = 3.14159265359d0

        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !---------------------------------------------------------------------
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("Z")), nm_mesh(nelem,:), Zo_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("R")), nm_mesh(nelem,:), Ro_loc )

        call copyArrayToLocalValues( TL(:,getVariableId("Vz")), nm_mesh(nelem,:), Vz_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("Vr")), nm_mesh(nelem,:), Vr_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("C")), nm_mesh(nelem,:) , C_loc )

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        TERM_1   = 0.D0
        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
            CALL BASIS_2d&
            ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDZ, DFDR, NGAUSS_1d )
            
            Ro = 0.d0 ; Zo = 0.d0 ; Vz = 0.d0 ; Vr = 0.d0 ; C = 0.d0 ; dCdZ = 0.d0 ; dCdR = 0.d0
            do ii = 1, nbf_2d
                    Ro    =  Ro    + Ro_loc(ii)  *  bfn   (ii,kk)

                    Zo    =  Zo    + Zo_loc(ii)  *  bfn   (ii,kk)

                    Vz    =  Vz    + Vz_loc(ii) *  bfn   (ii,kk)
                    Vr    =  Vr    + Vr_loc(ii) *  bfn   (ii,kk)

                    C    =  C    + C_loc(ii) *  bfn   (ii,kk)
                    dCdZ = dCdZ  + C_loc(ii) *  DFDZ(ii)
                    dCdR = dCdR  + C_loc(ii) *  DFDR(ii) 
            end do

            dZdt = (Z-Zo) / dt
            dRdt = (R-Ro) / dt

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n

            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,         &
            normalize = .true., forceOnObject = .true.) 
            WET = WO_1d(KK)*dL

            TERM_1  = TERM_1 + ( ( nZ*(velocity_bubble-dZdt)  ) * mol_bubble/volume_bubble ) * (R*WET)
            

        ENDDO LOOP_GAUSS


        TEMP_RES = 2.d0*pi*TERM_1

    End Function int_n_dot_UbubblemUmesh_z

! ********************************************************************

    Function int_n_dot_UbubblemUmesh_r( NELEM, NED ) Result(TEMP_RES)
        use time_integration, only:dt, time
        use PHYSICAL_MODULE, only: mol_bubble, volume_bubble, velocity_bubble
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  ARGUMENT
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer,      Intent(In)  :: NELEM, NED
        Real(8)                   :: TEMP_RES

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer                    :: KK, ii
        Real(8)                    :: WET
        Real(8)                    :: R, dRdC, dRdE
        Real(8)                    :: Z, dZdC, dZdE
        Real(8)                    :: CJAC, AJAC, dL, nR, nZ, Ro ,Zo ,Vz ,Vr ,C ,dCdZ ,dCdR, dRdt, dZdt

        REAL(8)                    :: TERM_1, term_2_, term_3_
        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc, Zo_loc, Ro_loc, Vz_loc, Vr_loc, C_loc

        Real(8), Parameter                   :: pi = 3.14159265359d0

        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !---------------------------------------------------------------------
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("Z")), nm_mesh(nelem,:), Zo_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("R")), nm_mesh(nelem,:), Ro_loc )

        call copyArrayToLocalValues( TL(:,getVariableId("Vz")), nm_mesh(nelem,:), Vz_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("Vr")), nm_mesh(nelem,:), Vr_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("C")), nm_mesh(nelem,:) , C_loc )

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        TERM_1   = 0.D0
        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
            CALL BASIS_2d&
            ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDZ, DFDR, NGAUSS_1d )
            
            Ro = 0.d0 ; Zo = 0.d0 ; Vz = 0.d0 ; Vr = 0.d0 ; C = 0.d0 ; dCdZ = 0.d0 ; dCdR = 0.d0
            do ii = 1, nbf_2d
                    Ro    =  Ro    + Ro_loc(ii)  *  bfn   (ii,kk)

                    Zo    =  Zo    + Zo_loc(ii)  *  bfn   (ii,kk)

                    Vz    =  Vz    + Vz_loc(ii) *  bfn   (ii,kk)
                    Vr    =  Vr    + Vr_loc(ii) *  bfn   (ii,kk)

                    C    =  C    + C_loc(ii) *  bfn   (ii,kk)
                    dCdZ = dCdZ  + C_loc(ii) *  DFDZ(ii)
                    dCdR = dCdR  + C_loc(ii) *  DFDR(ii) 
            end do

            dZdt = (Z-Zo) / dt
            dRdt = (R-Ro) / dt

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n

            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,         &
            normalize = .true., forceOnObject = .true.) 
            WET = WO_1d(KK)*dL

            TERM_1  = TERM_1 + ( ( nR*( -dRdt) ) * mol_bubble/volume_bubble ) * (R*WET)
            

        ENDDO LOOP_GAUSS


        TEMP_RES = 2.d0*pi*TERM_1

    End Function int_n_dot_UbubblemUmesh_r

! ********************************************************************



    Function int_n_dot_UmUmesh_z( NELEM, NED ) Result(TEMP_RES)
        use time_integration, only:dt, time
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  ARGUMENT
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer,      Intent(In)  :: NELEM, NED
        Real(8)                   :: TEMP_RES

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer                    :: KK, ii
        Real(8)                    :: WET
        Real(8)                    :: R, dRdC, dRdE
        Real(8)                    :: Z, dZdC, dZdE
        Real(8)                    :: CJAC, AJAC, dL, nR, nZ, Ro ,Zo ,Vz ,Vr ,C ,dCdZ ,dCdR, dRdt, dZdt

        REAL(8)                    :: TERM_1, term_2_, term_3_
        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc, Zo_loc, Ro_loc, Vz_loc, Vr_loc, C_loc

        Real(8), Parameter                   :: pi = 3.14159265359d0

        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !---------------------------------------------------------------------
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("Z")), nm_mesh(nelem,:), Zo_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("R")), nm_mesh(nelem,:), Ro_loc )

        call copyArrayToLocalValues( TL(:,getVariableId("Vz")), nm_mesh(nelem,:), Vz_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("Vr")), nm_mesh(nelem,:), Vr_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("C")), nm_mesh(nelem,:) , C_loc )

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        TERM_1   = 0.D0
        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
            CALL BASIS_2d&
            ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDZ, DFDR, NGAUSS_1d )
            
            Ro = 0.d0 ; Zo = 0.d0 ; Vz = 0.d0 ; Vr = 0.d0 ; C = 0.d0 ; dCdZ = 0.d0 ; dCdR = 0.d0
            do ii = 1, nbf_2d
                    Ro    =  Ro    + Ro_loc(ii)  *  bfn   (ii,kk)

                    Zo    =  Zo    + Zo_loc(ii)  *  bfn   (ii,kk)

                    Vz    =  Vz    + Vz_loc(ii) *  bfn   (ii,kk)
                    Vr    =  Vr    + Vr_loc(ii) *  bfn   (ii,kk)

                    C    =  C    + C_loc(ii) *  bfn   (ii,kk)
                    dCdZ = dCdZ  + C_loc(ii) *  DFDZ(ii)
                    dCdR = dCdR  + C_loc(ii) *  DFDR(ii) 
            end do

            dZdt = (Z-Zo) / dt
            dRdt = (R-Ro) / dt

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n

            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,         &
            normalize = .true., forceOnObject = .true.) 
            WET = WO_1d(KK)*dL

            TERM_1  = TERM_1 + ( (  + nZ*(Vz-dZdt) ) * C ) * (R*WET)
            

        ENDDO LOOP_GAUSS 


        TEMP_RES = 2.d0*pi*TERM_1

    End Function int_n_dot_UmUmesh_z

! ********************************************************************

    Function int_n_dot_UmUmesh_r( NELEM, NED ) Result(TEMP_RES)
        use time_integration, only:dt, time
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  ARGUMENT
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer,      Intent(In)  :: NELEM, NED
        Real(8)                   :: TEMP_RES

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer                    :: KK, ii
        Real(8)                    :: WET
        Real(8)                    :: R, dRdC, dRdE
        Real(8)                    :: Z, dZdC, dZdE
        Real(8)                    :: CJAC, AJAC, dL, nR, nZ, Ro ,Zo ,Vz ,Vr ,C ,dCdZ ,dCdR, dRdt, dZdt

        REAL(8)                    :: TERM_1, term_2_, term_3_
        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc, Zo_loc, Ro_loc, Vz_loc, Vr_loc, C_loc

        Real(8), Parameter                   :: pi = 3.14159265359d0

        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !---------------------------------------------------------------------
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("Z")), nm_mesh(nelem,:), Zo_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("R")), nm_mesh(nelem,:), Ro_loc )

        call copyArrayToLocalValues( TL(:,getVariableId("Vz")), nm_mesh(nelem,:), Vz_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("Vr")), nm_mesh(nelem,:), Vr_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("C")), nm_mesh(nelem,:) , C_loc )

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        TERM_1   = 0.D0
        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
            CALL BASIS_2d&
            ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDZ, DFDR, NGAUSS_1d )
            
            Ro = 0.d0 ; Zo = 0.d0 ; Vz = 0.d0 ; Vr = 0.d0 ; C = 0.d0 ; dCdZ = 0.d0 ; dCdR = 0.d0
            do ii = 1, nbf_2d
                    Ro    =  Ro    + Ro_loc(ii)  *  bfn   (ii,kk)

                    Zo    =  Zo    + Zo_loc(ii)  *  bfn   (ii,kk)

                    Vz    =  Vz    + Vz_loc(ii) *  bfn   (ii,kk)
                    Vr    =  Vr    + Vr_loc(ii) *  bfn   (ii,kk)

                    C    =  C    + C_loc(ii) *  bfn   (ii,kk)
                    dCdZ = dCdZ  + C_loc(ii) *  DFDZ(ii)
                    dCdR = dCdR  + C_loc(ii) *  DFDR(ii) 
            end do

            dZdt = (Z-Zo) / dt
            dRdt = (R-Ro) / dt

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n

            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,         &
            normalize = .true., forceOnObject = .true.) 
            WET = WO_1d(KK)*dL

            TERM_1  = TERM_1 + ( ( + nR*(Vr-dRdt) ) * C ) * (R*WET)
            

        ENDDO LOOP_GAUSS


        TEMP_RES = 2.d0*pi*TERM_1

    End Function int_n_dot_UmUmesh_r

! ********************************************************************

    Function int_n_dot_gradC_z( NELEM, NED ) Result(TEMP_RES)
        use time_integration, only:dt, time
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  ARGUMENT
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer,      Intent(In)  :: NELEM, NED
        Real(8)                   :: TEMP_RES

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer                    :: KK, ii
        Real(8)                    :: WET
        Real(8)                    :: R, dRdC, dRdE
        Real(8)                    :: Z, dZdC, dZdE
        Real(8)                    :: CJAC, AJAC, dL, nR, nZ, Ro ,Zo ,Vz ,Vr ,C ,dCdZ ,dCdR, dRdt, dZdt

        REAL(8)                    :: TERM_1, term_2_, term_3_
        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc, Zo_loc, Ro_loc, Vz_loc, Vr_loc, C_loc

        Real(8), Parameter                   :: pi = 3.14159265359d0

        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !---------------------------------------------------------------------
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("Z")), nm_mesh(nelem,:), Zo_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("R")), nm_mesh(nelem,:), Ro_loc )

        call copyArrayToLocalValues( TL(:,getVariableId("Vz")), nm_mesh(nelem,:), Vz_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("Vr")), nm_mesh(nelem,:), Vr_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("C")), nm_mesh(nelem,:) , C_loc )

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        TERM_1   = 0.D0
        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
            CALL BASIS_2d&
            ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDZ, DFDR, NGAUSS_1d )
            
            Ro = 0.d0 ; Zo = 0.d0 ; Vz = 0.d0 ; Vr = 0.d0 ; C = 0.d0 ; dCdZ = 0.d0 ; dCdR = 0.d0
            do ii = 1, nbf_2d
                    Ro    =  Ro    + Ro_loc(ii)  *  bfn   (ii,kk)

                    Zo    =  Zo    + Zo_loc(ii)  *  bfn   (ii,kk)

                    Vz    =  Vz    + Vz_loc(ii) *  bfn   (ii,kk)
                    Vr    =  Vr    + Vr_loc(ii) *  bfn   (ii,kk)

                    C    =  C    + C_loc(ii) *  bfn   (ii,kk)
                    ! dCdZ = dCdZ  + C_loc(ii) *  DFDZ(ii)
                    dCdZ = dCdZ  + 1.d0 *  DFDZ(ii)
                    dCdR = dCdR  + C_loc(ii) *  DFDR(ii) 
            end do

            dZdt = (Z-Zo) / dt
            dRdt = (R-Ro) / dt

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n

            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,         &
            normalize = .true., forceOnObject = .true.) 
            WET = WO_1d(KK)*dL

            ! TERM_1  = TERM_1 + ( (nZ*dCdZ)/PeN ) * (R*WET)
            TERM_1  = TERM_1 + ( (dCdZ) ) !* (R*WET)
            

        ENDDO LOOP_GAUSS


        TEMP_RES = 2.d0*pi*TERM_1

    End Function int_n_dot_gradC_z

! ********************************************************************

    Function int_n_dot_gradC_r( NELEM, NED ) Result(TEMP_RES)
        use time_integration, only:dt, time
        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  ARGUMENT
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer,      Intent(In)  :: NELEM, NED
        Real(8)                   :: TEMP_RES

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        !  LOCAL VARIABLES
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
        Integer                    :: KK, ii
        Real(8)                    :: WET
        Real(8)                    :: R, dRdC, dRdE
        Real(8)                    :: Z, dZdC, dZdE
        Real(8)                    :: CJAC, AJAC, dL, nR, nZ, Ro ,Zo ,Vz ,Vr ,C ,dCdZ ,dCdR, dRdt, dZdt

        REAL(8)                    :: TERM_1, term_2_, term_3_
        Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ

        Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
        Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc, Zo_loc, Ro_loc, Vz_loc, Vr_loc, C_loc

        Real(8), Parameter                   :: pi = 3.14159265359d0

        !---------------------------------------------------------------------
        !  COPY X VECTOR TO LOCAL VECTOR
        !---------------------------------------------------------------------
        call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("Z")), nm_mesh(nelem,:), Zo_loc )
        call copyArrayToLocalValues( TLo(:,getVariableId("R")), nm_mesh(nelem,:), Ro_loc )

        call copyArrayToLocalValues( TL(:,getVariableId("Vz")), nm_mesh(nelem,:), Vz_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("Vr")), nm_mesh(nelem,:), Vr_loc )
        call copyArrayToLocalValues( TL(:,getVariableId("C")), nm_mesh(nelem,:) , C_loc )

        !---------------------------------------------------------------------
        !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
        !---------------------------------------------------------------------
        call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        TERM_1   = 0.D0
        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------

        LOOP_GAUSS: DO KK = 1, NGAUSS_1d
            !---------------------------------------------------------------------
            !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
            !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
            !---------------------------------------------------------------------
            CALL BASIS_2d&
            ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDZ, DFDR,  NGAUSS_1d )
            
            Ro = 0.d0 ; Zo = 0.d0 ; Vz = 0.d0 ; Vr = 0.d0 ; C = 0.d0 ; dCdZ = 0.d0 ; dCdR = 0.d0
            do ii = 1, nbf_2d
                    Ro    =  Ro    + Ro_loc(ii)  *  bfn   (ii,kk)

                    Zo    =  Zo    + Zo_loc(ii)  *  bfn   (ii,kk)

                    Vz    =  Vz    + Vz_loc(ii) *  bfn   (ii,kk)
                    Vr    =  Vr    + Vr_loc(ii) *  bfn   (ii,kk)

                    C    =  C    + C_loc(ii) *  bfn   (ii,kk)
                    dCdZ = dCdZ  + C_loc(ii) *  DFDZ(ii)
                    ! dCdR = dCdR  + C_loc(ii) *  DFDR(ii) 
                    dCdR = dCdR  + 1.d0 *  DFDR(ii) 
            end do

            dZdt = (Z-Zo) / dt
            dRdt = (R-Ro) / dt

            !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n

            call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
            ned, nr, nz, dL,         &
            normalize = .true., forceOnObject = .true.) 
            WET = WO_1d(KK)*dL

            ! TERM_1  = TERM_1 + ( (nR*dCdR )/PeN ) * (R*WET)
            ! TERM_1  = TERM_1 + ( (dCdR )/PeN ) * (R*WET)
            TERM_1  = TERM_1 + ( (dCdR ) )! * (R*WET)
            ! write(*,*) DFDC
            ! pause
            

        ENDDO LOOP_GAUSS


        TEMP_RES = 2.d0*pi*TERM_1

    End Function int_n_dot_gradC_r

End Module ExtraEquations

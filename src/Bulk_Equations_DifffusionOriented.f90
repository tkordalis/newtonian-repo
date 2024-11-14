
Module BulkEquations
    use storage, only: MATRIX_STORAGE_RESIDUAL, MATRIX_STORAGE_JACOBIAN

    interface FEMinterpolation
        Module Procedure FEMinterpolationA1
        Module Procedure FEMinterpolationA2
    end interface FEMinterpolation

    Contains


    Subroutine FEMinterpolationA1(nodes, bfn, var)
        Implicit None 
        Real(8), Dimension(:), Intent(In) :: nodes
        Real(8), Dimension(:), Intent(In) :: bfn
        Real(8), Intent(Out)            :: var
        ! var = 0.d0
        var = dot_product(nodes, bfn)
    End Subroutine FEMinterpolationA1

    Subroutine FEMinterpolationA2(nodes, bfn, var)
        Implicit None 
        Real(8), Dimension(:,:), Intent(In)     :: nodes
        Real(8), Dimension(:),   Intent(In)     :: bfn
        Real(8), Dimension(:),   Intent(Out)  :: var
        integer :: i
            
        ! var = 0.d0
        do i=1,size(var)
            var(i) = dot_product(nodes(i,:), bfn)
        enddo
    End Subroutine FEMinterpolationA2


! *****************************************************************


    Subroutine DOMI_RESIDUAL_flowNgastransport( NELEM, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        Use basis_calculations
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d, NEQ_f, NUNKNOWNS_f, NCD
        Use GAUSS_MODULE,            Only: WO_2d, NGAUSS_2d, BFN_2d
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use MESH_MODULE,             Only: Xm, Ym
        Use TIME_INTEGRATION,        Only: Dt, increment
        use geometry,                only: distance, trace, secondInvariant

        Implicit None
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  ARGUMENTS
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer,                           Intent(In)  :: NELEM
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
        Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
        Logical,                           Intent(In)  :: STORE
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        !  LOCAL VARIABLES  
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
        Integer                          :: KK, II, JJ, IW
        Integer, Dimension(NBF_2d)       :: NM 

        Real(8)  :: BIFN, SBFN
        Real(8)  :: Helem, Uelem, Ha, tlsme, tlsic, tlsmt, hugn
        Real(8)  :: Gravity_Term, jac, Pgp, Cgp, Cogp, dCdt, dCdM
        real(8)  :: continuity_equation, mass_transfer

        Real(8), Dimension(NBF_2d)   :: BFN
        Real(8), Dimension(NCD)      :: Ugp, Uogp, dPgpdX, Xgp, Xogp, X0gp,  dCgpdX
        Real(8), dimension(NCD)      :: Q, S, dXdt, dUdt, dUdM, MomStr, tr_Ptot_d_gradW

        Real(8), Dimension(NCD)      :: gradq, gradk, gradm, momentum_equation, elliptic_grid, elliptic_grid_dummy, mtml

        Real(8), Dimension(NBF_2d)             :: P_n, C_n, Co_n
        Real(8), Dimension(NCD,NBF_2d)         :: X_n, Xo_n, dXdt_n, X0_n, U_n, Uo_n
        Real(8), Dimension(NCD,NBF_2d)         :: dBFNdX_, gradq_, gradk_, gradm_
        Real(8), Dimension(3,3,NBF_2d)         :: gradW_

        Real(8), Dimension(NEQ_f)            :: TERM_RES
        Real(8), Dimension(NBF_2d, NEQ_f)    :: TLo_loc


        Real(8), Dimension(NCD,NCD)          :: gradUgp, dXdX0gp, dX0dXgp
        Real(8), Dimension(3,3)              :: gradU, Gdot, Tnewt, Ptot, I1, gradW, Ptot_d_gradW



        TLo_loc(:,:) = TLo(NM_MESH(NELEM,:),:)


        P_n = TEMP_TL(:,getVariableId("P"))
        C_n = TEMP_TL(:,getVariableId("C"))       ;  Co_n = TLo_loc(:,getVariableId("C"))  

        U_n(1,:) = TEMP_TL(:,getVariableId("Vz")) ;  Uo_n(1,:) = TLo_loc(:,getVariableId("Vz"))
        U_n(2,:) = TEMP_TL(:,getVariableId("Vr")) ;  Uo_n(2,:) = TLo_loc(:,getVariableId("Vr"))

        X_n(1,:) = TEMP_TL(:,getVariableId("Z"))  ;  Xo_n(1,:) = TLo_loc(:,getVariableId("Z")) 
        X_n(2,:) = TEMP_TL(:,getVariableId("R"))  ;  Xo_n(2,:) = TLo_loc(:,getVariableId("R"))

        dXdt_n = (X_n - Xo_n)/dt

        X0_n(1,:) = Xm( NM_MESH(NELEM,:) )
        X0_n(2,:) = Ym( NM_MESH(NELEM,:) )
        
        
        Uelem = 0.d0
        ! Uelem = sum( [ (distance( U_n(:,ii), dXdt_n(:,ii) )/dble(NBF_2d) , ii=1, NBF_2d) ] )
        Uelem = sum( [ ( sqrt( (U_n(1,ii)- dXdt_n(1,ii))**2.d0 + (U_n(2,ii)- dXdt_n(2,ii))**2.d0 ), ii=1,NBF_2d ) ] )


        Gravity_Term =  -1.d0
        Gravity_Term = Gravity_Term*ratio_of_pressures

        ! unity tensor
        I1 = 0.d0 ; I1(1,1) = 1.d0 ; I1(2,2) = 1.d0 ; I1(3,3) = 1.d0 

        Xogp = 0.d0   ;    Uogp = 0.d0   ;  Xgp  = 0.d0   ;  dXdX0gp = 0.d0
        X0gp = 0.d0   ;  dX0dXgp = 0.d0  ;  Pgp  = 0.d0   ;  dPgpdX  = 0.d0
        Ugp  = 0.d0   ;  gradUgp = 0.d0  ;  Ptot = 0.d0   ;  Tnewt = 0.d0
        Cgp  = 0.d0   ;  Cogp  = 0.d0
        !---------------------------------------------------------------------
        !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
        !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
        !---------------------------------------------------------------------
        TEMP_RES = 0.D0
        !---------------------------------------------------------------------
        !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
        !---------------------------------------------------------------------
        LOOP_GAUSS: DO KK = 1, NGAUSS_2d


            BFN = BFN_2d(:,KK)
            
            call basis_spatial_derivs( KK  , X_n, dBFNdX_, Jac )
            

            call basis_interp_scalar ( P_n  , KK, X_n , Pgp , dPgpdX )
            
            call basis_interp_scalar ( C_n  , KK, X_n , Cgp , dCgpdX )
            
            call basis_interp_vector ( U_n  , KK, X_n , Ugp , gradUgp)
            
            call basis_interp_vector ( X_n  , KK, X0_n, Xgp , dXdX0gp)
            
            call basis_interp_vector ( X0_n , KK, X_n , X0gp, dX0dXgp)

            call FEMinterpolation(Xo_n, BFN(:), Xogp )  ;  call FEMinterpolation(Uo_n, BFN(:), Uogp )
            call FEMinterpolation(Co_n, BFN(:), Cogp )

            gradq_ = 0.d0 ; gradk_ = 0.d0 ; gradW_ = 0.d0 ; gradm_ = 0.d0
            do ii=1, NBF_2d
                gradq_(:, ii) = dBFNdX_(: ,ii)
                gradk_(:, ii) = dBFNdX_(: ,ii)
                gradm_(:, ii) = dBFNdX_(: ,ii)
                do jj=1,NCD
                    gradW_( 1:2,jj,ii) = dBFNdX_(: ,ii)
                enddo
                gradW_(3,3, ii) = BFN(ii)/Xgp(2)
            enddo

            Q = sqrt( [ (sum([ (dXdX0gp(ii,jj)**2, jj=1,NCD) ]),ii=1,NCD) ] )
            S(1) = Q(1)/Q(2)  ;  S(2) = Q(2)/Q(1)

            gradU = 0.d0
            gradU(1:2, 1:2) = gradUgp  ; gradU(3,3) = Ugp(2)/Xgp(2)

            Gdot = gradU+transpose(gradU)

            Tnewt = Gdot

            Ptot = - Pgp*I1 + Tnewt

            dXdt = (Xgp - Xogp)/dt  ;  dUdt = (Ugp - Uogp)/dt 
            dCdt = (Cgp - Cogp)/dt
            ! --------------------------------------------------
            !         Material derivative calculation
            ! --------------------------------------------------
            dUdM = dUdt + matmul( (Ugp-dXdt) , (gradU(1:2,1:2)) )

            dUdM = ArN*dUdM

            dCdM = dCdt + dot_product( (Ugp-dXdt) , dCgpdX )

            ! --------------------------------------------------
            ! --------------------------------------------------

            MomStr = dUdM + dPgpdX - [1,0]*Gravity_Term ! + stress terms for evp
            ! --------------------------------------------------
            !       Stabilizing parameters calculation
            ! --------------------------------------------------
            helem   = sum( [ ( sqrt(dBFNdX_(1,ii)**2.d0 + dBFNdX_(2,ii)**2.d0 ), ii=1,NBF_2d ) ] )
            helem   = 1.d0/helem


            hugn    = sum( [ (abs( dot_product((Ugp-dXdt),dBFNdX_(:,ii)) ), ii=1,NBF_2d) ] )
            hugn    = distance( Ugp, dXdt )/(hugn+1.d-8)

            Ha      = sqrt( 1.d0 + secondInvariant(Tnewt) )/sqrt( 1.d0 + secondInvariant(Gdot) )

            tlsme   = 1.d0/sqrt((ReN/Dt)**2 + (ReN*Uelem/Helem)**2 + (Ha/Helem**2)**2)
            ! tlsme   = 1.d0/sqrt((ArN/Dt)**2 + (ArN*Uelem/Helem)**2 + (Ha/Helem**2)**2)

            tlsic   = helem**2/tlsme

            tlsmt   = sqrt((2.d0/dt)**2 + (Uelem)/(hugn+1.d-8)**2)
            tlsmt = 1.d0/tlsmt
            ! --------------------------------------------------
            ! --------------------------------------------------

            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------
            LOOP_RESIDUALS_f:DO IW = 1, NBF_2d

                BIFN = BFN   (IW) ; gradW = gradW_(:,:,iw) ; gradq = gradq_(:,iw) ; gradk = gradk_(:,iw) ; gradm = gradm_(:,iw)

                SBFN = BIFN+tlsmt*dot_product((Ugp-dXdt),gradm)

                Ptot_d_gradW       = matmul(Ptot,gradW)

                tr_Ptot_d_gradW(1) = Ptot_d_gradW(1,1)  ; tr_Ptot_d_gradW(2) = trace(Ptot_d_gradW(2:3,2:3))

                ! Calculation of the bulk equations
                ! =====================================================================================
                momentum_equation   = ( dUdM*BIFN + tr_Ptot_d_gradW - [1,0]*Gravity_Term*BIFN + tlsic*[gradW(1,1),gradW(2,2)+gradW(3,3)]*trace(gradU) ) * Xgp(2)
                ! ------------------------------
                continuity_equation = ( trace(gradU)*BIFN + tlsme*dot_product(gradq,MomStr) ) * Xgp(2)
                ! ------------------------------
                elliptic_grid       = ( eo*S + (1.d0-eo) )*matmul(gradk, dX0dXgp)
                ! ------------------------------
                mass_transfer       = ( PeN*dCdM*SBFN + dot_product(gradm,dCgpdX) ) * Xgp(2) /PeN
                ! ------------------------------

                ! =====================================================================================
                TERM_RES                        = 0.d0
                TERM_RES(getVariableId("Vz" ))  = momentum_equation(1)
                TERM_RES(getVariableId("Vr" ))  = momentum_equation(2)

                TERM_RES(getVariableId("P"  ))  = continuity_equation

                TERM_RES(getVariableId("Z"  ))  = elliptic_grid(1) 
                TERM_RES(getVariableId("R"  ))  = elliptic_grid(2)
                
                TERM_RES(getVariableId("C"  ))  = mass_transfer

                ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
                TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_2d(KK) * Jac

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
    END SUBROUTINE DOMI_RESIDUAL_flowNgastransport


! *****************************************************************

    !------------------------------------------------
    !                  extra jacobian
    !------------------------------------------------
    Subroutine StoreToExtraJacobian(id, nm, temp_jac)
        Use CSR_STORAGE, Only: Ac_f
        Implicit None
        Integer,                 Intent(In) :: id
        Integer, Dimension(:)  , Intent(In) :: nm 
        Real(8), Dimension(:,:), Intent(In) :: temp_jac

        Integer                             :: NBF_2d
        Integer                             :: NEQ_f 
        Integer                             :: iw 
        Integer                             :: jw 
        Integer                             :: ieq 

        NBF_2d = size(temp_jac,1)
        NEQ_f  = size(temp_jac,2)

        do iw = 1, NBF_2d
            do ieq = 1, NEQ_f 
                jw = nm(iw) + ieq - 1
                Ac_f(jw,id) = Ac_f(jw,id) + temp_jac(iw,ieq)
            end do 
        end do 
    End Subroutine StoreToExtraJacobian


End Module BulkEquations



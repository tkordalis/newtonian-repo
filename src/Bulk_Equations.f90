
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



  !---------------------------------------------------------------------
  !                  SUBROUTINE   DOMI_RESIDUAL_fluid
  !---------------------------------------------------------------------
    Subroutine DOMI_RESIDUAL_fluid( NELEM, TEMP_TL, TEMP_RES, STORE )
        Use VariableMapping
        Use basis_calculations
        Use SR_REPRESENTATION
        Use CONTINUATION_MODULE,     Only: INCREMENT
        Use PHYSICAL_MODULE
        Use ELEMENTS_MODULE,         Only: NBF_2d, NEQ_f, NUNKNOWNS_f, NCD
        Use GAUSS_MODULE,            Only: WO_2d, NGAUSS_2d, BFN_2d
        Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
        Use GLOBAL_ARRAYS_MODULE,    Only: TLo
        Use FLOW_ARRAYS_MODULE,      Only: B_f
        Use MESH_MODULE,             Only: Xm, Ym
        Use TIME_INTEGRATION,        Only: Dt
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
        Real(8)  :: Helem, Uelem, Ha, tlsme
        Real(8)  :: Gravity_Term, jac, Pgp
        real(8)  :: continuity_equation

        Real(8), Dimension(NBF_2d)   :: BFN
        Real(8), Dimension(NCD)      :: Ugp, Uogp, dPgpdX, Xgp, Xogp, X0gp
        Real(8), dimension(NCD)      :: Q, S, dXdt, dUdt, dUdM, MomStr, tr_Ptot_d_gradW

        Real(8), Dimension(NCD)      :: gradq, gradk, momentum_equation, elliptic_grid

        Real(8), Dimension(NBF_2d)             :: P_
        Real(8), Dimension(NCD,NBF_2d)         :: X_, Xo_, dXdt_, X0_, U_, Uo_, dPdX_, dBFNdX_, gradq_, gradk_

        Real(8), Dimension(NEQ_f)            :: TERM_RES
        Real(8), Dimension(NBF_2d, NEQ_f)    :: TLo_loc


        Real(8), Dimension(NCD,NCD)          :: gradUgp, dXdX0gp, dX0dXgp
        Real(8), Dimension(3,3)              :: gradU, Gdot, Tnewt, Ptot, I1, gradW, Ptot_d_gradW
        Real(8), Dimension(3,3,NBF_2d)       :: gradW_



        TLo_loc(:,:) = TLo(NM_MESH(NELEM,:),:)


        P_ = TEMP_TL(:,getVariableId("P"))

        U_(1,:) = TEMP_TL(:,getVariableId("Vz")) ;  Uo_(1,:) = TLo_loc(:,getVariableId("Vz"))
        U_(2,:) = TEMP_TL(:,getVariableId("Vr")) ;  Uo_(2,:) = TLo_loc(:,getVariableId("Vr"))

        X_(1,:) = TEMP_TL(:,getVariableId("Z"))  ;  Xo_(1,:) = TLo_loc(:,getVariableId("Z")) 
        X_(2,:) = TEMP_TL(:,getVariableId("R"))  ;  Xo_(2,:) = TLo_loc(:,getVariableId("R"))

        dXdt_ = (X_ - Xo_)/dt

        X0_(1,:) = Xm( NM_MESH(NELEM,:) )
        X0_(2,:) = Ym( NM_MESH(NELEM,:) )

        Uelem = 0.d0
        Uelem = sum( [ (distance( U_(:,ii), dXdt_(:,ii) )/dble(NBF_2d) , ii=1, NBF_2d) ] )

        Gravity_Term = -1.d0
        ! Gravity_Term = Gravity_Term*ratio_of_pressures

        ! unity tensor
        I1 = 0.d0 ; I1(1,1) = 1.d0 ; I1(2,2) = 1.d0 ; I1(3,3) = 1.d0 

        Xogp = 0.d0   ;    Uogp = 0.d0   ;  Xgp  = 0.d0   ;  dXdX0gp = 0.d0
        X0gp = 0.d0   ;  dX0dXgp = 0.d0  ;  Pgp  = 0.d0   ;  dPgpdX  = 0.d0
        Ugp  = 0.d0   ;  gradUgp = 0.d0  ;  Ptot = 0.d0   ;  Tnewt = 0.d0
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

            call basis_spatial_derivs( KK  , X_, dBFNdX_, Jac )

            call basis_interp_scalar ( P_  , KK, X_ , Pgp , dPgpdX )
            call basis_interp_vector ( U_  , KK, X_ , Ugp , gradUgp)
            call basis_interp_vector ( X_  , KK, X0_, Xgp , dXdX0gp)
            call basis_interp_vector ( X0_ , KK, X_ , X0gp, dX0dXgp)

            call FEMinterpolation(Xo_, BFN(:), Xogp )  ;  call FEMinterpolation(Uo_, BFN(:), Uogp )

            gradq_ = 0.d0 ; gradk_ = 0.d0 ; gradW_ = 0.d0
            do ii=1, NBF_2d
                gradq_(:, ii) = dBFNdX_(: ,ii)
                gradk_(:, ii) = dBFNdX_(: ,ii)
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
            ! --------------------------------------------------
            !         Material derivative calculation
            ! --------------------------------------------------
            dUdM = dUdt + matmul( (Ugp-dXdt) , (gradU(1:2,1:2)) )

            dUdM = ArN*dUdM
            ! --------------------------------------------------
            ! --------------------------------------------------

            MomStr = dUdM + dPgpdX - [1,0]*Gravity_Term ! + stress terms for evp
            ! --------------------------------------------------
            !       Stabilizing parameters calculation
            ! --------------------------------------------------
            helem   = sum( [ ( sqrt(dBFNdX_(1,ii)**2.d0 + dBFNdX_(2,ii)**2.d0 ), ii=1,NBF_2d ) ] )
            helem   = 1.d0/helem

            Ha      = sqrt( 1.d0 + secondInvariant(Tnewt) )/sqrt( 1.d0 + secondInvariant(Gdot) )

            tlsme   = 1.d0/sqrt((ReN/Dt)**2 + (ReN*Uelem/Helem)**2 + (Ha/Helem**2)**2)
            ! --------------------------------------------------
            ! --------------------------------------------------

            !---------------------------------------------------------------------
            !    ITERATE OVER WEIGHTING FUNCTIONS
            !---------------------------------------------------------------------
            LOOP_RESIDUALS_f:DO IW = 1, NBF_2d

                BIFN = BFN   (IW) ; gradW = gradW_(:,:,iw) ; gradq = gradq_(:,iw) ; gradk = gradk_(:,iw)

                Ptot_d_gradW       = matmul(Ptot,gradW)

                tr_Ptot_d_gradW(1) = Ptot_d_gradW(1,1)  ; tr_Ptot_d_gradW(2) = trace(Ptot_d_gradW(2:3,2:3))

                ! Calculation of the bulk equations
                ! =====================================================================================
                momentum_equation   = ( dUdM*BIFN + tr_Ptot_d_gradW - [1,0]*Gravity_Term*BIFN ) * Xgp(2)

                continuity_equation = ( trace(gradU)*BIFN + tlsme*dot_product(gradq,MomStr) ) * Xgp(2)

                elliptic_grid       = ( eo*S + (1.d0-eo) )*matmul(gradk, dX0dXgp)

                ! =====================================================================================
                TERM_RES                        = 0.d0
                TERM_RES(getVariableId("Vz" ))  = momentum_equation(1)
                TERM_RES(getVariableId("Vr" ))  = momentum_equation(2)

                TERM_RES(getVariableId("P"  ))  = continuity_equation

                TERM_RES(getVariableId("Z"  ))  = elliptic_grid(1) ; TERM_RES(getVariableId("R"  ))  = elliptic_grid(2)

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
    END SUBROUTINE DOMI_RESIDUAL_fluid



! Subroutine DOMI_RESIDUAL_chemSpecies( NELEM, TEMP_TL, TEMP_RES, STORE )
!    use check_for_floating_point_exceptions
!    Use VariableMapping
!    Use basis_calculations
!    Use CONTINUATION_MODULE,     Only: INCREMENT
!    Use PHYSICAL_MODULE
!    Use ELEMENTS_MODULE,         Only: NBF_2d, NEQ_f, NUNKNOWNS_f
!    Use GAUSS_MODULE,            Only: WO_2d, NGAUSS_2d, BFN_2d, DFDC_2d, DFDE_2d
!    Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
!    Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb
!    Use FLOW_ARRAYS_MODULE,      Only: B_f
!    Use MESH_MODULE,             Only: Xm, Ym
!    Use TIME_INTEGRATION,        Only: Dt, TIME, DTo, dtb
!    Implicit None
!    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
!    !  ARGUMENTS
!    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
!    Integer,                           Intent(In)  :: NELEM
!    Real(8), Dimension(NBF_2d, NEQ_f), Intent(In)  :: TEMP_TL
!    Real(8), Dimension(NBF_2d, NEQ_f), Intent(Out) :: TEMP_RES
!    Logical,                           Intent(In)  :: STORE
!    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
!    !  LOCAL VARIABLES  
!    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
!    Integer  :: KK, II, JJ, IW

!    Real(8)  :: C, Vr ,Vz

!    Real(8)  :: JacT
!    Real(8)  :: dCdR , dCdZ, dCdt
!    Real(8)  :: dCdRo , dCdZo

!    Real(8)  :: Co
!    Real(8)  :: DCDM
   
!    Real(8)  :: BIFN, DBIR, DBIZ
!    Real(8)  :: E_TR, Helem, Uelem, Ha, tlsme, tlsic, tlsce, taudc
   
!    Integer, Dimension(NBF_2d)           :: NM 

   
!    Real(8), Dimension(NEQ_f)            :: TERM_RES
!    Real(8), Dimension(NBF_2d)           :: BFN, DFDX, DFDY, Ksi_loc, Eta_loc
!    Real(8), Dimension(NBF_2d, NEQ_f)    :: TLo_loc, TLb_loc
  
!    Real(8)                              :: X, dXdC, dXdE, Y, dYdC, dYdE, CJAC, AJAC


!     !---------------------------------------------------------------------
!     !  COPY X VECTOR TO LOCAL VECTOR
!     !---------------------------------------------------------------------
!    DO II = 1, NBF_2d
!      JJ = NM_MESH(NELEM,II)
     
!      ! Coordinates of the nodes at the computational domain 
!      Ksi_loc(II) = Xm(JJ) 
!      Eta_loc(II) = Ym(JJ)

!      TLo_loc(II,:) = TLo(JJ,:)
!      TLb_loc(II,:) = TLb(JJ,:)
!    ENDDO
  
   
!     !---------------------------------------------------------------------
!     !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
!     !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
!     !---------------------------------------------------------------------
!    TEMP_RES = 0.D0


!     !---------------------------------------------------------------------
!     !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
!     !---------------------------------------------------------------------
!    LOOP_GAUSS: DO KK = 1, NGAUSS_2d

!      BFN = BFN_2d(:,KK)

!      ! calculate the derivatives of basis functions with respect to physical coordinates at current timestep
!      call basis_spatial_derivatives( KK, Ksi_loc(:) , Eta_loc(:) , dBFNdZ , dBFNdR , JacT  )

    
!      ! give TEMP_TL to the subroutine and take back the variable calculated at the current gauss point
!      ! along with the derivatives of the variable with respect to the physical coordinates i.e. dVrdr
!      ! call basis_interpolation_chainrule( Z_nodes_elem(:), KK, Ksi_loc(:), Eta_loc(:),  Z, dZdKsi, dZdEta )

!      ! call basis_interpolation_chainrule( R_nodes_elem(:), KK, Ksi_loc(:), Eta_loc(:),  R, dRdKsi, dRdEta )

     
!      call basis_interpolation_chainrule( TEMP_TL(:,getVariableId("C")) , KK, Ksi_loc(:), Eta_loc(:),  C , dCdZ , dCdR  )

   

     

!      ! ************************************************************************
!      ! Local variables and its ferivatives calculated at the previous timestep
!      ! ************************************************************************
!      call basis_interpolation_chainrule( TLo_loc(:,getVariableId("C")) , KK, Ksi_loc(:), Eta_loc(:),  Co , dCdZo , dCdRo  )

             
!      dCdt = ( C - Co ) / dt
!     !---------------------------------------------------------------------     
!     !    DEFINE MATERIAL DERIVATIVES FOR MOMENTUM
!     !---------------------------------------------------------------------
!      ! dCdM  = dCdt  + Vr * dCdR   + Vz * dCdZ
!      ! dCdM  = dCdt  + 100.5d0* dCdZ
!      dCdM  = dCdt
     
!      ! dUrdM  = ReN*dUrdM


!     !---------------------------------------------------------------------
!     !    ITERATE OVER WEIGHTING FUNCTIONS
!     !---------------------------------------------------------------------
!     LOOP_RESIDUALS_f:DO IW = 1, NBF_2d
    
!         BIFN = BFN   (IW)  ;  DBIZ = dBFNdZ(IW)  ;  DBIR = dBFNdR(IW)
        
           
!         TERM_RES                        = 0.d0
!         TERM_RES(getVariableId("C" ))  = dCdM*BIFN  +  dCdZ*DBIZ  +  dCdR*DBIR

!         ! TERM_RES(getVariableId("C" ))  = TERM_RES(getVariableId("C" )) - 10.d0*C*BIFN
    
     
    
!         ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
!         TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_2d(KK) * JacT
!     ENDDO LOOP_RESIDUALS_f
         
   
!     ENDDO LOOP_GAUSS
       
    
!     !---------------------------------------------------------------------
!     !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
!     !---------------------------------------------------------------------
!     IF ( STORE ) THEN
       
!         NM = NM_f(NELEM,1:NBF_2d)
         
!         CALL MATRIX_STORAGE_RESIDUAL&
!         ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
         
!     ENDIF


! END SUBROUTINE DOMI_RESIDUAL_chemSpecies

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

  Real(8) Function Perturb (x)
    Implicit None
    !<><><><><><><><><><><><><><><>
    !  ARGUMENTS
    !<><><><><><><><><><><><><><><>
    Real(8), Intent(In)  :: X
    !<><><><><><><><><><><><><><><>
    !  LOCAL VARIABLES
    !<><><><><><><><><><><><><><><>
    Real(8)              :: EP_RES = 1.0D-8
      
    Perturb =  ep_res*dmax1(1.0d0,dabs(x))*sign(1.0d0,x)

  End Function Perturb

 

End Module BulkEquations



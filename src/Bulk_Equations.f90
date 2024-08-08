
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
   use check_for_floating_point_exceptions
   Use VariableMapping
   Use basis_calculations
   Use SR_REPRESENTATION
   Use CONTINUATION_MODULE,     Only: INCREMENT
   Use PHYSICAL_MODULE
   Use ELEMENTS_MODULE,         Only: NBF_2d, NEQ_f, NUNKNOWNS_f, NCD
   Use GAUSS_MODULE,            Only: WO_2d, NGAUSS_2d, BFN_2d, DFDC_2d, DFDE_2d
   Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
   Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb
   Use FLOW_ARRAYS_MODULE,      Only: B_f
   Use MESH_MODULE,             Only: Xm, Ym
   Use TIME_INTEGRATION,        Only: Dt, TIME, DTo, dtb
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
   Integer  :: KK, II, JJ, IW

   Real(8)  :: Ro, Rb, dRdt, QKsi, SKsi, Zo, Zb, dZdt, QEta, SEta
   Real(8)  :: QKsi2, SKsi2, QEta2, SEta2
   Real(8), dimension(NCD)  :: Q, S, dXdt, dUdt, dUdM, MomStr, momentum_equation, tr_Ptot_d_gradW, elliptic_grid

   Real(8)  :: Vr ,Vz , P
   Real(8)  :: Trr, Trz, Tzz, Ttt, Ksi, Eta, Z, R
   
   Real(8)  :: JacT
   Real(8)  :: dVrdR , dVrdZ, dVrdt, dVzdR , dVzdZ, dVzdt
   Real(8)  :: dPdR  , dPdZ

   Real(8)  :: dZdKsi, dZdEta, dRdKsi, dRdEta
   Real(8)  :: dKsidR, dKsidZ, dEtadR, dEtadZ
   Real(8)  :: Vzo, Vzb, dUzdM, dVrdZo, dVrdRo
   Real(8)  :: Vro, Vrb, dUrdM, dVzdZo, dVzdRo

   Real(8)  :: Continuity_Equ, Continuity_Equ2
   Real(8)  :: Drr , Drz , Dzr , Dzz , Dtt
   Real(8)  :: Gdot_rr , Gdot_rz , Gdot_zz , Gdot_tt
   Real(8)  :: Prr, Prz, Pzz, Ptt
   Real(8)  :: BIFN, DBIR, DBIZ, SBFN
   Real(8)  :: Helem, Helem2, Uelem, Ha, Ha2, Rmom, Zmom, tlsme, tlsme2

   Real(8), Dimension(NBF_2d)           :: Ksi_loc, Eta_loc
   Real(8), Dimension(NBF_2d)           :: BFN, dBFNdZ, dBFNdR
   
   Integer, Dimension(NBF_2d)           :: NM 
   Real(8), Dimension(NBF_2d)           :: Zo_nodes_elem, Ro_nodes_elem, Z_nodes_elem, R_nodes_elem,Uz_nodes_elem, Ur_nodes_elem, dZdt_nodes_elem, dRdt_nodes_elem
   ! if the code is 3d this needs to be dim(3,NBF)
   Real(8), Dimension(NBF_2d)             :: P_
   Real(8), Dimension(NCD,NBF_2d)         :: X_, Xo_, dXdt_, X0_, U_, Uo_, dPdX_, dBFNdX_, gradq_, gradk_
   
   Real(8), Dimension(NEQ_f)            :: TERM_RES
   Real(8), Dimension(NBF_2d, NEQ_f)    :: TLo_loc, TLb_loc  
   
   Real(8), Dimension(NCD)              :: Ugp, Uogp, dPgpdX, Xgp, Xogp, X0gp

   Real(8), Dimension(NCD,NCD)          :: gradUgp, dXdX0gp, dX0dXgp
   Real(8), Dimension(3,3)              :: gradU, Gdot, Tnewt, Ptot, I1, gradW, Ptot_d_gradW
   Real(8), Dimension(3,3,NBF_2d)       :: gradW_
   Real(8), Dimension(NCD)              :: gradq, gradk
   real(8) :: continuity_equation
   
   
   Real(8)                              :: Gravity_Term, e_tr, jac, Pgp, Uelem2
    !---------------------------------------------------------------------
    !  COPY X VECTOR TO LOCAL VECTOR
    !---------------------------------------------------------------------
   DO II = 1, NBF_2d
     JJ = NM_MESH(NELEM,II)
     
     ! Coordinates of the nodes at the computational domain 
     Ksi_loc(II) = Xm(JJ) 
     Eta_loc(II) = Ym(JJ)

     TLo_loc(II,:) = TLo(JJ,:)
     TLb_loc(II,:) = TLb(JJ,:)
   ENDDO

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


    !---------------------------------------------------------------------   
    !  DEFINE CHARACTERISTIC VELOCITY
    !---------------------------------------------------------------------   
   
   
   Z_nodes_elem = TEMP_TL(:,getVariableId("Z"))
   R_nodes_elem = TEMP_TL(:,getVariableId("R"))

   Zo_nodes_elem = TLo_loc(:,getVariableId("Z"))
   Ro_nodes_elem = TLo_loc(:,getVariableId("R"))
   Uz_nodes_elem = TEMP_TL(:,getVariableId("Vz")) ; Ur_nodes_elem = TEMP_TL(:,getVariableId("Vr"))

   dZdt_nodes_elem = ( Z_nodes_elem - Zo_nodes_elem ) / DT
   dRdt_nodes_elem = ( R_nodes_elem - Ro_nodes_elem ) / DT


    
    ! do II=1, NBF_2d
    !      Uelem = Uelem + SQRT( ( Uz_nodes_elem(II)-dZdt_nodes_elem(II) )**2 &
    !                              + ( Ur_nodes_elem(II)-dRdt_nodes_elem(II) )**2 ) / DBLE(NBF_2d)
    ! enddo
    !---------------------------------------------------------------------
    !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
    !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
    !---------------------------------------------------------------------
   TEMP_RES = 0.D0

    ! dBFNdX_ = 0.d0  ;  Jac = 0.d0
    Xogp = 0.d0   ;    Uogp = 0.d0
    Xgp  = 0.d0   ;  dXdX0gp = 0.d0
    X0gp = 0.d0   ;  dX0dXgp = 0.d0
    Pgp  = 0.d0   ;  dPgpdX  = 0.d0
    Ugp  = 0.d0   ;  gradUgp = 0.d0
    Ptot = 0.d0   ;  Tnewt = 0.d0

    !---------------------------------------------------------------------
    !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
    !---------------------------------------------------------------------
   LOOP_GAUSS: DO KK = 1, NGAUSS_2d

     BFN = BFN_2d(:,KK)

     ! calculate the derivatives of basis functions with respect to physical coordinates at current timestep
!      call basis_spatial_derivatives( KK,  Z_nodes_elem(:) ,  R_nodes_elem(:) , dBFNdZ , dBFNdR , JacT  )
    
!      ! give TEMP_TL to the subroutine and take back the variable calculated at the current gauss point
!      ! along with the derivatives of the variable with respect to the physical coordinates i.e. dVrdr
!      call basis_interpolation_chainrule( Z_nodes_elem(:), KK, Ksi_loc(:), Eta_loc(:),  Z, dZdKsi, dZdEta )

!      call basis_interpolation_chainrule( R_nodes_elem(:), KK, Ksi_loc(:), Eta_loc(:),  R, dRdKsi, dRdEta )

     


!      call basis_interpolation_chainrule( Ksi_loc(:), KK, Z_nodes_elem(:), R_nodes_elem(:),  Ksi, dKsidZ, dKsidR )

!      call basis_interpolation_chainrule( Eta_loc(:), KK, Z_nodes_elem(:), R_nodes_elem(:),  Eta, dEtadZ, dEtadR )




!      call basis_interpolation_chainrule( TEMP_TL(:,getVariableId("Vr")) , KK, Z_nodes_elem(:), R_nodes_elem(:),  Vr , dVrdZ , dVrdR  )

!      call basis_interpolation_chainrule( TEMP_TL(:,getVariableId("Vz")) , KK, Z_nodes_elem(:), R_nodes_elem(:),  Vz , dVzdZ , dVzdR  )

     


!      call basis_interpolation_chainrule( TEMP_TL(:,getVariableId("P"))  , KK, Z_nodes_elem(:), R_nodes_elem(:),  P  , dPdZ  , dPdR   )
!      ! ************************************************************************
!      ! Local variables and its ferivatives calculated at the previous timestep
!      ! ************************************************************************
!      call basis_interpolation_chainrule( TLo_loc(:,getVariableId("Vr")) , KK, Z_nodes_elem(:), R_nodes_elem(:),  Vro , dVrdZo , dVrdRo  )

!      call basis_interpolation_chainrule( TLo_loc(:,getVariableId("Vz")) , KK, Z_nodes_elem(:), R_nodes_elem(:),  Vzo , dVzdZo , dVzdRo  )

!      call FEMinterpolation(Zo_nodes_elem(:), BFN(:), Zo )
!      call FEMinterpolation(Zo_nodes_elem(:), BFN(:), Zb )

!      call FEMinterpolation(Ro_nodes_elem(:), BFN(:), Ro )
!      call FEMinterpolation(Ro_nodes_elem(:), BFN(:), Rb )
     
    
!      QKsi = DSQRT(dZdKsi**2+dRdKsi**2)  ;  QEta = DSQRT(dZdEta**2+dRdEta**2)

!      SKsi = QKsi/QEta  ;  SEta = QEta/QKsi
     
! ! 888888888888888888888888888888888888888888888888888888888888

!     !---------------------------------------------------------------------     
!     !    DEFINE DEFORMATION RATE TENSOR
!     !---------------------------------------------------------------------

!      Drr = dVrdR
!      Drz = dVzdR
!      Dzr = dVrdZ
!      Dzz = dVzdZ
!      Dtt = Vr/R


!      Continuity_Equ  = dVrdR + Vr/R + dVzdZ

     

!      Gdot_rr = 2.d0 * ( Drr       )
!      Gdot_rz =        ( Drz + Dzr )
!      Gdot_zz = 2.d0 * ( Dzz       )
!      Gdot_tt = 2.d0 * ( Dtt       )

        
!     !---------------------------------------------------------------------     
!     !    DEFINE TIME DERIVATIVES
!     !---------------------------------------------------------------------
!      dZdt   = (Z  - Zo)/DT    
!      dRdt   = (R  - Ro)/DT    

!      dVrdt  = (Vr - Vro)/DT   ;   dVzdt  = (Vz - Vzo)/DT
     

!     !---------------------------------------------------------------------     
!     !    DEFINE MATERIAL DERIVATIVES FOR MOMENTUM
!     !---------------------------------------------------------------------
!      dUzdM  = dVzdt + ( Vz - dZdt ) * dVzdZ  + ( Vr - dRdt ) * dVzdR   
!      dUrdM  = dVrdt + ( Vz - dZdt ) * dVrdZ  + ( Vr - dRdt ) * dVrdR   

          
!      dUrdM  = ArN*dUrdM
!      dUzdM  = ArN*dUzdM

!      Trr = Gdot_rr
!      Trz = Gdot_rz
!      Tzz = Gdot_zz
!      Ttt = Gdot_tt

!      Prr =Trr  
!      Prz =Trz  
!      Pzz =Tzz  
!      Ptt =Ttt  



!      Rmom = dUrdM + dPdR !- (dTrrdr + Trr/R - Ttt/R + dTrzdZ)               
!      Zmom = dUzdM + dPdZ - Gravity_Term!- (dTrzdr + Trz/R         + dTzzdZ) 


     


!      helem = sum([ (sqrt( dBFNdZ(ii)**2 + dBFNdR(ii)**2 ),ii=1,NBF_2d) ])
!      helem = 1.d0/helem


     
!     ! **************************************************************************************************
!     Ha = DSQRT( 1.d0 + 0.5D0*( Trr    **2 + 2.D0*Trz    **2 + Tzz    **2 + Ttt**2 ) ) / &
!         DSQRT(  1.D0    + 0.5D0*( Gdot_rr**2 + 2.D0*Gdot_rz**2 + Gdot_zz**2 + Gdot_tt**2) )
    

!     ! **************************************************************************************************
!     tlsme = 1.D0/DSQRT( (ReN/Dt)**2 + (ReN*Uelem/Helem)**2 + (Ha/Helem**2)**2 )

    
    ! **************************************************************************************************
    call basis_spatial_derivs( KK  , X_, dBFNdX_, Jac )

    call basis_interp_scalar ( P_  , KK, X_ , Pgp , dPgpdX )
    call basis_interp_vector ( U_  , KK, X_ , Ugp , gradUgp)
    call basis_interp_vector ( X_  , KK, X0_, Xgp , dXdX0gp)
    call basis_interp_vector ( X0_ , KK, X_ , X0gp, dX0dXgp)

    call FEMinterpolation(Xo_, BFN(:), Xogp )
    call FEMinterpolation(Uo_, BFN(:), Uogp )

    gradq_ = 0.d0 ; gradk_ = 0.d0 ; gradW_ = 0.d0
    do ii=1, NBF_2d
        gradq_(:, ii) = dBFNdX_(: ,ii)
        gradk_(:, ii) = dBFNdX_(: ,ii)
        do jj=1,NCD
            gradW_( 1:2,jj,ii) = dBFNdX_(: ,ii)
        enddo
        gradW_(3,3, ii) = BFN(ii)/Xgp(2)
    enddo

    Q=sqrt( [ (sum([ (dXdX0gp(ii,jj)**2, jj=1,NCD) ]),ii=1,NCD) ] )
    S(1) = Q(1)/Q(2)  ;  S(2) = Q(2)/Q(1)
    
    gradU = 0.d0
    gradU(1:2, 1:2) = gradUgp  ; gradU(3,3) = Ugp(2)/Xgp(2)

    Gdot = gradU+transpose(gradU)

    Tnewt = Gdot

    Ptot = - Pgp*I1 + Tnewt
    
    Continuity_Equ2 = trace(gradU)
    
    dXdt = ( Xgp - Xogp )/dt  ;  dUdt = ( Ugp - Uogp )/dt 
    
    dUdM = dUdt + matmul( (Ugp-dXdt) , (gradU(1:2,1:2)) )
    
    dUdM = ArN*dUdM
    
    MomStr = dUdM + dPgpdX - [1,0]*Gravity_Term ! + stress terms for evp
    
    helem2 = sum( [ ( sqrt(dBFNdX_(1,ii)**2.d0 + dBFNdX_(2,ii)**2.d0 ), ii=1,NBF_2d ) ] )
    helem2 = 1.d0/helem2
    Ha2 = sqrt( 1.d0 + secondInvariant(Tnewt) )/sqrt( 1.d0 + secondInvariant(Gdot) )
    tlsme2= 1.d0/sqrt((ReN/Dt)**2 + (ReN*Uelem/Helem2)**2 + (Ha2/Helem2**2)**2)

    !---------------------------------------------------------------------
    !    ITERATE OVER WEIGHTING FUNCTIONS
    !---------------------------------------------------------------------
    LOOP_RESIDUALS_f:DO IW = 1, NBF_2d
    
        ! BIFN = BFN   (IW)  ;  DBIZ = dBFNdZ(IW)  ;  DBIR = dBFNdR(IW)
        BIFN = BFN   (IW)  ;  gradW = gradW_(:,:,iw)  ;  gradq = gradq_(:,iw)  ; gradk = gradk_(:,iw)
        
        Ptot_d_gradW       = matmul(Ptot,gradW)
        
        tr_Ptot_d_gradW(1) = Ptot_d_gradW(1,1)  ; tr_Ptot_d_gradW(2) = trace(Ptot_d_gradW(2:3,2:3))

        ! Calculation of the bulk equations
        momentum_equation   = ( dUdM*BIFN + tr_Ptot_d_gradW - [1,0]*Gravity_Term*BIFN ) * Xgp(2)
        

        continuity_equation = ( trace(gradU)*BIFN + tlsme2*dot_product(gradq,MomStr) ) * Xgp(2)
   

        elliptic_grid       = ( eo*S + (1.d0-eo) )*matmul(gradk, dX0dXgp)
        
        
        TERM_RES                        = 0.d0
        TERM_RES(getVariableId("Vz" ))  = momentum_equation(1) ; TERM_RES(getVariableId("Vr" ))  = momentum_equation(2)
        
        TERM_RES(getVariableId("P"  ))  = continuity_equation
        
        TERM_RES(getVariableId("Z"  ))  = elliptic_grid(1) ; TERM_RES(getVariableId("R"  ))  = elliptic_grid(2)
        
        ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
        TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_2d(KK) * Jac
        
        ! TERM_RES(getVariableId("Vz" ))  = (dUzdM*BIFN  +  (Pzz-P)*DBIZ  +  (Prz)*DBIR  - Gravity_Term*BIFN )*R

        ! TERM_RES(getVariableId("Vr" ))  = (dUrdM*BIFN  +  (Prr -P)*DBIR  +  (Prz)*DBIZ +  (Ptt-P)*BIFN/R )*R
    
        ! TERM_RES(getVariableId("P"  ))  = (Continuity_Equ*BIFN + tlsme*(Rmom*DBIR+Zmom*DBIZ) )*R
    

       
        ! TERM_RES(getVariableId("Z"  ))  = ( eo1*SKsi + (1.d0 -eo1))*( dKsidR * DBIR + dKsidZ * DBIZ )
        ! TERM_RES(getVariableId("R"  ))  = ( eo2*SEta + (1.d0 -eo2))*( dEtadR * DBIR + dEtadZ * DBIZ )
        
        ! if ((KK .eq.1 ) .and. (nelem .eq.25)) then
        ! ! print*, dUrdM*BIFN, dUdM(2)*BIFN , dUrdM*BIFN- dUdM(2)*BIFN 
        ! write(404,*) TERM_RES(getVariableId("Z"  )), elliptic_grid(1), TERM_RES(getVariableId("Z"  ))- elliptic_grid(1)
        ! ! pause
        ! endif
    
        ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
        ! TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_2d(KK) * JacT
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



Subroutine DOMI_RESIDUAL_chemSpecies( NELEM, TEMP_TL, TEMP_RES, STORE )
   use check_for_floating_point_exceptions
   Use VariableMapping
   Use basis_calculations
   Use SR_REPRESENTATION
   Use CONTINUATION_MODULE,     Only: INCREMENT
   Use PHYSICAL_MODULE
   Use ELEMENTS_MODULE,         Only: NBF_2d, NEQ_f, NUNKNOWNS_f
   Use GAUSS_MODULE,            Only: WO_2d, NGAUSS_2d, BFN_2d, DFDC_2d, DFDE_2d
   Use ENUMERATION_MODULE,      Only: NM_MESH, NM_f
   Use GLOBAL_ARRAYS_MODULE,    Only: TLo, TLb
   Use FLOW_ARRAYS_MODULE,      Only: B_f
   Use MESH_MODULE,             Only: Xm, Ym
   Use TIME_INTEGRATION,        Only: Dt, TIME, DTo, dtb
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
   Integer  :: KK, II, JJ, IW

   Real(8)  :: C, Vr ,Vz

   Real(8)  :: JacT
   Real(8)  :: dCdR , dCdZ, dCdt
   Real(8)  :: dCdRo , dCdZo

   Real(8)  :: Co
   Real(8)  :: DCDM
   
   Real(8)  :: BIFN, DBIR, DBIZ
   Real(8)  :: E_TR, Helem, Uelem, Ha, Rmom, Zmom, Hugn, tlsme, tlsic, tlsce, taudc
   

   Real(8), Dimension(NBF_2d)           :: Ksi_loc, Eta_loc, dBFNdZ, dBFNdR
   
   Integer, Dimension(NBF_2d)           :: NM 

   
   Real(8), Dimension(NEQ_f)            :: TERM_RES
   Real(8), Dimension(NBF_2d)           :: BFN, DFDX, DFDY
   Real(8), Dimension(NBF_2d, NEQ_f)    :: TLo_loc, TLb_loc
  
   Real(8)                              :: X, dXdC, dXdE, Y, dYdC, dYdE, CJAC, AJAC


    !---------------------------------------------------------------------
    !  COPY X VECTOR TO LOCAL VECTOR
    !---------------------------------------------------------------------
   DO II = 1, NBF_2d
     JJ = NM_MESH(NELEM,II)
     
     ! Coordinates of the nodes at the computational domain 
     Ksi_loc(II) = Xm(JJ) 
     Eta_loc(II) = Ym(JJ)

     TLo_loc(II,:) = TLo(JJ,:)
     TLb_loc(II,:) = TLb(JJ,:)
   ENDDO
  
   
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

     ! calculate the derivatives of basis functions with respect to physical coordinates at current timestep
     call basis_spatial_derivatives( KK, Ksi_loc(:) , Eta_loc(:) , dBFNdZ , dBFNdR , JacT  )

    
     ! give TEMP_TL to the subroutine and take back the variable calculated at the current gauss point
     ! along with the derivatives of the variable with respect to the physical coordinates i.e. dVrdr
     ! call basis_interpolation_chainrule( Z_nodes_elem(:), KK, Ksi_loc(:), Eta_loc(:),  Z, dZdKsi, dZdEta )

     ! call basis_interpolation_chainrule( R_nodes_elem(:), KK, Ksi_loc(:), Eta_loc(:),  R, dRdKsi, dRdEta )

     
     call basis_interpolation_chainrule( TEMP_TL(:,getVariableId("C")) , KK, Ksi_loc(:), Eta_loc(:),  C , dCdZ , dCdR  )

   

     

     ! ************************************************************************
     ! Local variables and its ferivatives calculated at the previous timestep
     ! ************************************************************************
     call basis_interpolation_chainrule( TLo_loc(:,getVariableId("C")) , KK, Ksi_loc(:), Eta_loc(:),  Co , dCdZo , dCdRo  )

             
     dCdt = ( C - Co ) / dt
    !---------------------------------------------------------------------     
    !    DEFINE MATERIAL DERIVATIVES FOR MOMENTUM
    !---------------------------------------------------------------------
     ! dCdM  = dCdt  + Vr * dCdR   + Vz * dCdZ
     ! dCdM  = dCdt  + 100.5d0* dCdZ
     dCdM  = dCdt
     
     ! dUrdM  = ReN*dUrdM


    !---------------------------------------------------------------------
    !    ITERATE OVER WEIGHTING FUNCTIONS
    !---------------------------------------------------------------------
    LOOP_RESIDUALS_f:DO IW = 1, NBF_2d
    
        BIFN = BFN   (IW)  ;  DBIZ = dBFNdZ(IW)  ;  DBIR = dBFNdR(IW)
        
           
        TERM_RES                        = 0.d0
        TERM_RES(getVariableId("C" ))  = dCdM*BIFN  +  dCdZ*DBIZ  +  dCdR*DBIR

        ! TERM_RES(getVariableId("C" ))  = TERM_RES(getVariableId("C" )) - 10.d0*C*BIFN
    
     
    
        ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
        TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_2d(KK) * JacT
    ENDDO LOOP_RESIDUALS_f
         
      


    call check_fp_exceptions(dBFNdZ, "dBFNdZ")
    call check_fp_exceptions(dBFNdR, "dBFNdR")
    call check_fp_exceptions(JacT, "JacT")
    ! call check_fp_exceptions(Z, "Z")
    ! call check_fp_exceptions(dZdKsi, "dZdKsi")
    ! call check_fp_exceptions(dZdEta, "dZdEta")
    ! call check_fp_exceptions(R, "R")
    ! call check_fp_exceptions(dRdKsi, "dRdKsi")
    ! call check_fp_exceptions(dRdEta, "dRdEta")
    call check_fp_exceptions(C, "C")
    call check_fp_exceptions(dCdZ, "dCdZ")
    call check_fp_exceptions(dCdR, "dCdR")
    call check_fp_exceptions(Co, "Co")
    call check_fp_exceptions(dCdZo, "dCdZo")
    call check_fp_exceptions(dCdRo, "dCdRo")
    call check_fp_exceptions(dCdt, "dCdt")
    call check_fp_exceptions(dCdM, "dCdM")
    call check_fp_exceptions(TEMP_RES, "TEMP_RES")
    



    ENDDO LOOP_GAUSS
       
    
    !---------------------------------------------------------------------
    !  STORE THE ELEMENT RESIDUAL VECTOR IN THE GLOBAL VECTOR B
    !---------------------------------------------------------------------
    IF ( STORE ) THEN
       
        NM = NM_f(NELEM,1:NBF_2d)
         
        CALL MATRIX_STORAGE_RESIDUAL&
        ( TEMP_RES, NM, NBF_2d, NEQ_f, B_f, NUNKNOWNS_f )
         
    ENDIF


END SUBROUTINE DOMI_RESIDUAL_chemSpecies

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



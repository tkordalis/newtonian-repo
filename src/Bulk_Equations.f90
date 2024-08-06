
Module BulkEquations
    use storage, only: MATRIX_STORAGE_RESIDUAL, MATRIX_STORAGE_JACOBIAN

  Contains

  
  Subroutine FEMinterpolation(var, nodes, bfn)
        Implicit None 
        Real(8), Intent(InOut)            :: var
        Real(8), Dimension(:), Intent(In) :: nodes
        Real(8), Dimension(:), Intent(In) :: bfn
    
        var = dot_product(nodes, bfn)
  End Subroutine FEMinterpolation



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
   use geometry, only: distance
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

   Real(8)  :: Vr ,Vz , P
   Real(8)  :: Trr, Trz, Tzz, Ttt, Ksi, Eta, Z, R
   
   Real(8)  :: JacT
   Real(8)  :: dVrdR , dVrdZ, dVrdt, dVzdR , dVzdZ, dVzdt
   Real(8)  :: dPdR  , dPdZ

   Real(8)  :: dZdKsi, dZdEta, dRdKsi, dRdEta
   Real(8)  :: dKsidR, dKsidZ, dEtadR, dEtadZ
   Real(8)  :: Vzo, Vzb, dUzdM, dVrdZo, dVrdRo
   Real(8)  :: Vro, Vrb, dUrdM, dVzdZo, dVzdRo

   Real(8)  :: Continuity_Equ
   Real(8)  :: Drr , Drz , Dzr , Dzz , Dtt
   Real(8)  :: Gdot_rr , Gdot_rz , Gdot_zz , Gdot_tt
   Real(8)  :: Prr, Prz, Pzz, Ptt
   Real(8)  :: BIFN, DBIR, DBIZ, SBFN
   Real(8)  :: Helem, Uelem, Ha, Rmom, Zmom, tlsme

   Real(8), Dimension(NBF_2d)           :: Ksi_loc, Eta_loc
   Real(8), Dimension(NBF_2d)           :: BFN, dBFNdZ, dBFNdR
   
   Integer, Dimension(NBF_2d)           :: NM 
   Real(8), Dimension(NBF_2d)           :: Zo_nodes_elem, Ro_nodes_elem, Z_nodes_elem, R_nodes_elem
   Real(8), Dimension(NBF_2d)           :: Uz_nodes_elem, Ur_nodes_elem, dZdt_nodes_elem, dRdt_nodes_elem

   ! if the code is 3d this needs to be dim(3,NBF)
   Real(8), Dimension(NBF_2d)             :: P_
   Real(8), Dimension(NCD,NBF_2d)         :: X_, Xo_, dXdt_, X0_, U_, dPdX_, dbfndx_
   
   Real(8), Dimension(NEQ_f)            :: TERM_RES
   Real(8), Dimension(NBF_2d, NEQ_f)    :: TLo_loc, TLb_loc  
   
   Real(8), Dimension(3,3)              :: S_tensor_o, S_tensor, Stress_tensor_o, Stress_tensor
   Real(8), Dimension(NCD)              :: Ugp, dPgpdX
   Real(8), Dimension(NCD,NCD)          :: gradUgp
   
   
   Real(8)                              :: Gravity_Term, e_tr, jac2, Pgp, Uelem2
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

   X_(1,:) = TEMP_TL(:,getVariableId("Z"))  ;  Xo_(1,:) = TLo_loc(:,getVariableId("Z")) 
   X_(2,:) = TEMP_TL(:,getVariableId("R"))  ;  Xo_(2,:) = TLo_loc(:,getVariableId("R")) 

   dXdt_ = (X_ - Xo_)/dt

   U_(1,:) = TEMP_TL(:,getVariableId("Vz"))
   U_(2,:) = TEMP_TL(:,getVariableId("Vr"))

    Gravity_Term = -1.d0
   ! Gravity_Term = Gravity_Term*ratio_of_pressures

    !---------------------------------------------------------------------   
    !  DEFINE CHARACTERISTIC VELOCITY
    !---------------------------------------------------------------------   
   
   Uz_nodes_elem = TEMP_TL(:,getVariableId("Vz"))
   Ur_nodes_elem = TEMP_TL(:,getVariableId("Vr"))

   Z_nodes_elem = TEMP_TL(:,getVariableId("Z"))
   R_nodes_elem = TEMP_TL(:,getVariableId("R"))

   Zo_nodes_elem = TLo_loc(:,getVariableId("Z"))
   Ro_nodes_elem = TLo_loc(:,getVariableId("R"))

   
   dZdt_nodes_elem = ( Z_nodes_elem - Zo_nodes_elem ) / DT

   dRdt_nodes_elem = ( R_nodes_elem - Ro_nodes_elem ) / DT 


   Uelem = 0.d0
   do II=1, NBF_2d
        Uelem = Uelem + SQRT( ( Uz_nodes_elem(II)-dZdt_nodes_elem(II) )**2 &
                            + ( Ur_nodes_elem(II)-dRdt_nodes_elem(II) )**2   ) / DBLE(NBF_2d)
        ! Uelem = Uelem + SQRT( ( Uz_nodes_elem(II) )**2 &
        !                         + ( Ur_nodes_elem(II) )**2 ) / DBLE(NBF_2d)
   enddo

   Uelem2 = sum( [ (distance( U_(:,ii), dXdt_(:,ii) )/dble(NBF_2d) , ii=1, NBF_2d) ] )
  print*, Uelem, Uelem2, Uelem- Uelem2
  pause
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
     call basis_spatial_derivatives( KK,  Z_nodes_elem(:) ,  R_nodes_elem(:) , dBFNdZ , dBFNdR , JacT  )

    
     ! give TEMP_TL to the subroutine and take back the variable calculated at the current gauss point
     ! along with the derivatives of the variable with respect to the physical coordinates i.e. dVrdr
     call basis_interpolation_chainrule( Z_nodes_elem(:), KK, Ksi_loc(:), Eta_loc(:),  Z, dZdKsi, dZdEta )

     call basis_interpolation_chainrule( R_nodes_elem(:), KK, Ksi_loc(:), Eta_loc(:),  R, dRdKsi, dRdEta )

     call basis_interpolation_chainrule( Ksi_loc(:), KK, Z_nodes_elem(:), R_nodes_elem(:),  Ksi, dKsidZ, dKsidR )

     call basis_interpolation_chainrule( Eta_loc(:), KK, Z_nodes_elem(:), R_nodes_elem(:),  Eta, dEtadZ, dEtadR )



     call basis_interpolation_chainrule( TEMP_TL(:,getVariableId("Vr")) , KK, Z_nodes_elem(:), R_nodes_elem(:),  Vr , dVrdZ , dVrdR  )

     call basis_interpolation_chainrule( TEMP_TL(:,getVariableId("Vz")) , KK, Z_nodes_elem(:), R_nodes_elem(:),  Vz , dVzdZ , dVzdR  )

     ! call basis_interp_vector( U_ , kk , X_ , Ugp , gradUgp )
     call basis_interp_scalar( P_ , kk , X_ , Pgp , dPgpdX )
     


     call basis_interpolation_chainrule( TEMP_TL(:,getVariableId("P"))  , KK, Z_nodes_elem(:), R_nodes_elem(:),  P  , dPdZ  , dPdR   )

       ! print*, P, Pgp, P- Pgp
       ! print*," " 
       ! print*, dPdZ, dPgpdX(1), dPdZ - dPgpdX(1)
       ! print*," " 
       ! print*, dPdr, dPgpdX(2), dPdr - dPgpdX(2)
       ! print*," " 

       ! pause
     

     ! ************************************************************************
     ! Local variables and its ferivatives calculated at the previous timestep
     ! ************************************************************************
     call basis_interpolation_chainrule( TLo_loc(:,getVariableId("Vr")) , KK, Z_nodes_elem(:), R_nodes_elem(:),  Vro , dVrdZo , dVrdRo  )

     call basis_interpolation_chainrule( TLo_loc(:,getVariableId("Vz")) , KK, Z_nodes_elem(:), R_nodes_elem(:),  Vzo , dVzdZo , dVzdRo  )

     
     call FEMinterpolation(Zo  , Zo_nodes_elem(:), BFN(:))
     call FEMinterpolation(Zb  , Zo_nodes_elem(:), BFN(:))

     call FEMinterpolation(Ro  , Ro_nodes_elem(:), BFN(:))
     call FEMinterpolation(Rb  , Ro_nodes_elem(:), BFN(:))

    ! !-----------------------------------------------------------------------
    ! !    DEFINE SCALE FACTOR & RATIOS
    ! !-----------------------------------------------------------------------
     QKsi = DSQRT(dZdKsi**2+dRdKsi**2)  ;  QEta = DSQRT(dZdEta**2+dRdEta**2)

     SKsi = QKsi/QEta  ;  SEta = QEta/QKsi
    !---------------------------------------------------------------------     
    !    DEFINE DEFORMATION RATE TENSOR
    !---------------------------------------------------------------------
     Drr = dVrdR
     Drz = dVzdR
     Dzr = dVrdZ
     Dzz = dVzdZ
     Dtt = Vr/R

     Continuity_Equ = dVrdR + Vr/R + dVzdZ

     Gdot_rr = 2.d0 * ( Drr       )
     Gdot_rz =        ( Drz + Dzr )
     Gdot_zz = 2.d0 * ( Dzz       )
     Gdot_tt = 2.d0 * ( Dtt       )

        
    !---------------------------------------------------------------------     
    !    DEFINE TIME DERIVATIVES
    !---------------------------------------------------------------------
     dZdt   = (Z  - Zo)/DT    
     dRdt   = (R  - Ro)/DT    

     dVrdt  = (Vr - Vro)/DT   ;   dVzdt  = (Vz - Vzo)/DT
     
         

    !---------------------------------------------------------------------     
    !    DEFINE MATERIAL DERIVATIVES FOR MOMENTUM
    !---------------------------------------------------------------------
     dUrdM  = dVrdt  + ( Vr - dRdt ) * dVrdR   + ( Vz - dZdt ) * dVrdZ
     dUzdM  = dVzdt  + ( Vr - dRdt ) * dVzdR   + ( Vz - dZdt ) * dVzdZ

     ! dUrdM  = dVrdt  + Vr * dVrdR   + Vz * dVrdZ
     ! dUzdM  = dVzdt  + Vr * dVzdR   + Vz * dVzdZ 
     
     dUrdM  = ArN*dUrdM
     dUzdM  = ArN*dUzdM


    !-----------------------------------------------------------------------
    !     DEFINE DEVIATORIC TENSOR
    !-----------------------------------------------------------------------



    !---------------------------------------------------------------------     
    !    DEFINE EXTRA STRESS TENSOR (DEVSS-G FORMULATION)
    !---------------------------------------------------------------------
     Trr = Gdot_rr
     Trz = Gdot_rz
     Tzz = Gdot_zz
     Ttt = Gdot_tt

     Prr =Trr  
     Prz =Trz  
     Pzz =Tzz  
     Ptt =Ttt  



     Rmom = dUrdM + dPdR !- (dTrrdr + Trr/R - Ttt/R + dTrzdZ)               
     Zmom = dUzdM + dPdZ - Gravity_Term!- (dTrzdr + Trz/R         + dTzzdZ) 


     helem = sum([ (sqrt( dBFNdZ(ii)**2 + dBFNdR(ii)**2 ),ii=1,NBF_2d) ])
     helem = 1.d0/helem
     
    ! **************************************************************************************************
    Ha = DSQRT( 1.d0 + 0.5D0*( Trr    **2 + 2.D0*Trz    **2 + Tzz    **2 + Ttt**2 ) ) / &
        DSQRT(  1.D0    + 0.5D0*( Gdot_rr**2 + 2.D0*Gdot_rz**2 + Gdot_zz**2 + Gdot_tt**2) )
    
    ! **************************************************************************************************
    tlsme = 1.D0/DSQRT( (ReN/Dt)**2 + (ReN*Uelem/Helem)**2 + (Ha/Helem**2)**2 )
    ! if ( store) write(404,'(I5,2X,F12.8, 2X,F12.8, 2X,F12.8, 2X,F12.8, 2X,F12.8, 2X)') nelem, Uelem, Helem, Ha, DSQRT( (ReN/Dt)**2 + (ReN*Uelem/Helem)**2 + (Ha/Helem**2)**2 ), tlsme
    
    ! **************************************************************************************************
    
    ! if (kk==1 .and. store) then
    !   write(404,'(7(f16.8))'),  dUzdM, Pzz,P,Prz
    ! endif

    !---------------------------------------------------------------------
    !    ITERATE OVER WEIGHTING FUNCTIONS
    !---------------------------------------------------------------------
    LOOP_RESIDUALS_f:DO IW = 1, NBF_2d
    
        BIFN = BFN   (IW)  ;  DBIZ = dBFNdZ(IW)  ;  DBIR = dBFNdR(IW)
           
   
        TERM_RES                        = 0.d0

        TERM_RES(getVariableId("Vr" ))  = (dUrdM*BIFN  +  (Prr -P)*DBIR  +  (Prz)*DBIZ +  (Ptt-P)*BIFN/R )*R
    
        TERM_RES(getVariableId("Vz" ))  = (dUzdM*BIFN  +  (Pzz-P)*DBIZ  +  (Prz)*DBIR  - Gravity_Term*BIFN )*R
    
        TERM_RES(getVariableId("P"  ))  = (Continuity_Equ*BIFN + tlsme*(Rmom*DBIR+Zmom*DBIZ) )*R
       
        TERM_RES(getVariableId("Z"  ))  = ( eo1*SKsi + (1.d0 -eo1))*( dKsidR * DBIR + dKsidZ * DBIZ )
        TERM_RES(getVariableId("R"  ))  = ( eo2*SEta + (1.d0 -eo2))*( dEtadR * DBIR + dEtadZ * DBIZ )
    
        ! FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
        TEMP_RES(IW,1:NEQ_f) = TEMP_RES(IW,1:NEQ_f) + TERM_RES(1:NEQ_f)* WO_2d(KK) * JacT
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



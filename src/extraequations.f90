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
 Use PHYSICAL_MODULE,         Only: initial_position, position_o, &
                                    position, Pambient_o_Pchar, ratio_of_pressures



  Private
  Public :: CentroidOfTheBubble, DragForceCalculation, SurfaceIntegration, &
            quantitiesOfTheBubble,                     int_Z_dV, inflowVelocity

  
  Real(8) :: Bubble_volume

Contains


function inflowVelocity( NELEM, NED ) Result(TEMP_RES)
    Implicit None
    Integer,      Intent(In)  :: NELEM, NED
    Integer                   :: KK, ii
    Real(8)                   :: TEMP_RES

    Real(8)                    :: R, dRdC, dRdE
    Real(8)                    :: Z, dZdC, dZdE
    Real(8)                    :: CJAC, AJAC, JacT, dL, nR, nZ, nl

    REAL(8)                    :: TERM_1
    Real(8)                    :: WET

    Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ
    Real(8), Dimension(NEQ_f)  :: TERM_RES
       
    Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
    Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc
    Real(8), Dimension(:)  , Allocatable :: Vr_loc, Vz_loc
    Real(8)                              :: Vr, Vz
    
    Real(8), Parameter                   :: pi = 3.14159265359d0


    call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
    call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
    call copyArrayToLocalValues( TL(:,getVariableId("Vz")), nm_mesh(nelem,:), Vz_loc )
    call copyArrayToLocalValues( TL(:,getVariableId("Vr")), nm_mesh(nelem,:), Vr_loc )


    call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)


    TEMP_RES = 0.D0
    TERM_1   = 0.D0
    Vr = 0.d0
    Vz = 0.d0

  
    LOOP_GAUSS: DO KK = 1, NGAUSS_1d


        CALL BASIS_2d&
           ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
            CJAC, AJAC, DFDR, DFDZ, NGAUSS_1d )

        do ii = 1, nbf_2d
                
                Vr    =  Vr    + Vr_loc(ii) *  bfn   (ii,kk)
                
                Vz    =  Vz    + Vz_loc(ii) *  bfn   (ii,kk)
                
        end do

         
        call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
                                        ned, nr, nz, dL, normalize = .true., forceOnObject = .true.) 

        WET = WO_1d(KK)*dL
   
        TERM_1  = TERM_1 + (nr*Vr + nZ*Vz )*R*WET
       
    ENDDO LOOP_GAUSS
  

    TEMP_RES = 2.d0*TERM_1 


end function inflowVelocity




 Subroutine quantitiesOfTheBubble( Volume )
   Implicit None 
   Real(8), Intent(In) :: Volume

   Bubble_volume     = Volume

 End Subroutine quantitiesOfTheBubble


Function CentroidOfTheBubble( NELEM, NED ) Result(TEMP_RES)
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
     Real(8)                    :: CJAC, AJAC, JacT, dL, nR, nZ
    
     REAL(8)                    :: TERM_1
     Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ
     Real(8), Dimension(NEQ_f)  :: TERM_RES
     
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
         CJAC, AJAC, DFDR, DFDZ, NGAUSS_1d )

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
     !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

     TEMP_RES = 2.d0*pi*TERM_1 / (4.d0*Bubble_volume)

End Function CentroidOfTheBubble


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
   Real(8)                    :: CJAC, AJAC, JacT, dL, nR, nZ, nl
  
   Real(8)                    :: BIFN, DBIX, DBIY, DBIZ
   Real(8)                    :: BJFN, DBJX, DBJY, DBJZ
   
   Real(8)                    :: Ux, dUxdX, dUxdY
   Real(8)                    :: Uy, dUydX, dUydY

   REAL(8)                    :: TERM_1
   Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ
   Real(8), Dimension(NEQ_f)  :: TERM_RES
   
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
         CJAC, AJAC, DFDR, DFDZ, NGAUSS_1d )

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




Function DragForceCalculation( NELEM, NED ) Result(DragForce)
  use sr_representation, only: set_S_get_Stress
  implicit none
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
    !  ARGUMENT
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
    Integer,      Intent(In)  :: NELEM, NED
    Real(8)                   :: DragForce
        
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
    !  LOCAL VARIABLES
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
    Integer                    :: KK, II
    Integer                    :: IROW, JCOL
  
    Real(8)                    :: WET
    Real(8)                    :: R, dRdC, dRdE
    Real(8)                    :: Z, dZdC, dZdE
    Real(8)                    :: CJAC, AJAC, JacT, dL, nR, nZ, tr, tz
    Real(8)                    :: Prr, Prz, Pzz

    REAL(8)                    :: TERM_1
    Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ
    Real(8), Dimension(NEQ_f)  :: TERM_RES
     
    Real(8), Dimension(:,:), Allocatable :: BFN, DFDC, DFDE
    Real(8), Dimension(:)  , Allocatable :: Z_loc, R_loc

    Real(8), Dimension(:)  , Allocatable :: Szz_loc, Srz_loc, Srr_loc, Stt_loc
    Real(8), Dimension(:)  , Allocatable :: Prr_loc, Prz_loc, Pzz_loc

    Real(8), Dimension(:)  , Allocatable :: P_loc,  Pgeneralised_loc

    Real(8), Parameter                   :: pi = 3.14159265359d0

    Real(8), Dimension(3,3) :: S_tensor, T_tensor    
    REAL(8)                 :: Tzz, Trz, Trr

    
    !---------------------------------------------------------------------
  !  COPY X VECTOR TO LOCAL VECTOR
  !---------------------------------------------------------------------
     
    call copyArrayToLocalValues( TL(:,getVariableId("Z")), nm_mesh(nelem,:), Z_loc )
    call copyArrayToLocalValues( TL(:,getVariableId("R")), nm_mesh(nelem,:), R_loc )
    call copyArrayToLocalValues( TL(:,getVariableId("P")), nm_mesh(nelem,:), P_loc )
    
    call copyArrayToLocalValues( TL(:,getVariableId("Srr")),  nm_mesh(nelem,:), Srr_loc )
    call copyArrayToLocalValues( TL(:,getVariableId("Srz")),  nm_mesh(nelem,:), Srz_loc )
    call copyArrayToLocalValues( TL(:,getVariableId("Szz")),  nm_mesh(nelem,:), Szz_loc )
    call copyArrayToLocalValues( TL(:,getVariableId("Stt")),  nm_mesh(nelem,:), Stt_loc )

    allocate(Pgeneralised_loc(size(P_loc)))

    allocate( Prr_loc(size(Srr_loc)) )
    allocate( Prz_loc(size(Srz_loc)) )
    allocate( Pzz_loc(size(Szz_loc)) )


    do kk = 1, size(Prr_loc)
        ! At time t=0  the hydrostatic pressure is Po = Pambient + rho*g*ho - rho*g*z(node)
        ! At time t=t1 the coordinate system is at h=h1 because because it has risen by Dx
        ! h1 = ho - Dx, so the hydrostatic pressure is P1 = Pambient + rho*g*(ho-Dx) - rho*g*z(node)
        ! P1 = Pambient + rho*g*h1 - rho*g*z(node)

        S_tensor = 0.d0   ;   T_tensor = 0.d0

        S_Tensor(1,1) =  Srr_loc(kk)
        S_Tensor(1,2) =  Srz_loc(kk)
        S_Tensor(2,1) =  Srz_loc(kk)
        S_Tensor(2,2) =  Szz_loc(kk)
        S_Tensor(3,3) =  Stt_loc(kk)

        call set_S_get_Stress(S_tensor, T_tensor)

        Trr = T_tensor(1,1)
        Trz = T_tensor(1,2)
        Tzz = T_tensor(2,2)


        Pgeneralised_loc(kk)  =   P_loc(kk) - ( Pambient_o_Pchar + ratio_of_pressures*(position - Z_loc(kk)) )
        Prr_loc(kk)           = - Pgeneralised_loc(kk) + Trr
        Prz_loc(kk)           =   Trz
        Pzz_loc(kk)           = - Pgeneralised_loc(kk) + Tzz
    enddo
    !---------------------------------------------------------------------
  !  COPY BASIS FUNCTIONS & THEIR DERIVATIVES TO LOCAL VECTORS
  !---------------------------------------------------------------------
    call getBasisFunctionsAtFace(ned, bfn, dfdc, dfde)

    !---------------------------------------------------------------------
  !  INITIALIZE WORKING (TEMPORARY) AREAS FOR ELEMENT INTEGRATION
  !  BEFORE FORMING ELEMENTAL JACOBIAN AND RHS VECTOR
  !---------------------------------------------------------------------
    TERM_1   = 0.D0
    

    !---------------------------------------------------------------------
  !  ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
  !---------------------------------------------------------------------
   LOOP_GAUSS: DO KK = 1, NGAUSS_1d


  !---------------------------------------------------------------------
  !    CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANSFORMATION
  !    JACOBIAN AT THE GAUSS POINTS IN X,Y COORDINATES
  !---------------------------------------------------------------------

        Prr = 0.d0
        Prz = 0.d0
        Pzz = 0.d0  
    
        do II = 1, nbf_2d
    
          Prr = Prr + Prr_loc(II)*BFN(II,KK)
          Prz = Prz + Prz_loc(II)*BFN(II,KK)
          Pzz = Pzz + Pzz_loc(II)*BFN(II,KK)
    
        enddo
    
    
        CALL BASIS_2d&
        ( KK, Z_loc, R_loc, BFN, DFDC, DFDE, Z, dZdC, dZdE, R, dRdC, dRdE,&
          CJAC, AJAC, DFDR, DFDZ, NGAUSS_1d )
    
      !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n
         
          call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
                                      ned, nr, nz, dL,         &
                                      normalize = .true., forceOnObject = .True.) 
    
    
          WET = WO_1d(KK)*dL
    
    
          TERM_1  = TERM_1 + 2.d0*pi*(nr*(Prz) + (nZ)*(Pzz) )*R*WET
          
         
   ENDDO LOOP_GAUSS

  ! By adding the minus we acquire a positive value for the drag force.
  ! Since it is referred to as drag force, it is already implied that its
  ! direction is opposite to the direction of movement, hence the minus sign
  ! at the result is removed
  dragForce = -TERM_1

 
end Function DragForceCalculation





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
   Integer :: IROW, JCOL

   Real(8)  :: WET
   Real(8)  :: R, dRdC, dRdE
   Real(8)  :: Z, dZdC, dZdE
   Real(8)  :: CJAC, AJAC, JacT, dL, nR, nZ, nl

   Real(8), Dimension(NBF_2d) :: DFDR,  DFDZ
   Real(8), Dimension(NEQ_f)  :: TERM_RES
   
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


     call BASIS_2d&
     ( KK , Z_loc, R_loc, &
       BFN, DFDC , DFDE , &
       Z  , dZdC , dZdE , &
       R  , dRdC , dRdE , &
       CJAC, AJAC, DFDR, DFDZ, NGAUSS_1d)

     !    DEFINE DIFFERENTIAL ARCLENGTH dL & OUTWARD POINTING NORMAL VECTOR n
     call getNormalVectorAtFace( [dzdc, dzde, drdc, drde], &
                                   ned, nr, nz, dL,        &
                                 normalize = .true.)
     WET = WO_1d(KK)*dL
     
     TEMP_RES = TEMP_RES + (  nr * R  + nz * Z  )*R*WET
    
   end do loop_gauss
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   ! Reference: Karapetsas, Photeinos, Dimakopoulos, Tsamopoulos
   ! Dynamics and motion of a gas bubble in a viscoplastic medium under
   ! acoustic excitation, 2019
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   TEMP_RES = - 2.d0* PI * TEMP_RES / 3.d0

 End Function SurfaceIntegration

End Module ExtraEquations


Module FieldFunctions

contains

    function minimumAngleOfTriangle(Solution_, elements) Result(output)
        Use VariableMapping, only: getVariableId
        Use geometry
        Implicit None 
        Real(8), Dimension(:,:), Intent(In) :: Solution_
        Integer, Dimension(:,:), Intent(In) :: elements
                
        Real(8), Dimension(:)  , Allocatable:: output

        Integer                             :: ielem
        Integer, Dimension(:)  , Allocatable:: nm
        Real(8), Dimension(:)  , Allocatable:: x_, X
        Real(8), Dimension(:)  , Allocatable:: y_, Y
        Real(8), Dimension(:)  , Allocatable:: xm_
        Real(8), Dimension(:)  , Allocatable:: ym_

        Integer                             :: nbf_2d
        Integer                             :: nelem
        Integer                             :: nnodes




        if (allocated(output) ) deallocate(output)
        nelem  = size(elements,1)
        nbf_2d = size(elements,2)
        nnodes = size(Solution_,1)

        allocate( X( nnodes) )
        allocate( Y( nnodes) )
        X = Solution_(:, getVariableId("Z"))
        Y = Solution_(:, getVariableId("R"))


        allocate( nm (nbf_2d) )
        allocate( x_ (nbf_2d) )
        allocate( y_ (nbf_2d) )
        allocate( xm_(nbf_2d) )
        allocate( ym_(nbf_2d) )

        allocate( output (nelem) )

        do ielem = 1, nelem
        ! nm = (/1, 2, 3/)
            nm = elements(ielem,:)

            x_ = X     (nm(:))
            y_ = Y     (nm(:))
            

            output(ielem) =  minimumAngle(x_, y_)


        end do

        deallocate( nm  )
        deallocate(  x_ )
        deallocate(  y_ )
        deallocate( xm_ )
        deallocate( ym_ )

    End Function minimumAngleOfTriangle




    
    Function relativeElementArea(Solution_, elements, x_mesh, y_mesh) Result(output)
        Use VariableMapping, only: getVariableId
        Use geometry
        Implicit None 
        Real(8), Dimension(:,:), Intent(In) :: Solution_
        Integer, Dimension(:,:), Intent(In) :: elements
        Real(8), Dimension(:)  , Intent(In) :: x_mesh
        Real(8), Dimension(:)  , Intent(In) :: y_mesh
        
        Real(8), Dimension(:)  , Allocatable:: output

        Integer                             :: ielem
        Integer, Dimension(:)  , Allocatable:: nm
        Real(8), Dimension(:)  , Allocatable:: x_, X
        Real(8), Dimension(:)  , Allocatable:: y_, Y
        Real(8), Dimension(:)  , Allocatable:: xm_
        Real(8), Dimension(:)  , Allocatable:: ym_

        Integer                             :: nbf_2d
        Integer                             :: nelem
        Integer                             :: nnodes
        Real(8)                             :: Area_
        Real(8)                             :: Aream_
        Real(8)                             :: AreaEq_
        Real(8)                             :: AreaEqm_


        if (allocated(output) ) deallocate(output)
        nelem  = size(elements,1)
        nbf_2d = size(elements,2)
        nnodes = size(Solution_,1)

        allocate( X( nnodes) )
        allocate( Y( nnodes) )
        X = Solution_(:, getVariableId("Z"))
        Y = Solution_(:, getVariableId("R"))
        

        allocate( nm (nbf_2d) )
        allocate( x_ (nbf_2d) )
        allocate( y_ (nbf_2d) )
        allocate( xm_(nbf_2d) )
        allocate( ym_(nbf_2d) )

        allocate( output (nelem) )

        do ielem = 1, nelem
        ! nm = (/1, 2, 3/)
            nm = elements(ielem,:)

            x_ = x     (nm(:))
            y_ = y     (nm(:))
            xm_= x_mesh(nm(:))
            ym_= y_mesh(nm(:))

            Area_    = area( x_,  y_)
            Aream_   = area(xm_, ym_)

            AreaEq_  = areaEquilateral( x_,  y_)
            AreaEqm_ = areaEquilateral(xm_, ym_)

           
            output(ielem) =  ( Area_/AreaEq_ - Aream_/AreaEqm_ ) / ( Aream_/AreaEqm_ )


        end do

        deallocate( nm  )
        deallocate(  x_ )
        deallocate(  y_ )
        deallocate( xm_ )
        deallocate( ym_ )
    End Function relativeElementArea



   Function YieldedRegion(Stresses_nodes) Result(output)
        Use VariableMapping, only: getVariableId
        Use PHYSICAL_MODULE, Only: BnN, WiN
        Implicit None 
        Real(8), Dimension(:,:), Intent(In) :: Stresses_nodes
        Real(8), Dimension(:), Allocatable  :: output
        Integer                             :: nodtol_ 

        Integer                             :: j
        Real(8), Dimension(3,3)             :: Deviatoric_Stress_Tensor
        Real(8)                             :: Trace_Stress_Tensor
        Real(8)                             :: Magn_Stress_Dev
        Real(8)                             :: Trace_Deviatoric_Stress_Tensor
        Real(8)                             :: TraceStress
        Real(8)                             :: Magnitude
        Real(8)                             :: Trr, Trz, Tzz, Ttt
        Real(8), Dimension(3,3)             :: unity_tensor, SM, CM, Stress_Tensor, Stress_dot_Stress
        

        unity_tensor = 0.d0
        do j=1,3
            unity_tensor(j,j) = 1.d0
        enddo

        nodtol_ = size(Stresses_nodes,1)
        Allocate( output (nodtol_) )

        do j = 1, nodtol_

            Stress_Tensor = 0.d0

            Stress_Tensor(1,1) = Stresses_nodes(j,1)

            Stress_Tensor(1,2) = Stresses_nodes(j,2)
            Stress_Tensor(2,1) = Stresses_nodes(j,2)

            Stress_Tensor(2,2) = Stresses_nodes(j,3)

            Stress_Tensor(3,3) = Stresses_nodes(j,4)


            Trace_Stress_Tensor = Stress_Tensor(1,1) + Stress_Tensor(2,2) + Stress_Tensor(3,3)

            Deviatoric_Stress_Tensor = Stress_Tensor - (Trace_Stress_Tensor/3.0d0)*Unity_Tensor

            Stress_dot_Stress = matmul(Deviatoric_Stress_Tensor,Deviatoric_Stress_Tensor)

            Magnitude = sqrt( 0.5d0*(Stress_dot_Stress(1,1) + Stress_dot_Stress(2,2) + Stress_dot_Stress(3,3)) )

            output(J)  =  0.d0 !max(0.0d0, Magnitude - BnN)
        end do

    End Function YieldedRegion


    Function dynamicPressure(Solution_) Result(output)
        Use VariableMapping, only: getVariableId
        Use PHYSICAL_MODULE, Only: Pambient_o_Pchar, position, ratio_of_pressures
        Implicit None 
        Real(8), Dimension(:,:), Intent(In) :: Solution_
        Real(8), Dimension(:), Allocatable  :: output
        Integer                             :: nodtol_

        Real(8)                             :: Pressure, Dynamic_Pressure, Z_coord
        Integer                             :: j


        nodtol_ = size(Solution_,1)
        Allocate( output (nodtol_) )


        do j = 1, nodtol_

            Pressure = Solution_(j, getVariableId("P"))

            Z_coord  = Solution_(j, getVariableId("Z"))
            
            Dynamic_Pressure   =   Pressure - ( Pambient_o_Pchar + ratio_of_pressures*(position + Z_coord) )

            output(j)  =  Dynamic_Pressure

        enddo
    end function dynamicPressure


    function fromSgetStresses(Solution_) Result(output)
        Use VariableMapping, only: getVariableId
        Use PHYSICAL_MODULE, only: WiN 
        Implicit None
        Real(8), Dimension(:,:), Intent(In) :: Solution_
        Real(8), Dimension(:,:), Allocatable:: output
        Integer                             :: nodtol_


        Real(8)                             :: Pressure, Dynamic_Pressure, Z_coord
        Real(8), Dimension(3,3)             :: unity_tensor, SM, CM, Stress_Tensor
        Integer                             :: j

        nodtol_ = size(Solution_,1)
        Allocate( output (nodtol_,4) )

        unity_tensor = 0.d0
        do j=1,3
            unity_tensor(j,j) = 1.d0
        enddo

        do j = 1, nodtol_

            SM = 0.d0

            ! SM(1,1) = Solution_(j, getVariableId("Srr"))
            ! SM(1,2) = Solution_(j, getVariableId("Srz"))
            ! SM(2,1) = Solution_(j, getVariableId("Srz"))
            ! SM(2,2) = Solution_(j, getVariableId("Szz"))
            ! SM(3,3) = Solution_(j, getVariableId("Stt"))


            ! CM = matmul(SM,SM)

            ! Stress_Tensor = (CM - unity_tensor)/WiN

            ! output(j,1)  =  Stress_Tensor(1,1)
            ! output(j,2)  =  Stress_Tensor(1,2)
            ! output(j,3)  =  Stress_Tensor(2,2)
            ! output(j,4)  =  Stress_Tensor(3,3)
            output = 0.d0

        enddo
    end function fromSgetStresses

End Module FieldFunctions

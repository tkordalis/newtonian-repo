
Module FieldFunctions

contains

    function minimumAngleOfTriangle(Xnodes, Ynodes, elements) Result(output)
        Use VariableMapping, only: getVariableId
        Use geometry
        Implicit None 
        Real(8), Dimension(:), Intent(In) :: Xnodes
        Real(8), Dimension(:), Intent(In) :: Ynodes
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
        nnodes = size(Xnodes)

        allocate( X( nnodes) )
        allocate( Y( nnodes) )
        ! X = Solution_(:, getVariableId("Z"))
        ! Y = Solution_(:, getVariableId("R"))
        X = Xnodes
        Y = Ynodes



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
        use geometry,        only: trace, secondInvariant

        Implicit None 
        Real(8), Dimension(:,:), Intent(In) :: Stresses_nodes
        Real(8), Dimension(:), Allocatable  :: output
        Integer                             :: nodtol_ 

        Integer                             :: j
        Real(8), Dimension(3,3)             :: I1, Stress_Tensor

        I1 = 0.d0
        do j=1,3
            I1(j,j) = 1.d0
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

            output(J) = sqrt( secondInvariant( Stress_Tensor - (trace(Stress_Tensor)/3.d0)*I1 ) )
        end do

    End Function YieldedRegion


    Function dynamicPressure(time, Solution_) Result(output)
        Use VariableMapping, only: getVariableId
        Use PHYSICAL_MODULE, Only: Pambient_o_Pchar, initial_position, ratio_of_pressures
        use pressure_variation, only: PressureChamber
        Implicit None 
        Real(8), Intent(In)  :: time
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
            
            Dynamic_Pressure   =   Pressure - ( PressureChamber(time) + ratio_of_pressures*(initial_position -  Z_coord) )

            output(j)  =  Dynamic_Pressure

        enddo
    end function dynamicPressure

End Module FieldFunctions

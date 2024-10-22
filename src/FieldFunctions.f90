
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


    function CanalyticLandau(Solution_) Result(output)
        Use VariableMapping, only: getVariableId
        use BoundaryConditions
        Use PHYSICAL_MODULE, only: PeN, pi, KoN
        Use geometry
        use MESH_MODULE, only: Xm, Ym
        use check_for_floating_point_exceptions

        Implicit None 
        Real(8), Dimension(:,:), Intent(In) :: Solution_
                
        Real(8), Dimension(:)  , Allocatable:: output
        Integer                             :: nodtol_ 

        Real(8), Dimension(:)  , Allocatable:: X
        Real(8), Dimension(:)  , Allocatable:: Y
        Real(8)                             :: Svar, heta, ksi, Cao_o_Cinf

        Integer                             :: j
        Integer                             :: nnodes



         if (allocated(output) ) deallocate(output)
        nnodes = size(Solution_,1)

        allocate( X( nnodes) )
        allocate( Y( nnodes) )
        X = Solution_(:, getVariableId("Z"))
        Y = Solution_(:, getVariableId("R"))

        ! X = Xm
        ! Y = Ym

        ! X =  X - bubble%getCentroid()
        X = -X

        Cao_o_Cinf = KoN*bubble%getpressure()

        Allocate( output (nnodes) )


        nnodes = size(Solution_,1)
        do j = 1, nnodes

            ksi = sqrt(X(j)**2+Y(j)**2)
            heta = cos( atan2(Y(j),X(j)) )

            Svar = ( (ksi-1)*sqrt(PeN) * ( 1.d0- heta**2 +1.d-10) ) / ( sqrt(8.d0/3.d0) * (  ( sqrt(2.d0+3.d0*heta - heta**3.d0) +1.d-10) )  )

            output(j)  =  1.d0 - erf(Svar)

            ! call check_fp_exceptions(Svar, 'Svar')
        enddo

        output = Cao_o_Cinf + output*(Cao_o_Cinf-1.d0)


    end function CanalyticLandau


    Function dynamicPressure(Solution_) Result(output)
        Use VariableMapping, only: getVariableId
        Use PHYSICAL_MODULE, Only: Pambient_o_Pchar, initial_position, ratio_of_pressures
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
            
            Dynamic_Pressure   =   Pressure - ( Pambient_o_Pchar + ratio_of_pressures*(initial_position -  Z_coord) )

            output(j)  =  Dynamic_Pressure

        enddo
    end function dynamicPressure

End Module FieldFunctions

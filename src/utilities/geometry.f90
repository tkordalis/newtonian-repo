
module geometry
    private ::  RadiusEquilateral
    public  :: area, AreaEquilateral, distance, minimumAngle, &
                    trace, secondInvariant

    interface distance
        module procedure distance_array
    end interface distance

    interface area
        module procedure area_array
    end interface area

    interface RadiusEquilateral
        module procedure EquilateralRadius
    end interface RadiusEquilateral

    interface areaEquilateral
        module procedure EquilateralArea
    end interface AreaEquilateral


    contains

        function trace(array) Result(output)
            implicit none
            real(8), dimension(:,:), intent(in) :: array
            real(8)         :: output
            integer         :: i

            output = sum( [ (array(i,i),i=1,size(array,1)) ] )

        end function trace

        function secondInvariant(array) Result(output)
            implicit none
            real(8), dimension(:,:), intent(in) :: array
            real(8)                 :: output
            integer                 :: i
            real(8), dimension(size(array,1),size(array,2))  :: array2

            array2 = matmul(array,array) ! array**2

            output = 0.5d0*( (trace(array))**2 + trace(array2) )
        end function secondInvariant

        function minimumAngle(x_nodes, y_nodes) Result(output)
            implicit none
            real(8), dimension(:), Intent(in) :: x_nodes, y_nodes
            real(8) :: output

            Real(8), dimension(2)  :: v21, v31, v12, v32, v13, v23  
            Real(8)                :: v21_length, v31_length, v12_length, v32_length, v13_length, v23_length
            real(8), dimension(3)  :: theta_elem



            v21 = [ x_nodes(2)-x_nodes(1), y_nodes(2)-y_nodes(1) ]
            v31 = [ x_nodes(3)-x_nodes(1), y_nodes(3)-y_nodes(1) ]

            v21_length = sqrt( v21(1)**2 + v21(2)**2 )
            v31_length = sqrt( v31(1)**2 + v31(2)**2 )

            theta_elem(1) = acos( dot_product(v21,v31) / (v21_length*v31_length) )


            v12 = [ x_nodes(1)-x_nodes(2), y_nodes(1)-y_nodes(2) ]
            v32 = [ x_nodes(3)-x_nodes(2), y_nodes(3)-y_nodes(2) ]

            v12_length = sqrt( v12(1)**2 + v12(2)**2 )
            v32_length = sqrt( v32(1)**2 + v32(2)**2 )

            theta_elem(2) = acos( dot_product(v12,v32) / (v12_length*v32_length) )


            v13 = [ x_nodes(1)-x_nodes(3), y_nodes(1)-y_nodes(3) ]
            v23 = [ x_nodes(2)-x_nodes(3), y_nodes(2)-y_nodes(3) ]

            v13_length = sqrt( v13(1)**2 + v13(2)**2 )
            v23_length = sqrt( v23(1)**2 + v23(2)**2 )

            theta_elem(3) = acos( dot_product(v13,v23) / (v13_length*v23_length) )

            output = minval(theta_elem)


        end function minimumAngle




        function distance_array(point1, point2) result(output)
            Implicit None 
            Real(8), Dimension(:), Intent(In) :: point1 
            Real(8), Dimension(:), Intent(In) :: point2
            Real(8)                           :: dist
            Real(8)                           :: output

            Integer                           :: i
            Integer                           :: n

            n = size(point1)
            dist = 0.d0 
            do i = 1, n
                dist = dist + (point1(i) - point2(i))**2
            end do

            output = dsqrt(dist)
        End function distance_array


        Function area_array(x,y) Result(A)
            Implicit None 
            Real(8), Dimension(:), Intent(In) :: x
            Real(8), Dimension(:), Intent(In) :: y
            Real(8)                           :: A

            Real(8)                           :: AB
            Real(8)                           :: BC
            Real(8)                           :: CA
            Real(8)                           :: s


            If ( size(x) /= 3) Stop "[Error] area_array: the input is not triangle vertices."

            AB = distance( [x(1), y(1)], [x(2), y(2)])
            BC = distance( [x(2), y(2)], [x(3), y(3)])
            CA = distance( [x(3), y(3)], [x(1), y(1)])

            ! Herons Formula
            s  = (AB + BC + CA)/2.d0
            A  = dsqrt( s * (s-AB) * (s-BC) * (s-CA) )
        End Function area_array


        Function EquilateralRadius(x,y)   Result(Radius)
            Implicit None 
            Real(8), Dimension(:), Intent(In) :: x
            Real(8), Dimension(:), Intent(In) :: y
            Real(8)                           :: A
            Real(8)                           :: Radius

            Real(8)                           :: AB
            Real(8)                           :: BC
            Real(8)                           :: CA
            Real(8)                           :: s


            AB = distance( [x(1), y(1)], [x(2), y(2)])
            BC = distance( [x(2), y(2)], [x(3), y(3)])
            CA = distance( [x(3), y(3)], [x(1), y(1)])

            ! Herons Formula
            s  = (AB + BC + CA)/2.d0
            A  = dsqrt( s * (s-AB) * (s-BC) * (s-CA) )

            Radius = AB * BC * CA / ( 4.d0 * A )

        End Function EquilateralRadius


        Function EquilateralArea(x,y) Result(AreaEq)
            Implicit None
            Real(8), Dimension(:), Intent(In) :: x
            Real(8), Dimension(:), Intent(In) :: y
            Real(8)                           :: AreaEq
            Real(8)                           :: Radius

            Radius = EquilateralRadius(x,y)

            AreaEq = (3.d0 * sqrt(3.d0) / 4.d0) * Radius**2

        End Function EquilateralArea
        
end module geometry



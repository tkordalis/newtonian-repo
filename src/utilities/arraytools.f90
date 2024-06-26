
!*******************************************************************
! Useful functions for array manipulation
!*******************************************************************

Module ArrayTools
        Private 
        ! Join a sequence of arrays along an existing axis.
        Public :: Concatenate, copyArrayToLocalValues, copy

        Interface Concatenate
            module procedure Concatenate_Integer_dim_1_arrays_2
        End Interface Concatenate

        Interface Copy
            module procedure copy_integer_dim_1
        End Interface Copy

        Interface copyArrayToLocalValues
            module procedure copyArrayToLocalValues_dim1
            module procedure copyArrayToLocalValues_dim2
        End Interface copyArrayToLocalValues


    Contains 



    Subroutine copyArrayToLocalValues_dim1(Array, element, localArray)
        Implicit None 
        Real(8), Dimension(:), Intent(In) :: Array 
        Integer, Dimension(:), Intent(In) :: element 
        Real(8), Dimension(:), Allocatable:: localArray

        Integer                           :: i

        If ( Allocated(localArray) ) Deallocate(localArray) 
        Allocate( localArray ( size(element) ) )

        do i = 1, size(localArray)
        localArray(i) = Array(element(i))
        end do

    End Subroutine copyArrayToLocalValues_dim1


    Subroutine copyArrayToLocalValues_dim2(Array, element, cpid, localArray)
        Implicit None 
        Real(8), Dimension(:,:), Intent(In) :: Array 
        Integer, Dimension(:)  , Intent(In) :: element 
        Integer                , Intent(In) :: cpid
        Real(8), Dimension(:,:), Allocatable:: localArray
        
        Integer                             :: n1
        Integer                             :: n2
        Integer                             :: nel
        Integer                             :: i
        Integer                             :: j

        n1  = size(Array,1)
        n2  = size(Array,2)
        nel = size(element)

        If ( Allocated(localArray) ) Deallocate(localArray) 


        if ( cpid == 1 ) then 
            Allocate ( localArray( nel, n2 ) )
        
            do i = 1, nel
                do j = 1, n2
                    localArray(i,j) = Array(element(i), j )
                end do
            end do

        else 
            Print*, "The case has not implemented yes. copyArrayToLocalValues_dim2"
            Stop
        end if


    End Subroutine copyArrayToLocalValues_dim2

    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Private functions 
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Subroutine Concatenate_Integer_dim_1_arrays_2 (array1, array2, outputarray)
        Implicit None 
        Integer, Dimension(:), Intent(In)  :: array1 
        Integer, Dimension(:), Intent(In)  :: array2
        Integer, Dimension(:), Allocatable :: outputarray


        If ( Allocated(outputarray) ) Deallocate(outputarray)

        Allocate( outputarray( size(array1) + size(array2) ) )

        outputarray(1:size(array1))   = array1
        outputarray(size(array1)+1: ) = array2

    End Subroutine Concatenate_Integer_dim_1_arrays_2


    Subroutine copy_integer_dim_1( array1, array2)
        Implicit None
        Integer, Dimension(:), Intent(In)  :: array1
        Integer, Dimension(:), Allocatable :: array2


        if ( Allocated(array2) ) Deallocate(array2)
        Allocate(array2 (size(array1) ) )
        
        ! for the best results use the blas subs
        array2  = array1
    End Subroutine copy_integer_dim_1

    



End Module ArrayTools

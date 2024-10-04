    Function getdVtankdt(This)  Result(output)
        use time_integration, only: dt
        Implicit None 
        Class(BubbleDiffusionStaticCS)       :: This
        Real(8)             :: volume
        Real(8)             :: output


        output = ( this%volume - this%Volume_o ) / dt
    End Function getdVtankdt


    Function getDragForce(this) Result(output)
        Implicit none
        Class(BubbleDiffusionStaticCS)                               :: this
        Real(8)                                     :: output

        Real(8)                                     :: DragForce

        DragForce = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, DragForceCalculation)

        output = DragForce
    end Function getDragForce

    Function getAspectRatio(this)    Result(AR)
        Use GLOBAL_ARRAYS_MODULE,      Only: TL
        Implicit none
        Class(BubbleDiffusionStaticCS), Intent(In)            :: this
        Real(8)                              :: output

        Real(8), Dimension(:)  , Allocatable :: Z_coord, R_coord
        Real(8)                              :: Height, Width, AR
        
        call getBoundaryNodes(TL, This%elements, This%faces, "Z", Z_coord)
        call getBoundaryNodes(TL, This%elements, This%faces, "R", R_coord)
        
        Height = maxval(Z_coord) - minval(Z_coord)
        Width  = maxval(R_coord) ! - minval(R_coord) -> commented out since it is always zero

        AR = Height / (2.d0*Width)

        output = AR
    End Function getAspectRatio


    Function volumeConservation(this) Result(output)
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: Volume
        
        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        
        output = Volume - this%InitialVolume
        
        call loopOverElements(this%nelem, this%elements, this%faces, this%gidP, SurfaceIntegration   ) ! first  constrain
    end Function volumeConservation


    Function PressureVolumeConservation(this) Result(output)
        Implicit None
        Class(BubbleDiffusionStaticCS) :: this
        Real(8)       :: output

        Real(8)       :: Volume
        

        Volume = integrateOverAllElementsOfTheBoundary (this%elements, this%faces, SurfaceIntegration )
        
        output = this%pressure * Volume - this%InitialPressure * this%InitialVolume
     
        call loopOverElements(this%nelem, this%elements, this%faces, this%gidP, SurfaceIntegration, this%pressure ) 
        
    end Function PressureVolumeConservation
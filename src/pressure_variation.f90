module pressure_variation
    use physical_module
    use time_integration
    implicit none

    private :: transientPressureChange
    public ::  PressureChamber, areConditionsSteady
    Real(8), dimension(6), parameter :: crit_times = [4000.5d0, &
                                                    24.4d0,  &
                                                    4100.d0,  & 
                                                    6100.d0,  &
                                                    8100.d0,  &
                                                    10100.d0 ]
                                                     
 
    
    Real(8), dimension(7), parameter :: Pressure_before_crit_times = [Pambient/gravity_stress,    &
                                                                      Pambient/gravity_stress, &
                                                                      Pambient/gravity_stress, &
                                                                      20000.d0/gravity_stress, &
                                                                      17000.d0/gravity_stress, &
                                                                      14000.d0/gravity_stress, &
                                                                      12000.d0/gravity_stress]


    ! Real(8), parameter       :: deltat_transientPressure = 60.d0
    Real(8), parameter       :: deltat_transientPressure = 30.d0
    integer                  :: crit_times_counter, pressure_change_counter
    logical                  :: isThePressureTransient = .false.
    logical                  :: rampUpDt   = .false.
    logical                  :: rampDownDt = .false.
    Real(8)                  :: time_for_dto, time_for_dt1


    contains

    function PressureChamber(time) result(P_t)
        implicit none
        real(8) :: time, P_t
        real(8) :: time_previous_P, time_next_P, previous_P, next_P

        call transientPressureChange(time)

        if (isThePressureTransient .eqv. .true.) then
            ! crit_times_counter has already been given a value from function transientPressureChange
            ! it starts from value 1 and increases after each value has been reached
            ! it can be translated as: towards which crit_time is the problem heading to
            time_previous_P = crit_times(crit_times_counter - 1)
            time_next_P     = time_previous_P + deltat_transientPressure
            previous_P = Pressure_before_crit_times(pressure_change_counter - 1)
            next_P     = Pressure_before_crit_times(pressure_change_counter    )
            
            P_t = previous_P * ( time - time_next_P    )  / ( time_previous_P - time_next_P ) &
                + next_P     * ( time - time_previous_P ) / ( time_next_P - time_previous_P )
       elseif (isThePressureTransient .eqv. .false.) then
            P_t = Pressure_before_crit_times(pressure_change_counter)
        endif

    end function PressureChamber

! -----------------------------------------------------

    subroutine transientPressureChange(time)
        implicit none
        real(8) :: time, crit_time
        integer :: count_dummy


        do count_dummy = 1, size(crit_times)
            crit_time = crit_times(count_dummy)
            if ( time .gt. crit_times(size(crit_times)) ) then
                pressure_change_counter = size(crit_times) + 1
                crit_times_counter      = size(crit_times) + 1
                exit
            elseif ( time .lt. crit_time ) then
                crit_times_counter = count_dummy
                pressure_change_counter   = crit_times_counter
                exit
            endif
        enddo

        if (crit_times_counter .ne. 1) then
          crit_time = crit_times(crit_times_counter-1)
          if ( ( time .gt. crit_time ) .and. ( time .lt. (crit_time + deltat_transientPressure) ) ) then
                isThePressureTransient = .true.
          else
                isThePressureTransient = .false.
          endif
        endif

    end subroutine transientPressureChange

! -----------------------------------------------------

    subroutine areConditionsSteady(time)
        implicit none
        real(8), Intent(in):: time
        real(8)            :: fine_timestep_time
        real(8)            :: time_previous_dt, time_next_dt, previous_dt, next_dt
        real(8)            :: initial_time_period
        real(8)            :: before_critTime_time_period
        real(8)            :: rampUpDt_time_period
        real(8)            :: rampDownDt_time_period
        
        initial_time_period         = 0.68d0
        rampUpDt_time_period        = 5.d0 ! regulates the slope of the linear function increasing timestep
        before_critTime_time_period = 0.5d0
        rampDownDt_time_period      = 2.d0 ! regulates the slope of the linear function decreasing timestep
        fine_timestep_time          = 8.d0


        if ( crit_times_counter .eq. 1 ) then 
            if ( abs(time - initial_time_period) .lt. 0.499d0*dt ) then
                rampUpDt        = .true.
                time_for_dto = time
                time_for_dt1 = time + rampUpDt_time_period
            elseif ( abs(time - crit_times(crit_times_counter) + rampDownDt_time_period + before_critTime_time_period) .lt. 0.499d0*dt ) then
                rampDownDt     = .true.
                time_for_dto = time
                time_for_dt1 = time + rampDownDt_time_period
            endif
        else
            if ( abs(time - crit_times(crit_times_counter-1) - fine_timestep_time ) .lt. 0.499d0*dt ) then

                rampUpDt        = .true.
                time_for_dto = time
                time_for_dt1 = time + rampUpDt_time_period
            elseif ( abs(time - crit_times(crit_times_counter) + rampDownDt_time_period + before_critTime_time_period) .lt. 0.499d0*dt ) then

                rampDownDt     = .true.
                time_for_dto = time
                time_for_dt1 = time + rampDownDt_time_period
            endif
        endif
        
        if (rampUpDt) then
            time_previous_dt = time_for_dto
            time_next_dt     = time_for_dt1
            previous_dt      = Dt_constant
            next_dt          = Dt_max

            DTb = Dto
            Dto = DT
            ! DT = previous_dt * ( time - time_next_dt    )  / ( time_previous_dt - time_next_dt ) &
            ! + next_dt    * ( time - time_previous_dt ) / ( time_next_dt - time_previous_dt )
            DT = previous_dt * ( log10(time) - log10(time_next_dt)    )  / ( log10(time_previous_dt) - log10(time_next_dt) ) &
            + next_dt    * ( log10(time) - log10(time_previous_dt) ) / ( log10(time_next_dt) - log10(time_previous_dt) )

            if (dt .gt. dt_max) then
                dt = dt_max
                rampUpDt = .false.
            endif
        elseif (rampDownDt) then
            time_previous_dt = time_for_dto
            time_next_dt     = time_for_dt1
            previous_dt      = Dt_max
            next_dt          = Dt_constant

            DTb = Dto
            Dto = DT
            DT = previous_dt * ( log10(time) - log10(time_next_dt)    )  / ( log10(time_previous_dt) - log10(time_next_dt) ) &
            + next_dt    * ( log10(time) - log10(time_previous_dt) ) / ( log10(time_next_dt) - log10(time_previous_dt) )

            if (dt .lt. dt_constant) then
                dt = dt_constant
                rampDownDt = .false.
            endif
        endif
    end subroutine areConditionsSteady


end module pressure_variation
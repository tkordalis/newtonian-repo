module pressure_variation
	use physical_module
    use time_integration
	implicit none
	Real(8), dimension(6), parameter :: crit_times = [100.d0, &
                                                    2100.d0,  &
                                                    4100.d0,  & 
                                                    6100.d0,  &
                                                    8100.d0,  &
                                                    10100.d0 ]
                                                     
 
    
    Real(8), dimension(7), parameter :: Pressure_before_crit_times = [Patm/gravity_stress,    &
                                                                      32000.d0/gravity_stress, &
                                                                      27000.d0/gravity_stress, &
                                                                      19000.d0/gravity_stress, &
                                                                      16000.d0/gravity_stress, &
                                                                      13000.d0/gravity_stress, &
                                                                      11000.d0/gravity_stress]
    
    integer                  :: crit_times_counter, pressure_change_counter
    logical                  :: isThePressureTransient = .false.
    logical                  :: rampUpDt   = .false.
    logical                  :: rampDownDt = .false.
    Real(8)                  :: time_for_dto, time_for_dt1


    contains

    function PressureChamber(time) result(P_t)
        implicit none
        real(8) :: time, crit_time, P_t
        real(8) :: time_previous_P, time_next_P, previous_P, next_P
        integer :: count

        call transientPressureChange(time)

        if (isThePressureTransient .eqv. .true.) then
            ! crit_times_counter has already been given a value from function transientPressureChange
            ! it starts from value 1 and increases after each value has been reached
            ! it can be translated as: towards which crit_time is the problem heading to
            time_previous_P = crit_times(crit_times_counter - 1)
            time_next_P     = time_previous_P + 60.d0
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
          if ( ( time .gt. crit_time ) .and. ( time .lt. (crit_time + 60.d0) ) ) then
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
        
        initial_time_period         = 60.d0
        rampUpDt_time_period        = 60.d0 ! regulates the slope of the linear function increasing timestep
        before_critTime_time_period = 10.d0
        rampDownDt_time_period      = 30.d0 ! regulates the slope of the linear function decreasing timestep
        fine_timestep_time          = 70.d0

        if ( crit_times_counter .eq. 1 ) then 
            if ( abs(time - initial_time_period) .lt. 1.d-2 ) then
                rampUpDt        = .true.
                time_for_dto = time
                time_for_dt1 = time + rampUpDt_time_period
            elseif ( abs(time - crit_times(crit_times_counter) + rampDownDt_time_period + before_critTime_time_period) .lt. 0.4d0 ) then
                rampDownDt     = .true.
                time_for_dto = time
                time_for_dt1 = time + rampDownDt_time_period
            endif
        else
            if ( abs(time - crit_times(crit_times_counter-1) - fine_timestep_time ) .lt. 0.02d0 ) then
                rampUpDt        = .true.
                time_for_dto = time
                time_for_dt1 = time + rampUpDt_time_period
            elseif ( abs(time - crit_times(crit_times_counter) + rampDownDt_time_period + before_critTime_time_period) .lt. 0.5d0 ) then
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
              DT = previous_dt * ( time - time_next_dt    )  / ( time_previous_dt - time_next_dt ) &
                  + next_dt    * ( time - time_previous_dt ) / ( time_next_dt - time_previous_dt )
              
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
              DT = previous_dt * ( time - time_next_dt    )  / ( time_previous_dt - time_next_dt ) &
                  + next_dt    * ( time - time_previous_dt ) / ( time_next_dt - time_previous_dt )

              if (dt .lt. dt_constant) then
                  dt = dt_constant
                  rampDownDt = .false.
              endif
          endif

    end subroutine areConditionsSteady


end module pressure_variation
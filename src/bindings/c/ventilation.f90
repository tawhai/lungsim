module ventilation_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine evaluate_vent_c(gdirn, num_brths, num_itns, chest_wall_compliance, &
       dt, err_tol, press_in, volume_target, Tinsp, Texpn, &
       expiration_type, expiration_type_len) bind(C, name="evaluate_vent_c")

    use arrays,only: dp
    use iso_c_binding, only: c_ptr
    use other_consts, only: MAX_STRING_LEN
    use utils_c, only: strncpy
    use ventilation, only: evaluate_vent
    implicit none

    integer,intent(in) :: gdirn,num_brths,num_itns,expiration_type_len
    real(dp),intent(in) :: chest_wall_compliance,dt,err_tol,press_in,volume_target, &
         Tinsp,Texpn
    type(c_ptr), value, intent(in) :: expiration_type
    character(len=MAX_STRING_LEN) :: expiration_type_f

    call strncpy(expiration_type_f, expiration_type, expiration_type_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_evaluate_vent(gdirn, num_brths, num_itns, chest_wall_compliance, &
       dt, err_tol, press_in, volume_target, Tinsp, Texpn, expiration_type_f)
#else
    call evaluate_vent(gdirn, num_brths, num_itns, chest_wall_compliance, &
       dt, err_tol, press_in, volume_target, Tinsp, Texpn, expiration_type_f)
#endif

  end subroutine evaluate_vent_c


  !###################################################################################

  subroutine evaluate_uniform_flow_c() bind(C, name="evaluate_uniform_flow_c")

    use ventilation, only: evaluate_uniform_flow
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_evaluate_uniform_flow
#else
    call evaluate_uniform_flow
#endif

  end subroutine evaluate_uniform_flow_c


!###################################################################################

  subroutine two_unit_test_c() bind(C, name="two_unit_test_c")
    use ventilation, only: two_unit_test
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_two_unit_test
#else
    call two_unit_test
#endif

  end subroutine two_unit_test_c

!###################################################################################
end module ventilation_c

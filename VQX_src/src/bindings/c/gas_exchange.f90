module gas_exchange_c

  use precision
  
  implicit none
  
  private
  
contains
!!!######################################################################

  subroutine initial_gasexchange_c(initial_concentration, &
       inlet_concentration,p_ven_o2,shunt_fraction,cardiac_output, &
       surface_area,V_cap,species,species_len) bind(C, name="initial_gasexchange_c")
    use gas_exchange, only: initial_gasexchange
    use iso_c_binding, only: c_ptr
    use other_consts, only: MAX_FILENAME_LEN
    use utils_c, only: strncpy
    implicit none

    integer,intent(in) :: species_len
    real(dp),intent(in) :: initial_concentration,inlet_concentration,&
         p_ven_o2,shunt_fraction,cardiac_output
    real(dp),intent(in) :: surface_area,V_cap
    type(c_ptr),value,intent(in) :: species
    character(len=MAX_FILENAME_LEN) :: species_f

    call strncpy(species_f, species, species_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_initial_gasexchange(initial_concentration, &
       inlet_concentration,p_ven_o2,shunt_fraction,cardiac_output, &
       surface_area,V_cap,species_f)
#else
    call initial_gasexchange(initial_concentration, &
       inlet_concentration,p_ven_o2,shunt_fraction,cardiac_output, &
       surface_area,V_cap,species_f)
#endif

  end subroutine initial_gasexchange_c

!!! ######################################################################
  
  subroutine steadystate_gasexchange_c(cardiac_output,Vdot_alv,c_art_o2,c_ven_o2,&
       p_art_co2,p_art_o2,p_i_o2,p_ven_co2,p_ven_o2,shunt_fraction,&
       VCO2,VO2) bind(C, name="steadystate_gasexchange_c")
    use gas_exchange, only: steadystate_gasexchange
    implicit none
    
!!! Parameter List
    real(dp),intent(in) :: cardiac_output,Vdot_alv,p_i_o2,shunt_fraction,VCO2,VO2
    real(dp), intent(inout) :: c_art_o2,c_ven_o2,p_art_co2,p_art_o2,p_ven_o2,p_ven_co2
    
#if defined _WIN32 && defined __INTEL_COMPILER
    call so_steadystate_gasexchange(cardiac_output,Vdot_alv,c_art_o2,c_ven_o2,&
         p_art_co2,p_art_o2,p_i_o2,p_ven_co2,p_ven_o2,shunt_fraction,&
         VCO2,VO2)
#else
    call steadystate_gasexchange(cardiac_output,Vdot_alv,c_art_o2,c_ven_o2,&
         p_art_co2,p_art_o2,p_i_o2,p_ven_co2,p_ven_o2,shunt_fraction,&
         VCO2,VO2)
#endif
    
  end subroutine steadystate_gasexchange_c

!!! ######################################################################
  
  subroutine set_perfusion_gradient_c(Gdirn,cardiac_output,COV,Qmax,Qmin) &
       bind(C, name="set_perfusion_gradient_c")
    use gas_exchange, only: set_perfusion_gradient
    implicit none

    integer,intent(in) :: Gdirn
    real(dp),intent(in) :: cardiac_output,COV,Qmax,Qmin

#if defined _WIN32 && defined __INTEL_COMPILER
    call set_perfusion_gradient(Gdirn,cardiac_output,COV,Qmax,Qmin)
#else
    call set_perfusion_gradient(Gdirn,cardiac_output,COV,Qmax,Qmin)
#endif
    
  end subroutine set_perfusion_gradient_c

!!! ######################################################################
  
end module gas_exchange_c


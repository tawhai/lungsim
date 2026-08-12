module lymphatics
    !*Brief Description:* This module contains all lymphatic-specific subroutines
    !
  !*LICENSE:*
  !
  !Test Test Test
  !
  !*Full Description:*
  !
  !This module contains code for pulmonary fluid flux within the alveolo-capillary network,
  !and lymph transport through lymphatic collecting vessels.


  ! Supplementary equations for the original paper: https://www.protocols.io/view/supplementary-material-an-in-silico-approach-to-un-kxygx9zmwg8j/v1

  use arrays
  use diagnostics
  use indices
  use other_consts
  use precision ! sets dp for precision

  implicit none

  !Module parameters

  ! Baseline value parameters (eventually will be user-defined?)
  integer,protected :: sex !,n_timesteps
  !sex (0 = male, 1 = female) only determines the weight and therefore size of the lung. Should be updated based on CT

  integer :: printcount

  real(dp),protected :: lung_mass,capillary_volume_raw

  ! Capillary parameters
  real(dp),parameter :: capillary_conductivity = 4.41335e-8_dp*2.0_dp!*0.753!*2.0_dp! 0.753!5.98e-6!uL/s/mmHg/mm²   4.41335e-8 !mL.s-1.mmHg-1  obtained from Parker (6e-8 cm H2O)

  real(dp),parameter :: open_capillaries = 1.0_dp/6.0_dp !based on open capillaries at rest. Should be solved for by perfusion model?

  real(dp),protected :: interstitial_capacity
  real(dp),protected ::    interstitial_capacity_a  !arbitrarily sized - needs further studies on the capillary-lymph interface
  real(dp),protected ::   interstitial_capacity_b  !arbitrarily sized - needs further studies on the capillary-lymph interface
  real(dp),protected :: capillary_volume

  ! Simulation parameters
  real(dp),protected :: breathing_rate !constant but should be imported directly from ventilation model
  real(dp),protected :: breathing_function

  ! These two would presumably change in a geometrically consistent lymphatics model
  real(dp),parameter ::  int_max = -1.0_dp
  real(dp),parameter ::  int_min = -8.0_dp
  real(dp),parameter ::  int_diff = int_min - int_max! -8.00_dp - (-1.00_dp) ! intPmin - intPmax in mmHg
  real(dp),parameter ::  lymph_max = 1.0_dp! 1.00_dp - (-8.00_dp) ! lymphPmax-lymphPmin in mmHg
  real(dp),parameter ::  lymph_min = -8.0_dp! 1.00_dp - (-8.00_dp) ! lymphPmax-lymphPmin in mmHg
  real(dp),parameter ::  lymph_diff = lymph_max-lymph_min! 1.00_dp - (-8.00_dp) ! lymphPmax-lymphPmin in mmHg

  real(dp), allocatable ::  interstitial_pressure_a (:,:)
  real(dp), allocatable ::  interstitial_pressure_b (:,:)
  real(dp),  allocatable :: flux_a(:,:)
  real(dp),  allocatable :: flux_b(:,:)
  real(dp),  allocatable :: osm_flux(:,:)
  real(dp), allocatable :: lym_condition(:)

  real(dp), allocatable :: sats(:,:)

  real(dp), allocatable :: P_initial_lymphtix(:,:)
  real(dp), allocatable :: total_hydro_flux (:,:)
  real(dp), allocatable :: initial_lymph_flow(:,:)
  real(dp), allocatable :: initial_lymph_volume(:,:)
  real(dp), allocatable :: interstitial_volume(:,:)
  real(dp), allocatable :: interstitial_volume_a(:,:)
  real(dp), allocatable :: interstitial_volume_b(:,:)
  real(dp), allocatable :: interstitial_saturation(:,:)
  real(dp), allocatable :: int_osm_n(:,:)
  real(dp), allocatable :: alveolar_overflow(:)

  real(dp), allocatable ::  alveolar_volume(:,:)
  real(dp), allocatable :: initial_osm_n(:,:)
  real(dp), allocatable :: total_osm_flux(:,:)
  real(dp), dimension(:), allocatable :: unit_active_time

  ! protein parameters
  real(dp),parameter :: sigma = 0.6_dp*0.4!0.7!*0.4!0.2_dp !"Rat lung venules Lp=4.4×10⁻⁷; matches Safdar exactly" /Pulmonary σ=0.62 total protein; lung most consistent with lymph data"/Safdar/JCI (2003)/Parker/AJPLung (2006)
  real(dp),parameter :: Gp = 4.5e-7_dp!ul/s/mm2 with diameter 0.008 mm  1.13e-11_dp mu L/ms/mm
  real(dp),parameter :: R_contamination = 0.01555_dp
  real(dp),parameter :: c_plasma_baseline = 70.0_dp ! mg/ml = ug/ul
  real(dp),parameter :: c_interstitial_baseline = 45.0_dp
  real(dp),parameter :: volume_threshold = 1.0e-10_dp
  
  ! protein state variables (per unit)
  real(dp) :: V_plasma_unit
  real(dp), allocatable :: Q_plasma(:,:)
  real(dp), allocatable :: Q_int_a(:,:)
  real(dp), allocatable :: Q_int_b(:,:)
  real(dp), allocatable :: c_plasma(:,:)
  real(dp), allocatable :: c_int_a(:,:)
  real(dp), allocatable :: c_int_b(:,:)
  real(dp), allocatable :: osm_cap(:,:)
  real(dp), allocatable :: osm_int_a(:,:)
  real(dp), allocatable :: osm_int_b(:,:)
  real(dp), allocatable :: Jp_cap_a(:,:)
  real(dp), allocatable :: Jp_cap_b(:,:)
  real(dp), allocatable :: Jp_lymph_b(:,:)
  real(dp), allocatable :: Jp_diffusive(:,:)
  real(dp), allocatable :: Jp_convective(:,:)

  ! whether to printout the alveolar flux results
  logical,parameter :: write_out=.false.
  !Module types

  !Module variables

  integer :: n_active
  integer, allocatable :: active(:)
  
  !Interfaces
!  private
  public alveolar_volume
  public alveolar_flux
  public lymphatic_transport
  public alveolar_flux_dt

contains

!!!#############################################################################
  
  subroutine alveolar_flux_dt(dt, time, T_interval)

    real(dp), intent(in) :: dt,time, T_interval
    ! Local variables
    integer :: i,nunit, n_nunit, count,fluid_steps, new_n_active
    integer, allocatable :: active_copy(:)
    real(dp) :: capillary_osm_n, cap_osm_conc,diffusion,excess,flux_c, &
         initial_lymph_conc,interstitial_osmotic, int_osm_conc,lymph_conductivity,&
         net_flux,overflow,sin_breath,sumuptake,test_time, total_flux,transit_time,capillary_SA 
    real(dp) :: capillary_pressure, P_elastic, diff_Pe,fluctuation, fluid_dt !,P_elastic_pre
    ! protein variables:
    real(dp) :: c_overflow_a, Jq_diffusion, Q_overflow_alv, Q_overflow_b, hydrostatic_gradient ,&
         effective_gradient, sat, osmotic_reduction_factor, polynomial_factor
    logical :: cont
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    do nunit = 1, num_units
       call alveolar_flux_unit_dt(nunit, dt, time, T_interval)
    enddo

  end subroutine alveolar_flux_dt

!!!#############################################################################
  
  subroutine alveolar_flux_unit_dt(nunit, dt, time, T_interval)

    integer, intent(in) :: nunit
    real(dp), intent(in) :: dt,time, T_interval
    ! Local variables
    integer :: i,n_nunit, count,fluid_steps, new_n_active
    integer, allocatable :: active_copy(:)
    real(dp) :: capillary_osm_n, cap_osm_conc,diffusion,excess,flux_c, &
         initial_lymph_conc,interstitial_osmotic, int_osm_conc,lymph_conductivity,&
         net_flux,overflow,sin_breath,sumuptake,test_time, total_flux,transit_time,capillary_SA 
    real(dp) :: capillary_pressure, P_elastic, diff_Pe,fluctuation, fluid_dt !,P_elastic_pre
    ! protein variables:
    real(dp) :: c_overflow_a, Jq_diffusion, Q_overflow_alv, Q_overflow_b, hydrostatic_gradient ,&
         effective_gradient, sat, osmotic_reduction_factor, polynomial_factor
    logical :: cont
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    fluid_steps = 2
    fluid_dt =dt/fluid_steps
    
    ! only place that 'time' is used
    sin_breath = sin(2.0_dp*pi*0.25_dp*time)
    
    unit_active_time(nunit) = unit_active_time(nunit) + fluid_dt
    
    ! ============================================
    ! STEP 1: PRESSURES
    ! ============================================
    
    capillary_pressure = ((unit_field(nu_blood_press,nunit))/133.32239_dp)*2.0_dp
    fluctuation = ((unit_field(nu_Pe_max,nunit)/133.32239_dp)-(unit_field(nu_Pe_min,nunit)/133.32239_dp))
    
    interstitial_volume(nu_intsat, nunit)= interstitial_volume_a(nu_intsat, nunit)+ interstitial_volume_b(nu_intsat, nunit)
    interstitial_saturation(nu_intsat, nunit)= interstitial_volume(nu_intsat, nunit)/ interstitial_capacity
    
    interstitial_pressure_a(nu_Pe,nunit)= fluctuation/2.0_dp * sin_breath + &
         (int_diff +fluctuation) * (interstitial_volume_a (nu_intsat, nunit)/ interstitial_capacity_a)**2.0_dp + &
         (int_diff +fluctuation)*(-2.0_dp) * (interstitial_volume_a(nu_intsat, nunit) / interstitial_capacity_a) + &
         (int_min +fluctuation/2.0_dp)
    
    interstitial_pressure_b (nu_Pe,nunit)= fluctuation/2.0_dp * sin_breath + &
         (int_diff +fluctuation) * (interstitial_volume_b(nu_intsat, nunit) / interstitial_capacity_b)**2.0_dp + &
         (int_diff +fluctuation)*(-2.0_dp) * (interstitial_volume_b(nu_intsat, nunit) / interstitial_capacity_b) + &
         (int_min+fluctuation/2.0_dp)
    
    ! ============================================
    ! STEP 2: PROTEIN CONCENTRATIONS
    ! ============================================
    !       if(V_plasma_unit> volume_threshold) then
    !          c_plasma(nu_osmflux, nunit) = Q_plasma(nu_osmflux, nunit) / V_plasma_unit
    !       else
    c_plasma(nu_osmflux, nunit) = c_plasma_baseline
    !       endif
    if (interstitial_volume_a(nu_intsat, nunit) < volume_threshold) then
       interstitial_volume_a(nu_intsat, nunit) = volume_threshold
    endif
    
    if (interstitial_volume_b(nu_intsat, nunit) < volume_threshold) then
       interstitial_volume_b(nu_intsat, nunit) = volume_threshold
    endif

    if(interstitial_volume_a(nu_intsat, nunit) > volume_threshold) then
       c_int_a(nu_osmflux, nunit) = Q_int_a(nu_osmflux, nunit) / interstitial_volume_a(nu_intsat, nunit)
    else
       c_int_a(nu_osmflux, nunit) = c_interstitial_baseline
    endif
    
    if(interstitial_volume_b(nu_intsat, nunit) > volume_threshold) then
       c_int_b(nu_osmflux, nunit) = Q_int_b(nu_osmflux, nunit) / interstitial_volume_b(nu_intsat, nunit)
    else
       c_int_b(nu_osmflux, nunit) = c_interstitial_baseline
    endif
    
    ! ============================================
    ! STEP 3: OSMOTIC PRESSURES
    ! ============================================
    osm_cap(nu_osmflux, nunit) = 0.157_dp * c_plasma(nu_osmflux, nunit) + &
         0.0032_dp * c_plasma(nu_osmflux, nunit)**2
    osm_int_a(nu_osmflux, nunit) = 0.157_dp * c_int_a(nu_osmflux, nunit) + &
         0.0032_dp * c_int_a(nu_osmflux, nunit)**2
    osm_int_b(nu_osmflux, nunit) = 0.157_dp * c_int_b(nu_osmflux, nunit) + &
         0.0032_dp * c_int_b(nu_osmflux, nunit)**2
    
    !       ! ============================================
    !       ! STEP 4: OSMOTIC REDUCTION FACTOR (NOW IT'S SAFE!)
    !       ! ============================================
    !       hydrostatic_gradient = capillary_pressure - interstitial_pressure_b(nu_Pe,nunit)
    
    !
    !       if(abs(hydrostatic_gradient) > 0.1_dp) then
    !          osmotic_reduction_factor = (hydrostatic_gradient - osmotic_gradient) / hydrostatic_gradient
    !          osmotic_reduction_factor = max(0.01_dp, min(1.0_dp, osmotic_reduction_factor))
    !       else
    !          osmotic_reduction_factor = 0.5_dp
    !       endif
    
    ! ============================================
    ! STEP 5: FLUXES (using calculated osmotic pressures)
    ! ============================================
    if(capillary_pressure > interstitial_pressure_a(nu_Pe,nunit))then
       flux_a(nu_av_flux,nunit) = 0.5_dp * capillary_conductivity * unit_field(nu_sa,nunit) * &
            ((capillary_pressure - interstitial_pressure_a(nu_Pe,nunit)) - &
            sigma * (osm_cap(nu_osmflux,nunit) - osm_int_a(nu_osmflux,nunit))) * fluid_dt
    else
       flux_a(nu_av_flux,nunit) = 0.0_dp
    endif
    
    if(capillary_pressure > interstitial_pressure_b(nu_Pe,nunit))then
       flux_b(nu_av_flux,nunit) = 0.5_dp * capillary_conductivity * unit_field(nu_sa,nunit) * &
            ((capillary_pressure - interstitial_pressure_b(nu_Pe,nunit)) - &
            sigma * (osm_cap(nu_osmflux,nunit) - osm_int_b(nu_osmflux,nunit))) * fluid_dt
    else
       flux_b(nu_av_flux,nunit) = 0.0_dp
    endif
    
    flux_c = flux_a(nu_av_flux,  nunit) + flux_b(nu_av_flux,  nunit)
    total_hydro_flux (nu_flux,nunit) = total_hydro_flux (nu_flux,nunit) + flux_c
    unit_field(nu_osmflux,  nunit) = sigma * (osm_cap(nu_osmflux,nunit) - osm_int_a(nu_osmflux,nunit)) + &
         sigma * (osm_cap(nu_osmflux,nunit) - osm_int_b(nu_osmflux,nunit))
    
    ! Calculate protein fluxes to IntA
    Jp_diffusive(nu_osmflux,nunit) = Gp * (unit_field(nu_sa,nunit) * 0.5_dp) * &
         (c_plasma(nu_osmflux,nunit) - c_int_a(nu_osmflux,nunit)) * fluid_dt
    
    if(flux_a(nu_av_flux,nunit) > 0.0_dp) then
       Jp_convective(nu_osmflux,nunit) = flux_a(nu_av_flux,nunit) * &
            c_plasma(nu_osmflux,nunit) * R_contamination
    else
       Jp_convective(nu_osmflux,nunit) = flux_a(nu_av_flux,nunit) * &
            c_int_a(nu_osmflux,nunit) * R_contamination
    endif
    
    Jp_cap_a(nu_osmflux,nunit) = Jp_diffusive(nu_osmflux,nunit) + Jp_convective(nu_osmflux,nunit)
    
    ! Calculate protein fluxes to IntB
    Jp_diffusive(nu_osmflux,nunit) = Gp * (unit_field(nu_sa,nunit) * 0.5_dp) * &
         (c_plasma(nu_osmflux,nunit) - c_int_b(nu_osmflux,nunit)) * fluid_dt
    
    if(flux_b(nu_av_flux,nunit) > 0.0_dp) then
       Jp_convective(nu_osmflux,nunit) = flux_b(nu_av_flux,nunit) * &
            c_plasma(nu_osmflux,nunit) * R_contamination
    else
       Jp_convective(nu_osmflux,nunit) = flux_b(nu_av_flux,nunit) * &
            c_int_b(nu_osmflux, nunit) * R_contamination
    endif
    
    Jp_cap_b(nu_osmflux,nunit) = Jp_diffusive(nu_osmflux,nunit) + Jp_convective(nu_osmflux,nunit)
    
    ! Update protein amounts in interstitium
    Q_int_a(nu_osmflux,nunit) = Q_int_a(nu_osmflux,nunit) + Jp_cap_a(nu_osmflux,nunit)
    Q_int_b(nu_osmflux,nunit) = Q_int_b(nu_osmflux,nunit) + Jp_cap_b(nu_osmflux,nunit)
    
    if(interstitial_volume_a(nu_intsat, nunit) + flux_a(nu_av_flux,  nunit) > interstitial_capacity_a)then
       excess = flux_a(nu_av_flux,  nunit) - (interstitial_capacity_a - interstitial_volume_a(nu_intsat, nunit))
       interstitial_volume_a(nu_intsat, nunit) = interstitial_capacity_a
       alveolar_volume (nu_alvflow,nunit)  = alveolar_volume(nu_alvflow,nunit) + 0.5_dp*excess
       
       if((interstitial_volume_b(nu_intsat, nunit) + 0.5_dp*excess) > interstitial_capacity_b)then
          overflow = 0.5_dp * excess - (interstitial_capacity_b - interstitial_volume_b(nu_intsat, nunit))
          alveolar_volume(nu_alvflow,nunit) = alveolar_volume(nu_alvflow,nunit) + overflow
          interstitial_volume_b(nu_intsat, nunit) = interstitial_capacity_b
       else
          interstitial_volume_b(nu_intsat, nunit) = interstitial_volume_b(nu_intsat, nunit) + 0.5_dp*excess
       endif
    else
       interstitial_volume_a(nu_intsat, nunit) = interstitial_volume_a(nu_intsat, nunit) + flux_a(nu_av_flux,nunit)
    endif
    
    unit_field (nu_alvflow,nunit)=alveolar_volume(nu_alvflow,nunit)
    interstitial_volume_b(nu_intsat, nunit) = interstitial_volume_b(nu_intsat, nunit) + flux_b(nu_av_flux,nunit)
    
    ! Handle protein overflow (after fluid overflow)
    if(interstitial_volume_a(nu_intsat,nunit) >= interstitial_capacity_a .and. excess > 0.0_dp) then
       ! Calculate overflow concentration
       if(interstitial_volume_a(nu_intsat,nunit) > volume_threshold) then
          c_overflow_a = Q_int_a(nu_osmflux,nunit) / interstitial_capacity_a
       else
          c_overflow_a = c_int_a(nu_osmflux,nunit)
       endif
       
       ! Protein distribution
       if(overflow > 0.0_dp) then
          Q_overflow_alv = c_overflow_a * (0.5_dp * excess + overflow)
          Q_overflow_b = c_overflow_a * (0.5_dp * excess - overflow) !the amount of protein left in intB when dlooding
       else
          Q_overflow_b = c_overflow_a * 0.5_dp * excess
          Q_overflow_alv = c_overflow_a * 0.5_dp * excess ! the flooding from inta to alveolar
       endif
       
       ! Update amounts
       Q_int_a(nu_osmflux,nunit) = Q_int_a(nu_osmflux,nunit) - Q_overflow_alv - Q_overflow_b
       Q_int_b(nu_osmflux,nunit) = Q_int_b(nu_osmflux,nunit) + Q_overflow_b
       
       if(Q_int_a(nu_osmflux,nunit) < 0.0_dp) then
          Q_int_a(nu_osmflux,nunit) = 0.0_dp
       end if
    endif
    
!!!! DIMENSIONALLY INCONSISTENT??????? ==> doesn't reduce to mm3; 200 is presumably R_alv which is highly assumptive based on parameterisation
    
    diffusion = (((interstitial_volume_a(nu_intsat, nunit)/interstitial_capacity_a)- &
         (interstitial_volume_b(nu_intsat, nunit)/interstitial_capacity_b))/ &
         (200_dp)) * (fluid_dt)
    interstitial_volume_b(nu_intsat, nunit) = interstitial_volume_b(nu_intsat, nunit) + diffusion
    interstitial_volume_a(nu_intsat, nunit) = interstitial_volume_a (nu_intsat, nunit)- diffusion
    
    ! Protein diffusion follows fluid
    !if(abs(diffusion) > volume_threshold) then
    if(diffusion > 0.0_dp .and. interstitial_volume_a(nu_intsat,nunit) > volume_threshold) then
       Jq_diffusion = diffusion * (Q_int_a(nu_osmflux,nunit) / interstitial_volume_a(nu_intsat,nunit))
       Q_int_a(nu_osmflux,nunit) = Q_int_a(nu_osmflux,nunit) - Jq_diffusion
       Q_int_b(nu_osmflux,nunit) = Q_int_b(nu_osmflux,nunit) + Jq_diffusion
    else if(diffusion < 0.0_dp .and. interstitial_volume_b(nu_intsat,nunit) > volume_threshold) then
       Jq_diffusion = abs(diffusion) * (Q_int_b(nu_osmflux,nunit) / interstitial_volume_b(nu_intsat,nunit))
       Q_int_b(nu_osmflux,nunit) = Q_int_b(nu_osmflux,nunit) - Jq_diffusion
       Q_int_a(nu_osmflux,nunit) = Q_int_a(nu_osmflux,nunit) + Jq_diffusion
    endif
    !endif
    
    if(interstitial_volume_b(nu_intsat, nunit)/interstitial_capacity_b < 0.3_dp)then
       lymph_conductivity = 1.48_dp * capillary_conductivity !all calculated does as a function of capillary_conductivity
       !no information on the size of pores or similar for lympatic conductivity so assumed to be similar to capillary.
    else
       lymph_conductivity = ((845.87_dp * (interstitial_volume_b(nu_intsat, nunit) / interstitial_capacity_b)**5.0_dp) + &
            (-2416.7_dp * (interstitial_volume_b(nu_intsat, nunit) / interstitial_capacity_b)**4.0_dp) + (2388.5_dp * &
            (interstitial_volume_b(nu_intsat, nunit) / interstitial_capacity_b)**3.0_dp) + (-922.24_dp * &
            (interstitial_volume_b(nu_intsat, nunit) / interstitial_capacity_b)**2.0_dp) + &
            (125.85_dp * (interstitial_volume_b(nu_intsat, nunit) / interstitial_capacity_b)) - 0.0067_dp)* &
            capillary_conductivity !(capillary_conductivity)
    endif
    !! Calculate osmotic reduction dynamically
    !!hydrostatic_gradient = capillary_pressure - interstitial_pressure_b(nu_Pe,nunit)
    !!osmotic_gradient = sigma * (osm_cap(nu_osmflux,nunit) - osm_int_b(nu_osmflux,nunit))
    !!effective_gradient = hydrostatic_gradient - osmotic_gradient
    !!
!!! Safety checks
    !!if(hydrostatic_gradient > 0.1_dp) then
    !!   osmotic_reduction_factor = max(0.01_dp, min(1.0_dp, effective_gradient / hydrostatic_gradient))
    !!else
    !!   osmotic_reduction_factor = 0.01_dp
    !!endif
    !!
    !!
!!! Apply Ashworth's logic WITH scaling
    !!if(interstitial_volume_b(nu_intsat,nunit)/interstitial_capacity_b < 0.3_dp) then
    !!   ! Below 30%: constant baseline (SCALED)
    !!   lymph_conductivity = 1.48_dp * capillary_conductivity * osmotic_reduction_factor
    !!
    !!
    !!else
    !!   ! Above 30%: polynomial (SCALED)
    !!   sat = interstitial_volume_b(nu_intsat,nunit) / interstitial_capacity_b
    !!
    !!   polynomial_factor = 845.87_dp * sat**5 + &
    !!                      -2416.7_dp * sat**4 + &
    !!                       2388.5_dp * sat**3 + &
    !!                       -922.24_dp * sat**2 + &
    !!                        125.85_dp * sat - 0.0067_dp
    !!
    !!   lymph_conductivity = polynomial_factor * capillary_conductivity * osmotic_reduction_factor
    !!endif
    
    P_initial_lymphtix (nu_Pe,nunit) = fluctuation/2.0_dp * sin_breath + pi/2.0_dp + &
         (lymph_diff-fluctuation)* ((interstitial_volume_b(nu_intsat, nunit) / interstitial_capacity_b)**2.0_dp) + &
         (lymph_min+(fluctuation/2.0_dp))
    !
    !          !arbitrarily defined mathematical relationship to show that lymphatic pressure does not change much at low volumes with a
    !          !large volume change, but at high volumes only a small volume change is needed to cause a large change in pressure
    !          !write(*,'(''Plym: '',f8.4)')initial_lymphatic_pressure  (diff_Pe /(2.0_dp * pi * 0.25_dp))
    !
    !
    if(interstitial_volume(nu_intsat, nunit).le.0.0_dp)then
       initial_lymph_flow(nu_lymphflow,nunit) = 0.0_dp
       interstitial_volume(nu_intsat, nunit) = 0.0_dp
    elseif (interstitial_pressure_b (nu_Pe,nunit)>  P_initial_lymphtix (nu_Pe,nunit))then
       initial_lymph_flow(nu_lymphflow,nunit) = (lymph_conductivity * unit_field(nu_sa,nunit) * &
            (interstitial_pressure_b(nu_Pe,nunit)-P_initial_lymphtix (nu_Pe,nunit))) * (fluid_dt)
    else
       initial_lymph_flow(nu_lymphflow,nunit) = 0.0_dp
    endif
    !          iv_array(2) = iv_array(2) - initial_lymph_flow
    !          initial_lymph_volume = initial_lymph_volume + initial_lymph_flow
    interstitial_volume_b(nu_intsat, nunit) = interstitial_volume_b(nu_intsat, nunit) &
         - initial_lymph_flow(nu_lymphflow,nunit)
    initial_lymph_volume(nu_lymphflow,nunit) = initial_lymph_volume(nu_lymphflow,nunit) &
         + initial_lymph_flow(nu_lymphflow,nunit)
    
    !
    !          int_osm_conc = int_osm_n(nu_osmflux,nunit)/interstitial_volume_b(nu_intsat, nunit)
    !!          liflowcount = liflowcount + initial_lymphatic_flow
    !          initial_osm_n(nu_osmflux,nunit) = initial_osm_n(nu_osmflux,nunit) + &
    !                  (initial_lymph_flow(nu_lymphflow,nunit)*int_osm_conc)
    !          int_osm_n(nu_osmflux,nunit)  =int_osm_n(nu_osmflux,nunit) -&
    !                  (initial_lymph_flow(nu_lymphflow,nunit)*int_osm_conc)
    !          if (initial_lymph_volume(nu_lymphflow,nunit) > 0.0_dp)then
    !             initial_lymph_conc = initial_osm_n(nu_osmflux,nunit) /initial_lymph_volume(nu_lymphflow,nunit)
    !          else
    !             initial_lymph_conc = 0.0_dp
    !          endif
    !
    !! Lymphatic protein removal
    !if(initial_lymph_flow(nu_lymphflow,nunit) > 0.0_dp .and. &
    !   interstitial_volume_b(nu_intsat,nunit) > volume_threshold) then
    !
    Jp_lymph_b(nu_osmflux,nunit) = initial_lymph_flow(nu_lymphflow,nunit) * c_int_b(nu_osmflux,nunit)
    Q_int_b(nu_osmflux,nunit) = Q_int_b(nu_osmflux,nunit) - Jp_lymph_b(nu_osmflux,nunit)
    Q_int_b(nu_osmflux,nunit) = max(0.0_dp, Q_int_b(nu_osmflux,nunit))
    !endif
    !
    !! Update plasma protein
    Q_plasma(nu_osmflux,nunit) = Q_plasma(nu_osmflux,nunit) - &
         (Jp_cap_a(nu_osmflux,nunit) + Jp_cap_b(nu_osmflux,nunit))
    Q_plasma(nu_osmflux,nunit) = max(0.0_dp, Q_plasma(nu_osmflux,nunit))
    !!          total_flux = total_hydro_flux(nu_flux,nunit) ! +total_osm_flux

!!! new
    interstitial_volume(nu_intsat, nunit)= interstitial_volume_a(nu_intsat, nunit)+ interstitial_volume_b(nu_intsat, nunit)
    interstitial_saturation(nu_intsat, nunit)= interstitial_volume(nu_intsat, nunit)/ interstitial_capacity
    
    unit_field(nu_intsat,nunit) = interstitial_saturation(nu_intsat, nunit)
    unit_field(nu_time,nunit) = unit_active_time(nunit) ! why time here is global time not the transit time
    unit_field(nu_av_flux,nunit) = total_hydro_flux(nu_flux,nunit)/unit_active_time(nunit)!total_flux/time  !flux_c
    unit_field(nu_lymphflow,nunit) =  initial_lymph_volume(nu_lymphflow,nunit)/unit_active_time(nunit)!initial_lymph_flow(nu_lymphflow,nunit)!initial_lymph_volume(nu_lymphflow,nunit)

  end subroutine alveolar_flux_unit_dt
  
!!!#############################################################################
    
  subroutine alveolar_flux(dt, time, T_interval)
    !*alveolar_capillary_flux:* calculate fluid flux from blood to interstitium

    use parameter_types, only: solve_V_params
    
    real(dp), intent(in) :: dt,T_interval
    real(dp) :: time
    ! Local variables
    integer :: i,nunit,n_nunit, count,fluid_steps, new_n_active
    integer, allocatable :: active_copy(:)
    real(dp) :: capillary_osm_n, cap_osm_conc,diffusion,excess,flux_c, &
         initial_lymph_conc,interstitial_osmotic, int_osm_conc,lymph_conductivity,&
         net_flux,overflow,sin_breath,sumuptake,test_time, total_flux,transit_time,capillary_SA 
    real(dp) :: capillary_pressure, P_elastic, diff_Pe,fluctuation, fluid_dt !,P_elastic_pre
    ! protein variables:
    real(dp) :: c_overflow_a, Jq_diffusion, Q_overflow_alv, Q_overflow_b, hydrostatic_gradient ,&
         effective_gradient, sat, osmotic_reduction_factor, polynomial_factor, dvdt, old_int_vol
    logical :: cont
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'alveolar_flux'
    call enter_exit(sub_name,1)

    if(.not.allocated(active))then
       allocate(active(num_units))
       active = [(i, i=1,num_units)]
       n_active = num_units
    endif
    allocate(active_copy(num_units))
    active_copy = 0

    fluid_steps = 2
    fluid_dt =dt/fluid_steps
    sin_breath = sin(2.0_dp*pi*0.25_dp*time)
    
    count=1
    !do while (count <= fluid_steps)
    do while (n_active > 0)
       !do nunit = 1, num_units
       new_n_active = 0
       time = time + solve_V_params%dt
       do n_nunit = 1, n_active
          nunit = active(n_nunit)

          call alveolar_flux_unit_dt(nunit, dt, time, T_interval)

          ! check for convergence
          sats(5,nunit) = sats(4,nunit)
          sats(4,nunit) = sats(3,nunit)
          sats(3,nunit) = sats(2,nunit)
          sats(2,nunit) = sats(1,nunit)
          sats(1,nunit) = interstitial_saturation(nu_intsat, nunit)
          
          lym_condition(nunit) = abs(((sats(1,nunit) + sats(2,nunit) + sats(3,nunit) + &
               sats(4,nunit) + sats(5,nunit))/5.0_dp) - sats(1,nunit))

          if ((lym_condition(nunit) > 0.00001_dp .or. &
               unit_active_time(nunit) < 200.0_dp * unit_field(nu_tt,nunit)) .and. &
               unit_active_time(nunit) < 5000.0_dp * unit_field(nu_tt,nunit)) then
             new_n_active = new_n_active + 1
             active_copy(new_n_active) = nunit
          endif
       enddo
       count = count + 1
       n_active = new_n_active
       active(1:n_active) = active_copy(1:n_active)
    end do

    deallocate(active_copy)
    
    call enter_exit(sub_name,2)
    
  end subroutine alveolar_flux

!!!#############################################################################

  subroutine lymphatic_transport(filename)
    !*lymphatic_transport:* whole system transport

!!! Inputs
    character(len=MAX_FILENAME_LEN), intent(in) :: filename
!!! Locals
    real(dp) :: capillary_flow,capillary_osm_n, &
         cap_osm_conc,diffusion,excess,flux_c, &
         initial_lymph_flow,initial_lymph_pressure,initial_lymph_volume, &
         initial_lymph_conc,interstitial_osmotic,interstitial_saturation,interstitial_volume, &
         int_osm_conc,lymph_conductivity, &
         net_flux,fluctuation,osm_flux,overflow,sumuptake,test_time,time_sum, &
         time_variable,total_flux,total_hydro_flux,capillary_pressure, int_osm_n,osm_n_flux! ,transit_time,capillary_SA
    character(len=300) :: writefile
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'lymphatic_transport'
    call enter_exit(sub_name,1)
    call set_lymph_factors(0,213.00_dp) ! mass, breathing rate, capillary volume raw (not used), number of timesteps
    
    !if(index(filename, ".oplymph")> 0) then !full filename is given
    !   writefile = filename
    !else ! need to append the correct filename extension
    !   writefile = trim(filename)//'.oplymph'
    !endif
    !
    !open(10, file=writefile, status='replace')
    !Only used for the osmotic model at the moment (which isn't operational) can volume be obtained elsewhere?
    alveolar_volume = 0.0_dp !alveolar volume likely greater at rest, but is lost to respiration - further information needed to put in model
    ! Is this where the tidal volume importing could go????
    
    unit_active_time = 0.0_dp
    
    ! initial lymphatic values
    initial_lymph_volume = 0.0_dp ! in mL ===> dependent on capillary_conductivity volume units
    !
    !    ! Osmotic pressures
    !    capillary_osm_n = 0.0_dp ! This doesnt change????? The bleed on effect means the rest of the osmotic flux doesnt work 1.025_dp / 66.5!
    !    interstitial_osmotic = 0.0_dp
    !    int_osm_n = 0.0_dp
    !    initial_osm_n = 0.0_dp
    !    total_osm_flux = 0.0_dp
    !
    !    osm_n_flux = 0.0_dp
    !    osm_flux = 0.0_dp
    total_hydro_flux = 0.0_dp
    
    !    ! Initialize protein arrays
    !
    !c_plasma = c_plasma_baseline
    !c_int_a = c_interstitial_baseline
    !c_int_b = c_interstitial_baseline
    !osm_cap = 26.7_dp
    !osm_int_a = 12.2_dp
    !osm_int_b = 12.2_dp
    !Jp_cap_a = 0.0_dp
    !Jp_cap_b = 0.0_dp
    !Jp_lymph_b = 0.0_dp
    !Jp_diffusive = 0.0_dp
    !Jp_convective = 0.0_dp
    ! INITIALIZE PROTEIN AMOUNTS FOR ALL UNITS
    !     Q_plasma = c_plasma_baseline * V_plasma_unit
    Q_int_a = c_interstitial_baseline * interstitial_volume_a
    Q_int_b = c_interstitial_baseline * interstitial_volume_b
    
    printcount = 0
    !
    !    ! These two would presumably change in a geometrically consistent lymphatics model
    
    sats(5,:) = 5.0_dp
    sats(4,:) = 4.0_dp
    sats(3,:) = 3.0_dp
    sats(2,:) = 2.0_dp
    sats(1,:) = 1.0_dp
    lym_condition = 0.5_dp
    
    !close(10)

    call enter_exit(sub_name,2)

  end subroutine lymphatic_transport

!!!#############################################################################
  
  subroutine set_lymph_factors(mass,cvr)
    
    integer,intent(in) :: mass
    real(dp),intent(in) :: cvr
    
    sex = mass
    !  breathing_rate = br
    
    ! dt or n_timesteps should be controlled by the user
    !  n_timesteps = n_time
    
    ! lung_mass = mass ! Replace the calculated value with this one when implementing the CT update
    ! Calculated values
    lung_mass = abs(real((1-sex)*840.0_dp))+real(sex)*639.0_dp  ! in g;female lung weight of 639g and male of 840g - should be updated from CT
    !  breathing_function = (2.0_dp*pi)/(60.0_dp/breathing_rate)
    
    ! interstitial_capacity == maximal volume before spillover into alveolar in mm^3 - based on 30ml.100g of fluid (Drake 2002)
    interstitial_capacity = ((30.0_dp*(lung_mass/100.0_dp))/real(num_units))*1000.0_dp !based on lung mass which should be obtained from CT
    !  IGC_T = IGC*T
    ! this may need to be adaptable for a dynamic model
    !  capillary_osmotic = capillary_molar_conc*IGC_T  ! Van't Hoffs osmotic pressure reduction - van't hoff factor, 'i' [real(1)] has been reduced to 1 to save storage
    
    ! this value is unused
    capillary_volume_raw = cvr  ! in mL Gehr 1978 based on having a body mass of 74 kg - should be unique to each person
    capillary_volume = (capillary_volume_raw*open_capillaries)/real(num_units) !in mL  Ben: unit_field(nu_vol,nunit)/1000.0_dp!volume in mm3 (from perfusion model) converted to mL !  unit_field(nu_vol,nunit) from venti not works
    !       interstitial values (capacity, volume) | first index corresponds to A, second to B
    !  ic_array = (/ 0.005_dp*interstitial_capacity, 0.995_dp*interstitial_capacity /) ! in mm3 ! arbitrarily sized
    !  iv_array = (/ 0.0_dp, 0.48_dp*interstitial_capacity /) ! in mm3 ! assumption of 48% saturation at rest
    interstitial_capacity_a = 0.005_dp*interstitial_capacity
    interstitial_capacity_b = 0.995_dp*interstitial_capacity
    interstitial_volume_a = 0.000005_dp*interstitial_capacity
    interstitial_volume_b = 0.48_dp*interstitial_capacity
    
    ! Set the scalar V_plasma_unit
    V_plasma_unit = capillary_volume * 1000.0_dp  ! mL to mm³

  end subroutine set_lymph_factors

end module lymphatics


!FUTURE DIRECTIONS
!input a constant to account for difference between current values and expected values
     !Model appeared to be working within the range of the literature but is likely off by a factor of 1000 due to nl to ul conversion error.
     !Need to check that outputted units are correct - most things are in ml and mmHg
!Lymphatic network tree
     !currently all lymph is returned to the circulation immediately, in reality it moves up a tree of lymphatics against a pressure gradient
     !would require excessive modelling perhaps
!impairment of gas diffusion caused by high interstitial saturation
     !unclear at what level this would occur
!alveolar flooding changes
     !alveolar flooding should be able to move between adjacent compartments
     !alveolar fluid should be removed via respiration naturally and therefore should always occur naturally at some low level
!individuality needs to be added in line with the other modules
     !currently operates only on preset male/female values

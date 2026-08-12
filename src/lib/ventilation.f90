module ventilation
  !*Brief Description:* This module handles all code specific to
  ! simulating ventilation
  !
  !*LICENSE:*
  !TBC
  !
  !
  !*Full Description:*
  !
  ! This module handles all code specific to simulating ventilation

  use arrays
  use diagnostics
  use geometry
  use indices
  use other_consts
  use precision
  use surfactant
  use lymphatics

  implicit none

  !Module parameters

  real(dp), dimension(:,:), allocatable :: unit_field_pre
  real(dp), dimension(:,:), allocatable :: unit_field_current
  real(dp), dimension(:,:), allocatable :: alv_unit_field_pre
!  real(dp), dimension(:,:), allocatable :: alv_unit_field_current

  real(dp), dimension(:,:), allocatable :: alv_radii_pre
  real(dp), dimension(:,:), allocatable :: alv_radii_current
  real(dp), dimension(:,:), allocatable :: aci_radii_current
  real(dp), dimension(:,:), allocatable :: alv_area_pre
  real(dp), dimension(:,:), allocatable :: alv_area_current
  real(dp), dimension(:,:), allocatable :: alv_dA
  real(dp), dimension(:,:), allocatable :: surf_concentration
  real(dp), dimension(:,:), allocatable :: surface_tension
  real(dp), dimension(:,:), allocatable :: Pc
  real(dp), dimension(:,:), allocatable :: Pc_com
  real(dp), dimension(:,:), allocatable :: smoothed_Pc_com
  real(dp), dimension(:,:), allocatable :: surf_concentration_pre
  real(dp), dimension(:,:), allocatable :: surface_tension_pre
  real(dp), dimension(:,:), allocatable :: Pe_unit_field_pre
  real(dp), dimension(:,:), allocatable :: Pe

  real(dp), parameter :: Pa_cmH2O = 98.0665_dp ! Pa/cmH2O


  !Module types

  !Module variables

  public surf_concentration
  public surface_tension
  public Pc
  public alv_area_current

  real(dp) :: ds_scale_op, flow_scale_op, vol_scale_op               ! scaling to appropriate units
  character(len=8) :: units_ds, units_flow, units_vol   ! controls the units in output
  
  !Interfaces
  private
  public evaluate_vent
  public evaluate_uniform_flow
  public sum_elem_field_from_periphery
  public update_alveolar_info

  public tissue_compliance

  real(dp),parameter,private :: gravity = 9.81e3_dp         ! mm/s2
!!! for air
  real(dp),parameter,private :: gas_density =   1.146e-6_dp ! g.mm^-3
  real(dp),parameter,private :: gas_viscosity = 1.8e-5_dp   ! Pa.s

contains

!!!#############################################################################

  subroutine evaluate_vent(filename)
    !*evaluate_vent:* Sets up and solves dynamic ventilation model

    use parameter_types, only: lung_params, mech_params, solve_V_params, V_params

!!! Inputs
    character(len=MAX_FILENAME_LEN),intent(in) :: filename
!!! Locals
    integer :: iter_step,n,ne,nunit
    real(dp) :: chestwall_restvol     ! resting volume of chest wall
    real(dp) :: p_mus                 ! muscle (driving) pressure
    real(dp) :: pmus_factor_ex        ! pmus_factor (_in and _ex) used to scale
    real(dp) :: pmus_factor_in        ! modifies driving pressures to converge
    !                                   tidal volume and expired volume to the
    !                                   target volume.
    real(dp) :: press_in_total        ! dynamic pressure at entry to model (Pa)
    real(dp) :: sum_expid             ! sum of expired volume  (mm^3)
    real(dp) :: sum_tidal             ! sum of inspired volume  (mm^3)
    real(dp) :: Texpn                 ! time for expiration (s)
    real(dp) :: Tinsp                 ! time for inspiration (s)
    real(dp) :: undef                 ! the zero stress volume. undef < RV 

    real(dp) :: dpmus,endtime,err_est,init_vol,last_vol, &
         current_vol,Pcw,ppl_current,pptrans,prev_flow,ptrans_frc, &
         sum_dpmus,sum_dpmus_ei,time,Tpass,ttime,volume_tree,WOBe,WOBr, &
         WOBe_insp,WOBr_insp,WOB_insp,WOB

    integer :: n_time_steps, time_step, n_time, i
    integer, parameter :: n_brth_surf_equil = 10
    character(len=300) :: writefile
    character(len=MAX_FILENAME_LEN) :: lymphfile
    logical :: continue_solve, converged, op_litres
    logical, parameter :: lymphatics = .false.
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'evaluate_vent'
    call enter_exit(sub_name,1)

    allocate( unit_field_pre(nu_vol,num_units)) !n_time_steps
    allocate( unit_field_current(nu_vol,num_units))
!    allocate( alv_unit_field(nu_vol,num_units))
    allocate( alv_unit_field_pre(nu_vol,num_units))
    allocate( alv_unit_field_current(nu_vol,num_units))
!    allocate( alv_radii(nu_vol,num_units))
    allocate( alv_radii_pre(nu_vol,num_units))
    allocate( alv_radii_current(nu_vol,num_units))
    allocate( aci_radii_current(nu_vol,num_units))
!    allocate( alv_area(nu_vol,num_units))
    allocate( alv_area_pre(nu_vol,num_units))
    allocate( alv_area_current(nu_vol,num_units))
    allocate( alv_dA(nu_vol,num_units))
    allocate( surf_concentration(nu_vol,num_units))
    allocate( surface_tension(nu_vol,num_units))
    allocate( Pc(nu_vol,num_units))
    allocate( Pc_com(nu_vol,num_units))
    allocate( smoothed_Pc_com(nu_vol,num_units))
    allocate( surf_concentration_pre(nu_vol,num_units))
    allocate( surface_tension_pre(nu_vol,num_units))
    allocate(Pe_unit_field_pre(nu_pe,num_units))

    allocate(interstitial_pressure_a(nu_pe,num_units))
    allocate(interstitial_pressure_b(nu_pe,num_units))
    allocate(P_initial_lymphtix(nu_Pe,num_units))
    allocate(total_hydro_flux(nu_flux,num_units))
    allocate(interstitial_volume (nu_intsat,num_units))
    allocate(interstitial_volume_a (nu_intsat,num_units))
    allocate(interstitial_volume_b (nu_intsat,num_units))
    allocate(initial_lymph_flow(nu_lymphflow,num_units))
    allocate( initial_lymph_volume(nu_lymphflow,num_units))
    allocate( interstitial_saturation (nu_intsat,num_units))
    allocate(flux_a(nu_av_flux,num_units))
    allocate(flux_b(nu_av_flux,num_units))
!    allocate(time_sum(num_units))
!    allocate(lym_time(num_units))
    allocate(lym_condition(num_units))
    allocate(sats(5,num_units))
    allocate(bulk_c(nu_vol,num_units))
!    allocate(alveolar_volume (alv_flux,num_units))
!    allocate(unit_field(nu_time,num_units))
!    allocate(unit_field(nu_av_flux,num_units))
!    allocate(unit_field(nu_lymphflow,num_units))
!    allocate( Pe(nu_pe,num_units))
!    allocate(ip_s)
!!    allocate(capillary_pressure)
!    allocate(cap_pressure(num_units))
!    allocate(transit_ttime(num_units))
!    allocate(cap_SA(num_units))
!    allocate( unit_field(nu_pe,num_units))
!    allocate(ppl_current (nu_pe,num_units))
!    allocate(osm_flux(nu_osmflux,num_units))
!    allocate(int_osm_n(nu_osmflux,num_units))
!    allocate(osm_n_flux(nu_osmflux,num_units))
    allocate(alveolar_volume(nu_alvflow,num_units))
!    allocate(capillary_volume(nu_capflow,num_units))
!    allocate(initial_osm_n(nu_osmflux,num_units))
!    allocate(total_osm_flux(nu_osmflux,num_units))
    allocate(unit_active_time(num_units))
!    allocate(initial_osm_n(nu_osmflux,num_units))
! Allocate protein arrays
!  allocate(V_plasma_unit(nu_time, num_units))
  allocate(Q_plasma(nu_osmflux, num_units))
  allocate(Q_int_a(nu_osmflux, num_units))
  allocate(Q_int_b(nu_osmflux, num_units))
  allocate(c_plasma(nu_osmflux, num_units))
  allocate(c_int_a(nu_osmflux, num_units))
  allocate(c_int_b(nu_osmflux, num_units))
  allocate(osm_cap(nu_osmflux, num_units))
  allocate(osm_int_a(nu_osmflux, num_units))
  allocate(osm_int_b(nu_osmflux, num_units))
  allocate(Jp_cap_a(nu_osmflux, num_units))
  allocate(Jp_cap_b(nu_osmflux, num_units))
  allocate(Jp_lymph_b(nu_osmflux, num_units))
  allocate(Jp_diffusive(nu_osmflux, num_units))
  allocate(Jp_convective(nu_osmflux, num_units))


    if(index(filename, ".exflow")> 0) then !full filename is given
       writefile = filename
    else ! need to append the correct filename extension
       writefile = trim(filename)//'.exflow'
    endif
    open(10, file=writefile, status='replace')
    
    write(10,'(A)') '   Time   Inflow    V_t     Raw     Comp    Ppl     Ptp     VolL    Pmus    Pcw  Pmus-Pcw'
    write(10,'(A)') '   (s)    (mL/s)   (mL) (cmH/L.s) (L/cmH) (...cmH2O...)    (L)     (......cmH2O.......)'

    ! control the formatting of output to suit size of model
    if(lung_params%FRC < 1e+5_dp)then
       op_litres = .false.
       ds_scale_op = 1.0_dp ! to scale ds and Vt mm^3/s to uL/s as output
       flow_scale_op = 1.0_dp ! to scale model mm^3/s to uL/s as output
       vol_scale_op = 1.0e+3_dp ! to scale model mm^3 to mL as output
       units_ds  = 'uL'
       units_flow = 'uL/s'
       units_vol = 'mL'
    else
       op_litres = .true.
       ds_scale_op = 1.0e+3_dp ! to scale ds and Vt mm^3 to mL as output
       flow_scale_op = 1.0e+3_dp ! to scale model mm^3/s to mL/s as output
       vol_scale_op = 1.0e+6_dp ! to scale model mm^3 to litres as output
       units_ds  = 'mL'
       units_flow = 'mL/s'
       units_vol = 'L'
    endif
    
!!! Initialise variables:
    pmus_factor_in = 1.0_dp
    pmus_factor_ex = 1.0_dp
    time = 0.0_dp !initialise the simulation time.
    n = 0 !initialise the 'breath number'. incremented at start of each breath.
    sum_tidal = 0.0_dp ! initialise the inspired and expired volumes
    sum_expid = 0.0_dp
    last_vol = 0.0_dp

    if(mech_params%surfactant_mechanics)then
       surf_concentration = 0.3e-6_dp/2.0_dp ! gamma*/2
       alv_dA = 0.0_dp
       surface_tension = 0.0_dp
       Pc_com = 100.0_dp
       Pc=0.0_dp
       unit_field(nu_Pe_max,:) = -1.0e6_dp
       unit_field(nu_Pe_min,:) = 1.0e6_dp
       lymphfile = 'output/P2BRP268-H12816_terminal'
       call lymphatic_transport(lymphfile)
    endif

!!! set dynamic pressure at entry. only changes for the 'pressure' option
    press_in_total = V_params%press_in
    
!!! calculate key variables from the boundary conditions/problem parameters
    Texpn = V_params%T_interval / (1.0_dp + V_params%i_to_e_ratio)
    Tinsp = V_params%T_interval - Texpn

!!! store initial branch lengths, radii, resistance etc. in array 'elem_field'
    call update_elem_field(1.0_dp)
    call update_resistance
    call volume_of_mesh(init_vol,volume_tree)

!!! distribute the initial tissue unit volumes along the gravitational axis.
    call set_initial_volume(lung_params%FRC)
    undef = mech_params%ref_vol_ratio * (lung_params%FRC - volume_tree)/dble(elem_units_below(1))

!!! calculate the total model volume
    call volume_of_mesh(init_vol,volume_tree)

    write(*,'(" Anatomical deadspace = ",f8.3, 1x, a)') volume_tree / ds_scale_op, trim(units_ds)
    write(*,'(" Respiratory volume   = ",f8.3, 1x, a)') (init_vol - volume_tree) / vol_scale_op, trim(units_vol)
    write(*,'(" Total lung volume    = ",f8.3, 1x, a)') init_vol / vol_scale_op, trim(units_vol)

    unit_field(nu_dpdt,1:num_units) = 0.0_dp

!!! calculate the compliance of each tissue unit
    call tissue_compliance(undef)
    call update_pleural_pressure(ppl_current) !calculate new pleural pressure
    pptrans=SUM(unit_field(nu_pe,1:num_units))/num_units

    chestwall_restvol = init_vol + lung_params%chest_wall_compliance * (-ppl_current)
    Pcw = (chestwall_restvol - init_vol)/lung_params%chest_wall_compliance
    write(*,'('' Chest wall RV       = '',f8.3, 1x, a)') chestwall_restvol / vol_scale_op, trim(units_vol)
        
    call write_flow_step_results(init_vol, current_vol,ppl_current,pptrans,Pcw,p_mus,0.0_dp,0.0_dp)
    
    continue_solve = .true.
    do while (continue_solve)
       n = n + 1 ! increment the breath number
       call initialise_breath(n, endtime, p_mus, ptrans_frc, sum_dpmus, sum_dpmus_ei, &
            sum_expid, sum_tidal, ttime)

!!! solve for a single breath (for time up to endtime)
       do while (time.lt.endtime) 
          ttime = ttime + solve_V_params%dt ! increment the breath time
          time = time + solve_V_params%dt ! increment the whole simulation time
!!!.......calculate the flow and pressure distribution for one time-step
          call evaluate_vent_step(chestwall_restvol,init_vol,last_vol,current_vol, &
               Pcw,pmus_factor_ex,pmus_factor_in,p_mus, pptrans,press_in_total,prev_flow, &
               ptrans_frc,sum_expid,sum_tidal,texpn,tinsp,ttime,undef,WOBe,WOBr, &
               WOBe_insp,WOBr_insp,WOB_insp,dpmus,iter_step)

          if(lymphatics) &
               call alveolar_flux_dt(solve_V_params%dt,time, V_params%T_interval) ! lymphatics subroutine

!!!.......update the estimate of pleural pressure
          call update_pleural_pressure(ppl_current) ! new pleural pressure
           
          call write_flow_step_results(init_vol,current_vol,ppl_current,pptrans,Pcw,p_mus,time,ttime)

       enddo !while time<endtime
       
       ! write out 'end of breath' information
       call write_end_of_breath(init_vol,current_vol,pmus_factor_in, &
            sum_expid,sum_tidal,WOBe_insp,WOBr_insp,WOB_insp)
       
       if(abs(V_params%tidal_volume).gt.1.0e-5_dp)THEN
          ! modify driving muscle pressure by volume_target/sum_tidal
          ! this increases p_mus for volume_target>sum_tidal, and
          ! decreases p_mus for volume_target<sum_tidal
          pmus_factor_in = pmus_factor_in * abs(V_params%tidal_volume/sum_tidal)
          pmus_factor_ex = pmus_factor_ex * abs(V_params%tidal_volume/sum_expid)
       endif

!!!....check whether simulation continues
       continue_solve = ventilation_continue(n,sum_tidal)

    enddo !...WHILE(CONTINUE_SOLVE)

    if(lymphatics)then
       write(*,'('' ====== SOLVING FOR ALVEOLAR FLUX ====== '')')
       ! stopping condition is internal to alveolar_flux
       call alveolar_flux(solve_V_params%dt,time, V_params%T_interval)
!!!.......update the estimate of pleural pressure
!          call update_pleural_pressure(ppl_current) ! new pleural pressure

       surf_concentration = 0.3e-6_dp/2.0_dp ! gamma*/2
       surface_tension = 0.0_dp
       Pc_com = 100.0_dp
       Pc=0.0_dp

!!! run for surfactant equilibration
       write(*,'('' ====== SOLVING FOR SURFACTANT EQUILIBRATION ====== '')')
       time = 0.0_dp
       ttime = 0.0_dp
       do i = 1, n_brth_surf_equil 
          n = n + 1 ! increment breath number
          call initialise_breath(n, endtime, p_mus, ptrans_frc, sum_dpmus, sum_dpmus_ei, &
               sum_expid, sum_tidal, ttime)
          
          do while (time.lt.endtime)
             ttime = ttime + solve_V_params%dt
             time = time + solve_V_params%dt
             
             call evaluate_vent_step(chestwall_restvol,init_vol,last_vol,current_vol, &
                  Pcw,pmus_factor_ex,pmus_factor_in,p_mus, &
                  pptrans,press_in_total,prev_flow,ptrans_frc, &
                  sum_expid,sum_tidal,texpn,tinsp,ttime,undef,WOBe,WOBr, &
                  WOBe_insp,WOBr_insp,WOB_insp, &
                  dpmus,iter_step)
!!!.......update the estimate of pleural pressure
             call update_pleural_pressure(ppl_current) ! new pleural pressure
             
             call write_flow_step_results(init_vol,current_vol,ppl_current,pptrans,Pcw,p_mus,time,ttime)
          enddo
          
          ! write out 'end of breath' information
          call write_end_of_breath(init_vol,current_vol,pmus_factor_in, &
               sum_expid,sum_tidal,WOBe_insp,WOBr_insp,WOB_insp)
       
          if(abs(V_params%tidal_volume).gt.1.0e-5_dp)THEN
             ! modify driving muscle pressure by volume_target/sum_tidal
             ! this increases p_mus for volume_target>sum_tidal, and
             ! decreases p_mus for volume_target<sum_tidal
             pmus_factor_in = pmus_factor_in * abs(V_params%tidal_volume/sum_tidal)
             pmus_factor_ex = pmus_factor_ex * abs(V_params%tidal_volume/sum_expid)
          endif

          continue_solve = ventilation_continue(n, sum_tidal)
          
       enddo
       
       n = n + 1 ! increment breath number
       
       call initialise_breath(n, endtime, p_mus, ptrans_frc, sum_dpmus, sum_dpmus_ei, &
            sum_expid, sum_tidal, ttime)
       
       do while (time.lt.endtime)
          ttime = ttime + solve_V_params%dt
          time = time + solve_V_params%dt
          
          call evaluate_vent_step(chestwall_restvol,init_vol,last_vol,current_vol, &
               Pcw,pmus_factor_ex,pmus_factor_in,p_mus, &
               pptrans,press_in_total,prev_flow,ptrans_frc, &
               sum_expid,sum_tidal,texpn,tinsp,ttime,undef,WOBe,WOBr, &
               WOBe_insp,WOBr_insp,WOB_insp, &
               dpmus,iter_step)
          
          call alveolar_flux(solve_V_params%dt,time, V_params%T_interval)
!!!.......update the estimate of pleural pressure
          call update_pleural_pressure(ppl_current) ! new pleural pressure
          
          call write_flow_step_results(init_vol,current_vol,ppl_current,pptrans,Pcw,p_mus,time,ttime)
          
       enddo
    endif ! surfactant_mechanics
    
    call write_end_of_breath(init_vol,current_vol,pmus_factor_in, &
         sum_expid,sum_tidal,WOBe_insp,WOBr_insp,WOB_insp)
     
!!! Transfer the tidal volume for each elastic unit to the terminal branches,
!!! and sum up the tree. Divide by inlet flow. This gives the time-averaged and
!!! normalised flow field for the tree.
    do nunit = 1,num_units
       ne = units(nunit) !local element number
       elem_field(ne_Vdot,ne) = unit_field(nu_vt,nunit)
    enddo
    unit_field(nu_vent,:) = unit_field(nu_vt,:)/(Tinsp+Texpn)
    call sum_elem_field_from_periphery(ne_Vdot)
    elem_field(ne_Vdot,1:num_elems) = &
         elem_field(ne_Vdot,1:num_elems)/elem_field(ne_Vdot,1)

    close(10)

    call enter_exit(sub_name,2)
    
  end subroutine evaluate_vent

!!!#############################################################################

  subroutine initialise_breath(n, endtime, p_mus, ptrans_frc, sum_dpmus, sum_dpmus_ei, &
       sum_expid, sum_tidal, ttime)

    use parameter_types, only: V_params, solve_V_params
!!! Inputs
    integer :: n
    real(dp) :: endtime, p_mus, ptrans_frc, sum_dpmus, sum_dpmus_ei, &
         sum_expid, sum_tidal, ttime

    ttime = 0.0_dp ! each breath starts with ttime=0
    endtime = V_params%T_interval * n - 0.5_dp * solve_V_params%dt ! the end time of this breath
    p_mus = 0.0_dp 
    ptrans_frc = SUM(unit_field(nu_pe,1:num_units))/num_units !ptrans at frc
    sum_tidal = 0.0_dp !reset the tidal volume
    sum_expid = 0.0_dp !reset the expired volume
    unit_field(nu_vt,1:num_units) = 0.0_dp !reset acinar tidal volume
    sum_dpmus = 0.0_dp
    sum_dpmus_ei = 0.0_dp
    
  end subroutine initialise_breath
  
!!!#############################################################################

  subroutine evaluate_vent_step(chestwall_restvol,init_vol,last_vol,current_vol,Pcw, &
       pmus_factor_ex,pmus_factor_in,p_mus,pptrans, &
       press_in_total,prev_flow,ptrans_frc,sum_expid, &
       sum_tidal,texpn,tinsp,ttime,undef,WOBe,WOBr,WOBe_insp,WOBr_insp, &
       WOB_insp,dpmus,iter_step)

    use parameter_types, only: lung_params, mech_params, solve_V_params

    real(dp),intent(in) :: chestwall_restvol, &
         init_vol,pmus_factor_ex,pmus_factor_in,pptrans, &
         press_in_total,ptrans_frc,texpn,tinsp,ttime,undef
    real(dp) :: last_vol,current_vol,Pcw,prev_flow,p_mus, &
         sum_expid,sum_tidal,WOBe,WOB_insp,WOBe_insp, &
         WOBr,WOBr_insp
    ! Local variables
    integer :: iter_step
    real(dp) :: dpmus,err_est,volume_tree, fluctuation
    logical :: converged

    ! --------------------------------------------------------------------------
    
    

!!! Solve for a new flow and pressure field
!!! We will estimate the flow into each terminal lumped
!!! parameter unit (assumed to be an acinus), so we can calculate flow
!!! throughout the rest of the tree simply by summation. After summing
!!! the flows we can use the resistance equation (P0-P1=R1*Q1) to update
!!! the pressures throughout the tree.

    ! set the increment in driving (muscle) pressure
    call set_driving_pressures(dpmus,pmus_factor_ex,pmus_factor_in, &
         p_mus,Texpn,Tinsp,ttime)
    prev_flow = elem_field(ne_Vdot,1)

    !initialise Qinit to the previous flow
    elem_field(ne_Vdot0,1:num_elems) = elem_field(ne_Vdot,1:num_elems)
    converged = .FALSE.
    iter_step=0
    do while (.not.converged)
       iter_step = iter_step+1 !count the iterative steps
       call estimate_flow(dpmus,err_est) !analytic solution for Q
       if(iter_step.gt.1.and.err_est < solve_V_params%err_tolerance)then
          converged = .TRUE.
       else if(iter_step > solve_V_params%max_iterations)then
          converged = .TRUE.
          write(*,'('' Warning: lower convergence '// &
               'tolerance and time step - check values, Error='',D10.3)') &
               err_est
       endif
       call sum_elem_field_from_periphery(ne_Vdot) !sum flows UP tree
       call update_elem_field(1.0_dp)
       call update_resistance ! updates resistances
       call update_node_pressures(press_in_total) ! updates the pressures at nodes
       call update_unit_dpdt() ! update dP/dt at the terminal units
    enddo !converged
    
    call update_unit_volume() ! Update tissue unit volumes, unit tidal vols

    if(mech_params%surfactant_mechanics)then
       call update_alveolar_info(solve_V_params%dt, alv_dA)
       call update_surfactant_concentration(solve_V_params%dt, alv_area_current, alv_dA, surf_concentration)
       call update_surface_tension(surf_concentration, surface_tension)
       call update_Pc_compliance(unit_field,alv_radii_current,surface_tension, Pc, Pc_com, smoothed_Pc_com)
    endif
    
    call volume_of_mesh(current_vol,volume_tree) ! calculate mesh volume
    call update_elem_field(1.0_dp)
    call update_resistance  !update element lengths, volumes, resistances
    call tissue_compliance(undef) ! unit compliances
    call update_proximal_pressure ! pressure at proximal nodes of end branches
    !call calculate_work(current_vol-init_vol,current_vol-last_vol,WOBe,WOBr, &
    !     pptrans,WOB)!calculate work of breathing
    last_vol=current_vol
    Pcw = (chestwall_restvol - current_vol)/lung_params%chest_wall_compliance
    
    ! increment the tidal volume, or the volume expired
    if(elem_field(ne_Vdot,1).gt.0.0_dp)then
       sum_tidal = sum_tidal+elem_field(ne_Vdot,1) * solve_V_params%dt
    else
       sum_expid = sum_expid-elem_field(ne_Vdot,1) * solve_V_params%dt
       if(prev_flow.gt.0.0_dp)then
          WOBe_insp = (WOBe+sum_tidal*ptrans_frc*1.0e-9_dp)*(30.0_dp/Tinsp)
          WOBr_insp = WOBr*(30.0_dp/Tinsp)
          WOB_insp = WOBe_insp+WOBr_insp
          WOBe = 0.0_dp
          WOBr = 0.0_dp
       endif
    endif
!print*,'wobe',ptrans_frc
  end subroutine evaluate_vent_step

!!!#############################################################################

  subroutine evaluate_uniform_flow
    !*evaluate_uniform_flow:* Sets up and solves uniform ventilation model

    ! Local variables
    integer :: ne,nunit
    real(dp) :: init_vol,volume_tree

    ! --------------------------------------------------------------------------

!!! calculate the total model volume
    call volume_of_mesh(init_vol,volume_tree)

!!! initialise the flow field to zero
    elem_field(ne_Vdot,1:num_elems) = 0.0_dp

!!! For each elastic unit, calculate uniform ventilation
    do nunit = 1,num_units
       ne = units(nunit) !local element number
       unit_field(nu_Vdot0,nunit) = unit_field(nu_vol,nunit)/ &
            (init_vol-volume_tree)
       elem_field(ne_Vdot,ne) = unit_field(nu_Vdot0,nunit)
    enddo

    call sum_elem_field_from_periphery(ne_Vdot)

  end subroutine evaluate_uniform_flow


!!!#############################################################################

  subroutine set_driving_pressures(dpmus,pmus_factor_ex,pmus_factor_in, &
       p_mus,Texpn,Tinsp,ttime)

    use parameter_types, only: solve_V_params, V_params

    real(dp),intent(in) :: pmus_factor_ex,pmus_factor_in,Texpn, &
         Tinsp,ttime
    real(dp) :: dpmus,p_mus
    ! Local variables
    real(dp) :: sum_dpmus,sum_dpmus_ei,Tpass
    
    ! --------------------------------------------------------------------------

    select case(V_params%expiration_type)
       
    case("active")
       if(ttime.lt.Tinsp)then
          dpmus = V_params%insp_press_muscle * pmus_factor_in * PI * &
               sin(pi/Tinsp*ttime)/(2.0_dp*Tinsp)*solve_V_params%dt
       elseif(ttime.le.Tinsp+Texpn)then
          dpmus = V_params%insp_press_muscle * pmus_factor_ex * PI * &
               sin(2.0_dp*pi*(0.5_dp+(ttime-Tinsp)/(2.0_dp*Texpn)))/ &
               (2.0_dp*Texpn)*solve_V_params%dt
       endif

    case("passive")
       if(ttime.le.Tinsp+0.5_dp*solve_V_params%dt)then
          dpmus = V_params%insp_press_muscle * pmus_factor_in * PI * solve_V_params%dt * &
               sin(pi*ttime/Tinsp)/(2.0_dp*Tinsp)
          sum_dpmus = sum_dpmus+dpmus
          sum_dpmus_ei = sum_dpmus
       else
          Tpass = 0.1_dp
          dpmus = MIN(-sum_dpmus_ei/(Tpass*Texpn)*solve_V_params%dt,-sum_dpmus)
          sum_dpmus = sum_dpmus+dpmus
       endif

    end select

    p_mus = p_mus + dpmus !current value for muscle pressure

  end subroutine set_driving_pressures

!!!#############################################################################

  subroutine update_unit_dpdt()
    !*update_unit_dpdt:* updates the rate of change of pressure at the proximal
    ! end of element that supplies tissue unit. i.e. not the rate of change of
    ! pressure within the unit.

    use parameter_types, only: solve_V_params
    
    ! Local variables
    integer :: ne,np1,nunit
    real(dp) :: est

    ! --------------------------------------------------------------------------

    do nunit = 1,num_units
       ne = units(nunit)
       np1 = elem_nodes(1,ne)
       ! linear estimate
       est = (node_field(nj_aw_press,np1) &
            - unit_field(nu_air_press,nunit))/solve_V_params%dt
!!!    For stability, weight new estimate with the previous dP/dt
       unit_field(nu_dpdt,nunit) = 0.5_dp*(est+unit_field(nu_dpdt,nunit))
    enddo !nunit

  end subroutine update_unit_dpdt


!!!#############################################################################

  subroutine update_proximal_pressure
    !*update_proximal_pressure:* Update the pressure at the proximal node of
    ! the element that feeds an elastic unit

    ! Local variables
    integer :: ne,np1,nunit

    ! --------------------------------------------------------------------------

    do nunit = 1,num_units
       ne = units(nunit)
       np1 = elem_nodes(1,ne)
!!!    store the entry node pressure as an elastic unit air pressure
       unit_field(nu_air_press,nunit) = node_field(nj_aw_press,np1)
    enddo !noelem

  end subroutine update_proximal_pressure


!!!#############################################################################

  subroutine update_pleural_pressure(ppl_current)
    !*update_pleural_pressure:* Update the mean pleural pressure based on
    ! current Pel (=Ptp) and Palv, i.e. Ppl(unit) = -Pel(unit)+Palv(unit)

    real(dp),intent(out) :: ppl_current
!    real(dp):: Pi
    ! Local variables
    integer :: ne,np2,nunit

    ! --------------------------------------------------------------------------

    ppl_current = 0.0_dp
    do nunit = 1,num_units
       ne = units(nunit)
       np2 = elem_nodes(2,ne)

       !calculate Ppl(unit) = -Pel(unit)+Palv(unit)
       ppl_current  = ppl_current   - unit_field(nu_pe,nunit) + &
            node_field(nj_aw_press,np2)

       !calculate the interstitium pressure
!       Pi = node_field(nj_aw_press,np2) - Pc(nu_vol,nunit)
       !calculate Ppl(unit) = -Pel(unit)-Pc(unit)+Palv(unit)
!       ppl_current = ppl_current - unit_field(nu_pe,nunit) - Pc(nu_vol,nunit) + &
!            node_field(nj_aw_press,np2)
!       ppl_current = ppl_current - Pc(nu_vol,nunit) + &
!            node_field(nj_aw_press,np2)
    enddo !noelem
!    print*,'ppl pressure',ppl_current
    ppl_current = ppl_current/num_units
!    print*,'ppl pressure',ppl_current
!    print*,'alveolar pressure',node_field(nj_aw_press,1)

  end subroutine update_pleural_pressure


!!!#############################################################################

  subroutine update_node_pressures(press_in)
    !*update_node_pressures:* Use the known resistances and flows to calculate
    ! nodal pressures through whole tree

    real(dp),intent(in) :: press_in
    !Local parameters
    integer :: ne,np1,np2

    ! --------------------------------------------------------------------------

    ! set the initial node pressure to be the input pressure (usually zero)
    ne = 1 !element number at top of tree, usually = 1
    np1 = elem_nodes(1,ne) !first node in element
    node_field(nj_aw_press,np1) = press_in !set pressure at top of tree

    do ne = 1,num_elems !for each element
       np1 = elem_nodes(1,ne) !start node number
       np2 = elem_nodes(2,ne) !end node number
       !P(np2) = P(np1) - Resistance(ne)*Flow(ne)
       node_field(nj_aw_press,np2) = node_field(nj_aw_press,np1) &
            - (elem_field(ne_resist,ne)*  elem_field(ne_Vdot,ne))* &
           dble(elem_ordrs(no_type,ne))

       !noresistance
!       node_field(nj_aw_press,np2) = node_field(nj_aw_press,np1)
       !saline condition simulate
!            - ((0.8_dp)* elem_field(ne_resist,ne)*elem_field(ne_Vdot,ne))* &
!           dble(elem_ordrs(no_type,ne))
    enddo !noelem

  end subroutine update_node_pressures


!!!#############################################################################

  subroutine tissue_compliance(undef)

    use parameter_types, only: mech_params

    real(dp), intent(in) :: undef
    ! Local variables
    integer :: nalv, ne,nunit
    real(dp) :: ab_term, exp_term,lambda,ratio, time

    ! --------------------------------------------------------------------------

    !.....dV/dP=1/[(1/2h^2).c/2.(3a+b)exp().(4h(h^2-1)^2)+(h^2+1)/h^2)]

    ab_term = 3.0_dp * mech_params%a + mech_params%b
    
    do nunit = 1,num_units
       ne = units(nunit)
       !calculate a compliance for the tissue unit
       ratio = unit_field(nu_vol,nunit)/undef
       lambda = ratio**(1.0_dp/3.0_dp) !uniform extension ratio
       exp_term = exp(0.75_dp * ab_term * (lambda**2 - 1.0_dp)**2)

       unit_field(nu_comp,nunit) = mech_params%cc * exp_term/6.0_dp * (3.0_dp * ab_term**2 &
            * (lambda**2 - 1.0_dp)**2 / lambda**2 + ab_term &
            * (lambda**2 + 1.0_dp) / lambda**4)
       unit_field(nu_comp,nunit) = undef/unit_field(nu_comp,nunit) ! V/P

       if(mech_params%surfactant_mechanics)then
          !add the Pc compliance in parallel
          unit_field(nu_comp,nunit) =1.0_dp/(1.0_dp/unit_field(nu_comp,nunit)&
               +1.0_dp/Pc_com(nu_vol,nunit))
          Pe_unit_field_pre(nu_pe,nunit) =unit_field(nu_pe,nunit)
       endif
       
       !estimate an elastic recoil pressure for the unit
       unit_field(nu_pe,nunit) = mech_params%cc / 2.0_dp * ab_term * (lambda**2.0_dp &
            -1.0_dp) * exp_term/lambda

       if(mech_params%surfactant_mechanics)then
          ! set minimum and maximum Pe values
          unit_field(nu_Pe_max,nunit) = max(unit_field(nu_Pe_max,nunit), &
               unit_field(nu_pe,nunit))
          unit_field(nu_Pe_min,nunit) = min(unit_field(nu_Pe_min,nunit), &
               unit_field(nu_pe,nunit))
       endif
          
    enddo !nunit

  end subroutine tissue_compliance


!!!#############################################################################

  subroutine sum_elem_field_from_periphery(ne_field)

    integer,intent(in) :: ne_field
    !Local parameters
    real(dp) :: field_value
    integer :: i,ne,ne2

    ! --------------------------------------------------------------------------

    do ne = num_elems,1,-1
       if(elem_cnct(1,0,ne).gt.0)then !not terminal
          field_value = 0.0_dp
          do i = 1,elem_cnct(1,0,ne) !for each possible daughter branch (max 2)
             ne2 = elem_cnct(1,i,ne) !the daughter element number
             field_value = field_value+dble(elem_symmetry(ne2))* &
                  elem_field(ne_field,ne2) !sum daughter fields
          enddo !noelem2
          elem_field(ne_field,ne) = field_value
       endif
    enddo !noelem

  end subroutine sum_elem_field_from_periphery

!!!#############################################################################

  subroutine update_unit_volume()

    use parameter_types, only: mech_params, solve_V_params
    
    ! Local variables
    integer :: ne,np,nunit

    ! --------------------------------------------------------------------------

    do nunit = 1,num_units
       ne = units(nunit)
       np = elem_nodes(2,ne)

       if(mech_params%surfactant_mechanics)then
          ! collect the volume of the lumped tissue unit for previous time step
          unit_field_pre(nu_vol,nunit)=unit_field(nu_vol,nunit)
       endif
       
       ! update the volume of the lumped tissue unit
       unit_field(nu_vol,nunit) = unit_field(nu_vol,nunit)+solve_V_params%dt* &
            elem_field(ne_Vdot,ne)!in mm^3

       if(mech_params%surfactant_mechanics)then
          ! collect the volume of the lumped tissue unit for current time step
          unit_field_current(nu_vol,nunit)=unit_field(nu_vol,nunit)
       endif

       if(elem_field(ne_Vdot,1).gt.0.0_dp)then  !only store inspired volume
          unit_field(nu_vt,nunit) = unit_field(nu_vt,nunit)+solve_V_params%dt* &
               elem_field(ne_Vdot,ne)
       endif

    enddo !nunit

  end subroutine update_unit_volume

!!!#############################################################################

  subroutine update_alveolar_info(dt, alv_dA)


!    real(dp), dimension(:,:), intent(out) :: alv_area_current
    real(dp), dimension(:,:), intent(out) :: alv_dA

    real(dp),intent(in) :: dt
    real(dp) :: time

    integer :: nalv,nunit

    character(len=60) :: sub_name

  !--------------------------------------------------------------------------

    sub_name = 'update_alveolar_info'
    call enter_exit(sub_name,1)




    do nunit = 1,num_units
        !pre_step alveolus volume cm3
        !assume alveolus shape is hemisphere
        alv_unit_field_pre(nu_vol,nunit)=unit_field_pre(nu_vol,nunit)/ 37700000 !mm3 to cm3
        !r=sqrt[3]{{3V}{2*pi}}
        alv_radii_pre(nu_vol,nunit) = ((3.0_dp * alv_unit_field_pre(nu_vol,nunit)) / (2.0_dp * PI)) &
                **(1.0_dp/3.0_dp)

        alv_area_pre(nu_vol,nunit)=2.0_dp * PI * (alv_radii_pre(nu_vol,nunit)**2.0_dp)

        !current alveolus volume cm3
        alv_unit_field_current(nu_vol,nunit)=unit_field_current(nu_vol,nunit)/ 37700000 !mm3 to cm3

        alv_radii_current(nu_vol,nunit) = ((3.0_dp * alv_unit_field_current(nu_vol,nunit)) / (2.0_dp * PI)) &
                **(1.0_dp/3.0_dp)

        alv_area_current(nu_vol,nunit)=2.0_dp * PI * (alv_radii_current(nu_vol,nunit)**2.0_dp)

        !Area difference
        alv_dA(nu_vol,nunit)=(alv_area_current(nu_vol,nunit)-alv_area_pre(nu_vol,nunit))/dt

    enddo

    call enter_exit(sub_name,2)

  end subroutine update_alveolar_info

!!!#############################################################################

  subroutine update_elem_field(alpha)

    real(dp),intent(in) :: alpha   ! the factor by which the radius changes
    ! Local variables
    integer :: ne,np1,np2
    real(dp) :: gamma,resistance,reynolds,zeta
    real(dp) :: rad,le

    ! --------------------------------------------------------------------------

    do ne = 1,num_elems
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)

       ! element length
       elem_field(ne_length,ne) = sqrt((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)

       ! element radius
       elem_field(ne_radius,ne) = sqrt(alpha) * elem_field(ne_radius,ne)

       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
    enddo ! ne

  end subroutine update_elem_field

!!!#############################################################################

  subroutine update_resistance

    ! Local variables
    integer :: i,ne,ne2,np1,np2,nunit
    real(dp) :: ett_resistance,gamma,le,rad,resistance,reynolds,sum,zeta
    real(dp) :: tissue_resistance

    ! --------------------------------------------------------------------------

    elem_field(ne_t_resist,1:num_elems) = 0.0_dp

    tissue_resistance = 0.0_dp  ! 0.35_dp * 98.0665_dp/1.0e6_dp

    do nunit = 1,num_units
       ne = units(nunit)
       elem_field(ne_t_resist,ne) = tissue_resistance * dble(elem_units_below(1))
    enddo

    do ne = 1,num_elems
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)

       le = elem_field(ne_length,ne)
       rad = elem_field(ne_radius,ne)

       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3
       resistance = 8.0_dp*GAS_VISCOSITY*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius,ne)**4) !laminar resistance

       ! element turbulent resistance (flow in bifurcating tubes)
       gamma = 0.357_dp !inspiration
       if(elem_field(ne_Vdot,ne).lt.0.0_dp) gamma = 0.46_dp !expiration

       reynolds = abs(elem_field(ne_Vdot,ne)*2.0_dp*GAS_DENSITY/ &
            (pi*elem_field(ne_radius,ne)*GAS_VISCOSITY))
       zeta = MAX(1.0_dp,dsqrt(2.0_dp*elem_field(ne_radius,ne)* &
            reynolds/elem_field(ne_length,ne))*gamma)
       elem_field(ne_resist,ne) = resistance * zeta
       elem_field(ne_t_resist,ne) = elem_field(ne_resist,ne) + &
            elem_field(ne_t_resist,ne)
    enddo !noelem

    do ne = num_elems,1,-1
       sum = 0.0_dp
       if(elem_cnct(1,0,ne).gt.0)then !not terminal
          do i = 1,elem_cnct(1,0,ne) !for each possible daughter branch (max 2)
             ne2 = elem_cnct(1,i,ne) !the daughter element number
             ! line below is sum = sum + 1/R, where 1/R is multiplied by
             !  2 if this is a symmetric child branch
             sum = sum + dble(elem_symmetry(ne2))* &
                  dble(elem_ordrs(no_type,ne2))/elem_field(ne_t_resist,ne2)
          enddo
          if(sum.gt.zero_tol) elem_field(ne_t_resist,ne) = &
               elem_field(ne_t_resist,ne) + 1.0_dp/sum
       endif
    enddo

  end subroutine update_resistance

!!!#############################################################################

  subroutine estimate_flow(dp_external,err_est)

    use parameter_types, only: solve_V_params
    
    real(dp),intent(in) :: dp_external
    real(dp),intent(out) :: err_est
    ! Local variables
    integer :: ne,nunit
    real(dp) :: alpha,beta,flow_diff,flow_sum,Q,Qinit

    ! --------------------------------------------------------------------------

    err_est = 0.0_dp
    flow_sum = 0.0_dp

!!! For each elastic unit, calculate Qbar (equation 4.13 from Swan thesis)
    do nunit = 1,num_units !for each terminal only (with tissue units attached)
       ne = units(nunit) !local element number
       ! Calculate the mean flow into the unit in the time step
       ! alpha is rate of change of pressure at start node of terminal element
       alpha = unit_field(nu_dpdt,nunit) !dPaw/dt, updated each iter
       Qinit = elem_field(ne_Vdot0,ne) !terminal element flow, updated each dt
       ! beta is rate of change of 'external' pressure, incl muscle and entrance
       beta = dp_external/solve_V_params%dt ! == dPmus/dt (-ve for insp), updated each dt

!!!    Q = C*(alpha-beta)+(Qinit-C*(alpha-beta))*exp(-dt/(C*R))
       Q = unit_field(nu_comp,nunit)*(alpha-beta)+ &
            (Qinit-unit_field(nu_comp,nunit)*(alpha-beta))* &
            exp(-solve_V_params%dt/(unit_field(nu_comp,nunit)*elem_field(ne_t_resist,ne)))

       unit_field(nu_Vdot2,nunit) = unit_field(nu_Vdot1,nunit) !flow at iter-2
       unit_field(nu_Vdot1,nunit) = unit_field(nu_Vdot0,nunit) !flow at iter-1

!!!    for stability the flow estimate for current iteration
!!!    includes flow estimates from previous two iterations
       unit_field(nu_Vdot0,nunit) = 0.75_dp*unit_field(nu_Vdot2,nunit)+ &
            0.25_dp*(Q+unit_field(nu_Vdot1,nunit))*0.5_dp

       flow_diff = unit_field(nu_Vdot0,nunit) - elem_field(ne_Vdot,ne)
       if(abs(flow_diff).gt.zero_tol) &
            err_est = err_est+flow_diff**2 !sum up the error for all elements
       if(abs(unit_field(nu_Vdot0,nunit)).gt.zero_tol) &
            flow_sum = flow_sum+unit_field(nu_Vdot0,nunit)**2


!!! ARC: DO NOT CHANGE BELOW. THIS IS NEEDED FOR THE ITERATIVE STEP
!!! - SIMPLER OPTIONS JUST FORCE IT TO CONVERGE WHEN ITS NOT
       elem_field(ne_Vdot,ne) = (unit_field(nu_Vdot0,nunit)&
            +unit_field(nu_Vdot1,nunit))/2.0_dp
       unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
    enddo !nunit
!    print*,'Q',Q
    ! the estimate of error for the iterative solution
    if(abs(flow_sum*dble(num_units)).gt.zero_tol) then
       err_est = err_est/(flow_sum*dble(num_units))
    else
       err_est = err_est/dble(num_units)
    endif

  end subroutine estimate_flow

!!!#############################################################################

  subroutine calculate_work(breath_vol,dt_vol,WOBe,WOBr,pptrans,WOB)

    real(dp) :: breath_vol,dt_vol,WOBe,WOBr,pptrans,WOB
    ! Local variables
    integer :: ne,np1,nunit
    real(dp) :: p_resis,p_trans

    ! --------------------------------------------------------------------------

    p_resis = 0.0_dp
    !estimate elastic and resistive WOB for each dt (sum dP.V)

    p_trans = SUM(unit_field(nu_pe,1:num_units))/num_units

    do nunit = 1,num_units
       ne = units(nunit)
       np1 = elem_nodes(2,ne)
       p_resis = p_resis+node_field(nj_aw_press,1)-node_field(nj_aw_press,np1)
    enddo
    p_resis=p_resis/num_units
    ! vol in mm3 *1e-9=m3, pressure in Pa, hence *1d-9 = P.m3 (Joules)
!    Print*,'PP', p_trans
!    Print*,'PP2', pptrans
!     Print*,'v',breath_vol
    WOBe = WOBe+(p_trans-pptrans)*breath_vol*1.0e-9_dp
    WOBr = WOBr+p_resis*dt_vol*1.0e-9_dp
    WOB = WOBe+WOBr

    pptrans = p_trans

  end subroutine calculate_work

!!!#############################################################################

!!!#############################################################################

  subroutine write_end_of_breath(init_vol,current_vol,pmus_factor_in, &
       sum_expid,sum_tidal,WOBe_insp,WOBr_insp,WOB_insp)

    use parameter_types, only: V_params
    
    real(dp),intent(in) :: init_vol,current_vol,pmus_factor_in, &
         sum_expid,sum_tidal,WOBe_insp,WOBr_insp,WOB_insp

    ! --------------------------------------------------------------------------

    write(*,'('' End of breath, inspired = '',f10.2, 1x, a)') sum_tidal / vol_scale_op, trim(units_vol)
    write(*,'('' End of breath, expired  = '',f10.2, 1x, a)') sum_expid / vol_scale_op, trim(units_vol)
    write(*,'('' Peak muscle pressure    = '',f10.2,'' cmH2O'')') &
         V_params%insp_press_muscle * pmus_factor_in/Pa_cmH2O
    write(*,'('' Drift in FRC from start = '',f10.2,'' %'')') &
         100*(current_vol-init_vol)/init_vol
    write(*,'('' Difference from target Vt = '',f8.2,'' %'')') &
         100*(V_params%tidal_volume-sum_tidal) / V_params%tidal_volume
    write(*,'('' Total Work of Breathing ='',f7.3,''J/min'')')WOB_insp
    write(*,'('' elastic WOB ='',f7.3,''J/min'')')WOBe_insp
    write(*,'('' resistive WOB='',f7.3,''J/min'')')WOBr_insp

  end subroutine write_end_of_breath

!!!#############################################################################

  subroutine write_flow_step_results_surf(chest_wall_compliance,init_vol, &
       current_vol,ppl_current,pptrans,Pcw,p_mus,time,ttime, WOBe,WOBr,WOB)

    real(dp),intent(in) :: chest_wall_compliance,init_vol,current_vol, &
         ppl_current,pptrans,Pcw,p_mus,time,ttime, WOBe,WOBr,WOB
    ! Local variables
    real(dp) :: totalC,Precoil,Paw,totalV

    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'write_flow_step_results'
    call enter_exit(sub_name,1)

    !allocate( alv_unit_field(nu_vol,num_units))
    !allocate( alv_radii(nu_vol,num_units))
    !allocate( alv_area(nu_vol,num_units))
    !allocate( alv_dA(nu_vol,num_units))
    !allocate( unit_field_pre(nu_vol,num_units))
    !allocate(unit_field_current(nu_vol,num_units))

    !the total model compliance
    totalC = 1.0_dp/(1.0_dp/sum(unit_field(nu_comp,1:num_units))+ &
         1.0_dp/chest_wall_compliance)
    Precoil = sum(unit_field(nu_vol,1:num_units))/num_units
!    Precoil = sum(unit_field(nu_pe,1:num_units)*unit_field(nu_vol,1:num_units))/sum(unit_field(nu_vol,1:num_units))
    totalV = sum(unit_field(nu_vol,1:num_units))/num_units
    Paw = sum(unit_field(nu_air_press,1:num_units))/num_units!sum(node_field(nj_aw_press,1:num_elems))/num_elems
!    Paw = Precoil + ppl_current
    if(abs(time).lt.zero_tol)then
!!! write out the header information for run-time output
       write(*,'(2X,''Time'',3X,''Inflow'',4X,''V_t'',6X,''Raw'',6X,&
            &''Comp'',5X,''Ppl'',6X,''Ptp'',6X,''VolL'',5X,''Pmus'',&
            &5X,''Pcw'',4X,''Pmus-Pcw'',5X,''Aci_vol'',6X,''Alv_vol'',&
            &7X,''Sur_con'',7X,''alv_radi'',7X,''alv_area'',8X,''Pc'',&
            &9X,''Pe'',8X,''Pe_com'',9X,''Tension'',9X,''Pc_com'',9X,&
            &''smoothed_Pc_com'',7X,''Smoo_Pc_com'',7X,''aci_flow'')')
       write(*,'(3X,''(s)'',4X,''(mL/s)'',3X,''(mL)'',1X,''(cmH/L.s)'',&
            &1X,''(L/cmH)'',1X,''(...cmH2O...)'',&
            &4X,''(L)'',5X,''(......cmH2O.......)'')')

       write(*,'(F9.2,2(F8.1),8(F9.2),8(F18.8),10(F18.6))') &
            0.0_dp,0.0_dp,0.0_dp, &  !time, flow, tidal
            elem_field(ne_t_resist,1)*1.0e+6_dp/98.0665_dp, & !res (cmH2O/L.s)
            totalC*98.0665_dp/1.0e+6_dp, & !total model compliance
            ppl_current, & !Ppl (cmH2O)/98.0665_dp
            -ppl_current, & !mean Ptp (cmH2O) /98.0665_dp
            init_vol/1.0e+6_dp, & !total model volume (L)
            totalV, &
            0.0_dp, & !Pmuscle (cmH2O)
            Pcw/98.0665_dp, & !Pchest_wall (cmH2O)
            (-Pcw)/98.0665_dp, & !Pmuscle - Pchest_wall (cmH2O)
!            unit_field(nu_vol,19272), &
!            unit_field(nu_vol,2000), &
!            Precoil, &

!            surf_concentration(nu_vol,19272), & !g/cm^2
!            alv_area_pre(nu_vol,1),&
!            Pc(nu_vol,19272,&  !alv_collapse_pressure (cmH2O)
!            unit_field(nu_air_press,19272),&
!            Paw,&
            unit_field(nu_pe,19272), &
!            unit_field(nu_comp,1), & !unit_field(nu_vol,1), &
!            surface_tension(nu_vol,1),&
!            unit_field(nu_intsat,1),&
!            unit_field(nu_time, 1),&
!            lym_condition(2),&
            interstitial_pressure_a (nu_Pe,19272),&
            interstitial_pressure_b (nu_Pe,19272),&
            flux_b(nu_av_flux,19272),&
!            alv_radii_current(nu_vol,19272),&
!            alv_unit_field_current(nu_vol,19272), & !mm^3
!            alveolar_volume(nu_alvflow,19272),&
            P_initial_lymphtix (nu_Pe,19272),&
!            unit_field(nu_blood_press,19272),&
            unit_field(nu_lymphflow,19272),&
            unit_field(nu_av_flux,19272),&
            interstitial_saturation(nu_intsat, 19272), &
            unit_field(nu_time, 19272), &
            unit_field(nu_tt,19272),&
!            unit_active_time(19727), &
            alveolar_volume (nu_alvflow,19272), &
            Precoil/0.0042_dp, &
            surface_tension(nu_vol,19272),&
            osm_cap(nu_osmflux, 19272)

!            unit_field(nu_Vdot0,1)!,&
            !unit_field_pre(nu_vol,1),&
            !unit_field_current(nu_vol,1)
    else
       write(*,'(F9.2,2(F8.1),8(F9.2),8(F18.8),10(F18.6))') &
            time, & !time through breath (s)
            elem_field(ne_Vdot,1)/1.0e+3_dp, & !flow at the inlet (mL/s)
            (current_vol - init_vol)/1.0e+3_dp, & !current tidal volume (mL)
            elem_field(ne_t_resist,1)*1.0e+6_dp/98.0665_dp, & !res (cmH2O/L.s)
            totalC*98.0665_dp/1.0e+6_dp, & !total model compliance
            ppl_current, & !Ppl (cmH2O) /98.0665_dp
            pptrans, & !mean Ptp (cmH2O) /98.0665_dp
            current_vol/1.0e+6_dp, & !total model volume (L)
            totalV, &
            p_mus/98.0665_dp, & !Pmuscle (cmH2O)
            -Pcw/98.0665_dp, & !Pchest_wall (cmH2O)
            (p_mus+Pcw)/98.0665_dp, & !Pmuscle - Pchest_wall (cmH2O)
!            unit_field(nu_vol,19272), &
!            unit_field(nu_vol,2000), &
!            Precoil, &

!            surf_concentration(nu_vol,19272), & !g/cm^2
!            alv_area_pre(nu_vol,1),&
!            Pc(nu_vol,19272),&  !alv_collapse_pressure (cmH2O)
!            unit_field(nu_air_press,19272),&
!            Paw,&
            unit_field(nu_pe,19272), &
!            unit_field(nu_comp,1), & !unit_field(nu_vol,1), &
!            surface_tension(nu_vol,1),&
!            unit_field(nu_intsat,1),&
!            unit_field(nu_time, 1),&
!            lym_condition(2),&
            interstitial_pressure_a (nu_Pe,19272),&
            interstitial_pressure_b (nu_Pe,19272),&
            flux_b(nu_av_flux,19272),&
!            alv_radii_current(nu_vol,19272),&
!            alv_unit_field_current(nu_vol,19272), & !mm^3
!            alveolar_volume(nu_alvflow,19272),&
            P_initial_lymphtix (nu_Pe,19272),&
!            unit_field(nu_blood_press,19272),&
            unit_field(nu_lymphflow,19272),&
            unit_field(nu_av_flux,19272),&
            interstitial_saturation(nu_intsat, 19272), &
            unit_field(nu_time, 19272), &
            unit_field(nu_tt,19272),&
            alveolar_volume (nu_alvflow,19272), &
            Precoil/0.0042_dp, &
            surface_tension(nu_vol,19272),&
            osm_cap(nu_osmflux, 19272)
            !alv_radii(nu_vol,1), &
!            node_field(nj_aw_press,1), &
!            unit_field(nu_Vdot0,1)!,&
            !unit_field_pre(nu_vol,1),&
            !unit_field_current(nu_vol,1)

    endif

  end subroutine write_flow_step_results_surf

  subroutine write_flow_step_results(init_vol, &
       current_vol,ppl_current,pptrans,Pcw,p_mus,time,ttime)

    use parameter_types, only: lung_params

    real(dp),intent(in) :: init_vol,current_vol, &
         ppl_current,pptrans,Pcw,p_mus,time,ttime
    ! Local variables
    integer :: i, out_unit(2)
    real(dp) :: totalC,Precoil

    ! --------------------------------------------------------------------------

    !the total model compliance
    totalC = 1.0_dp/(1.0_dp/sum(unit_field(nu_comp,1:num_units))+ &
         1.0_dp/lung_params%chest_wall_compliance)
    Precoil = sum(unit_field(nu_pe,1:num_units))/num_units

    out_unit = [6, 10] ! 6 is stdout, 10 is my unit
    
    if(abs(time).lt.zero_tol)then
!!! write out the header information for run-time output
       write(*,'(2X,''Time'',3X,''Inflow'',4X,''V_t'',5X,''Raw'',5X,&
            &''Comp'',4X,''Ppl'',5X,''Ptp'',5X,''VolL'',4X,''Pmus'',&
            &4X,''Pcw'',2X,''Pmus-Pcw'')')
       write(*,'(3X,"(s)",4X,"(",A,")",3X,"(",A,")",1X,"(cmH/",A,".s)", &
            1X,"(",A,"/cmH)",1X,"(...cmH2O...)", &
            4X,"(",A,")",5X,"(......cmH2O.......)")') &
            trim(units_flow), trim(units_ds), trim(units_vol), &
            trim(units_vol), trim(units_vol)

       do i = 1,2
          write(out_unit(i),'(f7.3, 2(f8.1), 8(f8.2))') &
               0.0_dp, 0.0_dp, 0.0_dp, &  !time, flow, tidal
               elem_field(ne_t_resist,1) * vol_scale_op / Pa_cmH2O, & !res (cmH2O/L.s)
               totalC * Pa_cmH2O / vol_scale_op, & !total model compliance
               ppl_current / Pa_cmH2O, & !Ppl (cmH2O)
               -ppl_current / Pa_cmH2O, & !mean Ptp (cmH2O)
               init_vol / vol_scale_op, & !total model volume (L)
               0.0_dp, & !Pmuscle (cmH2O)
               Pcw / Pa_cmH2O, & !Pchest_wall (cmH2O)
               (-Pcw) / Pa_cmH2O !Pmuscle - Pchest_wall (cmH2O)
       enddo
    else
       do i = 1,2
          write(out_unit(i),'(F7.3,2(F8.1),8(F8.2))') &
               time, & !time through breath (s)
               elem_field(ne_Vdot,1) / flow_scale_op, & !flow at the inlet (mL/s)
               (current_vol - init_vol) / ds_scale_op, & !current tidal volume (mL)
               elem_field(ne_t_resist,1) * vol_scale_op / Pa_cmH2O, & !res (cmH2O/L.s)
               totalC * Pa_cmH2O / vol_scale_op, & !total model compliance
               ppl_current / Pa_cmH2O, & !Ppl (cmH2O)
               pptrans / Pa_cmH2O, & !mean Ptp (cmH2O)
               current_vol / vol_scale_op, & !total model volume (L)
               p_mus / Pa_cmH2O, & !Pmuscle (cmH2O)
               -Pcw / Pa_cmH2O, & !Pchest_wall (cmH2O)
               (p_mus+Pcw) / Pa_cmH2O !Pmuscle - Pchest_wall (cmH2O)
       enddo
    endif

  end subroutine write_flow_step_results

!!!#############################################################################

  function ventilation_continue(n,sum_tidal)

    use parameter_types, only: solve_V_params, V_params

    integer,intent(in) :: n
    real(dp),intent(in) :: sum_tidal
    ! Local variables
    logical :: ventilation_continue

    ! --------------------------------------------------------------------------

    ventilation_continue = .true.
    if(n >= solve_V_params%num_breaths)then
       ventilation_continue = .false.
    elseif(abs(V_params%tidal_volume).gt.1.0e-3_dp)then
       if(abs(100.0_dp*(V_params%tidal_volume-sum_tidal) &
            /V_params%tidal_volume).gt.0.1_dp.or.(n.lt.2))then
          ventilation_continue = .true.
       else
          ventilation_continue = .false.
       endif
    endif

  end function ventilation_continue
  
!!!#############################################################################
  
  logical function lymph_continue() result(keep_going)
    
    use parameter_types, only: solve_V_params, V_params

    integer :: nunit
    logical :: any_unit_running
    
    keep_going = .true.
    
    ! Check if any unit is still running
    any_unit_running = .false.
    do nunit = 1, num_units
       if ((lym_condition(nunit) > 0.00001_dp .or. &
            unit_active_time(nunit) < 200.0_dp * unit_field(nu_tt,nunit)) .and. &
            unit_active_time(nunit) < 5000.0_dp * unit_field(nu_tt,nunit)) then
          any_unit_running = .true.
          exit
       endif
    enddo

    ! If all units stopped
    if (.not. any_unit_running) then
       keep_going = .false.
    endif
    
  end function lymph_continue
  
!!!#############################################################################
  
end module ventilation

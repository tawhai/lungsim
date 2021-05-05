module geometry
  !*Brief Description:* This module handles all geometry read/write/generation.
  !
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !
  !This module handles all geometry read/write/generation.
  use arrays
  use diagnostics
  use exports
  use indices
  use mesh_utilities
  use other_consts ! currently has pi
  use precision ! sets dp for precision
  
  implicit none
  
  !Module parameters
  
  !Module types
  
  !Module variables
  
  !Interfaces
  private
  public add_mesh
  public add_matching_mesh
  public append_units
  public define_1d_elements
  public define_elem_geometry_2d
  public define_mesh_geometry_test
  public define_node_geometry
  public define_node_geometry_2d
  public define_data_geometry
  public enclosed_volume
  public define_rad_from_file
  public define_rad_from_geom
  public grow_tree_wrap
  public import_node_geometry_2d
  public get_final_real
  public make_data_grid
  public make_2d_vessel_from_1d
  public merge_2d_element
  public set_initial_volume
  public volume_of_mesh
  public write_geo_file
  public get_final_integer
  public get_four_nodes
  public write_elem_geometry_2d
  public write_node_geometry_2d
  
contains

!!!#############################################################################
  
  subroutine allocate_node_arrays(num_nodes)
    !*allocate_node_arrays:* allocate memory for arrays associated with 1D trees
    
    integer,intent(in) :: num_nodes
    ! Local variables
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'allocate_node_arrays'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(nodes)) allocate (nodes(num_nodes))
    if(.not.allocated(node_xyz)) allocate (node_xyz(3,num_nodes))
    if(.not.allocated(node_field)) allocate (node_field(num_nj,num_nodes))
    if(.not.allocated(elems_at_node)) allocate(elems_at_node(num_nodes,0:3))
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise
    
    call enter_exit(sub_name,2)
    
  end subroutine allocate_node_arrays
  
!!!#############################################################################

  subroutine add_mesh(AIRWAY_MESHFILE)
    !*add_mesh:* Reads in an ipmesh file and adds this mesh to the terminal
    ! branches of an existing tree geometry
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MESH" :: ADD_MESH

    character(len=MAX_FILENAME_LEN), intent(in) :: AIRWAY_MESHFILE
    ! Local parameters
    character(len=100) :: buffer
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line
    integer :: i,ibeg,iend,i_ss_end,j,ne,ne0,ne_global,ne_parent,ne_start, &
         ngen_parent,np,np0,np_global,&
         num_elems_new,num_elems_to_add,num_nodes_new,nunit,nlabel
    integer,dimension(1000) :: element_temp,generation, &
         parent_element,symmetry_temp
    real(dp),dimension(1000) :: length,radius,a_A
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'add_mesh'
    call enter_exit(sub_name,1)

    ios = 0
    line = 0
    open(fh, file=AIRWAY_MESHFILE)
    
    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    ios=0
    line=0
    i=0 ! count the number of elements read in
    
    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       ! line contains: element, parent element, generation,
       !                symmetry, length, outer radius, a/A ratio
       ! note that a/A ratio is always 1 for the conducting airways
       if (ios == 0) then
          line = line + 1
          i=i+1
          i_ss_end = len(buffer)
          
          do nlabel = 1,7
             ibeg = index(buffer," ") + 1 !get location of first integer beyond ws in string
             buffer = adjustl(buffer(ibeg:i_ss_end)) ! get info beyond ws, remove leading ws
             iend = index(buffer," ") !get location of next ws in sub-string
             select case (nlabel)
             case (1)
                read (buffer(1:iend-1), *, iostat=ios) element_temp(i) ! not used??
             case(2)
                read (buffer(1:iend-1), *, iostat=ios) parent_element(i)
             case(3)
                read (buffer(1:iend-1), *, iostat=ios) generation(i)
             case(4)
                read (buffer(1:iend-1), *, iostat=ios) symmetry_temp(i)
             case(5)
                read (buffer(1:iend-1), *, iostat=ios) length(i)
             case(6)
                read (buffer(1:iend-1), *, iostat=ios) radius(i)
             case(7)
                read (buffer(1:iend-1), *, iostat=ios) a_A(i)
             end select
          enddo
       endif
    enddo
    close(fh)
    
    num_elems_to_add = i
    
!!! increase the size of node and element arrays to accommodate the additional elements
    ! the number of nodes after adding mesh will be:
    num_nodes_new = num_nodes + num_units*num_elems_to_add
    ! the number of elems after adding mesh will be:
    num_elems_new = num_elems + num_units*num_elems_to_add
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    
    ne = num_elems ! the starting local element number
    ne_global = elems(ne) ! assumes this is the highest element number (!!!)
    np = num_nodes ! the starting local node number
    np_global = nodes(np) ! assumes this is the highest node number (!!!)
    
    do nunit = 1,num_units ! for all terminal branches, append the mesh
       
       ne_parent = units(nunit) ! local element number of terminal, to append to
       ngen_parent = elem_ordrs(1,ne_parent)
       ne_start = ne !starting element number for the unit
       
       do i=1,num_elems_to_add
          
          if(parent_element(i).eq.0)then
             ne_parent = units(nunit)
          else
             ne_parent = ne_start+parent_element(i)
          endif
          
          ne0 = ne_parent
          np0 = elem_nodes(2,ne0)
          
          ne_global = ne_global + 1 ! new global element number
          ne = ne + 1 ! new local element number
          np_global = np_global + 1 !new global node number
          np = np + 1 ! new local node number
          
          nodes(np) = np_global
          elems(ne) = ne_global
          
          elem_nodes(1,ne) = np0
          elem_nodes(2,ne) = np
          
          elem_ordrs(1,ne) = ngen_parent + generation(i)
          elem_ordrs(no_type,ne) = 1   ! ntype ! 0 for respiratory, 1 for conducting
          elem_symmetry(ne) = symmetry_temp(i)+1 ! uses 0/1 in file; 1/2 in code
          
          ! record the element connectivity
          elem_cnct(-1,0,ne) = 1 ! one parent branch
          elem_cnct(-1,1,ne) = ne0 ! store parent element
          elem_cnct(1,0,ne0) = elem_cnct(1,0,ne0) + 1
          elem_cnct(1,elem_cnct(1,0,ne0),ne0) = ne
          
          ! record the direction and location of the branch
          do j=1,3
             elem_direction(j,ne) = elem_direction(j,ne0)
             node_xyz(j,np) = node_xyz(j,np0) + &
                  elem_direction(j,ne)*length(i)
          enddo !j
          
          elem_field(ne_length,ne) = length(i)
          elem_field(ne_radius,ne) = radius(i)
          elem_field(ne_a_A,ne) = a_A(i)
          elem_field(ne_vol,ne) = PI*radius(i)**2*length(i)
          
       enddo !i
    enddo !nunit
    
    num_nodes = np
    num_elems = ne
    
    call element_connectivity_1d
    call evaluate_ordering ! calculate new ordering of tree
    
    call enter_exit(sub_name,2)

  end subroutine add_mesh

!!!#############################################################################

  subroutine add_matching_mesh()
    !*add_matching_mesh:* 
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MATCHING_MESH" :: ADD_MATCHING_MESH

    !Parameters to become inputs
    real(dp) :: offset(3)
    logical :: REVERSE=.TRUE.
    character(len=60) :: mesh_type='terminal'
    !local variables
    integer :: num_nodes_new,num_elems_new,ne,ne_global,np,np_global,np0, &
         nonode,np_m
    integer :: nj,ne_m,noelem,ne0,n,nindex,ne1,noelem0,nu,cap_conns, &
         cap_term,np1,np2
    integer, allocatable :: np_map(:)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'add_matching_mesh'
    call enter_exit(sub_name,1)
    !Ultimately offset should be an input argument
    offset(1)=0.0_dp
    offset(2)=1.0e-6_dp
    offset(3)=0.0_dp

    allocate(np_map(num_nodes))
!!! increase the size of node and element arrays to accommodate the additional elements
    ! the number of nodes after adding mesh will be:
    num_nodes_new = 2*num_nodes
    ! the number of elems after adding mesh will be:
    if(mesh_type.eq.'basic')then
       num_elems_new = 2*num_elems
    elseif(mesh_type.eq.'terminal')then
       num_elems_new = 2*num_elems + num_units
    elseif(mesh_type.eq.'ladder')then
       print *, "NOT YET IMPLEMENTED"
       call exit(0)
    endif
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    noelem0=0
    ne0 = num_elems ! the starting local element number
    ne_global = elems(ne0) ! assumes this is the highest element number (!!!)
    np0 = num_nodes ! the starting local node number
    np_global = nodes(np0) ! assumes this is the highest node number (!!!)
    
    do nonode=1,num_nodes
       np=np_global+nonode
       np_m=nodes(nonode)
       np_map(np_m)=np !maps new to old node numbering
       nodes(np0+nonode)=np
       do nj=1,3
          node_xyz(nj,np)=node_xyz(nj,np_m)+offset(nj)
       enddo
       elems_at_node(np,0)=0 !initialise
       !Doesnt map versions, would be added here
    enddo
    
    do noelem=1,num_elems
       ne=ne_global+noelem
       elem_field(ne_group,ne)=2.0_dp!VEIN
       ne_m=elems(noelem)
       elem_field(ne_group,ne_m)=0.0_dp!ARTERY
       elems(ne0+noelem)=ne
       if(.NOT.REVERSE)then
          elem_nodes(1,ne)=np_map(elem_nodes(1,ne_m))
          elem_nodes(2,ne)=np_map(elem_nodes(2,ne_m))
          elem_cnct(1,0,ne)=elem_cnct(1,0,ne_m)!The numberdownstream are the number downstream
          elem_cnct(-1,0,ne)=elem_cnct(-1,0,ne_m)
          do n=1,elem_cnct(1,0,ne)
             elem_cnct(1,n,ne)=elem_cnct(1,n,ne_m)+ne0
          enddo
          do n=1,elem_cnct(-1,0,ne)
             elem_cnct(-1,n,ne)=elem_cnct(-1,n,ne_m)+ne0
          enddo
       else
          elem_nodes(1,ne)=np_map(elem_nodes(2,ne_m))
          elem_nodes(2,ne)=np_map(elem_nodes(1,ne_m))
          elem_cnct(-1,0,ne)=elem_cnct(1,0,ne_m) !The number upstream are the number downstream
          elem_cnct(1,0,ne)=elem_cnct(-1,0,ne_m)!The number downstream are the number upstream
          do n=1,elem_cnct(1,0,ne)
             elem_cnct(1,n,ne)=elem_cnct(-1,n,ne_m)+ne0
          enddo
          do n=1,elem_cnct(-1,0,ne)
             elem_cnct(-1,n,ne)=elem_cnct(1,n,ne_m)+ne0
          enddo
       endif
       !if worrying about regions and versions do it here
       elems_at_node(elem_nodes(1,ne),0)=elems_at_node(elem_nodes(1,ne),0)+1
       elems_at_node(elem_nodes(1,ne),elems_at_node(elem_nodes(1,ne),0))=ne
       elems_at_node(elem_nodes(2,ne),0)=elems_at_node(elem_nodes(2,ne),0)+1
       elems_at_node(elem_nodes(2,ne),elems_at_node(elem_nodes(2,ne),0))=ne
       nindex=no_gen
       elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
       nindex=no_sord
       elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
       nindex=no_hord
       elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
    enddo
    
    !update current no of nodes and elements to determine connectivity
    np0=np !current highest node
    ne1=ne !current highest element
    noelem0=num_elems+noelem0
    if(mesh_type.eq.'ladder')then
       !To be implemented
    elseif(mesh_type.eq.'terminal')then
       cap_conns=0
       cap_term=0
       do nu=1,num_units
          ne=units(nu)
          cap_term=cap_term+1
          np1=elem_nodes(2,ne)
          np2=np_map(np1)
          noelem0=noelem0+1
          ne1=ne1+1
          elems(noelem0)=ne1
          elem_nodes(1,ne1)=np1
          elem_nodes(2,ne1)=np2
          elems_at_node(np1,0)=elems_at_node(np1,0)+1
          elems_at_node(np1,elems_at_node(np1,0))=ne1
          elems_at_node(np2,0)=elems_at_node(np2,0)+1
          elems_at_node(np2,elems_at_node(np2,0))=ne1
          elem_cnct(1,elem_cnct(1,0,ne)+1,ne)=ne1
          elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
          elem_cnct(-1,elem_cnct(-1,0,ne+ne_global)+1,ne+ne_global)=ne1
          elem_cnct(-1,0,ne+ne_global)=elem_cnct(-1,0,ne+ne_global)+1
          elem_cnct(-1,0,ne1)=1
          elem_cnct(1,0,ne1)=1
          elem_cnct(-1,1,ne1)=ne
          elem_cnct(1,1,ne1)=ne+ne0
          nindex=no_gen
          elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
          nindex=no_sord
          elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
          nindex=no_hord
          elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
          elem_field(ne_group,ne1)=1.0_dp!connection between meshes
       enddo
       print *, 'Number of connections', cap_term
    endif
    num_nodes=num_nodes_new
    num_elems=num_elems_new
    deallocate(np_map)

    call enter_exit(sub_name,2)
    
  end subroutine add_matching_mesh

!!!#############################################################################

  subroutine append_units()
    !*append_units:* Appends terminal units at the end of a tree structure
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_APPEND_UNITS" :: APPEND_UNITS

    ! Local parameters
    integer :: ne,ne0,nu
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'append_units'
    call enter_exit(sub_name,1)

    num_units = 0
    do ne=1,num_elems
       if(elem_cnct(1,0,ne).eq.0)THEN
          num_units=num_units+1
       endif
    enddo
    
    if(allocated(units))then !increasing the array size; just overwrite
       deallocate(units)
       deallocate(unit_field)
    endif
    allocate(units(num_units))
    allocate(unit_field(num_nu,num_units))
    
    unit_field=0.0_dp
    units=0
    elem_units_below(1:num_elems) = 0 !initialise the number of terminal units below a branch
    
    nu=0
    do ne=1,num_elems
       if(elem_cnct(1,0,ne).eq.0)THEN
          nu=nu+1
          units(nu)=ne     !Set up units array containing terminals
          elem_units_below(ne)=1
       endif
    enddo
    
    ! count the effective number of elements below each branch
    do ne=num_elems,2,-1
       ne0=elem_cnct(-1,1,ne)
       elem_units_below(ne0) = elem_units_below(ne0) &
            + elem_units_below(ne)*elem_symmetry(ne)
    enddo !ne
    
    call enter_exit(sub_name,2)

  end subroutine append_units

!!!#############################################################################

  subroutine define_1d_elements(ELEMFILE)
    !*define_1d_elements:* Reads in an 1D element ipelem file to define a geometry
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_1D_ELEMENTS" :: DEFINE_1D_ELEMENTS
    
    character(len=MAX_FILENAME_LEN), intent(in) :: ELEMFILE
    !     Local Variables
    integer :: ibeg,iend,ierror,i_ss_end,j,ne,ne_global,&
         nn,np,np1,np2,np_global
    character(LEN=132) :: ctemp1
    character(len=300) :: readfile
    character(LEN=40) :: sub_string
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_1d_elements'
    call enter_exit(sub_name,1)
    
    if(index(ELEMFILE, ".ipelem")> 0) then !full filename is given
       readfile = ELEMFILE
    else ! need to append the correct filename extension
       readfile = trim(ELEMFILE)//'.ipelem'
    endif
    
    open(10, file=readfile, status='old')
    
    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          num_elems = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_elements
       endif
    end do read_number_of_elements
    
!!! allocate memory for element arrays
    if(allocated(elems)) deallocate(elems)
    allocate(elems(num_elems))
    if(allocated(elem_cnct)) deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems))
    if(allocated(elem_nodes)) deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems))
    if(allocated(elem_ordrs)) deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems))
    if(allocated(elem_symmetry)) deallocate(elem_symmetry)
    allocate(elem_symmetry(num_elems))
    if(allocated(elem_units_below)) deallocate(elem_units_below)
    allocate(elem_units_below(num_elems))
    if(allocated(elems_at_node)) deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes,0:3))
    if(allocated(elem_field)) deallocate(elem_field)
    allocate(elem_field(num_ne,num_elems))
    if(allocated(elem_direction)) deallocate(elem_direction)
    allocate(elem_direction(3,num_elems))
    if(model_type.eq.'gas_mix')then
       if(allocated(expansile)) deallocate(expansile)
       allocate(expansile(num_elems))
    endif
    
!!! initialise element arrays
    elems = 0
    elem_nodes = 0
    elem_ordrs = 0  ! where the default is that 0==respiratory and 1==conducting
    elem_symmetry = 1
    elem_field = 0.0_dp
    if(model_type.eq.'gas_mix')expansile = .false.
    
    ne=0
    
    read_an_element : do
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          ne_global = get_final_integer(ctemp1) !return the final integer
          ne=ne+1
          elems(ne)=ne_global
          
          read_element_nodes : do
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "global")> 0) then !found the correct line
                iend=len(ctemp1)
                ibeg=index(ctemp1,":")+1 !get location of first integer in string
                sub_string = adjustl(ctemp1(ibeg:iend)) ! get the characters beyond : remove leading blanks
                i_ss_end=len(sub_string) !get the end location of the sub-string
                ibeg=1
                do nn=1,2
                   iend=index(sub_string," ") !get location of first blank in sub-string
                   read (sub_string(ibeg:iend-1), '(i7)' ) np_global
                   call get_local_node(np_global,np) ! get local node np for global node
                   elem_nodes(nn,ne)=np ! the local node number, not global
                   sub_string = adjustl(sub_string(iend:i_ss_end)) ! get chars beyond blank, remove leading blanks
                end do
                exit read_element_nodes
             endif !index
          end do read_element_nodes
          if(ne.ge.num_elems) exit read_an_element
       endif
       
    end do read_an_element
    
    close(10)
    
    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = sqrt((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)
       enddo !j
    enddo
    
    call element_connectivity_1d
    call evaluate_ordering

    elem_ordrs(no_type,:) = 1 ! 0 for respiratory, 1 for conducting
    
    call enter_exit(sub_name,2)

  end subroutine define_1d_elements

!!!#############################################################################

  subroutine define_elem_geometry_2d(ELEMFILE,sf_option)
    ! Reads in 2D ipelem file.
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_ELEM_GEOMETRY_2D" :: DEFINE_ELEM_GEOMETRY_2D

    character(len=*) :: ELEMFILE
    character(len=4) :: sf_option
    !     Local Variables
    integer :: ierror,ne,ne_global,nn,np,number_of_elements
    character(len=132) :: ctemp1,readfile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_elem_geometry_2d'
    call enter_exit(sub_name,1)
    
    if(index(ELEMFILE, ".ipelem")> 0) then !full filename is given
       readfile = ELEMFILE
    else ! need to append the correct filename extension
       readfile = trim(ELEMFILE)//'.ipelem'
    endif
    
    open(10, file=readfile, status='old')
    
    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          number_of_elements = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_elements
       endif
    end do read_number_of_elements
    
    num_elems_2d=number_of_elements
    if(allocated(elems_2d))then
       deallocate(elems_2d)
       deallocate(elem_nodes_2d)
       deallocate(elem_versn_2d)
    endif
    allocate(elems_2d(num_elems_2d))
    allocate(elem_nodes_2d(4,num_elems_2d))
    allocate(elem_versn_2d(4,num_elems_2d))
    
    ne = 0
    
    read_an_element : do 
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          ne_global = get_final_integer(ctemp1) !return the final integer
          ne = ne + 1
          elems_2d(ne) = ne_global
          
          read_element_nodes : do 
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "global")> 0) then !found the correct line
                call get_four_nodes(ne,ctemp1) !number of versions for node np
                ! note that only the ne'th data of elem_nodes_2d is passed to 'get_four_nodes'
                do nn=1,4
                   np=elem_nodes_2d(nn,ne)
                   if(node_versn_2d(np).gt.1)then
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      elem_versn_2d(nn,ne) = get_final_integer(ctemp1) !return the final integer
                   else
                      elem_versn_2d(nn,ne)= 1
                   endif !nversions
                enddo !nn
                exit read_element_nodes
             endif !index
          end do read_element_nodes
          
          if(ne.ge.number_of_elements) exit read_an_element
       endif
       
    end do read_an_element
    
    close(10)
    
    call element_connectivity_2d
    call line_segments_for_2d_mesh(sf_option)
    
    call enter_exit(sub_name,2)
    
  end subroutine define_elem_geometry_2d
  
!!!#############################################################################

  subroutine define_mesh_geometry_test()
    !*define_mesh_geometry_test:*
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_MESH_GEOMETRY_TEST" :: DEFINE_MESH_GEOMETRY_TEST

    !     Local Variables
    integer :: j,ne,np,np1,np2
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_mesh_geometry_test'
    call enter_exit(sub_name,1)
    
    num_elems = 400
    num_nodes = num_elems + 1
    
!!! allocate memory
    if(.not.allocated(nodes)) allocate (nodes(num_nodes))
    if(.not.allocated(node_xyz)) allocate (node_xyz(3,num_nodes))
    if(.not.allocated(node_field)) allocate (node_field(num_nj,num_nodes))
    if(.not.allocated(elems)) allocate(elems(num_elems))
    if(.not.allocated(elem_cnct)) allocate(elem_cnct(-1:1,0:2,0:num_elems))
    if(.not.allocated(elem_nodes)) allocate(elem_nodes(2,num_elems))
    if(.not.allocated(elem_ordrs)) allocate(elem_ordrs(num_ord,num_elems))
    if(.not.allocated(elem_symmetry)) allocate(elem_symmetry(num_elems))
    if(.not.allocated(elem_units_below)) allocate(elem_units_below(num_elems))
    if(.not.allocated(elems_at_node)) allocate(elems_at_node(num_nodes,0:3))
    if(.not.allocated(elem_field)) allocate(elem_field(num_ne,num_elems))
    if(.not.allocated(elem_direction)) allocate(elem_direction(3,num_elems))
    if(.not.allocated(expansile)) allocate(expansile(num_elems))
    
!!! initialise array values
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise
    elems=0
    elem_nodes=0
    elem_symmetry = 1
    elem_field = 0.0_dp
    expansile = .false.
    
!!! set up node arrays
    nodes(1) = 1
    do np=2,101
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 1.0_dp
    enddo
    
    np=102
    node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    do np=102,151
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo
    
    np=152
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,101) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,101) - 0.5_dp
    do np=153,201
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo
    
    np=202
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,151) - 0.5_dp
    node_xyz(3,np) = node_xyz(3,151) - 0.5_dp
    do np=203,251
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo
    
    np=252
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,151) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,151) - 0.5_dp
    do np=253,301
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo
    
    np=302
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,201) - 0.5_dp
    node_xyz(3,np) = node_xyz(3,201) - 0.5_dp
    do np=303,351
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo
    
    np=352
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,201) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,201) - 0.5_dp
    do np=353,401
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo
    
!!! set up element arrays
    do ne=1,num_elems
       elems(ne) = ne
       elem_nodes(1,ne) = ne
       elem_nodes(2,ne) = ne+1
    enddo
    
    elem_nodes(1,151) = 101
    elem_nodes(1,201) = 151
    elem_nodes(1,251) = 151
    elem_nodes(1,301) = 201
    elem_nodes(1,351) = 201
    
    elem_field(ne_radius,1:100) = 10.0_dp
    elem_field(ne_radius,101:200) = 5.0_dp
    elem_field(ne_radius,201:400) = sqrt(elem_field(ne_radius,101)**2/2.0_dp)
    
    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = sqrt((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)
       enddo !j
       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
       elem_field(ne_a_A,ne) = 1.0_dp ! set default for ratio a/A
    enddo
    
    call element_connectivity_1d
    call evaluate_ordering
    
    call enter_exit(sub_name,2)
    
  end subroutine define_mesh_geometry_test

!!!#############################################################################
  
  subroutine define_node_geometry(NODEFILE)
    !*define_node_geometry:* Reads in an ipnode file to define a tree geometry
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY" :: DEFINE_NODE_GEOMETRY
    
    character(len=MAX_FILENAME_LEN), intent(in) :: NODEFILE !Input nodefile
    !     Local Variables
    integer :: i,ierror,np,np_global,num_nodes_temp,num_versions,nv,NJT=0
    real(dp) :: point
    logical :: overwrite = .false. ! initialised
    character(len=300) :: ctemp1,readfile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_node_geometry'
    call enter_exit(sub_name,1)
    
    if(index(NODEFILE, ".ipnode")> 0) then !full filename is given
       readfile = NODEFILE
    else ! need to append the correct filename extension
       readfile = trim(NODEFILE)//'.ipnode'
    endif
    
    open(10, file=readfile, status='old')
    
    if(num_nodes.gt.0) overwrite = .true.
    
    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          num_nodes_temp = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_nodes !exit the named do loop
       endif
    end do read_number_of_nodes
    
    if(.not.overwrite) call allocate_node_arrays(num_nodes_temp) ! don't allocate if just overwriting
    
    !.....read in the number of coordinates
    read_number_of_coords : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "coordinates")> 0) then !keyword "coordinates" is found
          NJT = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_coords !exit the named do loop
       endif
    end do read_number_of_coords
    
    ! note that only the first version of coordinate is currently read in   
    
    !.....read the coordinate, derivative, and version information for each node. 
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          np_global = get_final_integer(ctemp1) !get node number
          
          np = np+1
          nodes(np) = np_global
          !.......read coordinates
          do i=1,3 ! for the x,y,z coordinates
             num_versions=1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "versions")> 0) then
                num_versions = get_final_integer(ctemp1)
                if(num_versions > 1)then
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   point = get_final_real(ctemp1)
                   do nv=2,num_versions
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   enddo
                else
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   point = get_final_real(ctemp1)
                endif
             else ! no prompting for versions
                point = get_final_real(ctemp1)
             endif
             node_xyz(i,np)=point
          end do !i
       endif !index
       if(np.ge.num_nodes_temp) exit read_a_node
    end do read_a_node
    
    if(.not.overwrite) num_nodes = num_nodes_temp
    
    close(10)
    
    call enter_exit(sub_name,2)
    
  end subroutine define_node_geometry

!!!#############################################################################

  subroutine define_node_geometry_2d(NODEFILE)
    !*define_node_geometry_2d:* Reads in an ipnode file to define surface nodes
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY_2D" :: DEFINE_NODE_GEOMETRY_2D
    
    character(len=*),intent(in) :: NODEFILE
    !     Local Variables
    integer :: i,ierror,np,np_global,&
         num_versions,nv
    integer,parameter :: num_derivs = 3
    character(len=132) :: ctemp1,readfile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_node_geometry_2d'
    call enter_exit(sub_name,1)
    
    if(index(NODEFILE, ".ipnode")> 0) then !full filename is given
       readfile = NODEFILE
    else ! need to append the correct filename extension
       readfile = trim(NODEFILE)//'.ipnode'
    endif
    
    open(10, file=readfile, status='old')
    
    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          num_nodes_2d = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_nodes !exit the named do loop
       endif
    end do read_number_of_nodes
    
!!!allocate memory to arrays that require node number
    if(allocated(nodes_2d))then ! deallocate
       deallocate(nodes_2d)
       deallocate(node_xyz_2d)
       deallocate(node_versn_2d)
    endif
    allocate(nodes_2d(num_nodes_2d))
    allocate(node_xyz_2d(4,10,3,num_nodes_2d))
    allocate(node_versn_2d(num_nodes_2d))
    
    !.....read the coordinate, derivative, and version information for each node. 
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          np_global = get_final_integer(ctemp1) !get node number
          
          np=np+1
          nodes_2d(np) = np_global
          
          !.......read coordinates
          do i=1,3 ! for the x,y,z coordinates
             num_versions = 0
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "versions")> 0) num_versions = &
                  get_final_integer(ctemp1)
             node_versn_2d(np) = max(1,num_versions) !number of versions for node np
             do nv=1,node_versn_2d(np)
                if(num_versions > 1)then
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 ! "For version number..."
                endif
                !...........coordinate          
                if(num_versions > 0) &
                     read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                node_xyz_2d(1,nv,i,np) = get_final_real(ctemp1)
                if(num_derivs.ge.1)then
                   !..........derivative 1
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(2,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.2)then
                   !..........derivative 2
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(3,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.3)then
                   !...........derivative 1&2
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(4,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.4)then
                   write(*,'(''This code is only valid for a surface geometry'')')
                   read(*,*)
                endif
             end do !nv
          end do !i
       endif !index
       if(np.ge.num_nodes_2d) exit read_a_node
    end do read_a_node
    
    close(10)
    
    call enter_exit(sub_name,2)
    
  end subroutine define_node_geometry_2d

!!!#############################################################################

  subroutine import_node_geometry_2d(NODEFILE)
    !*define_node_geometry_2d:* Reads in an exnode file to define surface nodes
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY_2D" :: DEFINE_NODE_GEOMETRY_2D
    
    character(len=*),intent(in) :: NODEFILE
    !     Local Variables
    integer :: i,ierror,index_location,np,np_global,num_versions,nv
    character(len=132) :: ctemp1,readfile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'import_node_geometry_2d'
    call enter_exit(sub_name,1)
    
    if(index(NODEFILE, ".exnode")> 0) then !full filename is given
       readfile = NODEFILE
    else ! need to append the correct filename extension
       readfile = trim(NODEFILE)//'.exnode'
    endif
    
    open(10, file=readfile, status='old')

    !.....get the total number of nodes.
    num_nodes_2d = 0
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(ierror<0) exit !ierror<0 means end of file
       if(index(ctemp1, "Node:")> 0) then !keyword "Node:" is found in ctemp1
          num_nodes_2d = num_nodes_2d+1
       endif
    end do read_number_of_nodes
    close(10)

!!!allocate memory to arrays that require node number
    if(.not.allocated(nodes_2d)) allocate(nodes_2d(num_nodes_2d))
    if(.not.allocated(node_xyz_2d)) allocate(node_xyz_2d(4,10,3,num_nodes_2d))
    if(.not.allocated(node_versn_2d)) allocate(node_versn_2d(num_nodes_2d))
    nodes_2d = 0
    node_xyz_2d = 0.0_dp
    node_versn_2d = 0
    
    !.....read the coordinate, derivative, and version information for each node. 
    open(10, file=readfile, status='old')
    np = 0
    num_versions = 1
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Derivatives") > 0)then
          index_location = index(ctemp1, "Versions")
          if(index_location > 0) then
             read(ctemp1(index_location+9:index_location+10), '(i2)', iostat=ierror) num_versions
          else
             num_versions = 1  ! the default
          endif
       endif
       if(index(ctemp1, "Node:")> 0) then
          np_global = get_final_integer(ctemp1) !get node number
          np = np+1
          nodes_2d(np) = np_global
          node_versn_2d(np) = num_versions
          
          !.......read coordinates
          do i =1,3 ! for the x,y,z coordinates
             do nv = 1,node_versn_2d(np)
                read(unit=10, fmt=*, iostat=ierror) node_xyz_2d(1:4,nv,i,np)
             end do !nv
          end do !i
          if((np_global.eq.80).or.(np_global.eq.56))then
             node_xyz_2d(3,5,:,np) = 0.0_dp
          elseif(np_global.eq.52)then
             node_xyz_2d(2,4,:,np) = 0.0_dp
          elseif(np_global.eq.94)then
             node_xyz_2d(2,3,:,np) = 0.0_dp
             node_xyz_2d(2,4,:,np) = 0.0_dp
          elseif(np_global.eq.96)then
             node_xyz_2d(2,5,:,np) = 0.0_dp
          endif
       endif !index
       if(np.ge.num_nodes_2d) exit read_a_node
    end do read_a_node
    
    close(10)
    write(*,*) 'warning: hardcoded zeroing of derivs at nodes 52,56,80,94,96'
    
    call enter_exit(sub_name,2)
    
  end subroutine import_node_geometry_2d

!!!#############################################################################

  subroutine define_data_geometry(datafile)
    !*define_data_geometry:* reads data points from a file
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_DATA_GEOMETRY" :: DEFINE_DATA_GEOMETRY

    character(len=*) :: datafile
    ! Local variables
    integer :: iend,ierror,length_string,ncount,nj,itemp
    character(len=132) :: buffer,readfile
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'define_data_geometry'
    call enter_exit(sub_name,1)
    
    if(index(datafile, ".ipdata")> 0) then !full filename is given
       readfile = datafile
    else ! need to append the correct filename extension
       readfile = trim(datafile)//'.ipdata'
    endif
    
    open(10, file=readfile, status='old')
    read(unit=10, fmt="(a)", iostat=ierror) buffer

    !set the counted number of data points to zero
    ncount = 0

!!! first run through to count the number of data points
    read_line_to_count : do
       read(unit=10, fmt="(a)", iostat=ierror) buffer
       if(ierror<0) exit !ierror<0 means end of file
       ncount = ncount + 1
    end do read_line_to_count
    num_data = ncount
    close (10)
    write(*,'('' Read'',I7,'' data points from file'')') num_data
    
!!! allocate arrays now that we know the size required
    if(allocated(data_xyz)) deallocate(data_xyz)
    if(allocated(data_weight)) deallocate(data_weight)
    allocate(data_xyz(3,num_data))
    allocate(data_weight(3,num_data))
    
!!! read the data point information
    open(10, file=readfile, status='old')
    read(unit=10, fmt="(a)", iostat=ierror) buffer
    
    !set the counted number of data points to zero
    ncount = 0
    read_line_of_data : do
       
       ! read the data #; z; y; z; wd1; wd2; wd3 for each data point
       read(unit=10, fmt="(a)", iostat=ierror) buffer
       if(ierror<0) exit !ierror<0 means end of file
       length_string = len_trim(buffer) !length of buffer, and removed trailing blanks
       
       ! read data number
       buffer=adjustl(buffer) !remove leading blanks
       iend=index(buffer," ",.false.)-1 !index returns location of first blank
       if(length_string == 0) exit
       ncount=ncount+1
       read (buffer(1:iend), '(i6)') itemp
       
       do nj=1,3
          ! read x,y,z coordinates
          buffer = adjustl(buffer(iend+1:length_string)) !remove data number from string
          buffer = adjustl(buffer) !remove leading blanks
          length_string = len(buffer) !new length of buffer
          iend=index(buffer," ",.false.)-1 !index returns location of first blank
          read (buffer(1:iend), '(D25.17)') data_xyz(nj,ncount)
       enddo !nj
       
       do nj=1,3
          !        ! read weightings
          !        buffer = adjustl(buffer(iend+1:length_string)) !remove data number from string
          !        buffer = adjustl(buffer) !remove leading blanks
          !        length_string = len(buffer) !new length of buffer
          !        iend=index(buffer," ",.false.)-1 !index returns location of first blank
          !        read (buffer(1:iend), '(D25.17)') data_weight(nj,ncount)
          data_weight(nj,ncount)=1.0_dp
       enddo !nj
       
    enddo read_line_of_data
    
    close(10)
    
    call enter_exit(sub_name,2)

  end subroutine define_data_geometry

!!!#############################################################################
  
  subroutine grow_tree_wrap(surface_elems,parent_ne,angle_max,angle_min,&
       branch_fraction,length_limit,shortest_length,rotation_limit)
    !temporary interface to the grow_tree subroutine until bindings sorted out 
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROW_TREE_WRAP" :: GROW_TREE_WRAP

    use growtree,only: grow_tree
    
    integer,intent(in)  :: surface_elems(:)         ! list of surface elements defining the host region
    integer,intent(in)  :: parent_ne                ! stem branch that supplies 'parents' to grow from
    real(dp),intent(in) :: angle_max                ! maximum branch angle with parent; in degrees
    real(dp),intent(in) :: angle_min                ! minimum branch angle with parent; in degrees
    real(dp),intent(in) :: branch_fraction          ! fraction of distance (to COFM) to branch
    real(dp),intent(in) :: length_limit             ! minimum length of a generated branch (shorter == terminal)
    real(dp),intent(in) :: shortest_length          ! length that short branches are reset to (shortest in model)
    real(dp),intent(in) :: rotation_limit           ! maximum angle of rotation of branching plane

    call grow_tree(surface_elems,parent_ne,angle_max,angle_min, &
         branch_fraction,length_limit,shortest_length,rotation_limit)
    
  end subroutine grow_tree_wrap

!!!#############################################################################
  
  !
  subroutine make_data_grid(surface_elems, offset, spacing, filename, groupname)
    !*make_data_grid:* makes a regularly-spaced 3D grid of data points to
    ! fill a bounding surface 
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_MAKE_DATA_GRID" :: MAKE_DATA_GRID
    
    integer,intent(in) :: surface_elems(:)
    real(dp),intent(in) :: offset, spacing
    logical :: to_export = .true.
    character(len=*),intent(in) :: filename
    character(len=*),intent(in) :: groupname
    ! Local Variables
    integer :: i,j,k,ne,nj,nline,nn,num_data_estimate,num_triangles,num_vertices
    integer,allocatable :: elem_list(:),triangle(:,:)
    real(dp) :: cofm1(3),cofm2(3),boxrange(3),max_bound(3),min_bound(3), &
         point_xyz(3),scale_mesh
    real(dp),allocatable :: data_temp(:,:),vertex_xyz(:,:)
    logical :: internal
    character(len=1) :: char1
    character(len=100) :: writefile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'make_data_grid'
    call enter_exit(sub_name,1)

    allocate(elem_list(count(surface_elems.ne.0)))
    do i = 1,count(surface_elems.ne.0)
       elem_list(i) = get_local_elem_2d(surface_elems(i))
    enddo

    call triangles_from_surface(num_triangles,num_vertices,elem_list, &
         triangle,vertex_xyz)

    if(to_export)then
!!! export vertices as nodes
       call export_triangle_nodes(num_vertices,vertex_xyz,filename,groupname)
!!! export the triangles as surface elements
       call export_triangle_elements(num_triangles,triangle,filename,groupname)
    endif
    
    scale_mesh = 1.0_dp-(offset/100.0_dp)
    cofm1 = sum(vertex_xyz,dim=2)/num_vertices
    forall (i = 1:num_vertices) vertex_xyz(1:3,i) = &
         vertex_xyz(1:3,i)*scale_mesh
    cofm2 = cofm1 * scale_mesh
    forall (i = 1:num_vertices) vertex_xyz(1:3,i) = &
         vertex_xyz(1:3,i) - (cofm2(1:3)-cofm1(1:3))

!!! find the bounding coordinates for the surface mesh
    
    min_bound = minval(vertex_xyz,2)
    max_bound = maxval(vertex_xyz,2)
    boxrange = max_bound - min_bound
    num_data_estimate = int(dble((boxrange(1)/spacing+1.0_dp) * &
         (boxrange(2)/spacing+1.0_dp) * (boxrange(3))/spacing+1.0_dp) * &
         volume_internal_to_surface(triangle,vertex_xyz)/ &
         (boxrange(1)*boxrange(2)*boxrange(3)))
    
!!! allocate arrays based on estimated number of data points
    
    if(allocated(data_xyz)) deallocate(data_xyz)
    allocate(data_xyz(3,num_data_estimate))
    i=0
    num_data = 0
    point_xyz = min_bound + 0.5_dp*spacing 
    do while(point_xyz(3).le.max_bound(3)) ! for z direction
       i=i+1
       j=0 
       do while(point_xyz(2).le.max_bound(2)) ! for y direction
          j=j+1
          k=0
          internal = .true.
          do while(point_xyz(1).le.max_bound(1)) ! for x direction
             k=k+1
             internal = point_internal_to_surface(num_vertices,triangle, &
                  point_xyz,vertex_xyz)
             if(internal)then
                num_data = num_data+1
                if(num_data.le.num_data_estimate)then
                   data_xyz(:,num_data) = point_xyz
                else
                   num_data_estimate = num_data_estimate + 1000
                   allocate(data_temp(3,num_data-1))
                   data_temp = data_xyz ! copy to temporary array
                   deallocate(data_xyz) !deallocate initially allocated memory
                   allocate(data_xyz(3,num_data_estimate))
                   data_xyz(:,1:num_data-1) = data_temp(:,1:num_data-1)
                   deallocate(data_temp) !deallocate the temporary array
                   
                   write(*,'('' WARNING: number of data is'',I6, &
                        &''; increased array size'')') num_data
                   data_xyz(:,num_data) = point_xyz
                endif
             endif
             ! increment the location in x direction by 'spacing'
             point_xyz(1) = point_xyz(1) + spacing
          enddo
          ! increment the location in y direction by 'spacing'
          point_xyz(1) = min_bound(1) + 0.5_dp*spacing
          point_xyz(2) = point_xyz(2) + spacing
       enddo
       ! increment the location in z direction by 'spacing'
       point_xyz(1:2) = min_bound(1:2) + 0.5_dp*spacing
       point_xyz(3) = point_xyz(3) + spacing
    enddo
    
    write(*,'('' Made'',I7,'' data points inside surface elements'')') num_data
    
    if(allocated(data_weight)) deallocate(data_weight)
    allocate(data_weight(3,num_data))
    data_weight(:,1:num_data) = 1.0_dp
    
    deallocate(triangle)
    deallocate(vertex_xyz)
    deallocate(elem_list)
    
    call enter_exit(sub_name,2)
    
  end subroutine make_data_grid
  
!!!#############################################################################

  subroutine make_2d_vessel_from_1d(elem_list)
    !*make_2d_vessel_from_1d:* create a surface mesh that aligns with the
    ! centrelines of a 1D tree, and located at distance 'radius' from the centre.
    ! a template for a set of 5 nodes (that together define a bifurcation) is
    ! scaled, rotated, translated to align with the 1d mesh and its radii. 
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_MAKE_2D_VESSEL_FROM_1D" :: MAKE_2D_VESSEL_FROM_1D

    integer,intent(in) :: elem_list(:)
    ! Local variables
    integer,allocatable :: elem_node_map(:,:,:) ! stores the surface nodes in ring around 1d nodes
    integer,allocatable :: short_elements(:)    ! stores short surface elements for removing.
    !                                             these will be at short successive bifurcations
    integer,allocatable :: elem_ring(:,:),node_ring(:,:)
    integer :: template_cnct(2,8)               ! the node numbering for the templated elements
    integer :: template_vrsn_map(2,8)           ! versions of nodes for the templated elements
    integer :: template_vrsns(5)                ! # of versions of derivatives for 'template' bifurcation
    integer :: i,j,k,ne_sib,np_side1,np_side2,ne,ne_child,ne_count,ne_global,ne_new, &
         ne0,nj,nk,nmax,nn,np_crux,np_new,np_now,np_prev,np0,np1,np2,np_close(2), &
         num_short,nv,nvb,max_num_nodes_2d,start_npth
    real(dp) :: new_coords_derivs(4,10,3,5)     ! coordinates of translated and rotated template
    real(dp) :: ring_coords(3,5)                ! the coordinates of nodes in a current 'ring'
    real(dp),allocatable :: ring_distance(:)    ! distance of new 'ring' of nodes from 1d start node
    real(dp) :: template_coords(4,6,3,5)        ! coordinates of 5 nodes in 'template' bifurcation
    real(dp) :: angle,cruxdist,distance,length,line_direction(3),l_to_d,min_angle,point1(3),point2(3),point3(3), &
         radius,radius_weighted,Rx(3,3),Ry(3,3),smallest,smallest_angle,Txyz(3),u(3),v(3),vector_1(3), &
         vector_2(3),zz(3)

    real(dp) :: ln_23
    logical :: bifn_parent,entry_branch,major_child,minor_child,single_child
    logical :: merge_tree = .true.

    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'make_2d_vessel_from_1d'
    call enter_exit(sub_name,1)

!!! adjust the 1D tree so that the elements are not shorter than they are wide
    do ne = 1,num_elems
       ne0 = elem_cnct(-1,1,ne)
       ne_sib = elem_cnct(1,1,ne0)
       if(ne_sib.eq.ne) ne_sib = elem_cnct(1,2,ne0)
       if(ne_sib.ne.0.and.elem_cnct(1,0,ne).eq.1)then ! only for branches with single child branch
          l_to_d = elem_field(ne_length,ne)/elem_field(ne_radius,ne)
          if(l_to_d.le.2.0_dp)then
             np1 = elem_nodes(1,ne) ! node at Xi1 = 0
             np2 = elem_nodes(2,ne) ! node at Xi1 = 1
             u(:) = elem_direction(:,ne) ! unit vector for direction of minor branch
             v(:) = elem_direction(:,ne_sib) ! unit vector for direction of sibling branch
             zz(:) = unit_vector(0.5_dp * (u(:) + v(:)))
             smallest_angle = acos(scalar_product_3(zz,v))
             distance = elem_field(ne_radius,ne)/sin(smallest_angle)
             node_xyz(:,np2) = node_xyz(:,np1) + elem_direction(:,ne) * distance
             elem_field(ne_length,ne) = distance
             ne_child = elem_cnct(1,1,ne)
             np1 = elem_nodes(1,ne_child) ! node at Xi1 = 0
             np2 = elem_nodes(2,ne_child) ! node at Xi1 = 1
             elem_field(ne_length,ne_child) = &
                  distance_between_points(node_xyz(:,np1),node_xyz(:,np2))
         endif
       endif
    enddo

!!! allocate memory for the 2d geometry arrays

    max_num_nodes_2d = num_nodes * 5 * 2 ! maximum possible
    num_elems_2d = num_elems * 4 * 2 !maximum possible
    allocate(nodes_2d(max_num_nodes_2d))
    allocate(node_xyz_2d(4,10,16,max_num_nodes_2d)) ! maximum possible
    allocate(node_versn_2d(max_num_nodes_2d))
    allocate(elem_node_map(3,5,max_num_nodes_2d))
    allocate(elems_2d(num_elems_2d))
    allocate(elem_nodes_2d(4,num_elems_2d))
    allocate(elem_versn_2d(4,num_elems_2d))
    allocate(short_elements(4*num_elems_2d))
    allocate(elem_ring(4,num_elems))
    allocate(node_ring(4,num_elems))

    node_versn_2d = 1 ! default to 1 version for node np
    elem_versn_2d = 1 ! default to 1 version for node nn in element ne
    elem_node_map = 0
    elem_ring = 0
    node_ring = 0
    ne_new = 0  ! initialise the surface mesh element numbering
    np_new = 0  ! initialise the surface mesh node numbering

!!! set up a generic structure (in template_coords) that will be rotated, scaled, and placed 
    ! the following arrays define the template bifurcation
    template_vrsns = (/2,6,2,6,2/)
    template_vrsn_map = reshape((/1,2,3,1,1,3,2,1,1,5,4,2,2,4,5,1/),shape(template_vrsn_map))
    template_cnct = reshape((/5,2,2,3,3,4,4,5,1,2,2,5,5,4,4,1/),shape(template_cnct))
    call mesh_2d_from_1d_generic(template_coords)

    ne_global = elem_list(1) ! the global stem element number for the 2d mesh
    ne = get_local_elem_1d(ne_global)  ! the local element number
    ne_count = 1 ! counter for the 1d elements in the list

    do while (ne /= 0)

       np1 = elem_nodes(1,ne)   ! start node of current element 
       np2 = elem_nodes(2,ne)   ! end node of current element
       radius = elem_field(ne_radius,ne)
       length = elem_field(ne_length,ne)  ! length of current element

       ! check which type of branch
       bifn_parent  = .false.
       entry_branch = .false.
       major_child  = .false.
       minor_child  = .false.
       single_child = .false.

       ne0 = elem_cnct(-1,1,ne) ! parent
       if(ne.eq.1)then
          entry_branch = .true.
       else
          if(elem_cnct(1,0,ne0).eq.1)then
             single_child = .true.
          else
             ne_sib = elem_cnct(1,1,ne0)
             if(ne_sib.eq.ne) ne_sib = elem_cnct(1,2,ne0)
             if(radius.gt.elem_field(ne_radius,ne_sib))then
                major_child = .true.
             else
                if(abs(radius - elem_field(ne_radius,ne_sib)).le.loose_tol)then
                   ! check branch angles
                   if(angle_btwn_vectors(elem_direction(:,ne0),elem_direction(:,ne)) &
                        .gt.angle_btwn_vectors(elem_direction(:,ne0),elem_direction(:,ne_sib)))then
                      minor_child = .true.
                   else
                      major_child = .true.
                   endif
                else
                   minor_child = .true.
                endif
             endif
          endif
       endif
       
       ! calculate the rotation angle and matrices for template mesh about 1d element
       !call calc_rotat_trans_mats(ne,np1,Rx,Ry,Txyz)
       !angle = 0.0_dp ! 
       call mesh_rotate_about_axis(ne,angle,Rx,Ry,Txyz,template_coords)
       
       ! calculate new locations and derivatives for scaling, translation, rotation of template
       call mesh_rotate_vector_about_axis(4,template_vrsns,angle,0.0_dp,radius,radius, &
            Rx,Ry,Txyz,template_coords,new_coords_derivs)

       ! create 4 new nodes at start and record their node numbers, coordinates
       if(entry_branch.or.minor_child)then
          do i = 1,4
             np_new = np_new+1
             nodes_2d(np_new) = np_new ! global element number stored in local index
             if(entry_branch)then
                elem_node_map(1,i,np1) = np_new !record nodes in ring around np1
                node_xyz_2d(1:3,1,:,np_new) = new_coords_derivs(1:3,1,:,i) ! coordinates and derivatives 1,2,1_2
             elseif(minor_child)then
                elem_node_map(2,i,np1) = np_new !record nodes in ring around np1
                node_xyz_2d(1,1,:,np_new) = new_coords_derivs(1,1,:,i) & ! coordinates and derivatives 1,2,1_2
                     + elem_direction(:,ne) * elem_field(ne_radius,ne0)
                node_xyz_2d(2:3,1,:,np_new) = new_coords_derivs(2:3,1,:,i)
             endif
          enddo !i
       endif

       ! create 4 new nodes at end and record their node numbers, coordinates
       do i = 1,4
          np_new = np_new+1
          nodes_2d(np_new) = np_new ! global element number stored in local index
          elem_node_map(1,i,np2) = np_new !record nodes in ring around np2
          node_xyz_2d(1:3,1,:,np_new) = new_coords_derivs(1:3,1,:,i) ! coordinates and derivatives 1,2,1_2
          forall (k = 1:3) node_xyz_2d(1,1,k,np_new) = node_xyz_2d(1,1,k,np_new) &
               + elem_direction(k,ne) * length! coordinates and derivatives 1,2,1_2
          ne_new = ne_new + 1
          elems_2d(ne_new) = ne_new  

          if(minor_child)then
             np_prev = elem_node_map(2,i,np1)
             elem_nodes_2d(1,ne_new) = np_prev
             elem_nodes_2d(3,ne_new) = np_new
             elem_ring(i,ne) = ne_new
             if(i.lt.4)then
                elem_nodes_2d(2,ne_new) = np_prev + 1
                elem_nodes_2d(4,ne_new) = np_new + 1
             else
                elem_nodes_2d(2,ne_new) = elem_node_map(2,1,np1)
                elem_nodes_2d(4,ne_new) = np_new - 3
             endif
          else
             min_angle = 1.0e6_dp
             if(i.eq.1)then ! first time through, figure out connectivity
                do j = 1,4
                   np_prev = elem_node_map(1,j,np1)
                   line_direction = direction_point_to_point(node_xyz_2d(1,1,:,np_prev),node_xyz_2d(1,1,:,np_new))
                   angle = angle_btwn_vectors(line_direction,elem_direction(:,ne))
                   if(angle.lt.min_angle)then
                      min_angle = angle
                      start_npth = j
                   endif
                enddo
             endif
             np_prev = elem_node_map(1,start_npth,np1)
             elem_nodes_2d(1,ne_new) = np_prev
             elem_nodes_2d(3,ne_new) = np_new
             elem_ring(i,ne) = ne_new
             if(i.lt.4)then
                elem_nodes_2d(4,ne_new) = np_new + 1
             else
                elem_nodes_2d(4,ne_new) = np_new - 3
             endif
             if(start_npth.lt.4)then
                elem_nodes_2d(2,ne_new) = np_prev + 1
                start_npth = start_npth + 1
             else
                elem_nodes_2d(2,ne_new) = elem_node_map(1,1,np1)
                start_npth = 1
             endif
          endif
       enddo !i

       ne_count = ne_count+1
       if(ne_count.gt.num_elems)then
          ne = 0
       else
          if(ne_count .gt. count(elem_list.ne.0))then
             ne = 0
          else
             ne_global = elem_list(ne_count)
             ne = get_local_elem_1d(ne_global)
          endif
       endif
    enddo ! ne

    num_nodes_2d = np_new
    num_elems_2d = ne_new
    call element_connectivity_2d
    call line_segments_for_2d_mesh('arcl')

    if(merge_tree)then
!!! run through a second time, merging the minor branches
    ne_global = elem_list(1) ! the global stem element number for the 2d mesh
    ne = get_local_elem_1d(ne_global)  ! the local element number
    ne_count = 1 ! counter for the 1d elements in the list

    do while (ne /= 0)
       np1 = elem_nodes(1,ne)   ! start node of current element 
       np2 = elem_nodes(2,ne)   ! end node of current element
       radius = elem_field(ne_radius,ne)
       length = elem_field(ne_length,ne)  ! length of current element

       ! check which type of branch
       bifn_parent  = .false.
       entry_branch = .false.
       major_child  = .false.
       minor_child  = .false.
       single_child = .false.

       ne0 = elem_cnct(-1,1,ne) ! parent
       if(ne.eq.1)then
          entry_branch = .true.
       else
          if(elem_cnct(1,0,ne0).eq.1)then
             single_child = .true.
          else
             ne_sib = elem_cnct(1,1,ne0)
             if(ne_sib.eq.ne) ne_sib = elem_cnct(1,2,ne0)
             if(radius.gt.elem_field(ne_radius,ne_sib))then
                ! major branch has largest radius
                major_child = .true.
             else
                if(abs(radius - elem_field(ne_radius,ne_sib)).le.loose_tol)then
                   ! if radii the same, check the angles
                   if(abs(angle_btwn_vectors(elem_direction(:,ne0),elem_direction(:,ne))) &
                        .gt.abs(angle_btwn_vectors(elem_direction(:,ne0),elem_direction(:,ne_sib))))then
                      ! minor branch has largest angle
                      minor_child = .true.
                   else
                      ! major branch has smallest angle
                      major_child = .true.
                   endif
                else
                   ! minor child has smallest radius
                   minor_child = .true.
                endif
             endif
          endif
       endif
       if(minor_child)then
          if(abs(angle_btwn_vectors(elem_direction(:,ne0), &
               elem_direction(:,ne))/pi*180.0_dp).gt.60.0_dp.or. &
               (elem_field(ne_radius,ne)/elem_field(ne_radius,ne0).lt.0.75_dp))then
             ! if the minor child is more than 75% radius of major, and the angles are similar
             !if((elem_field(ne_radius,ne)/elem_field(ne_radius,ne_sib).lt.0.75_dp).or. &
             !     (abs(angle_btwn_vectors(elem_direction(:,ne0),elem_direction(:,ne)) - &
             !     angle_btwn_vectors(elem_direction(:,ne0),elem_direction(:,ne_sib))).gt. &
             !     30.0_dp/180.0_dp*pi))then
             call create_minor_branch(elem_ring,ne,max_num_nodes_2d)
          else
             call create_bifurcation(elem_ring,ne)
          endif
       else if(single_child)then
!          write(*,*) 'single child at ne=',ne,' is ',elem_cnct(1,1,ne)
!          write(*,*) 'ring nodes',elem_ring(1:4,elem_cnct(1,1,ne))
       endif
       
       ne_count = ne_count+1
       if(ne_count.gt.num_elems)then
          ne = 0
       else
          if(ne_count .gt. count(elem_list.ne.0))then
             ne = 0
          else
             ne_global = elem_list(ne_count)
             ne = get_local_elem_1d(ne_global)
          endif
       endif
    enddo ! ne

    call element_connectivity_2d
    call line_segments_for_2d_mesh('arcl')

    endif
    
    deallocate(elem_node_map)
    deallocate(short_elements)
    deallocate(elem_ring)
    deallocate(node_ring)

  end subroutine make_2d_vessel_from_1d

!!!#############################################################################
  
  subroutine create_bifurcation(elem_ring,ne)

    use mesh_utilities,only: inlist
    
    integer :: elem_ring(:,:),ne
    
    integer :: elem_minor(2),i,j,ne0,ne_ring,ne_closest1,ne_closest2,ne_sib,n_elem, &
         n_node,nm_list,np,np1,np2,np_centre,np_close,np_crux,np_list(3)
    real(dp) :: centres(4,3,3),distance,dist_closest1,dist_closest2,mid_point(3), &
         radius,smallest_angle,u(3),v(3),w(3),zz(3)

!!! get the parent 1d element number
    ne0 = elem_cnct(-1,1,ne)
    
!!! get the 'sibling' element number
    ne_sib = elem_cnct(1,1,ne0)
    if(ne_sib.eq.ne) ne_sib = elem_cnct(1,2,ne0)

!!! DETERMINE THE CONNECTIVITY OF MAJOR AND MINOR ELEMENT RINGS
    ! calculate the centre of each element in the rings around the parent,
    ! minor, major branches. Will be used to determine which elements are
    ! to be edited to create the bifurcation.
    centres = 0.0_dp
    do n_elem = 1,4 ! for each element in the ring
       do n_node = 1,4 ! for each node in the element
          ! centres for the parent ring
          ne_ring = elem_ring(n_elem,ne0)
          np = elem_nodes_2d(n_node,ne_ring) ! parent
          centres(n_elem,1,:) = centres(n_elem,1,:) + 0.25_dp*node_xyz_2d(1,1,:,np)
          ! centres for minor branch ring
          ne_ring = elem_ring(n_elem,ne)
          np = elem_nodes_2d(n_node,ne_ring) ! minor branch
          centres(n_elem,2,:) = centres(n_elem,2,:) + 0.25_dp*node_xyz_2d(1,1,:,np)
          ! centres for major branch ring
          ne_ring = elem_ring(n_elem,ne_sib)
          np = elem_nodes_2d(n_node,ne_ring) ! major branch
          centres(n_elem,3,:) = centres(n_elem,3,:) + 0.25_dp*node_xyz_2d(1,1,:,np)
       enddo !n_node
    enddo !n_elem
    
    ! Calculate the mid-point coordinates of the minor branch centreline.
    ! Will be used to calculate which of the major branch ring elements
    ! are closest to the minor branch.
    np1 = elem_nodes(1,ne) ! node at Xi=0
    np2 = elem_nodes(2,ne) ! node at Xi=1
    mid_point(:) = 0.5_dp * (node_xyz(:,np1) + node_xyz(:,np2))

    ! determine which two major branch elements are closest to the minor branch mid-point
    ne_closest1 = 0
    ne_closest2 = 0
    dist_closest1 = 1.0e6_dp
    dist_closest2 = 1.0e6_dp
    ! closest element
    do n_elem = 1,4
       u(:) = centres(n_elem,3,:) ! coordinates for major branch ring
       if(distance_between_points(u,mid_point).lt.dist_closest1)then
          dist_closest1 = distance_between_points(u,mid_point)
          ne_closest1 = n_elem
       endif
    enddo ! n_elem
    ! second closest element
    do n_elem = 1,4
       if(n_elem.ne.ne_closest1)then
          u(:) = centres(n_elem,3,:) ! coordinates for major branch ring
          if(distance_between_points(u,mid_point).lt.dist_closest2)then
             dist_closest2 = distance_between_points(u,mid_point)
             ne_closest2 = n_elem
          endif
       endif
    enddo ! n_elem
    
    ! order the closest elements for convenience
    if(elem_ring(ne_closest1,ne_sib).gt.elem_ring(ne_closest2,ne_sib).and.(ne_closest1.ne.4) &
         .or. (ne_closest1.eq.1.and.ne_closest2.eq.4) &
         .or. (ne_closest1.eq.4.and.ne_closest2.eq.3)) &
         call swap_integers(ne_closest1,ne_closest2)

    ! get the node number that is shared by the two closest elements
    np_centre = elem_nodes_2d(2,elem_ring(ne_closest1,ne_sib))
    v(:) = node_xyz_2d(1,1,:,np_centre)
    
    ! find the closest node at Xi=0 from the minor branch ring. This will be
    ! replaced by np_centre in the ring of elements for the minor branch
    dist_closest1 = 1.0e6_dp
    do n_elem = 1,4
       ne_ring = elem_ring(n_elem,ne) ! for each element in the ring
       np = elem_nodes_2d(1,ne_ring)  ! for nodes at Xi=0
       u(:) = node_xyz_2d(1,1,:,np)
       if(distance_between_points(u,v).lt.dist_closest1)then
          dist_closest1 = distance_between_points(u,v)
          np_close = np
       endif
    enddo ! n_node

    ! replace np_close with np_centre in the minor branch ring
    ! np_list stores the nodes that have been 'joined up'
    np_list = 0
    nm_list = 0
    do n_elem = 1,4
       ne_ring = elem_ring(n_elem,ne) ! for each element in the ring
       if(elem_nodes_2d(1,ne_ring).eq.np_close)then
          np_list(3) = np_centre
          np_list(2) = elem_nodes_2d(2,elem_ring(ne_closest2,ne_sib))
          elem_nodes_2d(1,ne_ring) = np_centre
          elem_nodes_2d(2,ne_ring) = np_list(2)
          if(n_elem.ne.4)then !merge minor elements into a front node
             elem_nodes_2d(1,ne_ring+1) = np_list(2)
          else
             elem_nodes_2d(1,ne_ring-3) = np_list(2)
          endif
       else if(elem_nodes_2d(2,ne_ring).eq.np_close)then
          np_list(1) = elem_nodes_2d(1,elem_ring(ne_closest1,ne_sib))
          elem_nodes_2d(2,ne_ring) = np_centre
          elem_nodes_2d(1,ne_ring) = np_list(1)
          if(n_elem.ne.1)then !merge minor elements into a front node
             elem_nodes_2d(2,ne_ring-1) = np_list(1)
          else
             elem_nodes_2d(2,ne_ring+3) = np_list(1)
          endif
       else
          nm_list = nm_list + 1
          elem_minor(nm_list) = elem_ring(n_elem,ne)
       endif
    enddo ! n_elem

    ! get the node number for the node that is still 'hanging' (unjoined)
    do n_elem = 1,4
       ne_ring = elem_ring(n_elem,ne) ! for each element in the ring
       np = elem_nodes_2d(1,ne_ring)
       if(.not.inlist(np,np_list)) np_crux = np
    enddo

    ! replace np_centre in the sibling ring of elements with np_crux
    elem_nodes_2d(2,elem_ring(ne_closest1,ne_sib)) = np_crux
    elem_nodes_2d(1,elem_ring(ne_closest2,ne_sib)) = np_crux

!!! SET THE CRUX LOCATION AND DERIVATIVES
    ! set the crux node location to be in the direction of the average child branch 
    ! direction (ne and ne_sib) at a distance that is consistent with the radius of
    ! the major branch. That is:
    !     xyz_crux = xyz(at start of ne) + average direction (of ne and ne_sib)*distance
    !     where distance = radius(ne_sib)/angle(between ne_sib and average direction)
    !     and angle(between ne_sib and average direction) = acos(dot product(ne_sib, average))

    ! direction mid-way between the two child branches
    u(:) = elem_direction(:,ne) ! unit vector for direction of minor branch
    v(:) = elem_direction(:,ne_sib) ! unit vector for direction of sibling branch
    zz(:) = unit_vector(0.5_dp * (u(:) + v(:)))
    smallest_angle = acos(scalar_product_3(zz,v)) ! angle should be smallest for major branch
    distance = elem_field(ne_radius,ne_sib)/sin(smallest_angle)
    node_xyz_2d(1,1,:,np_crux) = node_xyz(:,elem_nodes(1,ne)) + zz(:) * distance

    ! increase the number of versions at the crux node
    node_versn_2d(np_crux) = 2
    node_xyz_2d(:,2,:,np_crux) = node_xyz_2d(:,1,:,np_crux) ! initialise to the 1st version values
    elem_versn_2d(2,elem_ring(ne_closest1,ne_sib)) = 2 ! minor branch uses version 1, major version 2
    elem_versn_2d(1,elem_ring(ne_closest2,ne_sib)) = 2

!!! SET THE FRONT AND BACK NODE LOCATION AND DERIVATIVES
    ! adjust the 'front' and 'back' bifurcation nodes so that the radial dimension
    ! is the weighted average of the child radii and parent radius
    radius = 0.5_dp * (elem_field(ne_radius,ne)+elem_field(ne_radius,ne_sib))
    radius = 0.5_dp * (elem_field(ne_radius,ne0) + radius)
    u(:) = node_xyz_2d(1,1,:,np_list(1))
    v(:) = node_xyz_2d(1,1,:,np_list(2))
    node_xyz_2d(1,1,:,np_list(1)) = 0.5_dp*(u(:)+v(:)) - radius * direction_point_to_point(u,v)
    node_xyz_2d(1,1,:,np_list(2)) = 0.5_dp*(u(:)+v(:)) + radius * direction_point_to_point(u,v)

    ! use the direction from front to back to define crux derivative
    node_xyz_2d(2,1,:,np_crux) = -direction_point_to_point(u,v)
    node_xyz_2d(2,2,:,np_crux) = -node_xyz_2d(2,1,:,np_crux)

    ! increase the number of versions at the front and back nodes (for different derivatives)
    do i = 1,2
       node_versn_2d(np_list(i)) = 3 
       forall(j=2:3) node_xyz_2d(1,j,:,np_list(i)) = node_xyz_2d(1,1,:,np_list(i))
       node_xyz_2d(2,2,:,np_list(i)) = node_xyz_2d(3,1,:,np_list(i)) ! 1st derivative
       node_xyz_2d(2,3,:,np_list(i)) = -node_xyz_2d(2,2,:,np_list(i))
       node_xyz_2d(3,2,:,np_list(i)) = 0.0_dp ! 2nd derivative
       node_xyz_2d(3,3,:,np_list(i)) = node_xyz_2d(3,2,:,np_list(i))
    enddo ! i
        
    elem_versn_2d(1,elem_ring(ne_closest1,ne_sib)) = 2 ! use 2nd version for 1st node in ne_closest1
    elem_versn_2d(2,elem_ring(ne_closest2,ne_sib)) = 3
    if(elem_nodes_2d(2,elem_minor(1)).eq.elem_nodes_2d(1,elem_ring(ne_closest1,ne_sib)))then
       elem_versn_2d(2,elem_minor(1)) = 3
       elem_versn_2d(1,elem_minor(2)) = 2
    else
       elem_versn_2d(1,elem_minor(1)) = 2
       elem_versn_2d(2,elem_minor(2)) = 3
    endif

!!! ADJUST THE 'SIDE' NODE LOCATIONS TO PROVIDE MORE CURVATURE INTO THE BIFURCATION
    np1 = elem_nodes_2d(2,elem_cnct_2d(-2,1,elem_ring(ne_closest1,ne_sib))) ! the Xi2=0 node in immediate distal element
    np2 = np_list(3)
    if(np1.ne.0) &
       node_xyz_2d(1,1,:,np2) = 0.8_dp * node_xyz_2d(1,1,:,np2) + 0.2_dp * node_xyz_2d(1,1,:,np1)
    np2 = elem_nodes_2d(2,elem_cnct_2d(1,1,elem_ring(ne_closest2,ne_sib))) ! the Xi1=1 node in immediate clockwise element
    np1 = elem_nodes_2d(2,elem_cnct_2d(-2,1,elem_cnct_2d(1,1,elem_ring(ne_closest2,ne_sib))))
    if(np1.ne.0.and.np2.ne.0) &
       node_xyz_2d(1,1,:,np2) = 0.8_dp * node_xyz_2d(1,1,:,np2) + 0.2_dp * node_xyz_2d(1,1,:,np1)
    
  end subroutine create_bifurcation
  
!!!#############################################################################

  subroutine swap_integers(n_1,n_2)

    integer :: n_1,n_2
    integer :: n_temp

    n_temp = n_1
    n_1 = n_2
    n_2 = n_temp

  end subroutine swap_integers
    
!!!#############################################################################
  
  subroutine create_minor_branch(elem_ring,ne,max_num_nodes_2d)

    integer :: elem_ring(:,:),ne,max_num_nodes_2d
    
    integer :: i,j,ne0,n_elem,ne_closest,ne_closest_parent,nelem_match(4,2),ne_new, &
         ne_minor,ne_next,ne_parent,ne_ring,np1,np2,np_centre,np_closest, &
         np_minor,np_new,np_ring,ntri,num_triangles,num_vertices,surface_elems(1)
    integer,allocatable :: triangle(:,:)
    real(dp) :: area,area_triangle,centre(3),denominator,dist_closest,distance, &
         gamma,keep_point(3),length,line_direction(3), &
         norm_v(3),point(3),P1(3),P2(3),P3(3),u(3),v(3)
    real(dp),allocatable :: vertex_xyz(:,:)
    character(len=MAX_FILENAME_LEN) :: EXNODEFILE
    character(len=MAX_STRING_LEN) :: groupname
    logical :: found

!!! get the parent 1d element number
    ne0 = elem_cnct(-1,1,ne)
    
!!! FIND CLOSEST NODE ON PARENT-MAJOR JUNCTION TO CENTRE OF MINOR BRANCH ELEMENT RING
    np1 = elem_nodes(1,ne) ! node at Xi1=0
    u(:) = node_xyz(:,np1) + elem_direction(:,ne) * elem_field(ne_radius,ne0)
    dist_closest = 1.0e6_dp
    do n_elem = 1,4
       ne_ring = elem_ring(n_elem,ne0) ! 2d element number for the n_elem'th ring around parent
       np_ring = elem_nodes_2d(3,ne_ring)
       v(:) = node_xyz_2d(1,1,:,np_ring)
       distance = distance_between_points(u,v)
       if(distance.lt.dist_closest)then
          dist_closest = distance
          ne_closest = ne_ring
          np_closest = np_ring
       endif
    enddo
    ne_closest_parent = ne_closest
    np_centre = np_closest

!!! FIND CLOSEST NODE ON MINOR RING TO THE 'CLOSEST' PARENT ELEMENT, AND SET UP
!!! VERSIONS OF DERIVATIVES FOR THE MINOR RING NODES
    u(:) = node_xyz_2d(1,1,:,elem_nodes_2d(4,ne_closest_parent)) 
    dist_closest = 1.0e6_dp
    do n_elem = 1,4
       ne_ring = elem_ring(n_elem,ne) ! minor branch elements
       np_ring = elem_nodes_2d(2,ne_ring)
       ! set up the versions of derivatives, to be used later
       node_versn_2d(np_ring) = 3 ! the number of versions
       forall (j=2:3) node_xyz_2d(1,j,:,np_ring) = node_xyz_2d(1,1,:,np_ring)
       v(:) = node_xyz_2d(1,1,:,np_ring)
       distance = distance_between_points(u,v)
       if(distance.lt.dist_closest)then
          dist_closest = distance
          ne_closest = ne_ring
          np_closest = np_ring
       endif
    enddo
    
    ! record the elements on the parent-major junction that align with the minor elements
    nelem_match(1,1) = ne_closest_parent
    nelem_match(1,2) = ne_closest
    nelem_match(2,1) = elem_cnct_2d(-1,1,nelem_match(1,1)) ! adjacent in -Xi1 direction
    nelem_match(2,2) = elem_cnct_2d(-1,1,nelem_match(1,2)) ! adjacent in -Xi2 direction
    nelem_match(3,1) = elem_cnct_2d(2,1,nelem_match(2,1))  ! adjacent in +Xi2 direction
    nelem_match(3,2) = elem_cnct_2d(-1,1,nelem_match(2,2)) ! adjacent in -Xi2 direction
    nelem_match(4,1) = elem_cnct_2d(1,1,nelem_match(3,1))  ! adjacent in +Xi1 direction
    nelem_match(4,2) = elem_cnct_2d(-1,1,nelem_match(3,2)) ! adjacent in -Xi2 direction
    
!!! SHIFT THE MINOR BRANCH NODES TO BE ON THE PARENT/MAJOR SURFACE
    do i = 1,4
       np1 = elem_nodes_2d(1,nelem_match(i,2)) ! minor branch node
       np2 = elem_nodes_2d(3,nelem_match(i,2)) ! adjacent node in +Xi2
       u(:) = node_xyz_2d(1,1,:,np1) ! location of minor branch node
       line_direction = direction_point_to_point(node_xyz_2d(1,1,:,np2),node_xyz_2d(1,1,:,np1))
       found = .false.
       do j = 1,4
          ne_parent = nelem_match(j,1) ! parent-major element
          surface_elems(1) = ne_parent
          call triangles_from_surface(num_triangles,num_vertices,surface_elems, &
               triangle,vertex_xyz)
          do ntri = 1,num_triangles
             ! get the normal to each triangle
             P1(1:3) = vertex_xyz(1:3,triangle(1,ntri))
             P2(1:3) = vertex_xyz(1:3,triangle(2,ntri))
             P3(1:3) = vertex_xyz(1:3,triangle(3,ntri))
             norm_v = unit_norm_to_three_points(P1,P2,P3) ! unit normal to triangle plane
             denominator = scalar_product_3(norm_v,line_direction) ! denominator is zero if line parallel to plane
             if(abs(denominator).gt.loose_tol)then ! not normal to the plane of the triangle
                ! calculate the distance of the surface point from point_xyz
                length = scalar_product_3(norm_v,P1-u)/denominator ! distance from np1 to the triangle plane
                point = u + line_direction * length ! point on the triangle plane that the line passes through
                area_triangle = area_between_two_vectors(P1-P2,P1-P3) ! area of triangle
                ! calculate summed area of triangles that join the triangle vertices to the
                ! point on the triangle plane. if the area is the same as the area of the triangle,
                ! then the point is within the triangle
                area = area_between_two_vectors(P1-point,P2-point)+ & 
                     area_between_two_vectors(P1-point,P3-point)+area_between_two_vectors(P2-point,P3-point)
                if(abs(area_triangle-area).lt.loose_tol)then
                   keep_point = point
                   found = .true.
                endif
             endif
          enddo ! ntri
       enddo !j
       node_xyz_2d(1,1,:,np1) = keep_point(:)
    enddo !i
    deallocate(triangle)
    deallocate(vertex_xyz)

!!! MOVE THE ADJACENT PARENT/MAJOR JUNCTION NODE TO BE CLOES TO ORTHOGONAL TO MINOR CENTRELINE
    do i = 3,4
       ne_parent = nelem_match(i,1)
       ne_minor = nelem_match(i,2)
       if(i.eq.4)then
          np1 = elem_nodes_2d(2,ne_parent) ! at Xi(1,0) on parent
          np2 = elem_nodes_2d(4,ne_parent) ! at Xi(1,1) on parent
          np_minor = elem_nodes_2d(1,ne_minor)
       else
          np1 = elem_nodes_2d(1,ne_parent) ! at Xi(1,0) on parent
          np2 = elem_nodes_2d(3,ne_parent) ! at Xi(1,1) on parent
          np_minor = elem_nodes_2d(2,ne_minor)
       endif
       u(:) = node_xyz_2d(1,1,:,np_minor) - node_xyz_2d(1,1,:,np1)
       v(:) = unit_vector(node_xyz_2d(1,1,:,np2) - node_xyz_2d(1,1,:,np1))
       node_xyz_2d(1,1,:,np1) = node_xyz_2d(1,1,:,np1) + v(:) * scalar_product_3(u,v)
    enddo
    
!!! MAKE A NEW NODE ON EACH RESPECTIVE PARENT/MAJOR ELEMENT
    do i = 1,4
       num_nodes_2d = num_nodes_2d + 1
       np_new = num_nodes_2d
       nodes_2d(np_new) = np_new
       ne_parent = nelem_match(i,1)
       ne_minor = nelem_match(i,2)
       if(i.eq.1.or.i.eq.4)then
          np1 = elem_nodes_2d(2,ne_parent) ! at Xi(1,0) on parent
          np2 = elem_nodes_2d(4,ne_parent) ! at Xi(1,1) on parent
       else
          np1 = elem_nodes_2d(1,ne_parent) ! at Xi(1,0) on parent
          np2 = elem_nodes_2d(3,ne_parent) ! at Xi(1,1) on parent
       endif
       
       if(i.eq.1.or.i.eq.3)then
          np_minor = elem_nodes_2d(1,ne_minor)
       else if(i.eq.2.or.i.eq.4)then
          np_minor = elem_nodes_2d(2,ne_minor)
       endif

       ! position the new node such that the minor branch node is orthogonal to the
       ! line joining np1 and np2
       u(:) = node_xyz_2d(1,1,:,np_minor) - node_xyz_2d(1,1,:,np1)
       v(:) = unit_vector(node_xyz_2d(1,1,:,np2) - node_xyz_2d(1,1,:,np1))
       node_xyz_2d(1,1,:,np_new) = node_xyz_2d(1,1,:,np1) + v(:) * scalar_product_3(u,v)

       ! check that the new node position isn't too close to the adjusted parent/major
       ! junction node. shift node np_new if it is too close
       gamma = distance_between_points(node_xyz_2d(1,1,:,np_new),node_xyz_2d(1,1,:,np2)) &
            /distance_between_points(node_xyz_2d(1,1,:,np1),node_xyz_2d(1,1,:,np2))
       if(i.le.2.and.gamma.lt.0.4_dp) then
          node_xyz_2d(1,1,:,np_new) = 0.4_dp * node_xyz_2d(1,1,:,np1) + 0.6_dp * node_xyz_2d(1,1,:,np2)
       else if(i.gt.2.and.gamma.gt.0.6_dp) then
          !node_xyz_2d(1,1,:,np_new) = 0.25_dp * node_xyz_2d(1,1,:,np1) + 0.75_dp * node_xyz_2d(1,1,:,np2)
       endif
!       if(i.le.2.and.gamma.lt.0.25_dp) then
!          node_xyz_2d(1,1,:,np_new) = 0.25_dp * node_xyz_2d(1,1,:,np1) + 0.75_dp * node_xyz_2d(1,1,:,np2)
!       else if(i.gt.2.and.gamma.gt.0.75_dp) then
!          !node_xyz_2d(1,1,:,np_new) = 0.25_dp * node_xyz_2d(1,1,:,np1) + 0.75_dp * node_xyz_2d(1,1,:,np2)
!       endif

       ! increase the number of versions for the minor branch element nodes, set deriv values
       ne_minor = nelem_match(i,2)
       np_minor = elem_nodes_2d(2,ne_minor)  ! at Xi(1,0) on minor
       if(i.eq.1.or.i.eq.3)then
          node_versn_2d(np_minor) = 2 ! the number of versions
          node_xyz_2d(1,2,:,np_minor) = node_xyz_2d(1,1,:,np_minor)
          node_xyz_2d(2,2,:,np_minor) = 0.4_dp * node_xyz_2d(2,1,:,np_centre) &
               - 0.6_dp * node_xyz_2d(3,1,:,np_minor)
!          node_xyz_2d(2,2,:,np_minor) = 0.25_dp * node_xyz_2d(2,1,:,np_centre) &
!               - 0.75_dp * node_xyz_2d(3,1,:,np_minor)
          if(i.eq.1)then
             node_xyz_2d(3,2,:,np_minor) = node_xyz_2d(2,1,:,np_minor) * 0.5_dp
          elseif(i.eq.3)then
             node_xyz_2d(3,2,:,np_minor) = -node_xyz_2d(2,1,:,np_minor) * 0.5_dp
          endif
       else ! 2 or 4
          node_versn_2d(np_minor) = 4 ! the number of versions
          forall (j = 2:4) node_xyz_2d(1,j,:,np_minor) = node_xyz_2d(1,1,:,np_minor)
          node_xyz_2d(2,2,:,np_minor) = node_xyz_2d(2,1,:,np_centre)
          node_xyz_2d(3,2,:,np_minor) = node_xyz_2d(3,1,:,np_centre)
          node_xyz_2d(2,3,:,np_minor) = node_xyz_2d(2,2,:,np_minor)
          if(i.eq.2)then
             node_xyz_2d(3,3,:,np_minor) = node_xyz_2d(2,1,:,np_minor) * 0.5_dp
             node_xyz_2d(3,4,:,np_minor) = -node_xyz_2d(2,1,:,np_minor) * 0.5_dp
          else if(i.eq.4)then
             node_xyz_2d(3,3,:,np_minor) = -node_xyz_2d(2,1,:,np_minor) * 0.5_dp
             node_xyz_2d(3,4,:,np_minor) = node_xyz_2d(2,1,:,np_minor) * 0.5_dp
          endif
          node_xyz_2d(2,4,:,np_minor) = node_xyz_2d(2,2,:,np_minor) !**
       endif
       ! make a new element that joins to the new parent/major node
       num_elems_2d = num_elems_2d + 1
       ne_new = num_elems_2d
       elems_2d(ne_new) = ne_new
       select case(i)
       case(1)
          elem_nodes_2d(1,ne_new) = elem_nodes_2d(1,ne_minor) ! 1,4
          elem_versn_2d(1,ne_new) = 3 !1
          elem_nodes_2d(2,ne_new) = num_nodes_2d !1
          elem_nodes_2d(3,ne_new) = np_minor !1,4
          elem_versn_2d(3,ne_new) = 2 !1
          elem_nodes_2d(4,ne_new) = elem_nodes_2d(4,ne_parent) !1
          ! replace the nodes at Xi2=1 in the original parent/major element
          elem_nodes_2d(3,ne_parent) = elem_nodes_2d(1,ne_minor) !1
          elem_nodes_2d(4,ne_parent) = num_nodes_2d !1
          elem_versn_2d(3,ne_parent) = 2 !1
       case(2)
          elem_nodes_2d(1,ne_new) = num_nodes_2d !2
          elem_nodes_2d(2,ne_new) = np_minor !2,3
          elem_versn_2d(2,ne_new) = 4 !2
          elem_nodes_2d(3,ne_new) = elem_nodes_2d(3,ne_parent) !2
          elem_nodes_2d(4,ne_new) = elem_nodes_2d(1,ne_minor) !2,3
          elem_versn_2d(4,ne_new) = 2 !2
          ! replace the nodes at Xi2=1 in the original parent/major element
          elem_nodes_2d(3,ne_parent) = num_nodes_2d !2
          elem_nodes_2d(4,ne_parent) = np_minor !2
          elem_versn_2d(4,ne_parent) = 2 !2
       case(3)
          elem_nodes_2d(1,ne_new) = elem_nodes_2d(1,ne_parent) !3
          elem_nodes_2d(2,ne_new) = np_minor !2,3
          elem_versn_2d(2,ne_new) = 2 !3
          elem_nodes_2d(3,ne_new) = num_nodes_2d !3
          elem_nodes_2d(4,ne_new) = elem_nodes_2d(1,ne_minor) !2,3
          elem_versn_2d(4,ne_new) = 3 !3
          ! replace the nodes at Xi2=0 in the original parent/major element
          elem_nodes_2d(1,ne_parent) = num_nodes_2d !3
          elem_nodes_2d(2,ne_parent) = elem_nodes_2d(1,ne_minor) !3
          elem_versn_2d(2,ne_parent) = 2 !3
       case(4)
          elem_nodes_2d(1,ne_new) = elem_nodes_2d(1,ne_minor) !1,4
          elem_versn_2d(1,ne_new) = 2 !4
          elem_nodes_2d(2,ne_new) = elem_nodes_2d(2,ne_parent) !4
          elem_nodes_2d(3,ne_new) = np_minor !1,4
          elem_versn_2d(3,ne_new) = 4 !4
          elem_nodes_2d(4,ne_new) = num_nodes_2d !4
          ! replace the nodes at Xi2=0 in the original parent/major element
          elem_nodes_2d(1,ne_parent) = np_minor !4
          elem_versn_2d(1,ne_parent) = 2 !4
          elem_nodes_2d(2,ne_parent) = num_nodes_2d !4
       end select
    enddo

  end subroutine create_minor_branch
  
!!!#############################################################################
  
  subroutine make_2d_vessel_from_1d0(elem_list)
    !*make_2d_vessel_from_1d:* create a surface mesh that aligns with the
    ! centrelines of a 1D tree, and located at distance 'radius' from the centre.
    ! a template for a set of 5 nodes (that together define a bifurcation) is
    ! scaled, rotated, translated to align with the 1d mesh and its radii. 
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_MAKE_2D_VESSEL_FROM_1D" :: MAKE_2D_VESSEL_FROM_1D

    integer,intent(in) :: elem_list(:)
    ! Local variables
    integer,allocatable :: elem_node_map(:,:,:) ! stores the surface nodes in ring around 1d nodes
    integer,allocatable :: short_elements(:)    ! stores short surface elements for removing.
    !                                             these will be at short successive bifurcations
    integer :: template_cnct(2,8)               ! the node numbering for the templated elements
    integer :: template_vrsn_map(2,8)           ! versions of nodes for the templated elements
    integer :: template_vrsns(5)                ! # of versions of derivatives for 'template' bifurcation
    integer :: i,j,k,np_side1,np_side2,ne,ne_child,ne_count,ne_global,ne_new, &
         ne0,nj,nk,nmax,nn,np_crux,np_new,np_now,np0,np1,np2,np_close(2), &
         num_short,nv,nvb
    real(dp) :: new_coords_derivs(4,10,3,5)     ! coordinates of translated and rotated template
    real(dp) :: ring_coords(3,5)                ! the coordinates of nodes in a current 'ring'
    real(dp),allocatable :: ring_distance(:)    ! distance of new 'ring' of nodes from 1d start node
    real(dp) :: template_coords(4,6,3,5)        ! coordinates of 5 nodes in 'template' bifurcation
    real(dp) :: angle,cruxdist,distance,length,point1(3),point2(3),point3(3), &
         radius,radius_weighted,Rx(3,3),Ry(3,3),smallest,Txyz(3),vector_1(3), &
         vector_2(3)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'make_2d_vessel_from_1d'
    call enter_exit(sub_name,1)

!!! allocate memory for the 2d geometry arrays

    num_nodes_2d = num_nodes * 5 * 2 ! maximum possible
    num_elems_2d = num_elems * 4 * 2 !maximum possible
    allocate(nodes_2d(num_nodes_2d))
    allocate(node_xyz_2d(4,10,16,num_nodes_2d)) ! maximum possible
    allocate(node_versn_2d(num_nodes_2d))
    allocate(elem_node_map(2,5,num_nodes_2d))
    allocate(ring_distance(num_nodes_2d))
    allocate(elems_2d(num_elems_2d))
    allocate(elem_nodes_2d(4,num_elems_2d))
    allocate(elem_versn_2d(4,num_elems_2d))
    allocate(short_elements(4*num_elems_2d))

    node_versn_2d = 1 ! default to 1 version for node np
    elem_versn_2d = 1 ! default to 1 version for node nn in element ne
    ring_distance = 0.0_dp
    short_elements = 0
    num_short = 0
    ne_new = 0  ! initialise the surface mesh element numbering
    np_new = 0  ! initialise the surface mesh node numbering

!!! set up a generic structure (in template_coords) that will be rotated, scaled, and placed 
    ! the following arrays define the template bifurcation
    template_vrsns = (/2,6,2,6,2/)
    template_vrsn_map = reshape((/1,2,3,1,1,3,2,1,1,5,4,2,2,4,5,1/),shape(template_vrsn_map))
    template_cnct = reshape((/5,2,2,3,3,4,4,5,1,2,2,5,5,4,4,1/),shape(template_cnct))
    call mesh_2d_from_1d_generic(template_coords)

!!! Place a first 'ring' of four nodes at the entrance to the model ---

    ne_global = elem_list(1) ! the global stem element number for the 2d mesh
    ne = get_local_elem_1d(ne_global)  ! the local element number
    np1 = elem_nodes(1,ne)  ! start node number of 1d element
    np2 = elem_nodes(2,ne)  ! end node number of 1d element
    radius = elem_field(ne_radius,ne)   ! radius for the 1D element
    ! calculate the rotation angle and matrices for template mesh about 1d element
    call mesh_rotate_about_axis(ne,angle,Rx,Ry,Txyz,template_coords)
    ! calculate new locations and derivatives for scaling, translation, rotation of template
    call mesh_rotate_vector_about_axis(4,template_vrsns,angle,0.0_dp,radius,radius, &
         Rx,Ry,Txyz,template_coords,new_coords_derivs)
    ! create 4 new nodes and record their node numbers, coordinates
    do i = 1,4
       np_new = np_new+1
       nodes_2d(np_new) = np_new ! global element number stored in local index
       elem_node_map(1,i,np1) = np_new !record nodes in ring around np1
       node_xyz_2d(:,1,:,np_new) = new_coords_derivs(:,1,:,i) ! coordinates and derivatives 1,2,1_2
    enddo !i
    elem_node_map(1,5,np1) = elem_node_map(1,1,np1) !dummy (5th) node to help with mapping
    
!!! work through each of the elements in the given list, creating a ring of new
!!! nodes and joining them to the previous ring of nodes. The assumption is that
!!! the elements comprise a continuous tree

    ne_count = 1 ! counter for the 1d elements in the list
    do while (ne /= 0)
       
       ne0 = elem_cnct(-1,1,ne) ! parent 1d element
       np1 = elem_nodes(1,ne)   ! start node of current element 
       np2 = elem_nodes(2,ne)   ! end node of current element
       
       radius = elem_field(ne_radius,ne)
       length = elem_field(ne_length,ne)  ! length of current element
       if(.not.bifurcation_element(ne).and.bifurcation_element(ne0)) &
            length = max(length,ring_distance(ne)*1.1_dp)

       ! check for the creation of 'short' surface elements that will later be removed
       if(length.le.radius)then
          forall (k = 1:4) short_elements(num_short + k) = ne_new + k
          num_short = num_short + 4
       endif
       
       ! calculate the rotation angle, translation (Txyz), and rotation matrices (Rx, Ry)
       call mesh_rotate_about_axis(ne,angle,Rx,Ry,Txyz,template_coords)

       nvb = 1 !version of node mapping to use
       if(.not.stem_element(ne))then     ! not the stem element (i.e. has a parent!)
          ! find out whether this is 1st or 2nd child of the parent. this indicates
          ! which 'ring' of parent nodes to join to when making new elements
          nvb = which_child(ne,ne0)  ! the 2nd child branch
       endif
       
       if(.not.bifurcation_element(ne))then
!!!    ---- for a single element that is the parent of a single element:
          ! make 4 new nodes using 'new_coords_derivs' at the end of the 1d element

          ! Apply rotation and translation to each new node from the generic template
          radius_weighted = 0.0_dp
          call mesh_rotate_vector_about_axis(4,template_vrsns,angle,length,radius, &
               radius_weighted,Rx,Ry,Txyz,template_coords,new_coords_derivs)
       
          ! as each new point is placed, check whether it is the closest to the first
          ! surface mesh node in ring of nodes surrounding the element, i.e. around
          ! node np1. This is used to determine the node numbering in the new element.
          smallest = 1.0e+6_dp
          point1(:) = node_xyz_2d(1,1,:,elem_node_map(nvb,1,np1)) 
          do i = 1,4 ! four nodes in the end 'ring'
             np_new = np_new + 1
             nodes_2d(np_new) = np_new 
             node_xyz_2d(:,1,:,np_new) = new_coords_derivs(:,1,:,i) ! coordinates and derivatives 1,2,1_2
             point2(:) = node_xyz_2d(1,1,:,np_new)
             distance = distance_between_points(point1,point2)
             if(distance.lt.smallest)then
                smallest = distance
                np_close(1) = np_new
             endif
          enddo ! i = 1,4
          ! the closest new node (in np_close(1)) to the first node in the ring around np1
          ! becomes the first node in the new ring. This aims to minimise distortion of
          ! new surface elements.
          elem_node_map(1,1,np2) = np_close(1) ! map the closest new node
          elem_node_map(1,5,np2) = np_close(1) ! and also store at the dummy node (5th)
          ! for the other three new nodes, map their order around the new ring and store
          ! this ring of nodes at the end node of the underlying 1d element
          do j = 2,4
             if(np_close(1)+j-1.gt.np_new)then
                elem_node_map(1,j,np2) = np_close(1)+j-5
             else
                elem_node_map(1,j,np2) = np_close(1)+j-1
             endif
          enddo ! j = 2,4
          
          ! make new elements, using the node mapping stored in elem_node_map.
          ! nodes 1 and 2 use mapping stored at 1d node np1, and nodes 3 and 4
          ! use mapping stored at 1d node np2. if the underlying 1d element is 
          ! a child branch in a bifurcation then we need to get the correct
          ! version numbers (for derivatives) from template_vrsn_map.
          do i = 1,4  ! 4 new elements
             ne_new = ne_new + 1
             elems_2d(ne_new) = ne_new  
             forall(nn = 1:2) elem_nodes_2d(nn,ne_new) = elem_node_map(nvb,nn+i-1,np1)
             forall(nn = 3:4) elem_nodes_2d(nn,ne_new) = elem_node_map(1,i+nn-3,np2)
             if(bifurcation_element(ne0)) & 
                forall(nn = 1:2) elem_versn_2d(nn,ne_new) = template_vrsn_map(nn,i+(nvb-1)*4)
          enddo ! i = 1,4
          forall (i=1:4) elem_node_map(1,i,np2) = np_new-4+i
          elem_node_map(1,5,np2) = elem_node_map(1,1,np2)
          
       else if(bifurcation_element(ne)) then
!!!    ---- for an element that is the parent of a bifurcation:
          ! create surface mesh nodes around the 1d element end node and at the crux.
          
          ! Apply rotation and translation to each new node from the generic template
          radius_weighted = 0.5_dp*radius+0.5_dp*max(elem_field(ne_radius,elem_cnct(1,1,ne)), &
               elem_field(ne_radius,elem_cnct(1,2,ne))) ! average of branch and largest child
          call mesh_rotate_vector_about_axis(5,template_vrsns,angle,length,radius, &
               radius_weighted,Rx,Ry,Txyz,template_coords,new_coords_derivs)
          new_coords_derivs(3,1,:,1) = new_coords_derivs(3,2,:,1)
          new_coords_derivs(3,1,:,3) = new_coords_derivs(3,2,:,3)
          
          ! adjust location of crux node using the location of the end node of
          ! the underlying 1d element, and end nodes of two child elements
          call mesh_2d_from_1d_crux(elem_nodes(2,ne),elem_cnct(1,1,ne), &
               elem_cnct(1,2,ne),cruxdist,new_coords_derivs)
          
          do i = 1,5  ! five nodes in the bifurcation
             np_new = np_new+1
             nodes_2d(np_new) = np_new 
             elem_node_map(2,i,np2) = np_new ! record nodes in ring around np2
             node_versn_2d(np_new) = template_vrsns(i)
             ring_coords(:,i) = new_coords_derivs(1,1,:,i)  ! store the coords for 'ring'
             node_xyz_2d(:,:,:,np_new) = new_coords_derivs(:,:,:,i) ! new node coords and derivs 1,2,1_2
          enddo ! i = 1,5
          elem_node_map(2,5,np2) = elem_node_map(2,1,np2) ! store mapping at dummy node

          np_side1 = np_new-4  ! first side node of bifurcation
          np_side2 = np_new-2  ! second side node of bifurcation

          ! find the node number in the start 'ring' (i.e. the ring of nodes at np1
          ! for the underlying 1d element) that is closest to the first side node (np_side1)
          smallest = 1.0e+6_dp
          point1(:) = node_xyz_2d(1,1,:,np_side1)
          do i = 1,4 ! four nodes in the end 'ring'
             point2(:) = node_xyz_2d(1,1,:,elem_node_map(nvb,i,np1)) &
                     + (node_xyz(:,np2)-node_xyz(:,np1))
             distance = distance_between_points(point1,point2)
             if(distance.lt.smallest)then
                smallest = distance
                np_close(1) = elem_node_map(nvb,i,np1)
             endif
          enddo

          ! adjust the node mapping around np1, so that starting from the node that is
          ! closest to the first side node
          NMAX = max(elem_node_map(nvb,1,np1),elem_node_map(nvb,2,np1), &
               elem_node_map(nvb,3,np1),elem_node_map(nvb,4,np1))
          elem_node_map(nvb,1,np1) = np_close(1)
          elem_node_map(nvb,5,np1) = elem_node_map(nvb,1,np1)
          do j = 2,4
             if(np_close(1)+j-1.gt.NMAX)then
                elem_node_map(nvb,j,np1) = np_close(1)+j-5
             else
                elem_node_map(nvb,j,np1) = np_close(1)+j-1
             endif
          enddo !j
          
          ! create new surface elements joining the np1 and np2 rings
          do i = 1,4
             ne_new = ne_new+1
             elems_2d(ne_new) = ne_new
             forall(nn = 1:2) elem_nodes_2d(nn,ne_new) = elem_node_map(nvb,nn+i-1,np1)
             forall(nn = 3:4) elem_nodes_2d(nn,ne_new) = np_new-5+(nn-3)+i
             elem_node_map(1:2,i,np2) = np_new-5+i !record nodes in ring 1 around np2
          enddo ! i
          elem_nodes_2d(2,ne_new) = elem_node_map(nvb,1,np1)
          elem_nodes_2d(4,ne_new) = np_new-4
          
          elem_node_map(1,1,np2) = np_new
          elem_node_map(2,3,np2) = np_new
          elem_node_map(1:2,5,np2) = elem_node_map(1:2,1,np2)

          np_crux = np_new
          
          ! for an element that is the parent of a bifurcation, also create
          ! nodes and elements along each of the two child branches
          np0 = np1  ! the start node of the parent element
          do k = 1,2 ! for each of two child branches
             ne_child = elem_cnct(1,k,ne) ! child element of the 1d branch
             np1 = elem_nodes(1,ne_child) ! child element start node
             np2 = elem_nodes(2,ne_child) ! child element end node
             length = max(radius*1.5_dp, cruxdist*1.25_dp, 0.5_dp*(cruxdist*1.25_dp+ &
                  (elem_field(ne_length,ne_child) - radius*0.5_dp)))
!             if(abs(angle_btwn_vectors(elem_direction(:,ne_child), &
!                  elem_direction(:,ne))/pi*180.0_dp).gt.60.0_dp)then
!                length = radius * 1.1_dp
!             elseif(elem_field(ne_radius,ne_child).lt.0.75_dp*radius)then
!                length = radius * 1.1_dp
!             endif
             ring_distance(ne_child) = length
             
             ! calculate the rotation angle for the branch
             call mesh_rotate_about_axis(ne_child,angle,Rx,Ry,Txyz,template_coords)
             
             radius = elem_field(ne_radius,ne_child)
             radius_weighted = radius
             call mesh_rotate_vector_about_axis(4,template_vrsns,angle,length,radius, &
                  radius_weighted,Rx,Ry,Txyz,template_coords,new_coords_derivs)
             
             do i = 1,4  ! four nodes in the new 'ring'
                np_new = np_new + 1
                nodes_2d(np_new) = np_new ! store new node number
                node_xyz_2d(:,1,:,np_new) = new_coords_derivs(:,1,:,i) ! coordinates and derivatives 1,2,1_2
             enddo ! i
             
             ! the bifurcation offers two potential rings of nodes to attach to.
             ! determine which is the correct ring using the angle between the
             ! direction of the child element and the direction from the start
             ! of the child element to the centre of each ring.
             point1(:) = node_xyz(:,np2)  ! coordinates at end of child branch
             point2(:) = node_xyz(:,np1)  ! coordinates at start of child branch
             point3 = 0.0_dp
             do i = 1,4
                point3(:) = point3(:) + ring_coords(:,template_cnct(1,i))/4.0_dp
             enddo !i
             angle = angle_btwn_points(point1,point2,point3)

             point3 = 0.0_dp
             do i = 1,4
                point3(:) = point3(:) + ring_coords(:,template_cnct(1,i+4))/4.0_dp
             enddo !i

             nvb = 1
             if(angle.gt.angle_btwn_points(point1,point2,point3)) nvb = 2
             
             ! determine which of the nodes in the bifurcation ring the new nodes
             ! should attach to based on which of the new nodes is the closest to
             ! the 5th (for nvb=1) or 1st (for nvb=2) ring node
             smallest = 1.0e+6_dp
             do i = 1,4
                vector_1(:) = new_coords_derivs(1,1,:,i) - ring_coords(:,2)
                vector_1 = unit_vector(vector_1)
                vector_2(:) = elem_direction(:,ne_child)
                angle = angle_btwn_vectors(vector_1,vector_2)
                if(angle.lt.smallest)then
                   smallest = angle
                   j = i - 4 - 1
                   if(j.eq.-4) j = 0
                   np_close(nvb) = np_new + j
                endif
             enddo !i
             
             ! record the node mapping in the ring
             elem_node_map(k,1,np1) = np_close(nvb)
             elem_node_map(k,5,np1) = elem_node_map(k,1,np1)
             do j = 2,4
                if(np_close(nvb)+j-1.gt.np_crux+4*k)then
                   elem_node_map(k,j,np1) = np_close(nvb)+j-5
                else
                   elem_node_map(k,j,np1) = np_close(nvb)+j-1
                endif
             enddo !j

             ! set up new elements
             do i = 1,4
                ne_new = ne_new+1
                elems_2d(ne_new) = ne_new
                do nn = 1,2
                   np_now = template_cnct(nn,i+(nvb-1)*4)+np_crux-5
                   elem_nodes_2d(nn,ne_new) = np_now
                   elem_versn_2d(nn,ne_new) = template_vrsn_map(nn,i+(nvb-1)*4)
                   if(np_now.eq.np_side1.or.np_now.eq.np_side2) elem_versn_2d(nn,ne_new) = 2
                enddo !nn
                forall(nn = 3:4) elem_nodes_2d(nn,ne_new) = elem_node_map(k,i+nn-3,np1)
             enddo !i
          enddo !k
          
       endif ! .not.bifurcation(ne)
       
       ne_count = ne_count+1
       if(ne_count.gt.num_elems)then
          ne = 0
       else
          if(ne_count .gt. count(elem_list.ne.0))then
             ne = 0
          else
             ne_global = elem_list(ne_count)
             ne = get_local_elem_1d(ne_global)
          endif
       endif
    enddo ! while (ne /= 0)
    
    num_nodes_2d = np_new
    num_elems_2d = ne_new

    call element_connectivity_2d
    call line_segments_for_2d_mesh('arcl')

    ! following is original for airways
    ! make sure that the nodes are nicely distributed along the branches
    call redistribute_mesh_nodes_2d_from_1d

    ! following is new because of vessel angles
    ! merge unnecessary elements. Most branches will end up with one element in Xi+2
    call merge_2d_from_1d_mesh
    ! remove very short branches where a trifurcation is more appropriate
    !call merge_trifurcations(short_elements)
    !call redistribute_mesh_nodes_2d_from_1d

    deallocate(elem_node_map)
    deallocate(ring_distance)
    deallocate(short_elements)
    
    call enter_exit(sub_name,2)
    
  end subroutine make_2d_vessel_from_1d0

!!!#############################################################################

  subroutine mesh_2d_from_1d_generic(template_coords)
    !*mesh_2d_from_1d_generic:* provide the coordinates of the templated nodes
    ! for a surface mesh creation about a 1d tree mesh

    real(dp) :: template_coords(:,:,:,:)
    !     Local Variables
    integer :: i,j,nj
    real(dp) :: angle,s1_dir(3)=0.0_dp,s2_dir(3)=0.0_dp
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'mesh_2d_from_1d_generic'
    call enter_exit(sub_name,1)
      
    template_coords = 0.0_dp

    !.....Set up a default ring of nodes. 
    do i = 1,4 ! step around the ring of nodes in four steps
       angle = 2.0_dp * pi * dble(i-1)/4.0_dp !=2*PI*step/number_of_steps
       !.......Derivatives
       s1_dir(1) = -sin(angle) !dxi_1/dx
       s1_dir(2) = cos(angle)  !dxi_1/dy
       s2_dir(3) = 0.3_dp      !dxi_2/dz
       template_coords(1,:,1,i) = s1_dir(2) !cos(angle)
       template_coords(1,:,2,i) = -s1_dir(1) !sin(angle)
       forall (nj = 1:3) template_coords(2,:,nj,i) = s1_dir(nj)
       forall (nj = 1:3) template_coords(3,:,nj,i) = s2_dir(nj)
    enddo !i

    !.....First side node
    template_coords(1,:,3,1) = -1.5_dp   ! side node
    template_coords(3,2,1,1) = 0.15_dp
    template_coords(3,2,2,1) = 0.0_dp
    template_coords(3,2,3,1) = 0.7_dp

    !.....Front bifurcation node
    template_coords(2,2,1:2,2) = 0.0_dp
    template_coords(2,2,3,2) = -1.0_dp
    template_coords(2,4,:,2) = -template_coords(2,2,:,2)
    template_coords(3,2:3,1,2) = -0.33_dp
    template_coords(3,2:5,2,2) = 0.0_dp
    template_coords(3,2:5,3,2)= 0.33_dp
    template_coords(3,4:5,1,2)= 0.33_dp
    
    !.....Second side node
    template_coords(1,:,3,3) = -1.5_dp   ! side node
    template_coords(3,2,1,3) = -0.15_dp
    template_coords(3,2,2:3,3) = template_coords(3,2,2:3,1)

    !.....Back bifurcation node
    template_coords(2,2,1:2,4) = 0.0_dp
    template_coords(2,2,3,4) = 1.0_dp
    template_coords(2,4,:,4) = -template_coords(2,2,:,4)
    template_coords(3,2:3,1,4) = -0.33_dp
    template_coords(3,2:5,2,4) = 0.0_dp
    template_coords(3,2:5,3,4)= 0.33_dp
    template_coords(3,4:5,1,4)= 0.33_dp
    
    !.....Crux node
    template_coords(1,:,3,5) = 1.5_dp
    template_coords(2,1,2,5) = 1.0_dp
    template_coords(2,2,:,5) = -template_coords(2,1,:,5)
    template_coords(3,1,1,5) = -0.3_dp
    template_coords(3,2,1,5) = 0.3_dp
    
    call enter_exit(sub_name,2)

  end subroutine mesh_2d_from_1d_generic

!!!#############################################################################

  subroutine mesh_rotate_about_axis(ne,angle_z,Rx,Ry,Txyz,template_coords)
    !*mesh_rotate_about_axis:* calculates the rotation matrices and z-angle for
    ! rotation of a vector about an arbitrary axis defined by the direction
    ! of element ne.
    
    integer,intent(in) :: ne
    real(dp) :: angle_z,Rx(:,:),Ry(:,:),Txyz(:)
    real(dp),intent(in) :: template_coords(:,:,:,:)
    ! Local Variables
    integer :: j,ne_parent,ne0,nj,np0,np1,np2,np3,np4
    real(dp) :: angle_z2,ln_23,NML(4),V(3),X1(3),X2(3),X3(3)
    logical :: found
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'mesh_rotate_about_axis'
    call enter_exit(sub_name,1)
      
    np1 = elem_nodes(1,ne)
    np2 = elem_nodes(2,ne)
    
    ! Find next bifurcation nodes, for calculating rotation about the z-axis
    if(elem_cnct(1,0,ne).ge.2)then !get adjacent nodes
       np0 = np1
       np3 = elem_nodes(2,elem_cnct(1,1,ne)) !end node of first child
       np4 = elem_nodes(2,elem_cnct(1,2,ne)) !end node of second child
    else !find next bifurcation
       found = .false.
       ne0 = ne
       do while (.not.found)
          if(elem_cnct(1,0,ne0).eq.1)then
             ne0 = elem_cnct(1,1,ne0)
             np0 = elem_nodes(1,ne0)
          elseif(elem_cnct(1,0,ne0).ge.2)then
             found = .true.
             np3 = elem_nodes(2,elem_cnct(1,1,ne0))
             np4 = elem_nodes(2,elem_cnct(1,2,ne0))
          elseif(elem_cnct(1,0,ne0).eq.0)then
             found = .true.
             np3 = 0
             np4 = 0
          endif
       enddo
    endif
    
    ! find the predecessor element that is first proximal to a bifurcation
    ne_parent = elem_cnct(-1,1,ne)
    do while (elem_cnct(1,0,ne_parent).ne.2.and.ne_parent.ne.0)
       ne_parent = elem_cnct(-1,1,ne_parent)
    enddo
    
!!! np1 == start node, np2 == end node, np3 == end of child 1, np4 == end of child 2
!!! .....Calculate the rotation and translation matrices      
    call calc_rotat_trans_mats(ne,np1,Rx,Ry,Txyz)

    !.....The angle for rotation about the z-axis is equal to the angle
    !.....between the normal to the plane containing bifurcation nodes and
    !.....the theta=0 direction
    if(np3.ne.0)then
!!!    get the normal to the plane (NML)
       X1(1:3) = node_xyz(1:3,np0)
       X2(1:3) = node_xyz(1:3,np3)
       X3(1:3) = node_xyz(1:3,np4)
       call make_plane_from_3points(NML,2,X1,X2,X3)

       V(:) = template_coords(1,1,:,2)   ! direction of 'front' bifurcation node == 0,1,0
       V = unit_vector(V)
       !.......Calculate location of V if rotation about x- and y- was applied
       X1 = 0.0_dp
       X2 = 0.0_dp
       do nj=1,3
          do j=1,3
             X1(nj) = X1(nj)+Ry(nj,j)*V(j)
          enddo !j
       enddo !nj
       do nj = 1,3
          do j = 1,3
             X2(nj) = X2(nj)+Rx(nj,j)*X1(j)
          enddo !j
       enddo !nj
       V = X2 !-Txyz(nj)
       V = unit_vector(V)   ! the direction of the vector rotated about x and y axes
       
       angle_z = min(scalar_product_3(V,NML),1.0_dp)
       angle_z = max(scalar_product_3(V,NML),-1.0_dp)

       if(angle_z.ge.1_dp)then
          angle_z = 0.0_dp
       elseif(angle_z.le.-1.0_dp)then
          angle_z = pi
       else
          angle_z = acos(angle_z) ! angle between normal and 2nd (front) node
       endif
       
       V(:) = template_coords(1,1,:,3)   ! direction to the 'side' node
       V = unit_vector(V)
       !.......Calculate location of V if rotation about x- and y- was applied
       X1 = 0.0_dp
       X2 = 0.0_dp
       do nj=1,3
          do j=1,3
             X1(nj) = X1(nj)+Ry(nj,j)*V(j)
          enddo !j
       enddo !nj
       do nj = 1,3
          do j = 1,3
             X2(nj) = X2(nj)+Rx(nj,j)*X1(j)
          enddo !j
       enddo !nj
       V = X2 !-Txyz(nj)
       V = unit_vector(V)
       angle_z2 = acos(scalar_product_3(V,NML)) ! angle between normal and 3rd (side) node
       if(angle_z2.lt.PI/2.0_dp)then
          angle_z = -angle_z
       endif
    else
       angle_z = 0.0_dp
    endif
    
    call enter_exit(sub_name,2)
    
  end subroutine mesh_rotate_about_axis
    
!!!#############################################################################

  subroutine calc_rotat_trans_mats(ne,np1,Rx,Ry,Txyz)
    !*calc_rotat_trans_mats:* calculates rotation and translation matrics for
    ! element ne with start node np1

    integer,intent(in) :: ne,np1
    real(dp) :: Rx(:,:),Ry(:,:),Txyz(:)
    real(dp) :: ln_23
    
    !!! translation = the displacement from (0,0,0) to np1 (= -node_xyz(np1))
    Txyz(1:3) = -node_xyz(1:3,np1)   ! translation of the generic nodes

!!! Rotation matrices: Rx(-1) = | 1  0    0  |   Ry(-1) = | d 0 -a |   Rz = | cos(t)  -sin(t) |
!!!                             | 0 c/d -b/d |            | 0 1  0 |        | sin(t)   cos(t) |
!!!                             | 0 b/d  c/d |            | a 0  d |        |   0        0    |
!!! where U(a,b,c) == branch direction (elem_direction(1:3,ne)) and d = sqrt(b2 + c2)
!!! see www.paulbourke.net/geometry/rotate for explanation
    
    Rx = 0.0_dp
    Ry = 0.0_dp
    Rx(1,1) = 1.0_dp !x-x = 1
    ln_23 = sqrt(elem_direction(2,ne)**2 + elem_direction(3,ne)**2)
    if(abs(ln_23).lt.zero_tol) ln_23 = 1.0_dp
    Rx(2,2) = elem_direction(3,ne)/ln_23
    Rx(2,3) = elem_direction(2,ne)/ln_23
    Rx(3,2) = -Rx(2,3)
    Rx(3,3) = Rx(2,2)
    
    Ry(2,2) = 1.0_dp !x-x = 1
    Ry(1,1) = ln_23
    Ry(1,3) = elem_direction(1,ne)
    Ry(3,1) = -Ry(1,3)
    Ry(3,3) = Ry(1,1)

  end subroutine calc_rotat_trans_mats
  
!!!#############################################################################

  subroutine mesh_rotate_about_axis_basic(np,axis,centre,theta)
    !*mesh_rotate_about_axis_basic:* calculates the rotation matrices and z-angle for
    ! rotation of a vector about an arbitrary axis
    ! of element ne.

    ! only set up for single versioned nodes
    
    integer,intent(in) :: np
    real(dp) :: axis(3),centre(3),theta
    ! Local Variables
    integer :: j,nj,nk
    real(dp) :: d,p1(3),q1(3),q2(3),u(3),Rx(3,3),Ry(3,3),x1(3),x2(3),ln_23
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'mesh_rotate_about_axis'
    call enter_exit(sub_name,1)
      
    !/* Step 1 */
    p1 = centre
    q1(1:3) = node_xyz_2d(1,1,1:3,np) - centre(1:3)
    u = unit_vector(axis)
    d = sqrt(u(2)*u(2) + u(3)*u(3))

    !/* Step 2 */
    if (d.gt.zero_tol)then
       q2(1) = q1(1)
       q2(2) = q1(2) * u(3) / d - q1(3) * u(2) / d
       q2(3) = q1(2) * u(2) / d + q1(3) * u(3) / d
    else
       q2 = q1
    endif

    !/* Step 3 */
    q1(1) = q2(1) * d - q2(3) * u(1)
    q1(2) = q2(2)
    q1(3) = q2(1) * u(1) + q2(3) * d

    !/* Step 4 */
    q2(1) = q1(1) * cos(theta) - q1(2) * sin(theta)
    q2(2) = q1(1) * sin(theta) + q1(2) * cos(theta)
    q2(3) = q1(3)

    !/* Inverse of step 3 */
    q1(1) = q2(1) * d + q2(3) * u(1)
    q1(2) = q2(2)
    q1(3) = - q2(1) * u(1) + q2(3) * d

    !/* Inverse of step 2 */
    if (abs(d).gt.zero_tol)then
       q2(1) = q1(1)
       q2(2) = q1(2) * u(3) / d + q1(3) * u(2) / d
       q2(3) = - q1(2) * u(2) / d + q1(3) * u(3) / d
    else
       q2 = q1
    endif

    !/* Inverse of step 1 */
    q1(1) = q2(1) + p1(1)
    q1(2) = q2(2) + p1(2)
    q1(3) = q2(3) + p1(3)

    node_xyz_2d(1,1,1:3,np) = q1(1:3)

    Rx = 0.0_dp
    Ry = 0.0_dp
    Rx(1,1) = 1.0_dp !x-x = 1
    ln_23 = sqrt(axis(2)**2 + axis(3)**2)
    if(abs(ln_23).lt.zero_tol) ln_23 = 1.0_dp
    Rx(2,2) = axis(3)/ln_23
    Rx(2,3) = axis(2)/ln_23
    Rx(3,2) = -Rx(2,3)
    Rx(3,3) = Rx(2,2)
    
    Ry(2,2) = 1.0_dp !x-x = 1
    Ry(1,1) = ln_23
    Ry(1,3) = axis(1)
    Ry(3,1) = -Ry(1,3)
    Ry(3,3) = Ry(1,1)
    do nk = 2,4 ! for the coordinates and derivatives
       X1(1) = cos(theta)*node_xyz_2d(nk,1,1,np) + sin(theta)*node_xyz_2d(nk,1,2,np)
       X1(2) = -sin(theta)*node_xyz_2d(nk,1,1,np) + cos(theta)*node_xyz_2d(nk,1,2,np)
       node_xyz_2d(nk,1,1,np) = X1(1)
       node_xyz_2d(nk,1,2,np) = X1(2)
    enddo ! nk
    do nk = 2,4
       X1 = 0.0_dp
       X2 = 0.0_dp
       do nj = 1,3
          do j = 1,3
             X1(nj) = X1(nj)+Ry(nj,j)*node_xyz_2d(nk,1,j,np)
          enddo !j
       enddo !nj
       do nj = 1,3
          do j = 1,3
             X2(nj) = X2(nj)+Rx(nj,j)*X1(j)
          enddo !j
       enddo !nj
       node_xyz_2d(nk,1,1:3,np) = X2(1:3)
    enddo ! nk

    call enter_exit(sub_name,2)
    
  end subroutine mesh_rotate_about_axis_basic
    
!!!#############################################################################

  subroutine mesh_rotate_vector_about_axis(N,template_vrsns,angle,length, &
       radius,radius_weighted,Rx,Ry,Txyz,template_coords,new_coords_derivs)
    !*mesh_rotate_vector_about_axis:* rotates a vector (starting at (0,0,0))
    ! about an arbitrary axis. Rotation matrices are in Rx and Ry.
    ! angle is the amount to rotate about z-axis. 

    integer,intent(in) :: N,template_vrsns(:)
    real(dp),intent(in) :: angle,length,radius,radius_weighted,Rx(3,3), &
         Ry(3,3),Txyz(3),template_coords(:,:,:,:)
    real(dp) :: new_coords_derivs(:,:,:,:)
    ! Local variables
    integer :: i,j,nj,nk,nv
    real(dp) :: X1(3),X2(3)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'mesh_rotate_vector_about_axis'
    call enter_exit(sub_name,1)
      
    new_coords_derivs = 0.0_dp
    
    do i = 1,N !each node in ring, plus crux if at bifurcation
       do nk = 1,4
          do nv = 1,template_vrsns(i)
             new_coords_derivs(nk,nv,:,i) = template_coords(nk,nv,:,i)
             X1(1) = cos(angle)*new_coords_derivs(nk,nv,1,i) + &
                  sin(angle)*new_coords_derivs(nk,nv,2,i)
             X1(2) = -sin(angle)*new_coords_derivs(nk,nv,1,i) + &
                  cos(angle)*new_coords_derivs(nk,nv,2,i)
             new_coords_derivs(nk,nv,1,i) = X1(1)
             new_coords_derivs(nk,nv,2,i) = X1(2)
          enddo !nv
       enddo !nk
       do nv = 1,template_vrsns(i)
          if(N.eq.5)then ! calculating for a bifurcation
             if(i.eq.2.or.i.eq.4)then
                new_coords_derivs(1,nv,1:2,i) = new_coords_derivs(1,nv,1:2,i)*radius_weighted
             else
                new_coords_derivs(1,nv,1:2,i) = new_coords_derivs(1,nv,1:2,i)*radius*1.05_dp
             endif
          else ! calculating for a non-bifurcation 'ring'
             new_coords_derivs(1,nv,1:2,i) = new_coords_derivs(1,nv,1:2,i)*radius
          endif
          if(N.lt.5)then
             new_coords_derivs(1,nv,3,i) = length
          elseif(N.EQ.5)then
             if(i.eq.3.or.i.eq.1)then
                new_coords_derivs(1,nv,3,i) = length-radius*0.5_dp
             else
                new_coords_derivs(1,nv,3,i) = length
             endif
          endif
       enddo !nv
       do nk = 1,4
          do nv = 1,template_vrsns(i)
             X1 = 0.0_dp
             X2 = 0.0_dp
             do nj = 1,3
                do j = 1,3
                   X1(nj) = X1(nj)+Ry(nj,j)*new_coords_derivs(nk,nv,j,i)
                enddo !j
             enddo !nj
             do nj = 1,3
                do j = 1,3
                   X2(nj) = X2(nj)+Rx(nj,j)*X1(j)
                enddo !j
             enddo !nj
             do nj = 1,3
                if(nk.eq.1)then !geometry
                   new_coords_derivs(nk,nv,nj,i) = X2(nj)-Txyz(nj)
                else !derivatives
                   new_coords_derivs(nk,nv,nj,i) = X2(nj)
                endif
             enddo !nj
          enddo !nv
       enddo !nk
    enddo !i

    call enter_exit(sub_name,2)
    
  end subroutine mesh_rotate_vector_about_axis

!!!#############################################################################

  subroutine mesh_2d_from_1d_crux(np0,ne_1,ne_2,distance,new_coords_derivs)
    !*mesh_2d_from_1d_crux:* adjusts the location of the crux node for deriving
    ! a 2d surface mesh over a 1d tree mesh.
    
    integer,intent(in) :: ne_1,ne_2,np0
    real(dp) :: distance,new_coords_derivs(:,:,:,:)
    ! Local variables
    real(dp) :: angle,angle_V,angle_U,gamma,length_U,length_V,matrix(3,3),&
         N(3),radius_U,radius_V,W(3),vector(3),weight,width
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'mesh_2d_from_1d_crux'
    call enter_exit(sub_name,1)
      
    !adjust location of crux node
    radius_U = elem_field(ne_radius,ne_1)
    radius_V = elem_field(ne_radius,ne_2)

    ! get the direction of the normal to the two branch directions (N)
    N = cross_product(elem_direction(1,ne_1),elem_direction(1,ne_2))
    N = unit_vector(N) ! normalise

    ! get the angle between the two branches
    angle = scalar_product_3(elem_direction(1,ne_1),elem_direction(1,ne_2))
    angle = acos(angle)

    width = sqrt(radius_U**2 +radius_V**2 - 2.0_dp*radius_U*radius_V*cos(pi-angle))

    gamma = asin(radius_U*sin(pi-angle)/width)
    length_U = width*sin(pi/2.0_dp-gamma)/sin(angle)
    gamma = asin(radius_V*sin(pi-angle)/width)
    length_V = width*sin(pi/2.0_dp-gamma)/sin(angle)

    ANGLE_U = atan(radius_U/length_U)
    ANGLE_V = atan(radius_V/length_V)

!!!..... n.w = 0
!!!..... u.w = cos(angle_u) * Lw
!!!..... v.w = cos(angle_v) * Lw
      
    matrix(1,:) = N(:)
    matrix(2,:) = elem_direction(:,ne_1)
    matrix(3,:) = elem_direction(:,ne_2)
    vector(1) = 0.0_dp
    vector(2) = cos(ANGLE_U)*sqrt(length_U**2+radius_U**2)
    vector(3) = cos(ANGLE_V)*sqrt(length_U**2+radius_U**2)

    w = mesh_a_x_eq_b(matrix,vector)
    w = unit_vector(w)

    distance = sqrt(length_U**2 +radius_U**2)
    if(angle.gt.pi)then
       write(*,'( ''WARNING - angle exceeds pi'')')
    else if(angle.gt.pi/2.0_dp)then
       weight = 1.1_dp
    else if(angle.gt.3.0_dp*pi/8.0_dp)then
       weight = 1.0_dp
    else if(angle.gt.pi/4.0_dp)then
       weight = 0.9_dp
    else
       weight = 0.75_dp
    endif
    
    new_coords_derivs(1,1,1:3,5) = node_xyz(1:3,np0) + W(1:3) * distance * weight
    new_coords_derivs(1,2,:,5) = new_coords_derivs(1,1,:,5)

    call enter_exit(sub_name,2)
    
  end subroutine mesh_2d_from_1d_crux

!!!###########################################################################

  subroutine merge_2d_from_1d_mesh
    !*merge_2d_from_1d_mesh:* used for deriving a 2d surface mesh over a 1d
    ! tree mesh. merges the elements in the 'rings' that are adjacent to the
    ! bifurcations in the +Xi2 direction

    ! Local variables
    integer :: j,ne,ne_parent,ne_parent2,non_zero,np,np1,np2
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'merge_2d_from_1d_mesh'
    call enter_exit(sub_name,1)

!    ne = 1
!    do while (ne.le.num_elems_2d)
!       if(elems_2d(ne).gt.0)then
!          if(elem_cnct_2d(2,0,ne).gt.0)then
!             np1 = elem_nodes_2d(3,ne)
!             np2 = elem_nodes_2d(4,ne)
!             if(node_versn_2d(np1).gt.1.and.node_versn_2d(np2).gt.1)then
!                ! merge with the previous ring
!                do j = 1,4
!                   ne_parent = elem_cnct_2d(-2,1,ne)  ! get element proximal to current
!                   elem_nodes_2d(1:2,ne) = elem_nodes_2d(1:2,ne_parent)
!                   elem_versn_2d(1:2,ne) = elem_versn_2d(1:2,ne_parent)
!                   elem_cnct_2d(-2,0:1,ne) = elem_cnct_2d(-2,0:1,ne_parent)
!                   if(elem_cnct_2d(-2,0,ne_parent).ne.0)then
!                      ne_parent2 = elem_cnct_2d(-2,1,ne_parent)
!                      elem_cnct_2d(2,0:1,ne_parent2) = elem_cnct_2d(2,0:1,ne_parent)
!                   endif
!                   nodes_2d(elem_nodes_2d(3,ne_parent)) = 0
!                   node_xyz_2d(:,:,:,elem_nodes_2d(3,ne_parent)) = 0.0_dp
!                   elems_2d(ne_parent) = 0
!                   elem_nodes_2d(:,ne_parent) = 0
!                   ne = ne + 1
!                enddo
!             else
!                ne = ne + 4 ! skip to next ring
!             endif
!          else
!             ne = ne + 4 ! skip to next ring
!          endif
!       else
!          ne = ne + 4 ! skip to next ring
!       endif
!    enddo

    ne = num_elems_2d
    do while (ne.gt.1)
       if(elems_2d(ne).gt.0)then
          if(elem_cnct_2d(2,0,ne).gt.0)then
             np1 = elem_nodes_2d(1,ne)
             np2 = elem_nodes_2d(2,ne)
             if(node_versn_2d(np1).gt.1.and.node_versn_2d(np2).gt.1)then
                ! merge with the next ring
                do j = 1,4
                   ne_parent = elem_cnct_2d(2,1,ne)  ! get element distal to current
                   elem_nodes_2d(3:4,ne) = elem_nodes_2d(3:4,ne_parent)
                   elem_versn_2d(3:4,ne) = elem_versn_2d(3:4,ne_parent)
                   elem_cnct_2d(2,0:1,ne) = elem_cnct_2d(2,0:1,ne_parent)
                   if(elem_cnct_2d(2,0,ne_parent).ne.0)then
                      ne_parent2 = elem_cnct_2d(2,1,ne_parent)
                      elem_cnct_2d(-2,0:1,ne_parent2) = elem_cnct_2d(-2,0:1,ne_parent)
                   endif
                   nodes_2d(elem_nodes_2d(1,ne_parent)) = 0
                   node_xyz_2d(:,:,:,elem_nodes_2d(1,ne_parent)) = 0.0_dp
                   elems_2d(ne_parent) = 0
                   elem_nodes_2d(:,ne_parent) = 0
                   ne = ne - 1
                enddo
             else
                ne = ne - 4 ! skip to next ring
             endif
          else
             ne = ne - 4 ! skip to next ring
          endif
       else
          ne = ne - 4 ! skip to next ring
       endif
    enddo
    
    non_zero = 0
    do np = 1,num_nodes_2d
       if(nodes_2d(np).ne.0)then
          non_zero = non_zero + 1
       endif
    enddo
    num_nodes_2d = non_zero
    
    non_zero = 0
    do ne = 1,num_elems_2d
       if(elems_2d(ne).ne.0)then
          non_zero = non_zero + 1
       endif
    enddo
    num_elems_2d = non_zero

    call enter_exit(sub_name,2)

  end subroutine merge_2d_from_1d_mesh

!!!#############################################################################

  subroutine merge_2d_element(ndirection,ne)
    !*merge_element:* Merges a 2d element with the neighbouring element in the 
    ! ndirection Xi direction
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_MERGE_2D_ELEMENT" :: MERGE_2D_ELEMENT
    
    integer :: ndirection,ne
    ! Local parameters
    integer :: j,ne_parent,ne_parent2
    
    do j = 1,4
       ne_parent = elem_cnct_2d(ndirection,1,ne)  ! get adjacent element in ndirection
       if(ndirection.eq.-2)then
          elem_nodes_2d(1:2,ne) = elem_nodes_2d(1:2,ne_parent)
          elem_versn_2d(1:2,ne) = elem_versn_2d(1:2,ne_parent)
          elem_cnct_2d(-2,0:1,ne) = elem_cnct_2d(-2,0:1,ne_parent)
          if(elem_cnct_2d(-2,0,ne_parent).ne.0)then
             ne_parent2 = elem_cnct_2d(-2,1,ne_parent)
             elem_cnct_2d(2,0:1,ne_parent2) = elem_cnct_2d(2,0:1,ne_parent)
          endif
          nodes_2d(elem_nodes_2d(3,ne_parent)) = 0
          node_xyz_2d(:,:,:,elem_nodes_2d(3,ne_parent)) = 0.0_dp
       endif
       elems_2d(ne_parent) = 0
       elem_nodes_2d(:,ne_parent) = 0
       ne = ne + 1
    enddo !j
    num_elems_2d = num_elems_2d - 4
    num_nodes_2d = num_nodes_2d - 4
    
  end subroutine merge_2d_element

!!!#############################################################################

  subroutine merge_trifurcations(short_elements)
    !*merge_trifurcations:* used when deriving a 2d surface mesh from a 1d tree
    ! mesh. merges short 2d elements with neighbouring elements, and removes.
    ! this is required for two bifurcations that occur with a very
    ! short distance between them. i.e. when approximately a trifurcation.

    integer,intent(in) :: short_elements(:)
    ! Local variables
    integer :: k,ne,ne_child,ne_next,ne_parent,np_current,np_next
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'merge_trifurcations'
    call enter_exit(sub_name,1)

    k = 1
    ne = short_elements(k)

    ! for each of the 'short' surface elements, delete them by joining the
    ! parent and child elements together. Work on a group of 4 elements at
    ! a time, where the four elements run in the Xi1 direction == a ring
    ! of elements in a cylinder . 
    
    do while (ne /= 0)
       ne_parent = elem_cnct_2d(-2,1,ne)  ! element proximal to the current one
       ne_child = elem_cnct_2d(2,1,ne)    ! element distal to the current one
       elem_nodes_2d(1,ne_child) = elem_nodes_2d(1,ne)
       elem_versn_2d(1,ne_child) = elem_versn_2d(1,ne)
       elem_nodes_2d(4,ne_parent) = elem_nodes_2d(4,ne)
       elem_versn_2d(4,ne_parent) = elem_versn_2d(4,ne)
       elems_2d(ne) = 0 ! remove the current element from the 2d element list

       k = k + 1
       ne = short_elements(k)
       ne_parent = elem_cnct_2d(-2,1,ne)
       ne_child = elem_cnct_2d(2,1,ne)
       elem_nodes_2d(2,ne_child) = elem_nodes_2d(2,ne)
       elem_versn_2d(2,ne_child) = elem_versn_2d(2,ne)
       elem_nodes_2d(3,ne_parent) = elem_nodes_2d(3,ne)
       elem_versn_2d(3,ne_parent) = elem_versn_2d(3,ne)
       elems_2d(ne) = 0
       
       k = k + 1
       ne = short_elements(k)
       ne_parent = elem_cnct_2d(-2,1,ne)
       ne_child = elem_cnct_2d(2,1,ne)
       elem_nodes_2d(1,ne_child) = elem_nodes_2d(1,ne)
       elem_versn_2d(1,ne_child) = elem_versn_2d(1,ne)
       elem_nodes_2d(2,ne_child) = elem_nodes_2d(2,ne)
       elem_versn_2d(2,ne_child) = elem_versn_2d(2,ne)

       ! update the 1st node in the neighbouring Xi+1 element
       ne_next = elem_cnct_2d(1,1,ne_child)
       np_next = elem_nodes_2d(2,ne)
       np_current = elem_nodes_2d(1,ne_next)
       elem_nodes_2d(1,ne_next) = np_next  ! replace np_current with np_next

       node_versn_2d(np_next) = node_versn_2d(np_next) + 1 ! add a version
       ! deriv in Xi1 == same as current node; deriv in Xi2 == same as new node
       node_xyz_2d(:,node_versn_2d(np_next),:,np_next) = node_xyz_2d(:,1,:,np_next)
       node_xyz_2d(2,node_versn_2d(np_next),:,np_next) = node_xyz_2d(2,elem_versn_2d(1,ne_next),:,np_current)
       node_xyz_2d(3,node_versn_2d(np_next),:,np_next) = node_xyz_2d(3,elem_versn_2d(2,ne_child),:,np_next)
       elem_versn_2d(1,ne_next) = node_versn_2d(np_next)
       elems_2d(ne) = 0

       ne_next = elem_cnct_2d(-1,1,ne_child)
       node_xyz_2d(3,elem_versn_2d(2,ne_next),:,np_next) = node_xyz_2d(3,node_versn_2d(np_next),:,np_next)

       k = k + 1
       ne = short_elements(k)
       ne_parent = elem_cnct_2d(-2,1,ne)
       ne_child = elem_cnct_2d(2,1,ne)
       elem_nodes_2d(1,ne_child) = elem_nodes_2d(1,ne)
       elem_versn_2d(1,ne_child) = elem_versn_2d(1,ne)
       elem_nodes_2d(2,ne_child) = elem_nodes_2d(2,ne)
       elem_versn_2d(2,ne_child) = elem_versn_2d(2,ne)

       ! update the 2nd node in the neighbouring Xi-1 direction
       ne_next = elem_cnct_2d(-1,1,ne_child)
       np_next = elem_nodes_2d(1,ne)
       np_current = elem_nodes_2d(2,ne_next)
       node_versn_2d(np_next) = node_versn_2d(np_next) + 1 ! add a version
       node_xyz_2d(:,node_versn_2d(np_next),:,np_next) = node_xyz_2d(:,1,:,np_next)
       node_xyz_2d(2,node_versn_2d(np_next),:,np_next) = -node_xyz_2d(2,node_versn_2d(np_next)-1,:,np_next)
       node_xyz_2d(3,node_versn_2d(np_next),:,np_next) = node_xyz_2d(3,elem_versn_2d(2,ne_next),:,np_current)
       elem_nodes_2d(2,ne_next) = np_next                    ! change the 2nd node number of adjacent element
       elem_versn_2d(2,ne_next) = node_versn_2d(np_next)
       elems_2d(ne) = 0

       ne_next = elem_cnct_2d(1,1,ne_next)
       node_xyz_2d(3,elem_versn_2d(1,ne_next),:,np_next) = node_xyz_2d(3,node_versn_2d(np_next),:,np_next)

       k = k + 1
       ne = short_elements(k)
    enddo
    
    num_nodes_2d = count(nodes_2d.ne.0)
    num_elems_2d = count(elems_2d.ne.0)
       
    call enter_exit(sub_name,2)
    
  end subroutine merge_trifurcations

!!!#############################################################################

  subroutine define_rad_from_file(FIELDFILE, radius_type_in)
    !*define_rad_from_file:* reads in a radius field associated with an 
    ! airway tree and assigns radius information to each element, also 
    ! calculates volume of each element
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_FILE" :: DEFINE_RAD_FROM_FILE

    character(len=MAX_FILENAME_LEN), intent(in) :: FIELDFILE
    character(len=MAX_STRING_LEN), optional ::  radius_type_in
    !     Local Variables
    integer :: ierror,ne,ne_counter,ne_global,np,np1,np2,np_global, &
         num_elem_rad,surround
    real(dp) :: radius
    logical :: node_based,versions
    character(len=MAX_STRING_LEN) ::  radius_type
    character(LEN=132) :: ctemp1
    character(len=250) :: readfile
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'define_rad_from_file'
    call enter_exit(sub_name,1)
    
    versions = .TRUE.
    if(present(radius_type_in))then
       radius_type = radius_type_in
    else
       radius_type = 'no_taper'
    endif
    
    if(index(FIELDFILE, ".ipfiel")> 0) then !full filename is given
       readfile = FIELDFILE
    else ! need to append the correct filename extension
       readfile = trim(FIELDFILE)//'.ipfiel'
    endif
    
    open(10, file=readfile, status='old')

!!! check whether reading in a node-based or element-based field
    check_type : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "of elem")> 0) then ! element-based field
          node_based = .false.
          num_elem_rad = get_final_integer(ctemp1) !get global element number
          exit check_type
       else if(index(ctemp1, "of node")> 0) then ! node-based field
          node_based = .true.
          exit check_type
       endif
    enddo check_type

    if(node_based)then
       !.....check whether versions are prompted (>1)
       read_versions : do !define a do loop name
          read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
          if(index(ctemp1, "different")> 0) then !keyword "different" is found
             if(index(ctemp1, " N")> 0) then !keyword " N" is found
                versions=.false.
             endif
             exit read_versions !exit the named do loop
          endif
       end do read_versions
       
       np = 0
       !.....read the coordinate, derivative, and version information for each node.
       read_a_node : do !define a do loop name
          !.......read node number
          read(unit=10, fmt="(a)", iostat=ierror) ctemp1
          if(index(ctemp1, "Node")> 0) then
             np_global = get_final_integer(ctemp1) !get global node number
             ! find the corresponding local node number
             call get_local_node(np_global,np) ! get local node np for global node
             surround=elems_at_node(np,0)         !get number of surrounding elems
             ne=elems_at_node(np,1)  !First element at this node
             if(surround==1)then !only one element at this node so either a terminal or inlet
                if(radius_type.eq.'taper')then !inlet radius needs to be defined
                   if(elem_cnct(-1,0,ne).eq.0)then!Inlet as it has no parent need to set up radius into this vessel
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                      if(index(ctemp1, "version number")>0) then
                         read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                      endif
                      if(index(ctemp1, "value")> 0) then
                         radius = get_final_real(ctemp1)
                         elem_field(ne_radius_in,ne) = radius
                      endif
                   endif
                endif
                if(elem_cnct(-1,0,ne).eq.0)cycle      !No parent therefore inlet. Skip and go to the next
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                if(index(ctemp1, "version number")>0) then
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                endif
                if(index(ctemp1, "value")> 0) then
                   radius = get_final_real(ctemp1)
                   if(radius_type.eq.'taper')then
                      elem_field(ne_radius_out,ne) = radius
                   else
                      elem_field(ne_radius,ne) = radius
                   endif
                endif
             elseif(surround.gt.1)then !Terminal airway - use first radius
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                if(index(ctemp1, "value")> 0) then
                   radius = get_final_real(ctemp1)
                   if(radius_type.eq.'taper')then
                      elem_field(ne_radius_out,ne) = radius
                   else
                      elem_field(ne_radius,ne) = radius
                   endif
                endif
             endif
          endif !index
          if(np.ge.num_nodes) exit read_a_node
       end do read_a_node

    else ! for element_based field file

       ne = 0
       ne_counter = 0
       
       read_an_elem : do !define a do loop name
          !.......read element number
          read(unit=10, fmt="(a)", iostat=ierror) ctemp1
          if(index(ctemp1, "Element number")> 0) then
             ne_global = get_final_integer(ctemp1) !get global element number
             ne = get_local_elem_1d(ne_global) ! get local elem ne for global elem
             if(ne.eq.0)then
                write(*,'('' WARNING! no local element for global element'',i6)') ne_global
                read(*,*)
             endif
             if(ne.gt.0)then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                if(index(ctemp1, "value")> 0) then
                   radius = get_final_real(ctemp1)
                   elem_field(ne_radius,ne) = radius
                   elem_field(ne_radius_in,ne) = radius
                   elem_field(ne_radius_out,ne) = radius
                endif
                ne_counter = ne_counter + 1
             endif
          endif !index
          if(ne_counter.ge.num_elem_rad) exit read_an_elem
       end do read_an_elem

    endif
    
!!! Calculate element volume
    do ne = 1,num_elems
       if(radius_type.eq.'taper')then
          if(elem_cnct(-1,0,ne).ne.0)then !radius in is radius of upstream vessel
             elem_field(ne_radius_in,ne)=elem_field(ne_radius_out, &
                  elem_cnct(-1,1,ne))
          endif
          elem_field(ne_radius,ne)=(elem_field(ne_radius_in,ne)+ &
               elem_field(ne_radius_out,ne))/2.0_dp
       endif
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
       elem_field(ne_a_A,ne) = 1.0_dp ! set default for ratio a/A
    enddo

    call enter_exit(sub_name,2)
    
  end subroutine define_rad_from_file

!!!#############################################################################

  subroutine define_rad_from_geom(ORDER_SYSTEM, CONTROL_PARAM, START_FROM, &
       USER_RAD, group_type_in, group_option_in)
    !*define_rad_from_geom:* Defines vessel or airway radius based on
    ! their geometric structure. For 'order_system' == 'strah' or 'horsf', uses a
    ! user-defined maximum radius and branching ratio; for == 'fit', uses pre-
    ! defined radii (read in) and a calculated branching ratio for each path so
    ! that the order 1 branches have radius = USER_RAD.
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_GEOM" :: DEFINE_RAD_FROM_GEOM

    real(dp), intent(in) :: CONTROL_PARAM
    real(dp), intent(in) :: USER_RAD   ! radius of largest branch when order_system
    !                                    =='strah' or 'horsf'; minimum radius when
    !                                    order_system = 'fit'
    character(LEN=*), intent(in) :: ORDER_SYSTEM,START_FROM
    character(LEN=*), optional :: group_type_in, group_option_in
    !Input options ORDER_SYSTEM=STRAHLER (CONTROL_PARAM=RDS), HORSFIELD (CONTROL_PARAM=RDH)
    ! Local variables
    integer :: inlet_count,n,ne,ne0,ne_max,ne_min,ne_start,nindex,norder,n_max_ord
    real(dp) :: max_radius,radius,ratio_diameter
    logical :: found
    character(LEN=100) :: group_type, group_options
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'define_rad_from_geom'
    call enter_exit(sub_name,1)

    !define list of elements you are going to operate on
    if(present(group_type_in))then
       group_type = group_type_in
    else!default to all
       group_type='all'
    endif
    if(group_type.eq.'all')then
       ne_min=1
       ne_max=num_elems
    elseif(group_type.eq.'efield')then
    elseif(group_type.eq.'list')then
       read (START_FROM,'(I10)') ne_min
       read (group_option_in,'(I10)') ne_max
    endif
    !Define start element
    if(START_FROM.eq.'inlet')then
       inlet_count=0
       do ne=ne_min,ne_max
          if(elem_cnct(-1,0,ne).eq.0)then
             inlet_count=inlet_count+1
             ne_start=ne
          endif
          if(inlet_count.gt.1)then
             write(*,'('' More than one inlet in this group, using last found, ne = '',i6)') ne
          endif
       enddo
    else!element number defined
       read (START_FROM,'(I10)') ne_start
    endif
    
    ne=ne_start

    if(ORDER_SYSTEM(1:3).eq.'fit')then
       nindex = no_hord ! default is Horsfield ordering; could be modified to either type
       if(elem_field(ne_radius,1).lt.loose_tol)then
          write(*,*) 'Have you defined some upper radii??'
          read(*,*)
       endif
       do ne = ne_min,ne_max
          if(elem_field(ne_radius,ne).lt.USER_RAD)then
             found = .false.
             norder = elem_ordrs(nindex,ne)
             ne0 = elem_cnct(-1,1,ne)
             do while(.not.found)
                if(elem_field(ne_radius,ne0).gt.USER_RAD)then
                   found = .true.
                   n_max_ord = elem_ordrs(nindex,ne0)
                   max_radius = elem_field(ne_radius,ne0)
                   ratio_diameter = 10.d0**(log10(USER_RAD/max_radius) &
                        /dble(1-n_max_ord))
                   radius = (10.0_dp**(log10(ratio_diameter)*dble(norder- &
                        n_max_ord)+log10(2.0_dp*max_radius)))*0.5_dp
                   elem_field(ne_radius,ne) = radius
                   if(ne_vol.gt.0)then
                      elem_field(ne_vol,ne) = pi*radius**2*elem_field(ne_length,ne)
                   endif
                else
                   ne0 = elem_cnct(-1,1,ne0)
                endif
             enddo
          endif
       enddo
    
    else

       !Strahler and Horsfield ordering system
       if(ORDER_SYSTEM(1:5).eq.'strah')THEN
          nindex = no_sord !for Strahler ordering
          elem_field(ne_radius,ne) = USER_RAD
       else if(ORDER_SYSTEM(1:5).eq.'horsf')then
          nindex = no_hord !for Horsfield ordering
          elem_field(ne_radius,ne) = USER_RAD
       endif
       n_max_ord=elem_ordrs(nindex,ne)
    
       do ne=ne_min,ne_max
          radius = 10.0_dp**(log10(CONTROL_PARAM)*dble(elem_ordrs(nindex,ne) &
               -n_max_ord)+log10(USER_RAD))
          elem_field(ne_radius,ne)=radius
          if(ne_vol.gt.0)then
            elem_field(ne_vol,ne) = pi*radius**2*elem_field(ne_length,ne)
          endif
          if(ne_radius_in.gt.0)then
             elem_field(ne_radius_in,ne)=radius
             elem_field(ne_radius_out,ne)=radius
          endif
       enddo
    endif
    
    call enter_exit(sub_name,2)
    
  end subroutine define_rad_from_geom

!!!#############################################################################

  subroutine line_segments_for_2d_mesh(sf_option)
    !*line_segments_for_2d_mesh:* set up the line segment arrays for a 2d mesh
    
    character(len=4),intent(in) :: sf_option
    ! Local variables
    integer :: ne,ne_adjacent,ni1,nj,nl,nl_adj,npn(2)
    logical :: MAKE
    logical :: based_on_elems = .true.
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'line_segments_for_2d_mesh'
    call enter_exit(sub_name,1)

    ! allocate elem_lines_2d, scale_factors_2d,lines_2d,line_versn_2d,lines_in_elem,nodes_in_line,arclength
    if(allocated(elem_lines_2d)) deallocate(elem_lines_2d)
    if(allocated(scale_factors_2d)) deallocate(scale_factors_2d)
    allocate(elem_lines_2d(4,num_elems_2d))
    allocate(scale_factors_2d(16,num_elems_2d))
    !if(.not.allocated(elem_lines_2d)) allocate(elem_lines_2d(4,num_elems_2d))
    !if(.not.allocated(scale_factors_2d)) allocate(scale_factors_2d(16,num_elems_2d))
    
    elem_lines_2d=0
    num_lines_2d = 0
    
    if(based_on_elems)then
!!! estimate number of lines, for allocating memory to arrays
!!! before setting up arrays, count the number of lines required
       do ne=1,num_elems_2d
          MAKE=.FALSE.
          if(elem_cnct_2d(-1,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-1,1,ne)
          if(ne_adjacent > 0)then
             if(elem_lines_2d(4,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          if(MAKE) num_lines_2d = num_lines_2d+1
          MAKE=.FALSE.
          if(elem_cnct_2d(-2,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-2,1,ne)
          if(ne_adjacent > 0)then
             if(elem_lines_2d(2,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          if(MAKE) num_lines_2d=num_lines_2d+1
          num_lines_2d = num_lines_2d+2
          elem_lines_2d(2,ne) = 1 ! at this stage just to tag it for conditional above
          elem_lines_2d(4,ne) = 1 ! at this stage just to tag it for conditional above
       enddo !ne
       
       elem_lines_2d = 0
       
       if(allocated(lines_2d)) deallocate(lines_2d)
       if(allocated(line_versn_2d)) deallocate(line_versn_2d)
       if(allocated(lines_in_elem)) deallocate(lines_in_elem)
       if(allocated(nodes_in_line)) deallocate(nodes_in_line)
       if(allocated(arclength)) deallocate(arclength)
       allocate(lines_2d(0:num_lines_2d))
       allocate(line_versn_2d(2,3,num_lines_2d))
       allocate(lines_in_elem(0:4,num_lines_2d))
       allocate(nodes_in_line(3,0:3,num_lines_2d))
       allocate(arclength(num_lines_2d)) 

       !if(.not.allocated(lines_2d)) allocate(lines_2d(0:num_lines_2d))
       !if(.not.allocated(line_versn_2d)) allocate(line_versn_2d(2,3,num_lines_2d))
       !if(.not.allocated(lines_in_elem)) allocate(lines_in_elem(0:4,num_lines_2d))
       !if(.not.allocated(nodes_in_line)) allocate(nodes_in_line(3,0:3,num_lines_2d))
       !if(.not.allocated(arclength)) allocate(arclength(num_lines_2d)) 
       lines_in_elem=0
       lines_2d=0
       nodes_in_line=0
       line_versn_2d=0
       num_lines_2d = 0 ! reset to zero for loop below
       
!!! Now run through the same as above, and set up the arrays
       do ne=1,num_elems_2d
          !check whether to make a line
          MAKE=.FALSE.
          if(elem_cnct_2d(-1,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-1,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(4,ne_adjacent) == 0) MAKE=.TRUE.
          endif

          if(MAKE)then
             num_lines_2d = num_lines_2d+1
             lines_2d(num_lines_2d) = num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d) = lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d) = ne !line num_lines_2d is in element ne
             elem_lines_2d(3,ne) = num_lines_2d !num_lines_2d is global line # corresponding to local line 3 of ne
             npn(1) = 1
             npn(2) = 3
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(1,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(3,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent=elem_cnct_2d(-1,1,ne)
             elem_lines_2d(3,ne)=elem_lines_2d(4,ne_adjacent)
          endif
          
          !check whether to make a line
          MAKE=.FALSE.
          if(elem_cnct_2d(-2,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-2,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(2,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          
          if(MAKE)then
             num_lines_2d=num_lines_2d+1
             lines_2d(num_lines_2d)=num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
             elem_lines_2d(1,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 1 of ne
             npn(1)=1
             npn(2)=2
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(1,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent = elem_cnct_2d(-2,1,ne)
             do nl_adj = 1,4
                nl = elem_lines_2d(nl_adj,ne_adjacent)
                if(nl /= 0)then
                   if(nodes_in_line(2,1,nl) == elem_nodes_2d(1,ne) .and. &
                        nodes_in_line(3,1,nl) == elem_nodes_2d(2,ne))then
                      elem_lines_2d(1,ne) = nl
                   elseif(nodes_in_line(2,1,nl) == elem_nodes_2d(2,ne) .and. &
                        nodes_in_line(3,1,nl) == elem_nodes_2d(1,ne))then
                      elem_lines_2d(1,ne) = nl
                   endif
                endif
             enddo
             !             elem_lines_2d(1,ne)=elem_lines_2d(2,ne_adjacent)
          endif
          
          MAKE=.TRUE.
          ne_adjacent=elem_cnct_2d(1,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(3,ne_adjacent) /= 0) MAKE=.FALSE.
          endif
          
          if(MAKE)then
             num_lines_2d=num_lines_2d+1
             lines_2d(num_lines_2d)=num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
             elem_lines_2d(4,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 4 of ne
             npn(1)=2
             npn(2)=4
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent=elem_cnct_2d(1,1,ne)
             elem_lines_2d(4,ne)=elem_lines_2d(3,ne_adjacent)
          endif
          
          MAKE=.TRUE.
          ne_adjacent=elem_cnct_2d(2,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(1,ne_adjacent) /= 0) MAKE=.FALSE.
          endif
          
          if(MAKE)then
             num_lines_2d = num_lines_2d+1
             lines_2d(num_lines_2d) = num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d) = lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d) = ne !line num_lines_2d is in element ne
             elem_lines_2d(2,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 2 of ne
             npn(1) = 3
             npn(2) = 4
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(3,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent=elem_cnct_2d(2,1,ne)
             elem_lines_2d(2,ne)=elem_lines_2d(1,ne_adjacent)
          endif
       enddo !ne
    endif
    
    call calc_scale_factors_2d(sf_option)
    
    call enter_exit(sub_name,2)
    
  end subroutine line_segments_for_2d_mesh

!!!#############################################################################

  subroutine set_initial_volume(Gdirn,COV,total_volume,Rmax,Rmin)
    !*set_initial_volume:* assigns a volume to terminal units appended on a
    ! tree structure based on an assumption of a linear gradient in the
    ! gravitational direction with max, min, and COV values defined.
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_INITIAL_VOLUME" :: SET_INITIAL_VOLUME
    
    integer,intent(in) :: Gdirn
    real(dp),intent(in) :: COV,total_volume,Rmax,Rmin
    !     Local parameters
    integer :: ne,np2,nunit
    real(dp) ::  factor_adjust,max_z,min_z,random_number,range_z,&
         volume_estimate,volume_of_tree,Vmax,Vmin,Xi
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'set_initial_volume'
    call enter_exit(sub_name,1)
    
    volume_estimate = 1.0_dp
    volume_of_tree = 0.0_dp
    
    call volume_of_mesh(volume_estimate,volume_of_tree)
    
    random_number=-1.1_dp
    
    Vmax = Rmax * (total_volume-volume_estimate)/elem_units_below(1)
    Vmin = Rmin * (total_volume-volume_estimate)/elem_units_below(1)
    
!!! for each elastic unit find the maximum and minimum coordinates in the Gdirn direction
    max_z=-1.0e+6_dp
    min_z=1.0e+6_dp
    do nunit=1,num_units
       ne=units(nunit)
       np2=elem_nodes(2,ne)
       max_z=MAX(max_z,node_xyz(Gdirn,np2))
       min_z=MIN(min_z,node_xyz(Gdirn,np2))
    enddo !nunit
    
    range_z=abs(max_z-min_z)
    if(abs(range_z).le.1.0e-5_dp) range_z=1.0_dp
    
!!! for each elastic unit allocate a size based on a gradient in the Gdirn direction, and
!!! perturb by a user-defined COV. This should be calling a random number generator.
    do nunit=1,num_units
       ne=units(nunit)
       np2=elem_nodes(2,ne) !end node
       Xi=(node_xyz(Gdirn,np2)-min_z)/range_z
       random_number=random_number+0.1_dp
       if(random_number.GT.1.0_dp) random_number=-1.1_dp
       unit_field(nu_vol,nunit)=(Vmax*Xi+Vmin*(1.0_dp-Xi))*(1.0_dp+COV*random_number)
       unit_field(nu_vt,nunit)=0.0_dp !initialise the tidal volume to a unit
    enddo !nunit
    
    ! correct unit volumes such that total volume is exactly as specified
    call volume_of_mesh(volume_estimate,volume_of_tree)
    factor_adjust = (total_volume-volume_of_tree)/(volume_estimate-volume_of_tree)
    do nunit=1,num_units
       unit_field(nu_vol,nunit) = unit_field(nu_vol,nunit)*factor_adjust
    enddo
    
    write(*,'('' Number of elements is '',I5)') num_elems
    write(*,'('' Initial volume is '',F6.2,'' L'')') total_volume/1.0e+6_dp
    write(*,'('' Deadspace volume is '',F6.1,'' mL'')') volume_of_tree/1.0e+3_dp
    
    call enter_exit(sub_name,2)
    
  end subroutine set_initial_volume

!!!#############################################################################

  subroutine volume_of_mesh(volume_model,volume_tree)
    !*volume_of_mesh:* calculates the volume of an airway mesh including
    ! conducting and respiratory airways
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VOLUME_OF_MESH" :: VOLUME_OF_MESH
    
    real(dp) :: volume_model,volume_tree
    !     Local Variables
    integer :: ne,ne0,nunit
    real(dp),allocatable :: vol_anat(:),vol_below(:)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'volume_of_mesh'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(vol_anat)) allocate(vol_anat(num_elems))
    if(.not.allocated(vol_below)) allocate(vol_below(num_elems))
    
    vol_anat = elem_field(ne_vol,1:num_elems) !initialise to branch volume
    vol_below = elem_field(ne_vol,1:num_elems) !initialise to branch volume

    do nunit = 1,num_units
       ne = units(nunit)
       if(ne.ne.0) vol_below(ne) = vol_below(ne) + unit_field(nu_vol,nunit) !add elastic unit volume
    enddo !nunit

    do ne = num_elems,2,-1
       ne0 = elem_cnct(-1,1,ne)
!!! don't include the respiratory airways (already included via lumped units). multiply
!!! by the element type (0 for respiratory) to account for this
       vol_anat(ne0) = vol_anat(ne0) + dble(elem_symmetry(ne))*dble(elem_ordrs(no_type,ne))*vol_anat(ne)
       vol_below(ne0) = vol_below(ne0) + dble(elem_symmetry(ne))*dble(elem_ordrs(no_type,ne))*vol_below(ne)
    enddo !noelem

    elem_field(ne_vd_bel,:) = vol_anat(:)
    elem_field(ne_vol_bel,:) = vol_below(:)
    volume_model = elem_field(ne_vol_bel,1)
    volume_tree = elem_field(ne_vd_bel,1)

    deallocate(vol_anat)
    deallocate(vol_below)
    
    call enter_exit(sub_name,2)
    
  end subroutine volume_of_mesh

!!!#############################################################################

  subroutine write_geo_file(type, filename)
    !*write_geo_file:* converts a surface mesh (created using make_2d_vessel_from_1d)
    ! into a gmsh formatted mesh and writes to file. 
    ! options on 'type': 1== single layered surface mesh of the vessel wall
    !                    2== double-layered thick-walled volume mesh of vessel wall
    !                    3== volume mesh of vessel lumen
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_WRITE_GEO_FILE" :: WRITE_GEO_FILE

    integer,intent(in) :: type
    character(len=*),intent(in) :: filename
    !     Local parameters
    integer :: j, ncount_loop = 0, ncount_point = 0, ncount_spline = 0, &
         nl_offset,np,np_offset
    integer,parameter :: ifile = 10
    integer,allocatable :: element_spline(:,:),elem_surfaces(:,:)
    real(dp),parameter :: lc0 = 1.0_dp, lc1 = 1.0_dp
    real(dp),allocatable :: node_xyz_offset(:,:)
    character(len=200) :: opfile
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'write_geo_file'
    call enter_exit(sub_name,1)

    opfile = trim(filename)//'.geo'
    open(10, file=opfile, status='replace')
      
    write(ifile,'(''/***********************'')')
    write(ifile,'(''*'')')
    write(ifile,'(''* Conversion of LungSim to GMSH'')')
    write(ifile,'(''*'')')
    write(ifile,'(''***********************/'')')
    
    write(ifile,'(/''lc ='',f8.4,'';'')') lc0
    write(ifile,'(/''sc ='',f8.4,'';'')') lc1
    write(ifile,'(/)')

    allocate(element_spline(4,num_elems_2d*2))
    allocate(elem_surfaces(5,num_elems_2d))
    element_spline = 0
    elem_surfaces = 0
    ncount_spline = 0 
    np_offset = 0

    if(type.eq.1)then
!!! write out a surface mesh that describes a structured vessel surface
       call write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline,np_offset)
       
    else if(type.eq.2)then
!!! write out a volume that encloses a thick-walled vessel tree. Make a gmsh .geo file
!!! for the surface of the tree, then copy, scale, and translate to get an 'outer shell'.
!!! Join the inner and outer shells at the entry and exits.

       allocate(node_xyz_offset(3,num_nodes_2d))
       node_xyz_offset = 0.0_dp
       call write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline,np_offset)
       call geo_node_offset(node_xyz_offset)

       do np = 1,num_nodes_2d
          forall (j = 1:3) node_xyz_2d(1,1,j,np) = node_xyz_2d(1,1,j,np) &
               + node_xyz_offset(j,np)
       enddo
       np_offset = ncount_point
       nl_offset = ncount_spline
       call write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline,np_offset)
       do np = 1,num_nodes_2d
          forall (j = 1:3) node_xyz_2d(1,1,j,np) = node_xyz_2d(1,1,j,np) &
               - node_xyz_offset(j,np)
       enddo
       ! cap the entry and exits
       call geo_entry_exit_cap(element_spline,ifile,ncount_loop, &
            ncount_spline,np_offset,nl_offset)
       deallocate(node_xyz_offset)

    else if(type.eq.3)then
!!! write out a volume mesh for the vessel lumen, where the vessel surface mesh is the
!!! exterior. Make a .gmsh file that includes the vessel surfaces, surfaces that join to a vessel
!!! centreline, and surfaces that 'cap' each vessel segment entry and exit.

       call write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline,np_offset)
       call write_3d_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline)
    endif

    deallocate(element_spline)
    deallocate(elem_surfaces)
    close(ifile)

    call enter_exit(sub_name,2)
    
  end subroutine write_geo_file
  
!!!#############################################################################
  
  function get_final_real(string)
    !*get_final_real:* gets the last real number on a string

    character,intent(in) :: string*(*)
    ! Local parameters
    integer :: ibeg,iend
    real(dp) :: rsign,rtemp,get_final_real
    character :: sub_string*(40)
    
    ! --------------------------------------------------------------------------
    
    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
   
    read (sub_string(ibeg:iend), * ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number
    
    get_final_real=rtemp !return the real value
    
  end function get_final_real

!!!#############################################################################

  subroutine get_final_string(string,rtemp)
    !*get_final_string:* gets the last set of characters surrounded
    !  by whitespace in a string

    character, intent(in) :: string*(*)
    ! Local parameters
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'get_final_string'
    call enter_exit(sub_name,1)
    
    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
    read (sub_string(ibeg:iend), '(D25.17)' ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number
    
    call enter_exit(sub_name,2)

  end subroutine get_final_string

!!!#############################################################################

  subroutine get_local_node(np_global,np_local)
    !*get_local_node:* gets the local node number for a given global node

    integer,intent(in) :: np_global
    integer,intent(out) :: np_local
    ! Local parameters
    integer :: np
    logical :: found
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'get_local_node'
    call enter_exit(sub_name,1)
    
    np=1
    found=.false.
    do while (.not.found)
       if(nodes(np).eq.np_global)then
          found=.true.
       elseif(np.gt.num_nodes)then
          found = .true.
          write(*,'('' Global node '',I6,'' not in node list'')') np_global
          read(*,*)
       else
          np=np+1
       endif
    enddo
    
    np_local = np

    call enter_exit(sub_name,2)

  end subroutine get_local_node
  
!!!#############################################################################

  subroutine geo_entry_exit_cap(element_spline,ifile,ncount_loop, &
       ncount_spline,np_offset,nl_offset)

    integer,intent(in) :: element_spline(:,:),ifile,np_offset,nl_offset
    integer :: ncount_loop,ncount_spline
    ! Local variables
    integer :: k,line1,line2,line3,line4,ne,np1,np2
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'geo_entry_exit_cap'
    call enter_exit(sub_name,1)
        
    ne = 1
    do while (ne.le.num_elems_2d)
       
       if(elem_cnct_2d(-2,0,ne).eq.0)then
          do k = 0,3
             np1 = elem_nodes_2d(1,ne+k)
             np2 = np1 + np_offset
             ncount_spline = ncount_spline + 1
             write(10,'(''Line('',I8,'') = {'',I8,'','',I8,''};'')') &
                  ncount_spline,np1,np2
          enddo
          do k = 0,3
             line1 = element_spline(1,ne+k) - nl_offset
             line3 = -element_spline(1,ne+k)
             if(k.lt.3)then
                line2 = ncount_spline + k - 2
                line4 = -(line2 - 1)
             else
                line2 = ncount_spline - 3 ! first new line
                line4 = -ncount_spline
             endif
             ncount_loop = ncount_loop + 1
             write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8, &
                  &'','',i8,'','',i8,''};'')') &
                  ncount_loop, line1, line2, line3, line4
             ncount_loop = ncount_loop + 1
             write(ifile,'(''Surface('',I8,'') = {'',I8,''};'')') &
                  ncount_loop, ncount_loop - 1
          enddo
       endif
       
       if(elem_cnct_2d(2,0,ne).eq.0)then
          do k = 0,3
             np1 = elem_nodes_2d(3,ne+k)
             np2 = np1 + np_offset
             ncount_spline = ncount_spline + 1
             write(10,'(''Line('',I8,'') = {'',I8,'','',I8,''};'')') &
                  ncount_spline,np1,np2
          enddo
          do k = 0,3
             line1 = -element_spline(3,ne+k) - nl_offset
             line3 = -element_spline(3,ne+k)
             if(k.lt.3)then
                line2 = ncount_spline + k - 2
                line4 = -(line2 - 1)
             else
                line2 = ncount_spline - 3 ! first new line
                line4 = -ncount_spline
             endif
             ncount_loop = ncount_loop + 1
             write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','' &
                  &,i8,'','',i8,''};'')') &
                  ncount_loop, line1, line2, line3, line4
             ncount_loop = ncount_loop + 1
             write(ifile,'(''Surface('',I8,'') = {'',I8,''};'')') &
                  ncount_loop, ncount_loop - 1
          enddo
       endif
       
       ne = ne + 4
    enddo

    call enter_exit(sub_name,2)
    
  end subroutine geo_entry_exit_cap

!!!#############################################################################

  subroutine geo_node_offset(node_xyz_offset)
    
    real(dp) :: node_xyz_offset(:,:)
    ! Local variables
    integer:: j,k,ne,np1
    real(dp) :: point_temp(3),point_xyz_centre(3),wall_thickness
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'geo_node_offset'
    call enter_exit(sub_name,1)
        
    ne = 1
    do while (ne.le.num_elems_2d)
       
!!!....nodes at model entry
       if(elem_cnct_2d(-2,0,ne).eq.0)then
          point_xyz_centre(:) = 0.0_dp
          do k = 0,3
             np1 = elem_nodes_2d(1,ne+k)
             forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) + &
                  0.25_dp * node_xyz_2d(1,1,j,np1)
          enddo
          point_temp(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(1,ne)) ! location of first ring node
          wall_thickness = &
               0.2_dp * distance_between_points(point_xyz_centre,point_temp)
          do k = 0,3
             np1 = elem_nodes_2d(1,ne+k)
             point_temp(1:3) = node_xyz_2d(1,1,1:3,np1)
             node_xyz_offset(:,np1) = wall_thickness * &
                  direction_point_to_point(point_xyz_centre,point_temp) 
          enddo  ! k
       endif  ! elem_cnct
       
!!!....nodes at Xi2=1 ends of 'rings'
       point_xyz_centre(:) = 0.0_dp
       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) &
               + 0.25_dp * node_xyz_2d(1,1,j,np1)
       enddo
       point_temp(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,ne)) ! location of first ring node
       wall_thickness = 0.2_dp * &
            distance_between_points(point_xyz_centre,point_temp)
       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          point_temp(1:3) = node_xyz_2d(1,1,1:3,np1)
          node_xyz_offset(:,np1) = wall_thickness * &
               direction_point_to_point(point_xyz_centre,point_temp)
       enddo  ! k

!!!....check for crux node
       if(node_versn_2d(elem_nodes_2d(3,ne)).eq.6.or. &
            node_versn_2d(elem_nodes_2d(4,ne)).eq.6)then
          np1 = elem_nodes_2d(3,ne+3) + 1   ! number of crux node
          point_temp(1:3) = node_xyz_2d(1,1,1:3,np1)
          node_xyz_offset(:,np1) =  wall_thickness * &
               direction_point_to_point(point_xyz_centre,point_temp)
       endif

       ne = ne + 4
       
    enddo ! while (ne.le.num_elems_2d)
    
    call enter_exit(sub_name,2)
    
  end subroutine geo_node_offset
  
!!!#############################################################################

  subroutine enclosed_volume(surface_elems)
    !*enclosed_volume:* estimates the volume that is inside a list of bounding
    ! surface elements
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ENCLOSED_VOLUME" :: ENCLOSED_VOLUME

    use mesh_utilities,only: triangles_from_surface,volume_internal_to_surface

    integer,intent(in) :: surface_elems(:)
    ! Local variables
    integer :: i,num_triangles,num_vertices
    integer,allocatable :: elem_list(:)
    integer,allocatable :: triangle(:,:)
    real(dp) :: volume
    real(dp),allocatable :: vertex_xyz(:,:)

    allocate(elem_list(count(surface_elems.ne.0)))
    do i = 1,count(surface_elems.ne.0)
       elem_list(i) = get_local_elem_2d(surface_elems(i))
    enddo

    call triangles_from_surface(num_triangles,num_vertices,elem_list, &
       triangle,vertex_xyz)
    volume = volume_internal_to_surface(triangle,vertex_xyz)
    write(*,'('' Enclosed volume = '',f9.2,'' mm^3'')') volume

    deallocate(elem_list)
    deallocate(triangle)
    deallocate(vertex_xyz)
    
  end subroutine enclosed_volume

!!!#############################################################################
  
  function get_final_integer(string)
    !*get_final_integer*
    
    character,intent(in) :: string*(*)
    ! Local parameters
    integer :: ibeg,iend,ierror,nsign,ntemp
    character :: sub_string*(40)
    integer :: get_final_integer
    
    ! --------------------------------------------------------------------------
    
    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of integer in string, follows ":"
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond ":"
    iend=len(sub_string) !length of the sub-string
    if(sub_string(1:1).eq.'-')then !check for negative sign
       nsign=-1
       ibeg=2
    else
       nsign=1
       ibeg=1
    endif
    
    read (sub_string(ibeg:iend), '(i10)', iostat=ierror ) ntemp !get integer values
    if(ierror.gt.0)then
       !... something wrong with data
       write(*,'(''Data read error'')')
       write(*,'(a)') sub_string(ibeg:iend)
    endif
    ntemp=ntemp*nsign !apply sign to number
    
    get_final_integer=ntemp !return the integer value
    
  end function get_final_integer
  
!!!#############################################################################
  
  subroutine get_four_nodes(ne,string)

    integer, intent(in) :: ne
    character(len=132), intent(in) :: string
    ! Local variables
    integer :: ibeg,iend,i_ss_end,nn,np_global
    character(len=40) :: sub_string
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'get_four_nodes'
    call enter_exit(sub_name,1)

    iend=len(string)
    ibeg=index(string,":")+1 !get location of first integer in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond : and remove the leading blanks
    i_ss_end=len(sub_string) !get the end location of the sub-string
    
    ibeg=1
    do nn=1,4
       iend=index(sub_string," ") !get location of first blank in sub-string
       read (sub_string(ibeg:iend-1), '(i7)' ) np_global
       sub_string = adjustl(sub_string(iend:i_ss_end)) ! get chars beyond blank, remove leading blanks
       elem_nodes_2d(nn,ne)=get_local_node_f(2,np_global)
    enddo ! nn

    call enter_exit(sub_name,2)
    
  end subroutine get_four_nodes
  
!!!#############################################################################

  subroutine redistribute_mesh_nodes_2d_from_1d_0

    integer :: i,j,k,ne,nelist(20),ne_adjacent,np,nplist(20),np_adjacent,np_last,num_list, &
         ring1_nodes(4)
    real(dp) :: centre(3),displace_length,distance_to_crux,distance_to_crux_last,line_length, &
         nedirection(3,20),point1(3),point2(3),point3(3),point4(3),theta,vector(3)
    logical :: continue
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'redistribute_mesh_nodes_2d_from_1d'
    call enter_exit(sub_name,1)

    nplist = 0
    nelist = 0

    np = 1
    do while (np.le.num_nodes_2d) ! run through all of the 2d mesh nodes
       ! use the 'front' bifurcation node to identify the crux node (based on the
       ! template structure that we used to set up the bifurcations). Crux node
       ! must be at front_node + 3.
       if(node_versn_2d(np).eq.6)then   ! this is the 'front' node at a bifurcation
          np = np + 3   ! this is the node number for the 'crux' node
          do i = 1,2 ! for each of the two branches that the bifurcation leads to
             ne = elems_at_node_2d(np,2*i-1)   ! get the first (and then third) element that node np (crux) is in
             num_list = 1                      ! count the number of 'rings' of nodes between bifurcations
             nelist(num_list) = ne             ! record the first (and third) element number ne
             np_adjacent = elem_nodes_2d(4,ne) ! node number in the +Xi2 direction (along the branch)
             nplist(num_list) = np_adjacent    ! the node number along the line from one bifurcation to the next
             ! get coordinates for three of the points on the first 'ring' that is in the direction of the
             ! branch. Not straightforward when at a bifurcation.
             point1(1:3) = node_xyz_2d(1,1,1:3,np_adjacent)
             point2(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,ne))
             point3(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,elem_cnct_2d(-1,1,ne)))
             point4(1:3) = node_xyz_2d(1,1,1:3,np)   ! the location of the crux node
             ! calculate the distance from the crux (point4) to the plane of first adjacent ring
             line_length = distance_from_plane_to_point(point1,point2,point3,point4)
             
             !point1(1:3) = node_xyz_2d(1,1,1:3,np-1)   ! the location of the 'back' node of bifurcation
!!             ! calculate the line length from back node to a point on the first ring
!!             distance_to_crux_last = distance_between_points(point1,point2)
             
             continue = .true.
             do while(continue)   ! keep going: note that bifurcation will have > 1 version at nodes
                if(elem_cnct_2d(2,0,ne).eq.0)then ! no adjacent 2d elements in Xi+2 direction
                   continue = .false.
                else
                   ne = elem_cnct_2d(2,1,ne)        ! get the next adjacent element in Xi+2 direction
                   num_list = num_list + 1 ! the number of 'rings' between bifurcations
                   nelist(num_list) = ne ! the element number of the adjacent element
                   np_last = np_adjacent   ! store the previous node number
                   np_adjacent = elem_nodes_2d(4,ne) ! the next node on the line
                   nplist(num_list) = np_adjacent    ! the node number on the line from one bifurcation to the next
                   ! calculate the distance between adjacent rings
                   point1(1:3) = node_xyz_2d(1,1,1:3,np_last)      ! coordinates of node on previous ring
                   point2(1:3) = node_xyz_2d(1,1,1:3,np_adjacent)  ! coordinate of node on current ring
                   line_length = line_length + &
                        distance_between_points(point1,point2)  ! sum distance between rings
                   ! calculate the direction between connected nodes on rings
                   vector(1:3) = node_xyz_2d(1,1,1:3,np_adjacent) - &
                        node_xyz_2d(1,1,1:3,np_last)
                   vector = unit_vector(vector)
                   nedirection(1:3,num_list-1) = vector(1:3)  ! store the direction
                   ! continue until the next bifurcation is detected (nodes with > 1 version)                   
                   if(node_versn_2d(np_adjacent).ne.1) continue = .false.
                endif
             enddo

             line_length = line_length/real(num_list) ! this is the length to redistribute rings to
             
!!!          adjust the location of the nodes in each 'ring'
             do j = 1,num_list - 1   ! only adjust the rings that are between bifns: last 'ring' is actually the next bifn
                
                ! first get the list of nodes in the ring
                ring1_nodes(1) = nplist(j)
                ne_adjacent = elem_cnct_2d(1,1,nelist(j)) ! get the next element in the +Xi1 direction
                do k = 2,4
                   ring1_nodes(k) = elem_nodes_2d(4,ne_adjacent)
                   ne_adjacent = elem_cnct_2d(1,1,ne_adjacent) ! get the next element in the +Xi1 direction
                enddo ! k
                
                ! assume that the direction for adjustment is defined by location of adjacent rings
                vector(1:3) = nedirection(1:3,j)
                
                ! calculate the ring displacement = j*line_length - (distance from crux)
                point1(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(1))
                point2(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(2))
                point3(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(3))
                point4(1:3) = node_xyz_2d(1,1,1:3,np)
                displace_length = real(j) * line_length - &
                     distance_from_plane_to_point(point1,point2,point3,point4)
                
                ! update the location of the four nodes in the current ring
                do k = 1,4
                   node_xyz_2d(1,1,1:3,ring1_nodes(k)) = &
                        node_xyz_2d(1,1,1:3,ring1_nodes(k)) + &
                        vector(1:3) * displace_length
                enddo

             enddo ! j
          enddo ! i
       endif
       np = np + 1 ! increment to check the next node
    enddo
    
    call enter_exit(sub_name,2)
    
  end subroutine redistribute_mesh_nodes_2d_from_1d_0

!!!#############################################################################

  subroutine redistribute_mesh_nodes_2d_from_1d

    integer :: i,j,k,ne,nelist(20),ne_adjacent,np,nplist(20),np_adjacent,np_last,num_list, &
         ring1_nodes(4)
    real(dp) :: centre(3),displace_length,distance_to_crux,distance_to_crux_last,line_length, &
         nedirection(3,20),point1(3),point2(3),point3(3),point4(3),ring_dist(10),theta,vector(3)
    logical :: continue
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'redistribute_mesh_nodes_2d_from_1d'
    call enter_exit(sub_name,1)

    nplist = 0
    nelist = 0

    np = 1
    do while (np.le.num_nodes_2d) ! run through all of the 2d mesh nodes
       ! use the 'front' bifurcation node to identify the crux node (based on the
       ! template structure that we used to set up the bifurcations). Crux node
       ! must be at front_node + 3.
       if(node_versn_2d(np).eq.6)then   ! this is the 'front' node at a bifurcation
          np = np + 3   ! this is the node number for the 'crux' node
          do i = 1,2 ! for each of the two branches that the bifurcation leads to
             ne = elems_at_node_2d(np,2*i-1)   ! get the first (and then third) element that node np (crux) is in
             num_list = 1                      ! count the number of 'rings' of nodes between bifurcations
             nelist(num_list) = ne             ! record the first (and third) element number ne
             np_adjacent = elem_nodes_2d(4,ne) ! node number in the +Xi2 direction (along the branch)
             nplist(num_list) = np_adjacent    ! the node number along the line from one bifurcation to the next
             ! get coordinates for three of the points on the first 'ring' that is in the direction of the
             ! branch. Not straightforward when at a bifurcation.
             point1(1:3) = node_xyz_2d(1,1,1:3,np_adjacent)
             point2(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,ne))
             point3(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,elem_cnct_2d(-1,1,ne)))
             point4(1:3) = node_xyz_2d(1,1,1:3,np)   ! the location of the crux node
             ! calculate the distance from the crux (point4) to the plane of first adjacent ring
             line_length = distance_from_plane_to_point(point1,point2,point3,point4)

             ! 2nd side node is at np-2 and 1st side node is at np-4
             if(i.eq.1)then
                point1(1:3) = node_xyz_2d(1,1,1:3,np-4)
                point2(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,ne))
             elseif(i.eq.2)then
                point1(1:3) = node_xyz_2d(1,1,1:3,np-2)
                point2(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(4,ne))
             endif
             line_length = distance_between_points(point1,point2)
             ring_dist(num_list) = line_length
             
             continue = .true.
             do while(continue)   ! keep going: note that bifurcation will have > 1 version at nodes
                if(elem_cnct_2d(2,0,ne).eq.0)then ! no adjacent 2d elements in Xi+2 direction
                   continue = .false.
                else
                   ne = elem_cnct_2d(2,1,ne)        ! get the next adjacent element in Xi+2 direction
                   num_list = num_list + 1 ! the number of 'rings' between bifurcations
                   nelist(num_list) = ne ! the element number of the adjacent element
                   np_last = np_adjacent   ! store the previous node number
                   np_adjacent = elem_nodes_2d(4,ne) ! the next node on the line
                   nplist(num_list) = np_adjacent    ! the node number on the line from one bifurcation to the next
                   ! calculate the distance between adjacent rings
                   point1(1:3) = node_xyz_2d(1,1,1:3,np_last)      ! coordinates of node on previous ring
                   point2(1:3) = node_xyz_2d(1,1,1:3,np_adjacent)  ! coordinate of node on current ring
                   line_length = line_length + &
                        distance_between_points(point1,point2)  ! sum distance between rings
                   ! calculate the direction between connected nodes on rings
                   vector(1:3) = node_xyz_2d(1,1,1:3,np_adjacent) - &
                        node_xyz_2d(1,1,1:3,np_last)
                   vector = unit_vector(vector)
                   nedirection(1:3,num_list-1) = vector(1:3)  ! store the direction
                   ! continue until the next bifurcation is detected (nodes with > 1 version)                   
                   if(node_versn_2d(np_adjacent).ne.1) continue = .false.
                   ring_dist(num_list) = line_length
                endif
             enddo

             line_length = line_length/real(num_list) ! this is the length to redistribute rings to
             
!!!          adjust the location of the nodes in each 'ring'
             do j = 1,num_list - 1   ! only adjust the rings that are between bifns: last 'ring' is actually the next bifn
                
                ! first get the list of nodes in the ring
                ring1_nodes(1) = nplist(j)
                ne_adjacent = elem_cnct_2d(1,1,nelist(j)) ! get the next element in the +Xi1 direction
                do k = 2,4
                   ring1_nodes(k) = elem_nodes_2d(4,ne_adjacent)
                   ne_adjacent = elem_cnct_2d(1,1,ne_adjacent) ! get the next element in the +Xi1 direction
                enddo ! k
                
                ! assume that the direction for adjustment is defined by location of adjacent rings
                vector(1:3) = nedirection(1:3,j)
                
                ! calculate the ring displacement = j*line_length - (distance from crux)
                !point1(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(1))
                !point2(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(2))
                !point3(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(3))
                !point4(1:3) = node_xyz_2d(1,1,1:3,np)
                !displace_length = real(j) * line_length - &
                !     distance_from_plane_to_point(point1,point2,point3,point4)
                displace_length = real(j) * line_length - ring_dist(j)
                
                ! update the location of the four nodes in the current ring
                do k = 1,4
                   node_xyz_2d(1,1,1:3,ring1_nodes(k)) = &
                        node_xyz_2d(1,1,1:3,ring1_nodes(k)) + &
                        vector(1:3) * displace_length
                enddo

!                ! rotate the ring a bit to test
!                theta = -pi*30.0_dp/180.0_dp ! 15 degrees
!                point4(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(4))
!                centre = (point1 + point2 + point3 + point4)/4.0_dp
!                do k = 1,4
!                   call mesh_rotate_about_axis_basic(ring1_nodes(k),vector,centre,theta)
!                enddo
                
             enddo ! j
          enddo ! i
       endif
       np = np + 1 ! increment to check the next node
    enddo
    
    call enter_exit(sub_name,2)
    
  end subroutine redistribute_mesh_nodes_2d_from_1d

!!!#############################################################################

  function get_local_elem_1d(ne_global)

    integer,intent(in) :: ne_global
    ! Local variables
    integer :: ne
    integer :: get_local_elem_1d

    ! --------------------------------------------------------------------------

    get_local_elem_1d = 0
    do ne=1,num_elems
       if(ne_global.eq.elems(ne)) get_local_elem_1d = ne
    enddo

  end function get_local_elem_1d

!!!###########################################################################

  subroutine write_elem_geometry_2d(elemfile)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_WRITE_ELEM_GEOMETRY_2D" :: WRITE_ELEM_GEOMETRY_2D

    character(len=*),intent(in) :: elemfile
    !     Local Variables
    integer :: ne,ne_count,nglobal_list(4),np,nv
    character(len=132) :: writefile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'write_elem_geometry_2d'
    call enter_exit(sub_name,1)

    writefile = trim(elemfile)//'.ipelem'
    open(10, file=writefile, status='replace')
    
    !.....write the total number of elems
    write(10,'('' CMISS Version 2.1  ipelem File Version 2'')')
    write(10,'('' Heading: 2D surface from 1D centreline'')')
    write(10,'()')
    write(10,'('' The number of elements is [1]:  '',i6)') num_elems_2d
    write(10,'()')

    !    do ne = 1,num_elems_2d
    ne_count = 1
    ne = 0
    do while (ne_count.le.num_elems_2d)
       ne = ne + 1
       if(elems_2d(ne).gt.0)then
          ne_count = ne_count + 1
          write(10,'('' Element number [    1]:  '',i6)')   elems_2d(ne)
          write(10,'('' The number of geometric Xj-coordinates is [3]: 3'')')
          write(10,'('' The basis function type for geometric variable 1 is [1]:  1'')')
          write(10,'('' The basis function type for geometric variable 2 is [1]:  1'')')
          write(10,'('' The basis function type for geometric variable 3 is [1]:  1'')')
          do np = 1,4
             nglobal_list(np) = nodes_2d(elem_nodes_2d(np,ne))
          enddo
          write(10,'('' Enter the 4 global numbers for basis 1: '',4(i6))') &
               nglobal_list(:)
          do np = 1,4
             if(node_versn_2d(elem_nodes_2d(np,ne)).gt.1)then ! has versions
                nv = elem_versn_2d(np,ne)
                write(10,'('' The version number for occurrence  1 of node'' &
                     &,i7,'', njj=1 is [ 1]: '',i3)') nglobal_list(np), nv
                write(10,'('' The version number for occurrence  1 of node'' &
                     &,i7,'', njj=2 is [ 1]: '',i3)') nglobal_list(np), nv
                write(10,'('' The version number for occurrence  1 of node'' &
                     &,i7,'', njj=3 is [ 1]: '',i3)') nglobal_list(np), nv
             endif
          enddo !np
          write(10,'()')
       endif
    enddo ! ne

  end subroutine write_elem_geometry_2d

!!!#############################################################################

  subroutine write_node_geometry_2d(NODEFILE)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_WRITE_NODE_GEOMETRY_2D" :: WRITE_NODE_GEOMETRY_2D

    character(len=*),intent(in) :: NODEFILE
    !     Local Variables
    integer :: i,np,np_count,nv
    character(len=132) :: writefile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'write_node_geometry_2d'
    call enter_exit(sub_name,1)

    writefile = trim(nodefile)//'.ipnode'
    open(10, file=writefile, status='replace')
    
    !.....write the total number of nodes
    write(10,'('' CMISS Version 1.21 ipnode File Version 2'')')
    write(10,'('' Heading: '')')
    write(10,'()')
    write(10,'('' The number of nodes is [    1]:  '',i6)') num_nodes_2d
    write(10,'('' Number of coordinates [ 3]:  3'')')
    do i=1,3
       write(10,'('' Do you want prompting for different versions of nj='',i1,'' [N]? Y'')') i
    enddo
    do i=1,3
       write(10,'('' The number of derivatives for coordinate '',i1,'' is [0]: 3'')') i
    enddo

    !    do np = 1,num_nodes_2d
    np_count = 1
    np = 0
    do while (np_count.le.num_nodes_2d)
       np = np + 1
       if(nodes_2d(np).gt.0)then
          np_count = np_count + 1
          write(10,'()')
          write(10,'('' Node number [    1]: '',i6)')  nodes_2d(np)
          do i=1,3
             write(10,'('' The number of versions for nj='',i1,'' is [1]:'',i2)')  i,node_versn_2d(np)
             do nv=1,node_versn_2d(np)
                if(node_versn_2d(np).gt.1) write(10,'('' For version number '',i1,'':'')') nv 
                !...........coordinate          
                write(10,'('' The Xj('',i1,'') coordinate is [ 0.00000E+00]: '',f12.5)') &
                     i,node_xyz_2d(1,nv,i,np)
                write(10,'('' The derivative wrt direction 1 is [ 0.00000E+00]: '',f12.5)') &
                     node_xyz_2d(2,nv,i,np)
                write(10,'('' The derivative wrt direction 2 is [ 0.00000E+00]: '',f12.5)') &
                     node_xyz_2d(3,nv,i,np)
                write(10,'('' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: '',f12.5)') &
                     node_xyz_2d(4,nv,i,np)
             enddo !nv
          end do !i
       endif
    enddo
    close(10)

     call enter_exit(sub_name,2)
 
  end subroutine write_node_geometry_2d

!!!#############################################################################

  subroutine write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
       ncount_loop,ncount_spline,np_offset)
    
    integer :: element_spline(:,:),elem_surfaces(:,:),ifile,ncount_point, &
         ncount_loop,ncount_spline,np_offset
    ! Local variables
    integer:: i,j,k,line1,line2,line3,line4,line_num(4),ne,ni1,ni2,nk,np,np_highest,np1,np2, &
         num_crux_lines,nv1,nv2
    integer,allocatable :: crux_lines(:,:)
    real(dp) :: phi_1_0,phi_1_1,phi_2_0,phi_2_1,point_xyz(3),xidivn(3)
    logical :: make
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'write_surface_geo'
    call enter_exit(sub_name,1)

    forall (i=1:3) xidivn(i) = 0.25_dp * i
    
!!! Make a gmsh 'point' at each of the surface mesh nodes
    write(ifile,'(''/* Points */'')')
    do np = 1,num_nodes_2d
       ncount_point = ncount_point + 1
       write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
            ncount_point,node_xyz_2d(1,1,1:3,np)
    enddo
    
    element_spline = 0
    num_crux_lines = 0
    ne = 1
    
    do while (ne.le.num_elems_2d)
       do k = 1,4
          if(k.eq.1)then
             ni1 = 1
             ni2 = 2
             nk = 2 ! calculate the location in the Xi+1 direction
             if(elem_cnct_2d(-2,0,ne).eq.0.or.ne.lt.elem_cnct_2d(-2,1,ne))then
                make = .true.
             else
                make = .false.
                if(elem_cnct_2d(-2,1,elem_cnct_2d(-2,1,ne)).eq.ne)then
                   ! the elements meet at Xi2=0 for both elements
                   line_num(k) = -element_spline(1,elem_cnct_2d(-2,1,ne))
                else
                   line_num(k) = element_spline(2,elem_cnct_2d(-2,1,ne))
                endif
             endif
          else if(k.eq.2)then
             ni1 = 3
             ni2 = 4
             nk = 2 ! calculate the location in the Xi+1 direction
             make = .true.
          else if(k.eq.3)then
             ni1 = 1
             ni2 = 3
             nk = 3 ! calculate the location in the Xi+2 direction
             if(elem_cnct_2d(-1,0,ne).eq.0.or.ne.lt.elem_cnct_2d(-1,1,ne))then
                make = .true.
             else
                make = .false.
                if(elem_cnct_2d(-2,1,elem_cnct_2d(-1,1,ne)).eq.ne)then
                   ! the elements change Xi direction where they meet
                   line_num(k) = element_spline(1,elem_cnct_2d(-1,1,ne))
                else
                   line_num(k) = element_spline(4,elem_cnct_2d(-1,1,ne))
                endif
             endif
          else if(k.eq.4)then
             ni1 = 2
             ni2 = 4
             nk = 3 ! calculate the location in the Xi+2 direction
             if(elem_cnct_2d(1,0,ne).eq.0.or.ne.lt.elem_cnct_2d(1,1,ne))then
                make = .true.
             else
                make = .false.
                if(elem_cnct_2d(-2,1,elem_cnct_2d(1,1,ne)).eq.ne)then
                   ! the elements change Xi direction where they meet
                   line_num(k) = -element_spline(1,elem_cnct_2d(1,1,ne))
                else
                   line_num(k) = element_spline(3,elem_cnct_2d(1,1,ne))
                endif
             endif
          endif

          if(make)then
             np1 = elem_nodes_2d(ni1,ne)
             np2 = elem_nodes_2d(ni2,ne)
             nv1 = elem_versn_2d(ni1,ne)
             nv2 = elem_versn_2d(ni2,ne)
             
             do i = 1,3
                phi_1_0 = 1.0_dp - 3.0_dp * xidivn(i)**2 + 2.0_dp * xidivn(i)**3
                phi_1_1 = xidivn(i) * (xidivn(i) - 1.0_dp)**2
                phi_2_0 = xidivn(i)**2 * (3.0_dp - 2.0_dp * xidivn(i))
                phi_2_1 = xidivn(i)**2 * (xidivn(i) - 1.0_dp)
                do j = 1,3
                   point_xyz(j) = phi_1_0 * node_xyz_2d(1,1,j,np1) + &
                        phi_2_0 * node_xyz_2d(1,1,j,np2) &
                        + phi_1_1 * node_xyz_2d(nk,nv1,j,np1) * scale_factors_2d((ni1-1)*4+nk,ne) &
                        + phi_2_1 * node_xyz_2d(nk,nv2,j,np2) * scale_factors_2d((ni2-1)*4+nk,ne)
                enddo
                ncount_point = ncount_point + 1
                write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
                     ncount_point, point_xyz(1:3)
             enddo ! i
             ncount_spline = ncount_spline + 1
             write(ifile,'(''Spline('',I8,'') = {'',I8,'','',I8,'','' &
                  &,I8,'','',I8,'','',I8,''};'')') ncount_spline, np1, &
                  ncount_point-2, ncount_point-1, ncount_point, np2
             element_spline(k,ne) = ncount_spline
             line_num(k) = ncount_spline
          endif
       enddo ! k
       
       ncount_loop = ncount_loop + 1
       write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,'','',i8,''};'')') &
            ncount_loop, line_num(1),line_num(4),-line_num(2),-line_num(3)
       ncount_loop = ncount_loop + 1
       write(ifile,'(''Surface('',I8,'') = {'',I8,''};'')') ncount_loop, ncount_loop - 1
       elem_surfaces(1,ne) = ncount_loop
       
       ne = ne + 1
       
    enddo ! while (ne.le.num_elems_2d)

    call enter_exit(sub_name,2)

  end subroutine write_surface_geo
  
!!!#############################################################################

  subroutine write_surface_geo0(element_spline,elem_surfaces,ifile,ncount_point, &
       ncount_loop,ncount_spline,np_offset)
    
    integer :: element_spline(:,:),elem_surfaces(:,:),ifile,ncount_point, &
         ncount_loop,ncount_spline,np_offset
    ! Local variables
    integer:: i,j,k,line1,line2,line3,line4,ne,nk,np,np_highest,np1,np2, &
         num_crux_lines,nv1,nv2
    integer,allocatable :: crux_lines(:,:)
    real(dp) :: phi_1_0,phi_1_1,phi_2_0,phi_2_1,point_xyz(3),xidivn(3)
    logical :: repeat
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'write_surface_geo'
    call enter_exit(sub_name,1)

    allocate(crux_lines(num_elems,3))
    
    forall (i=1:3) xidivn(i) = 0.25_dp * i
    
!!!    Make a gmsh 'point' at each of the surface mesh nodes
    write(ifile,'(''/* Points */'')')
    do np = 1,num_nodes_2d
       ncount_point = ncount_point + 1
       write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
            ncount_point,node_xyz_2d(1,1,1:3,np)
    enddo
    
!!!    Now work through each 'ring' of four adjoining surface elements. Points are created
!!!    between adjacent nodes in the Xi1 and Xi2 directions.

    element_spline = 0
    num_crux_lines = 0
    ne = 1
    
    do while (ne.le.num_elems_2d)
       
!!!....make intermediate points, lines, surfaces, at Xi2=0 for model entry .........
       if((elem_cnct_2d(-2,0,ne).eq.0) .or. &
            (node_versn_2d(elem_nodes_2d(1,ne)).eq.6.or.&
            node_versn_2d(elem_nodes_2d(2,ne)).eq.6))then
          ! location of points in the Xi+1 direction on the Xi2=0 boundary
          do k = 0,3
             if(elem_cnct_2d(-2,0,ne).eq.0 .or. element_spline(1,ne+k).eq.0)then  ! no line (or points) here already
                np1 = elem_nodes_2d(1,ne+k)
                np2 = elem_nodes_2d(2,ne+k)
                ! check that the line hasn't already been made. If it has, use it.
                repeat = .false.
                do i = 1,num_crux_lines
                   if((np1.eq.crux_lines(i,2).and.np2.eq.crux_lines(i,3)).or. &
                        (np1.eq.crux_lines(i,3).and.np2.eq.crux_lines(i,2)))then
                      repeat = .true.
                      element_spline(1,ne+k) = -crux_lines(i,1)
                   endif
                enddo

                if(.not.repeat)then
                   nv1 = elem_versn_2d(1,ne+k)
                   nv2 = elem_versn_2d(2,ne+k)
                   nk = 2 ! calculate the location in the Xi+1 direction
                   do i = 1,3
                      phi_1_0 = 1.0_dp - 3.0_dp * xidivn(i)**2 + 2.0_dp * xidivn(i)**3
                      phi_1_1 = xidivn(i) * (xidivn(i) - 1.0_dp)**2
                      phi_2_0 = xidivn(i)**2 * (3.0_dp - 2.0_dp * xidivn(i))
                      phi_2_1 = xidivn(i)**2 * (xidivn(i) - 1.0_dp)
                      do j = 1,3
                         point_xyz(j) = phi_1_0 * node_xyz_2d(1,1,j,np1) + &
                              phi_2_0 * node_xyz_2d(1,1,j,np2) &
                              + phi_1_1 * node_xyz_2d(nk,nv1,j,np1) * scale_factors_2d(2,ne+k) &
                              + phi_2_1 * node_xyz_2d(nk,nv2,j,np2) * scale_factors_2d(6,ne+k)
                      enddo
                      
                      ncount_point = ncount_point + 1
                      write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
                           ncount_point, point_xyz(1:3)
                   enddo ! i
                   ncount_spline = ncount_spline + 1
                   write(ifile,'(''Spline('',I8,'') = {'',I8,'','',I8,'','' &
                        &,I8,'','',I8,'','',I8,''};'')') ncount_spline, np1+np_offset, &
                        ncount_point-2, ncount_point-1, ncount_point, np2+np_offset
                   element_spline(1,ne+k) = ncount_spline
                   num_crux_lines = num_crux_lines + 1
                   crux_lines(num_crux_lines,1) = ncount_spline
                   crux_lines(num_crux_lines,2) = np1+np_offset
                   crux_lines(num_crux_lines,3) = np2+np_offset
!                   if(elem_cnct_2d(-2,0,ne).ne.0) element_spline(1,elem_cnct_2d(-2,1,ne)) = -ncount_spline
                endif
             endif
          enddo  ! k
       endif  ! elem_cnct etc
       
!!!.......make intermediate points, lines, surfaces, at Xi2=1 for each 'ring' .........
       
       ! location of points in the Xi+1 direction on the Xi2=1 boundary
       np_highest = elem_nodes_2d(3,ne)
       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          np2 = elem_nodes_2d(4,ne+k)
          nv1 = elem_versn_2d(3,ne+k)
          nv2 = elem_versn_2d(4,ne+k)
          nk = 2 
          do i = 1,3
             phi_1_0 = 1.0_dp - 3.0_dp * xidivn(i)**2 + 2.0_dp * xidivn(i)**3
             phi_1_1 = xidivn(i) * (xidivn(i) - 1.0_dp)**2
             phi_2_0 = xidivn(i)**2 * (3.0_dp - 2.0_dp * xidivn(i))
             phi_2_1 = xidivn(i)**2 * (xidivn(i) - 1.0_dp)
             forall (j=1:3) point_xyz(j) = phi_1_0 * node_xyz_2d(1,1,j,np1) &
                  + phi_2_0 * node_xyz_2d(1,1,j,np2) &
                  + phi_1_1 * node_xyz_2d(nk,nv1,j,np1) * scale_factors_2d(10,ne+k) &
                  + phi_2_1 * node_xyz_2d(nk,nv2,j,np2) * scale_factors_2d(14,ne+k) 
             
             ncount_point = ncount_point + 1
             write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
                  ncount_point,point_xyz(1:3)
          enddo ! i
          ncount_spline = ncount_spline + 1
          write(ifile,'(''Spline('',I8,'') = {'',I8,'','',I8,'','',I8,'','' &
               & ,I8,'','',I8,''};'')') ncount_spline, np1+np_offset, ncount_point-2, &
               ncount_point-1, ncount_point, np2+np_offset
          element_spline(3,ne+k) = ncount_spline
          if(elem_cnct_2d(2,0,ne+k).gt.0) element_spline(1,elem_cnct_2d(2,1,ne+k)) = ncount_spline
       enddo ! k
       
       ! location of points in the Xi+2 direction on the Xi1=0 boundary
       do k = 0,3
          np1 = elem_nodes_2d(1,ne+k)
          np2 = elem_nodes_2d(3,ne+k)
          nv1 = elem_versn_2d(1,ne+k)
          nv2 = elem_versn_2d(3,ne+k)
          nk = 3
          do i = 1,3
             phi_1_0 = 1.0_dp - 3.0_dp * xidivn(i)**2 + 2.0_dp * xidivn(i)**3
             phi_1_1 = xidivn(i) * (xidivn(i) - 1.0_dp)**2
             phi_2_0 = xidivn(i)**2 * (3.0_dp - 2.0_dp * xidivn(i))
             phi_2_1 = xidivn(i)**2 * (xidivn(i) - 1.0_dp)
             forall (j=1:3) point_xyz(j) = phi_1_0 * node_xyz_2d(1,1,j,np1) &
                  + phi_2_0 * node_xyz_2d(1,1,j,np2) &
                  + phi_1_1 * node_xyz_2d(nk,nv1,j,np1) * scale_factors_2d(3,ne+k) &
                  + phi_2_1 * node_xyz_2d(nk,nv2,j,np2) * scale_factors_2d(11,ne+k)
             
             ncount_point = ncount_point + 1
             write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
                  ncount_point,point_xyz(1:3)
          enddo ! i
          ncount_spline = ncount_spline + 1
          write(ifile,'(''Spline('',I8,'') = {'',I8,'','',I8,'','',I8,'','' &
               & ,I8,'','',I8,''};'')') ncount_spline, np1+np_offset, ncount_point-2, &
               ncount_point-1, ncount_point, np2+np_offset
          element_spline(4,ne+k) = ncount_spline
          element_spline(2,elem_cnct_2d(-1,1,ne+k)) = ncount_spline
       enddo ! k
       
       do k = 0,3
          line1 = element_spline(1,ne+k)
          line2 = element_spline(2,ne+k)
          line3 = -element_spline(3,ne+k)
          line4 = -element_spline(4,ne+k)
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3, line4
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',I8,'') = {'',I8,''};'')') ncount_loop, ncount_loop - 1
          elem_surfaces(3,ne+k) = ncount_loop
       enddo ! k
       
       ne = ne + 4
       
    enddo ! while (ne.le.num_elems_2d)

    deallocate(crux_lines)
    
    call enter_exit(sub_name,2)

  end subroutine write_surface_geo0
  
!!!#############################################################################

  subroutine write_3d_geo(element_spline,elem_surfaces,ifile,ncount_point, &
    ncount_loop,ncount_spline)

    integer,intent(in) :: ifile
    integer :: element_spline(:,:),elem_surfaces(:,:),ncount_point,ncount_loop, &
         ncount_spline
    ! Local variables
    integer :: i,j,k,line1,line2,line3,line4,ncount_cap_entry=0,ncount_cap_exit=0, &
         ncount_inner=0,ncount_centre=0,ncount_phys_vol=0,ncount_spline_0, &
         ncount_surface=0,ncount_volume=0,ncount_wall=0,ne,ne_next,np_highest,np1,np2
    integer,allocatable :: centre_points(:),ncap_entry(:),ncap_exit(:), &
         ncentre(:),ninner(:),nphys_vol(:),node_spoke(:,:),nwall(:)
    real(dp) :: point_xyz_centre(3), xidivn(3)
    logical :: at_bifurcation
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'write_3d_geo'
    call enter_exit(sub_name,1)
        
    ne = 1
    
    do while (ne.le.num_elems_2d)
!!!....the following works on four adjacent surface elements
       if(elem_cnct_2d(-2,0,ne).eq.0)then ! at the tree entry
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','' &
               &,i8,'','',i8,''};'')') ncount_loop, element_spline(1,ne:ne+3)
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          elem_surfaces(2,ne) = ncount_loop
       elseif(elem_cnct_2d(2,0,ne).eq.0)then ! at the tree exit
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','' &
               &,i8,'','',i8,''};'')') ncount_loop, element_spline(2,ne:ne+3)
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          elem_surfaces(2,ne) = ncount_loop
       endif
       
       ne = ne + 4

    enddo ! while (ne.le.num_elems_2d)

    !write(ifile,'(/''Surface Loop(1) = {'')', advance = "no")
    !ne = 1
    !do while (ne.le.num_elems_2d)
     !  do k = 0,3
     !     write(ifile,'(i6,'','')', advance = "no") elem_surfaces(1,ne+k)
     !  enddo
     !  if(elem_cnct_2d(-2,0,ne).eq.0)then ! at the tree entry
     !     write(ifile,'(i6,'','')', advance = "no") elem_surfaces(2,ne)
     !  elseif(elem_cnct_2d(2,0,ne).eq.0)then ! at the tree exit
     !     write(ifile,'(i6,'','')', advance = "no") elem_surfaces(2,ne)
     !  endif
     !  
     !  ne = ne + 4

    !enddo ! while (ne.le.num_elems_2d)
    !write(ifile,'(''};'')')
    
    !write(ifile,'(/''Volume(2) = {1};'')')

    write(ifile,'(/)')
    write(ifile,'(''Mesh.Algorithm = 3;'')') ! Anisotropic
    write(ifile,'(''Mesh.Smoothing = 4;'')')
    write(ifile,'(''Mesh.Algorithm3D = 2;'')') ! Netgen

    close(ifile)

    call enter_exit(sub_name,2)
    
  end subroutine write_3d_geo

!!!#############################################################################

  subroutine write_3d_geo0(element_spline,elem_surfaces,ifile,ncount_point, &
    ncount_loop,ncount_spline)

    integer,intent(in) :: ifile
    integer :: element_spline(:,:),elem_surfaces(:,:),ncount_point,ncount_loop, &
         ncount_spline
    ! Local variables
    integer :: i,j,k,line1,line2,line3,line4,ncount_cap_entry=0,ncount_cap_exit=0, &
         ncount_inner=0,ncount_centre=0,ncount_phys_vol=0,ncount_spline_0, &
         ncount_surface=0,ncount_volume=0,ncount_wall=0,ne,ne_next,np_highest,np1,np2
    integer,allocatable :: centre_points(:),ncap_entry(:),ncap_exit(:), &
         ncentre(:),ninner(:),nphys_vol(:),node_spoke(:,:),nwall(:)
    real(dp) :: point_xyz_centre(3), xidivn(3)
    logical :: at_bifurcation
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'write_3d_geo'
    call enter_exit(sub_name,1)
        
    allocate(centre_points(num_nodes_2d))
    allocate(ncap_entry(20))  ! assuming could have multiple inlets
    allocate(ncap_exit(num_elems_2d))  
    allocate(ninner(num_elems_2d*2))  
    allocate(ncentre(num_elems_2d))
    allocate(nphys_vol(num_elems_2d))
    allocate(nwall(num_elems_2d))  
    allocate(node_spoke(2,num_nodes_2d))
    node_spoke = 0

    forall (i=1:3) xidivn(i) = 0.25_dp * i

    ne = 1
    do while (ne.le.num_elems_2d)
!!!....the following works on four adjacent surface elements

!!!......... make intermediate points, lines, surfaces, at Xi2=0 for model entry .........
       
       if((elem_cnct_2d(-2,0,ne).eq.0) .or. &
            (node_versn_2d(elem_nodes_2d(1,ne)).eq.6.or.&
            node_versn_2d(elem_nodes_2d(2,ne)).eq.6))then
          ! location of points in the Xi+1 direction on the Xi2=0 boundary
          point_xyz_centre = 0.0_dp
          do k = 0,3
             np1 = elem_nodes_2d(1,ne+k)
             forall (j = 1:3) point_xyz_centre(j) = &
                  point_xyz_centre(j) + 0.25_dp * node_xyz_2d(1,1,j,np1)
          enddo  ! k

          if(elem_cnct_2d(-2,0,ne).eq.0)then
             ! make a point at the centre of the ring
             ncount_point = ncount_point + 1
             write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','' &
                  &,f12.7,'',lc};'')') ncount_point,point_xyz_centre(1:3)
             
             ! make a 'spoke' from centre of the ring to each surface node
             do k = 0,3
                np1 = elem_nodes_2d(1,ne+k)
                centre_points(np1) = ncount_point   ! record the centre point number for this 'ring'
                ncount_spline = ncount_spline + 1
                write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') &
                     ncount_spline,ncount_point,np1
                node_spoke(1,np1) = ncount_spline
             enddo  ! k

             ! make surfaces at the entry == a 'cap' of four surfaces
             do k = 0,3
                line1 = ncount_spline + k - 3 
                line2 = element_spline(1,ne+k)
                if(k.lt.3)then
                   line3 = -(line1 + 1)
                else
                   line3 = -(ncount_spline - 3)
                endif
                ncount_loop = ncount_loop + 1
                write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','' &
                     &,i8,''};'')') ncount_loop, line1, line2, line3
                ncount_loop = ncount_loop + 1
                write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
                     ncount_loop, ncount_loop - 1
                elem_surfaces(1,ne+k) = ncount_loop
                ncount_cap_entry = ncount_cap_entry + 1
                ncap_entry(ncount_cap_entry) = ncount_loop
             enddo  ! k
          endif
       endif  ! elem_cnct etc
       
!!!......... make intermediate points, lines, surfaces, at Xi2=1 for each 'ring' .........

       ! location of points in the Xi+1 direction on the Xi2=1 boundary
       point_xyz_centre = 0.0_dp
       np_highest = elem_nodes_2d(3,ne)
       ncount_spline_0 = ncount_spline + 1

       at_bifurcation = .false.
       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          np2 = elem_nodes_2d(4,ne+k)
          forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) &
               + node_xyz_2d(1,1,j,np1)
          if(node_versn_2d(np1).ge.6.or.node_versn_2d(np2).ge.6)  at_bifurcation = .true.
          if(np1.gt.np_highest) np_highest = np1
       enddo ! k
       ncount_point = ncount_point + 1
       if(at_bifurcation)then
          np_highest = np_highest + 1
          centre_points(np_highest) = ncount_point   ! record the centre point number for this 'ring'
          forall (j = 1:3) point_xyz_centre(j) = (point_xyz_centre(j) + node_xyz_2d(1,1,j,np_highest)) * 0.2_dp
       else
          forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) * 0.25_dp
       endif
       
       write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','' &
            &,f12.7,'',lc};'')') ncount_point,point_xyz_centre(1:3)

       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          centre_points(np1) = ncount_point   ! record the centre point number for this 'ring'
          ncount_spline = ncount_spline + 1
          write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') &
               ncount_spline,ncount_point,np1
          node_spoke(1,np1) = ncount_spline
       enddo  ! k

! not sure whether the following is needed. doesn't seem to be used       
!       if(at_bifurcation)then
!       !if(node_versn_2d(np1).ge.2.or.node_versn_2d(np2).ge.2)then  ! only for when there is a crux node
!          !np_highest = elem_nodes_2d(2,elem_cnct_2d(1,1,elem_cnct_2d(2,1,ne)))
          
!          if(elems_at_node_2d(np_highest,0).eq.6)then
!             do i = 1,elems_at_node_2d(np_highest,0)
!                ne_next = elems_at_node_2d(np_highest,i)
!                np1 = elem_nodes_2d(1,ne_next)
!                if(np1.eq.np_highest) np1 = &
!                     elem_nodes_2d(2,elems_at_node_2d(np_highest,i))
!                if(node_spoke(1,np1).eq.0)then
!                   ! make a spoke from the centre to this node
!                   ncount_spline = ncount_spline + 1
!                   write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') &
!                        ncount_spline,ncount_point,np1
!                   node_spoke(1,np1) = ncount_spline
!                   elem_surfaces(1,elem_cnct_2d(-2,1,ne_next)) = ncount_loop
!                endif
!             enddo
!          endif
!       endif

!!!....make surface elements at the Xi2=1 end
       do k = 0,3
          line1 = node_spoke(1,elem_nodes_2d(3,ne+k))
          line2 = element_spline(3,ne+k)
          line3 = -node_spoke(1,elem_nodes_2d(4,ne+k))
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          elem_surfaces(2,ne+k) = ncount_loop
          if(elem_cnct_2d(2,0,ne+k).ne.0)then
             elem_surfaces(1,elem_cnct_2d(2,1,ne+k)) = ncount_loop
          else
             ! this is an exit 'cap'. store to output as a physical surface
             ncount_cap_exit = ncount_cap_exit + 1
             ncap_exit(ncount_cap_exit) = ncount_loop
          endif
          ncount_inner = ncount_inner + 1
          ninner(ncount_inner) = ncount_loop
       enddo  ! k
             
       if(at_bifurcation)then
          ncount_spline = ncount_spline + 1
          write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') &
               ncount_spline,ncount_point,np_highest
          node_spoke(1,np_highest) = ncount_spline

          ne_next = elems_at_node_2d(np_highest,1)
          line1 = node_spoke(1,elem_nodes_2d(1,ne_next))
          line2 = element_spline(1,ne_next)
          line3 = -node_spoke(1,elem_nodes_2d(2,ne_next))

          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          ncount_inner = ncount_inner + 1
          ninner(ncount_inner) = ncount_loop

          elem_surfaces(1,ne_next) = ncount_loop
          elem_surfaces(1,elem_cnct_2d(-2,1,ne_next)) = ncount_loop

          ne_next = elems_at_node_2d(np_highest,2)
          line1 = node_spoke(1,elem_nodes_2d(1,ne_next))
          line2 = element_spline(1,ne_next)
          line3 = -node_spoke(1,elem_nodes_2d(2,ne_next))
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          ncount_inner = ncount_inner + 1
          ninner(ncount_inner) = ncount_loop

          elem_surfaces(1,ne_next) = ncount_loop
          ne_next = elem_cnct_2d(-2,1,ne_next)
          elem_surfaces(1,ne_next) = ncount_loop
          
          if(elems_at_node_2d(np_highest,0).eq.6)then
             do i = 1,elems_at_node_2d(np_highest,0)
                ne_next = elems_at_node_2d(np_highest,i)
                if(elem_surfaces(1,ne_next).eq.0)then
                   ! make an extra surface here
                   np1 = elem_nodes_2d(1,ne_next)
                   np2 = elem_nodes_2d(2,ne_next)
                   line1 = node_spoke(1,np1)
                   line2 = element_spline(1,ne_next)
                   line3 = -node_spoke(1,np2)
                   ncount_loop = ncount_loop + 1
                   write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','' &
                        &,i8,'','',i8,''};'')') &
                        ncount_loop, line1, line2, line3
                   ncount_loop = ncount_loop + 1
                   write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
                        ncount_loop, ncount_loop - 1
                   ncount_inner = ncount_inner + 1
                   ninner(ncount_inner) = ncount_loop
                   elem_surfaces(1,ne_next) = ncount_loop
                   ne_next = elem_cnct_2d(-2,1,ne_next)
                   elem_surfaces(1,ne_next) = ncount_loop
                   
                   ne_next = elem_cnct_2d(1,1,ne_next)
                   np1 = elem_nodes_2d(1,ne_next)
                   np2 = elem_nodes_2d(2,ne_next)
                   line1 = node_spoke(1,np1)
                   line2 = element_spline(1,ne_next)
                   line3 = -node_spoke(1,np2)
                   ncount_loop = ncount_loop + 1
                   write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','' &
                        &,i8,'','',i8,''};'')') &
                        ncount_loop, line1, line2, line3
                   ncount_loop = ncount_loop + 1
                   write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
                        ncount_loop, ncount_loop - 1
                   ncount_inner = ncount_inner + 1
                   ninner(ncount_inner) = ncount_loop
                   elem_surfaces(1,ne_next) = ncount_loop
                   ne_next = elem_cnct_2d(-2,1,ne_next)
                   elem_surfaces(1,ne_next) = ncount_loop
                   
                endif
             enddo
          endif
       endif

!!!.........Make line along the centre
       ncount_spline = ncount_spline + 1
       write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') ncount_spline, &
            centre_points(elem_nodes_2d(1,ne)),ncount_point
 
!!! Make surfaces from the centreline
       do k = 0,3
          line1 = node_spoke(1,elem_nodes_2d(1,ne+k))
          line2 = element_spline(4,ne+k)
          line3 = -node_spoke(1,elem_nodes_2d(3,ne+k))
          line4 = -ncount_spline ! the newest line
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','' &
               &,i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3, line4
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          ncount_centre = ncount_centre + 1
          ncentre(ncount_centre) = ncount_loop

          elem_surfaces(4,ne+k) = ncount_loop
          if(k.gt.0) elem_surfaces(5,ne+k-1) = ncount_loop
       enddo
       elem_surfaces(5,ne) = elem_surfaces(4,ne+1)
       elem_surfaces(5,ne+3) = elem_surfaces(4,ne)

!!! Make surface loops and volumes
       do k = 0,3
          ncount_volume = ncount_volume + 1
          write(ifile,'(''Surface Loop('',i8,'') = {'',i8,'','' &
               &,i8,'','',i8,'','',i8,'','',i8,''};'')') &
               ncount_volume,elem_surfaces(1:5,ne+k)
            
          ncount_volume = ncount_volume + 1
          write(ifile,'(''Volume('',i8,'') = {'',i8,''};'')') &
               ncount_volume,ncount_volume-1

          ncount_wall = ncount_wall + 1
          nwall(ncount_wall) = elem_surfaces(3,ne+k)
          ncount_phys_vol = ncount_phys_vol + 1
          nphys_vol(ncount_phys_vol) = ncount_volume
       enddo

       ne = ne + 4

    enddo ! while (ne.le.num_elems_2d)

    write(ifile,'(/''/* Physical surface for entry caps */'')')
    write(ifile,'(/''Physical Surface(1) = {'')', advance = "no")
    do i = 1,ncount_cap_entry-1
       write(ifile,'(i6,'','')', advance = "no") ncap_entry(i)
    enddo
    write(ifile,'(i6,''};'')') ncap_entry(ncount_cap_entry)
    
    write(ifile,'(/''/* Physical surface for exit caps */'')')
    write(ifile,'(/''Physical Surface(2) = {'')', advance = "no")
    do i = 1,ncount_cap_exit-1
       write(ifile,'(i6,'','')', advance = "no") ncap_exit(i)
    enddo
    write(ifile,'(i6,''};'')') ncap_exit(ncount_cap_exit)
    
    write(ifile,'(/''/* Physical surface for walls */'')')
    write(ifile,'(/''Physical Surface(3) = {'')', advance = "no")
    do i = 1,ncount_wall-1
       write(ifile,'(i6,'','')', advance = "no") nwall(i)
    enddo
    write(ifile,'(i6,''};'')') nwall(ncount_wall)
    
    write(ifile,'(/''/* Physical surface for centres */'')')
    write(ifile,'(/''Physical Surface(4) = {'')', advance = "no")
    do i = 1,ncount_centre-1
       write(ifile,'(i6,'','')', advance = "no") ncentre(i)
    enddo
    write(ifile,'(i6,''};'')') ncentre(ncount_centre)
    
    write(ifile,'(/''Physical Volume(1) = {'')', advance = "no")
    do i = 1,ncount_phys_vol-1
       write(ifile,'(i6,'','')', advance = "no") nphys_vol(i)
    enddo
    write(ifile,'(i6,''};'')') nphys_vol(ncount_phys_vol)

    write(ifile,'(/)')
    write(ifile,'(''Mesh.Algorithm = 3;'')') ! Anisotropic
    write(ifile,'(''Mesh.Smoothing = 4;'')')
    write(ifile,'(''Mesh.Algorithm3D = 2;'')') ! Netgen

    close(ifile)

    deallocate(centre_points)
    deallocate(ncap_entry)
    deallocate(ncap_exit)
    deallocate(ninner)
    deallocate(ncentre)
    deallocate(nphys_vol)
    deallocate(nwall)
    deallocate(node_spoke)

    call enter_exit(sub_name,2)
    
  end subroutine write_3d_geo0

!!!#############################################################################
  
end module geometry

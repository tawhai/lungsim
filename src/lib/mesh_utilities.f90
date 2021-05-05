module mesh_utilities

!!! Subroutines and functions for general calculations. Not specific to any
!!! particular application, although it is expected these will generally be used 
!!! for calculations to do with a mesh (not fields or solutions). 
!!! Any function that is used by more than one module should appear in here. 
!!! ALL subroutines and functions in this module are public.


  use arrays
  use diagnostics
  use other_consts
  use precision

  implicit none

  private

  public  area_between_three_points
  public  area_between_two_vectors
  public  angle_btwn_points
  public  angle_btwn_vectors
  public  bifurcation_element
  public  calc_arclengths
  public  calc_branch_direction
  public  calc_scale_factors_2d
  public  check_colinear_points
  public  coord_at_xi
  public  cross_product
  public  direction_point_to_point
  public  distance_between_points
  public  distance_from_plane_to_point
  public  element_connectivity_1d
  public  element_connectivity_2d
  public  evaluate_ordering
  public  get_local_elem
  public  get_local_elem_2d
  public  get_local_node_f
  public  group_elem_parent_term
  public  hermite
  public  inlist
  public  linear
  public  make_plane_from_3points
  public  mesh_a_x_eq_b
  public  point_internal_to_surface
  public  reallocate_node_elem_arrays
  public  scalar_product_3
  public  scalar_triple_product
  public  scale_mesh
  public  set_linear_derivatives
  public  stem_element
  public  terminal_element
  public  triangles_from_surface
  public  unit_norm_to_plane_two_vectors
  public  unit_norm_to_three_points
  public  unit_vector
  public  vector_length
  public  volume_internal_to_surface
  public  which_child

contains

!!! list of subroutines

  ! calc_branch_direction
  ! ..... calculates direction of branch and stores in elem_direction

  ! calc_scale_factors_2d
  ! ..... calculate the scale factors for a 2d mesh

  ! make_plane_from_3points
  ! ..... finds the equation of a plane in 3D and a vector normal to the plane from three
  ! ..... non-colinear points.

  ! scale_mesh
  ! ..... multiply mesh (coordinates and derivatives) by a constant

!!! list of functions

  ! angle_btwn_points
  ! .... returns the angle between three points

  ! angle_btwn_vectors
  ! .... returns the angle between two vectors
 
  ! check_colinear_points 
  ! .... returns true or false for whether 3 points are colinear

  ! cross_product ***
  ! .... returns the vector cross product of A*B

  ! distance_between_points ***
  ! .... calculates the distance between two arbitrary points

  ! mesh_a_x_eq_b ***
  ! .... solves a small matrix system

  ! scalar_product_3 ***
  ! .... dot product of two 3x1 vectors

  ! unit_vector ***
  ! .... Calculates the unit vector for an arbitrary 3x1 vector

  ! vector_length ***
  ! .... Calculates the length of a 3x1 vector

!!!#####################################################################

  subroutine calc_arclengths
    !*calc_arclengths*: estimates arclength for cubic Hermite using
    ! Gaussian quadrature with 4 points

    ! Local variables
    integer :: i,it,itmax,n,ng,nj,nl,nline,np,nv,xi_direction
    real(dp) :: est_length,incr_length,linear_est,line_xyz(2,3,2), &
         local_deriv(3),weight(4),xigg(4)

    xigg = [0.0694318442029_dp, 0.3300094782075_dp,&
         0.6699905217924_dp, 0.9305681557970_dp]    ! exact Gauss point locations
    weight = [0.1739274225687_dp, 0.3260725774313_dp,&
         0.3260725774313_dp, 0.1739274225687_dp]    ! exact Gauss point weightings

    do nline = 1,num_lines_2d                ! loop over all lines
       nl = lines_2d(nline)                  ! the line number 
       xi_direction = nodes_in_line(1,0,nl)  ! the Xi direction of the line
       
       ! get the nodal coordinates and derivatives for the line
       do n = 1,2                            ! for each node on the line
          np = nodes_in_line(n+1,1,nl)       ! the first and second node (np1,np2)
          nv = line_versn_2d(N,1,nl)         ! the version of the node for this line
          line_xyz(1,:,n) = node_xyz_2d(1,nv,:,np)  ! get the coordinates for the line
          if(xi_direction.eq.1) line_xyz(2,:,n) = node_xyz_2d(2,nv,:,np) ! dxi1
          if(xi_direction.eq.2) line_xyz(2,:,n) = node_xyz_2d(3,nv,:,np) ! dxi2
       enddo !n
       
       ! calculate the linear distance between start and end nodes. this should be
       ! used to check that the derivatives are appropriate when the start and 
       ! end nodes are coincident (i.e. a collapsed element). 
       est_length = 0.0_dp
       do ng = 1,4
          ! calculate the arclength derivatives
          do nj = 1,3
             ! local_deriv = phi_10' * xyz_1 + phi_20' * xyz_2
             local_deriv(nj) = linear(1,2,xigg(ng))*line_xyz(1,nj,1) &
                  + linear(2,2,xigg(ng))*line_xyz(1,nj,2)
          enddo
          incr_length = sqrt(scalar_product_3(local_deriv,local_deriv))
          est_length = est_length + weight(ng) * incr_length
       enddo !ng
       
       linear_est = est_length

       est_length = 0.0_dp
       do ng = 1,4
          ! calculate the arclength derivatives at Xi coordinates 
          do nj = 1,3
             !function' = phi_10'*xyz_1 + phi_11'*deriv_1 + phi_20'*xyz_2 + phi_21'*deriv_2
             local_deriv(nj) = hermite(1,1,2,xigg(ng))*line_xyz(1,nj,1) &
                  + hermite(1,2,2,xigg(ng))*line_xyz(2,nj,1) &
                  + hermite(2,1,2,xigg(ng))*line_xyz(1,nj,2) &
                  + hermite(2,2,2,xigg(ng))*line_xyz(2,nj,2)
          enddo
          incr_length = sqrt(scalar_product_3(local_deriv,local_deriv))
          est_length = est_length + weight(ng) * incr_length
       enddo !ng
       if(abs(linear_est).gt.zero_tol)then
          arclength(nl) = est_length
       else
          arclength(nl) = 0.0_dp
          np = nodes_in_line(2,1,nl)         ! the first node
          nv = line_versn_2d(1,1,nl)         ! the version of the node for this line
          node_xyz_2d(xi_direction+1,nv,:,np) = 0.0_dp
       endif
    enddo !loop over lines

  end subroutine calc_arclengths
  
!!!#####################################################################
  
  subroutine calc_branch_direction(ne)
    
!!! calculates the direction of element ne and stores in elem_direction    
    
    integer,intent(in) :: ne
    
    integer :: np_end,np_start
    real(dp) :: length
    
    np_start = elem_nodes(1,ne)
    np_end = elem_nodes(2,ne)
    
    length = distance_between_points(node_xyz(1,np_end),node_xyz(1,np_start))
    elem_direction(1:3,ne)=(node_xyz(1:3,np_end)-node_xyz(1:3,np_start))/length

  end subroutine calc_branch_direction

!!! ##########################################################################      

  subroutine calc_scale_factors_2d(sf_option)
    !*calc_scale_factors_2d*: calculates arclengths using Gaussian quadrature,
    ! and scale factors for 2d surface elements
  
    character(len=4),intent(in) :: sf_option
!!! local variables
    integer,parameter :: num_deriv = 4
    integer :: i,ido(num_deriv,2),it,ITMAX=20,k,N,NAE,ne,&
         ng,NGA=4,NI1(3),ni,ni2,nj,nk,nk2,nl,nline,nn,nn2,NNK,&
         np,ns,nv,NNL(2,4)
    real(dp) :: DA,SUM1,SUM2,SUM3,SUM4,W
    logical :: FOUND
    
    ido = reshape ([1,2,1,2,1,1,2,2],shape(ido))
    NI1 = [1,2,1]
    NNL = reshape([1,2,3,4,1,3,2,4],shape(NNL))


    if(.not.allocated(scale_factors_2d)) allocate(scale_factors_2d(16,num_elems_2d))

    select case (sf_option)
    case ('unit')
       scale_factors_2d = 1.0_dp
       
    case('arcl')

       call calc_arclengths
       
       ! calculate scale factors using the line derivatives 
       scale_factors_2d = 1.0_dp !initialise
       
       do ne = 1,num_elems_2d
          do NAE = 1,4      ! for each of the (up to) 4 lines in the element
             nl = elem_lines_2d(NAE,ne)     ! the line number
             if(nl /= 0)then                ! required for collapsed edges
                ni = nodes_in_line(1,0,nl)  ! the Xi direction of the line
                ni2 = NI1(ni+1)             ! 2,1 for Xi 1,2 resp.
                do N=1,2
                   nn=NNL(N,NAE)            ! 1,2,3,4 for n=1; 1,3,2,4 for n=2
                   ns=0                     
                   do nn2=1,nn-1
                      do nk2=1,num_deriv
                         ns=ns+1
                      enddo
                   enddo
                   do nk=2,num_deriv
                      if(IDO(nk,ni2).EQ.1) then
                         scale_factors_2d(nk+ns,ne) = arclength(nl)
                         if(abs(scale_factors_2d(nk+ns,ne)).LT.1.0e-6_dp) scale_factors_2d(nk+ns,ne) = 1.0_dp
                      endif
                   enddo !nk
                enddo !N=1,2
             endif
          enddo !NAE (nl)

          forall(i = 0:12:4) scale_factors_2d(i+4,ne) = scale_factors_2d(i+2,ne)* &
               scale_factors_2d(i+3,ne)
          
       enddo !noelem (ne)
    end select

  end subroutine calc_scale_factors_2d
  
!!!###############################################################
  
  function distance_from_plane_to_point(P1,P2,P3,P4)
    
    !###    calculates the distance from a plane defined by three points
    !###    and another arbitrary point
    
    real(dp),intent(in) :: P1(3),P2(3),P3(3),P4(3)
    real(dp) :: norml(4)
    real(dp) :: distance_from_plane_to_point

    call make_plane_from_3points(norml,2,P1,P2,P3)
    distance_from_plane_to_point = abs(scalar_product_3(norml,P4) + norml(4)) / &
         sqrt(scalar_product_3(norml,norml))
    
  end function distance_from_plane_to_point
  
!!!#############################################################################
  
  subroutine element_connectivity_1d()
    !*element_connectivity_1d:*  Calculates element connectivity in 1D and
    ! stores in array elem_cnct
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ELEMENT_CONNECTIVITY_1D" :: ELEMENT_CONNECTIVITY_1D

    !     Local Variables
    integer :: ne,ne2,nn,noelem,np,np2,np1
    integer,parameter :: NNT=2
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'element_connectivity_1d'
    call enter_exit(sub_name,1)
    
    elem_cnct = 0 !initialise

    ! calculate elems_at_node array: stores the elements that nodes are in
    elems_at_node(1:num_nodes,0) = 0 !initialise number of adjacent elements
    
    do ne=1,num_elems
       do nn=1,2
          np=elem_nodes(nn,ne)
          elems_at_node(np,0)=elems_at_node(np,0)+1
          elems_at_node(np,elems_at_node(np,0))=ne ! local element that np is in
       enddo !nn
    enddo !noelem
    
    ! calculate elem_cnct array: stores the connectivity of all elements
    
    elem_cnct=0 !initialise all elem_cnct
    
    do ne=1,num_elems
       !     ne_global=elems(noelem)
       if(NNT == 2) THEN !1d
          np1=elem_nodes(1,ne) !first local node
          np2=elem_nodes(2,ne) !second local node
          do noelem=1,elems_at_node(np2,0)
             ne2=elems_at_node(np2,noelem)
             if(ne2 /= ne)THEN
                elem_cnct(-1,0,ne2)=elem_cnct(-1,0,ne2)+1
                elem_cnct(-1,elem_cnct(-1,0,ne2),ne2)=ne !previous element
                elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
                elem_cnct(1,elem_cnct(1,0,ne),ne)=ne2
             endif !ne2
          enddo !noelem2
       endif
    enddo
    
    call enter_exit(sub_name,2)

  end subroutine element_connectivity_1d

!!!#############################################################################
  
  subroutine element_connectivity_2d
    !*element_connectivity_2d:*  Calculates element connectivity in 2D and
    ! stores in array elem_cnct_2d

    ! Local variables
    integer :: ne,ne2,nn,np,noelem2,np_list(4),np_list_2(4)
    integer,parameter :: num_elem_nodes = 4
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'element_connectivity_2d'
    call enter_exit(sub_name,1)
    
    if(allocated(elems_at_node_2d))then
       deallocate(elem_cnct_2d)
       deallocate(elems_at_node_2d)
    endif
    allocate(elem_cnct_2d(-2:2,0:10,num_elems_2d))
    allocate(elems_at_node_2d(num_nodes_2d,0:10))
    
!!! calculate elems_at_node_2d array: stores the elements that nodes are in
    
    elems_at_node_2d = 0 !initialise all
    do ne = 1,num_elems_2d
       do nn = 1,num_elem_nodes
          np = elem_nodes_2d(nn,ne)
          elems_at_node_2d(np,0) = elems_at_node_2d(np,0)+1
          elems_at_node_2d(np,elems_at_node_2d(np,0)) = ne !element that np is in
       enddo !nn
    enddo !noelem
!!! calculate elem_cnct_2d array: stores the connectivity of all elements
    
    elem_cnct_2d = 0 !initialise all elem_cnct_2d
    
    do ne = 1,num_elems_2d ! for each of the 2d elements
       np_list(1:4) = elem_nodes_2d(1:4,ne) ! the list of nodes in the element (including repeated)
!!! check the elements attached to the 1st node
       do noelem2 = 1,elems_at_node_2d(np_list(1),0) ! for each element attached to the 1st node
          ne2 = elems_at_node_2d(np_list(1),noelem2) ! attached element number
          if(ne2.ne.ne)then
             np_list_2(1:4) = elem_nodes_2d(1:4,ne2) !list of nodes in attached element
             if(np_list(2).ne.np_list(1))then ! only if first two nodes are not repeated
                if(inlist(np_list(2),np_list_2))then
                   elem_cnct_2d(-2,0,ne) = elem_cnct_2d(-2,0,ne)+1
                   elem_cnct_2d(-2,elem_cnct_2d(-2,0,ne),ne) = ne2 
                endif
             endif
             if(np_list(3).ne.np_list(1))then ! only if the two nodes are not repeated
                if(inlist(np_list(3),np_list_2))then
                   elem_cnct_2d(-1,0,ne) = elem_cnct_2d(-1,0,ne)+1
                   elem_cnct_2d(-1,elem_cnct_2d(-1,0,ne),ne) = ne2 
                endif
             endif
          endif
       enddo
!!! check the elements attached to the 4th node
       do noelem2 = 1,elems_at_node_2d(np_list(4),0) ! for each element attached to the 4th node
          ne2 = elems_at_node_2d(np_list(4),noelem2) ! attached element number
          if(ne2.ne.ne)then
             np_list_2(1:4) = elem_nodes_2d(1:4,ne2) !list of nodes in attached element
             if(np_list(2).ne.np_list(4))then ! only if two nodes are not repeated
                if(inlist(np_list(2),np_list_2))then
                   elem_cnct_2d(1,0,ne) = elem_cnct_2d(1,0,ne)+1
                   elem_cnct_2d(1,elem_cnct_2d(1,0,ne),ne) = ne2 
                endif
             endif
             if(np_list(3).ne.np_list(4))then ! only if the two nodes are not repeated
                if(inlist(np_list(3),np_list_2))then
                   elem_cnct_2d(2,0,ne) = elem_cnct_2d(2,0,ne)+1
                   elem_cnct_2d(2,elem_cnct_2d(2,0,ne),ne) = ne2 
                endif
             endif
          endif
       enddo !noelem2
    enddo
    
    call enter_exit(sub_name,2)
    
  end subroutine element_connectivity_2d

!!!#############################################################################

  subroutine evaluate_ordering()
    !*evaluate_ordering:* calculates generations, Horsfield orders,
    ! Strahler orders for a given tree
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_ORDERING" :: EVALUATE_ORDERING

    ! Local Variables
    integer :: INLETS,ne,ne0,ne2,noelem2,np,np2,nn,num_attach,n_children, &
         n_generation,n_horsfield,OUTLETS,STRAHLER,STRAHLER_ADD,temp1
    LOGICAL :: DISCONNECT,DUPLICATE
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'evaluate_ordering'
    call enter_exit(sub_name,1)
    
!!! Calculate branch generations
    elem_ordrs = 0
    maxgen=1
    do ne=1,num_elems
       ne0=elem_cnct(-1,1,ne) !parent
       if(ne0.NE.0)THEN
          n_generation=elem_ordrs(1,ne0) !parent generation
          if(elem_cnct(1,0,ne0).EQ.1)THEN !single daughter
             elem_ordrs(1,ne)=n_generation + (elem_symmetry(ne)-1)
          else if(elem_cnct(1,0,ne0).GE.2)THEN
             elem_ordrs(1,ne)=n_generation+1
          endif
       else
          elem_ordrs(1,ne)=1 !generation 1
       endif
       maxgen=max(maxgen,elem_ordrs(1,ne))
    enddo !noelem

!!! Calculate the branch orders
    do ne=num_elems,1,-1
       n_horsfield=MAX(elem_ordrs(2,ne),1)
       n_children=elem_cnct(1,0,ne) !number of child branches
       if(n_children.EQ.1)THEN
          if(elem_ordrs(1,elem_cnct(1,1,ne)).EQ.0)  n_children=0
       endif
       STRAHLER=0
       STRAHLER_ADD=1
       if(n_children.GE.2)THEN !branch has two or more daughters
          STRAHLER=elem_ordrs(3,elem_cnct(1,1,ne)) !first daughter
          do noelem2=1,n_children !for all daughters
             ne2=elem_cnct(1,noelem2,ne) !global element # of daughter
             temp1=elem_ordrs(2,ne2) !Horsfield order of daughter
             if(temp1.GT.n_horsfield)then
                n_horsfield=temp1
             endif
             if(elem_ordrs(3,ne2).LT.STRAHLER)THEN
                STRAHLER_ADD=0
             else if(elem_ordrs(3,ne2).GT.STRAHLER)THEN
                STRAHLER_ADD=0
                STRAHLER=elem_ordrs(3,ne2) !highest daughter
             endif
          enddo !noelem2 (ne2)
          n_horsfield=n_horsfield+1 !Horsfield ordering
       else if(n_children.EQ.1)THEN
          ne2=elem_cnct(1,1,ne) !local element # of daughter
          n_horsfield=elem_ordrs(2,ne2)+(elem_symmetry(ne)-1)
          STRAHLER_ADD=elem_ordrs(3,ne2)+(elem_symmetry(ne)-1)
       endif !elem_cnct
       elem_ordrs(2,ne)=n_horsfield !store the Horsfield order
       elem_ordrs(3,ne)=STRAHLER+STRAHLER_ADD !Strahler order
    enddo !noelem
    
!!! Check for disconnected nodes and number of inlets and outlets
    DUPLICATE=.FALSE.
    do ne=1,num_elems
       np=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       if(np.EQ.np2)THEN
          DUPLICATE=.TRUE.
       endif
    enddo
    
    DISCONNECT=.FALSE.
    INLETS=0
    OUTLETS=0
    do np=1,num_nodes
       num_attach=elems_at_node(np,0)
       if(num_attach.EQ.0)THEN
          DISCONNECT=.TRUE.
       elseif(num_attach.EQ.1)THEN
          ne=elems_at_node(np,1)
          if(elem_cnct(1,0,ne).EQ.0) OUTLETS=OUTLETS+1
          if(elem_cnct(-1,0,ne).EQ.0) INLETS=INLETS+1
       elseif(num_attach.GT.3)THEN
          WRITE(*,*) ' Node ',np,' attached to',num_attach,' elements'
       endif
    enddo
    
    call enter_exit(sub_name,2)
    
  end subroutine evaluate_ordering

!!!#############################################################################
  
  function inlist(item,ilist)
    
    integer :: item,ilist(:)
    ! Local variables
    integer :: n
    logical :: inlist

    ! --------------------------------------------------------------------------

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo
    
  end function inlist
  
!!!###############################################################
  
  subroutine make_plane_from_3points(NORML,NORMALTYPE,POINT1,POINT2,POINT3)
    
    !###    make_plane_from_3points finds the equation of a plane in three
    !###    dimensions and a vector normal to the plane from three
    !###    non-colinear points.
    !###    NORMALTYPE=1 for raw normal and plane equation
    !###    NORMALTYPE=2 for unit normal and plane equation
    !###    The coefficients represent aX + bY + cZ + d = 0
    !###    NORML(1)=a,NORML(2)=b,NORML(3)=c,NORML(4)=d
    
    
    !     Parameter list
    integer :: NORMALTYPE
    real(dp) :: POINT1(3),POINT2(3),POINT3(3),NORML(4)
    !     Local variables
    real(dp) :: DifF1(3),DifF2(3),NORMSIZE
    logical :: COLINEAR
    
    
    ! Check for colinearity
    COLINEAR=.FALSE.
    colinear = check_colinear_points(POINT1,POINT2,POINT3)
    if(.NOT.COLINEAR) then
       DifF1(1:3)=POINT2(1:3)-POINT1(1:3)
       DifF2(1:3)=POINT2(1:3)-POINT3(1:3)
       
       NORML(1)=(DifF1(2)*DifF2(3))-(DifF1(3)*DifF2(2))
       NORML(2)=(DifF1(3)*DifF2(1))-(DifF1(1)*DifF2(3))
       NORML(3)=(DifF1(1)*DifF2(2))-(DifF1(2)*DifF2(1))
       
       if(NORMALTYPE.EQ.2) then
          NORMSIZE = vector_length(NORML)
          NORML(1:3)=NORML(1:3)/NORMSIZE
       endif

       NORML(4) = -scalar_product_3(NORML,POINT1)
       
    else !Colinear
       
       WRITE(*,*) ' COLINEAR points in make_plane_from_3points '
       NORML = 0.0_dp
    endif
  end subroutine make_plane_from_3points
  
!!!##################################################

  subroutine scale_mesh(scaling,type)

    real(dp),intent(in) :: scaling
    character(len=2),intent(in) :: type

    select case(type)
    case('1d')
       node_xyz = node_xyz * scaling
    case('2d')
       node_xyz_2d = node_xyz_2d * scaling
       node_xyz_2d(4,:,:,:) = 0.0_dp
    end select
    
    scale_factors_2d = 1.0_dp

  end subroutine scale_mesh

!!!##################################################

  subroutine set_linear_derivatives

    ! Local variables
    integer :: ne,nk,nl,nline,np1,np2,nv1,nv2

!    do nline = 1,num_lines_2d                ! loop over all lines
!       nl = lines_2d(nline)                  ! the line number 
!       nk = nodes_in_line(1,0,nl)+1          ! the derivative = xi direction+1
!       np1 = nodes_in_line(2,1,nl)           ! the first node
!       nv1 = line_versn_2d(1,1,nl)           ! the version of the node for this line
!       np2 = nodes_in_line(3,1,nl)           ! the second node
!       node_xyz_2d(nk,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
!            - node_xyz_2d(1,1,:,np1)
!    enddo !nline
    do ne = 1,num_elems_2d
       np1 = elem_nodes_2d(1,ne)
       nv1 = elem_versn_2d(1,ne)
       np2 = elem_nodes_2d(2,ne)
       nv2 = elem_versn_2d(2,ne)
       node_xyz_2d(2,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
            - node_xyz_2d(1,1,:,np1)

!       write(*,'(''deriv1: ne='',i4,'' np1(nv1)='',i4,''('',i1,'') np2(nv2)='',i4,''('',i1,'')  length='',f5.1)') &
!            elems_2d(ne),nodes_2d(np1),nv1,nodes_2d(np2),nv2, &
!            node_xyz_2d(1,1,1,np2) - node_xyz_2d(1,1,1,np1)

       np2 = elem_nodes_2d(3,ne)
       nv2 = elem_versn_2d(3,ne)
       node_xyz_2d(3,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
            - node_xyz_2d(1,1,:,np1)

!       write(*,'(''deriv2: ne='',i4,'' np1(nv1)='',i4,''('',i1,'')  np2(nv2)='',i4,''('',i1,'')  length='',f5.1)') &
!            elems_2d(ne),nodes_2d(np1),nv1,nodes_2d(np2),nv2, &
!            node_xyz_2d(1,1,1,np2) - node_xyz_2d(1,1,1,np1)

       np1 = elem_nodes_2d(2,ne)
       nv1 = elem_versn_2d(2,ne)
       np2 = elem_nodes_2d(4,ne)
       nv2 = elem_versn_2d(4,ne)
       node_xyz_2d(3,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
            - node_xyz_2d(1,1,:,np1)

!       write(*,'(''deriv2: ne='',i4,'' np1(nv1)='',i4,''('',i1,'')  np2(nv2)='',i4,''('',i1,'')  length='',f5.1)') &
!            elems_2d(ne),nodes_2d(np1),nv1,nodes_2d(np2),nv2, &
!            node_xyz_2d(1,1,1,np2) - node_xyz_2d(1,1,1,np1)

       np1 = elem_nodes_2d(3,ne)
       nv1 = elem_versn_2d(3,ne)
       node_xyz_2d(2,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
            - node_xyz_2d(1,1,:,np1)

!       write(*,'(''deriv1: ne='',i4,'' np1(nv1)='',i4,''('',i1,'')  np2(nv2)='',i4,''('',i1,'')  length='',f5.1)') &
!            elems_2d(ne),nodes_2d(np1),nv1,nodes_2d(np2),nv2, &
!            node_xyz_2d(1,1,1,np2) - node_xyz_2d(1,1,1,np1)

       
    enddo !ne
    
  end subroutine set_linear_derivatives
  
!!!##################################################
  
  function area_between_two_vectors(vect_a,vect_b)
    
    !### 
    
    real(dp),intent(in) :: vect_a(3),vect_b(3)
    real(dp) :: cross(3)
    real(dp) :: area_between_two_vectors
    
    ! area = 1/2 x magnitude of the cross-product of vectors a and b
    
    cross = cross_product(vect_a,vect_b)
    area_between_two_vectors = 0.5_dp * sqrt(dot_product(cross,cross))
    
  end function area_between_two_vectors
  

!!!##################################################
  
  function area_between_three_points(point_a,point_b,point_c)
    
    !### 
    
    real(dp),intent(in) :: point_a(3),point_b(3),point_c(3)
    real(dp) :: norm(3),vect_a(3),vect_b(3)
    real(dp) :: area_between_three_points
    
    ! area = 1/2 x magnitude of the cross-product of vectors a and b
    
    vect_a(1:3) = point_a(1:3) - point_b(1:3)
    vect_b(1:3) = point_a(1:3) - point_c(1:3)
    norm = cross_product(vect_a,vect_b)
    area_between_three_points = 0.5_dp * sqrt(dot_product(norm,norm))
    
  end function area_between_three_points
  

!!!#############################################################################

  function get_local_elem(ne_global)

!!! dummy arguments
    integer,intent(in) :: ne_global
!!! local variables
    integer :: ne
    integer :: get_local_elem

    ! --------------------------------------------------------------------------

    do ne=1,num_elems_2d
       if(ne_global.eq.elems_2d(ne)) get_local_elem = ne
    enddo

  end function get_local_elem

!!!#############################################################################

  function get_local_elem_2d(ne_global)

    integer,intent(in) :: ne_global
    ! Local variables
    integer :: ne
    integer :: get_local_elem_2d

    ! --------------------------------------------------------------------------

    get_local_elem_2d = 0
    do ne = 1,num_elems_2d
       if(ne_global.eq.elems_2d(ne)) get_local_elem_2d = ne
    enddo

  end function get_local_elem_2d

!!!#############################################################################
  
  function get_local_node_f(ndimension,np_global) result(get_local_node)
    
    integer,intent(in) :: ndimension,np_global
    ! Local variables
    integer :: np
    integer :: get_local_node
    logical :: found

    ! --------------------------------------------------------------------------

    found = .false.
    np = 0
    
    select case (ndimension)
       
    case(1)
       do while (.not.found)
          np=np+1
          if(np.gt.num_nodes) then
             found = .true.
             write(*,'('' Warning: local node not found for global node'',I6)') np_global
          endif
          if(np_global.eq.nodes(np))then
             get_local_node = np
             found = .true.
          endif
       enddo
       
    case(2)
       do while (.not.found)
          np=np+1
          if(np.gt.num_nodes_2d) then
             found = .true.
             write(*,'('' Warning: local node not found for global node'',I6)') np_global
             read(*,*)
          endif
          if(np_global.eq.nodes_2d(np))then
             get_local_node = np
             found = .true.
          endif
       enddo
       
    end select
    
  end function get_local_node_f
  
!!!#############################################################################

  subroutine group_elem_by_parent(ne_parent,elemlist)
    !*group_elem_by_parent:* group elements that sit distal to a given
    ! parent element (ne_parent)

    integer,intent(in) :: ne_parent  ! the parent element number
    integer :: elemlist(:)
    ! Local Variables
    integer :: nt_bns,ne_count,num_nodes,m,n,ne0
    integer,allocatable :: ne_old(:),ne_temp(:)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name='group_elem_by_parent'
    call enter_exit(sub_name,1)
    
    elemlist = 0
    allocate(ne_old(size(elemlist)))
    allocate(ne_temp(size(elemlist)))
    
    nt_bns=1
    ne_old(1) = ne_parent
    ne_count = 1
    elemlist(ne_count)=ne_parent
    
    do while(nt_bns.ne.0)
       num_nodes=nt_bns
       nt_bns=0
       do m=1,num_nodes
          ne0=ne_old(m) !parent global element number
          do n=1,elem_cnct(1,0,ne0) !for each daughter branch
             nt_bns=nt_bns+1
             ne_temp(nt_bns)=elem_cnct(1,n,ne0)
          enddo !n
       enddo !m
       do n=1,nt_bns
          ne_old(n)=ne_temp(n) !updates list of previous generation element numbers
          ne_count=ne_count+1
          elemlist(ne_count)=ne_temp(n)
       enddo !n
    enddo !while
    
    deallocate(ne_old)
    deallocate(ne_temp)
    
    call enter_exit(sub_name,2)
    
  end subroutine group_elem_by_parent

!!!#############################################################################

  subroutine group_elem_parent_term(ne_parent,parent_term_list)
    !*group_elem_parent_term:* group the terminal elements that sit distal to
    !  a given parent element (ne_parent)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROUP_ELEM_PARENT_TERM" :: GROUP_ELEM_PARENT_TERM
    
    integer,intent(in) :: ne_parent  ! the parent element number
    integer :: parent_term_list(:)
    ! Local Variables
    integer :: ne,ne_count,noelem,num_parents
    integer,allocatable :: templist(:)
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'group_elem_parent_term'
    call enter_exit(sub_name,1)
    
    allocate(templist(num_elems))
    
    !reset the list of parent elements to zero
    call group_elem_by_parent(ne_parent,templist)
    
    ne_count=count(templist.ne.0)
    num_parents=0
    parent_term_list=0
    
    do noelem=1,ne_count
       ne = templist(noelem)
       if(elem_cnct(1,0,ne).eq.0)then
          num_parents=num_parents+1
          parent_term_list(num_parents)=ne
       endif !elem_cnct
    enddo !noelem
    
    deallocate(templist)
    
    call enter_exit(sub_name,2)
    
  end subroutine group_elem_parent_term
  
!!! ##########################################################################      

  function hermite(i,j,k,xi)
    !*hermite*: evaluates the cubic Hermite basis function at xi. the function
    ! is returned for values(k=1) or first(k=2) or second(k=3) derivatives, at
    ! xi=0 (i=1) or xi=1 (i=2). 
    
    integer,intent(in) :: i,j,k
    real(dp) :: xi
!!! local variables
    integer :: i_j_k
    real(dp) :: hermite
    
    ! K is 1,2, or 3; J is 1 or 2; I is 1 or 2
    
    i_j_k = 100*i + 10*j + k
    
    select case(i_j_k)
    case(111) !i=1,j=1,k=1
       hermite = 1.0_dp - 3.0_dp*xi**2 + 2.0_dp*xi**3  ! phi_10 = 1 -3xi^2 + 2xi^3
    case(121) !i=1,j=2,k=1
       hermite = xi*(xi - 1.0_dp)**2                   ! phi_11 = xi(xi - 1)^2
    case(211) !i=2,j=1,k=1
       hermite = xi**2 *(3.0_dp - 2.0_dp*xi)           ! phi_20 = xi^2(3 - 2xi)
    case(221) !i=2,j=2,k=1
       hermite = xi**2 *(xi - 1.0_dp)                  ! phi_21 = xi^2(xi - 1)
    case(112) !i=1,j=1,k=2
       hermite = 6.0_dp*xi *(xi - 1.0_dp)              ! phi_10' = 6xi(xi - 1)
    case(122) !i=1,j=2,k=2
       hermite = 3.0_dp*xi**2 - 4.0_dp*xi + 1.0_dp     ! phi_11' = 3xi^2-4xi+1
    case(212) !i=2,j=1,k=2
       hermite = 6.0_dp*xi*(1.0_dp - xi)               ! phi_20' = 6xi-6xi^2
    case(222) !i=2,j=2,k=2
       hermite = xi*(3.0_dp*xi - 2.0_dp)               ! phi_21' = 3xi^2-2xi
    case(113) !i=1,j=1,k=3
       hermite = 6.0_dp*(2.0_dp*xi - 1.0_dp)           ! phi_10'' = 12xi-6
    case(123) !i=1,j=2,k=3
       hermite = 6.0_dp*xi - 4.0_dp                    ! phi_11'' = 6xi-4
    case(213) !i=2,j=1,k=3
       hermite = 6.0_dp - 12.0_dp*xi                   ! phi_20'' = -12xi+6
    case(223) !i=2,j=2,k=3
       hermite = 6.0_dp*xi - 2.0_dp                    ! phi_21'' = 6xi-2
    end select
    
  end function hermite

!!! ##########################################################################      

  function linear(i,k,xi)
    
    integer,intent(in) :: i,k
    real(dp),intent(in) :: xi
!!! local variables
    integer :: i_k
    real(dp) :: linear
    
    i_k = 10*i + k
    
    select case(I_K)
    case(11) !i=1,k=1
       linear = 1.0_dp-xi                              ! phi_10 = 1-xi
    case(21) !i=2,k=1
       linear = xi                                     ! phi_20 = xi
    case(12) !i=1,k=2
       linear = -1.0_dp                                ! phi_10' = -1
    case(22) !i=2,k=2
       linear = 1.0_dp                                 ! phi_20' = 1
    case(30 :) !k=3
       linear = 0.0_dp                                 ! phi_10''= phi_20'' = 0
    end select
    
    return
  end function linear

!!!#############################################################################

  subroutine triangles_from_surface(num_triangles,num_vertices,surface_elems, &
       triangle,vertex_xyz)
    !*triangles_from_surface:* generates a linear surface mesh of triangles
    ! from an existing high order surface mesh. 
    
    integer :: num_triangles,num_vertices
    integer,intent(in) :: surface_elems(:)
    integer,allocatable :: triangle(:,:)
    real(dp),allocatable :: vertex_xyz(:,:)
    ! Local variables
    integer,parameter :: ndiv = 4
    integer :: i,index1,index2,j,ne,nelem,nmax_1,nmax_2,num_surfaces, &
         num_tri_vert,nvertex_row,step_1,step_2
    real(dp) :: X(3),xi(3)
    logical :: four_nodes
    character(len=3) :: repeat
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'triangles_from_surface'
    call enter_exit(sub_name,1)
    
    if(allocated(triangle)) deallocate(triangle)
    allocate(triangle(3,2*num_elems_2d*ndiv**2))
    if(allocated(vertex_xyz)) deallocate(vertex_xyz)
    allocate(vertex_xyz(3,num_elems_2d*(ndiv+1)**2))
    
    triangle = 0
    vertex_xyz = 0.0_dp
    num_surfaces = count(surface_elems.ne.0)
    num_triangles = 0
    num_vertices = 0
    num_tri_vert = 0 

    do nelem = 1,num_surfaces
       ne = surface_elems(nelem)
       four_nodes = .false.
       repeat = '0_0'
       if(elem_nodes_2d(1,ne).eq.elem_nodes_2d(2,ne)) repeat = '1_0'
       if(elem_nodes_2d(1,ne).eq.elem_nodes_2d(3,ne)) repeat = '2_0'
       if(elem_nodes_2d(2,ne).eq.elem_nodes_2d(4,ne)) repeat = '2_1'
       if(elem_nodes_2d(3,ne).eq.elem_nodes_2d(4,ne)) repeat = '1_1'

       select case(repeat)
       case ('0_0')
        
          nmax_1 = ndiv+1 ! ndiv+1 vertices in xi1 direction
          step_1 = 0      ! # of vertices in xi1 is constant
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
          step_2 = 0      ! # of vertices in xi2 is constant
          index1 = 1
          index2 = 2
          four_nodes = .true.
          
       case ('1_0')
          
          nmax_1 = 1      ! start with 1 vertex in xi1 direction
          step_1 = 1      ! increase # of vertices in xi1 with each step in xi2
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
          step_2 = 0      ! # of vertices in xi2 is constant
          index1 = 1
          index2 = 2
          
       case ('1_1')
          
          nmax_1 = ndiv+1 ! start with ndiv+1 vertices in xi1 direction
          step_1 = -1     ! decrease # of vertices in xi1 with each step in xi2
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
          step_2 = 0      ! # of vertices in xi2 is constant
          index1 = 1
          index2 = 2
          
       case ('2_0')
          
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi1 direction
          step_2 = 0      ! # of vertices in xi1 is constant
          nmax_1 = 1      ! start with 1 vertex in xi2 direction
          step_1 = 1      ! increase # of vertices in xi2 with each step in xi1          
          index2 = 1
          index1 = 2
          
       case ('2_1')
          
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi1 direction
          step_2 = 0      ! # of vertices in xi1 is constant
          nmax_1 = ndiv+1 ! start with ndiv+1 vertices in xi2 direction
          step_1 = -1     ! decrease # of vertices in xi2 with each step in xi1
          index2 = 1
          index1 = 2
       end select

       xi(index2) = 0.0_dp
       do i = 1,nmax_2
          xi(index1) = 0.0_dp
          do j = 1,nmax_1
             num_vertices = num_vertices + 1
             X = coord_at_xi(ne,xi,'hermite')
             vertex_xyz(1:3,num_vertices) = X(1:3)
             if(nmax_1.gt.1) xi(index1) = xi(index1) + 1.0_dp/(nmax_1-1)
             if(i.gt.1.and.j.gt.1)then
                num_triangles = num_triangles + 1
                triangle(1,num_triangles) = num_vertices
                triangle(2,num_triangles) = num_vertices-1
                triangle(3,num_triangles) = nvertex_row+j-1
                if(four_nodes.or.(.not.four_nodes.and.j.lt.nmax_1).or.step_1.eq.-1)then
                   num_triangles = num_triangles + 1
                   triangle(1,num_triangles) = num_vertices
                   triangle(2,num_triangles) = nvertex_row+j
                   triangle(3,num_triangles) = nvertex_row+j-1
                endif
                if(step_1.eq.-1.and.j.eq.nmax_1)then
                   num_triangles = num_triangles + 1
                   triangle(1,num_triangles) = num_vertices
                   triangle(2,num_triangles) = nvertex_row+j+1
                   triangle(3,num_triangles) = nvertex_row+j
                endif
             else if(step_1.eq.-1.and.i.eq.nmax_2.and.j.eq.1)then
                num_triangles = num_triangles + 1
                triangle(1,num_triangles) = num_vertices
                triangle(2,num_triangles) = num_vertices-1
                triangle(3,num_triangles) = num_vertices-2
             endif
          enddo !j
          nvertex_row = num_vertices-nmax_1 !record first vertex # in previous row
          if(nmax_2.gt.1) xi(index2) = xi(index2) + 1.0_dp/(nmax_2-1)
          nmax_1 = nmax_1 + step_1
       enddo !i
    enddo
    
!    write(*,'('' Made'',I8,'' triangles to cover'',I6,'' surface elements'')') &
!         num_triangles,num_surfaces !num_elems_2d
    
    call enter_exit(sub_name,2)
    
  end subroutine triangles_from_surface

!!!##################################################
  
  function unit_norm_to_plane_two_vectors(vect_a,vect_b)
    
    real(dp),intent(in) :: vect_a(3),vect_b(3)
    real(dp) :: magnitude,norm(3)
    real(dp) :: unit_norm_to_plane_two_vectors(3)
    
    norm = cross_product(vect_a,vect_b)
    magnitude = sqrt(dot_product(norm,norm))
    unit_norm_to_plane_two_vectors = norm/magnitude
    
  end function unit_norm_to_plane_two_vectors
  

!!!##################################################
  
  function unit_norm_to_three_points(point_a,point_b,point_c)

    real(dp),intent(in) :: point_a(3),point_b(3),point_c(3)
    real(dp) :: magnitude,norm(3),vect_a(3),vect_b(3)
    real(dp) :: unit_norm_to_three_points(3)

    vect_a(1:3) = point_a(1:3) - point_b(1:3)
    vect_b(1:3) = point_a(1:3) - point_c(1:3)
    norm = cross_product(vect_a,vect_b)
    magnitude = sqrt(dot_product(norm,norm))
    unit_norm_to_three_points = norm/magnitude

  end function unit_norm_to_three_points

!!!##################################################
  

  function angle_btwn_points(A,B,C)
    
    !###    calculates the angle between three points
    
    real(dp),intent(in) :: A(3),B(3),C(3)

    real(dp) :: U(3),V(3)
    real(dp) :: angle_btwn_points

    U = A - B
    V = C - B
    angle_btwn_points = angle_btwn_vectors(U,V)
        
  end function angle_btwn_points
  
!!!##################################################

  function angle_btwn_vectors(U,V)
    
    !###    ANGLE calculates the angle between two vectors
    
    real(dp),intent(in) :: U(3),V(3)

    real(dp) :: ANGLE,angle_btwn_vectors,N_U(3),N_V(3)
    
    N_U = unit_vector(U)
    N_V = unit_vector(V)

    ANGLE = scalar_product_3(N_U,N_V)
    ANGLE = max(-1.0_dp,ANGLE)
    ANGLE = min(1.0_dp,ANGLE)
    ANGLE = acos(ANGLE)

    angle_btwn_vectors=ANGLE
    
  end function angle_btwn_vectors
  
!!!###############################################################
  
  function check_colinear_points(POINT1,POINT2,POINT3)
    
    !###    check_colinear_points checks whether two vectors are colinear.
     
    !     Parameter list
    real(dp) :: POINT1(3),POINT2(3),POINT3(3)
    !     Local variables
    real(dp) :: ERR1(3),ERR2(3),LU,LV,U(3),V(3)
    logical :: check_colinear_points
    
    
    check_colinear_points =.FALSE.
    U(1:3)=POINT2(1:3)-POINT1(1:3)
    V(1:3)=POINT3(1:3)-POINT1(1:3)
    LU = vector_length(U)
    LV = vector_length(V)
    ! If 2 of the points are the same then LU and LV
    ! can be zero causing div by zero below and resulting in
    ! the wrong answer (on Linux) 
    if((abs(LU)>zero_tol).AND.(abs(LV)>zero_tol)) then
       ERR1(1:3)=abs(U(1:3)/LU-V(1:3)/LV)
       ERR2(1:3)=abs(U(1:3)/LU+V(1:3)/LV)
       if((ERR1(1).LE.ZERO_TOL.AND.ERR1(2).LE.ZERO_TOL.AND.ERR1(3).LE. &
            ZERO_TOL).OR.(ERR2(1).LE.ZERO_TOL.AND.ERR2(2).LE.ZERO_TOL.AND. &
            ERR2(3).LE.ZERO_TOL)) check_colinear_points=.TRUE.
    else
       check_colinear_points=.TRUE.
    endif
  end function check_colinear_points
  

!!!#############################################################################
  
  function coord_at_xi(ne,xi,basis)
    
    integer,intent(in) :: ne
    real(dp),intent(in) :: xi(:)
    character(len=*),intent(in) :: basis
    ! Local Variables
    integer :: nn,nv
    real(dp) :: phi(4),phi_10(2),phi_11(2),phi_20(2),phi_21(2),x(4,3,4)
    real(dp) :: coord_at_xi(3)

    ! --------------------------------------------------------------------------
    
    select case (basis)
    case('linear')
       forall (nn=1:4) x(1,1:3,nn) = node_xyz_2d(1,1,1:3,elem_nodes_2d(nn,ne))
       phi(1) = (1.0_dp - xi(1))*(1.0_dp - xi(2))
       phi(2) = xi(1)*(1.0_dp - xi(2))
       phi(3) = (1.0_dp - xi(1))*xi(2)
       phi(4) = xi(1)*xi(2)
       
       coord_at_xi(1:3) = phi(1)*x(1,1:3,1)+phi(2)*x(1,1:3,2)+phi(3)* &
            x(1,1:3,3)+phi(4)*x(1,1:3,4)
       
    case('hermite')
       do nn=1,4
          nv = elem_versn_2d(nn,ne)
          x(1:4,1:3,nn) = node_xyz_2d(1:4,nv,1:3,elem_nodes_2d(nn,ne))
       enddo
       phi_10(1) = (2.0_dp*xi(1)-3.0_dp)*xi(1)*xi(1)+1.0_dp  ! 2xi^3-3xi^2+1
       phi_10(2) = (2.0_dp*xi(2)-3.0_dp)*xi(2)*xi(2)+1.0_dp  ! 2xi^3-3xi^2+1
       phi_20(1) = xi(1)*xi(1)*(3.0_dp-2.0_dp*xi(1))         ! -2xi^3+3xi^2
       phi_20(2) = xi(2)*xi(2)*(3.0_dp-2.0_dp*xi(2))         ! -2xi^3+3xi^2
       phi_11(1) = ((xi(1)-2.0_dp)*xi(1)+1.0_dp)*xi(1)       ! xi^3-2xi^2+xi
       phi_11(2) = ((xi(2)-2.0_dp)*xi(2)+1.0_dp)*xi(2)       ! xi^3-2xi^2+xi
       phi_21(1) = xi(1)*xi(1)*(xi(1)-1.0_dp)                ! xi^3-xi^2
       phi_21(2) = xi(2)*xi(2)*(xi(2)-1.0_dp)                ! xi^3-xi^2
       coord_at_xi(1:3) = phi_10(1)*phi_10(2)*x(1,1:3,1) &
            + phi_20(1)*phi_10(2)*x(1,1:3,2) &
            + phi_10(1)*phi_20(2)*x(1,1:3,3) &
            + phi_20(1)*phi_20(2)*x(1,1:3,4) &
            + phi_11(1)*phi_10(2)*x(2,1:3,1) * scale_factors_2d(2,ne) &
            + phi_21(1)*phi_10(2)*x(2,1:3,2) * scale_factors_2d(6,ne) &
            + phi_11(1)*phi_20(2)*x(2,1:3,3) * scale_factors_2d(10,ne) &
            + phi_21(1)*phi_20(2)*x(2,1:3,4) * scale_factors_2d(14,ne) &
            + phi_10(1)*phi_11(2)*x(3,1:3,1) * scale_factors_2d(3,ne) &
            + phi_20(1)*phi_11(2)*x(3,1:3,2) * scale_factors_2d(7,ne) &
            + phi_10(1)*phi_21(2)*x(3,1:3,3) * scale_factors_2d(11,ne) &
            + phi_20(1)*phi_21(2)*x(3,1:3,4) * scale_factors_2d(15,ne) &
            + phi_11(1)*phi_11(2)*x(4,1:3,1) * scale_factors_2d(4,ne) &
            + phi_21(1)*phi_11(2)*x(4,1:3,2) * scale_factors_2d(8,ne) &
            + phi_11(1)*phi_21(2)*x(4,1:3,3) * scale_factors_2d(12,ne) &
            + phi_21(1)*phi_21(2)*x(4,1:3,4) * scale_factors_2d(16,ne)
    end select
    
  end function coord_at_xi
  
!!!###############################################################
  
  function cross_product(A,B)
    
    !###  cross_product returns the vector cross product of A*B in C.
    
    !     Parameter List
    real(dp),intent(in) :: A(3),B(3)

    real(dp) :: cross_product(3)
    
    cross_product(1) = A(2)*B(3)-A(3)*B(2)
    cross_product(2) = A(3)*B(1)-A(1)*B(3)
    cross_product(3) = A(1)*B(2)-A(2)*B(1)
    
  end function cross_product
  
!!! ##########################################################################   

  function direction_point_to_point(point_start,point_end)

    real(dp),intent(in) :: point_start(:),point_end(:)

    real(dp) :: vector(3)
    real(dp) :: direction_point_to_point(3)

    vector(1:3) = point_end(1:3) - point_start(1:3)
    vector(1:3) = unit_vector(vector)
    direction_point_to_point = vector

  end function direction_point_to_point

!!!###############################################################
  
  function scalar_triple_product(A,B,C)
    
    !###  scalar_triple_product returns A.(BxC)
    
    !     Parameter List
    real(dp),intent(in) :: A(3),B(3),C(3)

    real(dp) :: scalar_triple_product
    
    scalar_triple_product = A(1)*(B(2)*C(3)-B(3)*C(2)) + &
         A(2)*(B(3)*C(1)-B(1)*C(3)) + A(3)*(B(1)*C(2)-B(2)*C(1))
    
  end function scalar_triple_product
  
!!!###############################################################
  
  function distance_between_points(point1, point2)
    
    !###    calculates the distance between two arbitrary points
    
    real(dp),intent(in) :: point1(3),point2(3)
    integer :: i
    real(dp) :: distance_between_points
    
    distance_between_points = 0.0_dp
    do i=1,3
       distance_between_points = distance_between_points + (point1(i)-point2(i))**2
    enddo
    distance_between_points = sqrt(distance_between_points)
    
  end function distance_between_points
  
!!!###############################################################
  
  function mesh_a_x_eq_b(MATRIX,VECTOR)
    
    real(dp) :: MATRIX(3,3),VECTOR(3)
    !Local variables
    integer :: i,j,k,pivot_row
    real(dp) :: A(3,4),max,pivot_value,SOLUTION(3),TEMP(4)
    real(dp) :: mesh_a_x_eq_b(3)
    
    
    A(1:3,1:3) = MATRIX(1:3,1:3)
    A(1:3,4) = VECTOR(1:3)
    do k=1,2
       max=0.0_dp
       do i=k,3
          if(abs(A(i,k)).GT.max)then
             max=abs(A(i,k))
             pivot_row=i
          endif
       enddo !i
       if(pivot_row.ne.k)then
          do j=1,4
             TEMP(j)=A(k,j)
             A(k,j)=A(pivot_row,j)
             A(pivot_row,j)=TEMP(j)
          enddo !j
       endif
       pivot_value = A(k,k)
       A(k,1:4) = A(k,1:4)/pivot_value
       do i=k+1,3
          do j=k+1,4
             A(i,j) = A(i,j)-A(i,k)*A(k,j)
          enddo
          A(i,k) = 0.0_dp
       enddo
    enddo !N
    A(3,4) = A(3,4)/A(3,3)
    A(2,4) = A(2,4)-A(3,4)*A(2,3)
    A(1,4) = A(1,4)-A(3,4)*A(1,3)-A(2,4)*A(1,2)

    SOLUTION(1:3) = A(1:3,4)
    
    mesh_a_x_eq_b = solution

  end function mesh_a_x_eq_b
  
!!!#############################################################################
  
  subroutine reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    !*reallocate_node_elem_arrays:* Reallocates the size of geometric
    ! arrays when modifying geometries

    use indices

    integer,intent(in) :: num_elems_new,num_nodes_new
    ! Local variables
    integer,allocatable :: nodelem_temp(:),enodes_temp(:,:),enodes_temp2(:,:,:)
    real(dp),allocatable :: xyz_temp(:,:),rnodes_temp(:,:)
    logical,allocatable :: exp_temp(:)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'reallocate_node_elem_arrays'
    call enter_exit(sub_name,1)
    
    allocate(nodelem_temp(num_nodes))
    nodelem_temp = nodes ! copy to temporary array
    deallocate(nodes) !deallocate initially allocated memory
    allocate(nodes(num_nodes_new))
    nodes(1:num_nodes)=nodelem_temp(1:num_nodes)
    deallocate(nodelem_temp) !deallocate the temporary array
    
    allocate(xyz_temp(3,num_nodes))
    xyz_temp=node_xyz
    deallocate(node_xyz)
    allocate(node_xyz(3,num_nodes_new))
    node_xyz(1:3,1:num_nodes)=xyz_temp(1:3,1:num_nodes)
    
    allocate(nodelem_temp(num_elems))
    nodelem_temp = elems ! copy to temporary array
    deallocate(elems) !deallocate initially allocated memory
    allocate(elems(num_elems_new))
    elems(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array
    
    allocate(enodes_temp(2,num_elems))
    enodes_temp=elem_nodes
    deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems_new))
    elem_nodes(1:2,1:num_elems)=enodes_temp(1:2,1:num_elems)
    deallocate(enodes_temp)
    
    if(allocated(elem_field).and.num_ne.gt.0)then
       allocate(rnodes_temp(num_ne,num_elems))
       rnodes_temp=elem_field
       deallocate(elem_field)
       allocate(elem_field(num_ne,num_elems_new))
       elem_field(1:num_ne,1:num_elems)=rnodes_temp(1:num_ne,1:num_elems)
       deallocate(rnodes_temp)
       elem_field(1:num_ne,num_elems+1:num_elems_new) = 0.0_dp
    endif
    
    allocate(rnodes_temp(3,num_elems))
    rnodes_temp=elem_direction
    deallocate(elem_direction)
    allocate(elem_direction(3,num_elems_new))
    elem_direction(1:3,1:num_elems)=rnodes_temp(1:3,1:num_elems)
    deallocate(rnodes_temp)
    elem_direction(1:3,num_elems+1:num_elems_new) = 0.0_dp
    
    if(allocated(node_field).and.num_nj.gt.0)then
       allocate(rnodes_temp(num_nj,num_nodes))
       rnodes_temp=node_field
       deallocate(node_field)
       allocate(node_field(num_nj,num_nodes_new))
       node_field(1:num_nj,1:num_nodes)=rnodes_temp(1:num_nj,1:num_nodes)
       deallocate(rnodes_temp)
       node_field(1:num_nj,num_nodes+1:num_nodes_new)=0.0_dp
    endif
    
    allocate(nodelem_temp(num_elems))
    nodelem_temp = elem_symmetry ! copy to temporary array
    deallocate(elem_symmetry) !deallocate initially allocated memory
    allocate(elem_symmetry(num_elems_new))
    elem_symmetry(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array
    elem_symmetry(num_elems+1:num_elems_new)=1
    
    allocate(enodes_temp2(-1:1,0:2,0:num_elems))
    enodes_temp2=elem_cnct
    deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems_new))
    elem_cnct(-1:1,0:2,0:num_elems)=enodes_temp2(-1:1,0:2,0:num_elems)
    deallocate(enodes_temp2)
    elem_cnct(-1:1,0:2,num_elems+1:num_elems_new) = 0
    
    allocate(enodes_temp(num_ord,num_elems))
    enodes_temp=elem_ordrs
    deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems_new))
    elem_ordrs(1:num_ord,1:num_elems)=enodes_temp(1:num_ord,1:num_elems)
    deallocate(enodes_temp)
    elem_ordrs(1:num_ord,num_elems+1:num_elems_new) = 0
    
    if(allocated(elem_units_below).and.num_nu.gt.0)then
       allocate(nodelem_temp(num_elems))
       nodelem_temp=elem_units_below
       deallocate(elem_units_below)
       allocate(elem_units_below(num_elems_new))
       elem_units_below(1:num_elems)=nodelem_temp(1:num_elems)
       deallocate(nodelem_temp)
       elem_units_below(num_elems+1:num_elems_new)=0
    endif
    
    allocate(enodes_temp(num_nodes,0:3))
    enodes_temp=elems_at_node
    deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes_new,0:3))
    elems_at_node(1:num_nodes,0:3)=enodes_temp(1:num_nodes,0:3)
    deallocate(enodes_temp)
    elems_at_node(num_nodes+1:num_nodes_new,0:3)=0
    
    if(model_type.eq.'gas_mix')then
       allocate(exp_temp(num_elems))
       exp_temp = expansile
       deallocate(expansile)
       allocate(expansile(num_elems_new))
       expansile(1:num_elems)=exp_temp(1:num_elems)
       deallocate(exp_temp)
       expansile(num_elems+1:num_elems_new)=.false.
    endif
    
    call enter_exit(sub_name,2)
    
  end subroutine reallocate_node_elem_arrays
  
!!!##################################################
  
  function scalar_product_3(A,B)
    
    !### calculates scalar product of two vectors A,B of length 3.
    
    real(dp),intent(in) :: A(*),B(*)

    integer :: i
    real(dp) :: scalar_product_3
    
    scalar_product_3 = 0.0_dp
    do i=1,3
       scalar_product_3 = scalar_product_3 + A(i)*B(i)
    enddo
    
  end function scalar_product_3
  
!!!###############################################################
  
  function unit_vector(A)
    
    !###  Calculates the unit vector for an arbitrary 3x1 vector 
    
    real(dp),intent(in) :: A(*)
    real(dp) :: length_a,unit_vector(3)

    length_a = vector_length(A)
    if(length_a.gt.1.0e-6_dp)then
       unit_vector(1:3) = A(1:3)/length_a
    else
       WRITE(*,*) ' >>WARNING: Cannot normalise a zero length vector'
       WRITE(*,*) ' We recommend debugging, but hit enter to continue'
       read(*,*)
    endif

  end function unit_vector

!!!##################################################
  
  function vector_length(A)
    
    !###  Calculates the length of a 3x1 vector 
    
    real(dp),intent(in) :: A(*)
    real(dp) :: vector_length
    integer :: i
    
    vector_length = 0.0_dp
    do i=1,3
       vector_length = vector_length + A(i)*A(i)
    enddo
    vector_length = sqrt(vector_length)
    
  end function vector_length
  
!!!###############################################################

  function volume_internal_to_surface(triangles,vertex_xyz)

    ! calculates the volume enclosed by a list of surface elements

    integer,intent(in) :: triangles(:,:)
    real(dp),intent(in) :: vertex_xyz(:,:)
    real(dp) :: volume_internal_to_surface

!!! Local Variables
    integer :: ntri,num_triangles
    real(dp) :: volume,V1(3),V2(3),V3(3),P4(3)

    num_triangles = count(triangles(:,:).ne.0)/3

    P4 = sum(vertex_xyz,dim=2)/size(vertex_xyz,dim=2)

    volume = 0.0_dp

    do ntri = 1,num_triangles
       V1(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(1,ntri))
       V2(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(2,ntri))
       V3(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(3,ntri))
       volume = volume + abs(scalar_triple_product(V1,V2,V3))
    enddo

    volume_internal_to_surface = volume/6.0_dp

  end function volume_internal_to_surface

!!!###############################################################

  function point_internal_to_surface(num_vertices,triangles,point_xyz,vertex_xyz)
!!! Cast a line in positive x-direction from each data point and 
!!! then work out how many triangular elements it crosses. If even it is in the 
!!! shape and if odd it is outside the shape

    integer,intent(in) :: num_vertices,triangles(:,:)
    real(dp),intent(in) :: point_xyz(3),vertex_xyz(:,:)
    logical :: point_internal_to_surface

!!! Local Variables
    integer :: i,nclosest,ncrossed,ntri,num_triangles
    real(dp) :: area,area_triangle,cofm_surfaces(3),denominator,&
         dist,min_dist,norm_v(3),point(3),P1(3),P2(3),P3(3),u
    real(dp),parameter :: dist_tol = 1.0e-4_dp, user_tol = 1.0e-14_dp
    logical :: cross_any

    num_triangles = count(triangles(:,:).ne.0)/3

    forall (i=1:3) cofm_surfaces(i) = sum(vertex_xyz(i,1:num_vertices))/num_vertices

! check whether the line that joins the centre of mass of the surface mesh and the point
! in question crosses ANY face. If it does, then point not inside.

    cross_any = .false.
    ncrossed = 0

    do ntri = 1,num_triangles
       P1(1:3) = vertex_xyz(1:3,triangles(1,ntri))
       P2(1:3) = vertex_xyz(1:3,triangles(2,ntri))
       P3(1:3) = vertex_xyz(1:3,triangles(3,ntri))
       norm_v = unit_norm_to_three_points(P1,P2,P3) ! unit normal to triangle plane
       ! u = (a*x1+b*y1+c*z1+d)/(a*(x1-x2)+b*(y1-y2)+c*(z1-z2))
       denominator = norm_v(1)*(point_xyz(1)-cofm_surfaces(1)) + &
            norm_v(2)*(point_xyz(2)-cofm_surfaces(2)) + &
            norm_v(3)*(point_xyz(3)-cofm_surfaces(3))
       ! denominator is zero for line parallel to plane
       if(abs(denominator).gt.user_tol)then
          ! calculate the distance of the surface point from point_xyz
          u = (dot_product(norm_v,point_xyz)-dot_product(norm_v,P1))/denominator
          if(u.ge.0.0_dp.and.u.le.1.0_dp)then ! POTENTIALLY crosses. Test further (angle)
             point = point_xyz + u*(cofm_surfaces-point_xyz) ! projection to surface
             area = area_between_two_vectors(P1-point,P2-point)+ &
                  area_between_two_vectors(P1-point,P3-point)+area_between_two_vectors(P2-point,P3-point)
             area_triangle = area_between_two_vectors(P1-P2,P1-P3)
             if(abs(area_triangle-area).lt.dist_tol)then
                cross_any = .true.
                ncrossed = ncrossed + 1
             endif
          endif
       endif
    enddo
    
    if(.not.cross_any)then
       point_internal_to_surface = .true.
    else
       if(ncrossed.eq.2)then
          point_internal_to_surface = .true.
       else
          point_internal_to_surface = .false.
       endif
    endif

  end function point_internal_to_surface


!!!#############################################################################
  
  function terminal_element(ne)
    !*terminal element:* returns 'true' if a 1d element has no elements adjacent
    ! in the Xi+1 direction
    integer,intent(in) :: ne
    logical :: terminal_element

    ! --------------------------------------------------------------------------

    if(elem_cnct(1,0,ne).eq.0)then
       terminal_element = .true.
    else
       terminal_element = .false.
    endif
    
  end function terminal_element

!!!#############################################################################
  
  function stem_element(ne)
    !*stem element:* returns 'true' if a 1d element has no elements adjacent
    ! in the Xi-1 direction
    integer,intent(in) :: ne
    logical :: stem_element

    ! --------------------------------------------------------------------------

    if(elem_cnct(-1,0,ne).eq.0)then
       stem_element = .true.
    else
       stem_element = .false.
    endif
    
  end function stem_element

!!!#############################################################################
  
  function bifurcation_element(ne)
    !*bifurcation element:* returns 'true' if a 1d element has two elements 
    ! adjacent in the Xi+1 direction (i.e. parent of a bifurcation)
    integer,intent(in) :: ne
    logical :: bifurcation_element

    ! --------------------------------------------------------------------------

    if(ne.eq.0)then
       bifurcation_element = .false.
    elseif(elem_cnct(1,0,ne).eq.2)then
       bifurcation_element = .true.
    else
       bifurcation_element = .false.
    endif
    
  end function bifurcation_element

!!!#############################################################################
  
  function which_child(ne,ne0)
    !*which child:* returns '1' if ne is recorded as the first child element of
    ! element ne0, and '2' if the second
    integer,intent(in) :: ne,ne0
    integer :: which_child

    ! --------------------------------------------------------------------------

    if(elem_cnct(1,1,ne0).eq.ne) then
       which_child = 1
    elseif(elem_cnct(1,2,ne0).eq.ne) then
       which_child = 2
    else
       write(*,'('' Warning! element'',i6,'' is not a child element of'',i6)') &
            ne,ne0
       which_child = 0
    endif
    
  end function which_child

!!!#############################################################################
  
end module mesh_utilities

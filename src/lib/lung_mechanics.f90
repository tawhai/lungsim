module lung_mechanics

  ! this module is based on legacy code from CMISS (Continuum Mechanics, Image analysis, Signal processing and System identification).!
  ! [see www.cmiss.org]. it provides convenient and quick functionality for simulating soft tissue deformation of a hyper-elastic
  ! compressible body with an isotropic material law inside a contact surface. This is strictly 3D rectangular cartesian with 
  ! quadratic Lagrange basis functions. 
  ! simulation packages with more generalised capability should be used if additional functionality is needed.

  ! most of the original subroutine and array names are retained for ease of comparison with the legacy code
  
  use arrays
  use precision
  
  implicit none
  
!  private
  !  public
  
  integer,parameter :: n_basis = 2         ! # of basis functions for lung, cavity
  integer,parameter :: n_dirn = 3          ! max # of xi directions
  integer,parameter :: n_gauss = 3         ! # of Gauss points in a direction; maximum 3
  integer,parameter :: n_points = 27       ! max # of nodes in lung elements
  integer,parameter :: nb_lung = 1,nb_cavity = 2, nind_lung = 1, nind_cavity = 2

  integer :: inp(n_points,n_dirn,n_basis)
  real(dp) :: pg(n_points,11,n_gauss**3) ! for this problem
  real(dp) :: wg(n_gauss**3)
  real(dp) :: xig(n_dirn,n_gauss**3)
  
contains

!!!##############################################################################################

  subroutine problem_setup

    integer :: mm,nn
    integer :: i,j,k

    nit(nind_lung) = 3            ! 3 xi coordinates in lung region
    nut(nind_lung) = 11
    ngt(nind_lung) = n_gauss**3   ! # Gauss points in lung region
    nnt(nind_lung) = 27           ! # nodes in lung region
    nst(nind_lung) = 27

    nit(nind_cavity) = 2          ! 2 xi coordinates in cavity region
    nut(nind_cavity) = 6
    ngt(nind_cavity) = n_gauss**2 ! # Gauss points in cavity region
    nnt(nind_cavity) = 9          ! # nodes in cavity region
    nst(nind_cavity) = 9

    inp = 1
    nn = 0
    mm = 0
    do k = 1,3
       do j = 1,3
          do i = 1,3
             nn = nn + 1
             inp(nn,1,nind_lung) = i
             inp(nn,2,nind_lung) = j
             inp(nn,3,nind_lung) = k
          enddo
          mm = mm + 1
          inp(mm,1,nind_cavity) = j
          inp(mm,2,nind_cavity) = k
       enddo
    enddo
    
    call setup_gauss_quad(nind_lung)

  end subroutine problem_setup
  
!!!##############################################################################################
  
  subroutine setup_gauss_quad(nb)

    ! (CMISS GAUSS1) Define Gaussian quadrature coordinates xig and weights wg.
    ! Evaluate basis function Gauss point array pg
    
    !     Parameter List
    integer :: nb
    !     Local Variables
    integer :: I,J,K,ng,ng1,ng2,ng3,nitb,nn,nu
    real(dp) :: d(3,3),w(3,3),xi(3),xigg(3,3,3,3)

    data w/ 1.0_dp, 0.0_dp, 0.0_dp, 0.50_dp, 0.5_dp, 0.0_dp, 0.2777777777777780_dp, &
         0.4444444444444440_dp, 0.2777777777777780_dp/

    data d/ 0.0_dp, 0.0_dp, 0.0_dp, -0.2886751345948130_dp, 0.2886751345948130_dp, 0.0_dp, &
         -0.3872983346207410_dp, 0.0_dp, 0.3872983346207410_dp/
    
    nitb = nit(nb)

    do k = 1,n_gauss
       do j = 1,n_gauss
          do i = 1,n_gauss
             xigg(i,j,k,1) = 0.50_dp+d(i,n_gauss)
             xigg(i,j,k,2) = 0.50_dp+d(j,n_gauss)
             xigg(i,j,k,3) = 0.50_dp+d(k,n_gauss)
             ng = I+(J-1+(K-1)*n_gauss)*n_gauss
             wg(ng) = w(I,n_gauss)*w(J,n_gauss)*w(K,n_gauss)
             xi(1:nitb) = xigg(I,J,K,1:nitb)
             xig(1:nitb,ng) = xi(1:nitb)
             do nn = 1,nnt(nb)
                do nu = 1,nut(nb)
                   pg(nn,nu,ng) = psi1(nb,nu,nn,xi)
                enddo !nu
             enddo !nn
          enddo !i
       enddo !j
    enddo !k
    
  end subroutine setup_gauss_quad

!!!##############################################################################################

  subroutine define_node_tissue(filename)
    character :: filename*(*)
    integer :: nonode,np
    character :: readfile*(100)

    if(index(filename, ".lsnode")>0)then ! full filename is given
       readfile = trim(filename)
    else! append correct extension
       readfile = trim(filename)//'.lsnode'
    endif
    
    open(10, file = readfile, status = 'old')
    
    read(unit = 10, fmt =  *) tissue_num_nodes

    if(allocated(tissue_nodes)) deallocate(tissue_nodes)
    allocate(tissue_nodes(tissue_num_nodes))
    if(allocated(tissue_xyz)) deallocate(tissue_xyz)
    allocate(tissue_xyz(3,tissue_num_nodes))
    allocate(xp(njm,tissue_num_nodes))

    !npm = tissue_num_nodes
       
    do nonode = 1,tissue_num_nodes
       read(unit = 10, fmt =  *) np,xp(1:3,nonode)
       tissue_xyz(1:3,nonode) = xp(1:3,nonode) ! different storage
       tissue_nodes(nonode) = np
    enddo                     ! np
      
    close(10)
    
  end subroutine define_node_tissue

!!!##############################################################################################

  subroutine define_node_cavity(filename)
    character :: filename*(*)
    integer :: nonode,np
    real(dp) :: coords(3)
    character :: readfile*(100)

    if(index(filename, ".lsnode")>0)then ! full filename is given
       readfile = trim(filename)
    else! append correct extension
       readfile = trim(filename)//'.lsnode'
    endif
    
    open(10, file = readfile, status = 'old')
    
    read(unit = 10, fmt =  *) cavity_num_nodes

    if(allocated(cavity_nodes)) deallocate(cavity_nodes)
    allocate(cavity_nodes(cavity_num_nodes))
    allocate(xp_cavity(njm,cavity_num_nodes))
       
    do nonode = 1,cavity_num_nodes
       read(unit = 10, fmt =  *) np,coords(1:3) !xp(1,np),xp(2,np),xp(3,np)
       xp_cavity(1:3,nonode) = coords(1:3)
       cavity_nodes(nonode) = np
    enddo
      
    close(10)

  end subroutine define_node_cavity
  
!!!##############################################################################################

  subroutine define_elem_tissue(filename)
    character :: filename*(*)
    !     Local Variables
    integer :: centre_nodes(6),i,ne,nodelist(27),noelem,noelem2, &
         nn,np,npc,num_repeats
    integer,allocatable :: elems_in(:,:)
    character :: readfile*(100)

    centre_nodes = (/13,15,11,17,5,23/)
    if(index(filename, ".lselem")>0)then ! full filename is given
       readfile = trim(filename)
    else! append correct extension
       readfile = trim(filename)//'.lselem'
    endif
    
    open(10, file = readfile, status = 'old')
    read(unit = 10, fmt =  *) tissue_num_elems
    
    if(allocated(tissue_elems))then
       deallocate(tissue_elems)
       deallocate(tissue_elem_nodes)
       deallocate(num_adjacent)
    endif
    allocate(tissue_elems(tissue_num_elems))
    allocate(tissue_elem_nodes(27,tissue_num_elems))
    allocate(num_adjacent(tissue_num_nodes))
    allocate(elems_in(0:10,tissue_num_nodes))
    allocate(npne(nnm,nbm,tissue_num_elems))
    allocate(ce(nmm,tissue_num_elems))
    allocate(material_at_gp(4,ngm,tissue_num_elems))
    allocate(yg(2,ngm,tissue_num_elems))

    num_adjacent = 0
    elems_in = 0
       
    do noelem = 1,tissue_num_elems
       !     read the element and its nodes
       read(unit = 10, fmt =  *) ne,nodelist(1:27)
       npne(:,1,ne) = nodelist(:)
       tissue_elems(noelem) = ne
       do nn = 1,27
          np = local_tissue_node(nodelist(nn))
          tissue_elem_nodes(nn,noelem) = np
          num_adjacent(np) = num_adjacent(np) + 1
       enddo
    enddo                     ! noelem

    close(10)

    ! adjust the number of adjacent elements for nodes that are repeated
    do noelem = 1,tissue_num_elems
       do i = 1,6 ! for each face
          nn = centre_nodes(i)
          npc = tissue_elem_nodes(nn,noelem)
          if(count(tissue_elem_nodes(:,noelem).eq.npc).gt.1)then
             ! the centre node is repeated, must be a collapsed element
             num_repeats = 0
             do noelem2 = 1,tissue_num_elems
                if(count(tissue_elem_nodes(:,noelem2).eq.npc).gt.1) &
                     num_repeats = num_repeats + 1
             enddo
             if(num_repeats.ge.2) num_adjacent(npc) = 1
          endif
       enddo
    enddo

    deallocate(elems_in)
   
  end subroutine define_elem_tissue
  
!!!##############################################################################################

  subroutine define_elem_cavity(filename)
    character :: filename*(*)
    !     Local Variables
    integer :: i,ne,noelem,temp_npne(9)
    character :: readfile*(100)

    if(index(filename, ".lselem")>0)then ! full filename is given
       readfile = trim(filename)
    else! append correct extension
       readfile = trim(filename)//'.lselem'
    endif
    
    open(10, file = readfile, status = 'old')
    read(unit = 10, fmt =  *) cavity_num_elems
    
    if(allocated(cavity_elems)) deallocate(cavity_elems)
    allocate(cavity_elems(cavity_num_elems))
    allocate(npne_cavity(nnm,cavity_num_elems))
       
    do noelem = 1,cavity_num_elems
       read(unit = 10, fmt =  *) ne,temp_npne(:) 
       do i = 1,9
          npne_cavity(i,noelem) = local_cavity_node(temp_npne(i))
       enddo
       cavity_elems(noelem) = ne
    enddo ! noelem
    close(10)
    
  end subroutine define_elem_cavity
  
!!!##############################################################################################

  subroutine isort(n,idata)

!###    ISORT sorts N integer IDATA values into a non-decreasing
!###    sequence using ISHELLSORT if N<50 or IHEAPSORT if N>50 as
!###    recommended by numerical recipes.

    !     Parameter List
    integer :: IDATA(*),N

    if(N.LE.50) then
       call ISHELLSORT(N,IDATA)
    else
       call IHEAPSORT(N,IDATA)
    endif

  end subroutine isort

!!!##############################################################################################

  subroutine iheapsort(n,ra)

    !###    Numerical recipes integer heap sort algorithm. Use when N > 50

    !     Parameter List
    integer :: n,ra(n)
    !     Local Variables
    integer :: i,ir,j,l,rra

    if(N.LT.2) RETURN
    
    L = N/2+1
    IR = N
10  CONTINUE
    if(L.GT.1) then
       L = L-1
       RRA = RA(L)
    else
       RRA = RA(IR)
       RA(IR) = RA(1)
       IR = IR-1
       if(IR.EQ.1) then
          RA(1) = RRA
          RETURN
       endif
    endif
    I = L
    J = L+L
20  if(J.LE.IR) then
       if(J.LT.IR) then
          if(RA(J).LT.RA(J+1)) J = J+1
       endif
       if(RRA.LT.RA(J)) then
          RA(I) = RA(J)
          I = J
          J = J+J
       else
          J = IR+1
       endif
       GOTO 20
    endif
    RA(I) = RRA
    GOTO 10
  end subroutine iheapsort

!!!##############################################################################################

  subroutine ishellsort(n,a)

    !###    Numerical recipes integer shell sort algorithm. Use when N < 50

    !     Parameter List
    integer :: n,a(n)
    !     Local Variables
    integer :: i,j,inc,v

    INC = 1
1   INC = 3*INC+1
    if(INC.LE.N) GOTO 1
2   CONTINUE
    INC = INC/3
    do I = INC+1,N
       V = A(I)
       J = I
3      if(A(J-INC).GT.V)then
          A(J) = A(J-INC)
          J = J-INC
          if(J.LE.INC) GOTO 4
          GOTO 3
       endif
4      A(J) = V
    enddo
    if(INC.GT.1) GOTO 2

  end subroutine ishellsort

!!!##############################################################################################

  subroutine ilistrmdup(n,idata)

    !###    ILISTRMDUP sorts N integer IDATA values into a non-decreasing
    !###    sequence using IHEAPSORT (N>50) or ISHELLSORT (N>50) and then 
    !###    removes all duplicates efrom the list. On exit N contains the 
    !###    number of unique elements in the list.

    !     Parameter List
    integer :: idata(*),n
    !     Local Variables
    integer :: index,nolist

    call ISORT(N,IDATA)

    index = 0
    do nolist = 2,N
       if(IDATA(nolist).EQ.IDATA(nolist-1)) then
          index = index+1
       else
          IDATA(nolist-index) = IDATA(nolist)
       endif
    enddo !nolist
    
    N = N-index

  end subroutine ilistrmdup

!!!##############################################################################################

  subroutine exnode_3d(num_list,nplist,FILE,NODE_NAME)

    !     Parameter List
    integer :: num_list,nplist(:)
    character :: FILE*100,NODE_NAME*50
    !     Local Variables
    integer :: ifile = 10,nj,NOLIST,np
    character :: readfile*100

    if(index(file, ".exnode")>0)then ! full filename is given
       readfile = trim(file)
    else! append correct extension
       readfile = trim(file)//'.exnode'
    endif
    
    open(ifile, file = readfile, status = 'replace')
    
    !**   write the group name
    WRITE(IFILE,'( '' Group name: '',A)') trim(NODE_NAME)
    write(ifile,'('' #Fields=1'')')
    write(ifile,'('' 1) coordinates, coordinate,'',' &
         //''' rectangular cartesian, #Components=3'')')
    write(ifile,'(''   x. Value index=1, #Derivatives= 0'')')
    write(ifile,'(''   y. Value index=2, #Derivatives= 0'')')
    write(ifile,'(''   z. Value index=3, #Derivatives= 0'')')

    do NOLIST = 1,num_list !NPLIST(0)
       NP = NPLIST(NOLIST)
       !**   write the node
       WRITE(IFILE,'(1X,''Node: '',I12)') NP
       do nj = 1,3
          WRITE(IFILE,'(2X,1(1X,E24.16))') xp(nj,np) !tissue_xyz(nj,np) !xp(nj,np)
       enddo
    enddo                     !nolist (np)
    
    close(ifile)

  end subroutine exnode_3d

!!!##############################################################################################

  subroutine setup_mechanics(density,a,b,c,scale)

    real(dp),intent(in) :: a,b,c,density,scale
    integer :: nc,ne,nj_pressure,np

    call problem_setup
    
    nc = 1
    nyt = 0
    nym = tissue_num_nodes * 3 * 2 !nodes*nhm*nc
    
    allocate(nony(0:1,nym,2))
    allocate(npny(0:6,nym,0:2))  ! npny(3/4/5/6) = nh,np,nc,nr; use (3/4,ny,nrc1), (3,ny1,0)
    allocate(nynr(0:nym,0:2,2))
    allocate(yp(nym,5))
    allocate(fix_mech(nym,3))
    allocate(nynp(nhm,tissue_num_nodes,0:2,2))
    allocate(cg(nmm,ngm))
    allocate(zp(nhm,tissue_num_nodes,2))

    call setup_mechs_ny_dep
      
    nz_gk_m = nyt**2
    nz_gkk_m = nz_gk_m ! overly conservative but doesn't make much difference
    nom = 2*nyt ! note that if nom isn't big enough the the 'solution' just skips through
    
    allocate(gk(nz_gk_m))
    allocate(gkk(nz_gkk_m))
    allocate(gr(nyt))
    allocate(grr(nyt))
    allocate(nyno(0:1,nom,2))
    allocate(xo(nom))

    !---------------Calculate constitutive law params --------------------
    
    il_density = 4
    
    do ne = 1,tissue_num_elems
       ! define element material properties
       ce(1,ne) = a
       ce(2,ne) = b
       ce(3,ne) = c
       ce(il_density,ne) = density
       material_at_gp(1,:,ne) = a ! effective stiffness parameter
       material_at_gp(2,:,ne) = b
       material_at_gp(3,:,ne) = c
       material_at_gp(4,:,ne) = density
    enddo                     !noelem (ne)
    nj_pressure = 5           ! the first field after x,y,z,fibre
    xp(nj_pressure,:) = 0.0_dp
    
    do np = 1,tissue_num_nodes
       yp(nynp(1:3,np,0,1),1) = xp(1:3,np)   ! set initial solution same as read/fitted geometry
       xp(1:3,np) = xp(1:3,np) * scale       ! scale read geometry to reference size/shape
       zp(1:3,np,1) = yp(nynp(1:3,np,0,1),1)
    enddo 
    
    gravity(1:3) = 0.0_dp

  end subroutine setup_mechanics

!!!##############################################################################################

  subroutine setup_mechs_ny_dep

    !     Local Variables
    integer :: nc,nh,np,nrc,ny,ny_start(0:2,3),ny_max

    ny_start = 0
    nynp = 0
    npny = 0
    nynr = 0
    ny_max = 0

    !***  Set up mapping arrays for current region.
    do nrc = 0,2
       if(nrc.EQ.0) then
          ny = 0
          do nc = 1,2  !LHS and RHS and GD variables
             ny = ny+ny_start(nrc,nc)
          enddo
       endif
       do nc = 1,2  !LHS and RHS and GD variables
          if(nrc.NE.0) ny = ny_start(nrc,nc)
          nynr(0,nrc,nc)=0
          do np = 1,tissue_num_nodes
             do nh = 1,3
                ny = ny+1
                nynr(0,nrc,nc) = nynr(0,nrc,nc)+1
                if(ny.gt.ny_max) ny_max = ny
                if(nrc.NE.0.and.ny.GT.nyt) nyt = ny
                if(nynr(0,nrc,nc).le.nym) &
                     nynr(nynr(0,nrc,nc),nrc,nc) = ny
                if(ny.LE.NYM) then
                   nynp(nh,np,nrc,nc) = ny
                   if(nrc.NE.1.OR.nc.EQ.1) then
                      npny(0:6,ny,nrc) = 1
                      npny(3,ny,nrc) = nh
                      npny(4,ny,nrc) = np
                      npny(5,ny,nrc) = nc
                   endif
                endif         !ny
             enddo            !nh
          enddo               !np
       enddo                  !nc
    enddo                     !nrc
    
  end subroutine setup_mechs_ny_dep

!!!##############################################################################################

  subroutine fix_external(fix_group)
    
    !     Parameter List
    integer,intent(in) :: fix_group(:)
    !     Local Variables
    integer :: n,nc,n_fix,nh,np,ny

    n_fix = count(fix_group.ne.0)
    fix_mech(:,1:2) = .false.
    fix_mech(:,3) = .true.

    !     Apply bdry conditions to all nodes in group (fix displacement and forces)
    do n = 1,n_fix
       np = fix_group(n)
       do nh = 1,3      ! dep vars same as coords
          do nc = 1,2         ! displ / force
             ny = nynp(nh,np,0,nc)
             fix_mech(ny,1:2) = .true. ! displacements and forces are fixed at this node group
          enddo               ! nc
       enddo                  ! nh
    enddo                     ! n
    
    call ypzp(1)
    call ipsolv

  end subroutine fix_external

!!!##############################################################################################

  subroutine fix_disp_xyz(n_fix_xy,fix_xy,n_fix_z,fix_z)
    
    !     Parameter List
    integer,intent(in) :: n_fix_xy,n_fix_z,fix_xy(:),fix_z(:)
    !     Local Variables
    integer :: n,nc,nh,np,ny
    
    ! start by making sure nothing is fixed
    fix_mech(:,1:2) = .false.
    fix_mech(:,3) = .true.

    ! Fix nodes in x and y for group fix_xy
    do n = 1,n_fix_xy
       np = fix_xy(n)
       do nh = 1,2      ! dep vars same as coords
          do nc = 1,2
             ny = nynp(nh,np,0,nc)
             fix_mech(ny,1:2) = .true. ! displacements and forces are fixed at this node group
          enddo
       enddo                  ! nh
    enddo                     ! n
    
    ! Fix nodes in z for group fix_z
    nc = 1 ! displacement only
    nh = 3      ! dep vars same as coords
    do n = 1,n_fix_z
       np = fix_z(n)
       ny = nynp(nh,np,0,nc)
       fix_mech(ny,1:2) = .true. ! displacements are fixed at this node group
    enddo                     ! n
    
    call ypzp(1)
    call ipsolv

  end subroutine fix_disp_xyz

!!!##############################################################################################

  subroutine fix_xyz(fix_x,fix_y,fix_z,fix_force)

    !     Parameter List
    integer ::  n_fix_x,n_fix_y,n_fix_z,n_force,fix_x(:),fix_y(:),fix_z(:),fix_force(:)
    !     Local Variables
    integer :: n,nc,nh,np,ny

    nc = 1
    n_fix_x = count(fix_x.ne.0)
    n_fix_y = count(fix_y.ne.0)
    n_fix_z = count(fix_z.ne.0)
    n_force = count(fix_force.ne.0)

    fix_mech(:,1:2) = .false.
    fix_mech(:,3) = .true.

    !     Apply bdry conditions to all nodes in group
    do n = 1,n_fix_x
       np = fix_x(n)
       ny = nynp(1,np,0,nc)
       FIX_MECH(ny,1:2) = .true. ! displacements and forces fixed
    enddo                     ! n
    do n = 1,n_fix_y
       np = fix_y(n)
       ny = nynp(2,np,0,nc)
       FIX_MECH(ny,1:2) = .true.
    enddo                     ! n
    do n = 1,n_fix_z
       np = fix_z(n)
       ny = nynp(3,np,0,nc)
       FIX_MECH(ny,1:2) = .true.
    enddo                     ! n

    nc = 2 ! for forces
    do n = 1,n_force
       np = fix_force(n)
       do nh = 1,3      ! dep vars same as coords
          ny = nynp(nh,np,0,nc)
          fix_mech(ny,1:2) = .true.
       enddo
    enddo

    call ypzp(1)
    call ipsolv
    
  end subroutine fix_xyz
  
!!!##############################################################################################

  subroutine get_edge_tissue_nodes(edge_nodes,option)
    
    integer,allocatable :: edge_nodes(:)
    character(len = *) :: option
    
    integer,allocatable :: node_list(:)
    integer :: centre_nodes(6),face_nodes(9,6),i,j,k,nn,noelem,np1,np2,npc,num_edge, &
         nedge(4),index_fn(3,4)

    centre_nodes = (/13,15,11,17,5,23/)
    face_nodes = reshape((/1,4,7,10,13,16,19,22,25,3,6,9,12,15,18,21,24,27, &
         1,2,3,10,11,12,19,20,21,7,8,9,16,17,18,25,26,27, &
         1,2,3,4,5,6,7,8,9,19,20,21,22,23,24,25,26,27/), shape(face_nodes))
    index_fn = reshape((/1,2,3,1,4,7,3,6,9,7,8,9/), shape(index_fn))

    ! 13 --> node Xi1 = 0 face --> nodes 1,4,7,10,13,16,19,22,25
    ! 15 --> node Xi1 = 1 face --> nodes 3,6,9,12,15,18,21,24,27
    ! 11 --> node Xi2 = 0 face --> nodes 1,2,3,10,11,12,19,20,21
    ! 17 --> node Xi2 = 1 face --> nodes 7,8,9,16,17,18,25,26,27
    !  5 --> node Xi3 = 0 face --> nodes 1,2,3,4,5,6,7,8,9
    ! 23 --> node Xi3 = 1 face --> nodes 19,20,21,22,23,24,25,26,27

    allocate(node_list(tissue_num_nodes))
    node_list = 0
    if(allocated(edge_nodes)) deallocate(edge_nodes)
    num_edge = 0

    select case (trim(option))
    case ('xi1_0')
       i = 1
       nedge = (/4,10,16,22/) 
    case ('xi1_1')
       i = 2
       nedge = (/6,12,18,24/)
    case ('xi2_0')
       i = 3
       nedge = (/2,10,12,20/)
    case ('xi2_1')
       i = 4
       nedge = (/8,16,18,26/)
    case ('xi3_0')
       i = 5
       nedge = (/2,4,6,8/)
    case ('xi3_1')
       i = 6
       nedge = (/20,22,24,26/)
    end select
    
    do noelem = 1,tissue_num_elems
       nn = centre_nodes(i)
       npc = tissue_elem_nodes(nn,noelem)
       if(num_adjacent(npc).eq.1)then ! no adjacent elements
          do j = 1,4
             np1 = tissue_elem_nodes(nedge(j),noelem) ! nodes at centre of lines
             if(num_adjacent(np1).eq.1)then ! at an edge
                do k = 1,3
                   np2 = tissue_elem_nodes(face_nodes(index_fn(k,j),i),noelem)
                   if(.not.any(np2 == node_list))then
                      num_edge = num_edge + 1
                      node_list(num_edge) = np2
                   endif
                enddo !k
             endif
          enddo !j
       endif
    enddo !noelem
    
    allocate(edge_nodes(num_edge))
    edge_nodes(1:num_edge) = node_list(1:num_edge)
    deallocate(node_list)

  end subroutine get_edge_tissue_nodes
  
!!!##############################################################################################

  subroutine get_external_nodes_subset(subset,external_nodes,option)

    integer :: subset(:)
    integer,allocatable :: external_nodes(:)
    character(len=*) :: option
    
    integer,allocatable :: node_list(:)
    integer :: centre_nodes(6),face_nodes(9,6),i,i1,i2,j,ne,nn,noelem,np,npc,num_external

    centre_nodes = (/13,15,11,17,5,23/)
    face_nodes = reshape((/1,4,7,10,13,16,19,22,25,3,6,9,12,15,18,21,24,27, &
         1,2,3,10,11,12,19,20,21,7,8,9,16,17,18,25,26,27, &
         1,2,3,4,5,6,7,8,9,19,20,21,22,23,24,25,26,27/), shape(face_nodes))

    ! 13 --> node Xi1 = 0 face --> nodes 1,4,7,10,13,16,19,22,25
    ! 15 --> node Xi1 = 1 face --> nodes 3,6,9,12,15,18,21,24,27
    ! 11 --> node Xi2 = 0 face --> nodes 1,2,3,10,11,12,19,20,21
    ! 17 --> node Xi2 = 1 face --> nodes 7,8,9,16,17,18,25,26,27
    !  5 --> node Xi3 = 0 face --> nodes 1,2,3,4,5,6,7,8,9
    ! 23 --> node Xi3 = 1 face --> nodes 19,20,21,22,23,24,25,26,27

    allocate(node_list(tissue_num_nodes))
    node_list = 0
    if(allocated(external_nodes)) deallocate(external_nodes)
    num_external = 0

    select case (trim(option))
    case ('all')
       i1 = 1
       i2 = 6
    case ('xi1_0')
       i1 = 1
       i2 = 1
    case ('xi1_1')
       i1 = 2
       i2 = 2
    case ('xi2_0')
       i1 = 3
       i2 = 3
    case ('xi2_1')
       i1 = 4
       i2 = 4
    case ('xi3_0')
       i1 = 5
       i2 = 5
    case ('xi3_1')
       i1 = 6
       i2 = 6
    end select
    
    do noelem = 1,count(subset.ne.0)
       ne = subset(noelem)
       do i = i1,i2
          nn = centre_nodes(i)
          npc = tissue_elem_nodes(nn,ne)
          if(num_adjacent(npc).eq.1)then 
             do j = 1,9
                np = tissue_elem_nodes(face_nodes(j,i),ne)
                if(.not.any(np == node_list))then
                   num_external = num_external + 1
                   node_list(num_external) = np
                endif
             enddo !n
          endif
       enddo ! i
    enddo !noelem
    
    allocate(external_nodes(num_external))
    external_nodes(1:num_external) = node_list(1:num_external)
    deallocate(node_list)

  end subroutine get_external_nodes_subset
  
!!!##############################################################################################

  subroutine get_external_tissue_nodes(external_nodes,option)
    
    integer,allocatable :: external_nodes(:)
    character(len = *) :: option
    
    integer,allocatable :: node_list(:)
    integer :: centre_nodes(6),face_nodes(9,6),i,i1,i2,j,nn,noelem,np,npc,num_external

    centre_nodes = (/13,15,11,17,5,23/)
    face_nodes = reshape((/1,4,7,10,13,16,19,22,25,3,6,9,12,15,18,21,24,27, &
         1,2,3,10,11,12,19,20,21,7,8,9,16,17,18,25,26,27, &
         1,2,3,4,5,6,7,8,9,19,20,21,22,23,24,25,26,27/), shape(face_nodes))

    ! 13 --> node Xi1 = 0 face --> nodes 1,4,7,10,13,16,19,22,25
    ! 15 --> node Xi1 = 1 face --> nodes 3,6,9,12,15,18,21,24,27
    ! 11 --> node Xi2 = 0 face --> nodes 1,2,3,10,11,12,19,20,21
    ! 17 --> node Xi2 = 1 face --> nodes 7,8,9,16,17,18,25,26,27
    !  5 --> node Xi3 = 0 face --> nodes 1,2,3,4,5,6,7,8,9
    ! 23 --> node Xi3 = 1 face --> nodes 19,20,21,22,23,24,25,26,27

    allocate(node_list(tissue_num_nodes))
    node_list = 0
    if(allocated(external_nodes)) deallocate(external_nodes)
    num_external = 0

    select case (trim(option))
    case ('all')
       i1 = 1
       i2 = 6
    case ('xi1_0')
       i1 = 1
       i2 = 1
    case ('xi1_1')
       i1 = 2
       i2 = 2
    case ('xi2_0')
       i1 = 3
       i2 = 3
    case ('xi2_1')
       i1 = 4
       i2 = 4
    case ('xi3_0')
       i1 = 5
       i2 = 5
    case ('xi3_1')
       i1 = 6
       i2 = 6
    end select
    
    do noelem = 1,tissue_num_elems
       do i = i1,i2
          nn = centre_nodes(i)
          npc = tissue_elem_nodes(nn,noelem)
          if(num_adjacent(npc).eq.1)then 
             do j = 1,9
                np = tissue_elem_nodes(face_nodes(j,i),noelem)
                if(.not.any(np == node_list))then
                   num_external = num_external + 1
                   node_list(num_external) = np
                endif
             enddo !n
          endif
       enddo ! i
    enddo !noelem
    
    allocate(external_nodes(num_external))
    external_nodes(1:num_external) = node_list(1:num_external)
    deallocate(node_list)

  end subroutine get_external_tissue_nodes
  
!!!##############################################################################################

  subroutine get_external_tissue_elems(external_elems,option)
    
    integer,allocatable :: external_elems(:)
    character(len = *) :: option
    
    integer,allocatable :: elem_list(:)
    integer :: centre_nodes(6),i,i1,i2,nn,noelem,npc,num_external

    centre_nodes = (/13,15,11,17,5,23/)

    ! 13 --> node Xi1 = 0 face --> nodes 1,4,7,10,13,16,19,22,25
    ! 15 --> node Xi1 = 1 face --> nodes 3,6,9,12,15,18,21,24,27
    ! 11 --> node Xi2 = 0 face --> nodes 1,2,3,10,11,12,19,20,21
    ! 17 --> node Xi2 = 1 face --> nodes 7,8,9,16,17,18,25,26,27
    !  5 --> node Xi3 = 0 face --> nodes 1,2,3,4,5,6,7,8,9
    ! 23 --> node Xi3 = 1 face --> nodes 19,20,21,22,23,24,25,26,27

    allocate(elem_list(tissue_num_elems))
    elem_list = 0
    if(allocated(external_elems)) deallocate(external_elems)
    num_external = 0

    select case (trim(option))
    case ('all')
       i1 = 1
       i2 = 6
    case ('xi1_0')
       i1 = 1
       i2 = 1
    case ('xi1_1')
       i1 = 2
       i2 = 2
    case ('xi2_0')
       i1 = 3
       i2 = 3
    case ('xi2_1')
       i1 = 4
       i2 = 4
    case ('xi3_0')
       i1 = 5
       i2 = 5
    case ('xi3_1')
       i1 = 6
       i2 = 6
    end select
    
    do noelem = 1,tissue_num_elems
       do i = i1,i2
          nn = centre_nodes(i)
          npc = tissue_elem_nodes(nn,noelem)
          if(num_adjacent(npc).eq.1)then 
             if(.not.any(noelem == elem_list))then
                num_external = num_external + 1
                elem_list(num_external) = noelem
             endif
          endif
       enddo ! i
    enddo !noelem
    
    allocate(external_elems(num_external))
    external_elems(1:num_external) = elem_list(1:num_external)
    deallocate(elem_list)

  end subroutine get_external_tissue_elems
  
!!!##############################################################################################

  function get_cavity_elems(map_array,option)

    integer :: map_array(:,:)
    integer,allocatable :: get_cavity_elems(:)
    character(len = *) :: option
    
    integer,allocatable :: elem_list(:)
    integer :: i,i1,i2,noelem,num_external

    allocate(elem_list(tissue_num_elems))
    elem_list = 0
    if(allocated(get_cavity_elems)) deallocate(get_cavity_elems)
    num_external = 0

    select case (trim(option))
    case ('all')
       i1 = 1
       i2 = 6
    case ('xi1_0')
       i1 = 1
       i2 = 1
    case ('xi1_1')
       i1 = 2
       i2 = 2
    case ('xi2_0')
       i1 = 3
       i2 = 3
    case ('xi2_1')
       i1 = 4
       i2 = 4
    case ('xi3_0')
       i1 = 5
       i2 = 5
    case ('xi3_1')
       i1 = 6
       i2 = 6
    end select

    do noelem = 1,tissue_num_elems
       do i = i1,i2
          if(map_array(i,noelem).gt.0)then
             if(.not.any(map_array(i,noelem) == elem_list))then
                num_external = num_external + 1
                elem_list(num_external) = map_array(i,noelem)
             endif
          endif
       enddo ! i
    enddo !noelem
    
    allocate(get_cavity_elems(num_external))
    get_cavity_elems(1:num_external) = elem_list(1:num_external)
    deallocate(elem_list)

  end function get_cavity_elems

!!!##############################################################################################

  function map_tissue_to_cavity_elems()

    integer,allocatable :: map_tissue_cavity(:,:),map_tissue_to_cavity_elems(:,:)
    integer :: centre_nodes(6),i,j,ne_cavity,ne_map,ne_tissue,nn,npc
    real(dp) :: distance,min_distance,x1(3),x2(3)

    centre_nodes = (/13,15,11,17,5,23/)
    if(allocated(map_tissue_cavity)) deallocate(map_tissue_cavity)
    allocate(map_tissue_cavity(6,tissue_num_elems))
    allocate(map_tissue_to_cavity_elems(6,tissue_num_elems))
    map_tissue_cavity = 0
    map_tissue_to_cavity_elems = 0
    
    do ne_tissue = 1,tissue_num_elems
       do i = 1,6 ! for each face
          nn = centre_nodes(i) ! index of node at centre of face
          npc = npne(nn,1,ne_tissue) ! node number at centre of face
          x1(1:3) = xp(1:3,npc)
          if(num_adjacent(npc).eq.1)then  ! this is an external face
             ! find the closest cavity face
             min_distance = 1.0e+6_dp
             do j = 1,cavity_num_elems
                ne_cavity = j
                x2(1:3) = xp_cavity(1:3,npne_cavity(5,ne_cavity))
                distance = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)
                if(distance.lt.min_distance)then
                   min_distance = distance
                   ne_map = ne_cavity ! local cavity element number
                endif
             enddo ! j
             map_tissue_cavity(i,ne_tissue) = ne_map
          endif
       enddo ! i
    enddo !noelem
    map_tissue_to_cavity_elems = map_tissue_cavity

  end function map_tissue_to_cavity_elems
  
!!!##############################################################################################

  subroutine zpyp(iy)
    
    !###    ZPYP transfers global node parameters ZP(nh,np,nc) to global vector
    !###    YP(ny,iy).

    !     Parameter List
    integer :: iy
    !     Local Variables
    integer :: nc,nh,np,nrc,ny

    nrc = 0                     !only want the global variables
    do nc = 1,2
       do np = 1,tissue_num_nodes
          do nh = 1,3
             ny = nynp(nh,np,nrc,nc)
             yp(ny,iy) = zp(nh,np,nc)
          enddo !nh
       enddo !nonode (np)
    enddo !nc

  end subroutine ZPYP

!!!##############################################################################################

  subroutine YPZP(iy)

    !###    YPZP transfers global vector YP(ny,iy) to global node parameters
    !###    ZP(nh,np,nc).

    !     Parameter List
    integer :: iy
    !     Local Variables
    integer :: nc,nh,np,ny

    do nc = 1,2
       do np = 1,tissue_num_nodes
          do nh = 1,3
             ny = nynp(nh,np,0,nc)
             ZP(nh,np,nc) = YP(ny,iy)
          enddo                 !nh
       enddo                    !nonode (np)
    enddo !nc

  end subroutine YPZP

!!!##############################################################################################

  subroutine ipsolv

    !###    IPSOLV defines solution parameters

    !     Local Variables
    integer :: no_nynr,no_tot(2),nrc,ny1v,nyy(2)

    !***  Initialise mapping arrays above current region
    nyno = 0
    nony = 0
    not = 0

    do no_nynr = 1,nynr(0,0,1) !Loop over the global vars for nr
       ny1v = nynr(no_nynr,0,1) !global var#
       nyy(1) = getnyr(1,1,0,ny1v)      !row#
       nyy(2) = nynr(no_nynr,0,1)       !col#
       
       if(.not.FIX_MECH(ny1v,1))then  !add ny to list of solve variables
          do nrc = 1,2    !rows and columns
             not(nrc) = not(nrc)+1
             nony(0,nyy(nrc),nrc) = 1 ! indicates free dep variable
             nony(1,nyy(nrc),nrc) = not(nrc)
             if(no_tot(nrc).LE.NOM) then
                nyno(0,not(nrc),nrc) = 1
                nyno(1,not(nrc),nrc) = nyy(nrc)
             endif
          enddo               !nrc
          
       endif                  !.not.SPECIAL
       
    enddo                     !no_nynr   --------- end of no_nynr loop for FEM case ----------

    firsts = .true.
    
  end subroutine ipsolv

!!!##############################################################################################

  subroutine evpress

    !C###    EVPRESS evaluates effective pressures at nodes for a
    !C###    finite deformation elasticity problem. This has been implemented 
    !C###    for compressible mechanics, where an effective hydrostatic
    !C###    pressure is developed as the tissue changes in volume.

    !     Local Variables
    integer :: ne,ng,nm
    
    real(dp) :: AZ,AZL(3,3),AZU(3,3),RG2D,RGZ,RGZ2D,TC(3,3),TG(3,3),TN(3,3),xi(3)
    real(dp) :: xe(nsm,20),xg(20,num),ze(nsm,nhm),zg(nhm,num)

    yg = 0.0_dp
    do ne = 1,tissue_num_elems
       call xpxe(nb_lung,ne,xe)
       call zpze(nb_lung,ne,ze)
       forall(nm = 1:4) cg(nm,:) = ce(nm,ne) 
       do ng = 1,NGT(nb_lung)
          xi(1:3) = xig(1:3,ng)
          call ZETX50(nb_lung,ng,RG2D,RGZ,RGZ2D,TC,TG,TN,xe,xg,ze,zg) ! only call
          call ZGMG(nb_lung,AZ,AZL,AZU,zg)
          YG(1,ng,ne) = sqrt(DET(AZL)) !ratio deformed to undeformed
          yg(2,ng,ne) = (TC(1,1)+TC(2,2)+TC(3,3))/3.0_dp
          if(yg(1,ng,ne).lt.1.25_dp)then
             ! want to use positive value of pressure for xg in zgtg53
             !write(*,'('' >>Warning: Ratio at ng='',I5,'' ne='',I5,'' less than 1.25'')') ng,ne
             !write(*,'('' Ratio ='',D12.4)') yg(1,ng,ne)
             !write(*,'('' pressure'',f10.3)') yg(2,ng,ne)/98.0665_dp
          endif
       enddo                  !ng
    enddo                     !noelem

  end subroutine evpress

!!!##############################################################################################

  subroutine xpxe(nb,ne,xe)
    
    !     Parameter List
    integer :: nb,ne
    real(dp) :: xe(:,:)
    !     Local Variables
    integer :: nj,nn,np,ns
    
    do nj = 1,njm
       ns = 0
       do nn = 1,nnt(nb)
          np = npne(nn,nb,ne)
          ns = ns+1
          XE(ns,nj) = XP(nj,np)
       enddo
    enddo                     !njj

  end subroutine xpxe

!!!##############################################################################################

  subroutine xpxe_cavity(nb,ne,xe)
    
    !     Parameter List
    integer :: nb,ne
    real(dp) :: xe(:,:)
    !     Local Variables
    integer :: nj,nn,np,ns
    
    do nj = 1,njm
       ns = 0
       do nn = 1,nnt(nb)
          np = npne_cavity(nn,ne)
          ns = ns+1
          xe(ns,nj) = xp_cavity(nj,np)
       enddo
    enddo !nj

  end subroutine xpxe_cavity

!!!##############################################################################################

  subroutine zpze(nb,ne,ze)

    !###    zpze transfers global node parameters ZP(nh,np,nc) to element node
    !###    parameters ZE(ns,nhx) for nc.

    !     Parameter List
    integer :: nb,ne
    real(dp) :: ze(:,:)
    !     Local Variables
    integer :: nh,nn,np,ns

    do nh = 1,3
       ns = 0
       do nn = 1,nnt(nb)
          np = npne(nn,nb,ne)
          ns = ns+1
          ZE(ns,nh) = ZP(nh,np,1)
       enddo                  !nn
    enddo                     !nhx

  end subroutine zpze

!!!##############################################################################################

  subroutine CPCG(ne)

    !###    CPCG transfers element parameters CE(il,ne) in element ne to return Gauss pt values
    !###    CG(il,ng). 

    !     Parameter List
    integer :: ne
    !     Local Variables
    integer :: il

    cg = 0.0_dp
!    forall(il=1:4) CG(il,:) = CE(il,ne)
    forall(il=1:4) CG(il,:) = material_at_gp(il,:,ne)

  end subroutine CPCG

!!!##############################################################################################

  subroutine ZETX50(nb,ng,RGX2D,RGZ,RGZ2D,TC,TG,TN,xe,xg,ze,zg)

    !###    ZETX50 calculates 2nd Piola-Kirchhoff, Cauchy and  Nominal
    !###    stresses with respect to 'Reference' or 'Fibre' coords
    !###    (as specified by COORDS) at position XI in current element
    !###    if ng=0 else at Gauss point ng.  STRESSTYPE ('Total', 'Passive'
    !###    or 'Active') denotes the components of stress to be computed.

    !**** Since base vectors of theta coords are not orthonormal, cpts
    !**** of Cauchy and Nominal stress are converted to 'physical' values.
    !**** AZ,AZL,AZU  are deformed metric tensors wrt undeformed coords
    !**** RI1,RI2,RI3 are principal invariants of AZL
    !**** AXU  are contravariant cpts of undeformed metric tensor
    !**** XG   are undeformed theta coords and derivs wrt Xi
    !**** ZG   are deformed theta coords and derivs wrt undeformed coords
    !**** TG   are tensor cpts of 2nd Piola-Kirchhoff stresses
    !**** TN   are physical cpts of Nominal stresses
    !**** TC   are physical cpts of Cauchy stress
    !**** AZLZ are deformed Eulerian metrics (wrt deformed theta coords)

    !     Parameter List
    integer :: nb,ng
    real(dp) ::  RGX2D,RGZ,RGZ2D,TC(3,3),TG(3,3),TN(3,3),xe(:,:),xg(:,:),ze(:,:),zg(:,:)
    !    Local Variables
    integer :: k,mi,mix,mz,ni,nitb,nix,nz
    integer,parameter :: ncw=35 !CW must be dimen.d the same size as CE array
    real(dp) :: AZ,AZL(3,3),AZU(3,3),CW(NCW),DETERM,DET_DZDX, &
         DXDZ(3,3),DXIXJ(3,3),DXIXN(3,3),DXIZN(3,3),DZDX(3,3), &
         DZNXI(3,3),GXL(3,3),GXU(3,3),GZ,GZL(3,3), &
         GZU(3,3),RI3,RWX,SUM

    nitb = nit(nb)
    gzu = 0.0_dp
    CW(1:4) = CG(1:4,ng)   ! Put Gauss pt params into CW array
    
    !-stress referred to Xj--------------------------------------------------
    !     ***     Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
    call XEXG(nb,ng,xe,xg)
    !     ***   Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
    !     ***   derivatives of Xi wrt Xj (reference) coords, DXIXJ (IP=0)
    call XGMG(0,nitb,nb,DXIXJ,GXL,GXU,RWX,xg)
    !     ***   Calculate 2D Jacobian wrt undef coords for face integrals
    RGX2D = RWX* sqrt(GXU(3,3))
    !     ***   Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
    call DXIDXM(nitb,DXIXN,DETERM,xg)
    
    call ZEZG(0,nb,ng,DXIXJ,ze,zg)
    !     ***     Calculate deformed metric tensors wrt Xi (GZL,GZU)
    call ZGMG(nb,GZ,GZL,GZU,zg)
    RGZ =  sqrt(GZ)
    !     ***     Calculate 2D Jacobian wrt def coords for face integrals
    RGZ2D =  sqrt(GZ*GZU(3,3))
    !     Get derivs of Xi wrt deformed Nu coords, DXIZN
    call DXIDZM(nb,ng,DXIZN,DZNXI,'Fibre',xg,ze,zg)
    !     ***     Calculate derivs of deformed Nu wrt undeformed Nu (DZDX)
    do ni = 1,nitb
       do mi = 1,nitb
          SUM = 0.0_dp
          do k = 1,nitb
             SUM = SUM+DZNXI(ni,k)*DXIXN(k,mi)
          enddo               !k
          DZDX(ni,mi) = SUM
       enddo                  !mi
    enddo                     !ni

    call INVERT(nitb,DZDX,DXDZ,DET_DZDX)
    
    call ZEZG(1,nb,ng,DXIXN,ze,zg)
    !     ***   Calculate deformed metric tensors wrt Nu (AZL,AZU)
    call ZGMG(nb,AZ,AZL,AZU,zg)
    !     ***   Get contravariant cpts of 2nd Piola-Kirchhoff stress
    !     ***   tensor (TG) wrt undeformed Nu coordinates
    call ZGTG53(az,azl,azu,ri3,CW,TG,xg)
    
    !-----------------------------------------------------------------------
    
    do mz = 1,nitb
       do nz = 1,nitb
          SUM = 0.0_dp
          do mix = 1,nitb
             do nix = 1,nitb
                SUM = SUM+DZDX(mz,mix)*TG(mix,nix)*DZDX(nz,nix)
             enddo            !nix
          enddo               !mix
          TC(mz,nz) = SUM/sqrt(RI3)
       enddo                  !nz
    enddo                     !mz
    do nix = 1,nitb
       do nz = 1,nitb
          SUM = 0.0_dp
          do mix = 1,nitb
             SUM = SUM+TG(nix,mix)*DZDX(nz,mix)
          enddo               !mix
          TN(nix,nz) = SUM
       enddo                  !nz
    enddo                   !nix

  end subroutine ZETX50
      

!!!##############################################################################################

  subroutine DXIDXM(nitb,DXIXN,RG,xg)

    !C###    DXIDXM evaluates derivatives (DXIXN) of Xi- wrt
    !C###    undeformed Nu(fibre)-coords 
    !C###    RG returns the domain Jacobian.
    !C###    This routine assumes that XG contains
    !C###    the Gauss pt position, interpolated material axis
    !C###    orientations, and derivatives with respect to Xi coordinates.

    !     Parameter List
    integer :: nitb
    real(dp) :: DXIXN(3,3),RG,xg(:,:)
    !     Local Variables
    integer :: mjj,ni,njj
    real(dp) :: DXRCXN(3,3),DXRCXI(3,3),DXIXRC(3,3)

!!! Calculate dXrc/dXI
    dxrcxi(1:3,1) = xg(1:3,2)
    dxrcxi(1:3,2) = xg(1:3,4)
    dxrcxi(1:3,3) = xg(1:3,7)
      
!!! Compute undeformed anatomical fibre vectors wrt rc coordinates
    call MAT_VEC(DXRCXN(1,1),DXRCXN(1,2),DXRCXN(1,3),DXRCXI)
      
!!! Calculate dXI/dXrc
    call INVERT(nitb,DXRCXI,DXIXRC,RG)
      
!!! Calc derivatives of Xi wrt undeformed Nu/Wall
    dxixn = 0.0_dp
    do njj = 1,3
       do mjj = 1,3
          do ni = 1,3
             DXIXN(ni,njj) = DXIXN(ni,njj)+DXIXRC(ni,mjj)*DXRCXN(mjj,njj)
          enddo               !ni
       enddo                  !mjj
    enddo                     !njj
    
  end subroutine DXIDXM

!!!##############################################################################################

  subroutine XEXG(nb,ng,xe,xg)

!!! Evaluate Gauss point array XG from element node array XE at current Gauss point ng. 

    !     Parameter List
    integer :: nb,ng
    real(dp) :: xe(:,:),xg(:,:)
    !     Local Variables
    integer :: nj,nu
    real(dp) :: pg_temp(nsm),xe_temp(nsm)

    xg = 0.0_dp
    
    do nj = 1,njm
       xe_temp(:) = xe(:,nj)
       do nu = 1,NUT(nb)
          pg_temp(:) = pg(:,nu,ng)
          XG(nj,nu) = dot_product(pg_temp,xe_temp)
       enddo                  !nu
    enddo                     !njj2
    
  end subroutine XEXG


!!!##############################################################################################

  subroutine XGMG(IP,JAC,nb,DXIX,GL,GU,RGX,xg)

    !###    XGMG evaluates the covariant (GL) & contravariant (GU) metric
    !###    tensors wrt the Xi-coordinate system  and  the derivs  of
    !###    the Xi-coords wrt the Xj-coords (DXIX) at current Gauss pt.
    !**** If IP=0 DXIX contains derivatives of Xi wrt X(ref)-coords.
    !**** If IP=1 DXIX contains derivatives of Xi wrt Nu(fibre)-coords.
    !**** If IP=-1 DXIX is not touched.
    !**** The Jacobian RG for a length,area or volume integral is returned
    !****   if JAC=1,2 or 3, respec.

    !     Parameter List
    integer :: IP,JAC,nb
    real(dp) :: DXIX(3,3),GL(3,3),GU(3,3),RGX,xg(:,:)
    !     Local Variables
    integer :: mi,ni,nitb,njj,nu
    real(dp) :: D,DXXI(3,3),G


    !     Calculate derivatives of X wrt Xi
    nitb = nit(nb)
    do ni = 1,nitb
       nu = 1+ni*(1+ni)/2
       do njj = 1,3
          DXXI(njj,ni) = XG(njj,nu)
       enddo !njj
    enddo !ni

    !     Initialise metric tensors
    gl = 0.0_dp
    gu = 0.0_dp
    forall (mi = 1:3) gl(mi,mi) = 1.0_dp
    forall (mi = 1:3) gu(mi,mi) = 1.0_dp
    
    do mi = 1,nitb
       do ni = 1,nitb
          GL(mi,ni) = DXXI(1,mi)*DXXI(1,ni)
          do njj = 2,3
             GL(mi,ni) = GL(mi,ni)+DXXI(njj,mi)*DXXI(njj,ni)
          enddo !njj
       enddo !ni
    enddo !mi

    !     Calculate contravariant metric tensor GU(i,j)
    call INVERT(nitb,GL,GU,G)
    if(abs(G).LT.zero_tol) then
       RGX = 0.0_dp
       WRITE(*,'('' >>Warning: zero G in XGMG. G='',D12.5,' &
            //''' zero_tol='',D12.5)') G,zero_tol
       write(*,*) 'stopping here! check XGMG'
       read(*,*)
    endif

    !     Calculate derivs DXIX(i,j) of Xi wrt X (IP=0) or Nu (IP=1)
    if(IP.EQ.0) then          !DXIX is based on reference coords, X
       if(nitb.EQ.3) call INVERT(nitb,DXXI,DXIX,D)
    else if(IP.ge.1) then     !DXIX is based on material fibre coords, Nu
       call DXIDXM(nitb,DXIX,RGX,xg)
    endif
    
    !     Calculate Jacobian RG
    if(JAC.EQ.1) RGX = sqrt(abs(GL(1,1)))
    if(JAC.EQ.2) RGX = sqrt(abs(G*GU(3,3)))
    if(JAC.EQ.3) RGX = sqrt(abs(G))
    
  end subroutine XGMG


!!!##############################################################################################

  subroutine ZEZG(JP,nb,ng,DXIX,ze,zg)

    !     Parameter List
    integer :: JP,nb,ng
    real(dp) :: DXIX(:,:),ze(:,:),zg(:,:)
    !     Local Variables
    integer :: nh,ni,NSTNAT,NU1(0:3)
    real(dp) :: DZDXI(3),pg_temp(nsm),ze_temp(nsm)

    DATA NU1/1,2,4,7/

    zg = 0.0_dp
    do nh = 1,3
       NSTNAT = NST(nb)
       pg_temp(:) = pg(:,1,ng)
       ze_temp(:) = ze(:,nh)
       ZG(nh,1) = dot_product(pg_temp,ze_temp)
       
       if(JP.EQ.0) then      ! return 1st derivs wrt Xi
          do ni = 1,nit(nb)
             pg_temp(:) = pg(:,nu1(ni),ng)
             ZG(nh,NU1(ni)) =  dot_product(pg_temp,ze_temp) 
          enddo
          
       else if(JP.EQ.1) then !return 1st derivs multiplied by dXIX
          do ni = 1,nit(nb)     ! 1st derivatives wrt Xi
             pg_temp(:) = pg(:,nu1(ni),ng)
             dzdxi(ni) =  dot_product(pg_temp,ze_temp)
          enddo
          !         1st derivatives wrt X
          if(nit(nb).EQ.1) then
             ZG(nh,2) = DZDXI(1)*DXIX(1,1)
          else if(nit(nb).EQ.2) then
             ZG(nh,2) = DZDXI(1)*DXIX(1,1) + DZDXI(2)*DXIX(2,1)
             ZG(nh,4) = DZDXI(1)*DXIX(1,2) + DZDXI(2)*DXIX(2,2)
          else if(nit(nb).EQ.3) then
             ZG(nh,2) = DZDXI(1)*DXIX(1,1) + &
                  DZDXI(2)*DXIX(2,1) + DZDXI(3)*DXIX(3,1)
             ZG(nh,4) = DZDXI(1)*DXIX(1,2) + &
                  DZDXI(2)*DXIX(2,2) + DZDXI(3)*DXIX(3,2)
             ZG(nh,7) = DZDXI(1)*DXIX(1,3) + &
                  DZDXI(2)*DXIX(2,3) + DZDXI(3)*DXIX(3,3)
          endif !nit
          
       endif !jp
    enddo

  end subroutine ZEZG

!!!##############################################################################################

  subroutine ZGMG(nb,GZ,GZL,GZU,zg)

    !     Parameter List
    integer :: nb
    real(dp) :: GZ,GZL(3,3),GZU(3,3),zg(:,:)
    !     Local Variables
    integer :: mi,nhx,ni,nitb,NU1(0:3)
    real(dp) :: SUM
    
    DATA NU1/1,2,4,7/
    
    
    nitb = nit(nb)
    do mi = 1,nitb
       do ni = 1,nitb
          SUM = ZG(1,NU1(mi))*ZG(1,NU1(ni))
          do nhx = 2,3
             SUM = SUM+ZG(nhx,NU1(mi))*ZG(nhx,NU1(ni))
          enddo
          GZL(mi,ni) = SUM
       enddo
    enddo
    
    !     Calculate contravariant metric tensor GZU(i,j)
    call INVERT(nitb,GZL,GZU,GZ)
    if(abs(GZ).LT.zero_tol) then
       WRITE(*,'('' >>Warning: zero GZ in ZGMG. GZ='',D12.5,' &
            //''' zero_tol='',D12.5)') GZ,zero_tol
    endif
    
  end subroutine ZGMG

!!!##############################################################################################

  subroutine DXIDZM(nb,ng,DXIZN,DZNXI,COORDS,xg,ze,zg)

    !###    DXIDZM evaluates derivatives (DXIZN) of Xi- wrt
    !###    deformed Nu(fibre)-coords 

    !     Parameter List
    integer :: nb,ng
    real(dp) :: DXIZN(3,3),DZNXI(3,3),xg(:,:),ze(:,:),zg(:,:)
    character :: COORDS*(*)
    !     Local Variables
    integer :: mi,mhx,ni,ni2,nhx,nitb,NU1(0:3)
    real(dp) :: DETERM,DXIX(3,3),DZDNU(3,3),dZrc_dZref,GZ,GZL(3,3),GZU(3,3),SUM

    DATA NU1/1,2,4,7/
    
    nitb = nit(nb)
    dxix = 0.0_dp   
    dxizn = 0.0_dp
    dznxi = 0.0_dp
    forall(ni = 1:3) dxizn(ni,ni) = 1.0_dp
    forall(ni = 1:3) dznxi(ni,ni) = 1.0_dp
    dxizn = 1.0_dp
    dznxi = 1.0_dp
    
    if(COORDS(1:5).EQ.'Fibre') call mat_vec_def(nb,ng,DZDNU(1,1),DZDNU(1,2), &
         DZDNU(1,3),xg,ze,zg)     ! Compute deformed anatomical fibre vectors wrt rc coordinates
    call zezg(0,nb,ng,DXIX,ze,zg) ! Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
    call ZGMG(nb,GZ,GZL,GZU,zg)   ! Calculate deformed metric tensors wrt Xi (GZL,GZU)


    do ni = 1,nitb     !Calc derivs of Xi wrt deformed Nu/Wall coords
       do mi = 1,nitb
          SUM = 0.0_dp
          do ni2 = 1,nitb
             do nhx = 1,3
                do mhx = 1,3
                   dZrc_dZref = 0.0_dp
                   if(mhx.eq.nhx) dZrc_dZref = 1.0_dp
                   SUM = SUM+GZU(ni,ni2)*dZrc_dZref*ZG(nhx,NU1(ni2))* &
                        DZDNU(mhx,mi)
                enddo !mhx
             enddo !nhx
          enddo !ni2
          DXIZN(ni,mi) = SUM
       enddo !mi
    enddo !ni
    
    call INVERT(nitb,DXIZN,DZNXI,DETERM)
    
  end subroutine DXIDZM

!!!##############################################################################################

  subroutine INVERT(n,A,B,AA)

    !###    INVERT returns the inverse of matrix A as B and det(A) as AA.
    !###    Matrix A may be no larger than 3*3 (N=3). Note that in
    !###    both cases A and B are dimensioned to A(3,3) and B(3,3).

    !     Parameter List
    integer :: n
    real(dp) :: A(3,3),AA,B(3,3)
    !     Local Variables
    real(dp) :: MDET(3)

    B = 0.0_dp
    
    if(N.EQ.2) then !2*2 matrix
       AA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
       if(abs(AA).GT.ZERO_TOL) then
          B(1,1) =  A(2,2)/AA
          B(1,2) = -A(1,2)/AA
          B(2,1) = -A(2,1)/AA
          B(2,2) =  A(1,1)/AA
       else
          WRITE(*,'('' >>Warning: Zero determinant in 2*2 INVERT'')')
       endif
       
    else if(N.EQ.3) then !3*3 matrix
       MDET(1) = A(2,2)*A(3,3)-A(2,3)*A(3,2)
       MDET(2) = A(3,2)*A(1,3)-A(3,3)*A(1,2)
       MDET(3) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
       AA = MDET(1)*A(1,1)+MDET(2)*A(2,1)+MDET(3)*A(3,1)
       if(abs(AA).GT.ZERO_TOL) then
          B(1,:) = MDET(:)/AA
          B(2,1) = (A(2,3)*A(3,1)-A(3,3)*A(2,1))/AA
          B(2,2) = (A(3,3)*A(1,1)-A(1,3)*A(3,1))/AA
          B(2,3) = (A(1,3)*A(2,1)-A(2,3)*A(1,1))/AA
          B(3,1) = (A(2,1)*A(3,2)-A(3,1)*A(2,2))/AA
          B(3,2) = (A(3,1)*A(1,2)-A(1,1)*A(3,2))/AA
          B(3,3) = (A(1,1)*A(2,2)-A(2,1)*A(1,2))/AA
       else
          WRITE(*,'('' >>Warning: Zero determinant in 3*3 INVERT'')')
       endif
    else
       WRITE(*,'('' >>Warning: Matrix larger than 3x3 - cannot INVERT!'')')
    endif

  end subroutine INVERT
  
!!!##############################################################################################

  subroutine MAT_VEC(A_VECTOR,B_VECTOR,C_VECTOR,DXRCXI)

    !###    MAT_VEC calculates direction cosines of undeformed material vectors.

    !     Parameter List
    real(dp) :: A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),DXRCXI(3,3)
    !     Local Variables
    real(dp) :: FIBRE_ORIENT(3,3)

    ! Compute components of undeformed orthonormal fibre reference vectors at Gauss pt wrt rc coord system
    call FIBRE_REF_VECS(FIBRE_ORIENT(1,1),FIBRE_ORIENT(1,2),FIBRE_ORIENT(1,3),DXRCXI)
    A_VECTOR(1:3) = FIBRE_ORIENT(1:3,1)
    B_VECTOR(1:3) = FIBRE_ORIENT(1:3,2)
    C_VECTOR(1:3) = FIBRE_ORIENT(1:3,3)
    
  end subroutine MAT_VEC

!!!##############################################################################################

  subroutine MAT_VEC_DEF(nb,ng,AD_VECTOR,BD_VECTOR,CD_VECTOR,xg,ze,zg)

    !##    MAT_VEC_DEF calculates direction cosines of deformed
    !##    normalised orthogonal material vectors at Gauss point ng or at
    !##    XI if ng=0.
    !##    If ng=0 this routine assumes that XE/ZE contain element vertex
    !##    coordinates and microstructural orientations for the
    !##    undeformed/deformed state resp.
    !##    If ng>0 this routine assumes XG contains Gauss pt coordinates,
    !##    microstructural orientations, and derivatives wrt Xi.

    !     Parameter List
    integer :: ng
    real(dp) :: AD_VECTOR(3),BD_VECTOR(3),CD_VECTOR(3),xg(:,:),ze(:,:),zg(:,:)
    !     Local Variables
    integer :: mj,nb,ni,nitb,nj
    real(dp) :: DXDNU(3,3),DZDNU(3,3),DZDX(3,3),SUM
    real(dp) :: DXRCXI(3,3)

    nitb = nit(nb)

!!! Compute undeformed anatomical fibre vectors wrt rc coords at ng. 
!!! DXRCX(njj,ni) = dXRC(njj)/dXI(ni) from XG,
    dxrcxi(1:3,1) = xg(1:3,2)
    dxrcxi(1:3,2) = xg(1:3,4)
    dxrcxi(1:3,3) = xg(1:3,7)

    call MAT_VEC(DXDNU(1,1),DXDNU(1,2),DXDNU(1,3),DXRCXI)   ! Calculate direction cosines

!!! Initialise deformed material vectors with undef material vec dirns
    AD_VECTOR(1:3) = DXDNU(1:3,1)
    BD_VECTOR(1:3) = DXDNU(1:3,2)
    CD_VECTOR(1:3) = DXDNU(1:3,3)
    
!!! Compute the deformation gradient tensor wrt rc coords, i.e. want derivatives of deformed rc coordinates wrt
!!! undeformed rc coordinates.
    call DEFMGRADRC(nb,ng,DZDX,xg,ze,zg)
    
!!! Compute deformed material vectors wrt rc coordinates using the deformation gradient tensor and the undeformed material
!!! vectors wrt rc coordinates
    do nj = 1,3
       do ni = 1,nitb ! we don't actually need the third component
          SUM = 0.0_dp
          do mj = 1,3
             SUM = SUM+DZDX(nj,mj)*DXDNU(mj,ni)
          enddo !nj1
          DZDNU(nj,ni) = SUM
       enddo !ni
       AD_VECTOR(nj) = DZDNU(nj,1)
    enddo !nj
    
    !     Normalise the deformed vectors so that deformed material
    !     coordinates are measures of physical arc length
    call NORMALISE(3,AD_VECTOR)
    
!!! Make the second vector orthogonal to the first (still within the sheet)
    SUM = 0.0_dp
    do nj = 1,3
       SUM = SUM+AD_VECTOR(nj)*DZDNU(nj,2)
    enddo 
    do nj = 1,3
       BD_VECTOR(nj) = DZDNU(nj,2)-SUM*AD_VECTOR(nj)
    enddo !nj
    
    call NORMALISE(3,BD_VECTOR)
    
!!! The third vector is orthogonal to the others (normal to the sheet)
    if(nitb.EQ.3) call CROSS(AD_VECTOR,BD_VECTOR,CD_VECTOR)
    
  end subroutine MAT_VEC_DEF

!!!##############################################################################################

  subroutine CROSS(A,B,C)

    !###    CROSS returns the vector cross product of A*B in C.

    !     Parameter List
    real(dp) :: A(3),B(3),C(3)

    C(1) = A(2)*B(3)-A(3)*B(2)
    C(2) = A(3)*B(1)-A(1)*B(3)
    C(3) = A(1)*B(2)-A(2)*B(1)

  end subroutine CROSS
  
!!!##############################################################################################

  subroutine DEFMGRADRC(nb,ng,DZDX,xg,ze,zg)

    !###    DEFMGRADRC calculates components of the deformation
    !##    gradient tensor wrt rectangular cartesian coords.
    !##    This routines assumes that XG contains the Gauss pt position,
    !##    interpolated material axis orientations, and derivatives
    !##    with respect to Xi coordinates, and that ZE contains deformed
    !##    element vertex coordinates.
      
!     Parameter List
    integer ::  nb,ng
    real(dp) :: DZDX(3,3),xg(:,:),ze(:,:),zg(:,:)
    !     Local Variables
    integer :: mhx,mj,nhx,nitb,nj,NU1(0:3)
    real(dp) :: DETERM,DXIXJ(3,3),DXJXI(3,3),dXref_dXrc,dZrc_dZref,SUM

    DATA NU1/1,2,4,7/

    nitb = 3

!!! Calculate derivatives of Xi wrt Xj (reference) coords, DXIXJ
    dxjxi(1:3,1) = xg(1:3,2)
    dxjxi(1:3,2) = xg(1:3,4)
    dxjxi(1:3,3) = xg(1:3,7)
    call INVERT(nitb,DXJXI,DXIXJ,DETERM)
    
!!! Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
    call ZEZG(1,nb,ng,DXIXJ,ze,zg)
    do nhx = 1,3
       do nj = 1,3
          SUM = 0.0_dp
          do mhx = 1,3
             dZrc_dZref = 0.0_dp
             if(nhx.eq.mhx) dZrc_dZref = 1.0_dp
             do mj = 1,3
                dXref_dXrc = 0.0_dp
                if(mj.eq.nj) dXref_dXrc = 1.0_dp
                SUM = SUM+dZrc_dZref*ZG(mhx,NU1(mj))*dXref_dXrc
             enddo !mj
          enddo !mhx
          DZDX(nhx,nj) = SUM
       enddo !nj
    enddo !nhx

  end subroutine DEFMGRADRC

!!!##############################################################################################

  subroutine FIBRE_REF_VECS(F_VECTOR,G_VECTOR,H_VECTOR,dXRC_dXI)

    !     Parameter List
    real(dp) :: F_VECTOR(3),G_VECTOR(3),H_VECTOR(3),dXRC_dXI(3,3)
    !     Local Variables
    integer :: nj
    logical :: ZERO_F_VECTOR

    !     Initialise all vectors
    F_VECTOR = 0.0_dp
    G_VECTOR = 0.0_dp
    H_VECTOR = 0.0_dp
    
    !     F_VECTOR is the normalised undeformed Xi1 base vector
    F_VECTOR(1:3) = dXRC_dXI(1:3,1)
    call NORMALISE(3,F_VECTOR)
    
    ZERO_F_VECTOR = .TRUE.
    do nj = 1,3
       if(abs(F_VECTOR(nj)).GT.ZERO_TOL) ZERO_F_VECTOR = .FALSE.
    enddo !nj
    if(ZERO_F_VECTOR) then
       !         ...so set G_VECTOR to be the normalised undef Xi2 base vector
       G_VECTOR(1:3) = dXRC_dXI(1:3,2)
       call NORMALISE(3,G_VECTOR)
       !         ...then F_VECTOR is the undeformed Xi2-Xi3 plane normal
       call CROSS(dXRC_dXI(1,2),dXRC_dXI(1,3),F_VECTOR)
       call NORMALISE(3,F_VECTOR)
       !         ...and H_VECTOR lies in the undeformed Xi2-Xi3 plane and is
       !         normal to G_VECTOR
       call CROSS(F_VECTOR,G_VECTOR,H_VECTOR)
    else !F_VECTOR is not all zero
       !         H_VECTOR is the undeformed Xi1-Xi2 plane normal
       call CROSS(dXRC_dXI(1,1),dXRC_dXI(1,2),H_VECTOR)
       call NORMALISE(3,H_VECTOR)
       !         G_VECTOR lies in the undeformed Xi1-Xi2 plane and is
       !         normal to F_VECTOR
       call CROSS(H_VECTOR,F_VECTOR,G_VECTOR)
    endif !ZERO_F_VECTOR
    
  end subroutine FIBRE_REF_VECS

!!!##############################################################################################

  subroutine NORMALISE(NUMCMPTS,VECTOR)
    
    !##    NORMALISE divides the components of VECTOR by it's length
    
    !     Parameter List
    integer :: NUMCMPTS
    real(dp) :: VECTOR(*)
    !     Local Variables
    integer :: ni
    real(dp) :: VECTOR_LENGTH

    VECTOR_LENGTH = 0.0_dp
    do ni = 1,NUMCMPTS
       VECTOR_LENGTH = VECTOR_LENGTH+VECTOR(ni)*VECTOR(ni)
    enddo !ni
    VECTOR_LENGTH = sqrt(VECTOR_LENGTH)
    if(VECTOR_LENGTH.GT.ZERO_TOL) then
       do ni = 1,NUMCMPTS
          VECTOR(ni) = VECTOR(ni)/VECTOR_LENGTH
       enddo !ni
    else
       WRITE(*,'('' >>WARNING: Cannot normalise a zero length vector'')')
    endif
    
  end subroutine NORMALISE

!!!##############################################################################################

  subroutine ZGTG53(az,azl,azu,ri3,CG,TG,xg)
    
    !     Parameter List
    real(dp) :: AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),RI1,RI2,RI3,TG(3,3),xg(:,:)
    !     Local Variables
    integer :: i,mj,nj
    integer,parameter :: nj_pressure = 5
    real(dp) :: BG(3,3),DW(6),W3

    dw = 0.0_dp
    axu = 0.0_dp
    forall(i = 1:3) axu(i,i) = 1.0_dp
    
!!! stress referred to body/fibre Nu coords
    RI3 = AZ
    RI1 = AZL(1,1)+AZL(2,2)+AZL(3,3)
    RI2 =(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3

    bg = 0.0_dp
    forall(i = 1:3) bg(i,i) = ri1 - azl(i,i)
    forall(i = 2:3) bg(i,1) = -azl(i,1)
    bg(3,2) = -azl(3,2)
    
    call derivatives_sedf(CG,DW,RI1,RI2,RI3)

    W3 = RI3*DW(3)
    do mj = 1,3
       do nj = 1,mj
          TG(mj,nj) = 2.0_dp * (DW(1)*AXU(mj,nj) + DW(2)*BG(mj,nj) + W3*AZU(mj,nj))
       enddo
    enddo

    TG(1,2) = TG(2,1)
    TG(1,3) = TG(3,1)
    TG(2,3) = TG(3,2)

    tg = tg - xg(nj_pressure,1) * axu

  end subroutine ZGTG53

!!!##############################################################################################

  subroutine derivatives_sedf(CG,DW,P1,P2,P3)
    
    !##    calculates derivatives of strain energy function wrt
    !##     principal strain invariants.
    
    !***      P1 is the First  principal invariant RI1;
    !***      P2 is the Second principal invariant RI2;
    !***      P3 is the Third  principal invariant RI3;
    !***      DW(1..5) are dW/dI1,dW/dI2,dW/dI3,dW/dK1,dW/dK2.
    
    !     Parameter List
    real(dp) :: CG(NMM),DW(*),P1,P2,P3
    !     Local Variables
    real(dp) :: exp_term,R1,R2,R3

    R1 = P1-3.0_dp
    R2 = P2-3.0_dp
    R3 = P3-1.0_dp
    exp_term = CG(2)/4.0_dp*R1**2+CG(3)/4.0_dp*(R2-2.0_dp*R1)
    DW(1) = CG(1)/4.0_dp * exp(exp_term)*(CG(2)*R1-CG(3))
    DW(2) = CG(1)*CG(3)/8.0_dp*exp(exp_term)
    DW(3) = 0.0_dp
    
  end subroutine derivatives_sedf

!!!##############################################################################################

  subroutine solve_elasticity(niterate,errmax,factor)
    
!***  YP(ny,1) and ZP has current equilibrium solution
!***  YP(ny,2) has prescribed dep var/force increms set by FIX_MECH(ny,1)
!***  YP(ny,3) has prescribed initial equilibrium solution
!***  YP(ny,4) has current set of equilibrium equation residuals
!***  YP(ny,5) is temporary storage 1 used to store solution
!***           increments to add to YP(ny,1)
!***  YP(ny,10) is temporary storage 2 used to store current reference
!***           solution for FE50 cavity elements (set up in IPINI5/
!***           UPSOLU)
    !     Parameter List
    integer :: nb,niterate
    real(dp) :: errmax,factor
    !     Local Variables
    integer :: ntload
    integer :: i,ITER1,nc,no,noload,no_nynr,NWRIT,ny
    real(dp) :: RATIO(0:6),RSUM_SOLINCR
    real(dp) :: xe(nsm,20),xg(20,num),ze(nsm,nhm),zg(nhm,num)
    character :: CHAR1*100
    logical :: ADD_GRAVITY,CONVERgeD,FREE_VAR,OUTPUT,REITER,solve_system,update_matrix

    call evpress ! update the ratio of deformed to underformed at gauss points
    
    NTLOAD = 1
    nb = nb_lung
    
    FREE_VAR = .FALSE.
    !   ***   Solve equations (if at least one dof is not fixed)
    do no_nynr = 1,nynr(0,0,1) !loop over global variables
       ny = nynr(no_nynr,0,1) !is global variable number
       if(.not.FIX_MECH(ny,1)) FREE_VAR = .TRUE.
    enddo !no_nynr
    
    if(NTLOAD.EQ.0) then
       REITER = .TRUE.
       NTLOAD = 1
    else
       REITER = .FALSE.
    endif
    
    OUTPUT = .FALSE.
    NWRIT = 1
    do noload = 1,NTLOAD
       RATIO(0) = 0.0_dp
       OUTPUT = .TRUE.    !output initial residual
       NWRIT = 1
       ADD_GRAVITY = .TRUE.
       
       if(.not.REITER) then
          WRITE(*,'(/'' Load step'',I3/,1X,12(''=''))') noload
          !   ***       Apply increments to displacement and force b.c.s
          i = 0
          do nc = 1,2  !loop over RHS(displ) and LHS(force) vars
             do no_nynr = 1,nynr(0,0,nc) !loop over global variables
                ny = nynr(no_nynr,0,nc) !is global variable number
                if(FIX_MECH(ny,1)) then !ny has incremented essential bdry cond
                   YP(ny,1) = YP(ny,1) + YP(ny,2)*FACTOR
                endif
             enddo            !no_nynr (ny)
          enddo               !nc
       endif                  ! not.REITER
       
       CONVERgeD = .FALSE.

       solve_system = .false.
       if(.not.FREE_VAR) then
          WRITE(*,'(/'' no free variables'')')
       else
          if(niterate.ne.0) then
             solve_system = .TRUE.
          endif
       endif
       
       ITER1 = 0
       
       update_matrix = .TRUE.

       do WHILE (.not.CONVERgeD.and.solve_system)
          call ypzp(1)
          call zprp(nb,xe,xg,ze,zg)
          call calc_conv_ratio(ITER1,ERRMAX,RATIO)
          !   Calc sum of absolute solution vector increments
          RSUM_SOLINCR = 0.0_dp
          if(ITER1.GT.0) then ! after 1st Newton step
             do no = 1,not(2) !loop over global soln variables
                if(nyno(0,no,2).GT.0) then
                   ny = nyno(1,no,2) !is first global variable number
                   RSUM_SOLINCR = RSUM_SOLINCR + abs(YP(ny,5)) !coupled to no
                endif
             enddo            !no
          endif
          
          if(OUTPUT) then
             WRITE(*,'('' Sum of solution increments ='',D11.4)')  RSUM_SOLINCR
          endif

          solve_system = .false.
          converged = .false.
          if((RATIO(0).GT.ERRMAX.and.RSUM_SOLINCR.GT.ZERO_TOL).OR. &
               RSUM_SOLINCR.GT.ERRMAX) then
             if(ITER1.LT.niterate) then
                update_matrix = .TRUE.
                solve_system = .TRUE.
             endif            !ITER1
          else if(RATIO(0).GT.ERRMAX.and.RSUM_SOLINCR.LE.ZERO_TOL) then
             if(ITER1.EQ.0) then ! 1st iteration
                solve_system = .TRUE.
             else
                WRITE(CHAR1,'(D7.1)') ZERO_TOL
                WRITE(*,'(/'' Exiting since sum of solution ' &
                     //'increments < '//CHAR1(1:7)//''')')
             endif
          else
             converged = .TRUE.
          endif               !RATIO(0)
          
          !   ****    If not converged, then solve
          if (.not.converged.and.solve_system) then
             ITER1 = ITER1+1
             OUTPUT = .TRUE.
             
             if(update_matrix) then
                !   ***           Assemble global stiffness matrix gk
                call assemble_gk(nb,xe,xg,ze,zg)
                call solve5(update_matrix) !  Solve global system of equations
               
               !   Update current solution by adding increments
               do no_nynr = 1,nynr(0,0,1) !loop over global vars
                  ny = nynr(no_nynr,0,1) !is global variable number
                  YP(ny,1) = YP(ny,1) + YP(ny,5)
               enddo            !no_nynr
               
               WRITE(*,'(/'' --''/'' Completed iteration number '',i3)') iter1
            endif
         endif               ! converged
      enddo                  !END OF NEWTON STEP
      
      if(converged) then
         WRITE(*,'(/'' Convergence achieved after '',I3,'' iterations'')') ITER1
      else                   ! not converged
         WRITE(*,'(/'' Convergence has not been reached after '',I3,'' iterations'')') ITER1
      endif                  !converged
      
   enddo                     !noload (load step)

   ! update the ratio of deformed to underformed at gauss points
    call evpress
   
 end subroutine solve_elasticity

!!!##############################################################################################

 subroutine calc_conv_ratio(iter1,errmax,ratio)
   
   !     Parameter List
   integer :: ITER1
   real(dp) :: errmax, ratio(0:6)
   !     Local Variables
   integer :: nh,no,no_nynr,no_resid,noy,ny1,ny2,nyo
   real(dp) :: RSUM_CONSTRAINED(0:6),RSUM_UNCONSTRAIN(0:6)
   integer :: count_free,count_constrained
   
   count_free = 0
   count_constrained = 0

   !   reinitialise residual sums and ratios and sum of soln increments
   RSUM_CONSTRAINED(0:6) = 0.0_dp
   RSUM_UNCONSTRAIN(0:6) = 0.0_dp
   RATIO(0:6) = 0.0_dp
   if(ITER1.EQ.0) grr = 0.0_dp

   do no_nynr = 1,nynr(0,1,1) !loop over rows
      ny1 = nynr(no_nynr,1,1) !is row number
      if(nony(0,ny1,1).GT.0) then !free dependent variable
         if(ITER1.EQ.0) then
            do noy = 1,nony(0,ny1,1) !loop on rows assoc with ny1
               no = nony(noy,ny1,1) !is row number for ny1
               grr(no) = grr(no) + YP(ny1,4)
            enddo             !noy
         endif
      else                   !bdry cond applied to dependent variable
         nh = npny(3,ny1,0)
         RSUM_CONSTRAINED(1) = RSUM_CONSTRAINED(1) + abs(YP(ny1,4))
         count_constrained = count_constrained + 1
      endif                  ! free or fixed variable
   enddo                     !no_nynr
   do no = 1,not(1) !loop over global soln rows             
      do nyo = 1,nyno(0,no,1)
         ny1 = nyno(nyo,no,1) !is row number
         ny2 = nyno(nyo,no,2) !is global variable number
         nh = npny(3,ny1,0)
         RSUM_UNCONSTRAIN(1) = RSUM_UNCONSTRAIN(1)+ abs(grr(no))
         count_free = count_free + 1
      enddo                  !nyo
   enddo                     !no
   
   do no_resid = 1,6
      if(RSUM_CONSTRAINED(no_resid).GT.ZERO_TOL) then !only include the ratio of a particular residual type if some d.o.f are constrained
         RATIO(no_resid) =  RSUM_UNCONSTRAIN(no_resid)/RSUM_CONSTRAINED(no_resid)
      else if(RSUM_UNCONSTRAIN(no_resid).LT.ZERO_TOL) then !If both unconstrained and constrained are practically zero
         RATIO(no_resid) = 0.0_dp ! needs to be smaller than zero_tol as they get added up
      else                   !If unconstrained resid is high and constrained is zero
         RATIO(no_resid) = 2.0_dp*ERRMAX
      endif
      RATIO(0) = RATIO(0)+RATIO(no_resid)
   enddo
   
   WRITE(*,'('' Sum of (unconstr/constr) ratios of the degrees of freedom='',D11.4)') RATIO(0)   
   
 end subroutine calc_conv_ratio

!!!##############################################################################################

 subroutine ZPRP(nb,xe,xg,ze,zg)
   
   !##    ZPRP calculates global residual vector YP(ny,4) at current
   !##    solution.  RE(ns,nh) has been corrected with scaling factor
   !##    SE(ns,nb,ne)

   integer :: nb
   real(dp) :: xe(:,:),xg(:,:),ze(:,:),zg(:,:)
   !     Local Variables
   integer :: nc,ne,no_nynr,ny,ny1
   integer :: lge(192,2),nh,nhst(2),nhs,ns
   real(dp) :: re(nsm,nhm)

   nc = 1 !LHS
   yp(:,4) = 0.0_dp    !   Initialise residuals

    do ne = 1,tissue_num_elems
      call MELge(lge,nb,nc,ne,nhst)
      call xpxe(nb,ne,xe)  !     put XP into XE
      call zpze(nb,ne,ze)  !     put ZP into ZE
      call cpcg(ne)
      call ZERE50(nb,ne,re,xe,xg,ze,zg)  !     get element matrix from ZE
      nhs = 0
      do nh  =  1,3
         do ns = 1,NST(1)
            nhs = nhs+1
            ny = IABS(Lge(nhs,1)) !row number
            YP(ny,4) = YP(ny,4) + RE(ns,nh)
         enddo                  !ns
      enddo                     !nhx
   enddo                     !noelem (ne)
   
   !   *** Add in global loads 
   do no_nynr = 1,nynr(0,1,1) !loop over rows
      ny = nynr(no_nynr,1,1) !row number
      ny1 = geTNYR(2,0,1,ny) !is RHS var #
      if(npny(3,ny,0).LE.3) then !ny is a force/moment eqn
         YP(ny,4) = YP(ny,4) - YP(ny1,1)
      endif
   enddo !no_nynr

 end subroutine ZPRP


!!!##############################################################################################

 subroutine melge(lge,nb,nc,ne,nhst)

   !##    MELge calculates the row numbers (Lge(*,1)) and column numbers
   !##    (Lge(*,2)) in the global matrix nc for element variables nhs
   !##    in region nr. It also returns the total number of element
   !##    variables nhst(nrc).

   !     Parameter List
   integer :: lge(:,:),nb,nc,ne,nhst(2)
   !     Local Variables
   integer :: nh,nn,np,nrc

   do nrc = 1,2
      nhst(nrc) = 0
      do nh = 1,3
         !   !!!     Use the LHS (nc=1) basis to determine the # of equations
         do nn = 1,nnt(nb)     !nodal variables
            np = npne(nn,nb,ne)
            nhst(nrc) = nhst(nrc)+1
            Lge(nhst(nrc),nrc) = nynp(nh,np,nrc,nc)
         enddo               !nn
      enddo                  !nh
   enddo                     !nrc

 end subroutine melge

!!!##############################################################################################

 subroutine zere50(nb,ne,re,xe,xg,ze,zg)

   !     Parameter List
   integer :: nb,ne
   real(dp) :: re(:,:),xe(:,:),xg(:,:),ze(:,:),zg(:,:)
   !     Local Variables
   integer :: JP,ng,nh,nitb,ns,NU1(0:3)
   real(dp) :: Age,AZ,AZL(3,3),AZU(3,3),DXIX(3,3),DZDX(3,3),GXL(3,3), &
        GXU(3,3),PPGG(4),rgx,RI3,RWG,TG(3,3),Volume,ZG_temp(NHM,NUM)
   real(dp) :: yg_ne(2,ngm)
   CHARACTER STRESSTYPE*17
   DATA STRESSTYPE/' '/
   DATA NU1/1,2,4,7/

   nitb = nit(nb)
      
   RE = 0.0_dp
   
   !   *** Main Gauss point loop
   do ng = 1,NGT(nb)

      !   Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
      call XEXG(nb,ng,xe,xg)
      !   stresses referred to Nu in constitutive law.
      JP = 1
      !   Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
      !   derivs (DXIX) of Xi wrt Xj (JP=0) or Nu (JP=1) coords.
      call XGMG(JP,nitb,nb,DXIX,GXL,GXU,RGX,xg)
      !   Calculate the Jacobian for integration wrt undef coords:
      RWG = RGX*wg(ng)
               
      !   Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
      call ZEZG(1,nb,ng,DXIX,ze,zg)
      call DEFMGRADRC(nb,ng,DZDX,xg,ze,zg_temp)
      Volume = DET(DZDX)
      if(Volume .LT. 0.0_dp) then
         WRITE(*,'('' >>Warning: Volume at ng='',I5,'' ne='',I5,'' less than zero'')') ng,ne
         WRITE(*,'('' DET(DZDX)='',D12.4)') Volume
         !read(*,*)
      endif
         
      !   Calculate deformed metric tensors wrt Nu (AZL,AZU)
      call ZGMG(nb,AZ,AZL,AZU,zg)
      !   Get contravariant cpts of 2nd Piola-Kirchhoff stress
      !   tensor (TG) wrt undeformed Nu coordinates
      yg_ne(:,:) = yg(:,:,ne)
      call ZGTG53(az,azl,azu,ri3,CG,TG,xg)
      
      !   ***   Main element residual
      do nh = 1,3
         do ns = 1,NST(nb) !element variables
               
            PPGG(1) = PG(ns,1,ng)
            PPGG(2) = PG(ns,NU1(1),ng)*DXIX(1,1) + PG(ns,NU1(2),ng)*DXIX(2,1) + &
                 PG(ns,NU1(3),ng)*DXIX(3,1)
            PPGG(3) = PG(ns,NU1(1),ng)*DXIX(1,2) + PG(ns,NU1(2),ng)*DXIX(2,2) + &
                 PG(ns,NU1(3),ng)*DXIX(3,2)
            PPGG(4) = PG(ns,NU1(1),ng)*DXIX(1,3) + PG(ns,NU1(2),ng)*DXIX(2,3) + &
                 PG(ns,NU1(3),ng)*DXIX(3,3)
               
            Age = (TG(1,1)*ZG(nh,2)+TG(1,2)*ZG(nh,4)+TG(1,3)*ZG(nh,7))*PPGG(2) &
                 +(TG(2,1)*ZG(nh,2)+TG(2,2)*ZG(nh,4)+TG(2,3)*ZG(nh,7))*PPGG(3) &
                 +(TG(3,1)*ZG(nh,2)+TG(3,2)*ZG(nh,4)+TG(3,3)*ZG(nh,7))*PPGG(4)
            !   Main residual
            RE(ns,nh) = RE(ns,nh)+Age*RWG
         enddo               !ns
      enddo                  !nh
      
      !   ***   Gravity
      do nh = 1,3
         do ns = 1,NST(nb) !element variables
            Age = CG(IL_density,ng)*gravity(nh)*PG(ns,1,ng)
            RE(ns,nh) = RE(ns,nh)+Age*RWG !   Main residual
         enddo               !ns
      enddo                  !nh
      
   enddo                  !end of ng loop
   
 end subroutine zere50


!!!##############################################################################################

 subroutine assemble_gk(nb,xe,xg,ze,zg)

   !##    (CMISS ASSEMBLE: generates element stiffness matrix ES for nonlinear
   !##    problems by ZEES, and assembles into the global stiffness
   !##    matrix gk.

   !     Parameter List
   integer :: nb
   real(dp) :: xe(:,:),xg(:,:),ze(:,:),zg(:,:)
   !     Local Variables
   integer :: ne,nz
   integer lge(NHM*NSM,2)
   integer :: nhst(2),nhs1,nhs2,ny1,ny2
   real(dp) :: ES(NHM*NSM,NHM*NSM)

   gk = 0.0_dp
   call ypzp(1) ! transfer yp to zp
      
    do ne = 1,tissue_num_elems
      call melge(lge,nb,1,ne,nhst) ! gives lge and nhst
      es = 0.0_dp     
      call xpxe(nb,ne,xe)
      call zpze(nb,ne,ze)
      call zees(Lge,nb,ne,es,xe,xg,ze,zg) ! calculate element matrix es
      !**  Assemble element stiffness matrix into global system.
      do nhs1 = 1,nhst(1)
         ny1 = IABS(Lge(nhs1,1)) ! row
         do nhs2 = 1,nhst(2)
            ny2 = IABS(Lge(nhs2,2)) ! column
            call sparse(ny1,ny2,nyt,nz,nz_gk_m) ! nz for row,col
            gk(nz) = gk(nz)+ES(nhs1,nhs2)
         enddo               !nhs2
      enddo                  !nhs1
   enddo

 end subroutine assemble_gk

!!!##############################################################################################

 subroutine zees(lge,nb,ne,es,xe,xg,ze,zg)
   
   !##    ZEES calculates element tangent stiffness matrix ES from
   !##    current dependent variable array ZE.
   !**  Finite difference calculation of tangent stiffness matrix
   
   !     Parameter List
   integer :: Lge(:,:),nb,ne
   real(dp) :: ES(:,:),xe(:,:),xg(:,:),ze(:,:),zg(:,:)
   !     Local Variables
   integer :: mh,mhs,ms,my,my1,nh,nhs,ns,ny,ny1
   real(dp) :: ZE_STORE,re1(nsm,nhm),re2(nsm,nhm)
   real(dp),parameter :: delta = 1.0e-4_dp
   
   call cpcg(ne)
   
   call ZERE50(nb,ne,re1,xe,xg,ze,zg)
   
   nhs = 0
   do nh = 1,3
      do ns = 1,NST(nb)
         nhs = nhs+1
         ny = IABS(Lge(nhs,2)) !local variable number
         ny1 = geTNYR(1,0,2,ny) !global variable #
         if(.not.FIX_MECH(ny1,1)) then
            
            ZE_STORE = ZE(ns,nh)                 ! Store current solution
            ZE(ns,nh) = ZE(ns,nh)+DELTA          ! Normal perturbation for non-isochoric interpolation
            call ZERE50(nb,ne,re2,xe,xg,ze,zg) ! Evaluate perturbed residual
            ZE(ns,nh) = ZE_STORE                 ! Return ZE to original value
            !   ***         Assemble element stiffness matrix
            mhs = 0
            do mh = 1,3
               do ms = 1,NST(nb)
                  mhs = mhs+1
                  my = IABS(Lge(mhs,1)) !local variable number
                  my1 = geTNYR(1,0,2,my) !global var #
                  if(.not.FIX_MECH(my1,1)) then  ! Note: DELTA here is global delta
                     ES(mhs,nhs) = (RE2(ms,mh)-RE1(ms,mh))/DELTA
                  endif
               enddo
            enddo
         endif               !FIX_MECH
      enddo                  !ns
   enddo                     !nh
   
 end subroutine zees
 
 
!!!##############################################################################################
 
 subroutine sparse(I,J,N,nz,NZMAX)
   
   !     Parameter List
   integer :: I,J,N,nz,NZMAX
   
   !No sparsity
   nz = I+(J-1)*N
   if(nz.GT.NZMAX.OR.nz.LT.1) then
      nz = 0
      write(*,'(''>>Array coordinates outside range[0]'')')
      read(*,*)
   endif
   
 end subroutine sparse
   
 
!!!##############################################################################################
 
 subroutine solve5(update_matrix)
   
   !##    Solves the resulting system of linear equations for the
   !##    increments of the nonlinear solver in NONLIN.
   
   !     Parameter List
   logical :: update_matrix
   !     Local Variables
   integer :: no1,no2,no_nynr1,no_nynr2,noy1,noy2,ny1,ny2,ny3,nyo2,nz,nzz

   if(update_matrix) gkk = 0.0_dp
   grr = 0.0_dp
   
   !   *** Calculate global RHS vector
   
   do no_nynr1 = 1,nynr(0,1,1) !loop over rows
      ny1 = nynr(no_nynr1,1,1) !is row number
      GR(ny1) = YP(ny1,4)
   enddo

   do no_nynr1 = 1,nynr(0,1,1) !loop over rows of gk
      ny1 = nynr(no_nynr1,1,1) !is row number
      do noy1 = 1,nony(0,ny1,1)
         no1 = nony(noy1,ny1,1) !solution row # attached to ny1
         if(update_matrix) then
            do no_nynr2 = 1,nynr(0,2,1) !loop over local columns of gk
               ny2 = nynr(no_nynr2,0,1) !global column #
               ny3 = nynr(no_nynr2,2,1) !local column #
               call sparse(ny1,ny3,nyt,nz,nz_gk_m)
               if(nz.NE.0) then
                  do noy2 = 1,nony(0,ny2,2)
                     no2 = nony(noy2,ny2,2) !solution var # attached to ny2
                     call sparse(no1,no2,not(1),nzz,nz_gkk_m)
                     if(nzz.NE.0) gkk(nzz) = gkk(nzz)+gk(nz) 
                  enddo      !noy2
               endif
            enddo            !no_nynr2
         endif
         grr(no1) = grr(no1)+GR(ny1)
      enddo                  !noy1
   enddo                     !no_nynr1
   
   !** Solve reduced system

   call solve_system(not(1),not(2),gkk,grr,XO,firsts,update_matrix)
   
   yp(:,5) = 0.0_dp
   
   do no2 = 1,not(2)
      do nyo2 = 1,nyno(0,no2,2) !must be a sep loop since adding to YP
         ny2 = nyno(nyo2,no2,2)
         YP(ny2,5) = YP(ny2,5) - XO(no2) 
      enddo !nyo2
   enddo                     !no2  
   
 end subroutine solve5
 
!!!##############################################################################################
 
 subroutine solve_system(LDA,N,A,B,X,FIRST_A,UPDATE_A)
   
   !     Parameter List
   integer,intent(in) :: LDA,N
   real(dp) :: A(*),B(N),X(N)
   logical :: FIRST_A,UPDATE_A
   !     Local Variables
   integer :: ITER = 100, nres = 100
   real(dp) :: anorm, RESID = 0.10e-5_dp
   
   !   Free/Allocate memory as necessary
   if(FIRST_A) then
      if(allocated(rsolv1)) deallocate(rsolv1)
      allocate(rsolv1(N))
   else if(UPDATE_A) then
      !        call FREE_SOLVER(NX,FIRST_A)
      !deallocate(rsolv1)
   endif
   
   !   Factorise the system 
   call ITER_FACTOR_x(A,LDA,N,rsolv1,ANORM)
   FIRST_A = .FALSE.
   
   !   Solve the problem
   x(1:n) = 0.0_dp
   ITER = 100
   RESID = 0.1e-5_dp

   call gmres_x(A,LDA,N,X,B,rsolv1,ANORM,RESID,ITER,NRES)
   
 end subroutine solve_system
 
!!!##############################################################################################
 
 subroutine ITER_FACTOR_x(A,LDA,N,D,ANORM)
   
   !##    ITER_FACTOR factorises a system that is to be used in
   !##    one of the iterative solvers
   
   !     Parameter List
   integer :: LDA,N
   real(dp) :: A(*),D(*),ANORM
   !     Local variables
   integer :: I
   real*8 sum
   
   do I = 1,N
      if(abs(A(I+N*(I-1))).le.zero_tol) then
         D(I) = 1.0_dp          ! Kludge. We should die here
      else
         D(I) = 1.0_dp/A(I+LDA*(I-1))
      endif
   enddo
   
   !   Get the norm of A
   SUM = 0.0_dp
   do I = 1,N
      SUM = SUM+norm3(N,A(I),LDA)**2
   enddo
   ANORM = sqrt(SUM)
   
 end subroutine ITER_FACTOR_x
 
 
!!!##############################################################################################

 subroutine GMRES_x(A,LDA,N,X,B,D,ANORM,RESID,ITER,NRES)
   
   !##    GRMRES solves a system of linear equations using the preconditioned
   !##    Generalised Minimum Residual scheme.
   
   !     Parameter List
   integer :: LDA,N,ITER,NRES
   real(dp) :: A(*),X(N),B(N),D(N),ANORM,RESID
   !     Local Variables
   real(dp),allocatable :: R(:),V(:),S(:),T(:),U(:),Y(:),H(:)

   allocate(R(N))
   allocate(V(N*(NRES+1)))
   allocate(S(NRES+1))
   allocate(T(NRES+1))
   allocate(U(NRES+1))
   allocate(Y(NRES+1))
   allocate(H(NRES*(NRES+1)))

   call GMRES_SUB_x(A,LDA,N,X,B,D,R,V,S,T,U,Y,H,anorm,resid,iter,nres) 

   deallocate(H)
   deallocate(Y)
   deallocate(T)
   deallocate(S)
   deallocate(V)
   deallocate(U)
   deallocate(R)

 end subroutine GMRES_x

!!!##############################################################################################

 subroutine GMRES_SUB_x(A,LDA,N,X,B,D,R,V,S,T,U,Y,H,ANORM,RESID,ITER,NRES)
   
   !##  Solves a system of linear equations using the preconditioned
   !##  Generalised Minimum Residual scheme.

   !     Parameter List
   integer :: LDA,N,ITER,NRES
   real(dp) :: A(*),X(N),B(N),D(N),R(N),V(N,NRES+1)
   real(dp) :: S(NRES+1),T(NRES+1),U(NRES+1),Y(NRES+1),H(NRES+1,NRES)
   real(dp) :: ANORM,RESID
   !     Local Variables
   integer :: I,J,L,MAXIT
   real(dp) :: TOL,ABSV1,ALPHA,BETA,SIGMA,TMP,BNORM,SCALED_TOL,EPS
   real(dp) :: v1_temp(N),v2_temp(N)
   
   TOL = RESID
   MAXIT = ITER
   BNORM = norm2(b)
   eps = zero_tol **2
   SCALED_TOL = TOL*BNORM
   SIGMA = 1.0_dp
   r = 0.0_dp
   v = 0.0_dp
   
   !   Check for null system: if we have a zero RHS we exit.
   if(BNORM.LE.0.0_dp) then
      ITER = 0
      RESID = 0.0_dp
      x = 0.0_dp
      RETURN
   endif
   
   !   Solve the system
   ITER = 1
   do WHILE(ITER.LE.MAXIT)
      
      ! Calculate the residual r = b - A.x
      r(1:n) = b(1:n)
      call dgemv_x(N,N,-1.0_dp,A,LDA,X,1.0_dp,R)
      RESID = norm2(r) !DNRM2_x(N,R,1)
      
      forall(j = 1:n) v(j,1) = r(j)*d(j)   ! precondition
      
      ABSV1 = norm2(v(:,1)) !DNRM2_x(N,V(1,1),1)
      if(abs(ABSV1).le.abs(RESID*EPS)) then
         WRITE(*,'(''>>Sigma = Inf '')')
         read(*,*)
      endif
      if(ITER.EQ.1) then
         SIGMA = RESID/ABSV1
         SCALED_TOL = TOL*(ANORM* norm2(x) + bnorm) ! DNRM2_x(N,X,1) + BNORM)
      endif
      
      if(resid.lt.scaled_tol) goto 200  ! solution converged
      
      S(1) = ABSV1
      TMP = 1.0_dp/ABSV1
      v = v * tmp
      
      do I = 1,NRES !  Inner loop
         call dgemv_x(N,N,1.0_dp,A,LDA,v(1,I),0.0_dp,R) ! mat/vec multiplication of A and v(1,i)
         forall(j = 1:n) v(j,i+1) = r(j)*d(j)             ! equivalent to a preconditioning
         do J = 1,I
            v1_temp(1:n) = v(1:n,i+1)
            v2_temp(1:n) = v(1:n,j)
            H(J,I) = dot_product(v1_temp,v2_temp)
            v(1:n,i+1) = v(1:n,i+1) - h(j,i) * v(1:n,j)
         enddo
         H(I+1,I) =  norm2(v(:,i+1)) !DNRM2_x(N,V(1,I+1),1)
         if(abs(H(I+1,I)).LE.EPS) then
            WRITE(*,'(''>>H_i+1,i = 0 ('')')
            read(*,*)
         endif
         TMP = 1.0_dp/H(I+1,I)
         v(1:n,i+1) = v(1:n,i+1) * tmp
         do J = 1,I-1
            ALPHA = T(J)*H(  J,I)+U(J)*H(J+1,I)
            BETA  = T(J)*H(J+1,I)-U(J)*H(  J,I)
            H(J,I) = ALPHA
            H(J+1,I) = BETA
         enddo
         TMP = sqrt(H(I,I)**2+H(I+1,I)**2)
         if(abs(TMP).LE.EPS) then
            WRITE(*,'(''>>H_i,i^2 + H_i+1,i^2 = 0 '')')
            read(*,*)
         endif
         T(I) = H(I,I)/TMP
         U(I) = H(I+1,I)/TMP
         H(I,I) = TMP
         H(I+1,I) = 0.0_dp
         S(I+1) = -(U(I)*S(I))
         S(I) = S(I)*T(I)
         RESID = SIGMA * abs(S(I+1))
         if(RESID.LT.SCALED_TOL.OR.ITER.ge.MAXIT) GOTO 100
         ITER = ITER+1
      enddo
      I = I-1
100   CONTINUE

      do J = I,1,-1
         if(abs(H(J,J)).LE.abs(S(J)*EPS)) then
            WRITE(*,'(''>>Y(J) = Inf '')')
            read(*,*)
         endif
         Y(J) = S(J)/H(J,J)
         do L = J-1,1,-1
            S(L) = S(L)-Y(J)*H(L,J)
         enddo
      enddo
      
      do J = I,1,-1
         x(1:n) = x(1:n) + y(j) * v(1:n,j)
      enddo
      if(ITER.EQ.1) SCALED_TOL = TOL*(ANORM* norm2(x) + bnorm) !DNRM2_x(N,X,1) + BNORM)
      if(RESID.LT.SCALED_TOL.OR.ITER.ge.MAXIT) GOTO 200
      
   enddo
   ITER = MAXIT
   
200 CONTINUE
   
   RETURN
 end subroutine GMRES_SUB_x

!!!##############################################################################################

 subroutine dgemv_x (M, N, ALPHA, A, lda,X, BETA, Y)
   !  matrix-vector operation y := alpha*A*x + beta*y.

   real(dp) :: ALPHA, BETA
   integer :: lda,M, N
   real(dp) :: A(lda,*), X(*), Y(*)
   integer :: I, J

   if( abs(BETA - 1.0_dp).gt.zero_tol) y(1:n) = 0.0_dp
   
   !  y := alpha*A*x + y.
   do J = 1, N
      if( abs(X(j)).gt.zero_tol)then
         do I = 1, M
            Y( I ) = Y( I ) + alpha*x(j)*A( I, J )
         enddo
      end IF
   enddo

 end subroutine dgemv_x
      
!!!##############################################################################################

 subroutine project_lung_to_cavity(nelist,nplist,option)
   
   !     Parameter List
   integer :: num_ne,nelist(:),num_np,nplist(:)
   character(len = *) :: option
   !     local
   integer :: ne1,ne_project,nitb,nj,noelem,nonode,np2
   real(dp) :: point(3),sqnd,sq_dist,xe(nsm,njm),xe_nj(nsm),xi(2),xi_project(3)
   logical :: found

   num_ne = count(nelist.ne.0)
   num_np = count(nplist.ne.0)
   
   nitb = nit(nb_cavity)              ! the number of Xi coordinates in contact cavity (2)
   
   do nonode = 1,num_np               ! loop over lung surface nodes
      np2 = nplist(nonode)            ! lung surface node to project
      point(1:3) = ZP(1:3,np2,1)      ! coordinates of the surface node to project
      sq_dist = 1.0e+10_dp            ! initialise the minimum distance to surface
      noelem = 0
      do noelem = 1,num_ne            ! for a list of contact cavity elements
         ne1 = nelist(noelem)         ! element # for contact cavity
         xi = 0.5_dp ! initialise
         call xpxe_cavity(nb_cavity,ne1,xe)  ! surface node info into xe
         found = .true.               ! means that does not have to be orthogonal
         call project_closest(nb_cavity,sqnd,xe,xi,point,found)
         if(sqnd.lt.sq_dist)then
            sq_dist = sqnd
            ne_project = ne1
            xi_project(1:nitb) = xi(1:nitb) ! xi_project is xi on contact surface
         endif
      enddo
      
      select case (trim(option))
      case ('no_edge')
      case ('xi1_0_edge')
         xi_project(1) = 0.0_dp
      case ('xi1_1_edge')
         xi_project(1) = 1.0_dp
      case ('xi2_0_edge')
         xi_project(2) = 0.0_dp
      case ('xi2_1_edge')
         xi_project(2) = 1.0_dp
      end select
      
      call xpxe_cavity(nb_cavity,ne_project,xe)  ! surface node info into xe

      do nj = 1,3
         xe_nj(:) = xe(:,nj)
         ! update the deformed lung mesh coordinates from the xi location on contact element
         ZP(nj,np2,1) = pxi(nb_cavity,1,xi_project,xe_nj) ! pxi = function(xi, contact_element)
      enddo                  !nj
   enddo                     ! nonode
      
   call zpyp(1) ! update yp for lung mesh from zp
   call zpyp(5) ! update yp for lung mesh from zp
   
 end subroutine project_lung_to_cavity

      
!!!##############################################################################################
  
  subroutine project_closest(nb,SQ,XE,XI,point,INELEM)

    !##    Finds the XI-coordinates at the closest approach of a 2D
    !##    element to a data point with coordinates XD using a modified
    !##    Newton algorithm.
    !##    If INELEM is true then the closest point within the element is
    !##    returned; if false then the projection from the element to the
    !##    coordinate must be orthogonal.  If such a projection can't be
    !##    obtained the xi position returned is outside the element.
    
    !*** ITMAX is the maximum number of iterations.
    !*** TOL is the required tolerance of the solution.
    !*** VMAX denotes the maximum step length per iteration.
    !*** notE: Works only for 2-D elements in rectangular cartesian coords.

    !     Parameter List
    integer :: nb
    real(dp) :: SQ,XE(:,:),XI(:),point(:)
    logical :: INELEM
    !     Local Variables
    integer :: BOUND(2),it,it2,ni,nifix,nj
    integer,parameter :: itmax = 10
    real(dp) :: DELTA,DET,D2SQV2,D2SQVW2,D2SQXI(2,2),D2ZXI(3,2,2),DSQXI(2), &
         DSQXI1,DSQXI2,DSQV,DSQVW,DZ(3),DZXI(3,2),EVMIN,EVMAX,H(2), &
         MU,SQLIN,SQDIFF,SQdpRED,TEMP,TEMP1,TEMP2,TOL, &
         TOL2,V(2),V1,V2,W,xe_nj(nsm),XILIN(2),Z(3)
    logical :: converged,ENFORCE(2),FREE,NEWTON

    DELTA = 0.25_dp !VMAX/4.0_dp
    TOL = 5.0_dp*LOOSE_TOL !must be > sqrt(eps) or SQLIN<=SQ check may not work
    TOL2 = TOL**2
    SQ = 0.0_dp
    do nj = 1,3
       xe_nj(:) = xe(:,nj)
       Z(nj) =  pxi(nb,1,XI,xe_nj)
       DZ(nj) = Z(nj)- point(nj)
       SQ = SQ+DZ(nj)**2
    enddo

    IT = 0
    converged = .FALSE.
    do WHILE(.not.converged.and.IT.LT.ITMAX)
       dsqxi = 0.0_dp
       do nj = 1,3
          xe_nj(:) = xe(:,nj)
          DZXI(nj,1)= pxi(nb,2,XI,xe_nj)
          DZXI(nj,2) =  pxi(nb,4,XI,xe_nj)
          DSQXI(1) = DSQXI(1)+DZXI(nj,1)*DZ(nj)
          DSQXI(2) = DSQXI(2)+DZXI(nj,2)*DZ(nj)
       enddo
       do ni = 1,2
          if(abs(XI(ni)).le.zero_tol) then
             BOUND(ni) = 1
             ENFORCE(ni) =  DSQXI(ni).ge.0.0_dp
          else if(abs(XI(ni)-1.0_dp).le.zero_tol) then
             BOUND(ni) = -1
             ENFORCE(ni) =  DSQXI(ni).LE.0.0_dp
          else
             BOUND(ni) = 0
             ENFORCE(ni) = .FALSE.
          endif
       enddo
       if(ENFORCE(1).and.ENFORCE(2)) GO TO 9998

       d2sqxi = 0.0_dp
       do nj = 1,3
          xe_nj(:) = xe(:,nj)
          D2ZXI(nj,1,1) =  pxi(nb,3,XI,xe_nj)
          D2ZXI(nj,1,2) =  pxi(nb,6,XI,xe_nj)
          D2ZXI(nj,2,2) =  pxi(nb,5,XI,xe_nj)
          D2SQXI(1,1) =  D2SQXI(1,1)+DZXI(nj,1)*DZXI(nj,1)+D2ZXI(nj,1,1)*DZ(nj)
          D2SQXI(1,2) =  D2SQXI(1,2)+DZXI(nj,1)*DZXI(nj,2)+D2ZXI(nj,1,2)*DZ(nj)
          D2SQXI(2,2) =  D2SQXI(2,2)+DZXI(nj,2)*DZXI(nj,2)+D2ZXI(nj,2,2)*DZ(nj)
       enddo
       !     A Newton step is taken if the condition of the Hessian
       !     guarantees that the step will be within the trust region.
       !     Otherwise the Hessian is shifted towrds a diagonal matrix to
       !     shift the step towards steepest descent.  Usually it is a much
       !     better direction than steepest descent.  I think it is close to
       !     the best direction in the trust region.
       !**    Find the smallest eigen value of the Hessian.
       DSQXI2 = DSQXI(1)**2+DSQXI(2)**2
       DSQXI1 = sqrt(DSQXI2)
       TEMP1 = (D2SQXI(1,1)+D2SQXI(2,2))/2.0_dp
       TEMP2 = sqrt(((D2SQXI(1,1)-D2SQXI(2,2))/2.0_dp)**2+D2SQXI(1,2)**2)
       EVMIN = TEMP1-TEMP2
       EVMAX = TEMP1+TEMP2
       if(DSQXI1.LT.TOL2) GO TO 9998
       do it2 = 1,ITMAX
          if(abs(delta).gt.zero_tol)  TEMP = DSQXI1/DELTA
          NEWTON = EVMIN.ge.TEMP
          if(NEWTON) then !Newton is safe
             H(1) = D2SQXI(1,1)
             H(2) = D2SQXI(2,2)
             DET = EVMIN*EVMAX
          else
             !**        Shift eigenvalues to restrict step
             MU = TEMP-EVMIN
             H(1) = D2SQXI(1,1)+MU
             H(2) = D2SQXI(2,2)+MU
             DET = TEMP*(EVMAX+MU)
          endif
          V(1) = -(H(2)*DSQXI(1)-D2SQXI(1,2)*DSQXI(2))/DET
          V(2) = (D2SQXI(1,2)*DSQXI(1)-H(1)*DSQXI(2))/DET
          V2 = V(1)**2+V(2)**2
          DSQV = DSQXI(1)*V(1)+DSQXI(2)*V(2)
          !       This checks that numerical errors have not
          !       prevented the step direction being a descent direction.
          
          if(DSQV**2.LT.DSQXI2*V2*TOL2) then !try a smaller trust region
             DELTA = DELTA/10.0_dp
          else !step is good, check feasible and limit step size
             FREE = .TRUE.
             do ni = 1,2
                if(BOUND(ni).NE.0.and.(BOUND(ni).GT.0.EQV.V(ni).LT.0.0_dp)) then
                   FREE = .FALSE.
                   nifix = ni
                endif
             enddo
             W = 1.0_dp
             if(FREE) then
                V1 = sqrt(V2) !currently < DELTA
                D2SQV2 = V(1)*(V(1)*D2SQXI(1,1)+2.0_dp*V(2)*D2SQXI(1,2))+V(2)**2*D2SQXI(2,2)
                if(.not.NEWTON) then ! Try to step to estimate of minimum along line
                   if(V1.GT.0.0_dp) then
                      W = DELTA/V1
                      if(D2SQV2.GT.0.0_dp) then !minimum exists
                         W = DMIN1(W,-DSQV/D2SQV2) !minimum if within trust region
                      endif
                   endif
                endif !newton
             else
                if(ENFORCE(2)) then !gradient suggests must use ni=1
                   nifix = 2
                else if(ENFORCE(1)) then !gradient suggests must use ni=2
                   nifix = 1
                endif
                ni = 3-nifix
                if(.not.INELEM) then
                   !**            If stepping predominantly out of element then exit
                   nifix = 3-ni
                   if(abs(V(nifix)).GT.abs(DSQXI(ni)/H(ni))) then
                      XI(nifix) = XI(nifix)+V(nifix)
                      GO TO 9998
                   endif
                endif
                V(nifix) = 0.0_dp
                if(D2SQXI(ni,ni).GT.0.0_dp) then !minimum exists
                   V(ni) = -DSQXI(ni)/D2SQXI(ni,ni)
                   V1 = abs(V(ni))
                   NEWTON = V1.LE.DELTA
                endif
                if(.not.NEWTON) then
                   V(ni) = -DSIGN(DELTA,DSQXI(ni))
                   V1 = DELTA
                endif
                V2 = V1*V1
                DSQV = DSQXI(ni)*V(ni)
                D2SQV2 = V2*D2SQXI(ni,ni)
             endif !free
             !**        First half of convergence test.
             !         Should be before boundary colllision check
             converged = V1*W.LT.TOL
             !**        Try the step.  (Name: XILIN is historical)
             XILIN(1:2) = XI(1:2)+V(1:2)*W
             !**        Test for boundary collision
             do ni = 1,2
                if(XILIN(ni).LT.0.0_dp) then
                   XILIN(ni) = 0.0_dp
                   W = XI(ni)/(-V(ni))
                   XILIN(3-ni) = XI(3-ni)+V(3-ni)*W
                else if(XILIN(ni).GT.1.0_dp) then
                   XILIN(ni) = 1.0_dp
                   W = (1.0_dp-XI(ni))/V(ni)
                   XILIN(3-ni) = XI(3-ni)+V(3-ni)*W
                endif
             enddo !ni
             !**        Calculate new distance
             SQLIN = 0.0_dp
             do nj = 1,3
                xe_nj(:) = xe(:,nj)
                Z(nj) = pxi(nb,1,XILIN,xe_nj)
                DZ(nj) = Z(nj) - point(nj)
                SQLIN = SQLIN+DZ(nj)**2
             enddo
             !**        Second half of convergence test.
             converged = converged.and.abs(SQ-SQLIN)/(1.0_dp+SQ).LE.TOL
             if(converged) GO TO 5
             DSQVW = DSQV*W !<0
             D2SQVW2 = 0.5_dp*D2SQV2*W*W !1/2 for computational efficiency
             SQDIFF = 0.5_dp*(SQLIN-SQ) !1/2 because derivs are for SQ/2
             SQDPRED = SQDIFF-DSQVW-D2SQVW2
             !**        Exit loop if decrease is satisfactory
             if(SQDIFF.LE.0.25_dp*DSQVW) then
                if(NEWTON) then
                   DELTA = V1 !next step smaller unless this is increased
                else if(W.ge.1.0_dp) then
                   !             If the quadratic model is good increase trust region size
                   if(SQDPRED.LT.-0.1_dp*SQDIFF) then
                      DELTA = DMIN1(1.0_dp,DELTA*2.0_dp)
                   endif
                endif
                GO TO 5
             endif
             !**        Calculate new trust region size from an estimate of the
             !**        minimum along the step direction using a cubic approximation.
             TEMP = -3.0_dp*SQDPRED !<0
             DELTA =  W*V1*(D2SQVW2-sqrt(D2SQVW2**2+TEMP*DSQVW))/TEMP !>0
          endif !DSQV**2.LT.DSQXI2*V2*TOL2
       enddo !it2
       
5      SQ = SQLIN
       XI(1) = XILIN(1)
       XI(2) = XILIN(2)
       IT = IT+1
    enddo
    if(.not.converged) then
       WRITE(*,'('' >>WARNING!!! Projection iterations have not converged'')')
       WRITE(*,'(14X,''Estimate of error magnitude in xi:'',D9.2,''.'')') W*V1
    endif
    if(.not.inelem.and.XI(1).ge.0.0_dp.and.XI(1).LE.1.0_dp.and. &
         XI(2).ge.0.0_dp.and.XI(2).LE.1.0_dp) then
       inelem = .TRUE.
    endif

9998 RETURN
  end subroutine project_closest

!!!##############################################################################################

  subroutine exnode_deform(num_list,nplist,FILE,NODE_NAME)

    !     Parameter List
    integer :: nplist(:),num_list
    character :: FILE*100,NODE_NAME*50
    !     Local Variables
    integer :: ifile = 10,nc = 1,nh,NOLIST,np
    character :: readfile*100

    readfile = trim(file)//'.exnode'
    open(ifile, file = readfile, status = 'replace')
    !**   write the group name
    write(ifile,'( '' Group name: '',A)') trim(NODE_NAME)
    write(ifile,'('' #Fields=2'')')
    write(ifile,'('' 1) deformed, coordinate,'',' &
         //''' rectangular cartesian, #Components=3'')')
    write(ifile,'(''   x. Value index=1, #Derivatives= 0'')')
    write(ifile,'(''   y. Value index=2, #Derivatives= 0'')')
    write(ifile,'(''   z. Value index=3, #Derivatives= 0'')')
    write(ifile,'('' 2) deformed_fibres, anatomical, fibre, #Components=1'')')
    write(ifile,'(''   fibre angle.  Value index=4, #Derivatives= 0'')')

    do NOLIST = 1,num_list !NPLIST(0)
       NP = NPLIST(NOLIST)
       WRITE(IFILE,'(1X,''Node: '',I12)') NP
       do nh = 1,3
          WRITE(IFILE,'(2X,1(1X,E24.16))') zp(nh,np,nc) !ZP(1,1,nh,NP,nc)
       enddo
       nh = 4                 ! fibre
       WRITE(IFILE,'(2X,(1X,E24.16))') 0.0_dp
    enddo                     !nolist (np)

    close(ifile)
      
  end subroutine exnode_deform

!!!##############################################################################################

  subroutine output_mechanics_results(outfile)

    character(len = *) :: outfile
    !     Local Variables
    integer :: ncount,ne,ng,nj,nm
    real(dp) :: AZ,AZL(3,3),AZU(3,3),RG2D,RGZ,RGZ2D,TC(3,3),TG(3,3),TN(3,3),xi(3)
    real(dp) :: stress,xe(nsm,20),xg(20,num),ze(nsm,nhm),ze_nj(nsm),zg(nhm,num),z(3)
    character(len = 50) :: writefile

    writefile = trim(outfile)//'.test'
    open(10, file = writefile, status = 'replace')
    writefile = trim(outfile)//'.exnode'
    open(20, file = writefile, status = 'replace')
    
    !**   write the group name
    write(20,'( '' Group name: '',A)') 'gauss_points'
    write(20,'('' #Fields=3'')')
    write(20,'('' 1) coordinates, coordinate,'',' &
         //''' rectangular cartesian, #Components=3'')')
    write(20,'(''   x. Value index=1, #Derivatives= 0'')')
    write(20,'(''   y. Value index=2, #Derivatives= 0'')')
    write(20,'(''   z. Value index=3, #Derivatives= 0'')')
    write(20,'('' 2) expansion_ratio, field,'',' &
         //''' rectangular cartesian, #Components=1'')')
    write(20,'(''   1. Value index=4, #Derivatives= 0'')')
    write(20,'('' 3) hydro_stress, field,'',' &
         //''' rectangular cartesian, #Components=1'')')
    write(20,'(''   1. Value index=5, #Derivatives= 0'')')

    ncount = 0
    
    do ne = 1,tissue_num_elems
       call xpxe(nb_lung,ne,xe)
       call zpze(nb_lung,ne,ze)
       forall(nm = 1:4) cg(nm,:) = ce(nm,ne) 
       do ng = 1,NGT(nb_lung)
          xi(1:3) = xig(1:3,ng)
          do nj = 1,3
             ze_nj(:) = ze(:,nj)
             Z(nj) =  pxi(nb_lung,1,xi,ze_nj)
          enddo
          call ZETX50(nb_lung,ng,RG2D,RGZ,RGZ2D,TC,TG,TN,xe,xg,ze,zg) ! only call
          call ZGMG(nb_lung,AZ,AZL,AZU,zg)
          YG(1,ng,ne) = sqrt(DET(AZL)) !ratio deformed to undeformed
          ! Use Cauchy stress tensor to get pressure w.r.t. deformed geometry
          stress = (TC(1,1)+TC(2,2)+TC(3,3))/3.0_dp
          write(10,'(5(f10.3))') z(1:3),yg(1,ng,ne), stress / 98.0665_dp
          !           YG(11,ng,ne)=(TC(1,1)+TC(2,2)+TC(3,3))/3._dp
          ncount = ncount + 1
          write(20,'(1X,''Node: '',I12)') ncount
          write(20,'(2X,5(1X,E24.16))') z(1:3), yg(1,ng,ne), stress
       enddo                  !ng
    enddo                     !noelem
    close(10)    
    close(20)    

  end subroutine output_mechanics_results
  
!!!##############################################################################################
!!!                                      functions
!!!##############################################################################################

  function det(a)

    ! determinant of 3*3 matrix a.

    !     Parameter List
    real(dp) :: a(:,:),det

    DET = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) &
         +A(1,2)*(A(2,3)*A(3,1)-A(3,3)*A(2,1)) &
         +A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

  end function det
  
!!!##############################################################################################
 
 function norm3( N, X, INCX )
   ! returns the euclidean norm of a vector with steps (incx) > 1

   !     Parameter List
   integer :: N,INCX
   real(dp) :: X(*)
   !     Local Variables
   integer :: i,nhigh
   real(dp) :: sum,norm3

   nhigh = 1+( n-1 )*incx
   sum = 0.0_dp
   do i = 1,nhigh,incx
      sum = sum + x(i)**2
   enddo
   norm3 = sqrt(sum)
   
   return
 end function norm3

!!!##############################################################################################

  function getnyr(nc,nrc,nrc1,ny)

    !###    getnyr returns the corresponding variable nyr with a nrc and
    !###    nc for a region nr for a variable ny from nrc1, eg. the
    !###    equivalent flux variable.

    !     Parameter List
    integer :: nc,nrc,nrc1,ny
    !     Local Variables
    integer :: nh,np
    integer :: getnyr

    nh = npny(3,ny,nrc1)
    np = npny(4,ny,nrc1)
    getnyr = nynp(nh,np,nrc,nc)

  end function getnyr

!!!##############################################################################################

  function local_tissue_node(node)
    integer,intent(in) :: node
    integer :: local_tissue_node,np
    logical :: found

    found = .false.
    np = 1
    
    do while(.not.found.and.np.le.tissue_num_nodes)
       if(tissue_nodes(np).eq.node)then
          found = .true.
          local_tissue_node = np
       else
          np = np + 1
       endif
    enddo
    if(.not.found) local_tissue_node = 0
    
  end function local_tissue_node
  
!!!##############################################################################################

  function local_cavity_node(node)
    integer,intent(in) :: node
    integer :: local_cavity_node,np
    logical :: found

    found = .false.
    np = 1
    
    do while(.not.found.and.np.le.cavity_num_nodes)
       if(cavity_nodes(np).eq.node)then
          found = .true.
          local_cavity_node = np
       else
          np = np + 1
       endif
    enddo
    if(.not.found) local_cavity_node = 0
    
  end function local_cavity_node
  
!!!##############################################################################################

  function pl2(i,k,xi)

    ! Evaluate 1D quadratic Lagrange basis function at xi.
    integer,intent(in) :: i,k
    real(dp),intent(in) :: xi
    integer :: i_k
    real(dp) :: pl2

    i_k = 10*k + i

    select case(i_k)
    case(11)
       PL2 = 1.0_dp-3.0_dp*XI+2.0_dp*XI*XI
    case(12)
       PL2 = 4.0_dp*XI*(1.0_dp-XI)
    case(13)
       PL2 = XI*(XI+XI-1.0_dp)
    case(21)
       PL2 = 4.0_dp*XI-3.0_dp
    case(22)
       PL2 = 4.0_dp-8.0_dp*XI
    case(23)
       PL2 = 4.0_dp*XI-1.0_dp
    case(31)
       PL2 = 4.0_dp
    case(32)
       PL2 = -8.0_dp
    case(33)
       PL2 = 4.0_dp
    end select
    
  end function pl2

!!!##############################################################################################

  function psi1(nb,nu,nn,XI)

    ! psi1 evaluates tensor product Lagrange basis functions at XI.

    !**** IPU(nu,ni),nu=1,NUT(nb) identifies the complete set of partial
    !**** derivatives with respect to Xi(ni).

    !     Parameter List
    integer :: nb,nn,nu
    real(dp) :: XI(3)
    !     Local Variables
    integer :: i,IPU(11,3),k,ni
    real(dp) :: psi1

    DATA IPU/1,2,3,1,1,2,1,1,2,1,2, &
             1,1,1,2,3,2,1,1,1,2,2, &
             1,1,1,1,1,1,2,3,2,2,2/ 

    psi1 = 1.0_dp
    do ni = 1,nit(nb) ! 2 or 3
       i = inp(nn,ni,nb)
       k = ipu(nu,ni)
       psi1 = psi1*PL2(i,k,XI(ni))
    enddo

  end function psi1

!!!##############################################################################################

 function pxi(nb,nu,xi,xe)
   
   !##  pxi interpolates nodal array XE at XI
   !     Parameter List
   integer :: nb,nu
   real(dp) :: XE(NSM),XI(:)
   !     Local Variables
   integer :: nn,ns
   real(dp) :: pxi

   pxi = 0.0_dp
   ns = 0
   do nn = 1,nnt(nb)
      ns = ns+1
      pxi = pxi+psi1(nb,nu,nn,xi)*XE(ns)
   enddo
   
  end function pxi

!!!##############################################################################################

end module lung_mechanics


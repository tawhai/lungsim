    n_dr = 0
    diameter_ratio = 0.0_dp
    length_diameter = 0.0_dp
    nums = 0
    radii = 0.0_dp
    lengths = 0.0_dp
    nmax_gen = 0

    do ne = 1,num_elems
       ngen = elem_ordrs(1,ne) ! generation
       nhrd = elem_ordrs(2,ne) ! Horsfield order
       nstr = elem_ordrs(3,ne) ! Strahler order
       ! record by generation
       nums(ngen,1) = nums(ngen,1) + 1 ! count branches in generation
       radii(ngen,1) = radii(ngen,1) + elem_field(ne_radius,ne) ! sum the radii by gen
       lengths(ngen,1) = lengths(ngen,1) + elem_field(ne_length,ne) ! sum the lengths by gen
       ! record by Horsfield order
       nums(nhrd,2) = nums(nhrd,2) + 1
       radii(nhrd,2) = radii(nhrd,2) + elem_field(ne_radius,ne)
       lengths(nhrd,2) = lengths(nhrd,2) + elem_field(ne_length,ne)
       ! record by Strahler order
       nums(nstr,3) = nums(nstr,3) + 1
       radii(nstr,3) = radii(nstr,3) + elem_field(ne_radius,ne)
       lengths(nstr,3) = lengths(nstr,3) + elem_field(ne_length,ne)
       
       if (ngen.gt.nmax_gen) nmax_gen = ngen
       length_diameter = length_diameter + elem_field(ne_length,ne)/ &
            (2.0_dp*elem_field(ne_radius,ne))
       if(elem_cnct(-1,0,ne).ne.0)then
          ne0 = elem_cnct(-1,1,ne)
          if(elem_ordrs(1,ne0).ne.ngen)then
             diameter_ratio = diameter_ratio + elem_field(ne_radius,ne)/elem_field(ne_radius,ne0)
             n_dr = n_dr + 1
          endif
       endif
    enddo
    
    length_diameter = length_diameter/real(num_elems)
    diameter_ratio = diameter_ratio/real(n_dr)
    
    do ngen = 1,nmax_gen
       radii(ngen,1) = radii(ngen,1)/real(nums(ngen,1)) ! generations
       radii(ngen,2) = radii(ngen,2)/real(nums(ngen,2)) ! orders
       lengths(ngen,1) = lengths(ngen,1)/real(nums(ngen,1)) ! generations
       lengths(ngen,2) = lengths(ngen,2)/real(nums(ngen,2)) ! orders
    enddo
    write(*,'('' Model tree geometry:'')')
    write(*,'(''--gen#   #brns   avLen   avRad  avL/D--'')')
    do ngen = 1,nmax_gen
       write(*,'(i6, i8, 3(f10.2))') ngen,nums(ngen,1),lengths(ngen,1), &
            radii(ngen,1),lengths(ngen,1)/(2.0_dp*radii(ngen,1))
    enddo
    write(*,'(''--ord#   #brns   avLen   avRad  avL/D--'')')
    do nhrd = nmax_gen,1,-1
       write(*,'(i6, i8, 3(f10.2))') nhrd,nums(nhrd,2),lengths(nhrd,2), &
            radii(nhrd,2),lengths(nhrd,2)/(2.0_dp*radii(nhrd,2))
    enddo
    
    write(*,'('' L/D ='',f10.2,''; D/Dp ='',f10.2)') length_diameter, &
         diameter_ratio
    
!!!#############################################################################
 
  subroutine list_tree_statistics(filename)

    character(len=*),intent(in) :: filename
    !     Local Variables
    integer :: ne,ne0,ngen,nhrd,nmax_gen,nums(100,2),n_dr
    real(dp) :: diameter_ratio,radii(100,2),length_diameter,lengths(100,2)
    real(dp),allocatable :: diameters(:)
    character(len=300) :: treefile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'list_tree_statistics'
    call enter_exit(sub_name,1)

    if(len(filename) > 0)then
       if(index(filename, ".tree")> 0) then !full filename is given
          treefile = filename
       else ! need to append the correct filename extension
          treefile = trim(filename)//'.tree'
       endif
       open(10, file=treefile, status='replace')
    endif

    allocate(diameters(num_elems))
    diameters = 0.0_dp
    do ne = 1,num_elems
       diameters(ne) = elem_field(ne_radius,ne) * 2.0_dp
    enddo

!!! Initialise arrays
    means = 0.0_dp
    sdt = 0.0_dp
    ntotaln = 0
    sum_mean = 0.0_dp
    sd = 0.0_dp
    ntally = 0
    n_terminal = 0
    nbranches = 0
    branches = 0.0_dp
    nbins = 0
    bins = 0.0_dp

    STATS(2,ne)=undefined
    STATS(5,ne)=undefined
    STATS(6,ne)=undefined

!!! Initialise counters
    ntotal = 0
    num_ddp = 0
    num_llp = 0
    N = 0
        
    do ne = 1,num_elems
       ne0 = elem_cnct(-1,1,ne) ! parent element number
       ind(:) = elem_ords(:,ne) ! the gen, Hord, Sord for i=1,2,3

       ! to add to stats or not? only count as 'extra' branch if the element is at the start of a branch
       add = .false.
       if(ne0.eq.0)then ! this is a stem branch so add
          add = .true.
       else if(ne0.ne.0.and.elem_ordrs(1,ne0).ne.ind(1))then ! gen not same as parent so add
          add = .true.
       endif
       
       if(add)then
          N = N + 1
          nbranches(:,N) = ind(:) ! generation, H order, S order, D-D S order
          if(ne0.ne.0) nbranches(5,N) = elem_ordrs(3,ne0) ! Strahler order of parent

          ! Add length of all segments along branch, calculate their mean diameter
          n_segments=1
          mean_diam = diameters(ne)
          branches(1,N) = elem_field(ne_length,ne) 
          ne_next = ne
          do while(elem_cncts(1,0,ne_next).eq.1) ! while a line of elements
             ne_next = elem_cncts(1,1,ne_next) !next element
             branches(1,N) = branches(1,N) + elem_field(ne_length,ne_next) ! sum lengths
             mean_diam = mean_diam + diameters(ne_next) ! sum diameters
             n_segments=n_segments+1 !c ount number of segments in branch
          enddo
          stats(5,ne) = branches(1,N) ! the branch length
          stats(6,ne) = mean_diam/dble(n_segments) ! mean branch diameter
          branches(2,N) = mean_diam/dble(n_segments) ! record mean diameter
          branches(5,N) = branches(1,N)/branches(2,N) ! L/D

          ! Calculate branching angle to parent
          if(ind(1).gt.1)then ! only calculate angles for elements higher than stem element
             np0 = elem_nodes(1,ne0) ! start of parent
             np1 = elem_nodes(1,ne)  ! start node
             np2 = elem_nodes(2,ne)  ! end node
             xp0(:) = node_xyz(:,np0)
             xp1(:) = node_xyz(:,np1)
             xp2(:) = node_xyz(:,np2)
             angle = angle_btwn_points(xp0, xp1, xp2)
             stats(2,ne) = angle !temporary storage of angle
             branches(3,N) = angle*180.0_dp/pi !store the branching angle to parent
             ntotal = ntotal+1
             if(diameters(ne0).gt.0.0_dp.and.diameters(ne).gt.0.0_dp)then
                if(diameters(ne)/diameters(ne0).le.1.0_dp) num_ddp = num_ddp+1
             endif
             if(diameters(ne0).ge.4.0_dp)then
                nbins(1) = nbins(1) + 1
                bins(1) = bins(1) + angle
             else if(diameters(ne0).ge.2.0_dp)then
                nbins(2) = nbins(2) + 1
                bins(2) = bins(2) + angle
             else if(diameters(ne0).ge.1.0_dp)then
                nbins(3) = nbins(3) + 1
                bins(3) = bins(3) + angle
             else if(diameters(ne0).ge.0.7_dp)then
                nbins(4) = nbins(4) + 1
                bins(4) = bins(4) + angle
             endif
          else
             branches(3,N) = undefined ! no angle calculated
          endif
       endif ! end of add condition

       ! count the terminal branches in each generation
       if(elem_cncts(1,0,ne).eq.0)then ! this is a terminal element
          n_terminal(ind(1)) = n_terminal(ind(1)) + 1 ! ind(1) is element generation
       endif
       
       ! Calculate geometric properties of tree
       branches(4,N) = undefined !initialise to no rotation angle
       if(elem_cncts(-1,0,ne).gt.0.and.elem_cncts(1,0,ne).gt.1)then
          if(elem_cncts(1,0,ne0).gt.1)then
             ne1 = elem_cnct(1,1,ne0) ! first child of parent
             if(ne1.eq.ne) ne1 = elem_cnct(1,2,ne0) ! sibling element number
             np1 = elem_nodes(1,ne)  ! start node
             np2 = elem_nodes(2,ne)  ! end node
             np3 = elem_nodes(2,ne1) ! end node of sibling
             np4 = elem_nodes(2,elem_cncts(1,1,ne)) ! end node of fist child
             np5 = elem_nodes(2,elem_cncts(1,2,ne)) ! end node of second child
             xp1(:) = node_xyz(:,np1)
             xp2(:) = node_xyz(:,np2)
             xp3(:) = node_xyz(:,np3)
             xp4(:) = node_xyz(:,np4)
             xp5(:) = node_xyz(:,np5)
             call make_plane_from_3points(norm_1,2,xp1,xp2,xp3)) ! calculate unit normal and plane
             call make_plane_from_3points(norm_2,2,xp2,xp4,xp5)) ! calculate unit normal and plane
             angle = angle_btwn_vectors(norm_1,norm_2)
             branches(4,N) = angle*180.0_dp/pi ! rotation angle between branching planes
          endif
       endif
    enddo ! ne
        
    n_br = N
    do ne = 1,num_elems
       stats(11:21,ne) = undefined
       ne0 = elem_cncts(-1,1,ne) ! parent element
       if(ne0.ne.0)then ! not a stem
          if(elem_ordrs(1,ne0).ne.elem_ordrs(1,ne))then ! not an intermediate branch element
             if(stats(5,ne)/stats(5,ne0).le.1.0_dp) num_llp = num_llp + 1
             stats(19,ne) = stats(5,ne)/stats(5,ne0) ! L/Lparent
             if(diameters(ne0).gt.0.0_dp.and.diameters(ne).gt.0.0_dp)then
                stats(16,ne) = diameters(ne)/diameters(ne0) !D/Dparent
             endif
          endif
       endif
       if(elem_cncts(1,0,ne).ge.2)then !'bi'furcations only
          ne1 = elem_cncts(1,1,ne) !first child
          ne2 = elem_cncts(1,2,ne) !second child
            
!!!   Summary statistics
          if(stats(6,ne1).lt.undefined.and.stats(6,ne2).lt.undefined)then
             if(stats(6,ne1).ge.stats(6,ne2))then !diameter classification
                ne_major = ne1
                ne_minor = ne2
             else
                ne_major = ne2
                ne_minor = ne1
             endif
             if(stats(2,ne_minor).lt.undefined.and.stats(2,ne_major).lt.undefined)then
                stats(11,ne) = stats(2,ne_minor)*180.0_dp/pi
                stats(12,ne) = stats(2,ne_major)*180.0_dp/pi
             endif
              
             if(diameters(ne_minor).gt.0.0_dp.and.diameters(ne_major).gt.0.0_dp)then
                stats(13,ne) = stats(5,ne_minor)/diameters(ne_minor) !L/D minor
                stats(14,ne) = stats(5,ne_major)/diameters(ne_major) !L/D major
                stats(15,ne) = diameters(ne_minor)/diameters(ne_major) !minor D / major D
                stats(17,ne) = diameters(ne_minor)/diameters(ne) !minor D / D parent
                stats(18,ne) = diameters(ne_major)/diameters(ne) !major D / D parent
             endif
             if(stats(5,ne1).le.stats(5,ne2))then !length classification
                ne_major = ne1
                ne_minor = ne2
             else
                ne_major = ne2
                ne_minor = ne1
             endif !length criteria
             stats(20,ne) = stats(5,ne_major)/stats(5,ne_minor)
          endif
       endif ! elem_cncts
    enddo ! ne
        
!!! Calculate mean branching statistics from values in 'branches' (not elements!)
    do N = 1,N_BR
       do i = 1,3 !for generations, Horsfield orders, Strahler orders
          ind(i) = nbranches(i,N) ! the branch gen, Hord, Sord
!!!...... length and diameter            
          do j = 1,2
             sum_mean(i,j,ind(i)) = sum_mean(i,j,ind(i)) + branches(j,N)
             if(i.eq.3.and.j.eq.1)then
                if(ind(i).ne.nbranches(5,N))then !not same as parent 
                   ntally(i,j,ind(i))=ntally(i,j,ind(i))+1
                endif
             else
                ntally(i,j,ind(i))=ntally(i,j,ind(i))+1
             endif
          enddo !j
!!!...... branching angle and rotation angle            
          do j = 3,4
             if(branches(j,N).lt.undefined)then
                sum_mean(i,j,ind(i)) = sum_mean(i,j,ind(i)) + branches(j,N)
                ntally(i,j,ind(i)) = ntally(i,j,ind(i))+1
             endif
          enddo !j
!!!...... ratio of L:D            
          j = 5
          if(branches(j,N).lt.undefined)then
             sum_mean(i,j,ind(i)) = sum_mean(i,j,ind(i)) + branches(j,N)
             ntally(i,j,ind(i)) = ntally(i,j,ind(i))+1
          endif
       enddo !i
          
!!!... Summary statistics from branches
       do j = 3,5 !branching angle, rotation angle, L/D
          if(branches(j,N).lt.undefined)then
             means(j-2)=means(j-2)+branches(j,N)
          endif
       enddo !j
          
    enddo ! N (for all branches)

!!! Summary statistics by generation
    do N = 1,genm
       do j = 3,5
          ntotaln(j-2) = ntotaln(j-2) + ntally(1,j,N)
       enddo !j
    enddo !N
    do N = 1,genm
       do i = 1,3
          do j = 1,5
             if(ntally(i,j,N).gt.0)then
                sum_mean(i,j,N) = sum_mean(i,j,N)/dble(ntally(i,j,N))
                nmax_gen(i) = N
             else
                sum_mean(i,j,N) = 0.0_dp
             endif
          enddo ! j
       enddo ! i
    enddo ! N
    do N = 1,5
       if(nbins(N).ne.0)then
          bins(N) = bins(N)/dble(nbins(N))*180.0_dp/PI
       endif
    enddo ! N
        
!!!...... Summary statistics from branches
    do j = 3,5 !branching angle, rotation angle, L/D
       if(ntotaln(j-2).ne.0)then
          means(j-2) = means(j-2)/dble(ntotaln(j-2))
       else
          means(j-2) = 0.0_dp
       endif
    enddo !j
        
    i = 2 !Horsfield orders
    j = 6 !Nw/Nw-1
    do N = 1,genm-1
       if(ntally(i,1,N).gt.0.and.ntally(i,1,N+1).gt.0)then
          sum_mean(i,j,N)=dble(ntally(i,1,N))/dble(ntally(i,1,N+1))
       else
          sum_mean(i,j,N)=0.0_dp
       endif
    enddo !N
        
!!! Summary statistics from CE
    do ne = 1,num_elems
       do j = 11,21
          if(stats(j,ne).lt.undefined)then
             means(j-7) = means(j-7) + stats(j,ne)
             ntotaln(j-7) = ntotaln(j-7) + 1
          endif
       enddo !j
    enddo ! ne
        
    do j = 11,21
       if(ntotaln(j-7).gt.0)then
          means(j-7) = means(j-7)/dble(ntotaln(j-7))
       endif
    enddo !j
!!! End of mean calculations
        
!!! Calculate the standard deviations...... sum of (value-mean)^2
    do N = 1,n_br
       do i = 1,3 !for generations, Horsfield orders, Strahler orders
          ind(i) = nbranches(i,N)
          do j = 1,5 !length, diameter, branching angle, rotation angle, L/D
             if(branches(j,N).lt.undefined)then
                SD(i,j,ind(i)) = SD(i,j,ind(i)) + (branches(j,N)-sum_mean(i,j,ind(i)))**2.0_dp
             endif
          enddo !j
       enddo !i
       do j = 3,5 !branching angle, rotation angle, L/D
          if(branches(j,N).lt.undefined)then
             SDT(j-2) = SDT(j-2) + (branches(j,N)-means(j-2))**2.0_dp
          endif
       enddo !j
    enddo !N
    do ne = 1,num_elems
       do j = 11,21
          if(stats(j,ne).lt.undefined)then
             SDT(j-7) = SDT(j-7) + (stats(j,ne)-means(j-7))**2.0_dp
          endif
       enddo !j
    enddo !noelem
        
!!! SD = sqrt(1/(n-1)*sum)
    SD = 0.0_dp
    do N = 1,genm
       do i = 1,3 !for generations, Horsfield orders, Strahler orders
          do j = 1,5 !length, diameter, branching angle, rotation angle, L/D
             if(ntally(i,j,N).gt.1)then
                SD(i,j,N)=DSQRT(SD(i,j,N)/dble(ntally(i,j,N)-1))
             endif
          enddo !j
       enddo !i
    enddo !N
    do j = 1,13
       if(ntotaln(j).gt.1)then
          SDT(j) = sqrt((SDT(j)/dble(ntotaln(j)-1))
       endif
    enddo !j
!!! End of standard deviation calculation        
        
!!! Output tree statistics
    average_term_gen = 0.0_dp
    sum_term = 0
    write(*,'(/'' Generation  #branches  #terminal   Length'',10x,''Diameter&
         &        Branching        Rotation         ratio L:D'')')
    write(*,'(24x,''branches     (mm)'',13x,''(mm)'',11x,''angle(deg)&
         &      angle(deg)'')')
    write(*,'(115(''-''))')
        
    i = 1
    do N = 1,nmax_gen(i)
       write(*,'(3(i10),5(f8.2,'' ('',f6.2,'')''))') N,ntally(i,1,N),n_terminal(N), &
            sum_mean(i,1,N),SD(i,1,N),sum_mean(i,2,N),SD(i,2,N),sum_mean(i,3,N),SD(i,3,N), &
            sum_mean(i,4,N),SD(i,4,N),sum_mean(i,5,N),SD(i,5,N)
       average_term_gen = average_term_gen + n_terminal(N) * N
       sum_term = sum_term + n_terminal(N)
    enddo
    if(sum_term.gt.0)then
       average_term_gen = average_term_gen/dble(sum_term)
    else
       average_term_gen = 0.0_dp
    endif
        
    write(*,'(/'' Horsfield   #branches     Length'',11x,''Diameter&
         &       Branching        Rotation         ratio L:D      Nw/Nw-1'')')
    write(*,'(4x,''order'',20x,''(mm)'',14x,''(mm)'',9x,''angle(deg)&
         &     angle(deg)'')')
    write(*,'(115(''-''))')
        
    i = 2
    do N = 1,nmax_gen(i)
       write(*,'(2(i10),5(f8.2,'' ('',f6.2,'')''),f8.2)') N,ntally(2,1,N),sum_mean(i,1,N), &
            SD(i,1,N),sum_mean(i,2,N),SD(i,2,N),sum_mean(i,3,N),SD(i,3,N),sum_mean(i,4,N), &
            SD(i,4,N),sum_mean(i,5,N),SD(i,5,N),sum_mean(i,6,N)
    enddo
        
    write(*,'(/''   Strahler  #branches    Length'',10x,''Diameter'',8x,''Branching&
         &        Rotation'',10x,''ratio L:D'')')
    write(*,'(''    order'',18x,''(mm)'',13x,''(mm)'',10x,''angle(deg)      angle(deg)'')')
    write(*,'(115(''-''))')
    i = 3
    do N = 1,nmax_gen(i)
       write(*,'(2(i10),5(f8.2,'' ('',f6.2,'')''))') N,ntally(3,1,N),sum_mean(i,1,N),SD(i,1,N), &
            sum_mean(i,2,N),SD(i,2,N),sum_mean(i,3,N),SD(i,3,N),sum_mean(i,4,N),SD(i,4,N), &
            sum_mean(i,5,N),SD(i,5,N)
    enddo
        
    do i = 2,3 !Horsfield and Strahler orders
       do N = 1,nmax_gen(i)
          X(N) = N
          yregress(N,1) = dlog10(dble(ntally(i,1,N)))
          yregress(N,2) = dlog10(sum_mean(i,1,N))
          yregress(N,3) = dlog10(sum_mean(i,2,N))
       enddo !N
       do j = 1,3 !number of branches, length, diameter
          CALL LINREGRESS(nmax_gen(i),R_SQ(i,j),slope,X,yregress(1,j))
          RATIOS(i,j) = 10.0_dp**abs(slope)
       enddo !j
    enddo !i
    write(*,'(/''SUMMARY OF MEAN GEOMETRY STATISTICS'')')
    write(*,'(60(''-''))')
    write(*,'('' terminal generation  = '',f7.3, &
         &'' branching angle      = '',f7.3,'' ('',f6.3,'')'', &
         &'' rotation angle       = '',f7.3,'' ('',f6.3,'')'', &
         &'' minor angle          = '',f7.3,'' ('',f6.3,'')'', &
         &'' major angle          = '',f7.3,'' ('',f6.3,'')'', &
         &'' L/D                  = '',f7.3,'' ('',f6.3,'')'', &
         &'' L/D minor child      = '',f7.3,'' ('',f6.3,'')'', &
         &'' L/D major child      = '',f7.3,'' ('',f6.3,'')'', &
         &'' minor D/major D      = '',f7.3,'' ('',f6.3,'')'', &
         &'' D/Dparent            = '',f7.3,'' ('',f6.3,'')'', &
         &'' %D/Dparent  < 1      = '',f7.3, &
         &'' Dmin/Dparent         = '',f7.3,'' ('',f6.3,'')'', &
         &'' Dmaj/Dparent         = '',f7.3,'' ('',f6.3,'')'', &
         &'' L/Lp                 = '',f7.3,'' ('',f6.3,'')'', &
         &'' %L/Lp < 1            = '',f7.3, &
         &'' L1/L2 (L1 < L2)      = '',f7.3,'' ('',f6.3,'')'')') &
        
         average_term_gen,means(1),SDT(1),means(2),SDT(2),means(4),SDT(4),means(5), &
         SDT(5),means(3),SDT(3),means(6),SDT(6),means(7),SDT(7),means(8),SDT(8), &
         means(9),SDT(9),dble(num_ddp)/dble(ntotal)*100.0_dp,means(10),SDT(10), &
         means(11),SDT(11),means(12),SDT(12),dble(num_llp)/dble(num_elems-1)*100.0_dp, &
         means(13),SDT(13)
        
    write(*,'('' Rb Strahler          = '',f7.3, &
         &'' Rsq ='',f6.3, &
         &'' Rl Strahler          = '',f7.3,'' Rsq = '',f6.3, &
         &'' Rd Strahler          = '',f7.3,'' Rsq = '',f6.3, &
         &'' Rb Horsfield         = '',f7.3,'' Rsq = '',f6.3, &
         &'' Rl Horsfield         = '',f7.3,'' Rsq = '',f6.3, &
         &'' Rd Horsfield         = '',f7.3,'' Rsq = '',f6.3)') &
         RATIOS(3,1),R_SQ(3,1),RATIOS(3,2),R_SQ(3,2),RATIOS(3,3),R_SQ(3,3),RATIOS(2,1), &
         R_SQ(2,1),RATIOS(2,2),R_SQ(2,2),RATIOS(2,3),R_SQ(2,3)
    
    write(*,'('' mean angle Dp 4.0+   = '',f7.3, &
         &''  mean angle Dp 3.0+   = '',f7.3, &
         &''  mean angle Dp 2.0+   = '',f7.3, &
         &''  mean angle Dp 1.0+   = '',f7.3, &
         &''  mean angle Dp 0.7+   = '',f7.3)') (bins(j),j=1,5)

    deallocate(diameters)

  end subroutine list_tree_statistics

SUBROUTINE ibm_type(x,y,z)

! REQD 
! xnp,ynp,znp - Normals (global variables)

! GLOBAL VARS
! type_ibm (not defined)
! xbg,ybg,zbg,xnp,ynp,znp (body pts and normals to be read from a file)

! This subroutine assigns type_ibm(NI,NJ,NK,nblocks) either 0 (body) or 1 (fluid) 
! which can be used with Xgrid, Ygrid, Zgrid to classify points in the domain
	
	use declare_variables
	use kdtree2_module
	
	implicit none
	
	integer :: near_node_idx,near_node_dist,nbrhd_pts_ibm=1
	integer :: node
	
	real,dimension(nodes) :: x,y,z
	real xp,yp,zp,xbp,ybp,zbp,xbnp,ybnp,zbnp,dot_p
	real,allocatable :: global_ibm(:,:),qu_ibm(:)
	
	type(kdtree2_result),allocatable :: results_ibm(:)
	type(kdtree2),pointer :: tree_ibm
	
	allocate(global_ibm(n_dim,nodes))
	allocate(qu_ibm(n_dim))
	allocate(results_ibm(nbrhd_pts_ibm))

	
	do node = 1,nodes
		global_ibm(1,node) = x(node)
		global_ibm(2,node) = y(node)
		global_ibm(3,node) = z(node) ! comment this for n_dim = 2
	enddo
	
	tree_ibm => kdtree2_create(global_ibm,rearrange=.true.,sort=.true.)
	
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
		qu_ibm(1) = xgrid(i,j,k,nbl)
		qu_ibm(2) = ygrid(i,j,k,nbl)
		qu_ibm(3) = zgrid(i,j,k,nbl) ! comment this for n_dim = 2

		nbrhd_pts_ibm = 1

		call kdtree2_n_nearest(tree_ibm,qu_ibm,nbrhd_pts_ibm,results_ibm)
		
		near_node_idx = results_ibm(1)%idx ! reuslts_ibm is a struct with %idx and %dis
		near_node_dist = (results_ibm(1)%dis)**0.5
		
		xp = qu_ibm(1)
		yp = qu_ibm(2)
		zp = qu_ibm(3) ! comment this for n_dim = 2
		
		xbp = x(near_node_idx)
		ybp = y(near_node_idx)
		zbp = z(near_node_idx) ! comment this for n_dim = 2
		
		xbnp = xbn(near_node_idx) 
		ybnp = ybn(near_node_idx)
		zbnp = zbn(near_node_idx) ! comment this for n_dim = 2
		
		dot_p = (xp-xbp)*xbnp + (yp-ybp)*ybnp + (zp-zbp)*zbnp ! the z term is 0 for n_dim = 2
		
		if(dot_p.ge.0) type_ibm(i,j,k,nbl) = 1 ! 1:Fluid Point
		if(dot_p.le.0) type_ibm(i,j,k,nbl) = 0 ! 0:Body Point
	
	enddo
	enddo
	enddo
	enddo

END

!SUBROUTINE ghost_points()
!
!	use declare_variables
!	use kdtree2_module
!	
!	implicit none
!	
!	integer gridpts_iter,nbrhd_pts,nbrhd_pts
!	real,allocatable :: grid_pts(:,:),qv(:)
!	
!	allocatable(grid_pts(3,NImax*NJmax*NKmax*nblocks))
!	allocatable(qv(3))
!	
!	gridpts_iter = 0
!	ghost_pt_iter = 0
!	
!	if(n_dim.eq.3) then
!		nbrhd_pts = 6
!	elseif(n_dim.eq.2) then
!		nbrhd_pts = 6
!	endif
!	
!	do nbl = 1,nblocks
!	do k = NK(nbl)
!	do j = NJ(nbl)
!	do i = NI(nbl)
!	
!		gridpts_iter = gridpts_iter + 1
!		grid_pts(1,gridpts_iter) = Xgrid(i,j,k,nbl)
!		grid_pts(2,gridpts_iter) = Ygrid(i,j,k,nbl)
!		grid_pts(3,gridpts_iter) = Zgrid(i,j,k,nbl)
!	
!	enddo
!	enddo
!	enddo
!	enddo
!	
!	tree_grid_pts => kdtree2_create(grid_pts,rearrange=.true.,sort=.true.)
!	
!	do nbl = 1,nblocks
!	do k = 1,NK(nbl)
!	do j = 1,NJ(nbl)
!	do i = 1,NI(nbl)
!	
!		if(type_ibm(i,j,k,nbl).eq.0) then
!		
!			qv(1) = Xgrid(i,j,k,nbl)
!			qv(2) = Ygrid(i,j,k,nbl)
!			qv(3) = Zgrid(i,j,k,nbl)
!
!			call kdtree2_n_nearest(tree_grid_pts,qv,nbrhd_pts,results)
!		
!			do indx = 1,nbrhd_pts
!			
!				near_node_idx(indx) = results(indx)%idx
!				near_node_dist(indx) = results(indx)%dis
!			
!				call get_loc_index(near_node_idx(indx),nnbbll,ii,jj,kk)
!			
!				if(type_ibm(ii,jj,kk,nnbbll).eq.1) then
!				
!					ghost_pt_iter = ghost_pt_iter + 1
!					ghost_pt_idx(ghost_pt_iter) = near_node_idx(indx)
!					ghost_pt(1,ghost_pt_iter) = Xgrid(ii,jj,kk,nnbbll)
!					ghost_pt(2,ghost_pt_iter) = Ygrid(ii,jj,kk,nnbbll)
!					ghost_pt(3,ghost_pt_iter) = Zgrid(ii,jj,kk,nnbbll)
!					type_ibm(i,j,k,nbl) = -1 ! -1:Ghost Point
!					break
!				
!				endif
!			
!			
!			enddo
!		
!		endif
!	
!	enddo
!	enddo
!	enddo
!	enddo
!
!
!
!END

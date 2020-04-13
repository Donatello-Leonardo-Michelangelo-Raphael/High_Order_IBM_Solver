!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ IBM SOLVER ROUTINES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


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
	
	integer :: near_node_idx,near_node_dist,nbrhd_pts_ibm
	
	real,dimension(nodes) :: x,y,z
	real xp,yp,zp,xbp,ybp,zbp,xbnp,ybnp,zbnp,dot_p
	real,allocatable :: global_ibm(:,:),qu_ibm(:)
	
	type(kdtree2_result),allocatable :: results_ibm(:)
	type(kdtree2),pointer :: tree_ibm

	n_dim = 3
	
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



SUBROUTINE ghost_point()

! Global Vars
! ghost_pt_iter,

	use declare_variables
	use kdtree2_module
	
	implicit none
	
	integer n_dim,gridpts_iter,nbrhd_pts
	
	n_dim = 3
	gridpts_iter = 0
	ghost_pt_iter = 0
	
	if(n_dim.eq.3) then
		nbrhd_pts = 6
	elseif(n_dim.eq.2) then
		nbrhd_pts = 4
	endif
	
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
		gridpts_iter = gridpts_iter + 1
		grid_pts(1,gridpts_iter) = Xgrid(i,j,k,nbl)
		grid_pts(2,gridpts_iter) = Ygrid(i,j,k,nbl)
		grid_pts(3,gridpts_iter) = Zgrid(i,j,k,nbl)
	
	enddo
	enddo
	enddo
	enddo
	
	tree_grid_pts => kdtree2_create(grid_pts,rearrange=.true.,sort=.true.)
	
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
		if(type_ibm(i,j,k,nbl).eq.0) then
			
			qv(1) = Xgrid(i,j,k,nbl)
			qv(2) = Ygrid(i,j,k,nbl)
			qv(3) = Zgrid(i,j,k,nbl)
			
			call kdtree2_n_nearest(tree_grid_pts,qv,nbrhd_pts,results)
			
			do indx = 1,nbrhd_pts
			
				near_node_idx(indx) = results(indx)%idx
				near_node_dist(indx) = results(indx)%dis
				
				call get_loc_index(near_node_idx(indx),nnbbll,ii,jj,kk)
				
				if(type_ibm(ii,jj,kk,nnbbll).eq.1):
					
					ghost_pt_iter = ghost_pt_iter + 1
					ghost_pt_idx(ghost_pt_iter) = near_node_idx(indx)
					ghost_pt(1,ghost_pt_iter) = Xgrid(ii,jj,kk,nnbbll)
					ghost_pt(2,ghost_pt_iter) = Ygrid(ii,jj,kk,nnbbll)
					ghost_pt(3,ghost_pt_iter) = Zgrid(ii,jj,kk,nnbbll)
					type_ibm(i,j,k,nbl) = -1 ! -1:Ghost Point
					break
					
				endif
			
			enddo
			
		endif

	enddo
	enddo
	enddo
	enddo

END
	

SUBROUTINE boundary_intercept(x,y,z)

! x,y,z == xbg,ybg,zbg

! Global Vars
! pts
	
	use declare_variables
	use kdtree2_module
	
	implicit none
	
	integer nbrhd_pts 
	
	nbrhd_pts = 3
	
	do node = 1,nodes
		global_ibm(1,node) = x(node)
		global_ibm(2,node) = y(node)
		global_ibm(3,node) = z(node) ! comment this for n_dim = 2
	enddo
	
	tree => kdtree2_create(global_ibm,rearrange=.true.,sort=.true.)
	
	pts = 0
	
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
		if(type_ibm(i,j,k,nbl).eq.-1) then
			
			qv(1) = Xgrid(i,j,k,nbl)
			qv(2) = Ygrid(i,j,k,nbl)
			qv(3) = Zgrid(i,j,k,nbl)
			
			call kdtree2_n_nearest(tree,qv,nbrhd_pts,results)
			
		do indexing = 1,neighpts_ibm
			near_node = results(indexing)%idx
			dw = (results(indexing)%dis)**0.5		
		enddo
  
  	    ! Loop the elements sharing the nearest node
        dist2 = 1e5
   
		do pp = 1,num_share_elems(near_node)
			elem = ind_share_elems(near_node,pp)
			n1 = connect(elem,1)	  
			n2 = connect(elem,2)
			n3 = connect(elem,3)
			
			x1 = xbg(n1)
			y1 = ybg(n1)
			z1 = zbg(n1)
			
			x2 = xbg(n2)
			y2 = ybg(n2)
			z2 = zbg(n2)
		
			x3 = xbg(n3)
			y3 = ybg(n3)
			z3 = zbg(n3)
		
			!Ghost point coordinates
			xg = qu(1)	
			yg = qu(2)
			zg = qu(3)	  
			
			! Plane equation
			a = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1)
			b = -((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))
			c = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)	  
			d = a*x1 + b*y1 + c*z1
  
            tt = -(a*xg+b*yg+c*zg-d)/(a**2 + b**2 + c**2)
  
            ! Intersection point
            xi = xg + a*tt
            yi = yg + b*tt
            zi = zg + c*tt
  
            ! Verify if the point is in the triangle
           
            v1 = xi-x1
            v2 = yi-y1
            v3 = zi-z1
  
            v11 = x2-x1
            v12 = y2-y1
            v13 = z2-z1
  
            v21 = x3-x1
            v22 = y3-y1
            v23 = z3-z1           
  
            vdv1 = v1*v11 + v2*v12 + v3*v13
            vdv2 = v1*v21 + v2*v22 + v3*v23
            v1dv2 = v11*v21 + v12*v22 + v13*v23 
            v1dv1 = v11*v11 + v12*v12 + v13*v13
            v2dv2 = v21*v21 + v22*v22 + v23*v23
  
            alp = (vdv1*v2dv2 - v1dv2*vdv2)/(v2dv2*v1dv1-v1dv2*v1dv2)
            bet = (vdv2*v1dv1 - vdv1*v1dv2)/(v2dv2*v1dv1-v1dv2*v1dv2)
  
            !print*, 'point:', pts
            !print*, a,b,c,d,t
            !print*,xi,yi,zi
            !print*, alp,bet,alp+bet
  
            if(alp.ge.-1e-1.and.bet.ge.-1e-1.and.(alp+bet).le.1.1) then     ! NOTE - Not strictly less than 1 - Given some tolerance
              
			  dist2_new = (xi-xg)**2+(yi-yg)**2+(zi-zg)**2
              
			  if(dist2_new.le.dist2) then
				
				pts = pts + 1
				BII(1,pts) = xi
				BII(2,pts) = yi
				BII(3,pts) = zi
				dist2 = dist2_new
				flag = 1
              
			  endif        
            
			endif
  
		enddo
		
	enddo
	enddo
	enddo
	enddo
	
END


! ghost_pt_iter -> ghost_pt()
! pts -> BII()

SUBROUTINE nn_fluid(nbrhd_pts,x,y,z)

! x,y,z == BI
	
	use declare_variables
	use kdtree2_module
	
	implicit none
	
	integer nbrhd_pts
	
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
		gridpts_iter = gridpts_iter + 1
		grid_pts(1,gridpts_iter) = Xgrid(i,j,k,nbl)
		grid_pts(2,gridpts_iter) = Ygrid(i,j,k,nbl)
		grid_pts(3,gridpts_iter) = Zgrid(i,j,k,nbl)
	
	enddo
	enddo
	enddo
	enddo
	
	tree_grid_pts => kdtree2_create(grid_pts,rearrange=.true.,sort=.true.)
	
	b_incpt(1) = x
	b_incpt(2) = y
	b_incpt(3) = z
	
	call kdtree2_n_nearest(tree_grid_pts,b_incpt,nbrhd_pts,results)
	
	vicinity_pts_iter = 0
	
	do indx = 1,nbrhd_pts
		near_node = results(indx)%idx
		dist(indx) = results(indx)%dis
		call get_loc_index(near_node,nbl,i,j,k)
		if(type_ibm(i,j,k,nbl).eq.1) then
			vicinity_pts_iter = vicinity_pts_iter + 1
			vicinity_pts(1,vicinity_pts_iter) = Xgrid(i,j,k,nbl)
			vicinity_pts(2,vicinity_pts_iter) = Ygrid(i,j,k,nbl)
			vicinity_pts(3,vicinity_pts_iter) = Zgrid(i,j,k,nbl)
		endif
	enddo

END

! vicinity_pts,vicinity_pts_iter is global
! dist is global

SUBROUTINE matrix_calculations()

! W is global

	use declare_variables
	use kdtree2_n_nearest
	
	implicit none
	
	integer nbrhd_pts
	
	nbrhd_pts = 19
	total_comp_pts = nbrhd_pts+1 ! +1 ghost point
	N = 3
	n_dim = 3
		
	do pp = 1,ghost_pt_iter
	
		xbi = BII(1,pp)
		ybi = BII(2,pp)
		zbi = BII(3,pp)
		
		call nn_fluid(nbrhd_pts,xbi,ybi,zbi)
		
		R = max(dist)
		
		xgp = ghost_pt(1,pp)
		ygp = ghost_pt(2,pp)
		zgp = ghost_pt(3,pp)
		
		prime_coord(1,1) = xgp - xbi
		prime_coord(2,1) = ygp - ybi
		prime_coord(3,1) = zgp - zbi
		
		do qq = 2,total_comp_pts
			
			prime_coord(1,qq) = vicinity_pts(1,qq) - xbi
			prime_coord(2,qq) = vicinity_pts(2,qq) - ybi
			prime_coord(3,qq) = vicinity_pts(3,qq) - zbi
			W(qq,qq,pp) = 0.5*(1+cos(3.141*(prime_coord(1,qq)**2+prime_coord(2,qq)**2+prime_coord(3,qq)**2)**0.5/R))
			
		enddo
		
		call ibm_coeff_vandermonde(N,n_dim,pp,prime_coord,total_comp_pts)
		
		
		p_i = matmul(W(:,:,pp),V(:,:,pp))
		ps_inv = func(p_i) ! function to calculate psuedo inverse
		! code taken from here https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=5&t=2712 (not working)
		A(:,:,pp) = matmul(ps_inv,W(:,:,pp))
		
	
	enddo
	
END	


SUBROUTINE ibm_coeff_vandermonde(N,n_dim,ii,prime_coord,total_comp_pts)

! L_N is a global variable
! V is a global variable

	use declare_variables
	
	implicit none
	
	if(N.eq.1.and.n_dim.eq.2) then
		L_N = 3
		do jj = 1,total_comp_pts
			V(jj,1,ii) = 1
			V(jj,2,ii) = prime_coord(1,jj)
			V(jj,3,ii) = prime_coord(2,jj)
		enddo
	elseif(N.eq.1.and.n_dim.eq.3) then
		L_N = 4
		do jj = 1,total_comp_pts
			V(jj,1,ii) = 1
			V(jj,2,ii) = prime_coord(1,jj)
			V(jj,3,ii) = prime_coord(2,jj)
			V(jj,4,ii) = prime_coord(3,jj)
		enddo
	elseif(N.eq.2.and.n_dim.eq.2) then
		L_N = 6
		do jj = 1,total_comp_pts
			V(jj,1,ii) = 1
			V(jj,2,ii) = prime_coord(1,jj)
			V(jj,3,ii) = prime_coord(2,jj)
			V(jj,4,ii) = prime_coord(1,jj)*prime_coord(1,jj)
			V(jj,5,ii) = prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,6,ii) = prime_coord(2,jj)*prime_coord(1,jj)
		enddo
	elseif(N.eq.2.and.n_dim.eq.3) then
		L_N = 10
		do jj = 1,total_comp_pts
			V(jj,1,ii) = 1
			V(jj,2,ii) = prime_coord(1,jj)
			V(jj,3,ii) = prime_coord(2,jj)
			V(jj,4,ii) = prime_coord(3,jj)
			V(jj,5,ii) = prime_coord(1,jj)*prime_coord(1,jj)
			V(jj,6,ii) = prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,7,ii) = prime_coord(3,jj)*prime_coord(3,jj)
			V(jj,8,ii) = prime_coord(1,jj)*prime_coord(2,jj)
			V(jj,9,ii) = prime_coord(2,jj)*prime_coord(3,jj)
			V(jj,10,ii) = prime_coord(3,jj)*prime_coord(1,jj)
		enddo
	elseif(N.eq.3.and.n_dim.eq.2) then
		L_N = 10
		do jj = 1,total_comp_pts
			V(jj,1,ii) = 1
			V(jj,2,ii) = prime_coord(1,jj)
			V(jj,3,ii) = prime_coord(2,jj)
			V(jj,4,ii) = prime_coord(1,jj)*prime_coord(1,jj)
			V(jj,5,ii) = prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,6,ii) = prime_coord(2,jj)*prime_coord(1,jj)
			V(jj,7,ii) = prime_coord(1,jj)*prime_coord(1,jj)*prime_coord(1,jj)
			V(jj,8,ii) = prime_coord(2,jj)*prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,9,ii) = prime_coord(1,jj)*prime_coord(1,jj)*prime_coord(2,jj)
			V(jj,10,ii) = prime_coord(1,jj)*prime_coord(2,jj)*prime_coord(2,jj)
		enddo
	elseif(N.eq.3.and.n_dim.eq.3) then
		L_N = 20
		do jj = 1,total_comp_pts
			V(jj,1,ii) = 1
			V(jj,2,ii) = prime_coord(1,jj)
			V(jj,3,ii) = prime_coord(2,jj)
			V(jj,4,ii) = prime_coord(3,jj)
			V(jj,5,ii) = prime_coord(1,jj)*prime_coord(1,jj)
			V(jj,6,ii) = prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,7,ii) = prime_coord(3,jj)*prime_coord(3,jj)
			V(jj,8,ii) = prime_coord(1,jj)*prime_coord(2,jj)
			V(jj,9,ii) = prime_coord(2,jj)*prime_coord(3,jj)
			V(jj,10,ii) = prime_coord(3,jj)*prime_coord(1,jj)
			V(jj,11,ii) = prime_coord(1,jj)**3
			V(jj,12,ii) = prime_coord(2,jj)**3
			V(jj,13,ii) = prime_coord(3,jj)**3
			V(jj,14,ii) = prime_coord(1,jj)**2*prime_coord(2,jj)
			V(jj,15,ii) = prime_coord(2,jj)**2*prime_coord(3,jj)
			V(jj,16,ii) = prime_coord(3,jj)**2*prime_coord(1,jj)
			V(jj,17,ii) = prime_coord(1,jj)*prime_coord(2,jj)**2
			V(jj,18,ii) = prime_coord(2,jj)*prime_coord(3,jj)**2
			V(jj,19,ii) = prime_coord(3,jj)*prime_coord(1,jj)**2
			V(jj,20,ii) = prime_coord(1,jj)*prime_coord(2,jj)*prime_coord(3,jj)
		enddo
	elseif(N.eq.4.and.n_dim.eq.2) then
		L_N = 15
		do jj = 1,total_comp_pts
			V(jj,1,ii) = 1
			V(jj,2,ii) = prime_coord(1,jj)
			V(jj,3,ii) = prime_coord(2,jj)
			V(jj,4,ii) = prime_coord(1,jj)*prime_coord(1,jj)
			V(jj,5,ii) = prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,6,ii) = prime_coord(2,jj)*prime_coord(1,jj)
			V(jj,7,ii) = prime_coord(1,jj)*prime_coord(1,jj)*prime_coord(1,jj)
			V(jj,8,ii) = prime_coord(2,jj)*prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,9,ii) = prime_coord(1,jj)*prime_coord(1,jj)*prime_coord(2,jj)
			V(jj,10,ii) = prime_coord(1,jj)*prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,11,ii) = prime_coord(1,jj)**4
			V(jj,12,ii) = prime_coord(2,jj)**4
			V(jj,13,ii) = prime_coord(1,jj)**3*prime_coord(2,jj)
			V(jj,14,ii) = prime_coord(1,jj)**2*prime_coord(2,jj)**2
			V(jj,15,ii) = prime_coord(1,jj)*prime_coord(2,jj)**3
		enddo
	elseif(N.eq.4.and.n_dim.eq.3) then
		L_N = 35
		do jj = 1,total_comp_pts
			V(jj,1,ii) = 1
			V(jj,2,ii) = prime_coord(1,jj)
			V(jj,3,ii) = prime_coord(2,jj)
			V(jj,4,ii) = prime_coord(3,jj)
			V(jj,5,ii) = prime_coord(1,jj)*prime_coord(1,jj)
			V(jj,6,ii) = prime_coord(2,jj)*prime_coord(2,jj)
			V(jj,7,ii) = prime_coord(3,jj)*prime_coord(3,jj)
			V(jj,8,ii) = prime_coord(1,jj)*prime_coord(2,jj)
			V(jj,9,ii) = prime_coord(2,jj)*prime_coord(3,jj)
			V(jj,10,ii) = prime_coord(3,jj)*prime_coord(1,jj)
			V(jj,11,ii) = prime_coord(1,jj)**3
			V(jj,12,ii) = prime_coord(2,jj)**3
			V(jj,13,ii) = prime_coord(3,jj)**3
			V(jj,14,ii) = prime_coord(1,jj)**2*prime_coord(2,jj)
			V(jj,15,ii) = prime_coord(2,jj)**2*prime_coord(3,jj)
			V(jj,16,ii) = prime_coord(3,jj)**2*prime_coord(1,jj)
			V(jj,17,ii) = prime_coord(1,jj)*prime_coord(2,jj)**2
			V(jj,18,ii) = prime_coord(2,jj)*prime_coord(3,jj)**2
			V(jj,19,ii) = prime_coord(3,jj)*prime_coord(1,jj)**2
			V(jj,20,ii) = prime_coord(1,jj)*prime_coord(2,jj)*prime_coord(3,jj)
			V(jj,21,ii) = prime_coord(1,jj)**4
			V(jj,22,ii) = prime_coord(2,jj)**4
			V(jj,23,ii) = prime_coord(3,jj)**4
			V(jj,24,ii) = prime_coord(1,jj)**3*prime_coord(2,jj)
			V(jj,25,ii) = prime_coord(2,jj)**3*prime_coord(3,jj)
			V(jj,26,ii) = prime_coord(3,jj)**3*prime_coord(1,jj)
			V(jj,27,ii) = prime_coord(1,jj)*prime_coord(2,jj)**3
			V(jj,28,ii) = prime_coord(2,jj)*prime_coord(3,jj)**3
			V(jj,29,ii) = prime_coord(3,jj)*prime_coord(1,jj)**3
			V(jj,30,ii) = prime_coord(1,jj)**2*prime_coord(2,jj)**2
			V(jj,31,ii) = prime_coord(2,jj)**2*prime_coord(3,jj)**2
			V(jj,32,ii) = prime_coord(3,jj)**2*prime_coord(1,jj)**2
			V(jj,33,ii) = prime_coord(1,jj)**2*prime_coord(2,jj)*prime_coord(3,jj)
			V(jj,34,ii) = prime_coord(2,jj)**2*prime_coord(3,jj)*prime_coord(1,jj)
			V(jj,35,ii) = prime_coord(3,jj)**2*prime_coord(1,jj)*prime_coord(2,jj)
		enddo
	endif
	
END




SUBROUTINE c_coeff(PHI)

	use declare_variables
	
	implicit none
	
	do gg = 1,ghost_pt_iter
	do ii = 1,L_N
	do pp = 1,total_comp_pts
	
		call get_loc_index(ghost_pt_idx(gg),nbl,i,j,k)
		c(ii,gg) = c(ii,gg) + A(ii,pp,gg)*PHI(i,j,k,nbl) 
	
	enddo
	enddo


END



SUBROUTINE phi_gp(PHI_W,PHI)

	use declare_variables
	
	implicit none
	
	do gg = 1,ghost_pt_iter
	call get_loc_index(ghost_pt_idx(gg),nbl,i,j,k)
	do pp = 2,total_comp_pts
	
		call get_loc_index(vicinity_pts())
		PHI(i,j,k,nbl) = A(1,pp,gg)*PHI()
	
	enddo
	enddo

END




SUBROUTINE DISCRETIZATION_I_COMP_IBM(PHI,PHID,nvars)
	use declare_variables
	implicit none
	
	integer prim,nvars,coeff
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	real,dimension(NImax) :: RHS
	
	AC_COMP_IBM(1:ptsmax) = 1.d0
	AP_COMP_IBM(1:ptsmax) = alpha
	AM_COMP_IBM(1:ptsmax) = alpha
	
	do prim = 1,nvars
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
		
		
		do i = 3,NI(nbl)-2
			RHS(i) = (bdisc/4.d0)*(PHI(i+2,j,k,nbl,prim)-PHI(i-2,j,k,nbl,prim)) &
					+ (adisc/2.d0)*(PHI(i+1,j,k,nbl,prim)-PHI(i-1,j,k,nbl,prim))
			if(type_ibm(i,j,k,nbl).eq.0) then
				AM_COMP_IBM(i) = 0
				AP_COMP_IBM(i) = 0
			endif
		Enddo
					
		RHS(1) = (bdisc/4.d0)*(PHI(3,j,k,nbl,prim)-PHI(NI(nbl)-2,j,k,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(2,j,k,nbl,prim)-PHI(NI(nbl)-1,j,k,nbl,prim))
		RHS(2) = (bdisc/4.d0)*(PHI(4,j,k,nbl,prim)-PHI(NI(nbl)-1,j,k,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(3,j,k,nbl,prim)-PHI(1,j,k,nbl,prim))
		RHS(NI(nbl)-2) = (bdisc/4.d0)*(PHI(1,j,k,nbl,prim)-PHI(NI(nbl)-4,j,k,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(NI(nbl)-1,j,k,nbl,prim)-PHI(NI(nbl)-3,j,k,nbl,prim))
		RHS(NI(nbl)-1) = (bdisc/4.d0)*(PHI(2,j,k,nbl,prim)-PHI(NI(nbl)-3,j,k,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(1,j,k,nbl,prim)-PHI(NI(nbl)-2,j,k,nbl,prim))
						
		if(type_ibm(1,j,k,nbl).eq.0) then
			AM_COMP_IBM(1) = 0
			AP_COMP_IBM(1) = 0
		elseif(type_ibm(2,j,k,nbl).eq.0) then
			AM_COMP_IBM(2) = 0
			AP_COMP_IBM(2) = 0
		elseif(type_ibm(NI(nbl)-2,j,k,nbl).eq.0) then
			AM_COMP_IBM(NI(nbl)2) = 0
			AP_COMP_IBM(NI(nbl)-2) = 0
		elseif(type_ibm(NI(nbl)-1,j,k,nbl).eq.0) then
			AM_COMP_IBM(NI(nbl)-1) = 0
			AP_COMP_IBM(NI(nbl)-1) = 0
		endif
					
		call TDMAP(1,NI(nbl)-1,AP_COMP_IBM,AC_COMP_IBM,AM_COMP_IBM,RHS,NI(nbl)-1)
		
		PHID(1:NI(nbl)-1,j,k,nbl,prim) = RHS(1:NI(nbl)-1)
		PHID(NI(nbl),j,k,nbl,prim) = PHID(1,j,k,nbl,prim)
	  
	enddo
	enddo
	enddo
	enddo
	
END


SUBROUTINE DISCRETIZATION_J_COMP_IBM(PHI,PHID,nvars)
	use declare_variables
	implicit none
	
	integer prim,nvars
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	real,dimension(NJmax) :: RHS
	
	AC_COMP_IBM(1:ptsmax) = 1.d0
	AP_COMP_IBM(1:ptsmax) = alpha
	AM_COMP_IBM(1:ptsmax) = alpha
	
	do prim = 1,nvars
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do i = 1,NI(nbl)
	 
	 
		Do j = 3,NJ(nbl)-2
			RHS(j) = (bdisc/4.d0)*(PHI(i,j+2,k,nbl,prim)-PHI(i,j-2,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,j+1,k,nbl,prim)-PHI(i,j-1,k,nbl,prim))
			if(type_ibm(i,j,k,nbl).eq.0) then
				AM_COMP_IBM(j) = 0
				AP_COMP_IBM(j) = 0
			endif
							
		Enddo
		
		RHS(1) = (bdisc/4.d0)*(PHI(i,3,k,nbl,prim)-PHI(i,NJ(nbl)-2,k,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(i,2,k,nbl,prim)-PHI(i,NJ(nbl)-1,k,nbl,prim))
		RHS(2) = (bdisc/4.d0)*(PHI(i,4,k,nbl,prim)-PHI(i,NJ(nbl)-1,k,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(i,3,k,nbl,prim)-PHI(i,1,k,nbl,prim))
		RHS(NJ(nbl)-2) = (bdisc/4.d0)*(PHI(i,1,k,nbl,prim)-PHI(i,NJ(nbl)-4,k,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(i,NJ(nbl)-1,k,nbl,prim)-PHI(i,NJ(nbl)-3,k,nbl,prim))
		RHS(NJ(nbl)-1) = (bdisc/4.d0)*(PHI(i,2,k,nbl,prim)-PHI(i,NJ(nbl)-3,k,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(i,1,k,nbl,prim)-PHI(i,NJ(nbl)-2,k,nbl,prim))
						
		if(type_ibm(i,1,k,nbl).eq.0) then
			AM_COMP_IBM(1) = 0
			AP_COMP_IBM(1) = 0
		elseif(type_ibm(i,2,k,nbl).eq.0) then
			AM_COMP_IBM(2) = 0
			AP_COMP_IBM(2) = 0
		elseif(type_ibm(i,NJ(nbl)-2,k,nbl).eq.0) then
			AM_COMP_IBM(NJ(nbl)2) = 0
			AP_COMP_IBM(NJ(nbl)-2) = 0
		elseif(type_ibm(i,NJ(nbl)-1,k,nbl).eq.0) then
			AM_COMP_IBM(NJ(nbl)-1) = 0
			AP_COMP_IBM(NJ(nbl)-1) = 0
		endif

		call TDMAP(1,NJ(nbl)-1,AP_COMP_IBM,AC_COMP_IBM,AM_COMP_IBM,RHS,NJ(nbl)-1)
		
		PHID(i,1:NJ(nbl)-1,k,nbl,prim) = RHS(1:NJ(nbl)-1)
		PHID(i,NJ(nbl),k,nbl,prim) = PHID(i,1,k,nbl,prim)
	  
	    enddo
		enddo
		enddo
		enddo
	  
END


SUBROUTINE DISCRETIZATION_K_COMP_IBM(PHI,PHID,nvars)
	use declare_variables
	implicit none
	
	integer prim,nvars
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	real,dimension(NKmax) :: RHS
 
	AC_COMP_IBM(1:ptsmax) = 1.d0
	AP_COMP_IBM(1:ptsmax) = alpha
	AM_COMP_IBM(1:ptsmax) = alpha
	
	do prim = 1,nvars
	do nbl = 1,nblocks
	do j = 1,NJ(nbl)
	do i = 1,Ni(nbl)
		
		do k = 3,NK(nbl)-2
			RHS(k) = (bdisc/4.d0)*(PHI(i,j,k+2,nbl,prim)-PHI(i,j,k-2,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,j,k+1,nbl,prim)-PHI(i,j,k-1,nbl,prim))
			if(type_ibm(i,j,k,nbl).eq.0) then
				AM_COMP_IBM(k) = 0
				AP_COMP_IBM(k) = 0
			endif
			
		enddo
		
		RHS(1) = (bdisc/4.d0)*(PHI(i,j,3,nbl,prim)-PHI(i,j,NK(nbl)-2,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(i,j,2,nbl,prim)-PHI(i,j,NK(nbl)-1,nbl,prim))
		RHS(2) = (bdisc/4.d0)*(PHI(i,j,4,nbl,prim)-PHI(i,j,NK(nbl)-1,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(i,j,3,nbl,prim)-PHI(i,j,1,nbl,prim))
		RHS(NK(nbl)-2) = (bdisc/4.d0)*(PHI(i,j,1,nbl,prim)-PHI(i,j,NK(nbl)-4,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(i,j,NK(nbl)-1,nbl,prim)-PHI(i,j,NK(nbl)-3,nbl,prim))
		RHS(NK(nbl)-1) = (bdisc/4.d0)*(PHI(i,j,2,nbl,prim)-PHI(i,j,NK(nbl)-3,nbl,prim)) &
						+ (adisc/2.d0)*(PHI(i,j,1,nbl,prim)-PHI(i,j,NK(nbl)-2,nbl,prim))
			
		if(type_ibm(i,j,1,nbl).eq.0) then
			AM_COMP_IBM(1) = 0
			AP_COMP_IBM(1) = 0
		elseif(type_ibm(i,j,2,nbl).eq.0) then
			AM_COMP_IBM(2) = 0
			AP_COMP_IBM(2) = 0
		elseif(type_ibm(i,j,NK(nbl)-2,nbl).eq.0) then
			AM_COMP_IBM(NK(nbl)2) = 0
			AP_COMP_IBM(NK(nbl)-2) = 0
		elseif(type_ibm(i,j,NK(nbl)-1,nbl).eq.0) then
			AM_COMP_IBM(NK(nbl)-1) = 0
			AP_COMP_IBM(NK(nbl)-1) = 0
		endif
			
		call TDMAP(1,NK(nbl)-1,AP_COMP_IBM,AC_COMP_IBM,AM_COMP_IBM,RHS,NK(nbl)-1)
		
		PHID(i,j,1:NK(nbl)-1,nbl,prim) = RHS(1:NK(nbl)-1)
		PHID(i,j,NK(nbl),nbl,prim) = PHID(i,j,1,nbl,prim)
	 
	enddo
	enddo
	enddo
	enddo
 
END

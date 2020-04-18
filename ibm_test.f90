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
	integer :: flag1,flag2

	
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
	
	flag1=0
	flag2=0
	
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
		if(type_ibm(i,j,k,nbl).eq.1) then
			flag1 = flag1 + 1
		elseif(type_ibm(i,j,k,nbl).eq.0) then
			flag2 = flag2 + 1
		endif
	
	enddo
	enddo
	enddo
	enddo
	
	no_body_pts = flag2
	no_fluid_pts = flag1
	
	print*,'Fluid Points Assigned'
	print*,'Body Points Assigned'
	
	print*,'Number of Fluid Points : ', no_fluid_pts
	print*,'Number of Body Points : ', no_body_pts
	

END



SUBROUTINE ghost_points()

	use declare_variables	
	implicit none
	
	integer ghost_pt_iter,flag1
	
	ghost_pt_iter = 0
	flag1 = 0

	
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
		if(type_ibm(i,j,k,nbl).eq.0) then
			if(i+1.ne.NI(nbl)+1) then
				if(type_ibm(i+1,j,k,nbl).eq.1) then
					ghost_pt_iter = ghost_pt_iter + 1
				endif
			elseif(i-1.ne.0) then
				if(type_ibm(i-1,j,k,nbl).eq.1) then
					ghost_pt_iter = ghost_pt_iter + 1
				endif
			elseif(j+1.ne.NJ(nbl)+1) then
				if(type_ibm(i,j+1,k,nbl).eq.1) then
					ghost_pt_iter = ghost_pt_iter + 1
				endif
			elseif(j-1.ne.0) then
				if(type_ibm(i,j-1,k,nbl).eq.1) then
					ghost_pt_iter = ghost_pt_iter + 1
				endif
			elseif(k+1.ne.NK(nbl)+1) then
				if(type_ibm(i,j,k+1,nbl).eq.1) then
					ghost_pt_iter = ghost_pt_iter + 1
				endif
			elseif(k-1.ne.0) then
				if(type_ibm(i,j,k-1,nbl).eq.1) then
					ghost_pt_iter = ghost_pt_iter + 1
				endif
			endif
		endif
	
	enddo
	enddo
	enddo
	enddo
	
	no_ghost_pts = ghost_pt_iter
	
	allocate(ghost_pt(3,no_ghost_pts))
	allocate(ghost_pt_idx(no_ghost_pts))
	
	ghost_pt_iter = 0
	
	
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
		if(type_ibm(i,j,k,nbl).eq.0) then
			if(i+1.ne.NI(nbl)+1) then
				if(type_ibm(i+1,j,k,nbl).eq.1) then
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
				endif
			elseif(i-1.ne.0) then
				if(type_ibm(i-1,j,k,nbl).eq.1) then
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
				endif
			elseif(j+1.ne.NJ(nbl)+1) then
				if(type_ibm(i,j+1,k,nbl).eq.1) then
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
				endif
			elseif(j-1.ne.0) then
				if(type_ibm(i,j-1,k,nbl).eq.1) then
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
				endif
			elseif(k+1.ne.NK(nbl)+1) then
				if(type_ibm(i,j,k+1,nbl).eq.1) then
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
				endif
			elseif(k-1.ne.0) then
				if(type_ibm(i,j,k-1,nbl).eq.1) then
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
				endif
			endif

		endif
	
	enddo
	enddo
	enddo
	enddo
	
	print*,'Ghost Points Assigned'
	print*,'Number of Ghost Points :', no_ghost_pts
	
END





SUBROUTINE boundary_intercept(x,y,z)

! x,y,z == xbg,ybg,zbg

! Global Vars
! pts
	
	use declare_variables
	use kdtree2_module
	
	implicit none
	
	real,dimension(nodes) :: x,y,z
	integer neighpts_ibm,pts,flag,near_node,indexing,node
	real dw
	real x1,x2,x3,xi,xg
	real y1,y2,y3,yi,yg
	real z1,z2,z3,zi,zg
	real a,b,c,d,tt
	real v1,v2,v3,v11,v12,v13,v21,v22,v23
	real vdv1,vdv2,v1dv2,v1dv1,v2dv2,dist2,dist2_new
	real alp,bet
	integer pp,elem,n1,n2,n3
	type(kdtree2_result),allocatable :: res(:)
	type(kdtree2),pointer :: tree
	real,allocatable :: global_ibm(:,:),qu_ibm(:)
	
	neighpts_ibm = 1
	
	allocate(res(neighpts_ibm))
	allocate(BII(3,no_ghost_pts))
	allocate(global_ibm(n_dim,nodes))
	allocate(qu_ibm(n_dim))
	
	do node = 1,nodes
		global_ibm(1,node) = x(node)
		global_ibm(2,node) = y(node) ! already allocated in ibm_type
		global_ibm(3,node) = z(node) ! comment this for n_dim = 2
	enddo
	
	tree => kdtree2_create(global_ibm,rearrange=.true.,sort=.true.)

  	do pts = 1,no_ghost_pts
  	
        flag = 0
  
		!call get_loc_index(ghost_pt_idx(pts))
		!print*,ghost_pt_idx(pts),i_loc,j_loc,k_loc,nbl_loc
		qu_ibm(1) = ghost_pt(1,pts) !Xgrid(i_loc,j_loc,k_loc,nbl_loc)
		qu_ibm(2) = ghost_pt(2,pts) !Ygrid(i_loc,j_loc,k_loc,nbl_loc) ! already allocated in ibm_type
		qu_ibm(3) = ghost_pt(3,pts) !Zgrid(i_loc,j_loc,k_loc,nbl_loc)
		call kdtree2_n_nearest(tree,qu_ibm,neighpts_ibm,res)	
       
		do indexing = 1,neighpts_ibm
			near_node = res(indexing)%idx
			dw = (res(indexing)%dis)**0.5		
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
		
			xg = qu_ibm(1)	!Ghost point coordinates
			yg = qu_ibm(2)
			zg = qu_ibm(3)	  
			
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
				BII(1,pts) = xi
				BII(2,pts) = yi
				BII(3,pts) = zi
				dist2 = dist2_new
				flag = 1
              endif        
            endif
  
  	    enddo
  
        if(flag.eq.0) then
           print*, 'ISSUE:',nbl_loc,i_loc,j_loc,k_loc,BII(1,pts),BII(2,pts),BII(3,pts)
        endif
  
  	enddo
	
	print*,'Boundary Intercept Calculated'

END


SUBROUTINE nn_fluid(nbrhd_pts,x,y,z)

! x,y,z == BI
	
	use declare_variables
	use kdtree2_module
	implicit none
	
	integer nbrhd_pts
	real x,y,z
	integer gridpts_iter,near_node,vicinity_pts_iter,indx
	real dist
	real,allocatable :: grid_pts(:,:),b_incpt(:)
	type(kdtree2),pointer :: tree_grid_pts
	type(kdtree2_result),allocatable :: results(:)
	
	allocate(grid_pts(3,NImax*NJmax*NKmax*nblocks))
	allocate(b_incpt(3))
	allocate(results(nbrhd_pts))
	
	gridpts_iter = 0
	R = 0
	
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
		dist = results(indx)%dis
		call get_loc_index(near_node)
		if(type_ibm(i_loc,j_loc,k_loc,nbl_loc).eq.1) then
			vicinity_pts_iter = vicinity_pts_iter + 1
		endif
	enddo
	
	no_vicinity_pts = vicinity_pts_iter
	vicinity_pts_iter = 0

	if(flag_nn_alloc.eq.0) then
		allocate(vicinity_pts(3,no_vicinity_pts))
	endif

	do indx = 1,nbrhd_pts
		near_node = results(indx)%idx
		dist = results(indx)%dis
		call get_loc_index(near_node)
		if(type_ibm(i_loc,j_loc,k_loc,nbl_loc).eq.1) then
			vicinity_pts_iter = vicinity_pts_iter + 1
			vicinity_pts(1,vicinity_pts_iter) = Xgrid(i,j,k,nbl)
			vicinity_pts(2,vicinity_pts_iter) = Ygrid(i,j,k,nbl)
			vicinity_pts(3,vicinity_pts_iter) = Zgrid(i,j,k,nbl)
			if(R.le.dist) then
				R = dist
			endif
		endif
	enddo

	flag_nn_alloc = 1

END



SUBROUTINE matrix_calculations()

! W is global

	use declare_variables
	implicit none
	external dgesvd

	
	
	integer nbrhd_pts,total_comp_pts,pp,qq,LWORK,INFO
	real xbi,ybi,zbi,xgp,ygp,zgp
	real,allocatable :: prime_coord(:,:)
	real(dp),allocatable :: S(:),U(:,:),VT(:,:),WORK(:),p_inv(:,:)
	real(dp),allocatable :: sig(:),sig_plus(:,:)
	real(dp) dummy(1,1)
	character(len=1),allocatable :: JOBU,JOBVT
	
	nbrhd_pts = 59
	total_comp_pts = nbrhd_pts+1 ! +1 ghost point
	flag_nn_alloc = 0
	flag_vandermonde_alloc = 0
	flag_pi_alloc = 0
	flag_A_alloc = 0
	JOBU = 'S'
	JOBVT = 'S'
	
	allocate(prime_coord(3,total_comp_pts))
	allocate(W(total_comp_pts,total_comp_pts,no_ghost_pts))
	
	do pp = 1,no_ghost_pts

		xbi = BII(1,pp)
		ybi = BII(2,pp)
		zbi = BII(3,pp)
		
		call nn_fluid(nbrhd_pts,xbi,ybi,zbi)
		
		xgp = ghost_pt(1,pp)
		ygp = ghost_pt(2,pp)
		zgp = ghost_pt(3,pp)
		
		prime_coord(1,1) = xgp - xbi
		prime_coord(2,1) = ygp - ybi
		prime_coord(3,1) = zgp - zbi
		W(1,1,pp) = 0.5*(1+cos(3.141*(prime_coord(1,1)**2+prime_coord(2,1)**2+prime_coord(3,1)**2)**0.5/R))
		
		do qq = 2,total_comp_pts
			
			prime_coord(1,qq) = vicinity_pts(1,qq-1) - xbi
			prime_coord(2,qq) = vicinity_pts(2,qq-1) - ybi
			prime_coord(3,qq) = vicinity_pts(3,qq-1) - zbi
			W(qq,qq,pp) = 0.5*(1+cos(3.141*(prime_coord(1,qq)**2+prime_coord(2,qq)**2+prime_coord(3,qq)**2)**0.5/R))
			
		enddo
		
		call ibm_coeff_vandermonde(pp,prime_coord,total_comp_pts)
		
		LWORK = max(1,5*min(total_comp_pts,L_N),3*min(total_comp_pts,L_N)+max(total_comp_pts,L_N))

		if(flag_pi_alloc.eq.0) then ! deallocate these and reset flag for updated total_comp_pts
			allocate(p_i(total_comp_pts,L_N,no_ghost_pts))
			!allocate(p_inv(total_comp_pts,L_N))
			!allocate(S(min(L_N,total_comp_pts)))
			!allocate(U(total_comp_pts,min(L_N,total_comp_pts)))
			!allocate(VT(min(L_N,total_comp_pts),L_N))
			flag_pi_alloc = 1
		endif
		
		if(flag_A_alloc.eq.0) then
			allocate(A_matrix(L_N,total_comp_pts,no_ghost_pts))
			flag_A_alloc = 1
		endif
		
		p_i(:,:,pp) = matmul(W(:,:,pp),V(:,:,pp))
		p_inv = p_i(:,:,pp)

		allocate(S(min(total_comp_pts,L_N)))
		allocate(U(total_comp_pts,min(L_N,total_comp_pts)))
		allocate(VT(min(L_N,total_comp_pts),L_N))
		
		allocate(sig(min(total_comp_pts,L_N)))
		allocate(sig_plus(min(total_comp_pts,L_N),min(total_comp_pts,L_N)))
		
		sig_plus = 0

		LWORK = -1
		call dgesvd(JOBU,JOBVT,total_comp_pts,L_N,p_i(:,:,pp), &
			total_comp_pts,S,U,total_comp_pts,VT,min(L_N,total_comp_pts),dummy,LWORK,INFO)

		allocate(WORK(int(dummy(1,1))))
		LWORK = max(1,int(dummy(1,1)),5*min(total_comp_pts,L_N),3*min(total_comp_pts,L_N)+max(total_comp_pts,L_N))
		call dgesvd(JOBU,JOBVT,total_comp_pts,L_N,p_i(:,:,pp), &
			total_comp_pts,S,U,total_comp_pts,VT,min(L_N,total_comp_pts),WORK,LWORK,INFO)

		do i = 1,min(L_N,total_comp_pts)
			if(s(i).le.10e-6) then
				sig(i) = 0
				sig_plus(i,i) = sig(i)
			elseif(s(i).ge.10e-6) then
				sig(i) = s(i)
				sig_plus(i,i) = 1/sig(i)
			endif
		enddo
		
		A_matrix(:,:,pp) = matmul(matmul(matmul(transpose(VT),sig_plus),transpose(U)),W(:,:,pp))

		deallocate(WORK)
		deallocate(S)
		deallocate(U)
		deallocate(VT)
		
		deallocate(sig)
		deallocate(sig_plus)		
	
	enddo

	print*,'Matrix Calculations Done'
	
END	


SUBROUTINE ibm_coeff_vandermonde(pp,prime_coord,total_comp_pts)

! L_N is a global variable
! V is a global variable

	use declare_variables
	implicit none
	
	real,dimension(3,total_comp_pts) :: prime_coord
	integer qq,pp,total_comp_pts
	
	N_vandermonde = 3
	
	
	if(N_vandermonde.eq.1.and.n_dim.eq.2) then
		L_N = 3

	elseif(N_vandermonde.eq.1.and.n_dim.eq.3) then
		L_N = 4

	elseif(N_vandermonde.eq.2.and.n_dim.eq.2) then
		L_N = 6

	elseif(N_vandermonde.eq.2.and.n_dim.eq.3) then
		L_N = 10
		
	elseif(N_vandermonde.eq.3.and.n_dim.eq.2) then
		L_N = 10
		
	elseif(N_vandermonde.eq.3.and.n_dim.eq.3) then
		L_N = 20
		
	elseif(N_vandermonde.eq.4.and.n_dim.eq.2) then
		L_N = 15
		
	elseif(N_vandermonde.eq.4.and.n_dim.eq.3) then
		L_N = 35
		
	endif
	
	if(flag_vandermonde_alloc.eq.0) then
		allocate(V(total_comp_pts,L_N,no_ghost_pts))
	endif
	
	
	if(N_vandermonde.eq.1.and.n_dim.eq.2) then
		L_N = 3
		do qq = 1,total_comp_pts
			V(qq,1,pp) = 1
			V(qq,2,pp) = prime_coord(1,qq)
			V(qq,3,pp) = prime_coord(2,qq)
		enddo
	elseif(N_vandermonde.eq.1.and.n_dim.eq.3) then
		L_N = 4
		do qq = 1,total_comp_pts
			V(qq,1,pp) = 1
			V(qq,2,pp) = prime_coord(1,qq)
			V(qq,3,pp) = prime_coord(2,qq)
			V(qq,4,pp) = prime_coord(3,qq)
		enddo
	elseif(N_vandermonde.eq.2.and.n_dim.eq.2) then
		L_N = 6
		do qq = 1,total_comp_pts
			V(qq,1,pp) = 1
			V(qq,2,pp) = prime_coord(1,qq)
			V(qq,3,pp) = prime_coord(2,qq)
			V(qq,4,pp) = prime_coord(1,qq)*prime_coord(1,qq)
			V(qq,5,pp) = prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,6,pp) = prime_coord(2,qq)*prime_coord(1,qq)
		enddo
	elseif(N_vandermonde.eq.2.and.n_dim.eq.3) then
		L_N = 10
		do qq = 1,total_comp_pts
			V(qq,1,pp) = 1
			V(qq,2,pp) = prime_coord(1,qq)
			V(qq,3,pp) = prime_coord(2,qq)
			V(qq,4,pp) = prime_coord(3,qq)
			V(qq,5,pp) = prime_coord(1,qq)*prime_coord(1,qq)
			V(qq,6,pp) = prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,7,pp) = prime_coord(3,qq)*prime_coord(3,qq)
			V(qq,8,pp) = prime_coord(1,qq)*prime_coord(2,qq)
			V(qq,9,pp) = prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,10,pp) = prime_coord(3,qq)*prime_coord(1,qq)
		enddo
	elseif(N_vandermonde.eq.3.and.n_dim.eq.2) then
		L_N = 10
		do qq = 1,total_comp_pts
			V(qq,1,pp) = 1
			V(qq,2,pp) = prime_coord(1,qq)
			V(qq,3,pp) = prime_coord(2,qq)
			V(qq,4,pp) = prime_coord(1,qq)*prime_coord(1,qq)
			V(qq,5,pp) = prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,6,pp) = prime_coord(2,qq)*prime_coord(1,qq)
			V(qq,7,pp) = prime_coord(1,qq)*prime_coord(1,qq)*prime_coord(1,qq)
			V(qq,8,pp) = prime_coord(2,qq)*prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,9,pp) = prime_coord(1,qq)*prime_coord(1,qq)*prime_coord(2,qq)
			V(qq,10,pp) = prime_coord(1,qq)*prime_coord(2,qq)*prime_coord(2,qq)
		enddo
	elseif(N_vandermonde.eq.3.and.n_dim.eq.3) then
		L_N = 20
		do qq = 1,total_comp_pts
			V(qq,1,pp) = 1
			V(qq,2,pp) = prime_coord(1,qq)
			V(qq,3,pp) = prime_coord(2,qq)
			V(qq,4,pp) = prime_coord(3,qq)
			V(qq,5,pp) = prime_coord(1,qq)*prime_coord(1,qq)
			V(qq,6,pp) = prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,7,pp) = prime_coord(3,qq)*prime_coord(3,qq)
			V(qq,8,pp) = prime_coord(1,qq)*prime_coord(2,qq)
			V(qq,9,pp) = prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,10,pp) = prime_coord(3,qq)*prime_coord(1,qq)
			V(qq,11,pp) = prime_coord(1,qq)**3
			V(qq,12,pp) = prime_coord(2,qq)**3
			V(qq,13,pp) = prime_coord(3,qq)**3
			V(qq,14,pp) = prime_coord(1,qq)**2*prime_coord(2,qq)
			V(qq,15,pp) = prime_coord(2,qq)**2*prime_coord(3,qq)
			V(qq,16,pp) = prime_coord(3,qq)**2*prime_coord(1,qq)
			V(qq,17,pp) = prime_coord(1,qq)*prime_coord(2,qq)**2
			V(qq,18,pp) = prime_coord(2,qq)*prime_coord(3,qq)**2
			V(qq,19,pp) = prime_coord(3,qq)*prime_coord(1,qq)**2
			V(qq,20,pp) = prime_coord(1,qq)*prime_coord(2,qq)*prime_coord(3,qq)
		enddo
	elseif(N_vandermonde.eq.4.and.n_dim.eq.2) then
		L_N = 15
		do qq = 1,total_comp_pts
			V(qq,1,pp) = 1
			V(qq,2,pp) = prime_coord(1,qq)
			V(qq,3,pp) = prime_coord(2,qq)
			V(qq,4,pp) = prime_coord(1,qq)*prime_coord(1,qq)
			V(qq,5,pp) = prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,6,pp) = prime_coord(2,qq)*prime_coord(1,qq)
			V(qq,7,pp) = prime_coord(1,qq)*prime_coord(1,qq)*prime_coord(1,qq)
			V(qq,8,pp) = prime_coord(2,qq)*prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,9,pp) = prime_coord(1,qq)*prime_coord(1,qq)*prime_coord(2,qq)
			V(qq,10,pp) = prime_coord(1,qq)*prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,11,pp) = prime_coord(1,qq)**4
			V(qq,12,pp) = prime_coord(2,qq)**4
			V(qq,13,pp) = prime_coord(1,qq)**3*prime_coord(2,qq)
			V(qq,14,pp) = prime_coord(1,qq)**2*prime_coord(2,qq)**2
			V(qq,15,pp) = prime_coord(1,qq)*prime_coord(2,qq)**3
		enddo
	elseif(N_vandermonde.eq.4.and.n_dim.eq.3) then
		L_N = 35
		do qq = 1,total_comp_pts
			V(qq,1,pp) = 1
			V(qq,2,pp) = prime_coord(1,qq)
			V(qq,3,pp) = prime_coord(2,qq)
			V(qq,4,pp) = prime_coord(3,qq)
			V(qq,5,pp) = prime_coord(1,qq)*prime_coord(1,qq)
			V(qq,6,pp) = prime_coord(2,qq)*prime_coord(2,qq)
			V(qq,7,pp) = prime_coord(3,qq)*prime_coord(3,qq)
			V(qq,8,pp) = prime_coord(1,qq)*prime_coord(2,qq)
			V(qq,9,pp) = prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,10,ii) = prime_coord(3,qq)*prime_coord(1,qq)
			V(qq,11,pp) = prime_coord(1,qq)**3
			V(qq,12,pp) = prime_coord(2,qq)**3
			V(qq,13,pp) = prime_coord(3,qq)**3
			V(qq,14,pp) = prime_coord(1,qq)**2*prime_coord(2,qq)
			V(qq,15,pp) = prime_coord(2,qq)**2*prime_coord(3,qq)
			V(qq,16,pp) = prime_coord(3,qq)**2*prime_coord(1,qq)
			V(qq,17,pp) = prime_coord(1,qq)*prime_coord(2,qq)**2
			V(qq,18,pp) = prime_coord(2,qq)*prime_coord(3,qq)**2
			V(qq,19,pp) = prime_coord(3,qq)*prime_coord(1,qq)**2
			V(qq,20,pp) = prime_coord(1,qq)*prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,21,pp) = prime_coord(1,qq)**4
			V(qq,22,pp) = prime_coord(2,qq)**4
			V(qq,23,pp) = prime_coord(3,qq)**4
			V(qq,24,pp) = prime_coord(1,qq)**3*prime_coord(2,qq)
			V(qq,25,pp) = prime_coord(2,qq)**3*prime_coord(3,qq)
			V(qq,26,pp) = prime_coord(3,qq)**3*prime_coord(1,qq)
			V(qq,27,pp) = prime_coord(1,qq)*prime_coord(2,qq)**3
			V(qq,28,pp) = prime_coord(2,qq)*prime_coord(3,qq)**3
			V(qq,29,pp) = prime_coord(3,qq)*prime_coord(1,qq)**3
			V(qq,30,pp) = prime_coord(1,qq)**2*prime_coord(2,qq)**2
			V(qq,31,pp) = prime_coord(2,qq)**2*prime_coord(3,qq)**2
			V(qq,32,pp) = prime_coord(3,qq)**2*prime_coord(1,qq)**2
			V(qq,33,pp) = prime_coord(1,qq)**2*prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,34,pp) = prime_coord(2,qq)**2*prime_coord(3,qq)*prime_coord(1,qq)
			V(qq,35,pp) = prime_coord(3,qq)**2*prime_coord(1,qq)*prime_coord(2,qq)
		enddo
	endif
	
	flag_vandermonde_alloc = 1
	
END




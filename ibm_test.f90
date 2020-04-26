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
	
	deallocate(global_ibm)
	deallocate(qu_ibm)
	deallocate(results_ibm)
	

END


SUBROUTINE ghost_points()

	use declare_variables	
	implicit none
	
	integer ghost_pt_iter,flag1
	
	ghost_pt_iter = 0
	flag1 = 0

	
	do nbl = 1,nblocks
	do k = 1,NK(nbl) ! for 2 D because there are only 3 pts in z dirn
	do j = 2,NJ(nbl)-1
	do i = 2,NI(nbl)-1
	
		if(type_ibm(i,j,k,nbl).eq.0) then
		
			!if(i+1.ne.NI(nbl)+1) then
				if(type_ibm(i+1,j,k,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'i+'
					ghost_pt_iter = ghost_pt_iter + 1
					cycle
					print*,'should not get printed'
				!endif
			!elseif(i-1.ne.0) then
				elseif(type_ibm(i-1,j,k,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'i-'
					ghost_pt_iter = ghost_pt_iter + 1
					cycle
				!endif
			!elseif(j+1.ne.NJ(nbl)+1) then
				elseif(type_ibm(i,j+1,k,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'j+'
					ghost_pt_iter = ghost_pt_iter + 1
					cycle
				!endif
			!elseif(j-1.ne.0) then
				elseif(type_ibm(i,j-1,k,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'j-'
					ghost_pt_iter = ghost_pt_iter + 1
					cycle
				!endif
			!elseif(k+1.ne.NK(nbl)+1) then
				elseif(type_ibm(i,j,k+1,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'k+'
					ghost_pt_iter = ghost_pt_iter + 1
					cycle
				!endif
			!elseif(k-1.ne.0) then
				elseif(type_ibm(i,j,k-1,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'k-'
					ghost_pt_iter = ghost_pt_iter + 1
					cycle
				endif
			!endif
		endif
	
	enddo
	enddo
	enddo
	enddo
	
	no_ghost_pts = ghost_pt_iter
	
	allocate(ghost_pt(3,no_ghost_pts))
	allocate(ghost_pt_idx(no_ghost_pts))
	
	ghost_pt_iter = 0
	flag1 = 0
	
	
	do nbl = 1,nblocks
	do k = 1,NK(nbl) ! for 2 D because there are only 3 pts in z dirn
	do j = 2,NJ(nbl)-1
	do i = 2,NI(nbl)-1
		if(type_ibm(i,j,k,nbl).eq.0) then
			!if(i+1.ne.NI(nbl)+1) then
				if(type_ibm(i+1,j,k,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'i plus'
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
					cycle
					print*,'this shouldnt be printed'
				!endif
			!elseif(i-1.ne.0) then
				elseif(type_ibm(i-1,j,k,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'i minus'
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
					cycle
				!endif
			!elseif(j+1.ne.NJ(nbl)+1) then
				elseif(type_ibm(i,j+1,k,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'j plus'
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
					cycle
				!endif
			!elseif(j-1.ne.0) then
				elseif(type_ibm(i,j-1,k,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'j minus'
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
					cycle
				!endif
			!elseif(k+1.ne.NK(nbl)+1) then
				elseif(type_ibm(i,j,k+1,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'k plus'
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
					cycle
				!endif
			!elseif(k-1.ne.0) then
				elseif(type_ibm(i,j,k-1,nbl).eq.1) then
				!flag1 = flag1 + 1
				!print*,flag1,'k minus'
					type_ibm(i,j,k,nbl) = -1
					ghost_pt_iter = ghost_pt_iter + 1
					call get_global_index(i,j,k,nbl)
					ghost_pt_idx(ghost_pt_iter) = global_index
					ghost_pt(1,ghost_pt_iter) = Xgrid(i,j,k,nbl)
					ghost_pt(2,ghost_pt_iter) = Ygrid(i,j,k,nbl)
					ghost_pt(3,ghost_pt_iter) = Zgrid(i,j,k,nbl)
					cycle
				endif
			!endif

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
	
	deallocate(res)
	deallocate(global_ibm)
	deallocate(qu_ibm)

END


SUBROUTINE nn_fluid(nbrhd_pts,x,y,z)

! x,y,z == BI
	
	use declare_variables
	use kdtree2_module
	implicit none
	
	integer nbrhd_pts
	real x,y,z
	integer near_node,vicinity_pts_iter,indx,dist_dum
	real dist
	type(kdtree2_result),allocatable :: results(:)
	
	allocate(results(nbrhd_pts))

	R = 0
	dist_dum = 10e5	
	
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

	allocate(vicinity_pts(3,no_vicinity_pts))
	allocate(vicinity_pt_idx(no_vicinity_pts))
	
	do indx = 1,nbrhd_pts
		near_node = results(indx)%idx
		dist = results(indx)%dis
		call get_loc_index(near_node)
		if(type_ibm(i_loc,j_loc,k_loc,nbl_loc).eq.1) then
			vicinity_pts_iter = vicinity_pts_iter + 1
			vicinity_pt_idx(vicinity_pts_iter) = near_node
			vicinity_pts(1,vicinity_pts_iter) = Xgrid(i_loc,j_loc,k_loc,nbl_loc)
			vicinity_pts(2,vicinity_pts_iter) = Ygrid(i_loc,j_loc,k_loc,nbl_loc)
			vicinity_pts(3,vicinity_pts_iter) = Zgrid(i_loc,j_loc,k_loc,nbl_loc)
			if(dist<dist_dum) then
				bfp_idx = near_node
				bfp(1) = Xgrid(i_loc,j_loc,k_loc,nbl_loc)
				bfp(2) = Ygrid(i_loc,j_loc,k_loc,nbl_loc)
				bfp(3) = Zgrid(i_loc,j_loc,k_loc,nbl_loc)
				dist_dum = dist
			endif
			if(R.le.dist) then
				R = dist
			endif
		endif
	enddo

	flag_nn_alloc = 1

	deallocate(results)
	!print*,results(1)%idx
	
END


SUBROUTINE matrix_calculations()

! W is global

	use declare_variables
	implicit none
	external dgesvd	
	
	integer nbrhd_pts,pp,qq,LWORK,INFO
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
	
	allocate(dim_A(2,no_ghost_pts))
	allocate(A_matrix(nbrhd_pts,nbrhd_pts,no_ghost_pts))
	
	A_matrix = 0

	do pp = 1,no_ghost_pts

		xbi = BII(1,pp)
		ybi = BII(2,pp)
		zbi = BII(3,pp)

		call nn_fluid(nbrhd_pts,xbi,ybi,zbi)
	
		total_comp_pts = no_vicinity_pts + 1

		allocate(prime_coord(3,total_comp_pts))
		allocate(W(total_comp_pts,total_comp_pts))
		W = 0
		
		!if(pp.eq.1) then
		!open(5, form = 'formatted', file = 'vicinity_pts.txt')
		!write(5,*) 'x',',','y',',','z'
		!do i=1,no_vicinity_pts
		!write(5,'(*(G0.6,:,","))') vicinity_pts(1,i),vicinity_pts(2,i),vicinity_pts(3,i)
		!enddo
		!close(5)
		!endif
		
		xgp = ghost_pt(1,pp)
		ygp = ghost_pt(2,pp)
		zgp = ghost_pt(3,pp)
		
		prime_coord(1,1) = xgp - xbi
		prime_coord(2,1) = ygp - ybi
		prime_coord(3,1) = zgp - zbi
		!print*,prime_coord(3,1)
		W(1,1) = 0.5*(1+cos(3.141*(prime_coord(1,1)**2+prime_coord(2,1)**2+prime_coord(3,1)**2)**0.5/R))

		do qq = 2,total_comp_pts
			
			prime_coord(1,qq) = vicinity_pts(1,qq-1) - xbi
			prime_coord(2,qq) = vicinity_pts(2,qq-1) - ybi
			prime_coord(3,qq) = vicinity_pts(3,qq-1) - zbi
			W(qq,qq) = 0.5*(1+cos(3.141*(prime_coord(1,qq)**2+prime_coord(2,qq)**2+prime_coord(3,qq)**2)**0.5/R))

		enddo

		call ibm_coeff_vandermonde(prime_coord)
		
		allocate(p_i(total_comp_pts,L_N))

		LWORK = max(1,5*min(total_comp_pts,L_N),3*min(total_comp_pts,L_N)+max(total_comp_pts,L_N))

		p_i = matmul(W,V)

		allocate(S(min(total_comp_pts,L_N)))
		allocate(U(total_comp_pts,min(L_N,total_comp_pts)))
		allocate(VT(min(L_N,total_comp_pts),L_N))
		
		allocate(sig(min(total_comp_pts,L_N)))
		allocate(sig_plus(min(total_comp_pts,L_N),min(total_comp_pts,L_N)))
		
		dim_A(1,pp) = L_N
		dim_A(2,pp) = total_comp_pts
		
		sig_plus = 0

		LWORK = -1
		call dgesvd(JOBU,JOBVT,total_comp_pts,L_N,p_i, &
			total_comp_pts,S,U,total_comp_pts,VT,min(L_N,total_comp_pts),dummy,LWORK,INFO)

		allocate(WORK(int(dummy(1,1))))
		LWORK = max(1,int(dummy(1,1)),5*min(total_comp_pts,L_N),3*min(total_comp_pts,L_N)+max(total_comp_pts,L_N))
		call dgesvd(JOBU,JOBVT,total_comp_pts,L_N,p_i, &
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
		
		!print*,maxval(s)/maxval(sig_plus),info

		A_matrix(1:dim_A(1,pp),1:dim_A(2,pp),pp) = matmul(matmul(matmul(transpose(VT),sig_plus),transpose(U)),W)

		deallocate(WORK)
		deallocate(S)
		deallocate(U)
		deallocate(VT)
		
		deallocate(sig)
		deallocate(sig_plus)

		deallocate(p_i)
		deallocate(W)
		deallocate(V)
		deallocate(prime_coord)
		
		deallocate(vicinity_pts)
		deallocate(vicinity_pt_idx)
	
	enddo

	print*,'Matrix Calculations Done'
	
END	


SUBROUTINE ibm_coeff_vandermonde(prime_coord)

! L_N is a global variable
! V is a global variable

	use declare_variables
	implicit none
	real,dimension(3,total_comp_pts) :: prime_coord
	integer qq,pp,ndim
	
	ndim = 2
	N_vandermonde = 3
	
	if(N_vandermonde.eq.1.and.ndim.eq.2) then
		L_N = 3

	elseif(N_vandermonde.eq.1.and.ndim.eq.3) then
		L_N = 4

	elseif(N_vandermonde.eq.2.and.ndim.eq.2) then
		L_N = 6

	elseif(N_vandermonde.eq.2.and.ndim.eq.3) then
		L_N = 10
		
	elseif(N_vandermonde.eq.3.and.ndim.eq.2) then
		L_N = 10
		
	elseif(N_vandermonde.eq.3.and.ndim.eq.3) then
		L_N = 20
		
	elseif(N_vandermonde.eq.4.and.ndim.eq.2) then
		L_N = 15
		
	elseif(N_vandermonde.eq.4.and.ndim.eq.3) then
		L_N = 35
		
	endif
	
	allocate(V(total_comp_pts,L_N))
	V = 0
	
	
	if(N_vandermonde.eq.1.and.ndim.eq.2) then
		L_N = 3
		do qq = 1,total_comp_pts
			V(qq,1) = 1
			V(qq,2) = prime_coord(1,qq)
			V(qq,3) = prime_coord(2,qq)
		enddo
	elseif(N_vandermonde.eq.1.and.ndim.eq.3) then
		L_N = 4
		do qq = 1,total_comp_pts
			V(qq,1) = 1
			V(qq,2) = prime_coord(1,qq)
			V(qq,3) = prime_coord(2,qq)
			V(qq,4) = prime_coord(3,qq)
		enddo
	elseif(N_vandermonde.eq.2.and.ndim.eq.2) then
		L_N = 6
		do qq = 1,total_comp_pts
			V(qq,1) = 1
			V(qq,2) = prime_coord(1,qq)
			V(qq,3) = prime_coord(2,qq)
			V(qq,4) = prime_coord(1,qq)**2
			V(qq,5) = prime_coord(2,qq)**2
			V(qq,6) = prime_coord(2,qq)*prime_coord(1,qq)
		enddo
	elseif(N_vandermonde.eq.2.and.ndim.eq.3) then
		L_N = 10
		do qq = 1,total_comp_pts
			V(qq,1) = 1
			V(qq,2) = prime_coord(1,qq)
			V(qq,3) = prime_coord(2,qq)
			V(qq,4) = prime_coord(3,qq)
			V(qq,5) = prime_coord(1,qq)**2
			V(qq,6) = prime_coord(2,qq)**2
			V(qq,7) = prime_coord(3,qq)**2
			V(qq,8) = prime_coord(1,qq)*prime_coord(2,qq)
			V(qq,9) = prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,10) = prime_coord(3,qq)*prime_coord(1,qq)
		enddo
	elseif(N_vandermonde.eq.3.and.ndim.eq.2) then
		L_N = 10
		do qq = 1,total_comp_pts
			V(qq,1) = 1
			V(qq,2) = prime_coord(1,qq)
			V(qq,3) = prime_coord(2,qq)
			V(qq,4) = prime_coord(1,qq)**2
			V(qq,5) = prime_coord(2,qq)**2
			V(qq,6) = prime_coord(2,qq)*prime_coord(1,qq)
			V(qq,7) = prime_coord(1,qq)**3
			V(qq,8) = prime_coord(2,qq)**3
			V(qq,9) = prime_coord(1,qq)**2*prime_coord(2,qq)
			V(qq,10) = prime_coord(1,qq)*prime_coord(2,qq)*2
		enddo
	elseif(N_vandermonde.eq.3.and.ndim.eq.3) then
		L_N = 20
		do qq = 1,total_comp_pts
			V(qq,1) = 1
			V(qq,2) = prime_coord(1,qq)
			V(qq,3) = prime_coord(2,qq)
			V(qq,4) = prime_coord(3,qq)
			V(qq,5) = prime_coord(1,qq)**2
			V(qq,6) = prime_coord(2,qq)**2
			V(qq,7) = prime_coord(3,qq)**2
			V(qq,8) = prime_coord(1,qq)*prime_coord(2,qq)
			V(qq,9) = prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,10) = prime_coord(3,qq)*prime_coord(1,qq)
			V(qq,11) = prime_coord(1,qq)**3
			V(qq,12) = prime_coord(2,qq)**3
			V(qq,13) = prime_coord(3,qq)**3
			V(qq,14) = prime_coord(1,qq)**2*prime_coord(2,qq)
			V(qq,15) = prime_coord(2,qq)**2*prime_coord(3,qq)
			V(qq,16) = prime_coord(3,qq)**2*prime_coord(1,qq)
			V(qq,17) = prime_coord(1,qq)*prime_coord(2,qq)**2
			V(qq,18) = prime_coord(2,qq)*prime_coord(3,qq)**2
			V(qq,19) = prime_coord(3,qq)*prime_coord(1,qq)**2
			V(qq,20) = prime_coord(1,qq)*prime_coord(2,qq)*prime_coord(3,qq)
		enddo
	elseif(N_vandermonde.eq.4.and.ndim.eq.2) then
		L_N = 15
		do qq = 1,total_comp_pts
			V(qq,1) = 1
			V(qq,2) = prime_coord(1,qq)
			V(qq,3) = prime_coord(2,qq)
			V(qq,4) = prime_coord(1,qq)**2
			V(qq,5) = prime_coord(2,qq)**2
			V(qq,6) = prime_coord(2,qq)*prime_coord(1,qq)
			V(qq,7) = prime_coord(1,qq)**3
			V(qq,8) = prime_coord(2,qq)**3
			V(qq,9) = prime_coord(1,qq)**2*prime_coord(2,qq)
			V(qq,10) = prime_coord(1,qq)*prime_coord(2,qq)**2
			V(qq,11) = prime_coord(1,qq)**4
			V(qq,12) = prime_coord(2,qq)**4
			V(qq,13) = prime_coord(1,qq)**3*prime_coord(2,qq)
			V(qq,14) = prime_coord(1,qq)**2*prime_coord(2,qq)**2
			V(qq,15) = prime_coord(1,qq)*prime_coord(2,qq)**3
		enddo
	elseif(N_vandermonde.eq.4.and.ndim.eq.3) then
		L_N = 35
		do qq = 1,total_comp_pts
			V(qq,1) = 1
			V(qq,2) = prime_coord(1,qq)
			V(qq,3) = prime_coord(2,qq)
			V(qq,4) = prime_coord(3,qq)
			V(qq,5) = prime_coord(1,qq)**2
			V(qq,6) = prime_coord(2,qq)**2
			V(qq,7) = prime_coord(3,qq)**2
			V(qq,8) = prime_coord(1,qq)*prime_coord(2,qq)
			V(qq,9) = prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,10) = prime_coord(3,qq)*prime_coord(1,qq)
			V(qq,11) = prime_coord(1,qq)**3
			V(qq,12) = prime_coord(2,qq)**3
			V(qq,13) = prime_coord(3,qq)**3
			V(qq,14) = prime_coord(1,qq)**2*prime_coord(2,qq)
			V(qq,15) = prime_coord(2,qq)**2*prime_coord(3,qq)
			V(qq,16) = prime_coord(3,qq)**2*prime_coord(1,qq)
			V(qq,17) = prime_coord(1,qq)*prime_coord(2,qq)**2
			V(qq,18) = prime_coord(2,qq)*prime_coord(3,qq)**2
			V(qq,19) = prime_coord(3,qq)*prime_coord(1,qq)**2
			V(qq,20) = prime_coord(1,qq)*prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,21) = prime_coord(1,qq)**4
			V(qq,22) = prime_coord(2,qq)**4
			V(qq,23) = prime_coord(3,qq)**4
			V(qq,24) = prime_coord(1,qq)**3*prime_coord(2,qq)
			V(qq,25) = prime_coord(2,qq)**3*prime_coord(3,qq)
			V(qq,26) = prime_coord(3,qq)**3*prime_coord(1,qq)
			V(qq,27) = prime_coord(1,qq)*prime_coord(2,qq)**3
			V(qq,28) = prime_coord(2,qq)*prime_coord(3,qq)**3
			V(qq,29) = prime_coord(3,qq)*prime_coord(1,qq)**3
			V(qq,30) = prime_coord(1,qq)**2*prime_coord(2,qq)**2
			V(qq,31) = prime_coord(2,qq)**2*prime_coord(3,qq)**2
			V(qq,32) = prime_coord(3,qq)**2*prime_coord(1,qq)**2
			V(qq,33) = prime_coord(1,qq)**2*prime_coord(2,qq)*prime_coord(3,qq)
			V(qq,34) = prime_coord(2,qq)**2*prime_coord(3,qq)*prime_coord(1,qq)
			V(qq,35) = prime_coord(3,qq)**2*prime_coord(1,qq)*prime_coord(2,qq)
		enddo
	endif
	
	flag_vandermonde_alloc = 1
	
END


SUBROUTINE phi_gp(PHI_W,PHI,nvars,PHIC,nvarsc)

	use declare_variables
	use kdtree2_module
	implicit none
	
	integer nvars,pp,gg,i_l,j_l,k_l,nbl_l,nbrhd_pts,nvarsc
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
	real,dimension(nvars) :: PHI_W
	real,dimension(NImax,NJmax,NKmax,nblocks,nvarsc) :: PHIC 
	real xbi,ybi,zbi,x_normal,y_normal,z_normal,Etotal,norm,xgp,ygp,zgp

	
	nbrhd_pts = 59 ! same as that in matrix calculations
	
	do gg = 1,no_ghost_pts
		
		call get_loc_index(ghost_pt_idx(gg))
		i_l = i_loc
		j_l = j_loc
		k_l = k_loc
		nbl_l = nbl_loc
		xbi = BII(1,gg)
		ybi = BII(2,gg)
		zbi = BII(3,gg)
		xgp = ghost_pt(1,gg)
		ygp = ghost_pt(2,gg)
		zgp = ghost_pt(3,gg)
		norm = ((xbi-xgp)**2+(ybi-ygp)**2+(zbi-zgp)**2)**0.5
		
		call nn_fluid(nbrhd_pts,xbi,ybi,zbi) ! no of vicinity pts available for each BII

		x_normal = (xbi-xgp)/norm
		y_normal = (ybi-ygp)/norm
		z_normal = (zbi-zgp)/norm

		PHI(i_l,j_l,k_l,nbl_l,:) = 0
		PHIC(i_l,j_l,k_l,nbl_l,:) = 0
		!print*,no_vicinity_pts,dim_A(2,gg)

		do pp = 2,dim_A(2,gg) ! == no_vicinity_pts+1
			
			call get_loc_index(vicinity_pt_idx(pp-1)) ! all vicinity points from 1 to no_vicinity_pts
			PHI(i_l,j_l,k_l,nbl_l,2) = PHI(i_l,j_l,k_l,nbl_l,2) + A_matrix(1,pp,gg)*PHI(i_loc,j_loc,k_loc,nbl_loc,2)
			PHI(i_l,j_l,k_l,nbl_l,3) = PHI(i_l,j_l,k_l,nbl_l,3) + A_matrix(1,pp,gg)*PHI(i_loc,j_loc,k_loc,nbl_loc,3)
			PHI(i_l,j_l,k_l,nbl_l,4) = PHI(i_l,j_l,k_l,nbl_l,4) + A_matrix(1,pp,gg)*PHI(i_loc,j_loc,k_loc,nbl_loc,4)
			PHI(i_l,j_l,k_l,nbl_l,5) = PHI(i_l,j_l,k_l,nbl_l,5) + (x_normal*A_matrix(2,pp,gg)+y_normal*A_matrix(3,pp,gg) &
				+z_normal*A_matrix(4,pp,gg))*PHI(i_loc,j_loc,k_loc,nbl_loc,5)
			PHI(i_l,j_l,k_l,nbl_l,6) = PHI(i_l,j_l,k_l,nbl_l,6) + (x_normal*A_matrix(2,pp,gg)+y_normal*A_matrix(3,pp,gg) &
				+z_normal*A_matrix(4,pp,gg))*PHI(i_loc,j_loc,k_loc,nbl_loc,6)

		enddo

		PHI(i_l,j_l,k_l,nbl_l,2) = (PHI_W(2) - PHI(i_l,j_l,k_l,nbl_l,2))/A_matrix(1,1,gg)
		PHI(i_l,j_l,k_l,nbl_l,3) = (PHI_W(3) - PHI(i_l,j_l,k_l,nbl_l,3))/A_matrix(1,1,gg)
		PHI(i_l,j_l,k_l,nbl_l,4) = (PHI_W(4) - PHI(i_l,j_l,k_l,nbl_l,4))/A_matrix(1,1,gg)
		PHI(i_l,j_l,k_l,nbl_l,5) = (PHI_W(5) - PHI(i_l,j_l,k_l,nbl_l,5))/ &
			(x_normal*A_matrix(2,1,gg)+y_normal*A_matrix(3,1,gg)+z_normal*A_matrix(4,1,gg))
		PHI(i_l,j_l,k_l,nbl_l,6) = (PHI_W(6) - PHI(i_l,j_l,k_l,nbl_l,6))/ &
			(x_normal*A_matrix(2,1,gg)+y_normal*A_matrix(3,1,gg)+z_normal*A_matrix(4,1,gg))	
		PHI(i_l,j_l,k_l,nbl_l,1) = gamma*Mach**2.d0*PHI(i_l,j_l,k_l,nbl_l,5)/PHI(i_l,j_l,k_l,nbl_l,6)
		
		!PHI(i_l,j_l,k_l,nbl_l,2) = (0 - PHI(i_l,j_l,k_l,nbl_l,2))/A_matrix(1,1,gg)
		!PHI(i_l,j_l,k_l,nbl_l,3) = (0 - PHI(i_l,j_l,k_l,nbl_l,3))/A_matrix(1,1,gg)
		!PHI(i_l,j_l,k_l,nbl_l,4) = (0 - PHI(i_l,j_l,k_l,nbl_l,4))/A_matrix(1,1,gg)
		!PHI(i_l,j_l,k_l,nbl_l,5) = (0 - PHI(i_l,j_l,k_l,nbl_l,5))/ &
		!	(x_normal*A_matrix(2,1,gg)+y_normal*A_matrix(3,1,gg)+z_normal*A_matrix(4,1,gg))
		!PHI(i_l,j_l,k_l,nbl_l,6) = (0 - PHI(i_l,j_l,k_l,nbl_l,6))/ &
		!	(x_normal*A_matrix(2,1,gg)+y_normal*A_matrix(3,1,gg)+z_normal*A_matrix(4,1,gg))	
		!PHI(i_l,j_l,k_l,nbl_l,1) = gamma*Mach**2.d0*PHI(i_l,j_l,k_l,nbl_l,5)/PHI(i_l,j_l,k_l,nbl_l,6)

		PHIC(i_l,j_l,k_l,nbl_l,1) = PHI(i_l,j_l,k_l,nbl_l,1)
		PHIC(i_l,j_l,k_l,nbl_l,2) = PHI(i_l,j_l,k_l,nbl_l,1)*PHI(i_l,j_l,k_l,nbl_l,2)
		PHIC(i_l,j_l,k_l,nbl_l,3) = PHI(i_l,j_l,k_l,nbl_l,1)*PHI(i_l,j_l,k_l,nbl_l,3)
		PHIC(i_l,j_l,k_l,nbl_l,4) = PHI(i_l,j_l,k_l,nbl_l,1)*PHI(i_l,j_l,k_l,nbl_l,4)

		Etotal = PHI(i_l,j_l,k_l,nbl_l,5)/(PHI(i_l,j_l,k_l,nbl_l,1)*(gamma-1.d0)) + &
			0.5d0*(PHI(i_l,j_l,k_l,nbl_l,2)**2+PHI(i_l,j_l,k_l,nbl_l,3)**2+PHI(i_l,j_l,k_l,nbl_l,4)**2)

		PHIC(i_l,j_l,k_l,nbl_l,5) = PHI(i_l,j_l,k_l,nbl_l,1)*(Etotal)

		deallocate(vicinity_pts)
		deallocate(vicinity_pt_idx)

	enddo
	

END


SUBROUTINE boundary_fluid_points()

	use declare_variables
	implicit none

	integer node,bfp_iter,bfp_temp
	
	bfp_iter = 0

	do node=1,nodes
	
		call nn_fluid(5,xbg(node),ybg(node),zbg(node)) 
		
		!print*,no_vicinity_pts ! at no point should no_vicinity_pts go to 0
		!print*,bfp_idx
		if(node.ne.1) then
			if(bfp_temp.ne.bfp_idx) then
				bfp_iter = bfp_iter + 1
				bfp_temp = bfp_idx
			endif
		else
			bfp_iter = bfp_iter + 1
			bfp_temp = bfp_idx
		endif	
		
		deallocate(vicinity_pts)
		deallocate(vicinity_pt_idx)
		
	enddo
	
	no_bfp_pts = bfp_iter
	bfp_iter = 0
	
	allocate(boundary_fluid_pts(3,no_bfp_pts))
	allocate(boundary_fluid_pts_idx(no_bfp_pts))
	
	do node=1,nodes
	
		call nn_fluid(5,xbg(node),ybg(node),zbg(node)) 
		
		!print*,no_vicinity_pts ! at no point should no_vicinity_pts go to 0
		!print*,bfp_idx
		if(node.ne.1) then
			if(bfp_temp.ne.bfp_idx) then
				bfp_iter = bfp_iter + 1
				boundary_fluid_pts_idx(bfp_iter) = bfp_idx
				boundary_fluid_pts(1,bfp_iter) = bfp(1)
				boundary_fluid_pts(2,bfp_iter) = bfp(2)
				boundary_fluid_pts(3,bfp_iter) = bfp(3)
				bfp_temp = bfp_idx
			endif
		else
			bfp_iter = bfp_iter + 1
			boundary_fluid_pts_idx(bfp_iter) = bfp_idx
			boundary_fluid_pts(1,bfp_iter) = bfp(1)
			boundary_fluid_pts(2,bfp_iter) = bfp(2)
			boundary_fluid_pts(3,bfp_iter) = bfp(3)			
			bfp_temp = bfp_idx
		endif	
		
		deallocate(vicinity_pts)
		deallocate(vicinity_pt_idx)
		
	enddo
	


END


SUBROUTINE DISCRETIZATION_I_COMP_IBM(PHI,PHID,nvars)

	use declare_variables
	implicit none
	
	integer prim,nvars,node
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	real,dimension(NImax) :: RHS
	
	
	AC_COMP_IBM(1:ptsmax) = 1.d0
	AP_COMP_IBM(1:ptsmax) = 1.d0/3.d0
	AM_COMP_IBM(1:ptsmax) = 1.d0/3.d0
	
	do prim = 1,nvars
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
		
		
		do i = 3,NI(nbl)-2
		
			if(type_ibm(i,j,k,nbl).eq.0) then ! body points
				AM_COMP_IBM(i) = 0
				AC_COMP_IBM(i) = 1
				AP_COMP_IBM(i) = 0
				RHS(i) = 0
				cycle
			elseif(type_ibm(i,j,k,nbl).eq.-1) then ! ghost points
				AM_COMP_IBM(i) = 0
				AC_COMP_IBM(i) = 1
				AP_COMP_IBM(i) = 2 
				RHS(i) = -5.d0/2.d0*PHI(i,j,k,nbl,prim) + 2.d0*PHI(i+1,j,k,nbl,prim) + 1.d0/2.d0*PHI(i+2,j,k,nbl,prim)
				cycle
			elseif(type_ibm(i,j,k,nbl).eq.1) then
				RHS(i) = (bdisc/4.d0)*(PHI(i+2,j,k,nbl,prim)-PHI(i-2,j,k,nbl,prim)) &
					+ (adisc/2.d0)*(PHI(i+1,j,k,nbl,prim)-PHI(i-1,j,k,nbl,prim))
				do node = 1,no_bfp_pts
					call get_loc_index(boundary_fluid_pts_idx(node))
					if(i.eq.i_loc.and.j.eq.j_loc.and.k.eq.k_loc.and.nbl.eq.nbl_loc) then
						AM_COMP_IBM(i) = 1.d0/4.d0
						AC_COMP_IBM(i) = 1
						AP_COMP_IBM(i) = 1.d0/4.d0
						RHS(i) = -3.d0/4.d0*PHI(i-1,j,k,nbl,prim) + 3.d0/4.d0*PHI(i+1,j,k,nbl,prim) 
						exit	
					endif
				enddo 
			
			endif
			
		Enddo
	
		AM_COMP_IBM(1) = 0
		AC_COMP_IBM(1) = 1
		AP_COMP_IBM(1) = 2 
		RHS(1) = -5.d0/2.d0*PHI(1,j,k,nbl,prim) + 2.d0*PHI(2,j,k,nbl,prim) + 1.d0/2.d0*PHI(3,j,k,nbl,prim)
		
		
		AM_COMP_IBM(2) = 1.d0/4.d0
		AC_COMP_IBM(2) = 1
		AP_COMP_IBM(2) = 1.d0/4.d0
		RHS(2) = -3.d0/4.d0*PHI(1,j,k,nbl,prim) + 3.d0/4.d0*PHI(3,j,k,nbl,prim) 
			
		
		AM_COMP_IBM(NI(nbl)-1) = 1.d0/4.d0
		AC_COMP_IBM(NI(nbl)-1) = 1
		AP_COMP_IBM(NI(nbl)-1) = 1.d0/4.d0
		RHS(NI(nbl)-1) = 3.d0/4.d0*PHI(NI(nbl),j,k,nbl,prim) - 3.d0/4.d0*PHI(NI(nbl)-2,j,k,nbl,prim) 
			
		
		AM_COMP_IBM(NI(nbl)) = 2
		AC_COMP_IBM(NI(nbl)) = 1
		AP_COMP_IBM(NI(nbl)) = 0 
		RHS(NI(nbl)) = 5.d0/2.d0*PHI(NI(nbl),j,k,nbl,prim) - 2.d0*PHI(NI(nbl)-1,j,k,nbl,prim) - 1.d0/2.d0*PHI(NI(nbl)-2,j,k,nbl,prim) 
			
		! super diagonal AP_COMP_IBM(2:NI(nbl))
		! diagonal       AC_COMP_IBM(1:NI(nbl))
		! subdiagonal    AM_COMP_IBM(1:NI(nbl)-1)
		
		call TDMA(AM_COMP_IBM(1:NI(nbl)-1),AC_COMP_IBM(1:NI(nbl)),AP_COMP_IBM(2:NI(nbl)),RHS,PHID(1:NI(nbl),j,k,nbl,prim),NI(nbl))
		
	  
	enddo
	enddo
	enddo
	enddo
	

	
END


SUBROUTINE DISCRETIZATION_J_COMP_IBM(PHI,PHID,nvars)
	use declare_variables
	implicit none
	
	integer prim,nvars,node
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
			
			if(type_ibm(i,j,k,nbl).eq.0) then ! body points
				AM_COMP_IBM(j) = 0
				AC_COMP_IBM(j) = 1
				AP_COMP_IBM(j) = 0
				RHS(j) = 0
				cycle
			elseif(type_ibm(i,j,k,nbl).eq.-1) then ! ghost points
				AM_COMP_IBM(j) = 0
				AC_COMP_IBM(j) = 1
				AP_COMP_IBM(j) = 2 
				RHS(j) = -5.d0/2.d0*PHI(i,j,k,nbl,prim) + 2.d0*PHI(i,j+1,k,nbl,prim) + 1.d0/2.d0*PHI(i,j+2,k,nbl,prim) 
				cycle
			elseif(type_ibm(i,j,k,nbl).eq.1) then
				RHS(j) = (bdisc/4.d0)*(PHI(i,j+2,k,nbl,prim)-PHI(i,j-2,k,nbl,prim)) &
					+ (adisc/2.d0)*(PHI(i,j+1,k,nbl,prim)-PHI(i,j-1,k,nbl,prim))
				do node = 1,no_bfp_pts
					call get_loc_index(boundary_fluid_pts_idx(node))
					if(i.eq.i_loc.and.j.eq.j_loc.and.k.eq.k_loc.and.nbl.eq.nbl_loc) then
						AM_COMP_IBM(j) = 1.d0/4.d0
						AC_COMP_IBM(j) = 1
						AP_COMP_IBM(j) = 1.d0/4.d0
						RHS(j) = -3.d0/4.d0*PHI(i,j-1,k,nbl,prim) + 3.d0/4.d0*PHI(i,j+1,k,nbl,prim)
						exit	
					endif
				enddo 
			
			endif
			
		Enddo
		
		AM_COMP_IBM(1) = 0
		AC_COMP_IBM(1) = 1
		AP_COMP_IBM(1) = 2 
		RHS(1) = -5.d0/2.d0*PHI(i,1,k,nbl,prim) + 2.d0*PHI(i,2,k,nbl,prim) + 1.d0/2.d0*PHI(i,3,k,nbl,prim) 
		
		
		AM_COMP_IBM(2) = 1.d0/4.d0
		AC_COMP_IBM(2) = 1
		AP_COMP_IBM(2) = 1.d0/4.d0
		RHS(2) = -3.d0/4.d0*PHI(i,1,k,nbl,prim)+ 3.d0/4.d0*PHI(i,3,k,nbl,prim) 			
		
		AM_COMP_IBM(NJ(nbl)-1) = 1.d0/4.d0
		AC_COMP_IBM(NJ(nbl)-1) = 1
		AP_COMP_IBM(NJ(nbl)-1) = 1.d0/4.d0
		RHS(NJ(nbl)-1) = 3.d0/4.d0*PHI(i,NJ(nbl),k,nbl,prim) - 3.d0/4.d0*PHI(i,NJ(nbl)-2,k,nbl,prim)
		
		AM_COMP_IBM(NJ(nbl)) = 2
		AC_COMP_IBM(NJ(nbl)) = 1
		AP_COMP_IBM(NJ(nbl)) = 0 
		RHS(NJ(nbl)) = 5.d0/2.d0*PHI(i,NJ(nbl),k,nbl,prim) - 2.d0*PHI(i,NJ(nbl)-1,k,nbl,prim) - 1.d0/2.d0*PHI(i,NJ(nbl)-2,k,nbl,prim) 
		
		! super diagonal AP_COMP_IBM(2:NI(nbl))
		! diagonal       AC_COMP_IBM(1:NI(nbl))
		! subdiagonal    AM_COMP_IBM(1:NI(nbl)-1)

		call TDMA(AM_COMP_IBM(1:NJ(nbl)-1),AC_COMP_IBM(1:NJ(nbl)),AP_COMP_IBM(2:NJ(nbl)),RHS,PHID(i,1:NJ(nbl),k,nbl,prim),NJ(nbl))
		
		
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
			! insert here
			
		enddo
		
		! insert here
		
		! super diagonal AP_COMP_IBM(2:NI(nbl))
		! diagonal       AC_COMP_IBM(1:NI(nbl))
		! subdiagonal    AM_COMP_IBM(1:NI(nbl)-1)
		
		call TDMAP(1,NK(nbl)-1,AP_COMP_IBM(1:NK(nbl)-1),AC_COMP_IBM(1:NK(nbl)-1),AM_COMP_IBM(1:NK(nbl)-1),RHS,NK(nbl)-1)
		
		PHID(i,j,1:NK(nbl),nbl,prim) = RHS(1:NK(nbl))
	 
	enddo
	enddo
	enddo
	enddo
 
END


SUBROUTINE DISC_I_COMP_GRID_IBM(PHI,PHID)
		use declare_variables
		implicit none
		
		integer prim,nvars
		real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
		real,dimension(NImax) :: RHS
		real LL
		
		AC_COMP_IBM(1:ptsmax) = 1.d0
		AP_COMP_IBM(1:ptsmax) = 1.d0/3.d0
		AM_COMP_IBM(1:ptsmax) = 1.d0/3.d0
		
		do nbl = 1,nblocks
		do k = 1,NK(nbl)
		do j = 1,NJ(nbl)
		
		do i = 3,NI(nbl)-2
		RHS(i) = (bdisc/4.d0)*(PHI(i+2,j,k,nbl)-PHI(i-2,j,k,nbl)) &
					 + (adisc/2.d0)*(PHI(i+1,j,k,nbl)-PHI(i-1,j,k,nbl))
		enddo
		
		LL = PHI(NI(nbl),j,k,nbl) - PHI(1,j,k,nbl)
		
		AM_COMP_IBM(1) = 0
		AC_COMP_IBM(1) = 1
		AP_COMP_IBM(1) = 2 
		RHS(1) = -5.d0/2.d0*PHI(1,j,k,nbl) + 2.d0*PHI(2,j,k,nbl) + 1.d0/2.d0*PHI(3,j,k,nbl)
		
		AM_COMP_IBM(2) = 1.d0/4.d0
		AC_COMP_IBM(2) = 1
		AP_COMP_IBM(2) = 1.d0/4.d0
		RHS(2) = -3.d0/4.d0*PHI(1,j,k,nbl) + 3.d0/4.d0*PHI(3,j,k,nbl) 
		
		AM_COMP_IBM(NI(nbl)-1) = 1.d0/4.d0
		AC_COMP_IBM(NI(nbl)-1) = 1
		AP_COMP_IBM(NI(nbl)-1) = 1.d0/4.d0
		RHS(NI(nbl)-1) = 3.d0/4.d0*PHI(NI(nbl),j,k,nbl) - 3.d0/4.d0*PHI(NI(nbl)-2,j,k,nbl) 
		
		AM_COMP_IBM(NI(nbl)) = 2
		AC_COMP_IBM(NI(nbl)) = 1
		AP_COMP_IBM(NI(nbl)) = 0 
		RHS(NI(nbl)) = 5.d0/2.d0*PHI(NI(nbl),j,k,nbl) - 2.d0*PHI(NI(nbl)-1,j,k,nbl) - 1.d0/2.d0*PHI(NI(nbl)-2,j,k,nbl) 
		
		call TDMA(AM_COMP_IBM(1:NI(nbl)-1),AC_COMP_IBM(1:NI(nbl)),AP_COMP_IBM(2:NI(nbl)),RHS,PHID(1:NI(nbl),j,k,nbl),NI(nbl))
		
		
		enddo
		enddo
		enddo
	  
END


SUBROUTINE DISC_J_COMP_GRID_IBM(PHI,PHID)
		use declare_variables
		implicit none
		
		integer prim,nvars
		real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
		real,dimension(NJmax) :: RHS
		real LL
		
		AC_COMP_IBM(1:ptsmax) = 1.d0
		AP_COMP_IBM(1:ptsmax) = 1.d0/3.d0
		AM_COMP_IBM(1:ptsmax) = 1.d0/3.d0
		
		do nbl = 1,nblocks
		do k = 1,NK(nbl)
		do i = 1,NI(nbl)
		
		do j = 3,Nj(nbl)-2
		RHS(j) = (bdisc/4.d0)*(PHI(i,j+2,k,nbl)-PHI(i,j-2,k,nbl)) &
					 + (adisc/2.d0)*(PHI(i,j+1,k,nbl)-PHI(i,j-1,k,nbl))
		enddo
		
		LL = PHI(i,NJ(nbl),k,nbl) - PHI(i,1,k,nbl)
		
		AM_COMP_IBM(1) = 0
		AC_COMP_IBM(1) = 1
		AP_COMP_IBM(1) = 2 
		RHS(1) = -5.d0/2.d0*PHI(i,1,k,nbl) + 2.d0*PHI(i,2,k,nbl) + 1.d0/2.d0*PHI(i,3,k,nbl)
		
		AM_COMP_IBM(2) = 1.d0/4.d0
		AC_COMP_IBM(2) = 1
		AP_COMP_IBM(2) = 1.d0/4.d0
		RHS(2) = -3.d0/4.d0*PHI(i,1,k,nbl) + 3.d0/4.d0*PHI(i,3,k,nbl) 
		
		AM_COMP_IBM(NJ(nbl)-1) = 1.d0/4.d0
		AC_COMP_IBM(NJ(nbl)-1) = 1
		AP_COMP_IBM(NJ(nbl)-1) = 1.d0/4.d0
		RHS(NJ(nbl)-1) = 3.d0/4.d0*PHI(i,NJ(nbl),k,nbl) - 3.d0/4.d0*PHI(i,NJ(nbl)-2,k,nbl) 
		
		AM_COMP_IBM(NJ(nbl)) = 2
		AC_COMP_IBM(NJ(nbl)) = 1
		AP_COMP_IBM(NJ(nbl)) = 0 
		RHS(NJ(nbl)) = 5.d0/2.d0*PHI(i,NJ(nbl),k,nbl) - 2.d0*PHI(i,NJ(nbl)-1,k,nbl) - 1.d0/2.d0*PHI(i,NJ(nbl)-2,k,nbl) 
		
		call TDMA(AM_COMP_IBM(1:NJ(nbl)-1),AC_COMP_IBM(1:NJ(nbl)),AP_COMP_IBM(2:NJ(nbl)),RHS,PHID(i,1:NJ(nbl),k,nbl),NJ(nbl))
		
		
		enddo
		enddo
		enddo
	  
END


SUBROUTINE TDMA(a_arr,b_arr,c_arr,d_arr,x_arr,N_arr) 
! verified in scrap
! validation case https://matlabgeeks.weebly.com/uploads/8/0/4/8/8048228/thomas_algorithm_and_tridiagonal_matrix-_example.pdf

	use declare_variables
	implicit none
	
	integer :: N_arr,diag
	real :: a_arr(N_arr-1),c_arr(N_arr-1)
	real :: b_arr(N_arr),d_arr(N_arr),x_arr(N_arr)
	real :: w_arr
	
	do diag = 1,N_arr-1
		w_arr = a_arr(diag)/b_arr(diag)
		b_arr(diag+1) = b_arr(diag+1) - w_arr*c_arr(diag)
		d_arr(diag+1) = d_arr(diag+1) - w_arr*d_arr(diag)
	enddo
	
	x_arr(N_arr) = d_arr(N_arr)/b_arr(N_arr)
	
	do diag = N_arr-1,1,-1
		x_arr(diag) = (d_arr(diag) - c_arr(diag)*x_arr(diag+1))/b_arr(diag)
	enddo



END

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PREPROCESSING ROUTINES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!********************* READ INPUT ***********************************************************
      SUBROUTINE READ_INPUT()
		  use declare_variables
		  implicit none
		  
		  integer iii,node,elem
		  character(len=300) dummy
		  
		  OPEN(finput, file='input.dat', form = 'formatted')
		  
		  read(finput,*) Lx,Ly,Lz
		  read(finput,*) nblocks
		  
		  Allocate(NI(nblocks))			!Size of NI along the i directoin
		  Allocate(NJ(nblocks))
		  Allocate(NK(nblocks))
		  
		  do nbl = 1,nblocks
			read(finput,*) NI(nbl),NJ(nbl), NK(nbl)
		  enddo
		  
		  read(finput,*) nsteps, nprims,nconserv
		  read(finput,*) Re,Mach,gamma,Pr,T_ref
		  read(finput,*) dscheme, dschemek   
		  read(finput,*) fscheme, alphaf   
		  READ(finput,*) rk_steps,time_step
		  READ(finput,*) nsteps_time
		  READ(finput,*) restart
		  READ(finput,*) tgv_covo
		  READ(finput,*) mesh_type
		  READ(finput,*) perI,perJ,perK
		  READ(finput,*) bscheme1, bscheme2
		  READ(finput,*) exp_comp
		  
		  close(finput)
		  
		  print*, "The length of block ,nprims and nconserv are" ,nblocks,nprims,nconserv
		  
		  
		  ! IBM
		  
		  allocate(xbg(nodes))
		  allocate(ybg(nodes))
		  allocate(zbg(nodes))
		  allocate(xbn(nodes))
		  allocate(ybn(nodes))
		  allocate(zbn(nodes))
		  
  		  allocate(connect(elements,3))
		  
		  open(fcyl,file='cylinder-normals.dat',form='formatted')

		  do iii = 1,14
			read(fcyl,*) dummy
		  enddo	

		  do node = 1,nodes
			read(fcyl,*) xbg(node),ybg(node),zbg(node),xbn(node),ybn(node),zbn(node)
		  enddo

		  do elem = 1,elements
			read(fcyl,*) (connect(elem,i),i=1,3)
		  enddo
		  
		  close(fcyl)
		  
		  
	  
      END
!********************************************************************************************

!********************* ALLOCATE_ROUTINE *****************************************************
      SUBROUTINE ALLOCATE_ROUTINE()
		  use declare_variables
		  implicit none
		
		  NImax = maxval(NI)
		  NJmax = maxval(NJ)
		  NKmax = maxval(NK)
		
		  print*, "The max value of NI NJ NK are ",NImax, NJmax, NKmax
		
		  ALLOCATE(Xgrid(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(Ygrid(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(Zgrid(NImax,NJmax,NKmax,nblocks))
		  
		  ALLOCATE(Qp(NImax,NJmax,NKmax,nblocks,nprims))
		 
		  ALLOCATE(Qpi(NImax,NJmax,NKmax,nblocks,nprims))
		  ALLOCATE(Qpj(NImax,NJmax,NKmax,nblocks,nprims))
		  ALLOCATE(Qpk(NImax,NJmax,NKmax,nblocks,nprims))
		  
		  ALLOCATE(Qc(NImax,NJmax,NKmax,nblocks,nconserv))
		  ALLOCATE(Qcnew(NImax,NJmax,NKmax,nblocks,nconserv))
		  ALLOCATE(Qcini(NImax,NJmax,NKmax,nblocks,nconserv))
		 
		 
		  ALLOCATE(mu(NImax,NJmax,NKmax,nblocks))
		
		  ALLOCATE(xgrid_i(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(zgrid_i(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(ygrid_i(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(xgrid_j(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(zgrid_j(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(ygrid_j(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(xgrid_k(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(zgrid_k(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(ygrid_k(NImax,NJmax,NKmax,nblocks))
	
		  ALLOCATE(igrid_x(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(jgrid_x(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(kgrid_x(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(igrid_y(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(jgrid_y(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(kgrid_y(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(igrid_z(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(jgrid_z(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(kgrid_z(NImax,NJmax,NKmax,nblocks))
		
		  ALLOCATE(Fflux(NImax,NJmax,NKmax,nblocks,nconserv))
		  ALLOCATE(Gflux(NImax,NJmax,NKmax,nblocks,nconserv))
		  ALLOCATE(Hflux(NImax,NJmax,NKmax,nblocks,nconserv))
		  ALLOCATE(net_flux(NImax,NJmax,NKmax,nblocks,nconserv))
		  ALLOCATE(fluxD(NImax,NJmax,NKmax,nblocks,nconserv))
		
		  ALLOCATE(Jac(NImax,NJmax,NKmax,nblocks))
		
		  ALLOCATE(fluxDf(NImax,NJmax,NKmax,nblocks,nconserv))
		  ALLOCATE(fluxDg(NImax,NJmax,NKmax,nblocks,nconserv))
		  ALLOCATE(fluxDh(NImax,NJmax,NKmax,nblocks,nconserv))
		
		  ALLOCATE(facRK(rk_steps))
		  ALLOCATE(facqini(rk_steps))
		
		  ALLOCATE(vorticity_square(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(velocity_square(NImax,NJmax,NKmax,nblocks))
		  ALLOCATE(enstrophy(nsteps))
		  ALLOCATE(TKE(nsteps))
		  ALLOCATE(residual(nsteps))
		
		  ALLOCATE(fcoeff(6))
		  ALLOCATE(fdisc1(7))
		  ALLOCATE(fdisc2(7))
		
		  ptsmax = max(NImax,NJmax,NKmax)
		  ALLOCATE(AM(ptsmax))
		  ALLOCATE(AC(ptsmax))
		  ALLOCATE(AP(ptsmax))
		  ALLOCATE(AM_COMP(ptsmax))
		  ALLOCATE(AC_COMP(ptsmax))
		  ALLOCATE(AP_COMP(ptsmax))
		
		  ALLOCATE(swirl_vel_init(NImax,nbl))
		  ALLOCATE(swirl_vel_final(NImax,nbl))
		  
		  
		  ! IBM
		  
		  allocate(num_share_elems(nodes))
		  allocate(type_ibm(NImax,NJmax,NKmax,nblocks))
		  allocate(maxshare)
		  allocate(i_loc)
		  allocate(j_loc)
		  allocate(k_loc)
		  allocate(nbl_loc)
		  allocate(no_ghost_pts)
		  allocate(global_index)
		  allocate(no_body_pts)
		  allocate(no_fluid_pts)
		  allocate(no_vicinity_pts)
		  allocate(R)
		  allocate(N_vandermonde)
		  allocate(L_N)
		  allocate(flag_nn_alloc)
		  allocate(flag_vandermonde_alloc)
		  allocate(flag_pi_alloc)
		  allocate(flag_A_alloc)
		  allocate(total_comp_pts)
		  allocate(Qp_W(NImax,NJmax,NKmax,nblocks,nprims))
		  allocate(AC_COMP_IBM(ptsmax))
		  allocate(AP_COMP_IBM(ptsmax))
		  allocate(AM_COMP_IBM(ptsmax))
		  allocate(bfp(3))
		  allocate(bfp_idx)
		  allocate(no_bfp_pts)
		  allocate(grid_pts(3,NImax*NJmax*NKmax*nblocks))
		  allocate(b_incpt(3))

      END 
!********************************************************************************************

!********************* GENERATE_GRID ********************************************************
      SUBROUTINE GENERATE_GRID()
		  use declare_variables
		  implicit none
		  
		  if(mesh_type.eq.1) then
		  
			  do nbl = 1,nblocks
			  do k = 1,NK(nblocks)
			  do j = 1,NJ(nblocks)
			  do i = 1,NI(nblocks)
			  
					Xgrid(i,j,k,nbl) = Lx*(i-1.d0)/(NI(nbl)-1.d0) - Lx/2
					Ygrid(i,j,k,nbl) = Ly*(j-1.d0)/(NJ(nbl)-1.d0) - Ly/2
			 		Zgrid(i,j,k,nbl) = Lz*(k-1.d0)/(NK(nbl)-1.d0) 
			  
			  enddo
			  enddo
			  enddo
			  enddo
	  
		  elseif(mesh_type.eq.2) then
		  
			  do nbl = 1,nblocks
			  do j = 1,NJ(nblocks)
			  do i = 1,NI(nblocks)
		  
			  if(i.ge.10.and.i.le.NI(nbl)-10.and.j.ge.10.and.j.le.NJ(nbl)-10) then
					Xgrid(i,j,1,nbl) = Lx*(i-1.d0)/(NI(nbl)-1.d0)+0.2*rand(0)*Lx/(NI(nbl)-1.d0)
					Xgrid(i,j,2,nbl) = Xgrid(i,j,1,nbl)
					Xgrid(i,j,3,nbl) = Xgrid(i,j,1,nbl)
					Ygrid(i,j,1,nbl) = Ly*(j-1.d0)/(NJ(nbl)-1.d0)+0.2*rand(0)*Ly/(NJ(nbl)-1.d0)	  
					Ygrid(i,j,2,nbl) = Ygrid(i,j,1,nbl)
					Ygrid(i,j,3,nbl) = Ygrid(i,j,1,nbl)
			  else
					Xgrid(i,j,1,nbl) = Lx*(i-1.d0)/(NI(nbl)-1.d0)
					Xgrid(i,j,2,nbl) = Xgrid(i,j,1,nbl)
					Xgrid(i,j,3,nbl) = Xgrid(i,j,1,nbl)
					Ygrid(i,j,1,nbl) = Ly*(j-1.d0)/(NJ(nbl)-1.d0)
					Ygrid(i,j,2,nbl) = Ygrid(i,j,1,nbl)
					Ygrid(i,j,3,nbl) = Ygrid(i,j,1,nbl)
			  endif
	  
			  do k = 1,NK(nbl)
					Zgrid(i,j,1,nbl) = Lz*(k-1.d0)/(NK(nbl)-1.d0)
			  enddo
		  
			  enddo
			  enddo
			  enddo
	 
		  elseif(mesh_type.eq.3) then
		  
			  do nbl = 1,nblocks
			  do k = 1,NK(nblocks)
			  do j = 1,NJ(nblocks)
			  do i = 1,NI(nblocks)
	  
					Xgrid(i,j,k,nbl) = Lx/(NI(nbl)-1.d0)*(i-1+sin(3.14/2)*sin(6*3.14*(j-1)*Ly/(NJ(nbl)-1.d0)/Ly)) !-Lx/2.d0
					Ygrid(i,j,k,nbl) = Ly/(NJ(nbl)-1.d0)*(j-1+2*sin(3.14/2)*sin(6*3.14*(i-1)*Lx/(NI(nbl)-1.d0)/Lx)) !-Ly/2.d0
					Zgrid(i,j,k,nbl) = Lz*(k-1.d0)/(NK(nbl)-1.d0)
		
			  enddo
			  enddo
			  enddo
			  enddo
		  
		  endif
	  
      END 
!********************************************************************************************

!********************* IBM PREPROCESSING ********************************************************
	  SUBROUTINE ibm_preprocessing()
	  
		  use declare_variables
		  use kdtree2_module
		  implicit none
	  
		  integer elem,node,gridpts_iter
		  
		  num_share_elems = 0
		  
		  do elem = 1,elements
		  do i = 1,3
			node = connect(elem,i)
			num_share_elems(node) = num_share_elems(node) + 1
		  enddo
		  enddo

		  maxshare = maxval(num_share_elems)
		  allocate(ind_share_elems(nodes,maxshare))
		  
		  num_share_elems = 0
		  ind_share_elems = 0

		  do elem = 1,elements
		  do i = 1,3
			node = connect(elem,i)
			num_share_elems(node) = num_share_elems(node) + 1
			ind_share_elems(node,num_share_elems(node)) = elem
		  enddo
		  enddo		
		 
		  gridpts_iter = 0
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
	  
		  call ibm_type(xbg,ybg,zbg)
		  
		  call ghost_points()
		  
		  call boundary_intercept(xbg,ybg,zbg)
	  	  
		  call matrix_calculations()
		  
		  call boundary_fluid_points()

	  END
!********************************************************************************************

!********************* GLOBAL INDEX ****************************************
SUBROUTINE get_global_index(i_local,j_local,k_local,nbl_local)

	use declare_variables

	integer flag,i_local,j_local,k_local,nbl_local
	
	flag = 0
	
	do nnbbll = 1,nblocks
	do kk = 1,NK(nnbbll)
	do jj = 1,NJ(nnbbll)
	do ii = 1,NI(nnbbll)
	
		flag = flag + 1
	
		if(i_local.eq.ii.and.j_local.eq.jj.and.k_local.eq.kk.and.nbl_local.eq.nnbbll) then
		
			global_index = flag
		
		endif
	
	enddo
	enddo
	enddo
	enddo

END
!********************************************************************************************

!********************* LOCAL INDEX ****************************************
SUBROUTINE get_loc_index(indx)

	use declare_variables

	integer flag
	
	flag = 0

	do nnbbll = 1,nblocks
	do kk = 1,NK(nnbbll)
	do jj = 1,NJ(nnbbll)
	do ii = 1,NI(nnbbll)
	
		flag = flag + 1
		if(flag.eq.indx) then
		
			i_loc = ii
			j_loc = jj
			k_loc = kk
			nbl_loc = nnbbll
		
		endif
	
	enddo
	enddo
	enddo
	enddo

END
!********************************************************************************************

!********************* INITIALIZE_NON_DIMENSIONALIZE ****************************************
      SUBROUTINE INITIALIZE_NON_DIMENSIONALIZE()
		  use declare_variables
		  implicit none
		  
		  real xl,yl,zl,Etotal,x_c,y_c,vortex_strength,r_sq
		  integer nn
		  
		  if(restart.eq.0) then
		  
			  if(tgv_covo.eq.1) then
		  
				  do nbl = 1,nblocks
				  do k = 1,NK(nbl)
				  do j = 1,NJ(nbl)
				  do i = 1,NI(nbl)
	  
					xl = xgrid(i,j,k,nbl)
					yl = ygrid(i,j,k,nbl)
					zl = zgrid(i,j,k,nbl)
		 
	 
					Qp(i,j,k,nbl,2) =  sin(xl)*cos(yl)*sin(zl)
					Qp(i,j,k,nbl,3) = -cos(xl)*sin(yl)*sin(zl)
					Qp(i,j,k,nbl,4) = 0.d0
					Qp(i,j,k,nbl,5) = 1.d0/(Mach**2.d0*gamma) +(1.d0/16.d0)*(cos(2.d0*xl)+ cos(2.d0*yl))*(cos(2.d0*zl)+2.d0)
					Qp(i,j,k,nbl,6) = 1.d0           ! might be a logical error here.......temp=1k makes no sense
					Qp(i,j,k,nbl,1) = gamma*Mach**2.d0*Qp(i,j,k,nbl,5)/Qp(i,j,k,nbl,6)
					
					Qp_W(i,j,k,nbl,1) = 0.d0 ! not sure
					Qp_W(i,j,k,nbl,2) = 0.d0
					Qp_W(i,j,k,nbl,3) = 0.d0
					Qp_W(i,j,k,nbl,4) = 0.d0
					Qp_W(i,j,k,nbl,5) = 0.d0
					Qp_W(i,j,k,nbl,6) = 0.d0 ! not sure
				
					mu(i,j,k,nbl) = 1.d0
					
					Qc(i,j,k,nbl,1) = Qp(i,j,k,nbl,1)
					Qc(i,j,k,nbl,2) = Qp(i,j,k,nbl,1)*Qp(i,j,k,nbl,2)
					Qc(i,j,k,nbl,3) = Qp(i,j,k,nbl,1)*Qp(i,j,k,nbl,3)
					Qc(i,j,k,nbl,4) = Qp(i,j,k,nbl,1)*Qp(i,j,k,nbl,4)
					
					Etotal = Qp(i,j,k,nbl,5)/(Qp(i,j,k,nbl,1)*(gamma-1.d0)) + 0.5d0*(Qp(i,j,k,nbl,2)**2+Qp(i,j,k,nbl,3)**2+Qp(i,j,k,nbl,4)**2)
				
					Qc(i,j,k,nbl,5) = Qp(i,j,k,nbl,1)*(Etotal)
		
				  enddo
				  enddo
				  enddo
				  enddo
				  
				  call phi_gp(Qp_W,Qp,nprims,Qc,nconserv)
		  
			  elseif(tgv_covo.eq.2) then
	  
				  do nbl = 1,nblocks
				  do k = 1,NK(nbl)
				  do j = 1,NJ(nbl)
				  do i = 1,NI(nbl)
	  
					xl = xgrid(i,j,k,nbl)
					yl = ygrid(i,j,k,nbl)
					zl = zgrid(i,j,k,nbl)
					x_c = Lx/2.d0
					y_c = Ly/2.d0
					r_sq = ((Xgrid(i,j,k,nbl)-x_c)**2 + (Ygrid(i,j,k,nbl)-y_c)**2)
					vortex_strength = 0.02
		
	
					Qp(i,j,k,nbl,2) =  1 - vortex_strength*(Ygrid(i,j,k,nbl)-y_c)*exp(-r_sq/2.d0)
					Qp(i,j,k,nbl,3) = vortex_strength*(Xgrid(i,j,k,nbl)-x_c)*exp(-r_sq/2.d0)
					Qp(i,j,k,nbl,4) = 0.d0
					Qp(i,j,k,nbl,5) = 1.d0/(Mach**2.d0*gamma) - 1/2.d0*vortex_strength**2*exp(-r_sq)
					Qp(i,j,k,nbl,6) = 1.d0           ! might be a logical error here.......temp=1k makes no sense
					Qp(i,j,k,nbl,1) = gamma*Mach**2.d0*Qp(i,j,k,nbl,5)/Qp(i,j,k,nbl,6)
				
					Qp_W(i,j,k,nbl,1) = 0.d0 ! not sure
					Qp_W(i,j,k,nbl,2) = 0.d0
					Qp_W(i,j,k,nbl,3) = 0.d0
					Qp_W(i,j,k,nbl,4) = 0.d0
					Qp_W(i,j,k,nbl,5) = 0.d0
					Qp_W(i,j,k,nbl,6) = 0.d0 ! not sure
				
				
					mu(i,j,k,nbl) = 1.d0
				
					Qc(i,j,k,nbl,1) = Qp(i,j,k,nbl,1)
					Qc(i,j,k,nbl,2) = Qp(i,j,k,nbl,1)*Qp(i,j,k,nbl,2)
					Qc(i,j,k,nbl,3) = Qp(i,j,k,nbl,1)*Qp(i,j,k,nbl,3)
					Qc(i,j,k,nbl,4) = Qp(i,j,k,nbl,1)*Qp(i,j,k,nbl,4)
					
					Etotal = Qp(i,j,k,nbl,5)/(Qp(i,j,k,nbl,1)*(gamma-1.d0)) + 0.5d0*(Qp(i,j,k,nbl,2)**2+Qp(i,j,k,nbl,3)**2+Qp(i,j,k,nbl,4)**2)
					
					Qc(i,j,k,nbl,5) = Qp(i,j,k,nbl,1)*(Etotal)
		
				  enddo
				  enddo
				  enddo
				  enddo
				  
				  call phi_gp(Qp_W,Qp,nprims,Qc,nconserv)
		  
			  endif
			  
		  elseif(restart.eq.1) then

				open(fflow, form = 'unformatted',file='E2F10_flow.xyz')
				read(fflow) nblocks
				read(fflow) ( NI(n), NJ(n), NK(n), n= 1,nblocks )
				do  n = 1, nblocks
				read(fflow) (((( Qp(i,j,k,n,nn), i=1,NI(n)), j=1,NJ(n)), k=1,NK(n)), nn=1,nprims),   &
		
						(((( Qc(i,j,k,n,nn), i=1,NI(n)), j=1,NJ(n)), k=1,NK(n)), nn=1,nconserv),   &
			
						!(((( Qpi(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nprims),  &
						!(((( Qpj(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nprims),  &
						!(((( Qpk(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nprims),  &
						
						(((( net_flux(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv)!,  &
						!(((( Fflux(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv),  &
						!(((( Gflux(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv),  &
						!(((( Hflux(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv)!,  &
						!
						!((( xgrid_i(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
						!((( ygrid_i(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
						!((( zgrid_i(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
						!((( xgrid_j(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
						!((( ygrid_j(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
						!((( zgrid_j(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
						!((( xgrid_k(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
						!((( ygrid_k(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
						!((( zgrid_k(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))
			  enddo
			  close(fflow)

			  endif
	  
      END 
!******************************************************************************************** 

!********************* METRICS **************************************************************
      SUBROUTINE METRICS()
	  
		  use declare_variables
		  implicit none
		  
		  real xil,yil,zil,xjl,yjl,zjl,xkl,ykl,zkl,vol
		  
		  do nbl = 1,nblocks
		  do k = 1,NK(nblocks)
		  do j = 1,NJ(nblocks)
		  do i = 1,NI(nblocks)
		  
			xil = xgrid_i(i,j,k,nbl)
			xjl = xgrid_j(i,j,k,nbl)
			xkl = xgrid_k(i,j,k,nbl)
			yil = ygrid_i(i,j,k,nbl)
			yjl = ygrid_j(i,j,k,nbl)
			ykl = ygrid_k(i,j,k,nbl)
			zil = zgrid_i(i,j,k,nbl)
			zjl = zgrid_j(i,j,k,nbl)
			zkl = zgrid_k(i,j,k,nbl)
			
			vol = xil*(yjl*zkl-ykl*zjl) - xjl*(yil*zkl-ykl*zil) + xkl*(yil*zjl-zil*yjl)
			Jac(i,j,k,nbl) = 1.d0/vol
			
			igrid_x(i,j,k,nbl) = (yjl*zkl-ykl*zjl)*Jac(i,j,k,nbl)
			igrid_y(i,j,k,nbl) = (xkl*zjl-zkl*xjl)*Jac(i,j,k,nbl)
			igrid_z(i,j,k,nbl) = (xjl*ykl-xkl*yjl)*Jac(i,j,k,nbl) 
			jgrid_x(i,j,k,nbl) = (ykl*zil-yil*zkl)*Jac(i,j,k,nbl)
			jgrid_y(i,j,k,nbl) = (xil*zkl-zil*xkl)*Jac(i,j,k,nbl) 
			jgrid_z(i,j,k,nbl) = (yil*xkl-xil*ykl)*Jac(i,j,k,nbl)
			kgrid_x(i,j,k,nbl) = (yil*zjl-yjl*zil)*Jac(i,j,k,nbl)
			kgrid_y(i,j,k,nbl) = (xjl*zil-xil*zjl)*Jac(i,j,k,nbl)
			kgrid_z(i,j,k,nbl) = (yjl*xil-yil*xjl)*Jac(i,j,k,nbl)
		
			!Deallocate(xgrid_i,xgrid_j,xgrid_k,ygrid_i,ygrid_j,ygrid_k,zgrid_i,zgrid_j,zgrid_k)
	  
		  enddo
		  enddo
		  enddo
		  enddo

      END 
!********************************************************************************************

!********************* DISCRETIZATION_FILTER_RK_VALS_VALS ***********************************
      SUBROUTINE DISCRETIZATION_FILTER_RK_VALS()
		  use declare_variables
		  implicit none
	  
		  ! Setting Discretization Coeff I,J directions
		  if(dscheme.eq.1) then           !E2
			adisc=1.d0
			bdisc=0.d0
			alpha=0.d0
		  elseif(dscheme.eq.2) then       !E4
			adisc=4.d0/3.d0
			bdisc=-1.d0/3.d0
			alpha=0.d0
		  elseif(dscheme.eq.3) then       !C4
			adisc=3.d0/2.d0
			bdisc=0.d0
			alpha=1.d0/4.d0
		  elseif(dscheme.eq.4) then       !C6
			adisc=14.d0/9.d0
			bdisc=1.d0/9.d0
			alpha=1.d0/3.d0
		  endif
	  
		  !************** this has been implemented in j direction, which should be correct, code it in i and k
		  ! setting discretization coefficients for boundary points
		  ! 1,NI
		  if(dscheme.eq.1.or.dscheme.eq.2) then
			 if(bscheme1.eq.2) then   !E2
				fdisc1(1) = -1.5d0
				fdisc1(2) = 2.d0
				fdisc1(3) = -0.5d0
				fdisc1(4) = 0.d0
				fdisc1(5) = 0.d0
				fdisc1(6) = 0.d0
				fdisc1(7) = 0.d0
				alpha1 = 0.d0
			 elseif(bscheme1.eq.4) then !E4
				fdisc1(1) = -25.d0/12.d0
				fdisc1(2) = 4.d0
				fdisc1(3) = -3.d0
				fdisc1(4) = 4.d0/3.d0
				fdisc1(5) = -1.d0/4.d0
				fdisc1(6) = 0.d0
				fdisc1(7) = 0.d0
				alpha1 = 0.d0
			 endif
		  elseif(dscheme.eq.3.or.dscheme.eq.4) then
			 if(bscheme1.eq.2) then     !C2
				fdisc1(1) = -2.d0
				fdisc1(2) = 2d0
				fdisc1(3) = 0.d0
				fdisc1(4) = 0.d0
				fdisc1(5) = 0.d0
				fdisc1(6) = 0.d0
				fdisc1(7) = 0.d0
				alpha1 = 1.d0
			 elseif(bscheme1.eq.4) then    !C4
				fdisc1(1) = -17.d0/6.d0
				fdisc1(2) = 3.d0/2.d0
				fdisc1(3) = 3.d0/2.d0
				fdisc1(4) = -1.d0/6.d0
				fdisc1(5) = 0.d0
				fdisc1(6) = 0.d0
				fdisc1(7) = 0.d0
				alpha1 = 3.d0
			 elseif(bscheme1.eq.6) then  !C6 !sort the coefficients
				fdisc1(1) = -197.d0/60.d0
				fdisc1(2) = -5.d0/12.d0
				fdisc1(3) = 5.d0
				fdisc1(4) = -5.d0/3.d0
				fdisc1(5) = 5.d0/12.d0
				fdisc1(6) = 0.d0
				fdisc1(7) = 0.d0
				alpha1 = 5.d0
			 endif
		  endif
		  
	  
		  !2,NI-1 !sort this out ! fdisc2????
		  if(dscheme.eq.1.or.dscheme.eq.2) then
			  if(bscheme2.eq.2) then   !E2
				fdisc2(1) = -1.5d0
				fdisc2(2) = 2.d0
				fdisc2(3) = -0.5d0
				fdisc2(4) = 0.d0
				fdisc2(5) = 0.d0
				fdisc2(6) = 0.d0
				fdisc2(7) = 0.d0
				alpha2 = 0.d0
			  elseif(bscheme2.eq.4) then !E4
				fdisc2(1) = -25.d0/12.d0
				fdisc2(2) = 4.d0
				fdisc2(3) = -3.d0
				fdisc2(4) = 4.d0/3.d0
				fdisc2(5) = -1.d0/4.d0
				fdisc2(6) = 0.d0
				fdisc2(7) = 0.d0
				alpha2 = 0.d0
			  endif
		  elseif(dscheme.eq.3.or.dscheme.eq.4) then
		  ! C2 is not required for boundary points 2,NI-1. Interior schemes can be used at the boundary
			  if(bscheme2.eq.2) then     !C2
				!fdisc1(1) = -2.d0
				!fdisc1(2) = 2d0
				!fdisc1(3) = 0.d0
				!fdisc1(4) = 0.d0
				!fdisc1(5) = 0.d0
				!fdisc1(6) = 0.d0
				!fdisc1(7) = 0.d0
				!alpha1 = 1.d0
			  elseif(bscheme2.eq.4) then    !C4
				fdisc2(1) = -17.d0/6.d0
				fdisc2(2) = 3.d0/2.d0
				fdisc2(3) = 3.d0/2.d0
				fdisc2(4) = -1.d0/6.d0
				fdisc2(5) = 0.d0
				fdisc2(6) = 0.d0
				fdisc2(7) = 0.d0
				alpha2 = 3.d0
			  elseif(bscheme2.eq.6) then  !C6 !sort the coefficients
				fdisc2(1) = -197.d0/60.d0
				fdisc2(2) = -5.d0/12.d0
				fdisc2(3) = 5.d0
				fdisc2(4) = -5.d0/3.d0
				fdisc2(5) = 5.d0/12.d0
				fdisc2(6) = 0.d0
				fdisc2(7) = 0.d0
				alpha2 = 5.d0
			  endif
		  endif
	
	  ! for compact schemes for the boundary disc, we will have different AC_comp, AP_comp, AM_comp in different direction
	  
		  AC_COMP(1:ptsmax) = 1.d0
		  AP_COMP(1:ptsmax) = alpha
		  AM_COMP(1:ptsmax) = alpha
	  
		  ! Setting Discretization Coeff K direction
		  if(dschemek.eq.1) then           !E2
				adisck=1.d0
				bdisck=0.d0
				alphak=0.d0
		  elseif(dschemek.eq.2) then       !E4
				adisck=4.d0/3.d0
				bdisck=-1.d0/3.d0
				alphak=0.d0
		  elseif(dschemek.eq.3) then       !C4
				adisck=3.d0/2.d0
				bdisck=0.d0
				alphak=1.d0/4.d0
		  elseif(dschemek.eq.4) then       !C6
				adisck=14.d0/9.d0
				bdisck=1.d0/9.d0
				alphak=1.d0/3.d0
		  endif
	  
		  AC_COMP(1:ptsmax) = 1.d0
		  AP_COMP(1:ptsmax) = alpha
		  AM_COMP(1:ptsmax) = alpha
	  
		  !!!!!!!!!!!!!!!!!!!!!!1 RK Coeff !!!!!!!!!!!!!!!!!!!!11
		  facRK(1) = 1.d0/6.d0
		  facRK(2) = 2.d0/6.d0
		  facRK(3) = 2.d0/6.d0
		  facRK(4) = 1.d0/6.d0
		  
		  facqini(1) = 0.d0
		  facqini(2) = 0.5d0
		  facqini(3) = 0.5d0
		  facqini(4) = 1.d0
	  
		  !***************** Filter Coeff *************!
		  if(fscheme.eq.2) then          !F2  
				fcoeff(1) = 0.5 + alphaf
				fcoeff(2) = 0.5 + alphaf
				fcoeff(3) = 0.d0
				fcoeff(4) = 0.d0
				fcoeff(5) = 0.d0
				fcoeff(6) = 0.d0
		  elseif(fscheme.eq.4) then          !F4
				fcoeff(1) = 5.d0/8.d0 + 3.d0*alphaf/4.d0
				fcoeff(2) = 0.5d0 + alphaf
				fcoeff(3) = -1.d0/8.d0 + alphaf/4.d0 
				fcoeff(4) = 0.d0
				fcoeff(5) = 0.d0
				fcoeff(6) = 0.d0
		  elseif(fscheme.eq.6) then          !F6
				fcoeff(1) = 11.d0/16.d0 + 5.d0*alphaf/8.d0
				fcoeff(2) = 15.d0/32.d0 + 17.d0*alphaf/16.d0
				fcoeff(3) = -3.d0/16.d0 + 3.d0*alphaf/8.d0
				fcoeff(4) = 1.d0/32.d0 - alphaf/16.d0
				fcoeff(5) = 0.d0
				fcoeff(6) = 0.d0
		  elseif(fscheme.eq.8) then           !F8
				fcoeff(1) = 93.d0/128.d0 + 70.d0*alphaf/128.d0
				fcoeff(2) = 7.d0/16.d0 + 18.d0*alphaf/16.d0
				fcoeff(3) = -7.d0/32.d0 + 14.d0*alphaf/32.d0
				fcoeff(4) = 1.d0/16.d0 - alphaf/8.d0
				fcoeff(5) = -1.d0/128.d0 + alphaf/64.d0
				fcoeff(6) = 0.d0
		  elseif(fscheme.eq.10) then         !F10
				fcoeff(1) = 193.d0/256.d0 + 126.d0*alphaf/256.d0
				fcoeff(2) = 105.d0/256.d0 + 302.d0*alphaf/256.d0
				fcoeff(3) = -15.d0/64.d0 + 30.d0*alphaf/64.d0
				fcoeff(4) = 45.d0/512.d0 - 90.d0*alphaf/512.d0
				fcoeff(5) = -5.d0/256.d0 + 10.d0*alphaf/256.d0
				fcoeff(6) = 1.d0/512.d0 - 2.d0*alphaf/512.d0
		  endif
		
		  AM(1:ptsmax) = alphaf
		  AP(1:ptsmax) = alphaf
		  AC(1:ptsmax) = 1.d0

      END 
!********************************************************************************************

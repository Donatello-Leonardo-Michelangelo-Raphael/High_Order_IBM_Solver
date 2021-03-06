      MODULE declare_variables
	  
	  use kdtree2_module
	  implicit none
	  
	  
      integer, parameter :: fresidual = 1, finput = 2, fgrid = 3, fflow = 4, fenstrophy_tke = 5, &
							fswirl_error = 6
	  integer  nblocks, NImax, NJmax, NKmax
	  integer, allocatable :: NI(:), NJ(:), NK(:)
	  integer  nbl, i, j,k,n,ii,jj,kk,nnbbll
	  integer iter,nsteps, nprims, nconserv
	  
	  ! Discretization scheme coefficients
	  integer dscheme, dschemek          ! Discretization scheme in (i,j) and k->if in the k direction we don't have enough number of grid points to calculate the scheme on
	  real adisc, bdisc, alpha, adisck, bdisck, alphak
	  real alphaf, alpha1, alpha2
	  integer fscheme, ptsmax             ! Filtering scheme in i,j,k
	  
	  real time, time_step, Re, Mach, gamma, Pr, T_ref
	  real Lx,Ly,Lz
	  
	  integer rk_steps, nsteps_time
	  real timestep
	  
	  integer restart, tgv_covo, mesh_type, exp_comp
	  integer perI, perJ, perK
	  integer bscheme1, bscheme2        ! boundar scheme
	  
	  real swirl_error
	  
	  
	  real,allocatable :: xgrid(:,:,:,:),ygrid(:,:,:,:), zgrid(:,:,:,:)
	  
	  real,allocatable :: xgrid_i(:,:,:,:),ygrid_i(:,:,:,:),zgrid_i(:,:,:,:)
	  real,allocatable :: xgrid_j(:,:,:,:),ygrid_j(:,:,:,:),zgrid_j(:,:,:,:)
	  real,allocatable :: xgrid_k(:,:,:,:),ygrid_k(:,:,:,:),zgrid_k(:,:,:,:)
	  
	  real,allocatable :: igrid_x(:,:,:,:),igrid_y(:,:,:,:),igrid_z(:,:,:,:)
	  real,allocatable :: jgrid_x(:,:,:,:),jgrid_y(:,:,:,:),jgrid_z(:,:,:,:)
	  real,allocatable :: kgrid_x(:,:,:,:),kgrid_y(:,:,:,:),kgrid_z(:,:,:,:)
	  real,allocatable :: Jac(:,:,:,:)
	  
	  real,allocatable :: Qp(:,:,:,:,:) ! Primitive variables Rho,u,v,w, P,T Its 5 dimensional. 3 dimensions for the directions, one for block no, and one for the variable of interest
	  real,allocatable :: Qpi(:,:,:,:,:)
	  real,allocatable :: Qpj(:,:,:,:,:)
	  real,allocatable :: Qpk(:,:,:,:,:)
	  real,allocatable :: Qc(:,:,:,:,:),Qcnew(:,:,:,:,:),Qcini(:,:,:,:,:) ! Conservative variable Rhu,Rhv, RhW,RhE and Next Time Step
	  real,allocatable :: mu(:,:,:,:)   ! Non Dimensional Viscosity based on Sutherland's Law
	  
	  real,allocatable :: Fflux(:,:,:,:,:), Gflux(:,:,:,:,:), Hflux(:,:,:,:,:), net_flux(:,:,:,:,:), fluxD(:,:,:,:,:)  ! fluxD is the derivative obtained after calling the subroutine in respective dirn. 
	  real,allocatable :: fluxDf(:,:,:,:,:),fluxDg(:,:,:,:,:),fluxDh(:,:,:,:,:)
	  real,allocatable :: facQini(:),facRK(:)
	  
	  real,allocatable :: vorticity_square(:,:,:,:)
	  real,allocatable :: velocity_square(:,:,:,:)
	  real,allocatable :: enstrophy(:)
	  real,allocatable :: TKE(:)
	  real,allocatable :: fcoeff(:)                            ! Coefficients for filtering
	  
	  real,allocatable :: AM(:),AC(:),AP(:),AM_COMP(:),AC_COMP(:),AP_COMP(:)
	  
	  real,allocatable :: residual(:)
	  
	  real,allocatable :: swirl_vel_init(:,:)
	  real,allocatable :: swirl_vel_final(:,:)
	  
	  real,allocatable :: fdisc1(:),fdisc2(:)
	  
	  ! IBM 
	  
	  integer, parameter :: dp=kind(0.d0)
	  integer,parameter :: fcyl=10,nodes=1077,elements=1436,n_dim=3 ! declare_variables
	  integer,allocatable :: num_share_elems(:),ind_share_elems(:,:) ! declare_variables
	  real,allocatable :: xbg(:),ybg(:),zbg(:),xbn(:),ybn(:),zbn(:) ! declare_variables
	  real,allocatable :: connect(:,:) ! declare_variables
	  integer,allocatable :: type_ibm(:,:,:,:)
	  integer,allocatable :: maxshare,i_loc,j_loc,k_loc,nbl_loc,global_index
	  real,allocatable :: ghost_pt(:,:)
	  integer,allocatable :: ghost_pt_idx(:),vicinity_pt_idx(:)
	  integer,allocatable :: no_ghost_pts,no_body_pts,no_fluid_pts,no_vicinity_pts
	  real,allocatable :: BII(:,:)
	  real,allocatable :: vicinity_pts(:,:),R
	  real(dp),allocatable :: W(:,:),V(:,:),p_i(:,:),A_matrix(:,:,:)
	  integer,allocatable :: flag_nn_alloc,N_vandermonde,L_N,flag_vandermonde_alloc
	  integer,allocatable :: flag_pi_alloc,flag_A_alloc,total_comp_pts
	  integer,allocatable :: dim_A(:,:)
	  real,allocatable :: Qp_W(:)
	  real,allocatable :: AC_COMP_IBM(:),AP_COMP_IBM(:),AM_COMP_IBM(:)
	  real,allocatable :: bfp(:),boundary_fluid_pts(:,:)
	  integer,allocatable :: bfp_idx,boundary_fluid_pts_idx(:),no_bfp_pts
      real,allocatable :: grid_pts(:,:),b_incpt(:)
	  type(kdtree2),pointer :: tree_grid_pts
	  
	  
	  
	  
      END MODULE

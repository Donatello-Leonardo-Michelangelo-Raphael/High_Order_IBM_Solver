!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ POST PROCESSING ROUTINES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!******************************* OUTPUT GRID & FLOW *******************************************
      SUBROUTINE OUTPUT()
	  
		  use declare_variables
		  implicit none
		  
		  integer nn, iterl
		  
		  open(fresidual, form = 'formatted', file = 'E2F10_residual_covo.txt')
		  do iterl = 1,nsteps
			 write(fresidual,*) residual(iterl)
		  enddo
		  close(fresidual)
		  
		  open(fswirl_error, form = 'formatted', file = 'E2F10_swirlerror_covo.txt')
		  write(fswirl_error,*) swirl_error
		  do nbl = 1,nblocks
		 	 write(fswirl_error,*) (swirl_vel_final(i,nbl), i=1,NI(nbl))
		  enddo
		  close(fswirl_error)
		  
		  open(fgrid, form = 'unformatted', file = 'E2F10_grid_covo.xyz')
		  write(fgrid) nblocks
		  write(fgrid) (NI(nbl), NJ(nbl), NK(nbl), nbl = 1, nblocks)
		  
		  do nbl= 1, nblocks
		 	 write(fgrid) ((( xgrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
		 				 ((( ygrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), &
		 				 ((( zgrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))
		  enddo
		  close(fgrid)
	  	  	  
		  open(fflow, form = 'unformatted',file='E2F10_flow_covo.xyz')
		  write(fflow) nblocks
		  write(fflow) ( NI(nbl), NJ(nbl), NK(nbl), 2*nconserv+nprims, nbl= 1,nblocks )
		  do nbl = 1, nblocks
		  write(fflow) (((( Qp(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nprims),   &
		  
		 			 (((( Qc(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv),   &
		  
		 			 !(((( Qpi(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nprims),  &
		 			 !(((( Qpj(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nprims),  &
		 			 !(((( Qpk(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nprims),  &
		 			 
		 			 (((( net_flux(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv)!,  &
		 			 !(((( Fflux(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv),  &
		 			 !(((( Gflux(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv),  &
		 			 !(((( Hflux(i,j,k,nbl,nn), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), nn=1,nconserv)!,  &
		 			 
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
	  	 
      END
	  

!********************************* DEALLOCATE *************************************************	  
	  
      SUBROUTINE DEALLOCATE_ROUTINE()
		  
		  use declare_variables
		  implicit none
		  
		  deallocate(Xgrid,Ygrid,Zgrid)
	  END
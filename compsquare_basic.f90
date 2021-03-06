!***************************************************************************************
!****************** WRITTEN BY DR. NAGABHUSHANA RAO VADLAMANI **************************
!***************BASED ON THE HIGH ORDER COMPSQUARE SOLVER DEVELOPED BY DR. NRV**********
!***************DISTRIBUTED AS A PART OF AS6041 COURSE ON ADVANCED CFD *****************
!***************************************************************************************

      PROGRAM MISSION_AS6041

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
!@@@@@@@@@@@@@ PRE PROCESSING STEPS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	

!****** Declare Arrays, Integers, Real numbers *********************************	
      use declare_variables
	  integer step
	  real :: start_time, stop_time
	  call cpu_time(start_time)
	  
      print*, 'Declared Variables'

!******** READ INPUT ***********************************************************
      print*, 'Reading input file....'
      call READ_INPUT()
      print*, 'Reading input file Done'		  
!****** INITIALIZE ARRAYS ******************************************************
      call ALLOCATE_ROUTINE()
      print*, 'Allocated Integers, Real numbers and Arrays Done'
!****** GENERATE GRID ******************************************************
      call GENERATE_GRID()
      print*, 'Grid Generation Done'
!****** IBM PREPROCESSING ******************************************************
	  call ibm_preprocessing()
	  print*, 'IBM Preprocessing Done'
!******** Initialize and Non-dimensionalize arrays *****************************
      print*, 'Initializing and non-dimensionalizing arrays....'
      call INITIALIZE_NON_DIMENSIONALIZE()
      print*, 'Initializing and non-dimensionalizing arrays done'
!******** Read Discretization coeffs for flow !*********************************
      print*, 'Obtaining Discretization coefficients for flow....'
      call DISCRETIZATION_FILTER_RK_VALS()
      print*, 'Obtaining Discretization coefficients for flow Done'          
!******** Read Filter coeffs for flow !*****************************************		  
!******** Compute Metric terms *************************************************
      print*, 'Computing Metrics....'
	  
	  if(exp_comp.eq.1) then
		
		 call DISC_I_EXP_GRID(xgrid,xgrid_i)
		 call DISC_I_EXP_GRID(ygrid,ygrid_i)
		 call DISC_I_EXP_GRID(zgrid,zgrid_i)
					
		 call DISC_J_EXP_GRID(xgrid,xgrid_j)
		 call DISC_J_EXP_GRID(ygrid,ygrid_j)
		 call DISC_J_EXP_GRID(zgrid,zgrid_j)
		 	 
		 if(tgv_covo.eq.1) then
			 
			 call DISC_K_EXP_GRID(xgrid,xgrid_k)
			 call DISC_K_EXP_GRID(ygrid,ygrid_k)
			 call DISC_K_EXP_GRID(zgrid,zgrid_k)
			 
		 elseif(tgv_covo.eq.2) then
	  
			 call DISC_K2D_EXP_GRID(xgrid,xgrid_k)
			 call DISC_K2D_EXP_GRID(ygrid,ygrid_k)
			 call DISC_K2D_EXP_GRID(zgrid,zgrid_k)
	  
		 endif
	  
	  elseif(exp_comp.eq.2) then
	  
		 call DISC_I_COMP_GRID(xgrid,xgrid_i)
		 call DISC_I_COMP_GRID(ygrid,ygrid_i)
		 call DISC_I_COMP_GRID(zgrid,zgrid_i)
				  
		 call DISC_J_COMP_GRID(xgrid,xgrid_j)
		 call DISC_J_COMP_GRID(ygrid,ygrid_j)
		 call DISC_J_COMP_GRID(zgrid,zgrid_j)
				 
		 if(tgv_covo.eq.1) then
				 
		   	 call DISC_K_COMP_GRID(xgrid,xgrid_k)
			 call DISC_K_COMP_GRID(ygrid,ygrid_k)
			 call DISC_K_COMP_GRID(zgrid,zgrid_k)
	  
		 elseif(tgv_covo.eq.2) then
	  
			 call DISC_K2D_EXP_GRID(xgrid,xgrid_k)
			 call DISC_K2D_EXP_GRID(ygrid,ygrid_k)
			 call DISC_K2D_EXP_GRID(zgrid,zgrid_k)
	  
		 endif
	  
	  endif
	  
      call METRICS()
      print*, 'Computing Metrics Done'    

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ END OF PRE PROCESSING STEPS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ MAIN TIME LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      print*,''
	  
	  print*,'Qp Initial',maxval(Qp(:,:,:,:,1)),maxval(Qp(:,:,:,:,2)),maxval(Qc(:,:,:,:,3)) &
			,maxval(Qp(:,:,:,:,4)),maxval(Qp(:,:,:,:,5)),maxval(Qp(:,:,:,:,6))		
	  print*,'Qc Initial',maxval(Qc(:,:,:,:,1)),maxval(Qc(:,:,:,:,2)),maxval(Qc(:,:,:,:,3)) &
			,maxval(Qc(:,:,:,:,4)),maxval(Qc(:,:,:,:,5)) ,'NA'

	 
	  do iter = 1,nsteps   !EXPLICIT TIME STEPPING
		Qcini = Qc
		Qcnew = Qc
	  
		!if(tgv_covo.eq.2) then
		!	if(iter.eq.1) then
		!		do nbl = 1,nblocks
		!		do i = 1,NI(nbl)
		!			swirl_vel_init(i,nbl) = Qp(i,NJ(nbl)/2,2,nbl,3)
		!		enddo
		!		enddo
		!	endif
		!endif
		print*,''
		print*, "Time Step:", iter
	  
		do step = 1,rk_steps
			call UNSTEADY(step,iter)
			!if(step.eq.4) then
			!	call FILTERING_I(Qc,nconserv)
			!	call FILTERING_J(Qc,nconserv)
			!	if(tgv_covo.eq.1) then
			!		call FILTERING_K(Qc,nconserv)
			!	endif
			!endif
			call SET_PRIMITIVES()
		enddo
	  
		!if(tgv_covo.eq.2) then
		!	if(iter.eq.nsteps) then
		!		do nbl = 1,nblocks
		!		do i = 1,NI(nbl)
		!			swirl_vel_final(i,nbl) = Qp(i,NJ(nbl)/2,2,nbl,3)
		!		enddo
		!		enddo
		!	endif
		!endif
		
		!if(tgv_covo.eq.1) then
		!	call ENSTROPHY_TKE(iter)
		!	print*, 'Enstrophy : ', enstrophy(iter), 'TKE : ', tke(iter)
		!endif
		
		residual(iter) = maxval(abs(Qc-Qcini))
		print*, "Max Error:", residual(iter)
	!  	print*, "Maximum Flux:", maxval(net_flux)
	!	print*, "Max Net Flux", maxval(net_flux)
		print*,'Qp      ',maxval(Qp(:,:,:,:,1)),maxval(Qp(:,:,:,:,2)),maxval(Qc(:,:,:,:,3)) &
			,maxval(Qp(:,:,:,:,4)),maxval(Qp(:,:,:,:,5)),maxval(Qp(:,:,:,:,6))
		print*,'Qc      ',maxval(Qc(:,:,:,:,1)),maxval(Qc(:,:,:,:,2)),maxval(Qc(:,:,:,:,3)) &
			,maxval(Qc(:,:,:,:,4)),maxval(Qc(:,:,:,:,5)),'NA'
		print*,'Net Flux',maxval(net_flux(:,:,:,:,1)),maxval(net_flux(:,:,:,:,2)),maxval(net_flux(:,:,:,:,3)) &
			,maxval(net_flux(:,:,:,:,4)),maxval(net_flux(:,:,:,:,5)),'NA'
	  
		time = time + time_step
	  
	  enddo	
	  
	  print*, 'Time Integration Done'
	  
	 !if(tgv_covo.eq.2) then
	 !do nbl = 1,nblocks
	 !	swirl_error = maxval(abs(swirl_vel_final(:,nbl)-swirl_vel_init(:,nbl)))
	 !enddo
	 !endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!!@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN TIME LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!!@@@@@@@@@@@@@@@@@@@@@@@ POST PROCESSING STEPS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      
!	!****** Write the output ***************
!101	  print*, 'Writing Output....'
!      call SET_PRIMITIVES()
!      call OUTPUT(0)
!      print*, 'Writing output done'
! 
!	!****** Deallocate Arrays **************
!      call DEALLOCATE_ROUTINE()
!      print*, 'Deallocated Arrays'
!
!	  call cpu_time(stop_time)
!	  print*, "Run Time :", stop_time-start_time, "Seconds"
!	  print*, "Run Time :", (stop_time-start_time)/60, "Minutes"
!
      END PROGRAM
!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!!@@@@@@@@@@@@@@@@@@@@@@@ END OF POST PROCESSING STEPS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

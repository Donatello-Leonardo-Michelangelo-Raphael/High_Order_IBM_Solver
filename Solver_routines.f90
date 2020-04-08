!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SOLVER ROUTINES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!******************************************************************************
!************************** UNSTEADY ******************************************
!******************************************************************************
      SUBROUTINE UNSTEADY(stepl,iterl)
		  use declare_variables
		  implicit none
		  
		  real rhl,ul,vl,wl,pl,Tl,El
		  real Ucont,Vcont,Wcont,ixl,jxl,kxl,iyl,jyl,kyl,izl,jzl,kzl,vol
		  real uil,ujl,ukl,vil,vjl,vkl,wil,wjl,wkl,Til,Tjl,Tkl,mul
		  real Txx,Tyy,Tzz,Txy,Tyz,Tzx
		  real bx,by,bz
		  real div2b3,facPrM,net_fluxD,sterm
		  real T_x,T_y,T_z,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z
		  integer stepl,var,iterl
	  
		  !********** Inviscid Flux Estimation **********!
	  
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)	
		  do j = 1,NJ(nbl)
		  do i = 1,NI(nbl)
	  
			rhl = Qp(i,j,k,nbl,1)
			ul = Qp(i,j,k,nbl,2)
			vl = Qp(i,j,k,nbl,3)
			wl = Qp(i,j,k,nbl,4)
			pl = Qp(i,j,k,nbl,5)
			Tl = Qp(i,j,k,nbl,6)
			El = pl/(rhl*(gamma-1.d0)) + 0.5*(ul**2+vl**2+wl**2)
	  
			ixl = igrid_x(i,j,k,nbl)
			jxl = jgrid_x(i,j,k,nbl)
			kxl = kgrid_x(i,j,k,nbl)
			iyl = igrid_y(i,j,k,nbl)
			jyl = jgrid_y(i,j,k,nbl)
			kyl = kgrid_y(i,j,k,nbl)
			izl = igrid_z(i,j,k,nbl)
			jzl = jgrid_z(i,j,k,nbl)
			kzl = kgrid_z(i,j,k,nbl)
			
			Ucont = ixl*ul + iyl*vl + izl*wl
			Vcont = jxl*ul + jyl*vl + jzl*wl
			Wcont = kxl*ul + kyl*vl + kzl*wl
	  
			vol = 1.d0/Jac(i,j,k,nbl)
	  
			Fflux(i,j,k,nbl,1) = (rhl*Ucont)*vol
			Fflux(i,j,k,nbl,2) = (rhl*ul*Ucont + ixl*pl)*vol
			Fflux(i,j,k,nbl,3) = (rhl*vl*Ucont + iyl*pl)*vol
			Fflux(i,j,k,nbl,4) = (rhl*wl*Ucont + izl*pl)*vol
			Fflux(i,j,k,nbl,5) = (rhl*El*Ucont + pl*Ucont)*vol
							
			Gflux(i,j,k,nbl,1) = (rhl*Vcont)*vol
			Gflux(i,j,k,nbl,2) = (rhl*ul*Vcont + jxl*pl)*vol
			Gflux(i,j,k,nbl,3) = (rhl*vl*Vcont + jyl*pl)*vol
			Gflux(i,j,k,nbl,4) = (rhl*wl*Vcont + jzl*pl)*vol
			Gflux(i,j,k,nbl,5) = (rhl*El*Vcont + pl*Vcont)*vol
			
			Hflux(i,j,k,nbl,1) = (rhl*Wcont)*vol
			Hflux(i,j,k,nbl,2) = (rhl*ul*Wcont + kxl*pl)*vol
			Hflux(i,j,k,nbl,3) = (rhl*vl*Wcont + kyl*pl)*vol
			Hflux(i,j,k,nbl,4) = (rhl*wl*Wcont + kzl*pl)*vol
			Hflux(i,j,k,nbl,5) = (rhl*El*Wcont + pl*Wcont)*vol
	  
		  enddo
		  enddo
		  enddo
		  enddo
	  
		  !********** Viscous Flux Etimation *********!
	  
		  if(tgv_covo.eq.1) then
		  
			  if(exp_comp.eq.1) then
			  
				call DISCRETIZATION_I_EXP(Qp,Qpi,nprims)
				call DISCRETIZATION_J_EXP(Qp,Qpj,nprims)
				call DISCRETIZATION_K_EXP(Qp,Qpk,nprims)
				  
			  elseif(exp_comp.eq.2) then
			  
				call DISCRETIZATION_I_COMP(Qp,Qpi,nprims)
				call DISCRETIZATION_J_COMP(Qp,Qpj,nprims)
				call DISCRETIZATION_K_COMP(Qp,Qpk,nprims)
				  
			  endif
		  
			  do nbl = 1,nblocks
			  do k = 1,NK(nbl)	
			  do j = 1,NJ(nbl)
			  do i = 1,NI(nbl)
	  
				ixl = igrid_x(i,j,k,nbl)
				jxl = jgrid_x(i,j,k,nbl)
				kxl = kgrid_x(i,j,k,nbl)
				iyl = igrid_y(i,j,k,nbl)
				jyl = jgrid_y(i,j,k,nbl)
				kyl = kgrid_y(i,j,k,nbl)
				izl = igrid_z(i,j,k,nbl)
				jzl = jgrid_z(i,j,k,nbl)
				kzl = kgrid_z(i,j,k,nbl)
				
				rhl = Qp(i,j,k,nbl,1)
				ul = Qp(i,j,k,nbl,2)
				vl = Qp(i,j,k,nbl,3)
				wl = Qp(i,j,k,nbl,4)
				pl = Qp(i,j,k,nbl,5)
				Tl = Qp(i,j,k,nbl,6)
				El = pl/(rhl*(gamma-1.d0)) + 0.5*(ul**2+vl**2+wl**2)
				
				uil = qpi(i,j,k,nbl,2)
				vil = qpi(i,j,k,nbl,3)
				wil = qpi(i,j,k,nbl,4)
				Til = qpi(i,j,k,nbl,6)
				
				ujl = qpj(i,j,k,nbl,2)
				vjl = qpj(i,j,k,nbl,3)
				wjl = qpj(i,j,k,nbl,4)
				Tjl = qpj(i,j,k,nbl,6)
				
				ukl = qpk(i,j,k,nbl,2)
				vkl = qpk(i,j,k,nbl,3)
				wkl = qpk(i,j,k,nbl,4)
				Tkl = qpk(i,j,k,nbl,6)
			
				Sterm = 110.4/T_ref
				mu(i,j,k,nbl) = Tl**1.5*(1.d0+Sterm)/(Tl+Sterm)
				mul = mu(i,j,k,nbl)
				vol = 1.d0/Jac(i,j,k,nbl)
				
				
				u_x = ixl*uil + jxl*ujl + kxl*ukl
				v_x = ixl*vil + jxl*vjl + kxl*vkl
				w_x = ixl*wil + jxl*wjl + kxl*wkl
				T_x = ixl*Til + jxl*Tjl + kxl*Tkl
				
				u_y = iyl*uil + jyl*ujl + kyl*ukl
				v_y = iyl*vil + jyl*vjl + kyl*vkl
				w_y = iyl*wil + jyl*wjl + kyl*wkl
				T_y = iyl*Til + jyl*Tjl + kyl*Tkl
				
				u_z = izl*uil + jzl*ujl + kzl*ukl
				v_z = izl*vil + jzl*vjl + kzl*vkl
				w_z = izl*wil + jzl*wjl + kzl*wkl
				T_z = izl*Til + jzl*Tjl + kzl*Tkl
				
				div2b3 = (-2.d0/3.d0)*(u_x + v_y + w_z)
				facPrM = mul/((gamma-1.d0)*Pr*Re*Mach**2)
				
				Txx = (mul/Re)*(2.d0*u_x + div2b3)
				Tyy = (mul/Re)*(2.d0*v_y + div2b3)
				Tzz = (mul/Re)*(2.d0*w_z + div2b3)
				Txy = (mul/Re)*(u_y + v_x)
				Tyz = (mul/Re)*(v_z + w_y)
				Tzx = (mul/Re)*(w_x + u_z)
				
				bx = ul*Txx + vl*Txy + wl*Tzx + facPrM*T_x
				by = ul*Txy + vl*Tyy + wl*Tyz + facPrM*T_y
				bz = ul*Tzx + vl*Tyz + wl*Tzz + facPrM*T_z
			
			
				Fflux(i,j,k,nbl,1) = - Fflux(i,j,k,nbl,1) 
				Fflux(i,j,k,nbl,2) = - Fflux(i,j,k,nbl,2) + (ixl*Txx + iyl*Txy + izl*Tzx )*vol
				Fflux(i,j,k,nbl,3) = - Fflux(i,j,k,nbl,3) + (ixl*Txy + iyl*Tyy + izl*Tyz )*vol
				Fflux(i,j,k,nbl,4) = - Fflux(i,j,k,nbl,4) + (ixl*Tzx + iyl*Tyz + izl*Tzz )*vol
				Fflux(i,j,k,nbl,5) = - Fflux(i,j,k,nbl,5) + (ixl*bx  + iyl*by  + izl*bz  )*vol
				
				Gflux(i,j,k,nbl,1) = - Gflux(i,j,k,nbl,1)   
				Gflux(i,j,k,nbl,2) = - Gflux(i,j,k,nbl,2) + (jxl*Txx + jyl*Txy + jzl*Tzx )*vol
				Gflux(i,j,k,nbl,3) = - Gflux(i,j,k,nbl,3) + (jxl*Txy + jyl*Tyy + jzl*Tyz )*vol
				Gflux(i,j,k,nbl,4) = - Gflux(i,j,k,nbl,4) + (jxl*Tzx + jyl*Tyz + jzl*Tzz )*vol
				Gflux(i,j,k,nbl,5) = - Gflux(i,j,k,nbl,5) + (jxl*bx  + jyl*by  + jzl*bz  )*vol
				
				Hflux(i,j,k,nbl,1) = - Hflux(i,j,k,nbl,1) 
				Hflux(i,j,k,nbl,2) = - Hflux(i,j,k,nbl,2) + (kxl*Txx + kyl*Txy + kzl*Tzx )*vol
				Hflux(i,j,k,nbl,3) = - Hflux(i,j,k,nbl,3) + (kxl*Txy + kyl*Tyy + kzl*Tyz )*vol
				Hflux(i,j,k,nbl,4) = - Hflux(i,j,k,nbl,4) + (kxl*Tzx + kyl*Tyz + kzl*Tzz )*vol
				Hflux(i,j,k,nbl,5) = - Hflux(i,j,k,nbl,5) + (kxl*bx  + kyl*by  + kzl*bz  )*vol
			
			  enddo
			  enddo
			  enddo
			  enddo
	  
		  endif
	  
		  
		  if(exp_comp.eq.1) then
			  net_flux = 0.d0
			  call DISCRETIZATION_I_EXP(Fflux,fluxD,nconserv)
			  net_flux = net_flux + fluxD
			  fluxDf = fluxD
			  call DISCRETIZATION_J_EXP(Gflux,fluxD,nconserv)
			  net_flux = net_flux + fluxD
			  fluxDg = fluxD
			  if(tgv_covo.eq.1) then
			 		 call DISCRETIZATION_K_EXP(Hflux,fluxD,nconserv)
			 		 net_flux = net_flux + fluxD
			 		 fluxDh = fluxD
			  elseif(tgv_covo.eq.2) then
			 		 call DISCRETIZATION_K2D_EXP(Hflux,fluxD,nconserv)
			 		 net_flux = net_flux + fluxD
			 		 fluxDh = fluxD
			  endif
		  elseif(exp_comp.eq.2) then
			  net_flux = 0.d0
			  call DISCRETIZATION_I_COMP(Fflux,fluxD,nconserv)
			  net_flux = net_flux + fluxD
			  fluxDf = fluxD
			  call DISCRETIZATION_J_COMP(Gflux,fluxD,nconserv)
			  net_flux = net_flux + fluxD
			  fluxDg = fluxD
			  if(tgv_covo.eq.1) then
			 		 call DISCRETIZATION_K_COMP(Hflux,fluxD,nconserv)
			 		 net_flux = net_flux + fluxD
			 		 fluxDh = fluxD
			  elseif(tgv_covo.eq.2) then
			 		 call DISCRETIZATION_K2D_EXP(Hflux,fluxD,nconserv)
			 		 net_flux = net_flux + fluxD
			 		 fluxDh = fluxD
			  endif
		  endif
		  
	  
		  do var = 1,nconserv
		  do nbl = 1,nblocks 
		  do k = 1,NK(nbl)
		  do j = 1,NJ(nbl)
		  do i = 1,NI(nbl)
	  
			Qcnew(i,j,k,nbl,var) = Qcnew(i,j,k,nbl,var) + time_step*net_flux(i,j,k,nbl,var)*facRK(stepl)*Jac(i,j,k,nbl)
	  
			if(stepl.lt.4) then
				Qc(i,j,k,nbl,var) = Qcini(i,j,k,nbl,var)+ time_step*net_flux(i,j,k,nbl,var)*facqini(stepl+1)*Jac(i,j,k,nbl)
			else
				Qc(i,j,k,nbl,var)= Qcnew(i,j,k,nbl,var)
			endif
		  
		  enddo
		  enddo
		  enddo
		  enddo
		  enddo
	  
	  END

!******************************************************************************
!************************ END UNSTEADY ****************************************
!******************************************************************************


!******************************************************************************
!************************ ENSTROPHY TKE ***************************************
!******************************************************************************
	  
	  SUBROUTINE ENSTROPHY_TKE(iterl)
	  
		  use declare_variables
		  implicit none
		  
		  real ixl,jxl,kxl,iyl,jyl,kyl,izl,jzl,kzl,rhl
		  real uil,vil,wil,ujl,vjl,wjl,ukl,vkl,wkl
		  real u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,u1,v1,w1
		  integer iterl
		  real,dimension(NJmax,NKmax,nblocks) :: enstrophy_plane,TKE_plane
		  real,dimension(NKmax,nblocks) :: enstrophy_line,TKE_line
	 
		  if(exp_comp.eq.1) then
			  call DISCRETIZATION_I_EXP(Qp,Qpi,nprims)
			  call DISCRETIZATION_J_EXP(Qp,Qpj,nprims)
			  call DISCRETIZATION_K_EXP(Qp,Qpk,nprims)
		  elseif(exp_comp.eq.2) then
			  call DISCRETIZATION_I_COMP(Qp,Qpi,nprims)
			  call DISCRETIZATION_J_COMP(Qp,Qpj,nprims)
			  call DISCRETIZATION_K_COMP(Qp,Qpk,nprims)
		  endif
		  
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do j = 1,NJ(nbl)
		  do i = 1,NI(nbl)
		  
			u1 = Qp(i,j,k,nbl,2)
			v1 = Qp(i,j,k,nbl,3)
			w1 = Qp(i,j,k,nbl,4)
		
			ixl = igrid_x(i,j,k,nbl)
			jxl = jgrid_x(i,j,k,nbl)
			kxl = kgrid_x(i,j,k,nbl)
			iyl = igrid_y(i,j,k,nbl)
			jyl = jgrid_y(i,j,k,nbl)
			kyl = kgrid_y(i,j,k,nbl)
			izl = igrid_z(i,j,k,nbl)
			jzl = jgrid_z(i,j,k,nbl)
			kzl = kgrid_z(i,j,k,nbl)
		  
			rhl = Qp(i,j,k,nbl,1)
			
			uil = qpi(i,j,k,nbl,2)
			vil = qpi(i,j,k,nbl,3)
			wil = qpi(i,j,k,nbl,4)
			
			ujl = qpj(i,j,k,nbl,2)
			vjl = qpj(i,j,k,nbl,3)
			wjl = qpj(i,j,k,nbl,4)
			
			ukl = qpk(i,j,k,nbl,2)
			vkl = qpk(i,j,k,nbl,3)
			wkl = qpk(i,j,k,nbl,4)
		  
			u_x = ixl*uil + jxl*ujl + kxl*ukl
			v_x = ixl*vil + jxl*vjl + kxl*vkl
			w_x = ixl*wil + jxl*wjl + kxl*wkl
			
			u_y = iyl*uil + jyl*ujl + kyl*ukl
			v_y = iyl*vil + jyl*vjl + kyl*vkl
			w_y = iyl*wil + jyl*wjl + kyl*wkl
			
			u_z = izl*uil + jzl*ujl + kzl*ukl
			v_z = izl*vil + jzl*vjl + kzl*vkl
			w_z = izl*wil + jzl*wjl + kzl*wkl
			
			vorticity_square(i,j,k,nbl) = ((w_y-v_z)**2 + (w_x-u_z)**2 + (v_x-u_y)**2)*rhl
			velocity_square(i,j,k,nbl) = (u1**2 + v1**2 + w1**2)*rhl
		  
		  enddo
		  enddo
		  enddo
		  enddo
	 

		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do j = 1,NJ(nbl)
		  enstrophy_plane(j,k,nbl) = 0
		  TKE_plane(j,k,nbl) = 0.d0
		  do i = 1,NI(nbl)-1
		  
			enstrophy_plane(j,k,nbl) = enstrophy_plane(j,k,nbl) &
			+ 0.5*(Xgrid(i+1,j,k,nbl) - Xgrid(i,j,k,nbl))*(vorticity_square(i+1,j,k,nbl) + vorticity_square(i,j,k,nbl))
			
			TKE_plane(j,k,nbl) = TKE_plane(j,k,nbl) &
			+ 0.5*(Xgrid(i+1,j,k,nbl) - Xgrid(i,j,k,nbl))*(velocity_square(i+1,j,k,nbl) + velocity_square(i,j,k,nbl))
			
		  enddo
		  enddo
		  enddo
		  enddo
	 
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  enstrophy_line(k,nbl) = 0.d0
		  TKE_line(k,nbl) = 0.d0
		  do j = 1,NJ(nbl)-1
		  
			enstrophy_line(k,nbl) = enstrophy_line(k,nbl) &
			+ 0.5*(Ygrid(1,j+1,k,nbl) - Ygrid(1,j,k,nbl))*(enstrophy_plane(j+1,k,nbl) + enstrophy_plane(j,k,nbl))
			
			TKE_line(k,nbl) = TKE_line(k,nbl) &
			+ 0.5*(Ygrid(1,j+1,k,nbl) - Ygrid(1,j,k,nbl))*(TKE_plane(j+1,k,nbl) + TKE_plane(j,k,nbl))
		  
		  enddo
		  enddo
		  enddo
			
		  enstrophy(iterl) = 0.d0
		  TKE(iterl) = 0.d0
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)-1
		  
			enstrophy(iterl) = enstrophy(iterl) &
			+ 0.5*(Zgrid(1,1,k+1,nbl) - Zgrid(1,1,k,nbl))*(enstrophy_line(k+1,nbl) + enstrophy_line(k,nbl))
			
			TKE(iterl) = TKE(iterl) &
			+ 0.5*(Zgrid(1,1,k+1,nbl) - Zgrid(1,1,k,nbl))*(TKE_line(k+1,nbl) + TKE_line(k,nbl))
		  
		  enddo
		  enddo
		  
		  enstrophy(iterl) = enstrophy(iterl)/(2*Lx*Ly*Lz)
		  TKE(iterl) = TKE(iterl)/(2*Lx*Ly*Lz)
		  
		  !open(fenstrophy, file = 'enstrophy.txt', form = 'formatted')
		  !write(fenstrophy,*) enstrophy(iterl)
		  !close(fenstrophy)
	 
	  END
	 
!******************************************************************************
!************************ END ENSTROPHY TKE ***********************************
!******************************************************************************	  


!******************************************************************************
!************************ SET PRIMITIVES **************************************
!******************************************************************************	  
	 
	  SUBROUTINE SET_PRIMITIVES()
		  use declare_variables
		  IMPLICIT NONE
		  
		  real rhl,ul,vl,wl,pl,Tl
	  
		  do nbl = 1,nblocks                         
		  do k = 1,NK(nblocks)                       
		  do j = 1,NJ(nblocks)                       
		  do i = 1,NI(nblocks)  
		  
				rhl= Qc(i,j,k,nbl,1)
				ul= Qc(i,j,k,nbl,2)/Qc(i,j,k,nbl,1)
				vl= Qc(i,j,k,nbl,3)/Qc(i,j,k,nbl,1)
				wl= Qc(i,j,k,nbl,4)/Qc(i,j,k,nbl,1)
				pl= (Qc(i,j,k,nbl,5) - 0.5d0*Qc(i,j,k,nbl,1)*(ul**2 + vl**2 + wl**2))*(gamma-1)
				Tl= Mach**2*gamma*pl/rhl
		
				Qp(i,j,k,nbl,1) = rhl
				Qp(i,j,k,nbl,2) = ul
				Qp(i,j,k,nbl,3) = vl
				Qp(i,j,k,nbl,4) = wl
				Qp(i,j,k,nbl,5) = pl
				Qp(i,j,k,nbl,6) = Tl
		
		  enddo
		  enddo
		  enddo
		  enddo

      END
	  
!******************************************************************************
!************************ END SET PRIMITIVES **********************************
!******************************************************************************	  
	
	
!******************************************************************************
!************************ DESCREIZATION I EXPLICIT ****************************
!******************************************************************************	  
	  
      SUBROUTINE DISCRETIZATION_I_EXP(PHI,PHID,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars,coeff
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do j = 1,NJ(nbl)
		  do i = 3,NI(nbl)-2
		  
			PHID(i,j,k,nbl,prim) = (bdisc/4.d0)*(PHI(i+2,j,k,nbl,prim)-PHI(i-2,j,k,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(i+1,j,k,nbl,prim)-PHI(i-1,j,k,nbl,prim))
		  
		  enddo
		  enddo
		  enddo
		  enddo
		  enddo
		  
		  if(perI.eq.1) then
		  
			do prim = 1,nvars
			do nbl = 1,nblocks
			do k = 1,NK(nbl)
			do j = 1,NJ(nbl)
			
				PHID(1,j,k,nbl,prim) = (bdisc/4.d0)*(PHI(3,j,k,nbl,prim)-PHI(NI(nbl)-2,j,k,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(2,j,k,nbl,prim)-PHI(NI(nbl)-1,j,k,nbl,prim))
				PHID(2,j,k,nbl,prim) = (bdisc/4.d0)*(PHI(4,j,k,nbl,prim)-PHI(NI(nbl)-1,j,k,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(3,j,k,nbl,prim)-PHI(NI(nbl),j,k,nbl,prim))
				PHID(NI(nbl)-1,j,k,nbl,prim) = (bdisc/4.d0)*(PHI(2,j,k,nbl,prim)-PHI(NI(nbl)-3,j,k,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(NI(nbl),j,k,nbl,prim)-PHI(NI(nbl)-2,j,k,nbl,prim))
				PHID(NI(nbl),j,k,nbl,prim) = (bdisc/4.d0)*(PHI(3,j,k,nbl,prim)-PHI(NI(nbl)-2,j,k,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(2,j,k,nbl,prim)-PHI(NI(nbl)-1,j,k,nbl,prim))
				
			enddo
			enddo
			enddo
			enddo
			
		  else    ! boundary case
	  
			do prim = 1,nvars
			do nbl = 1,nblocks
			do k = 1,NK(nbl)
			do j = 1,NJ(nbl)
			
				PHID(1,j,k,nbl,prim) = 0.d0
				PHID(NI(nbl),j,k,nbl,prim) = 0.d0
				PHID(2,j,k,nbl,prim) = 0.d0
				PHID(NI(nbl)-1,j,k,nbl,prim) = 0.d0
				Do coeff = 1,7
				PHID(1,j,k,nbl,prim) = PHID(1,j,k,nbl,prim) + fdisc1(coeff)*PHI(coeff,j,k,nbl,prim)
				PHID(2,j,k,nbl,prim) = PHID(2,j,k,nbl,prim) + fdisc2(coeff)*PHI(coeff,j,k,nbl,prim)
				PHID(NI(nbl),j,k,nbl,prim) = PHID(NI(nbl),j,k,nbl,prim) - fdisc1(coeff)*PHI(NJ(nbl)-coeff+1,j,k,nbl,prim)
				PHID(NI(nbl)-1,j,k,nbl,prim) = PHID(NI(nbl)-1,j,k,nbl,prim) - fdisc2(coeff)*PHI(NI(nbl)-coeff+1,j,k,nbl,prim)
				Enddo
			 
			enddo
			enddo
			enddo
			enddo
	  
	      endif 
		  
	  
      END	  
	  
!******************************************************************************
!******************** END DESCREIZATION I EXPLICIT ****************************
!******************************************************************************


!******************************************************************************
!************************ DESCREIZATION J EXPLICIT ****************************
!******************************************************************************
	  
	  SUBROUTINE DISCRETIZATION_J_EXP(PHI,PHID,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars,coeff
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do j = 3,NJ(nbl)-2
		  do i = 1,NI(nbl)
		  
			PHID(i,j,k,nbl,prim) = (bdisc/4.d0)*(PHI(i,j+2,k,nbl,prim)-PHI(i,j-2,k,nbl,prim)) &
									 + (adisc/2.d0)*(PHI(i,j+1,k,nbl,prim)-PHI(i,j-1,k,nbl,prim))
		  
		  enddo
		  enddo
		  enddo
		  enddo
		  enddo
		  
		  if(perJ.eq.1) then
				  
			do prim = 1,nvars
			do nbl = 1,nblocks
			do k = 1,NK(nbl)
			do i = 1,NI(nbl)
				
				PHID(i,1,k,nbl,prim) = (bdisc/4.d0)*(PHI(i,3,k,nbl,prim)-PHI(i,NJ(nbl)-2,k,nbl,prim)) &
											+ (adisc/2.d0)*(PHI(i,2,k,nbl,prim)-PHI(i,NJ(nbl)-1,k,nbl,prim))
				PHID(i,2,k,nbl,prim) = (bdisc/4.d0)*(PHI(i,4,k,nbl,prim)-PHI(i,NJ(nbl)-1,k,nbl,prim)) &
											+ (adisc/2.d0)*(PHI(i,3,k,nbl,prim)-PHI(i,NJ(nbl),k,nbl,prim))
				PHID(i,NJ(nbl)-1,k,nbl,prim) = (bdisc/4.d0)*(PHI(i,2,k,nbl,prim)-PHI(i,NJ(nbl)-3,k,nbl,prim)) &
											+ (adisc/2.d0)*(PHI(i,NJ(nbl),k,nbl,prim)-PHI(i,NJ(nbl)-2,k,nbl,prim))
				PHID(i,NJ(nbl),k,nbl,prim) = (bdisc/4.d0)*(PHI(i,3,k,nbl,prim)-PHI(i,NJ(nbl)-2,k,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(i,2,k,nbl,prim)-PHI(i,NJ(nbl)-1,k,nbl,prim))
			
			enddo
			enddo
			enddo
			enddo
			
		  else    ! boundary case non periodic
	  
			do prim = 1,nvars
			do nbl = 1,nblocks
			do k = 1,NK(nbl)
			do i = 1,NI(nbl)
			
				PHID(i,1,k,nbl,prim) = 0.d0
				PHID(i,NJ(nbl),k,nbl,prim) = 0.d0
				PHID(i,2,k,nbl,prim) = 0.d0
				PHID(i,NJ(nbl)-1,k,nbl,prim) = 0.d0
				Do coeff = 1,7
				PHID(i,1,k,nbl,prim) = PHID(i,1,k,nbl,prim) + fdisc1(coeff)*PHI(i,coeff,k,nbl,prim)
				PHID(i,2,k,nbl,prim) = PHID(i,2,k,nbl,prim) + fdisc2(coeff)*PHI(i,coeff,k,nbl,prim)
				PHID(i,NJ(nbl),k,nbl,prim) = PHID(i,NJ(nbl),k,nbl,prim) - fdisc1(coeff)*PHI(i,NJ(nbl)-coeff+1,k,nbl,prim)
				PHID(i,NJ(nbl)-1,k,nbl,prim) = PHID(i,NJ(nbl)-1,k,nbl,prim) - fdisc2(coeff)*PHI(i,NJ(nbl)-coeff+1,k,nbl,prim)
				Enddo
			 
			enddo
			enddo
			enddo
			enddo
	  
	  endif
	  
	  
	  END	  
	  
!******************************************************************************
!******************** END DESCREIZATION J EXPLICIT ****************************
!******************************************************************************


!******************************************************************************
!************************ DESCREIZATION K EXPLICIT ****************************
!******************************************************************************	  
	  
	  SUBROUTINE DISCRETIZATION_K_EXP(PHI,PHID,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars,coeff
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do k = 3,NK(nbl)-2
		  do j = 1,NJ(nbl)
		  do i = 1,NI(nbl)
		  
			PHID(i,j,k,nbl,prim) = (bdisc/4.d0)*(PHI(i,j,k+2,nbl,prim)-PHI(i,j,k-2,nbl,prim)) &
									+ (adisc/2.d0)*(PHI(i,j,k+1,nbl,prim)-PHI(i,j,k-1,nbl,prim))
		  
		  enddo
		  enddo
		  enddo
		  enddo
		  enddo
	  	  
		  if(perK.eq.1) then
		  
			do prim = 1,nvars
			do nbl = 1,nblocks
			do j = 1,NJ(nbl)
			do i = 1,NI(nbl)
			
				PHID(i,j,1,nbl,prim) = (bdisc/4.d0)*(PHI(i,j,3,nbl,prim)-PHI(i,j,NK(nbl)-2,nbl,prim)) & 
										+ (adisc/2.d0)*(PHI(i,j,2,nbl,prim)-PHI(i,j,NK(nbl)-1,nbl,prim))
				PHID(i,j,2,nbl,prim) = (bdisc/4.d0)*(PHI(i,j,4,nbl,prim)-PHI(i,j,NK(nbl)-1,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(i,j,3,nbl,prim)-PHI(i,j,NK(nbl),nbl,prim))
				PHID(i,j,NK(nbl)-1,nbl,prim) = (bdisc/4.d0)*(PHI(i,j,2,nbl,prim)-PHI(i,j,NK(nbl)-3,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(i,j,NK(nbl),nbl,prim)-PHI(i,j,NK(nbl)-2,nbl,prim))
				PHID(i,j,NK(nbl),nbl,prim) = (bdisc/4.d0)*(PHI(i,j,3,nbl,prim)-PHI(i,j,NK(nbl)-2,nbl,prim)) &
										+ (adisc/2.d0)*(PHI(i,j,2,nbl,prim)-PHI(i,j,NK(nbl)-1,nbl,prim))
				
			enddo
			enddo
			enddo
			enddo
			
		  else    ! boundary case
	  
			do prim = 1,nvars
			do nbl = 1,nblocks
			do j = 1,NJ(nbl)
			do i = 1,NI(nbl)
			
				PHID(i,j,1,nbl,prim) = 0.d0
				PHID(i,j,NK(nbl),nbl,prim) = 0.d0
				PHID(i,j,2,nbl,prim) = 0.d0
				PHID(i,j,NK(nbl)-1,nbl,prim) = 0.d0
				Do coeff = 1,7
				PHID(i,j,1,nbl,prim) = PHID(i,j,1,nbl,prim) + fdisc1(coeff)*PHI(i,j,coeff,nbl,prim)
				PHID(i,j,2,nbl,prim) = PHID(i,j,2,nbl,prim) + fdisc2(coeff)*PHI(i,j,coeff,nbl,prim)
				PHID(i,j,NK(nbl),nbl,prim) = PHID(i,j,NK(nbl),nbl,prim) - fdisc1(coeff)*PHI(i,j,NK(nbl)-coeff+1,nbl,prim)
				PHID(i,j,NK(nbl)-1,nbl,prim) = PHID(i,j,NK(nbl)-1,nbl,prim) - fdisc2(coeff)*PHI(i,j,NK(nbl)-coeff+1,nbl,prim)
				Enddo
			 
			enddo
			enddo
			enddo
			enddo
	  
	      endif
	  
	  END	  
	  
!******************************************************************************
!******************** END DESCREIZATION K EXPLICIT ****************************
!******************************************************************************
	  
	  
!******************************************************************************
!************************ DESCREIZATION K2D EXPLICIT **************************
!******************************************************************************	

	  SUBROUTINE DISCRETIZATION_K2D_EXP(PHI,PHID,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do j = 1,NJ(nblocks)
		  do i = 1,NI(nblocks)
		  
		  
		  PHID(i,j,1,nbl,prim) =  (bdisc/4)*(PHI(i,j,3,nbl,prim)-PHI(i,j,1,nbl,prim)) &
		 						 + (adisc/2)*(PHI(i,j,2,nbl,prim)-PHI(i,j,2,nbl,prim))
		 						 
		  PHID(i,j,2,nbl,prim) =  (bdisc/4)*(PHI(i,j,2,nbl,prim)-PHI(i,j,2,nbl,prim)) &
								  + (adisc/2)*(PHI(i,j,3,nbl,prim)-PHI(i,j,1,nbl,prim)) ! might bean error here
							
		  PHID(i,j,3,nbl,prim) =  (bdisc/4)*(PHI(i,j,3,nbl,prim)-PHI(i,j,1,nbl,prim)) &
		 					 + (adisc/2)*(PHI(i,j,2,nbl,prim)-PHI(i,j,2,nbl,prim))
							
		  enddo
		  enddo
		  enddo
		  enddo
      
	  END	

!******************************************************************************
!******************** END DESCREIZATION K2D EXPLICIT **************************
!******************************************************************************


!******************************************************************************
!************************ DESCREIZATION I COMPACT *****************************
!******************************************************************************	

	  SUBROUTINE DISCRETIZATION_I_COMP(PHI,PHID,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars,coeff
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
		  real,dimension(NImax) :: RHS
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do j = 1,NJ(nbl)
		  
		  
			Do i = 3,NI(nbl)-2
			RHS(i) = (bdisc/4.d0)*(PHI(i+2,j,k,nbl,prim)-PHI(i-2,j,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i+1,j,k,nbl,prim)-PHI(i-1,j,k,nbl,prim))
			Enddo
						
			RHS(1) = (bdisc/4.d0)*(PHI(3,j,k,nbl,prim)-PHI(NI(nbl)-2,j,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(2,j,k,nbl,prim)-PHI(NI(nbl)-1,j,k,nbl,prim))
			RHS(2) = (bdisc/4.d0)*(PHI(4,j,k,nbl,prim)-PHI(NI(nbl)-1,j,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(3,j,k,nbl,prim)-PHI(1,j,k,nbl,prim))
			RHS(NI(nbl)-2) = (bdisc/4.d0)*(PHI(1,j,k,nbl,prim)-PHI(NI(nbl)-4,j,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(NI(nbl)-1,j,k,nbl,prim)-PHI(NI(nbl)-3,j,k,nbl,prim))
			RHS(NI(nbl)-1) = (bdisc/4.d0)*(PHI(2,j,k,nbl,prim)-PHI(NI(nbl)-3,j,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(1,j,k,nbl,prim)-PHI(NI(nbl)-2,j,k,nbl,prim))
						
			call TDMAP(1,NI(nbl)-1,AP_COMP,AC_COMP,AM_COMP,RHS,NI(nbl)-1)
			
			PHID(1:NI(nbl)-1,j,k,nbl,prim) = RHS(1:NI(nbl)-1)
			PHID(NI(nbl),j,k,nbl,prim) = PHID(1,j,k,nbl,prim)
		  
		  enddo
		  enddo
		  enddo
		  enddo
	  
	  END
	  
!******************************************************************************
!******************** END DESCREIZATION I COMPACT *****************************
!******************************************************************************
        
		
!******************************************************************************
!************************ DESCREIZATION J COMPACT *****************************
!******************************************************************************	

	  SUBROUTINE DISCRETIZATION_J_COMP(PHI,PHID,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
		  real,dimension(NJmax) :: RHS
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do i = 1,NI(nbl)
		  
		  
			Do j = 3,NJ(nbl)-2
			RHS(j) = (bdisc/4.d0)*(PHI(i,j+2,k,nbl,prim)-PHI(i,j-2,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,j+1,k,nbl,prim)-PHI(i,j-1,k,nbl,prim))
			Enddo
		
			RHS(1) = (bdisc/4.d0)*(PHI(i,3,k,nbl,prim)-PHI(i,NJ(nbl)-2,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,2,k,nbl,prim)-PHI(i,NJ(nbl)-1,k,nbl,prim))
			RHS(2) = (bdisc/4.d0)*(PHI(i,4,k,nbl,prim)-PHI(i,NJ(nbl)-1,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,3,k,nbl,prim)-PHI(i,1,k,nbl,prim))
			RHS(NJ(nbl)-2) = (bdisc/4.d0)*(PHI(i,1,k,nbl,prim)-PHI(i,NJ(nbl)-4,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,NJ(nbl)-1,k,nbl,prim)-PHI(i,NJ(nbl)-3,k,nbl,prim))
			RHS(NJ(nbl)-1) = (bdisc/4.d0)*(PHI(i,2,k,nbl,prim)-PHI(i,NJ(nbl)-3,k,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,1,k,nbl,prim)-PHI(i,NJ(nbl)-2,k,nbl,prim))
							
			call TDMAP(1,NJ(nbl)-1,AP_COMP,AC_COMP,AM_COMP,RHS,NJ(nbl)-1)
			
			PHID(i,1:NJ(nbl)-1,k,nbl,prim) = RHS(1:NJ(nbl)-1)
			PHID(i,NJ(nbl),k,nbl,prim) = PHID(i,1,k,nbl,prim)
		  
		  enddo
		  enddo
		  enddo
		  enddo
	  
	  END
	  
!******************************************************************************
!******************** END DESCREIZATION J COMPACT *****************************
!******************************************************************************


!******************************************************************************
!************************ DESCREIZATION K COMPACT *****************************
!******************************************************************************	
	  
	  SUBROUTINE DISCRETIZATION_K_COMP(PHI,PHID,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
		  real,dimension(NKmax) :: RHS
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do j = 1,NJ(nbl)
		  do i = 1,Ni(nbl)
		  
			do k = 3,NK(nbl)-2
			RHS(k) = (bdisc/4.d0)*(PHI(i,j,k+2,nbl,prim)-PHI(i,j,k-2,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,j,k+1,nbl,prim)-PHI(i,j,k-1,nbl,prim))
			enddo
		
			RHS(1) = (bdisc/4.d0)*(PHI(i,j,3,nbl,prim)-PHI(i,j,NK(nbl)-2,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,j,2,nbl,prim)-PHI(i,j,NK(nbl)-1,nbl,prim))
			RHS(2) = (bdisc/4.d0)*(PHI(i,j,4,nbl,prim)-PHI(i,j,NK(nbl)-1,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,j,3,nbl,prim)-PHI(i,j,1,nbl,prim))
			RHS(NK(nbl)-2) = (bdisc/4.d0)*(PHI(i,j,1,nbl,prim)-PHI(i,j,NK(nbl)-4,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,j,NK(nbl)-1,nbl,prim)-PHI(i,j,NK(nbl)-3,nbl,prim))
			RHS(NK(nbl)-1) = (bdisc/4.d0)*(PHI(i,j,2,nbl,prim)-PHI(i,j,NK(nbl)-3,nbl,prim)) &
							+ (adisc/2.d0)*(PHI(i,j,1,nbl,prim)-PHI(i,j,NK(nbl)-2,nbl,prim))
							
			call TDMAP(1,NK(nbl)-1,AP_COMP,AC_COMP,AM_COMP,RHS,NK(nbl)-1)
			
			PHID(i,j,1:NK(nbl)-1,nbl,prim) = RHS(1:NK(nbl)-1)
			PHID(i,j,NK(nbl),nbl,prim) = PHID(i,j,1,nbl,prim)
		  
		  enddo
		  enddo
		  enddo
		  enddo
	  
	  END
	  
!******************************************************************************
!******************** END DESCREIZATION K COMPACT *****************************
!******************************************************************************
 

!******************************************************************************
!************************ DESCREIZATION I EXPLICIT GRID ***********************
!******************************************************************************	
	  
	  SUBROUTINE DISC_I_EXP_GRID(PHI,PHID)
            use declare_variables
            implicit none

            REAL,DIMENSION(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
			real LL
			integer coeff

            do nbl = 1,nblocks
            do k = 1,NK(nblocks)
            do j = 1,NJ(nblocks)
            do i = 3,NI(nblocks)-2

                PHID(i,j,k,nbl)  = (bdisc/4)*(PHI(i+2,j,k,nbl) -PHI(i-2,j,k,nbl)) &
                                    + (adisc/2)*(PHI(i+1,j,k,nbl) - PHI(i-1,j,k,nbl))
            enddo 
            enddo
            enddo
            enddo

			if(perI.eq.1) then

				  do nbl = 1,nblocks
				  do k = 1,NK(nblocks)
				  do j = 1,NJ(nblocks)
				  
				  LL = PHI(NI(nbl),j,k,nbl) -PHI(1,j,k,nbl)
				 	 PHID(1,j,k,nbl)  = (bdisc/4)*(PHI(3,j,k,nbl) -PHI(NI(nbl)-2,j,k,nbl)+LL) &
				 							 + (adisc/2)*(PHI(2,j,k,nbl) - PHI(NI(nbl)-1,j,k,nbl)+LL)
				  
				 	 PHID(2,j,k,nbl)  = (bdisc/4)*(PHI(4,j,k,nbl) -PHI(NI(nbl)-1,j,k,nbl)+LL) &
				 							 + (adisc/2)*(PHI(3,j,k,nbl) - PHI(1,j,k,nbl)) ! might be an error here
				  
				 	 PHID(NI(nbl)-1,j,k,nbl)  = (bdisc/4)*(PHI(2,j,k,nbl) -PHI(NI(nbl)-3,j,k,nbl)+LL) &
				 								 + (adisc/2)*(PHI(NI(nbl),j,k,nbl) - PHI(NI(nbl)-2,j,k,nbl))
				  
				 	 PHID(NI(nbl),j,k,nbl)    =  (bdisc/4)*(PHI(3,j,k,nbl) -PHI(NI(nbl)-2,j,k,nbl)+LL) &
				 							 + (adisc/2)*(PHI(2,j,k,nbl) - PHI(NI(nbl)-1,j,k,nbl)+LL)
				 	 !print*,PHID(NI(nblocks),j,k,nbl)
				 	 !print*,PHID(NI(nblocks)-1,j,k,nbl)
				 	 !pause
				  
				  enddo
				  enddo
				  enddo
				  
			else

			  do nbl = 1,nblocks
			  do k = 1,NK(nbl)
			  do j = 1,NJ(nbl)
			  
				PHID(1,j,k,nbl) = 0.d0
				PHID(NI(nbl),j,k,nbl) = 0.d0
				PHID(2,j,k,nbl) = 0.d0
				PHID(NI(nbl)-1,j,k,nbl) = 0.d0
				do coeff = 1,7
				PHID(1,j,k,nbl) = PHID(1,j,k,nbl) + fdisc1(coeff)*PHI(coeff,j,k,nbl)
				PHID(2,j,k,nbl) = PHID(2,j,k,nbl) + fdisc2(coeff)*PHI(coeff,j,k,nbl)    !fdisc2 before???
				PHID(NI(nbl),j,k,nbl) = PHID(NI(nbl),j,k,nbl) - fdisc1(coeff)*PHI(NI(nbl)-coeff+1,j,k,nbl)
				PHID(NI(nbl)-1,j,k,nbl) = PHID(NI(nbl)-1,j,k,nbl) - fdisc2(coeff)*PHI(NI(nbl)-coeff+1,j,k,nbl)  !fdisc2 before?
				enddo
			  
			  enddo
			  enddo
			  enddo
		  
		    endif

      END
	  

!******************************************************************************
!******************** END DESCREIZATION I EXPLICIT GRID ***********************
!******************************************************************************
 

!******************************************************************************
!************************ DESCREIZATION I COMPACT GRID ************************
!******************************************************************************	
	  
      SUBROUTINE DISC_I_COMP_GRID(PHI,PHID)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars
		  real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
		  real,dimension(NImax) :: RHS
		  real LL
		  
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do j = 1,NJ(nbl)
		  
		  do i = 3,NI(nbl)-2
		  RHS(i) = (bdisc/4.d0)*(PHI(i+2,j,k,nbl)-PHI(i-2,j,k,nbl)) &
		 				 + (adisc/2.d0)*(PHI(i+1,j,k,nbl)-PHI(i-1,j,k,nbl))
		  enddo
		  
		  LL = PHI(NI(nbl),j,k,nbl) - PHI(1,j,k,nbl)
		  
		  RHS(1) = (bdisc/4.d0)*(PHI(3,j,k,nbl)-PHI(NI(nbl)-2,j,k,nbl)+LL) &
		 				 + (adisc/2.d0)*(PHI(2,j,k,nbl)-PHI(NI(nbl)-1,j,k,nbl)+LL)
		  RHS(2) = (bdisc/4.d0)*(PHI(4,j,k,nbl)-PHI(NI(nbl)-1,j,k,nbl)+LL) &
		 				 + (adisc/2.d0)*(PHI(3,j,k,nbl)-PHI(1,j,k,nbl))
		  RHS(NI(nbl)-2) = (bdisc/4.d0)*(PHI(1,j,k,nbl)-PHI(NI(nbl)-4,j,k,nbl)+LL) &
		 				 + (adisc/2.d0)*(PHI(NI(nbl)-1,j,k,nbl)-PHI(NI(nbl)-3,j,k,nbl))
		  RHS(NI(nbl)-1) = (bdisc/4.d0)*(PHI(2,j,k,nbl)-PHI(NI(nbl)-3,j,k,nbl)+LL) &
		 				 + (adisc/2.d0)*(PHI(1,j,k,nbl)-PHI(NI(nbl)-2,j,k,nbl)+LL)
		 				 
		  call TDMAP(1,NI(nbl)-1,AP_COMP,AC_COMP,AM_COMP,RHS,NI(nbl)-1)
		  
		  PHID(1:NI(nbl)-1,j,k,nbl) = RHS(1:NI(nbl)-1)
		  PHID(NI(nbl),j,k,nbl) = PHID(1,j,k,nbl)
		  
		  enddo
		  enddo
		  enddo
	  
	  End
	  
!******************************************************************************
!******************** END DESCREIZATION I COMPACT GRID ************************
!******************************************************************************



!******************************************************************************
!************************ DESCREIZATION J EXPLICIT GRID ***********************
!******************************************************************************		  

      SUBROUTINE DISC_J_EXP_GRID(PHI,PHID)
            use declare_variables
            implicit none

            REAL,DIMENSION(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
			real LL
			integer coeff

            do nbl = 1,nblocks
            do k = 1,NK(nblocks)
            do j = 3,NJ(nblocks)-2
            do i = 1,NI(nblocks)

                  PHID(i,j,k,nbl)  = (bdisc/4)*(PHI(i,j+2,k,nbl) -PHI(i,j-2,k,nbl)) &
                                          + (adisc/2)*(PHI(i,j+1,k,nbl) - PHI(i,j-1,k,nbl))
            enddo 
            enddo
            enddo
            enddo
			
			if(perJ.eq.1) then

				  do nbl = 1,nblocks
				  do k = 1,NK(nblocks)    
				  do i = 1,NI(nblocks)
				  LL = PHI(i,NJ(nbl),k,nbl) -PHI(i,1,k,nbl)
				 	 PHID(i,1,k,nbl)  = (bdisc/4.d0)*(PHI(i,3,k,nbl) -PHI(i,NJ(nbl)-2,k,nbl)+LL) &
				 						 + (adisc/2)*(PHI(i,2,k,nbl) - PHI(i,NJ(nbl)-1,k,nbl)+LL)
				  
				 	 PHID(i,2,k,nbl)  = (bdisc/4.d0)*(PHI(i,4,k,nbl) -PHI(i,NJ(nbl)-1,k,nbl)+LL) &
				 						 + (adisc/2)*(PHI(i,3,k,nbl) - PHI(i,1,k,nbl)) ! mihgt be an error here
				 																					 ! check where ll needs to be added	
				 	 PHID(i,NJ(nbl)-1,k,nbl)  = (bdisc/4.d0)*(PHI(i,2,k,nbl) -PHI(i,NJ(nbl)-3,k,nbl)+LL) &
				 						 + (adisc/2)*(PHI(i,NJ(nbl),k,nbl) - PHI(i,NJ(nbl)-2,k,nbl))
				  
				 	 PHID(i,NJ(nbl),k,nbl)    =  (bdisc/4.d0)*(PHI(i,3,k,nbl) -PHI(i,NJ(nbl)-2,k,nbl)+LL) &
				 						 + (adisc/2)*(PHI(i,2,k,nbl) - PHI(i,NJ(nbl)-1,k,nbl)+LL)
				  
				  enddo
				  enddo
				  enddo
				 
			else    ! boundary case
		  
			  do nbl = 1,nblocks
			  do k = 1,NK(nbl)
			  do i = 1,NI(nbl)
			  
				PHID(i,1,k,nbl) = 0.d0
				PHID(i,NJ(nbl),k,nbl) = 0.d0
				PHID(i,2,k,nbl) = 0.d0
				PHID(i,NJ(nbl)-1,k,nbl) = 0.d0
				do coeff = 1,7
				PHID(i,1,k,nbl) = PHID(i,1,k,nbl) + fdisc1(coeff)*PHI(i,coeff,k,nbl)
				PHID(i,2,k,nbl) = PHID(i,2,k,nbl) + fdisc2(coeff)*PHI(i,coeff,k,nbl)    !fdisc2 befoe?
				PHID(i,NJ(nbl),k,nbl) = PHID(i,NJ(nbl),k,nbl) - fdisc1(coeff)*PHI(i,NJ(nbl)-coeff+1,k,nbl)
				PHID(i,NJ(nbl)-1,k,nbl) = PHID(i,NJ(nbl)-1,k,nbl) - fdisc2(coeff)*PHI(i,NJ(nbl)-coeff+1,k,nbl)  !fdisc2 before??
				enddo
			  
			  enddo
			  enddo
			  enddo
		  
		    endif

      END
	  
!******************************************************************************
!******************** END DESCREIZATION J EXPLICIT GRID ***********************
!******************************************************************************


!******************************************************************************
!************************ DESCREIZATION J COMPACT GRID ************************
!******************************************************************************	

	  SUBROUTINE DISC_J_COMP_GRID(PHI,PHID)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars
		  real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
		  real,dimension(NJmax) :: RHS
		  real LL
		  
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do i = 1,NI(nbl)
		  
		  
			do j = 3,NJ(nbl)-2
			RHS(j) = (bdisc/4.d0)*(PHI(i,j+2,k,nbl)-PHI(i,j-2,k,nbl)) &
						+ (adisc/2.d0)*(PHI(i,j+1,k,nbl)-PHI(i,j-1,k,nbl))
			enddo
			
			LL = PHI(i,NJ(nbl),k,nbl) - PHI(i,1,k,nbl)
			
			RHS(1) = (bdisc/4.d0)*(PHI(i,3,k,nbl)-PHI(i,NJ(nbl)-2,k,nbl)+LL) &
						+ (adisc/2.d0)*(PHI(i,2,k,nbl)-PHI(i,NJ(nbl)-1,k,nbl)+LL)
			RHS(2) = (bdisc/4.d0)*(PHI(i,4,k,nbl)-PHI(i,NJ(nbl)-1,k,nbl)+LL) &
						+ (adisc/2.d0)*(PHI(i,3,k,nbl)-PHI(i,1,k,nbl))
			RHS(NJ(nbl)-2) = (bdisc/4.d0)*(PHI(i,1,k,nbl)-PHI(i,NJ(nbl)-4,k,nbl)+LL) &
						+ (adisc/2.d0)*(PHI(i,NJ(nbl)-1,k,nbl)-PHI(i,NJ(nbl)-3,k,nbl))
			RHS(NJ(nbl)-1) = (bdisc/4.d0)*(PHI(i,2,k,nbl)-PHI(i,NJ(nbl)-3,k,nbl)+LL) &
						+ (adisc/2.d0)*(PHI(i,1,k,nbl)-PHI(i,NJ(nbl)-2,k,nbl)+LL)
					
			call TDMAP(1,NJ(nbl)-1,AP_COMP,AC_COMP,AM_COMP,RHS,NJ(nbl)-1)
			
			PHID(i,1:NJ(nbl)-1,k,nbl) = RHS(1:NJ(nbl)-1)
			PHID(i,NJ(nbl),k,nbl) = PHID(i,1,k,nbl)
		  
		  enddo
		  enddo
		  enddo
		  
	  END
	  
!******************************************************************************
!******************** END DESCREIZATION J COMPACT GRID ************************
!******************************************************************************


!******************************************************************************
!************************ DESCREIZATION K EXPLICIT GRID ***********************
!******************************************************************************	

      SUBROUTINE DISC_K_EXP_GRID(PHI,PHID)
          use declare_variables
          implicit none

          REAL,DIMENSION(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
		  real LL
		  integer coeff

          do nbl = 1,nblocks
          do k = 3,NK(nblocks)-2
          do j = 1,NJ(nblocks)
          do i = 1,NI(nblocks)

                PHID(i,j,k,nbl)  = (bdisc/4)*(PHI(i,j,k+2,nbl) -PHI(i,j,k-2,nbl)) &
                                        + (adisc/2)*(PHI(i,j,k+1,nbl) - PHI(i,j,k-1,nbl))
          enddo 
          enddo
          enddo
          enddo
		  
		  if(perK.eq.1) then

			  do nbl = 1,nblocks
			  do j = 1,NJ(nblocks)    
			  do i = 1,NI(nblocks)
			  
			  LL = PHI(i,j,NK(nbl),nbl) -PHI(i,j,1,nbl)
			 	 
			 		 PHID(i,j,1,nbl)  = (bdisc/4.d0)*(PHI(i,j,3,nbl) -PHI(i,j,NK(nbl)-2,nbl)+LL) &
			 							 + (adisc/2.d0)*(PHI(i,j,2,nbl) - PHI(i,j,NK(nbl)-1,nbl)+LL)
			  
			 		 PHID(i,j,2,nbl)  = (bdisc/4.d0)*(PHI(i,j,4,nbl) -PHI(i,j,NK(nbl)-1,nbl)+LL) &
			 							 + (adisc/2.d0)*(PHI(i,j,3,nbl) - PHI(i,j,1,nbl)) ! might be an error here
			  
			 		 PHID(i,j,NK(nbl)-1,nbl)  = (bdisc/4.d0)*(PHI(i,j,2,nbl) -PHI(i,j,Nk(nbl)-3,nbl)+LL) &
			 								 + (adisc/2.d0)*(PHI(i,j,NK(nbl),nbl) - PHI(i,j,NK(nbl)-2,nbl))
			  
			 		 PHID(i,j,NK(nbl),nbl)    =  (bdisc/4.d0)*(PHI(i,j,3,nbl) -PHI(i,j,NK(nbl)-2,nbl)+LL) &
			 											 + (adisc/2.d0)*(PHI(i,j,2,nbl) - PHI(i,j,NK(nbl)-1,nbl)+LL)
			  
			  enddo
			  enddo
			  enddo
			 
		  else
		  
			  do nbl = 1,nblocks
			  do j = 1,NJ(nbl)
			  do i = 1,NI(nbl)
			  
				PHID(i,j,1,nbl) = 0.d0
				PHID(i,j,NK(nbl),nbl) = 0.d0
				PHID(i,j,2,nbl) = 0.d0
				PHID(i,j,NK(nbl)-1,nbl) = 0.d0
				do coeff = 1,7
				PHID(i,j,1,nbl) = PHID(i,j,1,nbl) + fdisc1(coeff)*PHI(i,j,coeff,nbl)
				PHID(i,j,2,nbl) = PHID(i,j,2,nbl) + fdisc2(coeff)*PHI(i,j,coeff,nbl)       !fdisc2 before??
				PHID(i,j,NK(nbl),nbl) = PHID(i,j,NK(nbl),nbl) - fdisc1(coeff)*PHI(i,j,NK(nbl)-coeff+1,nbl)
				PHID(i,j,NK(nbl)-1,nbl) = PHID(i,j,NK(nbl)-1,nbl) - fdisc2(coeff)*PHI(i,j,NK(nbl)-coeff+1,nbl)  !fdisc2 before??
				enddo
			  
			  enddo
			  enddo
			  enddo
		  
		  endif

      END
	  
!******************************************************************************
!******************** END DESCREIZATION K EXPLICIT GRID ***********************
!******************************************************************************


!******************************************************************************
!************************ DESCREIZATION K COMPACT GRID ************************
!******************************************************************************		  
	  
	  SUBROUTINE DISC_K_COMP_GRID(PHI,PHID)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars
		  real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
		  real,dimension(NKmax) :: RHS
		  real LL
		  
		  do nbl = 1,nblocks
		  do j = 1,NJ(nbl)
		  do i = 1,Ni(nbl)
		  
		  
			do k = 3,NK(nbl)-2
			RHS(k) = (bdisc/4.d0)*(PHI(i,j,k+2,nbl)-PHI(i,j,k-2,nbl)) &
							+ (adisc/2.d0)*(PHI(i,j,k+1,nbl)-PHI(i,j,k-1,nbl))
			enddo
			
			LL = PHI(i,j,NK(nbl),nbl) -PHI(i,j,1,nbl)
			
			RHS(1) = (bdisc/4.d0)*(PHI(i,j,3,nbl)-PHI(i,j,NK(nbl)-2,nbl)+LL) &
							+ (adisc/2.d0)*(PHI(i,j,2,nbl)-PHI(i,j,NK(nbl)-1,nbl)+LL)
			RHS(2) = (bdisc/4.d0)*(PHI(i,j,4,nbl)-PHI(i,j,NK(nbl)-1,nbl)+LL) &
							+ (adisc/2.d0)*(PHI(i,j,3,nbl)-PHI(i,j,1,nbl))
			RHS(NK(nbl)-2) = (bdisc/4.d0)*(PHI(i,j,1,nbl)-PHI(i,j,NK(nbl)-4,nbl)+LL) &
							+ (adisc/2.d0)*(PHI(i,j,NK(nbl)-1,nbl)-PHI(i,j,NK(nbl)-3,nbl))
			RHS(NK(nbl)-1) = (bdisc/4.d0)*(PHI(i,j,2,nbl)-PHI(i,j,NK(nbl)-3,nbl)+LL) &
							+ (adisc/2.d0)*(PHI(i,j,1,nbl)-PHI(i,j,NK(nbl)-2,nbl)+LL)
		 				 
			call TDMAP(1,NK(nbl)-1,AP_COMP,AC_COMP,AM_COMP,RHS,NK(nbl)-1)
			
			PHID(i,j,1:NK(nbl)-1,nbl) = RHS(1:NK(nbl)-1)
			PHID(i,j,NK(nbl),nbl) = PHID(i,j,1,nbl)
	  
		  enddo
		  enddo
		  enddo
	  
	  END
	  
!******************************************************************************
!******************** END DESCREIZATION K COMPACT GRID ************************
!******************************************************************************


!******************************************************************************
!************************ DESCREIZATION K2D EXPLICIT GRID *********************
!******************************************************************************	

      SUBROUTINE DISC_K2D_EXP_GRID(PHK2D,PHKD2D)
	    use declare_variables
	    implicit none

	    REAL,DIMENSION(NImax,NJmax,NKmax,nblocks) :: PHK2D,PHKD2D
		integer prim
		real LL
	
		do nbl = 1,nblocks
		do j = 1,NJ(nblocks)    
		do i = 1,NI(nblocks)
		
		LL = PHK2D(i,j,NK(nbl),nbl) -PHK2D(i,j,1,nbl)

			  PHKD2D(i,j,1,nbl)  = LL/2.d0
			  PHKD2D(i,j,2,nbl)  = LL/2.d0
			  PHKD2D(i,j,3,nbl)  = LL/2.d0
							
		enddo
		enddo
		enddo

      END
	  
!******************************************************************************
!******************** END DESCREIZATION K2D EXPLICIT GRID *********************
!******************************************************************************


!******************************************************************************
!****************************** FILTERING I ***********************************
!******************************************************************************

    SUBROUTINE FILTERING_I(PHI,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars,imc,ipc,coeff
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
		  real,dimension(NImax) :: RHS
		 	 
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do j = 1,NJ(nbl)
	  
			  do i = 6,NI(nbl)-5
			  RHS(i) = 0.d0
			  do coeff = 1,6
			  imc = i-coeff+1
			  ipc = i+coeff-1                     
			  RHS(i) = RHS(i) + &
			 				 0.5d0*fcoeff(coeff)*(PHI(imc,j,k,nbl,prim)+PHI(ipc,j,k,nbl,prim))
			  enddo
			  enddo
			  
			  do i = 1,5
			  RHS(i) = 0.d0
			  do coeff = 1,6
			  imc = i-coeff+1
			  if(imc.le.0) imc = imc + NI(nbl) - 1
			  ipc = i+coeff-1
			  RHS(i) = RHS(i) + &
			 				 0.5d0*fcoeff(coeff)*(PHI(imc,j,k,nbl,prim)+PHI(ipc,j,k,nbl,prim))	  
			  enddo
			  enddo
			  
			  do i = NI(nbl)-4,NI(nbl)
			  RHS(i) = 0.d0
			  do coeff = 1,6
			  imc = i-coeff+1
			  ipc = i+coeff-1
			  if(ipc.gt.NI(nbl)) ipc = ipc - NI(nbl) + 1
			  RHS(i) = RHS(i) + &
			 				 0.5d0*fcoeff(coeff)*(PHI(imc,j,k,nbl,prim)+PHI(ipc,j,k,nbl,prim))
			  enddo
			  enddo
			  
			  call TDMAP(1,NI(nbl)-1,AP,AC,AM,RHS,NI(nbl)-1)
			  !Do i = 1,NI(nbl)
			  !print*, PHI(i,j,k,nbl,prim),RHS(i)
			  !Enddo	  
			  PHI(1:NI(nbl)-1,j,k,nbl,prim) = RHS(1:NI(nbl)-1)
			  PHI(NI(nbl),j,k,nbl,prim) = PHI(1,j,k,nbl,prim)
	  
		  enddo
		  enddo
		  enddo
		  enddo
	  
	  END

!******************************************************************************
!************************** END FILTERING I ***********************************
!******************************************************************************	  


!******************************************************************************
!****************************** FILTERING J ***********************************
!******************************************************************************

    SUBROUTINE FILTERING_J(PHI,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars,imc,ipc,coeff
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
		  real,dimension(NJmax) :: RHS
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do k = 1,NK(nbl)
		  do i = 1,NI(nbl)
		  
			  do j = 6,NJ(nbl)-5
			  RHS(j) = 0.d0
			  do coeff = 1,6
			  imc = j-coeff+1
			  ipc = j+coeff-1                      ! some loop missing here 
			  RHS(j) = RHS(j) + &
			 		 0.5d0*fcoeff(coeff)*(PHI(i,imc,k,nbl,prim)+PHI(i,ipc,k,nbl,prim))  
			  enddo
			  enddo
		  
			  do j = 1,5
			  RHS(j) = 0.d0
			  do coeff = 1,6
			  imc = j-coeff+1
			  if(imc.le.0) imc = imc + NJ(nbl) - 1
			  ipc = j+coeff-1
			  RHS(j) = RHS(j) + &
			 				 0.5d0*fcoeff(coeff)*(PHI(i,imc,k,nbl,prim)+PHI(i,ipc,k,nbl,prim))
			  enddo
			  enddo
		  
		  
			  do j = NJ(nbl)-4,NJ(nbl)
			  RHS(j) = 0.d0
			  do coeff = 1,6
			  imc = j-coeff+1
			  ipc = j+coeff-1
			  if(ipc.gt.NJ(nbl)) ipc = ipc - NJ(nbl) + 1
			  RHS(j) = RHS(j) + &
			 				 0.5d0*fcoeff(coeff)*(PHI(i,imc,k,nbl,prim)+PHI(i,ipc,k,nbl,prim))
			  enddo
			  enddo
			  
			  call TDMAP(1,NJ(nbl)-1,AP,AC,AM,RHS,NJ(nbl)-1)
			  !Do j = 1,NJ(nbl)
			  !print*, PHI(i,j,k,nbl,prim),RHS(j)
			  !Enddo
			  PHI(i,1:NJ(nbl)-1,k,nbl,prim) = RHS(1:NJ(nbl)-1)
			  PHI(i,NJ(nbl),k,nbl,prim) = PHI(i,1,k,nbl,prim)
			  
		  enddo
		  enddo
		  enddo
		  enddo
	  
	  END

!******************************************************************************
!************************** END FILTERING J ***********************************
!******************************************************************************	


!******************************************************************************
!****************************** FILTERING K ***********************************
!******************************************************************************

    SUBROUTINE FILTERING_K(PHI,nvars)
		  use declare_variables
		  implicit none
		  
		  integer prim,nvars,imc,ipc,coeff
		  real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
		  real,dimension(NKmax) :: RHS
		  
		  do prim = 1,nvars
		  do nbl = 1,nblocks
		  do j = 1,NJ(nbl)
		  do i = 1,NI(nbl)
		  
			  do k = 6,NK(nbl)-5
			  RHS(k) = 0.d0	  
			  do coeff = 1,6
			  imc = k-coeff+1
			  ipc = k+coeff-1                      ! some loop missing here 		  
			  RHS(k) = RHS(k) + &
			 				 0.5d0*fcoeff(coeff)*(PHI(i,j,imc,nbl,prim)+PHI(i,j,ipc,nbl,prim))		  
			  enddo
			  enddo
			  
			  do k = 1,5
			  RHS(k) = 0.d0
			  do coeff = 1,6
			  imc = k-coeff+1
			  if(imc.le.0) imc = imc + NK(nbl) - 1
			  ipc = k+coeff-1		  
			  RHS(k) = RHS(k) + &
			 				 0.5d0*fcoeff(coeff)*(PHI(i,j,imc,nbl,prim)+PHI(i,j,ipc,nbl,prim))		  
			  enddo
			  enddo
			  
			  do k = NK(nbl)-4,NK(nbl)
			  RHS(k) = 0.d0
			  do coeff = 1,6
			  imc = k-coeff+1
			  ipc = k+coeff-1
			  if(ipc.gt.NI(nbl)) ipc = ipc - NK(nbl) + 1		  		  
			  RHS(k) = RHS(k) + &
			 				 0.5d0*fcoeff(coeff)*(PHI(i,j,imc,nbl,prim)+PHI(i,j,ipc,nbl,prim))		  
			  enddo
			  enddo
			  
			  call TDMAP(1,NK(nbl)-1,AP,AC,AM,RHS,NK(nbl)-1)
			  !Do k = 1,NK(nbl)
			  !print*, PHI(i,j,k,nbl,var),RHS(k)
			  !Enddo
			  PHI(i,j,1:NK(nbl)-1,nbl,prim) = RHS(1:NK(nbl)-1)
			  PHI(i,j,NK(nbl),nbl,prim) = PHI(i,j,1,nbl,prim)
			  
		  enddo
		  enddo
		  enddo
		  enddo
	  
	  END
	  
!******************************************************************************
!************************** END FILTERING K ***********************************
!******************************************************************************


!******************************************************************************
!******************************** CYCLIC TDMA *********************************
!******************************************************************************

	  SUBROUTINE TDMAP(ji,jf,ap,ac,am,fi,NMAXL)

		  !ap - super diagonal
		  !ac - diagonal
		  !am - sub diagonal
		  
		  implicit none
		  
		  ! -----------------------
		  !  Input/Output variables
		  ! -----------------------
		  integer:: ji, jf, NMAXL
		  real:: ap(NMAXL), ac(NMAXL), am(NMAXL), fi(NMAXL)
		  
		  ! -------------------
		  !  Internal variables
		  ! -------------------
		  integer:: i, j, ja, jj
		  real:: fnn, pp
		  real:: qq(NMAXL), ss(NMAXL), fei(NMAXL)
		  
		  ja=ji+1
		  jj=ji+jf
		  
		  qq(ji)=-ap(ji)/ac(ji)
		  ss(ji)=-am(ji)/ac(ji)
		  fnn=fi(jf)
		  fi(ji)=fi(ji)/ac(ji)
		  
		  !  forward elimination sweep
		  !----------------------------
		  do j=ja,jf
		 	 pp=1.0d0/(ac(j)+am(j)*qq(j-1))
		 		 qq(j)=-ap(j)*pp
		 		 ss(j)=-am(j)*ss(j-1)*pp
		 		 fi(j)=(fi(j)-am(j)*fi(j-1))*pp
		  enddo
		  
		  !  backward pass
		  !----------------
	  
		  ss(jf)=1.0d0
		  fei(jf)=0.0d0
		  
		  do i=ja,jf
		 	 j=jj-i
		 		 ss(j)=ss(j)+qq(j)*ss(j+1)
		 		 fei(j)=fi(j)+qq(j)*fei(j+1)
		  enddo
		  
		  fi(jf)=(fnn-ap(jf)*fei(ji)-am(jf)*fei(jf-1))/    &
		  &      (ap(jf)*ss(ji)+am(jf)*ss(jf-1)+ac(jf)) 
		  
		  !  backward substitution
		  !------------------------
		  do i=ja,jf
		 	 j=jj-i
		 		 fi(j)=fi(jf)*ss(j)+fei(j)
		  enddo
	  
	  END
	  
!******************************************************************************
!**************************** END CYCLIC TDMA *********************************
!******************************************************************************
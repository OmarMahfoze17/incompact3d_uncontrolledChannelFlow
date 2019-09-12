subroutine writePara(iStart,iFinish,iNAN,crrntTm,dt1)
USE variables
implicit none
real(8) ::dt1
integer :: nfil,i,iStart,iFinish,iNAN,crrntTm
character(len=90) :: lineInt,lineReal
!real(8) :: xlx
  print *, '*************** updat  incompact3d.prm *****************',xlx
nfil=41
10 format(I8,A)
20 format(f16.8,A)
  open(nfil,file='./incompact3d.prm')

  write(nfil,*) '#'
  write(nfil,*) '#             INCOMPACT 3D Flow parameters '
  write(nfil,*) '#' 
  write(nfil,20) xlx,'         #xlx  # Lx (Size of the box in x-direction)'
  write(nfil,20) yly,'          #yly  # Ly (Size of the box in y-direction) '
  write(nfil,20) zlz,'          #zlz  # Lz (Size of the box in z-direction) '
  write(nfil,20) 1/xnu,'        #re   # Reynolds number'
  write(nfil,20) 1.,'           #sc   # Schmidt number (if passive scalar)'
  write(nfil,20) 1.,'           #u1   # u1 (max velocity) (for inflow condition)'
  write(nfil,20) 1.,'           #u2   # u2 (min velocity) (for inflow condition)'
  write(nfil,20) 0.125,'          #noise# Turbulence intensity (1=100%) !! Initial condition'
  write(nfil,20) 0.0,'          #noise1# Turbulence intensity (1=100%) !! Inflow condition'
  write(nfil,20) dt1,'        #dt   # Time step'
  write(nfil,*) '#'
  write(nfil,*) '#             INCOMPACT3D Flow configuration'
  write(nfil,*) '#'
  write(nfil,10) nclx,'           #nclx   # nclx (BC)'
  write(nfil,10) ncly,'           #ncly   # ncly (BC)' 
  write(nfil,10) nclz,'           #nclz   # nclz (BC)'
  write(nfil,10) itype,'           #itype  # Type of flow'  
  write(nfil,10) 1,'           #iin    # Inflow condition (1: classic, 2: turbinit)'
  write(nfil,10) crrntTm+1,'           #ifirst # First iteration, attention add +1'   
  write(nfil,10) ilast,'           #ilast  # Last iteration '
  write(nfil,10) nscheme,'           #nscheme# Temporal scheme (1:AB2, 2: RK3, 3:RK4, 4:AB3)'
  write(nfil,10) istret,'          #istret # ymesh refinement(0:no, 1:center, 2:both sides, 3:bottom)'
  write(nfil,20) beta,'         #beta   # Refinement parameter (beta) '
  write(nfil,10) 1,'           #iskew  # (0:urotu, 1:skew, for the convective terms)'
  write(nfil,10) 0,'           #iscalar# (0: no scalar, 1:scalar)'
  write(nfil,*) '#'
  write(nfil,*) '#             INCOMPACT 3D File parameters '
  write(nfil,*) '#'
  write(nfil,10) 1,'           #ilit     # Read initial flow field ?'
  write(nfil,10) isave,'        #isave    # Frequency for writing backup file  '
  write(nfil,10) imodulo,'       #imodulo  # Frequency for visualization for VISU_INSTA'
  write(nfil,20) tStartSTAT,'        #tStartSTAT # collect Statistics after itime*dt passes this value'
  write(nfil,*) '#'
  write(nfil,*) '#             INCOMPACT 3D File parameters '
  write(nfil,*) '#'
  write(nfil,10) 0,'          #ivirt# IBM? (1: old school, 2: Lagrangian Poly)'
  write(nfil,20) 5.0,'          #cex  # X-centre position of the solid body'
  write(nfil,20) 5.0,'          #cey  # Y-centre position of the solid body'
  write(nfil,20) 0.0,'          #cez  # Z-centre position of the solid body'
  write(nfil,20) 0.5,'          #Re   # Radius of the solid body'
  write(nfil,*) '#'
 
  close(nfil)
  call system('nohup vim -c %s/\ // -c %s/\ // -c %s/\ // -c %s/\ // -c %s/\ // -c %s/\ // -c %s/\ // -c wq! ./incompact3d.prm &')
  call sleep(5)
end subroutine writePara 

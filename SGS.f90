!********************************************************************
!
subroutine SGS_stresses(ux1,uy1,uz1,smagx1,smagy1,smagz1)
! 
!********************************************************************
USE param
USE variables
USE decomp_2d


implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2 
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3

!WALE Arrays 
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1,syy1,szz1,&
     sxy1,sxz1,syz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxxd1,syyd1,szzd1,&
     sxyd1,sxzd1,syzd1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: nut1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: smagx1,smagy1,smagz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxx1,gyx1,gzx1,&
     gxy1,gyy1,gzy1,gxz1,gyz1,gzz1

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: sxx2,syy2,szz2,&
     sxy2,sxz2,syz2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: nut2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: smagx2,smagy2,smagz2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gxy2,gyy2,gzy2,&
     gxz2,gyz2,gzz2

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: sxx3,syy3,szz3,&
     sxy3,sxz3,syz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: nut3

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: smagx3,smagy3,smagz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: gxz3,gyz3,gzz3
! 

integer :: ijk,i,j,k
real(mytype) :: x,y,z,xcs

ta1=0.;tb1=0.;tc1=0.


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !FIRST COMPUTE gij,sij,sijd TENSOR TERMS FOR WALE LES   !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !WORK X-PENCILS
      !WALE Sxx, Sxz=Szx and Sxy=Syx terms
      call derx (sxx1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !TERME Sxx OK
      call derx (sxy1,uy1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
      call derx (sxz1,uz1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
      
      !Tenseur gradient vitesse
      gxx1(:,:,:)=sxx1(:,:,:) !TERME gxx OK
      gyx1(:,:,:)=sxy1(:,:,:) !TERME gyx OK
      gzx1(:,:,:)=sxz1(:,:,:) !TERME gzx OK

      !WORK Y-PENCILS
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_x_to_y(sxy1,ta2)
      !WALE Syx=Sxy, Syy and Syz=Szy terms
      call dery (sxy2,ux2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
      call dery (syy2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
      call dery (syz2,uz2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
      
      !Tenseur gradient vitesse
      gxy2(:,:,:)=sxy2(:,:,:)
      gyy2(:,:,:)=syy2(:,:,:)
      gzy2(:,:,:)=syz2(:,:,:)
      ! 
      sxy2(:,:,:)=0.5*(sxy2(:,:,:)+ta2(:,:,:)) 

      !WORK Z-PENCILS
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      call transpose_y_to_z(syz2,ta3)

      !WALE Sxz=Szx, Syz=Szy and Szz terms
      call derz(sxz3,ux3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
      call derz(syz3,uy3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
      call derz(szz3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

      !Tenseur gradient vitesse
      gxz3(:,:,:)=sxz3(:,:,:)
      gyz3(:,:,:)=syz3(:,:,:)
      gzz3(:,:,:)=szz3(:,:,:)
      ! 
      syz3(:,:,:)=0.5*(syz3(:,:,:)+ta3(:,:,:)) 

      !WORK Y-PENCILS
      call transpose_z_to_y(sxz3,sxz2)
      call transpose_z_to_y(syz3,syz2)
      call transpose_z_to_y(szz3,szz2)
      !Tenseur gradient vitesse
      call transpose_z_to_y(gxz3,gxz2)
      call transpose_z_to_y(gyz3,gyz2)
      call transpose_z_to_y(gzz3,gzz2)

      !WORK X-PENCILS
      call transpose_y_to_x(sxy2,sxy1) !TERME Sxy OK
      call transpose_y_to_x(syy2,syy1) !TERME Syy OK
      call transpose_y_to_x(syz2,syz1) !TERME Syz OK
      call transpose_y_to_x(szz2,szz1) !TERME Szz OK
      call transpose_y_to_x(sxz2,ta1)  !
      sxz1(:,:,:)=0.5*(sxz1(:,:,:)+ta1(:,:,:)) !TERME Sxz OK
   do j=1,ny-1
      dc(j)=(dx*dz*(yp(j+1)-yp(j)))**(1./3.)
   enddo
   dc(ny)=dc(ny-1)
   if (.true.) then 
		call SmogDynConst(ux1,uy1,uz1,sxx1,syy1,szz1,sxz1,sxy1,syz1,nut1)
   else 
		xcs=0.065
	   !!Then compute turbulent viscosity (SMAG OR WALE)
	   do k=1,xsize(3)
	      do j=1,xsize(2)
	         do i=1,xsize(1)
	            !SMAGORINSKY
	            nut1(i,j,k)=(xcs*dc(j+xstart(2)-1))**(2.)*&
	                 sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
	                 syy1(i,j,k)*syy1(i,j,k)+&
	                 szz1(i,j,k)*szz1(i,j,k)+&
	                 2.*sxy1(i,j,k)*sxy1(i,j,k)+&
	                 2.*sxz1(i,j,k)*sxz1(i,j,k)+&
	                 2.*syz1(i,j,k)*syz1(i,j,k)))
		         enddo
		      enddo
		   enddo
      endif

      ta1=0.
      ta2=0.
      ta3=0.
      !!And compute WALE extra terms
      !WORK X-PENCILS
      call derx (ta1,nut1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !ta1=d(nut1)/dx
      
      call derxx (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)!td1=d2(ux)/dx2
      call derxx (te1,uy1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)!te1=d2(uy)/dx2
      call derxx (tf1,uz1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)!tf1=d2(uz)/dx2

      smagx1(:,:,:)=td1(:,:,:)*nut1(:,:,:)+2.*sxx1(:,:,:)*ta1(:,:,:)
      smagy1(:,:,:)=te1(:,:,:)*nut1(:,:,:)+2.*sxy1(:,:,:)*ta1(:,:,:)
      smagz1(:,:,:)=tf1(:,:,:)*nut1(:,:,:)+2.*sxz1(:,:,:)*ta1(:,:,:)

      !WORK Y-PENCILS
      call transpose_x_to_y(smagx1,smagx2)
      call transpose_x_to_y(smagy1,smagy2)
      call transpose_x_to_y(smagz1,smagz2)
      call transpose_x_to_y(sxz1,sxz2)
      call transpose_x_to_y(nut1,nut2)

      call dery (ta2,nut2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)!ta2=d(nut2)/dy

      call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)!td2=d2(ux)/dy2
      call deryy (te2,uy2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)!te2=d2(uy)/dy2
      call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)!tf2=d2(uz)/dz2

      smagx2(:,:,:)=smagx2(:,:,:)+nut2(:,:,:)*td2(:,:,:)+2.*sxy2(:,:,:)*ta2(:,:,:)
      smagy2(:,:,:)=smagy2(:,:,:)+nut2(:,:,:)*te2(:,:,:)+2.*syy2(:,:,:)*ta2(:,:,:)
      smagz2(:,:,:)=smagz2(:,:,:)+nut2(:,:,:)*tf2(:,:,:)+2.*syz2(:,:,:)*ta2(:,:,:)

      !WORK Z-PENCILS
      call transpose_y_to_z(smagx2,smagx3)      
      call transpose_y_to_z(smagy2,smagy3)
      call transpose_y_to_z(smagz2,smagz3)
      call transpose_y_to_z(sxz2,sxz3)
      call transpose_y_to_z(nut2,nut3)
      
      call derz(ta3,nut3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)!ta3=d(nut3)/dz

      call derzz (td3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)!td3=d2(ux)/dz2
      call derzz (te3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)!te3=d2(uy)/dz2
      call derzz (tf3,uz3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)!tf3=d2(uz)/dz2

      smagx3(:,:,:)=smagx3(:,:,:)+nut3(:,:,:)*td3(:,:,:)+2.*sxz3(:,:,:)*ta3(:,:,:)
      smagy3(:,:,:)=smagy3(:,:,:)+nut3(:,:,:)*te3(:,:,:)+2.*syz3(:,:,:)*ta3(:,:,:)
      smagz3(:,:,:)=smagz3(:,:,:)+nut3(:,:,:)*tf3(:,:,:)+2.*szz3(:,:,:)*ta3(:,:,:)

      call transpose_z_to_y(smagx3,smagx2)      
      call transpose_z_to_y(smagy3,smagy2)      
      call transpose_z_to_y(smagz3,smagz2)      

      call transpose_y_to_x(smagx2,smagx1)      
      call transpose_y_to_x(smagy2,smagy1)      
      call transpose_y_to_x(smagz2,smagz1)      

end subroutine SGS_stresses


!********************************************************************
!
subroutine SmogDynConst(ux1,uy1,uz1,sxx1,syy1,szz1,sxz1,sxy1,syz1,nut1)
! 
!********************************************************************
USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE MPI

implicit none
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ux_test1,ux_test_f1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxx1,uyy1,uzz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxy1,uxz1,uyz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxy1f,uxz1f,uyz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1f,uy1f,uz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxx1f,uyy1f,uzz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: lxx1,lyy1,lzz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: lxy1,lxz1,lyz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1,syy1,szz1,sxy1,sxz1,syz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1f,syy1f,szz1f,sxy1f,sxz1f,syz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: bbxx1,bbyy1,bbzz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: bbxy1,bbxz1,bbyz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: axx1,ayy1,azz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: axy1,axz1,ayz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: axx1f,ayy1f,azz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: axy1f,axz1f,ayz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: mxx1,myy1,mzz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: mxy1,mxz1,myz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: nut1,smagC,smagC1f,Cs1
!-------------------------------------------------------------------------
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2,ux_test2,ux_test_f2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2f,uy2f,uz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uxx2f,uyy2f,uzz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uxy2f,uxz2f,uyz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: sxx2,syy2,szz2,sxy2,sxz2,syz2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: sxx2f,syy2f,szz2f,sxy2f,sxz2f,syz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: axx2,ayy2,azz2,axy2,axz2,ayz2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: axx2f,ayy2f,azz2f,axy2f,axz2f,ayz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: smagx2,smagy2,smagz2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: smagC2,smagC2f,Cs2
!-------------------------------------------------------------------------
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3,ux_test3,ux_test_f3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3f,uy3f,uz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uxx3f,uyy3f,uzz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uxy3f,uxz3f,uyz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: axx3,ayy3,azz3,axy3,axz3,ayz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: axx3f,ayy3f,azz3f,axy3f,axz3f,ayz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: sxx3,syy3,szz3,sxy3,sxz3,syz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: sxx3f,syy3f,szz3f,sxy3f,sxz3f,syz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: smagx3,smagy3,smagz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: smagC3,Cmz3,smagC3f,Cs3

real(mytype) :: Cs,Csm,Css,Cssm
integer :: ijk,nvect1,nvect2,nvect3,i,j,k,code
character(len=20) :: filename
!TEST FILTER for ui and uiuj 
   do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              uxx1(i,j,k)=ux1(i,j,k)*ux1(i,j,k)
              uyy1(i,j,k)=uy1(i,j,k)*uy1(i,j,k)
              uzz1(i,j,k)=uz1(i,j,k)*uz1(i,j,k)
              uxy1(i,j,k)=ux1(i,j,k)*uy1(i,j,k)
              uxz1(i,j,k)=ux1(i,j,k)*uz1(i,j,k)
              uyz1(i,j,k)=uy1(i,j,k)*uz1(i,j,k)
           enddo
        enddo
     enddo

     call filx(ux1f,ux1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!     

     call filx(uy1f,uy1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!  

     call filx(uz1f,uz1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!      

     call filx(uxx1f,uxx1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!        

     call filx(uyy1f,uyy1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!     

     call filx(uzz1f,uzz1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!     

     call filx(uxy1f,uxy1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!     

     call filx(uxz1f,uxz1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!     

     call filx(uyz1f,uyz1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!     
               

     call transpose_x_to_y(ux1f,ta2)
     call transpose_x_to_y(uy1f,tb2)
     call transpose_x_to_y(uz1f,tc2)
     call transpose_x_to_y(uxx1f,td2)
     call transpose_x_to_y(uyy1f,te2)
     call transpose_x_to_y(uzz1f,tf2)
     call transpose_x_to_y(uxy1f,tg2)
     call transpose_x_to_y(uxz1f,th2)
     call transpose_x_to_y(uyz1f,ti2)

     call fily(ux2f,ta2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(uy2f,tb2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(uz2f,tc2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(uxx2f,td2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(uyy2f,te2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(uzz2f,tf2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(uxy2f,tg2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(uxz2f,th2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(uyz2f,ti2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     


     ta2=0.;tb2=0.;tc2=0.
     td2=0.;te2=0.;tf2=0.
     tg2=0.;th2=0.;ti2=0.

     !
     call transpose_y_to_z(ux2f,ta3)
     call transpose_y_to_z(uy2f,tb3)
     call transpose_y_to_z(uz2f,tc3)
     call transpose_y_to_z(uxx2f,td3)
     call transpose_y_to_z(uyy2f,te3)
     call transpose_y_to_z(uzz2f,tf3)
     call transpose_y_to_z(uxy2f,tg3)
     call transpose_y_to_z(uxz2f,th3)
     call transpose_y_to_z(uyz2f,ti3)


     call filz(ux3f,ta3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(uy3f,tb3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(uz3f,tc3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(uxx3f,td3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(uyy3f,te3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(uzz3f,tf3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(uxy3f,tg3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(uxz3f,th3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(uyz3f,ti3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.
               

     ta3=0.;tb3=0.;tc3=0.
     td3=0.;te3=0.;tf3=0.
     tg3=0.;th3=0.;ti3=0.

     ux2f=0.;uy2f=0.;uz2f=0.
     call transpose_z_to_y(ux3f,ux2f)      
     call transpose_z_to_y(uy3f,uy2f)      
     call transpose_z_to_y(uz3f,uz2f)   
     call transpose_z_to_y(uxx3f,uxx2f)      
     call transpose_z_to_y(uyy3f,uyy2f)      
     call transpose_z_to_y(uzz3f,uzz2f)    
     call transpose_z_to_y(uxy3f,uxy2f)      
     call transpose_z_to_y(uxz3f,uxz2f)      
     call transpose_z_to_y(uyz3f,uyz2f)    

     ux1f=0.;uy1f=0.;uz1f=0.
     call transpose_y_to_x(ux2f,ux1f)    
     call transpose_y_to_x(uy2f,uy1f)     
     call transpose_y_to_x(uz2f,uz1f)   
     call transpose_y_to_x(uxx2f,uxx1f)    
     call transpose_y_to_x(uyy2f,uyy1f)     
     call transpose_y_to_x(uzz2f,uzz1f)   
     call transpose_y_to_x(uxy2f,uxy1f)    
     call transpose_y_to_x(uxz2f,uxz1f)     
     call transpose_y_to_x(uyz2f,uyz1f)   

     !Lij tensor OK
     lxx1(:,:,:)=uxx1f(:,:,:)-ux1f(:,:,:)*ux1f(:,:,:)
     lyy1(:,:,:)=uyy1f(:,:,:)-uy1f(:,:,:)*uy1f(:,:,:)
     lzz1(:,:,:)=uzz1f(:,:,:)-uz1f(:,:,:)*uz1f(:,:,:)
     lxy1(:,:,:)=uxy1f(:,:,:)-ux1f(:,:,:)*uy1f(:,:,:)
     lxz1(:,:,:)=uxz1f(:,:,:)-ux1f(:,:,:)*uz1f(:,:,:)
     lyz1(:,:,:)=uyz1f(:,:,:)-uy1f(:,:,:)*uz1f(:,:,:)
     
!======================================================================     
!================================= M_ij ===============================
!==========================================xcs============================     
!Sij and Sij^test

     !WORK X-PENCILS
     !LES Sxx, Sxz=Szx and Sxy=Syx terms
!     call derx (sxx1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !TERME Sxx OK !passed from SGS_stresses
!     call derx (sxy1,uy1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !passed from SGS_stresses
!     call derx (sxz1,uz1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !passed from SGS_stresses
     !Test filter
     call derx (sxx1f,ux1f,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !TERME Sxxf OK
     call derx (sxy1f,uy1f,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
     call derx (sxz1f,uz1f,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)

     !WORK Y-PENCILS
     call transpose_x_to_y(ux1,ux2)
     call transpose_x_to_y(uy1,uy2)
     call transpose_x_to_y(uz1,uz2)
!     call transpose_x_to_y(sxy1,ta2) !passed from SGS_stresses

     call transpose_x_to_y(ux1f,ux2f)
     call transpose_x_to_y(uy1f,uy2f)
     call transpose_x_to_y(uz1f,uz2f)
     call transpose_x_to_y(sxy1f,tb2)

     !Syx=Sxy, Syy and Syz=Szy terms
!     call dery (sxy2,ux2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) !passed from SGS_stresses
!     call dery (syy2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) !passed from SGS_stresses
!     call dery (syz2,uz2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) !passed from SGS_stresses

!     sxy2(:,:,:)=0.5*(sxy2(:,:,:)+ta2(:,:,:))  !passed from SGS_stresses
     !Test filter
     call dery (sxy2f,ux2f,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
     call dery (syy2f,uy2f,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
     call dery (syz2f,uz2f,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)

     sxy2f(:,:,:)=0.5*(sxy2f(:,:,:)+tb2(:,:,:)) 

     !WORK Z-PENCILS
     call transpose_y_to_z(ux2,ux3)
     call transpose_y_to_z(uy2,uy3)
     call transpose_y_to_z(uz2,uz3)
!     call transpose_y_to_z(syz2,ta3) !passed from SGS_stresses

     call transpose_y_to_z(ux2f,ux3f)
     call transpose_y_to_z(uy2f,uy3f)
     call transpose_y_to_z(uz2f,uz3f)
     call transpose_y_to_z(syz2f,tb3)

     !Sxz=Szx, Syz=Szy and Szz terms
!     call derz(sxz3,ux3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) !passed from SGS_stresses
!     call derz(syz3,uy3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) !passed from SGS_stresses
!     call derz(szz3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) !passed from SGS_stresses

!     syz3(:,:,:)=0.5*(syz3(:,:,:)+ta3(:,:,:))  !passed from SGS_stresses
     !Test filter
     call derz(sxz3f,ux3f,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
     call derz(syz3f,uy3f,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
     call derz(szz3f,uz3f,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

     syz3f(:,:,:)=0.5*(syz3f(:,:,:)+tb3(:,:,:)) 

     !WORK Y-PENCILS
!     call transpose_z_to_y(sxz3,sxz2) !passed from SGS_stresses
!     call transpose_z_to_y(syz3,syz2)  !passed from SGS_stresses
!     call transpose_z_to_y(szz3,szz2) !passed from SGS_stresses

     call transpose_z_to_y(sxz3f,sxz2f)
     call transpose_z_to_y(syz3f,syz2f)
     call transpose_z_to_y(szz3f,szz2f)

     !WORK X-PENCILS
!     call transpose_y_to_x(sxy2,sxy1) !TERME Sxy OK  !passed from SGS_stresses
!     call transpose_y_to_x(syy2,syy1) !TERME Syy OK  !passed from SGS_stresses
!     call transpose_y_to_x(syz2,syz1) !TERME Syz OK   !passed from SGS_stresses
!     call transpose_y_to_x(szz2,szz1) !TERME Szz OK   !passed from SGS_stresses
!     call transpose_y_to_x(sxz2,ta1)  ! !passed from SGS_stresses
!     sxz1(:,:,:)=0.5*(sxz1(:,:,:)+ta1(:,:,:)) !TERME Sxz OK  !passed from SGS_stresses

     call transpose_y_to_x(sxy2f,sxy1f) !TERME Sxyf OK
     call transpose_y_to_x(syy2f,syy1f) !TERME Syyf OK
     call transpose_y_to_x(syz2f,syz1f) !TERME Syzf OK
     call transpose_y_to_x(szz2f,szz1f) !TERME Szzf OK
     call transpose_y_to_x(sxz2f,tb1)  !
     sxz1f(:,:,:)=0.5*(sxz1f(:,:,:)+tb1(:,:,:)) !TERME Sxzf OK
     !

     !!Then compute B and A tensors
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)

              !Bij tensor with u test filtered OK
              bbxx1(i,j,k)=-2*(2.*dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
                   syy1f(i,j,k)*syy1f(i,j,k)+&
                   szz1f(i,j,k)*szz1f(i,j,k)+&
                   2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
                   2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
                   2.*syz1f(i,j,k)*syz1f(i,j,k)))&
                   *sxx1f(i,j,k)

              bbyy1(i,j,k)=-2*(2.*dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
                   syy1f(i,j,k)*syy1f(i,j,k)+&
                   szz1f(i,j,k)*szz1f(i,j,k)+&
                   2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
                   2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
                   2.*syz1f(i,j,k)*syz1f(i,j,k)))&
                   *syy1f(i,j,k)

              bbzz1(i,j,k)=-2*(2.*dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
                   syy1f(i,j,k)*syy1f(i,j,k)+&
                   szz1f(i,j,k)*szz1f(i,j,k)+&
                   2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
                   2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
                   2.*syz1f(i,j,k)*syz1f(i,j,k)))&
                   *szz1f(i,j,k)

              bbxy1(i,j,k)=-2*(2.*dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
                   syy1f(i,j,k)*syy1f(i,j,k)+&
                   szz1f(i,j,k)*szz1f(i,j,k)+&
                   2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
                   2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
                   2.*syz1f(i,j,k)*syz1f(i,j,k)))&
                   *sxy1f(i,j,k)

              bbxz1(i,j,k)=-2*(2.*dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
                   syy1f(i,j,k)*syy1f(i,j,k)+&
                   szz1f(i,j,k)*szz1f(i,j,k)+&
                   2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
                   2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
                   2.*syz1f(i,j,k)*syz1f(i,j,k)))&
                   *sxz1f(i,j,k)

              bbyz1(i,j,k)=-2*(2.*dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
                   syy1f(i,j,k)*syy1f(i,j,k)+&
                   szz1f(i,j,k)*szz1f(i,j,k)+&
                   2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
                   2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
                   2.*syz1f(i,j,k)*syz1f(i,j,k)))&
                   *syz1f(i,j,k)
              !

              !Aij tensor with u 
              axx1(i,j,k)=-2*(dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
                   syy1(i,j,k)*syy1(i,j,k)+&
                   szz1(i,j,k)*szz1(i,j,k)+&
                   2.*sxy1(i,j,k)*sxy1(i,j,k)+&
                   2.*sxz1(i,j,k)*sxz1(i,j,k)+&
                   2.*syz1(i,j,k)*syz1(i,j,k)))&
                   *sxx1(i,j,k)

              ayy1(i,j,k)=-2*(dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
                   syy1(i,j,k)*syy1(i,j,k)+&
                   szz1(i,j,k)*szz1(i,j,k)+&
                   2.*sxy1(i,j,k)*sxy1(i,j,k)+&
                   2.*sxz1(i,j,k)*sxz1(i,j,k)+&
                   2.*syz1(i,j,k)*syz1(i,j,k)))&
                   *syy1(i,j,k)

              azz1(i,j,k)=-2*(dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
                   syy1(i,j,k)*syy1(i,j,k)+&
                   szz1(i,j,k)*szz1(i,j,k)+&
                   2.*sxy1(i,j,k)*sxy1(i,j,k)+&
                   2.*sxz1(i,j,k)*sxz1(i,j,k)+&
                   2.*syz1(i,j,k)*syz1(i,j,k)))&
                   *szz1(i,j,k)

              axy1(i,j,k)=-2*(dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
                   syy1(i,j,k)*syy1(i,j,k)+&
                   szz1(i,j,k)*szz1(i,j,k)+&
                   2.*sxy1(i,j,k)*sxy1(i,j,k)+&
                   2.*sxz1(i,j,k)*sxz1(i,j,k)+&
                   2.*syz1(i,j,k)*syz1(i,j,k)))&
                   *sxy1(i,j,k)

              axz1(i,j,k)=-2*(dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
                   syy1(i,j,k)*syy1(i,j,k)+&
                   szz1(i,j,k)*szz1(i,j,k)+&
                   2.*sxy1(i,j,k)*sxy1(i,j,k)+&
                   2.*sxz1(i,j,k)*sxz1(i,j,k)+&
                   2.*syz1(i,j,k)*syz1(i,j,k)))&
                   *sxz1(i,j,k)

              ayz1(i,j,k)=-2*(dc(j+xstart(2)-1))**2. & 
                   *sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
                   syy1(i,j,k)*syy1(i,j,k)+&
                   szz1(i,j,k)*szz1(i,j,k)+&
                   2.*sxy1(i,j,k)*sxy1(i,j,k)+&
                   2.*sxz1(i,j,k)*sxz1(i,j,k)+&
                   2.*syz1(i,j,k)*syz1(i,j,k)))&
                   *syz1(i,j,k)
              !            
           enddo
        enddo
     enddo

     !Need to filter Aij components
     call filx(axx1f,axx1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!
              
     call filx(ayy1f,ayy1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!
               
     call filx(azz1f,azz1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!
               
     call filx(axy1f,axy1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!
               
     call filx(axz1f,axz1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!
               
     call filx(ayz1f,ayz1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!
     
     
     call transpose_x_to_y(axx1f,ta2)
     call transpose_x_to_y(ayy1f,tb2)
     call transpose_x_to_y(azz1f,tc2)
     call transpose_x_to_y(axy1f,td2)
     call transpose_x_to_y(axz1f,te2)
     call transpose_x_to_y(ayz1f,tf2)

     call fily(axx2f,ta2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(ayy2f,tb2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(azz2f,tc2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(axy2f,td2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(axz2f,te2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     

     call fily(ayz2f,tf2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!     


     ta2=0.;tb2=0.;tc2=0.
     td2=0.;te2=0.;tf2=0.

     call transpose_y_to_z(axx2f,ta3)
     call transpose_y_to_z(ayy2f,tb3)
     call transpose_y_to_z(azz2f,tc3)
     call transpose_y_to_z(axy2f,td3)
     call transpose_y_to_z(axz2f,te3)
     call transpose_y_to_z(ayz2f,tf3)

     call filz(axx3f,ta3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(ayy3f,tb3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(azz3f,tc3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(axy3f,td3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(axz3f,te3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.

     call filz(ayz3f,tf3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.
               

     ta3=0.;tb3=0.;tc3=0.
     td3=0.;te3=0.;tf3=0.

     call transpose_z_to_y(axx3f,axx2f)      
     call transpose_z_to_y(ayy3f,ayy2f)      
     call transpose_z_to_y(azz3f,azz2f)    
     call transpose_z_to_y(axy3f,axy2f)      
     call transpose_z_to_y(axz3f,axz2f)      
     call transpose_z_to_y(ayz3f,ayz2f)    

     call transpose_y_to_x(axx2f,axx1f)    
     call transpose_y_to_x(ayy2f,ayy1f)     
     call transpose_y_to_x(azz2f,azz1f)   
     call transpose_y_to_x(axy2f,axy1f)    
     call transpose_y_to_x(axz2f,axz1f)     
     call transpose_y_to_x(ayz2f,ayz1f)   

     !Mij tensor OK
     mxx1(:,:,:)=bbxx1(:,:,:)-axx1f(:,:,:)
     myy1(:,:,:)=bbyy1(:,:,:)-ayy1f(:,:,:)
     mzz1(:,:,:)=bbzz1(:,:,:)-azz1f(:,:,:)
     mxy1(:,:,:)=bbxy1(:,:,:)-axy1f(:,:,:)
     mxz1(:,:,:)=bbxz1(:,:,:)-axz1f(:,:,:)
     myz1(:,:,:)=bbyz1(:,:,:)-ayz1f(:,:,:)
     !

     !Lij deviator
     lxx1(:,:,:)=lxx1(:,:,:)-(lxx1(:,:,:)+lyy1(:,:,:)+lzz1(:,:,:))/3.
     lyy1(:,:,:)=lyy1(:,:,:)-(lxx1(:,:,:)+lyy1(:,:,:)+lzz1(:,:,:))/3.
     lzz1(:,:,:)=lzz1(:,:,:)-(lxx1(:,:,:)+lyy1(:,:,:)+lzz1(:,:,:))/3.

     !Compute the Smagorinsky constant dynamically
     !Original method without averaging
     if(itime==1) then
        smagC(:,:,:)=sqrt(0.18)
     else
        do k=1,xsize(3)
           do j=1,xsize(2)
              do i=1,xsize(1)
                 smagC(i,j,k)= (lxx1(i,j,k)*mxx1(i,j,k)+&
                      lyy1(i,j,k)*myy1(i,j,k)+&
                      lzz1(i,j,k)*mzz1(i,j,k)+&
                      2.*lxy1(i,j,k)*mxy1(i,j,k)+&
                      2.*lxz1(i,j,k)*mxz1(i,j,k)+&
                      2.*lyz1(i,j,k)*myz1(i,j,k))/&
                      (mxx1(i,j,k)*mxx1(i,j,k)+&
                      myy1(i,j,k)*myy1(i,j,k)+&
                      mzz1(i,j,k)*mzz1(i,j,k)+&
                      2.*mxy1(i,j,k)*mxy1(i,j,k)+&
                      2.*mxz1(i,j,k)*mxz1(i,j,k)+&
                      2.*myz1(i,j,k)*myz1(i,j,k))
                 
                 !ERIC LIMITEUR SI BESOIN
                 !if(smagC(i,j,k).lt.0) then 
                 !   smagC(i,j,k)=0.
                 !endif 
                 !
                 
              enddo
           enddo
        enddo
     endif

     !FILTERING THE NON-CONSTANT CONSTANT
	  call filx(smagC1f,smagC,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0) 
   !!!

     call transpose_x_to_y(smagC1f,ta2)

     call fily(smagC2f,ta2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)
  

     ta2=0.
     call transpose_y_to_z(smagC2f,ta3)

     call filz(smagC3f,ta3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)   !!!!!!!!!!!!.
               
     ta3=0.
     call transpose_z_to_y(smagC3f,smagC2f)
     call transpose_y_to_x(smagC2f,smagC)

     !ERIC : SPATIAL AVERAGING PROCEDURE (PLANS (x,z))
     !call transpose_x_to_y(smagC,smagC2)
     !do j=1,ysize(2)
     !   SmagCm(j)=0.
     !   do k=1,ysize(3)
     !      do i=1,ysize(1)
     !         SmagCm(j)=SmagCm(j)+smagC2(i,j,k)
     !      enddo
     !   enddo
     !   SmagCm(j)=SmagCm(j)/ysize(1)/ysize(3)
     !enddo
     !call MPI_ALLREDUCE(SmagCm,SmagCmm,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     !SmagCmm(:)=SmagCmm(:)/nproc
     
     !ERIC : UNCHECKED SPATIAL AVERAGING PROCEDURE (DOMAINE COMPLET)
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=2,xsize(1)
              Cs1(i,j,k)=DMAX1(0.0,smagC(i,j,k))
           enddo
        enddo
     enddo  
     
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=2,xsize(1)
              Cs1(1,j,k)=Cs1(1,j,k)+Cs1(i,j,k)
           enddo
        enddo
     enddo     
     do k=1,xsize(3)
        do j=1,xsize(2)
              Cs1(1,j,k)=Cs1(1,j,k)/float(xsize(1))
              Cs1(:,j,k)=Cs1(1,j,k)
        enddo
     enddo          
!     do k=1,xsize(3)
!        do j=1,xsize(2)
!           do i=2,xsize(1)
!              Cs1(i,j,k)=Cs1(1,j,k)
!           enddo
!        enddo
!     enddo
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
     Cs3(:,:,:)=0.
     call transpose_x_to_y(Cs1,Cs2)
     call transpose_y_to_z(Cs2,Cs3)
     do k=2,zsize(3)
        do j=1,zsize(2)
           do i=1,zsize(1)
              Cs3(i,j,1)=Cs3(i,j,1)+Cs3(i,j,k)
           enddo
        enddo
     enddo     
     do j=1,zsize(2)
        do i=1,zsize(1)
              Cs3(i,j,1)=Cs3(i,j,1) / float(zsize(3))
              Cs3(i, j, :) = Cs3(i, j, 1)
        enddo
     enddo          
!     do k=2,zsize(3)
!        do j=1,zsize(2)
!           do i=1,zsize(1)
!              Cs3(i,j,k)=Cs3(i,j,1)
!           enddo
!        enddo
!     enddo
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call transpose_z_to_y(Cs3,Cs2)
     call transpose_y_to_x(Cs2,Cs1)
     
     
     Cs=0.
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              Cs=Cs+smagC(i,j,k)
           enddo
        enddo
     enddo
     Cs=Cs/xsize(1)/xsize(2)/xsize(3)
     call MPI_ALLREDUCE(Cs,Csm,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     Csm=Csm/nproc
     if (nrank .eq. 0) then
		  print *, '------------------------'
		  print *, 'Smog Constant is ', Csm,Cssm
		  print *, '------------------------'
     endif
     
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)

              !Smagorinsky model
              nut1(i,j,k)=Cs1(i,j,k)*(dc(j+xstart(2)-1))**(2.)*&
                   sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
                   syy1(i,j,k)*syy1(i,j,k)+&
                   szz1(i,j,k)*szz1(i,j,k)+&
                   2.*sxy1(i,j,k)*sxy1(i,j,k)+&
                   2.*sxz1(i,j,k)*sxz1(i,j,k)+&
                   2.*syz1(i,j,k)*syz1(i,j,k)))  

              !Version moyenne dans les plans (x,z)
              !SmagCmm(j+xstart(2)-1)*(dc(j+xstart(2)-1))**(2.)*&
              !        sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
              !        syy1(i,j,k)*syy1(i,j,k)+&
              !        szz1(i,j,k)*szz1(i,j,k)+&
              !        2.*sxy1(i,j,k)*sxy1(i,j,k)+&
              !        2.*sxz1(i,j,k)*sxz1(i,j,k)+&
              !        2.*syz1(i,j,k)*syz1(i,j,k)))

              enddo
           enddo
        enddo
        
!!!!        do k=1,xsize(3)
!!!!        do j=1,xsize(2)
!!!!        do i=1,xsize(1)
!!!!       		ux_test1(i,j,k)=sin(dx*(i-1)/xlx*pi*2)+0.2*sin(dx*(i-1)/xlx*pi*200)
!!!!        		ux_test1(i,j,k)=ux_test1(i,j,k)+ sin(yp(j+xstart(2)-1)/yly*pi*2)+0.2*sin(yp(j+xstart(2)-1)/yly*pi*200)
!!!!        		ux_test1(i,j,k)=ux_test1(i,j,k)+sin(dz*(k-1+xstart(3)-1)/zlz*pi*2)+0.2*sin(dz*(k-1+xstart(3)-1)/zlz*pi*200)

!!!!        enddo
!!!!        enddo
!!!!        enddo

!!!!        call filx(ux_test_f1,ux_test1,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0)    !!!
!!!!        call transpose_x_to_y(ux_test_f1,ux_test2)
!!!!        call fily(ux_test_f2,ux_test2,fiffy,fifsy,fifwy,ppy,ysize(1),ysize(2),ysize(3),1)  !!!!
!!!!        call transpose_y_to_z(ux_test_f2,ux_test3)
!!!!        call filz(ux_test_f3,ux_test3,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),1)
!!!!        
!!!!        call transpose_z_to_y(ux_test_f3,ux_test_f2)
!!!!        call transpose_y_to_x(ux_test_f2,ux_test_f1)        
!!!!        
!!!!        
!!!!        if (mod(itime,2)==0) then
!!!!1001    format('ux_test',I3.3)
!!!!        write(filename, 1001) itime/2
!!!!        call decomp_2d_write_one(1,ux_test1,filename,2)
!!!!1002    format('ux_test_f',I3.3)
!!!!        write(filename, 1002) itime/2
!!!!        call decomp_2d_write_one(1,ux_test_f1,filename,2)
!!!!     end if
     
     
    if (mod(itime,imodulo)==0) then
1003    format('Q_',I3.3)
        write(filename, 1003) itime/imodulo
        call decomp_2d_write_one(1,smagC,filename,2)


1004    format('uz',I3.3)
        write(filename, 1004) itime/imodulo
        call decomp_2d_write_one(1,Cs1,filename,2)
     end if
     
end subroutine SmogDynConst

!***************************************************************************
!
subroutine matlab_data(ux,uy,uz,fz,fy) 
!
!***************************************************************************

USE param
USE variables


implicit none

integer  :: i,j,k
real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,epsi,p_x,p_y,p_z
real(8),dimension(nx,ny,nz) :: tx,ty,tzdi1,di2,px,py,pz
real(8),dimension(nx,ny,nz) :: uxm1,uym1,uzm1,uxm2,uym2,uzm2
real(8),dimension(nxm,nym,nzm) :: ppm
integer :: nxyz
real(8),dimension(nz,ny) :: fz,fy
character(len=20) :: filename


801 format('matlab',I4.4)
write(filename, 801) itime/imodulo

nxyz=nx*ny*nz
open(1,file='matlab_mesh.vrt')
!open(2,file='matlab.vrt')
!open(2,file=filename(1:11)//'.vtr')
!open(3,file='matlab_force.vrt')
!if (itime.lt.2) then
   write(1,*)nz
   write(1,*)ny
   do i=1,nz
      write(1,*) dz*(i-1)
   enddo

   do i=1,ny
      write(1,*) yp(i)
   enddo
!endif
!do i=1,nz
  !  do j=1,ny
	!write(3,*) Fz(i,j),Fy(i,j)
   ! enddo
!enddo

   !do i=1,nxyz 
      !write(2,*) ux(i,1,1),uy(i,1,1),uz(i,1,1)
   !enddo
close(1)
close(2)
close(3)
return
end subroutine matlab_data

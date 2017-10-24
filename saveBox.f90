!############################################################################
!
subroutine saveBox (ux1,uy1,uz1,pp3,nx1,nx2,ny1,ny2,nz1,nz2)
!
!############################################################################
USE param
USE variables
USE decomp_2d
USE decomp_2d_io
implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,pp3
integer :: code,icomplet,ip,jp,kp,nx1,nx2,ny1,ny2,nz1,nz2
character(len=20) nfichier,nfichier1
character(len=20) :: filename
if (itime==ifirst ) then
call system('rm -rf subSave ')
call system('mkdir subSave')
call system('mkdir ./subSave/ux')
call system('mkdir ./subSave/uy')
call system('mkdir ./subSave/uz')
call system('mkdir ./subSave/pp')
endif
!! ###########################################################################
!! ------------------- UX ---------------------------------------------------
1 format('./subSave/ux/ux',I5.5)

write(filename, 1) (itime)
   call decomp_2d_write_subdomain(1,ux1,nx1,nx2,ny1,ny2,nz1,nz2,filename) 
!! ###########################################################################
!! ------------------- UY ---------------------------------------------------
2 format('./subSave/uy/uy',I5.5)

write(filename, 2) (itime)
call decomp_2d_write_subdomain(1,uy1,nx1,nx2,ny1,ny2,nz1,nz2,filename) 

!! ###########################################################################
!! ------------------- UZ ---------------------------------------------------
3 format('./subSave/uz/uz',I5.5)

write(filename, 3) (itime)
call decomp_2d_write_subdomain(1,uz1,nx1,nx2,ny1,ny2,nz1,nz2,filename) 

!! ###########################################################################
!! ------------------- Prerssure ---------------------------------------------------
4 format('./subSave/pp/pp',I5.5)

write(filename, 4) (itime)
call decomp_2d_write_subdomain(1,pp3,nx1,nx2,ny1,ny2,nz1,nz2,filename) 
!############################################################################
end subroutine saveBox


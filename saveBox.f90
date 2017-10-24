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
integer :: code,icomplet,ip,jp,kp,ifile,nx1,nx2,ny1,ny2,nz1,nz2
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename
if (itime==ifirst ) then
call system('rm -rf boxSave ')
call system('mkdir boxSave')
call system('mkdir ./boxSave/ux')
call system('mkdir ./boxSave/uy')
call system('mkdir ./boxSave/uz')
call system('mkdir ./boxSave/pp')
endif
!! ###########################################################################
!! ------------------- UX ---------------------------------------------------
1 format('./boxSave/ux/ux',I5.5)
ifile=10+nrank
write(filename, 1) (10000*(nrank+1)+itime)
open(ifile,file=filename,ACCESS='STREAM',status='unknown')	
do kp=nz1,nz2
    if (kp.ge.xstart(3).and.kp.le.xend(3)) then
    do jp=ny1,ny2
        if (jp .ge.xstart(2).and.jp.le.xend(2)) then
	    write(ifile) (ux1(i,jp-xstart(2)+1,kp-xstart(3)+1),i=nx1,nx2)
	endif    
    enddo
    endif
enddo
 close(ifile)
!! ###########################################################################
!! ------------------- UY ---------------------------------------------------
2 format('./boxSave/uy/uy',I5.5)
ifile=10+nrank
write(filename, 2) (10000*(nrank+1)+itime)
open(ifile,file=filename,ACCESS='STREAM',status='unknown')	
do kp=nz1,nz2
    if (kp.ge.xstart(3).and.kp.le.xend(3)) then
    do jp=ny1,ny2
        if (jp .ge.xstart(2).and.jp.le.xend(2)) then
	    write(ifile) (uy1(i,jp-xstart(2)+1,kp-xstart(3)+1),i=nx1,nx2)
	endif    
    enddo
    endif
enddo
 close(ifile)

!! ###########################################################################
!! ------------------- UZ ---------------------------------------------------
3 format('./boxSave/uz/uz',I5.5)
ifile=10+nrank
write(filename, 3) (10000*(nrank+1)+itime)
open(ifile,file=filename,ACCESS='STREAM',status='unknown')	
do kp=nz1,nz2
    if (kp.ge.xstart(3).and.kp.le.xend(3)) then
    do jp=ny1,ny2
        if (jp .ge.xstart(2).and.jp.le.xend(2)) then
	    write(ifile) (uz1(i,jp-xstart(2)+1,kp-xstart(3)+1),i=nx1,nx2)
	endif    
    enddo
    endif
enddo
 close(ifile)

!! ###########################################################################
!! ------------------- Prerssure ---------------------------------------------------
4 format('./boxSave/pp/pp',I5.5)
ifile=10+nrank
write(filename, 4) (10000*(nrank+1)+itime)
open(ifile,file=filename,ACCESS='STREAM',status='unknown')	
do kp=nz1,nz2
    if (kp.ge.xstart(3).and.kp.le.xend(3)) then
    do jp=ny1,ny2
        if (jp .ge.xstart(2).and.jp.le.xend(2)) then
	    write(ifile) (pp3(i,jp-xstart(2)+1,kp-xstart(3)+1),i=nx1,nx2)
	endif    
    enddo
    endif
enddo
 close(ifile)
!############################################################################
end subroutine saveBox


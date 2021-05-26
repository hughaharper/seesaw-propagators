c
       subroutine getdipxy(d1,d2,dip,nsx,nsz,nel,x1,y1,x2,y2,F1,F2,F3,
     +			idble,sd1,sd2,sx1,sy1,sx2,sy2,hstepdx,hstepdy,
     +			dxystep3,sF1,sF2,sF3)
c

      parameter (nz=100,nexy=50)
      implicit real*8(a-h,o-z)
      real*8 L
      real*8 sd1(nz),sd2(nz)
      real*4 hstepdx(nexy,nz),hstepdy(nexy,nz),dxystep3(nexy,nz)
      real*4 sx1(nexy,nz),sx2(nexy,nz),sy1(nexy,nz),sy2(nexy,nz)
      real*4 sF1(nexy,nz),sF2(nexy,nz),sF3(nexy,nz)
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
c       create steps (xy and depth) for the force couples
c
c       input:
c
c	d1,d2	- lower and upper fault locking depth (km)
c	dip	- dip angle, defined as 0 for vertical, 
c		  positive counterclockise (deg)
c	nsx,nsz	- number of steps in x and z
c       x1,y1   - start position         (km)
c       x2,y2   - end position           (km)
c	F1,F2,F3 - Force magnitudes
c  	idble - force couple option (if 0, shift by Dtan(dip)
c
c       output:
c
c       sd1,d2	- vector of d1,d2 depths for steps
c	sx1,sx2 - vector of x1,x2 positions for steps
c       sy1,sy2 - vector of y1,y2 positions for steps
c	hstepdx	- horizontal step size (multiply by horizontal sheet)
c	sF1,sF2,sF3 - Force magnitudes (vectoized)
c


c convert dip angle to radians
c dip angle measured from vertical - 90 degree vertical
c  dip should give model dip of 0
	
	dipr=dip*(pi/180)

c compute total X-length of 2-D dipping fault
c need negative angle for dip east (45 deg) for positive xlength

       	Xlength=(d1-d2)*tan(-dipr)

c compute step size in x and z

	dxystep=Xlength/nsx
	dzstep=(d1-d2)/nsz 

c compute dxy step from dip angle
c need negative angle for dip east

	dxystep2=dzstep*tan(-dipr)
	deepslip_shift=d2*tan(-dipr)
	
	print*,'dip = ',dip,' =',dipr
	print*,'Xlength =',Xlength
	print*,'dz step size (km)=',dzstep  
	print*,'dxystep2 (km) = "', dxystep2
	print*,'idble = ',idble
	if(idble.eq.0) then
	print*,'shift by d2tan(dip) for deep slip=',deepslip_shift
	endif
	

c code incorrect for step size less than 2 km, so stop and warn

        if (dzstep.gt.-2) then
      	write(*,'(a)')' vert. step size smaller than 2km - exit!'
	stop
      endif    
 

c compute mapview length of segment

	Dx=(x2-x1)
	Dy=(y2-y1)
	L=sqrt(Dx*Dx + Dy*Dy)

c compute mapview orientation of fault

	rad=pi/180
	theta_0=atan(Dy/Dx)
	theta_d=(90*rad)-abs(theta_0)

	if (theta_0.gt.0) theta_d=-theta_d	

     
c compute step in x and y from dip and orientation

	Dyn=dxystep2*sin(theta_d)
	Dxn=dxystep2*cos(theta_d)
	deepslip_shifty=0
	deepslip_shiftx=0

	
	if(idble.eq.0) then
	deepslip_shifty=d2*tan(-dipr)*sin(theta_d)
	deepslip_shiftx=d2*tan(-dipr)*cos(theta_d)	
	print*,'shiftx by d2tan(dip) for deep slip=',deepslip_shiftx
	print*,'shifty by d2tan(dip) for deep slip=',deepslip_shifty
	endif
	

c loop over number of vertical steps
c nl2 parameter set at 100 - don't let this get bigger

      nsz1=nsz+1
      if (nsz1.ge.100) then
      	write(*,'(a)')' nz (steps) parameter exceeded'
	stop
      endif
      
      
      do 51 ii=1,nsz1
      
      	sF1(nel,ii)=F1
	sF2(nel,ii)=F2
	sF3(nel,ii)=F3

	if(ii.eq.1) then 
		sd1(ii)=d2+(dzstep/2);
		sd2(ii)=d2
		sx1(nel,ii)=x1+deepslip_shiftx
		sx2(nel,ii)=x2+deepslip_shiftx
		sy1(nel,ii)=y1+deepslip_shifty
		sy2(nel,ii)=y2+deepslip_shifty
c		sx1(nel,ii)=x1
c		sx2(nel,ii)=x2
c		sy1(nel,ii)=y1
c		sy2(nel,ii)=y2	
	else if(ii.eq.nsz1) then
		sd1(ii)=sd1(ii-1)+(dzstep/2)
		sd2(ii)=sd1(ii-1)
		sx1(nel,ii)=sx1(nel,ii-1)+Dxn
		sx2(nel,ii)=sx2(nel,ii-1)+Dxn
		sy1(nel,ii)=sy1(nel,ii-1)+Dyn
		sy2(nel,ii)=sy2(nel,ii-1)+Dyn
	else
		sd1(ii)=sd1(ii-1)+dzstep
		sd2(ii)=sd1(ii-1)
		sx1(nel,ii)=sx1(nel,ii-1)+Dxn
		sx2(nel,ii)=sx2(nel,ii-1)+Dxn
		sy1(nel,ii)=sy1(nel,ii-1)+Dyn
		sy2(nel,ii)=sy2(nel,ii-1)+Dyn
	endif

  51  continue


c next calculate horizontal step spacing

      lsxm1=ii-2
      do 52 jj=1,lsxm1
	hstepdx(nel,jj)=abs(sx1(nel,jj+1)-sx1(nel,jj))
	hstepdy(nel,jj)=abs(sy1(nel,jj+1)-sy1(nel,jj))
	dxystep3(nel,jj)=dxystep2
  52  continue
  
      hstepdx(nel,jj)=0.
      hstepdy(nel,jj)=0.






      return
      end




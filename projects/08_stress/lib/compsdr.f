c
      subroutine compsdr(v1,v2,v3,strike1,strike2,dip1,dip2,rake1,rake2)
c
c Compute the strike, dip and rake of two planes of maximum shear 
c stress given the three eigenvectors of the stress tensor.
c These outputs are in Aki and Richards format
c
c  input
c    v1 - 3-component eigenvector
c    v2 - 3-component eigenvector
c    v3 - 3-component eigenvector
c
c  output
c    strike1 - strike is degrees east of north, such that  
c    strike2   the plane is dipping down to the right [0,360]
c    dip1    - dip is degrees below horizontal, to the right
c    dip2      [0,90]
c    rake1   - rake is degrees within the plane up from strike  
c    rake2     direction [-180,180]
c
c       
      implicit real*8(a-h,o-z)
      real*8 v1(3),v2(3),v3(3),nA(3),nB(3),rA(3),rB(3)
      real*8 stvA(3),stvB(3)

      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi
c
      rad=pi/180.
      deg=180.0/pi

c
c find the normal to the two planes of maximum shear stress
c by finding the bisector between v1 and v3
c Then make sure normals are unit vectors
c      
c      nA = v1 + v3
      nA(1) = v1(1) + v3(1) 
      nA(2) = v1(2) + v3(2) 
      nA(3) = v1(3) + v3(3) 

c      nB = v1 - v3
      nB(1) = v1(1) - v3(1) 
      nB(2) = v1(2) - v3(2) 
      nB(3) = v1(3) - v3(3) 
c
c make sure normals are unit vectors
c (this norm should always be sqrt(2) )
c
      rnormA=(nA(1)**2+nA(2)**2+nA(3)**2)**0.5 
      rnormB=(nB(1)**2+nB(2)**2+nB(3)**2)**0.5
c      nA=nA/rnormA
      nA(1)=nA(1)/rnormA
      nA(2)=nA(2)/rnormA
      nA(3)=nA(3)/rnormA

c      nB=nB/rnormB
      nB(1)=nB(1)/rnormB
      nB(2)=nB(2)/rnormB
      nB(3)=nB(3)/rnormB
c
c Calculate mathematical strike, then subtract from 90
c to make it geologic. Then make sure it's between [0,360].
c Finally make sure plane is dipping to the right.
c      
      strikeA = atan2(-nA(1),nA(2))*deg
      strikeA = 90d0-strikeA
      strikeA = mod(strikeA+360d0,360d0)
      if(nA(2)*nA(3).gt.0)then
        strikeA=mod(strikeA+360d0,180d0)+180d0
      else
        strikeA=mod(strikeA+360d0,180d0)
      endif

      strikeB = atan2(-nB(1),nB(2))*deg
      strikeB = 90d0-strikeB
      strikeB = mod(strikeB+360d0,360d0)
      if(nB(2)*nB(3).gt.0)then
        strikeB=mod(strikeB+360d0,180d0)+180d0
      else
        strikeB=mod(strikeB+360d0,180d0)
      endif
c
c calculate dip. make dip always be 0(horizontal) to 90(vertical)
c      
      dipA = acos(abs(nA(3)))*deg
      dipB = acos(abs(nB(3)))*deg
c
c Ensure that normal vector is on the correct side of the plane.
c Ensure that v3 is on the same side as the now-correct normal.
c Make a rake vector (v2 cross normal).
c Ensure the rake vector is pointing in the proper direction (away from v3)
c Make the rake angle (positive if it points up-ish).
c If rake vector is horizontal, angle is 0 or 180 (same or opposite direction as s)
c
      if (sin(strikeA*rad)*nA(2)-cos(strikeA*rad)*nA(1).gt.0)then
        nA(1)=-nA(1)
        nA(2)=-nA(2)
        nA(3)=-nA(3)
      endif
      if (nA(1)*v3(1)+nA(2)*v3(2)+nA(3)*v3(3).lt.0)then
        v3(1)=-v3(1)
        v3(2)=-v3(2)
        v3(3)=-v3(3)
      endif
      rA(1) = nA(2)*v2(3)-nA(3)*v2(2)
      rA(2) = nA(3)*v2(1)-nA(1)*v2(3)
      rA(3) = nA(1)*v2(2)-nA(2)*v2(1)
      if (rA(1)*v3(1)+rA(2)*v3(2)+rA(3)*v3(3).gt.0)then
        rA(1)=-rA(1)
        rA(2)=-rA(2)
        rA(3)=-rA(3)
      endif
      arg=rA(1)*sin(strikeA*rad)+rA(2)*cos(strikeA*rad)
c this is necessary to avoid numerical silly nans
      if(arg.gt.1.0)arg=1.0
      if(arg.lt.-1.0)arg=-1.0
      rakeA=acos(arg)*deg*rA(3)/abs(rA(3))
      if(rA(3).eq.0)then
        rdots=rA(1)*sin(strikeA*rad)+rA(2)*cos(strikeA*rad)
c          print *,rA,strikeA,rdots
c	if(rdots.eq.1.)rakeA=0.
c	if(rdots.eq.-1.)rakeA=180.
	if(abs(rdots-1.).lt.1e-8)rakeA=0.
	if(abs(rdots+1.).lt.1e-8)rakeA=180.
      endif

      if (sin(strikeB*rad)*nB(2)-cos(strikeB*rad)*nB(1).gt.0)then
        nB(1)=-nB(1)
        nB(2)=-nB(2)
        nB(3)=-nB(3)
      endif
      if (nB(1)*v3(1)+nB(2)*v3(2)+nB(3)*v3(3).lt.0)then
        v3(1)=-v3(1)
        v3(2)=-v3(2)
        v3(3)=-v3(3)
      endif
      rB(1) = nB(2)*v2(3)-nB(3)*v2(2)
      rB(2) = nB(3)*v2(1)-nB(1)*v2(3)
      rB(3) = nB(1)*v2(2)-nB(2)*v2(1)
      if (rB(1)*v3(1)+rB(2)*v3(2)+rB(3)*v3(3).gt.0)then
        rB(1)=-rB(1)
        rB(2)=-rB(2)
        rB(3)=-rB(3)
      endif
      arg=rB(1)*sin(strikeB*rad)+rB(2)*cos(strikeB*rad)
      if(arg.gt.1.0)arg=1.0
      if(arg.lt.-1.0)arg=-1.0
      rakeB=acos(arg)*deg*rB(3)/abs(rB(3))
      if(rB(3).eq.0)then
        rdots=rB(1)*sin(strikeB*rad)+rB(2)*cos(strikeB*rad)
c	if(rdots.eq.1)rakeB=0.
c	if(rdots.eq.-1)rakeB=180.
	if(abs(rdots-1.).lt.1e-8)rakeB=0.
	if(abs(rdots+1.).lt.1e-8)rakeB=180.
      endif
c
c the smaller strike becomes strike1, the larger strike2
c      
      if (strikeA .lt. strikeB) then
        strike1 = strikeA
	strike2 = strikeB
	dip1 = dipA
	dip2 = dipB
	rake1 = rakeA
	rake2 = rakeB
      else
        strike1 = strikeB
	strike2 = strikeA
	dip1 = dipB
	dip2 = dipA
	rake1 = rakeB
	rake2 = rakeA
      endif

      return
      end
c

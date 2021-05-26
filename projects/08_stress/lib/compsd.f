c
      subroutine compsd(v1,v2,v3,strike1,strike2,dip1,dip2)
c
c Compute the strike and dip of two planes of maximum shear 
c stress given the three eigenvectors of the stress tensor
c
c  input
c    v1 - 3-component eigenvector
c    v2 - 3-component eigenvector
c    v3 - 3-component eigenvector
c
c  output
c    strike1 - strike of first plane of maximum shear stress
c    strike2 - strike of second plane of maximum shear stress
c    dip1    - dip of first plane of maximum shear stress
c    dip2    - dip of second plane of maximum shear stress
c       
      implicit real*8(a-h,o-z)
      real*8 v1(3),v2(3),v3(3),nA(3),nB(3)
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
c calculate strike. old strike will be [-90,90] (mathematical)
c new strike will be [0,180] degrees CW from North (geologic)
c      
      strikeA = atan2(-nA(1),nA(2))*deg
      strikeB = atan2(-nB(1),nB(2))*deg
      if(strikeA.lt.0)strikeA=strikeA+180.
      if(strikeB.lt.0)strikeB=strikeB+180.
      strikeA = 90d0-strikeA
      strikeB = 90d0-strikeB
      if(strikeA.lt.0)strikeA=strikeA+180.
      if(strikeB.lt.0)strikeB=strikeB+180.
c
c calculate dip. 
c vertical=90, dipping east=[0,90], dipping west=[90,180]
c      
      dipA = acos(nA(1)*nA(3)/abs(nA(1)))*deg
      dipB = acos(nB(1)*nB(3)/abs(nB(1)))*deg
c
c the smaller strike becomes strike1, the larger strike2
c      
      if (strikeA .lt. strikeB) then
        strike1 = strikeA
	strike2 = strikeB
	dip1 = dipA
	dip2 = dipB
      else
        strike1 = strikeB
	strike2 = strikeA
	dip1 = dipB
	dip2 = dipA
      endif

      return
      end
c

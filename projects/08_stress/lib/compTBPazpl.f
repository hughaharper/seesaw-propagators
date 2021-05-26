c
      subroutine compTBPazpl(v1,v2,v3,azimuth1,azimuth2,
     +                       azimuth3,plunge1,plunge2,plunge3)
c
c Compute the azimuth and plunge of three principal stress axes. 
c Plunge is degrees below horizontal [0,90]
c Azimuth is degrees EofN [0,360]
c
c  input
c    v1 - 3-component eigenvector
c    v2 - 3-component eigenvector
c    v3 - 3-component eigenvector
c
c  output
c    az1 - azimuth EofN
c    az2  
c    az3  
c    pl1 - plunge below horizontal
c    pl2  
c    pl3  
c
c       
      implicit real*8(a-h,o-z)
      real*8 v1(3),v2(3),v3(3)

      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi
c
      rad=pi/180.
      deg=180.0/pi
c
c make sure vectors are pointing down
c
      if (v1(3).gt.0)then
        v1(1)=-v1(1)
        v1(2)=-v1(2)
        v1(3)=-v1(3)
      endif
      if (v2(3).gt.0)then
        v2(1)=-v2(1)
        v2(2)=-v2(2)
        v2(3)=-v2(3)
      endif
      if (v3(3).gt.0)then
        v3(1)=-v3(1)
        v3(2)=-v3(2)
        v3(3)=-v3(3)
      endif
c
c azimuth is atan(x/y), convert to degrees, make between 0 and 360
c
      azimuth1 = mod(atan2(v1(1),v1(2))*deg+360d0,360d0)
      azimuth2 = mod(atan2(v2(1),v2(2))*deg+360d0,360d0)
      azimuth3 = mod(atan2(v3(1),v3(2))*deg+360d0,360d0)
c
c plunge is asin(abs(z), in degrees
c
      plunge1=asin(abs(v1(3)))*deg
      plunge2=asin(abs(v2(3)))*deg
      plunge3=asin(abs(v3(3)))*deg

      return
      end
c

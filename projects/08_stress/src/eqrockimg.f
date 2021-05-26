	program eqrockimg
c
c*****************   main program   ***************************************
c Program to transform  an img file by reading it in a strip at a time, and
c FFTing each strip.  Strips are overlapped by half their width and 
c cosine blended.  N and S edges are not obscured by tapering because 
c mirrors are used for first and last strip.
c
c These parameters must be set at compile time:
c
c nximg - the number of columns in the imgfile;  for 2' img use 10800 
c nyimg - the number of rows in the imgfile;     for 2' img use 6336
c nyfft - the number of rows in the FFT;         for 2' img use 576 or 704
c
c nyimg / nyfft must divide evenly with no remainder. 
c
c The shortest distance is the y-distance, and it is smallest in the
c first and last strips, where the scale is set for 72 degrees latitude.
c At this latitude, each pixel is 1.145 km across.  Thus the longest
c wavelength resolved in the FFT will be 1.145 * nyfft km, or:
c if nyfft=576, this is 660 km, or harmonic degree 60;
c if nyfft=704, this is 806 km, or harmonic degree 50.
c This isn't well-resolved; it is simply the 1 cycle wavelength in
c the y direction.  It isn't fully resolved because of the cosine window.
c
c Topography is input.  Stress at depth zobs resulting from the vertical
c loads of topography and airy-deformed moho is calculated.  Gravity is
c also calculated.  Cartesian stress components of full 3-D stress tensor
c are output as 7 img files.
c**************************************************************************
c
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
      real*8 kx,ky,kxm,kym
c
c  change nximg, nyimg, and nyfft as described above
c
c      parameter (nximg=10800, nyimg=6336, nyfft=704, drln=2./60.)
      parameter (nximg=21600, nyimg=17280, nyfft=1920, drln=1./60.)
      character*80 ifilnam,ofilnam1,crhoc
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
      real*4 fz(nximg,nyfft),ftopo(nximg,nyfft),grav(nximg,nyfft)
      integer*4 n(2)
      integer*4 ilandorsea
      logical exists
c
      pi=acos(-1.)
c
      n(1)=nximg
      n(2)=nyfft
      ni=nyfft
      ni2=nyfft/2+1
      nj=nximg
      nj2=nximg/2+1
      ntot=ni*nj
      nyfft2=nyfft/2
c
c  get values from command line
c
      narg = iargc()
      if(narg.lt.3) then
        write(*,'(a)')'  '
        write(*,'(a)')
     &  'Usage: eqrock topo.img topo.eqrk.img'
        write(*,'(a)')
        write(*,'(a)')
     &  '       topo.img - 2-byte integer mercator projected file'
        write(*,'(a)')
     &  '       topo.eqrk.img - output equivalent rock topography file'
        write(*,'(a)')
     &  '       rhoc     - crustal density [kg/m^3]'
        write(*,'(a)')
        stop
      else
        call getarg(1,ifilnam)
        nc=index(ifilnam,' ')
        ifilnam(nc:nc)=char(0)
	call getarg(2,ofilnam1)
        nc=index(ofilnam1,' ')
        ofilnam1(nc:nc)=char(0)
        call getarg(3,crhoc)
        read(crhoc,*)rhoc
c
c check if any output file exists
c        
  60    format('output file exists, remove it before continuing')
        inquire(file=ofilnam1, exist=exists)
        if (exists) then
          write(*,60)
          stop
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c set parameters 
c
c      rhoc=2800.
      rhom=3300.
      rhow=1025.
c
c Here is the main loop:
c ikillm=1 means subtract the mean upon read. =0 means don't.
c
      nstrips = nyimg/nyfft
      ikillm = 0
      jrowseek = 0
      jlatscl = 0
      do 901 istrip = 1, nstrips
c
c read one strip
c
        call fimgr(ftopo, nximg, nyfft, jrowseek, ikillm, ifilnam)
        jrowseek=jrowseek+nyfft
c
c what is the lat and lon of each point
c
        do 120 i=1,ni
          call getscal(jlatscl,nximg,nyimg,nyfft,width,height,rltrad)
	  rlt=rltrad*180./pi	    
	  jlatscl=jlatscl+1
        do 120 j=1,nj
          rln=(j-1)*drln
	  if (rln.gt.180.) rln=rln-360.
c
c ilandorsea = 1 for land, 0 for water
c
	  call readsrtmmask(rln,rlt,ilandorsea)
          if(ilandorsea.eq.0) then
	    feqrk=(rhoc-rhow)/rhoc
	  else
	    feqrk=1.
	  end if

	  ftopo(j,i)=ftopo(j,i)*feqrk
 120    continue
c
c
c write one strip
c
        call fimgw(ftopo,nximg,nyfft,ofilnam1)
c
c End of main loop
c
 901  continue
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getscal(jlatscl,nximg,nyimg,nyfft,xwide,ytall,drlat)
c
c Given nximg, the number of pixels in 360 degrees of longitude,
c and also the number of pixels in the horizontal direction of the FFT
c nyimg, such that nyimg/2 is the origin of the counter jlatscl,
c nyfft, the number of pixels in the vertical direction of the FFT,
c return xwide, ytall, and drlat, the width and height of the FFT in
c km and the double precision radian latitude for this scale factor.
c
      integer*4 jlatscl, nximg, nyimg, nyfft
      real*8 xwide, ytall
      real*8 drlat
c
c radius is the radius of the img mercator projection, that is, it is
c the number of pixels on the Equator per radian of longitude.
c twopi is 2 pi, and halfpi is pi/2; y is the argument to the 
c Gudermannian function, and pixelkm is the width of a pixel in km at
c the latitude drlat
c
      real*8 radius, twopi, halfpi, y, pixelkm
c
      halfpi = 2.0d+00 * datan(1.0d+00)
      twopi = 4.0d+00 * halfpi
      radius = nximg / twopi
      y = (nyimg/2 - jlatscl)/radius
      drlat = 2.0d+00 * datan(dexp(y)) - halfpi
      pixelkm = (360.0d+00 * 111.1949 * dcos(drlat))/nximg
c
      xwide = nximg * pixelkm
      ytall = nyfft * pixelkm
      return
      end

	program topo_stress
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
      character*100 ifilnam,ofilnam1,ofilnam2,ofilnam3
      character*80 ofilnam4,ofilnam5,ofilnam6
      character*80 ch,cTe,cz,crhoc
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
      real*4 savex(nximg, nyfft), blend(nximg, nyfft/2), window(nyfft)
      real*4 blend1(nximg, nyfft/2),blend2(nximg, nyfft/2)
      real*4 blend3(nximg, nyfft/2),blend4(nximg, nyfft/2)
      real*4 blend5(nximg, nyfft/2),blend6(nximg, nyfft/2)
      real*4 fz(nximg,nyfft)
      real*4 Txx(nximg,nyfft),Tyy(nximg,nyfft),Tzz(nximg,nyfft)
      real*4 Txy(nximg,nyfft),Txz(nximg,nyfft),Tyz(nximg,nyfft)
      complex*8 savek(nximg/2+1,nyfft), work(nximg)
      complex*8 fzk(nximg/2+1,nyfft)
      complex*8 Txxk(nximg/2+1,nyfft),Tyyk(nximg/2+1,nyfft)
      complex*8 Tzzk(nximg/2+1,nyfft),Txyk(nximg/2+1,nyfft)
      complex*8 Txzk(nximg/2+1,nyfft),Tyzk(nximg/2+1,nyfft)
      complex*8 fzval,gzval,tkval,grval
      equivalence (savex(1,1), savek(1,1))
      equivalence (fz(1,1),fzk(1,1))
      equivalence (Txx(1,1),Txxk(1,1))
      equivalence (Tyy(1,1),Tyyk(1,1))
      equivalence (Tzz(1,1),Tzzk(1,1))
      equivalence (Txy(1,1),Txyk(1,1))
      equivalence (Txz(1,1),Txzk(1,1))
      equivalence (Tyz(1,1),Tyzk(1,1))
      integer*4 n(2)
      integer*4 ilandorsea
      logical exists
c
      pi=acos(-1.)
c
      n(1)=nximg
      n(2)=nyfft
      ikillm = 1
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
      if(narg.lt.11) then
        write(*,'(a)')'  '
        write(*,'(a)')
     &  'Usage: topo_stress topo.img H Te rhoc Zobs outfiles.img (6)'
        write(*,'(a)')
        write(*,'(a)')
     &  '       topo.img - 2-byte integer mercator projected eqrk topo'
        write(*,'(a)')
     &  '       H        - depth to bottom of elastic layer [km]'
        write(*,'(a)')
     &  '       Te       - elastic thickness [km]'
        write(*,'(a)')
     &  '       rhoc     - crustal density [kg/m^3]'
        write(*,'(a)')
     &  '       Zobs     - stress observation depth < 0 [km]'
        write(*,'(a)')
     &  '       out.img  - 2-byte integer mercator projected files'
        write(*,'(a)')
     &  '                  Txx.img, Tyy.img, Tzz.img,'
        write(*,'(a)')
     &  '                  Txy.img, Txz.img, Tyz.img'
        write(*,'(a)')
        stop
      else
        call getarg(1,ifilnam)
        nc=index(ifilnam,' ')
        ifilnam(nc:nc)=char(0)
        call getarg(2,ch)
        call getarg(3,cTe)
        call getarg(4,crhoc)
        call getarg(5,cz)
        read(ch,*)thk
        read(cTe,*)Te
        read(crhoc,*)rhoc
        read(cz,*)zobs
	  call getarg(6,ofilnam1)
        nc=index(ofilnam1,' ')
        ofilnam1(nc:nc)=char(0)
	  call getarg(7,ofilnam2)
        nc=index(ofilnam2,' ')
        ofilnam2(nc:nc)=char(0)
	  call getarg(8,ofilnam3)
        nc=index(ofilnam3,' ')
        ofilnam3(nc:nc)=char(0)
	  call getarg(9,ofilnam4)
        nc=index(ofilnam4,' ')
        ofilnam4(nc:nc)=char(0)
	  call getarg(10,ofilnam5)
        nc=index(ofilnam5,' ')
        ofilnam5(nc:nc)=char(0)
	  call getarg(11,ofilnam6)
        nc=index(ofilnam6,' ')
        ofilnam6(nc:nc)=char(0)
c
c check if any output file exists
c        
  60    format('output file exists, remove it before continuing')
        inquire(file=ofilnam1, exist=exists)
        if (exists) then
          write(*,60)
          stop
        endif
        inquire(file=ofilnam2, exist=exists)
        if (exists) then
          write(*,60)
          stop
        endif
        inquire(file=ofilnam3, exist=exists)
        if (exists) then
          write(*,60)
          stop
        endif
        inquire(file=ofilnam4, exist=exists)
        if (exists) then
          write(*,60)
          stop
        endif
        inquire(file=ofilnam5, exist=exists)
        if (exists) then
          write(*,60)
          stop
        endif
        inquire(file=ofilnam6, exist=exists)
        if (exists) then
          write(*,60)
          stop
        endif
      endif
c
c check that nyimg and nyfft parameters are compatible, and if so,
c initialize nstrips, nyfft2, and window array:
c
      istrip = nyimg/nyfft
      if (istrip * nyfft .ne. nyimg) then
        write(*, 101)
 101    format ('topo_stress:  ERROR.  Incompatible parameters.')
        stop
      end if
      nstrips = 2 * istrip + 1
      winfact = pi/nyfft2
      do 201 j = 1, nyfft
        window(j)=(0.5 * (1.0 + dcos((j - nyfft2 - 0.5)*winfact)))
 201  continue
c
c set parameters for the common block and parameters to be passed
c (Young's modulus = 70 Gpa and Poisson's ratio = .25)
c
      grv=9.81
      Gcns=2.*pi*6.673e-11
c      rhoc=2800.
c      rhoc=2500.
      rhom=3300.
      rhow=1025.
      
      young=7.e10
c      rnu=0.499
      rnu=0.25
	  
      rlam1=young*rnu/(1+rnu)/(1-2*rnu)
      rmu1=young/2/(1+rnu)
      bulk=young/3/(1-2*rnu)
c      rlam1=2*young/5.
c      rmu1=rlam1
c      bulk=rlam1+2*rmu1/3.
      
      H=-abs(thk*1000.)
      zobs=-abs(zobs*1000.)
      Te=abs(Te*1000.)
c
c Here is the main loop:
c
      jrowseek = 0
      jlatscl = 0
      write (*,"(a,i2)") "number of strips =", nstrips
      do 901 istrip = 1, nstrips
c
        write (*,"(a,i2)") "do strip =", istrip
c Load the strip:
c
        if (istrip .gt. 2 .and. istrip .lt. nstrips) then
          jrowseek = jrowseek + nyfft2
        end if
        call fimgr(savex, nximg, nyfft, jrowseek, ikillm, ifilnam)
c
c Window the strip.  This reverses the window and mirrors if
c istrip .ew. 1 .or. istrip .eq. nstrips
c
        call winstrp(savex, nximg, nyfft, istrip, nstrips, window)
c
c Now figure out what latitude we are at and what the array size in km is.
c
        if (istrip .gt. 1) then
          jlatscl = jlatscl + nyfft2
        end if
        call getscal(jlatscl,nximg,nyimg,nyfft,width,height,drlat)
        write(*, 219) istrip, 180.0*drlat/pi, width, height
 219    format('FFT ',i2,' Lat ',f8.4,' is ',f7.1,' km by ',f7.1,' km') 
c
c define the topo and topo loads
c ilandorsea = 1 for land, 0 for water
c
        do 120 i=1,ni
	        if (istrip.eq.1) then
	          if (i.le.nyfft2) then
	            ipass=i-1
	          else
	            ipass=nyfft-i+1
	          end if
            call getscal(ipass,nximg,nyimg,nyfft,width,height,rltrad)
	          rlt=rltrad*180./pi	    
	        elseif (istrip.lt.nstrips) then
	          ipass=(istrip-2)*nyfft2+i-1
            call getscal(ipass,nximg,nyimg,nyfft,width,height,rltrad)
	          rlt=rltrad*180./pi
	        else
	          if (i.le.nyfft2) then
	            ipass=(istrip-1)*nyfft2-i+1
	          else
	            ipass=(istrip-3)*nyfft2+i-1
	          end if
            call getscal(ipass,nximg,nyimg,nyfft,width,height,rltrad)
	          rlt=rltrad*180./pi	    
	        end if
          do 121 j=1,nj
            rln=(j-1)*drln
	          if (rln.gt.180.) rln=rln-360.

c topo should be eqrk and highpassed upon input

	          feqrk=1.

	          fz(j,i)=savex(j,i)*feqrk*grv*rhoc/(ni*nj)
 121      continue       
 120    continue
c
c for debugging: print out the topography... does it make sense?
c
c        do 122 i=1,ni
c	        do 123 j=1,nj
c	          Txx(j,i)=savex(j,i)
c 123      continue
c 122    continue
c        
c	      scl=1
c        call scblwr(Txx,scl,ofilnam1,blend1,istrip,nstrips,nximg,nyfft)
c
c to skip the FFT part for debugging...
c
c        go to 700

c
c take the fourier transform
c 
        call fourt(fz,n,2,-1,0,work,nximg)
c
c Begin Wavnumber Loop: do something in the fourier domain
c
        do 255 i=1,ni
          ky=-(i-1)/height
          if(i.ge.ni2) ky= (ni-i+1)/height
          do 256 j=1,nj2
            kx=(j-1)/width

c make wavenumbers in meters for passing to subroutines
            kxm=kx/1000.
	          kym=ky/1000.

c compute moho topo
            fzval=fzk(j,i)
            call mohotopo(kxm,kym,Te,rhoc,rhom,fzval,gzval)

c compute the displacements      
            call airy2load(kxm,kym,zobs,H,fzval,gzval,cu,cv,cw
     +                                     ,cdwdz,cdudz,cdvdz)
c
c compute the cartesian stresses      
            call disp2stress(kxm,kym,cu,cv,cw,cdudz,cdvdz,cdwdz,
     +                       cTxx,cTyy,cTzz,cTxy,cTxz,cTyz)
c
c store the cartesian stresses for IFFT
c
	          Txxk(j,i)=cTxx
	          Tyyk(j,i)=cTyy
	          Tzzk(j,i)=cTzz
	          Txyk(j,i)=cTxy
	          Txzk(j,i)=cTxz
	          Tyzk(j,i)=cTyz

 256      continue
 255    continue
        call fourt(Txx,n,2,1,-1,work,nximg)
        call fourt(Tyy,n,2,1,-1,work,nximg)
        call fourt(Tzz,n,2,1,-1,work,nximg)
        call fourt(Txy,n,2,1,-1,work,nximg)
        call fourt(Txz,n,2,1,-1,work,nximg)
        call fourt(Tyz,n,2,1,-1,work,nximg)
c
c scale, blend, and write the files
c
        scl=.00001
        call scblwr(Txx,scl,ofilnam1,blend1,istrip,nstrips,nximg,nyfft)
        call scblwr(Tyy,scl,ofilnam2,blend2,istrip,nstrips,nximg,nyfft)
        call scblwr(Tzz,scl,ofilnam3,blend3,istrip,nstrips,nximg,nyfft)
        call scblwr(Txy,scl,ofilnam4,blend4,istrip,nstrips,nximg,nyfft)
        call scblwr(Txz,scl,ofilnam5,blend5,istrip,nstrips,nximg,nyfft)
        call scblwr(Tyz,scl,ofilnam6,blend6,istrip,nstrips,nximg,nyfft)

c 700    continue
c
c End of main loop
c
 901  continue
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine winstrp(savex, nximg, nyfft, istrip, nstrips, window)
c
c Window the data in savex with the cosine window stored in window
c Except when istrip .eq. 1 .or. istrip .eq. nstrips, then do a
c special windowing with a mirror.
c
      real*4 savex(nximg,1), window(nyfft)
      integer*4 istrip, nstrips
      integer*4 i, j, nyfft2, jw, jmirror
c
      nyfft2 = nyfft/2
c
      if (istrip .eq. 1) then
c
        do 702 j = 1, nyfft2
          jmirror = nyfft + 1 - j
          jw = j + nyfft2
          do 701 i = 1, nximg
            savex(i,j) = savex(i,j) * window(jw)
            savex(i,jmirror) = savex(i,j)
 701      continue
 702    continue
c
      else if (istrip .ne. nstrips) then
c
        do 752 j = 1, nyfft
          do 751 i = 1, nximg
            savex(i,j) = savex(i,j) * window(j)
 751      continue
 752    continue
c
      else
c
        do 792 j = 1, nyfft2
          jmirror = nyfft2 + 1 - j
          jw = j + nyfft2
          do 791 i = 1, nximg
            savex(i,jw) = savex(i,jw) * window(j)
            savex(i,jmirror) = savex(i,jw)
 791      continue
 792    continue
c
      end if
      return
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scblwr(savex,scl,ofilnam,blend,
     +                istrip,nstrips,nximg,nyfft)
c
c Scale, Blend, and Write.
c
      real*4 savex(nximg, nyfft), blend(nximg, nyfft/2)
      real*8 scl
      character*80 ofilnam
      integer*4 istrip,nstrips,nximg,nyfft
      integer*4 i,j,nyfft2


      nyfft2 = nyfft/2
c
c
        if (istrip .eq. 1) then
c
c copy the first half of the array to the blend strip.
c
          do 302 j = 1, nyfft2
            do 301 i = 1, nximg
              blend(i,j) = savex(i,j)*scl
 301        continue
 302      continue
c
        else if (istrip .ne. nstrips) then
c
c Add the first half to the blend strip; write out, then copy 2nd half.
c
          do 352 j = 1, nyfft2
            do 351 i = 1, nximg
              blend(i,j) = blend(i,j) + savex(i,j)*scl
 351        continue
 352      continue
c
          call fimgw(blend, nximg, nyfft2, ofilnam)
c
          do 362 j = 1, nyfft2
            do 361 i = 1, nximg
              blend(i,j) = savex(i,j+nyfft2)*scl
 361        continue
 362      continue
c
        else
c
c Add the second half to the blend strip and write out
c
          do 392 j = 1, nyfft2
            do 391 i = 1, nximg
              blend(i,j) = blend(i,j) + savex(i,j+nyfft2)*scl
 391        continue
 392      continue
c
          call fimgw(blend, nximg, nyfft2, ofilnam)
c
        end if
      return
      end

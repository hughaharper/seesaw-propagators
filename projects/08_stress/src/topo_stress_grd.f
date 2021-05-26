c
      program topo_stress_grd
c
c*****************   main program   ***************************************
c
c This program is intended to be a fast, local, grd version of the program
c /mgg/fftfault/src/topo_stress/topo_stress.f.  Should call all the same
c subroutines.  Only difference is the calculation is done over a small
c region instead of globally
c
c Topography is input.  Stress at depth zobs resulting from the vertical
c loads of topography and airy-deformed moho is calculated.  Gravity is
c also calculated.  Cartesian stress components of full 3-D stress tensor
c are output as 7 grd files.
c
c**************************************************************************
c
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
      real*8 kx,ky,kh2,kh,betz
      real*8 rln0,rlt0,dlt,dln,ddx,ddy,rland,rdum
      real*8 rln0g,rlt0g,dltg,dlng,trans
c      complex*8 gsurf,gmoho,grval
c
c  change ni and nj as needed
c
c      parameter(ni=1024,nj=1024,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
      parameter(ni=4096,nj=4096,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
c      parameter(ni=2000,nj=3000,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
      character*80 ifilnam,ofilnam1,ofilnam2,ofilnam3
      character*80 ofilnam4,ofilnam5,ofilnam6,ofilnam7
      character*80 ch,cTe,cz,crhoc,caobs,title
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
      real*4 fz(nj,ni),fpeek(nj,ni),topo(nj,ni),ftopo(nj,ni)
c      real*4 gz(nj,ni),gpeek(nj,ni),gtopo(nj,ni)
      real*4 Txx(nj,ni),Tyy(nj,ni),Tzz(nj,ni)
      real*4 Txy(nj,ni),Txz(nj,ni),Tyz(nj,ni)
      real*4 grav(nj,ni)
      real*4 xwind(nj),ywind(ni)
      complex*8 fzk(nj2,ni),ftopok(nj2,ni)
      complex*8 Txxk(nj2,ni),Tyyk(nj2,ni),Tzzk(nj2,ni)
      complex*8 Txyk(nj2,ni),Txzk(nj2,ni),Tyzk(nj2,ni)
      complex*8 gravk(nj2,ni)
      complex*8 fzval,gzval,tkval,grval
      dimension n(2)
      complex*8 work(nwork)
      equivalence (fz(1,1),fzk(1,1))
      equivalence (ftopo(1,1),ftopok(1,1))
      equivalence (grav(1,1),gravk(1,1))
      equivalence (Txx(1,1),Txxk(1,1))
      equivalence (Tyy(1,1),Tyyk(1,1))
      equivalence (Tzz(1,1),Tzzk(1,1))
      equivalence (Txy(1,1),Txyk(1,1))
      equivalence (Txz(1,1),Txzk(1,1))
      equivalence (Tyz(1,1),Tyzk(1,1))
c
      pi=acos(-1.)
c
c  set the dimensions for fourt
c
      n(1)=nj
      n(2)=ni
c
c   get values from command line
c
      narg = iargc()
      if(narg.lt.11) then
        write(*,'(a)')'  '
        write(*,'(a)')
     &  'Usage: topo_stress_grd topo.grd H Te rhoc Zobs out.grd (6)'
        write(*,'(a)')
        write(*,'(a)')
     &  '       topo.grd - high-passed eqivalent rock topo [m]'
        write(*,'(a)')
     &  '       H        - depth to bottom of elastic layer [km]'
        write(*,'(a)')
     &  '       Te       - elastic thickness [km]'
        write(*,'(a)')
     &  '       rhoc     - crustal density [kg/m^3]'
        write(*,'(a)')
     &  '       Zobs     - stress observation depth < 0 [km]'
        write(*,'(a)')
     &  '       out.grd  - stress component and gravity grd files'
        write(*,'(a)')
     &  '                  Txx.grd, Tyy.grd, Tzz.grd,'
        write(*,'(a)')
     &  '                  Txy.grd, Txz.grd, Tyz.grd'
        write(*,'(a)')
        stop
      else 
        call getarg(1,ifilnam)
c        nc=index(ifilnam,' ')
c        ifilnam(nc:nc)=char(0)
        call getarg(2,ch)
        call getarg(3,cTe)
        call getarg(4,crhoc)
        call getarg(5,cz)
        read(ch,*)thk
        read(cTe,*)Te
        read(crhoc,*)rhoc
        read(cz,*)zobs
	call getarg(6,ofilnam1)
c        nc=index(ofilnam1,' ')
c        ofilnam1(nc:nc)=char(0)
	call getarg(7,ofilnam2)
c        nc=index(ofilnam2,' ')
c        ofilnam2(nc:nc)=char(0)
	call getarg(8,ofilnam3)
c        nc=index(ofilnam3,' ')
c        ofilnam3(nc:nc)=char(0)
	call getarg(9,ofilnam4)
c        nc=index(ofilnam4,' ')
c        ofilnam4(nc:nc)=char(0)
	call getarg(10,ofilnam5)
c        nc=index(ofilnam5,' ')
c        ofilnam5(nc:nc)=char(0)
	call getarg(11,ofilnam6)
c        nc=index(ofilnam6,' ')
c        ofilnam6(nc:nc)=char(0)
      endif
c
c set Young's modulus, Poisson's ratio
c calculate the other parameters
c
      grv=9.81
      young=7.e10
c      rnu=0.499
      rnu=0.25
      
      rlam1=young*rnu/(1+rnu)/(1-2*rnu)
      rmu1=young/2/(1+rnu)
      bulk=young/3/(1-2*rnu)
      
      H=-abs(thk*1000.)
      Te=abs(Te*1000.)
      zobs=-abs(zobs*1000.)
      
c      rhoc=2800.0000000000000000000
      rhom=3300.0000000000000000000
      rhow=1025.0000000000000000000

c
c  define the Flexural Rigidity
c
      D=young*Te**3/12./(1-rnu**2)
c      write(*,*)rnu1,D
c
c  allow rmu2 to vary but don't let the bulk
c  modulus change
c
      rmu2=rmu1*shr
      rlam2=bulk-2.*rmu2/3.
c
c   read the grd files
c
      call readgrd(topo,nj1,ni1,rlt0,rln0,
     +            dlt,dln,rdum,title,trim(ifilnam)//char(0))
      if(ni1.ne.ni.or.nj1.ne.nj) then
        print *,ni1,ni,nj1,nj,rlt0,rln0,dlt,dln,rdum
        write(*,'(a)')' recompile program to match topo size'
        stop
      endif
c
c  compute the windows
c
      nsigy=ni/8
      do 70 i=1,ni
      if(i.lt.nsigy) then
       ywind(i)=0.5*(1.-cos(pi*(i-1)/nsigy))
      else if(i.gt.(ni-nsigy)) then
       ywind(i)=0.5*(1.-cos(pi*(ni-i)/nsigy))
      else
       ywind(i)=1.
      endif
   70 continue
      nsigx=nj/8
      do 80 j=1,nj
      if(j.lt.nsigx) then
       xwind(j)=0.5*(1.-cos(pi*(j-1)/nsigx))
      else if(j.gt.(nj-nsigx)) then
       xwind(j)=0.5*(1.-cos(pi*(nj-j)/nsigx))
      else
       xwind(j)=1.
      endif
   80 continue
c              
c  compute the height and width of the area in km
c
      rlat2=abs(rlt0+ni*dlt/2)*pi/180.
      xscl=cos(rlat2)
      dy=111000.*dlt
      dx=xscl*111000.*dln
      width=nj*dx
      height=abs(ni*dy)
c
c  apply window to input topo and define load
c
      do 120 i=1,ni
      do 120 j=1,nj
       fz(j,i)=topo(j,i)*grv*rhoc*xwind(j)*ywind(i)/(ni*nj)
c       ftopo(j,i)=topo(j,i)*xwind(j)*ywind(i)/(ni*nj)
c       fpeek(j,i)=fz(j,i)
  120 continue
c
c  take the fourier transform of the input
c
      call fourt(fz,n,2,-1,0,work,nwork)
c      call fourt(ftopo,n,2,-1,0,work,nximg)
c
c Begin Wavnumber Loop: do something in the fourier domain
c
      do 255 i=1,ni
      ky=-(i-1)/height
      if(i.ge.ni2) ky= (ni-i+1)/height
      do 255 j=1,nj2
      kx=(j-1)/width

c
c compute moho topo
c
          fzval=fzk(j,i)
          call mohotopo(kx,ky,Te,rhoc,rhom,fzval,gzval)
c
c compute the displacements
c      
          call airy2load(kx,ky,zobs,H,fzval,gzval,cu,cv,cw
     +                                     ,cdwdz,cdudz,cdvdz)
c
c compute the cartesian stresses
c      
          call disp2stress(kx,ky,cu,cv,cw,cdudz,cdvdz,cdwdz,
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

 255  continue
c
c  do the inverse fft's
c
      call fourt(Txx,n,2,1,-1,work,nwork)
      call fourt(Tyy,n,2,1,-1,work,nwork)
      call fourt(Tzz,n,2,1,-1,work,nwork)
      call fourt(Txy,n,2,1,-1,work,nwork)
      call fourt(Txz,n,2,1,-1,work,nwork)
      call fourt(Tyz,n,2,1,-1,work,nwork)
c
c  write grd file, some parameters must be real*8
c
      ddx=dlt
      ddy=dln
      rland=9998.
      rdum=9999.
      if(istr.ge.1) then
       rland=9998.d0*rmu1
       rdum=9999.d0*rmu1
      endif
      call writegrd(Txx,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ofilnam1)//char(0),trim(ofilnam1)//char(0))
      call writegrd(Tyy,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ofilnam2)//char(0),trim(ofilnam2)//char(0))
      call writegrd(Tzz,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ofilnam3)//char(0),trim(ofilnam3)//char(0))
      call writegrd(Txy,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ofilnam4)//char(0),trim(ofilnam4)//char(0))
      call writegrd(Txz,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ofilnam5)//char(0),trim(ofilnam5)//char(0))
      call writegrd(Tyz,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ofilnam6)//char(0),trim(ofilnam6)//char(0))
      stop
      end

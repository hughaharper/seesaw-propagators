c
      program predictgrav
c
c*****************   main program   ************************
c
c  Predicts free air gravity from a topography grid, or
c  upward continues a gravity file to the appropriate 
c  altitude.  Predicted gravity includes linear terms only 
c  of Parker's gravity expansion, for both contribution from
c  surface topography and contribution from Moho topography.
c  Moho topography is calculated assuming elastic plate 
c  flexure with given elastic thickness.
c
c***********************************************************
c
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
      real*8 kx,ky,kh2,kh,betz
      real*8 rln0,rlt0,dlt,dln,ddx,ddy,rland,rdum
      real*8 rln0g,rlt0g,dltg,dlng,trans
      complex*8 gsurf,gmoho,grval
c
c  change ni and nj as needed
c
c     parameter(ni=1024,nj=1024,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
      parameter(ni=1080,nj=1080,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
c      parameter(ni=3000,nj=1800,nwork=32768,nj2=nj/2+1,ni2=ni/2+1)
      character*80 fxdisp,fydisp,fzdisp,cTe
      character*80 ctopo,cgtopo,ch,cz,crho,cshr,ccomp,cstr,title
      character*80 czs,caobs,citg,crhoc
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
      real*4 fz(nj,ni),fpeek(nj,ni),topo(nj,ni)
      real*4 gz(nj,ni),gpeek(nj,ni),gtopo(nj,ni)
      real*4 u(nj,ni),v(nj,ni),w(nj,ni)
      real*4 xwind(nj),ywind(ni)
      complex*8 fkz(nj2,ni),gkz(nj2,ni)
      complex*8 uk(nj2,ni),vk(nj2,ni),wk(nj2,ni)
      complex*8 fzval,gzval
      dimension n(2)
      complex*8 work(nwork)
      equivalence (fz(1,1),fkz(1,1))
      equivalence (gz(1,1),gkz(1,1))
      equivalence (u(1,1),uk(1,1))
      equivalence (v(1,1),vk(1,1))
      equivalence (w(1,1),wk(1,1))
c
      pi=acos(-1.)
c
c  zero the arrays fx,fy,fz
c
      do 30 i=1,ni
      do 30 j=1,nj
      fz(j,i)=0.
      gz(j,i)=0.
  30  continue
c
c  set the dimensions for fourt
c
      n(1)=nj
      n(2)=ni
c
c   get values from command line
c
      narg = iargc()
      if(narg.lt.8) then
        write(*,'(a)')'  '
        write(*,'(a)')
     &  'Usage: predictgrav input.grd zs aobs H Te rhoc itg grav.grd'
        write(*,'(a)')
        write(*,'(a)')
     &  '       input.grd - either high passed surface topo file [m]'
        write(*,'(a)')
     &  '                   OR high passed gravity data file [mGal]'
        write(*,'(a)')
     &  '       zs        - mean elevation of original topo [m]'
        write(*,'(a)')
     &  '       aobs      - altitude of gravity observation [m]'
        write(*,'(a)')
     &  '       H         - depth to mean Moho [km]'
        write(*,'(a)')
     &  '       Te        - elastic thickness [km]'
        write(*,'(a)')
     &  '       rhoc      - crustal density [kg/m^3]'
        write(*,'(a)')
     &  '       itg       - input is topography or gravity'
        write(*,'(a)')
     &  '                   (0)-predict gravity from input topography'
        write(*,'(a)')
     &  '                   (1)-upward continue input gravity'
        write(*,'(a)')
     &  '       grav.grd  - output gravity [mGal]'
        write(*,'(a)')'  '
        stop
      else 
        call getarg(1,ctopo)
        nc=index(ctopo,' ')
        ctopo(nc:nc)=char(0)
        call getarg(2,czs)
        call getarg(3,caobs)
        call getarg(4,ch)
        call getarg(5,cTe)
        call getarg(6,crhoc)
        call getarg(7,citg)
        read(czs,*)zs
        read(caobs,*)aobs
        read(ch,*)thk
        read(cTe,*)Te
        read(crhoc,*)rhoc
        read(citg,*)itg
        call getarg(8,fxdisp)
        nc=index(fxdisp,' ')
        fxdisp(nc:nc)=char(0)
      endif
c
c   set Young's modulus to 70 Gpa and Poisson's ratio to .25
c
      grv=9.81
      young=7.e10
      rnu=0.499
      
      rlam1=young*rnu/(1+rnu)/(1-2*rnu)
      rmu1=young/2/(1+rnu)
      bulk=young/3/(1-2*rnu)
      
c      rlam1=2*young/5.
c      rmu1=rlam1
c      bulk=rlam1+2*rmu1/3.
      
      H=-abs(thk*1000.)
      Te=abs(Te*1000.)
      
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
     +            dlt,dln,rdum,title,ctopo)
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
c  apply window to input topo/grav
c
      do 120 i=1,ni
      do 120 j=1,nj
       fz(j,i)=topo(j,i)*xwind(j)*ywind(i)/(ni*nj)
       fpeek(j,i)=fz(j,i)
  120 continue
c
c  take the fourier transform of the input
c
      call fourt(fz,n,2,-1,0,work,nwork)
c
c  compute output gravity
c
      do 255 i=1,ni
      ky=-(i-1)/height
      if(i.ge.ni2) ky= (ni-i+1)/height
      do 255 j=1,nj2
      kx=(j-1)/width

      fzval=fkz(j,i)

      if(itg.eq.1) then
c
c upward continue input grav
c
        kh2=kx*kx+ky*ky
        kh=sqrt(kh2)
        beta=2*pi*kh
        uk(j,i)=fzval*exp(-beta*aobs)
      else
c
c caluclate grav from topo
c
        call gravcalc(kx,ky,H,Te,rhoc,rhom,aobs,fzval,zs,grval)
        uk(j,i)=grval
      endif

 255  continue
c
c  do the inverse fft's
c
      call fourt(u,n,2,1,-1,work,nwork)
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
      call writegrd(u,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              fxdisp,fxdisp)
      stop
      end

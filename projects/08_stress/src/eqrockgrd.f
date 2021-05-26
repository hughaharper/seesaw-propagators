c
      program eqrockgrd
c
c*****************   main program   ************************
c
c  Make equivalent rock topography for a specified rhoc value.
c  Wet/Dry areas are determined by srtm.bitmask file, pointed
c  to in readsrtmmask.c
c  In the oceans,  eqrk_topo=(rhoc-rhow)/rhoc*topo

c
c***********************************************************
c
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
      real*8 rln0,rlt0,dlt,dln,ddx,ddy,rland,rdum
c
c  change ni and nj as needed
c
      parameter(ni=2000,nj=3000)
c      parameter(ni=900,nj=1800)
      character*80 ofilnam1,ifilnam,title,crhoc
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi
c
      real*4 ftopo(nj,ni)
      integer*4 ilandorsea
c
      pi=acos(-1.)
c
c  zero the arrays fx,fy,fz
c
      do 30 i=1,ni
      do 30 j=1,nj
c      fz(j,i)=0.
c      gz(j,i)=0.
  30  continue
c
c   get values from command line
c
      narg = iargc()
      if(narg.lt.3) then
        write(*,'(a)')'  '
        write(*,'(a)')
     &  'Usage: eqrockgrd topo.grd topo.eqrk.grd rhoc'
        write(*,'(a)')
        write(*,'(a)')
     &  '       topo.grd - grd file'
        write(*,'(a)')
     &  '       topo.eqrk.grd - output equivalent rock topography file'
        write(*,'(a)')
     &  '       rhoc     - crustal density [kg/m^3]'
        write(*,'(a)')
        stop
      else 
        call getarg(1,ifilnam)
c        nc=index(ifilnam,' ')
c        ifilnam(nc:nc)='\0'
        call getarg(2,ofilnam1)
c        nc=index(ofilnam1,' ')
c        ofilnam1(nc:nc)='\0'
        call getarg(3,crhoc)
        read(crhoc,*)rhoc
      endif
c
c   read the grd files
c
      call readgrd(ftopo,nj1,ni1,rlt0,rln0,
     +            dlt,dln,rdum,title,trim(ifilnam)//char(0))
      if(ni1.ne.ni.or.nj1.ne.nj) then
        print *,ni1,ni,nj1,nj,rlt0,rln0,dlt,dln,rdum
        write(*,'(a)')' recompile program to match topo size'
        stop
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c set parameters 
c
c      rhoc=2800.
      rhom=3300.
      rhow=1025.
c
c main loop
c
      do 120 i=1,ni
      do 120 j=1,nj
c
c what is the lat and lon of each point
c
	rlt=(i-1)*dlt+rlt0
        rln=(j-1)*dln+rln0
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

  120 continue
c
c  write 3 grd files some parameters must be real*8
c
      ddx=dlt
      ddy=dln
      rland=9998.
      rdum=9999.
      if(istr.ge.1) then
       rland=9998.d0*rmu1
       rdum=9999.d0*rmu1
      endif
      call writegrd(ftopo,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ofilnam1)//char(0),trim(ofilnam1)//char(0))
  190 continue
  998 continue
      stop
      end

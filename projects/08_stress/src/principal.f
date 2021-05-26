c
      program principal
c
c*****************   main program   ************************
c
c  Calculates the principal stress components for a given
c  stress tensor.  Outputs magnitude and orientation information
c  in one of three formats: eigen vectors, optimal slip plane
c  orientations, or orientation of principal axes.
c
c***********************************************************
c
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
      real*8 rln0,rlt0,dlt,dln,ddx,ddy,rland,rdum
c
c  change ni and nj as needed
c
C       parameter(ni=334,nj=370,N=3)
c      parameter(ni=810,nj=1080,N=3)
c      parameter(ni=1024,nj=1024,N=3)
C       parameter(ni=4096,nj=4096,N=3)
      parameter(ni=1800,nj=4800,N=3)
      character*100 cs1,cs2,cs3,ct1,ct2,ct3,cp1,cp2,cp3
      character*100 cTxx,cTyy,cTzz,cTxy,cTxz,cTyz,title,cout
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi
c
      real*4 Txx(nj,ni),Tyy(nj,ni),Tzz(nj,ni)
      real*4 Txy(nj,ni),Txz(nj,ni),Tyz(nj,ni)
      real*4 sig1(nj,ni),sig2(nj,ni),sig3(nj,ni)
      real*4 theta1(nj,ni),theta2(nj,ni),theta3(nj,ni)
      real*4 phi1(nj,ni),phi2(nj,ni),phi3(nj,ni)
      real*4 H(N,N),eigen(N,N)
      real*8 v1(3),v2(3),v3(3)
      integer A,B,C
c
      pi=acos(-1.)
c
c   get values from command line
c
      narg = iargc()
      if(narg.lt.16) then
        write(*,'(a)')'  '
        write(*,'(a)')
     &  'Usage: principal Txx Tyy Tzz Txy Txz Tyz iout outfiles'
        write(*,'(a)')
        write(*,'(a)')
     &  '       Txx.grd - input cartesian stress (6)'
        write(*,'(a)')
     &  '       iout    - select output option'
        write(*,'(a)')
     &  '             0 = principal stress and eigenvector components'
        write(*,'(a)')
     &  '                 v1xyz & v3xyz (v2 is orth to these)(9 files)'
        write(*,'(a)')
     &  '             1 = principal stress and strike, dip, rake of two'
        write(*,'(a)')
     &  '                 planes of maximum shear stress in'
        write(*,'(a)')
     &  '                 Aki & Richards format (9 files)'
        write(*,'(a)')
     &  '             2 = principal stress and axis orientation'
        write(*,'(a)')
     &  '                 in spherical azimuth and plunge (9 files)'
        write(*,'(a)')
        write(*,'(a)')
     &  '       outfiles - sig1.grd sig2.grd sig3.grd (3)'
        write(*,'(a)')
     &  '             AND  v1x.grd v1y.grd v1z.grd (3)'
        write(*,'(a)')
     &  '                  v3x.grd v3y.grd v3z.grd (3)'
        write(*,'(a)')
     &  '              OR  strike.grd - deg CW from north [0,360] (2)'
        write(*,'(a)')
     &  '                  dip.grd - deg below horiz. RHS [0,90] (2)'
        write(*,'(a)')
     &  '                  rake.grd - deg up from strike [-180,180] (2)'
        write(*,'(a)')
     &  '              OR  azimuth.grd - deg CW from north [0,360] (3)'
        write(*,'(a)')
     &  '                  plunge.grd - deg below horizontal [0,90] (3)'
        write(*,'(a)')'  '
       stop
      else 
        call getarg(1,cTxx)
c         nc=index(cTxx,' ')
c         cTxx(nc:nc)='\0'
        call getarg(2,cTyy)
c         nc=index(cTyy,' ')
c         cTyy(nc:nc)='\0'
        call getarg(3,cTzz)
c         nc=index(cTzz,' ')
c         cTzz(nc:nc)='\0'
        call getarg(4,cTxy)
c         nc=index(cTxy,' ')
c         cTxy(nc:nc)='\0'
        call getarg(5,cTxz)
c         nc=index(cTxz,' ')
c         cTxz(nc:nc)='\0'
        call getarg(6,cTyz)
c         nc=index(cTyz,' ')
c         cTyz(nc:nc)='\0'
        call getarg(7,cout)
        read(cout,*)iout
        call getarg(8,cs1)
c         nc=index(cs1,' ')
c         cs1(nc:nc)='\0'
        call getarg(9,cs2)
c         nc=index(cs2,' ')
c         cs2(nc:nc)='\0'
        call getarg(10,cs3)
c         nc=index(cs3,' ')
c         cs3(nc:nc)='\0'
        call getarg(11,ct1)
c         nc=index(ct1,' ')
c         ct1(nc:nc)='\0'
        call getarg(12,ct2)
c         nc=index(ct2,' ')
c         ct2(nc:nc)='\0'
        call getarg(13,ct3)
c         nc=index(ct3,' ')
c         ct3(nc:nc)='\0'
        call getarg(14,cp1)
c         nc=index(cp1,' ')
c         cp1(nc:nc)='\0'
        call getarg(15,cp2)
c         nc=index(cp2,' ')
c         cp2(nc:nc)='\0'
        call getarg(16,cp3)
c         nc=index(cp3,' ')
c         cp3(nc:nc)='\0'
      endif
c
c   read the grd files
c
      call readgrd(Txx,nj1,ni1,rlt0,rln0,
     +            dlt,dln,rdum,title,trim(cTxx)//char(0))
      if(ni1.ne.ni.or.nj1.ne.nj) then
        print *,ni1,ni,nj1,nj,rlt0,rln0,dlt,dln,rdum
        write(*,'(a)')' recompile program to match topo size'
        stop
      endif

      call readgrd(Tyy,nj1,ni1,rlt0,rln0,
     +            dlt,dln,rdum,title,trim(cTyy)//char(0))
      call readgrd(Tzz,nj1,ni1,rlt0,rln0,
     +            dlt,dln,rdum,title,trim(cTzz)//char(0))
      call readgrd(Txy,nj1,ni1,rlt0,rln0,
     +            dlt,dln,rdum,title,trim(cTxy)//char(0))
      call readgrd(Txz,nj1,ni1,rlt0,rln0,
     +            dlt,dln,rdum,title,trim(cTxz)//char(0))
      call readgrd(Tyz,nj1,ni1,rlt0,rln0,
     +            dlt,dln,rdum,title,trim(cTyz)//char(0))
c
c  loop over values, fill H array with Tij, call jacobi
c
      do 120 i=1,ni
      do 120 j=1,nj
c       P=(Txx(j,i)+Tyy(j,i)+Tzz(j,i))/3
       P=0.0
       
       H(1,1)=Txx(j,i)-P
       H(2,2)=Tyy(j,i)-P
       H(3,3)=Tzz(j,i)-P
       H(1,2)=Txy(j,i)
       H(1,3)=Txz(j,i)
       H(2,3)=Tyz(j,i)
       H(2,1)=H(1,2)
       H(3,1)=H(1,3)
       H(3,2)=H(2,3)
       
c       print *,P
c       print *,H
       
       call jacobi(H,eigen,N,N,100)
c       print *,H
c       print *,eigen
c
c largest eigen value (sigma1) is in the A,A position
c smallest (sigma3) is in the C,C position
c

       if (H(1,1).ge.H(2,2).and.H(1,1).ge.H(3,3)) then
         A=1
	 if (H(2,2).gt.H(3,3)) then
           B=2
           C=3
	 else
	   C=2
	   B=3
	 endif
       elseif (H(2,2).ge.H(1,1).and.H(2,2).ge.H(3,3)) then
         A=2
	 if (H(1,1).gt.H(3,3)) then
           B=1
           C=3
	 else
	   B=3
	   C=1
	 endif
       else
         A=3
	 if (H(1,1).gt.H(2,2)) then
           B=1
           C=2
	 else
	   B=2
	   C=1
	 endif
       endif
c
c fill the ordered arrays
c  
       sig1(j,i)=H(A,A)
       sig2(j,i)=H(B,B)
       sig3(j,i)=H(C,C)
c
c if iout = 0 "theta" becomes v1 components, "phi" becomes v3 components
c
       if(iout.eq.0) then
         theta1(j,i)=eigen(1,A)
	 theta2(j,i)=eigen(2,A)
	 theta3(j,i)=eigen(3,A)
	 phi1(j,i)=eigen(1,C)
	 phi2(j,i)=eigen(2,C)
	 phi3(j,i)=eigen(3,C)
       endif
c
c if iout = 1, compute strike and dip and rake of maximum shear stress planes
c (this is in Aki & Richards format)
c
       if (iout.eq.1) then
c         v1 = (/ eigen(1,A), eigen(2,A), eigen(3,A) /)
c         v2 = (/ eigen(1,B), eigen(2,B), eigen(3,B) /)
c         v3 = (/ eigen(1,C), eigen(2,C), eigen(3,C) /)
         
	 v1(1) = eigen(1,A) 
	 v1(2) = eigen(2,A) 
	 v1(3) = eigen(3,A) 
         
	 v2(1) = eigen(1,B) 
	 v2(2) = eigen(2,B) 
	 v2(3) = eigen(3,B) 

	 v3(1) = eigen(1,C) 
	 v3(2) = eigen(2,C) 
	 v3(3) = eigen(3,C) 
      
         call compsdr(v1,v2,v3,strike1,strike2,dip1,dip2,rake1,rake2)
       
         theta1(j,i)=strike1
	 theta2(j,i)=strike2
	 theta3(j,i)=dip1
	 phi1(j,i)=dip2
	 phi2(j,i)=rake1
	 phi3(j,i)=rake2
       endif
c
c if iout = 2 compute the orientations of the principal axes
c
       if(iout.eq.2) then
c         v1 = (/ eigen(1,A), eigen(2,A), eigen(3,A) /)
c         v2 = (/ eigen(1,B), eigen(2,B), eigen(3,B) /)
c         v3 = (/ eigen(1,C), eigen(2,C), eigen(3,C) /)
         
	 v1(1) = eigen(1,A) 
	 v1(2) = eigen(2,A) 
	 v1(3) = eigen(3,A) 
         
	 v2(1) = eigen(1,B) 
	 v2(2) = eigen(2,B) 
	 v2(3) = eigen(3,B) 

	 v3(1) = eigen(1,C) 
	 v3(2) = eigen(2,C) 
	 v3(3) = eigen(3,C) 
        
	 call compTBPazpl(v1,v2,v3,azimuth1,azimuth2,
     +                       azimuth3,plunge1,plunge2,plunge3)         
	 
         theta1(j,i)=azimuth1
         theta2(j,i)=azimuth2
         theta3(j,i)=azimuth3
         phi1(j,i)=plunge1
         phi2(j,i)=plunge2
         phi3(j,i)=plunge3
       endif
c
c for debugging
c
c      write(*,*)Txx(j,i),Tyy(j,i),Tzz(j,i),Txy(j,i),Txz(j,i)
c     &,Tyz(j,i),sig1(j,i),sig2(j,i),sig3(j,i)
c      goto 998

  120 continue
c
c  write grd files some parameters must be real*8
c
      ddx=dlt
      ddy=dln
      rland=9998.
      rdum=9999.
      if(istr.ge.1) then
       rland=9998.d0*rmu1
       rdum=9999.d0*rmu1
      endif
      call writegrd(sig1,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(cs1)//char(0),trim(cs1)//char(0))
      call writegrd(sig2,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(cs2)//char(0),trim(cs2)//char(0))
      call writegrd(sig3,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(cs3)//char(0),trim(cs3)//char(0))
      call writegrd(theta1,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ct1)//char(0),trim(ct1)//char(0))
      call writegrd(theta2,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ct2)//char(0),trim(ct2)//char(0))
      call writegrd(theta3,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(ct3)//char(0),trim(ct3)//char(0))
      call writegrd(phi1,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(cp1)//char(0),trim(cp1)//char(0))
      call writegrd(phi2,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(cp2)//char(0),trim(cp2)//char(0))
      call writegrd(phi3,nj,ni,rlt0,rln0,ddx,ddy,rland,rdum,
     +              trim(cp3)//char(0),trim(cp3)//char(0))
  190 continue
  998 continue
      stop
      end

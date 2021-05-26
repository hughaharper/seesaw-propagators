c
      program Tij2principal
c
c*****************   main program   ************************
c
c  Reads six column ascii file containing cartesion tensor
c  components (usually of stress, but should be general).
c  Computes principal stresses and axes (eigen values and
c  vectors) and outputs them in any of a number of formats, 
c  depending on what is requested.  Output is a single ascii
c  file with columns defined in first line header and rows
c  with values following.  Non-defined values are represented
c  as a NaN.
c
c***********************************************************
c
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
      real*8 rln0,rlt0,dlt,dln,ddx,ddy,rland,rdum
c
c  change ni and nj as needed
c
      parameter(ni=1024,nj=1024,N=3,N2=2)
      character*80 filein,fileout,cdev,cxy
      character*80 cout1,cout2,cout3,cout4,cout5,cout6,cout7
      character*80 cH0,cH1,cH2,cH3,cH4,cH5,cH6,cH7,cHout
      data lui,luo/19,20/
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi
c
      real*4 H(N,N),eigen(N,N)
      real*4 H2(N2,N2),eigen2(N2,N2)
      real*8 v1(3),v2(3),v3(3)
      integer A,B,C
c
      pi=acos(-1.)
      fNaN=-1.0
      fNaN=sqrt(fNaN)
c
c   get values from command line
c
      narg = iargc()
      if(narg.lt.11) then
        write(*,'(a)')'  '
        write(*,'(a)')
     &  'Usage: Tij2principal infile.dat ixy iDev iout1 iout2 iout3 iout
     &4 iout5 iout6 iout7 outfile.dat'
        write(*,'(a)')
        write(*,'(a)')
     &  '       infile.dat - ascii file whose 6 columns are cartesian st
     &ress components: Txx Tyy Tzz Txy Txz Tyz'
        write(*,'(a)')
        write(*,'(a)')
     &  '       ixy        - infile contains xy columns [1=y/0=n]: if 1,
     & first two columns are ignored'
        write(*,'(a)')
        write(*,'(a)')
     &  '       iDev       - evaluate deviatoric tensor? [1=y/0=n]: if 1
     &, devTij = Tij - I*Tkk/3'
        write(*,'(a)')
        write(*,'(a)')
     &  '       iout#      - 7 flags [1/0] determining which values are 
     &output'
        write(*,'(a)')
     &  '              (1) - principal stresses (3 values) sig1 > sig2 >
     & sig3'
        write(*,'(a)')
     &  '                    "physics" sign convention (tension is posit
     &ive, compression is negative)'
        write(*,'(a)')
     &  '              (2) - iso dif rat (3 values)'
        write(*,'(a)')
     &  '                    isotropic stress: (sig1+sig2+sig3)/3'
        write(*,'(a)')
     &  '                    differential stress: sig1-sig3'
        write(*,'(a)')
     &  '                    stress ratio: (sig2-sig3)/(sig1-sig3) [0,1]
     & (also known as "R" or "phi")'
        write(*,'(a)')
     &  '              (3) - stress invariants (3 values) I_T II_T III_T
     &'
        write(*,'(a)')
     &  '                    roots of dev(T-kI)=0, such that k^3 - I_T*k
     &^2 - II_T*k - III_T = 0'
        write(*,'(a)')
     &  '                    I_T = Tkk = sig1 + sig2 + sig3'
        write(*,'(a)')
     &  '                    II_T = 0.5(T:T-I_T^2) = -(sig1*sig2 + sig2*
     &sig3 + sig1*sig3)'
        write(*,'(a)')
     &  '                    III_T = det(T) = sig1*sig2*sig3'
        write(*,'(a)')
     &  '              (4) - unit eigenvector components (9 values)'
        write(*,'(a)')
     &  '                    v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z
     &'
        write(*,'(a)')
     &  '              (5) - SHmax Aphi (2 values)'
        write(*,'(a)')
     &  '                    SHmax: azimuth of maximum horizontal compre
     &ssion'
        write(*,'(a)')
     &  '                    (physics sign convention) in degrees EofN [
     &-90,90]'
        write(*,'(a)')
     &  '                    Aphi: regime-modified stress ratio [0,3]'
        write(*,'(a)')
     &  '                          phi = (sig1-sig2)/(sig1-sig3)'
        write(*,'(a)')
     &  '                           { phi   if |v3z| is maximum vertical
     &'
        write(*,'(a)')
     &  '                    Aphi = { 2-phi if |v2z| is maximum vertical
     &'
        write(*,'(a)')
     &  '                           { 2+phi if |v1z| is maximum vertical
     &'
        write(*,'(a)')
     &  '                    note: for CSM convention, phi ~= stress rat
     &io'
        write(*,'(a)')
     &  '              (6) - strike dip rake (6 values) of two planes bi
     &secting v1 and v3'
        write(*,'(a)')
     &  '                    strike1 strike2 dip1 dip2 rake1 rake2 in Ak
     &i & Richards format'
        write(*,'(a)')
     &  '                    (strike - deg CW from north [0,360])'
        write(*,'(a)')
     &  '                    (dip- deg below horiz. RHS [0,90])'
        write(*,'(a)')
     &  '                    (rake - deg CCW (up) from strike [-180,180]
     &)'
        write(*,'(a)')
     &  '              (7) - azimuth plunge (6 values) of three principa
     &l stress axes'
        write(*,'(a)')
     &  '                    azim1 azim2 azim3 plunge1 plunge2 plunge3'
        write(*,'(a)')
     &  '                    (azimuth - deg CW from north [0,360])'
        write(*,'(a)')
     &  '                    (plunge - deg below horizontal [0,90])'
        write(*,'(a)')
        write(*,'(a)')
     &  '       outfile.dat - ascii file with nrows_out = nrowsin + 1'
        write(*,'(a)')
     &  '                    (first row is header line)'
        write(*,'(a)')
     &  '                    ncolumns depends on iout flags [0,32]'
        write(*,'(a)')'  '
       stop
      else 
        call getarg(1,filein)
        call getarg(2,cxy)
        call getarg(3,cdev)
        call getarg(4,cout1)
        call getarg(5,cout2)
        call getarg(6,cout3)
        call getarg(7,cout4)
        call getarg(8,cout5)
        call getarg(9,cout6)
        call getarg(10,cout7)
        call getarg(11,fileout)
        read(cxy,*)ixy
        read(cdev,*)idev
        read(cout1,*)iout1
        read(cout2,*)iout2
        read(cout3,*)iout3
        read(cout4,*)iout4
        read(cout5,*)iout5
        read(cout6,*)iout6
        read(cout7,*)iout7
      endif
c
c open the files
c
      open(unit=lui,file=filein,err=9000,status='old')
      open(unit=luo,file=fileout,err=9002,status='unknown')
c
c write the header line
c
      cH0 = '# '
      cxy = ' x y '
      cH1 = ' sig1 sig2 sig3 '
      cH2 = ' iso dif rat '
      cH3 = ' I_T II_T III_T '
      cH4 = ' v1x v1y v1z v2x v2y v2z v3x v3y v3z '
      cH5 = ' SHmax Aphi '
      cH6 = ' strike1 strike2 dip1 dip2 rake1 rake2 '
      cH7 = ' azim1 azim2 azim3 plunge1 plunge2 plunge3'
      
      
      write(luo,'(a,$)')trim(cH0)
      if (ixy.eq.1) write(luo,'(a,$)')trim(cxy)
      if (iout1.eq.1) write(luo,'(a,$)')trim(cH1)
      if (iout2.eq.1) write(luo,'(a,$)')trim(cH2)
      if (iout3.eq.1) write(luo,'(a,$)')trim(cH3)
      if (iout4.eq.1) write(luo,'(a,$)')trim(cH4)
      if (iout5.eq.1) write(luo,'(a,$)')trim(cH5)
      if (iout6.eq.1) write(luo,'(a,$)')trim(cH6)
      if (iout7.eq.1) write(luo,'(a,$)')trim(cH7)
      write(luo,'(a)')' '

c
c  loop over values, fill H array with Tij, call jacobi
c
  100 continue
      if (ixy.eq.0) read(lui,*,end=200)Txx,Tyy,Tzz,Txy,Txz,Tyz
      if (ixy.eq.1) then
        read(lui,*,end=200)rx,ry,Txx,Tyy,Tzz,Txy,Txz,Tyz
        write(luo,899)rx,ry
  899   format(2F11.6,$)     
      endif

c
c  make stress deviatoric, if requested
c
      if (idev.eq.1) P=(Txx+Tyy+Tzz)/3.
      if (idev.eq.0) P=0.0
c
c  fill the H array with Tij
c   
       H(1,1)=Txx-P
       H(2,2)=Tyy-P
       H(3,3)=Tzz-P
       H(1,2)=Txy
       H(1,3)=Txz
       H(2,3)=Tyz
       H(2,1)=H(1,2)
       H(3,1)=H(1,3)
       H(3,2)=H(2,3)
       
c       print *,P
c       print *,H

       
       call jacobi(H,eigen,N,N,100)
c
c special case: NaNs are present
c
       if(isnan(Txx).or.isnan(Tyy).or.isnan(Tzz).or.
     &    isnan(Txy).or.isnan(Txz).or.isnan(Tyz)) then
c         print *,'something is a nan'
	 eigen(1,1)=fNaN
	 eigen(2,2)=fNaN
	 eigen(3,3)=fNaN
	 eigen(1,2)=fNaN
	 eigen(1,3)=fNaN
	 eigen(2,3)=fNaN
	 eigen(2,1)=fNaN
	 eigen(3,1)=fNaN
	 eigen(3,2)=fNaN
       endif
c
c special case: all values are zero (sig123 is defined, but v123 isn't)
c
       if(Txx.eq.0.and.Tyy.eq.0.and.Tzz.eq.0.and.
     &    Txy.eq.0.and.Txz.eq.0.and.Tyz.eq.0) then
c         print *,'something is a nan'
	 eigen(1,1)=fNaN
	 eigen(2,2)=fNaN
	 eigen(3,3)=fNaN
	 eigen(1,2)=fNaN
	 eigen(1,3)=fNaN
	 eigen(2,3)=fNaN
	 eigen(2,1)=fNaN
	 eigen(3,1)=fNaN
	 eigen(3,2)=fNaN
       endif
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
       sig1=H(A,A)
       sig2=H(B,B)
       sig3=H(C,C)
       
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

c       print *,sig1,sig2,sig3
c       print *,v1
c       print *,v2
c       print *,v3

c
c (1) principal stresses
c
       if(iout1.eq.1) then
         write(luo,901)sig1,sig2,sig3
  901    format(3d15.6,$)     
       endif
c
c (2) iso, dif, rat
c
       if(iout2.eq.1) then
         fiso=(sig1+sig2+sig3)/3.0
	 fdif=sig1-sig3
	 frat=(sig2-sig3)/(sig1-sig3)
	 write(luo,902)fiso,fdif,frat
  902    format(3d15.6,$)     
       endif
c
c (3) invariants
c
       if(iout3.eq.1) then
         finv1=sig1+sig2+sig3
	 finv2=-(sig1*sig2+sig2*sig3+sig1*sig3)
	 finv3=sig1*sig2*sig3
	 write(luo,903)finv1,finv2,finv3
  903    format(3d15.6,$)     
       endif
c
c (4) unit eigenvector components
c
       if(iout4.eq.1) then
	 write(luo,904)v1,v2,v3
  904    format(9F10.5,$)     
       endif
c
c (5) SHmax and Aphi
c
       if(iout5.eq.1) then
	 v1z=abs(v1(3))
	 v2z=abs(v2(3))
	 v3z=abs(v3(3))

c	 frat=(sig2-sig3)/(sig1-sig3)
c	 if(v1z.gt.v2z.and.v1z.gt.v3z) Aphi=frat
c	 if(v2z.gt.v1z.and.v2z.gt.v3z) Aphi=2.0-frat
c	 if(v3z.gt.v1z.and.v3z.gt.v2z) Aphi=2.0+frat

	 frat=(sig1-sig2)/(sig1-sig3)
	 Aphi=fNaN
	 if(v1z.gt.v2z.and.v1z.gt.v3z) Aphi=2.0+frat
	 if(v2z.gt.v1z.and.v2z.gt.v3z) Aphi=2.0-frat
	 if(v3z.gt.v1z.and.v3z.gt.v2z) Aphi=frat

	 if (sig1.eq.sig3) Aphi=fNaN
	 
c	 print *,frat,Aphi

	 H2(1,1)=Txx-P
	 H2(2,2)=Tyy-P
	 H2(1,2)=Txy
	 H2(2,1)=H2(1,2)

	 call jacobi(H2,eigen2,N2,N2,100)

       if(isnan(Txx).or.isnan(Tyy).or.isnan(Txy)) then
c         print *,'something is a nan'
	 eigen2(1,1)=fNaN
	 eigen2(2,2)=fNaN
	 eigen2(1,2)=fNaN
	 eigen2(2,1)=fNaN
       endif

	 if (H2(1,1).gt.H2(2,2)) then
           A=1
           B=2
	 elseif (H2(2,2).gt.H2(1,1)) then
           A=2
	   B=1
	 endif
	 
	 SHmax=180.0/pi*atan2(eigen2(1,B),eigen2(2,B))
	 if (H2(1,1).eq.H2(2,2)) SHmax=fNaN
	 if (SHmax.lt.-90) SHmax=SHmax+180.
	 if (SHmax.gt.90) SHmax=SHmax-180.
	 
	 write(luo,905)SHmax,Aphi
  905    format(2F12.5,$)     
       endif
c
c (6) Strike, Dip, and Rake of 2 planes
c
       if(iout6.eq.1) then
	 call compsdr(v1,v2,v3,strike1,strike2,dip1,dip2,rake1,rake2)
	 
	 if (sig1.eq.sig3.or.sig2.eq.sig1.or.sig2.eq.sig3) then
	   strike1=fNaN
	   strike2=fNaN
	   dip1=fNaN
	   dip2=fNaN
	   rake1=fNaN
	   rake2=fNaN
	 endif
	 
	 write(luo,906)strike1,strike2,dip1,dip2,rake1,rake2
  906    format(6F12.5,$)     
       endif
c
c (7) Azimuth and Plunge of principal axes
c
       if(iout7.eq.1) then
	 call compTBPazpl(v1,v2,v3,azimuth1,azimuth2,
     +                       azimuth3,plunge1,plunge2,plunge3)         
	 write(luo,907)azimuth1,azimuth2,azimuth3,plunge1,plunge2,plunge3
  907    format(6F12.5,$)     
       endif
c
c  End of Loop
c
      write(luo,'(a)')' '
      go to 100
  200 continue
      stop
c
c  could not open file
c
 9000 continue
      write(*,*)' could not open ',filein
 9002 continue
      write(*,*)' could not open ',fileout
      end                

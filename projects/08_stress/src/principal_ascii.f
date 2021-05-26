c
      program principal_ascii
c
c*****************   main program   ************************
c
c  Same as principal, only input and output are ascii columns
c  in a single file instead of separate grdfiles.
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
      parameter(ni=1024,nj=1024,N=3)
      character*80 filein,fileout,cout
      data lui,luo/19,20/
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi
c
      real*4 H(N,N),eigen(N,N)
      real*8 v1(3),v2(3),v3(3)
      integer A,B,C
c
      pi=acos(-1.)
c
c   get values from command line
c
      narg = iargc()
      if(narg.lt.3) then
        write(*,'(a)')'  '
        write(*,'(a)')
     &  'Usage: principal_ascii Tij.dat iout outfile.dat'
        write(*,'(a)')
        write(*,'(a)')
     &  '       Tij.dat - input file whose 6 columns are cartesian'
        write(*,'(a)')
     &  '                 stress: Txx Tyy Tzz Txy Txz Tyz'
        write(*,'(a)')
     &  '       iout    - select output option'
        write(*,'(a)')
     &  '             0 = principal stress and eigenvector components'
        write(*,'(a)')
     &  '             1 = principal stress and strike, dip, rake of two'
        write(*,'(a)')
     &  '                 planes of maximum shear stress in'
        write(*,'(a)')
     &  '                 Aki & Richards format'
        write(*,'(a)')
     &  '             2 = principal stress and axis orientation'
        write(*,'(a)')
     &  '                 in spherical azimuth and plunge'
        write(*,'(a)')
        write(*,'(a)')
     &  '       outfile.dat - output file whose 9 columns are'
        write(*,'(a)')
     &  '                  sig1 sig2 sig3'
        write(*,'(a)')
     &  '             AND  v1x v1y v1z v3x v3y v3z'
        write(*,'(a)')
     &  '              OR  strike1 strike2 dip1 dip2 rake1 rake2'
        write(*,'(a)')
     &  '                  (strike - deg CW from north [0,360])'
        write(*,'(a)')
     &  '                  (dip- deg below horiz. RHS [0,90])'
        write(*,'(a)')
     &  '                  (rake - deg up from strike [-180,180])'
        write(*,'(a)')
     &  '              OR  azim1 azim2 azim3 plunge1 plunge2 plunge3'
        write(*,'(a)')
     &  '                  (azimuth - deg CW from north [0,360])'
        write(*,'(a)')
     &  '                  (plunge - deg below horizontal [0,90])'
        write(*,'(a)')'  '
       stop
      else 
        call getarg(1,filein)
        call getarg(2,cout)
        read(cout,*)iout
        call getarg(3,fileout)
      endif
c
c open the files
c
      open(unit=lui,file=filein,err=9000,status='old')
      open(unit=luo,file=fileout,err=9002,status='unknown')
c
c  loop over values, fill H array with Tij, call jacobi
c
  100 read(lui,*,end=200)Txx,Tyy,Tzz,Txy,Txz,Tyz

c       P=(Txx+Tyy+Tzz)/3
       P=0.0
       
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
c
c if iout=2 "theta" becomes v1 components, "phi" becomes v3 components
c
       if(iout.eq.0) then
         theta1=eigen(1,A)
	 theta2=eigen(2,A)
	 theta3=eigen(3,A)
	 phi1=eigen(1,C)
	 phi2=eigen(2,C)
	 phi3=eigen(3,C)
       endif
c
c if iout = 3, compute strike and dip and rake of maximum shear stress planes
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
       
         theta1=strike1
	 theta2=strike2
	 theta3=dip1
	 phi1=dip2
	 phi2=rake1
	 phi3=rake2
       endif
c
c if iout=0 compute the orientations of the principal axes
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
	 
         theta1=azimuth1
         theta2=azimuth2
         theta3=azimuth3
         phi1=plunge1
         phi2=plunge2
         phi3=plunge3
       endif
c
c  write outfile
c
       write(luo,901)sig1,sig2,sig3,theta1,theta2,theta3,phi1,phi2,phi3
  901  format(3d15.6,6F10.3)
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

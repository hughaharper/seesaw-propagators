c
      subroutine princstr(Txx,Tyy,Tzz,Txy,Txz,Tyz,
     +                    sig1,sig2,sig3,v1x,v1y,v1z,v3x,v3y,v3z)
c
c Computes principal stresses from a 3-D symmetric 
c cartesian stress tensor.
c
c  input 
c    Txx - x normal stress [Pa]
c    Tyy - y normal stress [Pa]
c    Tzz - z normal stress [Pa]
c    Txy - xy shear stress [Pa]
c    Txz - xz shear stress [Pa]
c    Tyz - yz shear stress [Pa]
c  output
c    sig1 - extentional principal stress  (+)   [Pa]
c    sig2 - intermediate principal stress (+/-) [Pa]
c    sig3 - compressional principal stress  (-) [Pa]
c
      parameter(N=3)
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      implicit real*8 (k)
      real*4 H(N,N),eigen(N,N)
      real*8 v1(3),v2(3),v3(3)
      integer A,B,C
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
c  subtract the pressure and shape the cartesian components into a single
c  3x3 matrix
c
      P=(Txx+Tyy+Tzz)/3
       
      H(1,1)=Txx-P
      H(2,2)=Tyy-P
      H(3,3)=Tzz-P
      H(1,2)=Txy
      H(1,3)=Txz
      H(2,3)=Tyz
      H(2,1)=H(1,2)
      H(3,1)=H(1,3)
      H(3,2)=H(2,3)

c      print *,P
c      print *,H
c
c get the eigenvalues and eigenvectors
c
      call jacobi(H,eigen,N,N,100)
c
c determine the order of the eigenvalues:
c largest eigen value (sigma1) is in the A,A position
c smallest (sigma3) is in the C,C position
c
       if (H(1,1).gt.H(2,2).and.H(1,1).gt.H(3,3)) then
         A=1
	 if (H(2,2).gt.H(3,3)) then
           B=2
           C=3
	 else
	   C=2
	   B=3
	 endif
       elseif (H(2,2).gt.H(1,1).and.H(2,2).gt.H(3,3)) then
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
c define the principal stresses (eigenvalues)
c  
      sig1=H(A,A)
      sig2=H(B,B)
      sig3=H(C,C)
c
c define the principal stress axes (eigenvectors)
c
c      v1 = (/ eigen(1,A), eigen(2,A), eigen(3,A) /)
c      v2 = (/ eigen(1,B), eigen(2,B), eigen(3,B) /)
c      v3 = (/ eigen(1,C), eigen(2,C), eigen(3,C) /)
c
c for g77 (ppc) fill arrays itteratively
c
      do 299 i=1,3
        v1(i) = eigen(i,A)
        v2(i) = eigen(i,B)
        v3(i) = eigen(i,C)
 299  continue
c
c output components of v1 and v3
c
      v1x = v1(1)
      v1y = v1(2)
      v1z = v1(3)
      v3x = v3(1)
      v3y = v3(2)
      v3z = v3(3)
      
c
c check that vectors are orthogonal and normal
c
c      rnrm1=v1x*v1x+v1y*v1y+v1z*v1z
c      rnrm2=v2x*v2x+v2y*v2y+v2z*v2z
c      rdot=v1x*v2x+v1y*v2y+v1z*v2z     
c
c print stuff out for debugging
c      
c      write(*,*),Txx,Tyy,Tzz,Txy,Txz,Tyz
      
c      write(*,*)A,B,C
c      write(*,*)eigen
c      write(*,*)v1
c      write(*,*)v2
c      write(*,*)rnrm1,rnrm2,rdot
      
      return
      end
c

c
      subroutine airycoef(kx,ky,z,H,fk,cub,cvb,cwb,cdwdz)
c
c  input 
c    kx  - x wavenumber             [1/m]
c    ky  - y wavenumber             [1/m]
c    z   - observation plane z < 0  [m]
c    H   - thickness of upper layer [m]
c    fk  - load on top of plate     [Pa m]
c
c  output
c    cub - x-displacement                      [m^2]
c    cvb - y-displacement                      [m^2]
c    cwb - z-displacement                      [m^2]
c    cdwdz - derivative of w with respect to z [m]
c
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      implicit real*8 (k)
      complex*8 fk,gk
      complex*8 A1,B1,C1,D1
      complex*8 arg1,arg2,arg3,arg4
c
c  routine to compute the greens function for point load on
c  top and bottom of a plate (top load = bottom load)
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
c  for the zero wavenumber, return zero
c
      if (kx.eq.0.and.ky.eq.0.) then
        cub=cmplx(0.,0.)
        cvb=cmplx(0.,0.)
        cwb=cmplx(0.,0.)
        cdwdz=cmplx(0.,0.)
        return
      endif

c
c  Compute the coefficients: A1 B1 C1 D1
c
      kh2=kx*kx+ky*ky
      kh=sqrt(kh2)
      beta=2*pi*kh

      alpha=(rlam1+rmu1)/(rlam1+2*rmu1)
      eta=(rlam1+rmu1)/(3*rlam1+4*rmu1)

      bH=beta*H
      if(bH.lt.-50.) then
        epH=0.
	ep2H=0.
      else
        epH=exp(bH)
	ep2H=exp(2*bH)
      endif

      denom=beta**3*(-eta-alpha*(-1+alpha+eta)
     +       +ep2H*(eta+alpha*(-1+alpha+eta))
     +       +2*alpha*beta*epH*(alpha-eta)*H)*rmu1
 
      A1=      eta*    (1-epH+alpha*(-2+epH*(2+beta*H)))/denom*fk
      B1=-alpha*eta*   (-1 + epH)			/denom*fk
      C1=-     eta*epH*(1-epH+alpha*(-2+2*epH+beta*H))  /denom*fk
      D1=-alpha*eta*epH*(-1 + epH)			/denom*fk

c
c  Compute the displacements: U V W dW/dz
c
      arg=beta*z
      if(arg.lt.-50.) then
        ep=0.
        en=exp(50.)
      else
        ep=exp(arg)
        en=exp(-arg)
      endif
      
      arg1=beta*ep*(A1+B1+B1*beta*z)-beta*en*(C1-D1+beta*D1*z)
      arg2=ep*(A1*alpha-2*B1+2*alpha*B1+alpha*B1*beta*z)
     +    +en*(C1*alpha+2*D1-2*alpha*D1+alpha*D1*beta*z)
      arg3=ep*(A1*alpha-2*B1+3*alpha*B1+alpha*B1*beta*z)
     +    -en*(C1*alpha+2*D1-3*alpha*D1+alpha*D1*beta*z)
      arg4=beta**2*ep*(A1+2*B1+B1*beta*z)
     +    +beta**2*en*(C1-2*D1+beta*D1*z)

      cub=cmplx(2.*pi*kx*alpha*aimag(arg1),-2.*pi*kx*alpha*real(arg1))
      cvb=cmplx(2.*pi*ky*alpha*aimag(arg1),-2.*pi*ky*alpha*real(arg1))
      cwb  =cmplx(-beta**2*real(arg2),-beta**2*aimag(arg2))
      cdwdz=cmplx(-beta**3*real(arg3),-beta**3*aimag(arg3))
      cdudz=cmplx(2.*pi*kx*alpha*aimag(arg4),-2.*pi*kx*alpha*real(arg4))
      cdvdz=cmplx(2.*pi*ky*alpha*aimag(arg4),-2.*pi*kx*alpha*real(arg4))

c      cub=cmplx(0.,-2.*pi*kx*alpha*arg1)
c      cvb=cmplx(0.,-2.*pi*ky*alpha*arg1)
c      cwb  =cmplx(-beta**2*arg2,0.)
c      cdwdz=cmplx(-beta**3*arg3,0.)

      return
      end
c

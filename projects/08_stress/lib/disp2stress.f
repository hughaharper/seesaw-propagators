c
      subroutine disp2stress(kx,ky,cu,cv,cw,cdudz,cdvdz,cdwdz,
     +                             cTxx,cTyy,cTzz,cTxy,cTxz,cTyz)
c
c  input 
c    kx    - x wavenumber                      [1/m]
c    ky    - y wavenumber                      [1/m]
c    cu    - x-displacement                    [m^2]
c    cv    - y-displacement                    [m^2]
c    cw    - z-displacement                    [m^2]
c    cdudz - derivative of u with respect to z [m]
c    cdvdz - derivative of v with respect to z [m]
c    cdwdz - derivative of w with respect to z [m]
c  output
c    cTxx - x normal stress [Pa m]
c    cTyy - y normal stress [Pa m]
c    cTzz - z normal stress [Pa m]
c    cTxy - xy shear stress [Pa m]
c    cTxz - xz shear stress [Pa m]
c    cTyz - yz shear stress [Pa m]
c
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      implicit real*8 (k)
c
c  routine to compute the greens function for point load on
c  top and bottom of a plate (top load = bottom load)
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
c  compute the horizontal derivatives, name the vertical ones
c  these derivatives in wavenumber domain have units [m]
c
      cux=cmplx(0.,2.*pi*kx)*cu
      cvy=cmplx(0.,2.*pi*ky)*cv
      cuy=cmplx(0.,2.*pi*ky)*cu
      cvx=cmplx(0.,2.*pi*kx)*cv
      cwx=cmplx(0.,2.*pi*kx)*cw
      cwy=cmplx(0.,2.*pi*ky)*cw
      cwz=cdwdz
      cuz=cdudz
      cvz=cdvdz
c
c  compute normal stresses
c
      cTxx=(rlam1+2*rmu1)*cux+rlam1*(cvy+cwz)
      cTyy=(rlam1+2*rmu1)*cvy+rlam1*(cux+cwz)
      cTzz=(rlam1+2*rmu1)*cwz+rlam1*(cux+cvy)
c
c  compute shear stresses
c
      cTxy=rmu1*(cuy+cvx)
      cTxz=rmu1*(cuz+cwx)
      cTyz=rmu1*(cvz+cwy)
c
c  write out values to debug
c
c      write(*,*),kx,ky,cTxx,cTyy,cTzz
      
      
      
      return
      end
c

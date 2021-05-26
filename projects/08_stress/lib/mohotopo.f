c
      subroutine mohotopo(kx,ky,Te,rhoc,rhom,fk,gk)
c
c  input 
c    kx   - x wavenumber              [1/m]
c    ky   - y wavenumber              [1/m]
c    Te   - elastic thickess of plate [m]
c    rhoc - crustal density           [kg/m^3]
c    rhom - mantle density            [kg/m^3]
c    g    - gravity                   [m/s^2]
c    fk   - load on top of plate      [Pa m]
c
c  output
c    gk   - load at base of plate     [Pa m]
c
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      implicit real*8 (k)
      complex*8 fk,gk
c
c  routine to compute the moho shape given topo shape,
c  elastic thickness, and material parameters
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
c  define the Flexural Rigidity and young's modulus
c
      young=rmu1*(3*rlam1+2*rmu1)/(rmu1+rlam1)
      rnu1=rlam1/2./(rlam1+rmu1)
      D=young*Te**3/12./(1-rnu1**2)
c
c  compute the transfer function
c
      kh2=kx*kx+ky*ky
      kh=sqrt(kh2)
      beta=2*pi*kh
      trans=1./(1.+ D*beta**4/grv/(rhom-rhoc))
c      rlamf=2*pi*(D/grv/(rhom-rhoc))**0.25
c
c  compute moho value
c
      gk=fk*trans
c
c check the material properties
c
c      write(*,*)rlam1,rmu1,rnu1,young,D,Te,rhoc,rhom,g
c      write(*,*),kh,beta,trans,gk
c      print *,Te,D,rlamf
      return
      end
c

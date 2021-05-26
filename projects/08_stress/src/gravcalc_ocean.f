c
      subroutine gravcalc_ocean(kx,ky,H,Te,rhoc,rhom,aobs,fk,zs,grval)
c
c  input 
c    kx   - x wavenumber              [1/m]
c    ky   - y wavenumber              [1/m]
c    H    - plate thickness           [m]
c    Te   - elastic thickess of plate [m]
c    rhoc - crustal density           [kg/m^3]
c    rhom - mantle density            [kg/m^3]
c    g    - gravity                   [m/s^2]
c    aobs - elevation of gravity obs. [m]
c    fk   - topo value               ??? [Pa m]
c    zs   - mean value of topo       ??? [Pa m]
c         this should be before it was filtered, but
c         for now this is the mean computed for the input, which should = 0
c
c  output
c    grval   - gravity in complex mGal     [Pa m]
c
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      implicit real*8 (k)
      complex*8 fk,gk,gsurf,gmoho,grval
c
c  routine to compute the moho shape given topo shape,
c  elastic thickness, and material parameters
c
      common/plate/rlam1,rlam2,rmu1,rmu2,rho,rc,rd,alph,pi,grv
c
c  define the Flexural Rigidity and young's modulus
c
      Gcns=2.*pi*6.673e-11
      young=rmu1*(3*rlam1+2*rmu1)/(rmu1+rlam1)
      rnu1=rlam1/2./(rlam1+rmu1)
      D=young*Te**3/12./(1-rnu1**2)
      rhow=1025.
c
c  define wavenumbers
c
      kh2=kx*kx+ky*ky
      kh=sqrt(kh2)
      beta=2*pi*kh
c
c  compute the transfer function
c       
c      trans=-rhoc/(rhom-rhoc)/(1.+ D*beta**4/grv/(rhom-rhoc))
      trans=-1./(1.+ D*beta**4/grv/(rhom-rhoc))
c
c  calculate gravity from surface and moho contributions
c       
c      gsurf=Gcns*fk*rhoc*exp(-beta*(aobs-zs))
cc      gmoho=Gcns*fk*trans*(rhom-rhoc)*exp(-beta*(aobs-zs-H))
c      gmoho=Gcns*fk*rhoc*exp(-beta*(aobs-zs-H))*trans
      gsurf=Gcns*fk*(rhoc-rhow)*exp(-beta*(aobs-zs))
      gmoho=Gcns*fk*(rhoc-rhow)*exp(-beta*(aobs-zs-H))*trans
c
c add surface and moho to get total gravity
c
      grval=(gsurf+gmoho)*1.e5
c
c check the material properties
c
c      write(*,*)kx,ky,H,Te,rhoc,rhom,g,aobs,fk,zs
c      write(*,*),kh,gsurf,gmoho,grval
c      write(*,*)kx,ky,fk,grval

      return
      end
c

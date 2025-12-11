      program example2
c nucleon-nucleon potential of Reid (T=1, soft core)
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16(c)
      parameter(nc=2,nmax=100,nc0=1)
      dimension cpot(nmax,nc,nc),lval(nc),eta(nc),qk(nc),zrma(nmax),
     1 cu(nc,nc),cf(nmax,nc,nc0),nvc(nc0)
      dimension fc(500),dfc(500),gc(500),dgc(500),cfx(nc)
      logical twf
      data pi/3.1415926535d0/,nvc(1)/1/
      data a1,a2,z1,z2/1,1,1,0/
c a1,a2,z1,z2=masses and charges of the colliding nuclei
      rmu=a1*a2/(a1+a2)
      ze=z1*z2*1.44d0
      hm=20.736d0/rmu
      eta0=ze/(2*hm)
c j=total spin of the proton-neutron system
c nr=number of Lagrange functions per interval
c ns=number of interval (ns=1 is propagation is not used)
c rmax=channel radius a
    1 read*,j,nr,ns,rmax
      if(j.lt.0)stop
      twf=rmax.lt.0
      rmax=abs(rmax)
      if(twf)open(1,file='wave_functions.txt')
      print1004,j,nr,ns,rmax
 1004 format(/,'total angular momentum=',i3,/,
     1 'number of basis functions per interval=',i3,/,
     2 'number of intervals=',i3,/,
     3 'channel radius=',f9.4)
      if(nr*ns.gt.nmax)then
      print*,'nmax too small'
      stop
      end if
      call rmat_ini(nr,ns,rmax,zrma)
      lval(1)=j-1
      lval(2)=j+1
      do 4 i=1,nr*ns
      r=zrma(i)
c Reid neutron-proton potential (T=1, soft core)
      xx=r*0.7d0
      zz=exp(-xx)
      vcen=(-10.463d0*zz+105.468d0*zz**2-3187.8d0*zz**4+9924.3d0*zz**6)
     1 /xx
      vtens=-10.463d0*((1+3/xx+3/xx**2)*zz-(12/xx+3/xx**2)*zz**4)/xx+
     1 351.77d0*zz**4/xx-1673.5d0*zz**6/xx
      vls=708.91d0*zz**4/xx-2713.1d0*zz**6/xx
      cpot(i,1,1)=vcen-2*(j-1)*vtens/(2*j+1)+(j-1)*vls
      cpot(i,1,2)=6*vtens*sqrt(j*(j+1.0d0))/(2*j+1)
      cpot(i,2,1)=cpot(i,1,2)
      cpot(i,2,2)=vcen-2*(j+2)*vtens/(2*j+1)-(j+2)*vls
    4 continue
      cpot=cpot/hm
    2 read*,ne,e0,estep
      print1005,ne,e0,estep
 1005 format('Number of energies:',i3,/,
     1 'Initial energy:',f10.4,/,
     2 'Energy step:',f10.4)
c ne=number of energies 
c e0=first energy (in MeV)
c estep=energy step (in MeV)
      if(ne.eq.0)go to 1
      do 3 ie=1,ne
      ecm=e0+(ie-1)*estep
      qk=sqrt(ecm/hm)
      eta=0
      call rmatrix(nc,lval,qk,eta,rmax,nr,ns,cpot,cu,nmax,nc,nopen,
     1 twf,cf,nmax,nc,nc0,nvc,0,cc)
      del1=atan2(aimag(cu(1,1)),real(cu(1,1)))/2
      del2=atan2(aimag(cu(2,2)),real(cu(2,2)))/2
      del12=abs(cu(1,2))
      print1000,ecm,del1,del2,del12
 1000 format('E (MeV)=',f5.1,' phase shift (rad.)=',2es12.4,
     1 ' eta_12=',es12.4)
      if(.not.twf)go to 3
c Write the wave function on file 'wave_function.txt'
      if(.not.twf)go to 3
      write(1,1001)ecm,nc0
 1001 format('Wave function, E=',f8.3,' entrance channel=',i2)
c Internal wave function (at the mesh points)
      write(1,*)'Internal wave functions'
      do 10 ir=1,nr*ns
      write(1,1002)zrma(ir),(cf(ir,1:nc,nvc(ic0)),ic0=1,nc0)
 1002 format(f8.3,10(2es12.4,2x))
   10 continue
c External wave function (from rmax to rmax+2 by step of 0.2 fm)
c This can be changed by the user
      write(1,*)'External wave functions'
      do 11 ir=0,11
      r=zrma(nr*ns)+0.2*ir
      nop=0
      nclo=nopen
      do 12 iv=1,nc
      xl=lval(iv)
      call coulfg(qk(iv)*r,eta(iv),xl,xl,fc,gc,dfc,dgc,1,0,ifail)
      if(qk(iv).gt.0)then
      ll=lval(iv)+1
      nop=nop+1
      co=dcmplx(gc(ll),fc(ll))*sqrt(qk(1)/qk(iv))
      cfx(iv)=-cu(nop,1)*co
      if(iv.eq.1)cfx(iv)=cfx(iv)+conjg(co)
      else
      nclo=nclo+1
      iw=0
      call whit(eta(iv),r,-qk(iv),-qk(iv)**2,ll,fc,dfc,iw)
      cfx(iv)=cu(nclo,1)*fc(ll+1)
      end if
   12 continue
      write(1,1002)r,cfx(1:nc)
   11 continue
    3 continue
      go to 2
      end

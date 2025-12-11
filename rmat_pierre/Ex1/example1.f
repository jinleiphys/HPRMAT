      program example1
c alpha-208Pb optical potential of Goldring et al. , Phys. Lett. B32 (1970) 465
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16(c)
      parameter(nc=1,nmax=100,nc0=1)
      dimension cpot(nmax,nc,nc),lval(nc),eta(nc),qk(nc),zrma(nmax),
     1 cu(nc,nc),cf(nmax,nc,nc0),nvc(nc0)
      dimension fc(500),dfc(500),gc(500),dgc(500),cfx(nc)
      logical twf
      data pi/3.1415926535d0/,nvc(1)/1/
      data a1,a2,z1,z2/208,4,82,2/
c a1,a2,z1,z2=masses and charges of the colliding nuclei
      rmu=a1*a2/(a1+a2)
      ze=z1*z2*1.44d0
      hm=20.736d0/rmu
      eta0=ze/(2*hm)
c l=relative angular momentum
c nr=number of Lagrange functions per interval
c ns=number of interval (ns=1 is the propagation is not used)
c rmax=channel radius a
    1 read*,l,nr,ns,rmax
      if(l.lt.0)stop
      twf=rmax.lt.0
      rmax=abs(rmax)
      if(twf)open(1,file='wave_functions.txt')
      print1004,l,nr,ns,rmax
 1004 format(/,'total angular momentum=',i3,/,
     1 'number of basis functions per interval=',i3,/,
     2 'number of intervals=',i3,/,
     3 'channel radius=',f9.4)
      if(nr*ns.gt.nmax)then
      print*,'nmax too small'
      stop
      end if
      call rmat_ini(nr,ns,rmax,zrma)
c Woods-Saxon potential with
c rn=range
c an=diffuseness
c amplitude=-(100+10*I) MeV
      an=0.5803d0
      rn=1.1132d0*(a1**(1.0d0/3)+a2**(1.0d0/3))
      do 4 i=1,nr*ns
      r=zrma(i)
      xx=1+exp((r-rn)/an)
      cvn=-dcmplx(100,10)/xx
      vc=ze/r
      cpot(i,1,1)=(cvn+vc)/hm
    4 continue
    2 read*,ne,e0,estep
      print1005,ne,e0,estep
 1005 format('Number of energies:',i3,/,
     1 'Initial energy:',f10.4,/,
     2 'Energy step:',f10.4)
c ne=number of energies 
c e0=first energy (in MeV)
c estep=energy step (in MeV)
      if(ne.eq.0)go to 1
      lval(1)=l
      do 3 ie=1,ne
      ecm=e0+(ie-1)*estep
      qk(1)=sqrt(ecm/hm)
      eta=eta0/qk
      call rmatrix(nc,lval,qk,eta,rmax,nr,ns,cpot,cu,nmax,nc,nopen,
     1 twf,cf,nmax,nc,nc0,nvc,0,cc)
      et=abs(cu(1,1))
      del=atan2(aimag(cu(1,1)),real(cu(1,1)))/2
c     print1000,ecm,del*180/pi,et
c1000 format('E (MeV)=',f8.3,' phase shift (deg.)=',2es12.4)
      print1000,ecm,cu(1,1)
 1000 format('E (MeV)=',f8.3,'   Collision matrix=',2es12.4)
 
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

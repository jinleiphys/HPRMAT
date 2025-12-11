      program example3
c 16O+44Ca potential of Rhoades-Brown et al., Phys. Rev. C21 (1980) 2417
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16(c)
      parameter(nc=4,nmax=100,nc0=1,ncha=2,npmax=2000)
      dimension ispi(ncha),th(ncha)
      dimension cpot(nmax,nc,nc),lval(nc),eta(nc),qk(nc),zrma(nmax),
     1 cu(nc,nc),cf(nmax,nc,nc0),nvc(nc0)
      dimension fc(500),dfc(500),gc(500),dgc(500),cfx(nc),xfac(nc,nc)
      dimension cwftab(npmax)
      logical twf
      data pi/3.1415926535d0/,nvc(1)/1/
      data a1,a2,z1,z2/44,16,20,8/
c a1,a2,z1,z2=masses and charges of the colliding nuclei
      data ispi,th/0,2,0,1.156d0/
      data v0,w0,r0,an/110,20,1.2d0,0.5d0/
c beta_N=0.4, beta_C=0
      data lam,ben,bec/2,0.40,0.0/
      pi4=4*acos(-1.0d0)
      rn=r0*(a1**(1.0d0/3)+a2**(1.0d0/3))
      rv=r0*a1**(1.0d0/3)
      rmu=a1*a2/(a1+a2)
      ze=z1*z2*1.44d0
      hm=20.736d0/rmu
      eta0=ze/(2*hm)
c j=total angular momentum
c nr=number of Lagrange functions per interval
c ns=number of interval (ns=1 is the propagation is not used)
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
      xfac=0
      nv=0
c Coefficients required for the coupling potentilas (PRC21, 1980, 2417)
      do 20 n=1,ncha
      ii=ispi(n)
      do 20 l=abs(j-ii),j+ii,2
      nv=nv+1
      nv2=0
      do 21 n2=1,ncha
      ii2=ispi(n2)
      do 21 l2=abs(j-ii2),j+ii2,2
      nv2=nv2+1
      fac=(2*l+1)*(2*l2+1)*(2*ii+1)*(2*ii2+1)*(2*lam+1)
      call troisj(2*ii2,2*lam,2*ii,0,0,0,t1)
      call troisj(2*l,2*lam,2*l2,0,0,0,t2)
      call sixj(2*ii,2*l,2*j,2*l2,2*ii2,2*lam,s6)
      xfac(nv,nv2)=t1*t2*s6*sqrt(fac/pi4)
      if(mod(abs(l-l2)/2,2).eq.1)xfac(nv,nv2)=-xfac(nv,nv2)
   21 continue
   20 continue
      if(mod(j+lam,2).eq.1)xfac=-xfac
      do 4 i=1,nr*ns
      r=zrma(i)
c Woods-Saxon potential with nuclear deformation
c rn=range
c an=diffuseness
c v0=real amplitude
c w0=imaginary amplitude
c ben=nuclear deformation parameter
c bec=Coulomb defromation parameter
      xx1=exp((r-rn)/an)
      xx=1+xx1
      cvn=-dcmplx(v0,w0)/xx
      if(r.le.rn)vc=(3-(r/rn)**2)/2/rn
      if(r.gt.rn)vc=1/r
      cpot(i,:,:)=0
      do 5 iv=1,nv
      cpot(i,iv,iv)=cvn+ze*vc
    5 continue
      cvn=-dcmplx(v0,w0)*ben*rv*xx1/(an*xx**2)
      if(r.le.rn)vc=r**lam/rn**(lam+lam+1)
      if(r.gt.rn)vc=1/r**(lam+1)
      vc=3*ze*vc*bec*rv**lam/(2*lam+1)
      cpot(i,1:nv,1:nv)=cpot(i,1:nv,1:nv)+xfac(1:nv,1:nv)*
     1 (cvn+vc)
    4 continue
      cpot=cpot/hm
      read*,npoin,h 
    2 read*,ne,e0,estep
c ne=number of energies 
c e0=first energy (in MeV)
c estep=energy step (in MeV)
      if(ne.eq.0)go to 1
      print1005,ne,e0,estep
 1005 format('Number of energies:',i3,/,
     1 'Initial energy:',f10.4,/,
     2 'Energy step:',f10.4)
      do 3 ie=1,ne
      ecm=e0+(ie-1)*estep
      nv=0
      do 22 n=1,ncha
      ii=ispi(n)
      do 22 l=abs(j-ii),j+ii,2
      nv=nv+1
      lval(nv)=l
      qk(nv)=sqrt((ecm-th(n))/hm)
      eta(nv)=eta0/qk(nv)
   22 continue
      call rmatrix(nv,lval,qk,eta,rmax,nr,ns,cpot,cu,nmax,nc,nopen,
     1 twf,cf,nmax,nc,nc0,nvc,0,cc)
      print1000,ecm,abs(cu(1:nopen,1))
      print1003,ecm,atan2(aimag(cu(1:nopen,1)),real(cu(1:nopen,1)))/2
 1000 format('E (MeV)=',f8.3,' Amplitude=         ',8es12.4)
 1003 format('E (MeV)=',f8.3,' phase shift (rad.)=',8es12.4)
c Write the wave function on file 'wave_function.txt'
      if(.not.twf)go to 3
      write(1,1001)ecm,nc0
 1001 format('Wave function, E=',f8.3,' entrance channel=',i2)
c Internal wave function (at the mesh points)
      write(1,*)'Internal wave functions'
      do 10 ir=1,nr*ns
      write(1,1002)zrma(ir),(cf(ir,1:nv,nvc(ic0)),ic0=1,nc0)
 1002 format(f8.3,10(2es12.4,2x))
   10 continue
c External wave function (from rmax to rmax+2 by step of 0.2 fm)
c This can be changed by the user
      write(1,*)'External wave functions'
      do 11 ir=0,11
      r=zrma(nr*ns)+0.2*ir
      nop=0
      nclo=nopen
      do 12 iv=1,nv
      ll=lval(iv)
      if(qk(iv).gt.0)then
      xl=lval(iv)
      call coulfg(qk(iv)*r,eta(iv),xl,xl,fc,gc,dfc,dgc,1,0,ifail)
      nop=nop+1
      co=dcmplx(gc(ll+1),fc(ll+1))*sqrt(qk(1)/qk(iv))
      cfx(iv)=-cu(nop,1)*co
      if(iv.eq.1)cfx(iv)=cfx(iv)+conjg(co)
      else
      nclo=nclo+1
      iw=0
      call whit(eta(iv),r,-qk(iv),-qk(iv)**2,ll,fc,dfc,iw)
      cfx(iv)=cu(nclo,1)*fc(ll+1)
      end if
   12 continue
      write(1,1002)r,cfx(1:nv)
   11 continue
      print*,'Wave function at a fixed step size'
      nom=1
      do 13 iv=1,nv
      print*,'Channel ',iv
      call wf_print(nv,lval,qk,eta,rmax,nr,ns,cu,nmax,nopen,cf,
     1 nmax,nc,zrma,iv,nom,npoin,h,cwftab)
      do 14 ir=1,npoin
      print1002,ir*h,cwftab(ir)
   14 continue
   13 continue
    3 continue
      go to 2
      end
*CLEBS
      SUBROUTINE CLEBS (L1,L2,L3,M1,M2,M3,Q)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION FT(0:1000)
      SAVE  
      DATALMEM,(FT(I),I=0,10)/9,2*1.0d0,2.0d0,6.0d0,24.0d0,
     1 120.0d0,720.0d0,5040.0d0,40320.0d0,362880.0d0,3628800.0d0/
      IS=(L1+L2+M1-M2)/2
      K3=-M3
      Q1=L3+1 
   12 Q=0  
      I1=L1 
      I2=L2 
      I3=L3 
      K1=M1 
      K2=M2 
      IF(K1+K2+K3.NE.0)RETURN 
      L=I1+I2+I3
      IF(MOD(L,2).EQ.1)GOTO8
      L=L/2 
      IF(L.LE.LMEM)GOTO6
      DO10I=LMEM,L
   10 FT(I+1)=(I+1)*FT(I) 
      LMEM=L
    6 J1=ABS(K1) 
      J2=ABS(K2) 
      J3=ABS(K3) 
      IF(I1.LT.ABS(I2-I3).OR.I1.GT.I2+I3)RETURN
      IF(J1+J2.EQ.0)GOTO11
      J1=I1-J1
      J2=I2-J2
      J3=I3-J3
      IF(J1)8,2,1 
    2 IF(I1.NE.0)GOTO1
      IF(J2.LT.0)GOTO8
   13 IF(J3.LT.0)GOTO8
    4 Q=SQRT(Q1/(I2+1)) 
      IS=IS+(I2-K2)/2 
      IF(MOD(IS,2).EQ.1)Q=-Q
      RETURN
    1 IF(J2.GT.J1)GOTO3 
      IF(J2.LT.0)GOTO8
      IS=IS+L 
      J1=J2 
      K1=K2 
      K2=M1 
      I1=I2 
      I2=L1 
      IF(I1.EQ.0)GOTO13 
    3 IF(J3.GT.J1)GOTO5 
      IF(J3.LT.0)GOTO8
      IS=IS+L 
      J1=K3 
      K3=K1 
      K1=J1 
      I3=I1 
      I1=L3 
      J1=J3 
      IF(I1.EQ.0)GOTO4
    5 IF(K1.GE.0)GOTO9
      K1=-K1
      K2=-K2
      K3=-K3
      IS=IS+L 
    9 CONTINUE
      Q1=Q1*FT(L-I3)/FT(L-I1)/FT(L-I2)/FT(L+1)
      I2=(I2+K2)/2
      I3=(I3+K3)/2
      K2=I2-K2
      K3=I3-K3
      J1=J1/2 
      I1=J1+K1
      J2=I3-K2
      J3=MAX(J2,0) 
      IS=IS+I1+K2 
      X=0  
      DO7I=J3,J1
    7 X=-X+FT(I1+I)*FT(I2+I3-I)/FT(J1-I)/FT(I3-I)/FT(I-J2)/FT(I)
      Q=X*       SQRT(Q1*FT(J1)*FT(K2)*FT(K3)*FT(I3)/FT(I1)/FT(I2)) 
      IF(MOD(IS,2).EQ.1)Q=-Q
      RETURN
    8 Q=0  
      PRINT1010,L1,L2,L3,M1,M2,M3 
 1010 FORMAT(10H ERREUR 3J,2(3X,3I3)) 
      RETURN
   11 IF(MOD(L,2).EQ.1)RETURN 
      I1=L-I1 
      I2=L-L2 
      I3=L-L3 
      Q=SQRT(FT(I1)*FT(I2)*FT(I3)/FT(L+1)*Q1) 
      I1=I1/2 
      I2=I2/2 
      I3=I3/2 
      L =L/2
      Q=Q*FT(L )/FT(I1)/FT(I2)/FT(I3) 
      IF(MOD(L +IS,2).EQ.1)Q=-Q 
      RETURN
      ENTRY TROISJI (L1,L2,L3,M1,M2,M3,Q) 
      ENTRY TROISJ (L1,L2,L3,M1,M2,M3,Q)  
      K3=M3 
      Q1=1 
      IS=0  
      GOTO12
      END 
*SIXJ 
      SUBROUTINE SIXJ(J1,J2,J3,L1,L2,L3,Q)
c Computes the 6j coefficient
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION M(7),M1(4),M2(4),M3(4),FT(0:1000)
      SAVE  
      DATA LMEM,(FT(I),I=0,10)/9,2*1.0d0,2.0d0,6.0d0,24.0d0,
     1 120.0d0,720.0d0,5040.0d0,40320.0d0,362880.0d0,3628800.0d0/  
      ENTRYSIXJI(J1,J2,J3,L1,L2,L3,Q) 
      I1=J1 
      I2=J2 
      I3=J3 
      K1=L1 
      K2=L2 
      K3=L3 
      IS=0  
   24 Q=0  
      M(1)=I1+I2+I3 
      M(2)=I1+K2+K3 
      M(3)=K1+I2+K3 
      M(4)=K1+K2+I3 
      DO 17 I=1,4 
      IF(MOD(M(I),2).EQ.1) GO TO 8
   17 CONTINUE
      L=MAX(I1+I2+K1+K2,I1+I3+K1+K3,I2+I3+K2+K3) 
      L=L/2 
      IF(L.LE.LMEM) GO TO 6 
      DO 10 I=LMEM,L  
   10 FT(I+1)=FT(I)*(I+1) 
      LMEM=L
    6 IF(I1.LT.ABS(I2-I3).OR.I1.GT.I2+I3) RETURN 
      IF(I1.LT.ABS(K2-K3).OR.I1.GT.K2+K3) RETURN 
      IF(K1.LT.ABS(I2-K3).OR.K1.GT.I2+K3) RETURN 
      IF(K1.LT.ABS(K2-I3).OR.K1.GT.K2+I3) RETURN 
      IF(I1) 8,2,1
    2 IF(I2.LT.0) GO TO 8 
    9 IF(I3.LT.0) GO TO 8 
   14 IF(K1.LT.0) GO TO 8 
   19 IF(K2.LT.0) GO TO 8 
   23 IF(K3.LT.0) GO TO 8 
   27 Q=SQRT(1.0d0/(I2+1)/(K2+1))  
      IS=(I2+K2+K1)/2+IS
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    1 IF(I1.GT.1) GO TO 3 
      IF(I2.LT.0) RETURN
   12 IF(I3.LT.0) RETURN
   16 IF(K1.LT.0) RETURN
   21 IF(K2.LT.0) RETURN
   25 IF(K3.LT.0) RETURN
   28 IF(I2.LT.I3) GO TO 4
      IC=I2 
      I2=I3 
      I3=IC 
      IC=K2 
      K2=K3 
      K3=IC 
    4 IF(K2.GT.K3) GO TO 5
      I11=I1+K1+I2-K2 
      I11=I11/2 
      I12=I11-I2+K2 
      Q=SQRT(I11*I12*1./I3/(I3+1)/K3/(K3+1))
      IS =I11+K2+IS 
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    5 I11=K3-K1+I2
      I11=I11/2+1 
      I12=I11+K1+1
      Q=SQRT(I11*I12*1.0d0/I3/(I3+1)/K2/(K2+1))
      IS =I12-1+IS
      IF(MOD(IS ,2).EQ.1) Q=-Q  
      RETURN
    3 IF(I2.GE.I1) GO TO 7
      IF(I2.LT.0) GO TO 8 
      IC=I2 
      I2=I1 
      I1=IC 
      IC=K1 
      K1=K2 
      K2=IC 
      IF(I1.EQ.0) GO TO 9 
      IF(I1.EQ.1) GO TO 12
    7 IF(I3.GE.I1) GO TO 13 
      IF(I3.LT.0) GO TO 8 
      IC=I3 
      I3=I1 
      I1=IC 
      IC=K3 
      K3=K1 
      K1=IC 
      IF(I1.EQ.0) GO TO 14
      IF(I1.EQ.1) GO TO 16
   13 IF(K1.GE.I1) GO TO 18 
      IF(K1.LT.0) GO TO 8 
      IC=K1 
      K1=I1 
      I1=IC 
      IC=K2 
      K2=I2 
      I2=IC 
      IF(I1.EQ.0) GO TO 19
      IF(I1.EQ.1) GO TO 21
   18 IF(K2.GE.I1) GO TO 22 
      IF(K2.LT.0) GO TO 8 
      IC=K2 
      K2=I1 
      I1=IC 
      IC=K1 
      K1=I2 
      I2=IC 
      IF (I1.EQ.0) GO TO 23 
      IF(I1.EQ.1) GO TO 25
   22 IF(K3.GE.I1) GO TO 26 
      IF(K3.LT.0) GO TO 8 
      IC=K3 
      K3=I1 
      I1=IC 
      IC=K1 
      K1=I3 
      I3=IC 
      IF(I1.EQ.0) GO TO 27
      IF(I1.EQ.1) GO TO 28
   26 M1(4)=I3
      M1(1)=I3
      M1(3)=K3
      M1(2)=K3
      M2(2)=I1
      M2(1)=I1
      M2(4)=K1
      M2(3)=K1
      M3(3)=I2
      M3(1)=I2
      M3(4)=K2
      M3(2)=K2
      M(1)=I1+I2+I3 
      M(2)=I1+K2+K3 
      M(3)=K1+I2+K3 
      M(4)=K1+K2+I3 
      Q1=1 
      DO 11 I=1,4 
      M(I)=M(I)/2 
   11 Q1=FT(M(I)-M1(I))*FT(M(I)-M2(I))*FT(M(I)-M3(I))*Q1/FT(M(I)+1) 
      Q1=SQRT(Q1) 
      M1(1)=I1+K1 
      M1(2)=I2+K2 
      M1(3)=I3+K3 
      IC=M1(1)+M1(2)  
      M(5)=IC/2 
      IC=M1(2)+M1(3)  
      M(6)=IC/2 
      IC=M1(1)+M1(3)  
      M(7)=IC/2 
      MAXZ=MIN(M(5),M(6),M(7)) 
      MINZ=MAX(M(1),M(2),M(3),M(4))
      X=0  
      DO 15 I=MINZ,MAXZ 
      Q2=1 
      DO 20 J=1,7 
      IJ=I-M(J) 
      IF(J.GT.4) IJ=-IJ 
   20 Q2=Q2*FT(IJ)
      Q2=FT(I+1)/Q2 
   15 X=-X+Q2 
      Q=X*Q1
      IS=MAXZ+IS
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    8 PRINT 1010,J1,J2,J3,L1,L2,L3
 1010 FORMAT(10H ERREUR 6J,2(3X,3I3)) 
      RETURN
      ENTRY RACAH(J1,J2,J3,L1,L2,L3,Q)
      IS=(J1+J2+J3+L1)/2
      I1=J1 
      I2=J2 
      I3=L2 
      K1=L1 
      K2=J3 
      K3=L3 
      GO TO 24
      END 

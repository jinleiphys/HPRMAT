!------------------------------------------------------------------------------
! example4_hp.f90 - High-Performance version of Example 4
! 12C+alpha scattering with coupled channels
!------------------------------------------------------------------------------
program example4_hp
  use rmat_hp_mod
  implicit real*8(a,b,d-h,o-z)
  implicit complex*16(c)

  parameter(nc=12, nmax=100, nc0=1, ncha=3)
  dimension ispi(ncha), th(ncha)
  dimension cpot(nmax,nc,nc), lval(nc), eta(nc), qk(nc), zrma(nmax)
  dimension cu(nc,nc), cf(nmax,nc,nc0), nvc(nc0)
  dimension fc(500), dfc(500), gc(500), dgc(500), cfx(nc), xfac(nc,nc)
  logical twf

  data pi/3.1415926535d0/, nvc(1)/1/
  data a1, a2, z1, z2 /12, 4, 6, 2/
  ! a1,a2,z1,z2=masses and charges of the colliding nuclei
  data ispi, th /0, 2, 4, 0, 4.44d0, 14.08d0/
  data v0, w0, r0, an /110, 20, 1.2d0, 0.5d0/
  data lam, ben, bec /2, 0.58d0, 0/

  ! HIGH-PERFORMANCE CONFIGURATION
  ! solver_type = 1: Dense LAPACK (default)
  ! solver_type = 2: Mixed Precision
  ! solver_type = 3: Woodbury-Kinetic
  ! solver_type = 4: GPU cuSOLVER
  solver_type = 1

  write(*,*) '================================================'
  write(*,*) 'HPRMAT Example 4: 12C+alpha Scattering'
  write(*,*) '================================================'
  write(*,'(A,I1)') ' Solver type: ', solver_type
  write(*,*)

  pi4 = 4 * acos(-1.0d0)
  rn = r0 * (a1**(1.0d0/3) + a2**(1.0d0/3))
  rv = r0 * a1**(1.0d0/3)
  rmu = a1 * a2 / (a1 + a2)
  ze = z1 * z2 * 1.44d0
  hm = 20.736d0 / rmu
  eta0 = ze / (2 * hm)

1 read(*, *, end=999) j, nr, ns, rmax
  if (j < 0) stop
  twf = rmax < 0
  rmax = abs(rmax)
  if (twf) open(1, file='wave_functions_hp.txt')

  write(*,1004) j, nr, ns, rmax
1004 format(/, 'Total angular momentum: ', i3, /, &
           'Lagrange functions per interval: ', i3, /, &
           'Number of intervals: ', i3, /, &
           'Channel radius: ', f9.4, ' fm')

  if (nr * ns > nmax) then
    write(*,*) 'ERROR: nmax too small'
    stop
  end if

  call rmat_ini(nr, ns, rmax, zrma)

  xfac = 0
  nv = 0
  ! Coefficients required for the coupling potentials (PRC21, 1980, 2417)
  do n = 1, ncha
    ii = ispi(n)
    do l = abs(j - ii), j + ii, 2
      nv = nv + 1
      nv2 = 0
      do n2 = 1, ncha
        ii2 = ispi(n2)
        do l2 = abs(j - ii2), j + ii2, 2
          nv2 = nv2 + 1
          fac = (2*l + 1) * (2*l2 + 1) * (2*ii + 1) * (2*ii2 + 1) * (2*lam + 1)
          call troisj(2*ii2, 2*lam, 2*ii, 0, 0, 0, t1)
          call troisj(2*l, 2*lam, 2*l2, 0, 0, 0, t2)
          call sixj(2*ii, 2*l, 2*j, 2*l2, 2*ii2, 2*lam, s6)
          xfac(nv, nv2) = t1 * t2 * s6 * sqrt(fac / pi4)
          if (mod(abs(l - l2) / 2, 2) == 1) xfac(nv, nv2) = -xfac(nv, nv2)
        end do
      end do
    end do
  end do
  if (mod(j + lam, 2) == 1) xfac = -xfac

  do i = 1, nr * ns
    r = zrma(i)
    ! Woods-Saxon potential with nuclear deformation
    xx1 = exp((r - rn) / an)
    xx = 1 + xx1
    cvn = -dcmplx(v0, w0) / xx
    if (r <= rn) vc = (3 - (r / rn)**2) / 2 / rn
    if (r > rn) vc = 1 / r
    cpot(i, :, :) = 0
    do iv = 1, nv
      cpot(i, iv, iv) = cvn + ze * vc
    end do
    cvn = -dcmplx(v0, w0) * ben * rv * xx1 / (an * xx**2)
    if (r <= rn) vc = r**lam / rn**(lam + lam + 1)
    if (r > rn) vc = 1 / r**(lam + 1)
    vc = 3 * ze * vc * bec * rv**lam / (2 * lam + 1)
    cpot(i, 1:nv, 1:nv) = cpot(i, 1:nv, 1:nv) + xfac(1:nv, 1:nv) * (cvn + vc)
  end do
  cpot = cpot / hm

2 read(*, *, end=1) ne, e0, estep
  if (ne == 0) go to 1

  write(*,1005) ne, e0, estep
1005 format('Number of energies: ', i3, /, &
           'Initial energy: ', f10.4, ' MeV', /, &
           'Energy step: ', f10.4, ' MeV')

  call cpu_time(t_start)

  do ie = 1, ne
    ecm = e0 + (ie - 1) * estep
    nv = 0
    do n = 1, ncha
      ii = ispi(n)
      do l = abs(j - ii), j + ii, 2
        nv = nv + 1
        lval(nv) = l
        ecm2 = ecm - th(n)
        qk(nv) = sqrt(abs(ecm2) / hm)
        eta(nv) = eta0 / qk(nv)
        if (ecm2 < 0) qk(nv) = -qk(nv)
      end do
    end do

    call rmatrix(nv, lval, qk, eta, rmax, nr, ns, cpot, cu, nmax, nc, nopen, &
                 twf, cf, nmax, nc, nc0, nvc, 0, cc, solver_type)

    write(*,1000) ecm, abs(cu(1:nopen, 1))
    write(*,1003) ecm, atan2(aimag(cu(1:nopen, 1)), real(cu(1:nopen, 1))) / 2
1000 format('E (MeV)=', f8.3, ' Amplitude=         ', 3(4es12.4, /, 36x))
1003 format('E (MeV)=', f8.3, ' phase shift (rad.)=', 3(4es12.4, /, 36x))

    if (.not. twf) cycle

    write(1,1001) ecm, nc0
1001 format('Wave function, E=', f8.3, ' entrance channel=', i2)

    write(1,*) 'Internal wave functions'
    do ir = 1, nr * ns
      write(1,1002) zrma(ir), (cf(ir, 1:nv, nvc(ic0)), ic0=1, nc0)
1002  format(f8.3, 10(2es12.4, 2x))
    end do

    write(1,*) 'External wave functions'
    do ir = 0, 11
      r = zrma(nr * ns) + 0.2d0 * ir
      nop = 0
      nclo = nopen
      do iv = 1, nv
        ll = lval(iv)
        if (qk(iv) > 0) then
          nop = nop + 1
          xl = lval(iv)
          call coulfg(qk(iv)*r, eta(iv), xl, xl, fc, gc, dfc, dgc, 1, 0, ifail)
          co = dcmplx(gc(ll+1), fc(ll+1)) * sqrt(qk(1) / qk(iv))
          cfx(iv) = -cu(nop, 1) * co
          if (iv == 1) cfx(iv) = cfx(iv) + conjg(co)
        else
          nclo = nclo + 1
          iw = 0
          call whit(eta(iv), r, -qk(iv), -qk(iv)**2, ll, fc, dfc, iw)
          cfx(iv) = cu(nclo, 1) * fc(ll+1)
        end if
      end do
      write(1,1002) r, cfx(1:nv)
    end do
  end do

  call cpu_time(t_end)
  write(*,'(A,F8.4,A)') 'Total time: ', t_end - t_start, ' seconds'
  write(*,*)

  go to 2

999 continue
  write(*,*) 'Calculation completed.'

end program example4_hp

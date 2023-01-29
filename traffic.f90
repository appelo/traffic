module FFT
  IMPLICIT NONE
  INCLUDE 'fftw3.f'                      ! needed for defining the plan
  INTEGER*8 :: forward,     backward     ! "plans" in FFTW
  ! work arrays
  double precision,    ALLOCATABLE, DIMENSION(:) :: u_r
  double complex,      ALLOCATABLE, DIMENSION(:) :: u_c
  ! wavenumber array
  double precision,    ALLOCATABLE, DIMENSION(:) :: k
contains

  SUBROUTINE setupfft(n,len)
    IMPLICIT NONE
    INCLUDE 'fftw3.f'           ! needed for defining the plan
    INTEGER i,n
    double precision pi,len
    pi = 4.D0*DATAN(1.D0)
    ! allocate work arrays
    ALLOCATE( u_r(n) )
    ALLOCATE( u_c(n/2+1) )
    ! allocate wavenumbers
    ALLOCATE( k(n/2+1))
    ! call initialization routine for FFTW 
    CALL dfftw_plan_dft_r2c_1d(forward ,n,u_r,u_c, FFTW_ESTIMATE)
    CALL dfftw_plan_dft_c2r_1d(backward,n,u_c,u_r, FFTW_ESTIMATE)
    ! compute wavenumbers
    DO i=1,(n/2+1)
       k(i) = DBLE(i-1)*2.d0*pi/len
    END DO
  END SUBROUTINE setupfft

  SUBROUTINE destroyfft
    CALL dfftw_destroy_plan(forward)
    CALL dfftw_destroy_plan(backward)
    DEALLOCATE(k,u_r,u_c)
  END SUBROUTINE destroyfft

  SUBROUTINE deriv(u,up,order,n)
    IMPLICIT NONE
    INTEGER :: n,order,i
    double precision, dimension(n) :: u,up
    u_r = u
    call dfftw_execute (forward)
    if (order .eq. 1) then
       DO i=1,(n/2+1)
          u_c(i) = (0.0D0,1.0D0)*k(i)*u_c(i)
       END DO
    END if
    if (order .eq. 2) then
       DO i=1,(n/2+1)
          u_c(i) = (-1.0D0,0.0D0)*k(i)*k(i)*u_c(i)
       END DO
    end if
    call dfftw_execute (backward)
    up = u_r/dble(n)
  end SUBROUTINE deriv

  SUBROUTINE filter(u,filter_style,n)
    IMPLICIT NONE
    INTEGER :: n,order,i,filter_style
    double precision, dimension(n) :: u
    double precision :: eta,filt
    if (filter_style.eq.0) then
       ! no filter
    else
       u_r = u
       call dfftw_execute (forward)
       if(filter_style .eq.1 ) then ! Anti Aliasing
          DO i=1,(n/2+1)
             if (i > floor(dble(n/2+1)*0.66d0)) u_c(i) = (0.0D0,0.d0)
          END DO
       ELSEif(filter_style .eq.2 ) then ! Cesaro filter
          DO i=1,(n/2+1)
             eta = k(i)/k(n/2+1)
             filt = 1-eta
             u_c(i) = filt*u_c(i)
          END DO
       ELSEif(filter_style .gt.2 ) then ! exp filter
          DO i=1,(n/2+1)
             eta = k(i)/k(n/2+1)
             filt = exp(-35.d0*eta**(2*filter_style))
             u_c(i) = filt*u_c(i)
          END DO
       end if
       call dfftw_execute (backward)
       u = u_r/dble(n)
    end if
  end SUBROUTINE filter

  Subroutine spectrum(u,up,n)
    IMPLICIT NONE
    INTEGER :: n,i
    double precision, dimension(n) :: u
    double precision, dimension(n) :: up
    u_r = u
    call dfftw_execute (forward)
    DO i=1,(n/2)
       up(i) = sqrt(real(u_c(i)*conjg(u_c(i))))
    END DO
    DO i=(n/2+1),n
       up(i) = k(i-n/2)
    END DO
  end SUBROUTINE spectrum

end module FFT

program TRAFFIC
  USE FFT
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: nmax = 512
  INTEGER, PARAMETER :: nmod = 20
 
  INTEGER, PARAMETER :: Filter_style = 0
  ! 0 : no filter
  ! 1 : de-aliasing (cut last 1/3)
  ! 2 : Cesaro filter
  ! > 2 : exponential of order 2*Filter_style
  double precision, parameter :: len = 10.d0
  double precision, PARAMETER :: tmax = 5.d0
  ! work arrays
  double precision,    ALLOCATABLE, DIMENSION(:) :: u,x,up,fh,us
  double precision,    ALLOCATABLE, DIMENSION(:,:) :: fu

  double precision :: a,pi,merr,t,dt
  double precision :: rkc(4),rnd
  integer :: i,n,it,irk,nt,iseed
  CHARACTER(6) charit

  rkc(1)=.5d0
  rkc(2)=.5d0
  rkc(3)=1.d0
  rkc(4)=1.d0/6.d0
  pi = 4.D0*DATAN(1.D0)
  do n=nmax,nmax,1
     dt = 0.1d0*(len / tmax / dble(n))
     nt = floor(1.d0 / dt)
     dt = tmax / dble(nt)

     call setupfft(n,len)
     ALLOCATE(u(n),x(n),up(n),fh(n),us(n))
     ALLOCATE(fu(n,3))

     ! Setup initial data
     DO i=1,n
        x(i) = dble(i-1)*len/n
        if (abs(x(i)-4.0d0) .lt. 1.d0) then 
           u(i) = 1.0d0
        else
           u(i) = 0.5d0
        end if
     END DO
     up = u
     !     now use the RK method
     do it=0,nt-1
        do irk=1,4
           if (irk.lt.2) then
              t = dble(it)*dt
           else
              t=(dble(it)+rkc(irk-1))*dt
           end if
           call rhs(fh,up,n,dt)
           up=u+rkc(irk)*fh
           if (irk.lt.4) then
              fu(:,irk)=fh
           else
              up=up+(fu(:,1)+2.d0*fu(:,2)+2.d0*fu(:,3))/6.d0
           end if
        end do
        !     now switch the time indices
        u=up
        ! filter
        call filter(u,filter_style,n)
        !
        if (mod(it,nmod).eq.0) then
           call spectrum(u,us,n)
           WRITE(charit,"(I6.6)") it
           open(unit=13,file="data" // charit  // ".txt",status="unknown")
           do i=1,n
              WRITE(13,FMT='(2E18.8)') x(i),u(i)!,us(i)
           end do
           close(13)
        end if
     end do

     call destroyfft
     DEALLOCATE(u,x,up,fu,fh)
  end do

contains
  subroutine rhs(fu,u,n,dt)
    IMPLICIT NONE
    ! work arrays
    integer :: n
    double precision  :: u(n),v(n)
    double precision  :: fu(n)
    double precision :: dt

    v = u*(1.d0-u)
    call deriv(v,fu,1,n)
    fu = -1.d0*dt*fu
    call deriv(u,v,2,n)
    fu = fu + 5.0d-3*dt*v
  end subroutine rhs
  
end program TRAFFIC

!! This module implements the generation of waves, using JONSWAP
!  spectrum and realistic weather conditions
!
MODULE WaveGeneration

  USE StateVariables
  USE PrecisionVar
  USE Solver
  USE iso_c_binding

  IMPLICIT NONE
  Include 'fftw3.f03'

  PRIVATE

 
  REAL(KIND=dp), PARAMETER :: tol = 1e-10 

  TYPE, PUBLIC :: TWaveGeneration
    ! inputs
    REAL(dp) :: length, width ! spatial domain dimensions
    INTEGER :: N,M ! grid 
    INTEGER  :: c1,c2 ! how many patches to extend such that the domain is covered 

    REAL(dp) :: windSpeed ! wind speed
    REAL(dp) :: fetch ! distance travelled by the wind
    REAL(dp) :: gammaVal ! characteristic of the wave
    ! calculated
    REAL(dp) :: omegaP ! frequency of the energy peak
    REAL(dp) :: thetaWind ! angle between the direction of the wind and the wave propagation dir

    REAL(dp) :: dkWaveX,dkWaveY ! grid for wavenumber mesh
    REAL(dp) :: dx,dy ! spatial grid -- uniform !!!!!!!!!!!!!!

    REAL(dp) :: dt, time

    REAL(dp), ALLOCATABLE :: x0(:), y0(:) ! domain grid
    REAL(dp), ALLOCATABLE :: kWave(:,:) ! wavenumber
    REAL(dp), ALLOCATABLE :: kWaveX(:), kWaveY(:) ! wavenumber vector components
    REAL(dp), ALLOCATABLE :: omega(:,:) ! frequency
    REAL(dp), ALLOCATABLE :: theta(:,:) ! direction of propagation

    REAL(dp), ALLOCATABLE :: freqSpectrum(:,:) ! JOSNWAP freq spectrum
    REAL(dp), ALLOCATABLE :: DirSpread(:,:) ! directional spreading function 
    REAL(dp), ALLOCATABLE :: EnergySpectrum(:,:) ! directional energy spectrum
    REAL(dp), ALLOCATABLE :: ampl(:,:) ! amplitude
    COMPLEX(dp), ALLOCATABLE :: complexAmpl(:,:,:) ! complex amplitude
    REAL(dp), ALLOCATABLE :: phase(:,:) ! complex amplitude
    COMPLEX(dp), ALLOCATABLE :: IFFT(:,:,:) ! complex amplitude
    REAL(dp), ALLOCATABLE :: Displ(:,:,:) ! displacement
    REAL(dp), ALLOCATABLE :: DisplWholeD(:,:,:) ! displacement

  CONTAINS
    !PROCEDURE, PASS(this), PUBLIC :: construct
    PROCEDURE, PASS(this), PUBLIC :: calcWaves
    PROCEDURE, PASS(this), PRIVATE :: writeDat
    PROCEDURE, PASS(this), PRIVATE :: genCharactWavesJONSWAP
    PROCEDURE, PASS(this), PRIVATE :: calcJONSWAPFreqSpectrum
    PROCEDURE, PASS(this), PRIVATE :: calcDirectSpreading
    !PROCEDURE, PASS(this), PRIVATE :: calcDispersionRel
    PROCEDURE, PASS(this), PRIVATE :: calcIFFT
  END TYPE TWaveGeneration    

  INTERFACE TWaveGeneration
    MODULE PROCEDURE construct
  END INTERFACE

CONTAINS

!--------------------------------------------------------------------------------------------------
  !! Constructor which allocates arrays,
  ! initialize parameters
  ! and domain characteristics
  !
  TYPE(TWaveGeneration) FUNCTION construct(N,M,c1,c2,U10,F,gammaVal,thetaW,length,width) RESULT( this )
  !SUBROUTINE construct(N,M,U10,F,gammaVal,thetaW,length,width, this) 
    !
    INTEGER, INTENT(in) :: N,M  ! size of the mesh
    INTEGER, INTENT(in) :: c1,c2  ! number of patches to fill the comp domain
    REAL(dp), INTENT(in) :: U10,F,gammaVal,thetaW,length,width
    !CLASS(TWaveGeneration), INTENT(out) :: this
    !

    ALLOCATE( this%x0(N), this%y0(M) )
    ALLOCATE( this%kWaveX(N), this%kWaveY(M) )
    ALLOCATE( this%kWave(N,M) )
    ALLOCATE( this%theta(N,M) )
    ALLOCATE( this%omega(N,M) )

    ALLOCATE( this%freqSpectrum(N,M) )
    ALLOCATE( this%DirSpread(N,M) )
    ALLOCATE( this%EnergySpectrum(N,M) )
    ALLOCATE( this%ampl(N,M) )
    ALLOCATE( this%Phase(N,M) )
    ALLOCATE( this%ComplexAmpl(1:3,N,M) )
    ALLOCATE( this%IFFT(1:3,N,M) )
    ALLOCATE( this%Displ(1:3,N,M) )
    ALLOCATE( this%DisplWholeD(1:3,N*c1,M*c2) )

    ! weather conditions
    this%windSpeed = U10
    this%fetch = F
    this%gammaVal = gammaVal
    this%thetaWind = thetaW

    ! domain characteristics
    this%length = length
    this%width = width
    this%N = N
    this%M = M
    this%c1 = c1
    this%c2 = c2


  END FUNCTION construct
  !END SUBROUTINE construct

!--------------------------------------------------------------------------------------------------
  !! Calculate the displacement of the water
  ! due to the waves
  SUBROUTINE calcWaves( this, time )
    !
    CLASS(TWaveGeneration), INTENT(inout) :: this
    REAL(dp), INTENT(in) :: time
    
    INTEGER :: i,j, end1, end2
    REAL(dp) :: r,c
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: intermVar 
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: intermVar2 

    ALLOCATE(intermVar(this%N, this%M))
    end1 = this%N*this%c1
    end2 = this%M*this%c2
    ALLOCATE(intermVar2(end1, end2))

    ! Generate the characteristics of the waves
    call genCharactWavesJONSWAP( this )

    ! calculate the complex temporal amplitudes of the waves
    do i = 1,this%N
      do j = 1,this%M
        r = dcos(this%omega(i,j)*Time+this%phase(i,j)) ! real part
        c = dsin(this%omega(i,j)*Time+this%phase(i,j)) ! imaginary part
        this%ComplexAmpl(1,i,j) = this%kWaveX(i)/(this%kWave(i,j)+tol)*this%ampl(i,j)*complex(r,c) ! for x
        this%ComplexAmpl(2,i,j) = this%kWaveY(j)/(this%kWave(i,j)+tol)*this%ampl(i,j)*complex(r,c) ! for y
        this%ComplexAmpl(3,i,j) = this%ampl(i,j)*complex(r,c) ! for z
      enddo
    enddo
    
    !! COMPUTE THE IFFT
    ! for x
    call this%CalcIFFT(1)
    do i =1, this%N
      do j=1, this%M
      this%Displ(1,i,j) = this%x0(i)+aimag(this%IFFT(1,i,j))
      enddo
    enddo

    ! for y
    call this%CalcIFFT(2)
    do i =1, this%N
      do j=1, this%M
      this%Displ(2,i,j) = this%y0(j)+aimag(this%IFFT(2,i,j))
      enddo
    enddo

    ! for z
    call this%CalcIFFT(3)
    do i =1, this%N
      do j=1, this%M
      this%Displ(3,i,j) = real(this%IFFT(3,i,j))
      enddo
    enddo

    !! extend the patch to larger domains
    !_x
    intermVar(:,:) = this%Displ(1,:,:)
    this%DisplWholeD(1,:,:) = extend(this%N, this%M, this%c1, this%c2, intermVar)
    ! readjust values 
    intermVar2(:,:) = this%DisplWholeD(1,:,:)
    this%DisplWholeD(1,:,:) = readjust(this%N, this%c1, this%M, this%c2, .True.,this%dx,intermVar2)

    !_y
    intermVar(:,:) = this%Displ(2,:,:)
    this%DisplWholeD(2,:,:) = extend(this%N, this%M, this%c1, this%c2, intermVar)
    ! readjust values 
    intermVar2(:,:) = this%DisplWholeD(2,:,:)
    this%DisplWholeD(2,:,:) = readjust(this%N, this%c1, this%M, this%c2, .False.,this%dy,intermVar2)

    !_z
    intermVar(:,:) = this%Displ(3,:,:)
    this%DisplWholeD(3,:,:) = extend(this%N, this%M, this%c1, this%c2, intermVar)


    ! check results
    call this%writeDat()

  END SUBROUTINE calcWaves

!--------------------------------------------------------------------------------------------------
  !! calculates the IFFT
  !
  SUBROUTINE calcIFFT( this, whichDir )

    CLASS(TWaveGeneration), INTENT(inout) :: this
    INTEGER, INTENT(in) :: whichDir

    INTEGER :: i,j
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: realPart, imagPart
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: inp,ou
    TYPE(C_PTR) :: plan_c2c
    !
    ALLOCATE(realPart(this%N,this%M))
    ALLOCATE(imagPart(this%N,this%M))
    ALLOCATE(inp(this%N,this%M))
    ALLOCATE(ou(this%N,this%M))

    inp(:,:) = this%ComplexAmpl(whichDir,:,:)
    ou(:,:) = this%IFFT(whichDir,:,:)

    ! define the plan
    plan_c2c = fftw_plan_dft_2d( this%N, this%M, inp, ou, FFTW_FORWARD,FFTW_ESTIMATE)

    ! determine real and imaginary part 
    RealPart = real(inp(:,:))
    ImagPart = aimag(inp(:,:))

    ! shift to have the terms for frequencies between 0 and n/2 and n/2 and n
    RealPart = shift2d(this%N, this%M, int(this%N/2), int(this%M/2), RealPart)
    ImagPart = shift2d(this%N, this%M, int(this%N/2), int(this%M/2), ImagPart)

    ! 
    do i=1,this%N
      do j=1,this%M
        inp(i,j) = complex(RealPart(i,j),ImagPart(i,j))
      enddo
    enddo

    ! calculate the ifft
    call fftw_execute_dft(plan_c2c, inp, ou)

    this%IFFT(whichDir,:,:)=ou(:,:)

    !! shift the result such that it is for the domain -n/2, n/2
    ! determine real and imaginary part
    RealPart = real(this%IFFT(whichDir,:,:))
    ImagPart = aimag(this%IFFT(whichDir,:,:))

    ! shift to have the terms for frequencies between 0 and n/2 and n/2 and n
    RealPart = shift2d(this%N, this%M, int(this%N/2), int(this%M/2), RealPart)
    ImagPart = shift2d(this%N, this%M, int(this%N/2), int(this%M/2), ImagPart)

    !
    do i=1,this%N
      do j=1,this%M
      this%IFFT(whichDir,i,j) = complex(RealPart(i,j),ImagPart(i,j))
      enddo
    enddo

    call fftw_destroy_plan(plan_c2c)

  END SUBROUTINE calcIFFT

!--------------------------------------------------------------------------------------------------
  !! Generate the characteristics of the waves
  !  using the JONSWAP spectrum. 
  !  The algo was proposed in J. Frechot: 
  !  " Realistic simulation of ocean surface usign wave spectra "
  !  in GRAPP 2006
  !
  SUBROUTINE genCharactWavesJONSWAP( this )
    !
    CLASS(TWaveGeneration), INTENT(inout) :: this

    INTEGER :: i,j

    !
    ! grid generation
    this%x0 = linspace(-this%N*0.5d0,this%N*0.5d0-1,this%N,.True.)*this%length/this%N
    this%y0 = linspace(-this%M*0.5d0,this%M*0.5d0-1,this%M,.True.)*this%width/this%M

    ! calculate the spatial grid
    this%dx = this%x0(2)-this%x0(1)
    this%dy = this%y0(2)-this%y0(1)
    print*, 'grid', this%dx,this%dy

    ! calculate the wavenumber component
    this%kWaveX = linspace(-this%N*0.5d0,this%N*0.5d0-1,this%N,.True.)*2.d0*pi/this%length 
    this%kWaveY = linspace(-this%M*0.5d0,this%M*0.5d0-1,this%M,.True.)*2.d0*pi/this%width

    ! calculate the wavenumber grid
    !! for the moment only uniform grid is considered
    this%dkWaveX = this%kWaveX(2)-this%kWaveX(1)
    this%dkWaveY = this%kWaveY(2)-this%kWaveY(1)

    do i = 1,this%N
      do j =1,this%M

        ! calculate the wavenumber norme
        this%kWave(i,j) = dsqrt(this%kWaveX(i)**2.d0+this%kWaveY(j)**2.d0)

        ! calculate the wave propagation direction
        this%theta(i,j) = atan2(this%kWaveY(j),this%kWaveX(i))

        ! calculate the frequency bandwidth through the dispersion relation
        this%omega(i,j) = calcDispersionRel(this%kWave(i,j))
      enddo
    enddo

    this%kWave = transpose(this%kWave)
    this%omega = transpose(this%omega)
    this%theta = transpose(this%theta)

    ! calculate then JOSNWAP frequency spectrum
    call this%calcJONSWAPFreqSpectrum() 

    ! calculate the directional spreading function
    call this%calcDirectSpreading()

    ! phase of the waves
    call random_seed()
    call random_number(this%phase)
    this%phase = this%phase*2.d0*pi

    do i = 1,this%N
      do j =1,this%M

        ! calculate the energy spectrum - and transform into wavenumber space
        this%EnergySpectrum(i,j) = this%freqSpectrum(i,j)*this%DirSpread(i,j)*0.5d0*sqrt(g/(this%kWave(i,j)+tol))/(this%kWave(i,j)+tol)

        ! calculate the amplitude
        this%ampl(i,j) = sqrt(2.d0*this%EnergySpectrum(i,j)*this%dkWaveX*this%dkWaveY)

      enddo
    enddo
    !
  END SUBROUTINE genCharactWavesJONSWAP

!--------------------------------------------------------------------------------------------------
  !! Calculate the JONSWAP frequency spectrum
  !
  SUBROUTINE calcJONSWAPFreqSpectrum( this)
    !
    CLASS(TWaveGeneration), INTENT(inout)  :: this
    !
    INTEGER(it4b) :: i,j
    REAL(dp) :: coef
    REAL(dp) :: alpha, sigma, r, sP, sPM

    ! Phillips constant
    alpha = 0.076d0*(this%windSpeed**2.d0/this%fetch/g)**0.22d0
    ! Frequency of the energy peak
    this%omegaP = 22.d0*(g/this%windSpeed)*(g*this%fetch/this%windSpeed**2.d0)**(-0.33d0)

    do i = 1, this%N
      do j= 1, this%M
        
        ! Width of the energy peak
        if (this%omega(i,j).le.this%omegaP) then
          sigma = 0.07d0 
        else
          sigma = 0.09d0
        endif
        ! Exposant constant
        coef = -(this%omega(i,j)-this%omegaP)**2.d0/(2.d0*sigma**2.d0*this%omegaP**2.d0)
        if (coef.le.-100d0) then
          r = 0.d0
        else
          r = dexp(coef)
        endif
        ! Phillips spectrum
        sP = alpha*g**2.d0/(this%omega(i,j)+tol)**5.d0
        ! Person-Moskowitz spectrum
        coef = -5.d0/4.d0*(this%omegaP/(this%omega(i,j)+tol))**4.d0
        if (coef.le.-100d0) then
          sPM = 0.d0
        else
          sPM = sP*dexp(coef)
        endif
        ! JONSWAP spectrum
        this%freqSpectrum(i,j) = sPM*this%gammaVal**r
      enddo  
    enddo
      !
  END SUBROUTINE calcJONSWAPFreqSpectrum

!--------------------------------------------------------------------------------------------------
  !! Calculate the directional spreading function
  !
  SUBROUTINE calcDirectSpreading( this )
    !
    CLASS(TWaveGeneration), INTENT(inout)  :: this
    !
    INTEGER(it4b) :: i,j
    REAL(dp) :: mu, sharpness, NormFactor
    REAL(dp) :: term

    do i = 1, this%N
      do j= 1, this%M
        ! Exposant constant
        if (this%omega(i,j).le.this%omegaP) then
          mu = 5.d0
        else
          mu = -2.5d0
        endif

        ! Sharpness of the directional spreading
        sharpness = 11.5d0*(g/this%omegaP/this%windSpeed)**2.5d0*(this%omega(i,j)/this%omegaP)**mu
        ! Normalization factor
        NormFactor = 0.5d0/sqrt(pi)*gamma(sharpness+1.d0)/gamma(sharpness+0.5d0) 
        ! Directional spreading
        term = dcos((this%theta(i,j)-this%thetaWind)/2.d0)
        if ( term .lt.1e-8) term = 0.d0
        this%DirSpread(i,j) = NormFactor*(term)**(2.d0*sharpness)
      enddo  
    enddo
    !
  END SUBROUTINE calcDirectSpreading

!--------------------------------------------------------------------------------------------------
  !! Calculate the wavenumber through the dispersion relation
  !
  FUNCTION calcDispersionRel( kWave ) RESULT( omega )
    !
    REAL(dp) :: kWave

    REAL(dp) :: omega
    
    omega = dsqrt(kWave*g)

  END FUNCTION calcDispersionRel 

!--------------------------------------------------------------------------------------------------
  !!  Returns the shifted matrix such that the 0 frequency
  !   is shifted at N/2 -- for N even number
  !
  FUNCTION shift2d(n1, n2, shift1, shift2, inp) result(ou)

    INTEGER,  INTENT(in) :: n1,n2,shift1,shift2
    REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(in) :: INP

    INTEGER :: i,j
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: OU

    allocate(ou(1:n1,1:n2))

    do i=1,shift1
      do j =1,shift2
        ou(i,j) = inp(i+shift1,j+shift1)
      enddo
    enddo
    do i=1,shift1
      do j =shift2+1,n2
        ou(i,j) = inp(i+shift1,j-shift1)
      enddo
    enddo
    do i=shift1+1,n1
      do j =1,shift2
        ou(i,j) = inp(i-shift1,j+shift1)
      enddo
    enddo
    do i=shift1+1,n1
      do j =shift2+1,n2
        ou(i,j) = inp(i-shift1,j-shift1)
      enddo
    enddo

  END FUNCTION shift2d

!--------------------------------------------------------------------------------------------------
  !!  readjust the values of the new grid 
  !   n1,c1, arr are related to the direction that needs to be readjust
  !   
  !
  FUNCTION readjust(n1,c1,n2,c2,is_x,delta,inp) result(ou)


    INTEGER, INTENT(in) :: n1,n2,c1,c2
    REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(in) :: inp
    REAL(dp), INTENT(in) :: delta
    LOGICAL, INTENT(in) :: is_x

    INTEGER :: i,j
    INTEGER :: end1, end2
    REAL(dp) :: endval
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ou

    end1 = n1*c1
    end2 = n2*c2

    ALLOCATE(ou(1:end1,1:end2))
    
    if (is_x.eqv..True.)then
      do i = 1,end1
        do j = 1, end2
          ou(i,j) = inp(i,j)+dabs(inp(1,j))*(1+2*int((i-1)/n1))
        enddo
      enddo 
      do i = 1,end1
        do j = 1, end2
          endval = ou(end1,j)+delta
          ou(i,j) = ou(i,j)-endval/2.d0
        enddo
      enddo
     else 
      do i = 1,end1
      do j = 1, end2
          ou(i,j) = inp(i,j)+dabs(inp(i,1))*(1+2*int((j-1)/n2))
        enddo
      enddo 
      do i = 1,end1
        do j = 1, end2
          endval = ou(i,end2)+delta
          ou(i,j) = ou(i,j)-endval/2.d0
        enddo
      enddo
    endif

  END FUNCTION readjust
!--------------------------------------------------------------------------------------------------
  !!  extend the patch to larger domains
  !   n1,n2 are the number of points on the patch
  !   c1, c2 is the number of patch to copy in each direction
  !
  FUNCTION extend(n1,n2,c1,c2,inp) result(ou)

    INTEGER, INTENT(in) :: n1,n2,c1,c2
    REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(in) :: inp

    INTEGER :: i,j,ii,jj
    INTEGER :: end1, end2
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ou

    end1 = n1*c1
    end2 = n2*c2

    ALLOCATE(ou(1:end1,1:end2))

      do i =1,end1
        ii = mod(i-1,n1)+1
        do j =1,end2
          jj=mod(j-1,n2)+1
          ou(i,j) = inp(ii,jj)       
        enddo
      enddo

  END FUNCTION extend

!--------------------------------------------------------------------------------------------------
  !!  Return evenly spaced numbers over a specified interval.
  !   Ported from the numpy routine.
  !   Author: Ivan Pribec
  !
  FUNCTION linspace(start,end,num,endpoint,step) result(samples)
    ! PARAMETERS
    real(dp), intent(in) :: start 
    !! The starting value of the sequence.
    real(dp), intent(in) :: end
    !! The end value of the sequence, unless `endpoint` is set to `.false.`. 
    !! In that case, the sequence consists of all but the last of `num + 1` 
    !! evenly spaced samples, so that `end` is excluded. Note that the 
    !! step size changes when `endpoint` is `.false.`.
    integer, intent(in), optional :: num
    !! Number of samples to generate. Default value is 50.
    logical, intent(in), optional :: endpoint
    !! If `.true.`, `end` is the last sample. Otherwise, it is not included. Default is `.true.`.
    real(dp), intent(out), optional :: step
    !! If present, `step` is the size of spacing between samples.
    ! RETURNS
    real(dp), allocatable :: samples(:)
    !! There are `num` equally spaced samples in the closed interval `[start, stop]` or 
    !! the half-open interval `[start, stop)` (depending on whether `endpoint` is `.true.` or `.false.`).

    integer :: num_, i
    logical :: endpoint_
    real(dp) :: step_

    num_ = 50
    if (present(num)) num_ = num

    endpoint_ = .true.
    if (present(endpoint)) endpoint_ = endpoint

    ! find step size
    if (endpoint_) then
      step_ = (end - start)/real(num_-1,dp)
    else
      step_ = (end - start)/real(num_,dp)
    end if

    if (present(step)) step = step_

    allocate(samples(num_))
    do i = 1, num_
       samples(i) = start + (i-1)*step_
    end do

  END FUNCTION linspace
!--------------------------------------------------------------------------------------------------
  !!  writes a dat file to check the results.
  !
  SUBROUTINE writeDat( this )

    CLASS(TWaveGeneration), INTENT(in) :: this

    INTEGER :: i,j
    INTEGER :: end1, end2

    end1 = this%N*this%c1
    end2 = this%M*this%c2

    open(unit=2, file='./Test_FFT.dat', status='replace')
      do i =1, this%N!end1
        do j=1, this%M!end2
          write(2, *) this%Displ(1,i,j), this%Displ(2,i,j), this%Displ(3,i,j)
        enddo
      enddo
    close(2)

    open(unit=3, file='./extTest_FFT.dat', status='replace')
      do i =1, end1
        do j=1,end2
          write(3, *) this%DisplWholeD(1,i,j), this%DisplWholeD(2,i,j), this%DisplWholeD(3,i,j)
        enddo
      enddo
    close(3)

  END SUBROUTINE writeDat

END MODULE 

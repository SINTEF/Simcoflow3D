!! This module initiates the level set and the vof fields
!  for the case under study
!
MODULE InitialVof

  USE StateVariables
  USE PrecisionVar
  USE Mesh
  USE CutCell
  USE WaveGeneration
  USE Clsvof

  IMPLICIT NONE

  PRIVATE

  REAL(dp), DIMENSION(:,:,:), POINTER :: vfl, phi
  REAL(dp), DIMENSION(:,:,:), POINTER :: nx, ny,nz
  REAL(dp), DIMENSION(:,:  ), POINTER :: etaZ

  REAL(dp), PARAMETER :: eta = 0.075d0

  PUBLIC :: InitialClsvofLiquidFieldWaves

  INTERFACE InitialClsvofLiquidFieldWaves
    MODULE PROCEDURE InitialClsvofLiquidFieldWaves
  END INTERFACE

  CONTAINS

  SUBROUTINE InitialClsvofLiquidFieldWaves(TGrid,TCell,TWave)

    IMPLICIT NONE

    TYPE(Grid), INTENT(in) :: TGrid
    TYPE(Cell), INTENT(inout), TARGET :: TCell
    TYPE(TWaveGeneration), INTENT(in), TARGET :: TWave

    INTEGER :: i,j,k
    REAL(dp) :: dis,vol,s
    REAL(dp) :: norm
    REAL(dp) :: tol

    tol = 1.d-20

    !NULLIFY(vfl)
    !if(associated(vfl).eqv..False.) then
    !  print*, 'Error using pointer for vfl'
    !else
      vfl => TCell%vofL
    !endif
    phi => TCell%phiL
    nx => TCell%nxL
    ny => TCell%nyL
    nz => TCell%nzL

    etaZ => Twave%DisplWholeD(3,:,:)

    ! wave

    do k =1,Kmax
      do j=1,Jmax
        do i=1,Imax
        ! for z between -lz/2 and lz/2
        ! level set
        phi(i,j,k) = TGrid%z(i,j,k)-etaZ(i,j)
        enddo
      enddo
    enddo
    do k =1,Kmax
      do j=1,Jmax
        do i=1,Imax
         
          ! normal vectors
          if (i>1.and.i<Imax.and.j>1.and.j<jmax.and.k>1.and.k<kmax) then
            call Normal_Vector_Irre(Tgrid,phi,i,j,k,nx(i,j,k),ny(i,j,k),nz(i,j,k))      
          else
            nx(i,j,k) = (phi(min(imax,i+1),j,k)-phi(max(1,i-1),j,k))/             &
                     (TGrid%x(min(imax,i+1),j,k)-TGrid%x(max(1,i-1),j,k))
            ny(i,j,k) =(phi(i,min(jmax,j+1),k)-phi(i,max(1,j-1),k))/             &
                     (TGrid%y(i,min(jmax,j+1),k)-TGrid%y(i,max(1,j-1),k))
            nz(i,j,k) =(phi(i,j,min(kmax,k+1))-phi(i,j,max(1,k-1)))/             &
                     (TGrid%z(i,j,min(kmax,k+1))-TGrid%z(i,j,max(1,k-1)))
          endif

          norm = dsqrt(nx(i,j,k)**2.d0+ny(i,j,k)**2.d0+nz(i,j,k)**2.d0)
          ! I do not know why exactly
          if (norm.lt.1d-14) then
            nx(i,j,k) = 0.d0
            ny(i,j,k) = 0.d0
            nz(i,j,k) = 1.d0
          else
            nx(i,j,k) = nx(i,j,k)/norm
            ny(i,j,k) = ny(i,j,k)/norm
            nz(i,j,k) = nz(i,j,k)/norm
          endif

          ! don't know why
          s = phi(i,j,k) + 0.5*(dabs(nx(i,j,k))*TGrid%dx(i,j,k)+                 &
                               dabs(ny(i,j,k))*TGrid%dy(i,j,k)+               &
                               dabs(nz(i,j,k))*TGrid%dz(i,j,k))
          ! compute volume of fraction
          call Volume_Fraction_Calc(TGrid%dx(i,j,k),TGrid%dy(i,j,k),         &
                 TGrid%dz(i,j,k),nx(i,j,k),ny(i,j,k),nz(i,j,k),s,vol)
          vfl(i,j,k) = 1.d0-vol/(TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k))

        enddo
      enddo
    enddo

    NULLIFY(vfl)
    NULLIFY(phi)
    NULLIFY(nx)
    NULLIFY(ny)
    NULLIFY(nz)

  END SUBROUTINE InitialClsvofLiquidFieldWaves

  !SUBROUTINE Normal_Vector_Irre(PGrid,phi,i,j,k,nxx,nyy,nzz)
  !  implicit none
  !  type(Grid),                        intent(in)  :: PGrid
  !  real(kind=dp),dimension(0:,0:,0:), intent(in)  :: phi
  !  integer,                           intent(in)  :: i,j,k
  !  real(dp),                          intent(out) :: nxx,nyy,nzz
  !  integer                                     :: case_deri,dx,dy,dz
  !  real(dp)                                    :: qi(-1:1),qj(-1:1),qk(-1:1)
  !  real(dp)                                    :: vx,vy,vz
  !     ! define qi for Dx
  !     if(i>=3) then
  !     ! for qi(-1)
  !       vx=(phi(i,j,k)-phi(i-2,j,k))/(PGrid%x(i,j,k)-PGrid%x(i-2,j,k))
  !     else
  !       vx=(phi(i,j,k)-phi(i-1,j,k))/(PGrid%x(i,j,k)-PGrid%x(i-1,j,k))
  !     end if
  !     vy=(phi(i-1,j+1,k)-phi(i-1,j-1,k))/(PGrid%y(i-1,j+1,k)-PGrid%y(i-1,j-1,k))
  !     vz=(phi(i-1,j,k+1)-phi(i-1,j,k-1))/(PGrid%z(i-1,j,k+1)-PGrid%z(i-1,j,k-1))
  !     qi(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
  !     ! for qi(1)
  !     if(i<=imax-2) then
  !       vx=(phi(i+2,j,k)-phi(i,j,k))/(PGrid%x(i+2,j,k)-PGrid%x(i,j,k))
  !     else
  !       vx=(phi(i+1,j,k)-phi(i,j,k))/(PGrid%x(i+1,j,k)-PGrid%x(i,j,k))
  !     end if
  !     vy=(phi(i+1,j+1,k)-phi(i+1,j-1,k))/(PGrid%y(i+1,j+1,k)-PGrid%y(i+1,j-1,k))
  !     vz=(phi(i+1,j,k+1)-phi(i+1,j,k-1))/(PGrid%z(i+1,j,k+1)-PGrid%z(i+1,j,k-1))
  !     qi(1) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
  !     if(qi(-1)<eta.and.qi(1)>=eta) then
  !       dx=-1
  !     elseif(qi(-1)>=eta.and.qi(1)<eta) then
  !       dx=1
  !     else
  !       dx=0
  !     end if
! !define qj for Dy
  !     if(j>=3) then
  !       vy=(phi(i,j,k)-phi(i,j-2,k))/(PGrid%y(i,j,k)-PGrid%y(i,j-2,k))
  !     else
  !       vy=(phi(i,j,k)-phi(i,j-1,k))/(PGrid%y(i,j,k)-PGrid%y(i,j-1,k))
  !     end if
  !     vx=(phi(i+1,j-1,k)-phi(i-1,j-1,k))/(PGrid%x(i+1,j-1,k)-PGrid%x(i-1,j-1,k))
  !     vz=(phi(i,j-1,k+1)-phi(i,j-1,k-1))/(PGrid%z(i,j-1,k+1)-PGrid%z(i,j-1,k-1))
  !     qj(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
  !     ! define qj(1)
  !     vx=(phi(i+1,j+1,k)-phi(i-1,j+1,k))/(PGrid%x(i+1,j+1,k)-PGrid%x(i-1,j+1,k))
  !     if(j<=jmax-2) then
  !       vy=(phi(i,j+2,k)-phi(i,j,k))/(PGrid%y(i,j+2,k)-PGrid%y(i,j,k))
  !     else
  !       vy=(phi(i,j+1,k)-phi(i,j,k))/(PGrid%y(i,j+1,k)-PGrid%y(i,j,k))
  !     end if
  !     vz=(phi(i,j+1,k+1)-phi(i,j+1,k-1))/(PGrid%z(i,j+1,k+1)-PGrid%z(i,j+1,k-1))
  !     qj(1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
  !     if(qj(-1)<eta.and.qj(1)>=eta) then
  !       dy=-1
  !     elseif(qj(-1)>=eta.and.qj(1)<eta) then
  !       dy=1
  !     else
  !       dy=0
  !     end if
! !define qk for dz
  !     ! qk(-1)
  !     vx=(phi(i+1,j,k-1)-phi(i-1,j,k-1))/(PGrid%x(i+1,j,k-1)-PGrid%x(i-1,j,k-1))
  !     vy=(phi(i,j+1,k-1)-phi(i,j-1,k-1))/(PGrid%y(i,j+1,k-1)-PGrid%y(i,j-1,k-1))
  !     if(k>=3) then
  !       vz=(phi(i,j,k)-phi(i,j,k-2))/(PGrid%z(i,j,k)-PGrid%z(i,j,k-2))
  !     else
  !       vz=(phi(i,j,k)-phi(i,j,k-1))/(PGrid%z(i,j,k)-PGrid%z(i,j,k-1))
  !     end if
  !     qk(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
  !     ! qk(1)
  !     vx=(phi(i+1,j,k+1)-phi(i-1,j,k+1))/(PGrid%x(i+1,j,k+1)-PGrid%x(i-1,j,k+1))
  !     vy=(phi(i,j+1,k+1)-phi(i,j-1,k+1))/(PGrid%y(i,j+1,k+1)-PGrid%y(i,j-1,k+1))
  !     if(k<=kmax-2) then
  !       vz=(phi(i,j,k+2)-phi(i,j,k))/(PGrid%z(i,j,k+2)-PGrid%z(i,j,k))
  !     else
  !       vz=(phi(i,j,k+1)-phi(i,j,k))/(PGrid%z(i,j,k+1)-PGrid%z(i,j,k))
  !     end if
  !     qk(1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
  !     if(qk(-1)<eta.and.qk(1)>=eta) then
  !       dz=-1
  !     elseif(qk(-1)>=eta.and.qk(1)<eta) then
  !       dz=1
  !     else
  !       dz=0
  !     end if
  !     if(dx==-1) nxx = (phi(i,j,k)-phi(i-1,j,k))/(PGrid%x(i,j,k)-PGrid%x(i-1,j,k))
  !     if(dx==1) nxx = (phi(i+1,j,k)-phi(i,j,k))/(PGrid%x(i+1,j,k)-PGrid%x(i,j,k))
  !     if(dx==0) nxx = (phi(i+1,j,k)-phi(i-1,j,k))/(PGrid%x(i+1,j,k)-PGrid%x(i-1,j,k))
  !     if(dy==-1) nyy = (phi(i,j,k)-phi(i,j-1,k))/(PGrid%y(i,j,k)-PGrid%y(i,j-1,k))
  !     if(dy==1) nyy = (phi(i,j+1,k)-phi(i,j,k))/(PGrid%y(i,j+1,k)-PGrid%y(i,j,k))
  !     if(dy==0) nyy = (phi(i,j+1,k)-phi(i,j-1,k))/(PGrid%y(i,j+1,k)-PGrid%y(i,j-1,k))
  !     if(dz==-1) nzz = (phi(i,j,k)-phi(i,j,k-1))/(PGrid%z(i,j,k)-PGrid%z(i,j,k-1))
  !     if(dz==1) nzz = (phi(i,j,k+1)-phi(i,j,k))/(PGrid%z(i,j,k+1)-PGrid%z(i,j,k))
  !     if(dz==0) nzz = (phi(i,j,k+1)-phi(i,j,k-1))/(PGrid%z(i,j,k+1)-PGrid%z(i,j,k-1))

  !  END SUBROUTINE ANormal_Vector_Irre

END MODULE InitialVof

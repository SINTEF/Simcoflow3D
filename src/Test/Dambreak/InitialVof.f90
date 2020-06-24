Module InitialVof
    USE ieee_arithmetic
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE Clsvof
    
    implicit none
    private
    
    real(kind=dp),dimension(:,:,:),pointer :: vfl,vflF    ! vfl represents the liquid volume fraction, vflF represents fluid volume fraction 
    real(kind=dp),dimension(:,:,:),pointer :: phi,phiF    ! phi represents the liquid level set function, phiF represents fluid level set function 
    real(kind=dp),dimension(:,:,:),pointer :: nxF,nyF,nzF ! nxF,nyF,nzF is fluid level set function
    
    public:: InitialClsvofFluidFieldDamBreak,InitialClsvofLiquidFieldDambreak,InitialGridDamBreak
    
    interface InitialClsvofFluidFieldDamBreak	
      module procedure InitialClsvofFluidFieldDamBreak
    end interface
    
    interface InitialClsvofLiquidFieldDamBreak
      module procedure InitialClsvofLiquidFieldDamBreak
    end interface
    
    interface InitialGridDamBreak
      module procedure InitialGridDamBreak
    end interface    
    
    contains
    
    subroutine InitialGridDamBreak(SPoint,EPoint,TGrid,Lref)
      implicit none
      type(Grid),intent(inout) :: TGrid
      type(Point),intent(in)   :: SPoint,EPoint
      real(kind=dp),intent(in) :: Lref
      integer		       :: i,j,k
        
      TGrid%Lref = Lref
      TGrid%dx(:,:,:) = (EPoint%x/Lref-SPoint%x/Lref)/dble(Imax)
      TGrid%dy(:,:,:) = (EPoint%y/Lref-SPoint%y/Lref)/dble(Jmax)
      TGrid%dz(:,:,:) = (EPoint%z/Lref-SPoint%z/Lref)/dble(Kmax)
      
      do i = 1,Imax
        do j = 1,Jmax
          do k = 1,Kmax
            TGrid%x(i,j,k) = SPoint%x/Lref+TGrid%dx(i,j,k)*(dble(i)-0.5d0)
            TGrid%y(i,j,k) = SPoint%y/Lref+TGrid%dy(i,j,k)*(dble(j)-0.5d0)
            TGrid%z(i,j,k) = SPoint%z/Lref+TGrid%dz(i,j,k)*(dble(k)-0.5d0)
          end do
        end do
      end do
    end subroutine InitialGridDamBreak
    
    subroutine InitialClsvofFluidFieldDambreak(TGrid,TCell)
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	      :: i,j,k
      real(kind=dp)		      :: dx,dy,dz,dis,vol,Radius,epsi,s
      real(kind=dp)		      :: tol
      
      tol = 1.d-20
      epsi = 1.d-40
      vflF => TCell%vof
      phiF => TCell%phi
      nxF => TCell%nx
      nyF => TCell%ny
      nzF => TCell%nz
      Radius = 0.5d0/TGrid%Lref
      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax
            dx=TGrid%x(i,j,k)-1000.d0
            dy=TGrid%y(i,j,k)-1000.d0
            dz=TGrid%z(i,j,k)-1000.d0
            dis=dsqrt(dx**2.d0+dy**2.d0+dz**2.d0)-Radius
            nxF(i,j,k)=dx/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            nyF(i,j,k)=dy/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            nzF(i,j,k)=dz/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            s=dis+0.5*(dabs(nxF(i,j,k))*TGrid%dx(i,j,k)+dabs(nyF(i,j,k))*      &
                       TGrid%dy(i,j,k)+dabs(nzF(i,j,k))*TGrid%dz(i,j,k))
            call Volume_Fraction_Calc(TGrid%dx(i,j,k),TGrid%dy(i,j,k),         &
                 TGrid%dz(i,j,k),nxF(i,j,k),nyF(i,j,k),nzF(i,j,k),s,vol)
            vflF(i,j,k)=vol/(TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k))
            phiF(i,j,k)=dis
          end do
        end do
      end do
      nullify(vflF)
      nullify(phiF)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
    end subroutine InitialClsvofFluidFieldDamBreak
    
    subroutine InitialClsvofLiquidFieldDamBreak(TGrid,TCell)
      implicit none
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	      :: i,j,k
      real(kind=dp)		      :: dx,dy,dz,dis,vol,epsi,s
      real(kind=dp)		      :: tol
      
      vfl => TCell%vofL
      phi => TCell%phiL
      nxF => TCell%nxL
      nyF => TCell%nyL
      nzF => TCell%nzL
    
    ! Dam break
      do i = 1,Imax
        do j = 1,Jmax
          do k = 1,Kmax
            dx = TGrid%x(i,j,k)-0.6d0/TGrid%Lref
            dz = TGrid%z(i,j,k)-0.3d0/TGrid%Lref
            if(dx<0.d0.and.dz<0.d0) then
              if(dabs(dx)<dabs(dz)) then
                phi(i,j,k) = dx
                nxF(i,j,k) = 1.d0
                nyF(i,j,k) = 0.d0
                nzF(i,j,k) = 0.d0
              else
                phi(i,j,k) = dz
                nxF(i,j,k) = 0.d0
                nyF(i,j,k) = 0.d0
                nzF(i,j,k) = 1.d0
              end if
            elseif(dx>=0.d0.and.dz<0.d0) then
              phi(i,j,k) = dx
              nxF(i,j,k) = 1.d0
              nyF(i,j,k) = 0.d0
              nzF(i,j,k) = 0.d0
            elseif(dx<0.d0.and.dz>=0.d0) then
              phi(i,j,k) = dz
              nxF(i,j,k) = 0.d0
              nyF(i,j,k) = 0.d0
              nzF(i,j,k) = 1.d0
            else
              phi(i,j,k) = dsqrt((dx**2.d0+dz**2.d0))
              nxF(i,j,k) = -dx/dsqrt(dx**2.d0+dz**2.d0)
              nyF(i,j,k) = 0.d0
              nzF(i,j,k) = -dz/dsqrt(dx**2.d0+dz**2.d0) 
            end if
            s=phi(i,j,k)+0.5*(dabs(nxF(i,j,k))*TGrid%dx(i,j,k)+	       		&
                               dabs(nyF(i,j,k))*TGrid%dy(i,j,k)+	       &
                               dabs(nzF(i,j,k))*TGrid%dz(i,j,k))
            call Volume_Fraction_Calc(TGrid%dx(i,j,k),TGrid%dy(i,j,k),         &
                 TGrid%dz(i,j,k),nxF(i,j,k),nyF(i,j,k),nzF(i,j,k),s,vol)
            vfl(i,j,k)=vol/(TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k))
          end do
        end do
      end do
      
      nullify(vfl)
      nullify(phi)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
    end subroutine InitialClsvofLiquidFieldDambreak
end module InitialVof

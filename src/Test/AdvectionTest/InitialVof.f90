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
    
    public:: InitialClsvofFluidFieldAdvectionTest,			       &
             InitialClsvofLiquidFieldAdvectionTest
    
    interface InitialClsvofFluidFieldAdvectionTest	
      module procedure InitialClsvofFluidFieldAdvectionTest
    end interface
    
    interface InitialClsvofLiquidFieldAdvectionTest
      module procedure InitialClsvofLiquidFieldAdvectionTest
    end interface
    
    contains
    
    
    subroutine InitialClsvofFluidFieldAdvectionTest(TGrid,TCell)
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	      :: i,j,k
      real(kind=dp)		      :: dx,dy,dz,dis,vol,Radius,epsi,s
      real(kind=dp)		      :: tol
      
      tol = 1.d-20
      epsi = 1.d-40
      if(associated(vflF).eqv..false.) then
        vflF => TCell%vof
      else
        print*, 'Error with using pointer for vflF InitialVof 67'
      end if
      if(associated(phiF).eqv..false.) then
        phiF => TCell%phi
      else
        print*, 'Error with using pointer for phiF InitialVof 72'
      end if
      if(associated(nxF).eqv..false.) then
        nxF => TCell%nx
      else
        print*, 'Error with using pointer for nxF Intial Vof 78'
      end if    
      if(associated(nyF).eqv..false.) then
        nyF => TCell%ny
      else
        print*, 'Error with using pointer for nyF InitialVof 83'
      end if
      if(associated(nzF).eqv..false.) then
        nzF => TCell%nz
      else
        print*, 'Error with using pointer for nzF InitialVof 88'
      end if  
      Radius = 0.5d0
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
    end subroutine InitialClsvofFluidFieldAdvectionTest
    
    subroutine InitialClsvofLiquidFieldAdvectionTest(TGrid,TCell)
      implicit none
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	      :: i,j,k
      real(kind=dp)		      :: dx,dy,dz,dis,vol,epsi,s
      real(kind=dp)		      :: tol,radius
      
      tol = 1.d-20
      epsi = 1.d-40 
      vfl => TCell%vofL
      phi => TCell%phiL
      nxF => TCell%nxL
      nyF => TCell%nyL
      nzF => TCell%nzL
    
    ! sphere     
      
      radius = 0.15d0
      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax
            dx=TGrid%x(i,j,k)-0.35d0
            dy=TGrid%y(i,j,k)-0.35d0
            dz=TGrid%z(i,j,k)-0.35d0
            phi(i,j,k)=dsqrt(dx**2.d0+dy**2.d0+dz**2.d0)-Radius
            nxF(i,j,k)=dx/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            nyF(i,j,k)=dy/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            nzF(i,j,k)=dz/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            s=phi(i,j,k)+0.5*(dabs(nxF(i,j,k))*TGrid%dx(i,j,k)+dabs(nyF(i,j,k))* &
                       TGrid%dy(i,j,k)+dabs(nzF(i,j,k))*TGrid%dz(i,j,k))
            call Volume_Fraction_Calc(TGrid%dx(i,j,k),TGrid%dy(i,j,k),           &
                 TGrid%dz(i,j,k),nxF(i,j,k),nyF(i,j,k),nzF(i,j,k),s,vol)
            vfl(i,j,k)=1.d0-vol/(TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k))
          end do
        end do
      end do
      nullify(vfl)
      nullify(phi)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
    end subroutine InitialClsvofLiquidFieldAdvectionTest
end module InitialVof

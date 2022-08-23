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
    
    public:: InitialClsvofLiquidFieldPersonal
    
    
    interface InitialClsvofLiquidFieldPersonal
      module procedure InitialClsvofLiquidFieldPersonal
    end interface
    
    contains
    
    subroutine InitialClsvofLiquidFieldPersonal(TGrid,TCell)
      implicit none
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	            :: i,j,k,ii,jj,kk
      real(kind=dp)		                :: dx,dy,dz,dis,vol,epsi,s
      real(kind=dp)		                :: tol    

      tol = 1.d-20
      epsi = 1.d-40 
      vfl => TCell%vofL
      phi => TCell%phiL
      nxF => TCell%nxL
      nyF => TCell%nyL
      nzF => TCell%nzL
    
      do i=0,Imax+1
        do j=0,Jmax+1
          do k=0,Kmax+1
            ii=max(1,min(Imax,i))
            jj=max(1,min(Jmax,j))
            kk=max(1,min(Kmax,k))
            dx=0.d0!TGrid%x(ii,jj,kk)-0.35d0
            dy=0.d0
            dz= TGrid%z(ii,jj,kk)
            phi(i,j,k)=dx+dy+dz
            nxF(i,j,k)=dx/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            nyF(i,j,k)=dy/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            nzF(i,j,k)=dz/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            s=phi(i,j,k)+0.5*(dabs(nxF(i,j,k))*TGrid%dx(ii,jj,kk)+                &
                              dabs(nyF(i,j,k))*TGrid%dy(ii,jj,kk)+                &
                              dabs(nzF(i,j,k))*TGrid%dz(ii,jj,kk))
            call Volume_Fraction_Calc(TGrid%dx(ii,jj,kk),TGrid%dy(ii,jj,kk),      &
                 TGrid%dz(ii,jj,kk),nxF(i,j,k),nyF(i,j,k),nzF(i,j,k),s,vol)
            vfl(i,j,k)=1.d0-vol/(TGrid%dx(ii,jj,kk)*TGrid%dy(ii,jj,kk)*TGrid%dz(ii,jj,kk))
          end do
        end do
      end do
      !
      nullify(vfl)
      nullify(phi)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
      !
    end subroutine InitialClsvofLiquidFieldPersonal
    !
end module InitialVof

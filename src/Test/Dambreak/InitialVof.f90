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
    
    public:: InitialClsvofFluidFieldDamBreak,				       &
             InitialClsvofLiquidFieldDambreak,				       &
             InitialVarDambreak,					       &
             InitialGridDamBreak
             
    interface InitialVarDambreak
      module procedure InitialVarDambreak
    end interface InitialVarDambreak
    
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
            vfl(i,j,k)=1.d0-vol/(TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k))
          end do
        end do
      end do
      
      nullify(vfl)
      nullify(phi)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
    end subroutine InitialClsvofLiquidFieldDambreak
    
    Subroutine InitialVarDambreak(Vari,Uint,Vint,Wint,Pint,Tint,Uref,Tref,Roref,Lref)
      Real(kind=dp),intent(in):: Uint,Vint,Wint,Pint,Tint,Uref,Tref,Roref,Lref
      Type(Variables),intent(inout):: Vari
      Integer(kind=it4b):: i,j,k
      Vari%Uint = Uint
      Vari%Vint = Vint
      Vari%Wint = Wint
      Vari%Pint = Pint
      Vari%Tint = Tint
      Vari%Uref = Uref
      Vari%Roref = Roref
      Vari%Pref = Roref*Uref**2.d0
      Vari%Tref = Tref
      nuref = nu
      Do i = 0,Imax+1
        Do j = 0,Jmax+1
          Do k = 0,Kmax+1
            Vari%u(i,j,k) = 0.d0 !Uint/Uref
            Vari%v(i,j,k) = 0.d0 !Vint/Uref
            Vari%w(i,j,k) = 0.d0
            Vari%p(i,j,k) = Pint/Vari%Pref
            Vari%Gpu(i,j,k) = 0.d0
            Vari%Gpv(i,j,k) = 0.d0
            Vari%t(i,j,k) = Tint/Tref
            Vari%mres(i,j,k) = 0.d0
            Vari%pres(i,j,k) = 0.d0
          End do
        End do
      End do
      Open(unit=5,file='Convergence.dat')
      close(5,status='delete')
      Rey = Uref*Lref/nuref
      Fr = Uref/dsqrt(g*Lref)
      Print*,"Reynolds number:",Rey
      Print*,"Froude number:",Fr
   !   Call BoundaryConditionVar(Vari)
    End subroutine InitialVarDambreak
end module InitialVof

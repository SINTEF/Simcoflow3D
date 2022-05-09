Module InitialVof
    USE ieee_arithmetic
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE Clsvof
    
    implicit none
    private
    
    real(kind=dp),dimension(:,:,:),pointer :: vfl,vflF    
    ! vfl represents the liquid volume fraction, vflF represents fluid volume fraction 
    real(kind=dp),dimension(:,:,:),pointer :: phi,phiF    
    ! phi represents the liquid level set function, phiF represents fluid level set function 
    real(kind=dp),dimension(:,:,:),pointer :: nxF,nyF,nzF 
    ! nxF,nyF,nzF is fluid level set function

    public:: InitialClsvofFluidFieldPersonal,                                  &
             InitialClsvofLiquidFieldPersonal,                                 &
             InitialGridPersonal,InitialVarPersonal
             
    interface InitialVarPersonal
      module procedure InitialVarPersonal
    end interface

    interface InitialClsvofFluidFieldPersonal
      module procedure InitialClsvofFluidFieldPersonal
    end interface
    
    interface InitialClsvofLiquidFieldPersonal
      module procedure InitialClsvofLiquidFieldPersonal
    end interface
    
    interface InitialGridPersonal
      module procedure InitialGridPersonal
    end interface    
    
    contains
    
    subroutine InitialGridPersonal(SPoint,EPoint,TGrid,Lref)
      !! The subroutine is used to compute the grid for dambreak test
      implicit none
      type(Grid),    intent(inout) :: TGrid
      !! The input grid
      type(Point),   intent(in)    :: SPoint,EPoint
      !! The position of computational domain corners
      real(kind=dp), intent(in)    :: Lref
      !! The reference length
      integer                      :: i,j,k
        
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
    end subroutine InitialGridPersonal
    
    subroutine InitialClsvofFluidFieldPersonal(TGrid,TCell)
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	      :: i,j,k
      real(kind=dp)		      :: dis
      
      nullify(vflF)
      nullify(phiF)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
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
      !
      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax
            !
            dis= -1d3
            !
            nxF(i,j,k)= 0.d0
            nyF(i,j,k)= 0.d0
            nzF(i,j,k)= 1.d0
            !
            vflF(i,j,k)= 1.d0
            !
            phiF(i,j,k)= dis
            !
          end do
        end do
      end do
      !
      nullify(vflF)
      nullify(phiF)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
      !
    end subroutine InitialClsvofFluidFieldPersonal
    !
    subroutine InitialClsvofLiquidFieldPersonal(TGrid,TCell)
      implicit none
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	      :: i, j, k
      real(kind=dp)		      :: dis
      
      vfl => TCell%vofL
      phi => TCell%phiL
      nxF => TCell%nxL
      nyF => TCell%nyL
      nzF => TCell%nzL
    
      do i = 1,Imax
        do j = 1,Jmax
          do k = 1,Kmax
            !
            dis= -1d3
            !
            nxF(i,j,k)= 0.d0
            nyF(i,j,k)= 0.d0
            nzF(i,j,k)= 1.d0
            !
            vfl(i,j,k)= 1.d0
            !
            phi(i,j,k)= dis
            !
          end do
        end do
      end do

      nullify(vfl)
      nullify(phi)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
    end subroutine InitialClsvofLiquidFieldPersonal
    !
    Subroutine InitialVarPersonal(Uint,Vint,Wint,Pint,Tint,Uref,Tref,Roref,Lref, Vari)
      !      
      Real(kind=dp),intent(in):: Uint,Vint,Wint,Pint,Tint,Uref,Tref,Roref,Lref
      Type(Variables),intent(inout):: Vari

      Integer(kind=it4b):: i,j,k
      !
      Vari%Uint = Uint
      Vari%Vint = Vint
      Vari%Wint = Wint
      Vari%Pint = Pint
      Vari%Tint = Tint
      Vari%Uref = Uref
      Vari%Roref = Roref
      Vari%Pref = Roref*Uref**2.d0
      Vari%Tref = Tref
      muref = mu
      !
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            !
            Vari%u(i,j,k) = Uint/Uref 
            Vari%v(i,j,k) = Vint/Uref
            Vari%w(i,j,k) = Wint/Uref
            Vari%p(i,j,k) = Pint/Vari%Pref
            Vari%Gpu(i,j,k) = 0.d0
            Vari%Gpv(i,j,k) = 0.d0
            Vari%Gpw(i,j,k) = 0.d0
            Vari%t(i,j,k) = Tint/Tref
            Vari%mres(i,j,k) = 0.d0
            Vari%pres(i,j,k) = 0.d0
            !
          End do
        End do
      End do
      !
      Rey = Uref*Lref*Roref/muref
      Fr = Uref/dsqrt(g*Lref)
      !
      Print*,"Reynolds number:",Rey
      Print*,"Froude number:",Fr
      !
    End subroutine InitialVarPersonal
   !
end module InitialVof

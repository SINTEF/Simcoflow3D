module Boundaryinterface
  use PrecisionVar

  implicit none
  private

  type :: BCBase
    integer(kind=it4b)                        :: flag(6)
    !
    real(kind=dp),dimension(:),allocatable   :: Const
    real(kind=dp),dimension(:,:),allocatable :: VarW,VarE,VarS,VarN,VarB,VarT 
    !
    procedure(SetDN),      pass(this),public,pointer :: SetDN=>null()
    procedure(SetConstant),pass(this),public,pointer :: SetConstant=>null()
    procedure(BCInterface),pass(this),public,pointer :: top=>null()
    procedure(BCInterface),pass(this),public,pointer :: bottom=>null() 
    procedure(BCInterface),pass(this),public,pointer :: north=>null()
    procedure(BCInterface),pass(this),public,pointer :: south=>null()
    procedure(BCInterface),pass(this),public,pointer :: east=>null()
    procedure(BCInterface),pass(this),public,pointer :: west=>null()
  end type BCBase
  
  abstract interface
     
     subroutine BCinterface(this,xin,yin,zin,dxin,dyin,dzin,pin,uin,vin,win,   &
                                             vofin,lvsin,time,nxin,nyin,nzin)
        import :: dp, BCBase
        class(BCBase),intent(inout)             :: this
        real(kind=dp),dimension(:,:),intent(in) :: xin,yin,zin,dxin,dyin,dzin
        real(kind=dp),dimension(:,:),intent(in) :: pin,uin,vin,win
        real(kind=dp),dimension(:,:),intent(in) :: vofin,lvsin
        real(kind=dp),               intent(in) :: time 
        real(kind=dp),dimension(:,:),intent(in),optional :: nxin,nyin,nzin
     
     end subroutine BCinterface
!
     subroutine SetDN(this,W,E,N,S,B,T)
       import it4b
       import BCBase
     !! Set the Dirichlet or Neumann boundary condition, 0 : Dirichlet, 1 : Neumann
       class(BCBase), intent(inout)  :: this
       integer(kind=it4b),intent(in) :: W,E,N,S,B,T
     
     end subroutine SetDN

     subroutine SetConstant(this, Constin)
        import :: dp, BCBase
        class(BCBase),                           intent(inout) :: this
        real(kind=dp), dimension(:), allocatable, intent(in)   :: Constin

     end subroutine SetConstant
  end interface

  public :: BCBase

end module Boundaryinterface

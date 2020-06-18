module BoundaryFunction
  use Boundaryinterface
  use PrecisionVar
  use Constants
  implicit none
  private

  interface BCBase
    module procedure construct
  end interface
  
  public :: BCBase
  
  public :: BCUW, BCUE, BCUN, BCUS, BCUB, BCUT, 					&
            BCVW, BCVE, BCVS, BCVN, BCVB, BCVT, 					&
            BCWW, BCWE, BCWS, BCWN, BCWB, BCWT, 					&
            BCPW, BCPE, BCPS, BCPN, BCPB, BCPT, 					&
            BCVofW, BCVofE, BCVofS, BCVofN, BCVofB, BCVofT, 				&
            BCLvsW, BCLvsE, BCLvsS, BCLvsN, BCLvsB, BCLvsT 		 
  
  contains
    
  type(BCBase) function construct(Isize, Jsize, Ksize) result(this)
     integer(kind=it4b), intent(in) :: Isize, Jsize, Ksize
     allocate(this%VarW(Jsize, Ksize))
     allocate(this%VarE(Jsize, Ksize))
     allocate(this%VarS(Isize, Ksize))
     allocate(this%VarN(Isize, Ksize))
     allocate(this%VarB(Isize, Jsize))
     allocate(this%VarT(Isize, Jsize))
     
     if(associated(this%SetConstant).eqv..false.) then
       this%SetConstant => SetConstant
     end if  
     if(associated(this%SetDN).eqv..false.) then
       this%SetDN => SetDN
     end if   
  end function
          
  subroutine SetDN(this, W, E, S, N, B, T)
  !! Set the Dirichlet or Neumann boundary condition, 0 : Dirichlet, 1 : Neumann    
    class(BCBase), intent(inout)   ::  this
    integer(kind=it4b), intent(in) :: W, E, S, N, B, T
    
    this%flag(1) = W
    this%flag(2) = E
    this%flag(3) = S
    this%flag(4) = N
    this%flag(5) = B
    this%flag(6) = T
  end subroutine SetDN
  
  subroutine SetConstant(this, Constin)
    class(BCBase), intent(inout)  		      :: this
    real(kind=dp), dimension(:), allocatable, intent(in) :: Constin	
    
    allocate(this%Const(sizeof(Constin))) 
    this%Const(:) = Constin(:) 
  end subroutine SetConstant
    
  subroutine BCUW(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,   &
  							 vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for u velocity
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarW(:,1))
    Array2D2 = sizeof(this%VarW(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarW, uin, dxin, this%flag(1), this%const(1), 		&
           					       Array2D1, Array2D2) 
    ! For user defined boundary condition.    
  end subroutine BCUW     
  
  subroutine BCUE(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for u velocity.
    class(BCBase), intent(inout)  	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarE(:,1))
    Array2D2 = sizeof(this%VarE(1,:))
    ! For simple boundary condition       
    call TypicalBC(this%VarE, uin, -dxin, this%flag(2), this%const(2), 		&
           					       Array2D1, Array2D2)
    ! For user defined boundary condition.   
  end subroutine BCUE
  
  subroutine BCUS(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for u velocity.
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarS(:,1))
    Array2D2 = sizeof(this%VarS(1,:))
    ! For simple boundary condition       
    call TypicalBC(this%VarS, uin, dyin, this%flag(3), this%const(3),  		&
           					       Array2D1, Array2D2)
    ! For user defined boundary condition.   
  end subroutine BCUS
  
  subroutine BCUN(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for u velocity.
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarN(:,1))
    Array2D2 = sizeof(this%VarN(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarN, uin, -dyin, this%flag(4), this%const(4),  	&
           					       Array2D1, Array2D2)
    ! For user defined boundary condition.   
  end subroutine BCUN 
  
  subroutine BCUB(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the bottom boundary for u velocity.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarB(:,1))
    Array2D2 = sizeof(this%VarB(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarB, uin, dzin, this%flag(5), this%const(5),  		&
           					       Array2D1, Array2D2)
  end subroutine BCUB
  
  subroutine BCUT(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the top boundary for u velocity.
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarT(:,1))
    Array2D2 = sizeof(this%VarT(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarT, uin, -dzin, this%flag(6), this%const(6),  	&
           					       Array2D1, Array2D2)   
  end subroutine BCUT
    
  subroutine BCVW(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for v velocity
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarW(:,1))
    Array2D2 = sizeof(this%VarW(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarW, vin, dxin, this%flag(1), this%const(1),  		&
           					       Array2D1, Array2D2) 
  end subroutine BCVW
  
  subroutine BCVE(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for v velocity
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarE(:,1))
    Array2D2 = sizeof(this%VarE(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarE, vin, -dxin, this%flag(2), this%const(2),  	&
           					        Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCVE 
  
  subroutine BCVS(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for v velocity
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarS(:,1))
    Array2D2 = sizeof(this%VarS(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarS, vin, dyin, this%flag(3), this%const(3), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCVS 
  
  subroutine BCVN(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for v velocity
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
  
    Array2D1 = sizeof(this%VarN(:,1))
    Array2D2 = sizeof(this%VarN(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarN, vin, -dyin, this%flag(4), this%const(4), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCVN
  
  subroutine BCVB(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
    !! Compute the boundary value at the bottom boundary for v velocity.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarB(:,1))
    Array2D2 = sizeof(this%VarB(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarB, vin, dzin, this%flag(5), this%const(5),  		&
           					       Array2D1, Array2D2)
  end subroutine BCVB  							
  
  subroutine BCVT(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
    !! Compute the boundary value at the tope boundary for v velocity.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarT(:,1))
    Array2D2 = sizeof(this%VarT(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarT, vin, -dzin, this%flag(6), this%const(6),  	&
           					       Array2D1, Array2D2) 
  end subroutine BCVT  
  
  subroutine BCWW(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for w velocity
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarW(:,1))
    Array2D2 = sizeof(this%VarW(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarW, win, dxin, this%flag(1), this%const(1),  		&
           					       Array2D1, Array2D2) 
  end subroutine BCWW
  
  subroutine BCWE(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for w velocity
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarE(:,1))
    Array2D2 = sizeof(this%VarE(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarE, win, -dxin, this%flag(2), this%const(2),  	&
           					        Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCWE 
  
  subroutine BCWS(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for w velocity
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarS(:,1))
    Array2D2 = sizeof(this%VarS(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarS, win, dyin, this%flag(3), this%const(3), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCWS 
  
  subroutine BCWN(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for w velocity
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
  
    Array2D1 = sizeof(this%VarN(:,1))
    Array2D2 = sizeof(this%VarN(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarN, win, -dyin, this%flag(4), this%const(4), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCWN
  
  subroutine BCWB(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
    !! Compute the boundary value at the bottom boundary for w velocity.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarB(:,1))
    Array2D2 = sizeof(this%VarB(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarB, win, dzin, this%flag(5), this%const(5),  		&
           					       Array2D1, Array2D2)
  end subroutine BCWB  							
  
  subroutine BCWT(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
    !! Compute the boundary value at the tope boundary for w velocity.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarT(:,1))
    Array2D2 = sizeof(this%VarT(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarT, win, -dzin, this%flag(6), this%const(6),  	&
           					       Array2D1, Array2D2) 
  end subroutine BCWT 
    
  subroutine BCPW(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for pressure
    class(BCBase), intent(inout)  	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarW(:,1))
    Array2D2 = sizeof(this%VarW(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarW, pin, dxin, this%flag(1), this%const(1), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCPW
  
  subroutine BCPE(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for pressure
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarE(:,1))
    Array2D2 = sizeof(this%VarE(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarE, pin, -dxin, this%flag(2), this%const(2), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCPE

  subroutine BCPS(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for pressure
    class(BCBase), intent(inout)              :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarS(:,1))
    Array2D2 = sizeof(this%VarS(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarS, pin, dyin, this%flag(3), this%const(3), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCPS

  subroutine BCPN(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for pressure
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarN(:,1))
    Array2D2 = sizeof(this%VarN(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarN, pin, -dyin, this%flag(4), this%const(4), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCPN
  
  subroutine BCPB(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
    !! Compute the boundary value at the bottom boundary for pressure.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarB(:,1))
    Array2D2 = sizeof(this%VarB(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarB, pin, dzin, this%flag(5), this%const(5),  		&
           					       Array2D1, Array2D2)
  end subroutine BCPB  							
  
  subroutine BCPT(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
    !! Compute the boundary value at the tope boundary for pressure .							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarT(:,1))
    Array2D2 = sizeof(this%VarT(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarT, pin, -dzin, this%flag(6), this%const(6),  	&
           					       Array2D1, Array2D2)
  end subroutine BCPT
  
  subroutine BCVofW(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,  &
  							vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for volume of fluid
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarW(:,1))
    Array2D2 = sizeof(this%VarW(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarW, vofin, dxin, this%flag(1), this%const(1), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCVofW
  
  subroutine BCVofE(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for volume of fluid
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarE(:,1))
    Array2D2 = sizeof(this%VarE(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarE, vofin, -dxin, this%flag(2), this%const(2), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCVofE

  subroutine BCVofS(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for volume of fluid
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarS(:,1))
    Array2D2 = sizeof(this%VarS(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarS, vofin, dyin, this%flag(3), this%const(3), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCVofS
  
  subroutine BCVofN(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for volume of fluid
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarN(:,1))
    Array2D2 = sizeof(this%VarN(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarN, vofin, -dyin, this%flag(4), this%const(4), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCVofN
  
  subroutine BCVofB(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win, &
  							vofin, lvsin, time)
    !! Compute the boundary value at the bottom boundary for volume of fluid.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarB(:,1))
    Array2D2 = sizeof(this%VarB(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarB, vofin, dzin, this%flag(5), this%const(5),        &	  
           					       Array2D1, Array2D2)
  end subroutine BCVofB  							
  
  subroutine BCVofT(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win, &
  							vofin, lvsin, time)
    !! Compute the boundary value at the tope boundary for volume of fluid .							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarT(:,1))
    Array2D2 = sizeof(this%VarT(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarT, vofin, -dzin, this%flag(6), this%const(6),  	&
           					       Array2D1, Array2D2)
  end subroutine BCVofT
  
  subroutine BCLvsW(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,  &
  							vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for level set function
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
  
    Array2D1 = sizeof(this%VarW(:,1))
    Array2D2 = sizeof(this%VarW(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarW, lvsin, dxin, this%flag(1), this%const(1), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCLvsW
  
  subroutine BCLvsE(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for level set function 
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarE(:,1))
    Array2D2 = sizeof(this%VarE(1,:)) 
    ! For simple boundary condition.   
    call TypicalBC(this%VarE, lvsin, -dxin, this%flag(2), this%const(2), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCLvsE

  subroutine BCLvsS(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for level set function 
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
    
    Array2D1 = sizeof(this%VarS(:,1))
    Array2D2 = sizeof(this%VarS(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarS, lvsin, dyin, this%flag(3), this%const(3), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCLvsS
  
  subroutine BCLvsN(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win,    &
  							vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for level set function 
    class(BCBase), intent(inout) 		   	 :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)	   	      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2
   
    Array2D1 = sizeof(this%VarN(:,1))
    Array2D2 = sizeof(this%VarN(1,:))
    ! For simple boundary condition   
    call TypicalBC(this%VarN, lvsin, -dyin, this%flag(4), this%const(4), 		&
    							Array2D1, Array2D2) 
    ! For user defined boundary condition.   
  end subroutine BCLvsN
  
  subroutine BCLvsB(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win, &
  							vofin, lvsin, time)
    !! Compute the boundary value at the bottom boundary for level set function.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarB(:,1))
    Array2D2 = sizeof(this%VarB(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarB, lvsin, dzin, this%flag(5), this%const(5),        &	  
           					       Array2D1, Array2D2)
  end subroutine BCLvsB  							
  
  subroutine BCLvsT(this, xin, yin, zin, dxin, dyin, dzin, pin, uin, vin, win, &
  							vofin, lvsin, time)
    !! Compute the boundary value at the tope boundary for level set function.							
    class(BCBase), intent(inout) 	      :: this
    real(kind=dp), dimension(:,:), intent(in) :: xin, yin, zin, dxin, dyin, dzin 
    real(kind=dp), dimension(:,:), intent(in) :: pin, uin, vin, win
    real(kind=dp), dimension(:,:), intent(in) :: vofin, lvsin
    real(kind=dp), intent(in)		      :: time 
    integer(kind=it4b)			      :: Array2D1, Array2D2

    Array2D1 = sizeof(this%VarT(:,1))
    Array2D2 = sizeof(this%VarT(1,:))
    ! For simple boundary condition
    call TypicalBC(this%VarT, lvsin, -dzin, this%flag(6), this%const(6),  	&
           					       Array2D1, Array2D2)
  end subroutine BCLvsT

  subroutine TypicalBC(Array, Varin, dxyz, flag, const, Array2D1, Array2D2)
    real(kind=dp), dimension(:,:), ALLOCATABLE, intent(inout) :: Array
    real(kind=dp), dimension(:,:), intent(in)    	      :: Varin, dxyz
    integer(kind=it4b), intent(in)    	                      :: flag
    integer(kind=it4b), intent(in)			      :: Array2D1, Array2D2
    real(kind=dp), intent(in)				      :: const
    integer(kind=it4b)	   	  			      :: i,j
     
    if(flag == 0) then ! Dirichlet BC
      do i = 1, Array2D1
        do j = 1, Array2D2
          Array(i,j) = const
        end do
      end do
    else ! Neumann BC
      do i = 1, Array2D1
        do j = 1, Array2D2
          Array(i,j) = Varin(i,j) - const*dxyz(i,j)/2.d0
        end do
      end do 
    end if
  end subroutine TypicalBC   
end module BoundaryFunction

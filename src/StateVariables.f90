Module StateVariables
    USE PrecisionVar
    USE Mesh
    Implicit none
    Integer(kind=it4b):: ight=1,jght=1,kght=1
    Integer(kind=it4b):: ite
    Real(kind=dp),parameter:: pi = 4.d0*datan(1.d0),Cp = 1.005d3
    Real(kind=dp),parameter:: nu = 1.0034d-6,kT = 0.0271d0,kTw = 0.0271d0
    Real(kind=dp),parameter:: nuw=1.0034d-6,nua=1.506d-5,                &
                              roa=1.225d0,row=998.2d0 !
    Real(kind=dp),public :: nuref
    Real(kind=dp):: Rey,wa,Ta,xc,yc,zc,R1 = 0.5d0, R2 = 4.d0

    Type :: Variables
      Real(kind=dp),dimension(:,:,:),allocatable:: u,v,w,p,t,Gpu,Gpv,Gpw,ures, &
                                                   vres,wres,pres,mres
      Real(kind=dp):: Uint,Vint,Wint,Pint,Tint,Uref,Roref,Pref,Tref
    End Type
    
    Public:: InitialVar,BoundaryConditionVar
    
    Interface InitialVar
      Module procedure InitialVar
    End interface InitialVar
    
    Interface BoundaryConditionVar
      Module procedure BoundaryConditionVar
    End interface BoundaryConditionVar
    
    Contains
    
    Subroutine InitialVar(Vari,Uint,Vint,Wint,Pint,Tint,Uref,Tref,Roref,Lref)
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
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
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
      Print*, Rey
      Call BoundaryConditionVar(Vari)
    End subroutine InitialVar
    !*******************************************************
    !             dp/dn = 0.d0, u,v = 0.d0
    !  wall________________________
    !
    !
    !
    ! => inlet: u = 5.d-3          => outflow:: du/dn = 0.d0
    !         : v = 0.d0                        dv/dn = 0.d0
    !         : dp/dn = 0.d0                    P = 0.d0
    !
    !  wall________________________
    !             dp/dn = 0.d0, u,v = 0.d0
    !*******************************************************
    
    subroutine BoundaryConditionVarNew(Vari)
      type(Grid), intent(in)  	      :: PGrid
      type(Cell), intent(in)  	      :: PCell
      type(TVariables), intent(inout) :: Vari
      type(BCBase), intent(inout)     :: BCu, BCv, BCw, BCp  	
      real(kind=dp), intent(in)       :: Time
      integer(kind=it4b)              :: i,j
    end subroutine BoundaryConditionVarNew
    
    Subroutine BoundaryConditionVar(Vari)
      Type(Variables),intent(inout):: Vari
      Integer(kind=it4b):: i,j,k
      Real(kind=dp),parameter:: Twall = 300.d0
      Do j = 1,Jmax
        Do k = 1,Kmax
      ! Inlet boundary
          Vari%p(1-ight,j,k) = Vari%p(1,j,k)        ! Vari%Pint/(Vari%Pref)
          Vari%t(1-ight,j,k) = Vari%t(1,j,k)
          Vari%u(1-ight,j,k) = Vari%Uint/Vari%Uref  ! -Vari%u(1,j,k)
          Vari%v(1-ight,j,k) = 0.d0 - Vari%v(1,j,k) ! Vari%Vint/Vari%Uref
          Vari%w(1-ight,j,k) = 0.d0 - Vari%w(1,j,k) !
      ! Outlet boundary
          Vari%p(Imax+ight,j,k) = 0.d0 - Vari%p(Imax,j,k)
          Vari%t(Imax+ight,j,k) = Vari%t(Imax,j,k)
          Vari%u(Imax+ight,j,k) = Vari%u(Imax,j,k)  ! Vari%Uint/Vari%Uref
          Vari%v(Imax+ight,j,k) = Vari%v(Imax,j,k)  ! Vari%Vint/Vari%Uref
          Vari%w(Imax+ight,j,k) = Vari%w(Imax,j,k)
        End do
      End do
      Do i = 1,Imax
        Do k = 1,Kmax
      ! slip wall boundary
          Vari%p(i,1-jght,k) = Vari%p(i,1,k)
          Vari%t(i,1-jght,k) = 2.d0*Twall/Vari%Tref-Vari%t(i,1,k)
          Vari%u(i,1-jght,k) = Vari%u(i,1,k)
          Vari%v(i,1-jght,k) = 0.d0 !- Vari%v(i,jbeg+jght) !VariVari%Vint/Vari%Uref
          Vari%w(i,1-jght,k) = Vari%w(i,1,k)
      ! slip wall boundary
          Vari%p(i,Jmax+jght,k) = Vari%p(i,Jmax,k)
          Vari%t(i,Jmax+jght,k) = 2.d0*Twall/Vari%Tref-Vari%t(i,Jmax,k)
          Vari%u(i,Jmax+jght,k) = Vari%u(i,Jmax,k)
          Vari%v(i,Jmax+jght,k) = 0.d0-Vari%v(i,Jmax-jght,k)
          Vari%v(i,Jmax,k) = 0.d0 !Vari%Uint/Vari%Uref
          Vari%w(i,Jmax+jght,k) = Vari%w(i,Jmax,k)
        End do
      End do
      Do i = 1,Imax
        Do j = 1,Jmax
      ! Slip wall boundary condition
          Vari%p(i,j,1-kght) = Vari%p(i,j,1)
          Vari%t(i,j,1-kght) = Vari%t(i,j,1)
          Vari%u(i,j,1-kght) = Vari%u(i,j,1)
          Vari%v(i,j,1-kght) = Vari%v(i,j,1)
          Vari%w(i,j,1-kght) = 0.d0
      ! Slip wall boundary condition
          Vari%p(i,j,Kmax+kght) = Vari%p(i,j,Kmax)
          Vari%t(i,j,Kmax+kght) = Vari%t(i,j,Kmax)
          Vari%u(i,j,Kmax+kght) = Vari%u(i,j,Kmax)
          Vari%v(i,j,Kmax+kght) = Vari%v(i,j,Kmax)
          Vari%w(i,j,Kmax+kght) = 0.d0-Vari%w(i,j,Kmax-1)
          Vari%w(i,j,Kmax) = 0.d0
        End do
      End do
    End subroutine BoundaryConditionVar
End Module StateVariables




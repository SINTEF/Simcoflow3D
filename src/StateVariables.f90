Module StateVariables
    USE PrecisionVar
    USE Mesh
    Implicit none
    Integer(kind=it4b):: ight=1,jght=1,kght=1
    Integer(kind=it4b):: ite
    Real(kind=dp),parameter   :: pi = 4.d0*datan(1.d0),Cp = 1.005d3
    Real(kind=dp),parameter   :: nu = 1.002d-3,kT = 0.0271d0,kTw = 0.0271d0
    Real(kind=dp),parameter   :: nuw=1.002d-3,nua=1.506d-5,                     &
                                 roa=1.205d0,row=998.3d0,g=9.80665d0 
    Real(kind=dp),public      :: nuref,gx,gy,gz
    Real(kind=dp) 	          :: Fr,Rey,wa,Ta,xc,yc,zc,R1=0.5d0,R2=4.d0
    Integer(kind=it4b),public :: TimeOrder,SpaceOrder

    Type :: Variables
      Real(kind=dp),dimension(:,:,:),allocatable:: u,v,w,p,t,Gpu,Gpv,Gpw,      &
                                                   ures,vres,wres,pres,mres
      Real(kind=dp):: Uint,Vint,Wint,Pint,Tint,Uref,Roref,Pref,Tref
    End Type
    
    Public:: InitialVar
    
    Interface InitialVar
      Module procedure InitialVar
    End interface InitialVar
    
    Contains
    !    
    Subroutine InitialVar(Uint,Vint,Wint,Pint,Tint,Uref,Tref,Roref,Lref, Vari)
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
      !
      nuref = nu
      !
      Do i = 0,Imax+1
        Do j = 0,Jmax+1
          Do k = 0,Kmax+1
            !
            Vari%u(i,j,k) = 0.d0 !Uint/Uref
            Vari%v(i,j,k) = 0.d0 !Vint/Uref
            Vari%w(i,j,k) = 0.d0
            Vari%p(i,j,k) = Pint/Vari%Pref
            Vari%Gpu(i,j,k) = 0.d0
            Vari%Gpv(i,j,k) = 0.d0
            Vari%t(i,j,k) = Tint/Tref
            Vari%mres(i,j,k) = 0.d0
            Vari%pres(i,j,k) = 0.d0
            !
          End do
        End do
      End do
      !
      Open(unit=5,file='Convergence.dat')
      close(5,status='delete')
      !
      Rey = Uref*Lref*Roref/nuref
      Fr = Uref/dsqrt(g*Lref)
      !
      Print*,"Reynolds number:",Rey,Uref,Lref,Roref,nuref
      Print*,"Froude number:",Fr
      !
    End subroutine InitialVar
    !
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
End Module StateVariables




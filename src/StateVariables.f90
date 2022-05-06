Module StateVariables
    USE PrecisionVar
    USE Mesh
    Implicit none
    Integer(kind=it4b):: ight=1,jght=1,kght=1
    Integer(kind=it4b):: ite
    Real(kind=dp),parameter   :: pi = 4.d0*datan(1.d0),Cp = 1.005d3
    Real(kind=dp),parameter   :: mu = 1.002d-3,kT = 0.0271d0,kTw = 0.0271d0
    Real(kind=dp),parameter   :: muw=1.002d-3,mua=1.506d-5,                     &
                                 roa=1.205d0,row=998.3d0,g=9.80665d0 
    Real(kind=dp),public      :: muref,gx,gy,gz
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
      muref = mu
      !
      !Do i = 0,Imax+1
      !  Do j = 0,Jmax+1
      !    Do k = 0,Kmax+1
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            !
            Vari%u(i,j,k) = Uint/Uref
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
      Rey = Uref*Lref*Roref/muref
      Fr = Uref/dsqrt(g*Lref)
      !
      Print*,"Reynolds number:",Rey,Uref,Lref,Roref,muref
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
    !
    subroutine get_max_min(var, maxvalue,minvalue,idxMaxvalue,idxMinvalue)
      !
      implicit none
      !
      real(dp), dimension(0:,0:,0:), intent(in) :: var
      real(dp), intent(out) :: maxvalue,minvalue
      integer(it4b), dimension(3),intent(out) :: idxMaxvalue,idxMinvalue

      integer(it4b) :: i,j,k

      maxvalue=10e-12
      minvalue=10e+12
      !
      i=-1
      j=-1
      k=-1
      !
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
          !print*, 'inside',i,j,k,var(i,j,k)
          if(var(i,j,k).gt.maxvalue)then
              maxvalue=var(i,j,k)
              idxMaxvalue(1)=i 
              idxMaxvalue(2)=j 
              idxMaxvalue(3)=k 
           endif
          if(var(i,j,k).lt.minvalue)then
              minvalue=var(i,j,k)
              idxMinvalue(1)=i 
              idxMinvalue(2)=j 
              idxMinvalue(3)=k 
           endif
      enddo
      enddo
      enddo
      !
      end subroutine get_max_min
      !
End Module StateVariables




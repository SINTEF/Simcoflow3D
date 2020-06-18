Module Solver
    USE PrecisionVar
    USE Mesh
    USE StateVariables
    USE CutCell
    USE Clsvof
    USE PrintResult
    USE ComputePUV
    USE MPI
    USE BoundaryInterface
    USE BoundaryFunction
    
    Implicit none
    Private
    Type,Public:: SolverTime
        Integer(kind=it8b):: iter
        Real(kind=dp):: cfl
        Real(kind=dp):: dt,PhysT,NondiT
    End Type SolverTime
    Type,Public:: SolverConvergence
        Real(kind=dp):: N1,N2,NInf,N1c,N2c,NInfc
    End Type SolverConvergence
    Real(kind=dp),parameter:: Twalls = 300.d0
    Real(kind=dp),parameter:: Twall = 400.d0
    Real(kind=dp):: Tref,Prn
    Real(dp),dimension(:,:,:),pointer:: Tem,u,v,w
    
    public:: IterationSolution
    
    interface IterationSolution
      module Procedure IterationSolution
    end Interface IterationSolution
    
    contains
    
    subroutine IterationSolution(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,    &
                                 WCell,BCu,BCv,BCw,BCp,BCVof,BCLvs,TVar,iprint)
        Implicit none
        Type(Grid),intent(in)         		   :: PGrid,UGrid,VGrid,WGrid
        Type(Cell),intent(inout)      		   :: PCell,UCell,VCell,WCell
        Type(Variables),intent(inout) 		   :: TVar
        type(BCBase),intent(inout)		   :: BCu,BCv,BCw,BCp,BCVof,BCLvs
        Integer(kind=it4b),intent(in) 		   :: iprint
        Real(kind=dp),dimension(:,:,:),allocatable :: GraP
        Type(SolverTime)			   :: Time
        Type(SolverConvergence)			   :: UConv,VConv,WConv,PConv
        Type(Variables)				   :: TVar_n
        Integer(kind=it8b)			   :: itt
        Allocate(TVar_n%p(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
        Allocate(TVar_n%u(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
        Allocate(TVar_n%v(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
        Allocate(TVar_n%w(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
        Allocate(GraP(Imax,Jmax,Kmax))
        
        Time%iter = 10**6
        Time%NondiT = 0.d0
        Time%Cfl = 0.5d0
        
        do itt = 1,Time%iter
          call AdamBasforthBDF2(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,     &
                     WCell,BCu,BCv,BCw,BCp,BCVof,BCLvs,TVar,TVar_n,UConv,      &
                     VConv,WConv,PConv,Time,itt)
          Time%NondiT = Time%NondiT+Time%dt
          Time%PhysT = Time%Nondit*PGrid%Lref/TVar%URef
          
          call PrintHistory(itt,Uconv)
          call PrintDragLiftCoef(TVar,PGrid,UGrid,VGrid,WGrid,PCell,UCell,     &
                                 VCell,WCell,itt,Time%NondiT)
          ite = itt
          ! print*, itt
          if(mod(itt,iprint)==0)then
            write(*,*), itt,Time%PhysT,Time%NondiT
            call PrintResultVTK(PGrid,TVar,PCell,itt)
            call PrintResultVTR3D(PGrid,TVar,PCell,itt)
          ! call PrintResultTecplotPCent(PGrid,TVar,PCell,itt)
          ! call PrintResultTecplotPCentXY(PGrid,TVar,PCell,itt)
          ! call PrintResultTecplotUCent(UGrid,TVar,UCell,itt)
          ! call PrintResultTecplotVCent(VGrid,TVar,VCell,itt)
          ! call PrintResultTecplotWCent(WGrid,TVar,WCell,itt)
          end if
        end do
        deallocate(GraP,TVar_n%p,TVar_n%u,Tvar_n%v)
    end subroutine IterationSolution

    Subroutine AdamBasforthBDF2(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,     &
                     WCell,BCu,BCv,BCw,BCp,BCVof,BCLvs,			       &
                     TVar,TVar_n,UConv,VConv,WConv,PConv,Time,itt)
        Implicit none
        Type(Grid),intent(in)               :: PGrid,UGrid,VGrid,WGrid
        Type(Cell),intent(inout)            :: PCell,UCell,VCell,WCell
        type(BCBase),intent(inout)	    :: BCu,BCv,BCw,BCp,BCVof,BCLvs
        Type(Variables),intent(inout)       :: TVar,TVar_n
        Type(SolverTime),intent(inout)      :: Time
        Type(SolverConvergence),intent(out) :: UConv,VConv,WConv,PConv
        Integer(kind=it8b),intent(in)       :: itt
        Integer(kind=it4b)   		    :: i,j,k
        Real(kind=dp)  			    :: dt,mres
        Call ComputeTimeStep(UGrid,VGrid,WGrid,TVar,Time)
        dt = Time%dt
        TVar_n%p(:,:,:) = TVar%p(:,:,:)
        TVar_n%u(:,:,:) = TVar%u(:,:,:)
        TVar_n%v(:,:,:) = TVar%v(:,:,:)
        TVar_n%w(:,:,:) = TVar%w(:,:,:)
        dt = Time%dt!/3.d0
     !  First Runge-Kutta substep
        Call UpdatePUV(UGrid,VGrid,WGrid,PGrid,UCell,VCell,WCell,PCell,	       &
                       BCu,BCv,BCw,BCp,BCVof,BCLvsTVar_n,TVar,dt,itt)
     !  Call UpdatePUV(UGrid,VGrid,WGrid,PGrid,UCell,VCell,WCell,PCell,TVar_n, &
     !                                                          TVar,dt,itt)
        Call VariablesInternalCellCondition(TVar,PCell,UCell,VCell,WCell)
        ! print*, TVar%u(20,jbeg)
        ! Second Runge-Kutta substep
        ! Calculate the three kind of norm for convergence
        Call ResidualNormCalculate(UCell,TVar%u,TVar_n%u,TVar%ures,UConv)
        Call ResidualNormCalculate(VCell,TVar%v,TVar_n%v,TVar%vres,VConv)
        Call ResidualNormCalculate(WCell,TVar%w,TVar_n%w,TVar%wres,WConv)
        Call ResidualNormCalculate(PCell,TVar%p,TVar_n%p,TVar%pres,PConv)
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              TVar%mres(i,j,k) = PGrid%dy(i,j,k)*PGrid%dz(i,j,k)*              &
                                (PCell%EEArea(i,j,k)*TVar%u(i,j,k)-            &
                                 PCell%EEArea(i-1,j,k)*TVar%u(i-1,j,k))+       &
                                 PGrid%dx(i,j,k)*PGrid%dz(i,j,k)*              &
                                (PCell%NEArea(i,j,k)*TVar%v(i,j,k)-            &
                                 PCell%NEArea(i,j-1,k)*TVar%v(i,j-1,k))+       &
                                 PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*              &
                                (PCell%TEArea(i,j,k)*TVar%w(i,j,k)-            &
                                 PCell%TEArea(i,j,k-1)*TVar%w(i,j,k-1))
            End do
          End do
        End do
    End Subroutine AdamBasforthBDF2

    Subroutine ComputeTimeStep(UGrid,VGrid,WGrid,TVar,Time)
        Implicit none
        Type(Grid),intent(in):: UGrid,VGrid,WGrid
        Type(Variables),intent(in):: TVar
        Type(SolverTime),intent(out):: Time
        Integer(kind=it4b):: i,j,k
        Real(kind=dp):: tol
        tol = 1.d-20
        Time%dt = 1.d0
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              Time%dt = dmin1(Time%dt,                                         &
                              Time%cfl*Ugrid%dx(i,j,k)/dabs(TVar%u(i,j,k)+tol),&
                              Time%cfl*VGrid%dy(i,j,k)/dabs(TVar%v(i,j,k)+tol),&
                              Time%cfl*WGrid%dz(i,j,k)/dabs(TVar%w(i,j,k)+tol),&
                              Time%cfl*Ugrid%dx(1,j,k)/(TVar%Uint/TVar%Uref))
            End do
          End do
        End do
    End Subroutine ComputeTimeStep

    Subroutine ResidualNormCalculate(TCell,Varn1,Varn,Tres,Conv)
        Implicit none
        Type(Cell),intent(in):: TCell
        Type(SolverConvergence),intent(inout):: Conv
        Real(kind=dp),dimension(:,:,:),intent(in),allocatable:: Varn1,Varn
        Real(kind=dp),dimension(:,:,:),intent(inout),allocatable:: Tres
        Real(kind=dp):: N1,N2,Ninf,N1c,N2c,Ninfc
        Integer(kind=it4b):: i,j,k,numb
        N1 = 0.d0
        N2 = 0.d0
        Ninf = 0.d0
        N1c = 0.d0
        N2c = 0.d0
        Ninfc = 0.d0
        numb = 0
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              Tres(i,j,k) = dabs(Varn1(i,j,k)-Varn(i,j,k))
              N1 = N1+dabs(Varn1(i,j,k)-Varn(i,j,k))
              N2 = N2+(Varn1(i,j,k)-Varn(i,j,k))**2.d0
              Ninf = dmax1(Ninf,dabs(Varn1(i,j,k)-Varn(i,j,k)))
              If(TCell%Cell_Type(i,j,k)==1) then
                numb = numb+1
                N1c = N1c+dabs(Varn1(i,j,k)-Varn(i,j,k))
                N2c = N2c+(Varn1(i,j,k)-Varn(i,j,k))**2.d0
                Ninfc = dmax1(Ninfc,dabs(Varn1(i,j,k)-Varn(i,j,k)))
              End if
            End do
          End do
        End do
        N1 = N1/(Imax*Jmax*Kmax)
        N2 = dsqrt(N2/(Imax*Jmax*Kmax))
        N1c = N1c/dble(numb)
        N2c = dsqrt(N2c/dble(numb))
        Conv%Ninf = Dlog10(Ninf)
        Conv%N1 = Dlog10(N1)
        Conv%N2 = Dlog10(N2)
        Conv%Ninfc = Dlog10(Ninfc)
        Conv%N1c = Dlog(N1c)
        Conv%N2c = Dlog(N2c)
    End subroutine ResidualNormCalculate

    Subroutine PrintHistory(itt,TNorm)
       Integer(kind=it8b),intent(in):: itt
       Type(SolverConvergence),intent(in):: TNorm
       Open(unit=5,file='Convergence.dat',access='append')
       Write(5,76) itt,TNorm%N1,TNorm%N2,TNorm%Ninf,TNorm%N1c,TNorm%N2c,       &
                                                                   TNorm%Ninfc
       Close(5)
 76	   Format(I10,f15.10,f15.10,f15.10,f15.10,f15.10,f15.10)
    End subroutine

    Subroutine VariablesInternalCellCondition(TVar,PCell,UCell,VCell,WCell)
       Implicit none
       Type(Variables),intent(inout):: TVar
       Type(Cell),intent(in):: PCell,UCell,VCell,WCell
       Integer(kind=it4b):: i,j,k
       Do i = 1,Imax
         Do j = 1,Jmax
           Do k = 1,Kmax
             If(PCell%Cell_Type(i,j,k)==2) then
               TVar%p(i,j,k) = 0.d0
               TVar%t(i,j,k) = Twall/TVar%Tref
             End if
             If(UCell%Cell_Type(i,j,k)==2) then
               TVar%u(i,j,k) = 0.d0
             End if
             If(VCell%Cell_Type(i,j,k)==2) then
               TVar%v(i,j,k) = 0.d0
             End if
             If(WCell%Cell_Type(i,j,k)==2) then
               TVar%w(i,j,k) = 0.d0
             End if
           End do
         End do
       End do
     End subroutine VariablesInternalCellCondition

     Subroutine PrintDragLiftCoef(TVar,PGrid,UGrid,VGrid,WGrid,PCell,UCell,    &
                                                          VCell,WCell,itt,time)
      Implicit none
      type(Grid),intent(in):: PGrid,UGrid,VGrid,WGrid
      type(Variables),intent(in):: TVar
      type(Cell),intent(in):: PCell,UCell,VCell,WCell
      Integer(kind=it8b),intent(in):: itt
      Real(kind=dp):: time
      Integer(kind=it4b) i,j,k
      Real(kind=dp):: nx,ny,nz,area,Cdp1,Cdf,Pr,dpx,dpy,dpz,pw
      Real(kind=dp):: Clpy1,Clfy,Clpz1,Clfz,tol,cdp2,Clpy2,Clpz2,Cl1,Cl2
      Character(15) curd
      Open(unit=10,file='DragLiftCoef.dat',access='append')
      tol = 1.d-8
      Pr = 0.d0 !TVar%p(1,Jm2)
      Cdp1 = 0.d0;Cdp2 = 0.d0
      Cdf = 0.d0
      Clpy1 = 0.d0;Clpy2 = 0.d0
      Clpz1 = 0.d0;Clpz2 = 0.d0
      Clfy = 0.d0
      Clfz = 0.d0
    ! Calculate the Cdp
    ! For upper half of cylinder
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            If(PCell%Cell_Type(i,j,k)==1) then
              If(PCell%Cell_Type(i-1,j,k)/=1) then
                dpx = (TVar%p(i,j,k)-Tvar%p(i-1,j,k))/PGrid%dx(i,j,k)
              Else
                dpx = (TVar%p(i+1,j,k)-TVar%p(i,j,k))/PGrid%dx(i,j,k)
              End if
              If(PCell%Cell_Type(i,j-1,k)/=1) then
                dpy = (TVar%p(i,j,k)-Tvar%p(i,j-1,k))/PGrid%dy(i,j,k)
              Else
                dpy = (TVar%p(i,j+1,k)-TVar%p(i,j,k))/PGrid%dy(i,j,k)
              End if
              If(PCell%Cell_Type(i,j,k-1)/=1) then
                dpz = (TVar%p(i,j,k)-TVar%p(i,j,k-1))/PGrid%dz(i,j,k)
              Else
                dpz = (TVar%p(i,j,k+1)-TVar%p(i,j,k))/PGrid%dz(i,j,k)
              End if
              pw = TVar%p(i,j,k)
              nx = -PCell%nx(i,j,k)
              Cdp1 = Cdp1+pw*nx*PCell%WlLh(i,j,k)
              ny = -PCell%ny(i,j,k)
              Clpy1 = Clpy1+pw*ny*PCell%WlLh(i,j,k)
              nz = -PCell%nz(i,j,k)
              Clpz1 = Clpz1+pw*nz*PCell%WlLh(i,j,k)
              pw = TVar%p(i,j,k)-(dpx*PCell%phi(i,j,k)*PCell%nx(i,j,k)+        &
                                  dpy*PCell%phi(i,j,k)*PCell%ny(i,j,k)+        &
                                  dpz*PCell%phi(i,j,k)*PCell%nz(i,j,k))
              Cdp2 = Cdp2+pw*nx*PCell%WlLh(i,j,k)
              Clpy2 = Clpy2+pw*ny*PCell%WlLh(i,j,k)
              Clpz2 = Clpz2+pw*nz*PCell%WlLh(i,j,k)
            End if
            If(UCell%Cell_Type(i,j,k)==1) then
              nx = -UCell%nx(i,j,k)
              Cdf = Cdf+TVar%u(i,j,k)*UCell%WlLh(i,j,k)/UCell%delh(i,j,k)/Rey
            End if
            If(VCell%Cell_Type(i,j,k)==1) then
              ny = -VCell%ny(i,j,k)
              Clfy = Clfy+TVar%v(i,j,k)*VCell%WlLh(i,j,k)/VCell%delh(i,j,k)/Rey
            End if
            If(WCell%Cell_Type(i,j,k)==1) then
              nz = -WCell%nz(i,j,k)
              Clfz = Clfz+TVar%w(i,j,k)*WCell%WlLh(i,j,k)/WCell%delh(i,j,k)/Rey
            End if
          End do
        End do
      End do
      Cl1 = dsqrt((8.d0*(Clpy1+Clfy)/pi)**2.d0+(8.d0*(Clpz1+Clfz)/pi)**2.d0)
      Cl2 = dsqrt((8.d0*(Clpy2+Clfy)/pi)**2.d0+(8.d0*(Clpz2+Clfz)/pi)**2.d0)
      Write(10,76) itt,time,8.d0*(Cdp1+Cdf)/pi,8.d0*(Cdp2+Cdf)/pi,Cl1,Cl2
      Close(10)
 76	  Format(I10,5(f15.10))
    End subroutine PrintDragLiftCoef
End module Solver

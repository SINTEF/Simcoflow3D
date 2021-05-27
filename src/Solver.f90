Module Solver
    USE PrecisionVar
    USE Mesh
    USE StateVariables
    USE CutCell
    USE Clsvof
    USE PrintResult
    USE ComputePUVW
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

  !  subroutine  initialiseModules(PGrid,PCell,TVar,BCLvs,BCvof)

  !      Implicit none 
  !      Type(Grid),intent(in)         		   :: PGrid
  !      Type(Cell),intent(inout)      		   :: PCell
  !      Type(Variables),intent(inout) 		   :: TVar
  !      type(BCBase),intent(inout)		   :: BCVof,BCLvs

  !      call initialiseClsvof(PGrid,PCell,TVar,BCLvs,BCvof)

  !  end subroutine initialiseModules
    !
    ! solve the Navier-Stokes equations
    !
    subroutine IterationSolution(iprint,PGrid,UGrid,VGrid,WGrid, PCell,UCell,VCell,    &
                                 WCell,BCu,BCv,BCw,BCp,BCVof,BCLvs,TVar)
        !                 
        Implicit none
        !
        Integer(kind=it4b),    intent(in)    :: iprint
        Type(Grid),            intent(in)    :: PGrid, UGrid, VGrid, WGrid
        Type(Cell),            intent(inout) :: PCell, UCell, VCell, WCell
        Type(Variables),       intent(inout) :: TVar
        type(BCBase),          intent(inout) :: BCu, BCv, BCw, BCp, BCVof, BCLvs

        Integer(kind=it8b)                   :: itt
        Real(kind=dp),dimension(:,:,:,:),allocatable :: FluxDivOld
        Real(kind=dp),dimension(:,:,:),allocatable :: GraP
        Type(SolverTime)                     :: Time
        Type(SolverConvergence)              :: UConv, VConv, WConv, PConv
        Type(Variables)                      :: TVar_n
        Type(Cell)                           :: PCellO, UCellO, VCellO, WCellO
        !
        
        allocate(TVar_n%p(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
        allocate(TVar_n%u(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
        allocate(TVar_n%v(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
        allocate(TVar_n%w(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
        allocate(GraP(Imax,Jmax,Kmax))
        allocate(FluxDivOld(Imax,Jmax,Kmax,3))
        !
        call allocateOldCell(PCellO, UCellO, VCellO, WCellO)

        ! Call initialiseModules(PGrid,PCell,TVar,BCLvs,BCvof)

        !
        FluxDivOld(:,:,:,:) = 0.d0
        !
        Time%iter = 10**6
        Time%NondiT = 0.d0
        Time%Cfl = 0.3d0
        !
        ! Print out the information about numerical method
        !
        print*, 'The accuracy order of time discretization    :',TimeOrder
        print*, 'The accuracy order of spatial discretization :',SpaceOrder
        !
        ! start time loop
        !
        do itt = 1,Time%iter
          !
          ! solve the equations using Adam Basfort scheme
          !
          call AdamBasforthBDF2(itt,PGrid,UGrid,VGrid,WGrid,                    &
                                Time, PCell, UCell, VCell, WCell,                    &
                                PCellO, UCellO, VCellO, WCellO,                &
                                BCu, BCv, BCw, BCp, BCVof, BCLvs,              &
                                TVar, TVar_n, FluxDivOld,                      &
                                UConv, VConv, WConv, PConv)
          !              
          Time%NondiT = Time%NondiT+Time%dt
          Time%PhysT  = Time%Nondit*PGrid%Lref/TVar%URef
          !
          ! ??
          !
          call PrintHistory(itt,Uconv)
          call PrintDragLiftCoef(TVar, PGrid, UGrid, VGrid, WGrid,             &
                                 PCell, UCell, VCell, WCell, itt, Time%NondiT)
          !               
          ite = itt
          ! print*, itt
          !
          ! save fields for visualization
          !
          if(mod(itt,iprint)==0)then 
            !      
            write(*,*), itt,Time%PhysT,Time%NondiT
            ! call PrintResultVTK(PGrid,TVar,PCell,itt)
            call PrintResultVTR3D(PGrid,TVar,PCell,"FlowFieldP",itt)
            call PrintResultVTR3D(UGrid,TVar,UCell,"FlowFieldU",itt)
            ! call PrintResultTecplotPCent(PGrid,TVar,PCell,itt)
            ! call PrintResultTecplotPCentXY(PGrid,TVar,PCell,itt)
            ! call PrintResultTecplotPCentXZ(PGrid,TVar,PCell,itt)
            ! call PrintResultTecplotUCent(UGrid,TVar,UCell,itt)
            ! call PrintResultTecplotVCent(VGrid,TVar,VCell,itt)
            ! call PrintResultTecplotWCent(WGrid,TVar,WCell,itt)
          end if
          !
        end do
        !
        ! deallocation
        !
        if(allocated(GraP)) then
          !      
          deallocate(GraP)
          !
        else
          !      
          print*,'Grap is unlocated Solver 87'
          !
        end if
        !
        if(allocated(TVar_n%p).and.allocated(TVar_n%u).and.                    &
           allocated(TVar_n%v).and.allocated(TVar_n%w)) then
          ! 
          deallocate(TVar_n%p,TVar_n%u,TVar_n%v,TVar_n%w)
          !
        else
          !      
          print*,'Tvar_n is unlocated Solver 93'
          !
        end if
        !
        if(allocated(FluxDivOld)) deallocate(FluxDivOld)
        !
    end subroutine IterationSolution
    !
    Subroutine AdamBasforthBDF2(itt, PGrid, UGrid, VGrid, WGrid,                    &
                                Time, PCell, UCell, VCell, WCell,                    &
                                PCellO, UCellO, VCellO, WCellO,                &
                                BCu, BCv, BCw, BCp, BCVof, BCLvs,              &
                                TVar, TVar_n, FluxDivOld,                      &
                                UConv, VConv, WConv, PConv)
        !                
        Implicit none
        !
        Integer(kind=it8b),               intent(in)    :: itt
        Type(Grid),                       intent(in)    :: PGrid, UGrid, VGrid, WGrid
        Type(Cell),                       intent(inout) :: PCell, UCell, VCell, WCell
        type(Cell),                       intent(inout) :: PCellO, UCellO, VCellO, WCellO
        type(BCBase),                     intent(inout) :: BCu, BCv, BCw, BCp, BCVof, BCLvs
        Type(Variables),                  intent(inout) :: TVar, TVar_n
        real(kind=dp),dimension(:,:,:,:),allocatable,intent(inout) :: FluxDivOld
        Type(SolverTime),                 intent(inout) :: Time
        Type(SolverConvergence),          intent(out)   :: UConv, VConv, WConv, PConv

        Integer(kind=it4b)                              :: i,j,k
        Real(kind=dp)                                   :: dt,mres
        !
        Call ComputeTimeStep(itt,UGrid,VGrid,WGrid,TVar, Time)
        !
        dt = Time%dt
        !
        if(itt==1) then
          TVar_n%p(:,:,:) = TVar%p(:,:,:)
          TVar_n%u(:,:,:) = TVar%u(:,:,:)
          TVar_n%v(:,:,:) = TVar%v(:,:,:)
          TVar_n%w(:,:,:) = TVar%w(:,:,:)
        end if  
        !
        print*, 'Solver.f90 151'
        print*, 'Time step size:', dt
        !
        ! Compute the previous cell configuration 
        !
        call CopyOldCellNewCell(PCell, PCellO)
        call CopyOldCellNewCell(Ucell, UCellO)
        call CopyOldCellNewCell(VCell, VCellO)
        call CopyOldCellNewCell(WCell, WCellO)
        !
        ! if(itt>1) then
        !   call Clsvof_Scheme(PGrid,PCell,TVar,BCu,BCv,BCw,BCLvs,BCvof,         &
        !                                                 Time%NondiT,dt,itt)
        !   call ComputeUVWLiquidField(PGrid,PCell,UCell,VCell,WCell,            &
        !                                          UGrid,VGrid,WGrid)
        ! end if
        !
        ! update the velocity and pressure
        !
        call UpdatePUVW(itt, Time%NondiT, dt, UGrid, VGrid, WGrid, PGrid,      &
                        UCell, VCell, WCell, PCell,                            &
                        UCellO, VCellO, WCellO, PCellO,                        &
                        BCu, BCv, BCw, BCp, BCVof, BCLvs,                      &
                        FluxDivOld, TVar_n, TVar)
        !
        ! correction for internal cells, zero fields
        !
        call VariablesInternalCellCondition(TVar, PCell, UCell, VCell, WCell)
        !
        ! Calculate the three kind of norm for convergence
        !
        ! call ResidualNormCalculate(UCell,TVar%u,TVar_n%u,TVar%ures,UConv)
        ! call ResidualNormCalculate(VCell,TVar%v,TVar_n%v,TVar%vres,VConv)
        ! call ResidualNormCalculate(WCell,TVar%w,TVar_n%w,TVar%wres,WConv)
        ! call ResidualNormCalculate(PCell,TVar%p,TVar_n%p,TVar%pres,PConv)
        !
        ! compute mass error
        !
        do i = 1,Imax
          do j = 1,Jmax
            do k = 1,Kmax
              TVar%mres(i,j,k) = PGrid%dy(i,j,k)*PGrid%dz(i,j,k)*              &
                                (PCell%EEArea(i,j,k)*TVar%u(i,j,k)-            &
                                 PCell%EEArea(i-1,j,k)*TVar%u(i-1,j,k))+       &
                                 PGrid%dx(i,j,k)*PGrid%dz(i,j,k)*              &
                                (PCell%NEArea(i,j,k)*TVar%v(i,j,k)-            &
                                 PCell%NEArea(i,j-1,k)*TVar%v(i,j-1,k))+       &
                                 PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*              &
                                (PCell%TEArea(i,j,k)*TVar%w(i,j,k)-            &
                                 PCell%TEArea(i,j,k-1)*TVar%w(i,j,k-1))
            end do
          end do
        end do
        !
    end Subroutine AdamBasforthBDF2
    !
    Subroutine ComputeTimeStep(itt,UGrid,VGrid,WGrid,TVar, Time)
        !    
        implicit none
        !
        integer(kind=it8b),intent(in) :: itt
        type(Grid),intent(in)	      :: UGrid,VGrid,WGrid
        type(Variables),intent(in)    :: TVar
        type(SolverTime),intent(out)  :: Time

        integer(kind=it4b)	      :: i,j,k
        real(kind=dp)		      :: tol
        !
        tol = 1.d-20
        Time%dt = 1.d0
        !
        do i = 1,Imax
          do j = 1,Jmax
            do k = 1,Kmax
              !
              ! Time%dt=dmin1(Time%dt,                                           &
              !               Time%cfl*Ugrid%dx(i,j,k)/dabs(TVar%u(i,j,k)+tol),&
              !               Time%cfl*VGrid%dy(i,j,k)/dabs(TVar%v(i,j,k)+tol),&
              !               Time%cfl*WGrid%dz(i,j,k)/dabs(TVar%w(i,j,k)+tol),&
              !               Time%cfl*Ugrid%dx(1,j,k)/(TVar%Uint/TVar%Uref))
              Time%dt=dmin1(Time%dt,2.d0*Time%cfl*VGrid%dx(i,j,k)/	       &
                     (dabs(TVar%u(i,j,k))+dsqrt(TVar%u(i,j,k)**2.d0+	       &
                      4.d0*VGrid%dx(i,j,k)*dabs(gx/g)/Fr)+tol))
              Time%dt=dmin1(Time%dt,2.d0*Time%cfl*VGrid%dy(i,j,k)/	       &
                     (dabs(TVar%v(i,j,k))+dsqrt(TVar%v(i,j,k)**2.d0+	       &
                      4.d0*VGrid%dy(i,j,k)*dabs(gy/g)/Fr)+tol))
              Time%dt=dmin1(Time%dt,2.d0*Time%cfl*VGrid%dz(i,j,k)/	       &
                     (dabs(TVar%w(i,j,k))+dsqrt(TVar%w(i,j,k)**2.d0+	       &
                      4.d0*VGrid%dz(i,j,k)*dabs(gz/g)/Fr)+tol))           
              !
            end do
          end do
        end do
        !
        if(itt==1) then
          !      
          do i=1,Imax
            do j=1,Jmax
              do k=1,Kmax
                !
                Time%dt=dmin1(Time%dt,Time%cfl*UGrid%dx(i,j,k))
                !
              end do
            end do
          end do
          !
        end if        
        !
    end Subroutine ComputeTimeStep
    !
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
    !
    Subroutine PrintHistory(itt,TNorm)
       !     
       Integer(kind=it8b),intent(in):: itt
       Type(SolverConvergence),intent(in):: TNorm
       !
       Open(unit=5,file='Convergence.dat',access='append')
       Write(5,76) itt,TNorm%N1,TNorm%N2,TNorm%Ninf,TNorm%N1c,TNorm%N2c,       &
                                                              TNorm%Ninfc
       Close(5)
       !
 76    Format(I10,f15.10,f15.10,f15.10,f15.10,f15.10,f15.10)
       !
    End subroutine
    !
    Subroutine VariablesInternalCellCondition(PCell,UCell,VCell,WCell, TVar)
       !     
       Implicit none
       !
       Type(Cell),intent(in):: PCell,UCell,VCell,WCell
       Type(Variables),intent(inout):: TVar

       Integer(kind=it4b):: i,j,k
       !
       Do i = 1,Imax
         Do j = 1,Jmax
           Do k = 1,Kmax
             !
             If(PCell%Cell_Type(i,j,k)==2) then
               TVar%p(i,j,k) = 0.d0
               TVar%t(i,j,k) = Twall/TVar%Tref
             End if
             !
             If(UCell%Cell_Type(i,j,k)==2) then
               TVar%u(i,j,k) = 0.d0
             End if
             !
             If(VCell%Cell_Type(i,j,k)==2) then
               TVar%v(i,j,k) = 0.d0
             End if
             !
             If(WCell%Cell_Type(i,j,k)==2) then
               TVar%w(i,j,k) = 0.d0
             End if
             !
           End do
         End do
       End do
       !
     End subroutine VariablesInternalCellCondition
     !
     subroutine AllocateOldCell(PCellO,UCellO,VCellO,WCellO)
       !      
       Implicit none
       !
       type(cell), intent(inout) :: PCellO, UCellO, VCellO, WCellO
       !
       allocate(PCellO%Cell_Type(Imax,Jmax,Kmax))
       allocate(PCellO%vof(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%phi(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%nx(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%ny(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%nz(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%nxL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%nyL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%nzL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
       allocate(PCellO%Cell_Cent(Imax,Jmax,Kmax,3))
       allocate(PCellO%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(PCellO%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))

       allocate(UCellO%Cell_Type(Imax,Jmax,Kmax))
       allocate(UCellO%vof(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%phi(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%nx(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%ny(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%nz(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%nxL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%nyL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%nzL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
       allocate(UCellO%Cell_Cent(Imax,Jmax,Kmax,3))
       allocate(UCellO%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%MoExCell(Imax,Jmax,Kmax))
       allocate(UCellO%EtaE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%EtaN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%EtaT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%AlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%AlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%AlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%SxE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%SyN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%SzT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(UCellO%FCE(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
       allocate(UCellO%FCN(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
       allocate(UCellO%FCT(0:Imax+1,0:Jmax+1,0:Kmax+1,3))

       allocate(VCellO%Cell_Type(Imax,Jmax,Kmax))
       allocate(VCellO%vof(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%phi(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%nx(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%ny(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%nz(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%nxL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%nyL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%nzL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
       allocate(VCellO%Cell_Cent(Imax,Jmax,Kmax,3))
       allocate(VCellO%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%MoExCell(Imax,Jmax,Kmax))
       allocate(VCellO%EtaE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%EtaN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%EtaT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%AlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%AlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%AlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%SxE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%SyN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%SzT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(VCellO%FCE(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
       allocate(VCellO%FCN(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
       allocate(VCellO%FCT(0:Imax+1,0:Jmax+1,0:Kmax+1,3))

       allocate(WCellO%Cell_Type(Imax,Jmax,Kmax))
       allocate(WCellO%vof(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%phi(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%nx(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%ny(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(wCellO%nz(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%nxL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%nyL(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(wCellO%nzL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
       allocate(WCellO%Cell_Cent(Imax,Jmax,Kmax,3))
       allocate(WCellO%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%MoExCell(Imax,Jmax,Kmax))
       allocate(WCellO%EtaE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%EtaN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%EtaT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%AlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%AlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%AlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%SxE(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%SyN(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(WCellO%SzT(0:Imax+1,0:Jmax+1,0:Kmax+1))
       allocate(wCellO%FCE(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
       allocate(WCellO%FCN(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
       allocate(WCellO%FCT(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
       !
     end subroutine AllocateOldCell
     !
     subroutine CopyOldCellNewCell(TCell, TCellO)
       !
       type(cell), intent(in)    :: TCell
       type(cell), intent(inout) :: TCellO
       !
       TCellO%Cell_Type(:,:,:)   = TCell%Cell_Type(:,:,:)
       TCellO%vof(:,:,:)         = TCell%vof(:,:,:)
       TCellO%phi(:,:,:)         = TCell%phi(:,:,:)
       TCellO%nx(:,:,:)          = TCell%nx(:,:,:)
       TCellO%ny(:,:,:)          = TCell%ny(:,:,:)
       TCellO%nz(:,:,:)          = TCell%nz(:,:,:)
       TCellO%vofL(:,:,:)        = TCell%vofL(:,:,:)
       TCellO%phiL(:,:,:)        = TCell%phiL(:,:,:)
       TCellO%nxL(:,:,:)         = TCell%nxL(:,:,:)
       TCellO%nyL(:,:,:)         = TCell%nyL(:,:,:)
       TCellO%nzL(:,:,:)         = TCell%nzL(:,:,:)
       TCellO%Cell_Cent(:,:,:,:) = TCell%Cell_Cent(:,:,:,:)
       TCellO%EEArea(:,:,:)      = TCell%EEArea(:,:,:)
       TCellO%NEArea(:,:,:)      = TCell%NEArea(:,:,:)
       TCellO%TEArea(:,:,:)      = TCell%TEArea(:,:,:)
       if(allocated(TCellO%FCE))      TCellO%FCE(:,:,:,:)    = TCell%FCE(:,:,:,:)
       if(allocated(TCellO%FCN))      TCellO%FCN(:,:,:,:)    = TCell%FCN(:,:,:,:)
       if(allocated(TCellO%FCT))      TCellO%FCT(:,:,:,:)    = TCell%FCT(:,:,:,:)
       if(allocated(TCellO%MoExCell)) TCellO%MoExCell(:,:,:) = TCell%MoExCell(:,:,:)
       if(allocated(TCellO%EtaE))     TCellO%EtaE(:,:,:)     = TCell%EtaE(:,:,:)
       if(allocated(TCellO%EtaN))     TCellO%EtaN(:,:,:)     = TCell%EtaN(:,:,:)
       if(allocated(TCellO%EtaT))     TCellO%EtaT(:,:,:)     = TCell%EtaT(:,:,:)
       if(allocated(TCellO%AlE))      TCellO%AlE(:,:,:)      = TCell%AlE(:,:,:)
       if(allocated(TCellO%AlN))      TCellO%AlN(:,:,:)      = TCell%AlN(:,:,:)
       if(allocated(TCellO%AlT))      TCellO%AlT(:,:,:)      = TCell%AlT(:,:,:)
       if(allocated(TCellO%SxE))      TCellO%SxE(:,:,:)      = TCell%SxE(:,:,:)
       if(allocated(TCellO%SyN))      TCellO%SyN(:,:,:)      = TCell%SyN(:,:,:)
       if(allocated(TCellO%SzT))      TCellO%SzT(:,:,:)      = TCell%SzT(:,:,:)
       !
     end subroutine CopyOldCellNewCell 
     !
     Subroutine PrintDragLiftCoef(TVar,PGrid,UGrid,VGrid,WGrid,PCell,UCell,    &
                                                          VCell,WCell,itt,time)
      !                                            
      Implicit none
      !
      type(Grid),         intent(in):: PGrid,UGrid,VGrid,WGrid
      type(Variables),    intent(in) :: TVar
      type(Cell),         intent(in):: PCell,UCell,VCell,WCell
      Integer(kind=it8b), intent(in):: itt

      Character(15) curd
      Integer(kind=it4b) :: i,j,k
      Real(kind=dp)      :: time
      Real(kind=dp)      :: nx,ny,nz,area,Cdp1,Cdf,Pr,dpx,dpy,dpz,pw
      Real(kind=dp)      :: Clpy1,Clfy,Clpz1,Clfz,tol,cdp2,Clpy2,Clpz2,Cl1,Cl2
      !
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
 76   Format(I10,5(f15.10))
      !
    End subroutine PrintDragLiftCoef
    !
End module Solver

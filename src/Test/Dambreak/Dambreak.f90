Program Main
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE Clsvof
    USE StateVariables
    USE PrintResult
    USE MPI
    USE Solver
    USE ComputePUVW
    USE BoundaryInterface
    USE BoundaryFunction
    USE InitialVof
    
    Implicit none
    
    Type(Grid)      :: UGrid,VGrid,WGrid,PGrid
    Type(Cell)      :: UCell,VCell,WCell,PCell
    Type(Point)     :: SPoint,EPoint,ReS,ReE
    Type(Variables) :: Var
    Type(BCBase)    :: BCp, BCu, BCv, BCw, BCVof, BCLvs, BCVofF, BCLvsF
    
    Integer(kind=it4b) :: Irec,Jrec,Krec,NI,NJ,NK,iprint
    Real(kind=dp)      :: Lref,vel
    real(kind=dp), dimension(:), allocatable :: Constin
    
    allocate(Constin(6))
    Open(unit=5,file='/home/sontd/code/CutCell3DGFMCLSVOF/src/Test/Dambreak/input.dat',action='read')
    Read(5,*),
    Read(5,*), Imax, Jmax, Kmax, Irec, Jrec, Krec, Rey, Lref, iprint
    Read(5,*),
    Read(5,*), TimeOrder, SpaceOrder
    close(5)
   
    
    NI = Imax+1
    NJ = Jmax+1
    NK = Kmax+1
    gx = 0.d0
    gy = 0.d0
    gz = -g
    vel = dsqrt(dabs(g*Lref))
  !  Imax = 98
  !  Jmax = 98
  !  Kmax = 98
    SPoint%x = 0.d0
    SPoint%y = 0.d0
    SPoint%z = 0.d0
    EPoint%x = 1.61d0
    EPoint%y = 0.15d0
    EPoint%z = 0.60d0
    ReS%x = -1.d0
    ReS%y = -1.d0
    ReS%z = -1.d0
    ReE%x = 1.d0
    ReE%y = 1.d0
    ReE%z = 1.d0
    
    BCp = BCBase(IMax,Jmax,Kmax)
    BCu = BCBase(Imax,Jmax,Kmax)
    BCv = BCBase(Imax,Jmax,Kmax)
    BCw = BCBase(Imax,Jmax,Kmax)
    BCvof = BCBase(Imax,Jmax,Kmax)
    BClvs = BCBase(Imax,Jmax,Kmax)
    BCVofF = BCBase(Imax,Jmax,Kmax)
    BCLvsF = BCBase(Imax,Jmax,Kmax)
    
    ! For dambreak flow 
    ! 0 for Dirichlet, 1 for Neumann Boundary Condition
    call BCp%SetDN(1,1,1,1,1,0)
    call BCu%SetDN(0,0,0,0,0,1)
    call BCv%SetDN(0,0,0,0,0,1)
    call BCw%SetDN(0,0,0,0,0,1)
    call BCVof%SetDN(1,1,1,1,1,1)
    call BCLvs%SetDN(1,1,1,1,1,1)
    call BCVofF%SetDN(1,1,1,1,1,1)
    call BCLvsF%SetDN(1,1,1,1,1,1)
    ! Set Constant for boundary condition
    Constin(:) = 0.d0
    call BCu%SetConstant(Constin)    
    Constin(:) = 0.d0
    call BCv%SetConstant(Constin)
    call BCw%SetConstant(Constin)
    call BCp%SetConstant(Constin)
    call BCVof%SetConstant(Constin)
    call BCLvs%SetConstant(Constin)
    call BCVofF%SetConstant(Constin)
    call BCLvsF%SetConstant(Constin)
    
    BCu%West   => BCUW
    BCu%East   => BCUE
    BCu%South  => BCUS
    BCu%North  => BCUN
    BCu%Bottom => BCUB
    BCu%Top    => BCUT
    
    BCv%West   => BCVW
    BCv%East   => BCVE
    BCv%South  => BCVS
    BCv%North  => BCVN
    BCv%Bottom => BCVB
    BCv%Top    => BCVT
    
    BCw%West   => BCWW
    BCw%East   => BCWE
    BCw%South  => BCWS
    BCw%North  => BCWN
    BCw%Bottom => BCWB
    BCw%Top    => BCWT
    
    BCp%West   => BCPW
    BCp%East   => BCPE
    BCp%South  => BCPS
    BCp%North  => BCPN
    BCp%Bottom => BCPB
    BCp%Top    => BCPT
    
    BCVof%West   => BCVofW
    BCVof%East   => BCVofE
    BCVof%South  => BCVofS
    BCVof%North  => BCVofN
    BCVof%Bottom => BCVofB
    BCVof%Top    => BCVofT
    
    BCLvs%West   => BCLvsW
    BCLvs%East   => BCLvsE
    BCLvs%South  => BCLvsS
    BCLvs%North  => BCLvsN
    BCLvs%Bottom => BCLvsB
    BCLvs%Top    => BCLvsT

    BCVofF%West   => BCVofW
    BCVofF%East   => BCVofE
    BCVofF%South  => BCVofS
    BCVofF%North  => BCVofN
    BCVofF%Bottom => BCVofB
    BCVofF%Top    => BCVofT

    BCLvsF%West   => BCLvsW
    BCLvsF%East   => BCLvsE
    BCLvsF%South  => BCLvsS
    BCLvsF%North  => BCLvsN
    BCLvsF%Bottom => BCLvsB
    BCLvsF%Top    => BCLvsT
    deallocate(Constin)
    Call AllocateVar(Pgrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,Var)
    Call InitialGridDamBreak(SPoint,EPoint,PGrid,Lref)
    Call InitialUVGrid(PGrid,UGrid,0,Lref)
    Call InitialUVGrid(PGrid,VGrid,1,Lref)
    Call InitialUVGrid(PGrid,WGrid,2,Lref)
    Call MPI_Initial
    Call HYPRE_CreateGrid(PGrid)
    Call InitialClsvofFluidFieldDamBreak(PGrid,PCell)
    Call InitialClsvofFluidFieldDamBreak(UGrid,UCell)
    Call InitialClsvofFluidFieldDamBreak(VGrid,VCell)
    Call InitialClsvofFluidFieldDamBreak(WGrid,WCell)
    
    Call InitialClsvofLiquidFieldDamBreak(PGrid,PCell)
    call ComputeUVWLiquidField(PGrid,PCell,UCell,VCell,WCell,            &
                                                 UGrid,VGrid,WGrid)
 !   Call InitialClsvofLiquidFieldDamBreak(UGrid,UCell)
 !   Call InitialClsvofLiquidFieldDamBreak(VGrid,VCell)
 !   Call InitialClsvofLiquidFieldDamBreak(WGrid,WCell)

    call BoundaryConditionLvsVof(PGrid, PCell, Var, BCLvs, BCVof, 0.d0)
    call BoundaryConditionLvsVofFluid(PGrid, PCell, Var, BCLvsF,BCVofF, 0.d0)
 !   Call PrintResultTecplotPCent(PGrid,Var,PCell,INT8(0))
 !   Call PrintResultTecplotUCent(UGrid,Var,UCell,INT8(0))
 !   Call PrintResultTecplotVCent(VGrid,Var,VCell,INT8(0))
 !   Call PrintResultTecplotWCent(WGrid,Var,WCell,INT8(0))
 !   Call PrintResultVTK(PGrid,Var,PCell,INT8(0)) 
    Call PrintResultVTR3D(PGrid,Var,PCell,"FlowField",INT8(0))
    Call PrintResultVTR3D(UGrid,Var,UCell,"FlowFieldU",INT8(0))
    Call PrintResultVTR3D(VGrid,Var,VCell,"FlowFieldV",INT8(0))
    Call PrintResultVTR3D(WGrid,Var,WCell,"FlowFieldW",INT8(0))
    
    Call GridPreProcess(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,int8(1))
    Call DefineMomentumExchangeCell(PCell,UCell,VCell,WCell)
    
    Call NumberExternalCell(PCell,0,0,0)
    Call NumberExternalCell(UCell,1,0,0)
    Call NumberExternalCell(VCell,0,1,0)
    Call NumberExternalCell(WCell,0,0,1)
    Call NewCellFace(PCell,UCell,VCell,WCell,PGrid,UGrid,VGrid,WGrid)
    
    Call InitialVarDambreak(Var,vel,0.d0,0.d0,0.d0,300.d0,vel,300.d0,row,Lref)
    call BoundaryConditionVarNew(PGrid, PCell, Var, BCp, BCu, BCv, BCw, 0.d0)
    Call IterationSolution(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,    &
                           BCu,BCv,BCw,BCp,BCVof,BCLvs,Var,1)
    Pause
End program main

Subroutine AllocateVar(Pgrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,Var)
    USE Mesh
    USE Cutcell
    USE StateVariables
    Implicit none
    Type(Grid),intent(inout):: PGrid,UGrid,VGrid,WGrid
    Type(Cell),intent(inout):: PCell,UCell,VCell,WCell
    Type(Variables),intent(inout):: Var
    Allocate(UGrid%x(Imax,Jmax,Kmax))
    Allocate(UGrid%y(Imax,Jmax,Kmax))
    Allocate(UGrid%z(Imax,Jmax,Kmax))
    Allocate(VGrid%x(Imax,Jmax,Kmax))
    Allocate(VGrid%y(Imax,Jmax,Kmax))
    Allocate(VGrid%z(Imax,Jmax,Kmax))
    Allocate(WGrid%x(Imax,Jmax,Kmax))
    Allocate(WGrid%y(Imax,Jmax,Kmax))
    Allocate(WGrid%z(Imax,Jmax,Kmax))
    Allocate(PGrid%x(Imax,Jmax,Kmax))
    Allocate(PGrid%y(Imax,Jmax,Kmax))
    Allocate(PGrid%z(Imax,Jmax,Kmax))

    Allocate(UGrid%dx(Imax,Jmax,Kmax))
    Allocate(UGrid%dy(Imax,Jmax,Kmax))
    Allocate(UGrid%dz(Imax,Jmax,Kmax))
    Allocate(VGrid%dx(Imax,Jmax,Kmax))
    Allocate(VGrid%dy(Imax,Jmax,Kmax))
    Allocate(VGrid%dz(Imax,Jmax,Kmax))
    Allocate(WGrid%dx(Imax,Jmax,Kmax))
    Allocate(WGrid%dy(Imax,Jmax,Kmax))
    Allocate(WGrid%dz(Imax,Jmax,Kmax))
    Allocate(PGrid%dx(Imax,Jmax,Kmax))
    Allocate(PGrid%dy(Imax,Jmax,Kmax))
    Allocate(PGrid%dz(Imax,Jmax,Kmax))

    allocate(UCell%Cell_Type(Imax,Jmax,Kmax))
    allocate(UCell%vof(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%phi(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nx(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%ny(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nz(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nxL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nyL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nzL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
    Allocate(UCell%Cell_Cent(Imax,Jmax,Kmax,3))
    Allocate(UCell%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%MoExCell(Imax,Jmax,Kmax))
    Allocate(UCell%EtaE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%EtaN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%EtaT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%DAlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%DAlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%DAlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%AlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%AlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%AlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%SxE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%SyN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%SzT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(UCell%FCE(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(UCell%FCN(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(UCell%FCT(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(UCell%PosNu(Imax-1,Jmax,Kmax))
    Allocate(UCell%WlLh(Imax,Jmax,Kmax))
    Allocate(UCell%delh(Imax,Jmax,Kmax))
    Allocate(UCell%MsCe(Imax,Jmax,Kmax,3))
    Allocate(UCell%WePr(Imax-1,Jmax,Kmax))

    allocate(VCell%Cell_Type(Imax,Jmax,Kmax))
    allocate(VCell%vof(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%phi(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nx(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%ny(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nz(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nxL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nyL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nzL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
    allocate(VCell%Cell_Cent(Imax,Jmax,Kmax,3))
    allocate(VCell%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%MoExCell(Imax,Jmax,Kmax))
    allocate(VCell%EtaE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%EtaN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%EtaT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%AlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%AlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%AlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%DAlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%DAlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%DAlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%SxE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%SyN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%SzT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%FCE(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    allocate(VCell%FCN(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    allocate(VCell%FCT(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    allocate(VCell%PosNu(Imax,Jmax-1,Kmax))
    allocate(VCell%WlLh(Imax,Jmax,Kmax))
    allocate(VCell%delh(Imax,Jmax,Kmax))
    allocate(VCell%MsCe(Imax,Jmax,Kmax,3))
    allocate(VCell%WePr(Imax,Jmax,Kmax))

    Allocate(WCell%Cell_Type(Imax,Jmax,Kmax))
    Allocate(WCell%vof(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%phi(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%nx(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%ny(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(wCell%nz(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%nxL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%nyL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(wCell%nzL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
    Allocate(WCell%Cell_Cent(Imax,Jmax,Kmax,3))
    Allocate(WCell%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%MoExCell(Imax,Jmax,Kmax))
    Allocate(WCell%EtaE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%EtaN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%EtaT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%AlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%AlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%AlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%DAlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%DAlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%DAlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%SxE(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%SyN(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%SzT(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(wCell%FCE(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(WCell%FCN(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(WCell%FCT(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(WCell%PosNu(Imax,Jmax,Kmax-1))
    Allocate(WCell%WlLh(Imax,Jmax,Kmax))
    Allocate(WCell%delh(Imax,Jmax,Kmax))
    Allocate(WCell%MsCe(Imax,Jmax,Kmax,3))
    Allocate(WCell%WePr(Imax,Jmax,Kmax))

    Allocate(PCell%Cell_Type(Imax,Jmax,Kmax))
    Allocate(PCell%vof(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%phi(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nx(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%ny(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nz(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nxL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nyL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nzL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
    Allocate(PCell%Cell_Cent(Imax,Jmax,Kmax,3))
    Allocate(PCell%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
    
    Allocate(PCell%FCE(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(PCell%FCN(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(PCell%FCT(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
    Allocate(PCell%PosNu(Imax,Jmax,Kmax))
    Allocate(PCell%WlLh(Imax,Jmax,Kmax))
    Allocate(PCell%delh(Imax,Jmax,Kmax))

    Allocate(Var%u(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%v(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%w(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%p(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%Gpu(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%Gpv(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%Gpw(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%ures(Imax,Jmax,Kmax))
    Allocate(Var%vres(Imax,Jmax,Kmax))
    Allocate(Var%wres(Imax,Jmax,Kmax))
    Allocate(Var%pres(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%mres(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%t(0:Imax+1,0:Jmax+1,0:Kmax+1))
 End subroutine

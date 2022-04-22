    !*******************************************************
    !              dp/dn = 0.d0,du/dn=0, v = 0
    !    slip wall________________________
    !           --->                     --->
    !           --->                     --->
    !           --->                     --->
    ! => inlet : u = 5.d-3          => outflow :: du/dn = 0.d0
    !          : v = 0.d0                 --->    dv/dn = 0.d0
    !          : dp/dn = 0.d0             --->    P = 0.d0
    !           --->                      --->
    !    slip wall________________________
    !             dp/dn = 0.d0, u,v = 0.d0
    !*******************************************************
Program Main
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE Clsvof
    USE StateVariables
    !USE PrintResult
    USE MPI
    USE Solver
    USE ComputePUVW
    USE BoundaryInterface
    USE BoundaryFunction
    
    Implicit none
    
    Type(Grid)      :: UGrid,VGrid,WGrid,PGrid
    Type(Cell)      :: UCell,VCell,WCell,PCell
    Type(Point)     :: SPoint,EPoint,ReS,ReE
    Type(Variables) :: Var
    Type(BCBase)    :: BCp, BCu, BCv, BCw, BCVof, BCLvs, BCVofF, BCLvsF
    
    integer(it4b) :: i,j,k
    Integer(kind=it4b) :: Irec,Jrec,Krec,NI,NJ,NK,iprint
    Real(kind=dp)      :: Lref,vel
    real(kind=dp), dimension(:), allocatable :: Constin
    
    allocate(Constin(6))
    Open(unit=5,file='/home/elena-roxanap/Documents/Iceload/simco3d/src/Test/FlowSphere/input.dat',action='read')
    Read(5,*)
    Read(5,*) Imax, Jmax, Kmax, Irec, Jrec, Krec, Rey, Lref, iprint
    Read(5,*)
    Read(5,*) TimeOrder, SpaceOrder
    close(5)

    gx = 0.d0
    gy = 0.d0
    gz = -g
    print*, 'gx,gy,gz', gx,gy,gz
    
    Ta = 1000.d0
    wa = dsqrt(2.d0*Ta*nu**2.d0/((R1+R2)*(R2-R1)**3.d0))
    
    xc = 0.d0
    yc = 0.d0
    zc = 0.d0
    
    vel = dble(Rey*nuw/Lref/row)
    print*, "vel", vel, Rey, nuw,Lref, row
    
    NI = Imax+1
    NJ = Jmax+1
    NK = Kmax+1
    
    SPoint%x = -20.d0
    SPoint%y = -20.d0
    SPoint%z = -20.d0
    EPoint%x = 40.d0
    EPoint%y = 20.d0
    EPoint%z = 20.d0
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
    BCvofF = BCBase(Imax,Jmax,Kmax)
    BClvsF = BCBase(Imax,Jmax,Kmax)


    ! For flow over sphere
    call BCp%SetDN(1,0,1,1,1,1)
    call BCu%SetDN(0,1,0,0,0,0)
    call BCv%SetDN(0,1,0,0,0,0)
    call BCw%SetDN(0,1,0,0,0,0)
    call BCVof%SetDN(0,0,1,1,1,1)
    call BCLvs%SetDN(0,0,1,1,1,1)
    ! Set Constant for boundary condition
    Constin(1) = vel
    Constin(2:6) = 0.d0
    call BCu%SetConstant(Constin)    
    Constin(:) = 0.d0
    call BCv%SetConstant(Constin)
    call BCw%SetConstant(Constin)
    call BCp%SetConstant(Constin)
    Constin(1)=1.d0
    Constin(2)=1.d0
    call BCVof%SetConstant(Constin)
    call BCVofF%SetConstant(Constin)
    Constin(1)=-1.d3
    Constin(2)=-1.d3
    call BCLvs%SetConstant(Constin)
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

    Call InitialGrid(NI,NJ,NK,Irec,Jrec,Krec,Lref,SPoint,EPoint,ReS,ReE, PGrid)
    !Call InitialGrid(SPoint,EPoint,ReS,ReE,NI,NJ,NK,Irec,Jrec,Krec,PGrid,Lref)

    Call InitialUVGrid(0,Lref,PGrid, UGrid)
    Call InitialUVGrid(1,Lref,PGrid, VGrid)
    Call InitialUVGrid(2,Lref,PGrid, WGrid)

    Call MPI_Initial

    Call HYPRE_CreateGrid(PGrid)

    Call InitialClsvofFluidField(PGrid,PCell)
    Call InitialClsvofFluidField(UGrid,UCell)
    Call InitialClsvofFluidField(VGrid,VCell)
    Call InitialClsvofFluidField(WGrid,WCell)
    
    Call InitialClsvofLiquidField(PGrid,PCell)
    Call InitialClsvofLiquidField(UGrid,UCell)
    Call InitialClsvofLiquidField(VGrid,VCell)
    Call InitialClsvofLiquidField(WGrid,WCell)

    call BoundaryConditionLvsVof(PGrid, PCell, Var, BCLvs, BCVof, 0.d0)
    call BoundaryConditionLvsVofFluid(PGrid, PCell, Var, BCLvsF,BCVofF, 0.d0)

    Call InitialVar(vel,0.d0,0.d0,0.d0,300.d0,vel,300.d0,row,Lref, Var)
    !Call PrintResultTecplotPCent(PGrid,Var,PCell,INT8(0))
    !Call PrintResultTecplotUCent(UGrid,Var,UCell,INT8(0))
    !Call PrintResultTecplotVCent(VGrid,Var,VCell,INT8(0))
    !Call PrintResultTecplotWCent(WGrid,Var,WCell,INT8(0))
    !<per -nag
    !Call PrintResultVTK(PGrid,Var,PCell,INT8(0)) 
    !Call PrintResultVTR3D(PGrid,Var,PCell,"FlowField",INT8(0))
    !<per -nag

    Call GridPreProcess(int(1,it8b),PGrid,UGrid,VGrid,WGrid, PCell,UCell,VCell,WCell)

    Call DefineMomentumExchangeCell(PCell, UCell,VCell,WCell)
    
    Call NumberExternalCell(0,0,0, PCell)
    Call NumberExternalCell(1,0,0, UCell)
    Call NumberExternalCell(0,1,0, VCell)
    Call NumberExternalCell(0,0,1, WCell)

    Call NewCellFace(PGrid,UGrid,VGrid,WGrid, PCell,UCell,VCell,WCell)

    call BoundaryConditionVarNew(0.d0,PGrid, PCell, Var, BCp, BCu, BCv, BCw)

    Call IterationSolution(1,PGrid,UGrid,VGrid,WGrid, PCell,UCell,VCell,WCell,    &
                           BCu,BCv,BCw,BCp,BCVof,BCLvs,Var)
 
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
    allocate(UCell%vof( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%phi( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nx(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%ny(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nz(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nxL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nyL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(UCell%nzL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    
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
    allocate(VCell%vof( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%phi( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nx(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%ny(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nz(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nxL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nyL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(VCell%nzL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    
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
    Allocate(WCell%vof( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%phi( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%nx(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%ny(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(wCell%nz(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%nxL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(WCell%nyL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(wCell%nzL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    
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
    Allocate(PCell%phi( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nx(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%ny(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nz(  0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nxL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nyL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(PCell%nzL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
    
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
    Allocate(Var%Gpu(Imax,Jmax,Kmax))
    Allocate(Var%Gpv(Imax,Jmax,Kmax))
    Allocate(Var%Gpw(Imax,Jmax,Kmax))
    Allocate(Var%ures(Imax,Jmax,Kmax))
    Allocate(Var%vres(Imax,Jmax,Kmax))
    Allocate(Var%wres(Imax,Jmax,Kmax))
    Allocate(Var%pres(Imax,Jmax,Kmax))
    Allocate(Var%mres(Imax,Jmax,Kmax))
    Allocate(Var%t(0:Imax+1,0:Jmax+1,0:Kmax+1))
 End subroutine

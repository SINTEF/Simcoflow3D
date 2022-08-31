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
PROGRAM Main

    USE PrecisionVar
    USE StateVariables
    USE WaveGeneration
    USE ClsVof
    USE Mesh
    USE Cutcell
    USE Solver
    USE PrintResult
    USE InitialVof

    Implicit none
    
    TYPE(SolverTime) :: time
    TYPE(Grid) :: Ugrid, VGrid, WGrid, Pgrid
    TYPE(Cell) :: Ucell, VCell, WCell, PCell
    TYPE(Variables) :: Var
    TYPE(Point) :: Spoint,EPoint,ReS,ReE
    TYPE(TWaveGeneration) :: Wave 

    INTEGER(it4b) :: i,j,k 
    INTEGER(it4b) :: Irec,Jrec,Krec,iprint
    INTEGER(it4b) :: N,M
    INTEGER(it4b) :: c1, c2
    REAL(dp) :: Lref
    REAL(dp) :: length, width, U10, F, gammaVal, thetaW

    ! Read input file
    OPEN(UNIT=5,FILE='/home/elena-roxanap/Documents/Iceload/simco3d/src/Test/WavesSimulation/input.dat',ACTION='read')
    READ(5,*)
    READ(5,*) Imax, Jmax, Kmax, Irec, Jrec, Krec, Rey, Lref, iprint
    READ(5,*)
    READ(5,*) TimeOrder, SpaceOrder
    CLOSE(5)
   
    ! Read input file for waves simulation
    OPEN(UNIT=10, FILE='/home/elena-roxanap/Documents/Iceload/simco3d/src/Test/WavesSimulation/InputFileWavesSimulation.dat', ACTION='read')
    READ(10,*)
    READ(10,*)
    READ(10,*)
    READ(10,*)
    READ(10,*)
    READ(10,*) length,width
    READ(10,*) 
    READ(10,*) N,M
    READ(10,*) 
    READ(10,*) 
    READ(10,*) U10, F, gammaVal, thetaW
    F = F * 1000.d0 ! [m]
    CLOSE(10)
    c1 = 4
    c2 = 4
   
    !! WAVES
    ! constructor for the Waves Simulatiton
    ! allocation of arrays and initialization of parameters
    Wave = TWaveGeneration(N,M,c1,c2,U10,F,gammaVal,thetaW,length,width) 

    ! generate characteristics of waves
    call Wave%calcWaves(0.d0)

    ! allocate variables
    Imax = Wave%N
    Jmax = Wave%M
    Kmax = Jmax
    call AllocateVar(PGrid,UGrid,VGrid,Wgrid,Pcell,UCell,VCell,WCell,Var)


    !! GRID
    ! initialize grid
    !  domain length

    SPoint%x = Wave%x0(1)
    EPoint%x = Wave%x0(Wave%N)
    SPoint%y = Wave%y0(1)
    EPoint%y = Wave%y0(Wave%M)
    SPoint%z = -3.d0
    EPoint%z = 3.d0

    PGrid%dx(:,:,:) = (EPoint%x-Spoint%x)/dble(Imax-1)
    PGrid%dy(:,:,:) = (EPoint%y-Spoint%y)/dble(Jmax-1)
    PGrid%dz(:,:,:) = (EPoint%z-Spoint%z)/dble(Kmax-1)

    do i = 1,Imax
      do j = 1,Jmax
        do k = 1,Kmax
          PGrid%x(i,j,k)=Wave%x0(i)!SPoint%x+PGrid%dx(i,j,k)*(dble(i)-0.5d0)
          PGrid%y(i,j,k)=Wave%y0(j)!SPoint%y+PGrid%dy(i,j,k)*(dble(j)-0.5d0)
          PGrid%z(i,j,k)=SPoint%z+PGrid%dz(i,j,k)*(dble(k)-0.5d0)
        end do
      end do
    end do


    print*, minval(Pgrid%x), maxval(Pgrid%x)
    !do i = 1,Imax
    !  print*, 'wave grid', Wave%x0(i)
    !  print*, ' pgrid', PGrid%x(i,1,1)
    !enddo

    ! mpi and hypre intialization
    !call MPI_Initial

    !call HYPRE_CreateGrid(Pgrid)

    ! initiate level set and vof field
    Wave%time = 0.d0
    Wave%dt = 0.25d0*2.d0*pi/maxval(Wave%omega)
    print*, 'dt', Wave%dt
    i=0
    do while (Wave%time.lt.10.d0)
      call Wave%calcWaves(Wave%time)
      call InitialClsvofLiquidFieldWaves(PGrid,PCell,Wave)
      call PrintResultVTK(Pgrid,Var,Pcell,int8(i)) 
      Wave%time = Wave%time+Wave%dt
      i=i+1
    enddo


END PROGRAM Main

SUBROUTINE AllocateVar(PGrid,UGrid,VGrid,Wgrid,Pcell,UCell,VCell,WCell,Var)
  USE Mesh
  USE Cutcell
  USE StateVariables

  IMPLICIT NONE

  TYPE(Grid), INTENT(inout) :: Pgrid,UGrid,VGrid,WGrid
  TYPE(Cell), INTENT(inout) :: PCell,UCell,VCell,WCell
  TYPE(Variables), INTENT(inout) :: Var


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

  Allocate(UCell%Cell_Type(Imax,Jmax,Kmax))
  Allocate(UCell%vof( 0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%phi( 0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%nx(  0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%ny(  0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%nz(  0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%nxL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%nyL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(UCell%nzL( 0:Imax+1,0:Jmax+1,0:Kmax+1))

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

  Allocate(VCell%Cell_Type(Imax,Jmax,Kmax))
  Allocate(VCell%vof( 0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%phi( 0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%nx(  0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%ny(  0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%nz(  0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%vofL(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%phiL(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%nxL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%nyL( 0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%nzL( 0:Imax+1,0:Jmax+1,0:Kmax+1))

  Allocate(VCell%Cell_Cent(Imax,Jmax,Kmax,3))
  Allocate(VCell%EEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%NEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%TEArea(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%MoExCell(Imax,Jmax,Kmax))
  Allocate(VCell%EtaE(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%EtaN(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%EtaT(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%AlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%AlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%AlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%DAlE(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%DAlN(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%DAlT(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%SxE(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%SyN(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%SzT(0:Imax+1,0:Jmax+1,0:Kmax+1))
  Allocate(VCell%FCE(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
  Allocate(VCell%FCN(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
  Allocate(VCell%FCT(0:Imax+1,0:Jmax+1,0:Kmax+1,3))
  Allocate(VCell%PosNu(Imax,Jmax-1,Kmax))
  Allocate(VCell%WlLh(Imax,Jmax,Kmax))
  Allocate(VCell%delh(Imax,Jmax,Kmax))
  Allocate(VCell%MsCe(Imax,Jmax,Kmax,3))
  Allocate(VCell%WePr(Imax,Jmax,Kmax))

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

  END SUBROUTINE AllocateVar


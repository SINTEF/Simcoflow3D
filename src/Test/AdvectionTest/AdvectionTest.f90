program main
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE Clsvof
    USE StateVariables
    USE PrintResult
    USE MPI
    USE Solver
    USE BoundaryInterface
    USE BoundaryFunction
    USE ComputePUV
    USE InitialVof
    
    implicit none
    
    Type(Grid)      :: PGrid
    Type(Cell)      :: PCell
    Type(Point)     :: SPoint,EPoint,ReS,ReE
    Type(Variables) :: Var
    type(BCBase)    :: BCp,BCu,BCv,BCw,BCVof,BCLvs
    
    real(kind=dp)   :: dt,t,tp,eta,nv,cfl
    integer(kind=it4b) :: i,j,k,iprint,iprint1
    integer(kind=it8b) :: itt
    real(kind=dp), dimension(:), allocatable :: Constin
    
    allocate(Constin(6))
    
    open(unit=1,file='/home/sontd/code/CutCell3DGFMCLSVOF/src/Test/AdvectionTest/input.txt',status='old',action='read')
    read(1,*)
    read(1,*)
    read(1,*) imax,jmax,kmax,tp,iprint,iprint1,eta,nv,cfl
    
    SPoint%x = 0.d0
    SPoint%y = 0.d0
    SPoint%z = 0.d0
    EPoint%x = 1.d0
    EPoint%y = 1.d0
    EPoint%z = 1.d0
    
    Allocate(PGrid%x(Imax,Jmax,Kmax))
    Allocate(PGrid%y(Imax,Jmax,Kmax))
    Allocate(PGrid%z(Imax,Jmax,Kmax))
    
    Allocate(PGrid%dx(Imax,Jmax,Kmax))
    Allocate(PGrid%dy(Imax,Jmax,Kmax))
    Allocate(PGrid%dz(Imax,Jmax,Kmax))
    Allocate(PCell%vof(Imax,Jmax,Kmax))
    Allocate(PCell%phi(Imax,Jmax,Kmax))
    Allocate(PCell%nx(Imax,Jmax,Kmax))
    Allocate(PCell%ny(Imax,Jmax,Kmax))
    Allocate(PCell%nz(Imax,Jmax,Kmax))
    Allocate(PCell%vofL(Imax,Jmax,Kmax))
    Allocate(PCell%phiL(Imax,Jmax,Kmax))
    Allocate(PCell%nxL(Imax,Jmax,Kmax))
    Allocate(PCell%nyL(Imax,Jmax,Kmax))
    Allocate(PCell%nzL(Imax,Jmax,Kmax))
    Allocate(Var%u(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%v(0:Imax+1,0:Jmax+1,0:Kmax+1))
    Allocate(Var%w(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(Var%p(0:Imax+1,0:Jmax+1,0:Kmax+1))
    allocate(Var%mres(Imax,Jmax,Kmax))
    
    pgrid%dx(:,:,:) = (EPoint%x-SPoint%x)/dble(Imax)
    pgrid%dy(:,:,:) = (EPoint%y-SPoint%y)/dble(Jmax)
    pgrid%dz(:,:,:) = (EPoint%z-SPoint%z)/dble(Kmax)
  
    do i = 1,Imax
      do j = 1,Jmax
        do k = 1,Kmax
          pgrid%x(i,j,k)=SPoint%x+pGrid%dx(i,j,k)*(dble(i)-0.5d0)
          pgrid%y(i,j,k)=SPoint%y+pGrid%dy(i,j,k)*(dble(j)-0.5d0)
          pgrid%z(i,j,k)=SPoint%z+pGrid%dz(i,j,k)*(dble(k)-0.5d0)
        end do
      end do
    end do
    BCu = BCBase(Imax,Jmax,Kmax)
    BCv = BCBase(Imax,Jmax,Kmax)
    BCw = BCBase(Imax,Jmax,Kmax)
    BCvof = BCBase(Imax,Jmax,Kmax)
    BClvs = BCBase(Imax,Jmax,Kmax)
    
    ! Advection Test
    call BCu%SetDN(0,0,0,0,0,0)
    call BCv%SetDN(0,0,0,0,0,0)
    call BCw%SetDN(0,0,0,0,0,0)
    call BCVof%SetDN(0,0,0,0,0,0)
    call BCLvs%SetDN(0,0,0,0,0,0)
    ! Set Constant for boundary condition  
    Constin(:) = 0.d0
    call BCu%SetConstant(Constin)   
    call BCv%SetConstant(Constin)
    call BCw%SetConstant(Constin)
    call BCVof%SetConstant(Constin)
    call BCLvs%SetConstant(Constin)
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
    
    deallocate(Constin)
    call InitialClsvofFluidFieldAdvectionTest(PGrid,PCell)
    call InitialClsvofLiquidFieldAdvectionTest(PGrid,PCell)
    call PrintResultVTR3D(PGrid,Var,PCell,'AdvectionLvsVof',INT8(0))
    print*, 'After print result at initial stage'
    do itt = 1,INT8(100000)
      dt = 1.d0
      do i = 1,imax
        do j = 1,jmax
          do k = 1,kmax
            Var%u(i,j,k)=2.d0*(dsin(pi*PGrid%x(i,j,k)))**2.d0*dsin(2.d0*pi*    &
                      PGrid%y(i,j,k))*dsin(2.d0*pi*PGrid%z(i,j,k))*dcos(pi*t/tp) 
            Var%v(i,j,k)=-dsin(2.d0*pi*PGrid%x(i,j,k))*(dsin(pi*               &
               PGrid%y(i,j,k)))**2.d0*dsin(2.d0*pi*PGrid%z(i,j,k))*dcos(pi*t/tp)
            Var%w(i,j,k)=-dsin(2.d0*pi*PGrid%x(i,j,k))*dsin(2.d0*pi*           &
               PGrid%y(i,j,k))*(dsin(pi*PGrid%z(i,j,k)))**2.d0*dcos(pi*t/tp)
            dt=dmin1(dt,cfl*PGrid%dx(1,1,1)/dabs(Var%u(i,j,k)+1.d-30),         &
                     cfl*PGrid%dy(1,1,1)/dabs(Var%v(i,j,k)+1.d-30),  	         &
                     cfl*PGrid%dz(1,1,1)/dabs(Var%w(i,j,k)+1.d-30))
             end do
          end do
       end do  
       call BoundaryConditionVarNew(PGrid, PCell, Var, BCp, BCu, BCv, BCw, t)
       if(dt>=cfl*PGrid%dx(1,1,1)/0.2d0) dt = cfl*PGrid%dx(1,1,1)/0.2d0
       if(t<tp.and.(t+dt)>=tp) then
         pause ' do you want to continue '
         dt = tp - t
         call Clsvof_Scheme(PGrid,PCell,Var,BCu,BCv,BCw,BCLvs,BCvof,t,dt,itt)
         iprint = iprint1
       else
         call Clsvof_Scheme(PGrid,PCell,Var,BCu,BCv,BCw,BCLvs,BCvof,t,dt,itt)
       end if
       t = t+dt
       if(mod(itt,iprint)==0) then
         call PrintResultVTR3D(PGrid,Var,PCell,'AdvectionLvsVof',itt)
     !     call Print_Result_3d(Grid,Sta,Inte)
     !     call Print_Result_2D_YZ(Grid,Sta,Inte,88)
         print*,'Iteration number and time:', itt , t
       end if
       if(t>=tp/2.d0.and.t-dt<tp/2.d0) then
         call PrintResultVTR3D(PGrid,Var,PCell,'AdvectionLvsVof',itt)
     !     call Print_Result_3d(Grid,Sta,Inte)
     !     call Print_Result_2D_YZ(Grid,Sta,Inte,88)
         print*,'Iteration number, time and time step:0', itt,t,dt,t-dt
       end if 
       if(t>1.2d0*tp) exit
    end do       
    deallocate(Var%u,Var%v,Var%w)
    deallocate(PGrid%x,PGrid%y,PGrid%z)
    deallocate(PGrid%dx,PGrid%dy,PGrid%dz)
    deallocate(PCell%phi)
    deallocate(PCell%vof)
    deallocate(PCell%nx)
    deallocate(PCell%ny)
    deallocate(PCell%nz)
    deallocate(PCell%vofL)
    deallocate(PCell%phiL)
    deallocate(PCell%nxL)
    deallocate(PCell%nyL)
    deallocate(PCell%nzL)
end program

Module ComputePUVW
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE PredictorUVW
    USE ProjectionP
    use BoundaryInterface
    use BoundaryFunction
    use MPI

    Implicit none
    Private
    Public:: UpdatePUVW,BoundaryConditionVarNew

    Interface UpdatePUVW
      Module Procedure UpdatePUVW
    End interface

    Interface BoundaryConditionVarNew
      module procedure BoundaryConditionVarNew
    End interface

    Contains
    Subroutine UpdatePUVW(UGrid, VGrid, WGrid, PGrid,                           &
                         UCell, VCell, WCell, PCell,                           &
                         UCellO, VCellO, WCellO, PCellO,                       &
                         BCu, BCv, BCw, BCp, BCVof, BCLvs,                     &
                         FluxDivOld, TVar_n, TVar, Time, dt, itt)
      !! The subroutine is used to update the velocity and pressure
      Implicit none
      Type(Grid),      intent(in)    :: UGrid, VGrid, WGrid, PGrid
      !! The input grid
      Type(Cell),      intent(in)    :: UCell, VCell, WCell, PCell
      !! The present cell configuration 
      Type(Cell),      intent(in)    :: UCellO, VCellO, WCellO, PCellO
      !! The previous cell configuration
      type(BCBase),    intent(inout) :: BCu, BCv, BCw, BCp, BCVof, BCLvs
      !! The boundary parameter
      type(Variables), intent(inout) :: TVar_n
      !! The state variables at n-1
      type(Variables), intent(inout) :: TVar
      !! The state variables at n for 'in' (intent(inout)) and n+1 for 'out' (intent(inout))
      real(kind=dp),dimension(:,:,:,:),allocatable,intent(inout) :: FluxDivOld
      !! The previous flux difference
      real(kind=dp),     intent(in)  :: Time
      !! The time for computing mpi
      Real(dp),          intent(in)  :: dt
      !! The time step size
      Integer(kind=it8b),intent(in)  :: itt
      !! The iteration number
      Type(PoissonCoefficient)       :: PU, PV, PW
      Type(Predictor)                :: Pred
      Type(Projection)               :: Proj
      Integer(kind=it4b) 	           :: i,j,k,ii,jj,kk,iu,iv,iw
      Real(kind=dp)                  :: dps,ute,GradPUVW,maxPoCoef
      real(kind=dp),dimension(:,:,:,:),allocatable :: PoCoef

      allocate(Pred%u(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(Pred%v(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(Pred%w(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(Proj%Pp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(PU%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(PV%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(PW%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(PoCoef(Imax,JMax,KMax,6)) ! The order of the face : W-1,S-2,B-3,E-4,N-5,T-6
      ! Set the predictor velocity
      Pred%u(:,:,:) = TVar%u(:,:,:)
      Pred%v(:,:,:) = TVar%v(:,:,:)
      Pred%w(:,:,:) = TVar%w(:,:,:)
      Proj%Pp(:,:,:) = 0.d0
      ! The predictor step to compute the predicting velocities
      call PredictingUVW(PGrid, UGrid, VGrid, WGrid,                           &
                        PCell, UCell, VCell, WCell,                            &
                        PCellO, UCellO, VCellO, WCellO,                        &
                        FluxDivOld, TVar_n, TVar,                              &
                        BCu, BCv, BCw, PU, PV, PW, Pred, dt, itt)
      ! Solve the Poisson equation
      call PoissonEquationSolver(PGrid, UGrid, VGrid, WGrid,                   &
                                PCell, UCell, VCell, WCell,                    &
                                TVar, Pred, PU, PV, PW, BCp, PoCoef, Proj, dt)
      ! Set the velocity at time step n-1
      if(itt>=2) then
        TVar_n%u(:,:,:)=TVar%u(:,:,:)
        TVar_n%v(:,:,:)=TVar%v(:,:,:)
        TVar_n%w(:,:,:)=TVar%w(:,:,:)
        TVar_n%p(:,:,:)=Tvar%p(:,:,:)
      end if  
      maxPoCoef=0.d0
      do i = 1,Imax
        do j = 1,Jmax
          do k = 1,Kmax
            if(UCell%MoExCell(i,j,k)/=1.and.UCell%Cell_Type(i,j,k)/=2) then
              GradPUVW=(Proj%Pp(i+1,j,k)-Proj%Pp(i,j,k))*PoCoef(i,j,k,4)
              TVar%u(i,j,k)=Pred%u(i,j,k)-GradPUVW
            else
              ii=UCell%MsCe(i,j,k,1)
              jj=UCell%MsCe(i,j,k,2)
              kk=UCell%MsCe(i,j,k,3)
              GradPUVW=(Proj%Pp(ii+1,jj,kk)-Proj%Pp(ii,jj,kk))*PoCoef(i,j,k,4)
           !   TVar%u(i,j,k)=Pred%u(i,j,k)-GradPUVW
              TVar%u(i,j,k)=0.d0
            end if
            if(VCell%MoExCell(i,j,k)/=1.and.VCell%Cell_Type(i,j,k)/=2) then
              GradPUVW=(Proj%Pp(i,j+1,k)-Proj%Pp(i,j,k))*PoCoef(i,j,k,5)
              TVar%v(i,j,k)=Pred%v(i,j,k)-GradPUVW
            else
              ii = VCell%MsCe(i,j,k,1)
              jj = VCell%MsCe(i,j,k,2)
              kk = VCell%MsCe(i,j,k,3)
              GradPUVW=(Proj%Pp(ii,jj+1,kk)-Proj%Pp(ii,jj,kk))*PoCoef(i,j,k,5)
           !   TVar%v(i,j,k)=Pred%v(i,j,k)-GradPUVW
              TVar%v(i,j,k)=0.d0
            end if
            if(WCell%MoExCell(i,j,k)/=1.and.WCell%Cell_Type(i,j,k)/=2) then
              GradPUVW=(Proj%Pp(i,j,k+1)-Proj%Pp(i,j,k))*PoCoef(i,j,k,6)
              TVar%w(i,j,k)=Pred%w(i,j,k)-GradPUVW
            else
              ii = WCell%MsCe(i,j,k,1)
              jj = WCell%MsCe(i,j,k,2)
              kk = WCell%MsCe(i,j,k,3)
              GradPUVW=(Proj%Pp(ii,jj,kk+1)-Proj%Pp(ii,jj,kk))*PoCoef(i,j,k,6)
           !   TVar%w(i,j,k)=Pred%w(i,j,k)-GradPUVW
              TVar%w(i,j,k)=0.d0
            end if
            TVar%p(i,j,k)=Proj%Pp(i,j,k)
            if(dabs(Tvar%u(i,j,k))>MaxPocoef) then
              MaxPocoef=dabs(TVar%u(i,j,k))
              iu=i
              iv=j
              iw=k
            endif 
          end do
        end do
      end do
  !    print*,PoCoef(iu,iv,iw,:)
  !    print*,'Check predicted velocity'
  !    print*,Pred%u(iu,iv,iw),Pred%u(iu-1,iv,iw)
  !    print*,Pred%v(iu,iv,iw),Pred%v(iu,iv-1,iw)
  !    print*,Pred%w(iu,iv,iw),Pred%w(iu,iv,iw-1)
  !    print*,'================================='
  !    call BoundaryConditionVar(TVar)
      call BoundaryConditionVarNew(PGrid, PCell, TVar, BCp, BCu, BCv, BCw, Time)
      ! Correcting velocity at the domain boundary to assure the global mass conservation
      ! At the western boundary
      if(BCp%flag(1)==0) then
        i = 1
        do j = 1,Jmax
          do k = 1,Kmax
            TVar%u(i-1,j,k)=(PGrid%dx(i,j,k)*PGrid%dz(i,j,k)*                  &
                            (PCell%NEArea(i,j,k)*TVar%v(i,j,k)-                &
                             PCell%NEArea(i,j,k)*TVar%v(i,j-1,k))+             &
                             PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*                  &
                            (PCell%TEArea(i,j,k)*TVar%w(i,j,k)-                &
                             PCell%TEArea(i,j,k)*TVar%w(i,j,k-1))+             &
                             PGrid%dy(i,j,k)*PGrid%dz(i,j,k)*                  &
                             PCell%EEArea(i,j,k)*TVar%u(i,j,k))/               &
                             PCell%EEArea(i,j,k)/PGrid%dy(i,j,k)/PGrid%dz(i,j,k)
            if(isnan(TVar%u(i-1,j,k))) then
              print*,PCell%EEArea(i,j,k)
              pause 'very small face East fraction computepuvw 125'
            end if  
          end do
        end do
      ! At the eastern boundary  
      elseif(BCp%flag(2)==0) then
        i = Imax
        do j = 1,Jmax
          do k = 1,Kmax
            TVar%u(i,j,k)=(PGrid%dx(i,j,k)*PGrid%dz(i,j,k)*                    &
                         (-PCell%NEArea(i,j,k)*TVar%v(i,j,k)+                  &
                           PCell%NEArea(i,j,k)*TVar%v(i,j-1,k))+               &
                           PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*                    &
                         (-PCell%TEArea(i,j,k)*TVar%w(i,j,k)+                  &
                           PCell%TEArea(i,j,k)*TVar%w(i,j,k-1))+               &
                           PGrid%dy(i,j,k)*PGrid%dz(i,j,k)*                    &
                           PCell%EEArea(i,j,k)*TVar%u(i-1,j,k))/               &
                           PCell%EEArea(i,j,k)/PGrid%dy(i,j,k)/PGrid%dz(i,j,k)
            if(isnan(TVar%u(i,j,k))) then
              print*,PCell%EEArea(i,j,k)
              pause 'very small face East fraction computepuvw 144'
            end if 
          end do
        end do
      ! At the southern boundary  
      elseif(BCp%flag(3)==0) then
        j=1
        do i=1,Imax
          do k=1,Kmax
            TVar%v(i,j-1,k)=(PGrid%dy(i,j,k)*PGrid%dz(i,j,k)*                  &
                           (PCell%EEArea(i,j,k)*TVar%u(i,j,k)-                 &
                            PCell%EEArea(i,j,k)*TVar%u(i-1,j,k))+              &
                            PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*                   &
                           (PCell%TEArea(i,j,k)*TVar%w(i,j,k)-                 &
                            PCell%TEArea(i,j,k)*TVar%w(i,j,k-1))+              &
                            PGrid%dx(i,j,k)*PGrid%dz(i,j,k)*                   &
                            PCell%NEArea(i,j,k)*TVar%v(i,j,k))/                &
                            PCell%NEArea(i,j,k)/PGrid%dx(i,j,k)/PGrid%dz(i,j,k)
            if(isnan(TVar%v(i,j-1,k))) then
              print*,PCell%NEArea(i,j,k)
              pause 'very small face North fraction computepuvw 163'
            end if                 
          end do
        end do
      ! At the northern boundary  
      elseif(BCp%flag(4)==0) then
        j=Jmax
        do i=1,Imax
          do j=1,Jmax
            TVar%v(i,j,k) =(-PGrid%dy(i,j,k)*PGrid%dz(i,j,k)*                  &
                           (PCell%EEArea(i,j,k)*TVar%u(i,j,k)-                 &
                            PCell%EEArea(i,j,k)*TVar%u(i-1,j,k))               &
                            -PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*                  &
                           (PCell%TEArea(i,j,k)*TVar%w(i,j,k)-                 &
                            PCell%TEArea(i,j,k)*TVar%w(i,j,k-1))+              &
                            PGrid%dx(i,j,k)*PGrid%dz(i,j,k)*                   &
                            PCell%NEArea(i,j,k)*TVar%v(i,j-1,k))/              &
                            PCell%NEArea(i,j,k)/PGrid%dx(i,j,k)/PGrid%dz(i,j,k)
            if(isnan(TVar%v(i,j,k))) then
              print*,PCell%NEArea(i,j,k)
              pause 'very small face North fraction computepuvw 182'
            end if
          end do
        end do
      ! At the bottom boundary  
      elseif(BCp%flag(5)==0) then
        k=1
        do i=1,Imax
          do j=1,Jmax
            TVar%w(i,j,k-1)=(PGrid%dy(i,j,k)*PGrid%dz(i,j,k)*                  &
                           (PCell%EEArea(i,j,k)*TVar%u(i,j,k)-                 &
                            PCell%EEArea(i,j,k)*TVar%u(i-1,j,k))+              &
                            PGrid%dx(i,j,k)*PGrid%dz(i,j,k)*                   &
                           (PCell%NEArea(i,j,k)*TVar%v(i,j,k)-                 &
                            PCell%NEArea(i,j,k)*TVar%v(i,j-1,k))+              &
                            PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*                   &
                            PCell%TEArea(i,j,k)*TVar%w(i,j,k))/                &
                            PCell%TEArea(i,j,k)/PGrid%dx(i,j,k)/PGrid%dy(i,j,k)
            if(isnan(TVar%w(i,j,k-1))) then
              print*,PCell%TEArea(i,j,k)
              pause 'very small face Top fraction computepuvw 201'
            end if
          end do
        end do
      ! At the top boundary  
      elseif(BCp%flag(6)==0) then
        k=kmax
        do i=1,Imax
          do j=1,Jmax
            TVar%w(i,j,k)=(-PGrid%dy(i,j,k)*PGrid%dz(i,j,k)*                   &
                          (PCell%EEArea(i,j,k)*TVar%u(i,j,k)-                  &
                           PCell%EEArea(i,j,k)*TVar%u(i-1,j,k))                &
                           -PGrid%dx(i,j,k)*PGrid%dz(i,j,k)*                   &
                          (PCell%NEArea(i,j,k)*TVar%v(i,j,k)-                  &
                           PCell%NEArea(i,j,k)*TVar%v(i,j-1,k))+               &
                           PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*                    &
                           PCell%TEArea(i,j,k)*TVar%w(i,j,k-1))/               &
                           PCell%TEArea(i,j,k)/PGrid%dx(i,j,k)/PGrid%dy(i,j,k)
            if(isnan(TVar%w(i,j,k))) then
              print*,PCell%TEArea(i,j,k)
              pause 'very small face Top fraction computepuvw 220'
            end if
          end do
        end do
      end if
      
      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax
            if(Pred%v(i,j,k)>MaxPoCoef) then
              MaxPocoef = dabs(TVar%v(i,j,k))
              ii = i
              jj = j
              kk = k
            end if
          end do
        end do
      end do     
      print*, 'Set up velocity v to 0 ComputePUVW.f90 278'
      TVar%v=0.d0
      deallocate(Pred%u)
      deallocate(Pred%v)
      deallocate(Pred%w)
      deallocate(Proj%Pp)
      deallocate(PU%Dp)
      deallocate(PV%Dp)
      deallocate(PW%Dp)
      deallocate(PoCoef)
    end subroutine UpdatePUVW

    subroutine BoundaryConditionVarNew(PGrid, PCell, Vari, BCp, BCu, BCv, BCw, Time)
      type(Grid),      intent(in)    :: PGrid
      type(Cell),      intent(in)    :: PCell
      type(Variables), intent(inout) :: Vari
      type(BCBase),    intent(inout) :: BCu, BCv, BCw, BCp
      real(kind=dp),   intent(in)    :: Time
      integer(kind=it4b)             :: i,j,k

      ! For the western boundary
      call BCu%west(PGrid%x(1,:,:)-PGrid%dx(1,:,:)/2.d0, PGrid%y(1,:,:),       &
                    PGrid%z(1,:,:), PGrid%dx(1,:,:), PGrid%dy(1,:,:), 	       &
                    PGrid%dz(1,:,:), Vari%p(1,:,:), Vari%u(1,:,:), 	       &
                    Vari%v(1,:,:), Vari%w(1,:,:), PCell%vof(1,:,:), 	       &
                    PCell%phi(1,:,:), Time)
      Vari%u(1-ight,1:Jmax,1:Kmax) = BCu%VarW(1:Jmax,1:Kmax)

      ! For the eastern boundary
      call BCu%east(PGrid%x(Imax,:,:)+PGrid%dx(Imax,:,:)/2.d0, PGrid%y(Imax,:,:), &
                    PGrid%z(Imax,:,:), PGrid%dx(Imax,:,:), PGrid%dy(Imax,:,:), &
                    PGrid%dz(Imax,:,:), Vari%p(Imax,:,:), Vari%u(Imax-1,:,:),  &
                    Vari%v(Imax,:,:), Vari%w(Imax,:,:), PCell%vof(Imax,:,:),   &
                    PCell%phi(Imax,:,:), Time)
      Vari%u(Imax,1:Jmax,1:Kmax) = BCu%VarE(1:Jmax,1:Kmax)
      call ComputeGhostVarBoundary(Vari%u(Imax,1:Jmax,1:Kmax),BCu%VarE,BCu%flag(2),      &
                                   Vari%u(Imax+ight,1:Jmax,1:Kmax))

      ! For the southern boundary
      call BCu%South(PGrid%x(:,1,:)+PGrid%dx(:,1,:)/2.d0, 		       &
                     PGrid%y(:,1,:)-PGrid%dy(:,1,:)/2.d0, PGrid%z(:,1,:),      &
                     PGrid%dx(:,1,:), PGrid%dy(:,1,:), PGrid%dz(:,1,:),        &
                     Vari%p(:,1,:), Vari%u(:,1,:), Vari%v(:,1,:), 	           &
                     Vari%w(:,1,:), PCell%vof(:,1,:), PCell%phi(:,1,:), Time)
      call ComputeGhostVarBoundary(Vari%u(1:Imax,1,1:Kmax),BCu%VarS,BCu%flag(3),         &
                                                 Vari%u(1:Imax,1-jght,1:Kmax))

      ! For the northern boundary
      call BCu%North(PGrid%x(:,Jmax,:)+PGrid%dx(:,Jmax,:)/2.d0, 	             &
                     PGrid%y(:,Jmax,:)+PGrid%dy(:,Jmax,:)/2.d0, 	             &
                     PGrid%z(:,Jmax,:), PGrid%dx(:,Jmax,:), PGrid%dy(:,Jmax,:),&
                     PGrid%dz(:,Jmax,:), Vari%p(:,Jmax,:), Vari%u(:,Jmax,:),   &
                     Vari%v(:,Jmax,:), Vari%w(:,Jmax,:), PCell%vof(:,Jmax,:),  &
                     PCell%phi(:,Jmax,:), Time)
      call ComputeGhostVarBoundary(Vari%u(1:Imax,Jmax,1:Kmax),BCu%VarN,        &
                                   BCu%flag(4),Vari%u(1:Imax,Jmax+jght,1:Kmax))

      ! For the bottom boundary
      call BCu%Bottom(PGrid%x(:,:,1)+PGrid%dx(:,:,1)/2.d0, PGrid%y(:,:,1),     &
                      PGrid%z(:,:,1)-PGrid%dz(:,:,1)/2.d0, PGrid%dx(:,:,1),    &
                      PGrid%dy(:,:,1), PGrid%dz(:,:,1), Vari%p(:,:,1),         &
                      Vari%u(:,:,1), Vari%v(:,:,1), Vari%w(:,:,1),             &
                      PCell%vof(:,:,1), PCell%phi(:,:,1), Time)
      call ComputeGhostVarBoundary(Vari%u(1:Imax,1:Jmax,1),BCu%VarB,BCu%flag(5),        &
                                   Vari%u(1:Imax,1:Jmax,1-kght))

      ! For the top boundary
      call BCu%Top(PGrid%x(:,:,Kmax)+PGrid%dx(:,:,Kmax)/2.d0,                  &
                  PGrid%y(:,:,Kmax), PGrid%z(:,:,Kmax)+PGrid%dz(:,:,Kmax)/2.d0,&
                  PGrid%dx(:,:,Kmax), PGrid%dy(:,:,Kmax), PGrid%dz(:,:,Kmax),  &
                  Vari%p(:,:,Kmax), Vari%u(:,:,Kmax), Vari%v(:,:,Kmax),        &
                  Vari%w(:,:,Kmax), PCell%vof(:,:,Kmax), PCell%phi(:,:,Kmax), Time)
      call ComputeGhostVarBoundary(Vari%u(1:Imax,1:Jmax,Kmax),BCu%VarT,BCu%flag(6),      &
                                   Vari%u(1:Imax,1:Jmax,Kmax+kght))

      ! For the western boundary
      call BCv%west(PGrid%x(1,:,:)-PGrid%dx(1,:,:)/2.d0,                       &
                    PGrid%y(1,:,:)+PGrid%dy(1,:,:)/2.d0,                       &
                    PGrid%z(1,:,:), PGrid%dx(1,:,:), PGrid%dy(1,:,:),          &
                    PGrid%dz(1,:,:), Vari%p(1,:,:), Vari%u(1,:,:),             &
                    Vari%v(1,:,:), Vari%w(1,:,:), PCell%vof(1,:,:),            &
                    PCell%phi(1,:,:), Time)
      call ComputeGhostVarBoundary(Vari%v(1,1:Jmax,1:Kmax),BCv%VarW,BCv%flag(1),         &
                                   Vari%v(1-ight,1:Jmax,1:Kmax))

      ! For the eastern boundary
      call BCv%east(PGrid%x(Imax,:,:)+PGrid%dx(Imax,:,:)/2.d0, 		       &
                    PGrid%y(Imax,:,:)+PGrid%dy(Imax,:,:)/2.d0, 		       &
                    PGrid%z(Imax,:,:), PGrid%dx(Imax,:,:), PGrid%dy(Imax,:,:), &
                    PGrid%dz(Imax,:,:), Vari%p(Imax,:,:), Vari%u(Imax-1,:,:),  &
                    Vari%v(Imax,:,:), Vari%w(Imax,:,:), PCell%vof(Imax,:,:),   &
                    PCell%phi(Imax,:,:), Time)
      call ComputeGhostVarBoundary(Vari%v(Imax,1:Jmax,1:Kmax),BCv%VarE,        &
                                   BCv%flag(2),Vari%v(Imax+ight,1:Jmax,1:Kmax))

      ! For the southern boundary
      call BCv%South(PGrid%x(:,1,:), PGrid%y(:,1,:)-PGrid%dy(:,1,:)/2.d0,      &
                     PGrid%z(:,1,:), PGrid%dx(:,1,:), PGrid%dy(:,1,:), 	       &
                     PGrid%dz(:,1,:), Vari%p(:,1,:), Vari%u(:,1,:), 	         &
                     Vari%v(:,1,:), Vari%w(:,1,:), PCell%vof(:,1,:), 	         &
                     PCell%phi(:,1,:), Time)
      Vari%v(1:Imax,1-jght,1:Kmax) = BCv%VarS(1:Imax,1:Kmax)   	    
      ! For the northern boundary
      call BCv%North(PGrid%x(:,Jmax,:), PGrid%y(:,Jmax,:)+PGrid%dy(:,Jmax,:)/2.d0, &
                     PGrid%z(:,Jmax,:), PGrid%dx(:,Jmax,:), PGrid%dy(:,Jmax,:),&
                     PGrid%dz(:,Jmax,:), Vari%p(:,Jmax,:), Vari%u(:,Jmax,:),   &
                     Vari%v(:,Jmax-1,:), Vari%w(:,Jmax,:), PCell%vof(:,Jmax,:),&
                     PCell%phi(:,Jmax,:), Time)
      Vari%v(1:Imax,Jmax,1:Kmax) = BCv%VarN(1:Imax,1:Kmax)
      call ComputeGhostVarBoundary(Vari%v(1:Imax,Jmax,1:Kmax),BCv%VarN,        &
                                   BCv%flag(4),Vari%v(1:Imax,Jmax+jght,1:Kmax))
      ! For the bottom boundary
      call BCv%Bottom(PGrid%x(:,:,1), PGrid%y(:,:,1)+PGrid%dy(:,:,1)/2.d0,     &
                      PGrid%z(:,:,1)-PGrid%dz(:,:,1)/2.d0, PGrid%dx(:,:,1),    &
                      PGrid%dy(:,:,1), PGrid%dz(:,:,1), Vari%p(:,:,1), 	       &
                      Vari%u(:,:,1), Vari%v(:,:,1), Vari%w(:,:,1), 	           &
                      PCell%vof(:,:,1), PCell%phi(:,:,1), Time)
      call ComputeGhostVarBoundary(Vari%v(1:Imax,1:Jmax,1),BCv%VarB,BCv%flag(5),&
                                   Vari%v(1:Imax,1:Jmax,1-kght))
      ! For the top boundary
      call BCv%Top(PGrid%x(:,:,Kmax), PGrid%y(:,:,Kmax)+PGrid%dy(:,:,Kmax)/2.d0,&
                   PGrid%z(:,:,Kmax)+PGrid%dz(:,:,Kmax)/2.d0,                   &
                   PGrid%dx(:,:,Kmax), PGrid%dy(:,:,Kmax), PGrid%dz(:,:,Kmax),  &
                   Vari%p(:,:,Kmax), Vari%u(:,:,Kmax), Vari%v(:,:,Kmax),        &
                   Vari%w(:,:,Kmax), PCell%vof(:,:,Kmax), PCell%phi(:,:,Kmax), Time)
      call ComputeGhostVarBoundary(Vari%v(1:Imax,1:Jmax,Kmax),BCv%VarT,        &
                                   BCv%flag(6),Vari%v(1:Imax,1:Jmax,Kmax+kght))

      ! For the western boundary
      call BCw%west(PGrid%x(1,:,:)-PGrid%dx(1,:,:)/2.d0, PGrid%y(1,:,:),       &
                    PGrid%z(1,:,:)+PGrid%dz(1,:,:)/2.d0, PGrid%dx(1,:,:),      &
                    PGrid%dy(1,:,:), PGrid%dz(1,:,:), Vari%p(1,:,:), 	         &
                    Vari%u(1,:,:), Vari%v(1,:,:), Vari%w(1,:,:), 	             &
                    PCell%vof(1,:,:), PCell%phi(1,:,:), Time)
      call ComputeGhostVarBoundary(Vari%w(1,1:Jmax,1:Kmax),BCw%VarW,           &
                                   BCw%flag(1),Vari%w(1-ight,1:Jmax,1:Kmax))
      ! For the eastern boundary
      call BCw%east(PGrid%x(Imax,:,:)+PGrid%dx(Imax,:,:)/2.d0, 		             &
                    PGrid%y(Imax,:,:),PGrid%z(Imax,:,:)+PGrid%dz(Imax,:,:)/2.d0,&
                    PGrid%dx(Imax,:,:), PGrid%dy(Imax,:,:), 		               &
                    PGrid%dz(Imax,:,:), Vari%p(Imax,:,:), Vari%u(Imax-1,:,:),  &
                    Vari%v(Imax,:,:), Vari%w(Imax,:,:), PCell%vof(Imax,:,:),   &
                    PCell%phi(Imax,:,:), Time)
      call ComputeGhostVarBoundary(Vari%w(Imax,1:Jmax,1:Kmax),BCw%VarE,BCw%flag(2),      &
                                                 Vari%w(Imax+ight,1:Jmax,1:Kmax))
      ! For the southern boundary
      call BCw%South(PGrid%x(:,1,:), PGrid%y(:,1,:)-PGrid%dy(:,1,:)/2.d0,      &
                     PGrid%z(:,1,:)+PGrid%dz(:,1,:)/2.d0, PGrid%dx(:,1,:),     &
                     PGrid%dy(:,1,:), PGrid%dz(:,1,:), Vari%p(:,1,:),          &
                     Vari%u(:,1,:), Vari%v(:,1,:), Vari%w(:,1,:),              &
                     PCell%vof(:,1,:), PCell%phi(:,1,:), Time)
      call ComputeGhostVarBoundary(Vari%w(1:Imax,1,1:Kmax),BCw%VarS,BCw%flag(3),         &
                                                 Vari%w(1:Imax,1-jght,1:Kmax))
      ! For the northern boundary
      call BCw%North(PGrid%x(:,Jmax,:), PGrid%y(:,Jmax,:)+PGrid%dy(:,Jmax,:)/2.d0, &
                     PGrid%z(:,Jmax,:)+PGrid%dz(:,Jmax,:)/2.d0,                & 
                     PGrid%dx(:,Jmax,:), PGrid%dy(:,Jmax,:),                   &
                     PGrid%dz(:,Jmax,:), Vari%p(:,Jmax,:), Vari%u(:,Jmax,:),   &
                     Vari%v(:,Jmax,:), Vari%w(:,Jmax,:), PCell%vof(:,Jmax,:),  &
                     PCell%phi(:,Jmax,:), Time)
      call ComputeGhostVarBoundary(Vari%w(1:Imax,Jmax,1:Kmax),BCw%VarN,BCw%flag(4),      &
                                   Vari%w(1:Imax,Jmax+jght,1:Kmax))
      ! For the bottom boundary
      call BCw%Bottom(PGrid%x(:,:,1), PGrid%y(:,:,1),                          &
                      PGrid%z(:,:,1)-PGrid%dz(:,:,1)/2.d0, PGrid%dx(:,:,1),    &
                      PGrid%dy(:,:,1), PGrid%dz(:,:,1), Vari%p(:,:,1),         &
                      Vari%u(:,:,1), Vari%v(:,:,1), Vari%w(:,:,1),             &
                      PCell%vof(:,:,1), PCell%phi(:,:,1), Time)
      Vari%w(1:Imax,1:Jmax,1-kght) = BCw%VarB(1:Imax,1:Jmax)
      ! For the top boundary
      call BCw%Top(PGrid%x(:,:,Kmax), PGrid%y(:,:,Kmax),                       &
                   PGrid%z(:,:,Kmax)+PGrid%dz(:,:,Kmax)/2.d0,                  &
                   PGrid%dx(:,:,Kmax), PGrid%dy(:,:,Kmax), PGrid%dz(:,:,Kmax), &
                   Vari%p(:,:,Kmax), Vari%u(:,:,Kmax), Vari%v(:,:,Kmax),       &
                   Vari%w(:,:,Kmax-1), PCell%vof(:,:,Kmax), PCell%phi(:,:,Kmax), Time)
      Vari%w(1:Imax,1:Jmax,Kmax)=BCw%VarT(1:Imax,1:Jmax)
      call ComputeGhostVarBoundary(Vari%w(1:Imax,1:Jmax,Kmax),BCw%VarT,        &
                                   BCw%flag(6),Vari%w(1:Imax,1:Jmax,Kmax+kght))
      ! For the western boundary
      call BCp%west(PGrid%x(1,:,:)-PGrid%dx(1,:,:)/2.d0, PGrid%y(1,:,:),       &
                    PGrid%z(1,:,:), PGrid%dx(1,:,:),                           &
                    PGrid%dy(1,:,:), PGrid%dz(1,:,:), Vari%p(1,:,:),           &
                    Vari%u(1,:,:), Vari%v(1,:,:), Vari%w(1,:,:),               &
                    PCell%vof(1,:,:), PCell%phi(1,:,:), Time)
      call ComputeGhostVarBoundary(Vari%p(1,1:Jmax,1:Kmax),BCp%VarW,BCp%flag(1),&
                                                 Vari%p(1-ight,1:Jmax,1:Kmax))
      ! For the eastern boundary
      call BCp%east(PGrid%x(Imax,:,:)+PGrid%dx(Imax,:,:)/2.d0,                 &
                   PGrid%y(Imax,:,:), PGrid%z(Imax,:,:),                       &
                   PGrid%dx(Imax,:,:), PGrid%dy(Imax,:,:),                     &
                   PGrid%dz(Imax,:,:), Vari%p(Imax,:,:), Vari%u(Imax-1,:,:),   &
                   Vari%v(Imax,:,:), Vari%w(Imax,:,:), PCell%vof(Imax,:,:),    &
                   PCell%phi(Imax,:,:), Time)
      call ComputeGhostVarBoundary(Vari%p(Imax,1:Jmax,1:Kmax),BCp%VarE,        &
                                   BCp%flag(2),Vari%p(Imax+ight,1:Jmax,1:Kmax))
      ! For the southern boundary
      call BCp%South(PGrid%x(:,1,:), PGrid%y(:,1,:)-PGrid%dy(:,1,:)/2.d0,      &
                     PGrid%z(:,1,:), PGrid%dx(:,1,:),                          &
                     PGrid%dy(:,1,:), PGrid%dz(:,1,:), Vari%p(:,1,:),          &
                     Vari%u(:,1,:), Vari%v(:,1,:), Vari%w(:,1,:),              &
                     PCell%vof(:,1,:), PCell%phi(:,1,:), Time)
      call ComputeGhostVarBoundary(Vari%p(1:Imax,1,1:Kmax),BCp%VarS,           &
                                   BCp%flag(3),Vari%p(1:Imax,1-jght,1:Kmax))
      ! For the northern boundary
      call BCp%North(PGrid%x(:,Jmax,:), PGrid%y(:,Jmax,:)+PGrid%dy(:,Jmax,:)/2.d0, &
                     PGrid%z(:,Jmax,:), PGrid%dx(:,Jmax,:), PGrid%dy(:,Jmax,:),&
                     PGrid%dz(:,Jmax,:), Vari%p(:,Jmax,:), Vari%u(:,Jmax,:),   &
                     Vari%v(:,Jmax,:), Vari%w(:,Jmax,:), PCell%vof(:,Jmax,:),  &
                     PCell%phi(:,Jmax,:), Time)
      call ComputeGhostVarBoundary(Vari%p(1:Imax,Jmax,1:Kmax),BCp%VarN,        &
                                    BCp%flag(4),Vari%p(1:Imax,Jmax+jght,1:Kmax))
      ! For the bottom boundary
      call BCp%Bottom(PGrid%x(:,:,1), PGrid%y(:,:,1),                          &
                      PGrid%z(:,:,1)-PGrid%dz(:,:,1)/2.d0, PGrid%dx(:,:,1),    &
                      PGrid%dy(:,:,1), PGrid%dz(:,:,1), Vari%p(:,:,1),         &
                      Vari%u(:,:,1), Vari%v(:,:,1), Vari%w(:,:,1),             &
                      PCell%vof(:,:,1), PCell%phi(:,:,1), Time)
      call ComputeGhostVarBoundary(Vari%p(1:Imax,1:Jmax,1),BCp%VarB,           &
                      BCp%flag(5),Vari%p(1:Imax,1:Jmax,1-ight))
      ! For the top boundary
      call BCp%Top(PGrid%x(:,:,Kmax), PGrid%y(:,:,Kmax),                       &
                   PGrid%z(:,:,Kmax)+PGrid%dz(:,:,Kmax)/2.d0,                  &
                   PGrid%dx(:,:,Kmax), PGrid%dy(:,:,Kmax), PGrid%dz(:,:,Kmax), &
                   Vari%p(:,:,Kmax), Vari%u(:,:,Kmax), Vari%v(:,:,Kmax),       &
            Vari%w(:,:,Kmax-1), PCell%vof(:,:,Kmax), PCell%phi(:,:,Kmax), Time)
      Vari%p(1:Imax,1:Jmax,Kmax)=BCp%VarT(1:Imax,1:Jmax)
      call ComputeGhostVarBoundary(Vari%p(1:Imax,1:Jmax,Kmax),BCp%VarT,        &
                                   BCp%flag(6),Vari%p(1:Imax,1:Jmax,Kmax+kght))                                                                               
    end subroutine BoundaryConditionVarNew

    subroutine ComputeGhostVarBoundary(Varin,BCin,flag,Varout)
      real(kind=dp), dimension(:,:), intent(in)  :: Varin,BCin
      real(kind=dp), dimension(:,:), intent(out) :: Varout
      integer(kind=it4b), intent(in) 		 :: flag

      if(flag==1) then
        Varout(:,:) = BCin(:,:)-Varin
      else
        Varout(:,:) = 2.d0*BCin(:,:)-Varin(:,:)
      end if
    end subroutine ComputeGhostVarBoundary

    subroutine VariablesInternalCellCondition(TVar,PCell,UCell,VCell,WCell)
      implicit none
      type(Variables),intent(inout):: TVar
      type(Cell),intent(in):: PCell,UCell,VCell,WCell
      integer(kind=it4b):: i,j,k

      do i = 1,Imax
        do j = 1,Jmax
          do k = 1,Kmax
            if(PCell%Cell_Type(i,j,k)==2) then
              TVar%p(i,j,k) = 0.d0
            end if
            if(UCell%Cell_Type(i,j,k)==2) then
              TVar%u(i,j,k) = 0.d0
            end if
            if(VCell%Cell_Type(i,j,k)==2) then
              TVar%v(i,j,k) = 0.d0
            end if
          end do
        end do
      end do
    end subroutine VariablesInternalCellCondition
end module ComputePUVW

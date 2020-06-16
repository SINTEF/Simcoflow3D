Module ComputePUV
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE PredictorUV
    USE ProjectionP
    Implicit none
    Private
    Public:: UpdatePUV
    Interface UpdatePUV
      Module Procedure UpdatePUV
    End interface
    Contains
    Subroutine UpdatePUV(UGrid,VGrid,WGrid,PGrid,UCell,VCell,WCell,PCell,      &
                                                       TVar_n,TVar,dt,itt)
      Implicit none
      Type(Grid),intent(in)         :: UGrid,VGrid,WGrid,PGrid
      Type(Cell),intent(inout)      :: UCell,VCell,WCell,PCell
      Type(Variables),intent(in)    :: TVar_n
      Type(Variables),intent(inout) :: TVar
      Real(dp),intent(in)           :: dt
      Integer(kind=it8b),intent(in) :: itt
      Type(PoissonCoefficient)      :: PU,PV,PW
      Type(Predictor)		    :: Pred
      Type(Projection)              :: Proj
      Integer(kind=it4b) 	    :: i,j,k,ii,jj,kk
      Real(kind=dp)                 :: dps,ute,GradPUVW,maxPoCoef
      real(kind=dp),dimension(:,:,:,:),allocatable :: PoCoef 
      
      allocate(Pred%u(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(Pred%v(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(Pred%w(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(Proj%Pp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(PU%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(PV%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(PW%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      allocate(PoCoef(Imax,JMax,KMax,6)) ! The order of the face : W-1,S-2,B-3,E-4,N-5,T-6
      
      Pred%u(:,:,:) = TVar%u(:,:,:)
      Pred%v(:,:,:) = TVar%v(:,:,:)
      Pred%w(:,:,:) = TVar%w(:,:,:)
      Proj%Pp(:,:,:) = 0.d0
      
      call PredictorUVW(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,             &
                        WCell,TVar_n,TVar,PU,PV,PW,Pred,dt)
      call PoissonEquationSolver(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,    &
                        WCell,TVar,Pred,PU,PV,PW,PoCoef,Proj,dt)
      
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
          end do
        end do
      end do
  !    print*, 'Compute the Poisson Coefficient'        
  !    print*, MaxPoCoef
  !    print*, Pred%u(ii,jj,kk)
  !    print*, '------------------------------------'
  !    print*, TVar%u(63,41,30),TVar%u(64,41,30)
  !    print*, TVar%v(63,41,30),TVar%v(64,41,30)
  !    print*, TVar%w(63,41,30),TVar%w(64,41,30)
  !    print*, ii,jj,kk
  !    print*,'Print out the maximum coefficient Poisson Coefficient'
  !    print*,MaxPoCoef
  !    print*,PoCoef(ii,jj,kk,4)
  !    print*,PoCoef(ii,jj,kk,1)
  !    print*,PCell%vof(ii,jj,kk),UCell%vof(ii,jj,kk)
  !    print*,'test 222222222222222222222222'
  !    print*,PU%dp(59,48,48)/UGrid%dx(59,48,48)
  !    print*,UCell%MoExCell(59,48,48)	
  !    print*,UCell%vof(59,48,48)
  !    print*,ii,jj,kk
  !    print*,'1111111111111111'
      call BoundaryConditionVar(TVar)
      i = Imax
      do j = 1,Jmax
        do k = 1,Kmax
          TVar%u(i,j,k) = (PGrid%dx(Imax,j,k)*PGrid%dz(Imax,j,k)*              &
                PCell%NEArea(Imax,j,k)*(-TVar%v(Imax,j,k)+TVar%v(Imax,j-1,k))+ &
                              PGrid%dx(Imax,j,k)*PGrid%dy(Imax,j,k)*           &
                PCell%TEArea(Imax,j,k)*(-TVar%w(Imax,j,k)+TVar%w(Imax,j,k-1))+ &
                              PGrid%dy(Imax,j,k)*PGrid%dz(Imax,j,k)*           &
                PCell%EEArea(Imax,j,k)*TVar%u(Imax-1,j,k))/                    &
                PCell%EEArea(Imax,j,k)/PGrid%dy(Imax,j,k)/PGrid%dz(Imax,j,k)
        end do
      end do
      
      deallocate(Pred%u)
      deallocate(Pred%v)
      deallocate(Pred%w)
      deallocate(Proj%Pp)
      deallocate(PU%Dp)
      deallocate(PV%Dp)
      deallocate(PW%Dp)
      deallocate(PoCoef)
    end subroutine UpdatePUV

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
end module ComputePUV

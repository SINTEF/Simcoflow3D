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
      Type(Grid),intent(in):: UGrid,VGrid,WGrid,PGrid
      Type(Cell),intent(inout):: UCell,VCell,WCell,PCell
      Type(Variables),intent(in):: TVar_n
      Type(Variables),intent(inout):: TVar
      Real(dp),intent(in):: dt
      Integer(kind=it8b),intent(in):: itt
      Type(PoissonCoefficient):: PU,PV,PW
      Type(Predictor):: Pred
      Type(Projection):: Proj
      Integer(kind=it4b):: i,j,k,ii,jj,kk
      Real(kind=dp):: dps,ute
      Allocate(Pred%u(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      Allocate(Pred%v(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      Allocate(Pred%w(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      Allocate(Proj%Pp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      Allocate(PU%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      Allocate(PV%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      Allocate(PW%Dp(1-ight:Imax+ight,1-jght:Jmax+jght,1-kght:Kmax+kght))
      Pred%u(:,:,:) = TVar%u(:,:,:)
      Pred%v(:,:,:) = TVar%v(:,:,:)
      Pred%w(:,:,:) = TVar%w(:,:,:)
   !   TVar%p(:,:,:) = 0.d0
      Proj%Pp(:,:,:) = 0.d0 !TVar%p(:,:)
      Call PredictorUVW(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,             &
                        WCell,TVar_n,TVar,PU,PV,PW,Pred,dt)
      Call PoissonEquationSolver(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,    &
                        WCell,TVar,Pred,PU,PV,PW,Proj,dt)
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            If(UCell%MoExCell(i,j,k)/=1.and.UCell%Cell_Type(i,j,k)/=2) then
              TVar%u(i,j,k) = Pred%u(i,j,k)-(Proj%Pp(i+1,j,k)-Proj%Pp(i,j,k))* &
                                                                PU%Dp(i,j,k)
              TVar%Gpu(i,j,k)=(TVar%p(i+1,j,k)+Proj%Pp(i+1,j,k)-Proj%Pp(i,j,k)-&
                                               TVar%p(i,j,k))/UGrid%dx(i,j,k)
            Elseif(UCell%MoExCell(i,j,k)==1) then
              ii = UCell%MsCe(i,j,k,1)
              jj = UCell%MsCe(i,j,k,2)
              kk = UCell%MsCe(i,j,k,3)
              dps = (Proj%Pp(ii+1,jj,kk)-Proj%Pp(ii,jj,kk))
              TVar%u(i,j,k) = Pred%u(i,j,k)-dps*PU%Dp(i,j,k)
              TVar%Gpu(i,j,k) = (Proj%Pp(ii+1,jj,kk)+TVar%p(ii+1,jj,kk)-       &
                            Proj%Pp(ii,jj,kk)-TVar%p(ii,jj,kk))/UGrid%dx(i,j,k)
            End if
            If(VCell%MoExCell(i,j,k)/=1.and.VCell%Cell_Type(i,j,k)/=2) then
              TVar%v(i,j,k) = Pred%v(i,j,k)-(Proj%Pp(i,j+1,k)-Proj%Pp(i,j,k))* &
                                                                  PV%Dp(i,j,k)
              TVar%Gpv(i,j,k)=(TVar%p(i,j+1,k)+Proj%Pp(i,j+1,k)-Proj%Pp(i,j,k)-&
                                                  TVar%p(i,j,k))/VGrid%dy(i,j,k)
            ElseIf(VCell%MoExCell(i,j,k)==1) then
              ii = VCell%MsCe(i,j,k,1)
              jj = VCell%MsCe(i,j,k,2)
              kk = VCell%MsCe(i,j,k,3)
              dps = (Proj%Pp(ii,jj+1,kk)-Proj%Pp(ii,jj,kk))
              TVar%v(i,j,k) = Pred%v(i,j,k)-dps*PV%Dp(i,j,k)
              TVar%Gpv(i,j,k) = (Proj%Pp(ii,jj+1,kk)+TVar%p(ii,jj+1,kk)-       &
                               Proj%Pp(ii,jj,kk)-TVar%p(ii,jj,kk))
            End if
            If(WCell%MoExCell(i,j,k)/=1.and.WCell%Cell_Type(i,j,k)/=2) then
              TVar%w(i,j,k) = Pred%w(i,j,k)-(Proj%Pp(i,j,k+1)-Proj%Pp(i,j,k))* &
                                                                  PW%Dp(i,j,k)
              TVar%Gpw(i,j,k)=(TVar%p(i,j,k+1)+Proj%Pp(i,j,k+1)-Proj%Pp(i,j,k)-&
                                                  TVar%p(i,j,k))/WGrid%dz(i,j,k)
            Elseif(WCell%MoExCell(i,j,k)==1) then
              ii = WCell%MsCe(i,j,k,1)
              jj = WCell%MsCe(i,j,k,2)
              kk = WCell%MsCe(i,j,k,3)
              dps = (Proj%Pp(ii,jj,kk+1)-Proj%Pp(ii,jj,kk))
              TVar%w(i,j,k) = Pred%w(i,j,k)-dps*PW%Dp(i,j,k)
              TVar%Gpw(i,j,k) = (Proj%Pp(ii,jj,kk+1)+TVar%p(ii,jj,kk+1)-       &
                                 Proj%Pp(ii,jj,kk)-TVar%p(ii,jj,kk))
            End if
            TVar%p(i,j,k) = TVar%p(i,j,k)+Proj%Pp(i,j,k)
          End do
        End do
      End do
      Call BoundaryConditionVar(TVar)
      i = Imax
      Do j = 1,Jmax
        Do k = 1,Kmax
          TVar%u(i,j,k) = (PGrid%dx(Imax,j,k)*PGrid%dz(Imax,j,k)*              &
                PCell%NEArea(Imax,j,k)*(-TVar%v(Imax,j,k)+TVar%v(Imax,j-1,k))+ &
                              PGrid%dx(Imax,j,k)*PGrid%dy(Imax,j,k)*           &
                PCell%TEArea(Imax,j,k)*(-TVar%w(Imax,j,k)+TVar%w(Imax,j,k-1))+ &
                              PGrid%dy(Imax,j,k)*PGrid%dz(Imax,j,k)*           &
                PCell%EEArea(Imax,j,k)*TVar%u(Imax-1,j,k))/                    &
                PCell%EEArea(Imax,j,k)/PGrid%dy(Imax,j,k)/PGrid%dz(Imax,j,k)
        End do
      End do
      Deallocate(Pred%u)
      Deallocate(Pred%v)
      Deallocate(Pred%w)
      Deallocate(Proj%Pp)
      Deallocate(PU%Dp)
      Deallocate(PV%Dp)
      Deallocate(PW%Dp)
    End Subroutine UpdatePUV

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
            End if
            If(UCell%Cell_Type(i,j,k)==2) then
              TVar%u(i,j,k) = 0.d0
            End if
            If(VCell%Cell_Type(i,j,k)==2) then
              TVar%v(i,j,k) = 0.d0
            End if
          End do
        End do
      End do
    End subroutine VariablesInternalCellCondition
End Module ComputePUV

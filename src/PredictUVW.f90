Module PredictorUV
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE Matrix
    USE Printresult
    USE MPI
    Implicit none
    Private
    Real(kind=dp),dimension(:,:,:),pointer:: u
    Real(kind=dp),dimension(:,:,:),pointer:: v
    Real(kind=dp),dimension(:,:,:),pointer:: w
    Real(dp),dimension(:,:,:),pointer:: p
    Real(dp),parameter:: alpha=0.d0,beta = 1.d0
    Type,public:: Predictor
        Real(kind=dp),dimension(:,:,:),allocatable:: u,v,w
    End Type
    Type,public:: PoissonCoefficient
        Real(kind=dp),dimension(:,:,:),allocatable:: Dp
    End Type
    Public:: PredictorUVW
    Interface PredictorUVW
        Module Procedure PredictorUVW
    End interface
    Contains
    Subroutine PredictorUVW(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,  &
                                              TVar_n,TVar,PU,PV,PW,Pred,dt)
      Implicit none
      Real(kind=dp),intent(in):: dt
      Type(Grid),intent(in):: PGrid,UGrid,VGrid,WGrid
      Type(Cell),intent(inout):: PCell,UCell,VCell,WCell
      Type(Variables),intent(in),target:: TVar
      Type(Variables),intent(in):: TVar_n
      Type(Predictor),intent(inout):: Pred
      Type(PoissonCoefficient),intent(inout):: PU,PV,PW
      Type(Variables):: TVart
      Integer(kind=it4b):: i,j,k,ii,jj,kk
      Integer*8:: A,parcsr_A,b,par_b,x,par_x,solver,precond

      ! Convective flux and Diffusive flux
      Real(kind=dp),dimension(:,:,:,:),allocatable:: CFEW,CFNS,CFTB,           &
                                     DFEW,DFNS,DFTB,EDFEW,EDFNS,EDFTB,FluxDiv
      Real(kind=dp),dimension(:,:,:),allocatable:: UFric,VFric,WFric,          &
                                           UWE,USN,UBT,VWE,VSN,VBT,WWE,WSN,WBT
      Integer(kind=it4b):: num_iterations
      Real(kind=dp):: final_res_norm,MaxU
      Real(kind=dp):: Fe,Fw,Fn,Fs,Ft,Fb
      Allocate(CFEW(0:Imax+1,Jmax,Kmax,3))
      Allocate(CFNS(Imax,0:Jmax+1,Kmax,3))
      Allocate(CFTB(Imax,Jmax,0:Kmax+1,3))
      Allocate(DFEW(0:Imax+1,Jmax,Kmax,3))
      Allocate(DFNS(Imax,0:Jmax+1,Kmax,3))
      Allocate(DFTB(Imax,Jmax,0:Kmax+1,3))
      Allocate(EDFEW(0:Imax+1,Jmax,Kmax,3))
      Allocate(EDFNS(Imax,0:Jmax+1,Kmax,3))
      Allocate(EDFTB(Imax,Jmax,0:Kmax+1,3))
      Allocate(FluxDiv(Imax,Jmax,Kmax,3))
      Allocate(UFric(Imax,Jmax,Kmax))
      Allocate(VFric(Imax,Jmax,Kmax))
      Allocate(WFric(Imax,Jmax,Kmax))
      Allocate(UWE(Jmax,Kmax,2))
      Allocate(USN(Imax,Kmax,2))
      Allocate(UBT(Imax,Jmax,2))
      Allocate(VWE(Jmax,Kmax,2))
      Allocate(VSN(Imax,Kmax,2))
      Allocate(VBT(Imax,Jmax,2))
      Allocate(WWE(Jmax,Kmax,2))
      Allocate(WSN(Imax,Kmax,2))
      Allocate(WBT(Imax,Jmax,2))
      u => TVar%u
      v => TVar%v
      w => TVar%w
      p => TVar%p

    ! Step 1: Calculate the convective coefficient

      Call ModifiedConvectiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,   &
                                        WCell,CFEW,1,0,0)
      Call ModifiedConvectiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,   &
                                        WCell,CFNS,0,1,0)
      Call ModifiedConvectiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,   &
                                        WCell,CFTB,0,0,1)
    ! Step 2: Calculate the diffusive coefficient
      EDFEW = 0.d0
      EDFNS = 0.d0
      EDFTB = 0.d0
      FluxDiv = 0.d0
      Call DiffusiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,      &
                                         DFEW,EDFEW,1,0,0)
      Call DiffusiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,      &
                                         DFNS,EDFNS,0,1,0)
      Call DiffusiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,      &
                                         DFTB,EDFTB,0,0,1)
    ! Step 3: Calculate source term coefficient as wall function
      do i = 1,Imax-1
        do j = 1,Jmax
          do k = 1,Kmax
          ! U Cell
            if(UCell%Cell_Type(i,j,k)/=2) then
              Fe = CFEW(i+1,j,k,1)
              Fw = CFEW(i,j,k,1)
              Fn = CFNS(i,j+1,k,1)
              Fs = CFNS(i,j,k,1)
              Ft = CFTB(i,j,k+1,1)
              Fb = CFTB(i,j,k,1)
              FluxDiv(i,j,k,1)=(Fe-Fw+Fn-Fs+Ft-Fb)/UCell%vof(i,j,k)
            end if
            if(UCell%Cell_Type(i,j,k)/=2) then
              UFric(i,j,k)=(UCell%vofL(i,j,k)/UCell%vof(i,j,k)*nuw/nuref+      &
                     (1.d0-UCell%vofL(i,j,k)/UCell%vof(i,j,k))*nua/nuref)/Rey* &
                     UCell%WlLh(i,j,k)/UCell%delh(i,j,k)/UCell%vof(i,j,k)
            else
              UFric(i,j,k) = 0.d0
            end if
            if(i==68.and.j==46.and.k==46) then
              print*, 'Test the flux for computing the flux difference'
              print*, Fe,Fw
              print*, Fn,Fs
              print*, Ft,Fb
              print*, FluxDiv(i,j,k,1),UCell%vof(i,j,k)
              print*, (Fe-Fw+Fn-Fs+Ft-Fb)/UCell%vof(i,j,k)
              print*, '+++++++++++++++++++++++++++++++++++++++++++++++'
            end if
            if(isnan(Fluxdiv(i,j,k,1)).or.dabs(Fluxdiv(i,j,k,1))>1.d5) then
              print*, i,j,k,Pred%u(i,j,k)
              print*, fluxdiv(i,j,k,1),i,j,k
              print*, fe,fw,fn,fs,ft,fb
              pause 'predictoruv 114'
            end if
          end do
        end do
      end do
      print*, 'Flux after computing the totall flux'
      print*, FluxDiv(68,48,46,1)
      Do i = 1,Imax
        Do j = 1,Jmax-1
          Do k = 1,Kmax
            If(VCell%Cell_Type(i,j,k)/=2) then ! for V Cell
              Fe = CFEW(i+1,j,k,2)
              Fw = CFEW(i,j,k,2)
              Fn = CFNS(i,j+1,k,2)
              Fs = CFNS(i,j,k,2)
              Ft = CFTB(i,j,k+1,2)
              Fb = CFTB(i,j,k,2)
              FluxDiv(i,j,k,2)=(Fe-Fw+Fn-Fs+Ft-Fb)/VCell%vof(i,j,k)
            End if
          ! V Cell
            If(VCell%Cell_Type(i,j,k)/=2) then
              VFric(i,j,k)=(VCell%vofL(i,j,k)/VCell%vof(i,j,k)*nuw/nuref+      &
                     (1.d0-VCell%vofL(i,j,k)/VCell%vof(i,j,k))*nua/nuref)/Rey* &
                     VCell%WlLh(i,j,k)/VCell%delh(i,j,k)/VCell%vof(i,j,k)
            Else
              VFric(i,j,k)=0.d0
            End if
            if(i==66.and.j==48.and.k==46) then
              print*, 'Test the flux for computing the flux difference'
              print*, Fe,Fw
              print*, Fn,Fs
              print*, Ft,Fb
              print*, FluxDiv(i,j,k,2),UCell%vof(i,j,k)
              print*, (Fe-Fw+Fn-Fs+Ft-Fb)/UCell%vof(i,j,k)
              print*, '+++++++++++++++++++++++++++++++++++++++++++++++'
            end if
            If(isnan(Fluxdiv(i,j,k,2)).or.dabs(Fluxdiv(i,j,k,2))>1.d5) then
              print*,i,j,k,pred%v(i,j,k)
              print*,Fe,Fw,Fn,Fs,Ft,Fb
              pause 'predictoruv 143'
            End if
          End do
        End do
      End do
      print*, 'Flux after computing the total flux'
      print*, FluxDiv(68,46,46,2)
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax-1
            If(WCell%Cell_Type(i,j,k)/=2) then ! for WCell
              Fe = CFEW(i+1,j,k,3)
              Fw = CFEW(i,j,k,3)
              Fn = CFNS(i,j+1,k,3)
              Fs = CFNS(i,j,k,3)
              Ft = CFTB(i,j,k+1,3)
              Fb = CFTB(i,j,k,3)
              FluxDiv(i,j,k,3)=(Fe-Fw+Fn-Fs+Ft-Fb)/WCell%vof(i,j,k)
            End if
          ! W Cell
            If(WCell%Cell_Type(i,j,k)/=2) then
              WFric(i,j,k)=(WCell%vofL(i,j,k)/WCell%vof(i,j,k)*nuw/nuref+      &
                   (1.d0-WCell%vofL(i,j,k)/WCell%vof(i,j,k))*nua/nuref)/Rey*   &
                     WCell%WlLh(i,j,k)/WCell%delh(i,j,k)/WCell%vof(i,j,k)
            Else
              WFric(i,j,k) = 0.d0
            End if
            if(i==68.and.j==46.and.k==33) then
              print*, 'Test the flux for computing the flux difference'
              print*, Fe,Fw
              print*, Fn,Fs
              print*, Ft,Fb
              print*, FluxDiv(i,j,k,3)
              print*, (Fe-Fw+Fn-Fs+Ft-Fb)/UCell%vof(i,j,k)
              print*, '+++++++++++++++++++++++++++++++++++++++++++++++'
            end if
            If(isnan(Fluxdiv(i,j,k,3)).or.dabs(Fluxdiv(i,j,k,3))>1.d5) then
              print*,i,j,k,pred%w(i,j,k)
              print*,Fe,Fw,Fn,Fs,Ft,Fb
              pause 'predictoruv 167'
            End if
          End do
        End do
      End do
      
      
    !  print*,Fluxdiv(1,20,20,1),UGrid%dz(1,20,20)
    !  print*,Fluxdiv(1,20,20,1)/UGrid%dz(1,20,20)
    !  print*,
    ! Solving for UCell
      Call SetBasicSolver(solver,precond)
    ! call SetBasicSolver(solver=solver,ierr=ierr)
      Call SetMatrix(A,parcsr_A,UGrid,UCell,DFEW,DFNS,DFTB,EDFEW,EDFNS,EDFTB,  &
                                            UFric,PU,UWE,USN,UBT,dt,1,0,0)
      Call SetVectors(b,x,par_b,par_x,UGrid,UCell,PU,UWE,USN,UBT,              &
                                            FluxDiv(:,:,:,1),TVar_n,dt,1,0,0)
      Call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
      Call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    ! Run info - needed logging turned on
      Call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      Call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      Call DeltaGetValues(x,UCell,Pred%u,1,0,0)
      Call HYPRE_IJMatrixDestroy(A,ierr)
      Call HYPRE_IJVectorDestroy(b,ierr)
      Call HYPRE_IJVectorDestroy(x,ierr)
      Call HYPRE_BoomerAMGDestroy(precond,ierr)
      Call HYPRE_ParCSRPCGDestroy(solver,ierr)
      print*, 'After solving equations'
      print*, Pred%u(59,48,48)
    ! For VCell
      Call SetBasicSolver(solver,precond)
      Call SetMatrix(A,parcsr_A,VGrid,VCell,DFEW,DFNS,DFTB,EDFEW,EDFNS,EDFTB,  &
                                            VFric,PV,VWE,VSN,VBT,dt,0,1,0)
      Call SetVectors(b,x,par_b,par_x,VGrid,VCell,PV,VWE,VSN,VBT,              &
                                            FluxDiv(:,:,:,2),TVar_n,dt,0,1,0)
      Call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
      Call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    ! Run info - needed logging turned on
      Call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      Call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      Call DeltaGetValues(x,VCell,Pred%v,0,1,0)
      Call HYPRE_IJMatrixDestroy(A,ierr)
      Call HYPRE_IJVectorDestroy(b,ierr)
      Call HYPRE_IJVectorDestroy(x,ierr)
      Call HYPRE_BoomerAMGDestroy(precond,ierr)
      Call HYPRE_ParCSRPCGDestroy(solver,ierr)
      print*, Pred%v(59,48,48)
    ! For WCell
      Call SetBasicSolver(solver,precond)
      Call SetMatrix(A,parcsr_A,WGrid,WCell,DFEW,DFNS,DFTB,EDFEW,EDFNS,EDFTB,  &
                                            WFric,PW,WWE,WSN,WBT,dt,0,0,1)
      Call SetVectors(b,x,par_b,par_x,WGrid,WCell,PW,WWE,WSN,WBT,              &
                                            FluxDiv(:,:,:,3),TVar_n,dt,0,0,1)
      Call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
      Call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    ! Run info - needed logging turned on
      Call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      Call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      Call DeltaGetValues(x,WCell,Pred%w,0,0,1)
      Call HYPRE_IJMatrixDestroy(A,ierr)
      Call HYPRE_IJVectorDestroy(b,ierr)
      Call HYPRE_IJVectorDestroy(x,ierr)
      Call HYPRE_BoomerAMGDestroy(precond,ierr)
      Call HYPRE_ParCSRPCGDestroy(solver,ierr)

      PU%Dp(1-ight,:,:) = PU%Dp(1,:,:)
      PU%Dp(Imax,:,:) = PU%Dp(Imax-1,:,:)
      PU%Dp(Imax+ight,:,:) = PU%Dp(Imax,:,:)
      PU%Dp(:,1-jght,:) = PU%Dp(:,1,:)
      PU%Dp(:,Jmax+jght,:) = PU%Dp(:,Jmax,:)
      PU%Dp(:,:,1-kght) = PU%Dp(:,:,1)
      PU%Dp(:,:,Kmax+kght) = PU%Dp(:,:,Kmax)

      PV%Dp(1-ight,:,:) = PV%Dp(1,:,:)
      PV%Dp(Imax+ight,:,:) = PV%Dp(Imax,:,:)
      PV%Dp(:,1-jght,:) = PV%Dp(:,1,:)
      PV%Dp(:,Jmax,:) = PV%Dp(:,Jmax-1,:)
      PV%Dp(:,Jmax+jght,:) = PV%Dp(:,Jmax,:)
      PV%Dp(:,:,1-kght) = PV%Dp(:,:,1)
      PV%Dp(:,:,Kmax+kght) = PV%Dp(:,:,Kmax)

      PW%Dp(1-ight,:,:) = PW%Dp(1,:,:)
      PW%Dp(Imax+ight,:,:) = PW%Dp(Imax,:,:)
      PW%Dp(:,1-jght,:) = PW%Dp(:,1,:)
      PW%Dp(:,Jmax+jght,:) = PW%Dp(:,Jmax,:)
      PW%Dp(:,:,1-kght) = PW%Dp(:,:,1)
      PW%Dp(:,:,Kmax) = PW%Dp(:,:,Kmax-1)
      PW%Dp(:,:,Kmax+kght) = PW%Dp(:,:,Kmax)
      Call PredictorVelocityBoundaryCondition(Pred,TVar)
      Call PredictorVelocityInternalCellCondition(Pred,UCell,VCell,WCell)
      Deallocate(CFEW,CFNS,CFTB)
      Deallocate(DFEW,DFNS,DFTB)
      Deallocate(EDFEW,EDFNS,EDFTB)
      Deallocate(UFric,VFric,WFric)
      Deallocate(UWE,USN,UBT)
      Deallocate(VWE,VSN,VBT)
      Deallocate(WWE,WSN,WBT)
      Deallocate(FluxDiv)
      Nullify(u,v,w,p)
    End subroutine PredictorUVW
! calculate normal flux for both normal cell and cut cell
    Subroutine SetBasicSolver(solver,precond)
        Implicit none
        Integer*8,intent(inout):: solver
        Integer*8,intent(inout),optional:: precond
!       Set up and use a solver
        Call HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD,solver,ierr)
!       Set some parameters
        Call HYPRE_ParCSRPCGSetMaxIter(solver,50,ierr)
        Call HYPRE_ParCSRPCGSetTol(solver,1.0d-10,ierr)
        Call HYPRE_ParCSRPCGSetTwoNorm(solver,1,ierr)
        Call HYPRE_ParCSRPCGSetPrintLevel(solver,2,ierr)
        Call HYPRE_ParCSRPCGSetLogging(solver,1,ierr)
!       Now set up the AMG preconditioner and specify any parameters
        If(present(precond)) then
          Call HYPRE_BoomerAMGCreate(precond,ierr)
!        Set some parameters
!        Print less solver info since a preconditioner
!          Call HYPRE_BoomerAMGSetPrintLevel(precond,1,ierr);
!        Falgout coarsening
          Call HYPRE_BoomerAMGSetCoarsenType(precond,6,ierr)
!        SYMMETRIC G-S/Jacobi hybrid relaxation
          Call HYPRE_BoomerAMGSetRelaxType(precond,6,ierr)
!        Sweeeps on each level
          Call HYPRE_BoomerAMGSetNumSweeps(precond,1,ierr)
!        conv. tolerance
          Call HYPRE_BoomerAMGSetTol(precond,0.0d0,ierr)
!        do only one iteration!
          Call HYPRE_BoomerAMGSetMaxIter(precond,1,ierr)
!        set amg as the pcg preconditioner
!         precond_id = 2
          Call HYPRE_ParCSRPCGSetPrecond(solver,2,precond,ierr)
        End if
    End subroutine SetBasicSolver

    Subroutine SetMatrix(A,parcsr_A,TGrid,TCell,DFEW,DFNS,DFTB,EDFEW,EDFNS,    &
                                       EDFTB,Fric,PUVW,CWE,CSN,CBT,dt,iu,iv,iw)
        Implicit none
        Integer*8,intent(inout):: A,parcsr_A
        Type(Grid),intent(in):: TGrid
        Type(Cell),intent(in):: TCell
        Real(kind=dp),intent(in):: dt
        Real(kind=dp),dimension(:,:,:,:),allocatable,intent(in):: DFEW,DFNS,   &
                                                     DFTB,EDFEW,EDFNS,EDFTB
        Real(kind=dp),dimension(:,:,:),allocatable,intent(in):: Fric
        Type(PoissonCoefficient),intent(inout):: PUVW
        Real(kind=dp),dimension(:,:,:),allocatable,intent(inout):: CWE,CSN,CBT
        Integer,intent(in):: iu,iv,iw
        Integer(kind=it4b):: nnz,ictr,ilower,iupper,cols(0:6)
        Integer(kind=it4b):: i,j,k,ii,jj,kk
        Real(kind=dp):: Fep,Fem,Fwp,Fwm,Fnp,Fnm,Fsp,Fsm,Ftp,Ftm,Fbp,Fbm
        Real(kind=dp):: aP,aE,aW,aN,aS,aT,aB,Dfe,Dfw,Dfn,Dfs,Dft,Dfb,Sp
        Real(kind=dp):: values(0:6),MaxDiff
        ilower = 0
        iupper = TCell%ExtCell
      ! Create and Set up matrix
        Call HYPRE_IJMatrixCreate(MPI_COMM_WORLD,ilower,iupper,ilower,iupper,  &
                                                                        A,ierr)
        Call HYPRE_IJMatrixSetObjectType(A,HYPRE_PARCSR,ierr)
        Call HYPRE_IJMatrixInitialize(A,ierr)
        MaxDiff=0.d0
        do i = 1,Imax-iu
          do j = 1,Jmax-iv
            do k = 1,Kmax-iw
              PUVW%Dp(i,j,k)=dt ! It is used to avoid the bug causing by vof=0
              ! If(i==1) j = Jmax-iv-j
              if(TCell%Cell_Type(i,j,k)/=2) then
                Dfe=0.d0;Dfw=0.d0;Dfn=0.d0;Dfs=0.d0;Dft=0.d0;Dfb=0.d0
                Fep=0.d0;Fem=0.d0;Fwp=0.d0;Fwm=0.d0;Fnp=0.d0;Fnm=0.d0
                Fsp=0.d0;Fsm=0.d0;Ftp=0.d0;Ftm=0.d0;Fbp=0.d0;Fbm=0.d0
                
              ! For UCell,VCell,WCell
                Dfe = DFEW(i+1,j,k,iu+2*iv+3*iw)/TCell%Vof(i,j,k)
                Dfw = DFEW(i,j,k,iu+2*iv+3*iw)/TCell%Vof(i,j,k)
                Dfn = DFNS(i,j+1,k,iu+2*iv+3*iw)/TCell%Vof(i,j,k)
                Dfs = DFNS(i,j,k,iu+2*iv+3*iw)/TCell%Vof(i,j,k)
                Dft = DFTB(i,j,k+1,iu+2*iv+3*iw)/TCell%Vof(i,j,k)
                Dfb = DFTB(i,j,k,iu+2*iv+3*iw)/TCell%Vof(i,j,k)
                
                if(dabs(MaxDiff)<dabs(DFEW(i+1,j,k,iu+2*iv+3*iw))) then
                  MaxDiff = DFEW(i+1,j,k,iu+2*iv+3*iw)
                  ii = i+1
                  jj = j
                  kk = k
                end if 
                
                Fep = EDFEW(i+1,j,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaE(i,j,k))/    &
                      TCell%Vof(i,j,k)
                Fem = EDFEW(i+1,j,k,iu+2*iv+3*iw)*TCell%EtaE(i,j,k)/	       &
                      TCell%Vof(i,j,k)		
                Fwp = EDFEW(i,j,k,iu+2*iv+3*iw)*TCell%EtaE(i-1,j,k)/	       &
                      TCell%Vof(i,j,k) 	
                Fwm = EDFEW(i,j,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaE(i-1,j,k))/    &
                      TCell%Vof(i,j,k) 	
                Fnp = EDFNS(i,j+1,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaN(i,j,k))/    &
                      TCell%Vof(i,j,k)
                Fnm = EDFNS(i,j+1,k,iu+2*iv+3*iw)*TCell%EtaN(i,j,k)/           &
                      TCell%Vof(i,j,k)
                Fsp = EDFNS(i,j,k,iu+2*iv+3*iw)*TCell%EtaN(i,j-1,k)/           &
                      TCell%Vof(i,j,k) 
                Fsm = EDFNS(i,j,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaN(i,j-1,k))/    &
                      TCell%Vof(i,j,k)
                Ftp = EDFTB(i,j,k+1,iu+2*iv+3*iw)*(1.d0-TCell%EtaT(i,j,k))/    &
                      TCell%Vof(i,j,k)
                Ftm = EDFTB(i,j,k+1,iu+2*iv+3*iw)*TCell%EtaT(i,j,k)/           &
                      TCell%Vof(i,j,k)
                Fbp = EDFTB(i,j,k,iu+2*iv+3*iw)*TCell%EtaT(i,j,k-1)/           &
                      TCell%Vof(i,j,k)
                Fbm = EDFTB(i,j,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaT(i,j,k-1))/    &
                      TCell%Vof(i,j,k)

                if(i==Imax-iu) Dfe = 0.d0
                if(j==iu.or.j==iw) Dfs = 0.d0
                if(j==Jmax+1-iu.or.j==Jmax+1-iw) Dfn = 0.d0
                if(k==iu.or.k==iv) Dfb = 0.d0
                if(k==Kmax+1-iu.or.k==Kmax+1-iv) Dft = 0.d0
                aE = Dfe;aW = Dfw;aN = Dfn;aS = Dfs;aT = Dft;aB = Dfb
		
                if(i==1)CWE(j,k,1) = aW
                if(i==Imax-iu)CWE(j,k,2) = aE
                if(j==1)CSN(i,k,1) = aS
                if(j==Jmax-iv)CSN(i,k,2) = aN
                if(k==1)CBT(i,j,1) = aB
                if(k==Kmax-iw)CBT(i,j,2) = aT
                Sp = Fric(i,j,k)
                aP = TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*         		       &
                                      TGrid%dz(i,j,k)/dt+aE+aW+aN+aS+aT+aB+Sp                  
                aP = aP+Fep-Fwp+Fnp-Fsp+Ftp-Fbp
                aE = aE-Fem
                aW = aW+Fwm
                aN = aN-Fnm
                aS = aS+Fsm
                aT = aT-Ftm
                aB = aB+Fbm
                PUVW%Dp(i,j,k)=TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k)/&
                               (aP-aE-aW-aN-aS-aT-aB)
                if(aP<0.d0) then
                  print*, 'Negative central coefficient'
                  print*, aE,aW
                  print*, aN,aS
                  print*, aT,aB
                  print*, '++++++++++++++++++++++++++++'
                  print*, Fep,Fwp
                  print*, Fnp,Fsp
                  print*, Ftp,Fbp
                end if                 
                if(isnan(PUVW%Dp(i,j,k))) then
                  print*,'Test compute the velocity velocity'
                  print*,PUVW%Dp(i,j,k)
                  print*,aE,aW
                  print*,aN,aS
                  print*,aT,aB
                  print*,'----------------------------------'
                end if
                If(dabs(aw)+dabs(ae)+dabs(as)+dabs(an)+dabs(ab)+dabs(at)>dabs(ap)) then
                  print*, 'undominant coefficient'
                  print*,i,j,k
                  print*,iu,iv,iw
                  print*, '------------------------------------'
                  print*,TCell%vof(i,j,k)
                  print*,ap,aw,ae,as,an,ab,sp
                  print*,fep,fem,TCell%EEArea(i,j,k),TCell%EtaE(i,j,k)
                  print*,fwp,fwm,TCell%EEArea(i-1,j,k),TCell%EtaE(i-1,j,k)
                  print*,fnp,fnm,TCell%NEArea(i,j,k),TCell%EtaN(i,j,k)
                  print*,fsp,fsm,TCell%NEArea(i,j-1,k),TCell%EtaN(i,j-1,k)
                  print*,ftp,ftm,TCell%TEArea(i,j,k),TCell%EtaT(i,j,k)
                  print*,fbp,fbm,TCEll%TEArea(i,j,k-1),TCell%EtaT(i,j,k-1)
                  print*,i,j,k,iu,iv,iw
                  pause
                  print*,'+++++++++++++++++++++++++++++++++++++'
                End if
                ictr=Tcell%Posnu(i,j,k)
                nnz=0
                values=0.d0
                cols=0
              ! West of current cell
                if(i>1) then
                  if(TCell%Posnu(i-1,j,k)/=-1) then
                    cols(nnz) = TCell%Posnu(i-1,j,k)
                    values(nnz) = -aW
                    nnz = nnz+1
                  end if
                end if
              ! South of current cell
                if(j>1) then
                  if(TCell%Posnu(i,j-1,k)/=-1) then
                    cols(nnz) = TCell%Posnu(i,j-1,k)
                    values(nnz) = -aS
                    nnz = nnz+1
                  end if
                end if
              ! Bottom of current
                if(k>1) then
                  if(TCell%Posnu(i,j,k-1)/=-1) then
                    cols(nnz) = TCell%Posnu(i,j,k-1)
                    values(nnz) = -aB
                    nnz = nnz+1
                  end if
                end if
              ! Set the diagonal cell
                cols(nnz) = TCell%Posnu(i,j,k)
                values(nnz) = aP
                nnz = nnz+1
              ! East of current cell
                if(i<Imax-iu) then
                  if(TCell%Posnu(i+1,j,k)/=-1) then
                    cols(nnz) = TCell%Posnu(i+1,j,k)
                    values(nnz) = -aE
                    nnz = nnz+1
                  end if
                end if
              ! North of current cell
                if(j<Jmax-iv) then
                  if(TCell%Posnu(i,j+1,k)/=-1) then
                    cols(nnz) = TCell%Posnu(i,j+1,k)
                    values(nnz) = -aN
                    nnz = nnz+1
                  end if
                end if
              ! Top of current cell
                if(k<Kmax-iw) then
                  if(TCell%Posnu(i,j,k+1)/=-1) then
                    cols(nnz) = TCell%Posnu(i,j,k+1)
                    values(nnz) = -aT
                    nnz = nnz+1
                  end if
                end if
                call HYPRE_IJMatrixSetValues(A,1,nnz,ictr,cols,values,ierr)
              end if
            end do
          end do
        end do
     !   print*,'Test maximum diffusive coefficient'
     !   print*, MaxDiff
     !   print*, iu,iv,iw
     !   print*, ii,jj,kk
        Call HYPRE_IJMatrixAssemble(A,ierr)
        Call HYPRE_IJMatrixGetObject(A,parcsr_A,ierr)
    End subroutine SetMatrix

    Subroutine SetVectors(b,x,par_b,par_x,TGrid,TCell,PUVW,CWE,CSN,CBT,IJKFlux,&
                                                            TVar_n,dt,iu,iv,iw)
        Integer*8:: b,x,par_b,par_x
        Integer,intent(in):: iu,iv,iw
        Type(Grid),intent(in):: TGrid
        Type(Cell),intent(in):: TCell
        Type(PoissonCoefficient),intent(inout):: PUVW
        Type(Variables),intent(in):: TVar_n
        Real(dp),dimension(:,:,:),allocatable,intent(in):: CWE,CSN,CBT
        Real(dp),dimension(:,:,:),intent(in):: IJKFlux
        Real(dp),intent(in):: dt
        Integer(kind=it4b):: i,j,k,ii,jj,kk,imp,jmp,kmp
        Integer:: ilower,iupper,ictr,local_size
        Integer(kind=it4b),dimension(:),allocatable:: rows
        Real(kind=dp),dimension(:),allocatable:: rhs,xval
        real(kind=dp) :: MaxVect
        ilower = 0
        iupper = TCell%ExtCell
        local_size = iupper-ilower+1 ! the number of rows
        ! In here, we apply boundary condition for deltaP with its values is 0 at
        ! all boundary. therefore, we do not need to set boundary in vector b
        Allocate(rhs(0:TCell%ExtCell))
        Allocate(xval(0:TCell%ExtCell))
        Allocate(rows(0:TCell%ExtCell))
        Call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,b,ierr)
        Call HYPRE_IJVectorSetObjectType(b,HYPRE_PARCSR,ierr)
        Call HYPRE_IJVectorInitialize(b,ierr)
        Call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,x,ierr)
        Call HYPRE_IJVectorSetObjectType(x,HYPRE_PARCSR,ierr)
        Call HYPRE_IJVectorInitialize(x,ierr)
        rhs(:) = 0.d0
        imp=1
        jmp=1
        kmp=1
        Do i = 1,Imax-iu
          Do j = 1,Jmax-iv
            Do k = 1,Kmax-iw
              If(TCell%Cell_Type(i,j,k)/=2) then
                ictr=TCell%PosNu(i,j,k)
                rhs(ictr)=(dble(iu)*TVar_n%u(i,j,k)+dble(iv)*TVar_n%v(i,j,k)+  &
                           dble(iw)*TVar_n%w(i,j,k))*TGrid%dx(i,j,k)*	       &
                           TGrid%dy(i,j,k)*TGrid%dz(i,j,k)/dt
                           
                rhs(ictr)=rhs(ictr)-IJKFlux(i,j,k)
                
                if(MaxVect<dabs(rhs(ictr))) then
                  MaxVect=rhs(ictr)
                  imp=i
                  jmp=j
                  kmp=k
                end if
                if(i==120.and.j==1.and.k==1) then
                  print*, rhs(ictr)
                  print*, (iu*TVar_n%u(i,j,k)+iv*TVar_n%v(i,j,k)+              &
                           iw*TVar_n%w(i,j,k))*TGrid%dx(i,j,k)*		       &
                           TGrid%dy(i,j,k)*TGrid%dz(i,j,k)/dt
                  print*, iu,iv,iw, dt
                  print*, TGrid%dx(i,j,k),TGrid%dy(i,j,k),TGrid%dz(i,j,k)
                  print*, TVar_n%u(i,j,k),TVar_n%v(i,j,k),TVar_n%w(i,j,k)
                end if
                
                If(i==1)rhs(ictr)=rhs(ictr)+CWE(j,k,1)*(iu*TVar_n%u(i-1,j,k)+  &
                             iv*TVar_n%v(i-1,j,k)+iw*TVar_n%w(i-1,j,k))
                If(i==Imax-iu)rhs(ictr)=rhs(ictr)+CWE(j,k,2)*                  &
                  (iu*TVar_n%u(i+1,j,k)+iv*TVar_n%v(i+1,j,k)+iw*TVar_n%w(i+1,j,k))
                If(j==1)rhs(ictr)=rhs(ictr)+CSN(i,k,1)*                        &
                  (iu*TVar_n%u(i,j-1,k)+iv*TVar_n%v(i,j-1,k)+iw*TVar_n%w(i,j-1,k))
                If(j==Jmax-iv)rhs(ictr)=rhs(ictr)+CSN(i,k,2)*                  &
                  (iu*TVar_n%u(i,j+1,k)+iv*TVar_n%v(i,j+1,k)+iw*TVar_n%w(i,j+1,k))
                If(k==1)rhs(ictr)=rhs(ictr)+CBT(i,j,1)*                        &
                  (iu*TVar_n%u(i,j,k-1)+iv*TVar_n%v(i,j,k-1)+iw*TVar_n%w(i,j,k-1))
                If(k==Kmax-iw)rhs(ictr)=rhs(ictr)+CBT(i,j,2)*                  &
                  (iu*TVar_n%u(i,j,k+1)+iv*TVar_n%v(i,j,k+1)+iw*TVar_n%w(i,j,k+1))
                xval(ictr) = 0.d0
                rows(ictr) = ilower+ictr
              End if
            End do
          End do
        End do
        print*,'Test maximum vector'
        print*, MaxVect
        print*, imp,jmp,kmp
        print*, TVar_n%u(imp,jmp,kmp),TVar_n%v(imp,jmp,kmp),TVar_n%w(imp,jmp,kmp)
        print*, IJKFlux(imp,jmp,kmp),CWE(jmp,kmp,2)
        print*, '==============================================================='
        Call HYPRE_IJVectorSetValues(b,local_size,rows,rhs,ierr)
        Call HYPRE_IJVectorSetValues(x,local_size,rows,xval,ierr)
        Call HYPRE_IJVectorAssemble(b,ierr)
        Call HYPRE_IJVectorAssemble(x,ierr)
	! get the x and b objects
        Call HYPRE_IJVectorGetObject(b,par_b,ierr)
        Call HYPRE_IJVectorGetObject(x,par_x,ierr)
        Deallocate(rhs)
        Deallocate(xval)
        Deallocate(rows)
    End subroutine SetVectors

    Subroutine DeltaGetValues(x,TCell,Var,iu,iv,iw)
        Integer*8,intent(in):: x
        Integer,intent(in):: iu,iv,iw
        Type(Cell),intent(in):: TCell
        Real(kind=dp),dimension(:,:,:),allocatable,intent(inout):: Var
        Integer(kind=it4b):: i,j,k
        Integer(kind=it4b):: ilower,iupper,local_size,ctr
        Integer(kind=it4b),dimension(:),allocatable:: rows
        Real(kind=dp),dimension(:),allocatable:: values
        ilower = 0
        iupper = TCell%ExtCell
        local_size = TCell%ExtCell+1 ! number of element
        Allocate(values(ilower:iupper))
        Allocate(rows(ilower:iupper))
        Do i = 1,Imax-iu
          Do j = 1,Jmax-iv
            Do k = 1,Kmax-iw
              If(TCell%Cell_Type(i,j,k)/=2) then
                rows(TCell%PosNu(i,j,k)) = TCell%PosNu(i,j,k)+ilower
              End if
            End do
          End do
        End do
        Call HYPRE_IJVectorGetValues(x,local_size,rows,values,ierr)
        ctr = 0
        Do i = 1,Imax-iu
          Do j = 1,Jmax-iv
            Do k = 1,Kmax-iw
              If(TCell%PosNu(i,j,k)==ctr) then
                Var(i,j,k) = values(ctr)
                ctr = ctr+1
              Else
                Var(i,j,k) = 0.d0
              End if
            End do
          End do
        End do
        Deallocate(values,rows)
    End subroutine DeltaGetValues
! Face flux with MUSCL scheme
    Subroutine ModifiedConvectiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,     &
                                      VCell,WCell,flux,iu,iv,iw)
      Implicit none
      Type(Grid),intent(in):: PGrid,UGrid,VGrid,WGrid
      Type(Cell),intent(in):: PCell,UCell,VCell,WCell
      Integer(kind=it4b),intent(in):: iu,iv,iw
      Real(kind=dp),dimension(:,:,:,:),allocatable,intent(inout):: flux
      Integer(kind=it4b):: i,j,k
      Real(kind=dp):: epsi,delh,delhec,delhec1,delhec2
      Real(kind=dp):: eta,sx,sy,sz,uw,vs,wb,uwp,uwn,vsp,vsn,wbp,wbn
      Real(kind=dp):: vw,vb,ww,ws,us,ub
      epsi = 1.d-20
      Do i = 1,Imax+iu
        Do j = 1,Jmax+iv
          Do k = 1,Kmax+iw
            If(iu==1) then
              uw = 0.5d0*(u(i,j,k)+u(i-1,j,k))
            ! For UCell both for convective velocity and scalar velocity
              If(i==1) then
                Flux(i,j,k,1) = uw**2.d0*UCell%EEArea(i,j,k)*UGrid%dy(i,j,k)*  &
                                         UGrid%dz(i,j,k)
              Elseif(i>=Imax) then
                Flux(i,j,k,1) = uw**2.d0*UCell%EEArea(i-1,j,k)*                &
                             UGrid%dy(i-1,j,k)*UGrid%dz(i-1,j,k)
              Else
                eta = UCell%EtaE(i-1,j,k)
                uw = (1.d0-eta)*u(i-1,j,k)+eta*u(i,j,k)
              ! central second order
                Flux(i,j,k,1) = (uw*UCell%AlE(i-1,j,k))**2.d0                  &
                      *UCell%EEArea(i-1,j,k)*UGrid%dy(i-1,j,k)*UGrid%dz(i-1,j,k)
              ! Central Second order
              !  Flux(i,j,k,1) = uw**2.d0*UCell%EEArea(i-1,j,k)*                &
              !                           UGrid%dy(i-1,j,k)
              End if
              ! Convective velocity: u, scalar advective : v
              uw = 0.5d0*(u(i-1,j+1,k)+u(i-1,j,k))
              If(i>Imax) then
                Flux(i,j,k,2) = uw*0.5d0*(v(i-1,j,k)+v(i,j,k))*                &
                       VCell%EEArea(i-1,j,k)*VGrid%dy(i-1,j,k)*VGrid%dz(i-1,j,k)
              Elseif(i==1) then
                Flux(i,j,k,2) = uw*0.5d0*(v(i-1,j,k)+v(i,j,k))*                &
                       VCell%EEArea(i,j,k)*VGrid%dy(i,j,k)*VGrid%dz(i,j,k)
              Else
                Flux(i,j,k,2) = 0.d0
                If(PCell%EEArea(i-1,j,k)<=1.d-5.or.PCell%EEArea(i-1,j+1,k)     &
                                                                   <=1.d-5)then
                  delhec1 = VCell%FCE(i-1,j,k,1)*VCell%nx(i-1,j,k)+            &
                    VCell%FCE(i-1,j,k,2)*VCell%ny(i-1,j,k)+VCell%FCE(i-1,j,k,3)&
                                        *VCell%nz(i-1,j,k)+VCell%phi(i-1,j,k)
                  delhec2 = (VCell%FCE(i-1,j,k,1)-UGrid%dx(i-1,j,k))*          &
                    VCell%nx(i,j,k)+VCell%FCE(i-1,j,k,2)*VCell%ny(i,j,k)+      &
                      VCell%FCE(i-1,j,k,3)*VCell%nz(i,j,k)+VCell%phi(i,j,k)
                  delhec = 0.5d0*(delhec1+delhec2)
                  delh = delhec
                  If(UCell%MoExCell(i-1,j+1,k)/=1.and.UCell%Cell_Type(i-1,j+1,k)&
                          /=2.and.UCell%vof(i-1,j+1,k)>=UCell%vof(i-1,j,k))then
                    delh = dabs(UCell%Cell_Cent(i-1,j+1,k,1)*UCell%nx(i-1,j+1,k)&
                        +UCell%Cell_Cent(i-1,j+1,k,2)*UCell%ny(i-1,j+1,k)+     &
                         UCell%Cell_Cent(i-1,j+1,k,3)*UCell%nz(i-1,j+1,k)+     &
                         UCell%phi(i-1,j+1,k))
                    uw = u(i-1,j+1,k)
                  End if
                  If(UCell%MoExCell(i-1,j,k)/=1.and.UCell%Cell_Type(i-1,j,k)   &
                           /=2.and.UCell%vof(i-1,j,k)>UCell%vof(i-1,j+1,k))then
                    delh = dabs(UCell%Cell_Cent(i-1,j,k,1)*UCell%nx(i-1,j,k)+  &
                        UCell%Cell_Cent(i-1,j,k,2)*UCell%ny(i-1,j,k)+          &
                        UCell%Cell_Cent(i-1,j,k,3)*UCell%nz(i-1,j,k)+          &
                        UCell%phi(i-1,j,k))
                    uw = u(i-1,j,k)
                  End if
                  uw = uw*delhec/(delh+epsi)
                Else
                  Sy = UCell%SyN(i-1,j,k)
                  eta = dabs(VCell%FCE(i-1,j,k,2)+VGrid%dy(i-1,j,k)/2.d0-      &
                             UCell%Cell_Cent(i-1,j,k,2))/Sy
                  uw = (1.d0-eta)*u(i-1,j,k)+eta*u(i-1,j+1,k)
                End if
                eta = VCell%EtaE(i-1,j,k)
                vw = (1.d0-eta)*v(i-1,j,k)+eta*v(i,j,k)
                Flux(i,j,k,2) = uw*vw*VCell%AlE(i-1,j,k)*VCell%EEArea(i-1,j,k)*&
                                      VGrid%dy(i-1,j,k)*VGrid%dz(i-1,j,k)
              End if
            ! Convective velocity: u, scalar advective: w
              uw = 0.5d0*(u(i-1,j,k+1)+u(i-1,j,k))
              If(i>Imax) then
                Flux(i,j,k,3) = uw*0.5d0*(w(i-1,j,k)+w(i,j,k))*                &
                       WCell%EEArea(i-1,j,k)*WGrid%dy(i-1,j,k)*WGrid%dz(i-1,j,k)
              Elseif(i==1) then
                Flux(i,j,k,3) = uw*0.5d0*(w(i-1,j,k)+w(i,j,k))*                &
                       WCell%EEArea(i,j,k)*WGrid%dy(i,j,k)*WGrid%dz(i,j,k)
              Else
                Flux(i,j,k,3) = 0.d0
                If(PCell%EEArea(i-1,j,k)<=1.d-5.or.PCell%EEArea(i-1,j,k+1)     &
                                                                   <=1.d-5)then
                  delhec1 = WCell%FCE(i-1,j,k,1)*WCell%nx(i-1,j,k)+            &
                    WCell%FCE(i-1,j,k,2)*WCell%ny(i-1,j,k)+WCell%FCE(i-1,j,k,3)&
                                        *WCell%nz(i-1,j,k)+WCell%phi(i-1,j,k)
                  delhec2 = (WCell%FCE(i-1,j,k,1)-UGrid%dx(i-1,j,k))*          &
                    WCell%nx(i,j,k)+WCell%FCE(i-1,j,k,2)*WCell%ny(i,j,k)+      &
                      WCell%FCE(i-1,j,k,3)*WCell%nz(i,j,k)+WCell%phi(i,j,k)
                  delhec = 0.5d0*(delhec1+delhec2)
                  If(UCell%MoExCell(i-1,j,k+1)/=1.and.UCell%Cell_Type(i-1,j,k+1)&
                          /=2.and.UCell%vof(i-1,j,k+1)>=UCell%vof(i-1,j,k))then
                    delh = dabs(UCell%Cell_Cent(i-1,j,k+1,1)*UCell%nx(i-1,j,k+1)&
                        +UCell%Cell_Cent(i-1,j,k+1,2)*UCell%ny(i-1,j,k+1)+     &
                         UCell%Cell_Cent(i-1,j,k+1,3)*UCell%nz(i-1,j,k+1)+     &
                         UCell%phi(i-1,j,k+1))
                    uw = u(i-1,j,k+1)
                  End if
                  If(UCell%MoExCell(i-1,j,k)/=1.and.UCell%Cell_Type(i-1,j,k)   &
                           /=2.and.UCell%vof(i-1,j,k)>UCell%vof(i-1,j,k+1))then
                    delh = dabs(UCell%Cell_Cent(i-1,j,k,1)*UCell%nx(i-1,j,k)+  &
                        UCell%Cell_Cent(i-1,j,k,2)*UCell%ny(i-1,j,k)+          &
                        UCell%Cell_Cent(i-1,j,k,3)*UCell%nz(i-1,j,k)+          &
                        UCell%phi(i-1,j,k))
                    uw = u(i-1,j,k)
                  End if
                  uw = uw*delhec/(delh+epsi)
                Else
                  Sz = UCell%SzT(i-1,j,k)
                  eta = dabs(WCell%FCE(i-1,j,k,3)+WGrid%dz(i-1,j,k)/2.d0-      &
                             UCell%Cell_Cent(i-1,j,k,3))/Sz
                  uw = (1.d0-eta)*u(i-1,j,k)+eta*u(i-1,j,k+1)
                End if
                eta = WCell%EtaE(i-1,j,k)
                ww = (1.d0-eta)*w(i-1,j,k)+eta*w(i,j,k)
                Flux(i,j,k,3) = uw*ww*WCell%AlE(i-1,j,k)*WCell%EEArea(i-1,j,k)*&
                                      WGrid%dy(i-1,j,k)*WGrid%dz(i-1,j,k)
              End if
            End if
            If(iv==1) then ! Jflux
          ! Convective velocity: v, scalar advective: u
              vs = 0.5d0*(v(i,j-1,k)+v(i+1,j-1,k))
              If(j>Jmax) then
                Flux(i,j,k,1) = vs*0.5d0*(u(i,j,k)+u(i,j-1,k))*                &
                       UCell%NEArea(i,j-1,k)*UGrid%dx(i,j-1,k)*UGrid%dz(i,j-1,k)
              Elseif(j==1) then
                Flux(i,j,k,1) = vs*0.5d0*(u(i,j,k)+u(i,j-1,k))*                &
                       UCell%NEArea(i,j,k)*UGrid%dx(i,j,k)*UGrid%dz(i,j,k)
              Else
                Flux(i,j,k,1) = 0.d0
                If(PCell%NEArea(i,j-1,k)<=1.d-5.or.PCell%NEArea(i+1,j-1,k)    &
                                                                   <=1.d-5)then
                  delhec1 = UCell%FCN(i,j-1,k,1)*UCell%nx(i,j-1,k)+            &
                    UCell%FCN(i,j-1,k,2)*UCell%ny(i,j-1,k)+UCell%FCN(i,j-1,k,3)&
                                        *UCell%nz(i,j-1,k)+UCell%phi(i,j-1,k)
                  delhec2 = UCell%FCN(i,j-1,k,1)*UCell%nx(i,j,k)+              &
                   (UCell%FCN(i,j-1,k,2)-VGrid%dy(i,j-1,k))*UCell%ny(i,j,k)+   &
                      UCell%FCN(i,j-1,k,3)*UCell%nz(i,j,k)+UCell%phi(i,j,k)
                  delhec = 0.5d0*(delhec1+delhec2)
                  If(VCell%MoExCell(i,j-1,k)/=1.and.VCell%Cell_Type(i,j-1,k)   &
                         /=2.and.VCell%vof(i,j-1,k)>=VCell%vof(i+1,j-1,k)) then
                    delh = dabs(VCell%Cell_Cent(i,j-1,k,1)*VCell%nx(i,j-1,k)+  &
                      VCell%Cell_Cent(i,j-1,k,2)*VCell%ny(i,j-1,k)+            &
                      VCell%Cell_Cent(i,j-1,k,3)*VCell%nz(i,j-1,k)+            &
                      VCell%phi(i,j-1,k))
                    vs = v(i,j-1,k)
                  End if
                  If(VCell%MoExCell(i+1,j-1,k)/=1.and.VCell%Cell_Type(i+1,j-1,k)&
                          /=2.and.VCell%vof(i+1,j-1,k)>VCell%vof(i,j-1,k)) then
                    delh = dabs(VCell%Cell_Cent(i+1,j-1,k,1)*VCell%nx(i+1,j-1,k)&
                          +VCell%Cell_Cent(i+1,j-1,k,2)*VCell%ny(i+1,j-1,k)+   &
                           VCell%Cell_Cent(i+1,j-1,k,3)*VCell%nz(i+1,j-1,k)+   &
                           VCell%phi(i+1,j-1,k))
                    vs = v(i+1,j-1,k)
                  End if
                  vs = vs*delhec/(delh+epsi)
                Else
                  Sx = VCell%SxE(i,j-1,k)
                  eta = dabs(UCell%FCN(i,j-1,k,1)+0.5d0*UGrid%dx(i,j-1,k)-     &
                                            VCell%Cell_Cent(i,j-1,k,1))/Sx
                  vs = (1.d0-eta)*v(i,j-1,k)+eta*v(i+1,j-1,k)
                End if
                eta = UCell%EtaN(i,j-1,k)
                us = (1.d0-eta)*u(i,j-1,k)+eta*u(i,j,k)
                Flux(i,j,k,1) = us*vs*UCell%AlN(i,j-1,k)*UCell%NEArea(i,j-1,k)*&
                                      UGrid%dx(i,j-1,k)*UGrid%dz(i,j-1,k)
              End if
            ! Convective velocity: v, scalar convective: v
              vs = 0.5d0*(v(i,j,k)+v(i,j-1,k))
              If(j>=Jmax) then
                Flux(i,j,k,2) = vs**2.d0*VCell%NEArea(i,j-1,k)*                &
                                             VGrid%dx(i,j-1,k)*VGrid%dz(i,j-1,k)
              Elseif(j==1) then
                Flux(i,j,k,2) = vs**2.d0*VCell%NEArea(i,j,k)*VGrid%dx(i,j,k)*  &
                                                             VGrid%dz(i,j,k)
              Else
                eta = VCell%EtaN(i,j-1,k)
                vs = (1.d0-eta)*v(i,j-1,k)+eta*v(i,j,k)
                Flux(i,j,k,2) = (vs*VCell%AlN(i,j-1,k))**2.d0*                 &
                       VCell%NEArea(i,j-1,k)*VGrid%dx(i,j-1,k)*VGrid%dz(i,j-1,k)
              End if
            ! Convective velocity: v, scalar advective: w
              vs = 0.5d0*(v(i,j-1,k)+v(i,j-1,k+1))
              If(j>Jmax) then
                Flux(i,j,k,3) = vs*0.5d0*(w(i,j,k)+w(i,j-1,k))*                &
                       WCell%NEArea(i,j-1,k)*WGrid%dx(i,j-1,k)*WGrid%dz(i,j-1,k)
              Elseif(j==1) then
                Flux(i,j,k,3) = vs*0.5d0*(w(i,j,k)+w(i,j-1,k))*                &
                       WCell%NEArea(i,j,k)*WGrid%dx(i,j,k)*wGrid%dz(i,j,k)
              Else
                Flux(i,j,k,3) = 0.d0
                If(PCell%NEArea(i,j-1,k)<=1.d-5.or.PCell%NEArea(i,j-1,k+1)    &
                                                                   <=1.d-5)then
                  delhec1 = WCell%FCN(i,j-1,k,1)*WCell%nx(i,j-1,k)+            &
                    WCell%FCN(i,j-1,k,2)*WCell%ny(i,j-1,k)+WCell%FCN(i,j-1,k,3)&
                                        *WCell%nz(i,j-1,k)+WCell%phi(i,j-1,k)
                  delhec2 = WCell%FCN(i,j-1,k,1)*WCell%nx(i,j,k)+              &
                   (WCell%FCN(i,j-1,k,2)-WGrid%dy(i,j-1,k))*WCell%ny(i,j,k)+   &
                      WCell%FCN(i,j-1,k,3)*WCell%nz(i,j,k)+WCell%phi(i,j,k)
                  delhec = 0.5d0*(delhec1+delhec2)
                  If(VCell%MoExCell(i,j-1,k)/=1.and.VCell%Cell_Type(i,j-1,k)   &
                         /=2.and.VCell%vof(i,j-1,k)>=VCell%vof(i,j-1,k+1)) then
                    delh = dabs(VCell%Cell_Cent(i,j-1,k,1)*VCell%nx(i,j-1,k)+  &
                      VCell%Cell_Cent(i,j-1,k,2)*VCell%ny(i,j-1,k)+            &
                      VCell%Cell_Cent(i,j-1,k,3)*VCell%nz(i,j-1,k)+            &
                      VCell%phi(i,j-1,k))
                    vs = v(i,j-1,k)
                  End if
                  If(VCell%MoExCell(i,j-1,k+1)/=1.and.VCell%Cell_Type(i,j-1,k+1)&
                          /=2.and.VCell%vof(i,j-1,k+1)>VCell%vof(i,j-1,k)) then
                    delh = dabs(VCell%Cell_Cent(i,j-1,k+1,1)*VCell%nx(i,j-1,k+1)&
                          +VCell%Cell_Cent(i,j-1,k+1,2)*VCell%ny(i,j-1,k+1)+   &
                           VCell%Cell_Cent(i,j-1,k+1,3)*VCell%nz(i,j-1,k+1)+   &
                           VCell%phi(i,j-1,k+1))
                    vs = v(i,j-1,k+1)
                  End if
                  vs = vs*delhec/(delh+epsi)
                Else
                  Sz = VCell%SzT(i,j-1,k)
                  eta = dabs(WCell%FCN(i,j-1,k,3)+0.5d0*WGrid%dz(i,j-1,k)-     &
                                            WCell%Cell_Cent(i,j-1,k,3))/Sz
                  vs = (1.d0-eta)*v(i,j-1,k)+eta*v(i,j-1,k+1)
                End if
                eta = WCell%EtaN(i,j-1,k)
                ws = (1.d0-eta)*w(i,j-1,k)+eta*w(i,j,k)
                Flux(i,j,k,3) = vs*ws*WCell%AlN(i,j-1,k)*WCell%NEArea(i,j-1,k)*&
                                      WGrid%dx(i,j-1,k)*WGrid%dz(i,j-1,k)
              End if
            End if
            If(iw==1) then
          ! Convective velocity: w, scalar advective: u
              wb = 0.5d0*(w(i,j,k-1)+v(i+1,j,k-1))
              If(k>Kmax) then
                Flux(i,j,k,1) = wb*0.5d0*(u(i,j,k)+u(i,j,k-1))*                &
                       UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k-1)*UGrid%dy(i,j,k-1)
              Elseif(k==1) then
                Flux(i,j,k,1) = wb*0.5d0*(u(i,j,k)+u(i,j,k-1))*                &
                       UCell%TEArea(i,j,k)*UGrid%dx(i,j,k)*UGrid%dy(i,j,k)
              Else
                Flux(i,j,k,1) = 0.d0
                If(PCell%TEArea(i,j,k-1)<=1.d-5.or.PCell%TEArea(i+1,j,k-1)    &
                                                                   <=1.d-5)then
                  delhec1 = UCell%FCT(i,j,k-1,1)*UCell%nx(i,j,k-1)+            &
                    UCell%FCT(i,j,k-1,2)*UCell%ny(i,j,k-1)+UCell%FCT(i,j,k-1,3)&
                                        *UCell%nz(i,j,k-1)+UCell%phi(i,j,k-1)
                  delhec2 = UCell%FCT(i,j,k-1,1)*UCell%nx(i,j,k)+              &
                   (UCell%FCT(i,j,k-1,2)-WGrid%dz(i,j,k-1))*UCell%ny(i,j,k)+   &
                      UCell%FCT(i,j,k-1,3)*UCell%nz(i,j,k)+UCell%phi(i,j,k)
                  delhec = 0.5d0*(delhec1+delhec2)
                  If(WCell%MoExCell(i,j,k-1)/=1.and.WCell%Cell_Type(i,j,k-1)   &
                         /=2.and.WCell%vof(i,j,k-1)>=WCell%vof(i+1,j,k-1)) then
                    delh = dabs(WCell%Cell_Cent(i,j,k-1,1)*WCell%nx(i,j,k-1)+  &
                                WCell%Cell_Cent(i,j,k-1,2)*WCell%ny(i,j,k-1)+  &
                                WCell%Cell_Cent(i,j,k-1,3)*WCell%nz(i,j,k-1)+  &
                                WCell%phi(i,j,k-1))
                    wb = w(i,j,k-1)
                  End if
                  If(WCell%MoExCell(i+1,j,k-1)/=1.and.WCell%Cell_Type(i+1,j,k-1)&
                          /=2.and.WCell%vof(i+1,j,k-1)>WCell%vof(i,j,k-1)) then
                    delh = dabs(WCell%Cell_Cent(i+1,j,k-1,1)*WCell%nx(i+1,j,k-1)&
                          +WCell%Cell_Cent(i+1,j,k-1,2)*WCell%ny(i+1,j,k-1)+   &
                           WCell%Cell_Cent(i+1,j,k-1,3)*WCell%nz(i+1,j,k-1)+   &
                           WCell%phi(i+1,j,k-1))
                    wb = w(i+1,j,k-1)
                  End if
                  wb = wb*delhec/(delh+epsi)
                Else
                  Sx = WCell%SxE(i,j,k-1)
                  eta = dabs(UCell%FCT(i,j,k-1,1)+0.5d0*UGrid%dx(i,j,k-1)-     &
                                            WCell%Cell_Cent(i,j,k-1,1))/Sx
                  wb = (1.d0-eta)*w(i,j,k-1)+eta*w(i+1,j,k-1)
                End if
                eta = UCell%EtaT(i,j,k-1)
                ub = (1.d0-eta)*u(i,j,k-1)+eta*u(i,j,k)
                Flux(i,j,k,1) = wb*ub*UCell%AlT(i,j,k-1)*UCell%TEArea(i,j,k-1)*&
                                      UGrid%dx(i,j,k-1)*UGrid%dy(i,j,k-1)
                If(isnan(flux(i,j,k,1))) then
                  print*, UCell%FCT(i,j,k-1,1),UCell%FCT(i,j,k-1,2),UCell%FCT(i,j,k-1,3)
                  print*, PCell%TEArea(i,j,k-1),PCell%TEArea(i+1,j,k-1)
                  print*,delhec1,delhec2,delh,wb,delhec/delh,i,j,k
                  print*,wbp,wbn,UCell%TEArea(i,j,k-1),WCell%SxE(i,j,k-1)
                  print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
                End if
             ! Convective velocity: w, scalar advective: v
              End if
              wb = 0.5d0*(w(i,j+1,k-1)+w(i,j,k-1))
              If(k>Kmax) then
                Flux(i,j,k,2) = wb*0.5d0*(v(i,j,k-1)+v(i,j,k))*                &
                       VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k-1)*VGrid%dy(i,j,k-1)
              Elseif(k==1) then
                Flux(i,j,k,2) = wb*0.5d0*(v(i,j,k-1)+v(i,j,k))*                &
                       VCell%TEArea(i,j,k)*VGrid%dx(i,j,k)*VGrid%dy(i,j,k)
              Else
                Flux(i,j,k,2) = 0.d0
                If(PCell%TEArea(i,j,k-1)<=1.d-5.or.PCell%TEArea(i,j+1,k-1)     &
                                                                   <=1.d-5)then
                  delhec1 = VCell%FCT(i,j,k-1,1)*VCell%nx(i,j,k-1)+            &
                    VCell%FCT(i,j,k-1,2)*VCell%ny(i,j,k-1)+VCell%FCT(i,j,k-1,3)&
                                        *VCell%nz(i,j,k-1)+VCell%phi(i,j,k-1)
                  delhec2 = (VCell%FCT(i,j,k-1,1)-WGrid%dz(i,j,k-1))*          &
                      VCell%nx(i,j,k)+VCell%FCT(i,j,k-1,2)*VCell%ny(i,j,k)+    &
                      VCell%FCT(i,j,k-1,3)*VCell%nz(i,j,k)+VCell%phi(i,j,k)
                  delhec = 0.5d0*(delhec1+delhec2)
                  delh = delhec
                  If(WCell%MoExCell(i,j+1,k-1)/=1.and.WCell%Cell_Type(i,j+1,k-1)&
                          /=2.and.WCell%vof(i,j+1,k-1)>=WCell%vof(i,j,k-1))then
                    delh = dabs(WCell%Cell_Cent(i,j+1,k-1,1)*WCell%nx(i,j+1,k-1)&
                        +WCell%Cell_Cent(i,j+1,k-1,2)*WCell%ny(i,j+1,k-1)+     &
                         WCell%Cell_Cent(i,j+1,k-1,3)*WCell%nz(i,j+1,k-1)+     &
                         WCell%phi(i,j+1,k-1))
                    wb = w(i,j+1,k-1)
                  End if
                  If(WCell%MoExCell(i,j,k-1)/=1.and.WCell%Cell_Type(i,j,k-1)   &
                           /=2.and.WCell%vof(i,j,k-1)>WCell%vof(i,j+1,k-1))then
                    delh = dabs(WCell%Cell_Cent(i,j,k-1,1)*WCell%nx(i,j,k-1)+  &
                        WCell%Cell_Cent(i,j,k-1,2)*WCell%ny(i,j,k-1)+          &
                        WCell%Cell_Cent(i,j,k-1,3)*WCell%nz(i,j,k-1)+          &
                        WCell%phi(i,j,k-1))
                    wb = w(i,j,k-1)
                  End if
                  wb = wb*delhec/(delh+epsi)
                Else
                  Sy = WCell%SyN(i,j,k-1)
                  eta = dabs(VCell%FCT(i,j,k-1,2)+0.5d0*VGrid%dy(i,j,k-1)-     &
                                                WCell%Cell_Cent(i,j,k-1,2))/Sy
                  wb = (1.d0-eta)*w(i,j,k-1)+eta*w(i,j+1,k-1)
                End if
                eta = VCell%EtaT(i,j,k-1)
                vb = (1.d0-eta)*v(i,j,k-1)+eta*v(i,j,k)
                Flux(i,j,k,2) = wb*vb*VCell%AlT(i,j,k-1)*VCell%TEArea(i,j,k-1)*&
                                      VGrid%dx(i,j,k-1)*VGrid%dy(i,j,k-1)
             ! Convective velocity: w, scalar advective: w
              End if
              wb = 0.5d0*(w(i,j,k)+w(i,j,k-1))
              If(k>=Kmax) then
                Flux(i,j,k,3) = wb**2.d0*WCell%TEArea(i,j,k-1)*                &
                                         WGrid%dx(i,j,k-1)*WGrid%dy(i,j,k-1)
              Elseif(k==1) then
                Flux(i,j,k,3) = wb**2.d0*WCell%TEArea(i,j,k)*WGrid%dx(i,j,k)*  &
                                                             WGrid%dy(i,j,k)
              Else
                eta = WCell%EtaT(i,j,k-1)
                wb = (1.d0-eta)*w(i,j,k-1)+eta*w(i,j,k)
                Flux(i,j,k,3) = (wb*WCell%AlT(i,j,k-1))**2.d0*                 &
                      WCell%TEArea(i,j,k-1)*WGrid%dx(i,j,k-1)*WGrid%dy(i,j,k-1)
              End if
            End if
          End do
        End do
      End do
    End subroutine ModifiedConvectiveFlux

    Subroutine DiffusiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,  &
                                                           flux,EFlux,iu,iv,iw)
      Implicit none
      Integer(kind=it4b),intent(in):: iu,iv,iw
      Type(Cell),intent(in):: PCell,UCell,VCell,WCell
      Type(Grid),intent(in):: PGrid,UGrid,VGrid,WGrid
      Integer(kind=it4b):: i,j,k
      Real(kind=dp),dimension(:,:,:,:),allocatable,intent(inout):: flux,EFlux
      Real(kind=dp):: Sx,Sy,Sz,tol
      tol = 1.d-14
      Do i = 1,Imax+iu
        Do j = 1,Jmax+iv
          Do k = 1,Kmax+iw
            If(iu==1) then
              If(i==Imax+iu) then
                flux(i,j,k,1) = UCell%EEArea(i-1,j,k)*UGrid%dy(i-1,j,k)*       &
                                UGrid%dz(i-1,j,k)/PGrid%dx(i-1,j,k)/Rey
                flux(i,j,k,2) = VCell%EEArea(i-1,j,k)*VGrid%dy(i-1,j,k)*       &
                                     VGrid%dz(i-1,j,k)/UGrid%dx(i-1,j,k)/Rey
                flux(i,j,k,3) = WCell%EEArea(i-1,j,k)*WGrid%dy(i-1,j,k)*       &
                                     WGrid%dz(i-1,j,k)/UGrid%dx(i-1,j,k)/Rey
              Elseif(i==1) then
                flux(i,j,k,1) = UCell%EEArea(i,j,k)*UGrid%dy(i,j,k)*           &
                              UGrid%dz(i,j,k)/PGrid%dx(i,j,k)/Rey
                flux(i,j,k,2) = VCell%EEArea(i,j,k)*VGrid%dy(i,j,k)*           &
                                   VGrid%dz(i,j,k)/UGrid%dx(i,j,k)/Rey
                flux(i,j,k,3) = WCell%EEArea(i,j,k)*WGrid%dy(i,j,k)*           &
                                   WGrid%dz(i,j,k)/UGrid%dx(i,j,k)/Rey
              Else
                Sx = UCell%SxE(i-1,j,k)
                Sy = UCell%Cell_Cent(i,j,k,2)-UCell%Cell_Cent(i-1,j,k,2)
                Sz = UCell%Cell_Cent(i,j,k,3)-UCell%Cell_Cent(i-1,j,k,3)
                flux(i,j,k,1) = UCell%EEArea(i-1,j,k)*UGrid%dy(i,j,k)*         &
                                                      UGrid%dz(i,j,k)/Sx/Rey
                If(dabs(Sy)>1.d-4*UGrid%dy(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         UGrid%dz(i,j,k)) then
                  EFlux(i,j,k,1) = UCell%DAlE(i-1,j,k)*UCell%EEArea(i-1,j,k)*  &
                                        UGrid%dy(i-1,j,k)*UGrid%dz(i-1,j,k)/Rey
                End if
                Sx = VCell%SxE(i-1,j,k)
                Sy = VCell%Cell_Cent(i,j,k,2)-VCell%Cell_Cent(i-1,j,k,2)
                Sz = VCell%Cell_Cent(i,j,k,3)-VCell%Cell_Cent(i-1,j,k,3)
                flux(i,j,k,2) = VCell%EEArea(i-1,j,k)*VGrid%dy(i,j,k)*         &
                                                      VGrid%dz(i,j,k)/Sx/Rey
                If(dabs(Sy)>1.d-4*VGrid%dy(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         VGrid%dz(i,j,k)) then
                  EFlux(i,j,k,2) = VCell%DAlE(i-1,j,k)*VCell%EEArea(i-1,j,k)*  &
                                        VGrid%dy(i-1,j,k)*VGrid%dz(i-1,j,k)/Rey
                End if
                Sx = WCell%SxE(i-1,j,k)
                Sy = WCell%Cell_Cent(i,j,k,2)-WCell%Cell_Cent(i-1,j,k,2)
                Sz = WCell%Cell_Cent(i,j,k,3)-WCell%Cell_Cent(i-1,j,k,3)
                flux(i,j,k,3) = WCell%EEArea(i-1,j,k)*WGrid%dy(i,j,k)*         &
                                                      WGrid%dz(i,j,k)/Sx/Rey
                If(dabs(Sy)>1.d-4*WGrid%dy(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         WGrid%dz(i,j,k)) then
                  EFlux(i,j,k,3) = WCell%DAlE(i-1,j,k)*WCell%EEArea(i-1,j,k)*  &
                                        WGrid%dy(i-1,j,k)*WGrid%dz(i-1,j,k)/Rey
                End if
              End if
            End if
            If(iv==1) then
              If(j==Jmax+iv) then
                flux(i,j,k,1) = UCell%NEArea(i,j-1,k)*UGrid%dx(i,j-1,k)*       &
                                     UGrid%dz(i,j-1,k)/VGrid%dy(i,j-1,k)/Rey
                flux(i,j,k,2) = VCell%NEArea(i,j-1,k)*VGrid%dx(i,j-1,k)*       &
                                     VGrid%dz(i,j-1,k)/PGrid%dy(i,j-1,k)/Rey
                flux(i,j,k,3) = WCell%NEArea(i,j-1,k)*WGrid%dx(i,j-1,k)*       &
                                     WGrid%dz(i,j-1,k)/VGrid%dy(i,j-1,k)/Rey
              Elseif(j==1) then
                flux(i,j,k,1) = UCell%NEArea(i,j,k)*UGrid%dx(i,j,k)*           &
                                     UGrid%dz(i,j,k)/VGrid%dy(i,j,k)/Rey
                flux(i,j,k,2) = VCell%NEArea(i,j,k)*VGrid%dx(i,j,k)*           &
                                     VGrid%dz(i,j,k)/PGrid%dy(i,j,k)/Rey
                flux(i,j,k,3) = WCell%NEArea(i,j,k)*WGrid%dx(i,j,k)*           &
                                     WGrid%dz(i,j,k)/VGrid%dy(i,j,k)/Rey
              Else
                Sx = UCell%Cell_Cent(i,j,k,1)-UCell%Cell_Cent(i,j-1,k,1)
                Sy = UCell%SyN(i,j-1,k)
                Sz = UCell%Cell_Cent(i,j,k,3)-UCell%Cell_Cent(i,j-1,k,3)
                flux(i,j,k,1) = UCell%NEArea(i,j-1,k)*UGrid%dx(i,j,k)*         &
                                                      UGrid%dz(i,j,k)/Sy/Rey
                If(dabs(Sx)>1.d-4*UGrid%dx(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         UGrid%dz(i,j,k)) then
                  EFlux(i,j,k,1) = UCell%DAlN(i,j-1,k)*UCell%NEArea(i,j-1,k)*  &
                                   UGrid%dx(i,j-1,k)*UGrid%dz(i,j-1,k)/Rey
                End if
                Sx = VCell%Cell_Cent(i,j,k,1)-VCell%Cell_Cent(i,j-1,k,1)
                Sy = VCell%SyN(i,j-1,k)
                Sz = VCell%Cell_Cent(i,j,k,3)-VCell%Cell_Cent(i,j-1,k,3)
                flux(i,j,k,2) = VCell%NEArea(i,j-1,k)*VGrid%dx(i,j,k)*         &
                                                      VGrid%dz(i,j,k)/Sy/Rey
                If(dabs(Sx)>1.d-4*VGrid%dx(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         VGrid%dz(i,j,k)) then
                  EFlux(i,j,k,2) = VCell%DAlN(i,j-1,k)*VCell%NEArea(i,j-1,k)*  &
                                   VGrid%dx(i,j-1,k)*VGrid%dz(i,j-1,k)/Rey
                End if
                Sx = WCell%Cell_Cent(i,j,k,1)-WCell%Cell_Cent(i,j-1,k,1)
                Sy = WCell%SyN(i,j-1,k)
                Sz = WCell%Cell_Cent(i,j,k,3)-WCell%Cell_Cent(i,j-1,k,3)
                flux(i,j,k,3) = WCell%NEArea(i,j-1,k)*WGrid%dx(i,j,k)*         &
                                                      WGrid%dz(i,j,k)/Sy/Rey
                If(dabs(Sx)>1.d-4*WGrid%dx(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         WGrid%dz(i,j,k)) then
                  EFlux(i,j,k,3) = WCell%DAlN(i,j-1,k)*WCell%NEArea(i,j-1,k)*  &
                                   WGrid%dx(i,j-1,k)*WGrid%dz(i,j-1,k)/Rey
                End if
              End if
            End if
            If(iw==1) then
              If(k==Kmax+iw) then
                flux(i,j,k,1) = UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k-1)*       &
                                     UGrid%dy(i,j,k-1)/WGrid%dz(i,j,k-1)/Rey
                flux(i,j,k,2) = VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k-1)*       &
                                     VGrid%dy(i,j,k-1)/WGrid%dz(i,j,k-1)/Rey
                flux(i,j,k,3) = WCell%TEArea(i,j,k-1)*WGrid%dx(i,j,k-1)*       &
                                     WGrid%dy(i,j,k-1)/PGrid%dz(i,j,k-1)/Rey
              Elseif(k==1) then
                flux(i,j,k,1) = UCell%TEArea(i,j,k)*UGrid%dx(i,j,k)*           &
                                     UGrid%dy(i,j,k)/WGrid%dz(i,j,k)/Rey
                flux(i,j,k,2) = VCell%TEArea(i,j,k)*VGrid%dx(i,j,k)*           &
                                     VGrid%dy(i,j,k)/WGrid%dz(i,j,k)/Rey
                flux(i,j,k,3) = WCell%TEArea(i,j,k)*WGrid%dx(i,j,k)*           &
                                     WGrid%dy(i,j,k)/PGrid%dz(i,j,k)/Rey
              Else
                Sx = UCell%Cell_Cent(i,j,k,1)-UCell%Cell_Cent(i,j,k-1,1)
                Sy = UCell%Cell_Cent(i,j,k,2)-UCell%Cell_Cent(i,j,k-1,2)
                Sz = UCell%SzT(i,j,k-1)
                flux(i,j,k,1) = UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k)*         &
                                                      UGrid%dy(i,j,k)/Sz/Rey
                If(dabs(Sx)>1.d-4*UGrid%dx(i,j,k).or.dabs(Sy)>1.d-4*           &
                                                         UGrid%dy(i,j,k)) then
                  EFlux(i,j,k,1) = UCell%DAlT(i,j,k-1)*UCell%TEArea(i,j,k-1)*  &
                                   UGrid%dx(i,j,k-1)*UGrid%dy(i,j,k-1)/Rey
                End if

                Sx = VCell%Cell_Cent(i,j,k,1)-VCell%Cell_Cent(i,j,k-1,1)
                Sy = VCell%Cell_Cent(i,j,k,2)-VCell%Cell_Cent(i,j,k-1,2)
                Sz = VCell%SzT(i,j,k-1)
                flux(i,j,k,2) = VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k)*         &
                                                      VGrid%dy(i,j,k)/Sz/Rey
                If(dabs(Sx)>1.d-4*VGrid%dx(i,j,k).or.dabs(Sy)>1.d-4*           &
                                                         VGrid%dy(i,j,k)) then
                  EFlux(i,j,k,2) = VCell%DAlT(i,j,k-1)*VCell%TEArea(i,j,k-1)*  &
                                   VGrid%dx(i,j,k-1)*VGrid%dy(i,j,k-1)/Rey
                End if
                Sx = WCell%Cell_Cent(i,j,k,1)-WCell%Cell_Cent(i,j,k-1,1)
                Sy = WCell%Cell_Cent(i,j,k,2)-WCell%Cell_Cent(i,j,k-1,2)
                Sz = WCell%SzT(i,j,k-1)
                flux(i,j,k,3) = WCell%TEArea(i,j,k-1)*WGrid%dx(i,j,k)*         &
                                                      WGrid%dy(i,j,k)/Sz/Rey
                If(dabs(Sx)>1.d-4*WGrid%dx(i,j,k).or.dabs(Sy)>1.d-4*           &
                                                         WGrid%dy(i,j,k)) then
                  EFlux(i,j,k,3) = WCell%DAlT(i,j,k-1)*WCell%TEArea(i,j,k-1)*  &
                                   WGrid%dx(i,j,k-1)*WGrid%dy(i,j,k-1)/Rey
                End if
              End if
            End if
          End do
        End do
      End do
    End subroutine DiffusiveFlux
    Subroutine PredictorVelocityBoundaryCondition(Pred,TVar)
      Implicit none
      Type(Predictor),intent(inout):: Pred
      Type(Variables),intent(in):: TVar
      Integer(kind=it4b):: i,j,k
      Do j = 1,Jmax
        Do k = 1,Kmax
          Pred%u(1-ight,j,k) = TVar%Uint/TVar%Uref ! 0.d0
          Pred%v(1-ight,j,k) = 0.d0-Pred%v(1,j,k)
          Pred%w(1-ight,j,k) = 0.d0-Pred%w(1,j,k)
          Pred%u(Imax,j,k) = Pred%u(Imax-1,j,k)
          Pred%u(Imax+ight,j,k) = Pred%u(Imax,j,k)
          Pred%v(Imax+ight,j,k) = Pred%v(Imax,j,k)
          Pred%w(Imax+ight,j,k) = Pred%w(Imax,j,k)
        End do
      End do
      Do i = 1,Imax
        Do k = 1,Kmax
      ! wall boundary
          Pred%v(i,1-jght,k) = 0.d0 ! - Pred%v(i,jbeg+jght)
          Pred%u(i,1-jght,k) = Pred%u(i,1,k)
          Pred%w(i,1-jght,k) = Pred%w(i,1,k)
      ! wall boundary
          Pred%v(i,Jmax,k) = 0.d0 !TVar%Uint/TVar%Uref
          Pred%v(i,Jmax+jght,k) = 0.d0-Pred%v(i,Jmax-1,k)
          Pred%u(i,Jmax+jght,k) = Pred%u(i,Jmax,k)
          Pred%w(i,Jmax+jght,k) = Pred%w(i,Jmax,k)
        End do
      End do
      Do i = 1,Imax
        Do j = 1,Jmax
      ! Wall boundary
          Pred%u(i,j,1-kght) = Pred%u(i,j,1)
          Pred%v(i,j,1-kght) = Pred%v(i,j,1)
          Pred%w(i,j,1-kght) = 0.d0
      ! Wall boundary
          Pred%u(i,j,Kmax+kght) = Pred%u(i,j,Kmax)
          Pred%v(i,j,Kmax+kght) = Pred%v(i,j,Kmax)
          Pred%w(i,j,Kmax) = 0.d0
          Pred%w(i,j,Kmax+kght) = 0.d0-Pred%w(i,j,Kmax-1)
        End do
      End do
    End subroutine PredictorVelocityBoundaryCondition

    Subroutine PredictorVelocityInternalCellCondition(Pred,UCell,VCell,WCell)
      Implicit none
      Type(Predictor),intent(inout):: Pred
      Type(Cell),intent(in):: UCell,VCell,WCell
      Integer(kind=it4b):: i,j,k
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            If(UCell%Cell_Type(i,j,k)==2) then
              Pred%u(i,j,k) = 0.d0
            End if
            If(VCell%Cell_Type(i,j,k)==2) then
              Pred%v(i,j,k) = 0.d0
            End if
            If(WCell%Cell_Type(i,j,k)==2) then
              Pred%w(i,j,k) = 0.d0
            End if
          End do
        End do
      End do
    End subroutine
End module PredictorUV

Module PredictorUV
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE Matrix
    USE Printresult
    USE MPI
    use BoundaryInterface
    use BoundaryFunction
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

    contains

    subroutine PredictorUVW(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,  &
                    FluxDivOld,TVar_n,TVar,BCu,BCv,BCw,PU,PV,PW,Pred,dt,itt)
      implicit none
      integer(kind=it8b),intent(in)        :: itt
      real(kind=dp),intent(in)         :: dt
      type(Grid),intent(in)        :: PGrid,UGrid,VGrid,WGrid
      type(Cell),intent(inout)         :: PCell,UCell,VCell,WCell
      real(kind=dp),dimension(:,:,:,:),allocatable,intent(inout) :: FluxDivOld
      type(Variables),intent(in),target      :: TVar
      type(Variables),intent(in)             :: TVar_n
      type(BCBase),intent(inout)             :: BCu,BCv,BCw
      type(Predictor),intent(inout)          :: Pred
      type(PoissonCoefficient),intent(inout) :: PU,PV,PW
      type(Variables)          :: TVart
      integer(kind=it4b)         :: i,j,k,ii,jj,kk
      integer*8 :: A,parcsr_A,b,par_b,x,par_x,solver,precond

      ! Convective flux and Diffusive flux
      real(kind=dp),dimension(:,:,:,:),allocatable :: CFEW,CFNS,CFTB,          &
                                     DFEW,DFNS,DFTB,EDFEW,EDFNS,EDFTB,FluxDiv
      real(kind=dp),dimension(:,:,:),allocatable   :: UFric,VFric,WFric,       &
                                           UWE,USN,UBT,VWE,VSN,VBT,WWE,WSN,WBT
      integer(kind=it4b)         :: num_iterations
      real(kind=dp)          :: final_res_norm,MaxU
      real(kind=dp)                :: Fe,Fw,Fn,Fs,Ft,Fb,Fluxn0

      allocate(CFEW(0:Imax+1,Jmax,Kmax,3))
      allocate(CFNS(Imax,0:Jmax+1,Kmax,3))
      allocate(CFTB(Imax,Jmax,0:Kmax+1,3))
      allocate(DFEW(0:Imax+1,Jmax,Kmax,3))
      allocate(DFNS(Imax,0:Jmax+1,Kmax,3))
      allocate(DFTB(Imax,Jmax,0:Kmax+1,3))
      allocate(EDFEW(0:Imax+1,Jmax,Kmax,3))
      allocate(EDFNS(Imax,0:Jmax+1,Kmax,3))
      allocate(EDFTB(Imax,Jmax,0:Kmax+1,3))
      allocate(FluxDiv(Imax,Jmax,Kmax,3))
      allocate(UFric(Imax,Jmax,Kmax))
      allocate(VFric(Imax,Jmax,Kmax))
      allocate(WFric(Imax,Jmax,Kmax))
      allocate(UWE(Jmax,Kmax,2))
      allocate(USN(Imax,Kmax,2))
      allocate(UBT(Imax,Jmax,2))
      allocate(VWE(Jmax,Kmax,2))
      allocate(VSN(Imax,Kmax,2))
      allocate(VBT(Imax,Jmax,2))
      allocate(WWE(Jmax,Kmax,2))
      allocate(WSN(Imax,Kmax,2))
      allocate(WBT(Imax,Jmax,2))
      u => TVar%u
      v => TVar%v
      w => TVar%w
      p => TVar%p

    ! Step 1: Calculate the convective coefficient
    !  print*, 'Boundary value here predictUVW 85'
    !  print*, BCv%VarS(1,13)
      Call ModifiedConvectiveFluxFirstOrder(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,   &
                                        WCell,CFEW,1,0,0)
      Call ModifiedConvectiveFluxFirstOrder(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,   &
                                        WCell,CFNS,0,1,0)
      Call ModifiedConvectiveFluxFirstOrder(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,   &
                                        WCell,CFTB,0,0,1)

    ! Step 2: Calculate the diffusive coefficient

      EDFEW = 0.d0
      EDFNS = 0.d0
      EDFTB = 0.d0
      FluxDiv = 0.d0

      Call DiffusiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,      &
                                         DFEW,EDFEW,nuw,1,0,0)
      Call DiffusiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,      &
                                         DFNS,EDFNS,nuw,0,1,0)
      Call DiffusiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,      &
                                         DFTB,EDFTB,nuw,0,0,1)

    ! Step 3: Calculate source term coefficient as wall function

      Do i = 1,Imax-1
        Do j = 1,Jmax
          Do k = 1,Kmax
          ! U Cell
            if(UCell%Cell_Type(i,j,k)/=2) then
              Fe = CFEW(i+1,j,k,1)
              Fw = CFEW(i,j,k,1)
              Fn = CFNS(i,j+1,k,1)
              Fs = CFNS(i,j,k,1)
              Ft = CFTB(i,j,k+1,1)
              Fb = CFTB(i,j,k,1)
              Fluxn0 = Fe-Fw+Fn-Fs+Ft-Fb
           !   if(itt==1) then
                FluxDiv(i,j,k,1)=Fluxn0
           !   else
           !     FluxDiv(i,j,k,1)=1.5d0*Fluxn0-0.5d0*FluxDivOld(i,j,k,1)
           !     FluxDivOld(i,j,k,1)=Fluxn0
           !   end if
            end if
            if(UCell%Cell_Type(i,j,k)/=2) then
              UFric(i,j,k)=(UCell%vofL(i,j,k)/UCell%vof(i,j,k)*nuw/nuref+      &
                     (1.d0-UCell%vofL(i,j,k)/UCell%vof(i,j,k))*nua/nuref)/Rey* &
                     UCell%WlLh(i,j,k)/UCell%delh(i,j,k)
            else
              UFric(i,j,k) = 0.d0
            end if
            if(isnan(Fluxdiv(i,j,k,1)).or.dabs(Fluxdiv(i,j,k,1))>1.d5) then
              print*,i,j,k,Pred%u(i,j,k)
              print*, fluxdiv(i,j,k,1),i,j,k
              print*, fe,fw,fn,fs,ft,fb
              pause 'predictoruv 114'
            end if
          end do
        end do
      end do
      do i = 1,Imax
        do j = 1,Jmax-1
          do k = 1,Kmax
            if(VCell%Cell_Type(i,j,k)/=2) then ! for VCell
              Fe = CFEW(i+1,j,k,2)
              Fw = CFEW(i,j,k,2)
              Fn = CFNS(i,j+1,k,2)
              Fs = CFNS(i,j,k,2)
              Ft = CFTB(i,j,k+1,2)
              Fb = CFTB(i,j,k,2)
              Fluxn0=Fe-Fw+Fn-Fs+Ft-Fb
           !   if(itt==1) then
                FluxDiv(i,j,k,2)=Fluxn0
           !   else
           !     FluxDiv(i,j,k,2)=1.5d0*Fluxn0-0.5d0*FluxDivOld(i,j,k,2)
           !     FluxDivOld(i,j,k,2)=Fluxn0
           !   end if
            end if
            ! V Cell
            if(VCell%Cell_Type(i,j,k)/=2) then
              VFric(i,j,k)=(VCell%vofL(i,j,k)/VCell%vof(i,j,k)*nuw/nuref+       &
                     (1.d0-VCell%vofL(i,j,k)/VCell%vof(i,j,k))*nua/nuref)/Rey* &
                     VCell%WlLh(i,j,k)/VCell%delh(i,j,k)
            else
              VFric(i,j,k)=0.d0
            end if
            if(isnan(Fluxdiv(i,j,k,2)).or.dabs(Fluxdiv(i,j,k,2))>1.d5) then
              print*,i,j,k,pred%v(i,j,k)
              print*,Fe,Fw,Fn,Fs,Ft,Fb
              pause 'predictoruv 143'
            end if
          end do
        end do
      end do

      do i = 1,Imax
        do j = 1,Jmax
          do k = 1,Kmax-1
            if(WCell%Cell_Type(i,j,k)/=2) then ! for WCell
              Fe=CFEW(i+1,j,k,3)
              Fw=CFEW(i,j,k,3)
              Fn=CFNS(i,j+1,k,3)
              Fs=CFNS(i,j,k,3)
              Ft=CFTB(i,j,k+1,3)
              Fb=CFTB(i,j,k,3)
              Fluxn0=Fe-Fw+Fn-Fs+Ft-Fb
            !  if(itt==1) then
                FluxDiv(i,j,k,3)=Fluxn0
            !  else
            !    FluxDiv(i,j,k,3)=1.5d0*Fluxn0-0.5d0*FluxDivOld(i,j,k,3)
             ! end if
            end if
          ! W Cell
            if(WCell%Cell_Type(i,j,k)/=2) then
              WFric(i,j,k)=(WCell%vofL(i,j,k)/WCell%vof(i,j,k)*nuw/nuref+      &
                   (1.d0-WCell%vofL(i,j,k)/WCell%vof(i,j,k))*nua/nuref)/Rey*   &
                     WCell%WlLh(i,j,k)/WCell%delh(i,j,k)
            else
              WFric(i,j,k) = 0.d0
            end if
            if(isnan(Fluxdiv(i,j,k,3)).or.dabs(Fluxdiv(i,j,k,3))>1.d5) then
              print*,i,j,k,pred%w(i,j,k)
              print*,Fe,Fw,Fn,Fs,Ft,Fb
              pause 'predictoruv 167'
            end if
          end do
        end do
      end do
    ! Solving for UCell
      call SetBasicSolver(solver,precond)
    ! call SetBasicSolver(solver=solver,ierr=ierr)
      call SetMatrix(A,parcsr_A,UGrid,UCell,BCu,DFEW,DFNS,DFTB,            &
                       EDFEW,EDFNS,EDFTB,UFric,PU,UWE,USN,UBT,dt,1,0,0)
      call SetVectors(b,x,par_b,par_x,UGrid,UCell,BCu,PU,UWE,USN,UBT,           &
                                            FluxDiv(:,:,:,1),TVar_n,dt,1,0,0)
      call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
      call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    ! Run info - needed logging turned on
      call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      call DeltaGetValues(x,UCell,Pred%u,1,0,0)
      call HYPRE_IJMatrixDestroy(A,ierr)
      call HYPRE_IJVectorDestroy(b,ierr)
      call HYPRE_IJVectorDestroy(x,ierr)
      call HYPRE_BoomerAMGDestroy(precond,ierr)
      call HYPRE_ParCSRPCGDestroy(solver,ierr)
    ! For VCell
      call SetBasicSolver(solver,precond)
      call SetMatrix(A,parcsr_A,VGrid,VCell,BCv,DFEW,DFNS,DFTB,          &
                       EDFEW,EDFNS,EDFTB,VFric,PV,VWE,VSN,VBT,dt,0,1,0)
      call SetVectors(b,x,par_b,par_x,VGrid,VCell,BCv,PV,VWE,VSN,VBT,           &
                                      FluxDiv(:,:,:,2),TVar_n,dt,0,1,0)
      call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
      call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    ! Run info - needed logging turned on
      call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      call DeltaGetValues(x,VCell,Pred%v,0,1,0)
      call HYPRE_IJMatrixDestroy(A,ierr)
      call HYPRE_IJVectorDestroy(b,ierr)
      call HYPRE_IJVectorDestroy(x,ierr)
      call HYPRE_BoomerAMGDestroy(precond,ierr)
      call HYPRE_ParCSRPCGDestroy(solver,ierr)
    ! For WCell
      call SetBasicSolver(solver,precond)
      call SetMatrix(A,parcsr_A,WGrid,WCell,BCw,DFEW,DFNS,DFTB,          &
                       EDFEW,EDFNS,EDFTB,WFric,PW,WWE,WSN,WBT,dt,0,0,1)
      call SetVectors(b,x,par_b,par_x,WGrid,WCell,BCw,PW,WWE,WSN,WBT,           &
                                      FluxDiv(:,:,:,3),TVar_n,dt,0,0,1)
      call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr) 
      call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    ! Run info - needed logging turned on
      call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      call DeltaGetValues(x,WCell,Pred%w,0,0,1)
      call HYPRE_IJMatrixDestroy(A,ierr)
      call HYPRE_IJVectorDestroy(b,ierr)
      call HYPRE_IJVectorDestroy(x,ierr)
      call HYPRE_BoomerAMGDestroy(precond,ierr)
      call HYPRE_ParCSRPCGDestroy(solver,ierr)

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
      Call PredictorVelocityBoundaryCondition(Pred,BCu,BCv,BCw)
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
        Call HYPRE_ParCSRPCGSetTol(solver,1.0d-30,ierr)
        Call HYPRE_ParCSRPCGSetTwoNorm(solver,0,ierr)
!        Call HYPRE_ParCSRPCGSetPrintLevel(solver,2,ierr)
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
          Call HYPRE_BoomerAMGSetMaxIter(precond,10,ierr)
!        set amg as the pcg preconditioner
!         precond_id = 2
          Call HYPRE_ParCSRPCGSetPrecond(solver,2,precond,ierr)
        End if
    End subroutine SetBasicSolver

    Subroutine SetMatrix(A,parcsr_A,TGrid,TCell,BC,DFEW,DFNS,DFTB,EDFEW,EDFNS,    &
                                    EDFTB,Fric,PUVW,CWE,CSN,CBT,dt,iu,iv,iw)
        Implicit none
        Integer*8,intent(inout)  :: A,parcsr_A
        Type(Grid),intent(in)  :: TGrid
        Type(Cell),intent(in)  :: TCell
        type(BCBase),intent(in)  :: BC
        Real(kind=dp),intent(in) :: dt
        Real(kind=dp),dimension(:,:,:,:),allocatable,intent(in):: DFEW,DFNS,   &
                                                     DFTB,EDFEW,EDFNS,EDFTB
        Real(kind=dp),dimension(:,:,:),allocatable,intent(in):: Fric
        Type(PoissonCoefficient),intent(inout):: PUVW
        Real(kind=dp),dimension(:,:,:),allocatable,intent(inout):: CWE,CSN,CBT
        Integer,intent(in):: iu,iv,iw
        Integer(kind=it4b):: nnz,ictr,ilower,iupper,cols(0:6)
        Integer(kind=it4b):: i,j,k,ii,jj,kk
        Real(kind=dp) :: Fep,Fem,Fwp,Fwm,Fnp,Fnm,Fsp,Fsm,Ftp,Ftm,Fbp,Fbm
        Real(kind=dp) :: aP,aE,aW,aN,aS,aT,aB,Dfe,Dfw,Dfn,Dfs,Dft,Dfb,Sp
        real(kind=dp) :: CoefEW,CoefNS,CoefTB
        Real(kind=dp):: values(0:6),MaxDiff
        ilower = 0
        iupper = TCell%ExtCell
        CoefEW=dble(iu+iv/2.d0+iw/2.d0)
        CoefNS=dble(iu/2.d0+iv+iw/2.d0)
        CoefTB=dble(iu/2.d0+iv/2.d0+iw)
      ! Create and Set up matrix
        Call HYPRE_IJMatrixCreate(MPI_COMM_WORLD,ilower,iupper,ilower,iupper,  &
                                                                        A,ierr)
        Call HYPRE_IJMatrixSetObjectType(A,HYPRE_PARCSR,ierr)
        Call HYPRE_IJMatrixInitialize(A,ierr)
        MaxDiff=0.d0

        Do i = 1,Imax-iu
          Do j = 1,Jmax-iv
            Do k = 1,Kmax-iw
            !  If(i==1) j = Jmax-iv-j
              PUVW%dp(i,j,k)=0.d0
              If(TCell%Cell_Type(i,j,k)/=2) then
              !  Dfe=0.d0;Dfw=0.d0;Dfn=0.d0;Dfs=0.d0;Dft=0.d0;Dfb=0.d0
              !  Fep=0.d0;Fem=0.d0;Fwp=0.d0;Fwm=0.d0;Fnp=0.d0;Fnm=0.d0
              !  Fsp=0.d0;Fsm=0.d0;Ftp=0.d0;Ftm=0.d0;Fbp=0.d0;Fbm=0.d0

              ! For UCell,VCell,WCell
                Dfe = DFEW(i+1,j,k,iu+2*iv+3*iw)
                Dfw = DFEW(i,j,k,iu+2*iv+3*iw)
                Dfn = DFNS(i,j+1,k,iu+2*iv+3*iw)
                Dfs = DFNS(i,j,k,iu+2*iv+3*iw)
                Dft = DFTB(i,j,k+1,iu+2*iv+3*iw)
                Dfb = DFTB(i,j,k,iu+2*iv+3*iw)

                if(dabs(MaxDiff)<dabs(Dfe)) then
                  MaxDiff = DFEW(i+1,j,k,iu+2*iv+3*iw)
                  ii = i
                  jj = j
                  kk = k
                end if

                Fep = EDFEW(i+1,j,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaE(i,j,k))
                Fem = EDFEW(i+1,j,k,iu+2*iv+3*iw)*TCell%EtaE(i,j,k)
                Fwp = EDFEW(i,j,k,iu+2*iv+3*iw)*TCell%EtaE(i-1,j,k)
                Fwm = EDFEW(i,j,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaE(i-1,j,k))
                Fnp = EDFNS(i,j+1,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaN(i,j,k))
                Fnm = EDFNS(i,j+1,k,iu+2*iv+3*iw)*TCell%EtaN(i,j,k)
                Fsp = EDFNS(i,j,k,iu+2*iv+3*iw)*TCell%EtaN(i,j-1,k)
                Fsm = EDFNS(i,j,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaN(i,j-1,k))
                Ftp = EDFTB(i,j,k+1,iu+2*iv+3*iw)*(1.d0-TCell%EtaT(i,j,k))
                Ftm = EDFTB(i,j,k+1,iu+2*iv+3*iw)*TCell%EtaT(i,j,k)
                Fbp = EDFTB(i,j,k,iu+2*iv+3*iw)*TCell%EtaT(i,j,k-1)
                Fbm = EDFTB(i,j,k,iu+2*iv+3*iw)*(1.d0-TCell%EtaT(i,j,k-1))
                ! At the boundary
                if(i==1) then
                  if(BC%flag(1)==1) then
                    Dfw=0.d0
                    fwp=0.d0
                    fwm=0.d0
                  else
                    Dfw=Dfw*CoefEW
                    fwp=fwp*CoefEW
                    fwm=fwm*CoefEW
                  end if
                end if

                if(i==Imax-iu) then
                  if(BC%flag(2)==1) then
                    Dfe=0.d0
                    fep=0.d0
                    fem=0.d0
                  else
                    Dfe=Dfe*CoefEW
                    fep=fep*CoefEW
                    fem=fem*CoefEW
                  end if
                end if

                if(j==1) then
                  if(BC%flag(3)==1) then
                    Dfs=0.d0
                    fsp=0.d0
                    fsm=0.d0
                  else
                    Dfs=Dfs*CoefNS
                    fsp=fsp*CoefNS
                    fsm=fsm*CoefNS
                  end if
                end if

                if(j==Jmax-iv) then
                  if(BC%flag(4)==1) then
                    Dfn=0.d0
                    fnp=0.d0
                    fnm=0.d0
                  else
                    Dfn=Dfn*CoefNS
                    fnp=fnp*CoefNS
                    fnm=fnm*CoefNS
                  end if
                end if

                if(k==1) then
                  if(BC%flag(5)==1) then
                    Dfb=0.d0
                    fbp=0.d0
                    fbm=0.d0
                  else
                    Dfb=Dfb*CoefTB
                    fbp=fbp*CoefTB
                    fbm=fbm*CoefTB
                  end if
                end if

                if(k==Kmax-iw) then
                  if(BC%flag(6)==1) then
                    Dft=0.d0
                    ftp=0.d0
                    ftm=0.d0
                  else
                    Dft=Dft*CoefTB
                    ftp=ftp*CoefTB
                    ftm=ftm*CoefTB
                  end if
                end if

                aE = Dfe;aW = Dfw;aN = Dfn;aS = Dfs;aT = Dft;aB = Dfb

                if(i==1)CWE(j,k,1) = aW
                if(i==Imax-iu)CWE(j,k,2) = aE
                if(j==1)CSN(i,k,1) = aS
                if(j==Jmax-iv)CSN(i,k,2) = aN
                if(k==1)CBT(i,j,1) = aB
                if(k==Kmax-iw)CBT(i,j,2) = aT

                Sp = Fric(i,j,k)
                aP = TCell%vof(i,j,k)*TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*         &
                                      TGrid%dz(i,j,k)/dt+aE+aW+aN+aS+aT+aB+Sp

                aP = aP+Fep-Fwp+Fnp-Fsp+Ftp-Fbp

                aE = aE-Fem
                aW = aW+Fwm
                aN = aN-Fnm
                aS = aS+Fsm
                aT = aT-Ftm
                aB = aB+Fbm

                PUVW%Dp(i,j,k) = 1.d0/(aP-aE-aW-aN-aS-aT-aB)

                ictr = Tcell%Posnu(i,j,k)
                nnz = 0
                values = 0.d0
                cols = 0

                if(isnan(PUVW%Dp(i,j,k)).or.isnan(ap).or.dabs(ap)>1.d10.or.dabs(PUVW%Dp(i,j,k))>1.d10) then
                  print*, aP
                  print*, aE
                  print*, aW
                  print*, aN
                  print*, aS
                  print*, aT
                  print*, aB
                  pause 'Test problem with nan in velocity coefficient'
                end if
                if(dabs(aw)>1.d10.or.dabs(ae)>1.d10.or.dabs(as)>1.d10.or. &
                   dabs(an)>1.d10.or.dabs(ab)>1.d10.or.dabs(at)>1.d10)then
                  print*, 'Test problem with coefficient'
                  print*, aw,ae
                  print*, as,an
                  print*, ab,at
                end if   
                if(dabs(aw)+dabs(ae)+dabs(as)+dabs(an)+dabs(ab)+dabs(at)>dabs(ap)) then
                  print*,ap,aw,ae,as,an,ab,sp
                  print*,fep,fem,TCell%EEArea(i,j,k),TCell%EtaE(i,j,k)
                  print*,fwp,fwm,TCell%EEArea(i-1,j,k),TCell%EtaE(i-1,j,k)
                  print*,fnp,fnm,TCell%NEArea(i,j,k),TCell%EtaN(i,j,k)
                  print*,fsp,fsm,TCell%NEArea(i,j-1,k),TCell%EtaN(i,j-1,k)
                  print*,ftp,ftm,TCell%TEArea(i,j,k),TCell%EtaT(i,j,k)
                  print*,fbp,fbm,TCEll%TEArea(i,j,k-1),TCell%EtaT(i,j,k-1)
                  print*,i,j,k,iu,iv,iw
                  pause 'Predictuvw 416'
                end if
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
                do ii=0,nnz-1
                  if(isnan(values(ii)).or.dabs(values(ii))>1.d5)then
                    print*, 'something wroing with computing coefficient'
                    print*, values(nnz),cols(ii)
                    print*, ii
                  end if  
                end do
                call HYPRE_IJMatrixSetValues(A,1,nnz,ictr,cols,values,ierr)
              end if
            end do
          end do
        end do
     !   print*,'Test the maximum diffusive coefficient'
     !   print*,MaxDiff
     !   print*, ii,jj,kk
        call HYPRE_IJMatrixAssemble(A,ierr)
        call HYPRE_IJMatrixGetObject(A,parcsr_A,ierr)
    end subroutine SetMatrix

    Subroutine SetVectors(b,x,par_b,par_x,TGrid,TCell,BC,PUVW,CWE,CSN,CBT,IJKFlux,&
                                                            TVar_n,dt,iu,iv,iw)
        Integer*8:: b,x,par_b,par_x
        Integer,intent(in)               :: iu,iv,iw
        Type(Grid),intent(in)            :: TGrid
        Type(Cell),intent(in)            :: TCell
        type(BCBase),intent(in)          :: BC
        Type(PoissonCoefficient),intent(inout) :: PUVW
        Type(Variables),intent(in)         :: TVar_n
        Real(dp),dimension(:,:,:),allocatable,intent(in) :: CWE,CSN,CBT
        Real(dp),dimension(:,:,:),intent(in) :: IJKFlux
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
        MaxVect=0.d0
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
        Do i = 1,Imax-iu
          Do j = 1,Jmax-iv
            Do k = 1,Kmax-iw
              If(TCell%Cell_Type(i,j,k)/=2) then
                ictr = TCell%PosNu(i,j,k)
                rhs(ictr)=(iu*TVar_n%u(i,j,k)+iv*TVar_n%v(i,j,k)+              &
                           iw*TVar_n%w(i,j,k))*TCell%vof(i,j,k)*               &
                           TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k)/dt+ &
                           1.d0/Fr**2.d0*TCell%vof(i,j,k)*           &
                           TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k)*    &
                           (gx*dble(iu)+gy*dble(iv)+gz*dble(iw))/g
                If(TCell%MoExCell(i,j,k)/=1) then
                !  rhs(ictr) = rhs(ictr)-TCell%vof(i,j,k)*TGrid%dx(i,j,k)*      &
                !             TGrid%dy(i,j,k)*TGrid%dz(i,j,k)*(p(i+iu,j+iv,k+iw)&
                !             -p(i,j,k))/(dble(iu)*TGrid%dx(i,j,k)+dble(iv)*    &
                !                      TGrid%dy(i,j,k)+dble(iw)*TGrid%dz(i,j,k))
                  PUVW%Dp(i,j,k) = TCell%vof(i,j,k)*TGrid%dx(i,j,k)*           &
                                   TGrid%dy(i,j,k)*TGrid%dz(i,j,k)*PUVW%Dp(i,j,k)
                Else
                  ii = TCell%MsCe(i,j,k,1)
                  jj = TCell%MsCe(i,j,k,2)
                  kk = TCell%MsCe(i,j,k,3)
              !    rhs(ictr) = rhs(ictr)-TCell%vof(i,j,k)*TGrid%dx(i,j,k)*      &
              !            TGrid%dy(i,j,k)*TGrid%dz(i,j,k)*(p(ii+iu,jj+iv,kk+iw)&
              !            -p(ii,jj,kk))/(dble(iu)*TGrid%dx(i,j,k)+dble(iv)*    &
              !                       TGrid%dy(i,j,k)+dble(iw)*TGrid%dz(i,j,k))
                  PUVW%Dp(i,j,k) = TCell%vof(i,j,k)*TGrid%dx(i,j,k)*           &
                                   TGrid%dy(i,j,k)*TGrid%dz(i,j,k)*PUVW%Dp(i,j,k)
                End if
                
                rhs(ictr) = rhs(ictr)-IJKFlux(i,j,k)
                if(i==1)rhs(ictr)=rhs(ictr)+CWE(j,k,1)*BC%VarW(j,k)
                if(i==Imax-iu)rhs(ictr)=rhs(ictr)+CWE(j,k,2)*BC%VarE(j,k)
                if(j==1)rhs(ictr)=rhs(ictr)+CSN(i,k,1)*BC%VarS(i,k)
                if(j==Jmax-iv)rhs(ictr)=rhs(ictr)+CSN(i,k,2)*BC%VarN(i,k)
                if(k==1)rhs(ictr)=rhs(ictr)+CBT(i,j,1)*BC%VarB(i,j)
                if(k==Kmax-iw)rhs(ictr)=rhs(ictr)+CBT(i,j,2)*BC%VarT(i,j)
                xval(ictr) = 0.d0
                if(isnan(PUVW%Dp(i,j,k)).or.isnan(rhs(ictr)).or.dabs(rhs(ictr))>1.d10) then
                  print*,TCell%vof(i,j,k)*TGrid%dx(i,j,k)*           &
                         TGrid%dy(i,j,k)*TGrid%dz(i,j,k)/(dble(iu)*  &
                         TGrid%dx(i,j,k)+dble(iv)*TGrid%dy(i,j,k)+   &
                         dble(iw)*TGrid%dz(i,j,k))
                  print*,i,j,k
                  print*,iu,iv,iw
                  print*, CWE(j,k,1),BC%VarW(j,k)
                  print*, CWE(j,k,2),BC%VarE(j,k)
                  print*, CSN(i,k,1),BC%VarS(i,k)
                  print*, CSN(i,k,2),BC%VarN(i,k)
                  print*, CBT(i,j,1),BC%VarB(i,j)
                  print*, CBT(i,j,2),BC%VarT(i,j)
                  PRINT*,IJKFlux(i,j,k)
                  pause 'Set Vector 557'
                end if
                rows(ictr) = ilower+ictr
              end if
            end do
          end do
        end do
        call HYPRE_IJVectorSetValues(b,local_size,rows,rhs,ierr)
        call HYPRE_IJVectorSetValues(x,local_size,rows,xval,ierr)
        call HYPRE_IJVectorAssemble(b,ierr)
        call HYPRE_IJVectorAssemble(x,ierr)
  ! get the x and b objects
        call HYPRE_IJVectorGetObject(b,par_b,ierr)
        call HYPRE_IJVectorGetObject(x,par_x,ierr)
        deallocate(rhs)
        deallocate(xval)
        deallocate(rows)
    end subroutine SetVectors

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
    Subroutine ModifiedConvectiveFluxFirstOrder(PGrid,UGrid,VGrid,WGrid,PCell,UCell,     &
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
              uwp=0.5d0*(uw+dabs(uw))
              uwn=0.5d0*(uw-dabs(uw))
            ! For UCell both for convective velocity and scalar velocity
              If(i==1) then
                Flux(i,j,k,1)=(uwn*u(i,j,k)+uwp*u(i-1,j,k))*           &
                             UCell%EEArea(i,j,k)*UGrid%dy(i,j,k)*UGrid%dz(i,j,k)
              Elseif(i>=Imax) then
                Flux(i,j,k,1)=(uwn*u(i,j,k)+uwp*u(i-1,j,k))*           &
                       UCell%EEArea(i-1,j,k)*UGrid%dy(i-1,j,k)*UGrid%dz(i-1,j,k)
              Else
                eta=UCell%EtaE(i-1,j,k)
                uw=(1.d0-eta)*u(i-1,j,k)+eta*u(i,j,k)
              ! central second order
                uwp=0.5d0*(uw+dabs(uw))
                uwn=0.5d0*(uw-dabs(uw))
                
                Flux(i,j,k,1)=(uwn*u(i,j,k)+uwp*u(i-1,j,k))*UCell%AlE(i-1,j,k)*&
                      UCell%EEArea(i-1,j,k)*UGrid%dy(i-1,j,k)*UGrid%dz(i-1,j,k)
              End if
              ! Convective velocity: u, scalar advective : v
              uw=0.5d0*(u(i-1,j+1,k)+u(i-1,j,k))
              uwp=0.5d0*(uw+dabs(uw))
              uwn=0.5d0*(uw-dabs(uw))
              If(i>Imax) then
                Flux(i,j,k,2)=(uwp*v(i-1,j,k)+uwn*v(i,j,k))*                   &
                       VCell%EEArea(i-1,j,k)*VGrid%dy(i-1,j,k)*VGrid%dz(i-1,j,k)
              Elseif(i==1) then
                Flux(i,j,k,2)=(uwp*v(i-1,j,k)+uwn*v(i,j,k))*                   &
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
                uwp=0.5d0*(uw+dabs(uw))
                uwn=0.5d0*(uw-dabs(uw))
                vw = (1.d0-eta)*v(i-1,j,k)+eta*v(i,j,k)
                Flux(i,j,k,2)=(uwp*v(i-1,j,k)+uwn*v(i,j,k))*           &
                     VCell%EEArea(i-1,j,k)*VGrid%dy(i-1,j,k)*VGrid%dz(i-1,j,k)
              End if
            ! Convective velocity: u, scalar advective: w
              uw = 0.5d0*(u(i-1,j,k+1)+u(i-1,j,k))
              uwp=0.5d0*(uw+dabs(uw))
              uwn=0.5d0*(uw-dabs(uw))
              If(i>Imax) then
                Flux(i,j,k,3)=(uwp*w(i-1,j,k)+uwn*w(i,j,k))*                   &
                       WCell%EEArea(i-1,j,k)*WGrid%dy(i-1,j,k)*WGrid%dz(i-1,j,k)
              Elseif(i==1) then
                Flux(i,j,k,3)=(uwp*w(i-1,j,k)+uwn*w(i,j,k))*                   &
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
                uwp=0.5d0*(uw+dabs(uw))
                uwn=0.5d0*(uw-dabs(uw))
                Flux(i,j,k,3)=(uwp*w(i-1,j,k)+uwn*w(i,j,k))*           &
                               WCell%EEArea(i-1,j,k)*                  &
                               WGrid%dy(i-1,j,k)*WGrid%dz(i-1,j,k)
              End if
            End if
            If(iv==1) then ! Jflux
          ! Convective velocity: v, scalar advective: u
              vs = 0.5d0*(v(i,j-1,k)+v(i+1,j-1,k))
              vsp = 0.5d0*(vs+dabs(vs))
              vsn = 0.5d0*(vs-dabs(vs))
              If(j>Jmax) then
                Flux(i,j,k,1)=(vsn*u(i,j,k)+vsp*u(i,j-1,k))*                   &
                       UCell%NEArea(i,j-1,k)*UGrid%dx(i,j-1,k)*UGrid%dz(i,j-1,k)
              Elseif(j==1) then
                Flux(i,j,k,1)=(vsn*u(i,j,k)+vsp*u(i,j-1,k))*                   &
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
                    delh=dabs(VCell%Cell_Cent(i+1,j-1,k,1)*VCell%nx(i+1,j-1,k) &
                          +VCell%Cell_Cent(i+1,j-1,k,2)*VCell%ny(i+1,j-1,k)+   &
                           VCell%Cell_Cent(i+1,j-1,k,3)*VCell%nz(i+1,j-1,k)+   &
                           VCell%phi(i+1,j-1,k))
                    vs=v(i+1,j-1,k)
                  End if
                  vs=vs*delhec/(delh+epsi)
                Else
                  Sx=VCell%SxE(i,j-1,k)
                  eta=dabs(UCell%FCN(i,j-1,k,1)+0.5d0*UGrid%dx(i,j-1,k)-       &
                                            VCell%Cell_Cent(i,j-1,k,1))/Sx
                  vs=(1.d0-eta)*v(i,j-1,k)+eta*v(i+1,j-1,k)
                End if
                eta=UCell%EtaN(i,j-1,k)
                vsp=0.5d0*(vs+dabs(vs))
                vsn=0.5d0*(vs-dabs(vs))
              !  us=(1.d0-eta)*u(i,j-1,k)+eta*u(i,j,k)
                Flux(i,j,k,1)=(vsn*u(i,j,k)+vsp*u(i,j-1,k))*          &
                     UCell%NEArea(i,j-1,k)*UGrid%dx(i,j-1,k)*UGrid%dz(i,j-1,k)
              End if
              ! Convective velocity: v, scalar convective: v
              vs = 0.5d0*(v(i,j,k)+v(i,j-1,k))
              vsp = 0.5d0*(vs+dabs(vs))
              vsn = 0.5d0*(vs-dabs(vs))
              If(j>=Jmax) then
                Flux(i,j,k,2)=(vsn*v(i,j,k)+vsp*v(i,j-1,k))*           &
                       VCell%NEArea(i,j-1,k)*VGrid%dx(i,j-1,k)*VGrid%dz(i,j-1,k)
              Elseif(j==1) then
                Flux(i,j,k,2)=(vsn*v(i,j,k)+vsp*v(i,j-1,k))*           &
                       VCell%NEArea(i,j,k)*VGrid%dx(i,j,k)*VGrid%dz(i,j,k)
              Else
                eta=VCell%EtaN(i,j-1,k)
                vs=(1.d0-eta)*v(i,j-1,k)+eta*v(i,j,k)
                Flux(i,j,k,2)=(vsn*v(i,j,k)+vsp*v(i,j-1,k))*VCell%AlN(i,j-1,k)*&
                       VCell%NEArea(i,j-1,k)*VGrid%dx(i,j-1,k)*VGrid%dz(i,j-1,k)
              End if
            ! Convective velocity: v, scalar advective: w
              vs=0.5d0*(v(i,j-1,k)+v(i,j-1,k+1))
              If(j>Jmax) then
                Flux(i,j,k,3)=vs*0.5d0*(w(i,j,k)+w(i,j-1,k))*                  &
                       WCell%NEArea(i,j-1,k)*WGrid%dx(i,j-1,k)*WGrid%dz(i,j-1,k)
              Elseif(j==1) then
                Flux(i,j,k,3)=vs*0.5d0*(w(i,j,k)+w(i,j-1,k))*                  &
                       WCell%NEArea(i,j,k)*WGrid%dx(i,j,k)*wGrid%dz(i,j,k)
              Else
                Flux(i,j,k,3)=0.d0
                If(PCell%NEArea(i,j-1,k)<=1.d-5.or.PCell%NEArea(i,j-1,k+1)     &
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
                vsp = 0.5d0*(vs+dabs(vs))
                vsn = 0.5d0*(vs-dabs(vs))
                eta = WCell%EtaN(i,j-1,k)
                ws = (1.d0-eta)*w(i,j-1,k)+eta*w(i,j,k)
                Flux(i,j,k,3)=(vsp*w(i,j-1,k)+vsn*w(i,j,k))*      &
                   WCell%NEArea(i,j-1,k)*WGrid%dx(i,j-1,k)*WGrid%dz(i,j-1,k)
              End if
            End if
            If(iw==1) then
          ! Convective velocity: w, scalar advective: u
              wb = 0.5d0*(w(i,j,k-1)+v(i+1,j,k-1))
              wbp=0.5d0*(wb+dabs(wb))
              wbn=0.5d0*(wb-dabs(wb))
              If(k>Kmax) then
                Flux(i,j,k,1)=(u(i,j,k)*wbn+u(i,j,k-1)*wbp)*                   &
                       UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k-1)*UGrid%dy(i,j,k-1)
              Elseif(k==1) then
                Flux(i,j,k,1)=(u(i,j,k)*wbn+u(i,j,k-1)*wbp)*                   &
                       UCell%TEArea(i,j,k)*UGrid%dx(i,j,k)*UGrid%dy(i,j,k)
              Else
                Flux(i,j,k,1) = 0.d0
                If(PCell%TEArea(i,j,k-1)<=1.d-5.or.PCell%TEArea(i+1,j,k-1)     &
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
                wbp=0.5d0*(wb+dabs(wb))
                wbn=0.5d0*(wb-dabs(wb))
                eta = UCell%EtaT(i,j,k-1)
                ub = (1.d0-eta)*u(i,j,k-1)+eta*u(i,j,k)
                Flux(i,j,k,1)=(wbp*u(i,j,k-1)+wbn*u(i,j,k))*UCell%TEArea(i,j,k-1)*&
                                      UGrid%dx(i,j,k-1)*UGrid%dy(i,j,k-1)
                If(isnan(flux(i,j,k,1))) then
                  print*, UCell%FCT(i,j,k-1,1),UCell%FCT(i,j,k-1,2),UCell%FCT(i,j,k-1,3)
                  print*, PCell%TEArea(i,j,k-1),PCell%TEArea(i+1,j,k-1)
                  print*,delhec1,delhec2,delh,wb,delhec/delh,i,j,k
                  print*,wbp,wbn,UCell%TEArea(i,j,k-1),WCell%SxE(i,j,k-1)
                  pause '++++++++++++++Predictuvw 922'
                End if
             ! Convective velocity: w, scalar advective: v
              End if
              wb = 0.5d0*(w(i,j+1,k-1)+w(i,j,k-1))
              wbp=0.5d0*(wb+dabs(wb))
              wbn=0.5d0*(wb-dabs(wb))
              If(k>Kmax) then
                Flux(i,j,k,2)=(wbp*v(i,j,k-1)+wbn*v(i,j,k))*                   &
                       VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k-1)*VGrid%dy(i,j,k-1)
              Elseif(k==1) then
                Flux(i,j,k,2)=(wbp*v(i,j,k-1)+wbn*v(i,j,k))*                   &
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
                wbp=0.5d0*(wb+dabs(wb))
                wbn=0.5d0*(wb-dabs(wb))
                eta = VCell%EtaT(i,j,k-1)
                vb = (1.d0-eta)*v(i,j,k-1)+eta*v(i,j,k)
                Flux(i,j,k,2)=(wbp*v(i,j,k-1)+wbn*v(i,j,k))*           &
                    VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k-1)*VGrid%dy(i,j,k-1)
             ! Convective velocity: w, scalar advective: w
              End if
              wb = 0.5d0*(w(i,j,k)+w(i,j,k-1))
              wbp=0.5d0*(wb+dabs(wb))
              wbn=0.5d0*(wb-dabs(wb))
              If(k>=Kmax) then
                Flux(i,j,k,3)=(wbp*w(i,j,k-1)+wbn*w(i,j,k))*      &
                    WCell%TEArea(i,j,k-1)*WGrid%dx(i,j,k-1)*WGrid%dy(i,j,k-1)
              Elseif(k==1) then
                Flux(i,j,k,3)=(wbp*w(i,j,k-1)+wbn*w(i,j,k))*      &
                    WCell%TEArea(i,j,k)*WGrid%dx(i,j,k)*WGrid%dy(i,j,k)
              Else
                eta = WCell%EtaT(i,j,k-1)
                wb = (1.d0-eta)*w(i,j,k-1)+eta*w(i,j,k)
                wbp=0.5d0*(wb+dabs(wb))
                wbn=0.5d0*(wb-dabs(wb))
                Flux(i,j,k,3)=((wbp*w(i,j,k-1)+wbn*w(i,j,k))*     &
                      WCell%AlT(i,j,k-1))*WCell%TEArea(i,j,k-1)*    &
                      WGrid%dx(i,j,k-1)*WGrid%dy(i,j,k-1)
              End if
            End if
          End do
        End do
      End do
    End subroutine ModifiedConvectiveFluxFirstOrder
        
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
                  pause '++++++++++++++Predictuvw 922'
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
    
    SUBROUTINE HighOrderConvectiveFluxForXDir(PGrid,UGrid,VGrid,WGrid,UCell,         &
                                              VCell,WCell,vb,flux)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid,WGrid
      TYPE(Cell),INTENT(IN):: UCell,VCell,WCell
      REAL(KIND=dp),INTENT(IN):: vb
      REAL(KIND=dp),DIMENSION(:,:,:,:),allocatable,INTENT(INOUT)::flux
      INTEGER(kind=it4b):: i,j,k,Lim
      REAL(KIND=dp):: ul,ur,vl,vr,wl,wr,alr,uwp,uwn,sx,sy,sz,uw,delhec,delh,eta
      REAL(KIND=dp):: omei,omei1,ri,ri1,tolim,tol
      tol=1.d-24
      Lim=2
      do j=1,Jmax
        do k=1,Kmax
          ul=u(0,j,k);ur=u(1,j,k)
          alr=(ur+ul)
          alr=dmax1(dabs(ur),dabs(ul))
          flux(1,j,k,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*           &
                               UCell%EEArea(1,j,k)*UGrid%dy(1,j,k)*UGrid%dz(1,j,k)
          vl=v(0,j,k);vr=v(1,j,k)
          alr=0.5d0*(u(0,j,k)+u(0,j+1,k))
          flux(1,j,k,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*               &
                               VCell%EEArea(1,j,k)*VGrid%dy(1,j,k)
          wl=w(0,j,k);wr=w(1,j,k)
          alr=0.5d0*(w(0,j,k)+w(0,j+1,k))
          flux(1,j,k,3)=0.5d0*(alr*wr+alr*wl-dabs(alr)*(wr-wl))*         &
                     WCell%EEArea(1,j,k)*WGrid%dy(1,j,k)  
          do i=2,Imax
        ! Calculate threshold for MUSCL
        ! from 'A MUSCL scheme on staggered grids with kinetic-like fluxes
        ! for the barotropic Euler system', Thierry Goundon, Julie Llobell
            tolim=dmin1(2.d0*UGrid%dx(i-1,j,k)/PGrid%dx(i,j,k),                &
                        2.d0*UGrid%dx(i,j,k)/PGrid%dx(i,j,k))
            ri=((u(i-1,j,k)-u(i-2,j,k))/PGrid%dx(i-1,j,k))/                    &
               ((u(i,j,k)-u(i-1,j,k))/PGrid%dx(i,j,k))
            omei=MUSCLLimiter(ri,Lim,tolim)*(u(i,j,k)-u(i-1,j,k))/PGrid%dx(i,j,k)
            ri1=((u(i,j,k)-u(i-1,j,k))/PGrid%dx(i,j,k))/                       &
                ((u(i+1,j,k)-u(i,j,k))/PGrid%dx(min(Imax,i+1),j,k))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(u(i+1,j,k)-u(i,j,k))/           &
                                       PGrid%dx(min(Imax,i+1),j,k)
            ul=u(i-1,j,k)+0.5d0*(PGrid%dx(i,j,k))*omei
            ur=u(i,j,k)-0.5d0*(PGrid%dx(i,j,k))*omei1
            alr=(ur+ul)
            alr=dmax1(dabs(ur),dabs(ul))
            flux(i,j,k,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*         &
                                 UCell%EEArea(i-1,j,k)*UGrid%dy(i,j,k)
            if(i>2.and.i<Imax) then
              if(UCell%vof(i-2,j,k)<1.d0-epsi.or.            &
                 UCell%vof(i-1,j,k)<1.d0-epsi.or.                      &
                 UCell%Vof(i,j,k)<1.d0-epsi.or.              &
                 UCell%vof(i+1,j,k)<1.d0-epsi) then
                eta=UCell%EtaE(i-1,j,k)
                uw=(1-eta)*u(i-1,j,k)+eta*u(i,j,k)
              ! first order for x-direction
              !   if(UCell%vofS(i-1,j)>epsi.or.UCell%VofS(i,j)>epsi) then
                uwp=0.5d0*(uw+dabs(uw))
                uwn=0.5d0*(uw-dabs(uw))
                Flux(i,j,k,1)=(uwp*u(i-1,j,k)+uwn*u(i,j,k))*           &
                        UCell%AlE(i-1,j,k)*UGrid%dy(i,j,k)*UCell%EEArea(i-1,j,k)
              end if
            end if
            
            tolim=dmin1(2.d0*VGrid%dx(i-1,j,k)/UGrid%dx(i-1,j,k),              &
                        2.d0*VGrid%dx(i,j,k)/UGrid%dx(i-1,j,k))
            ri=((v(i-1,j,k)-v(i-2,j,k))/UGrid%dx(max(1,i-2),j,k))/             &
               ((v(i,j,k)-v(i-1,j,k))/UGrid%dx(i-1,j,k))
            omei=MUSCLlimiter(ri,Lim,tolim)*(v(i,j,k)-v(i-1,j,k))/UGrid%dx(i-1,j,k)
            ri1=((v(i,j,k)-v(i-1,j,k))/UGrid%dx(i-1,j,k))/                     &
                ((v(i+1,j,k)-v(i,j,k))/UGrid%dx(i,j,k))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(v(i+1,j,k)-v(i,j,k))/UGrid%dx(i,j,k)
            vl=v(i-1,j,k)+0.5d0*(PGrid%dx(i-1,j,k))*omei
            vr=v(i,j,k)-0.50d0*(PGrid%dx(i,j,k))*omei1
            alr=0.5d0*(u(i-1,j,k)+u(i-1,j+1,k))
        !  alr=dmax1(dabs(u(i-1,j)),dabs(u(i-1,j+1)))
            flux(i,j,k,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*             &
                               VCell%EEArea(i-1,j,k)*VGrid%dy(i,j,k)
            if(i>2.and.i<Imax) then
              if(VCell%vof(i-2,j,k)<1.d0-epsi.or.            &
                 VCell%vof(i-1,j,k)<1.d0-epsi.or.                    &
                 VCell%vof(i,j,k)<1.d0-epsi.or.              &
                 VCell%vof(i+1,j,k)<1.d0-epsi) then
                
                if(VCell%EEArea(i,j,k)<0.5d0) then
                  delhec=dabs(VCell%FCE(i-1,j,k,1)*VCell%nx(i-1,j,k)+          &
                              VCell%FCE(i-1,j,k,2)*VCell%ny(i-1,j,k)+        &
                              VCell%FCE(i-1,j,k,3)*VCell%nz(i-1,j,k)+        &
                              VCell%phi(i-1,j,k))
                  if(UCell%MoExCell(i-1,j+1,k)/=1.and.             &
                     UCell%Vof(i-1,j+1,k)>epsi)then
                    delh=dabs(UCell%Cell_Cent(i-1,j+1,k,1)*UCell%nx(i-1,j+1,k)+&
                              UCell%Cell_Cent(i-1,j+1,k,2)*UCell%ny(i-1,j+1,k)+&
                              UCell%Cell_Cent(i-1,j+1,k,3)*UCell%nz(i-1,j+1,k)+&
                              UCell%phi(i-1,j+1,k))+tol
                    uw=u(i-1,j+1,k)
                  elseif(UCell%MoExCell(i-1,j,k)/=1.and.UCell%Vof(i-1,j,k)>epsi)then
                    delh=dabs(UCell%Cell_Cent(i-1,j,k,1)*UCell%nx(i-1,j,k)+    &
                              UCell%Cell_Cent(i-1,j,k,2)*UCell%ny(i-1,j,k)+    &
                              UCell%Cell_Cent(i-1,j,k,3)*UCell%nz(i-1,j,k)+    &
                              UCell%phi(i-1,j,k))+tol
                     uw=u(i-1,j,k)
                  else
                    delh=delhec+tol
                    uw=0.d0
                  end if
                  uwp=0.5d0*(uw+dabs(uw))
                  uwn=0.5d0*(uw-dabs(uw))
                  Flux(i,j,k,2)=(uwp*v(i-1,j,k)+uwn*v(i,j,k))*delhec/delh*     &
                                 VCell%EEArea(i-1,j,k)*VGrid%dy(i,j,k)
                else
                  Sy=UCell%SyN(i-1,j,k)
                  eta=dabs(VCell%FCE(i-1,j,k,2)+VGrid%dy(i-1,j,k)/2.d0-        &
                           UCell%Cell_Cent(i-1,j,k,2))/Sy
                  if(dabs(eta)>=1.d0) eta=0.5d0
                  uw=(1.d0-eta)*u(i-1,j,k)+eta*u(i-1,j+1,k)
                  uwp=0.5d0*(uw+dabs(uw))
                  uwn=0.5d0*(uw-dabs(uw))
                  Flux(i,j,k,2)=(uwp*v(i-1,j,k)+uwn*v(i,j,k))*           &
                                 VCell%EEArea(i-1,j,k)*VGrid%dy(i,j,k)

                end if
              end if
            end if
            
            tolim=dmin1(2.d0*WGrid%dx(i-1,j,k)/WGrid%dx(i-1,j,k),              &
                        2.d0*WGrid%dx(i,j,k)/WGrid%dx(i-1,j,k))
            ri=((w(i-1,j,k)-w(i-2,j,k))/UGrid%dx(max(1,i-2),j,k))/             &
               ((w(i,j,k)-w(i-1,j,k))/UGrid%dx(i-1,j,k))
            omei=MUSCLlimiter(ri,Lim,tolim)*(w(i,j,k)-w(i-1,j,k))/UGrid%dx(i-1,j,k)
            ri1=((w(i,j,k)-w(i-1,j,k))/UGrid%dx(i-1,j,k))/                     &
                ((w(i+1,j,k)-w(i,j,k))/UGrid%dx(i,j,k))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(w(i+1,j,k)-w(i,j,k))/UGrid%dx(i,j,k)
            wl=w(i-1,j,k)+0.5d0*(PGrid%dx(i-1,j,k))*omei
            wr=w(i,j,k)-0.50d0*(PGrid%dx(i,j,k))*omei1
            alr=0.5d0*(u(i-1,j,k)+u(i-1,j,k+1))
            flux(i,j,k,3)=0.5d0*(alr*wr+alr*wl-dabs(alr)*(wr-wl))*               &
                               WCell%EEArea(i-1,j,k)*WGrid%dy(i,j,k)
            if(i>2.and.i<Imax) then
              if(WCell%vof(i-2,j,k)<1.d0-epsi.or.            &
                 WCell%vof(i-1,j,k)<1.d0-epsi.or.                    &
                 WCell%vof(i,j,k)<1.d0-epsi.or.              &
                 WCell%vof(i+1,j,k)<1.d0-epsi) then
                
                if(WCell%EEArea(i,j,k)<0.5d0) then
                  delhec=dabs(WCell%FCE(i-1,j,k,1)*WCell%nx(i-1,j,k)+          &
                              WCell%FCE(i-1,j,k,2)*WCell%ny(i-1,j,k)+        &
                              WCell%FCE(i-1,j,k,3)*WCell%nz(i-1,j,k)+        &
                              WCell%phi(i-1,j,k))
                  if(UCell%MoExCell(i-1,j,k+1)/=1.and.             &
                     UCell%Vof(i-1,j,k+1)>epsi)then
                    delh=dabs(UCell%Cell_Cent(i-1,j,k+1,1)*UCell%nx(i-1,j,k+1)+&
                              UCell%Cell_Cent(i-1,j,k+1,2)*UCell%ny(i-1,j,k+1)+&
                              UCell%Cell_Cent(i-1,j,k+1,3)*UCell%nz(i-1,j,k+1)+&
                              UCell%phi(i-1,j,k+1))+tol
                    uw=u(i-1,j,k+1)
                  elseif(UCell%MoExCell(i-1,j,k)/=1.and.UCell%Vof(i-1,j,k)>epsi)then
                    delh=dabs(UCell%Cell_Cent(i-1,j,k,1)*UCell%nx(i-1,j,k)+    &
                              UCell%Cell_Cent(i-1,j,k,2)*UCell%ny(i-1,j,k)+    &
                              UCell%Cell_Cent(i-1,j,k,3)*UCell%nz(i-1,j,k)+    &
                              UCell%phi(i-1,j,k))+tol
                     uw=u(i-1,j,k)
                  else
                    delh=delhec+tol
                    uw=0.d0
                  end if
                  uwp=0.5d0*(uw+dabs(uw))
                  uwn=0.5d0*(uw-dabs(uw))
                  Flux(i,j,k,3)=(uwp*w(i-1,j,k)+uwn*w(i,j,k))*delhec/delh*     &
                                 WCell%EEArea(i-1,j,k)*WGrid%dy(i,j,k)
                else
                  Sz=UCell%SzT(i-1,j,k)
                  eta=dabs(WCell%FCE(i-1,j,k,3)+WGrid%dz(i-1,j,k)/2.d0-        &
                           WCell%Cell_Cent(i-1,j,k,3))/Sz
                  if(dabs(eta)>=1.d0) eta=0.5d0
                  uw=(1.d0-eta)*u(i-1,j,k)+eta*u(i-1,j,k+1)
                  uwp=0.5d0*(uw+dabs(uw))
                  uwn=0.5d0*(uw-dabs(uw))
                  Flux(i,j,k,3)=(uwp*w(i-1,j,k)+uwn*w(i,j,k))*           &
                                 WCell%EEArea(i-1,j,k)*WGrid%dy(i,j,k)

                end if
              end if
            end if
          end do  
        end do
        ul=u(Imax,j,k);ur=u(Imax+1,j,k)
        alr=(ur+ul)
        alr=dmax1(dabs(ur),dabs(ul))
        flux(Imax+1,j,k,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*       &
                                     UCell%EEArea(Imax,j,k)*UGrid%dy(Imax,j,k)
        vl=v(Imax,j,k);vr=v(Imax+1,j,k)
        alr=0.5d0*(u(Imax,j,k)+u(Imax,j+1,k))
        flux(Imax+1,j,k,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*           &
                                     VCell%EEArea(Imax,j,k)*VGrid%dy(Imax,j,k)
        wl=w(Imax,j,k);wr=w(Imax+1,j,k)
        alr=0.5d0*(u(Imax,j,k)+u(Imax,j,k+1))
        flux(Imax+1,j,k,3)=0.5d0*(alr*wr+alr*wl-dabs(alr)*(vr-vl))*        &
                                     WCell%EEArea(Imax,j,k)*WGrid%dy(Imax,j,k)  
      end do
    end subroutine HighOrderConvectiveFluxForXDir
    
    subroutine HighOrderConvectiveFluxForYDir(PGrid,UGrid,VGrid,WGrid,UCell,   &
                                                    VCell,WCell,vb,flux)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid,WGrid
      TYPE(Cell),INTENT(IN):: UCell,VCell,WCell
      REAL(KIND=dp),INTENT(IN):: vb
      REAL(KIND=dp),DIMENSION(:,:,:,:),allocatable,INTENT(INOUT)::flux
      INTEGER(kind=it4b):: i,j,k,Lim
      REAL(KIND=dp):: ul,ur,vl,vr,wl,wr,alr,delhec,delh,tol
      REAL(KIND=dp):: omei,omei1,ri,ri1,tolim,Sx,Sy,Sz,vs,vsp,vsn,eta
      Lim=2
      tol=1.d-24
      do i=1,Imax
        do k=1,Kmax
          ul=u(i,0,k);ur=u(i,1,k)
          alr=0.5d0*(v(i,0,k)+v(i+1,0,k))
          flux(i,1,k,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*               &
                        UCell%NEArea(i,1,k)*UGrid%dx(i,1,k)*UGrid%dz(i,1,k)
          vl=v(i,0,k);vr=v(i,1,k)
          alr=(vr+vl)
          alr=dmax1(dabs(vl),dabs(vr))
          flux(i,1,k,2)=0.5d0*(vr**2.d0+vl**2.d0-dabs(alr)*(vr-vl))*           &
                        VCell%NEArea(i,1,k)*VGrid%dx(i,1,k)*VGrid%dz(i,1,k)
          wl=w(i,0,k);wr=w(i,1,k)
          alr=0.5d0*(v(i,0,k)+v(i,0,k+1))
          flux(i,1,k,3)=0.5d0*(alr*wr+alr*wl-dabs(alr)*(wr-wl))*         &
              WCell%NEArea(i,1,k)*WGrid%dx(i,1,k)*WGrid%dz(i,1,k) 
          do j=2,Jmax
            tolim=dmin1(2.d0*UGrid%dy(i,j-1,k)/VGrid%dy(i,j-1,k),              &
                        2.d0*UGrid%dy(i,j,k)/VGrid%dy(i,j-1,k))
            ri=((u(i,j-1,k)-u(i,j-2,k))/VGrid%dy(i,max(1,j-2),k))/             &
               ((u(i,j,k)-u(i,j-1,k))/VGrid%dy(i,j-1,k))
            omei=MUSCLLimiter(ri,Lim,tolim)*(u(i,j,k)-u(i,j-1,k))/VGrid%dy(i,j-1,k)
            ri1=((u(i,j,k)-u(i,j-1,k))/VGrid%dy(i,j-1,k))/                     &
                ((u(i,j+1,k)-u(i,j,k))/VGrid%dy(i,j,k))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(u(i,j+1,k)-u(i,j,k))/VGrid%dy(i,j,k)
            ul=u(i,j-1,k)+0.5d0*(PGrid%dy(i,j-1,k))*omei
            ur=u(i,j,k)-0.5d0*(PGrid%dy(i,j,k))*omei1
            alr=0.5d0*(v(i,j-1,k)+v(i+1,j-1,k))
       !   alr=dmax1(dabs(v(i,j-1)),dabs(v(i+1,j-1)))
            flux(i,j,k,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*             &
                          UCell%NEArea(i,j-1,k)*UGrid%dx(i,j,k)*UGrid%dz(i,j,k)
            if(j>2.and.j<Jmax) then
              if(UCell%vof(i,j-2,k)<1.d0-epsi.or.            &
                 UCell%vof(i,j-1,k)<1.d0-epsi.or.                      &
                 UCell%Vof(i,j,k)<1.d0-epsi.or.              &
                 UCell%vof(i,j+1,k)<1.d0-epsi) then
                if(UCell%NEArea(i,j-1,k)<0.5d0) then
                  delhec=dabs(UCell%FCN(i,j-1,k,1)*UCell%nx(i,j-1,k)+          &
                              UCell%FCN(i,j-1,k,2)*UCell%ny(i,j-1,k)+        &
                              UCell%FCN(i,j-1,k,3)*UCell%nz(i,j-1,k)+          &
                              UCell%phi(i,j-1,k))
                  if(VCell%MoExCell(i,j-1,k)/=1.and.             &
                     VCell%Vof(i,j-1,k)>epsi)then
                    delh=dabs(VCell%Cell_Cent(i,j-1,k,1)*VCell%nx(i,j-1,k)+    &
                              VCell%Cell_Cent(i,j-1,k,2)*VCell%ny(i,j-1,k)+    &
                              VCell%Cell_Cent(i,j-1,k,3)*VCell%nz(i,j-1,k)+    &  
                              VCell%phi(i,j-1,k))+tol
                    vs=v(i,j-1,k)
                  elseif(VCell%MoExCell(i+1,j-1,k)/=1.and.                     &
                         VCell%Vof(i+1,j-1,k)>epsi) then
                    delh=dabs(VCell%Cell_Cent(i+1,j-1,k,1)*VCell%nx(i+1,j-1,k)+&
                            VCell%Cell_Cent(i+1,j-1,k,2)*VCell%ny(i+1,j-1,k)+&
                            VCell%Cell_Cent(i+1,j-1,k,3)*VCell%nz(i+1,j-1,k)+&
                              VCell%phi(i+1,j-1,k))+tol
                    vs=v(i+1,j-1,k)
                  else
                    delh=delhec+tol
                    vs=0.d0
                  end if
                  vs=vb+(vs-vb)*delhec/delh
                  vsp=0.5d0*(vs+dabs(vs))
                  vsn=0.5d0*(vs-dabs(vs))
                  Flux(i,j,k,1)=(vsp*u(i,j-1,k)+vsn*u(i,j,k))*delhec/delh*     &
                         UCell%NEArea(i,j-1,k)*UGrid%dx(i,j,k)*UGrid%dz(i,j,k)
                else
                  Sx=VCell%SxE(i,j-1,k)
                  eta=dabs(UCell%FCN(i,j-1,k,1)+0.5d0*UGrid%dx(i,j-1,k)-       &
                           VCell%Cell_Cent(i,j-1,k,1))/Sx
                  if(dabs(eta)>=1.d0) eta=0.5d0
                  vs=(1.d0-eta)*v(i,j-1,k)+eta*v(i+1,j-1,k)
                  vsp=0.5d0*(vs+dabs(vs))
                  vsn=0.5d0*(vs-dabs(vs))
                  Flux(i,j,k,1)=(vsp*u(i,j-1,k)+vsn*u(i,j,k))*           &
                         UCell%NEArea(i,j-1,k)*UGrid%dx(i,j,k)*UGrid%dz(i,j,k)
                end if
              end if
            end if

            tolim=dmin1(2.d0*VGrid%dy(i,j-1,k)/PGrid%dy(i,j,k),                &
                        2.d0*VGrid%dy(i,j,k)/PGrid%dy(i,j,k))
            ri=((v(i,j-1,k)-v(i,j-2,k))/PGrid%dy(i,j-1,k))/                    &
               ((v(i,j,k)-v(i,j-1,k))/PGrid%dy(i,j,k))
            omei=MUSCLlimiter(ri,Lim,tolim)*(v(i,j,k)-v(i,j-1,k))/PGrid%dy(i,j-1,k)
            ri1=((v(i,j,k)-v(i,j-1,k))/PGrid%dy(i,j,k))/                       &
                ((v(i,j+1,k)-v(i,j,k))/PGrid%dy(i,min(Jmax,j+1),k))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(v(i,j+1,k)-v(i,j,k))/           &
                                              PGrid%dy(i,min(Jmax,j+1),k)
            vl=v(i,j-1,k)+0.5d0*(PGrid%dy(i,j,k))*omei
            vr=v(i,j,k)-0.5d0*(PGrid%dy(i,j,k))*omei1
            alr=(vr+vl)
            alr=dmax1(dabs(vl),dabs(vr))
            flux(i,j,k,2)=0.5d0*(vr*vr+vl*vl-dabs(alr)*(vr-vl))*               &
                           VCell%NEArea(i,j-1,k)*VGrid%dx(i,j,k)*VGrid%dz(i,j,k)
            if(j>2.and.j<Jmax) then
              if(VCell%vof(i,j-2,k)<1.d0-epsi.or.            &
                 VCell%vof(i,j-1,k)<1.d0-epsi.or.                      &
                 VCell%Vof(i,j,k)<1.d0-epsi.or.              &
                 VCell%vof(i,j+1,k)<1.d0-epsi) then
                Sx=VCell%Cell_Cent(i,j,k,1)-VCell%Cell_Cent(i,j-1,k,1)
                Sy=VCell%SyN(i,j-1,k)
                eta=VCell%EtaN(i,j-1,k)
                vs=(1.d0-eta)*v(i,j-1,k)+eta*v(i,j,k)
                vs=(vb+(vs-vb)*VCell%AlN(i,j-1,k))
                vsp=0.5d0*(vs+dabs(vs))
                vsn=0.5d0*(vs-dabs(vs))
                Flux(i,j,k,2)=(vsp*v(i,j-1,k)+vsn*v(i,j,k))*                   &
                           VCell%NEArea(i,j-1,k)*VGrid%dx(i,j,k)*VGrid%dz(i,j,k)
              end if
            end if
            
            tolim=dmin1(2.d0*WGrid%dy(i,j-1,k)/VGrid%dy(i,j-1,k),              &
                        2.d0*WGrid%dy(i,j,k)/VGrid%dy(i,j-1,k))
            ri=((w(i,j-1,k)-w(i,j-2,k))/VGrid%dy(i,max(1,j-2),k))/             &
               ((w(i,j,k)-w(i,j-1,k))/VGrid%dy(i,j-1,k))
            omei=MUSCLLimiter(ri,Lim,tolim)*(w(i,j,k)-w(i,j-1,k))/VGrid%dy(i,j-1,k)
            ri1=((w(i,j,k)-w(i,j-1,k))/VGrid%dy(i,j-1,k))/                     &
                ((w(i,j+1,k)-w(i,j,k))/VGrid%dy(i,j,k))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(w(i,j+1,k)-w(i,j,k))/WGrid%dy(i,j,k)
            ul=w(i,j-1,k)+0.5d0*(PGrid%dy(i,j-1,k))*omei
            ur=w(i,j,k)-0.5d0*(PGrid%dy(i,j,k))*omei1
            alr=0.5d0*(v(i,j-1,k)+v(i,j-1,k+1))
            flux(i,j,k,3)=0.5d0*(alr*wr+alr*wl-dabs(alr)*(wr-wl))*             &
                          WCell%NEArea(i,j-1,k)*WGrid%dx(i,j,k)*WGrid%dz(i,j,k)
            if(j>2.and.j<Jmax) then
              if(WCell%vof(i,j-2,k)<1.d0-epsi.or.            &
                 WCell%vof(i,j-1,k)<1.d0-epsi.or.                      &
                 WCell%Vof(i,j,k)<1.d0-epsi.or.              &
                 WCell%vof(i,j+1,k)<1.d0-epsi) then
                if(WCell%NEArea(i,j-1,k)<0.5d0) then
                  delhec=dabs(WCell%FCN(i,j-1,k,1)*WCell%nx(i,j-1,k)+          &
                              WCell%FCN(i,j-1,k,2)*WCell%ny(i,j-1,k)+        &
                              WCell%FCN(i,j-1,k,3)*WCell%nz(i,j-1,k)+          &
                              WCell%phi(i,j-1,k))
                  if(VCell%MoExCell(i,j-1,k)/=1.and.             &
                     VCell%Vof(i,j-1,k)>epsi)then
                    delh=dabs(VCell%Cell_Cent(i,j-1,k,1)*VCell%nx(i,j-1,k)+    &
                              VCell%Cell_Cent(i,j-1,k,2)*VCell%ny(i,j-1,k)+    &
                              VCell%Cell_Cent(i,j-1,k,3)*VCell%nz(i,j-1,k)+    &  
                              VCell%phi(i,j-1,k))+tol
                    vs=v(i,j-1,k)
                  elseif(VCell%MoExCell(i,j-1,k+1)/=1.and.                     &
                         VCell%Vof(i,j-1,k+1)>epsi) then
                    delh=dabs(VCell%Cell_Cent(i,j-1,k+1,1)*VCell%nx(i,j-1,k+1)+&
                            VCell%Cell_Cent(i,j-1,k+1,2)*VCell%ny(i,j-1,k+1)+&
                            VCell%Cell_Cent(i,j-1,k+1,3)*VCell%nz(i,j-1,k+1)+&
                              VCell%phi(i,j-1,k+1))+tol
                    vs=v(i,j-1,k+1)
                  else
                    delh=delhec+tol
                    vs=0.d0
                  end if
                  vs=vb+(vs-vb)*delhec/delh
                  vsp=0.5d0*(vs+dabs(vs))
                  vsn=0.5d0*(vs-dabs(vs))
                  Flux(i,j,k,3)=(vsp*w(i,j-1,k)+vsn*w(i,j,k))*delhec/delh*     &
                         WCell%NEArea(i,j-1,k)*WGrid%dx(i,j,k)*WGrid%dz(i,j,k)
                else
                  Sz=VCell%SzT(i,j-1,k)
                  eta=dabs(WCell%FCT(i,j-1,k,3)+0.5d0*WGrid%dz(i,j-1,k)-       &
                           VCell%Cell_Cent(i,j-1,k,3))/Sz
                  if(dabs(eta)>=1.d0) eta=0.5d0
                  vs=(1.d0-eta)*v(i,j-1,k)+eta*v(i,j-1,k+1)
                  vsp=0.5d0*(vs+dabs(vs))
                  vsn=0.5d0*(vs-dabs(vs))
                  Flux(i,j,k,3)=(vsp*w(i,j-1,k)+vsn*w(i,j,k))*           &
                         WCell%NEArea(i,j-1,k)*WGrid%dx(i,j,k)*WGrid%dz(i,j,k)
                end if
              end if
            end if
          end do
        
        
          ul=u(i,Jmax,k);ur=u(i,Jmax+1,k)
          alr=0.5d0*(v(i,Jmax,k)+v(i+1,Jmax,k))
          flux(i,Jmax+1,k,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*            &
                   UCell%NEArea(i,Jmax,k)*UGrid%dx(i,Jmax,k)*UGrid%dz(i,Jmax,k)
          vl=v(i,Jmax,k);vr=v(i,Jmax+1,k)
          alr=(vr+vl)
          alr=dmax1(dabs(vl),dabs(vr))
          flux(i,Jmax+1,k,2)=0.5d0*(vr**2.d0+vl**2.d0-dabs(alr)*(vr-vl))*        &
                    VCell%NEArea(i,Jmax,k)*VGrid%dx(i,Jmax,k)*VGrid%dz(i,Jmax,k)
          wl=w(i,Jmax,k);wr=w(i,Jmax+1,k)
          alr=0.5d0*(v(i,Jmax,k)+v(i,Jmax,k+1))
          flux(i,Jmax+1,k,3)=0.5d0*(alr*wr+alr*wl-dabs(alr)*(wr-wl))*        &
                     WCell%NEArea(i,Jmax,k)*WGrid%dx(i,Jmax,k)*WGrid%dz(i,Jmax,k)        
        end do
      end do
    end subroutine HighOrderConvectiveFluxForYDir
    
    subroutine HighOrderConvectiveFluxForZDir(PGrid,UGrid,VGrid,WGrid,UCell,   &
                                                    VCell,WCell,vb,flux)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN)    :: PGrid,UGrid,VGrid,WGrid
      TYPE(Cell),INTENT(IN)    :: UCell,VCell,WCell
      REAL(KIND=dp),INTENT(IN) :: vb
      REAL(KIND=dp),DIMENSION(:,:,:,:),allocatable,INTENT(INOUT)::flux
      INTEGER(kind=it4b)       :: i,j,k,Lim
      REAL(KIND=dp)            :: ul,ur,vl,vr,wl,wr,alr,delhec,delh,tol
      REAL(KIND=dp)            :: omei,omei1,ri,ri1,tolim,Sx,Sy,Sz,wb,wbp,wbn,eta
      Lim=2
      tol=1.d-24
      do i=1,Imax
        do j=1,Jmax
          ul=u(i,j,0);ur=u(i,j,1)
          alr=0.5d0*(w(i,j,0)+w(i+1,j,0))
          flux(i,j,1,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*               &
                        UCell%TEArea(i,j,1)*UGrid%dx(i,j,1)*UGrid%dy(i,j,1)
          vl=v(i,j,0);vr=v(i,j,1)
          alr=0.5d0*(w(i,j,0)+w(i,j+1,0))
          flux(i,j,1,2)=0.5d0*(vr**2.d0+vl**2.d0-dabs(alr)*(vr-vl))*           &
                        VCell%TEArea(i,j,1)*VGrid%dx(i,j,1)*VGrid%dy(i,j,1)
          wl=w(i,j,0);wr=w(i,j,1)
          alr=dmax1(dabs(wl),dabs(wr))
          flux(i,j,1,3)=0.5d0*(alr*wr+alr*wl-dabs(alr)*(wr-wl))*         &
              WCell%TEArea(i,j,1)*WGrid%dx(i,j,1)*WGrid%dy(i,j,1) 
          do k=2,Kmax
            tolim=dmin1(2.d0*UGrid%dz(i,j,k-1)/WGrid%dz(i,j,k-1),              &
                        2.d0*UGrid%dz(i,j,k)/WGrid%dz(i,j,k-1))
            ri=((u(i,j,k-1)-u(i,j,k-2))/WGrid%dz(i,j,max(1,k-2)))/             &
               ((u(i,j,k)-u(i,j,k-1))/WGrid%dz(i,j,k-1))
            omei=MUSCLLimiter(ri,Lim,tolim)*(u(i,j,k)-u(i,j,k-1))/WGrid%dz(i,j,k-1)
            ri1=((u(i,j,k)-u(i,j,k-1))/WGrid%dz(i,j,k-1))/                     &
                ((u(i,j,k+1)-u(i,j,k))/WGrid%dz(i,j,k))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(u(i,j,k+1)-u(i,j,k))/WGrid%dz(i,j,k)
            ul=u(i,j,k-1)+0.5d0*(PGrid%dz(i,j,k-1))*omei
            ur=u(i,j,k)-0.5d0*(PGrid%dz(i,j,k))*omei1
            alr=0.5d0*(v(i,j,k-1)+v(i+1,j,k-1))
       !   alr=dmax1(dabs(v(i,j-1)),dabs(v(i+1,j-1)))
            flux(i,j,k,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*             &
                          UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k)*UGrid%dy(i,j,k)
            if(k>2.and.k<Kmax) then
              if(UCell%vof(i,j,k-2)<1.d0-epsi.or.            &
                 UCell%vof(i,j,k-1)<1.d0-epsi.or.                      &
                 UCell%Vof(i,j,k)<1.d0-epsi.or.              &
                 UCell%vof(i,j,k+1)<1.d0-epsi) then
                if(UCell%TEArea(i,j,k-1)<0.5d0) then
                  delhec=dabs(UCell%FCT(i,j,k-1,1)*UCell%nx(i,j,k-1)+          &
                              UCell%FCT(i,j,k-1,2)*UCell%ny(i,j,k-1)+        &
                              UCell%FCT(i,j,k-1,3)*UCell%nz(i,j,k-1)+          &
                              UCell%phi(i,j,k-1))
                  if(WCell%MoExCell(i,j,k-1)/=1.and.             &
                     WCell%Vof(i,j,k-1)>epsi)then
                    delh=dabs(WCell%Cell_Cent(i,j,k-1,1)*WCell%nx(i,j,k-1)+    &
                              WCell%Cell_Cent(i,j,k-1,2)*WCell%ny(i,j,k-1)+    &
                              WCell%Cell_Cent(i,j,k-1,3)*WCell%nz(i,j,k-1)+    &  
                              WCell%phi(i,j,k-1))+tol
                    wb=w(i,j,k-1)
                  elseif(WCell%MoExCell(i+1,j,k-1)/=1.and.                     &
                         WCell%Vof(i-1,j,k+1)>epsi) then
                    delh=dabs(WCell%Cell_Cent(i+1,j,k-1,1)*WCell%nx(i+1,j,k-1)+&
                            WCell%Cell_Cent(i+1,j,k-1,2)*WCell%ny(i+1,j,k-1)+&
                            WCell%Cell_Cent(i+1,j,k-1,3)*WCell%nz(i+1,j,k-1)+&
                              WCell%phi(i+1,j,k-1))+tol
                    wb=w(i+1,j,k-1)
                  else
                    delh=delhec+tol
                    wb=0.d0
                  end if
                 ! wb=vb+(wb-vb)*delhec/delh
                  wbp=0.5d0*(wb+dabs(wb))
                  wbn=0.5d0*(wb-dabs(wb))
                  Flux(i,j,k,1)=(wbp*u(i,j,k-1)+wbn*u(i,j,k))*delhec/delh*     &
                           UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k)*UGrid%dy(i,j,k)
                else
                  Sx=WCell%SxE(i,j,k-1)
                  eta=dabs(UCell%FCT(i,j,k-1,1)+0.5d0*UGrid%dx(i,j,k-1)-       &
                           WCell%Cell_Cent(i,j,k-1,1))/Sx
                  if(dabs(eta)>=1.d0) eta=0.5d0
                  wb=(1.d0-eta)*w(i,j,k-1)+eta*w(i+1,j,k-1)
                  wbp=0.5d0*(wb+dabs(wb))
                  wbn=0.5d0*(wb-dabs(wb))
                  Flux(i,j,k,1)=(wbp*u(i,j,k-1)+wbn*u(i,j,k))*           &
                         UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k)*UGrid%dy(i,j,k)
                end if
              end if
            end if
            
            ! For VCell
            
            tolim=dmin1(2.d0*VGrid%dz(i,j,k-1)/WGrid%dz(i,j,k-1),              &
                        2.d0*VGrid%dz(i,j,k)/WGrid%dz(i,j,k-1))
            ri=((v(i,j,k-1)-v(i,j,k-2))/WGrid%dz(i,j,max(1,k-2)))/             &
               ((v(i,j,k)-v(i,j,k-1))/WGrid%dz(i,j,k-1))
            omei=MUSCLLimiter(ri,Lim,tolim)*(v(i,j,k)-v(i,j,k-1))/WGrid%dz(i,j,k-1)
            ri1=((v(i,j,k)-v(i,j,k-1))/WGrid%dz(i,j,k-1))/                     &
                ((v(i,j,k+1)-v(i,j,k))/WGrid%dz(i,j,k))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(v(i,j,k+1)-v(i,j,k))/WGrid%dz(i,j,k)
            vl=v(i,j,k-1)+0.5d0*(PGrid%dz(i,j,k-1))*omei
            vr=v(i,j,k)-0.5d0*(PGrid%dz(i,j,k))*omei1
            alr=0.5d0*(v(i,j,k-1)+v(i,j+1,k-1))
       !   alr=dmax1(dabs(v(i,j-1)),dabs(v(i+1,j-1)))
            flux(i,j,k,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*             &
                          VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k)*VGrid%dy(i,j,k)
            if(k>2.and.k<Kmax) then
              if(VCell%vof(i,j,k-2)<1.d0-epsi.or.            &
                 VCell%vof(i,j,k-1)<1.d0-epsi.or.                      &
                 VCell%Vof(i,j,k)<1.d0-epsi.or.              &
                 VCell%vof(i,j,k+1)<1.d0-epsi) then
                if(VCell%TEArea(i,j,k-1)<0.5d0) then
                  delhec=dabs(VCell%FCT(i,j,k-1,1)*VCell%nx(i,j,k-1)+          &
                              VCell%FCT(i,j,k-1,2)*VCell%ny(i,j,k-1)+        &
                              VCell%FCT(i,j,k-1,3)*VCell%nz(i,j,k-1)+          &
                              VCell%phi(i,j,k-1))
                  if(WCell%MoExCell(i,j,k-1)/=1.and.             &
                     WCell%Vof(i,j,k-1)>epsi)then
                    delh=dabs(WCell%Cell_Cent(i,j,k-1,1)*WCell%nx(i,j,k-1)+    &
                              WCell%Cell_Cent(i,j,k-1,2)*WCell%ny(i,j,k-1)+    &
                              WCell%Cell_Cent(i,j,k-1,3)*WCell%nz(i,j,k-1)+    &  
                              WCell%phi(i,j,k-1))+tol
                    wb=w(i,j,k-1)
                  elseif(WCell%MoExCell(i,j+1,k-1)/=1.and.                     &
                         WCell%Vof(i,j+1,k-1)>epsi) then
                    delh=dabs(WCell%Cell_Cent(i,j+1,k-1,1)*VCell%nx(i,j+1,k-1)+&
                            WCell%Cell_Cent(i,j+1,k-1,2)*VCell%ny(i,j+1,k-1)+&
                            WCell%Cell_Cent(i,j+1,k-1,3)*VCell%nz(i,j+1,k-1)+&
                              WCell%phi(i,j+1,k-1))+tol
                    wb=w(i,j+1,k-1)
                  else
                    delh=delhec+tol
                    wb=0.d0
                  end if
                  wbp=0.5d0*(wb+dabs(wb))
                  wbn=0.5d0*(wb-dabs(wb))
                  Flux(i,j,k,2)=(wbp*v(i,j,k-1)+wbn*v(i,j,k))*delhec/delh*     &
                           VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k)*VGrid%dy(i,j,k)
                else
                  Sy=WCell%SyN(i,j,k-1)
                  eta=dabs(VCell%FCT(i,j,k-1,2)+0.5d0*VGrid%dy(i,j,k-1)-       &
                           WCell%Cell_Cent(i,j,k-1,2))/Sy
                  if(dabs(eta)>=1.d0) eta=0.5d0
                  wb=(1.d0-eta)*w(i,j,k-1)+eta*w(i,j+1,k-1)
                  wbp=0.5d0*(wb+dabs(wb))
                  wbn=0.5d0*(wb-dabs(wb))
                  Flux(i,j,k,2)=(wbp*v(i,j,k-1)+wbn*v(i,j,k))*           &
                         VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k)*VGrid%dy(i,j,k)
                end if
              end if
            end if
            
            ! For WCell 

            tolim=dmin1(2.d0*WGrid%dz(i,j,k-1)/PGrid%dz(i,j,k),                &
                        2.d0*WGrid%dz(i,j,k)/PGrid%dz(i,j,k))
            ri=((w(i,j,k-1)-w(i,j,k-2))/PGrid%dz(i,j,k-1))/                    &
               ((w(i,j,k)-w(i,j,k-1))/PGrid%dz(i,j,k))
            omei=MUSCLlimiter(ri,Lim,tolim)*(w(i,j,k)-w(i,j,k-1))/PGrid%dz(i,j,k-1)
            ri1=((v(i,j,k)-v(i,j,k-1))/PGrid%dz(i,j,k))/                       &
                ((v(i,j,k+1)-v(i,j,k))/PGrid%dz(i,j,min(Kmax,k+1)))
            omei1=MUSCLLimiter(ri1,Lim,tolim)*(w(i,j,k+1)-w(i,j,k))/           &
                                              PGrid%dz(i,j,min(Kmax,k+1))
            wl=w(i,j,k-1)+0.5d0*(PGrid%dz(i,j,k))*omei
            wr=w(i,j,k)-0.5d0*(PGrid%dz(i,j,k))*omei1
            alr=(wr+wl)
            alr=dmax1(dabs(wl),dabs(wr))
            flux(i,j,k,3)=0.5d0*(wr*wr+wl*wl-dabs(alr)*(wr-wl))*               &
                           WCell%TEArea(i,j,k-1)*WGrid%dx(i,j,k)*WGrid%dy(i,j,k)
            if(k>2.and.k<Kmax) then
              if(WCell%vof(i,j,k-2)<1.d0-epsi.or.            &
                 WCell%vof(i,j,k-1)<1.d0-epsi.or.                      &
                 WCell%Vof(i,j,k)<1.d0-epsi.or.              &
                 WCell%vof(i,j,k+1)<1.d0-epsi) then
                eta=WCell%EtaT(i,j,k-1)
                wb=(1.d0-eta)*w(i,j,k-1)+eta*w(i,j,k)
                wbp=0.5d0*(wb+dabs(wb))
                wbn=0.5d0*(wb-dabs(wb))
                Flux(i,j,k,3)=(wbp*w(i,j,k-1)+wbn*w(i,j,k))*                   &
                           WCell%TEArea(i,j,k-1)*WGrid%dx(i,j,k)*WGrid%dy(i,j,k)
              end if
            end if
          end do
        
          ul=u(i,j,kmax);ur=u(i,j,kmax+1)
          alr=0.5d0*(w(i,j,kmax)+w(i+1,j,kmax))
          flux(i,j,kmax+1,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*          &
                   UCell%TEArea(i,j,kmax)*UGrid%dx(i,j,kmax)*UGrid%dy(i,j,kmax)
          vl=v(i,j,kmax);vr=v(i,j,kmax+1)
          alr=0.5d0*(w(i,j,kmax)+w(i,j+1,kmax))
          flux(i,j,kmax+1,2)=0.5d0*(alr*ur+alr*vl-dabs(alr)*(vr-vl))*          &
                   VCell%TEArea(i,j,kmax)*VGrid%dx(i,j,kmax)*VGrid%dy(i,j,kmax)
          wl=w(i,Jmax,k);wr=w(i,Jmax+1,k)
          alr=dmax1(dabs(wl),dabs(wr))
          flux(i,j,kmax+1,3)=0.5d0*(wr**2.d0+wl**2.d0-dabs(alr)*(wr-wl))*      &
                   WCell%TEArea(i,j,kmax)*WGrid%dx(i,j,kmax)*WGrid%dy(i,j,kmax)        
        end do
      end do
    end subroutine HighOrderConvectiveFluxForZDir
    
    Subroutine DiffusiveFlux(PGrid,UGrid,VGrid,WGrid,PCell,UCell,VCell,WCell,  &
                                                     flux,EFlux,nuref,iu,iv,iw)
      !! The subroutine is used to compute the coefficient for diffusive flux.
      Implicit none
      Integer(kind=it4b),intent(in)        :: iu,iv,iw
      Type(Cell),intent(in)          :: PCell,UCell,VCell,WCell
      Type(Grid),intent(in)          :: PGrid,UGrid,VGrid,WGrid
      Real(kind=dp),dimension(:,:,:,:),allocatable,intent(inout) :: flux,EFlux
      Real(kind=dp),intent(in)           :: nuref
      Integer(kind=it4b)           :: i,j,k
      Real(kind=dp)            :: Sx,Sy,Sz,tol
      real(kind=dp)            :: NuF,VflF
      tol = 1.d-30
      Do i = 1,Imax+iu
        Do j = 1,Jmax+iv
          Do k = 1,Kmax+iw
            If(iu==1) then
              If(i==Imax+iu) then
                ! Compute the mixing viscosity for U-Cell
                NuF=UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,1) = UCell%EEArea(i-1,j,k)*UGrid%dy(i-1,j,k)*NuF*   &
                                UGrid%dz(i-1,j,k)/PGrid%dx(i-1,j,k)/Rey
                ! Compute the mixing viscosity for V-Cell
                NuF=VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,2) = VCell%EEArea(i-1,j,k)*VGrid%dy(i-1,j,k)*NuF*   &
                                VGrid%dz(i-1,j,k)/UGrid%dx(i-1,j,k)/Rey
                ! Compute the mixing viscosity for W-Celll
                NuF=WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,3) = WCell%EEArea(i-1,j,k)*WGrid%dy(i-1,j,k)*NuF*   &
                                WGrid%dz(i-1,j,k)/UGrid%dx(i-1,j,k)/Rey
              Elseif(i==1) then
                ! Compute the mixing viscosity for U-Cell
                NuF=UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,1)=UCell%EEArea(i,j,k)*UGrid%dy(i,j,k)*NuF*         &
                              UGrid%dz(i,j,k)/PGrid%dx(i,j,k)/Rey
                ! Compute the mixing viscosity for V-Cell
                NuF=VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,2)=VCell%EEArea(i,j,k)*VGrid%dy(i,j,k)*nuF*         &
                                   VGrid%dz(i,j,k)/UGrid%dx(i,j,k)/Rey
                ! Compute the mixing viscosity for W-Celll
                NuF=WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,3)=WCell%EEArea(i,j,k)*WGrid%dy(i,j,k)*nuF*         &
                                   WGrid%dz(i,j,k)/UGrid%dx(i,j,k)/Rey
              Else
                Sx=UCell%SxE(i-1,j,k)
                Sy=UCell%Cell_Cent(i,j,k,2)-UCell%Cell_Cent(i-1,j,k,2)
                Sz=UCell%Cell_Cent(i,j,k,3)-UCell%Cell_Cent(i-1,j,k,3)
                if(UCell%vof(i,j,k)>epsi.and.UCell%vof(i+1,j,k)>epsi) then
                  VflF=(1.d0-UCell%EtaE(i,j,k))*             &
                     UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)+           &
                     UCell%EtaE(i,j,k)*UCell%vofL(i+1,j,k)/(UCell%vof(i+1,j,k)+tol)
                else
                  VflF=dmax1(UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol),         &
                             UCell%vofL(i+1,j,k)/(UCell%vof(i+1,j,k)+tol))
                end if
                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
              !  if(dabs(NuF-1.d0)>1.d-5) then
              !    print*, vflF
              !    print*, UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)
              !    print*, UCell%vofL(i+1,j,k),UCell%vof(i+1,j,k)
              !    print*, UCell%vofL(i+1,j,k)/(UCell%vof(i+1,j,k)+tol)
              !    pause 'PredicUVW 1052'
              !  end if
                flux(i,j,k,1)=UCell%EEArea(i-1,j,k)*UGrid%dy(i,j,k)*NuF*       &
                                                      UGrid%dz(i,j,k)/Sx/Rey
                If(dabs(Sy)>1.d-4*UGrid%dy(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         UGrid%dz(i,j,k)) then
                  EFlux(i,j,k,1)=UCell%DAlE(i-1,j,k)*UCell%EEArea(i-1,j,k)*    &
                                 NuF*UGrid%dy(i-1,j,k)*UGrid%dz(i-1,j,k)/Rey
                End if
                Sx=VCell%SxE(i-1,j,k)
                Sy=VCell%Cell_Cent(i,j,k,2)-VCell%Cell_Cent(i-1,j,k,2)
                Sz=VCell%Cell_Cent(i,j,k,3)-VCell%Cell_Cent(i-1,j,k,3)

                if(VCell%vof(i,j,k)>epsi.and.VCell%vof(i+1,j,k)>epsi) then
                  VflF=(1.d0-VCell%EtaE(i,j,k))*             &
                     VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)+           &
                     VCell%EtaE(i,j,k)*VCell%vofL(i+1,j,k)/(VCell%vof(i+1,j,k)+tol)
                else
                  VflF=dmax1(VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol),         &
                             VCell%vofL(i+1,j,k)/(VCell%vof(i+1,j,k)+tol))
                end if

                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
              !  if(dabs(NuF-1.d0)>1.d-5) then
              !    pause 'PredicUVW 1069'
              !  end if
                flux(i,j,k,2)=VCell%EEArea(i-1,j,k)*VGrid%dy(i,j,k)*NuF*       &
                                                      VGrid%dz(i,j,k)/Sx/Rey
                If(dabs(Sy)>1.d-4*VGrid%dy(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         VGrid%dz(i,j,k)) then
                  EFlux(i,j,k,2)=VCell%DAlE(i-1,j,k)*VCell%EEArea(i-1,j,k)*    &
                                 NuF*VGrid%dy(i-1,j,k)*VGrid%dz(i-1,j,k)/Rey
                End if
                Sx=WCell%SxE(i-1,j,k)
                Sy=WCell%Cell_Cent(i,j,k,2)-WCell%Cell_Cent(i-1,j,k,2)
                Sz=WCell%Cell_Cent(i,j,k,3)-WCell%Cell_Cent(i-1,j,k,3)

                if(WCell%vof(i,j,k)>epsi.and.WCell%vof(i+1,j,k)>epsi) then
                  VflF=(1.d0-WCell%EtaE(i,j,k))*             &
                     WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)+           &
                     WCell%EtaE(i,j,k)*WCell%vofL(i+1,j,k)/(WCell%vof(i+1,j,k)+tol)
                else
                  VflF=dmax1(WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol),         &
                             WCell%vofL(i+1,j,k)/(WCell%vof(i+1,j,k)+tol))
                end if
                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
             !   if(dabs(NuF-1.d0)>1.d-5) then
             !     pause 'PredicUVW 1088'
             !   end if
                flux(i,j,k,3)=WCell%EEArea(i-1,j,k)*WGrid%dy(i,j,k)*NuF*       &
                                                      WGrid%dz(i,j,k)/Sx/Rey
                If(dabs(Sy)>1.d-4*WGrid%dy(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         WGrid%dz(i,j,k)) then
                  EFlux(i,j,k,3)=WCell%DAlE(i-1,j,k)*WCell%EEArea(i-1,j,k)*    &
                                   NuF*WGrid%dy(i-1,j,k)*WGrid%dz(i-1,j,k)/Rey
                End if
              End if
            End if
            If(iv==1) then
              If(j==Jmax+iv) then
                ! Compute the mixing viscosity for U-Cell
                NuF=UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,1) = UCell%NEArea(i,j-1,k)*UGrid%dx(i,j-1,k)*NuF*   &
                                     UGrid%dz(i,j-1,k)/VGrid%dy(i,j-1,k)/Rey
                ! Compute the mixing viscosity for V-Cell
                NuF=VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,2) = VCell%NEArea(i,j-1,k)*VGrid%dx(i,j-1,k)*NuF*   &
                                     VGrid%dz(i,j-1,k)/PGrid%dy(i,j-1,k)/Rey
                ! Compute the mixing viscosity for W-Cell
                NuF=WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,3) = WCell%NEArea(i,j-1,k)*WGrid%dx(i,j-1,k)*NuF*   &
                                     WGrid%dz(i,j-1,k)/VGrid%dy(i,j-1,k)/Rey
              Elseif(j==1) then
                ! Compute the mixing viscosity for U-Cell
                NuF=UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,1) = UCell%NEArea(i,j,k)*UGrid%dx(i,j,k)*NuF*        &
                                     UGrid%dz(i,j,k)/VGrid%dy(i,j,k)/Rey
                ! Compute the mixing viscosity for V-Cell
                NuF=VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,2) = VCell%NEArea(i,j,k)*VGrid%dx(i,j,k)*NuF*       &
                                     VGrid%dz(i,j,k)/PGrid%dy(i,j,k)/Rey
                ! Compute the mixing viscosity for W-Cell
                NuF=WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,3) = WCell%NEArea(i,j,k)*WGrid%dx(i,j,k)*NuF*       &
                                     WGrid%dz(i,j,k)/VGrid%dy(i,j,k)/Rey
              Else
                Sx = UCell%Cell_Cent(i,j,k,1)-UCell%Cell_Cent(i,j-1,k,1)
                Sy = UCell%SyN(i,j-1,k)
                Sz = UCell%Cell_Cent(i,j,k,3)-UCell%Cell_Cent(i,j-1,k,3)
                if(UCell%vof(i,j,k)>epsi.and.UCell%vof(i,j+1,k)>epsi) then
                  VflF=(1.d0-UCell%EtaN(i,j,k))*             &
                     UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)+           &
                     UCell%EtaN(i,j,k)*UCell%vofL(i,j+1,k)/(UCell%vof(i,j+1,k)+tol)
                else
                  VflF=dmax1(UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol),         &
                             UCell%vofL(i,j+1,k)/(UCell%vof(i,j+1,k)+tol))
                end if
                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
              !  if(dabs(NuF-1.d0)>1.d-5) then
              !    pause 'PredicUVW 1141'
              !  end if
                flux(i,j,k,1)=UCell%NEArea(i,j-1,k)*UGrid%dx(i,j,k)*NuF*       &
                                                      UGrid%dz(i,j,k)/Sy/Rey
                If(dabs(Sx)>1.d-4*UGrid%dx(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         UGrid%dz(i,j,k)) then
                  EFlux(i,j,k,1)=UCell%DAlN(i,j-1,k)*UCell%NEArea(i,j-1,k)*NuF*&
                                 UGrid%dx(i,j-1,k)*UGrid%dz(i,j-1,k)/Rey
                End if
                Sx = VCell%Cell_Cent(i,j,k,1)-VCell%Cell_Cent(i,j-1,k,1)
                Sy = VCell%SyN(i,j-1,k)
                Sz = VCell%Cell_Cent(i,j,k,3)-VCell%Cell_Cent(i,j-1,k,3)
                if(VCell%vof(i,j,k)>epsi.and.VCell%vof(i,j+1,k)>epsi) then
                  VflF=(1.d0-VCell%EtaN(i,j,k))*             &
                     VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)+           &
                     VCell%EtaN(i,j,k)*VCell%vofL(i,j+1,k)/(VCell%vof(i,j+1,k)+tol)
                else
                  VflF=dmax1(VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol),         &
                             VCell%vofL(i,j+1,k)/(VCell%vof(i,j+1,k)+tol))
                end if
                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
              !  if(dabs(NuF-1.d0)>1.d-5) then
              !    pause 'PredicUVW 1158'
              !  end if
                flux(i,j,k,2)=VCell%NEArea(i,j-1,k)*VGrid%dx(i,j,k)*NuF*       &
                                                      VGrid%dz(i,j,k)/Sy/Rey
                If(dabs(Sx)>1.d-4*VGrid%dx(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         VGrid%dz(i,j,k)) then
                  EFlux(i,j,k,2)=VCell%DAlN(i,j-1,k)*VCell%NEArea(i,j-1,k)*NuF*&
                                   VGrid%dx(i,j-1,k)*VGrid%dz(i,j-1,k)/Rey
                End if
                Sx=WCell%Cell_Cent(i,j,k,1)-WCell%Cell_Cent(i,j-1,k,1)
                Sy=WCell%SyN(i,j-1,k)
                Sz=WCell%Cell_Cent(i,j,k,3)-WCell%Cell_Cent(i,j-1,k,3)
                if(WCell%vof(i,j,k)>epsi.and.WCell%vof(i,j+1,k)>epsi) then
                  VflF=(1.d0-WCell%EtaN(i,j,k))*             &
                     WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)+           &
                     WCell%EtaN(i,j,k)*WCell%vofL(i,j+1,k)/(WCell%vof(i,j+1,k)+tol)
                else
                  VflF=dmax1(WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol),         &
                             WCell%vofL(i,j+1,k)/(WCell%vof(i,j+1,k)+tol))
                end if
                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
             !   if(dabs(NuF-1.d0)>1.d-5) then
             !     pause 'PredicUVW 1175'
             !   end if
                flux(i,j,k,3)=WCell%NEArea(i,j-1,k)*WGrid%dx(i,j,k)*NuF*       &
                                                      WGrid%dz(i,j,k)/Sy/Rey
                If(dabs(Sx)>1.d-4*WGrid%dx(i,j,k).or.dabs(Sz)>1.d-4*           &
                                                         WGrid%dz(i,j,k)) then
                  EFlux(i,j,k,3)=WCell%DAlN(i,j-1,k)*WCell%NEArea(i,j-1,k)*NuF*&
                                   WGrid%dx(i,j-1,k)*WGrid%dz(i,j-1,k)/Rey
                End if
              End if
            End if
            If(iw==1) then
              If(k==Kmax+iw) then
                ! Compute the mixing viscosity for U-Cell
                NuF=UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,1) = UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k-1)*NuF*   &
                                     UGrid%dy(i,j,k-1)/WGrid%dz(i,j,k-1)/Rey
                ! Compute the mixing viscosity for V-Cell
                NuF=VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,2) = VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k-1)*NuF*   &
                                     VGrid%dy(i,j,k-1)/WGrid%dz(i,j,k-1)/Rey
                ! Compute the mixing viscosity for W-Cell
                NuF=WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,3) = WCell%TEArea(i,j,k-1)*WGrid%dx(i,j,k-1)*NuF*   &
                                     WGrid%dy(i,j,k-1)/PGrid%dz(i,j,k-1)/Rey
              Elseif(k==1) then
                ! Compute the mixing viscosity for U-Cell
                NuF=UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,1) = UCell%TEArea(i,j,k)*UGrid%dx(i,j,k)*NuF*       &
                                     UGrid%dy(i,j,k)/WGrid%dz(i,j,k)/Rey
                ! Compute the mixing viscosity for V-Cell
                NuF=VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,2) = VCell%TEArea(i,j,k)*VGrid%dx(i,j,k)*NuF*       &
                                     VGrid%dy(i,j,k)/WGrid%dz(i,j,k)/Rey
                ! Compute the mixing viscosity for W-Cell
                NuF=WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)*nuw/nuref+        &
                    (1.d0-WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol))*nua/nuref
                flux(i,j,k,3) = WCell%TEArea(i,j,k)*WGrid%dx(i,j,k)*NuF*       &
                                     WGrid%dy(i,j,k)/PGrid%dz(i,j,k)/Rey
              Else
                Sx = UCell%Cell_Cent(i,j,k,1)-UCell%Cell_Cent(i,j,k-1,1)
                Sy = UCell%Cell_Cent(i,j,k,2)-UCell%Cell_Cent(i,j,k-1,2)
                Sz = UCell%SzT(i,j,k-1)
                if(UCell%vof(i,j,k)>epsi.and.UCell%vof(i,j,k+1)>epsi) then
                  VflF=(1.d0-UCell%EtaT(i,j,k))*             &
                     UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol)+           &
                     UCell%EtaT(i,j,k)*UCell%vofL(i,j,k+1)/(UCell%vof(i,j,k+1)+tol)
                else
                  VflF=dmax1(UCell%vofL(i,j,k)/(UCell%vof(i,j,k)+tol),         &
                             UCell%vofL(i,j,k+1)/(UCell%vof(i,j,k+1)+tol))
                end if
                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
             !   if(dabs(NuF-1.d0)>1.d-5) then
             !     pause 'PredicUVW 1228'
             !   end if
                flux(i,j,k,1) = UCell%TEArea(i,j,k-1)*UGrid%dx(i,j,k)*NuF*     &
                                                      UGrid%dy(i,j,k)/Sz/Rey
                If(dabs(Sx)>1.d-4*UGrid%dx(i,j,k).or.dabs(Sy)>1.d-4*           &
                                                         UGrid%dy(i,j,k)) then
                  EFlux(i,j,k,1)=UCell%DAlT(i,j,k-1)*UCell%TEArea(i,j,k-1)*NuF*&
                                   UGrid%dx(i,j,k-1)*UGrid%dy(i,j,k-1)/Rey
                End if

                Sx = VCell%Cell_Cent(i,j,k,1)-VCell%Cell_Cent(i,j,k-1,1)
                Sy = VCell%Cell_Cent(i,j,k,2)-VCell%Cell_Cent(i,j,k-1,2)
                Sz = VCell%SzT(i,j,k-1)
                if(VCell%vof(i,j,k)>epsi.and.VCell%vof(i,j,k+1)>epsi) then
                  VflF=(1.d0-VCell%EtaT(i,j,k))*             &
                     VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol)+           &
                     VCell%EtaT(i,j,k)*VCell%vofL(i,j,k+1)/(VCell%vof(i,j,k+1)+tol)
                else
                  VflF=dmax1(VCell%vofL(i,j,k)/(VCell%vof(i,j,k)+tol),         &
                             VCell%vofL(i,j,k+1)/(VCell%vof(i,j,k+1)+tol))
                end if
                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
             !   if(dabs(NuF-1.d0)>1.d-5) then
             !     pause 'PredicUVW 1246'
             !   end if
                flux(i,j,k,2) = VCell%TEArea(i,j,k-1)*VGrid%dx(i,j,k)*NuF*     &
                                                      VGrid%dy(i,j,k)/Sz/Rey
                If(dabs(Sx)>1.d-4*VGrid%dx(i,j,k).or.dabs(Sy)>1.d-4*           &
                                                         VGrid%dy(i,j,k)) then
                  EFlux(i,j,k,2)=VCell%DAlT(i,j,k-1)*VCell%TEArea(i,j,k-1)*NuF*&
                                 VGrid%dx(i,j,k-1)*VGrid%dy(i,j,k-1)/Rey
                End if
                Sx = WCell%Cell_Cent(i,j,k,1)-WCell%Cell_Cent(i,j,k-1,1)
                Sy = WCell%Cell_Cent(i,j,k,2)-WCell%Cell_Cent(i,j,k-1,2)
                Sz = WCell%SzT(i,j,k-1)
                if(WCell%vof(i,j,k)>epsi.and.WCell%vof(i,j,k+1)>epsi) then
                  VflF=(1.d0-WCell%EtaT(i,j,k))*             &
                     WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol)+           &
                     WCell%EtaT(i,j,k)*WCell%vofL(i,j,k+1)/(WCell%vof(i,j,k+1)+tol)
                else
                  VflF=dmax1(WCell%vofL(i,j,k)/(WCell%vof(i,j,k)+tol),         &
                             WCell%vofL(i,j,k+1)/(WCell%vof(i,j,k+1)+tol))
                end if
                NuF=VflF*nuw/nuref+(1.d0-VflF)*nua/nuref
              !  if(dabs(NuF-1.d0)>1.d-5) then
              !    pause 'PredicUVW 1263'
              !  end if
                flux(i,j,k,3)=WCell%TEArea(i,j,k-1)*WGrid%dx(i,j,k)*NuF*       &
                                                      WGrid%dy(i,j,k)/Sz/Rey

                If(dabs(Sx)>1.d-4*WGrid%dx(i,j,k).or.dabs(Sy)>1.d-4*           &
                                                         WGrid%dy(i,j,k)) then
                  EFlux(i,j,k,3)=WCell%DAlT(i,j,k-1)*WCell%TEArea(i,j,k-1)*NuF*&
                                   WGrid%dx(i,j,k-1)*WGrid%dy(i,j,k-1)/Rey
                End if
              End if
            End if
          End do
        End do
      End do
    End subroutine DiffusiveFlux

    Subroutine PredictorVelocityBoundaryCondition(Pred,BCu,BCv,BCw)
      Implicit none
      Type(Predictor),intent(inout):: Pred
      Type(BCBase),intent(in):: BCu,BCv,BCw
      Integer(kind=it4b):: i,j,k
      ! The boundary condition for u velocity
      ! At the Western boundary
      if(BCu%flag(1)==1) then
        do j = 1,Jmax
          do k = 1,Kmax
            Pred%u(1-ight,j,k)=Pred%u(1,j,k)
          end do
        end do
      else
        do j = 1,Jmax
          do k = 1,Kmax
            Pred%u(1-ight,j,k)=BCu%VarW(j,k)
          end do
        end do
      end if
      ! At the Eastern boundary
      if(BCu%flag(2)==1) then
        do j = 1,Jmax
          do k = 1,Kmax
            Pred%u(Imax,j,k)=Pred%u(Imax-1,j,k)
            Pred%u(Imax+1,j,k)=2.d0*BCu%VarE(j,k)-Pred%u(Imax-1,j,k)
          end do
        end do
      else
        do j = 1,Jmax
          do k = 1,Kmax
            Pred%u(Imax,j,k)=BCu%VarE(j,k)
            Pred%u(Imax+1,j,k)=2.d0*BCu%VarE(j,k)-Pred%u(Imax-1,j,k)
          end do
        end do  
      end if
      ! The boundary condition for v velocity
      ! At the Souhthern boundary
      if(BCv%flag(3)==1) then ! Neumann boundary condition
        do i = 1,Imax
          do k = 1,Kmax
            Pred%v(i,1-jght,k)=Pred%v(i,1,k)
          end do
        end do
      else ! Dirichlet boundary condition
        do i = 1,Imax
          do k = 1,Kmax
            Pred%v(i,1-jght,k)=BCv%VarW(i,k)
          end do
        end do 
      end if
      ! At the Western boundary
      if(BCv%flag(4)==1) then
        do i=1,Imax
          do k=1,Kmax
            Pred%v(i,Jmax,k)=Pred%v(i,Jmax-1,k)
            Pred%v(i,Jmax+1,k)=2.d0*BCv%VarE(i,j)-Pred%v(i,Jmax-1,k)
          end do
        end do
      else
        do i=1,Imax
          do k=1,Kmax
            Pred%v(i,Jmax,k)=BCv%VarE(i,k)
            Pred%v(i,Jmax+1,k)=2.d0*BCv%VarE(i,j)-Pred%v(i,Jmax-1,k)
          end do
        end do
      end if
      ! The boundary condition for w velocity
      ! At the bottom boundary
      if(BCw%flag(5)==1) then ! Neumann boundary condition
        do i=1,Imax
          do j=1,Jmax
            Pred%w(i,j,1-kght)=Pred%w(i,j,1)
          end do
        end do
      else
        do i=1,Imax
          do j=1,Jmax
            Pred%w(i,j,1-kght)=BCw%VarB(i,j)
          end do
        end do  
      end if
      ! At the top boundary 
      if(BCw%flag(6)==1) then
        do i=1,Imax
          do j=1,Jmax
            Pred%w(i,j,Kmax)=Pred%w(i,j,Kmax-1)
            Pred%w(i,j,Kmax+1)=2.d0*BCw%VarT(i,j)-Pred%w(i,j,Kmax-1)
          end do
        end do
      else
        do i=1,Imax
          do j=1,Jmax
            Pred%w(i,j,Kmax)=BCw%VarT(i,j)
            Pred%w(i,j,Kmax+1)=2.d0*BCw%VarT(i,j)-Pred%w(i,j,Kmax-1)
          end do
        end do  
      end if
    end subroutine PredictorVelocityBoundaryCondition

    subroutine PredictorVelocityInternalCellCondition(Pred,UCell,VCell,WCell)
      implicit none
      type(Predictor),intent(inout):: Pred
      type(Cell),intent(in):: UCell,VCell,WCell
      integer(kind=it4b):: i,j,k
      do i = 1,Imax
        do j = 1,Jmax
          do k = 1,Kmax
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

    real(kind=dp) function DensityMixture(Lvof,Fvof,Roref) result (mixrho)
      real(kind=dp),intent(in)      :: Lvof,Fvof,Roref

      mixrho = (Fvof-Lvof)*roa/Roref+Lvof*row/Roref
    end function DensityMixture

    real(kind=dp) function ViscosityMixture(Lvof,Fvof,Roref) result (mixnuy)
      real(kind=dp),intent(in)      :: Lvof,Fvof,Roref

      mixnuy = (Fvof-Lvof)*nua/nuref+Lvof*nuw/nuref
    end function ViscosityMixture

    real(kind=dp) function MUSCLLimiter(x,opt,tolim) result (y)
      REAL(KIND=dp),INTENT(IN)      :: x
      INTEGER(kind=it4b),INTENT(IN) :: opt
      REAL(KIND=dp),INTENT(IN)      :: tolim
      select case(opt)
      ! MinMod limiter
        case (1)
          y=dmax1(0.d0,dmin1(tolim,x))
      ! SuperBee limiter
        case (2)
          y=dmax1(0.d0,dmin1(tolim*x,1.d0),dmin1(x,tolim))
      end select
    end function MUSCLLimiter
End module PredictorUV

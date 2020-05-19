!*********************************************************************
!* For cut-cell method
!* To avoid the singularity we also solve the internal cell and set up
!* the value of internal cell afterwards.
Module ProjectionP
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE PredictorUV
    USE MPI
    Implicit none
    Private
    Real(kind=dp),dimension(:,:,:),pointer :: p,u,v,w
    Real(kind=dp),parameter :: amp = 1.d0, tol = 1.d-30
    Real(kind=dp),parameter :: alp = 0.4d0, bet = 0.5d0
    Type,public :: Projection
      Real(dp),dimension(:,:,:),allocatable :: Pp
    End type
    Public :: PoissonEquationSolver
    Interface PoissonEquationSolver
      Module procedure PoissonEquationSolver
    End interface PoissonEquationSolver
    Contains
    subroutine PoissonEquationSolver(PGrid,UGrid,VGrid,WGrid,PCell,UCell,      &
                                     VCell,WCell,TVar,TPred,PU,PV,PW,Proj,dt)
        implicit none
        type(Grid),intent(in)		      :: PGrid,UGrid,VGrid,WGrid
        type(Cell),intent(in)		      :: PCell,UCell,VCell,WCell
        type(Variables),intent(in),target     :: TVar
        type(Predictor),intent(in),target     :: TPred
        type(PoissonCoefficient),intent(in)   :: PU,PV,PW
        type(Projection),intent(inout),target :: Proj
      
        real(kind=dp),intent(in)	      :: dt
        integer*8			      :: A,parcsr_A,b,par_b,x,par_x
        integer*8			      :: solver,precond
        integer(kind=it4b)		      :: num_iterations
        real(kind=dp)			      :: final_res_norm,tol
        real(kind=dp),dimension(:,:,:,:),allocatable :: PoCoef,HJump   ! the coefficient for Poisson solving
        
        allocate(PoCoef(ISize,JSize,KSize,6)) ! The order of the face : W-1,S-2,B-3,E-4,N-5,T-6
        allocate(HJump(ISize,JSize,KSize,6))  
        call Compute1DGFMCoefficient(PGrid,PCell,UGrid,UCell,PU,1,0,0,0,0,0,PoCoef(:,:,:,1))
        call Compute1DGFMCoefficient(PGrid,PCell,VGrid,VCell,PV,0,1,0,0,0,0,PoCoef(:,:,:,2))
        call Compute1DGFMCoefficient(PGrid,PCell,WGrid,WCell,PW,0,0,1,0,0,0,PoCoef(:,:,:,3))
        call Compute1DGFMCoefficient(PGrid,PCell,UGrid,UCell,PU,0,0,0,1,0,0,PoCoef(:,:,:,4))
        call Compute1DGFMCoefficient(PGrid,PCell,VGrid,VCell,PV,0,0,0,0,1,0,PoCoef(:,:,:,5))
        call Compute1DGFMCoefficient(PGrid,PCell,WGrid,WCell,PW,0,0,0,0,0,1,PoCoef(:,:,:,6))
        
        p => TVar%p
        u => TPred%u
        v => TPred%v
        w => TPred%w
        Call SetBasicSolver(solver,precond)
    !   call SetBasicSolver(solver=solver,ierr=ierr)
        Call SetPoissonMatrix(A,parcsr_A,PGrid,PCell,PoCoef,1,0,1,1,1,1)
        Call SetPoissonVectors(b,x,par_b,par_x,PGrid,PCell,dt)
        Call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
        Call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
!        Run info - needed logging turned on
        Call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
        call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
   !     call HYPRE_ParCSRPCGGetResidual(solver,tol,ierr)
        Call DeltaPressureGetValues(x,PCell,Proj)
        Call DeltaPressureBoundaryCondition(Proj,1,0,1,1,1,1)
        Call HYPRE_IJMatrixDestroy(A,ierr)
        Call HYPRE_IJVectorDestroy(b,ierr)
        Call HYPRE_IJVectorDestroy(x,ierr)
        call HYPRE_BoomerAMGDestroy(precond,ierr)
        call HYPRE_ParCSRPCGDestroy(solver,ierr)
        deallocate(PoCoef,HJump)
        Nullify(p)
        Nullify(u)
        Nullify(v)
        Nullify(w)
    end subroutine PoissonEquationSolver
    
    subroutine Compute1DGFMCoefficient(PGrid,PCell,TGrid,TCell,PVel,ium,jvm,  &
                                                    kwm,iup,jvp,kwp,TPoCoef)
      type(Grid),intent(in)		  :: PGrid,TGrid
      type(Cell),intent(in)		  :: PCell,TCell
      type(PoissonCoefficient),intent(in) :: PVel
      integer(kind=it4b),intent(in)	  :: ium,jvm,kwm,iup,jvp,kwp
      real(kind=dp)			  :: Lamda,BetaP,BetaM,BetaW
      integer(kind=it4b)		  :: ii,jj,kk,im,jm,km
      real(kind=dp)			  :: SArea,DGrid,BetaP,BetaM,BetaW
      
      BetaP = 1.d0/(row/Roref)
      BetaM = 1.d0/(roa/Roref)
      do i=1+ium,Isize-iup
        do j=1+jvm,JSize-jvp
          do k=1+kwm,KSize-kwp
            ii = i-ium+iup
            jj = j-jvm+jvp
            kk = k-kwm+kwp
            im = i-ium
            jm = j-jvm
            km = k-kwm
            SArea = dble(ium+iup)*PGrid%dy(i,j,k)*PGrid%dz(i,j,k)+           &
                    dble(jvm+jvp)*PGrid%dx(i,j,k)*PGrid%dz(i,j,k)+           &
                    dble(kwm+kwp)*PGrid%dx(i,j,k)*PGrid%dy(i,j,k)
            DGrid = dble(ium+iup)*TGrid%dx(im,jm,km)+                        &
                    dble(jvm+jvp)*TGrid%dy(im,jm,km)+		             &
                    dble(kwm+kwp)*TGrid%dz(im,jm,km)
                    
            Lamda=dabs(PCell%phi(i,j,k))/(dabs(PCell%phi(i,j,k))+            &
                        dabs(PCell%phi(ii,jj,kk))+tol)                 
            if((PCell%vof(i,j,k)>=0.5d0.and.PCell%vofS(i,j,k)<epsi).or.      &
               (PCell%phi(i,j,k)<0.d0.and.PCell%vofS(i,j,k)>=epsi)) then     ! The cell is in liquid and it is assigned wet cell 
            
              if((PCell%vof(ii,jj,kk)<0.5d0.and.        		     &
                  PCell%vofS(ii,jj,kk)<epsi).or. 	     		     &
                 (PCell%phi(ii,jj,kk)>=0.d0.and.        		     &
                  PCell%vofS(ii,jj,kk)>=epsi))then                           ! The neighbour cell is in gas and assigned dry cell            
                BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP           
              elseif((PCell%vof(ii,jj,kk)>=0.5d0.and.   		     &
                      PCell%vofS(ii,jj,kk)<epsi).or.    		     &
                     (PCell%phi(ii,jj,kk)<0.d0.and.     		     &
                      PCell%vofS(ii,jj,kk)>=epsi))then                       ! The neighbour cell is in liquid and assigned wet cell
                BetaW=BetaM
              end if
              TPoCoef(i,j,k)=PVel%Dp(im,jm,km)*SArea/DGrid*BetaP*BetaM/BetaW
            elseif((PCell%vof(i,j,k)<0.5d0.and.PCell%vofS(i,j,k)<epsi).or.   &
                   (PCell%phi(i,j,k)>=0.d0.and.PCell%vofS(i,j,k)>=epsi)) then ! The cell is in gas and it is assigned dry cell
              if((PCell%vof(ii,jj,kk)>=0.5d0.and.        		     &
                  PCell%vofS(ii,jj,kk)<epsi).or. 	     		     &
                 (PCell%phi(ii,jj,kk)<0.d0.and.        		             &
                  PCell%vofS(ii,jj,kk)>=epsi))then      
                BetaD=Lamda*BetaP+(1.d0-Lamda)*BetaM
              elseif((PCell%vof(ii,jj,kk)<0.5d0.and.        		     &
                      PCell%vofS(ii,jj,kk)<epsi).or. 	     		     &
                     (PCell%phi(ii,jj,kk)>=0.d0.and.        		     &
                      PCell%vofS(ii,jj,kk)>=epsi))then  
                BetaD=BetaP
              end if     
              TPoCoef(i,j,k)=PVel%Dp(im,jm,km)*SArea/DGrid*BetaP*BetaM/BetaD  
            end if
          end do
        end do
      end do    
      ! Set up value for variables at boundary cells 
      PVel%Dp(ium*1+iup*ISize,:,:)=PVel%Dp(ium*2+iup*(ISize-1),:,:)
      PVel%Dp(:,jvm*1+jvp*JSize,:)=PVel%Dp(:,jvm*2+jvp*(JSize-1),:)
      PVel%Dp(:,:,kwm*1+kwp*KSize)=PVel%Dp(:,:,kwm*2+kwp*(KSize-1))
    end subroutine Compute1DGFMCoefficient
    
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
          Call HYPRE_BoomerAMGSetMaxIter(precond,1,ierr)
!        set amg as the pcg preconditioner
!         precond_id = 2
          Call HYPRE_ParCSRPCGSetPrecond(solver,2,precond,ierr)
        End if
    End subroutine SetBasicSolver

    Subroutine SetPoissonMatrix(A,parcsr_A,PGrid,PCell,PoCoef,WB,EB,SB,NB,BB,TB)
        implicit none
        Integer*8,intent(inout)		   	  :: A,parcsr_A
        Type(Grid),intent(in)		   	  :: PGrid
        Type(Cell),intent(in)		   	  :: PCell
        real(kind=dp),dimension(:,:,:),intent(in) :: PU,PV,PW
        Integer(kind=it4b),intent(in)	   	  :: WB,EB,SB,NB,BB,TB ! boundary condition for west face, east face,
                                                          ! south face, north face, 0: Dirichlet, 1: Neumann
        Integer(kind=it4b)                 	  :: nnz,ictr,ilower,iupper,cols(0:6)
        Integer(kind=it4b)		   	  :: i,j,k
        Real(kind=dp)			          :: values(0:6)
        Real(kind=dp)			          :: dx,dy,dz
        ilower = 0
        iupper = PCell%ExtCell
      ! Create and Set up matrix
        Call HYPRE_IJMatrixCreate(MPI_COMM_WORLD,ilower,iupper,ilower,iupper,  &
                                                                        A,ierr)
        Call HYPRE_IJMatrixSetObjectType(A,HYPRE_PARCSR,ierr)
        Call HYPRE_IJMatrixInitialize(A,ierr)
!     Now go through my local rows and set the matrix entries.
!     Each row has at most 5 entries. For example, if n=3:
!
!      A = [M -I 0; -I M -I; 0 -I M]
!      M = [4 -1 0; -1 4 -1; 0 -1 4]
!
!     Note that here we are setting one row at a time, though
!     one could set all the rows together (see the User's Manual).
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              If(PCell%Cell_Type(i,j,k)/=2) then
                dx = PGrid%dx(i,j,k)
                dy = PGrid%dy(i,j,k)
                dz = PGrid%dz(i,j,k)
                ictr = Pcell%Posnu(i,j,k)
                nnz = 0
                values = 0.d0
                cols = 0
              ! West of current cell
                If(i>1) then
                  If(PCell%Posnu(i-1,j,k)/=-1) then
                    cols(nnz) = PCell%Posnu(i-1,j,k)
                    values(nnz) = -PoCoef(i,j,1)*PCell%EEArea(i-1,j,k)
                    nnz = nnz+1
                  End if
                End if
              ! South of current cell
                If(j>1) then
                  If(PCell%Posnu(i,j-1,k)/=-1) then
                    cols(nnz) = PCell%Posnu(i,j-1,k)
                    values(nnz) = -PoCoef(i,j,2)*PCell%NEArea(i,j-1,k)
                    nnz = nnz+1
                  End if
                End if
              ! Bottom of current cell
                If(k>1) then
                  If(PCell%Posnu(i,j,k-1)/=-1) then
                    cols(nnz) = PCell%Posnu(i,j,k-1)
                    values(nnz) = -PoCoef(i,j,3)*PCell%TEArea(i,j,k-1)
                    nnz = nnz+1
                  End if
                End if
              ! Set the diagonal cell
                cols(nnz) = PCell%Posnu(i,j,k)
                values(nnz) =  PoCoef(i,j,k,1)*PCell%EEArea(i-1,j,k)+          &
                               PoCoef(i,j,k,2)*PCell%NEArea(i,j-1,k)+          &
                               PoCoef(i,j,k,3)*PCell%TEArea(i,j,k-1)+          &
			       PoCoef(i,j,k,4)*PCell%EEArea(i,j,k)+            &
                               PoCoef(i,j,k,5)*PCell%NEArea(i,j,k)+            &
                               PoCoef(i,j,k,6)*PCell%TEArea(i,j,k) 
                If(isnan(values(nnz)).or.dabs(values(nnz))>1.d10) then
                  print*,values(nnz),i,j,k
                  pause 'project_167'
                End if
              !  If(i==Imax) values(nnz) = values(nnz)+1.d30
              ! Apply boundary condition for matrix
                If(WB==1.and.i==1) then
                  values(nnz) = values(nnz)-PoCoef(i,j,k,1)*PCell%EEArea(i-1,j,k)
                End if
                If(SB==1.and.j==1) then
                  values(nnz) = values(nnz)-PoCoef(i,j,k,2)*PCell%NEArea(i,j-1,k)
                End if
                If(BB==1.and.k==1) then
                  values(nnz) = values(nnz)-PoCoef(i,j,k,3)*PCell%TEArea(i,j,k-1)
                End if
                If(EB==1.and.i==Imax) then
                  values(nnz) = values(nnz)-PoCoef(i,j,k,4)*PCell%EEArea(i,j,k)
                End if
                If(NB==1.and.j==Jmax) then
                  values(nnz) = values(nnz)-PoCoef(i,j,k,5)*PCell%NEArea(i,j,k)
                End if
                If(TB==1.and.k==Kmax) then
                  values(nnz) = values(nnz)-PoCoef(i,j,k,6)*PCell%TEArea(i,j,k)
                End if        
                
                nnz = nnz+1
              ! East of current cell
                If(i<Imax) then
                  If(PCell%Posnu(i+1,j,k)/=-1) then
                    cols(nnz) = PCell%Posnu(i+1,j,k)
                    values(nnz) = -PoCoef(i,j,k,4)*PCell%EEArea(i,j,k)
                    nnz = nnz+1
                  End if
                End if
              ! North of current cell
                If(j<Jmax) then
                  If(PCell%Posnu(i,j+1,k)/=-1) then
                    cols(nnz) = PCell%Posnu(i,j+1,k)
                    values(nnz) = -PoCoef(i,j,k,5)*PCell%NEArea(i,j,k)
                    nnz = nnz+1
                  End if
                End if
              ! Top of current cell
                If(k<Kmax) then
                  If(PCell%Posnu(i,j,k+1)/=-1) then
                    cols(nnz) = PCell%Posnu(i,j,k+1)
                    values(nnz) = -PoCoef(i,j,k,6)*PCell%TEArea(i,j,k)
                    nnz = nnz+1
                  End if
                End if
                Call HYPRE_IJMatrixSetValues(A,1,nnz,ictr,cols,values,ierr)
              End if
            End do
          End do
        End do
        Call HYPRE_IJMatrixAssemble(A,ierr)
        Call HYPRE_IJMatrixGetObject(A,parcsr_A,ierr)
    End subroutine SetPoissonMatrix

    Subroutine SetPoissonVectors(b,x,par_b,par_x,PGrid,PCell,dt)
        Integer*8:: b,x,par_b,par_x
        Type(Grid),intent(in):: PGrid
        Type(Cell),intent(in):: PCell
        Real(kind=dp),intent(in):: dt
        Integer(kind=it4b):: i,j,k
        Integer:: ilower,iupper,ictr,local_size
        Real(kind=dp):: dx,dy,dz
        Integer(kind=it4b),dimension(:),allocatable:: rows
        Real(kind=dp),dimension(:),allocatable:: rhs,xval
        Real(kind=dp),dimension(:,:),allocatable:: ExtFlux
        ilower = 0
        iupper = PCell%ExtCell
        local_size = iupper-ilower+1 ! the number of rows
        ! In here, we apply boundary condition for deltaP with its values is 0 at
        ! all boundary. therefore, we do not need to set boundary in vector b
        Allocate(rhs(0:PCell%ExtCell))
        Allocate(xval(0:PCell%ExtCell))
        Allocate(rows(0:PCell%ExtCell))
        Call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,b,ierr)
        Call HYPRE_IJVectorSetObjectType(b,HYPRE_PARCSR,ierr)
        Call HYPRE_IJVectorInitialize(b,ierr)

        Call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,x,ierr)
        Call HYPRE_IJVectorSetObjectType(x,HYPRE_PARCSR,ierr)
        Call HYPRE_IJVectorInitialize(x,ierr)
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              If(PCell%Cell_Type(i,j,k)/=2) then
                dx = PGrid%dx(i,j,k)
                dy = PGrid%dy(i,j,k)
                dz = PGrid%dz(i,j,k)
                ictr = PCell%PosNu(i,j,k)
                rhs(ictr) = -dy*dz*(PCell%EEArea(i,j,k)*u(i,j,k)-             &
                                    PCell%EEArea(i-1,j,k)*u(i-1,j,k))         &
                            -dx*dz*(PCell%NEArea(i,j,k)*v(i,j,k)-             &
                                    PCell%NEArea(i,j-1,k)*v(i,j-1,k))         &
                            -dx*dy*(PCell%TEArea(i,j,k)*w(i,j,k)-             &
                                    PCell%TEArea(i,j,k-1)*w(i,j,k-1))
                xval(ictr) = 0.d0
                rows(ictr) = ilower+ictr
                If(isnan(rhs(ictr)).or.dabs(rhs(ictr))>1.d10) then
                  print*,i,j
                  pause 'fuck you 350'
                End if
              End if
            End do
          End do
        End do

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
    End subroutine SetPoissonVectors

    Subroutine DeltaPressureGetValues(x,PCell,Projp)
        Integer*8,intent(in):: x
        Type(Cell),intent(in):: PCell
        Type(Projection),intent(inout):: Projp
        Integer(kind=it4b):: i,j,k
        Integer(kind=it4b):: ilower,iupper,local_size,ctr
        Integer(kind=it4b),dimension(:),allocatable:: rows
        Real(kind=dp),dimension(:),allocatable:: values
        ilower = 0
        iupper = PCell%ExtCell
        local_size = PCell%ExtCell+1 ! number of element
        Allocate(values(ilower:iupper))
        Allocate(rows(ilower:iupper))
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              If(PCell%Cell_Type(i,j,k)/=2) then
                rows(PCell%PosNu(i,j,k)) = PCell%PosNu(i,j,k)+ilower
              End if
            End do
          End do
        End do
        Call HYPRE_IJVectorGetValues(x,local_size,rows,values,ierr)
        ctr = 0
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              If(PCell%PosNu(i,j,k)==ctr) then
                Projp%Pp(i,j,k) = values(ctr)
                ctr = ctr+1
                If(isnan(Projp%Pp(i,j,k)).or.projp%Pp(i,j,k)>1.d10) then
                  pause 'ProjectionP_Mod'
                End if
              Else
                Projp%Pp(i,j,k) = 0.d0
              End if
            End do
          End do
        End do
        Deallocate(values,rows)
    End subroutine DeltaPressureGetValues

    Subroutine DeltaPressureBoundaryCondition(Proj,WB,EB,SB,NB,BB,TB)
        Implicit none
        Type(Projection),intent(inout):: Proj
        Integer(kind=it4b),intent(in):: WB,EB,SB,NB,BB,TB
        Integer(kind=it4b):: i,j,k
        Do j = 1,Jmax
          Do k = 1,Kmax
            Proj%Pp(0,j,k) = dble(WB)*Proj%Pp(1,j,k)
            Proj%Pp(Imax+1,j,k) = dble(EB)*Proj%Pp(Imax,j,k)
          End do
        End do
        Do i = 1,Imax
          Do k = 1,Kmax
            Proj%Pp(i,0,k) = dble(SB)*Proj%Pp(i,1,k)
            Proj%Pp(i,Jmax+1,k) = dble(NB)*Proj%Pp(i,Jmax,k)
          End do
        End do
        Do i = 1,Imax
          Do j = 1,Jmax
            Proj%Pp(i,j,0) = dble(BB)*Proj%Pp(i,j,1)
            Proj%Pp(i,j,Kmax+1) = dble(TB)*Proj%Pp(i,j,Kmax)
          End do
        End do
    End subroutine DeltaPressureBoundaryCondition
End module ProjectionP

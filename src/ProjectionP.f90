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
    USE ieee_arithmetic
    use BoundaryInterface
    use BoundaryFunction
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
                               VCell,WCell,TVar,TPred,PU,PV,PW,BCp,PoCoef,Proj,dt)
        implicit none
        type(Grid),intent(in)		      :: PGrid,UGrid,VGrid,WGrid
        type(Cell),intent(in)		      :: PCell,UCell,VCell,WCell
        type(Variables),intent(in),target     :: TVar
        type(Predictor),intent(in),target     :: TPred
        type(PoissonCoefficient),intent(in)   :: PU,PV,PW
        type(BCBase),intent(in)		      :: BCp
        type(Projection),intent(inout),target :: Proj
      
        real(kind=dp),intent(in)	      :: dt
        integer*8			      :: A,parcsr_A,b,par_b,x,par_x
        integer*8			      :: solver,precond
        integer(kind=it4b)		      :: num_iterations
        real(kind=dp)			      :: final_res_norm,tol
        real(kind=dp),dimension(:,:,:,:),allocatable,intent(inout) :: PoCoef   ! the coefficient for Poisson solving
        
        call Compute1DGFMCoefficient(PGrid,PCell,UGrid,PU,row,	               &
                                                   1,0,0,0,0,0,PoCoef(:,:,:,1))
        call Compute1DGFMCoefficient(PGrid,PCell,VGrid,PV,row,	               &
                                                   0,1,0,0,0,0,PoCoef(:,:,:,2))
        call Compute1DGFMCoefficient(PGrid,PCell,WGrid,PW,row,	               &
                                                   0,0,1,0,0,0,PoCoef(:,:,:,3))
        call Compute1DGFMCoefficient(PGrid,PCell,UGrid,PU,row,	               &
                                                   0,0,0,1,0,0,PoCoef(:,:,:,4))
        call Compute1DGFMCoefficient(PGrid,PCell,VGrid,PV,row,                 &
                                                   0,0,0,0,1,0,PoCoef(:,:,:,5))
        call Compute1DGFMCoefficient(PGrid,PCell,WGrid,PW,row,                 &
                                                   0,0,0,0,0,1,PoCoef(:,:,:,6))
        p => TVar%p
        u => TPred%u
        v => TPred%v
        w => TPred%w  
        Call SetBasicSolver(solver,precond)
    !    Call SetBasicSolver(solver=solver,ierr=ierr)
        Call SetPoissonMatrix(A,parcsr_A,PGrid,PCell,PoCoef,		       &
                                                 BCp%flag(1),BCp%flag(2),      &
                                                 BCp%flag(3),BCp%flag(4),      &
                                                 BCp%flag(5),BCp%flag(6))
        Call SetPoissonVectors(b,x,par_b,par_x,PGrid,PCell,BCp,PoCoef,dt)
        Call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
        Call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    !    Run info - Needed logging turned on
        Call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
        call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
    !    Call HYPRE_ParCSRPCGGetResidual(solver,tol,ierr)
        Call DeltaPressureGetValues(x,PCell,Proj)
        Call DeltaPressureBoundaryCondition(Proj,BCp%flag(1),BCp%flag(2),      &
                                                 BCp%flag(3),BCp%flag(4),      &
                                                 BCp%flag(5),BCp%flag(6))
        Call HYPRE_IJMatrixDestroy(A,ierr)
        Call HYPRE_IJVectorDestroy(b,ierr)
        Call HYPRE_IJVectorDestroy(x,ierr)
        call HYPRE_BoomerAMGDestroy(precond,ierr)
        call HYPRE_ParCSRPCGDestroy(solver,ierr)
        Nullify(p)
        Nullify(u)
        Nullify(v)
        Nullify(w)
    end subroutine PoissonEquationSolver
    
    subroutine Compute1DGFMCoefficient(PGrid,PCell,TGrid,PVel,Roref,     &
                                               ium,jvm,kwm,iup,jvp,kwp,TPoCoef)
      type(Grid),intent(in)		  	   :: PGrid,TGrid
      type(Cell),intent(in)		  	   :: PCell
      type(PoissonCoefficient),intent(in) 	   :: PVel
      real(kind=dp),intent(in)			   :: Roref
      integer(kind=it4b),intent(in)	  	   :: ium,jvm,kwm,iup,jvp,kwp
      real(kind=dp),dimension(:,:,:)		   :: TPoCoef
      real(kind=dp)			  	   :: Lamda,BetaP,BetaM
      real(kind=dp)				   :: BetaW,BetaD
      integer(kind=it4b)			   :: i,j,k
      integer(kind=it4b)		  	   :: ii,jj,kk,im,jm,km
      real(kind=dp)			  	   :: DGrid
      
      BetaP = 1.d0/(row/Roref)
      BetaM = 1.d0/(roa/Roref)
      do i=1,Imax
        do j=1,JMax
          do k=1,KMax
            ii = i-ium+iup
            jj = j-jvm+jvp
            kk = k-kwm+kwp
            im = i-ium
            jm = j-jvm
            km = k-kwm
            DGrid = dble(ium+iup)*TGrid%dx(max(im,1),max(jm,1),max(km,1))+                          &
                    dble(jvm+jvp)*TGrid%dy(max(im,1),max(jm,1),max(km,1))+		               &
                    dble(kwm+kwp)*TGrid%dz(max(im,1),max(jm,1),max(km,1))
                    
            Lamda=dabs(PCell%phiL(i,j,k))/(dabs(PCell%phiL(i,j,k))+            &
                                           dabs(PCell%phiL(ii,jj,kk))+tol)                  
            if((PCell%vofL(i,j,k)>=0.5d0.and.                                  &
                PCell%vof(i,j,k)>1.d0-epsi).or.                                &
               (PCell%phiL(i,j,k)<0.d0.and.                                    &
                PCell%vof(i,j,k)<=1.d0-epsi.and.PCell%vof(i,j,k)>epsi)) then 	 ! The cell is in liquid and it is assigned wet cell 
              if((PCell%vofL(ii,jj,kk)<0.5d0.and.                              &
                  PCell%vof(ii,jj,kk)>1.d0-epsi).or.                           &
                 (PCell%phiL(ii,jj,kk)>=0.d0.and.                              &
                  PCell%vof(ii,jj,kk)<=1.d0-epsi.and.PCell%vof(ii,jj,kk)>epsi)) then  ! The neighbour cell is in gas and assigned dry cell            
                BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP           
              elseif((PCell%vofL(ii,jj,kk)>=0.5d0.and.   		                   &
                      PCell%vof(ii,jj,kk)>1.d0-epsi).or.    		               &
                     (PCell%phiL(ii,jj,kk)<0.d0.and.     		                   &
                  PCell%vof(ii,jj,kk)<=1.d0-epsi.and.PCell%vof(ii,jj,kk)>epsi)) then ! The neighbour cell is in liquid and assigned wet cell
                BetaW=BetaM
              end if
              TPoCoef(i,j,k)=PVel%Dp(im,jm,km)/DGrid*BetaP*BetaM/BetaW
              if(isnan(betaW).or.isnan(TPoCoef(i,j,k))) then
                print*,'print out values'
                print*,Lamda*BetaM,(1.d0-Lamda)*BetaP 
                print*,BetaM 
                print*,Lamda
                print*,i,j,k
                print*,ii,jj,kk
                print*,'++++++++'
                print*,BetaP,BetaM,BetaW,DGrid
                print*,TPoCoef(i,j,k),PVel%Dp(im,jm,km)
                print*,PCell%phiL(ii,jj,kk),PCell%vofL(ii,jj,kk)
                print*,PCell%phiL(ii,jj,kk),PCell%vof(ii,jj,kk)
                print*, 'Level set test'
                print*, PCell%phiL(i,j,k),PCell%phiL(ii,jj,kk)
                pause 'Test Lamda 142'
              end if  
            elseif((PCell%vofL(i,j,k)<0.5d0.and.			                         &
                    PCell%vof(i,j,k)>1.d0-epsi).or.   			                   &
                   (PCell%phiL(i,j,k)>=0.d0.and.			                         &
                    PCell%vof(i,j,k)<=1.d0-epsi.and.PCell%vof(i,j,k)>epsi)) then  ! The cell is in gas and it is assigned dry cell
              if((PCell%vofL(ii,jj,kk)>=0.5d0.and.        		       &
                  PCell%vof(ii,jj,kk)>1.d0-epsi).or. 	     		       &
                 (PCell%phiL(ii,jj,kk)<0.d0.and.        		       &
                  PCell%vof(ii,jj,kk)<=1.d0-epsi.and.PCell%vof(ii,jj,kk)>epsi)) then  ! The neighbour cell is in liquid and assigned wet cell     
                BetaD=Lamda*BetaP+(1.d0-Lamda)*BetaM
              elseif((PCell%vofL(ii,jj,kk)<0.5d0.and.        		       &
                      PCell%vof(ii,jj,kk)>1.d0-epsi).or. 	     	       &
                     (PCell%phiL(ii,jj,kk)>=0.d0.and.        		       &
                  PCell%vof(ii,jj,kk)<=1.d0-epsi.and.PCell%vof(ii,jj,kk)>epsi)) then  ! The neighbour cell is in gas and assigned dry cell    
                BetaD=BetaP
              end if     
              TPoCoef(i,j,k)=PVel%Dp(im,jm,km)/DGrid*BetaP*BetaM/BetaD 
              if(isnan(betaD).or.isnan(TPoCoef(i,j,k))) then
                print*,Lamda
                print*,i,j,k
                print*,ii,jj,kk
                print*,'++++++++'
                print*,BetaD,TPoCoef(i,j,k)
                print*,PCell%phiL(ii,jj,kk),PCell%vofL(ii,jj,kk)
                print*,PCell%phi(ii,jj,kk),PCell%vof(ii,jj,kk)
                pause 'Test Lamda 168'
              end if  
            end if   
          end do
        end do
      end do    
      ! Set up value for poisson coefficients at boundary cells 
      ! We need a better method to compute poisson coefficients at boundary
      ! Now, we assume the poisson coefficient at boudnary is equal to the neighbour cell
      ! The better approach is calculate these values based on geometry shape at pressure cell. 
      if(ium==1.or.iup==1) then
        TPoCoef(ium*1+iup*Imax,:,:)=TPoCoef(ium*2+iup*(Imax-1),:,:)
      end if
      if(jvm==1.or.jvp==1) then  
        TPoCoef(:,jvm*1+jvp*JMax,:)=TPoCoef(:,jvm*2+jvp*(JMax-1),:)
      end if
      if(kwm==1.or.kwp==1) then
        TPoCoef(:,:,kwm*1+kwp*KMax)=TPoCoef(:,:,kwm*2+kwp*(KMax-1))
      end if  
    end subroutine Compute1DGFMCoefficient
    
    Subroutine SetBasicSolver(solver,precond)
        Implicit none
        Integer*8,intent(inout):: solver
        Integer*8,intent(inout),optional:: precond
!       Set up and use a solver
        Call HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD,solver,ierr)
!       Set some parameters
        Call HYPRE_ParCSRPCGSetMaxIter(solver,50,ierr)
        Call HYPRE_ParCSRPCGSetTol(solver,1.d-30,ierr)
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

    Subroutine SetPoissonMatrix(A,parcsr_A,PGrid,PCell,PoCoef,WB,EB,SB,NB,BB,TB)
        implicit none
        Integer*8,intent(inout)		   	  :: A,parcsr_A
        Type(Grid),intent(in)		   	  :: PGrid
        Type(Cell),intent(in)		   	  :: PCell
        Integer(kind=it4b),intent(in)	   	  :: WB,EB,SB,NB,BB,TB ! boundary condition for west face, east face,
                                                                       ! south face, north face, 0: Dirichlet, 1: Neumann
        real(kind=dp),dimension(:,:,:,:),allocatable :: PoCoef
        Integer(kind=it4b)                 	  :: nnz,ictr,ilower,iupper,cols(0:6)
        Integer(kind=it4b)		   	  :: i,j,k
        Real(kind=dp)			          :: values(0:6)
        Real(kind=dp)			          :: dx,dy,dz,mp,mp1,mp2,mp3,mp4,mp5,mp6
        
        ilower = 0
        iupper = PCell%ExtCell
      ! Create and Set up matrix
        Call HYPRE_IJMatrixCreate(MPI_COMM_WORLD,ilower,iupper,ilower,iupper,  &
                                                                        A,ierr)
        Call HYPRE_IJMatrixSetObjectType(A,HYPRE_PARCSR,ierr)
        Call HYPRE_IJMatrixInitialize(A,ierr)
      ! Now go through my local rows and set the matrix entries. 
      ! Each row has at most 5 entries. For example, if n=3:
      !
      !      A = [M -I 0; -I M -I; 0 -I M]
      !      M = [4 -1 0; -1 4 -1; 0 -1 4]
      !
      ! Note that here we are setting one row at a time, though
      ! one could set all the rows together (see the User's Manual).
        
        do i = 1,Imax
          do j = 1,Jmax
            do k = 1,Kmax
              if(PCell%Cell_Type(i,j,k)/=2) then
                mp1=0.d0
        	mp2=0.d0
       		mp3=0.d0
        	mp4=0.d0
        	mp5=0.d0
        	mp6=0.d0
                dx = PGrid%dx(i,j,k)
                dy = PGrid%dy(i,j,k)
                dz = PGrid%dz(i,j,k)
                ictr = Pcell%Posnu(i,j,k)
                nnz = 0
                values = 0.d0
                cols = 0
              ! West of current cell
                if(i>1) then
                  if(PCell%Posnu(i-1,j,k)/=-1) then
                    cols(nnz)=PCell%Posnu(i-1,j,k)
                    values(nnz)=-PoCoef(i,j,k,1)*PCell%EEArea(i-1,j,k)*        &
                                 PGrid%dy(i,j,k)*PGrid%dz(i,j,k)
                    mp1=dabs(values(nnz))
                    nnz=nnz+1
                  end if
                end if
              ! South of current cell
                if(j>1) then
                  if(PCell%Posnu(i,j-1,k)/=-1) then
                    cols(nnz)=PCell%Posnu(i,j-1,k)
                    values(nnz)=-PoCoef(i,j,k,2)*PCell%NEArea(i,j-1,k)*        &
                                 PGrid%dx(i,j,k)*PGrid%dz(i,j,k)
                    mp2=dabs(values(nnz))             
                    nnz=nnz+1
                  end if
                end if
              ! Bottom of current cell
                if(k>1) then
                  if(PCell%Posnu(i,j,k-1)/=-1) then
                    cols(nnz)=PCell%Posnu(i,j,k-1)
                    values(nnz)=-PoCoef(i,j,k,3)*PCell%TEArea(i,j,k-1)*	       &
                                 PGrid%dx(i,j,k)*PGrid%dy(i,j,k)
                    mp3=dabs(values(nnz))	
                    nnz=nnz+1
                  end if
                end if
              ! Set the diagonal cell
                cols(nnz)=PCell%Posnu(i,j,k)
                values(nnz)=PoCoef(i,j,k,1)*PCell%EEArea(i-1,j,k)*             &
                            PGrid%dy(i,j,k)*PGrid%dz(i,j,k)+                   &
                            PoCoef(i,j,k,2)*PCell%NEArea(i,j-1,k)*             &
                            PGrid%dx(i,j,k)*PGrid%dz(i,j,k)+                   &
                            PoCoef(i,j,k,3)*PCell%TEArea(i,j,k-1)*             &
                            PGrid%dx(i,j,k)*PGrid%dy(i,j,k)+                   &
                            PoCoef(i,j,k,4)*PCell%EEArea(i,j,k)*               &
                            PGrid%dy(i,j,k)*PGrid%dz(i,j,k)+                   &
                            PoCoef(i,j,k,5)*PCell%NEArea(i,j,k)*               &
                            PGrid%dx(i,j,k)*PGrid%dz(i,j,k)+                   & 
                            PoCoef(i,j,k,6)*PCell%TEArea(i,j,k)*               &
                            PGrid%dx(i,j,k)*PGrid%dy(i,j,k) 
                mp=dabs(values(nnz))
                If(isnan(values(nnz)).or.dabs(values(nnz))>1.d10) then
                  print*,'000000000000000000000000000000000000000000000000000000'
                  print*,values(nnz)
                  print*,i,j,k
                  print*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                  print*,PoCoef(i,j,k,1)
                  print*,PoCoef(i,j,k,2)
                  print*,PoCoef(i,j,k,3)
                  print*,PoCoef(i,j,k,4)
                  print*,PoCoef(i,j,k,5)
                  print*,PoCoef(i,j,k,6)
                  pause 
                End if
              ! If(i==Imax) values(nnz) = values(nnz)+1.d30
              ! Apply boundary condition for matrix
                If(WB==1.and.i==1) then
                  values(nnz)=values(nnz)-PoCoef(i,j,k,1)*PCell%EEArea(i-1,j,k)*&
                              PGrid%dy(i,j,k)*PGrid%dz(i,j,k)
                End if
                If(SB==1.and.j==1) then
                  values(nnz)=values(nnz)-PoCoef(i,j,k,2)*PCell%NEArea(i,j-1,k)*&
                              PGrid%dx(i,j,k)*PGrid%dz(i,j,k)
                End if
                If(BB==1.and.k==1) then
                  values(nnz)=values(nnz)-PoCoef(i,j,k,3)*PCell%TEArea(i,j,k-1)*&
                              PGrid%dx(i,j,k)*PGrid%dy(i,j,k)
                End if
                If(EB==1.and.i==Imax) then
                  values(nnz)=values(nnz)-PoCoef(i,j,k,4)*PCell%EEArea(i,j,k)*  &
                              PGrid%dy(i,j,k)*PGrid%dz(i,j,k)
                End if
                If(NB==1.and.j==Jmax) then
                  values(nnz)=values(nnz)-PoCoef(i,j,k,5)*PCell%NEArea(i,j,k)*  &
                              PGrid%dx(i,j,k)*PGrid%dz(i,j,k) 
                End if
                If(TB==1.and.k==Kmax) then
                  values(nnz)=values(nnz)-PoCoef(i,j,k,6)*PCell%TEArea(i,j,k)*  &
                              PGrid%dx(i,j,k)*PGrid%dy(i,j,k)
                End if                 
                nnz = nnz+1
              ! East of current cell
                If(i<Imax) then
                  If(PCell%Posnu(i+1,j,k)/=-1) then
                    cols(nnz)=PCell%Posnu(i+1,j,k)
                    values(nnz)=-PoCoef(i,j,k,4)*PCell%EEArea(i,j,k)*	       &
                                 PGrid%dy(i,j,k)*PGrid%dz(i,j,k)
                    mp4=dabs(values(nnz))             	
                    nnz = nnz+1
                  End if
                End if
              ! North of current cell
                If(j<Jmax) then
                  If(PCell%Posnu(i,j+1,k)/=-1) then
                    cols(nnz)=PCell%Posnu(i,j+1,k)
                    values(nnz)=-PoCoef(i,j,k,5)*PCell%NEArea(i,j,k)*	       &
                                 PGrid%dx(i,j,k)*PGrid%dz(i,j,k)
                    mp5=dabs(values(nnz))              	
                    nnz = nnz+1
                  End if
                End if
              ! Top of current cell
                If(k<Kmax) then
                  If(PCell%Posnu(i,j,k+1)/=-1) then
                    cols(nnz)=PCell%Posnu(i,j,k+1)
                    values(nnz)=-PoCoef(i,j,k,6)*PCell%TEArea(i,j,k)*	       &
                                 PGrid%dx(i,j,k)*PGrid%dy(i,j,k)	
                    mp6=dabs(values(nnz))
                    nnz=nnz+1
                  End if
                End if
                Call HYPRE_IJMatrixSetValues(A,1,nnz,ictr,cols,values,ierr)
                if(mp<mp1+mp2+mp3+mp4+mp5+mp6.and.i>1.and.j>1.and.k>1) then
                  print*,'undominant coefficient'
                  print*, mp1,mp2,mp3
                  print*, mp4,mp5,mp6
                  print*, mp1+mp2+mp3+mp4+mp5+mp6
                  print*, mp
                  print*, i,j,k
                  print*, 'tttttttttttttttttttt'
                  print*, PoCoef(i,j,k,1),PoCoef(i,j,k,2),PoCoef(i,j,k,3)
                  print*, PoCoef(i,j,k,4),PoCoef(i,j,k,5),PoCoef(i,j,k,6)
                  pause 
                end if
              End if
            End do
          End do
        End do
        
        Call HYPRE_IJMatrixAssemble(A,ierr)
        Call HYPRE_IJMatrixGetObject(A,parcsr_A,ierr)
    End subroutine SetPoissonMatrix

    Subroutine SetPoissonVectors(b,x,par_b,par_x,PGrid,PCell,BCp,PoCoef,dt)
        Integer*8:: b,x,par_b,par_x
        Type(Grid),intent(in)   				:: PGrid
        Type(Cell),intent(in)   				:: PCell
        type(BCBase),intent(in) 				:: BCp
        real(kind=dp),dimension(:,:,:,:),allocatable,intent(in) :: PoCoef
        Real(kind=dp),intent(in)				:: dt
        Integer(kind=it4b)					:: i,j,k
        Integer							:: ilower,iupper,ictr,local_size
        Real(kind=dp)						:: dx,dy,dz
        Integer(kind=it4b),dimension(:),allocatable		:: rows
        Real(kind=dp),dimension(:),allocatable			:: rhs,xval
        Real(kind=dp),dimension(:,:),allocatable		:: ExtFlux
        ilower = 0
        iupper = PCell%ExtCell
        local_size = iupper-ilower+1 ! the number of rows
        ! In here, we apply boundary condition for deltaP with its values is 0 at
        ! all boundary. therefore, we do not need to set boundary in vector b
        allocate(rhs(0:PCell%ExtCell))
        allocate(xval(0:PCell%ExtCell))
        allocate(rows(0:PCell%ExtCell))
        call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,b,ierr)
        call HYPRE_IJVectorSetObjectType(b,HYPRE_PARCSR,ierr)
        call HYPRE_IJVectorInitialize(b,ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,x,ierr)
        call HYPRE_IJVectorSetObjectType(x,HYPRE_PARCSR,ierr)
        call HYPRE_IJVectorInitialize(x,ierr)
        do i = 1,Imax
          do j = 1,Jmax
            do k = 1,Kmax
              if(PCell%Cell_Type(i,j,k)/=2) then
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
                                    
                if(BCp%flag(1)==0.and.i==1) then
                  rhs(ictr)=rhs(ictr)+PoCoef(i,j,k,1)*BCp%VarW(j,k)
                end if
                if(BCp%flag(2)==0.and.i==Imax) then
                  rhs(ictr)=rhs(ictr)+PoCoef(i,j,k,4)*BCp%VarE(j,k)  
                end if
                if(BCp%flag(3)==0.and.j==1) then
                  rhs(ictr)=rhs(ictr)+PoCoef(i,j,k,2)*BCp%VarS(i,k)
                end if
                if(BCp%flag(4)==0.and.j==Jmax) then
                  rhs(ictr)=rhs(ictr)+PoCoef(i,j,k,5)*BCp%VarN(i,k)
                end if
                if(BCp%flag(5)==0.and.k==1) then
                  rhs(ictr)=rhs(ictr)+PoCoef(i,j,k,3)*BCp%VarB(i,j)
                end if
                if(BCp%flag(6)==0.and.k==Kmax) then
                  rhs(ictr)=rhs(ictr)+PoCoef(i,j,k,6)*BCp%VarT(i,j)
                end if
                        
                xval(ictr) = 0.d0
                rows(ictr) = ilower+ictr
                
                if(isnan(rhs(ictr)).or.dabs(rhs(ictr))>1.d10) then
                  print*,i,j
                  print*,u(i,j,k)
                  print*,u(i-1,j,k)
                  print*,v(i,j,k)
                  print*,v(i,j-1,k)
                  print*,w(i,j,k)
                  print*,w(i,j,k-1)
                  pause 'fuck you 350'
                end if
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
    end subroutine SetPoissonVectors

    Subroutine DeltaPressureGetValues(x,PCell,Projp)
        Integer*8,intent(in)			    :: x
        Type(Cell),intent(in)			    :: PCell
        Type(Projection),intent(inout)		    :: Projp
        Integer(kind=it4b)			    :: i,j,k
        Integer(kind=it4b)			    :: ilower,iupper,local_size,ctr
        Integer(kind=it4b),dimension(:),allocatable :: rows
        Real(kind=dp),dimension(:),allocatable	    :: values
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

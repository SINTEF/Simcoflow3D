Module Clsvof
   !! Description:
   !! The module tracks the interface movement. The transport equation for volume of fluid
   !! and level set function is solved to get their value
   !! Method:
   !! The spilting operator is applied for discretizing the transport equations.
   !! It means that equations will be solved in each direction seperately.
   ! Current Code Owner: SIMCOFlow
   ! Code Description:
   ! Language: Fortran 90.
   ! Software Standards: "European Standards for Writing and
   ! Documenting Exchangeable Fortran 90 Code".
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Author : Son Tung Dang
   !        : NTNU,SINTEF
   ! Date : 20.09.2019
    USE ieee_arithmetic
    USE Matrix
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE BoundaryFunction

    implicit none
    private

    integer,parameter                      :: band_width = 4,nv = 6,nl = 6
    real(dp),parameter                     :: eta=0.075d0,vofeps=1.d-14,toldeno=1.d-14
    ! Parameters for setting normal vector to avoid underflow problem
    real(kind=dp),parameter                :: n1=1.d0-1.d-14, n2=dsqrt((2.d0*1.d-14-1.d-28)/2.d0)
    real(kind=dp),parameter                :: tolpar=1.d-14

    public:: InitialClsvofFluidField,InitialClsvofLiquidField,Clsvof_Scheme,   &
             ComputeUVWLiquidField,volume_fraction_calc,ConvertTheCoordinateSystem, &
             BoundaryConditionLvsVof, BoundaryConditionLvsVofFluid,InitialiseClsvof

    interface InitialClsvofFluidField
      module procedure InitialClsvofFluidField
    end interface

    interface InitialClsvofLiquidField
      module procedure InitialClsvofLiquidField
    end interface

    interface ClsVof_Scheme
      module procedure ClsVof_Scheme
    end interface

    interface ComputeUVWLiquidField
      module procedure ComputeUVWLiquidField
    end interface

    interface volume_fraction_calc
      module procedure volume_fraction_calc
    end interface

    interface ConvertTheCoordinateSystem
      module procedure ConvertTheCoordinateSystem
    end interface
    interface BoundaryConditionLvsVof
      module procedure BoundaryConditionLvsVof
    end interface    

    interface BoundaryConditionLvsVofFluid
      module procedure BoundaryConditionLvsVofFluid
    end interface 
 
    contains

    subroutine initialiseClsvof(PGrid,PCell,TVar,BCLvs,BCvof)

        Implicit none 
        Type(Grid),intent(in)         		   :: PGrid
        Type(Cell),intent(inout)      		   :: PCell
        Type(Variables),intent(inout) 		   :: TVar
        type(BCBase),intent(inout)		       :: BCVof,BCLvs

      !  nullify(vfl,vflF)
      !  nullify(phi,phiF)
      !  nullify(nxF,nyF,nzF)

      !  call BoundaryConditionLvsVof(PGrid,PCell,TVar,BCLvs,BCVof,0.d0)

    end subroutine initialiseClsvof

    subroutine InitialClsvofFluidField(TGrid,TCell)
      type(Grid),       intent(in)           :: TGrid
      type(Cell),       intent(inout),target :: TCell
      integer(kind=it4b)	                   :: i,j,k
      real(kind=dp)		                       :: dx,dy,dz,dis,vol,Radius,epsi,s
      real(kind=dp)		                       :: tol
      real(kind=dp),dimension(:,:,:),pointer :: vflF    ! vfl represents the liquid volume fraction, vflF represents fluid volume fraction
      real(kind=dp),dimension(:,:,:),pointer :: phiF    ! phi represents the liquid level set function, phiF represents fluid level set function
      real(kind=dp),dimension(:,:,:),pointer :: nxF,nyF,nzF ! nxF,nyF,nzF is fluid normal vector

      tol = 1.d-20
      epsi = 1.d-40
      vflF => TCell%vof
      phiF => TCell%phi
      nxF => TCell%nx
      nyF => TCell%ny
      nzF => TCell%nz
      Radius = 0.5d0/TGrid%Lref
      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax
            dx=TGrid%x(i,j,k)!-100.d0
            dy=TGrid%y(i,j,k)
            dz=TGrid%z(i,j,k)
            dis=dsqrt(dx**2.d0+dy**2.d0+dz**2.d0)-Radius
            nxF(i,j,k)=dx/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            nyF(i,j,k)=dy/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            nzF(i,j,k)=dz/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
            s=dis+0.5*(dabs(nxF(i,j,k))*TGrid%dx(i,j,k)+dabs(nyF(i,j,k))*      &
                       TGrid%dy(i,j,k)+dabs(nzF(i,j,k))*TGrid%dz(i,j,k))
            call Volume_Fraction_Calc(TGrid%dx(i,j,k),TGrid%dy(i,j,k),         &
                 TGrid%dz(i,j,k),nxF(i,j,k),nyF(i,j,k),nzF(i,j,k),s,vol)
            vflF(i,j,k)=vol/(TGrid%dx(i,j,k)*TGrid%dy(i,j,k)*TGrid%dz(i,j,k))
            phiF(i,j,k)=dis
          end do
        end do
      end do
      nullify(vflF)
      nullify(phiF)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
    end subroutine InitialClsvofFluidField

    subroutine InitialClsvofLiquidField(TGrid,TCell)
      implicit none
      type(Grid),       intent(in)           :: TGrid
      type(Cell),       intent(inout),target :: TCell
      integer(kind=it4b)	                   :: i,j,k
      real(kind=dp)		                       :: dx,dy,dz,dis,vol,Radius,epsi,s
      real(kind=dp)		                       :: tol
      real(kind=dp),dimension(:,:,:),pointer :: vfl    ! vfl represents the liquid volume fraction, vflF represents fluid volume fraction
      real(kind=dp),dimension(:,:,:),pointer :: phi    ! phi represents the liquid level set function, phiF represents fluid level set function
      real(kind=dp),dimension(:,:,:),pointer :: nxL,nyL,nzL ! nxF,nyF,nzF is fluid normal vector

      vfl => TCell%vofL
      phi => TCell%phiL
      nxL => TCell%nxL
      nyL => TCell%nyL
      nzL => TCell%nzL

      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax
            vfl(i,j,k)=1.d0*TCell%vof(i,j,k)
            phi(i,j,k)=-1.d3
            nxL(i,j,k)=0.d0
            nyL(i,j,k)=0.d0
            nzL(i,j,k)=1.d0
          end do
        end do
      end do

      nullify(vfl)
      nullify(phi)
      nullify(nxL)
      nullify(nyL)
      nullify(nzL)
    end subroutine InitialClsvofLiquidField

    subroutine ComputeUVWLiquidField(PGrid,PCell,UCell,VCell,WCell,            &
                                                 UGrid,VGrid,WGrid)
      !! The subroutine is used to compute the liquid volume fraction for UCell,VCell,WCell
      !! based on the volume fraction from PCell.
      implicit none
      type(Grid),intent(in)	:: PGrid
      type(Cell),intent(in)     :: PCell
      type(Cell),intent(inout)  :: UCell,VCell,WCell
      type(Grid),intent(in)     :: UGrid,VGrid,WGrid
      integer(kind=it4b)        :: i,j,k
      real(kind=dp)           :: VofLeft,VofRight,Dist

      call ComputeVelocityCellLiquidField(PGrid%dx,PCell,UCell,1,0,0)
      ! Compute the liquid volume fraction for UCell
      do j=1,Jmax
        do k=1,Kmax
          do i=1,Imax-1
            ! Compute the left volume fraction
            Dist=ConvertTheCoordinateSystem(PGrid%dx(i,j,k),PGrid%dy(i,j,k),   &
                 PGrid%dz(i,j,k),PCell%nxL(i,j,k),PCell%nyL(i,j,k),            &
                 PCell%nzL(i,j,k),PCell%phiL(i,j,k))
            call Volume_Fraction_Calc(PGrid%dx(i,j,k)/2.d0,PGrid%dy(i,j,k),    &
                                      PGrid%dz(i,j,k),PCell%nxL(i,j,k),        &
                   PCell%nyL(i,j,k),PCell%nzL(i,j,k),Dist,VofLeft)
            if(PCell%nxL(i,j,k)>=0.d0) then
              VofLeft=(PCell%vofL(i,j,k)*PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*      &
                                         PGrid%dz(i,j,k)-VofLeft)
            end if
            ! Compute the right volume fraction
            Dist=ConvertTheCoordinateSystem(PGrid%dx(i+1,j,k),                 &
                 PGrid%dy(i+1,j,k),PGrid%dz(i+1,j,k),PCell%nxL(i+1,j,k),       &
                 PCell%nyL(i+1,j,k),PCell%nzL(i+1,j,k),PCell%phiL(i+1,j,k))
            call Volume_Fraction_Calc(PGrid%dx(i+1,j,k)/2.d0,PGrid%dy(i+1,j,k),&
                                      PGrid%dz(i+1,j,k),PCell%nxL(i+1,j,k),    &
                PCell%nyL(i+1,j,k),PCell%nzL(i+1,j,k),Dist,VofRight)
            if(PCell%nxL(i+1,j,k)<=0.d0) then
              VofRight=(PCell%VofL(i+1,j,k)*PGrid%dx(i+1,j,k)*                 &
                        PGrid%dy(i+1,j,k)*PGrid%dz(i+1,j,k)-VofRight)
            end if
            UCell%vofL(i,j,k)=(VofLeft+VofRight)/UGrid%dx(i,j,k)/              &
                                                 UGrid%dy(i,j,k)/UGrid%dz(i,j,k)  
          end do
          UCell%vofL(Imax,j,k)=PCell%vofL(Imax,j,k)
        end do
      end do    
      call ComputeVelocityCellLiquidField(PGrid%dy,PCell,VCell,0,1,0)
      ! Compute the liquid volume fraction for VCell
      do i=1,Imax
        do k=1,Kmax
          do j=1,Jmax-1
            ! Compute the left volume fraction
            Dist=ConvertTheCoordinateSystem(PGrid%dx(i,j,k),PGrid%dy(i,j,k),   &
                 PGrid%dz(i,j,k),PCell%nxL(i,j,k),PCell%nyL(i,j,k),            &
                 PCell%nzL(i,j,k),PCell%phiL(i,j,k))
            call Volume_Fraction_Calc(PGrid%dx(i,j,k),PGrid%dy(i,j,k)/2.d0,    &
                                      PGrid%dz(i,j,k),PCell%nxL(i,j,k),        &
                       PCell%nyL(i,j,k),PCell%nzL(i,j,k),Dist,VofLeft)
            if(PCell%nyL(i,j,k)>=0.d0) then
              VofLeft=(PCell%vofL(i,j,k)*PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*      &
                                         PGrid%dz(i,j,k)-VofLeft) 
            end if  
            ! Compute the right volume fraction
            Dist=ConvertTheCoordinateSystem(PGrid%dx(i,j+1,k),                 &
                 PGrid%dy(i,j+1,k),PGrid%dz(i,j+1,k),PCell%nxL(i,j+1,k),       &
                 PCell%nyL(i,j+1,k),PCell%nzL(i,j+1,k),PCell%phiL(i,j+1,k))
            call Volume_Fraction_Calc(PGrid%dx(i,j+1,k),PGrid%dy(i,j+1,k)/2.d0,&
                                      PGrid%dz(i,j+1,k),PCell%nxL(i,j+1,k),    &
                       PCell%nyL(i,j+1,k),PCell%nzL(i,j+1,k),Dist,VofRight) 
            if(PCell%nyL(i,j+1,k)<=0.d0) then
              VofRight=(PCell%vofL(i,j+1,k)*PGrid%dx(i,j+1,k)*                 &
                                 PGrid%dy(i,j+1,k)*PGrid%dz(i,j+1,k)-VofRight)
            endif  
            VCell%vofL(i,j,k)=(VofLeft+VofRight)/VGrid%dx(i,j,k)/              &
                                   VGrid%dy(i,j,k)/VGrid%dz(i,j,k) 
          end do
          VCell%vofL(i,Jmax,k)=PCell%vofL(i,Jmax,k)
        end do
      end do      
      call ComputeVelocityCellLiquidField(PGrid%dz,PCell,WCell,0,0,1)
      ! Compute the liquid volume fraction for WCell
      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax-1
          ! Compute the left volume fraction
            Dist=ConvertTheCoordinateSystem(PGrid%dx(i,j,k),PGrid%dy(i,j,k),   &
                 PGrid%dz(i,j,k),PCell%nxL(i,j,k),PCell%nyL(i,j,k),            &
                 PCell%nzL(i,j,k),PCell%phiL(i,j,k))
            call Volume_Fraction_Calc(PGrid%dx(i,j,k),PGrid%dy(i,j,k),         &
                                      PGrid%dz(i,j,k)/2.d0,PCell%nxL(i,j,k),   &
                       PCell%nyL(i,j,k),PCell%nzL(i,j,k),Dist,VofLeft)
            if(PCell%nzL(i,j,k)>=0.d0) then
              VofLeft=(PCell%VofL(i,j,k)*PGrid%dx(i,j,k)*PGrid%dy(i,j,k)*      &
                                         PGrid%dz(i,j,k)-VofLeft)
            end if
          ! Compute the right volume fraction
            Dist=ConvertTheCoordinateSystem(PGrid%dx(i,j,k+1),                 &
                 PGrid%dy(i,j,k+1),PGrid%dz(i,j,k+1),PCell%nxL(i,j,k+1),       &
                 PCell%nyL(i,j,k+1),PCell%nzL(i,j,k+1),PCell%phiL(i,j,k+1))
            call Volume_Fraction_Calc(PGrid%dx(i,j,k+1),PGrid%dy(i,j,k+1),     &
                                  PGrid%dz(i,j,k+1)/2.d0,PCell%nxL(i,j,k+1),   &
                PCell%nyL(i,j,k+1),PCell%nzL(i,j,k+1),Dist,VofRight) 
            if(PCell%nzL(i,j,k)<0.d0) then
              VofRight=(PCell%vofL(i,j,k+1)*PGrid%dx(i,j,k+1)*                 &
                        PGrid%dy(i,j,k+1)*PGrid%dz(i,j,k+1)-VofRight)                                                
            end if
            WCell%vofL(i,j,k)=(VofLeft+VofRight)/WGrid%dx(i,j,k)/              &
                                   WGrid%dy(i,j,k)/WGrid%dz(i,j,k)
          end do
          WCell%vofl(i,j,Kmax)=PCell%vofL(i,j,Kmax)
        end do
      end do      
    end subroutine ComputeUVWLiquidField

    subroutine ComputeVelocityCellLiquidField(dxyz,PCell,VelCell,iu,iv,iw)
      !! The subroutine is used to compute the velocity cell properties such as
      !! volume of fluid, level set function, normal vector. The boundary conditions
      !! are applied later.
      implicit none
      real(kind=dp),dimension(:,:,:), intent(in)    :: dxyz
      type(Cell),                     intent(in)    :: PCell
      type(Cell),                     intent(inout) :: VelCell
      integer(kind=it4b),             intent(in)    :: iu,iv,iw
      integer(kind=it4b)			:: i,j,k

    !  call DirectionAverageArray(dxyz,PCell%vofL,VelCell%vofL,iu,iv,iw)
      call DirectionAverageArray(dxyz,PCell%phiL,VelCell%phiL,iu,iv,iw)
      call DirectionAverageArray(dxyz,PCell%nxL,VelCell%nxL,iu,iv,iw)
      call DirectionAverageArray(dxyz,PCell%nyL,VelCell%nyL,iu,iv,iw)
      call DirectionAverageArray(dxyz,PCell%nzL,VelCell%nzL,iu,iv,iw)

      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax
            if(VelCell%vof(i,j,k)+VelCell%vofL(i,j,k)>1.d0) then
              VelCell%vofL(i,j,k)=1.d0-VelCell%vof(i,j,k)
            end if
          end do
        end do
      end do
    end subroutine ComputeVelocityCellLiquidField

    subroutine DirectionAverageArray(dxyz,Varin,Varout,iu,iv,iw)
      !! The subroutine is used to compute the variables based on centre average.
      !! More complex interpolation technique can be used later.
      implicit none
      real(kind=dp),dimension(:,:,:),intent(in)  :: dxyz,Varin
      real(kind=dp),dimension(:,:,:),intent(out) :: Varout
      integer(kind=it4b),            intent(in)  :: iu,iv,iw
      real(kind=dp)				     :: lamda
      integer(kind=it4b)			 :: i,j,k

      do i=1,Imax-iu
        do j=1,Jmax-iv
          do k=1,Kmax-iw
            lamda=dxyz(i+iu,j+iv,k+iw)/(dxyz(i,j,k)+dxyz(i+iu,j+iv,k+iw))
            Varout(i,j,k)=lamda*Varin(i,j,k)+(1.d0-lamda)*Varin(i+iu,j+iv,k+iw)
          end do
        end do
      end do

      if(iu==1) Varout(Imax,1:Jmax,1:Kmax)=Varout(Imax-1,1:Jmax,1:Kmax)
      if(iv==1) Varout(1:Imax,Jmax,1:Kmax)=Varout(1:Imax,Jmax-1,1:Kmax)
      if(iw==1) Varout(1:Imax,1:Jmax,Kmax)=Varout(1:Imax,1:Jmax,Kmax-1)
    end subroutine DirectionAverageArray

    subroutine Clsvof_Scheme(PGrid,PCell,TVar,BCu,BCv,BCw,BCLvs,BCvof,Time,dt,itt)
      !! The subroutine is used to track the interface using the couple leve set volume of fluid
      implicit none
      type(Grid),intent(in)           :: PGrid
      !! The pressure grid
      type(Cell),intent(inout),target :: PCell
      !! The pressure cell
      type(Variables),intent(in)      :: TVar
      !! The StateVariables
      type(BCBase),intent(in)         :: BCu,BCv,BCw
      !! The velocity boundary condition
      type(BCBase),intent(inout)      :: BCLvs,BCvof
      !! The interface boundary condition
      real(kind=dp),intent(in)        :: Time
      !! The time
      integer(it8b),intent(in)        :: itt
      !! The total number of iterations
      real(dp),intent(in)             :: dt
      !! The time step size
      integer(it4b)			              :: i,j,k,kk
      real(dp)				                :: dtv
      real(kind=dp),dimension(:,:,:),pointer :: vfl,vflF    ! vfl represents the liquid volume fraction, vflF represents fluid volume fraction
      real(kind=dp),dimension(:,:,:),pointer :: phi,phiF    ! phi represents the liquid level set function, phiF represents fluid level set function
      real(kind=dp),dimension(:,:,:),pointer :: nxF,nyF,nzF ! nxF,nyF,nzF is fluid normal vector
      real(dp),dimension(:,:,:),allocatable  :: ue,ve,we,nx,ny,nz,dis
      real(dp),dimension(:,:,:),allocatable  :: temvfx,temvfy,temvfz ! The sub-step volume of fluid in x,y,z direction
      real(dp),dimension(:,:,:),allocatable  :: temlsx,temlsy,temlsz ! The sub-step level-set function in x,y,z direction
       
      ! The normal vector of interface 
      allocate(nx(0:imax+1,0:jmax+1,0:kmax+1))
      allocate(ny(0:imax+1,0:jmax+1,0:kmax+1))
      allocate(nz(0:imax+1,0:jmax+1,0:kmax+1))
      ! The distance from the cell centre to the interface
      allocate(dis(0:imax+1,0:jmax+1,0:kmax+1))
      ! The velocity at the cell face
      allocate(ue(0:imax+1,0:jmax+1,0:kmax+1))
      allocate(ve(0:imax+1,0:jmax+1,0:kmax+1))
      allocate(we(0:imax+1,0:jmax+1,0:kmax+1))
      ! The temporary volume fraction
      allocate(temvfx(0:imax+1,0:jmax+1,0:kmax+1))
      allocate(temvfy(0:imax+1,0:jmax+1,0:kmax+1))
      allocate(temvfz(0:imax+1,0:jmax+1,0:kmax+1))
      !  The temporary level set function
      allocate(temlsx(0:imax+1,0:jmax+1,0:kmax+1))
      allocate(temlsy(0:imax+1,0:jmax+1,0:kmax+1))
      allocate(temlsz(0:imax+1,0:jmax+1,0:kmax+1))

      if(associated(vfl).eqv..false.) then
        vfl(0:,0:,0:)=>PCell%vofL
      else
        nullify(vfl)
        vfl(0:,0:,0:)=>PCell%vofL
      end if  

      if(associated(phi).eqv..false.) then
        phi(0:,0:,0:)=>PCell%phiL
      else
        nullify(phi)
        phi(0:,0:,0:)=>PCell%phiL
      end if 

      if(associated(vflF).eqv..false.) then
        vflF(0:,0:,0:)=>PCell%Vof
      else
        nullify(vflF)
        vflF(0:,0:,0:)=>PCell%Vof
      end if
      if(associated(phiF).eqv..false.) then
        phiF(0:,0:,0:)=>PCell%phi
      else
        nullify(phiF)
        phiF(0:,0:,0:)=>PCell%phi
      end if

      if(associated(nxF).eqv..false.) then
        nxF(0:,0:,0:)=>PCell%nx
      else
        nullify(nxF)
        nxF(0:,0:,0:)=>PCell%nx
      end if  

      if(associated(nyF).eqv..false.) then
        nyF(0:,0:,0:)=>PCell%ny
      else
        nullify(nyF)
        nyF(0:,0:,0:)=>PCell%ny
      end if
        
      if(associated(nzF).eqv..false.) then
        nzF(0:,0:,0:)=>PCell%nz
      else
        nullify(nzF)
        nzF(0:,0:,0:)=>PCell%nz
      end if  
      ! Compute the sub time step size for clsvof
      dtv=dt/dble(nv)
      ! Set the cell face velocity to the state velocities 
      do i = 0,imax+1
        do j = 0,jmax+1
          do k = 0,kmax+1
            ue(i,j,k) = TVar%u(i,j,k)
            ve(i,j,k) = Tvar%v(i,j,k)
            we(i,j,k) = TVar%w(i,j,k)
          end do
        end do
      end do
    !  print*, ue(32,32,32),ve(10,10,10),we(10,10,10)
      do kk=1,nv
        if(mod(kk,3)==0) then
          temvfx(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsx(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
          ! print*,'X_Sweep first 000000000000000000000000000000000000000000000'
          ! print*,vfl(24,32,32)
          ! Advance the interface in the x-direction
          call X_Sweep(PGrid,vfl,phi,temvfx,temlsx,ue,ve,we,BCu,BCVof,BCLvs,   &
                                                              nx,ny,nz,dis,dtv)
          
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                ! Compute the temporary liquid volume fraction and level set function
                temvfx(i,j,k)=temvfx(i,j,k)/(1.d0-dtv/PGrid%dx(i,j,k)*         &
                                     (ue(i,j,k)-ue(i-1,j,k)))
                temlsx(i,j,k)=temlsx(i,j,k)/(1.d0-dtv/PGrid%dx(i,j,k)*         &
                                     (ue(i,j,k)-ue(i-1,j,k)))
                ! Check for undercomputed or overcomputed volume fraction 
                if(temvfx(i,j,k)<=vofeps) temvfx(i,j,k) = 0.d0
                if(temvfx(i,j,k)>=(1.d0-vofeps)) temvfx(i,j,k) = 1.d0
                ! Check for nan volume fraction or level set function
                if(isnan(temvfx(i,j,k)).or.isnan(temlsx(i,j,k)).or.            &
                                           dabs(temlsx(i,j,k))>1.d10) then
                   print*, i,j,k
                   print*, ue(i,j,k)-ue(i-1,j,k)
                !  pause 'X_Sweep 277'
                end if  
                vfl(i,j,k) = temvfx(i,j,k)
                phi(i,j,k) = temlsx(i,j,k)
              end do
            end do
          end do
          ! print*,'X_Sweep first after X_Sweep 00000000000000000000000000000000'
          ! print*,vfl(24,32,32),temvfx(24,32,32)
          ! print*,vfl(23,32,32),temvfx(23,32,32)
          ! print*,vfl(25,32,32),temvfx(25,32,32)
          temvfy(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsy(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
          ! Advance the interface in the y-direction
          call Y_Sweep(PGrid,vfl,phi,temvfy,temlsy,ue,ve,we,BCv,BCVof,BCLvs,   &
                                                              nx,ny,nz,dis,dtv)
          
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                temvfy(i,j,k)=temvfy(i,j,k)+dtv/PGrid%dy(i,j,k)*temvfx(i,j,k)* &
                                                (ve(i,j,k)-ve(i,j-1,k))
                temlsy(i,j,k)=temlsy(i,j,k)+dtv/PGrid%dy(i,j,k)*temlsx(i,j,k)* &
                                                (ve(i,j,k)-ve(i,j-1,k))
                if(temvfy(i,j,k)<=vofeps) temvfy(i,j,k) = 0.d0
                if(temvfy(i,j,k)>=(1.d0-vofeps)) temvfy(i,j,k) = 1.d0
                if(isnan(temvfy(i,j,k)).or.isnan(temlsy(i,j,k)).or.            &
                                           dabs(temlsy(i,j,k))>1.d10) then
                  print*, i,j,k
                  print*, ve(i,j,k)-ve(i,j-1,k)
                  pause 'Y_Sweep 305'
                end if  
                vfl(i,j,k) = temvfy(i,j,k)
                phi(i,j,k) = temlsy(i,j,k)
              end do
            end do
          end do
          ! print*,'X_Sweep first after Y_Sweep 000000000000000000000000000000000'
          ! print*,vfl(24,32,32),temvfy(24,32,32)
          ! print*,vfl(23,32,32),temvfy(23,32,32)
          ! print*,vfl(25,32,32),temvfy(25,32,32)
          temvfz(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsz(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
          call Z_Sweep(PGrid,vfl,phi,temvfz,temlsz,ue,ve,we,BCw,BCVof,BCLvs,   &
                                                              nx,ny,nz,dis,dtv)
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                vfl(i,j,k) = temvfz(i,j,k)+dtv/PGrid%dz(i,j,k)*temvfx(i,j,k)*  &
                                                        (we(i,j,k)-we(i,j,k-1))
                phi(i,j,k) = temlsz(i,j,k)+dtv/PGrid%dz(i,j,k)*temlsx(i,j,k)*  &
                                                      (we(i,j,k)-we(i,j,k-1))
                if(isnan(vfl(i,j,k)).or.isnan(phi(i,j,k)).or.                  &
                                        dabs(phi(i,j,k))>1.d10) then
                  print*, i,j,k
                  print*, we(i,j,k)-we(i,j,k-1)
                  pause 'X_Sweep 329'
                end if  
                if(vfl(i,j,k)<=vofeps) vfl(i,j,k) = 0.d0
                if(vfl(i,j,k)>=(1.d0-vofeps)) vfl(i,j,k) = 1.d0  
              end do
            end do
          end do
         ! print*, 'X_Sweep first after Z_Sweep 000000000000000000000000000000000'
         ! print*, vfl(24,32,32), temvfx(24,32,32), temvfz(24,32,32)
         ! print*, vfl(23,32,32), vfl(25,32,32), we(23,32,32), we(24,32,32)
          call BoundaryConditionLvsVof(PGrid,PCell,TVar,BCLvs,BCVof,Time)
        elseif(mod(kk,3)==1) then
          temvfy(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsy(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
          call Y_Sweep(PGrid,vfl,phi,temvfy,temlsy,ue,ve,we,BCv,BCVof,BCLvs,   &
                                                              nx,ny,nz,dis,dtv)
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                temvfy(i,j,k) = temvfy(i,j,k)/(1.d0-dtv/PGrid%dy(i,j,k)*       &
                                                        (ve(i,j,k)-ve(i,j-1,k)))
                temlsy(i,j,k) = temlsy(i,j,k)/(1.d0-dtv/PGrid%dy(i,j,k)*       &
                                                        (ve(i,j,k)-ve(i,j-1,k)))
                if(temvfy(i,j,k)<=vofeps) temvfy(i,j,k) = 0.d0
                if(temvfy(i,j,k)>=(1.d0-vofeps)) temvfy(i,j,k) = 1.d0
                if(isnan(temvfy(i,j,k)).or.isnan(temlsy(i,j,k)).or.            &
                                           dabs(temlsy(i,j,k))>1.d10) then
                  print*, temvfy(i,j,k), temlsy(i,j,k)
                  print*,i,j,k
                  print*,ve(i,j,k)-ve(i,j-1,k)
                  print*,ve(i,j-1,k)
                  pause 'Y_Sweep 354'
                end if  
                vfl(i,j,k) = temvfy(i,j,k)
                phi(i,j,k) = temlsy(i,j,k)
              end do
            end do
          end do
           ! print*, 'Y first after Y_Sweep 111111111111111111111111111111111111'
           ! print*, vfl(24,32,32)
           ! print*, temvfy(24,32,32)
           ! print*, 'End test Y_Sweep'
          temvfz(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsz(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
          call Z_Sweep(PGrid,vfl,phi,temvfz,temlsz,ue,ve,we,BCw,BCVof,BCLvs,   &
                                                               nx,ny,nz,dis,dtv)
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                temvfz(i,j,k)=temvfz(i,j,k)+dtv/PGrid%dz(i,j,k)*temvfy(i,j,k)* &
                                                      (we(i,j,k)-we(i,j,k-1))
                temlsz(i,j,k)=temlsz(i,j,k)+dtv/PGrid%dz(i,j,k)*temlsy(i,j,k)* &
                                                      (we(i,j,k)-we(i,j,k-1))
                if(temvfz(i,j,k)<=vofeps) temvfz(i,j,k)=0.d0
                if(temvfz(i,j,k)>=(1.d0-vofeps)) temvfz(i,j,k)=1.d0
                if(isnan(temvfz(i,j,k)).or.isnan(temlsz(i,j,k)).or.            &
                                           dabs(temlsz(i,j,k))>1.d10) then
                  print*, temlsz(i,j,k)
                  print*,i,j,k
                  print*,we(i,j,k)-we(i,j,k-1)
                  pause 'Z_Sweep 381'
                end if  
                vfl(i,j,k) = temvfz(i,j,k)
                phi(i,j,k) = temlsz(i,j,k)
              end do
            end do
          end do
          ! print*, 'Y_Sweep first after Z_Sweep 11111111111111111111111111111111'
          ! print*, vfl(24,32,32)
          ! print*, temvfy(24,32,32), temvfz(24,32,32)
          ! print*, '==============================='
          ! print*, 'End test Z_Sweep'
          temvfx(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsx(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
          call X_Sweep(PGrid,vfl,phi,temvfx,temlsx,ue,ve,we,BCu,BCVof,BCLvs,   &
                                                               nx,ny,nz,dis,dtv)
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                vfl(i,j,k)=temvfx(i,j,k)+dtv/PGrid%dx(i,j,k)*temvfy(i,j,k)*    &
                                                   (ue(i,j,k)-ue(i-1,j,k))
                phi(i,j,k)=temlsx(i,j,k)+dtv/PGrid%dx(i,j,k)*temlsy(i,j,k)*    &
                                                   (ue(i,j,k)-ue(i-1,j,k))
                if(isnan(vfl(i,j,k)).or.isnan(phi(i,j,k)).or.                  &
                                                  dabs(phi(i,j,k))>1.d10) then
                  print*, vfl(i,j,k), temvfx(i,j,k), temvfy(i,j,k)
                  print*, phi(i,j,k), temlsx(i,j,k), temlsy(i,j,k)
                  print*,i,j,k
                  print*,ue(i,j,k)-ue(i-1,j,k)
                  pause 'X_Sweep 404'
                end if  
                if(vfl(i,j,k)<=vofeps) vfl(i,j,k)=0.d0
                if(vfl(i,j,k)>=(1.d0-vofeps)) vfl(i,j,k) = 1.d0
              end do
            end do
          end do

            ! print*, 'Y_Sweep first after X Sweep 1111111111111111111111111111111'
            ! print*, '1111'
            ! print*, temvfx(24,32,32), temvfy(24,32,32)
            ! print*, temvfz(24,32,32), vfl(24,32,32)
            ! print*, 'End test volume of fluid compuation'
          call BoundaryConditionLvsVof(PGrid,PCell,TVar,BCLvs,BCVof,Time)
            ! print*, 'After using boundary condition'
            ! print*, PCell%vofL(24,32,32), vfl(24,32,32)
            ! print*, 'oooooooooooooooooooooooooooo'
        else
          temvfz(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsz(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
            ! print*, 'Z_SWeep before doing the Z_Sweep 2222222222222222222222222'
            ! print*, vfl(24,32,32)
            ! print*, 'size of vfl'
            ! print*, size(vfl,2)
            ! print*, temvfz(24,32,32)
            ! print*, '==================='
          call Z_Sweep(PGrid,vfl,phi,temvfz,temlsz,ue,ve,we,BCw,BCVof,BCLvs,   &
                                                               nx,ny,nz,dis,dtv)
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                temvfz(i,j,k)=temvfz(i,j,k)/(1.d0-dtv/PGrid%dz(i,j,k)*         &
                                                   (we(i,j,k)-we(i,j,k-1)))
                temlsz(i,j,k) = temlsz(i,j,k)/(1.d0-dtv/PGrid%dz(i,j,k)*       &
                                                   (we(i,j,k)-we(i,j,k-1)))
                if(temvfz(i,j,k)<=vofeps) temvfz(i,j,k)=0.d0
                if(temvfz(i,j,k)>=(1.d0-vofeps)) temvfz(i,j,k)=1.d0
                if(isnan(temvfz(i,j,k)).or.isnan(temlsz(i,j,k)).or.            &
                                           dabs(temlsz(i,j,k))>1.d10) then
                  print*,we(i,j,k)-we(i,j,k-1)
                  pause 'Z_Sweep 428'  
                endif  
                vfl(i,j,k) = temvfz(i,j,k)
                phi(i,j,k) = temlsz(i,j,k)
              end do
            end do
          end do
            ! print*, 'Z_Sweep first after Z_Sweep 2222222222222222222222222222'
            ! print*, 'kk '
            ! print*, temvfz(24,32,32), vfl(24,32,32)
            ! print*, '++++++++++++++++++++++++++'
          temvfx(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsx(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
          call X_Sweep(PGrid,vfl,phi,temvfx,temlsx,ue,ve,we,BCu,BCVof,BCLvs,   &
                                                               nx,ny,nz,dis,dtv)
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                temvfx(i,j,k)=temvfx(i,j,k)+dtv/PGrid%dx(i,j,k)*temvfz(i,j,k)* &
                                                         (ue(i,j,k)-ue(i-1,j,k))
                temlsx(i,j,k)=temlsx(i,j,k)+dtv/PGrid%dx(i,j,k)*temlsz(i,j,k)* &
                                                         (ue(i,j,k)-ue(i-1,j,k))
                if(temvfx(i,j,k)<=vofeps) temvfx(i,j,k)=0.d0
                if(temvfx(i,j,k)>=(1.d0-vofeps)) temvfx(i,j,k) = 1.d0
                if(isnan(temvfx(i,j,k)).or.isnan(temlsx(i,j,k)).or.            &
                                           dabs(temlsx(i,j,k))>1.d10) then
                  print*, ue(i,j,k)-ue(i-1,j,k)
                  pause 'X_Sweep 454'
                end if  
                vfl(i,j,k) = temvfx(i,j,k)
                phi(i,j,k) = temlsx(i,j,k)
              end do
            end do
          end do
          temvfy(0:Imax+1,0:Jmax+1,0:Kmax+1) = vfl(0:Imax+1,0:Jmax+1,0:Kmax+1)
          temlsy(0:Imax+1,0:Jmax+1,0:Kmax+1) = phi(0:Imax+1,0:Jmax+1,0:Kmax+1)
            ! print*, 'Z first 000000000000000 Test volume for X_Sweep'
            ! print*, temvfx(24,32,32), vfl(24,32,32)
            ! print*, 'End test X_Sweep'
          call Y_Sweep(PGrid,vfl,phi,temvfy,temlsy,ue,ve,we,BCv,BCVof,BCLvs,   &
                                                               nx,ny,nz,dis,dtv)
          do i = 1,imax
            do j = 1,jmax
              do k = 1,kmax
                vfl(i,j,k)=temvfy(i,j,k)+dtv/PGrid%dy(i,j,k)*temvfz(i,j,k)*    &
                                                      (ve(i,j,k)-ve(i,j-1,k))
                phi(i,j,k)=temlsy(i,j,k)+dtv/PGrid%dy(i,j,k)*temlsz(i,j,k)*    &
                                                      (ve(i,j,k)-ve(i,j-1,k))
                if(isnan(vfl(i,j,k)).or.isnan(phi(i,j,k)).or.                  &
                                        dabs(phi(i,j,k))>1.d10) then
                  print*,temlsy(i,j,k),temlsz(i,j,k)
                  print*, ve(i,j,k)-ve(i,j-1,k)
                  pause 'Y_Sweep 476'
                end if                                      
                if(vfl(i,j,k)<=vofeps) vfl(i,j,k) = 0.d0
                if(vfl(i,j,k)>=(1.d0-vofeps)) vfl(i,j,k) = 1.d0
              end do
            end do
          end do
            ! print*, 'Z first 000000000000000 Test volume for X_Sweep'
            ! print*, temvfy(24,32,32), vfl(24,32,32)
            ! print*, 'End test y_Sweep'
          call BoundaryConditionLvsVof(PGrid,PCell,TVar,BCLvs,BCVof,Time)
        end if
        call Interface_Reconstruct(PGrid,vfl,phi,nx,ny,nz,dis)
        do i = 1,imax
          do j = 1,jmax
            do k = 1,kmax
              PCell%nxL(i,j,k) = nx(i,j,k)
              PCell%nyL(i,j,k) = ny(i,j,k)
              PCell%nzL(i,j,k) = nz(i,j,k)
              dis(i,j,k) = (0.5d0*PGrid%dx(i,j,k)*dabs(nx(i,j,k))+             &
                            0.5d0*PGrid%dy(i,j,k)*dabs(ny(i,j,k))+             &
                            0.5d0*PGrid%dz(i,j,k)*dabs(nz(i,j,k)))-dis(i,j,k) 
            end do
          end do
        end do
        call Redistance_Levelset(PGrid,vfl,phi,nx,ny,nz,dis)
      end do
        ! print*, 'Test volume of fluid'
        ! print*, vfl(24,32,32)
        ! print*, 'End test volume of fluid compuation'
      if(associated(vfl).eqv..true.) nullify(vfl)
      if(associated(phi).eqv..true.) nullify(phi)
      if(associated(vflF).eqv..true.) nullify(vflF)
      if(associated(phiF).eqv..true.) nullify(phiF)
      if(associated(nxF).eqv..true.) nullify(nxF)
      if(associated(nyF).eqv..true.) nullify(nyF)
      if(associated(nzF).eqv..true.) nullify(nzF)

      deallocate(nx)
      deallocate(ny)
      deallocate(nz)
      deallocate(dis)
      deallocate(ue)
      deallocate(ve)
      deallocate(we)
      deallocate(temvfx)
      deallocate(temvfy)
      deallocate(temvfz)
      deallocate(temlsx)
      deallocate(temlsy)
      deallocate(temlsz)
    end subroutine Clsvof_Scheme

    ! build-up interface
    subroutine Isinterface(vfl,i,j,k,flag)
      implicit none
      real(kind=dp),dimension(0:,0:,0:), intent(in)  :: vfl
      integer(kind=it4b),                intent(in)  :: i,j,k
      logical,                           intent(out) :: flag
      flag = .false.
      if(vfl(i,j,k)>vofeps.and.vfl(i,j,k)<(1.d0-vofeps))      flag = .true.
      if(dabs(vfl(i-1,j,k)-vfl(i,j,k))>=(1.0d0-2.1d0*vofeps)) flag = .true.
      if(dabs(vfl(i+1,j,k)-vfl(i,j,k))>=(1.0d0-2.1d0*vofeps)) flag = .true.
      if(dabs(vfl(i,j-1,k)-vfl(i,j,k))>=(1.0d0-2.1d0*vofeps)) flag = .true.
      if(dabs(vfl(i,j+1,k)-vfl(i,j,k))>=(1.0d0-2.1d0*vofeps)) flag = .true.
      if(dabs(vfl(i,j,k-1)-vfl(i,j,k))>=(1.0d0-2.1d0*vofeps)) flag = .true.
      if(dabs(vfl(i,j,k+1)-vfl(i,j,k))>=(1.0d0-2.1d0*vofeps)) flag = .true.
      return
    end subroutine Isinterface

    subroutine X_Sweep(PGrid,vfl,phi,temvf,temls,ue,ve,we,BCu,BCVof,BCLvs,     &
                                                      nxx,nyy,nzz,diss,dtv)
       !! The subroutine is use to compute the liquid volume fraction in the 
       !! x-direction after advancing the interface.  
       implicit none
       TYPE(Grid),                           intent(in)    :: PGrid
       real(dp),dimension(0:,0:,0:),         intent(in)    :: vfl,phi
       real(dp),                             intent(in)    :: dtv
       real(dp),dimension(:,:,:),allocatable,intent(in)    :: ue,ve,we
       type(BCBase),intent(in)                             :: BCu,BCVof,BCLvs
       real(dp),dimension(:,:,:),allocatable,intent(inout) :: nxx,nyy,nzz
       real(dp),dimension(:,:,:),allocatable,intent(inout) :: diss
       real(dp),dimension(:,:,:),allocatable,intent(inout) :: temvf,temls
       integer                                             :: i,j,k
       real(dp)                                            :: flux,lse
       call Interface_Reconstruct(PGrid,vfl,phi,nxx,nyy,nzz,diss)
       flux = 0.d0
       ! volume of fluid
       do j = 1,jmax
         do k = 1,kmax
           do i = 1,imax-1
             if(ue(i,j,k)>0.d0) then
               if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps) then
                 flux = vfl(i,j,k)*ue(i,j,k)*dtv/PGrid%dx(i,j,k)
               else
                 call East_Flux(nxx(i,j,k),nyy(i,j,k),nzz(i,j,k),              &
                         diss(i,j,k),vfl(i,j,k),PGrid%dx(i,j,k),               &
                         PGrid%dy(i,j,k),PGrid%dz(i,j,k),ue(i,j,k)*dtv,flux)
               end if
             else
               if(vfl(i+1,j,k)>=(1.d0-vofeps).or.vfl(i+1,j,k)<=vofeps) then
                 flux=vfl(i+1,j,k)*ue(i,j,k)*dtv/PGrid%dx(i,j,k)
               else
                 call West_Flux(nxx(i+1,j,k),nyy(i+1,j,k),nzz(i+1,j,k),        &
                       diss(i+1,j,k),vfl(i+1,j,k),PGrid%dx(i+1,j,k),           &
                       PGrid%dy(i+1,j,k),PGrid%dz(i+1,j,k),-ue(i,j,k)*dtv,flux)
                 flux = -flux
               end if
             end if
         
             temvf(i,j,k) = temvf(i,j,k)-flux
             if(i<imax) temvf(i+1,j,k) = temvf(i+1,j,k)+flux
             ! if((i==24.or.i+1==24).and.j==32.and.k==32) then
             !    print*, 'Inside do loop test the flux X_Sweep xxxxxxx'
             !    print*, flux
             !    print*, i,j,k
             !    print*, nxx(i,j,k), nyy(i,j,k), nzz(i,j,k)
             !    print*, '========'
             !    print*, vfl(i,j,k)
             !    print*, ve(i,j,k), dtv, Pgrid%dy(i,j,k)
             !    print*,'End test X_Sweep xxxxxx'
             !  end if  
           end do
           ! for i=1
           ! flux=BCVof%VarW(j,k)*BCu%VarW(j,k)*dtv/PGrid%dx(1,j,k)
           flux=vfl(1,j,k)*ue(0,j,k)*dtv/PGrid%dx(1,j,k)
           temvf(1,j,k)=temvf(1,j,k)+flux
           ! for i=imax
           ! flux=BCVof%VarE(j,k)*BCu%VarE(j,k)*dtv/PGrid%dx(imax,j,k)
           flux=vfl(imax,j,k)*ue(imax,j,k)*dtv/PGrid%dx(imax,j,k)
           temvf(imax,j,k)=temvf(imax,j,k)-flux
         end do
       end do
       ! level set
       lse = 0.d0
       flux = 0.d0
       do j = 1,jmax
         do k = 1,kmax
           do i = 2,imax-1
             if(ue(i,j,k)>=0.d0) then
               lse=phi(i,j,k)+PGrid%dx(i,j,k)/2.d0*(1.d0-ue(i,j,k)*            &
                   dtv/PGrid%dx(i,j,k))*(phi(i+1,j,k)-phi(i-1,j,k))/           &
                   (2.d0*PGrid%dx(i,j,k))
             else
               if(i<=imax-2) then
                 lse=phi(i+1,j,k)-PGrid%dx(i,j,k)/2.d0*(1.d0+ue(i,j,k)*        &
                     dtv/PGrid%dx(i,j,k))*(phi(i+2,j,k)-phi(i,j,k))/           &
                     (2.d0*PGrid%dx(i,j,k))
               else
                 lse=phi(i+1,j,k)-PGrid%dx(i,j,k)/2.d0*(1.d0+ue(i,j,k)*        &
                     dtv/PGrid%dx(i,j,k))*(phi(i+1,j,k)-phi(i,j,k))/           &
                     (PGrid%dx(i,j,k))
               end if
             end if
             flux = ue(i,j,k)*lse*dtv/PGrid%dx(i,j,k)
             if(i>1) temls(i,j,k) = temls(i,j,k)-flux
             if(i<imax) temls(i+1,j,k) = temls(i+1,j,k)+flux
             if(isnan(temls(i,j,k)).or.isnan(temls(i+1,j,k))) then
               print*, PGrid%dx(i,j,k)/2.d0*(1.d0-ue(i,j,k)*            &
                   dtv/PGrid%dx(i,j,k))*(phi(i+1,j,k)-phi(i-1,j,k))/           &
                   (2.d0*PGrid%dx(i,j,k))
               print*, PGrid%dx(i,j,k)/2.d0*(1.d0+ue(i,j,k)*        &
                     dtv/PGrid%dx(i,j,k))*(phi(i+2,j,k)-phi(i,j,k))/           &
                     (2.d0*PGrid%dx(i,j,k))    
               print*, flux
               print*, lse
               print*, phi(i,j,k)
               print*, i,j,k
               pause 'Flux inside the X_Sweep with flux'
             end if 
           end do
           ! Reduce to the first order
           if(ue(1,j,k)>=0.d0) then
             lse=phi(1,j,k)+PGrid%dx(1,j,k)/2.d0*(1.d0-ue(1,j,k)*dtv/          &
                            PGrid%dx(1,j,k))*(phi(2,j,k)-phi(1,j,k))/PGrid%dx(1,j,k)
           else
             lse=phi(2,j,k)-PGrid%dx(2,j,k)/2.d0*(1.d0+ue(1,j,k)*dtv/          &
                            PGrid%dx(2,j,k))*(phi(3,j,k)-phi(1,j,k))/2.d0*PGrid%dx(2,j,k)
           end if
           flux=lse*ue(1,j,k)*dtv
           temls(1,j,k)=temls(1,j,k)-flux/PGrid%dx(1,j,k)
           temls(2,j,k)=temls(2,j,k)+flux/PGRid%dx(2,j,k)
           ! For boundary value
           flux=BCu%VarW(j,k)*BCLvs%VarW(j,k)*dtv/PGrid%dx(1,j,k)
           temls(1,j,k)=temls(1,j,k)+flux
           ! For i=kmax
           if(ue(Imax,j,k)>0.d0) then
             lse=phi(Imax,j,k)+PGrid%dx(Imax,j,k)/2.d0*                        &
                 (1.d0-ue(Imax,j,k)*dtv/PGrid%dx(Imax,j,k))*                   &
                 (phi(Imax,j,k)-phi(Imax-1,j,k))/                              &
                 (PGrid%x(Imax,j,k)-PGrid%x(Imax-1,j,k))
           else
             lse=phi(Imax,j,k)-PGrid%dx(Imax-1,j,k)/2.d0*                      &
                 (1.d0+ue(Imax,j,k)*dtv/PGrid%dx(Imax-1,j,k))*                 &
                 (phi(Imax,j,k)-phi(Imax-1,j,k))/                              &
                 (PGrid%x(Imax,j,k)-PGrid%x(Imax-1,j,k))
           end if
           flux=lse*ue(Imax,j,k)*dtv/PGrid%dx(Imax,j,k)
           temls(Imax,j,k)=temls(Imax,j,k)-flux
         end do
       end do
    end subroutine X_Sweep

    subroutine Y_Sweep(PGrid,vfl,phi,temvf,temls,ue,ve,we,BCv,BCVof,BCLvs,     &
                                                        nxx,nyy,nzz,diss,dtv)
       implicit none
       TYPE(Grid),                           intent(in)    :: PGrid
       real(dp),dimension(0:,0:,0:),         intent(in)    :: vfl,phi
       real(dp),dimension(:,:,:),allocatable,intent(in)    :: ue,ve,we
       type(BCBase),                         intent(in)    :: BCv,BCVof,BCLvs
       real(dp),                             intent(in)    :: dtv
       real(dp),dimension(:,:,:),allocatable,intent(inout) :: nxx,nyy,nzz,diss
       real(dp),dimension(:,:,:),allocatable,intent(inout) :: temvf,temls
       integer                                             :: i,j,k
       real(dp)                                            :: flux,lsn
       call Interface_Reconstruct(PGrid,vfl,phi,nxx,nyy,nzz,diss)
       flux = 0.d0
       ! volume of fluid
       do i = 1,imax
         do k = 1,kmax
           do j = 1,jmax-1
             if(ve(i,j,k)>0.d0) then
               ! The condition for applying flux computation
               if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps) then
                 flux = vfl(i,j,k)*ve(i,j,k)*dtv/PGrid%dy(i,j,k)
               else
                 ! For interface cell
                 call North_Flux(nxx(i,j,k),nyy(i,j,k),nzz(i,j,k),diss(i,j,k), &
                                vfl(i,j,k),PGrid%dx(i,j,k),PGrid%dy(i,j,k),    &
                                PGrid%dz(i,j,k),ve(i,j,k)*dtv,flux)
               end if
             else
               if(vfl(i,j+1,k)>=(1.d0-vofeps).or.vfl(i,j+1,k)<=vofeps) then
                 flux = vfl(i,j+1,k)*ve(i,j,k)*dtv/PGrid%dy(i,j,k)
               else
                 call South_Flux(nxx(i,j+1,k),nyy(i,j+1,k),nzz(i,j+1,k),       &
                               diss(i,j+1,k),vfl(i,j+1,k),PGrid%dx(i,j+1,k),   &
                               PGrid%dy(i,j+1,k),PGrid%dz(i,j+1,k),-ve(i,j,k)*dtv,flux)
                 flux = -flux
               end if
             end if
             temvf(i,j,k) = temvf(i,j,k)-flux
             if(j<jmax) temvf(i,j+1,k) = temvf(i,j+1,k)+flux
              ! if((i==24.or.i+1==24).and.j==32.and.k==32) then
              !   print*, 'Inside do loop test the flux'
              !   print*, flux
              !   print*, i,j,k
              !   print*, '========'
              !   print*, vfl(i,j,k), vfl(32,32,10)
              !   print*, ve(i,j,k), dtv, Pgrid%dy(i,j,k)
              !   print*, temvf(1,1,2)
              ! end if  
           end do
           ! for j = 1
           ! flux=BCVof%VarS(i,k)*BCv%VarS(i,k)*dtv/PGrid%dy(i,1,k)
           flux=vfl(i,1,k)*ve(i,0,k)*dtv/PGrid%dy(i,1,k)
           temvf(i,1,k)=temvf(i,1,k)+flux
           ! for j = jmax
           ! flux=BCVof%VarN(i,k)*BCv%VarN(i,k)*dtv/PGrid%dy(i,jmax,k)
           flux=vfl(i,jmax,k)*ve(i,jmax,k)*dtv/PGrid%dy(i,jmax,k)
           temvf(i,jmax,k)=temvf(i,jmax,k)-flux 
         end do
       end do
         ! print*, 'Inside Y_Sweep Clsvof 895'
         ! print*, flux
         ! print*, temvf(1,1,2)
         ! print*, 'End test'
       lsn = 0.d0
       flux = 0.d0
    ! level set
       do i = 1,imax
         do k = 1,kmax
           do j = 2,jmax-1 
             if(ve(i,j,k)>=0.d0) then
               lsn=phi(i,j,k)+PGrid%dy(i,j,k)/2.d0*(1.d0-ve(i,j,k)*dtv/        &
                   PGrid%dy(i,j,k))*(phi(i,j+1,k)-phi(i,j-1,k))/               &
                   (2.d0*PGrid%dy(i,j,k))
             else
               if(j<=jmax-2) then
                 lsn=phi(i,j+1,k)-PGrid%dy(i,j,k)/2.d0*(1.d0+ve(i,j,k)*dtv/    &
                     PGrid%dy(i,j,k))*(phi(i,j+2,k)-phi(i,j,k))/               &
                     (2.d0*PGrid%dy(i,j,k))
               else
                 lsn=phi(i,j+1,k)-PGrid%dy(i,j,k)/2.d0*(1.d0+ve(i,j,k)*dtv/    &
                     PGrid%dy(i,j,k))*(phi(i,j+1,k)-phi(i,j,k))/PGrid%dy(i,j,k)
               end if
             end if
             flux=lsn*ve(i,j,k)*dtv/PGrid%dy(i,j,k)
             if(j>=2) temls(i,j,k) = temls(i,j,k)-flux
             if(j<jmax) temls(i,j+1,k) = temls(i,j+1,k)+flux
           end do
           if(ve(i,1,k)>=0.d0) then
             lsn=phi(i,1,k)+PGrid%dy(i,1,k)/2.d0*(1.d0-ve(i,1,k)*              &
                 dtv/PGrid%dy(i,1,k))*(phi(i,2,k)-phi(i,1,k))/                 &
                (PGrid%y(i,2,k)-PGrid%y(i,1,k))
           else
             lsn=phi(i,2,k)-PGrid%dy(i,2,k)/2.d0*(1.d0+ve(i,1,k)*              &
                 dtv/PGrid%dy(i,2,k))*(phi(i,3,k)-phi(i,1,k))/                 &
                (PGrid%y(i,3,k)-PGrid%y(i,1,k))
           end if
           flux=lsn*ve(i,1,k)*dtv
           temls(i,1,k)=temls(i,1,k)-flux/PGrid%dy(i,1,k)
           temls(i,2,k)=temls(i,2,k)+flux/PGrid%dy(i,2,k)
           ! for i=1
           flux=BCLvs%VarS(i,k)*BCv%VarS(i,k)*dtv/PGrid%dy(i,1,k)
           temls(i,1,k)=temls(i,1,k)+flux
           if(ve(i,jmax,k)>0.d0) then
             lsn=phi(i,jmax,k)+PGrid%dy(i,jmax,k)/2.d0*(1.d0-ve(i,jmax,k)*dtv/ &
                 PGrid%dy(i,jmax,k))*(phi(i,jmax,k)-phi(i,jmax-1,k))/          &
                (PGrid%y(i,jmax,k)-PGrid%y(i,jmax-1,k))
           else
             lsn=phi(i,jmax,k)-PGrid%dy(i,jmax-1,k)/2.d0*(1.d0+ve(i,jmax,k)*dtv/&
                 PGrid%dy(i,jmax-1,k))*(phi(i,jmax,k)-phi(i,jmax-1,k))/        &
                 (PGrid%y(i,jmax,k)-PGrid%y(i,jmax-1,k))
           end if
           flux=lsn*ve(i,jmax,k)*dtv
           temls(i,jmax,k)=temls(i,jmax,k)-flux/PGrid%dy(i,jmax,k) 
         end do
       end do
    end subroutine Y_Sweep

    subroutine Z_Sweep(PGrid,vfl,phi,temvf,temls,ue,ve,we,BCw,BCVof,BCLvs,     &
                                                          nxx,nyy,nzz,diss,dtv)
       implicit none
       type(Grid),                           intent(in)    :: PGrid
       real(dp),dimension(0:,0:,0:),         intent(in)    :: vfl,phi
       real(dp),                             intent(in)    :: dtv
       real(dp),dimension(:,:,:),allocatable,intent(in)    :: ue,ve,we
       type(BCBase),                         intent(in)    :: BCw,BCVof,BCLvs
       real(dp),dimension(:,:,:),allocatable,intent(inout) :: nxx,nyy,nzz,diss
       real(dp),dimension(:,:,:),allocatable,intent(inout) :: temvf,temls
       integer                                             :: i,j,k
       real(dp)						                                 :: flux,lst
       call Interface_Reconstruct(PGrid,vfl,phi,nxx,nyy,nzz,diss)
       flux = 0.d0
       ! volume of fluid
       do i = 1,imax
         do j = 1,jmax
           do k = 1,kmax-1
             if(we(i,j,k)>=0.d0) then
               if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps)then
                 flux = vfl(i,j,k)*we(i,j,k)*dtv/PGrid%dz(i,j,k)
               else
                 call Top_Flux(nxx(i,j,k),nyy(i,j,k),nzz(i,j,k),diss(i,j,k),   &
                               vfl(i,j,k),PGrid%dx(i,j,k),PGrid%dy(i,j,k),     &
                               PGrid%dz(i,j,k),we(i,j,k)*dtv,flux)
               end if
             else
               if(vfl(i,j,k+1)>=(1.d0-vofeps).or.vfl(i,j,k+1)<=vofeps) then
                 flux = vfl(i,j,k+1)*we(i,j,k)*dtv/PGrid%dz(i,j,k)
               else
                 call Bottom_Flux(nxx(i,j,k+1),nyy(i,j,k+1),nzz(i,j,k+1),      &
                            diss(i,j,k+1),vfl(i,j,k+1),PGrid%dx(i,j,k+1),      &
                            PGrid%dy(i,j,k+1),PGrid%dz(i,j,k+1),-we(i,j,k)*dtv,flux)
                 flux = -flux
               end if
             end if
             temvf(i,j,k) = temvf(i,j,k)-flux
             if(k<kmax) temvf(i,j,k+1) = temvf(i,j,k+1)+flux
             if(isnan(flux)) then 
               print*, i,j,k
               print*, flux, vfl(i,j,k)
               pause 'Vof-scheme 361'
             end if  
             ! if(i==32.and.j==32.and.(k==28.or.k+1==28)) then
             !   print*, 'Inside the Z_Sweep'
             !   print*, vfl(i,j,k), we(i,j,k)
             !   print*, 'normal vector'
             !   print*, nxx(i,j,k), nyy(i,j,k), nzz(i,j,k)
             !   print*, flux, temvf(24,32,32)
             !   print*, 'End test inside Z_Sweep====='
             ! end if  
           end do
           ! for k=1 at bottom boundary
           ! flux = BCVof%VarB(i,j)*BCw%VarB(i,j)*dtv/PGrid%dz(i,j,1)
           flux=vfl(i,j,1)*we(i,j,0)*dtv/PGrid%dz(i,j,1)
           temvf(i,j,1)=temvf(i,j,1)+flux
           ! fro k=kmax at top boundary
           ! flux = BCVof%VarT(i,j)*BCw%VarT(i,j)*dtv/PGrid%dz(i,j,kmax)
           flux=vfl(i,j,kmax)*we(i,j,kmax)*dtv/PGrid%dz(i,j,kmax)
           temvf(i,j,kmax)=temvf(i,j,kmax)-flux
         end do
       end do
       ! level set
       do j = 1,jmax
         do i = 1,imax
           do k = 2,kmax-1
             if(we(i,j,k)>=0.d0) then
               lst=phi(i,j,k)+PGrid%dz(i,j,k)/2.d0*(1.d0-we(i,j,k)*dtv/        &
                              PGrid%dz(i,j,k))*(phi(i,j,k+1)-phi(i,j,k-1))/    &
                              (2.d0*PGrid%dz(i,j,k))
             else
               if(k<=kmax-2) then
                 lst=phi(i,j,k+1)-PGrid%dz(i,j,k)/2.d0*(1.d0+we(i,j,k)*dtv/    &
                     PGrid%dz(i,j,k))*(phi(i,j,k+2)-phi(i,j,k))/               &
                     (2.d0*PGrid%dz(i,j,k))
               else
                 lst=phi(i,j,k+1)-PGrid%dz(i,j,k)/2.d0*(1.d0+we(i,j,k)*dtv/    &
                     PGrid%dz(i,j,k))*(phi(i,j,k+1)-phi(i,j,k))/PGrid%dz(i,j,k)
               end if
             end if
             flux = lst*we(i,j,k)*dtv/PGrid%dz(i,j,k)
             if(k>=2) temls(i,j,k) = temls(i,j,k)-flux
             if(k<=kmax-1) temls(i,j,k+1) = temls(i,j,k+1)+flux
           end do
           ! For k=1 and k=2
           if(we(i,j,1)>=0.d0) then
             lst=phi(i,j,1)+PGrid%dz(i,j,1)/2.d0*(1.d0-we(i,j,1)*              &
                 dtv/PGrid%dz(i,j,1))*(phi(i,j,2)-phi(i,j,1))/                 &
                (PGrid%z(i,j,2)-PGrid%z(i,j,1))
           else
             lst=phi(i,j,2)+PGrid%dz(i,j,2)/2.d0*(1.d0+we(i,j,1)*              &
                 dtv/PGrid%dz(i,j,2))*(phi(i,j,3)-phi(i,j,1))/                 &
                (PGrid%z(i,j,3)-PGrid%z(i,j,1))
           end if
           flux=lst*we(i,j,1)*dtv 
           temls(i,j,1)=temls(i,j,1)-flux/PGrid%dz(i,j,1)
           temls(i,j,2)=temls(i,j,2)+flux/PGrid%dz(i,j,2)
           ! For boundary cell
           flux=BCw%VarB(i,j)*BCLvs%VarB(i,j)*dtv/PGrid%dz(i,j,1)
           temls(i,j,1)=temls(i,j,1)+flux
           ! For boundary cell k=kmax
           if(we(i,j,kmax)>0.d0) then
             lst=phi(i,j,kmax)+PGrid%dz(i,j,kmax)/2.d0*(1.d0-we(i,j,kmax)*dtv/ &
                 PGrid%dz(i,j,kmax))*(phi(i,j,kmax)-phi(i,j,kmax-1))/          &
                (PGrid%z(i,j,kmax)-PGrid%z(i,j,kmax-1))
           else
             lst=phi(i,j,kmax)-PGrid%dz(i,j,kmax-1)/2.d0*                      &
                 (1.d0+we(i,j,kmax)*dtv/PGrid%dy(i,j,kmax-1))*                 &
                 (phi(i,j,kmax)-phi(i,j,kmax-1))/                              &
                 (PGrid%z(i,j,kmax)-PGrid%z(i,j,kmax-1))
           end if
           flux=lst*we(i,j,kmax)*dtv
           temls(i,j,kmax)=temls(i,j,kmax)-flux/PGrid%dz(i,j,kmax)
           if(dabs(temls(i,j,kmax))>1.d10.or.isnan(temls(i,j,kmax))) then
             print*, 'Too large level set function Z_Sweep'
             print*, lst
             print*, '====='
             print*, PGrid%dz(i,j,kmax)/2.d0*(1.d0-we(i,j,kmax)*dtv/           &
                 PGrid%dz(i,j,kmax))*(phi(i,j,kmax)-phi(i,j,kmax-1))/          &
                (PGrid%z(i,j,kmax)-PGrid%z(i,j,kmax-1))
             print*, PGrid%dz(i,j,kmax-1)/2.d0*                                &
                 (1.d0+we(i,j,kmax)*dtv/PGrid%dy(i,j,kmax-1))*                 &
                 (phi(i,j,kmax)-phi(i,j,kmax-1))/                              &
                 (PGrid%z(i,j,kmax)-PGrid%z(i,j,kmax-1))
             print*, flux, phi(i,j,k), phi(i,j,k+1)
             pause 'Z_Sweep' 
           end if 
         end do
       end do
    end subroutine Z_Sweep

    subroutine Interface_Reconstruct(PGrid,vfl,phi,nxx,nyy,nzz,diss)
      !! The subroutine is used to compute the interface normal vector and 
      !! the distance from the cell centre to the interface.
      implicit none
      type(Grid),                                 intent(in)   :: PGrid
      !! The pressure grid
      real(kind=dp),dimension(0:,0:,0:),          intent(in)   :: vfl
      !! The volume of fluid
      real(kind=dp),dimension(0:,0:,0:),          intent(in)   :: phi
      !! The level-set function
      real(kind=dp),dimension(:,:,:),allocatable, intent(inout) :: nxx,nyy,nzz
      !! The normal vector
      real(kind=dp),dimension(:,:,:),allocatable, intent(inout) :: diss
      !! The distance from cell centre to the interface
      integer                                    :: i,j,k
      real(kind=dp)						                   :: nxx1,nyy1,nzz1,diss1
      real(kind=dp)                              :: temp
      logical                                    :: flag
      nxx(:,:,:)=0.d0
      nyy(:,:,:)=0.d0
      nzz(:,:,:)=0.d0
      do i = 1,imax
        do j = 1,jmax
          do k = 1,kmax
            ! Check whether this current cell contains the interface or not
            call Isinterface(vfl,i,j,k,flag)
            if(flag.eqv..true.) then
              if(i>1.and.i<Imax.and.j>1.and.j<jmax.and.k>1.and.k<kmax) then
                ! Compute the normal vector using the scheme in my thesis 
                ! "A Cartesian Cut-Cell Methodology, applied to large scale interface 
                ! dynamics and wave impacts spray formation"
                call Normal_Vector_Irre(PGrid,vfl,phi,i,j,k,nxx1,nyy1,nzz1)
              else
                ! For boundary cell
                nxx1=(phi(min(imax,i+1),j,k)-phi(max(1,i-1),j,k))/             &
                     (PGrid%x(min(imax,i+1),j,k)-PGrid%x(max(1,i-1),j,k))  
                nyy1=(phi(i,min(jmax,j+1),k)-phi(i,max(1,j-1),k))/             &
                     (PGrid%y(i,min(jmax,j+1),k)-PGrid%y(i,max(1,j-1),k))
                nzz1=(phi(i,j,min(kmax,k+1))-phi(i,j,max(1,k-1)))/             &
                     (PGrid%z(i,j,min(kmax,k+1))-PGrid%z(i,j,max(1,k-1)))
              end if  
              ! Recompute the normal vector such that summation of the square of all its component is 1
              temp = dsqrt(nxx1**2.d0+nyy1**2.d0+nzz1**2.d0)
              if(isnan(temp)) then
                print*, i,j,k
                print*, nxx1,nyy1,nzz1
                print*, phi(i,j,k+1),phi(i,j,k),phi(i,j,k-1)
                pause 'vof 903'
              end if
              if(temp<1.d-14) then
                nxx1 = 0.d0
                nyy1 = 0.d0
                nzz1 = 1.d0
              else
                nxx1 = nxx1/temp
                nyy1 = nyy1/temp
                nzz1 = nzz1/temp
              end if
              ! Compute the distance from cell centre to the interface
              call Find_Distance(PGrid%dx(i,j,k),PGrid%dy(i,j,k),             &
                   PGrid%dz(i,j,k),nxx1,nyy1,nzz1,vfl(i,j,k),diss1)
              if(isnan(diss1)) then
                pause 'reconstruct 245'
              end if
              if(j==32.and.k==32.and.vfl(i,j,k)>0.d0.and.i>20) then
               print*, 'Test is interface'
               print*, i,j,k
               print*, vfl(i,j,k)
               print*, nxx1,nyy1,nzz1,temp
               print*, flag
             end if
            else
            ! Set the normal vector for cells containing only liquid or gas  
             if(j==32.and.k==32.and.vfl(i,j,k)>0.d0.and.i>20) then
               print*, 'Test is interface'
               print*, i,j,k
               print*, vfl(i,j,k)
               print*, flag
             end if  
              nxx1 = 0.d0
              nyy1 = 0.d0
              nzz1 = 0.d0
              diss1 = 0.d0
            ! Set the distance for these cells  
            ! diss1 = (0.5d0-vfl(i,j,k))*4.d0*PGrid%dz(i,j,k)
            end if
            nxx(i,j,k)=nxx1
            nyy(i,j,k)=nyy1
            nzz(i,j,k)=nzz1
            diss(i,j,k)=diss1
          end do
        end do
      end do
    end subroutine Interface_Reconstruct

    subroutine Find_Distance(dx,dy,dz,nxx,nyy,nzz,f,s)
       real(dp),intent(in):: dx,dy,dz,nxx,nyy,nzz,f
       real(dp),intent(out):: s
       real(dp):: nx1,ny1,nz1,dx1,dy1,dz1,case_choose
       real(dp):: sc,sm,fc,v1,v2,v3,v31,v32
       integer:: case2
       nx1 = dabs(nxx)
       ny1 = dabs(nyy)
       nz1 = dabs(nzz)
       dx1 = dmax1(nx1*dx,ny1*dy,nz1*dz)
       dz1 = dmin1(nx1*dx,ny1*dy,nz1*dz)
       if(dabs(dz1)<tolpar) dz1=0.d0
       dy1 = (nx1*dz+ny1*dy+nz1*dx)-dx1-dz1
       if(dabs(dy1)<tolpar) dy1=0.d0
       fc = dmin1(f,1.d0-f)
       sm = dx1+dy1+dz1
       sc = (6.d0*fc*dx1*dy1*dz1)**(1.d0/3.d0)
       if(isnan(sc)) then
         print*, 'Something wrong with sc Find_Distance clsvof 1183'
         print*, fc, dx1, dy1, dz1
         print*, nx1, ny1, nz1
         pause 'End test'
       end if  
       if(sc<dz1) then
          call Final_Distance(f,sc,sm,s)
       else
          sc = 0.5d0*dz1+dsqrt(2.d0*fc*dx1*dy1-dz1**2.d0/12.d0)
          if(sc<dy1) then
             call Final_Distance(f,sc,sm,s)
             if(isnan(s)) then
                pause
             end if
          else
             if(dx1>=(dy1+dz1)) then
                sc = fc*dx1+(dy1+dz1)/2.d0
                if(sc>=(dy1+dz1)) then
                   call Final_Distance(f,sc,sm,s)
                else
                   call Cubic_Equation_Solve(fc,dx1,dy1,dz1,sc,case2)
                   if(sc>=(dy1+dz1)) then
                      pause 'final_distance 546'
                   end if
                   call Final_Distance(f,sc,sm,s)
                end if
             else
                sc = fc*dx1+(dy1+dz1)/2.d0
                if(sc>=dx1) then
                   call Cubic_Equation_Solve(fc,dx1,dy1,dz1,sc)
                   if(sc<dx1) then
                   ! pause 'final_distance 555'
                      call Cubic_Equation_Solve(fc,dx1,dy1,dz1,sc,case2)
                      if(sc>dx1) then
                         pause 'final_distance 558'
                       end if
                   end if
                   call Final_Distance(f,sc,sm,s)
                else
                   call Cubic_Equation_Solve(fc,dx1,dy1,dz1,sc,case2)
                   if(sc>dx1) then
                      pause 'final_distance 561'
                   end if
                   call Final_Distance(f,sc,sm,s)
                end if
             end if
          end if
       end if
    end subroutine Find_Distance

    subroutine Cubic_Equation_Solve(fc,dx1,dy1,dz1,sc,case2)
       real(dp),intent(in):: fc,dx1,dy1,dz1
       real(dp),intent(inout):: sc
       integer,optional:: case2
       real(dp):: a1,a2,a3,a0,p0,q0,delta,teta
       if(present(case2)) then
          a3 = 1.d0
          a2 = -3.d0*(dy1+dz1)
          a1 = 3.d0*(dy1**2.d0+dz1**2.d0)
          a0 = 6.d0*fc*dx1*dy1*dz1-(dy1**3.d0+dz1**3.d0)
          call Newton_Raphson(a3,a2,a1,a0,sc)
       else
          a3 = 2.d0
          a2 = -3.d0*(dx1+dy1+dz1)
          a1 = 3.d0*(dx1**2.d0+dy1**2.d0+dz1**2.d0)
          a0 = 6.d0*fc*dx1*dy1*dz1-(dx1**3.d0+dy1**3.d0+dz1**3.d0)
          call Newton_Raphson(a3,a2,a1,a0,sc)
       end if
    end subroutine Cubic_Equation_Solve

    subroutine Final_Distance(f,sc,sm,s)
       real(dp),intent(in):: f,sc,sm
       real(dp),intent(out):: s
       if(f<=0.5d0) then
          s = sc
       else
          s = sm-sc
       end if
    end subroutine Final_Distance

    subroutine Volume_Calc(dx,dy,dz,nxx,nyy,nzz,alpha,vol)
       implicit none
       real(dp),intent(in):: dx,dy,dz,nxx,nyy,nzz,alpha
       real(dp),intent(out):: vol
       real(dp):: nx1,ny1,nz1,alphamax,f1,f2,f3,fm1,fm2,fm3
       real(dp):: vol1,vol2,slop_eps = 1.d-12
       nx1 = dabs(nxx)
       ny1 = dabs(nyy)
       nz1 = dabs(nzz)
       if(alpha<=0.d0) then
          vol = 0.d0
          return
       end if
       if(alpha-(nx1*dx+ny1*dy+nz1*dz)>=0.d0) then
          vol = 1.d0*dx*dy*dz
          return
       end if
       if(alpha-nx1*dx>0) then
          f1 = (alpha-nx1*dx)**3.d0
       else
          f1 = 0.d0
       end if
       if(alpha-ny1*dy>0) then
          f2 = (alpha-ny1*dy)**3.d0
       else
          f2 = 0.d0
       end if
       if(alpha-nz1*dz>0) then
          f3 = (alpha-nz1*dz)**3.d0
       else
          f3 = 0.d0
       end if
       alphamax = nx1*dx+ny1*dy+nz1*dz
       if(alpha-alphamax+nx1*dx>0) then
          fm1 = (alpha-alphamax+nx1*dx)**3.d0
       else
          fm1 = 0.d0
       end if
       if(alpha-alphamax+ny1*dy>0) then
          fm2 = (alpha-alphamax+ny1*dy)**3.d0
       else
          fm2 = 0.d0
       end if
       if(alpha-alphamax+nz1*dz>0) then
          fm3 = (alpha-alphamax+nz1*dz)**3.d0
       else
          fm3 = 0.d0
       end if
       vol = 1.d0/(6.d0*nx1*ny1*nz1)*(alpha**3.d0-f1-f2-f3+fm1+fm2+fm3)
       return
    end subroutine Volume_Calc

    subroutine Newton_Raphson(a3,a2,a1,a0,sc)
       implicit none
       real(dp),intent(in):: a3,a2,a1,a0
       real(dp),intent(inout):: sc
       real(dp),parameter:: tol = 1.d-6
       real(dp):: delvar
       delvar = 1.d6
       do while (dabs(delvar)>tol*dabs(sc))
          delvar = (a3*sc**3.d0+a2*sc**2.d0+a1*sc+a0)/(3.d0*a3*sc**2.d0+2.d0*a2*sc+a1)
          sc = sc-delvar
       end do
    end subroutine Newton_Raphson

    subroutine Volume_Fraction_Calc(dx,dy,dz,nxx,nyy,nzz,s,vol)
      !! The subroutine is used to compute the volume inside the box.
      !! The inputs are box size: dx,dy,dz and the interface information: nxx,nyy,nzz,s.
      !! The method is from "CLSVOF as a fast and mass-conserving extension
      !! level-set method for the simulation of two phases flow" Griebel and Klitz.
      !! Section 5 in page 14  
       real(dp),intent(in)  :: nxx,nyy,nzz,dx,dy,dz,s
       real(dp),intent(out) :: vol
       real(dp)             :: nx1,ny1,nz1,dx1,dy1,dz1
       real(dp)             :: sc,sm,fc,vol1,vol2
       nx1 = dabs(nxx)
       ny1 = dabs(nyy)
       nz1 = dabs(nzz)

       dx1 = dmax1(nx1*dx,ny1*dy,nz1*dz)
       dz1 = dmin1(nx1*dx,ny1*dy,nz1*dz)
       dy1 = (nx1*dx+ny1*dy+nz1*dz)-dx1-dz1
       sm = dx1+dy1+dz1
       sc = dmin1(s,sm-s)
       if(s<=0.d0) then
          vol =  0.d0
          return
       end if
       if(s>=sm) then     ! be careful with this condition. it affects to region contains fluid
          vol =  1.d0*dx*dy*dz
          return
       end if
       if(sc<=dz1) then
          fc = sc**3.d0/(6.d0*dx1*dy1*dz1)
       elseif(sc<=dy1) then
          fc = (sc**2.d0-dz1*sc+dz1**2.d0/3.d0)/(2.d0*dx1*dy1)
       elseif(dx1>=dy1+dz1) then
          if(sc>=dy1+dz1) then
             fc = (2.d0*sc-dy1-dz1)/(2.d0*dx1)
          else
             fc = (sc**3.d0-(sc-dz1)**3.d0-(sc-dy1)**3.d0)
             fc = fc/(6.d0*dx1*dy1*dz1)
          end if
       else
          if(sc>dx1) then
             fc = (sc**3.d0-(sc-dz1)**3.d0-(sc-dy1)**3.d0-(sc-dx1)**3.d0)
          else
             fc = (sc**3.d0-(sc-dz1)**3.d0-(sc-dy1)**3.d0)
          end if
          fc = fc/(6.d0*dx1*dy1*dz1)
       end if
       if(s<=0.5d0*sm) then      ! be careful with this step. the chosen of
          vol = fc*dx*dy*dz      ! the region contains fluid
          return
       else
          vol = (1.d0-fc)*dx*dy*dz
          return
       end if
       if(isnan(vol)) then
          pause 'Volume Fraction Calculate'
       end if
       if(ieee_is_finite(fc)) then
          pause 'infinity in volume fraction calculate 1155'
       end if
    end subroutine Volume_Fraction_Calc
    ! calculate flux throught the east of cell
    subroutine East_Flux(nxx,nyy,nzz,diss,volf,dx,dy,dz,udt,flux)
       implicit none
       real(dp),intent(in):: nxx,nyy,nzz,diss,volf,dx,dy,dz,udt
       real(dp),intent(out):: flux
       real(dp):: eps,vol
       eps = vofeps
       if(udt==0.d0) then
          flux = 0.d0
          return
       end if
       if(nxx<=0) then
          call Volume_Fraction_Calc(udt,dy,dz,nxx,nyy,nzz,diss,vol)
          flux = vol/(udt*dy*dz)
       else
          call Volume_Fraction_Calc(dx-udt,dy,dz,nxx,nyy,nzz,diss,vol)
          flux = (volf*dx*dy*dz-vol)/(udt*dy*dz)
       end if
       if(isnan(flux)) then
          pause 'east flux 1178'
       end if
       flux  = flux*udt/dx
    end subroutine East_Flux

    ! calculate flux throught the west of cell
    subroutine West_Flux(nxx,nyy,nzz,diss,volf,dx,dy,dz,udt,flux)
       implicit none
       real(dp),intent(in):: nxx,nyy,nzz,diss,volf,dx,dy,dz,udt
       real(dp),intent(out):: flux
       real(dp):: eps,vol
       eps = vofeps
       if(udt==0.d0) then
          flux = 0.d0
          return
       end if
       if(nxx>=0) then
          call Volume_Fraction_Calc(udt,dy,dz,nxx,nyy,nzz,diss,vol)
          flux = vol/(udt*dy*dz)
       else
          call Volume_Fraction_Calc(dx-udt,dy,dz,nxx,nyy,nzz,diss,vol)
          flux = (volf*dx*dy*dz-vol)/(udt*dy*dz)
       end if
       if(isnan(flux)) then
          pause 'west flux 778'
       end if
       flux  = flux*udt/dx
    end subroutine West_Flux

    ! calculate flux through the north side of cell
    subroutine North_Flux(nxx,nyy,nzz,diss,volf,dx,dy,dz,vdt,flux)
       implicit none
       real(dp),intent(in):: nxx,nyy,nzz,diss,volf,dx,dy,dz,vdt
       real(dp),intent(out):: flux
       real(dp):: eps,vol
       eps = vofeps
       if(vdt==0.d0) then
          flux = 0.d0
          return
       end if
       if(nyy<=0.d0) then
          call Volume_Fraction_Calc(dx,vdt,dz,nxx,nyy,nzz,diss,vol)
          flux = vol/(dx*vdt*dz)
       else
          call Volume_Fraction_Calc(dx,dy-vdt,dz,nxx,nyy,nzz,diss,vol)
          flux = (volf*dx*dy*dz-vol)/(dx*vdt*dz)
       end if
       flux = flux*vdt/dy
    end subroutine North_Flux

    ! calculate flux throght the south side of cell
    subroutine South_Flux(nxx,nyy,nzz,diss,volf,dx,dy,dz,vdt,flux)
       implicit none
       real(dp),intent(in):: nxx,nyy,nzz,diss,volf,dx,dy,dz,vdt
       real(dp),intent(out):: flux
       real(dp):: vol,eps 
       eps = vofeps
       if(vdt == 0.d0) then
          flux = 0.d0
          return
       end if
       if(vdt==0.d0) then
          flux = 0.d0
          return
       end if
       if(nyy>=0.d0) then
          call Volume_Fraction_Calc(dx,vdt,dz,nxx,nyy,nzz,diss,vol)
          flux = vol/(dx*vdt*dz)
       else
          call Volume_Fraction_Calc(dx,dy-vdt,dz,nxx,nyy,nzz,diss,vol)
          flux = (volf*dx*dy*dz-vol)/(dx*vdt*dz)
       end if
       flux = flux*vdt/dy
    end subroutine South_Flux

    ! calculate the flux through top of cell
    subroutine Top_Flux(nxx,nyy,nzz,diss,volf,dx,dy,dz,wdt,flux)
       implicit none
       real(dp),intent(in):: nxx,nyy,nzz,diss,volf,dx,dy,dz,wdt
       real(dp),intent(out):: flux
       real(dp):: eps,vol
       eps = vofeps
       if(wdt==0.d0) then
          flux = 0.d0
          return
       end if
       if(nzz<=0) then
          call Volume_Fraction_Calc(dx,dy,wdt,nxx,nyy,nzz,diss,vol)
          flux = vol/(dx*dy*wdt)
       else
          call Volume_Fraction_Calc(dx,dy,dz-wdt,nxx,nyy,nzz,diss,vol)
          flux = (volf*dx*dy*dz-vol)/(dx*dy*wdt)
       end if
       flux = flux*wdt/dz
    end subroutine Top_Flux

    ! calculate the flux through bottom of cell
    subroutine Bottom_Flux(nxx,nyy,nzz,diss,volf,dx,dy,dz,wdt,flux)
       implicit none
       real(dp),intent(in):: nxx,nyy,nzz,diss,volf,dx,dy,dz,wdt
       real(dp),intent(out):: flux
       real(dp):: eps,vol
       eps = vofeps
       if(wdt==0.d0) then
          flux = 0.d0
          return
       end if
       if(nzz>=0) then
          call Volume_Fraction_Calc(dx,dy,wdt,nxx,nyy,nzz,diss,vol)
          flux = vol/(dx*dy*wdt)
       else
          call Volume_Fraction_Calc(dx,dy,dz-wdt,nxx,nyy,nzz,diss,vol)
          flux = (volf*dx*dy*dz-vol)/(dx*dy*wdt)
       end if
       flux = flux*wdt/dz
    end subroutine Bottom_Flux

    ! Redistance from vof to level set
    subroutine Redistance_Levelset(PGrid,vfl,phi,Tnx,Tny,Tnz,Tdis)
      !! The subroutine is used to compute the level set value from the interface 
      !! information. The methodology for parallization is given in the page 22 of 
      !! "CLSVOF as a fast and mass-conserving extension
      !! of the level-set method for the simulation of twophase flow problems"
      !! M. Griebel & M. Klitz
      implicit none
      type(Grid),                                 intent(in)    :: PGrid
      !! The pressure grid
      real(kind=dp),dimension(0:,0:,0:),          intent(in)    :: vfl
      !! The volume of fluid
      real(kind=dp), dimension(0:,0:,0:),         intent(inout) :: phi
      !! The level set function
      real(kind=dp),dimension(:,:,:),allocatable, intent(in)    :: Tnx,Tny,Tnz,Tdis
      integer                   :: i,j,k,ii,jj,kk
      integer                   :: iii,jjj,kkk,ij
      integer                   :: l,m,n,cas
      real(kind=dp)				      :: xv,yv,zv,dv
      real(kind=dp)					    :: deuc,dij1
      real(kind=dp)					    :: xp,yp,zp
      real(kind=dp)					    :: xoff,yoff,zoff
      real(kind=dp) 				    :: xfc,yfc,zfc
      real(kind=dp)					    :: xs,ys,zs
      real(kind=dp)					    :: diss1,diss,nxx,nyy,nzz
      logical                   :: pointp,face
      real(kind=dp),dimension(:,:,:),allocatable      :: phiaux
      integer,dimension(:,:,:),allocatable            :: fix
      real(kind=dp)					    :: epsi

      epsi = 1.d-15
      allocate(fix(imax,jmax,kmax))
      allocate(phiaux(imax,jmax,kmax))
      fix(:,:,:) = 0
      phiaux(:,:,:) = 1.d4
      do i = 1,imax
        do j = 1,jmax
          do k = 1,kmax
            if((0.d0+epsi)<vfl(i,j,k).and.vfl(i,j,k)<(1.d0-epsi)) then
              phiaux(i,j,k) = dsign(1.d0,0.5d0-vfl(i,j,k))*dabs(Tdis(i,j,k))
              fix(i,j,k) = 1
              do ii = -band_width,band_width
                do jj = -band_width,band_width
                  do kk = -band_width,band_width
                  ! determine the point xv on the boundary cell(i,j,k) with the
                  ! shortest distance to the cell center of (i+-ii,j+-jj)
                    if(1<=i+ii.and.i+ii<=imax.and.1<=j+jj.and.j+jj<=jmax.and.  &
                                                  1<=k+kk.and.k+kk<=kmax)then
                      if(vfl(i+ii,j+jj,k+kk)<0.d0+epsi.or.             	       &
                         vfl(i+ii,j+jj,k+kk)>1.d0-epsi) then
                        if(fix(i+ii,j+jj,k+kk)==0) fix(i+ii,j+jj,k+kk) = 1

                        nxx=Tnx(i,j,k)
                        nyy=Tny(i,j,k)
                        nzz=Tnz(i,j,k)
                        diss=Tdis(i,j,k)
                        l=max(-1,min(1,ii))
                        m=max(-1,min(1,jj))
                        n=max(-1,min(1,kk))
                        xv=PGrid%dx(i,j,k)*dble(l)/2.d0
                        yv=PGrid%dy(i,j,k)*dble(m)/2.d0
                        zv=PGrid%dz(i,j,k)*dble(n)/2.d0
                        dv=nxx*xv+nyy*yv+nzz*zv+diss
                        if(dv*dsign(1.d0,phi(i+ii,j+jj,k+kk))<=0.d0) then
                          deuc=dsqrt((PGrid%dx(i,j,k)*dble(ii)-xv)**2.d0+      &
                                     (PGrid%dy(i,j,k)*dble(jj)-yv)**2.d0+      &
                                     (PGrid%dz(i,j,k)*dble(kk)-zv)**2.d0)
                          phiaux(i+ii,j+jj,k+kk)=dsign(1.d0,0.5d0-             &
                                                 vfl(i+ii,j+jj,k+kk))*         &
                                         dmin1(deuc,dabs(phiaux(i+ii,j+jj,k+kk)))
                     ! third step: find the projection of x' onto the interface
                        else
                          dij1=nxx*PGrid%dx(i,j,k)*dble(ii)+       	       &
                               nyy*PGrid%dy(i,j,k)*dble(jj)+                   &
                               nzz*PGrid%dz(i,j,k)*dble(kk)+diss
                          xp=PGrid%dx(i,j,k)*dble(ii)-dij1*nxx
                          yp=PGrid%dy(i,j,k)*dble(jj)-dij1*nyy
                          zp=PGrid%dz(i,j,k)*dble(kk)-dij1*nzz
                          pointp=Isinsidecell(xp,yp,zp,PGrid%dx(i,j,k),        &
                                              PGrid%dy(i,j,k),PGrid%dz(i,j,k))
                          if(pointp.eqv..true.) then
                            fix(i+ii,j+jj,k+kk) = 1     !  set to "2" to reduce computational time
                            deuc=dsqrt((PGrid%dx(i,j,k)*dble(ii)-xp)**2.d0+    &
                                       (PGrid%dy(i,j,k)*dble(jj)-yp)**2.d0+    &
                                       (PGrid%dz(i,j,k)*dble(kk)-zp)**2.d0)
                            phiaux(i+ii,j+jj,k+kk)=dsign(1.d0,0.5d0-           &
                                            vfl(i+ii,j+jj,k+kk))*dmin1(deuc,   &
                                            dabs(phiaux(i+ii,j+jj,k+kk)))
                          else
                            ! forth step: find the corner of calculate
                            xoff=dmax1(dabs(xp)-0.5d0*PGrid%dx(i,j,k),0.d0)
                            yoff=dmax1(dabs(yp)-0.5d0*PGrid%dy(i,j,k),0.d0)
                            zoff=dmax1(dabs(zp)-0.5d0*PGrid%dz(i,j,k),0.d0)
                            xfc=dsign(1.d0,xp)*0.5d0*PGrid%dx(i,j,k)
                            yfc=dsign(1.d0,yp)*0.5d0*PGrid%dy(i,j,k)
                            zfc=sign(1.d0,zp)*0.5d0*PGrid%dz(i,j,k)
                            if(xoff*dabs(nxx)>=yoff*dabs(nyy).and.             &
                               xoff*dabs(nxx)>=zoff*dabs(nzz)) then
                              face=iscutface(xfc,nxx,nyy,nzz,diss,	       &
                                   PGrid%dx(i,j,k),PGrid%dy(i,j,k),            &
                                   PGrid%dz(i,j,k),1)
                              if(face.eqv..true.) then
                                xs=xfc
                                diss1=nxx*xs+diss  
                                cas=1
                              else
                                if(yoff*dabs(nyy)>zoff*dabs(nzz)) then
                                  ys=yfc
                                  diss1=nyy*ys+diss
                                  cas=2
                                else
                                  zs=zfc
                                  diss1=nzz*zs+diss
                                  cas = 3
                                endif
                              end if
                              call Shortest_Point_2d(nxx,nyy,nzz,diss1,        &
                                ii,jj,kk,cas,PGrid%dx(i,j,k),PGrid%dy(i,j,k),  &
                                             PGrid%dz(i,j,k),xs,ys,zs)
                            end if
                            if(yoff*dabs(nyy)>=xoff*dabs(nxx).and.             &
                               yoff*dabs(nyy)>=zoff*dabs(nzz)) then
                              face=Iscutface(yfc,nxx,nyy,nzz,diss,             &
                                   PGrid%dx(i,j,k),PGrid%dy(i,j,k),            &
                                   PGrid%dz(i,j,k),2)
                              if(face.eqv..true.) then
                                ys=yfc
                                diss1=nyy*ys+diss
                                cas = 2
                              else
                                if(xoff*dabs(nxx)>zoff*dabs(nzz)) then
                                  xs=xfc
                                  diss1=nxx*xs+diss
                                  cas=1
                                else
                                  zs=zfc
                                  diss1=nzz*zs+diss
                                  cas=3
                                endif
                              end if
                              call Shortest_Point_2d(nxx,nyy,nzz,diss1,        &
                                 ii,jj,kk,cas,PGrid%dx(i,j,k),PGrid%dy(i,j,k), &
                                              PGrid%dz(i,j,k),xs,ys,zs)
                            end if
                            if(zoff*dabs(nzz)>=xoff*dabs(nxx).and.             &
                               zoff*dabs(nzz)>=yoff*dabs(nyy)) then
                              face=Iscutface(zfc,nxx,nyy,nzz,diss,             &
                                   PGrid%dx(i,j,k),PGrid%dy(i,j,k),            &
                                   PGrid%dz(i,j,k),3)
                              if(face.eqv..true.) then
                                zs=zfc
                                diss1=nzz*zs+diss
                                cas=3
                              else
                                if(xoff*dabs(nxx)>yoff*dabs(nyy)) then
                                  xs=xfc
                                  diss1=nxx*xs+diss
                                  cas=1
                                else
                                  ys=yfc
                                  diss1=nyy*ys+diss
                                  cas=2
                                endif
                              end if
                              call Shortest_Point_2d(nxx,nyy,nzz,diss1,        &
                                 ii,jj,kk,cas,PGrid%dx(i,j,k),PGrid%dy(i,j,k), &
                                              PGrid%dz(i,j,k),xs,ys,zs)
                            end if
                            deuc=dsqrt((PGrid%dx(i,j,k)*dble(ii)-xs)**2.d0+    &
                                       (PGrid%dy(i,j,k)*dble(jj)-ys)**2.d0+    &
                                       (PGrid%dz(i,j,k)*dble(kk)-zs)**2.d0)
                            phiaux(i+ii,j+jj,k+kk)=dsign(1.d0,0.5d0-           &
                                        vfl(i+ii,j+jj,k+kk))*                  &
                                        dmin1(deuc,dabs(phiaux(i+ii,j+jj,k+kk)))
                          end if
                        end if
                      end if
                    end if
                  end do
                end do
              end do
            end if
          end do
        end do
      end do
      do i = 1,imax
        do j = 1,jmax
          do k = 1,kmax
            if(fix(i,j,k)==1) then
              phi(i,j,k) = phiaux(i,j,k)
            elseif(fix(i,j,k)==2) then
              phi(i,j,k) = phiaux(i,j,k)
            else
              phi(i,j,k) = (band_width+1)*PGrid%dx(i,j,k)*dsqrt(3.d0)*         &
                                        dsign(1.d0,0.5-vfl(i,j,k))
            end if
          end do
        end do
      end do
      deallocate(fix)
      deallocate(phiaux)
    end subroutine Redistance_Levelset

 !   function Heavisideph(phir) result(dhp)
 !      real(dp):: phir,dhp
 !      real(dp):: epsil
 !      epsil = 1.d-5
 !      if(dabs(phir)<=epsilon) then
 !         dhp = -phir/epsilon**2.d0*(1.d0+dcos(pi*phir/epsilon))
 !      else
 !         dhp = 0.d0
 !      end if
 !   end function heavisideph

    function Isinsidecell(xp,yp,zp,dx,dy,dz) result(logic)
      implicit none
      real(dp),intent(in) :: xp,yp,zp
      real(dp),intent(in) :: dx,dy,dz
      logical:: logic
      if(-0.5d0*dx<=xp.and.xp<=0.5d0*dx.and.	       		               &
         -0.5d0*dy<=yp.and.yp<=0.5d0*dy.and.                                   &
         -0.5d0*dz<=zp.and.zp<=0.5d0*dz) then
        logic = .true.
      else
        logic = .false.
      end if
    end function Isinsidecell

    function Isinsiderect(dire1,dire2,rang10,rang11,rang20,rang21) result(logic)
       implicit none
       real(dp),intent(in):: dire1,dire2,rang10,rang11,rang20,rang21
       logical:: logic
       if(rang10<=dire1.and.dire1<=rang11.and.rang20<=dire2.and.dire2<=rang21)then
          logic = .true.
       else
          logic = .false.
       end if
    end function Isinsiderect

    function Iscutface(face_point,nxx,nyy,nzz,diss,dx,dy,dz,cas) result(logic)
       implicit none
       real(dp),intent(in) :: face_point,nxx,nyy,nzz,diss
       real(dp),intent(in) :: dx,dy,dz
       integer,intent(in)  :: cas
       logical:: logic
       real(dp):: point1,point2,point3,point4
       select case(cas)
          case(1)   ! in x surface
             point1=face_point*nxx-0.5d0*dy*nyy-0.5d0*dz*nzz+diss
             point2=face_point*nxx-0.5d0*dy*nyy+0.5d0*dz*nzz+diss
             point3=face_point*nxx+0.5d0*dy*nyy+0.5d0*dz*nzz+diss
             point4=face_point*nxx+0.5d0*dy*nyy-0.5d0*dz*nzz+diss
          case(2)   ! in y surface
             point1=-0.5d0*dx*nxx+face_point*nyy-0.5d0*dz*nzz+diss
             point2=-0.5d0*dx*nxx+face_point*nyy+0.5d0*dz*nzz+diss
             point3=0.5d0*dx*nxx+face_point*nyy+0.5d0*dz*nzz+diss
             point4=0.5d0*dx*nxx+face_point*nyy-0.5d0*dz*nzz+diss
          case(3)   ! in z surface
             point1=-0.5d0*dx*nxx-0.5d0*dy*nyy+face_point*nzz+diss
             point2=-0.5d0*dx*nxx+0.5d0*dy*nyy+face_point*nzz+diss
             point3=0.5d0*dx*nxx+0.5d0*dy*nyy+face_point*nzz+diss
             point4=0.5d0*dx*nxx-0.5d0*dy*nyy+face_point*nzz+diss
       end select
       if(point1*point2>=0.d0.and.point2*point3>=0.and.point3*point4>=0.d0     &
          .and.point4*point1>=0) then
          logic = .false.
          return
       else
          logic = .true.
          return
       end if
    end function Iscutface

    subroutine Shortest_Point_2d(nx1,ny1,nz1,dis1,ii,jj,kk,cas,dx,dy,dz,xs,ys,zs)
       implicit none
       real(dp),intent(in)    :: nx1,ny1,nz1,dis1
       integer, intent(in)    :: ii,jj,kk,cas
       real(dp),intent(in)    :: dx,dy,dz
       real(dp),intent(inout) :: xs,ys,zs
       real(dp)               :: temp,nxx,nyy,nzz,diss
       real(dp)               :: dij1,xp,yp,zp,xoff,yoff,zoff,xfc,yfc,zfc
       logical                :: pointp
       real(dp)		      :: epsi=1.d-30
       select case(cas)
         case(1) ! for x face
           temp=dsqrt(ny1**2.d0+nz1**2.d0+epsi**2.d0)
           nyy=ny1/temp
           nzz=nz1/temp
           diss=dis1/temp
           dij1=nyy*dy*dble(jj)+nzz*dz*dble(kk)+diss
           yp=dy*dble(jj)-dij1*nyy
           zp=dz*dble(kk)-dij1*nzz
           pointp=Isinsiderect(yp,zp,-0.5d0*dy,0.5d0*dy,-0.5d0*dz,0.5d0*dz)
           if(pointp.eqv..true.) then
             ys = yp
             zs = zp
             return
           else
           ! forth step: find the corner of
             yoff = dmax1(dabs(yp)-0.5d0*dy,0.d0)
             zoff = dmax1(dabs(zp)-0.5d0*dz,0.d0)
             yfc = dsign(1.d0,yp)*0.5d0*dy
             zfc = dsign(1.d0,zp)*0.5d0*dz
             if(yoff*dabs(nyy)>zoff*dabs(nzz)) then
               ys = yfc
               zs = (diss+nyy*ys)/(dsign(1.d0,-nzz)*dmax1(dabs(nzz),tolpar))  
               return
             else
               zs = zfc
               ys = (diss+nzz*zs)/(dsign(1.d0,-nyy)*dmax1(dabs(nyy),tolpar))  
               return
             end if
           end if
         case(2) ! for y face
           temp = dsqrt(nx1**2.d0+nz1**2.d0+epsi**2.d0)
           nxx = nx1/temp
           nzz = nz1/temp
           diss = dis1/temp
           dij1 = nxx*dx*dble(ii)+nzz*dz*dble(kk)+diss
           xp = dx*dble(ii)-dij1*nxx
           zp = dz*dble(kk)-dij1*nzz
           pointp=Isinsiderect(xp,zp,-0.5d0*dx,0.5d0*dx,-0.5d0*dz,0.5d0*dz)
           if(pointp.eqv..true.) then
             xs=xp
             zs=zp
             return
           else
           ! forth step: find the corner of
             xoff = dmax1(dabs(xp)-0.5d0*dx,0.d0)
             zoff = dmax1(dabs(zp)-0.5d0*dz,0.d0)
             xfc = dsign(1.d0,xp)*0.5d0*dx
             zfc = dsign(1.d0,zp)*0.5d0*dz
             if(xoff*dabs(nxx)>zoff*dabs(nzz)) then
               xs = xfc
               zs = (diss+nxx*xs)/(dsign(1.d0,-nzz)*dmax1(dabs(nzz),tolpar))
               return
             else
               zs = zfc
               xs = (diss+nzz*zs)/(dsign(1.d0,-nxx)*dmax1(dabs(nzz),tolpar))
               return
             end if
           end if
         case(3) ! for z face
           temp=dsqrt(nx1**2.d0+ny1**2.d0+epsi**2.d0)
           nxx=nx1/temp
           nyy=ny1/temp
           diss=dis1/temp
           dij1=nxx*dx*dble(ii)+nyy*dy*dble(jj)+diss
           xp=dx*dble(ii)-dij1*nxx
           yp=dy*dble(jj)-dij1*nyy
           pointp=Isinsiderect(xp,yp,-0.5d0*dx,0.5d0*dx,-0.5d0*dy,0.5d0*dy)
           if(pointp.eqv..true.) then
             xs=xp
             ys=yp
             return
           else
           ! forth step: find the corner of
             xoff=dmax1(dabs(xp)-0.5d0*dx,0.d0)
             yoff=dmax1(dabs(yp)-0.5d0*dy,0.d0)
             xfc=dsign(1.d0,xp)*0.5d0*dx
             yfc=dsign(1.d0,yp)*0.5d0*dy
             if(xoff*dabs(nxx)>yoff*dabs(nyy)) then
               xs=xfc
               ys=(diss+nxx*xs)/(dsign(1.d0,-nyy)*dmax1(dabs(nyy),tolpar))
               return
             else
               ys=yfc
               xs=(diss+nyy*ys)/(dsign(1.d0,-nxx)*dmax1(dabs(nxx),tolpar))
               return
             end if
           end if
         end select
    end subroutine Shortest_Point_2d

    subroutine Boundary_Condition(vari)
       implicit none
       integer(kind=it4b) 		       	    :: i,j,k
       real(kind=dp),dimension(:,:,:),intent(inout) :: vari
       do j = 1,jmax
          do k = 1, kmax
             vari(1,j,k) = vari(2,j,k)
             vari(imax,j,k) = vari(imax-1,j,k) 
          end do
       end do
       do i = 1,imax
          do k = 1,kmax
             vari(i,1,k) = vari(i,2,k)
             vari(i,jmax,k) = vari(i,jmax-1,k)
          end do
       end do
       do i = 1,imax
          do j = 1,jmax
             vari(i,j,1) = vari(i,j,2)
             vari(i,j,kmax) = vari(i,j,kmax-1)
          end do
       end do
       return
    end subroutine Boundary_Condition

    subroutine BoundaryConditionLvsVof(PGrid, PCell, Vari, BCLvs, BCVof, Time)
      type(Grid), intent(in)         :: PGrid
      type(Cell), intent(inout)      :: PCell
      type(Variables), intent(in)    :: Vari
      type(BCBase), intent(inout)    :: BCLvs,BCVof
      real(kind=dp), intent(in)      :: Time
      integer(kind=it4b)             :: i,j,k

      ! For the western boundary
      call BCLvs%West(PGrid%x(1,1:Jmax,1:Kmax)-PGrid%dx(1,1:Jmax,1:Kmax)/2.d0,&
                      PGrid%y(1,1:Jmax,1:Kmax), PGrid%z(1,1:Jmax,1:Kmax),     &
                      PGrid%dx(1,1:Jmax,1:Kmax), PGrid%dy(1,1:Jmax,1:Kmax),   &
                      PGrid%dz(1,1:Jmax,1:Kmax), Vari%p(1,1:Jmax,1:Kmax),     &
                      Vari%u(1,1:Jmax,1:Kmax), Vari%v(1,1:Jmax,1:Kmax),       &
                      Vari%w(1,1:Jmax,1:Kmax), PCell%vofL(1,1:Jmax,1:Kmax),   &
                      PCell%phiL(1,1:Jmax,1:Kmax), Time,                     &
                      PCell%nxL(1,1:Jmax,1:Kmax),PCell%nyL(1,1:Jmax,1:Kmax),  &
                      PCell%nzL(1,1:Jmax,1:Kmax))
      PCell%phiL(0,1:Jmax,1:Kmax)=PCell%phiL(1,1:Jmax,1:Kmax)-                &
                          PCell%nx(1,1:Jmax,1:Kmax)*PGrid%dx(1,1:Jmax,1:Kmax)
      ! For the eastern boundary
      call BCLvs%East(PGrid%x(Imax,1:Jmax,1:Kmax)+PGrid%dx(Imax,1:Jmax,1:Kmax)/2.d0, &
                    PGrid%y(Imax,1:Jmax,1:Kmax),PGrid%z(Imax,1:Jmax,1:Kmax),  &
                    PGrid%dx(Imax,1:Jmax,1:Kmax),PGrid%dy(Imax,1:Jmax,1:Kmax),&
                    PGrid%dz(Imax,1:Jmax,1:Kmax),Vari%p(Imax,1:Jmax,1:Kmax),  &
                    Vari%u(Imax-1,1:Jmax,1:Kmax),Vari%v(Imax,1:Jmax,1:Kmax),  &
                    Vari%w(Imax,1:Jmax,1:Kmax),PCell%vofL(Imax,1:Jmax,1:Kmax),&
                    PCell%phiL(Imax,1:Jmax,1:Kmax),Time,                      &
                    PCell%nxL(Imax,1:Jmax,1:Kmax),                            &
                    PCell%nyL(Imax,1:Jmax,1:Kmax),PCell%nzL(Imax,1:Jmax,1:Kmax))
      PCell%phiL(Imax+1,1:Jmax,1:Kmax)=PCell%phiL(Imax,1:Jmax,1:Kmax)+        &
                          PCell%nx(Imax,1:Jmax,1:Kmax)*PGrid%dx(Imax,1:Jmax,1:Kmax)
      ! For the southern boundary
      call BCLvs%South(PGrid%x(1:Imax,1,1:Kmax),                              &
                  PGrid%y(1:Imax,1,1:Kmax)-PGrid%dy(1:Imax,1,1:Kmax)/2.d0,    &
                  PGrid%z(1:Imax,1,1:Kmax), PGrid%dx(1:Imax,1,1:Kmax),        &
                  PGrid%dy(1:Imax,1,1:Kmax), PGrid%dz(1:Imax,1,1:Kmax),       &
                  Vari%p(1:Imax,1,1:Kmax), Vari%u(1:Imax,1,1:Kmax),           &
                  Vari%v(1:Imax,1,1:Kmax), Vari%w(1:Imax,1,1:Kmax),           &
                  PCell%vofL(1:Imax,1,1:Kmax), PCell%phiL(1:Imax,1,1:Kmax), Time)
      PCell%phiL(1:Imax,0,1:Kmax)=PCell%phiL(1:Imax,1,1:Kmax)-                &
                            PCell%ny(1:Imax,1,1:Kmax)*PGrid%dy(1:Imax,1,1:Kmax)
      ! For the northern boundary
      call BCLvs%North(PGrid%x(1:Imax,Jmax,1:Kmax),                           &
               PGrid%y(1:Imax,Jmax,1:Kmax)+PGrid%dy(1:Imax,Jmax,1:Kmax)/2.d0, &
               PGrid%z(1:Imax,Jmax,1:Kmax), PGrid%dx(1:Imax,Jmax,1:Kmax),     &
               PGrid%dy(1:Imax,Jmax,1:Kmax), PGrid%dz(1:Imax,Jmax,1:Kmax),    &
               Vari%p(1:Imax,Jmax,1:Kmax), Vari%u(1:Imax,Jmax,1:Kmax),        &
               Vari%v(1:Imax,Jmax,1:Kmax), Vari%w(1:Imax,Jmax,1:Kmax),        &
               PCell%vofL(1:Imax,Jmax,1:Kmax), PCell%phiL(1:Imax,Jmax,1:Kmax), Time) 
      PCell%phiL(1:Imax,Jmax+1,1:Kmax)=PCell%phiL(1:Imax,Jmax,1:Kmax)-        &
                      PCell%ny(1:Imax,Jmax,1:Kmax)*PGrid%dy(1:Imax,Jmax,1:Kmax)
      ! For the bottom boundary
      call BCLvs%Bottom(PGrid%x(1:Imax,1:Jmax,1), PGrid%y(1:Imax,1:Jmax,1),   &
               PGrid%z(1:Imax,1:Jmax,1)-PGrid%dz(1:Imax,1:Jmax,1)/2.d0,       &
               PGrid%dx(1:Imax,1:Jmax,1), PGrid%dy(1:Imax,1:Jmax,1),          &
               PGrid%dz(1:Imax,1:Jmax,1), Vari%p(1:Imax,1:Jmax,1),            &
               Vari%u(1:Imax,1:Jmax,1), Vari%v(1:Imax,1:Jmax,1),              &
               Vari%w(1:Imax,1:Jmax,1), PCell%vofL(1:Imax,1:Jmax,1),          &
               PCell%phiL(1:Imax,1:Jmax,1), Time)
      PCell%phiL(1:Imax,1:Jmax,0)=PCell%phiL(1:Imax,1:Jmax,1)-                &
                            PCell%nz(1:Imax,1:Jmax,1)*PGrid%dz(1:Imax,1:Jmax,1)
      ! For the top boundary
      call BCLvs%Top(PGrid%x(1:Imax,1:Jmax,Kmax),PGrid%y(1:Imax,1:Jmax,Kmax), &
               PGrid%z(1:Imax,1:Jmax,Kmax)+PGrid%dz(1:Imax,1:Jmax,Kmax)/2.d0, &
               PGrid%dx(1:Imax,1:Jmax,Kmax), PGrid%dy(1:Imax,1:Jmax,Kmax),    &
               PGrid%dz(1:Imax,1:Jmax,Kmax), Vari%p(1:Imax,1:Jmax,Kmax),      &
               Vari%u(1:Imax,1:Jmax,Kmax), Vari%v(1:Imax,1:Jmax,Kmax),        &
               Vari%w(1:Imax,1:Jmax,Kmax), PCell%vofL(1:Imax,1:Jmax,Kmax),    &
               PCell%phiL(1:Imax,1:Jmax,Kmax), Time)
      PCell%phiL(1:Imax,1:Jmax,Kmax+1)=PCell%phiL(1:Imax,1:Jmax,Kmax)+         &
                      PCell%nz(1:Imax,1:Jmax,Kmax)*PGrid%dz(1:Imax,1:Jmax,Kmax)
      

      call BCVof%West(PGrid%x(1,1:Jmax,1:Kmax)-PGrid%dx(1,1:Jmax,1:Kmax)/2.d0, &
               PGrid%y(1,1:Jmax,1:Kmax), PGrid%z(1,1:Jmax,1:Kmax),             &
               PGrid%dx(1,1:Jmax,1:Kmax), PGrid%dy(1,1:Jmax,1:Kmax),           &
               PGrid%dz(1,1:Jmax,1:Kmax), Vari%p(1,1:Jmax,1:Kmax),             &
               Vari%u(1,1:Jmax,1:Kmax), Vari%v(1,1:Jmax,1:Kmax),               &
               Vari%w(1,1:Jmax,1:Kmax), PCell%vofL(1,1:Jmax,1:Kmax),           &
               PCell%phiL(1,1:Jmax,1:Kmax), Time)
      PCell%vofL(0,1:Jmax,1:Kmax)=PCell%vofL(1,1:Jmax,1:Kmax)
      ! For the eastern boundary
      call BCVof%East(PGrid%x(Imax,1:Jmax,1:Kmax)+PGrid%dx(Imax,1:Jmax,1:Kmax)/2.d0, &
               PGrid%y(Imax,1:Jmax,1:Kmax), PGrid%z(Imax,1:Jmax,1:Kmax),       &
               PGrid%dx(Imax,1:Jmax,1:Kmax), PGrid%dy(Imax,1:Jmax,1:Kmax),     &
               PGrid%dz(Imax,1:Jmax,1:Kmax), Vari%p(Imax,1:Jmax,1:Kmax),       &
               Vari%u(Imax-1,:,:), Vari%v(Imax,1:Jmax,1:Kmax),                 &
               Vari%w(Imax,1:Jmax,1:Kmax), PCell%vofL(Imax,1:Jmax,1:Kmax),     &
               PCell%phiL(Imax,1:Jmax,1:Kmax), Time)
      PCell%vofL(Imax+1,1:Jmax,1:Kmax)=PCell%vofL(Imax,1:Jmax,1:Kmax)      
      ! For the southern boundary
      call BCVof%South(PGrid%x(1:Imax,1,1:Kmax),                               &
               PGrid%y(1:Imax,1,1:Kmax)-PGrid%dy(1:Imax,1,1:Kmax)/2.d0,        &
               PGrid%z(1:Imax,1,1:Kmax), PGrid%dx(1:Imax,1,1:Kmax),            &
               PGrid%dy(1:Imax,1,1:Kmax), PGrid%dz(1:Imax,1,1:Kmax),           &
               Vari%p(1:Imax,1,1:Kmax), Vari%u(1:Imax,1,1:Kmax),               &
               Vari%v(1:Imax,1,1:Kmax), Vari%w(1:Imax,1,1:Kmax),               &
               PCell%vofL(1:Imax,1,1:Kmax), PCell%phiL(1:Imax,1,1:Kmax), Time)
      PCell%vofL(1:Imax,0,1:Kmax)=PCell%vofL(1:Imax,1,1:Kmax)      
      ! For the northern boundary
      call BCVof%North(PGrid%x(1:Imax,Jmax,1:Kmax),                            &
               PGrid%y(1:Imax,Jmax,1:Kmax)+PGrid%dy(1:Imax,Jmax,1:Kmax)/2.d0,  &
               PGrid%z(1:Imax,Jmax,1:Kmax), PGrid%dx(1:Imax,Jmax,1:Kmax),      &
               PGrid%dy(1:Imax,Jmax,1:Kmax), PGrid%dz(1:Imax,Jmax,1:Kmax),     &
               Vari%p(1:Imax,Jmax,1:Kmax), Vari%u(1:Imax,Jmax,1:Kmax),         &
               Vari%v(1:Imax,Jmax,1:Kmax), Vari%w(1:Imax,Jmax,1:Kmax),         &
               PCell%vofL(1:Imax,Jmax,1:Kmax), PCell%phiL(1:Imax,Jmax,1:Kmax), Time) 
      PCell%vofL(1:Imax,Jmax+1,1:Kmax)=PCell%vofL(1:Imax,Jmax,1:Kmax)
      ! For the bottom boundary
      call BCVof%Bottom(PGrid%x(1:Imax,1:Jmax,1), PGrid%y(1:Imax,1:Jmax,1),    &
               PGrid%z(1:Imax,1:Jmax,1)-PGrid%dz(1:Imax,1:Jmax,1)/2.d0,        &
               PGrid%dx(1:Imax,1:Jmax,1), PGrid%dy(1:Imax,1:Jmax,1),           &
               PGrid%dz(1:Imax,1:Jmax,1), Vari%p(1:Imax,1:Jmax,1),             &
               Vari%u(1:Imax,1:Jmax,1), Vari%v(1:Imax,1:Jmax,1),               &
               Vari%w(1:Imax,1:Jmax,1), PCell%vofL(1:Imax,1:Jmax,1),           &
               PCell%phiL(1:Imax,1:Jmax,1), Time)
      PCell%vofL(1:Imax,1:Jmax,0)=PCell%vofL(1:Imax,1:Jmax,1)
      ! For the top boundary
      call BCVof%Top(PGrid%x(1:Imax,1:Jmax,Kmax),PGrid%y(1:Imax,1:Jmax,Kmax),  &
               PGrid%z(1:Imax,1:Jmax,Kmax)+PGrid%dz(1:Imax,1:Jmax,Kmax)/2.d0,  &
               PGrid%dx(1:Imax,1:Jmax,Kmax), PGrid%dy(1:Imax,1:Jmax,Kmax),     &
               PGrid%dz(1:Imax,1:Jmax,Kmax), Vari%p(1:Imax,1:Jmax,Kmax),       &
               Vari%u(1:Imax,1:Jmax,Kmax), Vari%v(1:Imax,1:Jmax,Kmax),         &
               Vari%w(1:Imax,1:Jmax,Kmax), PCell%vofL(1:Imax,1:Jmax,Kmax), PCell%phiL(1:Imax,1:Jmax,Kmax), Time)
      PCell%vofL(1:Imax,1:Jmax,Kmax+1)=PCell%vofL(1:Imax,1:Jmax,Kmax)
    end subroutine BoundaryConditionLvsVof 
  ! Boundary condition for fluid field
    subroutine BoundaryConditionLvsVofFluid(PGrid, PCell, Vari, BCLvsF, BCVofF, Time)
      type(Grid), intent(in)         :: PGrid
      type(Cell), intent(inout)      :: PCell
      type(Variables), intent(in)    :: Vari
      type(BCBase), intent(inout)    :: BCLvsF,BCVofF
      real(kind=dp), intent(in)      :: Time
      integer(kind=it4b)             :: i,j,k

      ! For the western boundary
      call BCLvsF%West(PGrid%x(1,1:Jmax,1:Kmax)-PGrid%dx(1,1:Jmax,1:Kmax)/2.d0,&
                      PGrid%y(1,1:Jmax,1:Kmax), PGrid%z(1,1:Jmax,1:Kmax),     &
                      PGrid%dx(1,1:Jmax,1:Kmax), PGrid%dy(1,1:Jmax,1:Kmax),   &
                      PGrid%dz(1,1:Jmax,1:Kmax), Vari%p(1,1:Jmax,1:Kmax),     &
                      Vari%u(1,1:Jmax,1:Kmax), Vari%v(1,1:Jmax,1:Kmax),       &
                      Vari%w(1,1:Jmax,1:Kmax), PCell%vof(1,1:Jmax,1:Kmax),    &
                      PCell%phi(1,1:Jmax,1:Kmax), Time,                       &
                      PCell%nx(1,1:Jmax,1:Kmax),PCell%ny(1,1:Jmax,1:Kmax),    &
                      PCell%nz(1,1:Jmax,1:Kmax))
      PCell%phi(0,1:Jmax,1:Kmax)=PCell%phi(1,1:Jmax,1:Kmax)-                  &
                          PCell%nx(1,1:Jmax,1:Kmax)*PGrid%dx(1,1:Jmax,1:Kmax)
      ! For the eastern boundary
      call BCLvsF%East(PGrid%x(Imax,1:Jmax,1:Kmax)+PGrid%dx(Imax,1:Jmax,1:Kmax)/2.d0, &
                    PGrid%y(Imax,1:Jmax,1:Kmax),PGrid%z(Imax,1:Jmax,1:Kmax),  &
                    PGrid%dx(Imax,1:Jmax,1:Kmax),PGrid%dy(Imax,1:Jmax,1:Kmax),&
                    PGrid%dz(Imax,1:Jmax,1:Kmax),Vari%p(Imax,1:Jmax,1:Kmax),  &
                    Vari%u(Imax-1,1:Jmax,1:Kmax),Vari%v(Imax,1:Jmax,1:Kmax),  &
                    Vari%w(Imax,1:Jmax,1:Kmax),PCell%vof(Imax,1:Jmax,1:Kmax), &
                    PCell%phi(Imax,1:Jmax,1:Kmax),Time,                       &
                    PCell%nx(Imax,1:Jmax,1:Kmax),                             &
                    PCell%ny(Imax,1:Jmax,1:Kmax),PCell%nz(Imax,1:Jmax,1:Kmax))
      PCell%phi(Imax+1,1:Jmax,1:Kmax)=PCell%phi(Imax,1:Jmax,1:Kmax)+          &
                          PCell%nx(Imax,1:Jmax,1:Kmax)*PGrid%dx(Imax,1:Jmax,1:Kmax)
      ! For the southern boundary
      call BCLvsF%South(PGrid%x(1:Imax,1,1:Kmax),                             &
                  PGrid%y(1:Imax,1,1:Kmax)-PGrid%dy(1:Imax,1,1:Kmax)/2.d0,    &
                  PGrid%z(1:Imax,1,1:Kmax), PGrid%dx(1:Imax,1,1:Kmax),        &
                  PGrid%dy(1:Imax,1,1:Kmax), PGrid%dz(1:Imax,1,1:Kmax),       &
                  Vari%p(1:Imax,1,1:Kmax), Vari%u(1:Imax,1,1:Kmax),           &
                  Vari%v(1:Imax,1,1:Kmax), Vari%w(1:Imax,1,1:Kmax),           &
                  PCell%vof(1:Imax,1,1:Kmax), PCell%phi(1:Imax,1,1:Kmax), Time)
      PCell%phi(1:Imax,0,1:Kmax)=PCell%phi(1:Imax,1,1:Kmax)-PCell%ny(1:Imax,1,1:Kmax)*PGrid%dy(1:Imax,1,1:Kmax)
      ! For the northern boundary
      call BCLvsF%North(PGrid%x(1:Imax,Jmax,1:Kmax),                          &
               PGrid%y(1:Imax,Jmax,1:Kmax)+PGrid%dy(1:Imax,Jmax,1:Kmax)/2.d0, &
               PGrid%z(1:Imax,Jmax,1:Kmax), PGrid%dx(1:Imax,Jmax,1:Kmax),     &
               PGrid%dy(1:Imax,Jmax,1:Kmax), PGrid%dz(1:Imax,Jmax,1:Kmax),    &
               Vari%p(1:Imax,Jmax,1:Kmax), Vari%u(1:Imax,Jmax,1:Kmax),        &
               Vari%v(1:Imax,Jmax,1:Kmax), Vari%w(1:Imax,Jmax,1:Kmax),        &
               PCell%vof(1:Imax,Jmax,1:Kmax), PCell%phi(1:Imax,Jmax,1:Kmax), Time) 
      PCell%phi(1:Imax,Jmax+1,1:Kmax)=PCell%phi(1:Imax,Jmax,1:Kmax)-PCell%ny(1:Imax,Jmax,1:Kmax)*PGrid%dy(1:Imax,Jmax,1:Kmax)
      ! For the bottom boundary
      call BCLvsF%Bottom(PGrid%x(1:Imax,1:Jmax,1), PGrid%y(1:Imax,1:Jmax,1),   &
               PGrid%z(1:Imax,1:Jmax,1)-PGrid%dz(1:Imax,1:Jmax,1)/2.d0,        &
               PGrid%dx(1:Imax,1:Jmax,1), PGrid%dy(1:Imax,1:Jmax,1),           &
               PGrid%dz(1:Imax,1:Jmax,1), Vari%p(1:Imax,1:Jmax,1),             &
               Vari%u(1:Imax,1:Jmax,1), Vari%v(1:Imax,1:Jmax,1),               &
               Vari%w(1:Imax,1:Jmax,1), PCell%vof(1:Imax,1:Jmax,1),            &
               PCell%phi(1:Imax,1:Jmax,1), Time)
      PCell%phi(1:Imax,1:Jmax,0)=PCell%phi(1:Imax,1:Jmax,1)-PCell%nz(1:Imax,1:Jmax,1)*PGrid%dz(1:Imax,1:Jmax,1)
      ! For the top boundary
      call BCLvsF%Top(PGrid%x(1:Imax,1:Jmax,Kmax),PGrid%y(1:Imax,1:Jmax,Kmax), &
               PGrid%z(1:Imax,1:Jmax,Kmax)+PGrid%dz(1:Imax,1:Jmax,Kmax)/2.d0,  &
               PGrid%dx(1:Imax,1:Jmax,Kmax), PGrid%dy(1:Imax,1:Jmax,Kmax),     &
               PGrid%dz(1:Imax,1:Jmax,Kmax), Vari%p(1:Imax,1:Jmax,Kmax),       &
               Vari%u(1:Imax,1:Jmax,Kmax), Vari%v(1:Imax,1:Jmax,Kmax),         &
               Vari%w(1:Imax,1:Jmax,Kmax), PCell%vof(1:Imax,1:Jmax,Kmax),      &
               PCell%phi(1:Imax,1:Jmax,Kmax), Time)
      PCell%phi(1:Imax,1:Jmax,Kmax+1)=PCell%phi(1:Imax,1:Jmax,Kmax)+PCell%nz(1:Imax,1:Jmax,Kmax)*PGrid%dz(1:Imax,1:Jmax,Kmax)
      

      call BCVofF%West(PGrid%x(1,1:Jmax,1:Kmax)-PGrid%dx(1,1:Jmax,1:Kmax)/2.d0, &
               PGrid%y(1,1:Jmax,1:Kmax), PGrid%z(1,1:Jmax,1:Kmax),             &
               PGrid%dx(1,1:Jmax,1:Kmax), PGrid%dy(1,1:Jmax,1:Kmax),           &
               PGrid%dz(1,1:Jmax,1:Kmax), Vari%p(1,1:Jmax,1:Kmax),             &
               Vari%u(1,1:Jmax,1:Kmax), Vari%v(1,1:Jmax,1:Kmax),               &
               Vari%w(1,1:Jmax,1:Kmax), PCell%vof(1,1:Jmax,1:Kmax),           &
               PCell%phi(1,1:Jmax,1:Kmax), Time)
      PCell%vof(0,1:Jmax,1:Kmax)=PCell%vof(1,1:Jmax,1:Kmax)
      ! For the eastern boundary
      call BCVofF%East(PGrid%x(Imax,1:Jmax,1:Kmax)+PGrid%dx(Imax,1:Jmax,1:Kmax)/2.d0, &
               PGrid%y(Imax,1:Jmax,1:Kmax), PGrid%z(Imax,1:Jmax,1:Kmax),       &
               PGrid%dx(Imax,1:Jmax,1:Kmax), PGrid%dy(Imax,1:Jmax,1:Kmax),     &
               PGrid%dz(Imax,1:Jmax,1:Kmax), Vari%p(Imax,1:Jmax,1:Kmax),       &
               Vari%u(Imax-1,:,:), Vari%v(Imax,1:Jmax,1:Kmax),                 &
               Vari%w(Imax,1:Jmax,1:Kmax), PCell%vof(Imax,1:Jmax,1:Kmax),     &
               PCell%phi(Imax,1:Jmax,1:Kmax), Time)
      PCell%vof(Imax+1,1:Jmax,1:Kmax)=PCell%vof(Imax,1:Jmax,1:Kmax)      
      ! For the southern boundary
      call BCVofF%South(PGrid%x(1:Imax,1,1:Kmax),                               &
               PGrid%y(1:Imax,1,1:Kmax)-PGrid%dy(1:Imax,1,1:Kmax)/2.d0,        &
               PGrid%z(1:Imax,1,1:Kmax), PGrid%dx(1:Imax,1,1:Kmax),            &
               PGrid%dy(1:Imax,1,1:Kmax), PGrid%dz(1:Imax,1,1:Kmax),           &
               Vari%p(1:Imax,1,1:Kmax), Vari%u(1:Imax,1,1:Kmax),               &
               Vari%v(1:Imax,1,1:Kmax), Vari%w(1:Imax,1,1:Kmax), PCell%vof(1:Imax,1,1:Kmax),         &
               PCell%phi(1:Imax,1,1:Kmax), Time)
      PCell%vof(1:Imax,0,1:Kmax)=PCell%vof(1:Imax,1,1:Kmax)      
      ! For the northern boundary
      call BCVofF%North(PGrid%x(1:Imax,Jmax,1:Kmax),                                      &
               PGrid%y(1:Imax,Jmax,1:Kmax)+PGrid%dy(1:Imax,Jmax,1:Kmax)/2.d0,                &
               PGrid%z(1:Imax,Jmax,1:Kmax), PGrid%dx(1:Imax,Jmax,1:Kmax),      &
               PGrid%dy(1:Imax,Jmax,1:Kmax), PGrid%dz(1:Imax,Jmax,1:Kmax),     &
               Vari%p(1:Imax,Jmax,1:Kmax), Vari%u(1:Imax,Jmax,1:Kmax),         &
               Vari%v(1:Imax,Jmax,1:Kmax), Vari%w(1:Imax,Jmax,1:Kmax),         &
               PCell%vof(1:Imax,Jmax,1:Kmax), PCell%phi(1:Imax,Jmax,1:Kmax), Time) 
      PCell%vof(1:Imax,Jmax+1,1:Kmax)=PCell%vof(1:Imax,Jmax,1:Kmax)
      ! For the bottom boundary
      call BCVofF%Bottom(PGrid%x(1:Imax,1:Jmax,1), PGrid%y(1:Imax,1:Jmax,1),                         &
               PGrid%z(1:Imax,1:Jmax,1)-PGrid%dz(1:Imax,1:Jmax,1)/2.d0,        &
               PGrid%dx(1:Imax,1:Jmax,1), PGrid%dy(1:Imax,1:Jmax,1),           &
               PGrid%dz(1:Imax,1:Jmax,1), Vari%p(1:Imax,1:Jmax,1),             &
               Vari%u(1:Imax,1:Jmax,1), Vari%v(1:Imax,1:Jmax,1),               &
               Vari%w(1:Imax,1:Jmax,1), PCell%vof(1:Imax,1:Jmax,1),           &
               PCell%phi(1:Imax,1:Jmax,1), Time)
      PCell%vof(1:Imax,1:Jmax,0)=PCell%vof(1:Imax,1:Jmax,1)
      ! For the top boundary
      call BCVofF%Top(PGrid%x(1:Imax,1:Jmax,Kmax),PGrid%y(1:Imax,1:Jmax,Kmax),                      &
               PGrid%z(1:Imax,1:Jmax,Kmax)+PGrid%dz(1:Imax,1:Jmax,Kmax)/2.d0,                   &
               PGrid%dx(1:Imax,1:Jmax,Kmax), PGrid%dy(1:Imax,1:Jmax,Kmax), PGrid%dz(1:Imax,1:Jmax,Kmax),  &
               Vari%p(1:Imax,1:Jmax,Kmax), Vari%u(1:Imax,1:Jmax,Kmax), Vari%v(1:Imax,1:Jmax,Kmax),        &
               Vari%w(1:Imax,1:Jmax,Kmax), PCell%vof(1:Imax,1:Jmax,Kmax), PCell%phi(1:Imax,1:Jmax,Kmax), Time)
      PCell%vof(1:Imax,1:Jmax,Kmax+1)=PCell%vof(1:Imax,1:Jmax,Kmax)
    end subroutine BoundaryConditionLvsVofFluid

    subroutine Normal_Vector_Irre(PGrid,vfl,phi,i,j,k,nxx,nyy,nzz)
       implicit none
       type(Grid),                        intent(in)  :: PGrid
       real(kind=dp),dimension(0:,0:,0:), intent(in)  :: vfl,phi
       integer,                           intent(in)  :: i,j,k
       real(dp),                          intent(out) :: nxx,nyy,nzz
       integer                                     :: case_deri,dx,dy,dz
       real(dp)                                    :: qi(-1:1),qj(-1:1),qk(-1:1)
       real(dp)                                    :: vx,vy,vz
       ! define qi for Dx
       if(i>=3) then
       ! for qi(-1)
         vx=(phi(i,j,k)-phi(i-2,j,k))/(PGrid%x(i,j,k)-PGrid%x(i-2,j,k))
       else
         vx=(phi(i,j,k)-phi(i-1,j,k))/(PGrid%x(i,j,k)-PGrid%x(i-1,j,k))
       end if
       vy=(phi(i-1,j+1,k)-phi(i-1,j-1,k))/(PGrid%y(i-1,j+1,k)-PGrid%y(i-1,j-1,k))
       vz=(phi(i-1,j,k+1)-phi(i-1,j,k-1))/(PGrid%z(i-1,j,k+1)-PGrid%z(i-1,j,k-1))
       qi(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! for qi(1)
       if(i<=imax-2) then
         vx=(phi(i+2,j,k)-phi(i,j,k))/(PGrid%x(i+2,j,k)-PGrid%x(i,j,k))
       else
         vx=(phi(i+1,j,k)-phi(i,j,k))/(PGrid%x(i+1,j,k)-PGrid%x(i,j,k))
       end if
       vy=(phi(i+1,j+1,k)-phi(i+1,j-1,k))/(PGrid%y(i+1,j+1,k)-PGrid%y(i+1,j-1,k))
       vz=(phi(i+1,j,k+1)-phi(i+1,j,k-1))/(PGrid%z(i+1,j,k+1)-PGrid%z(i+1,j,k-1))
       qi(1) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       if(qi(-1)<eta.and.qi(1)>=eta) then
         dx=-1
       elseif(qi(-1)>=eta.and.qi(1)<eta) then
         dx=1
       else
         dx=0
       end if
       ! define qj for Dy
       if(j>=3) then
         vy=(phi(i,j,k)-phi(i,j-2,k))/(PGrid%y(i,j,k)-PGrid%y(i,j-2,k))
       else
         vy=(phi(i,j,k)-phi(i,j-1,k))/(PGrid%y(i,j,k)-PGrid%y(i,j-1,k))
       end if
       vx=(phi(i+1,j-1,k)-phi(i-1,j-1,k))/(PGrid%x(i+1,j-1,k)-PGrid%x(i-1,j-1,k))
       vz=(phi(i,j-1,k+1)-phi(i,j-1,k-1))/(PGrid%z(i,j-1,k+1)-PGrid%z(i,j-1,k-1))
       qj(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! define qj(1)
       vx=(phi(i+1,j+1,k)-phi(i-1,j+1,k))/(PGrid%x(i+1,j+1,k)-PGrid%x(i-1,j+1,k))
       if(j<=jmax-2) then
         vy=(phi(i,j+2,k)-phi(i,j,k))/(PGrid%y(i,j+2,k)-PGrid%y(i,j,k))
       else
         vy=(phi(i,j+1,k)-phi(i,j,k))/(PGrid%y(i,j+1,k)-PGrid%y(i,j,k))
       end if
       vz=(phi(i,j+1,k+1)-phi(i,j+1,k-1))/(PGrid%z(i,j+1,k+1)-PGrid%z(i,j+1,k-1))
       qj(1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       if(qj(-1)<eta.and.qj(1)>=eta) then
         dy=-1
       elseif(qj(-1)>=eta.and.qj(1)<eta) then
         dy=1
       else
         dy=0
       end if
       ! define qk for dz
       ! qk(-1)
       vx=(phi(i+1,j,k-1)-phi(i-1,j,k-1))/(PGrid%x(i+1,j,k-1)-PGrid%x(i-1,j,k-1))
       vy=(phi(i,j+1,k-1)-phi(i,j-1,k-1))/(PGrid%y(i,j+1,k-1)-PGrid%y(i,j-1,k-1))
       if(k>=3) then
         vz=(phi(i,j,k)-phi(i,j,k-2))/(PGrid%z(i,j,k)-PGrid%z(i,j,k-2))
       else
         vz=(phi(i,j,k)-phi(i,j,k-1))/(PGrid%z(i,j,k)-PGrid%z(i,j,k-1))
       end if
       qk(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! qk(1)
       vx=(phi(i+1,j,k+1)-phi(i-1,j,k+1))/(PGrid%x(i+1,j,k+1)-PGrid%x(i-1,j,k+1))
       vy=(phi(i,j+1,k+1)-phi(i,j-1,k+1))/(PGrid%y(i,j+1,k+1)-PGrid%y(i,j-1,k+1))
       if(k<=kmax-2) then
         vz=(phi(i,j,k+2)-phi(i,j,k))/(PGrid%z(i,j,k+2)-PGrid%z(i,j,k))
       else
         vz=(phi(i,j,k+1)-phi(i,j,k))/(PGrid%z(i,j,k+1)-PGrid%z(i,j,k))
       end if
       qk(1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       if(qk(-1)<eta.and.qk(1)>=eta) then
         dz=-1
       elseif(qk(-1)>=eta.and.qk(1)<eta) then
         dz=1
       else
         dz=0
       end if
       if(dx==-1) nxx = (phi(i,j,k)-phi(i-1,j,k))/(PGrid%x(i,j,k)-PGrid%x(i-1,j,k))
       if(dx==1) nxx = (phi(i+1,j,k)-phi(i,j,k))/(PGrid%x(i+1,j,k)-PGrid%x(i,j,k))
       if(dx==0) nxx = (phi(i+1,j,k)-phi(i-1,j,k))/(PGrid%x(i+1,j,k)-PGrid%x(i-1,j,k))
       if(dy==-1) nyy = (phi(i,j,k)-phi(i,j-1,k))/(PGrid%y(i,j,k)-PGrid%y(i,j-1,k))
       if(dy==1) nyy = (phi(i,j+1,k)-phi(i,j,k))/(PGrid%y(i,j+1,k)-PGrid%y(i,j,k))
       if(dy==0) nyy = (phi(i,j+1,k)-phi(i,j-1,k))/(PGrid%y(i,j+1,k)-PGrid%y(i,j-1,k))
       if(dz==-1) nzz = (phi(i,j,k)-phi(i,j,k-1))/(PGrid%z(i,j,k)-PGrid%z(i,j,k-1))
       if(dz==1) nzz = (phi(i,j,k+1)-phi(i,j,k))/(PGrid%z(i,j,k+1)-PGrid%z(i,j,k))
       if(dz==0) nzz = (phi(i,j,k+1)-phi(i,j,k-1))/(PGrid%z(i,j,k+1)-PGrid%z(i,j,k-1))
    end subroutine Normal_Vector_Irre

    real(kind=dp) Function ConvertTheCoordinateSystem(dx,dy,dz,NxIn,NyIn,NzIn,phiIn) result(dout)
      !! The subroutine is used to change the coordinate system from cell centre origin to corner origin
      !! The formulation is derived from the work "CLSVOF as a fast and mass-conserving extension
      !! of the level-set method for the simulation of twophase flow problems" Griebel and Klitz.
      real(kind=dp),intent(in) :: dx,dy,dz
      real(kind=dp),intent(in) :: NxIn,NyIn,NzIn,phiIn
      
      dout=phiIn+0.5d0*(dabs(NxIn)*dx+dabs(NyIn)*dy+dabs(NzIn)*dz)
    end Function  

    real(kind=dp) Function DivideError(nume,deno,opt) result(rs)
      real(kind=dp),intent(in)      :: nume,deno
      integer(kind=it4b),intent(in) :: opt
      real(kind=dp)                 :: tol=1.d-15
      
      select case(opt)
      case(1)
        rs=nume/deno
        if(isnan(rs).or.rs>1.d100) then
          print*,'Overflow with division!!!'
          print*,nume,deno
          print*,'result will be set to 1.d15'
          rs=1.d15
        elseif(dabs(rs)<1.d-20.and.dabs(rs)>1.d-100) then
          print*,'Underflow with division!!!'
          print*,nume,deno
          print*,'result will be set to 1.d-15'
          rs=1.d-15
        end if  
      case(2)
        rs=(nume+tol)/(deno+tol)
      end select    
    end function  

    subroutine CorrectNormalVector(nxx,nyy,nzz)
      real(kind=dp),intent(inout) :: nxx,nyy,nzz
      real(kind=dp)               :: temp
      
      temp = dsqrt(nxx**2.d0+nyy**2.d0+nzz**2.d0)
      if(isnan(temp)) then
        print*, nxx,nyy,nzz
        pause 'clsvof 2138'
      end if
      if(temp<1.d-14) then
        nxx=dsign(1.d0,nxx)*n2
        nyy=dsign(1.d0,nyy)*n2
        nzz=dsign(1.d0,nzz)*n1
      else
        nxx = nxx/temp
        nyy = nyy/temp
        nzz = nzz/temp
        if(nxx>=n1) then
          nxx=dsign(1.d0,nxx)*n1
          nyy=dsign(1.d0,nyy)*n2
          nzz=dsign(1.d0,nzz)*n2
        end if
        if(nyy>=n1) then
          nxx=dsign(1.d0,nxx)*n2
          nyy=dsign(1.d0,nyy)*n1
          nzz=dsign(1.d0,nzz)*n2
        end if
        if(nzz>=n1) then
          nxx=dsign(1.d0,nxx)*n2
          nyy=dsign(1.d0,nyy)*n2
          nzz=dsign(1.d0,nzz)*n1
        end if
      end if  
    end subroutine  
end module Clsvof

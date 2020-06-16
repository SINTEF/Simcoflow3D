Module Clsvof
    USE ieee_arithmetic
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    
    implicit none
    private
    
    integer,parameter			   :: band_width = 4,nv = 6,nl = 6
    real(kind=dp),dimension(:,:,:),pointer :: vfl,vflF    ! vfl represents the liquid volume fraction, vflF represents fluid volume fraction 
    real(kind=dp),dimension(:,:,:),pointer :: phi,phiF    ! phi represents the liquid level set function, phiF represents fluid level set function 
    real(kind=dp),dimension(:,:,:),pointer :: nxF,nyF,nzF ! nxF,nyF,nzF is fluid level set function
    real(dp),parameter                     :: eta=0.075d0,vofeps=1.d-14
    
    public:: InitialClsvofFluidField,InitialClsvofLiquidField,Clsvof_Scheme
    
    interface InitialClsvofFluidFiel
      module procedure InitialClsvofFluidField
    end interface
    
    interface InitialClsvofLiquidField
      module procedure InitialClsvofLiquidField
    end interface
    
    interface ClsVof_Scheme
      module procedure ClsVof_Scheme
    end interface
    
    contains
    
    subroutine InitialClsvofFluidField(TGrid,TCell)
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	      :: i,j,k
      real(kind=dp)		      :: dx,dy,dz,dis,vol,Radius,epsi,s
      real(kind=dp)		      :: tol
      
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
      nullify(vfl)
      nullify(phi)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
    end subroutine InitialClsvofFluidField
    
    subroutine InitialClsvofLiquidField(TGrid,TCell)
      implicit none
      type(Grid),intent(in)           :: TGrid
      type(Cell),intent(inout),target :: TCell
      integer(kind=it4b)	      :: i,j,k
      real(kind=dp)		      :: dx,dy,dz,dis,vol,Radius,epsi,s
      real(kind=dp)		      :: tol
      
      vfl => TCell%vofL
      phi => TCell%phiL
      nxF => TCell%nxL
      nyF => TCell%nyL
      nzF => TCell%nzL
      
      do i=1,Imax
        do j=1,Jmax
          do k=1,Kmax
            vfl(i,j,k)=1.d0*TCell%vof(i,j,k)
            phi(i,j,k)=-1.d3
            nxF(i,j,k)=0.d0
            nyF(i,j,k)=0.d0
            nzF(i,j,k)=1.d0
          end do
        end do
      end do    
      
      nullify(vfl)
      nullify(phi)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
    end subroutine InitialClsvofLiquidField

    subroutine Clsvof_Scheme(PGrid,PCell,UCell,VCell,WCell,TVar,dt,itt)
      implicit none
      type(Grid),intent(in) 		    :: PGrid
      type(Cell),intent(inout),target       :: PCell
      type(Cell),intent(in)		    :: UCell,VCell,WCell
      type(Variables),intent(in)	    :: TVar
      integer(it8b),intent(in)	    	    :: itt
      real(dp),intent(in)		    :: dt
      integer(it4b)			    :: i,j,k,kk
      real(dp)				    :: dtv
      real(dp),dimension(:,:,:),allocatable :: ue,ve,we,nx,ny,nz,dis
      real(dp),dimension(:,:,:),allocatable :: temvfx,temvfy,temvfz
      real(dp),dimension(:,:,:),allocatable :: temlsx,temlsy,temlsz
      
      allocate(nx(imax,jmax,kmax))
      allocate(ny(imax,jmax,kmax))
      allocate(nz(imax,jmax,kmax))
      allocate(dis(imax,jmax,kmax))
      allocate(ue(imax,jmax,kmax))
      allocate(ve(imax,jmax,kmax))
      allocate(we(imax,jmax,kmax))
      allocate(temvfx(imax,jmax,kmax))
      allocate(temvfy(imax,jmax,kmax))
      allocate(temvfz(imax,jmax,kmax))
      allocate(temlsx(imax,jmax,kmax))
      allocate(temlsy(imax,jmax,kmax))
      allocate(temlsz(imax,jmax,kmax))
      
      vfl=>PCell%vofL
      phi=>PCell%phiL
      
      vflF=>PCell%Vof
      phiF=>PCell%phi
      nxF=>PCell%nx
      nyF=>PCell%ny
      nzF=>PCell%nz
      dtv=dt/dble(nv)
      
      do i = 1,imax
        do j = 1,jmax
          do k = 1,kmax
            ue(i,j,k) = TVar%u(i,j,k)
          end do
        end do
      end do
      
      do j = 1,jmax
        do i = 1,imax
          do k = 1,kmax
            ve(i,j,k) = Tvar%v(i,j,k)
          end do
        end do
      end do
      
      do k = 1,kmax
        do i = 1,imax
          do j = 1,jmax
            we(i,j,k) = TVar%w(i,j,k)
          end do
        end do
      end do
      
      do kk=1,nv
        if(mod(kk,3)==0) then
          temvfx(:,:,:) = vfl(:,:,:)
          temlsx(:,:,:) = phi(:,:,:)
          call X_Sweep(PGrid,temvfx,temlsx,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                temvfx(i,j,k) = temvfx(i,j,k)/(1.d0-dtv/PGrid%dx(i,j,k)*(ue(i,j,k)-ue(i-1,j,k)))
                temlsx(i,j,k) = temlsx(i,j,k)/(1.d0-dtv/PGrid%dx(i,j,k)*(ue(i,j,k)-ue(i-1,j,k))) 
                if(temvfx(i,j,k)<=vofeps) temvfx(i,j,k) = 0.d0
                if(temvfx(i,j,k)>=(1.d0-vofeps)) temvfx(i,j,k) = 1.d0
                vfl(i,j,k) = temvfx(i,j,k)
                phi(i,j,k) = temlsx(i,j,k)
              end do
            end do
          end do
          call Boundary_Condition(vfl)
          call Boundary_Condition(phi)
          call Boundary_Condition(temvfx)
          call Boundary_Condition(temlsx)
        
          temvfy(:,:,:) = vfl(:,:,:)
          temlsy(:,:,:) = phi(:,:,:)
          call Y_Sweep(PGrid,temvfy,temlsy,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                temvfy(i,j,k)=temvfy(i,j,k)+dtv/PGrid%dy(i,j,k)*temvfx(i,j,k)* &
                                                (ve(i,j,k)-ve(i,j-1,k))
                temlsy(i,j,k)=temlsy(i,j,k)+dtv/PGrid%dy(i,j,k)*temlsx(i,j,k)* &
                                                (ve(i,j,k)-ve(i,j-1,k))
                if(temvfy(i,j,k)<=vofeps) temvfy(i,j,k) = 0.d0
                if(temvfy(i,j,k)>=(1.d0-vofeps)) temvfy(i,j,k) = 1.d0
                vfl(i,j,k) = temvfy(i,j,k)
                phi(i,j,k) = temlsy(i,j,k)
              end do
            end do
          end do
          call Boundary_Condition(vfl)
          call Boundary_Condition(phi)
        
          temvfz(:,:,:) = vfl(:,:,:)
          temlsz(:,:,:) = phi(:,:,:)
          call Z_Sweep(PGrid,temvfz,temlsz,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                vfl(i,j,k) = temvfz(i,j,k)+dtv/PGrid%dz(i,j,k)*temvfx(i,j,k)*  &
                                                        (we(i,j,k)-we(i,j,k-1))
                phi(i,j,k) = temlsz(i,j,k)+dtv/PGrid%dz(i,j,k)*temlsx(i,j,k)*  &
                                                      (we(i,j,k)-we(i,j,k-1))
                if(vfl(i,j,k)<=vofeps) vfl(i,j,k) = 0.d0
                if(vfl(i,j,k)>=(1.d0-vofeps)) vfl(i,j,k) = 1.d0
              end do
            end do
          end do
        elseif(mod(kk,3)==1) then
          temvfy(:,:,:) = vfl(:,:,:)
          temlsy(:,:,:) = phi(:,:,:)
          call Y_Sweep(PGrid,temvfy,temlsy,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                temvfy(i,j,k) = temvfy(i,j,k)/(1.d0-dtv/PGrid%dy(i,j,k)*       &
                                                        (ve(i,j,k)-ve(i,j-1,k)))
                temlsy(i,j,k) = temlsy(i,j,k)/(1.d0-dtv/PGrid%dy(i,j,k)*       &
                                                        (ve(i,j,k)-ve(i,j-1,k))) 
                if(temvfy(i,j,k)<=vofeps) temvfy(i,j,k) = 0.d0
                if(temvfy(i,j,k)>=(1.d0-vofeps)) temvfy(i,j,k) = 1.d0
                vfl(i,j,k) = temvfy(i,j,k)
                phi(i,j,k) = temlsy(i,j,k)
              end do
            end do
          end do
          call Boundary_Condition(vfl)
          call Boundary_Condition(phi)
          call Boundary_Condition(temvfy)
          call Boundary_Condition(temlsy)
          temvfz(:,:,:) = vfl(:,:,:)
          temlsz(:,:,:) = phi(:,:,:)
          call Z_Sweep(PGrid,temvfz,temlsz,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                temvfz(i,j,k)=temvfz(i,j,k)+dtv/PGrid%dz(i,j,k)*temvfy(i,j,k)* &
                                                      (we(i,j,k)-we(i,j,k-1))
                temlsz(i,j,k)=temlsz(i,j,k)+dtv/PGrid%dz(i,j,k)*temlsy(i,j,k)* &
                                                      (we(i,j,k)-we(i,j,k-1))
                if(temvfz(i,j,k)<=vofeps) temvfz(i,j,k)=0.d0
                if(temvfz(i,j,k)>=(1.d0-vofeps)) temvfz(i,j,k)=1.d0
                vfl(i,j,k) = temvfz(i,j,k)
                phi(i,j,k) = temlsz(i,j,k)
              end do
            end do
          end do
          call Boundary_Condition(vfl)
          call Boundary_Condition(phi)
          temvfx(:,:,:) = vfl(:,:,:)
          temlsx(:,:,:) = phi(:,:,:)
          call X_Sweep(PGrid,temvfx,temlsx,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                vfl(i,j,k)=temvfx(i,j,k)+dtv/PGrid%dx(i,j,k)*temvfy(i,j,k)*    &
                                                   (ue(i,j,k)-ue(i-1,j,k))
                phi(i,j,k)=temlsx(i,j,k)+dtv/PGrid%dx(i,j,k)*temlsy(i,j,k)*    &
                                                   (ue(i,j,k)-ue(i-1,j,k))
                if(vfl(i,j,k)<=vofeps) vfl(i,j,k)=0.d0
                if(vfl(i,j,k)>=(1.d0-vofeps)) vfl(i,j,k) = 1.d0
              end do
            end do
          end do
        else 
          temvfz(:,:,:) = vfl(:,:,:)
          temlsz(:,:,:) = phi(:,:,:)
          call Z_Sweep(PGrid,temvfz,temlsz,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                temvfz(i,j,k)=temvfz(i,j,k)/(1.d0-dtv/PGrid%dz(i,j,k)*         &
                                                   (we(i,j,k)-we(i,j,k-1)))
                temlsz(i,j,k) = temlsz(i,j,k)/(1.d0-dtv/PGrid%dz(i,j,k)*       &
                                                   (we(i,j,k)-we(i,j,k-1)))
                if(temvfz(i,j,k)<=vofeps) temvfz(i,j,k)=0.d0
                if(temvfz(i,j,k)>=(1.d0-vofeps)) temvfz(i,j,k)=1.d0
                vfl(i,j,k) = temvfz(i,j,k)
                phi(i,j,k) = temlsz(i,j,k)
              end do
            end do
          end do
          call Boundary_Condition(vfl)
          call Boundary_Condition(phi)
          call Boundary_Condition(temvfz)
          call Boundary_Condition(temlsz)
          temvfx(:,:,:) = vfl(:,:,:)
          temlsx(:,:,:) = phi(:,:,:)
          call X_Sweep(PGrid,temvfx,temlsx,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                temvfx(i,j,k)=temvfx(i,j,k)+dtv/PGrid%dx(i,j,k)*temvfz(i,j,k)* &
                                                         (ue(i,j,k)-ue(i-1,j,k))
                temlsx(i,j,k)=temlsx(i,j,k)+dtv/PGrid%dx(i,j,k)*temlsz(i,j,k)* &
                                                         (ue(i,j,k)-ue(i-1,j,k))
                if(temvfx(i,j,k)<=vofeps) temvfx(i,j,k)=0.d0
                if(temvfx(i,j,k)>=(1.d0-vofeps)) temvfx(i,j,k) = 1.d0 
                vfl(i,j,k) = temvfx(i,j,k)
                phi(i,j,k) = temlsx(i,j,k)
              end do
            end do
          end do
          call Boundary_Condition(vfl)
          call Boundary_Condition(phi)  
          temvfy(:,:,:) = vfl(:,:,:)
          temlsy(:,:,:) = phi(:,:,:)
          call Y_Sweep(PGrid,temvfy,temlsy,ue,ve,we,nx,ny,nz,dis,dtv)
          do i = 2,imax-1
            do j = 2,jmax-1
              do k = 2,kmax-1
                vfl(i,j,k)=temvfy(i,j,k)+dtv/PGrid%dy(i,j,k)*temvfz(i,j,k)*    &
                                                      (ve(i,j,k)-ve(i,j-1,k))
                phi(i,j,k)=temlsy(i,j,k)+dtv/PGrid%dy(i,j,k)*temlsz(i,j,k)*    &
                                                      (ve(i,j,k)-ve(i,j-1,k))
                if(vfl(i,j,k)<=vofeps) vfl(i,j,k) = 0.d0
                if(vfl(i,j,k)>=(1.d0-vofeps)) vfl(i,j,k) = 1.d0
              end do
            end do
          end do
        end if
        call Boundary_Condition(vfl)
        call Boundary_Condition(phi)
        call Interface_Reconstruct(PGrid,nx,ny,nz,dis)
        do i = 1,imax
          do j = 1,jmax
            do k = 1,kmax
              PCell%nxL(i,j,k) = nx(i,j,k)
              PCell%nyL(i,j,k) = ny(i,j,k)
              PCell%nzL(i,j,k) = nz(i,j,k)
            end do
          end do
        end do 
        call Redistance_Levelset(PGrid,nx,ny,nz,dis)
      end do
      nullify(vfl)
      nullify(phi)
      nullify(vflF)
      nullify(phiF)
      nullify(nxF)
      nullify(nyF)
      nullify(nzF)
      
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
    subroutine Isinterface(i,j,k,flag)
       implicit none
       real(dp):: esp
       integer i,j,k
       logical flag
       esp = 1.d-14
       flag = .false.
       if(vfl(i,j,k)>esp.and.vfl(i,j,k)<(1.d0-esp)) flag = .true.
       if(dabs(vfl(i-1,j,k)-vfl(i,j,k))>=(1.0d0-2.1d0*esp)) flag = .true.
       if(dabs(vfl(i+1,j,k)-vfl(i,j,k))>=(1.0d0-2.1d0*esp)) flag = .true.
       if(dabs(vfl(i,j-1,k)-vfl(i,j,k))>=(1.0d0-2.1d0*esp)) flag = .true.
       if(dabs(vfl(i,j+1,k)-vfl(i,j,k))>=(1.0d0-2.1d0*esp)) flag = .true.
       if(dabs(vfl(i,j,k-1)-vfl(i,j,k))>=(1.0d0-2.1d0*esp)) flag = .true.
       if(dabs(vfl(i,j,k+1)-vfl(i,j,k))>=(1.0d0-2.1d0*esp)) flag = .true.
       return
    end subroutine Isinterface    
    
    subroutine X_Sweep(PGrid,temvf,temls,ue,ve,we,nxx,nyy,nzz,diss,dtv)
       implicit none
       TYPE(Grid),intent(in)                               :: PGrid
       real(dp),intent(in)				   :: dtv
       real(dp),dimension(:,:,:),intent(in),allocatable    :: ue,ve,we
       real(dp),dimension(:,:,:),intent(inout) 		   :: nxx,nyy,nzz
       real(dp),dimension(:,:,:),intent(inout),allocatable :: diss
       real(dp),dimension(:,:,:),intent(inout),allocatable :: temvf,temls
       integer i,j,k
       real(dp):: flux,lse
       call Interface_Reconstruct(PGrid,nxx,nyy,nzz,diss)
       flux = 0.d0 
       ! volume of fluid 
       do j = 2,jmax-1
         do k = 2,kmax-1
           do i = 2,imax-1
             if(ue(i,j,k)>0.d0) then
               if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps.or.          &
                 (nxx(i,j,k)==0.d0.and.nyy(i,j,k)==0.d0.and.nzz(i,j,k)==0.d0))then
                 flux = vfl(i,j,k)*ue(i,j,k)*dtv/PGrid%dx(i,j,k)
               else
                 call East_Flux(nxx(i,j,k),nyy(i,j,k),nzz(i,j,k),              &
                         diss(i,j,k),vfl(i,j,k),PGrid%dx(i,j,k),               &
                         PGrid%dy(i,j,k),PGrid%dz(i,j,k),ue(i,j,k)*dtv,flux)
               end if
             else
               if(vfl(i+1,j,k)>=(1.d0-vofeps).or.vfl(i+1,j,k)<=vofeps.or.      &
                 (nxx(i+1,j,k)==0.d0.and.nyy(i+1,j,k)==0.d0                    &
                                    .and.nzz(i+1,j,k)==0.d0)) then
                 flux=vfl(i+1,j,k)*ue(i,j,k)*dtv/PGrid%dx(i,j,k)
               else
                 call West_Flux(nxx(i+1,j,k),nyy(i+1,j,k),nzz(i+1,j,k),        &
                           diss(i+1,j,k),vfl(i+1,j,k),PGrid%dx(i,j,k),         &
                           PGrid%dy(i,j,k),PGrid%dz(i,j,k),-ue(i,j,k)*dtv,flux)
                 flux = -flux
               end if
             end if
             if(i==1) then
               flux = vfl(i,j,k)*ue(i,j,k)*dtv/PGrid%dx(i,j,k)
             end if
             if(i>=2) temvf(i,j,k) = temvf(i,j,k)-flux
             if(i<=imax-1) temvf(i+1,j,k) = temvf(i+1,j,k)+flux
           end do
         end do
       end do
       ! level set 
       lse = 0.d0
       flux = 0.d0
       do j = 2,jmax-1
         do i = 2,imax-1
           do k = 2,kmax-1
             if(ue(i,j,k)>=0.d0) then
               lse=phi(i,j,k)+PGrid%dx(i,j,k)/2.d0*(1.d0-ue(i,j,k)*         &
                   dtv/PGrid%dx(i,j,k))*(phi(i+1,j,k)-phi(i-1,j,k))/        &
                   (2.d0*PGrid%dx(i,j,k))
             else
               if(i<=imax-2) then
                 lse=phi(i+1,j,k)-PGrid%dx(i,j,k)/2.d0*(1.d0+ue(i,j,k)*     &
                     dtv/PGrid%dx(i,j,k))*(phi(i+2,j,k)-phi(i,j,k))/        &
                     (2.d0*PGrid%dx(i,j,k))
               else
                 lse=phi(i+1,j,k)-PGrid%dx(i,j,k)/2.d0*(1.d0+ue(i,j,k)*     &
                     dtv/PGrid%dx(i,j,k))*(phi(i+1,j,k)-phi(i,j,k))/        &
                     (PGrid%dx(i,j,k))
               end if
             end if
             flux = ue(i,j,k)*lse*dtv/PGrid%dx(i,j,k)
             if(i>=2) temls(i,j,k) = temls(i,j,k)-flux
             if(i<=imax-1) temls(i+1,j,k) = temls(i+1,j,k)+flux
           end do
         end do
       end do
    end subroutine X_Sweep    
    
    subroutine Y_Sweep(PGrid,temvf,temls,ue,ve,we,nxx,nyy,nzz,diss,dtv)
       implicit none
       TYPE(Grid),intent(in)                               :: PGrid
       real(dp),dimension(:,:,:),intent(in),allocatable    :: ue,ve,we
       real(dp),intent(in)				   :: dtv
       real(dp),dimension(:,:,:),intent(inout),allocatable :: nxx,nyy,nzz,diss
       real(dp),dimension(:,:,:),intent(inout),allocatable :: temvf,temls
       integer :: i,j,k
       real(dp):: flux,lsn 
       call Interface_Reconstruct(PGrid,nxx,nyy,nzz,diss)
       flux = 0.d0 
       ! volume of fluid 
       do i = 2,imax-1
         do j = 2,jmax-1
           do k = 2,kmax-1
             if(ve(i,j,k)>0.d0) then
               if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps.or.          &
                 (nxx(i,j,k)==0.d0.and.nyy(i,j,k)==0.d0.and.nzz(i,j,k)==0.d0)) then
                 flux = vfl(i,j,k)*ve(i,j,k)*dtv/PGrid%dy(i,j,k)
               else
                 call North_Flux(nxx(i,j,k),nyy(i,j,k),nzz(i,j,k),diss(i,j,k), &
                                vfl(i,j,k),PGrid%dx(i,j,k),PGrid%dy(i,j,k),    &
                                PGrid%dz(i,j,k),ve(i,j,k)*dtv,flux)
               end if
             else
               if(vfl(i,j+1,k)>=(1.d0-vofeps).or.vfl(i,j+1,k)<=vofeps.or.      &
                            (nxx(i,j+1,k)==0.d0.and.nyy(i,j+1,k)==0.d0.and.    &
                                                nzz(i,j+1,k)==0.d0)) then
                 flux = vfl(i,j+1,k)*ve(i,j,k)*dtv/PGrid%dy(i,j,k)
               else
                 call South_Flux(nxx(i,j+1,k),nyy(i,j+1,k),nzz(i,j+1,k),       &
                               diss(i,j+1,k),vfl(i,j+1,k),PGrid%dx(i,j,k),     &
                               PGrid%dy(i,j,k),PGrid%dz(i,j,k),-ve(i,j,k)*dtv,flux)
                 flux = -flux
               end if
             end if
             if(j>=2) temvf(i,j,k) = temvf(i,j,k)-flux
             if(j<=jmax-1) temvf(i,j+1,k) = temvf(i,j+1,k)+flux 
           end do
         end do
       end do
       lsn = 0.d0
       flux = 0.d0
    ! level set 
       do i = 2,imax-1
         do j = 2,jmax-1
           do k = 2,kmax-1
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
             if(j<=jmax-1) temls(i,j+1,k) = temls(i,j+1,k)+flux
           end do
         end do
       end do
    end subroutine Y_Sweep
    
    subroutine Z_Sweep(PGrid,temvf,temls,ue,ve,we,nxx,nyy,nzz,diss,dtv)
       implicit none
       type(Grid),intent(in)                               :: PGrid
       real(dp),intent(in)				   :: dtv
       real(dp),dimension(:,:,:),intent(in),allocatable    :: ue,ve,we
       real(dp),dimension(:,:,:),intent(inout),allocatable :: nxx,nyy,nzz,diss
       real(dp),dimension(:,:,:),intent(inout),allocatable :: temvf,temls
       integer 						   :: i,j,k
       real(dp)						   :: flux,lst 
       call Interface_Reconstruct(PGrid,nxx,nyy,nzz,diss)
       flux = 0.d0 
       ! volume of fluid 
       do i = 2,imax-1
         do j = 2,jmax-1
           do k = 2,kmax-1
             if(we(i,j,k)>=0.d0) then
               if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps.or.          &
                 (nxx(i,j,k)==0.d0.and.nyy(i,j,k)==0.d0.and.nzz(i,j,k)==0.d0))then
                 flux = vfl(i,j,k)*we(i,j,k)*dtv/PGrid%dz(i,j,k)
               else
                 call Top_Flux(nxx(i,j,k),nyy(i,j,k),nzz(i,j,k),diss(i,j,k),   &
                               vfl(i,j,k),PGrid%dx(i,j,k),PGrid%dy(i,j,k),     &
                               PGrid%dz(i,j,k),we(i,j,k)*dtv,flux)
               end if
             else
               if(vfl(i,j,k+1)>=(1.d0-vofeps).or.vfl(i,j,k+1)<=vofeps.or.      &
                 (nxx(i,j,k+1)==0.d0.and.nyy(i,j,k+1)==0.d0                    &
                                    .and.nzz(i,j,k+1)==0.d0)) then
                 flux = vfl(i,j,k+1)*we(i,j,k)*dtv/PGrid%dz(i,j,k)
               else
                 call Bottom_Flux(nxx(i,j,k+1),nyy(i,j,k+1),nzz(i,j,k+1),      &
                            diss(i,j,k+1),vfl(i,j,k+1),PGrid%dx(i,j,k),        &
                            PGrid%dy(i,j,k),PGrid%dz(i,j,k),-we(i,j,k)*dtv,flux)
                 flux = -flux
               end if
             end if
             if(k>=2) temvf(i,j,k) = temvf(i,j,k)-flux
             if(k<=kmax-1) temvf(i,j,k+1) = temvf(i,j,k+1)+flux
             if(isnan(flux)) pause 'Vof-scheme 361'
           end do
         end do
       end do
       ! level set 
       do j = 2,jmax-1
         do i = 2,imax-1
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
         end do
       end do
    end subroutine Z_Sweep
    
    subroutine Interface_Reconstruct(PGrid,nxx,nyy,nzz,diss)
       implicit none
       type(Grid),intent(in)                               :: PGrid
       real(dp),dimension(:,:,:),intent(inout)             :: nxx,nyy,nzz
       real(dp),dimension(:,:,:),allocatable,intent(inout) :: diss
       integer						   :: i,j,k
       real(dp)						   :: dx,dy,dz
       real(dp)						   :: nxx1,nyy1,nzz1,diss1
       do i = 2,imax-1
         do j = 2,jmax-1
           do k = 2,kmax-1
             call Interface_Reconstruct_ijk(i,j,k,PGrid%dx(i,j,k),             &
                  PGrid%dy(i,j,k),PGrid%dz(i,j,k),nxx1,nyy1,nzz1,diss1)
             nxx(i,j,k)=nxx1
             nyy(i,j,k)=nyy1
             nzz(i,j,k)=nzz1
             diss(i,j,k)=diss1
           end do
         end do
       end do
       do i = 1,imax
         do k = 1,kmax
           nxx(i,jmax,k)=nxx(i,jmax-1,k)
           nyy(i,jmax,k)=nyy(i,jmax-1,k)
           nzz(i,jmax,k)=nzz(i,jmax-1,k)
           diss(i,jmax,k)=diss(i,jmax-1,k)
           nxx(i,1,k)=nxx(i,2,k)
           nyy(i,1,k)=nyy(i,2,k)
           nzz(i,1,k)=nzz(i,2,k)
           diss(i,1,k)=diss(i,2,k)
         end do
       end do
       do j = 1,jmax
         do k = 1,kmax
           nxx(1,j,k)=nxx(2,j,k)
           nyy(1,j,k)=nyy(2,j,k)
           nzz(1,j,k)=nzz(2,j,k)
           diss(1,j,k)=diss(2,j,k)
           nxx(imax,j,k)=nxx(imax-1,j,k)
           nyy(imax,j,k)=nyy(imax-1,j,k)
           nzz(imax,j,k)=nzz(imax-1,j,k)
           diss(imax,j,k)=diss(imax-1,j,k)
         end do
       end do
       do i = 1,imax
         do j = 1,jmax
           nxx(i,j,1)=nxx(i,j,2)
           nyy(i,j,1)=nyy(i,j,2)
           nzz(i,j,1)=nzz(i,j,2)
           diss(i,j,1)=diss(i,j,2)
           nxx(i,j,kmax)=nxx(i,j,kmax-1)
           nyy(i,j,kmax)=nyy(i,j,kmax-1)
           nzz(i,j,kmax)=nzz(i,j,kmax-1)
           diss(i,j,kmax)=diss(i,j,kmax-1)
         end do
       end do
    end subroutine Interface_Reconstruct
    
    subroutine Interface_Reconstruct_ijk(i,j,k,dx,dy,dz,nxx,nyy,nzz,diss) 
       implicit none
       integer:: i,j,k,ii,jj,kk
       real(dp):: dx,dy,dz,nxx,nyy,nzz,diss,temp,vofeps
       real(dp):: nx(3),ny(3),nz(3),delt
       real(dp):: dxp,dxn,dyp,dyn,dzp,dzn,nxx1,nyy1,nzz1,nxy,nyz,nxz
       logical:: flag
       vofeps = 1.d-14
       delt = 2.d0
       nxx = 0.d0
       nyy = 0.d0
       nzz = 0.d0
       diss = 0.d0
       call Isinterface(i,j,k,flag)
       if(flag.eqv..true.) then
          call Normal_Vector_Irre(i,j,k,dx,dy,dz,nxx,nyy,nzz)
          temp = dsqrt(nxx**2.d0+nyy**2.d0+nzz**2.d0)
          if(isnan(temp)) then
            pause 'vof 903'
          end if
          if(temp<1.d-14) then
             nxx = 0.d0
             nyy = 0.d0
             nzz = 1.d0
          else
             nxx = nxx/temp
             nyy = nyy/temp
             nzz = nzz/temp
          end if
          call Find_Distance(dx,dy,dz,nxx,nyy,nzz,vfl(i,j,k),diss)         
          if(isnan(diss)) then
             pause 'reconstruct 245'
          end if
       end if
    end subroutine Interface_Reconstruct_ijk
    
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
       dy1 = (nx1*dz+ny1*dy+nz1*dx)-dx1-dz1
       fc = dmin1(f,1.d0-f)
       sm = dx1+dy1+dz1
       sc = (6.d0*fc*dx1*dy1*dz1)**(1.d0/3.d0)
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
       real(dp),intent(in):: nxx,nyy,nzz,dx,dy,dz,s
       real(dp),intent(out):: vol
       real(dp):: nx1,ny1,nz1,dx1,dy1,dz1
       real(dp):: sc,sm,fc,vol1,vol2
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
       if(s<=0.5d0*sm) then                 ! be careful with this step. the chosen of 
          vol = fc*dx*dy*dz          ! the region contains fluid
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
       eps = 1.d-14
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
       eps = 1.d-14
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
       eps = 1.d-14
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
       real(dp):: vol,eps = 1.d-14
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
       eps = 1.d-14
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
       eps = 1.d-14
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
    subroutine Redistance_Levelset(PGrid,Tnx,Tny,Tnz,Tdis)
      implicit none
      type(Grid),intent(in)				    :: PGrid
      real(kind=dp),dimension(:,:,:),allocatable,intent(in) :: Tnx,Tny,Tnz,Tdis
      integer						    :: i,j,k,ii,jj,kk
      integer						    :: iii,jjj,kkk,ij
      integer						    :: l,m,n,cas
      real(kind=dp)				  	    :: xv,yv,zv,dv
      real(kind=dp)					    :: deuc,dij1
      real(kind=dp)					    :: xp,yp,zp
      real(kind=dp)					    :: xoff,yoff,zoff
      real(kind=dp) 					    :: xfc,yfc,zfc
      real(kind=dp)					    :: xs,ys,zs
      real(kind=dp)					    :: diss1,diss,nxx,nyy,nzz
      logical						    :: pointp,face
      real(kind=dp),dimension(:,:,:),allocatable	    :: phiaux
      integer,dimension(:,:,:),allocatable		    :: fix
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
                                diss1=nxx*xs+diss  ! diss for dis1 and diss1 for dis
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
          case(2) ! in y surface
             point1=-0.5d0*dx*nxx+face_point*nyy-0.5d0*dz*nzz+diss
             point2=-0.5d0*dx*nxx+face_point*nyy+0.5d0*dz*nzz+diss
             point3=0.5d0*dx*nxx+face_point*nyy+0.5d0*dz*nzz+diss
             point4=0.5d0*dx*nxx+face_point*nyy-0.5d0*dz*nzz+diss
          case(3) ! in z surface
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
       integer,intent(in)     :: ii,jj,kk,cas
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
               zs = (diss+nyy*ys)/(-nzz)
               return 
             else
               zs = zfc
               ys = (diss+nzz*zs)/(-nyy)
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
               zs = (diss+nxx*xs)/(-nzz)
               return 
             else
               zs = zfc
               xs = (diss+nzz*zs)/(-nxx)
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
               ys=(diss+nxx*xs)/(-nyy)
               return 
             else
               ys=yfc
               xs=(diss+nyy*ys)/(-nxx)
               return
             end if
           end if
         end select 
    end subroutine Shortest_Point_2d
    
    subroutine Boundary_Condition(vari)
       implicit none
       integer i,j,k
       real(dp),dimension(:,:,:),intent(inout) :: vari
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
    
    subroutine Normal_Vector_Irre(i,j,k,delx,dely,delz,nxx,nyy,nzz)
       implicit none
       integer,intent(in)   :: i,j,k
       real(dp),intent(in)  :: delx,dely,delz
       real(dp),intent(out) :: nxx,nyy,nzz
       integer              :: case_deri,dx,dy,dz
       real(dp)             :: qi(-1:1),qj(-1:1),qk(-1:1)
       real(dp)             :: vx,vy,vz
       ! define qi for Dx
       if(i>=3) then
       ! for qi(-1)
          vx = 0.5d0*(phi(i,j,k)-phi(i-2,j,k))/delx
       else
          vx = (phi(i,j,k)-phi(i-1,j,k))/delx
       end if
       vy = 0.5d0*(phi(i-1,j+1,k)-phi(i-1,j-1,k))/dely
       vz = 0.5d0*(phi(i-1,j,k+1)-phi(i-1,j,k-1))/delz
       qi(-1) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! for qi(1) 
       if(i<=imax-2) then
          vx = 0.5d0*(phi(i+2,j,k)-phi(i,j,k))/delx
       else
          vx = (phi(i+1,j,k)-phi(i,j,k))/delx
       end if
       vy = 0.5d0*(phi(i+1,j+1,k)-phi(i+1,j-1,k))/dely
       vz = 0.5d0*(phi(i+1,j,k+1)-phi(i+1,j,k-1))/delz
       qi(1) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0)) 
       if(qi(-1)<eta.and.qi(1)>=eta) then
          dx = -1
       elseif(qi(-1)>=eta.and.qi(1)<eta) then
          dx = 1
       else
          dx = 0
       end if
       ! define qj for Dy
       if(j>=3) then
          vy = 0.5d0*(phi(i,j,k)-phi(i,j-2,k))/dely
       else
          vy = (phi(i,j,k)-phi(i,j-1,k))/dely
       end if
       vx = 0.5d0*(phi(i+1,j-1,k)-phi(i-1,j-1,k))/delx
       vz = 0.5d0*(phi(i,j-1,k+1)-phi(i,j-1,k-1))/delz
       qj(-1) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! define qj(1)
       vx = 0.5d0*(phi(i+1,j+1,k)-phi(i-1,j+1,k))/delx
       if(j<=jmax-2) then
          vy = 0.5d0*(phi(i,j+2,k)-phi(i,j,k))/dely
       else
          vy = (phi(i,j+1,k)-phi(i,j,k))/dely
       end if
       vz = 0.5d0*(phi(i,j+1,k+1)-phi(i,j+1,k-1))/delz
       qj(1) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       if(qj(-1)<eta.and.qj(1)>=eta) then
          dy = -1
       elseif(qj(-1)>=eta.and.qj(1)<eta) then
          dy = 1
       else
          dy = 0
       end if
       ! define qk for dz
       ! qk(-1)
       vx = 0.5d0*(phi(i+1,j,k-1)-phi(i-1,j,k-1))/delx
       vy = 0.5d0*(phi(i,j+1,k-1)-phi(i,j-1,k-1))/dely
       if(k>=3) then
          vz = 0.5d0*(phi(i,j,k)-phi(i,j,k-2))/delz
       else
          vz = (phi(i,j,k)-phi(i,j,k-1))/delz
       end if
       qk(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! qk(1)
       vx = 0.5d0*(phi(i+1,j,k+1)-phi(i-1,j,k+1))/delx
       vy = 0.5d0*(phi(i,j+1,k+1)-phi(i,j-1,k+1))/dely
       if(k<=kmax-2) then
          vz = 0.5d0*(phi(i,j,k+2)-phi(i,j,k))/delz
       else
          vz = (phi(i,j,k+1)-phi(i,j,k))/delz
       end if
       qk(1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       if(qk(-1)<eta.and.qk(1)>=eta) then
          dz = -1
       elseif(qk(-1)>=eta.and.qk(1)<eta) then
          dz = 1
       else
          dz = 0
       end if
       if(dx==-1) nxx = (phi(i,j,k)-phi(i-1,j,k))/delx
       if(dx==1) nxx = (phi(i+1,j,k)-phi(i,j,k))/delx
       if(dx==0) nxx = 0.5d0*(phi(i+1,j,k)-phi(i-1,j,k))/delx
       if(dy==-1) nyy = (phi(i,j,k)-phi(i,j-1,k))/dely
       if(dy==1) nyy = (phi(i,j+1,k)-phi(i,j,k))/dely
       if(dy==0) nyy = 0.5d0*(phi(i,j+1,k)-phi(i,j-1,k))/dely
       if(dz==-1) nzz = (phi(i,j,k)-phi(i,j,k-1))/delz
       if(dz==1) nzz = (phi(i,j,k+1)-phi(i,j,k))/delz
       if(dz==0) nzz = 0.5d0*(phi(i,j,k+1)-phi(i,j,k-1))/delz
    end subroutine Normal_Vector_Irre
end module Clsvof

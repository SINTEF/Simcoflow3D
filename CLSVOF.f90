module Clsvof
    use ieee_arithmetic
    use variables
    implicit none
    integer:: nv
    integer,parameter:: band_width = 4
    real(dp),dimension(:,:,:),allocatable:: vfl
    real(dp),dimension(:,:,:),allocatable:: phi
    real(dp),parameter:: dzero = 0.d0,done = 1.d0,vofeps=1.d-14
    real(dp):: dtv,eta,epsilon
    contains
     
    subroutine Initial_Clsvof(Intv,Iinfo,Gri)
       implicit none
       type(Interface_Var),dimension(:,:,:),intent(inout),allocatable:: Intv
       type(Grid_Var),dimension(:,:,:),intent(in),allocatable:: Gri
       type(Interface_Information),dimension(:,:,:),intent(inout),allocatable:: Iinfo
       integer:: i,j,k
       real(dp):: nx,ny,nz,vol,dis,dx,dy,dz,s
       do i = 1,imax
          do j = 1,jmax 
             do k = 1,kmax
                dx = Gri(i,j,k)%x-0.35d0
                dy = Gri(i,j,k)%y-0.35d0  
                dz = Gri(i,j,k)%z-0.35d0
                dis = dsqrt(dx**2.d0+dy**2.d0+dz**2.d0)-ra
                nx = dx/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
                ny = dy/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
                nz = dz/dsqrt(dx**2.d0+dy**2.d0+dz**2.d0+epsi)
                s = dis+0.5*(dabs(nx)*del_x+dabs(ny)*del_y+dabs(nz)*del_z)
                Iinfo(i,j,k)%nx = 0.d0  
                Iinfo(i,j,k)%ny = 0.d0
                Iinfo(i,j,k)%nz = 0.d0
                Iinfo(i,j,k)%dis = 0.d0
                call Volume_Fraction_Calc(del_x,del_y,del_z,nx,ny,nz,s,vol)
                Intv(i,j,k)%vfl = 1.d0-vol/(del_x*del_y*del_z)
                Intv(i,j,k)%phi = dis 
             end do
          end do
       end do
    end subroutine   
!***********************************************************    
! coupled level set and volume of fluid 
!***********************************************************
    subroutine Clsvof_Scheme(Intv,Iinfo,Grid)
       implicit none
       type(Interface_Var),dimension(:,:,:),allocatable:: Intv
       type(Interface_Information),dimension(:,:,:),allocatable:: Iinfo
       type(Grid_Var),dimension(:,:,:),allocatable:: Grid
       integer i,j,k,kk
       real(dp),dimension(:,:,:),allocatable:: u,v,w,ue,ve,we,nx,ny,nz,dis
       real(dp),dimension(:,:,:),allocatable:: temvfx,temvfy,temvfz
       real(dp),dimension(:,:,:),allocatable:: temlsx,temlsy,temlsz
       allocate(vfl(imax,jmax,kmax))
       allocate(phi(imax,jmax,kmax))
       allocate(nx(imax,jmax,kmax))
       allocate(ny(imax,jmax,kmax))
       allocate(nz(imax,jmax,kmax))
       allocate(u(imax,jmax,kmax))
       allocate(v(imax,jmax,kmax))
       allocate(w(imax,jmax,kmax))
       allocate(ue(imax,jmax,kmax))
       allocate(ve(imax,jmax,kmax))
       allocate(we(imax,jmax,kmax))
       allocate(dis(imax,jmax,kmax))
       allocate(temvfx(imax,jmax,kmax))
       allocate(temvfy(imax,jmax,kmax))
       allocate(temvfz(imax,jmax,kmax))
       allocate(temlsx(imax,jmax,kmax))
       allocate(temlsy(imax,jmax,kmax))
       allocate(temlsz(imax,jmax,kmax))
       vfl(:,:,:) = Intv(:,:,:)%vfl
       phi(:,:,:) = Intv(:,:,:)%phi
       dtv = dt/dble(nv)
       do kk = 0,nv-1
          do i = 1,imax
             do j = 1,jmax
                do k = 1,kmax
                   u(i,j,k) = 2.d0*(dsin(pi*Grid(i,j,k)%x))**2.d0*dsin(2.d0*pi*     &
                                  Grid(i,j,k)%y)*dsin(2.d0*pi*Grid(i,j,k)%z)*       &
                                                 dcos(pi*(t+dtv*dble(kk))/tp) 
                   v(i,j,k) = -dsin(2.d0*pi*Grid(i,j,k)%x)*(dsin(pi*                &
                                  Grid(i,j,k)%y))**2.d0*dsin(2.d0*pi*Grid(i,j,k)%z)*&
                                                 dcos(pi*(t+dtv*dble(kk))/tp)
                   w(i,j,k) = -dsin(2.d0*pi*Grid(i,j,k)%x)*dsin(2.d0*pi*            &
                                   Grid(i,j,k)%y)*(dsin(pi*Grid(i,j,k)%z))**2.d0*   &
                                                 dcos(pi*(t+dtv*dble(kk))/tp)
                end do
             end do
          end do  
          do i = 1,imax-1
             do j = 1,jmax
                do k = 1,kmax
                   ue(i,j,k) = 0.5d0*(u(i,j,k)+u(i+1,j,k))
                end do
             end do
          end do
          ue(imax,:,:) = ue(imax-1,:,:)
          do j = 1,jmax-1
             do i = 1,imax
                do k = 1,kmax
                   ve(i,j,k) = 0.5d0*(v(i,j,k)+v(i,j+1,k))
                end do
             end do
          end do
          ve(:,jmax,:) = ve(:,jmax-1,:)
          do k = 1,kmax-1
             do i = 1,imax
                do j = 1,jmax
                   we(i,j,k) = 0.5d0*(w(i,j,k)+w(i,j,k+1))
                end do
             end do
          end do
          we(:,:,kmax) = we(:,:,kmax-1)
          if(mod(kk,3)==0) then
             temvfx(:,:,:) = vfl(:,:,:)
             temlsx(:,:,:) = phi(:,:,:)
             call X_Sweep(temvfx,temlsx,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
               do j = 2,jmax-1
                  do k = 2,kmax-1
                     temvfx(i,j,k) = temvfx(i,j,k)/(1.d0-dtv/del_x*(ue(i,j,k)-      &
                                                                    ue(i-1,j,k)))
                     temlsx(i,j,k) = temlsx(i,j,k)/(1.d0-dtv/del_x*(ue(i,j,k)-      &
                                                                    ue(i-1,j,k))) 
                     if(temvfx(i,j,k)<=vofeps.or.temlsx(i,j,k)>=epsilon)            &
                                                            temvfx(i,j,k) = 0.d0
                     if(temvfx(i,j,k)>=(1.d0-vofeps).or.temlsx(i,j,k)<=-epsilon)    &
                                                            temvfx(i,j,k) = 1.d0
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
             call Y_Sweep(temvfy,temlsy,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
                do j = 2,jmax-1
                   do k = 2,kmax-1
                      temvfy(i,j,k) = temvfy(i,j,k) + dtv/del_y*temvfx(i,j,k)*      &
                                                      (ve(i,j,k)-ve(i,j-1,k))
                      temlsy(i,j,k) = temlsy(i,j,k) + dtv/del_y*temlsx(i,j,k)*      &
                                                      (ve(i,j,k)-ve(i,j-1,k))
                      if(temvfy(i,j,k)<=vofeps.or.temlsy(i,j,k)>=epsilon)           &
                                                         temvfy(i,j,k) = 0.d0
                      if(temvfy(i,j,k)>=(1.d0-vofeps).or.temlsy(i,j,k)<=-epsilon)   &
                                                         temvfy(i,j,k) = 1.d0
                      vfl(i,j,k) = temvfy(i,j,k)
                      phi(i,j,k) = temlsy(i,j,k)
                   end do
                end do
             end do
             call Boundary_Condition(vfl)
             call Boundary_Condition(phi)
             temvfz(:,:,:) = vfl(:,:,:)
             temlsz(:,:,:) = phi(:,:,:)
             call Z_Sweep(temvfz,temlsz,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
                do j = 2,jmax-1
                   do k = 2,kmax-1
                      vfl(i,j,k) = temvfz(i,j,k) + dtv/del_z*temvfx(i,j,k)*         &
                                                      (we(i,j,k)-we(i,j,k-1))
                      phi(i,j,k) = temlsz(i,j,k) + dtv/del_z*temlsx(i,j,k)*         &
                                                      (we(i,j,k)-we(i,j,k-1))
                      if(vfl(i,j,k)<=vofeps.or.temlsz(i,j,k)>=epsilon)           &
                                                             vfl(i,j,k) = 0.d0
                      if(vfl(i,j,k)>=(1.d0-vofeps).or.temlsz(i,j,k)<=-epsilon)   &
                                                             vfl(i,j,k) = 1.d0
                   end do
                end do
             end do
          elseif(mod(kk,3)==1) then
             temvfy(:,:,:) = vfl(:,:,:)
             temlsy(:,:,:) = phi(:,:,:)
             call Y_Sweep(temvfy,temlsy,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
               do j = 2,jmax-1
                  do k = 2,kmax-1
                     temvfy(i,j,k) = temvfy(i,j,k)/(1.d0-dtv/del_y*(ve(i,j,k)-      &
                                                                    ve(i,j-1,k)))
                     temlsy(i,j,k) = temlsy(i,j,k)/(1.d0-dtv/del_y*(ve(i,j,k)-      &
                                                                    ve(i,j-1,k))) 
                     if(temvfy(i,j,k)<=vofeps.or.temlsy(i,j,k)>=epsilon)            &
                                                             temvfy(i,j,k) = 0.d0
                     if(temvfy(i,j,k)>=(1.d0-vofeps).or.temlsy(i,j,k)<=-epsilon)    &
                                                             temvfy(i,j,k) = 1.d0
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
             call Z_Sweep(temvfz,temlsz,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
                do j = 2,jmax-1
                   do k = 2,kmax-1
                      temvfz(i,j,k) = temvfz(i,j,k) + dtv/del_z*temvfy(i,j,k)*     &
                                                          (we(i,j,k)-we(i,j,k-1))
                      temlsz(i,j,k) = temlsz(i,j,k) + dtv/del_z*temlsy(i,j,k)*     &
                                                          (we(i,j,k)-we(i,j,k-1))
                      if(temvfz(i,j,k)<=vofeps.or.temlsz(i,j,k)>=epsilon)          &
                                                            temvfz(i,j,k)=0.d0
                      if(temvfz(i,j,k)>=(1.d0-vofeps).or.temlsz(i,j,k)<=-epsilon)  &
                                                            temvfz(i,j,k)=1.d0
                      vfl(i,j,k) = temvfz(i,j,k)
                      phi(i,j,k) = temlsz(i,j,k)
                   end do
                end do
             end do
             call Boundary_Condition(vfl)
             call Boundary_Condition(phi)
             temvfx(:,:,:) = vfl(:,:,:)
             temlsx(:,:,:) = phi(:,:,:)
             call X_Sweep(temvfx,temlsx,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
                do j = 2,jmax-1
                   do k = 2,kmax-1
                      vfl(i,j,k) = temvfx(i,j,k) + dtv/del_x*temvfy(i,j,k)*         &
                                                         (ue(i,j,k)-ue(i-1,j,k))
                      phi(i,j,k) = temlsx(i,j,k) + dtv/del_x*temlsy(i,j,k)*         &
                                                         (ue(i,j,k)-ue(i-1,j,k))
                      if(vfl(i,j,k)<=vofeps.or.temlsx(i,j,k)>=epsilon)           &
                                                             vfl(i,j,k)=0.d0
                      if(vfl(i,j,k)>=(1.d0-vofeps).or.temlsx(i,j,k)<=-epsilon)   &
                                                             vfl(i,j,k) = 1.d0
                   end do
                end do
             end do
          else 
             temvfz(:,:,:) = vfl(:,:,:)
             temlsz(:,:,:) = phi(:,:,:)
             call Z_Sweep(temvfz,temlsz,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
                do j = 2,jmax-1
                   do k = 2,kmax-1
                      temvfz(i,j,k) = temvfz(i,j,k)/(1.d0-dtv/del_z*               &
                                                   (we(i,j,k)-we(i,j,k-1)))
                      temlsz(i,j,k) = temlsz(i,j,k)/(1.d0-dtv/del_z*               &
                                                   (we(i,j,k)-we(i,j,k-1)))
                      if(temvfz(i,j,k)<=vofeps.or.temlsz(i,j,k)>=epsilon)          &
                                                    temvfz(i,j,k)=0.d0
                      if(temvfz(i,j,k)>=(1.d0-vofeps).or.temlsz(i,j,k)<=-epsilon)  &
                                                    temvfz(i,j,k)=1.d0
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
             call X_Sweep(temvfx,temlsx,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
                do j = 2,jmax-1
                   do k = 2,kmax-1
                      temvfx(i,j,k) = temvfx(i,j,k) + dtv/del_x*temvfz(i,j,k)*      &
                                                         (ue(i,j,k)-ue(i-1,j,k))
                      temlsx(i,j,k) = temlsx(i,j,k) + dtv/del_x*temlsz(i,j,k)*      &
                                                         (ue(i,j,k)-ue(i-1,j,k))
                      if(temvfx(i,j,k)<=vofeps.or.temlsx(i,j,k)>=epsilon)           &
                                                             temvfx(i,j,k)=0.d0
                      if(temvfx(i,j,k)>=(1.d0-vofeps).or.temlsx(i,j,k)<=-epsilon)   &
                                                             temvfx(i,j,k) = 1.d0 
                      vfl(i,j,k) = temvfx(i,j,k)
                      phi(i,j,k) = temlsx(i,j,k)
                   end do
                end do
             end do
             call Boundary_Condition(vfl)
             call Boundary_Condition(phi)  
             temvfy(:,:,:) = vfl(:,:,:)
             temlsy(:,:,:) = phi(:,:,:)
             call Y_Sweep(temvfy,temlsy,ue,ve,we,nx,ny,nz,dis)
             do i = 2,imax-1
                do j = 2,jmax-1
                   do k = 2,kmax-1
                      vfl(i,j,k) = temvfy(i,j,k) + dtv/del_y*temvfz(i,j,k)*         &
                                                      (ve(i,j,k)-ve(i,j-1,k))
                      phi(i,j,k) = temlsy(i,j,k) + dtv/del_y*temlsz(i,j,k)*         &
                                                      (ve(i,j,k)-ve(i,j-1,k))
                      if(vfl(i,j,k)<=vofeps.or.temlsy(i,j,k)>=epsilon)           &
                                                        vfl(i,j,k) = 0.d0
                      if(vfl(i,j,k)>=(1.d0-vofeps).or.temlsy(i,j,k)<=-epsilon)   &
                                                        vfl(i,j,k) = 1.d0
                   end do
                end do
             end do
          end if
          call Boundary_Condition(vfl)
          call Boundary_Condition(phi)
          call Interface_Reconstruct(del_x,del_y,del_z,nx,ny,nz,dis)
          do i = 1,imax
             do j = 1,jmax
                do k = 1,kmax
                   Iinfo(i,j,k)%nx = nx(i,j,k)
                   Iinfo(i,j,k)%ny = ny(i,j,k)
                   Iinfo(i,j,k)%nz = nz(i,j,k)
                   Iinfo(i,j,k)%dis = (0.5d0*del_x*dabs(nx(i,j,k))+0.5d0*del_y*       &
                               dabs(ny(i,j,k))+0.5d0*del_z*dabs(nz(i,j,k)))-dis(i,j,k) 
                end do
             end do
          end do 
          call Redistance_Levelset(Iinfo)
       end do
       Intv(:,:,:)%vfl = vfl(:,:,:)
       Intv(:,:,:)%phi = phi(:,:,:)
       deallocate(vfl)
       deallocate(phi)
       deallocate(nx)
       deallocate(ny)
       deallocate(nz)
       deallocate(dis)
       deallocate(u)
       deallocate(v)
       deallocate(w)
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
    subroutine X_Sweep(temvf,temls,ue,ve,we,nx,ny,nz,dis)
       implicit none
       real(dp),dimension(:,:,:),intent(in),allocatable:: ue,ve,we
       real(dp),dimension(:,:,:),intent(inout),allocatable:: nx,ny,nz,dis
       real(dp),dimension(:,:,:),intent(inout),allocatable:: temvf,temls
       integer i,j,k
       real(dp):: flux,lse
       call Interface_Reconstruct(del_x,del_y,del_z,nx,ny,nz,dis)
       flux = 0.d0 
       ! volume of fluid 
       do j = 2,jmax-1
          do k = 2,kmax-1
             do i = 2,imax-1
                if(ue(i,j,k)>0.d0) then
                   if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps.or.           &
                    (nx(i,j,k)==0.d0.and.ny(i,j,k)==0.d0.and.nz(i,j,k)==0.d0))then
                      flux = vfl(i,j,k)*ue(i,j,k)*dtv/del_x
                   else
                      call East_Flux(nx(i,j,k),ny(i,j,k),nz(i,j,k),dis(i,j,k),      &
                                vfl(i,j,k),del_x,del_y,del_z,ue(i,j,k)*dtv,flux)
                   end if
                else
                   if(vfl(i+1,j,k)>=(1.d0-vofeps).or.vfl(i+1,j,k)<=vofeps.or.       &
                     (nx(i+1,j,k)==0.d0.and.ny(i+1,j,k)==0.d0                       &
                                       .and.nz(i+1,j,k)==0.d0)) then
                      flux = vfl(i+1,j,k)*ue(i,j,k)*dtv/del_x
                   else
                      call West_Flux(nx(i+1,j,k),ny(i+1,j,k),nz(i+1,j,k),           &
                                     dis(i+1,j,k),vfl(i+1,j,k),del_x,del_y,del_z,   &
                                               -ue(i,j,k)*dtv,flux)
                      flux = -flux
                   end if
                end if
                if(i==1) then
                   flux = vfl(i,j,k)*ue(i,j,k)*dtv/del_x
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
                   lse =  phi(i,j,k)+del_x/2.d0*(1.d0-ue(i,j,k)*dtv/del_x)*         &
                          (phi(i+1,j,k)-phi(i-1,j,k))/(2.d0*del_x)
                else
                   if(i<=imax-2) then
                      lse = phi(i+1,j,k)-del_x/2.d0*(1.d0+ue(i,j,k)*dtv/del_x)*     &
                            (phi(i+2,j,k)-phi(i,j,k))/(2.d0*del_x)
                   else
                      lse = phi(i+1,j,k)-del_x/2.d0*(1.d0+ue(i,j,k)*dtv/del_x)*     &
                            (phi(i+1,j,k)-phi(i,j,k))/(del_x)
                   end if
                end if
                flux = ue(i,j,k)*lse*dtv/del_x
                if(i>=2) temls(i,j,k) = temls(i,j,k)-flux
                if(i<=imax-1) temls(i+1,j,k) = temls(i+1,j,k)+flux
             end do
          end do
       end do
    end subroutine X_Sweep    
    subroutine Y_Sweep(temvf,temls,ue,ve,we,nx,ny,nz,dis)
       implicit none
       real(dp),dimension(:,:,:),intent(in),allocatable:: ue,ve,we
       real(dp),dimension(:,:,:),intent(inout),allocatable:: nx,ny,nz,dis
       real(dp),dimension(:,:,:),intent(inout),allocatable:: temvf,temls
       integer i,j,k
       real(dp):: flux,lsn 
       call Interface_Reconstruct(del_x,del_y,del_z,nx,ny,nz,dis)
       flux = 0.d0 
       ! volume of fluid 
       do i = 2,imax-1
          do j = 2,jmax-1
             do k = 2,kmax-1
                if(ve(i,j,k)>0.d0) then
                   if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps.or.           &
                     (nx(i,j,k)==0.d0.and.ny(i,j,k)==0.d0.and.nz(i,j,k)==0.d0)) then
                      flux = vfl(i,j,k)*ve(i,j,k)*dtv/del_y
                   else
                      call North_Flux(nx(i,j,k),ny(i,j,k),nz(i,j,k),dis(i,j,k),     &
                               vfl(i,j,k),del_x,del_y,del_z,ve(i,j,k)*dtv,flux)
                   end if
                else
                   if(vfl(i,j+1,k)>=(1.d0-vofeps).or.vfl(i,j+1,k)<=vofeps.or.       &
                            (nx(i,j+1,k)==0.d0.and.ny(i,j+1,k)==0.d0.and.           &
                                                nz(i,j+1,k)==0.d0)) then
                   flux = vfl(i,j+1,k)*ve(i,j,k)*dtv/del_y
                else
                   call South_Flux(nx(i,j+1,k),ny(i,j+1,k),nz(i,j+1,k),             &
                                      dis(i,j+1,k),vfl(i,j+1,k),del_x,del_y,del_z,  &
                                            -ve(i,j,k)*dtv,flux)
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
                   lsn =  phi(i,j,k)+del_y/2.d0*(1.d0-ve(i,j,k)*dtv/del_y)*         &
                          (phi(i,j+1,k)-phi(i,j-1,k))/(2.d0*del_y)
                else
                   if(j<=jmax-2) then
                      lsn = phi(i,j+1,k)-del_y/2.d0*(1.d0+ve(i,j,k)*dtv/del_y)*     &
                            (phi(i,j+2,k)-phi(i,j,k))/(2.d0*del_y)
                   else
                      lsn = phi(i,j+1,k)-del_y/2.d0*(1.d0+ve(i,j,k)*dtv/del_y)*     &
                            (phi(i,j+1,k)-phi(i,j,k))/del_y
                   end if
                end if
                flux = lsn*ve(i,j,k)*dt/del_y
                if(j>=2) temls(i,j,k) = temls(i,j,k)-flux
                if(j<=jmax-1) temls(i,j+1,k) = temls(i,j+1,k)+flux
             end do
          end do
       end do
    end subroutine Y_Sweep
    
    subroutine Z_Sweep(temvf,temls,ue,ve,we,nx,ny,nz,dis)
       implicit none
       real(dp),dimension(:,:,:),intent(in),allocatable:: ue,ve,we
       real(dp),dimension(:,:,:),intent(inout),allocatable:: nx,ny,nz,dis
       real(dp),dimension(:,:,:),intent(inout),allocatable:: temvf,temls
       integer i,j,k
       real(dp):: flux,lst 
       call Interface_Reconstruct(del_x,del_y,del_z,nx,ny,nz,dis)
       flux = 0.d0 
       ! volume of fluid 
       do i = 2,imax-1
          do j = 2,jmax-1
             do k = 2,kmax-1
                if(we(i,j,k)>=0.d0) then
                   if(vfl(i,j,k)>=(1.d0-vofeps).or.vfl(i,j,k)<=vofeps.or.           &
                     (nx(i,j,k)==0.d0.and.ny(i,j,k)==0.d0.and.nz(i,j,k)==0.d0))then
                      flux = vfl(i,j,k)*we(i,j,k)*dtv/del_z
                   else
                      call Top_Flux(nx(i,j,k),ny(i,j,k),nz(i,j,k),dis(i,j,k),       &
                                    vfl(i,j,k),del_x,del_y,del_z,we(i,j,k)*dtv,flux)
                   end if
                else
                   if(vfl(i,j,k+1)>=(1.d0-vofeps).or.vfl(i,j,k+1)<=vofeps.or.       &
                     (nx(i,j,k+1)==0.d0.and.ny(i,j,k+1)==0.d0                       &
                                       .and.nz(i,j,k+1)==0.d0)) then
                      flux = vfl(i,j,k+1)*we(i,j,k)*dtv/del_z
                   else
                      call Bottom_Flux(nx(i,j,k+1),ny(i,j,k+1),nz(i,j,k+1),         &
                                       dis(i,j,k+1),vfl(i,j,k+1),del_x,del_y,del_z, &
                                       -we(i,j,k)*dtv,flux)
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
                   lst =  phi(i,j,k)+del_z/2.d0*(1.d0-we(i,j,k)*dtv/del_z)*         &
                          (phi(i,j,k+1)-phi(i,j,k-1))/(2.d0*del_z)
                else
                   if(k<=kmax-2) then
                      lst = phi(i,j,k+1)-del_z/2.d0*(1.d0+we(i,j,k)*dtv/del_z)*     &
                            (phi(i,j,k+2)-phi(i,j,k))/(2.d0*del_z)
                   else
                      lst = phi(i,j,k+1)-del_z/2.d0*(1.d0+we(i,j,k)*dtv/del_z)*     &
                            (phi(i,j,k+1)-phi(i,j,k))/del_z
                   end if
                end if
                flux = lst*we(i,j,k)*dt/del_z
                if(k>=2) temls(i,j,k) = temls(i,j,k)-flux
                if(k<=kmax-1) temls(i,j,k+1) = temls(i,j,k+1)+flux
             end do
          end do
       end do
    end subroutine Z_Sweep
    
    subroutine Interface_Reconstruct(dx,dy,dz,nx,ny,nz,dis)
       implicit none
       real(dp),dimension(:,:,:),intent(inout),allocatable:: nx,ny,nz,dis
       integer:: i,j,k
       real(dp):: dx,dy,dz,nxx,nyy,nzz,diss
       do i = 2,imax-1
          do j = 2,jmax-1
             do k = 2,kmax-1
                call Interface_Reconstruct_ijk(i,j,k,dx,dy,dz,nxx,nyy,nzz,diss)
                nx(i,j,k) = nxx
                ny(i,j,k) = nyy
                nz(i,j,k) = nzz
                dis(i,j,k) = diss
             end do
          end do
       end do
       do i = 1,imax
          do k = 1,kmax
             nx(i,jmax,k) = nx(i,jmax-1,k)
             ny(i,jmax,k) = ny(i,jmax-1,k)
             nz(i,jmax,k) = ny(i,jmax-1,k)
             dis(i,jmax,k) = dis(i,jmax-1,k)
             nx(i,1,k) = nx(i,2,k)
             ny(i,1,k) = ny(i,2,k)
             nz(i,1,k) = nz(i,2,k)
             dis(i,1,k) = dis(i,2,k)
          end do
       end do
       do j = 1,jmax
          do k = 1,kmax
             nx(1,j,k) = nx(2,j,k)
             ny(1,j,k) = ny(2,j,k)
             nz(1,j,k) = nz(2,j,k)
             dis(1,j,k) = dis(2,j,k)
             nx(imax,j,k) = nx(imax-1,j,k)
             ny(imax,j,k) = ny(imax-1,j,k)
             nz(imax,j,k) = nz(imax-1,j,k)
             dis(imax,j,k) = dis(imax-1,j,k)
          end do
       end do
       do i = 1,imax
          do j = 1,jmax
             nx(i,j,1) = nx(i,j,2)
             ny(i,j,1) = ny(i,j,2)
             nz(i,j,1) = nz(i,j,2)
             dis(i,j,1) = dis(i,j,2)
             nx(i,j,kmax) = nx(i,j,kmax-1)
             ny(i,j,kmax) = ny(i,j,kmax-1)
             nz(i,j,kmax) = nz(i,j,kmax-1)
             dis(i,j,kmax) = dis(i,j,kmax-1)
          end do
       end do
    end subroutine Interface_Reconstruct
    
    subroutine Interface_Reconstruct_ijk(i,j,k,dx,dy,dz,nxx,nyy,nzz,diss) 
       implicit none
       integer:: i,j,k,ii,jj,kk
       real(dp):: dx,dy,dz,nxx,nyy,nzz,diss,temp,vofeps
       real(dp):: nx(3),ny(3),nz(3),delt
       real(dp):: delxp,delxn,delyp,delyn,delzp,delzn,nxx1,nyy1,nzz1,nxy,nyz,nxz
       logical:: flag
       vofeps = 1.d-14
       delt = 2.d0
       nxx = 0.d0
       nyy = 0.d0
       nzz = 0.d0
       diss = 0.d0
       call Isinterface(i,j,k,flag)
       if(flag == .true.) then
          call Normal_Vector_Irre(i,j,k,nxx,nyy,nzz)
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
    
    subroutine Find_Distance(delx,dely,delz,nx,ny,nz,f,s)
       real(dp),intent(in):: delx,dely,delz,nx,ny,nz,f
       real(dp),intent(out):: s
       real(dp):: nx1,ny1,nz1,delx1,dely1,delz1,case_choose
       real(dp):: sc,sm,fc,v1,v2,v3,v31,v32
       integer:: case2
       nx1 = dabs(nx)
       ny1 = dabs(ny)
       nz1 = dabs(nz)
       delx1 = dmax1(nx1*delx,ny1*dely,nz1*delz)
       delz1 = dmin1(nx1*delx,ny1*dely,nz1*delz)
       dely1 = (nx1*delz+ny1*dely+nz1*delx)-delx1-delz1
       fc = dmin1(f,1.d0-f)
       sm = delx1+dely1+delz1
       sc = (6.d0*fc*delx1*dely1*delz1)**(1.d0/3.d0)
       if(sc<delz1) then
          call Final_Distance(f,sc,sm,s)
       else
          sc = 0.5d0*delz1+dsqrt(2.d0*fc*delx1*dely1-delz1**2.d0/12.d0)
          if(sc<dely1) then
             call Final_Distance(f,sc,sm,s)
             if(isnan(s)) then 
                pause
             end if
          else
             if(delx1>=(dely1+delz1)) then
                sc = fc*delx1+(dely1+delz1)/2.d0
                if(sc>=(dely1+delz1)) then
                   call Final_Distance(f,sc,sm,s)
                else
                   call Cubic_Equation_Solve(fc,delx1,dely1,delz1,sc,case2)
                   if(sc>=(dely1+delz1)) then
                      pause 'final_distance 546'
                   end if
                   call Final_Distance(f,sc,sm,s)
                end if
             else
                sc = fc*delx1+(dely1+delz1)/2.d0
                if(sc>=delx1) then
                   call Cubic_Equation_Solve(fc,delx1,dely1,delz1,sc)
                   if(sc<delx1) then
      !                pause 'final_distance 555'
                      call Cubic_Equation_Solve(fc,delx1,dely1,delz1,sc,case2)
                      if(sc>delx1) then
                         pause 'final_distance 558'
                       end if
                   end if
                   call Final_Distance(f,sc,sm,s)
                else
                   call Cubic_Equation_Solve(fc,delx1,dely1,delz1,sc,case2)
                   if(sc>delx1) then
                      pause 'final_distance 561'
                   end if
                   call Final_Distance(f,sc,sm,s)
                end if
             end if
          end if
       end if
    end subroutine Find_Distance
    subroutine Cubic_Equation_Solve(fc,delx1,dely1,delz1,sc,case2)
       real(dp),intent(in):: fc,delx1,dely1,delz1
       real(dp),intent(inout):: sc
       integer,optional:: case2
       real(dp):: a1,a2,a3,a0,p0,q0,delta,teta
       if(present(case2)) then
          a3 = 1.d0
          a2 = -3.d0*(dely1+delz1)
          a1 = 3.d0*(dely1**2.d0+delz1**2.d0)
          a0 = 6.d0*fc*delx1*dely1*delz1-(dely1**3.d0+delz1**3.d0)
          call Newton_Raphson(a3,a2,a1,a0,sc)
       else
          a3 = 2.d0
          a2 = -3.d0*(delx1+dely1+delz1)
          a1 = 3.d0*(delx1**2.d0+dely1**2.d0+delz1**2.d0)
          a0 = 6.d0*fc*delx1*dely1*delz1-(delx1**3.d0+dely1**3.d0+delz1**3.d0)
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
    subroutine Volume_Calc(delx,dely,delz,nx,ny,nz,alpha,vol)
       implicit none
       real(dp),intent(in):: delx,dely,delz,nx,ny,nz,alpha
       real(dp),intent(out):: vol
       real(dp):: nx1,ny1,nz1,alphamax,f1,f2,f3,fm1,fm2,fm3
       real(dp):: vol1,vol2,slop_eps = 1.d-12
       nx1 = dabs(nx)
       ny1 = dabs(ny)
       nz1 = dabs(nz)
       if(alpha<=0.d0) then
          vol = 0.d0
          return
       end if
       if(alpha-(nx1*delx+ny1*dely+nz1*delz)>=0.d0) then
          vol = 1.d0*delx*dely*delz
          return
       end if
       if(alpha-nx1*delx>0) then
          f1 = (alpha-nx1*delx)**3.d0
       else
          f1 = 0.d0
       end if
       if(alpha-ny1*dely>0) then
          f2 = (alpha-ny1*dely)**3.d0
       else
          f2 = 0.d0
       end if
       if(alpha-nz1*delz>0) then
          f3 = (alpha-nz1*delz)**3.d0
       else
          f3 = 0.d0
       end if
       alphamax = nx1*delx+ny1*dely+nz1*delz
       if(alpha-alphamax+nx1*delx>0) then
          fm1 = (alpha-alphamax+nx1*delx)**3.d0
       else
          fm1 = 0.d0
       end if
       if(alpha-alphamax+ny1*dely>0) then
          fm2 = (alpha-alphamax+ny1*dely)**3.d0
       else
          fm2 = 0.d0
       end if
       if(alpha-alphamax+nz1*delz>0) then
          fm3 = (alpha-alphamax+nz1*delz)**3.d0
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
    
    subroutine Volume_Fraction_Calc(delx,dely,delz,nx,ny,nz,s,vol)
       real(dp),intent(in):: nx,ny,nz,delx,dely,delz,s
       real(dp),intent(out):: vol
       real(dp):: nx1,ny1,nz1,delx1,dely1,delz1
       real(dp):: sc,sm,fc,vol1,vol2
       nx1 = dabs(nx)
       ny1 = dabs(ny)
       nz1 = dabs(nz)

       delx1 = dmax1(nx1*delx,ny1*dely,nz1*delz)
       delz1 = dmin1(nx1*delx,ny1*dely,nz1*delz)
       dely1 = (nx1*delx+ny1*dely+nz1*delz)-delx1-delz1
       sm = delx1+dely1+delz1
       sc = dmin1(s,sm-s)
       if(s<=0.d0) then
          vol =  0.d0 
          return
       end if
       if(s>=sm) then     ! be careful with this condition. it affects to region contains fluid
          vol =  1.d0*delx*dely*delz
          return
       end if
       if(sc<=delz1) then
          fc = sc**3.d0/(6.d0*delx1*dely1*delz1)
       elseif(sc<=dely1) then
          fc = (sc**2.d0-delz1*sc+delz1**2.d0/3.d0)/(2.d0*delx1*dely1)
       elseif(delx1>=dely1+delz1) then
          if(sc>=dely1+delz1) then
             fc = (2.d0*sc-dely1-delz1)/(2.d0*delx1)
          else
             fc = (sc**3.d0-(sc-delz1)**3.d0-(sc-dely1)**3.d0)
             fc = fc/(6.d0*delx1*dely1*delz1)
          end if
       else
          if(sc>delx1) then
             fc = (sc**3.d0-(sc-delz1)**3.d0-(sc-dely1)**3.d0-(sc-delx1)**3.d0)
          else
             fc = (sc**3.d0-(sc-delz1)**3.d0-(sc-dely1)**3.d0)
          end if
          fc = fc/(6.d0*delx1*dely1*delz1)
       end if
       if(s<=0.5d0*sm) then                 ! be careful with this step. the chosen of 
          vol = fc*delx*dely*delz          ! the region contains fluid
          return
       else
          vol = (1.d0-fc)*delx*dely*delz
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
! redistance from vof to level set
    subroutine Redistance_Levelset(Iinfo)
       implicit none
       type(Interface_Information),dimension(:,:,:),allocatable:: Iinfo
       integer:: i,j,k,ii,jj,kk,iii,jjj,kkk,ij,l,m,n,cas
       real(dp):: xv,yv,zv,dv,deuc,dij1
       real(dp):: xp,yp,zp,xoff,yoff,zoff,xfc,yfc,zfc,xs,ys,zs
       real(dp):: dis,dis1,nx,ny,nz
       logical:: pointp,face
       real(dp),dimension(:,:,:),allocatable:: phiaux
       integer,dimension(:,:,:),allocatable:: fix
       allocate(fix(imax,jmax,kmax))
       allocate(phiaux(imax,jmax,kmax))
       fix(:,:,:) = 0
       phiaux(:,:,:) = 1.d4
       do i = 1,imax
          do j = 1,jmax
             do k = 1,kmax
                if((0.d0+epsi)<vfl(i,j,k).and.vfl(i,j,k)<(1.d0-epsi)) then               
                   phiaux(i,j,k) = dsign(1.d0,0.5d0-vfl(i,j,k))*dabs(Iinfo(i,j,k)%dis)
                   fix(i,j,k) = 1
                   do ii = -band_width,band_width
                      do jj = -band_width,band_width
                         do kk = -band_width,band_width
                      ! determine the point xv on the boundary cell(i,j,k) with the 
                      ! shortest distance to the cell center of (i+-ii,j+-jj)
                           if(1<=i+ii.and.i+ii<=imax.and.1<=j+jj.and.j+jj<=jmax.and.&
                                                         1<=k+kk.and.k+kk<=kmax)then
                              if(vfl(i+ii,j+jj,k+kk)<0.d0+epsi.or.                  &
                                 vfl(i+ii,j+jj,k+kk)>1.d0-epsi) then
                                 if(fix(i+ii,j+jj,k+kk)==0) fix(i+ii,j+jj,k+kk) = 1
                                  
                                  nx = Iinfo(i,j,k)%nx
                                  ny = Iinfo(i,j,k)%ny
                                  nz = Iinfo(i,j,k)%nz
                                  dis1 = Iinfo(i,j,k)%dis
                                  l = max(-1,min(1,ii))
                                  m = max(-1,min(1,jj))
                                  n = max(-1,min(1,kk))
                                  xv = del_x*dble(l)/2.d0
                                  yv = del_y*dble(m)/2.d0
                                  zv = del_z*dble(n)/2.d0
                                  dv = nx*xv+ny*yv+nz*zv+dis1
                                  if(dv*dsign(1.d0,phi(i+ii,j+jj,k+kk))<=0.d0) then
                                      deuc = dsqrt((del_x*dble(ii)-xv)**2.d0+       &
                                                   (del_y*dble(jj)-yv)**2.d0+       &
                                                   (del_z*dble(kk)-zv)**2.d0)
                                      phiaux(i+ii,j+jj,k+kk) = dsign(1.d0,0.5d0-    &
                                                     vfl(i+ii,j+jj,k+kk))*          &
                                            dmin1(deuc,dabs(phiaux(i+ii,j+jj,k+kk)))
                     ! third step: find the projection of x' onto the interface  
                                  else
                                   dij1 = nx*del_x*dble(ii)+ny*del_y*dble(jj)+       &
                                          nz*del_z*dble(kk)+dis1
                                   xp = del_x*dble(ii)-dij1*nx
                                   yp = del_y*dble(jj)-dij1*ny
                                   zp = del_z*dble(kk)-dij1*nz
                                   pointp = Isinsidecell(xp,yp,zp)   
                                   if(pointp == .true.) then
                                     fix(i+ii,j+jj,k+kk) = 1     !  set to 2 to reduce time calculate thing
                                     deuc = dsqrt((del_x*dble(ii)-xp)**2.d0+        &
                                                  (del_y*dble(jj)-yp)**2.d0+        &
                                                  (del_z*dble(kk)-zp)**2.d0)
                                     phiaux(i+ii,j+jj,k+kk) = dsign(1.d0,0.5d0-     &
                                                  vfl(i+ii,j+jj,k+kk))*dmin1(deuc,  &
                                                  dabs(phiaux(i+ii,j+jj,k+kk)))
                                   else
                            ! forth step: find the corner of calculate 
                                    xoff = dmax1(dabs(xp)-0.5d0*del_x,0.d0)
                                    yoff = dmax1(dabs(yp)-0.5d0*del_y,0.d0) 
                                    zoff = dmax1(dabs(zp)-0.5d0*del_z,0.d0)
                                    xfc = dsign(1.d0,xp)*0.5d0*del_x
                                    yfc = dsign(1.d0,yp)*0.5d0*del_y 
                                    zfc = dsign(1.d0,zp)*0.5d0*del_z
                                    if(xoff*dabs(nx)>=yoff*dabs(ny).and.            &
                                       xoff*dabs(nx)>=zoff*dabs(nz))then
                                      face = Iscutface(xfc,nx,ny,nz,dis1,1)
                                      if(face==.true.) then                                                                            
                                         xs = xfc
                                         dis = nx*xs+dis1
                                         cas = 1
                                      else
                                       if(yoff*dabs(ny)>zoff*dabs(nz)) then
                                         ys = yfc
                                         dis = ny*ys+dis1
                                         cas = 2
                                       else
                                         zs = zfc
                                         dis = nz*zs+dis1
                                         cas = 3
                                       endif
                                      end if
                                      call Shortest_Point_2d(nx,ny,nz,dis,ii,jj,kk, &
                                              cas,xs,ys,zs)
                                    end if
                                    if(yoff*dabs(ny)>=xoff*dabs(nx).and.             &
                                       yoff*dabs(ny)>=zoff*dabs(nz))then
                                       face = Iscutface(yfc,nx,ny,nz,dis1,2)
                                      if(face==.true.) then                                                                            
                                         ys = yfc
                                         dis = ny*ys+dis1
                                         cas = 2
                                      else
                                       if(xoff*dabs(nx)>zoff*dabs(nz)) then
                                         xs = xfc
                                         dis = nx*xs+dis1
                                         cas = 1
                                       else
                                         zs = zfc
                                         dis = nz*zs+dis1
                                         cas = 3
                                       endif
                                      end if
                                      call Shortest_Point_2d(nx,ny,nz,dis,ii,jj,kk, &
                                              cas,xs,ys,zs)
                                    end if 
                                    if(zoff*dabs(nz)>=xoff*dabs(nx).and.            &
                                       zoff*dabs(nz)>=yoff*dabs(ny))then
                                      face = Iscutface(zfc,nx,ny,nz,dis1,3)
                                      if(face==.true.) then                                                                            
                                         zs = zfc
                                         dis = nz*zs+dis1
                                         cas = 3
                                      else
                                       if(xoff*dabs(nx)>yoff*dabs(ny)) then
                                         xs = xfc
                                         dis = nx*xs+dis1
                                         cas = 1
                                       else
                                         ys = yfc
                                         dis = ny*ys+dis1
                                         cas = 2
                                       endif
                                      end if
                                      
                                      call Shortest_Point_2d(nx,ny,nz,dis,ii,jj,kk, &
                                                            cas,xs,ys,zs)
                                    end if 
                                    deuc = dsqrt((del_x*dble(ii)-xs)**2.d0+         &
                                                 (del_y*dble(jj)-ys)**2.d0+         &
                                                 (del_z*dble(kk)-zs)**2.d0)
                                    phiaux(i+ii,j+jj,k+kk) = dsign(1.d0,0.5d0-      &
                                           vfl(i+ii,j+jj,k+kk))*                    &
                                            dmin1(deuc,dabs(phiaux(i+ii,j+jj,k+kk)))
                                     
                                   end if          
                                  end if
                               endif
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
                   phi(i,j,k) = (band_width+1)*del_x*dsqrt(3.d0)*                   &
                                        dsign(1.d0,0.5-vfl(i,j,k))
                end if
             end do
          end do
       end do
       deallocate(fix)
       deallocate(phiaux)
    end subroutine Redistance_Levelset
    
    function Heavisideph(phir) result(dhp)
       real(dp):: phir,dhp
       if(dabs(phir)<=epsilon) then
          dhp = -phir/epsilon**2.d0*(1.d0+dcos(pi*phir/epsilon))
       else
          dhp = 0.d0
       end if
    end function heavisideph
    
    function Isinsidecell(xp,yp,zp) result(logic)
       implicit none
       real(dp),intent(in):: xp,yp,zp
       logical:: logic
       if(-0.5d0*del_x<=xp.and.xp<=0.5d0*del_x.and.-0.5d0*del_y<=yp.and.yp<=0.5d0*  &
          del_y.and.-0.5d0*del_z<=zp.and.zp<=0.5d0*del_z) then
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
    
    function Iscutface(face_point,nx,ny,nz,dis,cas) result(logic)
       implicit none
       real(dp),intent(in):: face_point,nx,ny,nz,dis
       integer,intent(in):: cas
       logical:: logic
       real(dp):: point1,point2,point3,point4
       select case(cas)
          case(1)   ! in x surface
             point1 = face_point*nx-0.5d0*del_y*ny-0.5d0*del_z*nz+dis
             point2 = face_point*nx-0.5d0*del_y*ny+0.5d0*del_z*nz+dis
             point3 = face_point*nx+0.5d0*del_y*ny+0.5d0*del_z*nz+dis
             point4 = face_point*nx+0.5d0*del_y*ny-0.5d0*del_z*nz+dis
          case(2) ! in y surface
             point1 = -0.5d0*del_x*nx+face_point*ny-0.5d0*del_z*nz+dis
             point2 = -0.5d0*del_x*nx+face_point*ny+0.5d0*del_z*nz+dis
             point3 = 0.5d0*del_x*nx+face_point*ny+0.5d0*del_z*nz+dis
             point4 = 0.5d0*del_x*nx+face_point*ny-0.5d0*del_z*nz+dis
          case(3) ! in z surface
             point1 = -0.5d0*del_x*nx-0.5d0*del_y*ny+face_point*nz+dis
             point2 = -0.5d0*del_x*nx+0.5d0*del_y*ny+face_point*nz+dis
             point3 = 0.5d0*del_x*nx+0.5d0*del_y*ny+face_point*nz+dis
             point4 = 0.5d0*del_x*nx-0.5d0*del_y*ny+face_point*nz+dis
       end select
       if(point1*point2>=0.d0.and.point2*point3>=0.and.point3*point4>=0.d0          &
          .and.point4*point1>=0) then
          logic = .false.
          return
       else
          logic = .true.
          return 
       end if
    end function Iscutface 
    subroutine Shortest_Point_2d(nx1,ny1,nz1,dis1,ii,jj,kk,cas,xs,ys,zs)
       implicit none
       real(dp),intent(in):: nx1,ny1,nz1,dis1
       integer,intent(in):: ii,jj,kk,cas
       real(dp),intent(inout):: xs,ys,zs
       real(dp):: temp,nx,ny,nz,dis
       real(dp):: dij1,xp,yp,zp,xoff,yoff,zoff,xfc,yfc,zfc
       logical:: pointp
       select case(cas)
          case(1) ! for x face
             temp = dsqrt(ny1**2.d0+nz1**2.d0+epsi**2.d0)
             ny = ny1/temp
             nz = nz1/temp
             dis = dis1/temp
             dij1 = ny*del_y*dble(jj)+nz*del_z*dble(kk)+dis
             yp = del_y*dble(jj)-dij1*ny
             zp = del_z*dble(kk)-dij1*nz                         
             pointp = Isinsiderect(yp,zp,-0.5d0*del_y,0.5d0*del_y,                  &
                                   -0.5d0*del_z,0.5d0*del_z)   
             if(pointp == .true.) then
                ys = yp
                zs = zp
                return
             else
             ! forth step: find the corner of 
                yoff = dmax1(dabs(yp)-0.5d0*del_y,0.d0)
                zoff = dmax1(dabs(zp)-0.5d0*del_z,0.d0) 
                yfc = dsign(1.d0,yp)*0.5d0*del_y
                zfc = dsign(1.d0,zp)*0.5d0*del_z
                if(yoff*dabs(ny)>zoff*dabs(nz)) then
                   ys = yfc
                   zs = (dis+ny*ys)/(-nz)
                   return 
                else
                   zs = zfc
                   ys = (dis+nz*zs)/(-ny)
                   return
                end if
             end if
          case(2) ! for y face     
             temp = dsqrt(nx1**2.d0+nz1**2.d0+epsi**2.d0)
             nx = nx1/temp
             nz = nz1/temp
             dis = dis1/temp
             dij1 = nx*del_x*dble(ii)+nz*del_z*dble(kk)+dis
             xp = del_x*dble(ii)-dij1*nx
             zp = del_z*dble(kk)-dij1*nz
             pointp = Isinsiderect(xp,zp,-0.5d0*del_x,0.5d0*del_x,                  &
                                   -0.5d0*del_z,0.5d0*del_z)   
             if(pointp == .true.) then
                xs = xp
                zs = zp
                return
             else
             ! forth step: find the corner of 
                xoff = dmax1(dabs(xp)-0.5d0*del_x,0.d0)
                zoff = dmax1(dabs(zp)-0.5d0*del_z,0.d0) 
                xfc = dsign(1.d0,xp)*0.5d0*del_x
                zfc = dsign(1.d0,zp)*0.5d0*del_z
                if(xoff*dabs(nx)>zoff*dabs(nz)) then
                   xs = xfc
                   zs = (dis+nx*xs)/(-nz)
                   return 
                else
                   zs = zfc
                   xs = (dis+nz*zs)/(-nx)
                   return
                end if
             end if
          case(3) ! for z face  
             temp = dsqrt(nx1**2.d0+ny1**2.d0+epsi**2.d0)
             nx = nx1/temp
             ny = ny1/temp
             dis = dis1/temp
             dij1 = nx*del_x*dble(ii)+ny*del_y*dble(jj)+dis
             xp = del_x*dble(ii)-dij1*nx
             yp = del_y*dble(jj)-dij1*ny
             pointp = Isinsiderect(xp,yp,-0.5d0*del_x,0.5d0*del_x,                  &
                                   -0.5d0*del_y,0.5d0*del_y)   
             if(pointp == .true.) then
                xs = xp
                ys = yp
                return
             else
             ! forth step: find the corner of 
                xoff = dmax1(dabs(xp)-0.5d0*del_x,0.d0)
                yoff = dmax1(dabs(yp)-0.5d0*del_y,0.d0) 
                xfc = dsign(1.d0,xp)*0.5d0*del_x
                yfc = dsign(1.d0,yp)*0.5d0*del_y
                if(xoff*dabs(nx)>yoff*dabs(ny)) then
                   xs = xfc
                   ys = (dis+nx*xs)/(-ny)
                   return 
                else
                   ys = yfc
                   xs = (dis+ny*ys)/(-nx)
                   return
                end if
             end if
          end select 
    end subroutine Shortest_Point_2d
    
    subroutine Boundary_Condition(vari)
       implicit none
       integer i,j,k
       real(dp),dimension(:,:,:),intent(inout),allocatable:: vari
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
    
    subroutine Normal_Vector_Irre(i,j,k,nx,ny,nz)
       implicit none
       integer,intent(in):: i,j,k
       real(dp),intent(out):: nx,ny,nz
       integer:: case_deri,dx,dy,dz
       real(dp):: qi(-1:1),qj(-1:1),qk(-1:1)
       real(dp):: vx,vy,vz
       ! define qi for Dx
       if(i>=3) then
       ! for qi(-1)
          vx = 0.5d0*(phi(i,j,k)-phi(i-2,j,k))/del_x
       else
          vx = (phi(i,j,k)-phi(i-1,j,k))/del_x
       end if
       vy = 0.5d0*(phi(i-1,j+1,k)-phi(i-1,j-1,k))/del_y
       vz = 0.5d0*(phi(i-1,j,k+1)-phi(i-1,j,k-1))/del_z
       qi(-1) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! for qi(1) 
       if(i<=imax-2) then
          vx = 0.5d0*(phi(i+2,j,k)-phi(i,j,k))/del_x
       else
          vx = (phi(i+1,j,k)-phi(i,j,k))/del_x
       end if
       vy = 0.5d0*(phi(i+1,j+1,k)-phi(i+1,j-1,k))/del_y
       vz = 0.5d0*(phi(i+1,j,k+1)-phi(i+1,j,k-1))/del_z
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
          vy = 0.5d0*(phi(i,j,k)-phi(i,j-2,k))/del_y
       else
          vy = (phi(i,j,k)-phi(i,j-1,k))/del_y
       end if
       vx = 0.5d0*(phi(i+1,j-1,k)-phi(i-1,j-1,k))/del_x
       vz = 0.5d0*(phi(i,j-1,k+1)-phi(i,j-1,k-1))/del_z
       qj(-1) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! define qj(1)
       vx = 0.5d0*(phi(i+1,j+1,k)-phi(i-1,j+1,k))/del_x
       if(j<=jmax-2) then
          vy = 0.5d0*(phi(i,j+2,k)-phi(i,j,k))/del_y
       else
          vy = (phi(i,j+1,k)-phi(i,j,k))/del_y
       end if
       vz = 0.5d0*(phi(i,j+1,k+1)-phi(i,j+1,k-1))/del_z
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
       vx = 0.5d0*(phi(i+1,j,k-1)-phi(i-1,j,k-1))/del_x
       vy = 0.5d0*(phi(i,j+1,k-1)-phi(i,j-1,k-1))/del_y
       if(k>=3) then
          vz = 0.5d0*(phi(i,j,k)-phi(i,j,k-2))/del_z
       else
          vz = (phi(i,j,k)-phi(i,j,k-1))/del_z
       end if
       qk(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       ! qk(1)
       vx = 0.5d0*(phi(i+1,j,k+1)-phi(i-1,j,k+1))/del_x
       vy = 0.5d0*(phi(i,j+1,k+1)-phi(i,j-1,k+1))/del_y
       if(k<=kmax-2) then
          vz = 0.5d0*(phi(i,j,k+2)-phi(i,j,k))/del_z
       else
          vz = (phi(i,j,k+1)-phi(i,j,k))/del_z
       end if
       qk(1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0+vz**2.d0))
       if(qk(-1)<eta.and.qk(1)>=eta) then
          dz = -1
       elseif(qk(-1)>=eta.and.qk(1)<eta) then
          dz = 1
       else
          dz = 0
       end if
       if(dx==-1) nx = (phi(i,j,k)-phi(i-1,j,k))/del_x
       if(dx==1) nx = (phi(i+1,j,k)-phi(i,j,k))/del_x
       if(dx==0) nx = 0.5d0*(phi(i+1,j,k)-phi(i-1,j,k))/del_x
       if(dy==-1) ny = (phi(i,j,k)-phi(i,j-1,k))/del_y
       if(dy==1) ny = (phi(i,j+1,k)-phi(i,j,k))/del_y
       if(dy==0) ny = 0.5d0*(phi(i,j+1,k)-phi(i,j-1,k))/del_y
       if(dz==-1) nz = (phi(i,j,k)-phi(i,j,k-1))/del_z
       if(dz==1) nz = (phi(i,j,k+1)-phi(i,j,k))/del_z
       if(dz==0) nz = 0.5d0*(phi(i,j,k+1)-phi(i,j,k-1))/del_z
    end subroutine Normal_Vector_Irre
end module clsvof     
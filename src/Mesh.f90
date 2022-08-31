Module Mesh
    !    
    USE PrecisionVar
    USE MPI
    !
    Private
    !
    Integer(kind=it4b),public :: Imax,Jmax,Kmax
    !
    Type,public :: Point
      Real(kind=dp) ::x,y,z
    End type Point
    !
    Type,public :: Grid
      Integer*8 :: Grid_Id
      Real(kind=dp),dimension(:,:,:),allocatable :: x,y,z,dx,dy,dz
      Real(kind=dp) :: Lref
    End type Grid
    !
    Public :: InitialGrid,HYPRE_CreateGrid,InitialUVGrid
    !
    Interface InitialGrid
      Module procedure InitialGrid
    End interface
    Interface HYPRE_CreateGrid
      Module procedure HYPRE_CreateGrid
    End interface
    Interface InitialUVGrid
      Module procedure InitialUVGrid
    End interface
    !
    Contains
      !      
      Subroutine InitialGrid(NI,NJ,NK,Irec,Jrec,Krec,Lref, &
                             SPoint,EPoint,ReS,ReE, TGrid)
        !                                                   
        Implicit none
        !
        Integer(kind=it4b),intent(in) :: Irec,Jrec,Krec,NI,NJ,NK
        Real(kind=dp),intent(in)      :: Lref
        Type(Point),intent(in)        :: SPoint,EPoint,ReS,ReE
        Type(Grid),intent(inout)      :: TGrid

        Real(kind=dp),dimension(:),allocatable:: x,y,z
        Integer:: i,j,k,IJKsize
        Real(kind=dp):: beta,dx,dy,dz,dl
        !
        Allocate(x(NI))
        allocate(y(NJ))
        Allocate(z(NK))
        !
        IJKsize = (NJ-Jrec)/2
        !
        dy = (ReE%y-ReS%y)/dble(Jrec-1)
        dl = EPoint%y-ReE%y
        !
        call NewtonRaphson(IJKsize,dl,dy, beta)
        !
        y(1) = SPoint%y
        !
        do j = 2,2+IJKsize
          y(j) = y(j-1)+dy*beta**(IJKsize-(j-2))
        end do
        !
        do j = 2+IJKsize,IJKsize+Jrec
          y(j) = y(j-1)+dy
        end do
        !
        do j = 1+IJKsize+Jrec,NJ
          y(j) = y(j-1)+dy*beta**(j-IJKsize-Jrec)
        end do
        !
        IJKsize = (NK-Krec)/2
        !
        dz = (ReE%z-ReS%z)/dble(Krec-1)
        dl = ReS%z-SPoint%z
        !
        call NewtonRaphson(IJKsize,dl,dz, beta)
        !
        z(1) = SPoint%z
        !
        do k = 2,2+IJKsize
          z(k) = z(k-1)+dz*beta**(IJKsize-(k-2))
        end do
        !
        do k = 2+IJKsize,IJKsize+Krec
          z(k) = z(k-1)+dz
        end do
        !
        dl = EPoint%z-ReE%z
        !
        call NewtonRaphson(IJKsize,dl,dz, beta)
        !
        do k = 1+IJKsize+Krec,NK
          z(k) = z(k-1)+dz*beta**(k-IJKsize-Krec)
        end do
        !
        IJKsize = (NI-Irec)/2
        !
        dx = (ReE%x-ReS%x)/dble(Irec-1)
        dl = ReS%x-SPoint%x
        !
        call NewtonRaphson(IJKsize,dl,dx, beta)
        !
        x(1) = SPoint%x
        !
        do i = 2,1+IJKsize
          x(i) = x(i-1)+dx*beta**(IJKsize-(i-2))
        end do
        !
        do i = 2+IJKsize,IJKsize+Irec
          x(i) = x(i-1)+dx
        end do
        !
        dl = EPoint%x-ReE%x
        !
        call NewtonRaphson(IJKsize,dl,dx, beta)
        !
        do i = 1+IJKsize+Irec,NI
          x(i) = x(i-1)+dx*beta**(i-IJKsize-Irec)
        end do
        !
        Do i = 1,NI-1
          Do j = 1,NJ-1
            Do k = 1,NK-1
              TGrid%x(i,j,k) = 0.5d0*(x(i+1)+x(i))
              TGrid%dx(i,j,k) = x(i+1)-x(i)
              TGrid%y(i,j,k) = 0.5d0*(y(j+1)+y(j))
              TGrid%dy(i,j,k) = y(j+1)-y(j)
              TGrid%z(i,j,k) = 0.5d0*(z(k+1)+z(k))
              TGrid%dz(i,j,k) = z(k+1)-z(k)
            End do
          End do
        End do
        !
        TGrid%Lref = Lref
        !
        !  TGrid%dx(:,:,:) = (EPoint%x/Lref-SPoint%x/Lref)/dble(Imax-1)
        !  TGrid%dy(:,:,:) = (EPoint%y/Lref-SPoint%y/Lref)/dble(Jmax-1)
        !  TGrid%dz(:,:,:) = (EPoint%z/Lref-SPoint%z/Lref)/dble(Kmax-1)
        ! Do i = 1,Imax
        !   Do j = 1,Jmax
        !     Do k = 1,Kmax
        !       TGrid%x(i,j,k) = SPoint%x/Lref+TGrid%dx(i,j,k)*dble(i-1)
        !       TGrid%y(i,j,k) = SPoint%y/Lref+TGrid%dy(i,j,k)*dble(j-1)
        !       TGrid%z(i,j,k) = SPoint%z/Lref+TGrid%dz(i,j,k)*dble(k-1)
        !       If(TGrid%x(i,j,k)>1.d5.or.TGrid%dy(i,j,k)>1.d5.or.TGrid%dz(i,j,k)>1.d5) then
        !         print*, TGrid%dx(i,j,k),TGrid%dy(i,j,k),TGrid%dz(i,j,k)
        !         print*, i,j,k
        !         print*,
        !       End if
        !     End do
        !   End do
        ! End do
        !
        Deallocate(x,y,z)
        !
      End subroutine InitialGrid
      !
      Subroutine InitialUVGrid(UVW,Lref,PGrid, TGrid)
        !      
        Implicit none
        !
        Integer(kind=it4b),intent(in):: UVW
        Real(kind=dp),intent(in):: Lref
        Type(Grid),intent(in):: PGrid
        Type(Grid),intent(inout):: TGrid

        Integer:: i,j,k
        !
        TGrid%Lref = Lref
        !
        If(UVW==0) then ! for UGrid
          !      
          Do i = 1,Imax
            Do j = 1,Jmax
              Do k = 1,Kmax
                !
                TGrid%x(i,j,k) = PGrid%x(i,j,k)+0.5d0*PGrid%dx(i,j,k)
                TGrid%y(i,j,k) = PGrid%y(i,j,k)
                TGrid%z(i,j,k) = PGrid%z(i,j,k)
                TGrid%dy(i,j,k) = PGrid%dy(i,j,k)
                TGrid%dz(i,j,k) = PGrid%dz(i,j,k)
                !
                If(i<Imax) then
                  TGrid%dx(i,j,k) = PGrid%x(i+1,j,k)-PGrid%x(i,j,k)
                Else
                  TGrid%dx(i,j,k) = PGrid%dx(i,j,k)
                End if
                !
                If(TGrid%x(i,j,k)>=45.d0) then
                  print*, TGrid%x(i,j,k)
                  print*, PGrid%dx(i,j,k)
                  print*, PGrid%x(i,j,k)
                End if
                !
              End do
            End do
          End do
          !
        Elseif(UVW==1) then ! for V-Cell
          !      
          Do i = 1,Imax
            Do j = 1,Jmax
              Do k = 1,Kmax
                !
                TGrid%x(i,j,k) = Pgrid%x(i,j,k)
                TGrid%z(i,j,k) = PGrid%z(i,j,k)
                TGrid%y(i,j,k) = PGrid%y(i,j,k)+0.5d0*PGrid%dy(i,j,k)
                TGrid%dx(i,j,k) = PGrid%dx(i,j,k)
                TGrid%dz(i,j,k) = PGrid%dz(i,j,k)
                !
                If(j<Jmax) then
                  TGrid%dy(i,j,k) = PGrid%y(i,j+1,k)-PGrid%y(i,j,k)
                Else
                  TGrid%dy(i,j,k) = PGrid%dy(i,j,k)
                End if
                !
              End do
            End do
          End do
          !  
        Else ! for WCell
          !      
          Do i = 1,Imax
            Do j = 1,Jmax
              Do k = 1,Kmax
                !
                TGrid%x(i,j,k) = Pgrid%x(i,j,k)
                TGrid%y(i,j,k) = PGrid%y(i,j,k)
                TGrid%z(i,j,k) = PGrid%z(i,j,k)+0.5d0*PGrid%dz(i,j,k)
                TGrid%dx(i,j,k) = PGrid%dx(i,j,k)
                TGrid%dy(i,j,k) = PGrid%dy(i,j,k)
                !
                If(k<Kmax) then
                  TGrid%dz(i,j,k) = PGrid%z(i,j,k+1)-PGrid%z(i,j,k)
                Else
                  TGrid%dz(i,j,k) = PGrid%dz(i,j,k)
                End if
                !
              End do
            End do
          End do
          !
        End if
        !
      End subroutine
      !
      Subroutine HYPRE_CreateGrid(TGrid)
        !      
        Implicit none
        !
        Type(Grid),intent(inout):: TGrid

        Integer:: ilower(0:2),iupper(0:2)
        !
        ilower(0) = 1
        ilower(1) = 1
        ilower(2) = 1
        iupper(0) = Imax
        iupper(1) = Jmax
        iupper(2) = Kmax
        !
        Call HYPRE_StructGridCreate(MPI_COMM_WORLD,3,TGrid%Grid_Id,ierr)
        Call HYPRE_StructGridSetExtents(TGrid%Grid_Id,ilower,iupper,ierr)
        Call HYPRE_StructGridAssemble(TGrid%Grid_Id,ierr)
        !
      End subroutine HYPRE_CreateGrid

      !
      Subroutine NewtonRaphson(IJsize,dl,dx, beta)
        !
        Real(kind=dp),intent(in):: dl,dx
        Integer(kind=it4b),intent(in):: IJsize
        Real(kind=dp),intent(inout):: beta

        Real(kind=dp):: tol,fx,dfx
        !
        tol = 1.d0
        beta = 1.001d0
        !
        do while(tol>1.d-14)
          fx = dx*(beta**(IJsize+1)-beta)/(beta-1.d0)-dl
          dfx = dx*(IJsize*beta**(IJsize+1)-(IJsize+1)*beta**IJsize+1)/(beta-1)**2.d0
          tol = dabs(fx/dfx)
          beta = beta-fx/dfx
        end do
        !
      End subroutine
      !
End Module Mesh

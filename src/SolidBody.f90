!! This module implements routines for keeping track of a solid body moving in
!  the domain and for interpolating its position on the main grid
!
MODULE SolidBody

  USE PrecisionVar
  USE Geometry,        ONLY : TGpoint => Gpoint, TVector => vector, CenterPoint, CheckPointTriangle, Distance, CheckPointLine, &
                              RayCutTriangle
  USE Mesh,            ONLY : TGrid => Grid
  USE CutCell,         ONLY : TCell => Cell, GridPreProcess, NumberExternalCell, DefineMomentumExchangeCell, &
                              NewCellFace
  USE StateVariables,  ONLY : TSolverTime => SolverTime
  USE Clsvof,          ONLY : Volume_Fraction_Calc, Normal_Vector_Irre
  USE STL,             ONLY : TsimcoSTL => simcoSTL

  IMPLICIT NONE

  PRIVATE

  REAL(KIND=dp), PARAMETER :: tolpar = 1.d-14, pi=4.0_dp*DATAN(1.0_dp)

  !! Solid body type. It has its own uniform mesh, with cells containing level set, normal and vof.
  !  TODO: This could use generic grid and cell objects instead of redefining everything
  !
  TYPE, PUBLIC :: TSolidBody
    REAL(dp)              :: CMpoint(3)       ! Coordinates of the center of mass
    REAL(dp)              :: CMpoint0(3)      ! Value at start of time step
    REAL(dp)              :: CMorient(3)      ! Orientation angles of solid at the center of mass
                                              ! (Roll, Pitch, Yaw), (p, q, r), (Phi, Theta, Psi)
    REAL(dp)              :: CMorient0(3)     ! Value at start of time step
    REAL(dp)              :: rotMatLtoG(3,3)  ! Rotation matrix to express a vector defined in main grid coordinates
                                              ! in terms of local grid coordinates
    REAL(dp)              :: rotMatLtoG0(3,3) ! Matrix at start of time step
    INTEGER               :: Imax, Jmax, Kmax ! Size of local mesh
    REAL(dp)              :: dx, dy, dz       ! Sizes of the cells (Uniform mesh)
    REAL(dp)              :: dt               ! Time step
    REAL(dp), ALLOCATABLE :: x(:), y(:), z(:) ! Coordinates of the cell centers
    REAL(dp), ALLOCATABLE :: phi(:,:,:)       ! Level set function on local mesh
    REAL(dp), ALLOCATABLE :: nx (:,:,:)       ! x-component of normal to surface on local mesh
    REAL(dp), ALLOCATABLE :: ny (:,:,:)       ! y-component of normal to surface on local mesh
    REAL(dp), ALLOCATABLE :: nz (:,:,:)       ! z-component of normal to surface on local mesh
    REAL(dp), ALLOCATABLE :: vof(:,:,:)       ! Volume of fluid on local mesh
    TYPE(TsimcoSTL)       :: ComSTL           ! Communicator for STL file
  CONTAINS
    PROCEDURE, PASS(this), PUBLIC  :: setInitialPosition
    PROCEDURE, PASS(this), PUBLIC  :: setUpSolid
    PROCEDURE, PASS(this), PUBLIC  :: solidVelocity
    PROCEDURE, PASS(this), PUBLIC  :: advance
    PROCEDURE, PASS(this), PRIVATE :: interpolateLvS
    PROCEDURE, PASS(this), PRIVATE :: LvsObject
    PROCEDURE, PASS(this), PRIVATE :: calcRotationMatrix
    PROCEDURE, PASS(this), PRIVATE :: transferMeshGtoL
    PROCEDURE, PASS(this), PRIVATE :: transferMeshLtoG
  END TYPE TSolidBody

  INTERFACE TSolidBody
     MODULE PROCEDURE construct1
  END INTERFACE TSolidBody

CONTAINS

  !! Contructor which allocates arrays with given mesh size
  !  They are not initialised
  !
  TYPE(TSolidBody) FUNCTION construct1(Imax, Jmax, Kmax) RESULT( this )
    !
    INTEGER, INTENT(in) :: Imax, Jmax, Kmax        ! Size of the local mesh
    !
    IF( MOD(Imax,2) == 1 .OR. MOD(Jmax,2) == 1 .OR. MOD(Kmax,2) == 1 ) THEN
       WRITE(*,*) "SolidBody was written for even grid sizes"
       WRITE(*,*) "This restriction can be removed when center of mass is properly calculated"
       STOP
    END IF

    this%Imax = Imax
    this%Jmax = Jmax
    this%Kmax = Kmax

    ALLOCATE( this%x(Imax), this%y(Jmax), this%z(Kmax) )
    ALLOCATE( this%phi(Imax, Jmax, Kmax), this%vof(Imax, Jmax, Kmax) )
    ALLOCATE( this%nx(Imax, Jmax, Kmax), this%ny(Imax, Jmax, Kmax), this%nz(Imax, Jmax, Kmax) )

  END FUNCTION construct1


  SUBROUTINE purge(this)
    !
    CLASS(TSolidBody), INTENT(inout) :: this

    IF( ALLOCATED( this%x   ) ) DEALLOCATE( this%x   )
    IF( ALLOCATED( this%y   ) ) DEALLOCATE( this%y   )
    IF( ALLOCATED( this%z   ) ) DEALLOCATE( this%z   )
    IF( ALLOCATED( this%phi ) ) DEALLOCATE( this%phi )
    IF( ALLOCATED( this%vof ) ) DEALLOCATE( this%vof )
    IF( ALLOCATED( this%nx  ) ) DEALLOCATE( this%nx  )
    IF( ALLOCATED( this%ny  ) ) DEALLOCATE( this%ny  )
    IF( ALLOCATED( this%nz  ) ) DEALLOCATE( this%nz  )

  END SUBROUTINE purge


  SUBROUTINE setUpSolid( this, ComSTL )
    !
    CLASS(TSolidBody), INTENT(inout) :: this
    TYPE(TsimcoSTL),   INTENT(in)    :: ComSTL
    !
    INTEGER  :: i, id
    REAL(dp) :: minDim, maxDim
    
    ! TODO: Find center of mass of solid and set it as origin of local grid 
    ! For now, set origin as mid-point of shape

    ! Find dx/dy/dz to contain the whole solid given the mesh size
    DO id = 1, 3
       minDim =  HUGE(0.0_dp)
       maxDim = -HUGE(0.0_dp)
       DO i = 1, 3
          minDim = MIN(minDim,MINVAL(ComSTL%tri(:)%ptr(i)%p(id)))
          maxDim = MAX(maxDim,MAXVAL(ComSTL%tri(:)%ptr(i)%p(id)))
       END DO
       SELECT CASE( id )
       CASE( 1 ) 
          this%dx = (maxDim - minDim)/(this%Imax-8)  ! Keep 4 cells margin around shape
          this%x(this%Imax/2) = -0.5_dp * this%dx    ! Assume even grid size
       CASE( 2 ) 
          this%dy = (maxDim - minDim)/(this%Jmax-8)  ! Keep 4 cells margin around shape
          this%y(this%Jmax/2) = -0.5_dp * this%dy    ! Assume even grid size
       CASE( 3 ) 
          this%dz = (maxDim - minDim)/(this%Kmax-8)  ! Keep 4 cells margin around shape
          this%z(this%Kmax/2) = -0.5_dp * this%dz    ! Assume even grid size
       END SELECT
    END DO

    DO i = this%Imax/2-1, 1, -1
       this%x(i) = this%x(i+1) - this%dx
    END DO
    DO i = this%Imax/2+1, this%Imax
       this%x(i) = this%x(i-1) + this%dx
    END DO

    DO i = this%Jmax/2-1, 1, -1
       this%y(i) = this%y(i+1) - this%dy
    END DO
    DO i = this%Jmax/2+1, this%Jmax
       this%y(i) = this%y(i-1) + this%dy
    END DO

    DO i = this%Kmax/2-1, 1, -1
       this%z(i) = this%z(i+1) - this%dz
    END DO
    DO i = this%Kmax/2+1, this%Kmax
       this%z(i) = this%z(i-1) + this%dz
    END DO

    ! Calculate level set, vof and normal on the above set-up mesh
    CALL this%LvsObject(ComSTL)

  END SUBROUTINE setUpSolid


  !! The subroutine is used to compute the distance from the grid points to
  !! the object surface. Only several layers from the body surface contains
  !! the actual value of level set function. For remaining points, the level
  !! set is set up to a constant value.
  !
  SUBROUTINE LvsObject(this, ComSTL)
    !
    CLASS(TSolidBody), INTENT(inout)  :: this
    TYPE(TsimcoSTL),   INTENT(in)     :: ComSTL
    !
    INTEGER(KIND=it4b)                :: i,j,k,ii,jj,kk,i1,j1,k1,bw
    INTEGER(KIND=it4b), ALLOCATABLE   :: GridIndex(:,:)
    INTEGER(KIND=it4b), ALLOCATABLE   :: BoxSizePos(:,:), BoxSizeNeg(:,:)
    REAL(KIND=dp)                     :: tempd,dSur,d1,d2,d3
    REAL(KIND=dp)                     :: nx, ny, nz, s, vol, maxLvs
    INTEGER(KIND=it4b)                :: iboxL, iboxU, jboxU, jboxL,&
                                         kboxU, kboxL
    REAL(KIND=dp)                     :: MaxObjX, MinObjX, MaxObjY, &
                                         MinObjY, MaxObjZ, MinObjZ
    TYPE(TGpoint)                     :: pt,ProjPt
    TYPE(Tvector)                     :: vec,v12,v1t,npoint
    LOGICAL                           :: Isin, IsinE1, IsinE2, IsinE3
    TYPE(TGpoint)                     :: TriCenter
    REAL(dp) :: start,finish

    bw = 4
    maxLvs=1.d4

    ! Compute the level set function for a box with certain size
    ! The size of the box is determined based on the minimum and maximum values of 
    ! the ship's location

    MaxObjX = -1.d6; MaxObjY = -1.d6; MaxObjZ = -1.d6
    MinObjX = 1.d6; MinObjY = 1.d6; MinObjZ = 1.d6

    this%phi = maxLvs

    !do i = 1, ComSTL%ntri
    !  ! Compute maximum object position in the x direction
    !  if(ComSTL%tri(i)%pTr(1)%p(1)>MaxObjX) MaxObjX = ComSTL%tri(i)%PTr(1)%p(1)
    !  if(ComSTL%tri(i)%pTr(2)%p(1)>MaxObjX) MaxObjX = ComSTL%tri(i)%PTr(2)%p(1)
    !  if(ComSTL%tri(i)%pTr(3)%p(1)>MaxObjX) MaxObjX = ComSTL%tri(i)%PTr(3)%p(1)
    !  ! Compute maximum object position in the y direction
    !  if(ComSTL%tri(i)%pTr(1)%p(2)>MaxObjY) MaxObjY = ComSTL%tri(i)%PTr(1)%p(2)
    !  if(ComSTL%tri(i)%pTr(2)%p(2)>MaxObjY) MaxObjY = ComSTL%tri(i)%PTr(2)%p(2)
    !  if(ComSTL%tri(i)%pTr(3)%p(2)>MaxObjY) MaxObjY = ComSTL%tri(i)%PTr(3)%p(2)
    !  ! Compute maximum object position in the z direction
    !  if(ComSTL%tri(i)%pTr(1)%p(3)>MaxObjZ) MaxObjZ = ComSTL%tri(i)%PTr(1)%p(3)
    !  if(ComSTL%tri(i)%pTr(2)%p(3)>MaxObjZ) MaxObjZ = ComSTL%tri(i)%PTr(2)%p(3)
    !  if(ComSTL%tri(i)%pTr(3)%p(3)>MaxObjZ) MaxObjZ = ComSTL%tri(i)%PTr(3)%p(3)
    !  ! Compute minimum object position in the x direction
    !  if(ComSTL%tri(i)%pTr(1)%p(1)<MinObjX) MinObjX = ComSTL%tri(i)%PTr(1)%p(1)
    !  if(ComSTL%tri(i)%pTr(2)%p(1)<MinObjX) MinObjX = ComSTL%tri(i)%PTr(2)%p(1)
    !  if(ComSTL%tri(i)%pTr(3)%p(1)<MinObjX) MinObjX = ComSTL%tri(i)%PTr(3)%p(1)
    !  ! Compute minimum object position in the y direction
    !  if(ComSTL%tri(i)%pTr(1)%p(2)<MinObjY) MinObjY = ComSTL%tri(i)%PTr(1)%p(2)
    !  if(ComSTL%tri(i)%pTr(2)%p(2)<MinObjY) MinObjY = ComSTL%tri(i)%PTr(2)%p(2)
    !  if(ComSTL%tri(i)%pTr(3)%p(2)<MinObjY) MinObjY = ComSTL%tri(i)%PTr(3)%p(2)
    !  ! Compute minimum object position in the z direction
    !  if(ComSTL%tri(i)%pTr(1)%p(3)<MinObjZ) MinObjZ = ComSTL%tri(i)%PTr(1)%p(3)
    !  if(ComSTL%tri(i)%pTr(2)%p(3)<MinObjZ) MinObjZ = ComSTL%tri(i)%PTr(2)%p(3)
    !  if(ComSTL%tri(i)%pTr(3)%p(3)<MinObjZ) MinObjZ = ComSTL%tri(i)%PTr(3)%p(3)
    !end do
    !! Determine the box which contains the object
    !do i = 2, Imax-1
    !  if(this%x(i,1,1)>=MinObjX.and.this%x(i-1,1,1)<MinObjX) iboxL = i-1
    !  if(this%x(i,1,1)>=MaxObjX.and.this%x(i-1,1,1)<MaxObjX) iboxU = i
    !end do
    !do j = 2, Jmax-1
    !  if(this%y(1,j,1)>=MinObjY.and.this%y(1,j-1,1)<MinObjY) jboxL = j-1
    !  if(this%y(1,j,1)>=MaxObjY.and.this%y(1,j-1,1)<MaxObjY) jboxU = j
    !end do
    !do k = 2, Kmax-1
    !  if(this%z(1,1,k)>=MinObjZ.and.this%z(1,1,k-1)<MinObjZ) kboxL = k-1
    !  if(this%z(1,1,k)>=MaxObjZ.and.this%z(1,1,k-1)<MaxObjZ) kboxU = k
    !end do
    !iboxL = max(1,iboxL-bw)
    !iboxU = min(Imax,iboxU+bw)
    !jboxL = max(1,jboxL-bw)
    !jboxU = min(Jmax,jboxU+bw)
    !kboxL = max(1,kboxL-bw)
    !kboxU = min(Kmax,kboxU+bw)
    ! print*, 'Clsvof.f90 155'
    ! print*, this%dx(1,1,1), this%dy(1,1,1), this%dz(1,1,1)
    ! print*, iboxL, MinObjX, iboxU, MaxObjX
    ! print*, jboxL, MinObjY, jboxU, MaxObjY
    ! print*, kboxL, MinObjZ, kboxU, MaxObjZ

    allocate(GridIndex(ComSTL%ntri,3))
    allocate(BoxSizePos(ComSTL%ntri,3))
    allocate(BoxSizeNeg(ComSTL%ntri,3))


    call cpu_time(start)
    ! Determine local box size for each triangles
    do i=1,ComSTL%ntri
      call CenterPoint(ComSTL%tri(i),TriCenter)
      MaxObjX = dmax1((ComSTL%tri(i)%pTr(1)%p(1)-TriCenter%p(1)),        &
                      (ComSTL%tri(i)%pTr(2)%p(1)-TriCenter%p(1)),        &
                      (ComSTL%tri(i)%pTr(3)%p(1)-TriCenter%p(1)))
      MaxObjY = dmax1((ComSTL%tri(i)%pTr(1)%p(2)-TriCenter%p(2)),        &
                      (ComSTL%tri(i)%pTr(2)%p(2)-TriCenter%p(2)),        &
                      (ComSTL%tri(i)%pTr(3)%p(2)-TriCenter%p(2)))
      MaxObjZ = dmax1((ComSTL%tri(i)%pTr(1)%p(3)-TriCenter%p(3)),        &
                      (ComSTL%tri(i)%pTr(2)%p(3)-TriCenter%p(3)),        &
                      (ComSTL%tri(i)%pTr(3)%p(3)-TriCenter%p(3)))

      MinObjX = dmin1((ComSTL%tri(i)%pTr(1)%p(1)-TriCenter%p(1)),        &
                      (ComSTL%tri(i)%pTr(2)%p(1)-TriCenter%p(1)),        &
                      (ComSTL%tri(i)%pTr(3)%p(1)-TriCenter%p(1)))
      MinObjY = dmin1((ComSTL%tri(i)%pTr(1)%p(2)-TriCenter%p(2)),        &
                      (ComSTL%tri(i)%pTr(2)%p(2)-TriCenter%p(2)),        &
                      (ComSTL%tri(i)%pTr(3)%p(2)-TriCenter%p(2)))
      MinObjZ = dmin1((ComSTL%tri(i)%pTr(1)%p(3)-TriCenter%p(3)),        &
                      (ComSTL%tri(i)%pTr(2)%p(3)-TriCenter%p(3)),        &
                      (ComSTL%tri(i)%pTr(3)%p(3)-TriCenter%p(3)))

      BoxSizePos(i,1) = max(4,int(MaxObjX/this%dx)+4)    !(Imax/2,1,1)
      BoxSizePos(i,2) = max(4,int(MaxObjY/this%dy)+4)    !(1,Jmax/2,1)
      BoxSizePos(i,3) = max(4,int(MaxObjZ/this%dz)+4)    !(1,1,Kmax/2)

      BoxSizeNeg(i,1) = min(-4,int(MinObjX/this%dx)-4)   !(Imax/2,1,1)
      BoxSizeNeg(i,2) = min(-4,int(MinObjY/this%dy)-4)   !(1,Jmax/2,1)
      BoxSizeNeg(i,3) = min(-4,int(MinObjZ/this%dz)-4)   !(1,1,Kmax/2)
    end do

    call cpu_time(finish)
    print*, 'cycle_ntri_time spent here_4 is: ', finish-start

    call ComputeGridIndexObject(this,ComSTL,GridIndex)

    call cpu_time(start)
    do i=1,ComSTL%ntri
      do ii=BoxSizeNeg(i,1),BoxSizePos(i,1)
        do jj=BoxSizeNeg(i,2),BoxSizePos(i,2)
          do kk=BoxSizeNeg(i,3),BoxSizePos(i,3)
            ! Set a point to the cell center
            i1 = max(1,GridIndex(i,1)+ii)
            j1 = max(1,GridIndex(i,2)+jj)
            k1 = max(1,GridIndex(i,3)+kk)
            i1 = min(i1,this%Imax)
            j1 = min(j1,this%Jmax)
            k1 = min(k1,this%Kmax)
            pt%p(1) = this%x(i1)
            pt%p(2) = this%y(j1)
            pt%p(3) = this%z(k1)
            ! Check whether the projection of a point lie inside the triangle.
            Isin    = CheckPointTriangle(ComSTL%tri(i),ComSTL%nt(i),pt)
            if(Isin.eqv..TRUE.) then
              ! Compute the distance from a point to a surface
              dSur  = -dot_product(ComSTL%tri(i)%pTr(1)%p,ComSTL%nt(i)%v)
              ! Surface = nx*x+ny*y+nz*z+dSur
              tempd = Distance(ComSTL%nt(i),dSur,pt)
              npoint%v = ComSTL%nt(i)%v
            else
              ! Check the projected point is inside the edge of triangle, the edge connect point p1 and p2
              IsinE1 = CheckPointLine(ComSTL%tri(i)%pTr(1),ComSTL%tri(i)%pTr(2),pt)
              ! The edge connects point 2 and 3
              IsinE2 = CheckPointLine(ComSTL%tri(i)%pTr(2),ComSTL%tri(i)%pTr(3),pt)
              ! The edge connects point 3 and 1
              IsinE3 = CheckPointLine(ComSTL%tri(i)%pTr(3),ComSTL%tri(i)%pTr(1),pt)
              ! Compute the distance from point to edges(d1,d2,d3)
              if(IsinE1.eqv..TRUE.) then
                d1 = Distance(ComSTL%tri(i)%pTr(1),ComSTL%tri(i)%pTr(2),pt)
              else
                d1 = maxLvs
              end if
              if(IsinE2.eqv..TRUE.) then
                d2 = Distance(ComSTL%tri(i)%pTr(2),ComSTL%tri(i)%pTr(3),pt)
              else
                d2 = maxLvs
              end if
              if(IsinE3.eqv..TRUE.) then
                d3 = Distance(ComSTL%tri(i)%pTr(3),ComSTL%tri(i)%pTr(1),pt)
              else
                d3 = maxLvs
              end if
              ! When the projected point is inside one of edges
              if((IsinE1.eqv..TRUE.).or.(IsinE2.eqv..TRUE.).or.(IsinE3.eqv..TRUE.)) then
                if(dabs(d1)<dabs(d2).and.dabs(d1)<dabs(d3).and.(IsinE1.eqv..TRUE.)) then
                  v1t%v = pt%p-ComSTL%tri(i)%pTr(1)%p
                  v12%v = ComSTL%tri(i)%pTr(2)%p-ComSTL%tri(i)%pTr(1)%p
                  ! Compute the projection point of P into line connecting p1 and p2
                  projpt%p = ComSTL%tri(i)%pTr(1)%p+                         &
                      dot_product(v1t%v,v12%v)/dot_product(v12%v,v12%v)*v12%v
                  ! Vector connecting the point and its projected point
                  npoint%v = pt%p-projpt%p
                !  tempd = dabs(d1)*dsign(1.d0,dot_product(ComSTL%nee(i,1)%v,npoint%v))
                end if
                if(dabs(d2)<dabs(d1).and.dabs(d2)<dabs(d3).and.(IsinE2.eqv..TRUE.)) then
                  v1t%v = pt%p-ComSTL%tri(i)%pTr(2)%p
                  v12%v = ComSTL%tri(i)%pTr(3)%p-ComSTL%tri(i)%pTr(2)%p
                  ! Compute the projection point of P into line connecting p2 and p3
                  projpt%p = ComSTL%tri(i)%pTr(2)%p+                         &
                      dot_product(v1t%v,v12%v)/dot_product(v12%v,v12%v)*v12%v
                  npoint%v = pt%p-projpt%p
                !  tempd = dabs(d2)*dsign(1.d0,dot_product(ComSTL%nee(i,2)%v,npoint%v))
                end if
                if(dabs(d3)<dabs(d1).and.dabs(d3)<dabs(d2).and.(IsinE3.eqv..TRUE.)) then
                  v1t%v = pt%p-ComSTL%tri(i)%pTr(3)%p
                  v12%v = ComSTL%tri(i)%pTr(1)%p-ComSTL%tri(i)%pTr(3)%p
                  ! Compute the projection point of P into line connecting p3 and p1
                  projpt%p = ComSTL%tri(i)%pTr(3)%p+                         &
                      dot_product(v1t%v,v12%v)/dot_product(v12%v,v12%v)*v12%v
                  npoint%v = pt%p-projpt%p
                !  tempd = dabs(d3)*dsign(1.d0,dot_product(ComSTL%nee(i,3)%v,npoint%v))
                end if
                tempd=min(d1,d2,d3)
              else
                tempd = maxLvs
              end if
              ! Compute the distance from point to triangels' points
              d1=Distance(pt,ComSTL%tri(i)%pTr(1))
              d2=Distance(pt,ComSTL%tri(i)%pTr(2))
              d3=Distance(pt,ComSTL%tri(i)%pTr(3))
              ! Compare distances and chose the smallest one
              if(dabs(d1)<dabs(d2).and.dabs(d1)<dabs(d3).and.dabs(d1)<dabs(tempd)) then
                npoint%v=pt%p-ComSTL%tri(i)%pTr(1)%p
              !  tempd=dabs(d1)*dsign(1.d0,dot_product(ComSTL%npt(i,1)%v,npoint%v))
              end if
              if(dabs(d2)<dabs(d1).and.dabs(d2)<dabs(d3).and.dabs(d2)<dabs(tempd)) then
                npoint%v=pt%p-ComSTL%tri(i)%pTr(2)%p
              !  tempd=dabs(d2)*dsign(1.d0,dot_product(ComSTL%npt(i,2)%v,npoint%v))
              end if
              if(dabs(d3)<dabs(d1).and.dabs(d3)<dabs(d2).and.dabs(d3)<dabs(tempd)) then
                npoint%v=pt%p-ComSTL%tri(i)%pTr(3)%p
              !  tempd=dabs(d3)*dsign(1.d0,dot_product(ComSTL%npt(i,3)%v,npoint%v))
              end if
              tempd=min(d1,d2,d3)
            end if

            if(dabs(tempd)<dabs(this%phi(i1,j1,k1))) then
                    this%phi(i1,j1,k1) = tempd
              call npoint%NormalizeVector
              this%nx(i1,j1,k1)  = npoint%v(1)
              this%ny(i1,j1,k1)  = npoint%v(2)
              this%nz(i1,j1,k1)  = npoint%v(3)
              this%vof(i1,j1,k1) = 1.d0
            end if
          end do
        end do
      end do
    end do
    call cpu_time(finish)
    print*, 'cycle_ntri_time spent here_5 is: ', finish-start

    call cpu_time(start)
    ! Using ray casting method to check whether the point inside or outside object
    do i = 1, this%Imax
      do j = 1, this%Jmax
        do k = 1, this%Kmax
          if(dabs(this%phi(i,j,k))<0.9d0*maxLvs) then
            pt%p(1) = this%x(i)
            pt%p(2) = this%y(j)
            pt%p(3) = this%z(k)
            Isin = RayCasting(pt,ComSTL)
            if(Isin.eqv..true.) then
              this%phi(i,j,k) = -dabs(this%phi(i,j,k))
            else
              this%phi(i,j,k) = dabs(this%phi(i,j,k))
            end if
          end if
        end do
      end do
    end do
    call cpu_time(finish)
    print*, 'ijk_raycasting_time spent here_6 is: ', finish-start

    deallocate(GridIndex)
    deallocate(BoxSizeNeg, BoxSizePos)

    call cpu_time(start)
    do i=1,this%Imax
      do j=1,this%Jmax
        do k=1,this%Kmax
          nx = this%nx(i,j,k)
          ny = this%ny(i,j,k)
          nz = this%nz(i,j,k)
          s=this%phi(i,j,k)+0.5*(dabs(nx)*this%dx+               &     !(i,j,k)
                                 dabs(ny)*this%dy+               &     !(i,j,k)
                                 dabs(nz)*this%dz)                     !(i,j,k)
          call Volume_Fraction_Calc(this%dx,this%dy,      &            !(i,j,k)
                                    this%dz,nx,ny,nz,s,vol)
          this%vof(i,j,k)=vol/(this%dx*this%dy*this%dz)
          ! per debug 13.05
          if(this%vof(i,j,k)<tolpar) this%vof(i,j,k)      = 0.d0
          if(this%vof(i,j,k)>1.d0-tolpar) this%vof(i,j,k) = 1.d0
          !if(isnan(this%vof(i,j,k))) then
          !  print*, 'Clsvof.f90 341'
          !  print*, this%phi(i,j,k)
          !  print*, vol
          !  print*, i,j,k
          !  print*, this%phi(i,j,k)
          !end if
          !write(5,*) this%vof(i,j,k)
        end do
      end do
    end do
    call cpu_time(finish)
    print*, 'ijk_vof_time spent here_7 is: ', finish-start

    call cpu_time(start)
    !per debug
    do i=1,this%Imax-1
      do j=1,this%Jmax
        do k=1,this%Kmax
        if(this%phi(i,j,k).lt.0.d0.and.this%phi(i+1,j,k).eq.1d4)then
                this%phi(i+1,j,k)=-1d4
                this%vof(i+1,j,k)=0.0_dp
        endif
        end do
      end do
    end do
    call cpu_time(finish)
    print*, 'ijk_vof_time spent here_8 is: ', finish-start

    close(5)

  END SUBROUTINE LvsObject


  LOGICAL FUNCTION RayCasting(pt,ComSTL) RESULT(Isin)
    !
    TYPE(TGpoint),   INTENT(in) :: pt
    TYPE(TsimcoSTL), INTENT(in) :: ComSTL
    !
    TYPE(Tvector)               :: nray
    INTEGER(KIND=it4b)          :: i, ctr
    LOGICAL                     :: Iscut
    INTEGER(KIND=it4b)          :: CutPoint, CutEdge, CutPlane, CutType

    ctr = 0
    call RANDOM_NUMBER(nray%v(1))
    call RANDOM_NUMBER(nray%v(2))
    call RANDOM_NUMBER(nray%v(3))
    call nray%NormalizeVector

    CutPoint = 0
    CutEdge = 0
    CutPlane = 0
    do i=1,ComSTL%ntri
      call RayCutTriangle(pt,nray,ComSTL%tri(i),ComSTL%nt(i), Iscut, CutType)
      if(Iscut.eqv..TRUE.) then
        ctr=ctr+1
      end if
      if(CutType == 0) CutPoint = CutPoint+1
      if(CutType == 1) CutEdge = CutEdge+1
      if(CutType == 2) CutPlane = CutPlane+1
    end do
    if(mod(ctr,2)==0) then
      Isin = .FALSE.
      return
    else
      Isin = .TRUE.
      return
    end if

  END FUNCTION RayCasting


  !! Determine the cell where the element triangle of STL locates at
  !
  PURE SUBROUTINE ComputeGridIndexObject(this, ComSTL, GridIndex)
    !
    CLASS(TSolidBody),               INTENT(inout) :: this
    TYPE(TsimcoSTL),                 INTENT(in)    :: ComSTL
    INTEGER(KIND=it4b), ALLOCATABLE, INTENT(inout) :: GridIndex(:,:)
    !
    INTEGER(KIND=it4b) :: i,j
    TYPE(TGpoint)      :: TriCenter

    do j=1,ComSTL%ntri
      call CenterPoint(ComSTL%tri(j),TriCenter)
      do i=1,this%Imax-1
        if(this%x(i)<=TriCenter%p(1).and.this%x(i+1)>TriCenter%p(1)) then
          GridIndex(j,1)=i
        end if
      end do
      do i=1,this%Jmax-1
        if(this%y(i)<=TriCenter%p(2).and.this%y(i+1)>TriCenter%p(2)) then
          GridIndex(j,2)=i
        end if
      end do
      do i=1,this%Kmax-1
        if(this%z(i)<=TriCenter%p(3).and.this%z(i+1)>TriCenter%p(3)) then
          GridIndex(j,3)=i
        end if
      end do
    end do
    !
  END SUBROUTINE ComputeGridIndexObject


  SUBROUTINE setInitialPosition( this, PGrid, UGrid, VGrid, WGrid, PCell, UCell, VCell, WCell, xc, yc, zc, p, q, r )
    !
    CLASS(TSolidBody), INTENT(inout) :: this
    TYPE(TGrid),       INTENT(in)    :: PGrid, UGrid, VGrid, WGrid
    TYPE(TCell),       INTENT(inout) :: PCell, UCell, VCell, WCell
    REAL(dp),          INTENT(in)    :: xc, yc, zc, p, q, r

    this%CMpoint (1) = xc
    this%CMpoint (2) = yc
    this%CMpoint (3) = zc
    this%CMorient(1) = p
    this%CMorient(2) = q
    this%CMorient(3) = r
    this%rotMatLtoG = this%calcRotationMatrix(this%CMorient)

    ! Interpolate to global grids
    CALL this%interpolateLvS( PGrid, PCell )
    CALL this%interpolateLvS( UGrid, UCell )
    CALL this%interpolateLvS( VGrid, VCell )
    CALL this%interpolateLvS( WGrid, WCell )

  END SUBROUTINE setInitialPosition


  !! Calculate the rotation matrix to express a vector defined in main grid coordinates
  !! in terms of local grid coordinates
  !
  FUNCTION calcRotationMatrix( this, angles ) RESULT( Mat )
    !
    CLASS(TSolidBody), INTENT(in)    :: this
    REAL(dp),          INTENT(in)    :: angles(3)
    REAL(dp)                         :: Mat(3,3)
    !
    REAL(dp) :: cosPsi, cosThe, cosPhi, sinPsi, sinThe, sinPhi

    cosPsi = COS(angles(1))
    cosThe = COS(angles(2))
    cosPhi = COS(angles(3))
    sinPsi = SIN(angles(1))
    sinThe = SIN(angles(2))
    sinPhi = SIN(angles(3))

    ! Wikipedia, Tait-Bryan angles, 
    ! Order: Roll-Pitch-Yaw from ship's to global referential (Phi-Theta-Psi, X-Y-Z)
    !        Yaw-Picth-Roll from global to ship's referential
    Mat(1,:) = (/cosPhi*cosThe, cosPhi*sinThe*sinPsi - cosPsi*sinPhi, sinPhi*sinPsi + cosPhi*cosPsi*sinThe/)
    Mat(2,:) = (/cosThe*sinPhi, cosPhi*cosPsi + sinPhi*sinThe*sinPsi, cosPsi*sinPhi*sinThe - cosPhi*sinPsi/)
    Mat(3,:) = (/-sinThe      , cosThe*sinPsi                       , cosThe*cosPsi                       /)


  END FUNCTION calcRotationMatrix


  !! Find coordinates in solid mesh of a point in the main mesh
  !
  PURE SUBROUTINE transferMeshGtoL( this, pointG, pointL )
    !
    CLASS(TSolidBody), INTENT(in)    :: this
    REAL(dp),          INTENT(in)    :: pointG(3)      ! Coordinates of the point in global reference frame
    REAL(dp),          INTENT(inout) :: pointL(3)      ! Coordinates of the point in local  reference frame

    ! Translate point to origin of local mesh, then rotate
    pointL = MATMUL(TRANSPOSE(this%rotMatLtoG), pointG - this%CMpoint)

  END SUBROUTINE transferMeshGtoL


  !! Find coordinates in solid mesh of a point in the main mesh
  !
  PURE SUBROUTINE transferMeshLtoG( this, pointL, pointG )
    !
    CLASS(TSolidBody), INTENT(in)    :: this
    REAL(dp),          INTENT(in)    :: pointL(3)      ! Coordinates of the point in local  reference frame
    REAL(dp),          INTENT(inout) :: pointG(3)      ! Coordinates of the point in global reference frame

    ! Rotate back to global mesh orientation, the translate point to origin of global mesh
    pointG = MATMUL(this%rotMatLtoG, pointL) + this%CMpoint

  END SUBROUTINE transferMeshLtoG


  !! Interpolate the level set function onto global grid
  !
  SUBROUTINE interpolateLvS( this, grid, cell )
    !
    CLASS(TSolidBody), INTENT(in)    :: this
    TYPE(TGrid),       INTENT(in)    :: grid    ! Global grid
    TYPE(TCell),       INTENT(inout) :: cell    ! Cell values on global grid
    !
    INTEGER  :: i, j, k, iL, jL, kL, ii, jj, kk, ImaxG, JmaxG, KmaxG
    REAL(dp) :: coeff, pointL(3), vol, s, temp

    ImaxG = SIZE(grid%x, 1)
    JmaxG = SIZE(grid%x, 2)
    KmaxG = SIZE(grid%x, 3)

    ! Remove any solid from grid
    cell%phi = 1.0e4_dp
    cell%vof = 1.0_dp

    DO i = 1, ImaxG
       DO j = 1, JmaxG
          DO k = 1, KmaxG
             ! Find coordinates of global cell center in local grid
             CALL this%transferMeshGtoL((/grid%x(i,j,k),grid%y(i,j,k),grid%z(i,j,k)/), pointL)

             ! Find a first node (i,j,k) in local grid in the cube of nodes around the point
             iL = CEILING((pointL(1) - this%x(1))/this%dx)
             jL = CEILING((pointL(2) - this%y(1))/this%dy)
             kL = CEILING((pointL(3) - this%z(1))/this%dz)

             ! If it ends up outside of the local mesh, there is no solid
             IF( iL < 1 .OR. iL > this%Imax-1 .OR. jL < 1 .OR. jL > this%Jmax-1 .OR. kL < 1 .OR. kL > this%Kmax-1 ) CYCLE

             ! Interpolate in local grid and store value in global cell i,j,k
             cell%phi(i,j,k) = 0.0_dp
             DO ii = 0, 1
                DO jj = 0, 1
                   DO kk = 0, 1
                      coeff = (this%dx-ABS(pointL(1)-this%x(iL+ii)))/this%dx *              &
                        &     (this%dy-ABS(pointL(2)-this%y(jL+jj)))/this%dy *              &
                        &     (this%dz-ABS(pointL(3)-this%z(kL+kk)))/this%dz
                      cell%phi(i,j,k) = cell%phi(i,j,k) + coeff * this%phi(iL+ii,jL+jj,kL+kk)
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    ! We need to have all interpolated phi to move on to calculate normal vectors
    ! After calculating normal vectors, we can recalculate more accurate vof. vof far from surface
    ! (either inside or outside of solid) should be ok from interpolation in previous loop.
    DO i = 1, ImaxG
       DO j = 1, JmaxG
          DO k = 1, KmaxG
             IF( ABS(cell%phi(i,j,k)) > 4.0_dp ) CYCLE

             ! Calculate normal vector from interpolated phi
             ! Note that vof is sent as argument, but actually unused
             IF(i>1.AND.i<ImaxG.AND.j>1.AND.j<jmaxG.AND.k>1.AND.k<kmaxG) THEN
                CALL Normal_Vector_Irre(grid, cell%vof, cell%phi, i, j, k,                 &
                                        cell%nx(i,j,k), cell%ny(i,j,k), cell%nz(i,j,k))
             ELSE
               ! For boundary cell                
               cell%nx(i,j,k)=(cell%phi(min(ImaxG,i+1),j,k)-cell%phi(max(1,i-1),j,k))/             &
                         (grid%x(min(ImaxG,i+1),j,k)-grid%x(max(1,i-1),j,k))
               cell%ny(i,j,k)=(cell%phi(i,min(JmaxG,j+1),k)-cell%phi(i,max(1,j-1),k))/             &
                         (grid%y(i,min(JmaxG,j+1),k)-grid%y(i,max(1,j-1),k))
               cell%nz(i,j,k)=(cell%phi(i,j,min(KmaxG,k+1))-cell%phi(i,j,max(1,k-1)))/             &
                         (grid%z(i,j,min(KmaxG,k+1))-grid%z(i,j,max(1,k-1)))
             END IF

             ! Recompute the normal vector such that summation of the square of all its component is 1
             temp = SQRT(cell%nx(i,j,k)**2+cell%ny(i,j,k)**2+cell%nz(i,j,k)**2)
             IF(temp<1.d-14) THEN
               cell%nx(i,j,k) = 0.0_dp
               cell%ny(i,j,k) = 0.0_dp
               cell%nz(i,j,k) = 1.0_dp
             ELSE
               cell%nx(i,j,k) = cell%nx(i,j,k)/temp
               cell%ny(i,j,k) = cell%ny(i,j,k)/temp
               cell%nz(i,j,k) = cell%nz(i,j,k)/temp
             END IF

             ! Recalculate accurate vof based on interpolated phi and normal
             s = cell%phi(i,j,k)+0.5_dp*(DABS(cell%nx(i,j,k))*grid%dx(i,j,k) + &
                                         DABS(cell%ny(i,j,k))*grid%dy(i,j,k) + &
                                         DABS(cell%nz(i,j,k))*grid%dz(i,j,k)   )

             CALL Volume_Fraction_Calc(grid%dx(i,j,k),grid%dy(i,j,k),         &
                  grid%dz(i,j,k),cell%nx(i,j,k),cell%ny(i,j,k),cell%nz(i,j,k),s,vol)

             cell%vof(i,j,k)=vol/(grid%dx(i,j,k)*grid%dy(i,j,k)*grid%dz(i,j,k))

          END DO
       END DO
    END DO

  END SUBROUTINE interpolateLvS


  !! Calculate one component of solid velocity on global grid
  !
  SUBROUTINE solidVelocity( this, pointG, vels, us, vs, ws )
    !
    CLASS(TSolidBody),  INTENT(in)    :: this
    REAL(dp),           INTENT(in)    :: pointG(3)  ! Coordinates of point where to calculate velocity
    REAL(dp), OPTIONAL, INTENT(out)   :: vels(3)    ! Solid velocity
    REAL(dp), OPTIONAL, INTENT(out)   :: us         ! Solid velocity - u component
    REAL(dp), OPTIONAL, INTENT(out)   :: vs         ! Solid velocity - v component
    REAL(dp), OPTIONAL, INTENT(out)   :: ws         ! Solid velocity - w component
    !
    REAL(dp) :: tmp(3), pointL0(3)

    ! Find coordinates of global point on local grid
    CALL this%transferMeshGtoL(pointG, tmp)

    ! Find coordinates of global point on local grid at the start of the time step
    ! Translate point to origin of local mesh, then rotate
    pointL0 = MATMUL(TRANSPOSE(this%rotMatLtoG0), pointG - this%CMpoint0)

    tmp = MATMUL(this%rotMatLtoG, (pointL0 - tmp)/this%dt)

    IF( PRESENT( vels ) ) THEN
       vels(1) = tmp(1)
       vels(2) = tmp(2)
       vels(3) = tmp(3)
    END IF
    
    IF( PRESENT( us ) ) us = tmp(1)
    IF( PRESENT( vs ) ) vs = tmp(2)
    IF( PRESENT( ws ) ) ws = tmp(3)

  END SUBROUTINE solidVelocity


  SUBROUTINE advance( this, PGrid, UGrid, VGrid, WGrid, PCell, UCell, VCell, WCell, Time )
    !
    CLASS(TSolidBody), INTENT(inout) :: this
    TYPE(TGrid),       INTENT(in)    :: PGrid, UGrid, VGrid, WGrid
    TYPE(TCell),       INTENT(inout) :: PCell, UCell, VCell, WCell
    TYPE(TSolverTime), INTENT(inout) :: Time

    this%dt          = Time%dt
    this%CMpoint0    = this%CMpoint
    this%CMorient0   = this%CMorient
    this%rotMatLtoG0 = this%rotMatLtoG

    ! Collect forces on the solid
    !TODO

    ! Move the solid
    !TODO
    ! Apply given movement for now
    this%CMpoint(1) = 0.01_dp * SIN(MAX(Time%PhysT-3600.0,0.0)/3600.0_dp*2.0_dp*pi)
    this%CMorient(3) = pi/6.0_dp * SIN(MAX(Time%PhysT-3600.0,0.0)/3600.0_dp*2.0_dp*pi)
    
    this%rotMatLtoG = this%calcRotationMatrix(this%CMorient)

    ! Interpolate to global grids
    CALL this%interpolateLvS( PGrid, PCell )
    CALL this%interpolateLvS( UGrid, UCell )
    CALL this%interpolateLvS( VGrid, VCell )
    CALL this%interpolateLvS( WGrid, WCell )

    CALL GridPreProcess(PGrid,UGrid,VGrid,WGrid, PCell,UCell,VCell,WCell)

    CALL DefineMomentumExchangeCell(PCell, UCell,VCell,WCell)

    CALL NumberExternalCell(0,0,0, PCell)
    CALL NumberExternalCell(1,0,0, UCell)
    CALL NumberExternalCell(0,1,0, VCell)
    CALL NumberExternalCell(0,0,1, WCell)

    CALL NewCellFace(PGrid,UGrid,VGrid,WGrid, PCell,UCell,VCell,WCell)

  END SUBROUTINE advance

END MODULE SolidBody

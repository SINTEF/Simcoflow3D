!! This module is design to read object information. The object reads
!!
module STL

  use PrecisionVar,only : sp, dp, it4b, it1b
  use Geometry,    only : GPoint, triangle, vector, compare, VecCpr

  implicit none
  private
  !! Define data for object
  type, public :: simcoSTL
    integer(kind=it4b)         :: ntri    ! Number of triangles
    type(triangle),allocatable :: tri(:)  ! Triangle position
    type(vector),allocatable   :: nt(:)   ! Normal vector of triangle
    type(vector),allocatable   :: npt(:,:) ! Normal vector of point
    type(vector),allocatable   :: nee(:,:) ! Normal vector of edge
    contains
      procedure,pass(this),public  :: ReadSTLFile
      procedure,pass(this),public  :: ModifySTLFile
      final                        :: destructor
  end type simcoSTL

  contains

  subroutine ReadSTLFile(filename,                                             &
                         this)
    !! The subroutine is used to read data from STL file
    class(simcoSTL),intent(inout)  :: this
    character(len=80),intent(in)   :: filename
    type(vector)                   :: e1, e2, n0
    integer*2                      :: padding
    integer(kind=it4b)             :: iunit,i,k,ctrp1,ctrp2,ctrp3
    character(len=80)              :: title
    real                           :: arr(3)
    real                           :: maxarr1, maxarr2, maxarr3
    real                              meanarr1, meanarr2, meanarr3
    real                              minarr1, minarr2, minarr3
    logical                        :: isExist,Edge1,Edge2,Edge3
    integer(kind=it4b)             :: stat,im,ierr 

    iunit=13
    inquire(file=filename,exist=isExist)
    if (.NOT.isExist) then
      print*, "File .STL does not exist, please check gain"
    end if

    open(unit=iunit,file=filename,status='old',form='unformatted',             &
                                                  access='stream')
    k=1
    read(iunit,POS=k) title
    k=k+sizeof(title)
    read(iunit,POS=k) this%ntri
    k=k+sizeof(this%ntri)
    allocate(this%nt(this%ntri))
    allocate(this%tri(this%ntri))
    allocate(this%npt(this%ntri,3))
    allocate(this%nee(this%ntri,3))
    meanarr1=0.d0
    meanarr2=0.d0
    meanarr3=0.d0
    maxarr1=-1.d6
    maxarr2=-1.d6
    maxarr3=-1.d6
    minarr1=1.d6
    minarr2=1.d6
    minarr3=1.d6
    do i=1,this%ntri
      ! Read normal vector
      call ReadArray(iunit,k,arr)
      call this%nt(i)%SetVector(dble(arr(1)),dble(arr(2)),dble(arr(3)))
      call this%nt(i)%NormalizeVector
      ! Read the position of triangle

      ! First point
      call ReadArray(iunit,k,arr)
      call this%tri(i)%pTr(1)%SetPoint(dble(arr(1)),dble(arr(2)),dble(arr(3)))
      meanarr1=meanarr1+arr(1)
      meanarr2=meanarr2+arr(2)
      meanarr3=meanarr3+arr(3)

      ! Second point
      call ReadArray(iunit,k,arr)
      call this%tri(i)%pTr(2)%SetPoint(dble(arr(1)),dble(arr(2)),dble(arr(3)))
      meanarr1=meanarr1+arr(1)
      meanarr2=meanarr2+arr(2)
      meanarr3=meanarr3+arr(3)

      ! Third point
      call ReadArray(iunit,k,arr)
      call this%tri(i)%pTr(3)%SetPoint(dble(arr(1)),dble(arr(2)),dble(arr(3)))

      read(iunit,POS=k) padding
      k=k+sizeof(padding)

      ! Compute the mean value of array
      meanarr1=meanarr1+arr(1)
      meanarr2=meanarr2+arr(2)
      meanarr3=meanarr3+arr(3)
      if(arr(1)>maxarr1) maxarr1 = arr(1)
      if(arr(2)>maxarr2) maxarr2 = arr(2)
      if(arr(3)>maxarr3) maxarr3 = arr(3)
      if(arr(1)<minarr1) minarr1 = arr(1)
      if(arr(2)<minarr2) minarr2 = arr(2)
      if(arr(3)<minarr3) minarr3 = arr(3)
    end do
    write(*,*) trim(filename),' has this title ',trim(title),' and has',       &
                                                 this%ntri,' triangles'
    close(iunit)
    meanarr1=meanarr1/this%ntri/3.d0
    meanarr2=meanarr2/this%ntri/3.d0
    meanarr3=meanarr3/this%ntri/3.d0
    ! Compute the normal vector at edges and vertices of triangles
    ! Recompute coordinate of points and normal vector
    ! The origin of coordinate is (0,0,0)
    ! The centre point of object is set at the origin
    do i=1,this%ntri
      this%tri(i)%PTr(:)%p(1)=(this%tri(i)%PTr(:)%p(1)-meanarr1)
      this%tri(i)%PTr(:)%p(2)=(this%tri(i)%PTr(:)%p(2)-meanarr2)
      this%tri(i)%PTr(:)%p(3)=(this%tri(i)%PTr(:)%p(3)-meanarr3)
    !  Compute vector connecting two points
      e1%v = this%tri(i)%pTr(3)%p-this%tri(i)%pTr(1)%p
      e2%v = this%tri(i)%pTr(3)%p-this%tri(i)%pTr(2)%p
      ! Cross product of the two vectors
      n0 = VecCpr(e1,e2)
      n0%v(:) = abs(n0%v(:))*dsign(1.d0,this%nt(i)%v(:))
      call n0%NormalizeVector
      this%nt(i)%v(:) = n0%v(:)
    end do  

    do i=1,this%ntri
      this%npt(i,1)%v=this%nt(i)%v
      this%npt(i,2)%v=this%nt(i)%v
      this%npt(i,3)%v=this%nt(i)%v

      ctrp1=1
      ctrp2=1
      ctrp3=1
      do k=1,this%ntri
        if(i/=k) then
      ! Compute the normal vector for shared points
          if((compare(this%tri(i)%pTr(1),this%tri(k)%pTr(1)).eqv..TRUE.).or.   &
             (compare(this%tri(i)%pTr(1),this%tri(k)%pTr(2)).eqv..TRUE.).or.   &
             (compare(this%tri(i)%pTr(1),this%tri(k)%pTr(3)).eqv..TRUE.)) then 
            ctrp1=ctrp1+1
      ! the formulation is based on the work of J.Andreas Barentzen and Henrik Aanes
      ! "Generating Signed Distance Fields from Triangle Meshes"  
            this%npt(i,1)%v=this%npt(i,1)%v+this%nt(k)%v
            Edge1=.TRUE.
          else
          Edge1 = .false.
          end if
          if((compare(this%tri(i)%pTr(2),this%tri(k)%pTr(1)).eqv..TRUE.).or.   &
             (compare(this%tri(i)%pTr(2),this%tri(k)%pTr(2)).eqv..TRUE.).or.   &
             (compare(this%tri(i)%pTr(2),this%tri(k)%pTr(3)).eqv..TRUE.)) then
            ctrp2=ctrp2+1
            this%npt(i,2)%v=this%npt(i,2)%v+this%nt(k)%v
            Edge2=.TRUE.
          else
            Edge2=.false.
          end if
          if((compare(this%tri(i)%pTr(3),this%tri(k)%pTr(1)).eqv..TRUE.).or.   &
             (compare(this%tri(i)%pTr(3),this%tri(k)%pTr(2)).eqv..TRUE.).or.   &
             (compare(this%tri(i)%pTr(3),this%tri(k)%pTr(3)).eqv..TRUE.)) then
            ctrp3=ctrp3+1
            this%npt(i,3)%v=this%npt(i,3)%v+this%nt(k)%v
            Edge3=.TRUE.
          else
            Edge3=.false.
          end if
          ! Compute the normal vector for triangle edges
          if((Edge1.eqv..TRUE.).and.(Edge2.eqv..TRUE.)) then
            this%nee(i,1)%v=0.5d0*(this%nt(i)%v+this%nt(k)%v)
            call this%nee(i,1)%NormalizeVector
          end if
          if((Edge2.eqv..TRUE.).and.(Edge3.eqv..TRUE.)) then
            this%nee(i,2)%v=0.5d0*(this%nt(i)%v+this%nt(k)%v)
            call this%nee(i,2)%NormalizeVector
          end if
          if((Edge3.eqv..TRUE.).and.(Edge1.eqv..TRUE.)) then
            this%nee(i,3)%v=0.5d0*(this%nt(i)%v+this%nt(k)%v)
            call this%nee(i,3)%NormalizeVector
          end if      
        end if
      end do
      if(ctrp1>0) this%npt(i,1)%v=this%npt(i,1)%v/dble(ctrp1)
      call this%npt(i,1)%NormalizeVector
      if(ctrp2>0) this%npt(i,2)%v=this%npt(i,2)%v/dble(ctrp2)
      call this%npt(i,2)%NormalizeVector
      if(ctrp3>0) this%npt(i,3)%v=this%npt(i,3)%v/dble(ctrp3)
      call this%npt(i,3)%NormalizeVector
    end do
  end subroutine ReadSTLFile

  subroutine ModifySTLFile(this, ScaleFactor, CentPos, Lref)
    !! The subroutine is used to scale the object and put it at the specific position
    class(simcoSTL)          , intent(inout)  :: this
    real(kind=dp)            , intent(in)     :: ScaleFactor 
    type(GPoint)             , intent(in)     :: CentPos
    real(kind=dp)  , optional, intent(in)     :: Lref
    integer(kind=it4b)  :: i
    ! Scale the object
    do i = 1, this%ntri
      this%tri(i)%PTr(1)%p(:) = this%tri(i)%PTr(1)%p(:)/ScaleFactor
      this%tri(i)%PTr(2)%p(:) = this%tri(i)%PTr(2)%p(:)/ScaleFactor
      this%tri(i)%PTr(3)%p(:) = this%tri(i)%PTr(3)%p(:)/ScaleFactor
    end do  
    ! Set new the origin of coordinate
    do i = 1, this%ntri
      this%tri(i)%PTr(1)%p(:) = this%tri(i)%PTr(1)%p(:) + CentPos%p(:)
      this%tri(i)%PTr(2)%p(:) = this%tri(i)%PTr(2)%p(:) + CentPos%p(:)
      this%tri(i)%PTr(3)%p(:) = this%tri(i)%PTr(3)%p(:) + CentPos%p(:)
      if(present(Lref)) then
        this%tri(i)%PTr(1)%p(:) = this%tri(i)%PTr(1)%p(:)/Lref
        this%tri(i)%PTr(2)%p(:) = this%tri(i)%PTr(2)%p(:)/Lref
        this%tri(i)%PTr(3)%p(:) = this%tri(i)%PTr(3)%p(:)/Lref
      end if    
    end do  
  end subroutine ModifySTLFile  

  subroutine destructor(this)
    type(simcoSTL),intent(inout) :: this

    if(allocated(this%tri))  deallocate(this%tri)
    if(allocated(this%nt))   deallocate(this%nt)
  end subroutine destructor

  subroutine ReadArray(iunit,k,arr)
    integer(kind=it4b),intent(in)    :: iunit
    integer(kind=it4b),intent(inout) :: k
    real,intent(out)                 :: arr(3)

    read(iunit,POS=k) arr(1),arr(2),arr(3)
    k=k+sizeof(arr)
  end subroutine ReadArray

end module STL

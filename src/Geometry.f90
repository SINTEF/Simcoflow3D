 !! Description:
 !! The module declare several entities which are used for geometrical operations.
 !
 !! Method:
 !! It involves the operation for point, edge, and triangles
 !
 ! Current Code Owner: SIMCOFlow
 !
 ! Code Description:
 ! Language: Fortran 90.
 ! Software Standards: "European Standards for Writing and
 ! Documenting Exchangeable Fortran 90 Code".
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! Author : Son Tung Dang
 !        : EPT, NTNU
 ! Date : 20.09.2019
module Geometry

  use PrecisionVar,only : dp,it4b
  implicit none
  private
  
  real(kind=dp), parameter :: tol=1.d-14

  type, public :: Gpoint
    real(kind=dp)   :: p(3)
    contains
      procedure,pass(this),public  :: SetPoint
  end type Gpoint

  type,public :: edge
    type(Gpoint)     :: pE(2)
    contains
      procedure,pass(this),public  :: SetEdge
  end type edge

  type,public :: triangle
    type(Gpoint)     :: pTr(3)
    contains
      procedure,pass(this),public  :: SetTriangle
  end type triangle

  type,public :: vector
    real(kind=dp)   :: v(3)
    contains
      procedure,pass(this),public :: SetVector
      procedure,pass(this),public :: NormalizeVector
  end type vector

  interface Distance
    module procedure DistancetoSurface, DistancetoLine, DistanceTwoPoints
  end interface Distance

  interface VecCPr
    module procedure VectorCrossProduct
  end interface VecCPr
  
  interface Compare
    module procedure ComparePoints,CompareEdges
  end interface Compare    

  public :: VecCPr, Compare, Distance, CenterPoint, CheckPointTriangle,        &
            CheckPointLine, ProjectionPoint, VolumeFractionCalculation,        &
            AbsValueVector, RayCutTriangle, PointInsideTriangle

  contains

  pure subroutine SetPoint(this,x,y,z)
    class(Gpoint),intent(inout)     :: this
    real(kind=dp),intent(in)       :: x,y,z

    this%p(1)=x
    this%p(2)=y
    this%p(3)=z
  end subroutine SetPoint

  pure subroutine SetEdge(this,p1,p2)
    class(Edge),intent(inout)      :: this
    type(Gpoint),intent(in)         :: p1,p2

    this%pE(1)=p1
    this%pE(2)=p2
  end subroutine SetEdge

  pure subroutine SetTriangle(this,p1,p2,p3)
    class(Triangle),intent(inout)  :: this
    type(Gpoint),intent(in)         :: p1,p2,p3

    this%pTr(1)=p1
    this%pTr(2)=p2
    this%pTr(3)=p3
  end subroutine SetTriangle

  pure subroutine SetVector(this,x,y,z)
    class(Vector),intent(inout)    :: this
    real(kind=dp),intent(in)       :: x,y,z

    this%v(1)=x
    this%v(2)=y
    this%v(3)=z
  end subroutine SetVector

  pure subroutine CenterPoint(this,Cpoint)
    !! The subroutine is designed to compute the mass centre point of triangle
    type(triangle),intent(in)  :: this
    type(Gpoint),intent(out)   :: Cpoint

    Cpoint%p(1)=1.d0/3.d0*(this%pTr(1)%p(1)+this%pTr(2)%p(1)+this%pTr(3)%p(1))
    Cpoint%p(2)=1.d0/3.d0*(this%pTr(1)%p(2)+this%pTr(2)%p(2)+this%pTr(3)%p(2))
    Cpoint%p(3)=1.d0/3.d0*(this%pTr(1)%p(3)+this%pTr(2)%p(3)+this%pTr(3)%p(3))
  end subroutine CenterPoint

  pure logical function ComparePoints(p1,p2) result(Iseq)
    !! The function is designed to compare two points. 
    !! True if they are identical and false if not. 
    type(Gpoint),intent(in)    :: p1,p2
    real(kind=dp)              :: p1p2

    p1p2 = dsqrt((p1%p(1)-p2%p(1))**2.d0+                                      &
                 (p1%p(2)-p2%p(2))**2.d0+                                      &
                 (p1%p(3)-p2%p(3))**2.d0)
    if(p1p2<tol) then
      Iseq = .TRUE.
    else
      Iseq = .FALSE.
    end if    
  end function ComparePoints  

  pure logical function CompareEdges(pa1,pa2,pb1,pb2) result(Iseq)
    !! This function is designed to compare two edges.
    !! True if they are identical and false if not.
    type(Gpoint),intent(in)    :: pa1, pa2, pb1, pb2
    logical                    :: ab1, ab2
    if(compare(pa1,pb1).eqv..true..or.compare(pa1,pb2).eqv..true.) then
      ab1=.true. 
    end if
    if(compare(pa2,pb1).eqv..true..or.compare(pa2,pb2).eqv..true.) then
      ab2=.true.
    end if
    
    if(ab1.eqv..true..and.ab2.eqv..true.) then
      Iseq = .TRUE.
    else
      Iseq = .FALSE.
    end if        
  end function CompareEdges
  
  pure logical function CheckPointTriangle(tri,nt,pt) result(Isin)
    !! The subroutine is used to check whether the projection of a point lie inside
    !! the triangle.
    !! There is an elegant solution to this given by W. Heidrich,
    !! Journal of Graphics, GPU, and Game Tools,Volume 10, Issue 3, 2005.
    !! The formulation is given in :
    !! https://math.stackexchange.com/questions/544946/
    type(triangle),intent(in) :: tri
    type(vector),intent(in)   :: nt
    type(Gpoint),intent(in)   :: pt
    real(kind=dp)             :: a,b,c
    type(vector)              :: u,v,w,uw,wv,n

    u%v=tri%pTr(2)%p-tri%pTr(1)%p
    v%v=tri%pTr(3)%p-tri%pTr(1)%p
    w%v=pt%p-tri%pTr(1)%p
    uw=VecCPr(u,w)
    wv=VecCpr(w,v)
    n=VecCPr(u,v)
    a=dot_product(uw%v,n%v)/(n%v(1)**2.d0+n%v(2)**2.d0+n%v(3)**2.d0)
    b=dot_product(wv%v,n%v)/(n%v(1)**2.d0+n%v(2)**2.d0+n%v(3)**2.d0)
    c=1.d0-a-b
    if((a>=0.d0.and.a<=1.d0).and.                                              &
       (b>=0.d0.and.b<=1.d0).and.                                              &
       (c>=0.d0.and.c<=1.d0)) then
      Isin=.TRUE.
      return
    else
      Isin=.FALSE.
      return
    end if
  end function CheckPointTriangle
  
  logical function CheckPointLine(p1,p2,pt) result(Isin)
    !! The subroutine is used to check whether projection of a point lie inside 
    !! The line segment
    type(Gpoint),intent(in) :: p1,p2,pt
    type(Gpoint)            :: projpt
    type(vector)            :: v1t,v12
    real(kind=dp)           :: ang

    v1t%v = pt%p-p1%p
    v12%v = p2%p-p1%p
    ! Compute the projection point of P into line connecting p1 and p2
    ! P'(P) on AB = A+AP*AB/(AB*AB)*AB
    projpt%p = p1%p+dot_product(v1t%v,v12%v)/dot_product(v12%v,v12%v)*v12%v
    v1t%v = p1%p-projpt%p
    v12%v = p2%p-projpt%p
    ang = AngleVectors(v1t,v12)
    if(ang>0.d0) then
      Isin = .FALSE.
    else
      Isin = .TRUE.
    end if      
  end function CheckPointLine
  
  subroutine RayCutTriangle(pt, nray, tri, ntri, Iscut, CutType) 
    type(Gpoint),   intent(in) :: pt
    type(vector),   intent(in) :: nray
    type(triangle), intent(in) :: tri
    type(vector),   intent(in) :: ntri
    logical,        intent(out) :: Iscut
    integer,        intent(out) :: CutType
    type(vector)               :: vptr
    type(Gpoint)               :: IPtr
    real(kind=dp)              :: a, b, r
    logical                    :: Isin

    CutType = -1
    ! CutType = 0 : Ray cut triangle at corner's points
    ! CutType = 1 : Ray cut triangle at edges
    ! CutType = 2 : Ray is parallel to triangle
    vptr%v = pt%p-tri%pTr(1)%p
    a = -dot_product(ntri%v,vptr%v)
    b = dot_product(ntri%v,nray%v)
    
    if(dabs(b)<tol) then
      if(a==0) then
        ! The ray is inside the triangle
        Iscut = .TRUE.
        CutType = 2
        return
      else
        ! The ray disjoint from triangle
        Iscut = .FALSE.  
        ! print*, 'The ray disjoint from triangle'
        return
      end if
    end if  
    r = a/b
    ! ray goes away from triangle
    if(r<0.d0) then
      Iscut = .FALSE.
      ! print*, 'The ray goes away from triangle'
      return     
    end if  
    IPtr%p=pt%p+r*nray%v
    call PointInsideTriangle(IPtr, tri, Isin, CutType)
    if(Isin.eqv..TRUE.) then
      Iscut = .TRUE.
      return
    else
      Iscut = .FALSE.
      return
    end if    
  end subroutine RayCutTriangle

  ! pure logical function PointInsideTriangle(pt,tri) result(Isin)
  !   !! The function is used to check whether the point is inside the triangle
  !   type(Gpoint),   intent(in) :: pt
  !   type(triangle), intent(in) :: tri
  !   real(kind=dp)              :: D,s,t,uu,uv,vv,wu,wv
  !   type(vector)               :: u,v,w

  !   u%v=tri%pTr(2)%p-tri%pTr(1)%p
  !   v%v=tri%pTr(3)%p-tri%pTr(1)%p
  !   w%v=pt%p-tri%pTr(1)%p
  !   uu = dot_product(u%v,u%v)
  !   uv = dot_product(u%v,v%v)
  !   vv = dot_product(v%v,v%v)
  !   wu = dot_product(w%v,u%v)
  !   wv = dot_product(w%v,v%v)
  !   D = uv*uv-uu*vv
  !   s = (uv*wv-vv*wu)/dmax1(D,tol)
  !   if(s<0.d0.or.s>1.d0) then
  !     Isin = .TRUE.
  !     return
  !   end if  
  !   t = (uv*wu-uu*wv)/dmax1(D,tol)
  !   if(t<0.d0.or.(s+t)>0.d0) then
  !     Isin = .TRUE.
  !     return
  !   end if
  !   Isin = .FALSE.  
  ! end function PointInsideTriangle 
  
  pure subroutine PointInsideTriangle(pt, tri, Isin, CutType)
    !! The function is used to check whether the point is inside the triangle
    !! It based on https://math.stackexchange.com/questions/4322/check-whether-a-point-is-within-a-3d-triangle
    !! The barycentric coordinates
    type(Gpoint),   intent(in)  :: pt
    type(triangle), intent(in)  :: tri
    logical,        intent(out) :: Isin
    integer,        intent(out) :: CutType

    real(kind=dp)               :: D,s,t,Area,Area1,Area2,Area3, alpha, beta, gamma
    type(vector)                :: u,v,w,uv,uw,vw
     
    ! Compute vector representing two edges of triangle 
    u%v = tri%pTr(1)%p-tri%pTr(2)%p
    v%v = tri%pTr(1)%p-tri%pTr(3)%p
    uv = VecCPr(u,v)
    ! Compute the Area of triangle
    Area = AbsValueVector(uv)/2.d0
    
    ! Compute the vector connecting points to triangle's corners
    u%v=tri%pTr(1)%p-pt%p
    v%v=tri%pTr(2)%p-pt%p
    w%v=tri%pTr(3)%p-pt%p
    
    vw = VecCPr(v,w)
    uv = VecCPr(u,v)
    uw = VecCPr(u,w)
    ! Compute the area of triangle connecting point p and other triangle's cornerss
    Area1 = AbsValueVector(vw)/2.d0
    Area2 = AbsValueVector(uv)/2.d0
    Area3 = AbsValueVector(uw)/2.d0
    gamma = Area1/Area
    alpha = Area2/Area
    beta = Area3/Area
    
    ! The condition
    if((alpha>=0.d0.and.alpha<=1.d0).and.(beta>=0.d0.and.beta<=1.d0).and.      &
       (gamma>=0.d0.and.gamma<=1.d0).and.dabs(alpha+beta+gamma-1.d0)<tol) then
      Isin = .TRUE.
      ! Cut the triangle at corner points
      if(dabs(alpha-1.d0)<tol.or.dabs(beta-1.d0)<tol.or.dabs(gamma-1.d0)<tol) then
        CutType = 0
        return
      end if  
      ! Cut the triangle at edges
      if(dabs(alpha)<tol.or.dabs(beta)<tol.or.dabs(gamma)<tol) then
        CutType = 1
        return
      end if  
    else
      Isin = .FALSE.
    end if
  end subroutine PointInsideTriangle

  pure subroutine NormalizeVector(this)
    class(vector),intent(inout) :: this
    real(kind=dp)               :: sq
  
    sq=dsqrt(this%v(1)**2.d0+this%v(2)**2.d0+this%v(3)**2.d0)
    this%v(1)=this%v(1)/sq
    this%v(2)=this%v(2)/sq
    this%v(3)=this%v(3)/sq
    ! Set normal vector when all components are close to 0
    if(dabs(this%v(1))<tol.and.dabs(this%v(2))<tol.and.dabs(this%v(3))<tol) then
      this%v(1) = 0.d0
      this%v(2) = 0.d0
      this%v(3) = 1.d0
      return
    end if  
    ! Renormalize vector if all components are too small
    sq = dsqrt(this%v(1)**2.d0+this%v(2)**2.d0+this%v(3)**2.d0)
    if(dabs(sq-1.d0)>tol) then
      this%v(1) = this%v(1)/sq
      this%v(2) = this%v(2)/sq
      this%v(3) = this%v(3)/sq
    end if  

  end subroutine NormalizeVector

  pure type(vector) function VectorCrossProduct(a,b) result(c)
    type(vector),intent(in)   :: a,b

    c%v(1)=(a%v(2)*b%v(3)-a%v(3)*b%v(2))
    c%v(2)=(a%v(3)*b%v(1)-a%v(1)*b%v(3))
    c%v(3)=(a%v(1)*b%v(2)-a%v(2)*b%v(1))
  end function VectorCrossProduct

  pure real(kind=dp) function AbsValueVector(a) result(abs)
    type(vector),intent(in)   :: a

    abs = dsqrt(a%v(1)**2.d0+a%v(2)**2.d0+a%v(3)**2.d0)
  end function AbsValueVector

  pure real(kind=dp) function AngleVectors(a,b) result(ang)
    type(vector),intent(in)   :: a,b

    ang = dot_product(a%v,b%v)/dmax1(AbsValueVector(a),tol)/                   &
                               dmax1(AbsValueVector(b),tol)
  end function AngleVectors
    
  pure real(kind=dp) function DistancetoSurface(nSur,dSur,pTar) result(d)
    type(vector),intent(in)   :: nSur
    type(Gpoint),intent(in)   :: pTar
    real(kind=dp),intent(in)  :: dSur

    d=dabs(pTar%p(1)*nSur%v(1)+pTar%p(2)*nSur%v(2)+pTar%p(3)*nSur%v(3)+dSur)/dsqrt(nSur%v(1)**2.d0+nSur%v(2)**2.d0+nSur%v(3)**2.d0)
  end function DistancetoSurface

  pure real(kind=dp) function DistanceTwoPoints(p1,p2) result(d)
    type(Gpoint),intent(in)    :: p1,p2

    d=dsqrt((p2%p(1)-p1%p(1))**2.d0+(p2%p(2)-p1%p(2))**2.d0+                   &
            (p2%p(3)-p1%p(3))**2.d0)
  end function DistanceTwoPoints
  
  pure real(kind=dp) function DistancetoLine(p1,p2,pTar) result(d)
    !! The function is used to compute the distance from a point to a line 
    !! connecting point p1 and p2
    type(Gpoint),intent(in)    :: p1,p2,pTar
    type(vector)               :: v1t,v2t,v12,v12t 
    
    v1t%v = PTar%p - p1%p
    v2t%v = pTar%p - p2%p
    v12%v = p2%p - p1%p
    v12t = VecCPr(v1t,v2t)
    d = AbsValueVector(v12t)/dmax1(AbsValueVector(v12),tol)
  end function DistancetoLine 

  pure subroutine ProjectionPoint(nSur,dSur,PTar,pt)
    !! Method to compute projection point based on its location and distance from
    !! a point to the plane.
    type(vector),intent(in)   :: nSur
    type(Gpoint),intent(in)    :: pTar
    real(kind=dp),intent(in)  :: dSur
    type(Gpoint),intent(out)   :: pt
    real(kind=dp)             :: d

    d=DistancetoSurface(nSur,dSur,pTar)
    pt%p(:)=PTar%p(:)+d*nSur%v(:)
  end subroutine ProjectionPoint
  
  subroutine VolumeFractionCalculation(delx,dely,delz,n,s,vol)
    !! The subroutine is used to compute the fluid (solid) volume fraction inside the box.
    !! The interface is given by normal vector (n) and the distance (s) from the box centre to the interface
    implicit none
    real(kind=dp),intent(in)  :: delx,dely,delz
    type(vector),intent(in)   :: n
    real(kind=dp),intent(in)  :: s
    real(kind=dp),intent(out) :: vol
    real(kind=dp)             :: nx1,ny1,nz1
    real(kind=dp)	      :: delx1,dely1,delz1
    real(kind=dp)	      :: sc,sm,fc
    
    nx1 = dabs(n%v(1))
    ny1 = dabs(n%v(2))
    nz1 = dabs(n%v(3))
    
    ! Change the coordinate system. 
    delx1 = dmax1(nx1*delx,ny1*dely,nz1*delz)
    delz1 = dmin1(nx1*delx,ny1*dely,nz1*delz)
    dely1 = (nx1*delx+ny1*dely+nz1*delz)-delx1-delz1
    sm    = delx1+dely1+delz1
    sc    = dmin1(s,sm-s)
    if(s<=0.d0) then
      vol =  0.d0
      return
    end if
    if(s>=sm) then     ! be careful with this condition. it affects to region containing fluid
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
      vol = fc*delx*dely*delz            ! region containing fluid
      return
    else
      vol = (1.d0-fc)*delx*dely*delz
      return
    end if
    if(isnan(vol)) then
      pause 'Nan in Volume Fraction Calculation'
    end if
    if(fc>1.d10) then
      pause 'Infinity in Volume Fraction Calculation'
    end if
  end subroutine VolumeFractionCalculation
end module geometry

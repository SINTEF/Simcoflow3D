Module Matrix
    use precisionvar
    implicit none
    public

    interface SetArraytoZero
      module procedure SetArraytoZero1DReal,SetArraytoZero2DReal,SetArraytoZero3DReal
    end interface
    contains

    subroutine inverse(a,c,n)
!***********************************************************
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer n
real(kind=dp):: a(n,n), c(n,n)
real(kind=dp):: L(n,n), U(n,n), b(n), d(n), x(n)
real(kind=dp):: coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.d0
U=0.d0
b=0.d0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.d0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.d0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.d0
end do
end subroutine inverse

subroutine Matrixinv(a,b,n)
 ! subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
 ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
  integer :: i,j,k,l,m,n,irow
  double precision:: big,a(n,n),b(n,n),dum

 !build the identity matrix
  do i = 1,n
    do j = 1,n
      b(i,j) = 0.d0
    end do
    b(i,i) = 1.d0
  end do

  do i = 1,n ! this is the big loop over all the columns of a(n,n)
 ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
 ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
    big = a(i,i)
    do j = i,n
      if(dabs(a(j,i)).gt.big)then
        big = dabs(a(j,i))
        irow = j
      end if
    end do
 ! interchange lines i with irow for both a() and b() matrices
    if (big.gt.a(i,i)) then
      do k = 1,n
        dum = a(i,k) ! matrix a()
        a(i,k) = a(irow,k)
        a(irow,k) = dum
        dum = b(i,k) ! matrix b()
        b(i,k) = b(irow,k)
        b(irow,k) = dum
      end do
    end if
 ! divide all entries in line i from a(i,j) by the value a(i,i);
 ! same operation for the identity matrix
    dum = a(i,i)
    if(dum==0.d0) then
      pause 'matrixinv 126'
    end if
    do j = 1,n
      a(i,j) = a(i,j)/dum
      b(i,j) = b(i,j)/dum
    end do
 ! make zero all entries in the column a(j,i); same operation for indent()
    do j = i+1,n
      dum = a(j,i)
      do k = 1,n
        a(j,k) = a(j,k) - dum*a(i,k)
        b(j,k) = b(j,k) - dum*b(i,k)
      end do
    end do
  end do

 ! substract appropiate multiple of row j from row j-1
  do i = 1,n-1
    do j = i+1,n
      dum = a(i,j)
      do l = 1,n
        a(i,l) = a(i,l)-dum*a(j,l)
        b(i,l) = b(i,l)-dum*b(j,l)
      end do
    end do
  end do
end

real(kind=dp) function LangrangePolynomial(Pos,Var,PosTar) result(VarTar)
  !! Function for computing value based on 1D Langrange polynomial interpolation
  implicit none
  real(kind=dp),dimension(:),intent(in) :: Pos
  !! The position of variables
  real(kind=dp),dimension(:),intent(in) :: Var
  !! The value of variables
  real(kind=dp),intent(in)              :: PosTar
  !! The position of interpolated variable
  integer(kind=it4b)                    :: i,j,sizeArr
  real(kind=dp)                         :: Mult
  sizeArr = size(Pos)
  VarTar = 0.d0
  do i=1,sizeArr
    Mult = Var(i)
    do j=1,sizeArr
      if(j/=i) then
        Mult=Mult*(PosTar-Pos(j))/(Pos(i)-Pos(j))
      end if
    end do
    VarTar = VarTar+Mult
  end do      
end function LangrangePolynomial  

subroutine SetArraytoZero3DReal(Arr)
  implicit none
  real(kind=dp),dimension(:,:,:),intent(inout) :: Arr
  integer(kind=it4b)                           :: i,j,k,ArrD1,ArrD2,ArrD3
  ArrD1=size(Arr,1)
  ArrD2=size(Arr,2)
  ArrD3=size(Arr,3)
  do i=1,ArrD1
    do j=1,ArrD2
      do k=1,ArrD3
        Arr(i,j,k)=0.d0
      end do
    end do
  end do   
end subroutine SetArraytoZero3DReal

subroutine SetArraytoZero2DReal(Arr)
  implicit none
  real(kind=dp),dimension(:,:),intent(inout) :: Arr
  integer(kind=it4b)                         :: i,j,ArrD1,ArrD2
  ArrD1=size(Arr,1)
  ArrD2=size(Arr,2)
  do i=1,ArrD1
    do j=1,ArrD2
      Arr(i,j)=0.d0
    end do
  end do   
end subroutine SetArraytoZero2DReal  

subroutine SetArraytoZero1DReal(Arr)
  implicit none
  real(kind=dp),dimension(:),intent(inout) :: Arr
  integer(kind=it4b)                         :: i,ArrD1
  ArrD1=size(Arr,1)
  do i=1,ArrD1
    Arr(i)=0.d0
  end do   
end subroutine SetArraytoZero1DReal    
End module Matrix

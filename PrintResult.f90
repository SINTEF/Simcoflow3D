Module PrintResult
    USE PrecisionVar
    USE Mesh
    USE StateVariables,ONLY : ight,jght,kght,variables,pi,Ktw,rey,xc,yc
    USE Clsvof
    USE Cutcell
    Implicit none
    Private
    character(len=1), parameter :: newline=achar(10)
    Public:: PrintResultTecplotPCent,PrintResultTecplotUCent,                  &
             PrintResultTecplotVCent,PrintResultTecplotWCent,                  &
             PrintResultTecplotPCentXY
    Interface PrintResultTecplotPCent
       Module procedure PrintResultTecplotPCent
    End Interface PrintResultTecplotPCent
    Interface PrintResultTecplotUCent
       Module procedure PrintResultTecplotUCent
    End Interface PrintResultTecplotUCent
    Interface PrintResultTecplotVCent
       Module procedure PrintResultTecplotVCent
    End Interface PrintResultTecplotVCent
    Interface PrintResultTecplotWCent
       Module procedure PrintResultTecplotWCent
    End Interface PrintResultTecplotWCent
    Interface PrintResultTecplotPCentXY
       Module procedure PrintResultTecplotPCentXY
    End Interface
    Contains
    Subroutine PrintResultTecplotPCent(TGrid,TVar,TCell,iter)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(inout):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: iter
      Integer(kind=it4b) i,j,k
      Real(kind=dp),dimension(:,:,:),allocatable:: p
      Character(15) curd
      Allocate(p(Imax,Jmax,Kmax))
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            p(i,j,k) = TVar%p(i,j,k)
          End do
        End do
      End do
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) iter
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'Pressure.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","z","p","u","v","w","phi","vof","Pres","Mres"'
      Write(5,1110) Imax,Jmax,Kmax
      Write(5,"(f18.10)")(((TGrid%x(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%y(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%z(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((p(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      write(5,"(f18.10)")(((0.5d0*(TVar%u(i-1,j,k)+TVar%u(i,j,k)),i=1,Imax),   &
                                                            j=1,Jmax),k=1,Kmax)
      write(5,"(f18.10)")(((0.5d0*(TVar%v(i,j-1,k)+TVar%v(i,j,k)),i=1,Imax),   &
                                                            j=1,Jmax),k=1,Kmax)
      write(5,"(f18.10)")(((0.5d0*(TVar%w(i,j,k-1)+TVar%w(i,j,k)),i=1,Imax),   &
                                                            j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%phi(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vof(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%pres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%mres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
      Deallocate(p)
    End Subroutine PrintResultTecplotPCent

    Subroutine PrintResultTecplotPCentXY(TGrid,TVar,TCell,iter)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(inout):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: iter
      Integer(kind=it4b) i,j,k,kp
      Real(kind=dp),dimension(:,:,:),allocatable:: p
      Character(15) curd
      Allocate(p(Imax,Jmax,Kmax))
      kp = int(Kmax/2)
      kp = kp
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            p(i,j,k) = TVar%p(i,j,k)
          End do
        End do
      End do
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) iter
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'PresXY.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","p","u","v","phi","vof","Pres","Mres"'
      Write(5,1110) Imax,Jmax
      Write(5,"(f18.10)")((TGrid%x(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((TGrid%y(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((p(i,j,kp),i=1,Imax),j=1,Jmax)
      write(5,"(f18.10)")((0.5d0*(TVar%u(i-1,j,kp)+TVar%u(i,j,kp)),i=1,Imax),  &
                                                            j=1,Jmax)
      write(5,"(f18.10)")((0.5d0*(TVar%v(i,j-1,kp)+TVar%v(i,j,kp)),i=1,Imax),  &
                                                            j=1,Jmax)
      Write(5,"(f18.10)")((TCell%phi(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((TCell%vof(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((TVar%pres(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((TVar%mres(i,j,kp),i=1,Imax),j=1,Jmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,'f=block')
      Deallocate(p)
    End Subroutine PrintResultTecplotPCentXY

    Subroutine PrintResultTecplotUCent(TGrid,TVar,TCell,itt)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(in):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: itt
      Integer(kind=it4b) i,j,k
      Character(15) curd
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) itt
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'Uvelocity.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","z","u","Gp","phi","vof","ures"'
      Write(5,1110) Imax,Jmax,Kmax
      Write(5,"(f18.10)")(((TGrid%x(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%y(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%z(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%u(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%Gpu(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%phi(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vof(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%ures(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
    End Subroutine PrintResultTecplotUCent

    Subroutine PrintResultTecplotVCent(TGrid,TVar,TCell,itt)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(in):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: itt
      Integer(kind=it4b) i,j,k
      Character(15) curd
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) itt
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'Vvelocity.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","z","v","Gp","phi","vof","vres"'
      Write(5,1110) Imax,Jmax,Kmax
      Write(5,"(f18.10)")(((TGrid%x(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%y(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%z(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%v(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%Gpv(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%phi(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vof(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%vres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
    End Subroutine PrintResultTecplotVCent

    Subroutine PrintResultTecplotWCent(TGrid,TVar,TCell,itt)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(in):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: itt
      Integer(kind=it4b) i,j,k
      Character(15) curd
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) itt
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'Wvelocity.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","z","w","Gp","phi","vof","wres"'
      Write(5,1110) Imax,Jmax,Kmax
      Write(5,"(f18.10)")(((TGrid%x(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%y(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%z(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%w(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%Gpw(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%phi(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vof(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%wres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
    End Subroutine PrintResultTecplotWCent
End module PrintResult
